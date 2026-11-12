"""Batch replacement of parcel boundaries on the parcel ID map.

The workflow is:

1. Validate the base parcel map (single-band, north-up, UInt16/UInt32,
   NoData = 0) and the WGS84 GeoJSON payload (one Polygon per existing
   parcel ID, at most 32 features).
2. Reproject every polygon to the raster CRS and express it in raster
   pixel coordinates.
3. Dry pass over the raster in its own (bounded) blocks: rasterize each
   polygon to find candidate pixels, then keep exactly the pixels whose
   *center* lies inside the polygon or on its boundary.  While doing so,
   collect every conflict and every polygon that hits no pixel.
4. Only when the whole batch is clean, write a second, block-wise pass to
   a new GeoTIFF that preserves the original size, CRS, grid and dtype.

Nothing is written when any validation rule fails.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.warp import transform_geom
from rasterio.windows import Window

from app.analysis import _split_window

_DEFAULT_BLOCK_SIZE = 1024
# Tolerance, in pixels, for deciding that a pixel center sits on an edge.
# Loose enough to absorb round-off from the WGS84 round trip, yet far below
# one pixel (0.1 mm at a 1 m resolution, ~1 mm at 10 m) so nearby pixels are
# never pulled in by mistake.
_EDGE_EPS = 1e-4
_MAX_PARCELS = 32
_IDENTITY = Affine.identity()


class BoundaryError(Exception):
    """Rejection of a boundary-update request, with structured details."""

    def __init__(
        self,
        message: str,
        *,
        geojson_errors: list[str] | None = None,
        unknown_parcel_ids: list[int] | None = None,
        no_hit_parcel_ids: list[int] | None = None,
        conflicts: list[dict[str, Any]] | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.geojson_errors = geojson_errors or []
        self.unknown_parcel_ids = unknown_parcel_ids or []
        self.no_hit_parcel_ids = no_hit_parcel_ids or []
        self.conflicts = conflicts or []

    def as_detail(self) -> dict[str, Any]:
        return {
            "message": self.message,
            "geojson_errors": self.geojson_errors,
            "unknown_parcel_ids": self.unknown_parcel_ids,
            "no_hit_parcel_ids": self.no_hit_parcel_ids,
            "conflicts": self.conflicts,
        }


@dataclass(frozen=True)
class RasterSpec:
    path: str
    width: int
    height: int
    crs: CRS
    transform: Affine
    dtype: str
    nodata: float | int | None
    id_max: int


@dataclass(frozen=True)
class BoundaryUpdate:
    """One corrected boundary, already mapped into raster pixel space."""

    parcel_id: int
    # List of rings; ring 0 is the exterior, the rest are holes.  Each ring
    # is an (N, 2) array of (column, row) pixel coordinates, closed.
    pixel_rings: list[np.ndarray]


@dataclass
class _ScanFindings:
    hit_counts: dict[int, int] = field(default_factory=dict)
    existing_ids: set[int] = field(default_factory=set)
    # Keyed by (conflict_type, sorted id pair) so each pair is reported once.
    conflicts: dict[tuple[str, tuple[int, int]], dict[str, Any]] = field(
        default_factory=dict
    )

    def remember_conflict(
        self, conflict_type: str, first: int, second: int, row: int, col: int
    ) -> None:
        key = (conflict_type, tuple(sorted((first, second))))
        if key not in self.conflicts:
            self.conflicts[key] = {
                "type": conflict_type,
                "parcel_ids": list(key[1]),
                "row": int(row),
                "column": int(col),
            }


# ---------------------------------------------------------------------------
# Raster validation
# ---------------------------------------------------------------------------


def open_parcel_raster(path: str) -> RasterSpec:
    """Validate the structural contract of the parcel ID map."""

    try:
        ds = rasterio.open(path)
    except Exception as exc:  # driver-specific errors share no common type
        raise BoundaryError(
            f"地块编号图无法作为 GeoTIFF 打开：{exc}",
        ) from exc

    with ds:
        if ds.driver.lower() != "gtiff":
            raise BoundaryError(f"地块编号图格式不是 GeoTIFF（driver={ds.driver}）。")
        if ds.count != 1:
            raise BoundaryError(
                f"地块编号图必须为单波段，实际有 {ds.count} 个波段。"
            )
        dtype = ds.dtypes[0]
        if dtype not in ("uint16", "uint32"):
            raise BoundaryError(
                f"地块编号图数据类型必须是 UInt16 或 UInt32，实际为 {dtype}。"
            )
        if ds.nodata != 0:
            raise BoundaryError(
                f"地块编号图的 NoData 必须为 0，实际为 {ds.nodata!r}。"
            )
        if ds.crs is None:
            raise BoundaryError("地块编号图缺少坐标系（CRS）。")
        t = ds.transform
        if abs(t.b) > 1e-12 or abs(t.d) > 1e-12:
            raise BoundaryError("地块编号图不是北向上（仿射变换含旋转/剪切）。")
        if t.a <= 0 or t.e >= 0:
            raise BoundaryError("地块编号图不是北向上（像元宽应为正、高应为负）。")
        if not math.isfinite(t.a) or not math.isfinite(t.e):
            raise BoundaryError("地块编号图的像元分辨率无效。")

        return RasterSpec(
            path=path,
            width=ds.width,
            height=ds.height,
            crs=ds.crs,
            transform=t,
            dtype=dtype,
            nodata=ds.nodata,
            id_max=int(np.iinfo(np.dtype(dtype)).max),
        )


# ---------------------------------------------------------------------------
# GeoJSON parsing and reprojection
# ---------------------------------------------------------------------------


def _feature_id(feature: dict[str, Any], index: int) -> tuple[int | None, str | None]:
    properties = feature.get("properties") or {}
    raw = properties.get("parcel_id", feature.get("id"))
    if raw is None:
        return None, f"第 {index + 1} 个要素缺少地块编号（properties.parcel_id）。"
    # Reject bool (a bool is an int in Python) and non-integral values.
    if isinstance(raw, bool) or not isinstance(raw, int):
        return None, (
            f"第 {index + 1} 个要素的地块编号必须是整数，实际为 {raw!r}。"
        )
    if raw <= 0:
        return None, f"第 {index + 1} 个要素的地块编号必须为正整数，实际为 {raw}。"
    return int(raw), None


def _check_named_crs(payload: dict[str, Any]) -> str | None:
    """Honour the legacy GeoJSON ``crs`` member; only WGS84 is accepted."""

    crs_member = payload.get("crs")
    if not crs_member:
        return None
    name = ""
    try:
        name = crs_member["properties"]["name"]
    except (KeyError, TypeError):
        return "GeoJSON 的 crs 成员无法识别；边界必须使用 WGS84（EPSG:4326）。"
    normalized = name.replace(" ", "").lower()
    if normalized in ("urn:ogc:def:crs:ogc:1.3:crs84", "urn:ogc:def:crs:epsg::4326"):
        return None
    if "4326" in normalized or normalized.endswith(":crs84"):
        return None
    return f"GeoJSON 坐标系必须是 WGS84（EPSG:4326），实际声明为 {name!r}。"


def _polygon_rings(geometry: dict[str, Any]) -> list[list[tuple[float, float]]] | str:
    if geometry.get("type") != "Polygon":
        return (
            f"几何类型必须是 Polygon，实际为 {geometry.get('type')!r}；"
            "每个地块编号对应一个 Polygon。"
        )
    coordinates = geometry.get("coordinates")
    if not isinstance(coordinates, list) or not coordinates:
        return "Polygon 缺少坐标环。"
    rings: list[list[tuple[float, float]]] = []
    for ring in coordinates:
        if not isinstance(ring, list) or len(ring) < 4:
            return "Polygon 的坐标环至少需要 4 个点（含闭合点）。"
        cleaned: list[tuple[float, float]] = []
        for point in ring:
            if not isinstance(point, list) or len(point) < 2:
                return "坐标点必须至少包含经度、纬度两个数值。"
            lon, lat = point[0], point[1]
            if not isinstance(lon, (int, float)) or not isinstance(lat, (int, float)):
                return "坐标点的经度、纬度必须是数值。"
            if not math.isfinite(lon) or not math.isfinite(lat):
                return "坐标点包含非有限数值。"
            # Drop Z/M dimensions; the service works in 2D.
            cleaned.append((float(lon), float(lat)))
        if cleaned[0] != cleaned[-1]:
            cleaned.append(cleaned[0])
        rings.append(cleaned)
    return rings


def parse_boundary_geojson(
    text: str | bytes, spec: RasterSpec
) -> list[BoundaryUpdate]:
    """Parse and reproject the WGS84 GeoJSON into raster pixel space."""

    try:
        payload = json.loads(text)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise BoundaryError(f"GeoJSON 无法解析：{exc}") from exc
    if not isinstance(payload, dict):
        raise BoundaryError("GeoJSON 根节点必须是对象（FeatureCollection 或 Feature）。")

    errors: list[str] = []
    crs_error = _check_named_crs(payload)
    if crs_error:
        errors.append(crs_error)

    if payload.get("type") == "Feature":
        features = [payload]
    elif payload.get("type") == "FeatureCollection":
        features = payload.get("features")
        if not isinstance(features, list):
            raise BoundaryError("FeatureCollection 缺少 features 数组。")
    else:
        raise BoundaryError(
            "GeoJSON 必须是 FeatureCollection（或单个 Feature），"
            f"实际类型为 {payload.get('type')!r}。"
        )

    if not features:
        errors.append("GeoJSON 不包含任何要素。")
    if len(features) > _MAX_PARCELS:
        errors.append(f"一次最多更正 {_MAX_PARCELS} 块地，实际提交了 {len(features)} 个。")

    parsed: list[tuple[int, list[list[tuple[float, float]]]]] = []
    seen_ids: set[int] = set()
    for index, feature in enumerate(features[:_MAX_PARCELS]):
        if not isinstance(feature, dict) or feature.get("type") != "Feature":
            errors.append(f"第 {index + 1} 个要素不是 Feature 对象。")
            continue
        parcel_id, id_error = _feature_id(feature, index)
        if id_error:
            errors.append(id_error)
        geometry = feature.get("geometry")
        if not isinstance(geometry, dict):
            errors.append(f"第 {index + 1} 个要素缺少 Polygon 几何。")
            continue
        rings = _polygon_rings(geometry)
        if isinstance(rings, str):
            errors.append(f"第 {index + 1} 个要素：{rings}")
            continue
        if parcel_id is None:
            continue
        if parcel_id in seen_ids:
            errors.append(f"地块编号 {parcel_id} 在 GeoJSON 中重复出现。")
            continue
        if parcel_id > spec.id_max:
            errors.append(
                f"地块编号 {parcel_id} 超出地块图数据类型 {spec.dtype} 的上限"
                f"（{spec.id_max}）。"
            )
            continue
        seen_ids.add(parcel_id)
        parsed.append((parcel_id, rings))

    if errors:
        raise BoundaryError("GeoJSON 校验未通过，未执行边界更正。", geojson_errors=errors)

    source_crs = CRS.from_epsg(4326)
    inverse = ~spec.transform
    updates: list[BoundaryUpdate] = []
    for parcel_id, rings in parsed:
        # transform_geom densifies long edges while keeping every ring
        # (exterior plus holes) in order.
        wgs84_geom = {
            "type": "Polygon",
            "coordinates": [[list(point) for point in ring] for ring in rings],
        }
        projected = transform_geom(source_crs, spec.crs, wgs84_geom)

        # Normalize winding on the projected (Cartesian) coordinates: exterior
        # counter-clockwise, holes clockwise.  Inputs may use either winding.
        normalized_rings: list[list[list[float]]] = []
        for index, ring in enumerate(projected["coordinates"]):
            ring_arr = np.asarray(ring, dtype=np.float64)
            area = _signed_area(ring_arr)
            should_be_ccw = index == 0
            needs_flip = (area < 0) if should_be_ccw else (area > 0)
            normalized_rings.append(ring[::-1] if needs_flip else [list(p) for p in ring])
        projected["coordinates"] = normalized_rings

        pixel_rings: list[np.ndarray] = []
        for ring in projected["coordinates"]:
            xs = [point[0] for point in ring]
            ys = [point[1] for point in ring]
            cols, rows = inverse @ (np.asarray(xs), np.asarray(ys))
            pixel_rings.append(np.column_stack([np.asarray(cols), np.asarray(rows)]))
        updates.append(BoundaryUpdate(parcel_id=parcel_id, pixel_rings=pixel_rings))
    return updates


# ---------------------------------------------------------------------------
# Pixel-center hit testing
# ---------------------------------------------------------------------------


def _on_ring(
    cols: np.ndarray, rows: np.ndarray, ring: np.ndarray
) -> np.ndarray:
    """Vectorized point-to-segment distance test against one closed ring."""

    hit = np.zeros(cols.shape, dtype=bool)
    for start, end in zip(ring[:-1], ring[1:]):
        ax, ay = start
        bx, by = end
        dx, dy = bx - ax, by - ay
        length_sq = dx * dx + dy * dy
        if length_sq == 0.0:
            hit |= (np.abs(cols - ax) <= _EDGE_EPS) & (
                np.abs(rows - ay) <= _EDGE_EPS
            )
            continue
        t = ((cols - ax) * dx + (rows - ay) * dy) / length_sq
        t = np.clip(t, 0.0, 1.0)
        nearest_x = ax + t * dx
        nearest_y = ay + t * dy
        hit |= (cols - nearest_x) ** 2 + (rows - nearest_y) ** 2 <= _EDGE_EPS**2
    return hit


def _strictly_inside(
    cols: np.ndarray, rows: np.ndarray, ring: np.ndarray
) -> np.ndarray:
    """Even-odd ray casting with the half-open edge rule.

    An edge straddles the horizontal ray through a point when one endpoint
    is at or above the point and the other is strictly below it.  This makes
    polygons behave like half-open pixel regions: a point on the top or
    right boundary is not counted as strictly inside (boundary points are
    captured separately by :func:`_on_ring`), while bottom/left edges are
    included.
    """

    inside = np.zeros(cols.shape, dtype=bool)
    for start, end in zip(ring[:-1], ring[1:]):
        ax, ay = start
        bx, by = end
        if ay == by:
            continue
        y_min, y_max = min(ay, by), max(ay, by)
        crosses = (rows >= y_min) & (rows < y_max)
        if not crosses.any():
            continue
        x_cross = ax + (rows - ay) * (bx - ax) / (by - ay)
        inside ^= crosses & (x_cross > cols)
    return inside


def _signed_area(ring: np.ndarray) -> float:
    """Shoelace signed area; positive for counter-clockwise rings."""

    xs, ys = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.sum(xs[:-1] * ys[1:] - xs[1:] * ys[:-1]))


def _hit_pixels(
    candidate_cols: np.ndarray,
    candidate_rows: np.ndarray,
    rings: list[np.ndarray],
) -> np.ndarray:
    """Center inside or on the polygon boundary, holes excluded."""

    on_boundary = np.zeros(candidate_cols.shape, dtype=bool)
    for ring in rings:
        on_boundary |= _on_ring(candidate_cols, candidate_rows, ring)

    inside = _strictly_inside(candidate_cols, candidate_rows, rings[0])
    for hole in rings[1:]:
        inside &= ~_strictly_inside(candidate_cols, candidate_rows, hole)
    return on_boundary | inside


def _rasterize_candidates(
    rings: list[np.ndarray], shape: tuple[int, int]
) -> np.ndarray:
    geometry = {
        "type": "Polygon",
        "coordinates": [ring.tolist() for ring in rings],
    }
    return rasterize(
        shapes=[(geometry, 1)],
        out_shape=shape,
        transform=_IDENTITY,
        dtype="uint8",
        fill=0,
        all_touched=True,
    )


def _window_hit(
    update: BoundaryUpdate, window: Window
) -> tuple[slice, slice, np.ndarray] | None:
    """Boolean hit mask of one polygon within its intersection with a window."""

    exterior = update.pixel_rings[0]
    min_col, min_row = exterior.min(axis=0)
    max_col, max_row = exterior.max(axis=0)

    # Generous margin around the polygon bounds: candidate selection is only
    # a coarse pre-filter, while the exact center-in/on test runs on every
    # candidate.  A two-pixel margin guarantees that centers lying on an edge
    # (up to half a pixel outside the geometry bounds) are always considered,
    # regardless of how rasterize() handles corner-only pixel contact.
    col0 = max(int(window.col_off), int(math.floor(min_col)) - 2)
    row0 = max(int(window.row_off), int(math.floor(min_row)) - 2)
    col1 = min(int(window.col_off + window.width), int(math.ceil(max_col)) + 2)
    row1 = min(int(window.row_off + window.height), int(math.ceil(max_row)) + 2)
    if col1 <= col0 or row1 <= row0:
        return None

    local_rings = []
    for ring in update.pixel_rings:
        shifted = ring.copy()
        shifted[:, 0] -= col0
        shifted[:, 1] -= row0
        local_rings.append(shifted)

    sub_height, sub_width = row1 - row0, col1 - col0
    candidates = _rasterize_candidates(local_rings, (sub_height, sub_width))
    cand_rows, cand_cols = np.nonzero(candidates)
    hit = np.zeros((sub_height, sub_width), dtype=bool)
    if cand_rows.size:
        hit[cand_rows, cand_cols] = _hit_pixels(
            cand_cols.astype(np.float64),
            cand_rows.astype(np.float64),
            local_rings,
        )
    if not hit.any():
        return None

    row_slice = slice(row0 - int(window.row_off), row1 - int(window.row_off))
    col_slice = slice(col0 - int(window.col_off), col1 - int(window.col_off))
    return row_slice, col_slice, hit


# ---------------------------------------------------------------------------
# Dry pass and write pass
# ---------------------------------------------------------------------------


def _block_windows(ds: rasterio.DatasetReader, block_size: int):
    for _, block_window in ds.block_windows(1):
        yield from _split_window(block_window, max(1, block_size))


def _first_pixel(
    mask: np.ndarray, window: Window, row_slice: slice, col_slice: slice
) -> tuple[int, int]:
    rows, cols = np.nonzero(mask)
    return (
        int(window.row_off) + row_slice.start + int(rows[0]),
        int(window.col_off) + col_slice.start + int(cols[0]),
    )


def _scan(spec: RasterSpec, updates: list[BoundaryUpdate], block_size: int) -> _ScanFindings:
    findings = _ScanFindings()
    findings.hit_counts = {update.parcel_id: 0 for update in updates}
    batch_id_array = np.fromiter(
        (update.parcel_id for update in updates), dtype=np.int64, count=len(updates)
    )
    batch_ids = set(batch_id_array.tolist())

    with rasterio.open(spec.path) as ds:
        for window in _block_windows(ds, block_size):
            tile = ds.read(1, window=window)
            positive = tile[tile > 0]
            if positive.size:
                findings.existing_ids.update(
                    int(value) for value in np.unique(positive)
                )

            owner = np.full(tile.shape, -1, dtype=np.int64)
            for update in updates:
                result = _window_hit(update, window)
                if result is None:
                    continue
                row_slice, col_slice, hit = result
                findings.hit_counts[update.parcel_id] += int(hit.sum())

                region_owner = owner[row_slice, col_slice]
                region_old = tile[row_slice, col_slice].astype(np.int64)

                # Two new boundaries claiming the same pixel.
                overlap = hit & (region_owner >= 0)
                if overlap.any():
                    for other_id in np.unique(region_owner[overlap]).tolist():
                        pair_mask = overlap & (region_owner == other_id)
                        row, col = _first_pixel(
                            pair_mask, window, row_slice, col_slice
                        )
                        findings.remember_conflict(
                            "new_new", int(other_id), update.parcel_id, row, col
                        )

                # A new boundary occupying a pixel of a parcel not in the
                # batch; old pixels of batch parcels are reusable.
                invaded = (
                    hit
                    & (region_owner < 0)
                    & (region_old > 0)
                    & (region_old != update.parcel_id)
                )
                invaded &= ~np.isin(region_old, batch_id_array)
                if invaded.any():
                    for other_id in np.unique(region_old[invaded]).tolist():
                        pair_mask = invaded & (region_old == other_id)
                        row, col = _first_pixel(
                            pair_mask, window, row_slice, col_slice
                        )
                        findings.remember_conflict(
                            "new_old", update.parcel_id, int(other_id), row, col
                        )

                region_owner[hit] = update.parcel_id
    return findings


def _assert_clean(spec: RasterSpec, updates: list[BoundaryUpdate], findings: _ScanFindings) -> None:
    unknown = sorted(
        update.parcel_id
        for update in updates
        if update.parcel_id not in findings.existing_ids
    )
    no_hit = sorted(
        parcel_id
        for parcel_id, count in findings.hit_counts.items()
        if count == 0
    )
    conflicts = [findings.conflicts[key] for key in sorted(findings.conflicts)]
    if unknown or no_hit or conflicts:
        raise BoundaryError(
            "边界核对未通过，未生成新的地块编号图。",
            unknown_parcel_ids=unknown,
            no_hit_parcel_ids=no_hit,
            conflicts=conflicts,
        )


def apply_boundary_updates(
    spec: RasterSpec,
    updates: list[BoundaryUpdate],
    output_path: str,
    block_size: int = _DEFAULT_BLOCK_SIZE,
) -> None:
    """Validate the whole batch, then rewrite the map block by block."""

    findings = _scan(spec, updates, block_size)
    _assert_clean(spec, updates, findings)

    batch_ids = np.array(
        sorted(update.parcel_id for update in updates), dtype=np.int64
    )
    with rasterio.open(spec.path) as src:
        profile = src.profile.copy()
        with rasterio.open(output_path, "w", **profile) as dst:
            for window in _block_windows(src, block_size):
                tile = src.read(1, window=window)
                result_tile = tile.copy()

                # Old pixels of every corrected parcel are released first;
                # new boundaries then stamp their own parcel IDs.
                result_tile[np.isin(tile.astype(np.int64), batch_ids)] = 0

                owner = np.full(tile.shape, -1, dtype=np.int64)
                for update in updates:
                    hit = _window_hit(update, window)
                    if hit is None:
                        continue
                    row_slice, col_slice, hit_mask = hit
                    owner[row_slice, col_slice][hit_mask] = update.parcel_id

                stamped = owner >= 0
                result_tile[stamped] = owner[stamped].astype(result_tile.dtype)
                dst.write(result_tile, 1, window=window)
