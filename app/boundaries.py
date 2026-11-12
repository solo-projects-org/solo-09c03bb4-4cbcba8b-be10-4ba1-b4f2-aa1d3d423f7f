"""Windowed parcel-boundary correction from WGS84 GeoJSON polygons.

The caller uploads the existing parcel map (a north-up, single-band UInt16/
UInt32 GeoTIFF with NoData 0) together with a WGS84 GeoJSON whose features
map parcel IDs to replacement polygons.  This module:

* validates the raster and the GeoJSON,
* reprojects the polygons to the raster CRS and rasterizes them by pixel
  centres (a centre on the polygon edge still counts, via ``covers``),
* checks the whole batch for conflicts in bounded tiles,
* writes the corrected map only when every check passed.

Nothing is returned as an image when any check fails.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import rasterio
import shapely
from affine import Affine
from rasterio.crs import CRS
from rasterio.warp import transform_geom
from rasterio.windows import Window
from shapely.geometry import shape as shapely_shape
from shapely.prepared import prep as shapely_prep

from app.analysis import RasterMeta, UploadedRaster, _inspect, _split_window

DEFAULT_BLOCK_SIZE = 1024
MAX_PARCELS = 32
# Maximum segment length, in geographic degrees, before reprojecting a ring.
# Densification keeps the reprojected polygon faithful to the boundary
# (0.01 degree is roughly 1.1 km or less, reprojected error well under a pixel
# for the 10 m-class parcel maps this service handles).
DENSIFY_MAX_SEGMENT_DEG = 0.01
_GEOJSON_CRS = CRS.from_epsg(4326)
_ALLOWED_DTYPES = ("uint16", "uint32")


class BoundaryError(Exception):
    """A structured, batch-level validation error.

    ``kind`` selects the JSON field layout the API layer emits.  Conflicts
    carry the two parcel IDs and the row/column of one conflicting pixel.
    """

    def __init__(
        self,
        kind: str,
        message: str,
        *,
        parcel_id: int | None = None,
        parcel_ids: tuple[int, int] | None = None,
        row: int | None = None,
        column: int | None = None,
        related: Sequence["BoundaryError"] | None = None,
    ) -> None:
        super().__init__(message)
        self.kind = kind
        self.message = message
        self.parcel_id = parcel_id
        self.parcel_ids = parcel_ids
        self.row = row
        self.column = column
        self.related: list[BoundaryError] = list(related or [])

    def as_dict(self) -> dict[str, object]:
        data: dict[str, object] = {"kind": self.kind, "message": self.message}
        if self.parcel_id is not None:
            data["parcel_id"] = self.parcel_id
        if self.parcel_ids is not None:
            data["parcel_ids"] = list(self.parcel_ids)
        if self.column is not None:
            data["column"] = self.column
        if self.row is not None:
            data["row"] = self.row
        if self.related:
            data["related_errors"] = [e.as_dict() for e in self.related]
        return data


@dataclass(frozen=True)
class ParcelPolygon:
    parcel_id: int
    geometry: object  # shapely Polygon in the raster CRS
    prepared: object  # shapely PreparedGeometry


@dataclass
class CorrectionPlan:
    raster_meta: RasterMeta
    parcels: list[ParcelPolygon]
    block_size: int = DEFAULT_BLOCK_SIZE
    # Pixel-space bounding boxes clipped to the raster extent, or None when the
    # polygon cannot cover any pixel centre inside the map.
    boxes: dict[int, tuple[int, int, int, int] | None] = field(default_factory=dict)

    def box(self, parcel: ParcelPolygon) -> tuple[int, int, int, int] | None:
        return self.boxes[parcel.parcel_id]


def inspect_parcel_raster(upload: UploadedRaster) -> RasterMeta:
    """Open the parcel map and enforce the boundary-correction contract."""

    meta, error = _inspect(upload)
    if error is not None:
        raise BoundaryError("invalid_raster", f"地块图不兼容：{error.reason}")
    assert meta is not None

    if meta.dtype not in _ALLOWED_DTYPES:
        raise BoundaryError(
            "invalid_raster",
            f"地块图数据类型必须是 UInt16 或 UInt32，实际为 {meta.dtype}。",
        )
    if meta.nodata is None or int(meta.nodata) != 0:
        actual = "未设置" if meta.nodata is None else str(meta.nodata)
        raise BoundaryError(
            "invalid_raster",
            f"地块图 NoData 必须为 0，实际为 {actual}。",
        )
    return meta


def _feature_id(feature: dict, index: int) -> int:
    props = feature.get("properties") or {}
    raw = props.get("parcel_id", props.get("id"))
    if raw is None and feature.get("id") is not None:
        raw = feature.get("id")
    if raw is None:
        raise BoundaryError(
            "invalid_geojson",
            f"第 {index + 1} 个 Feature 缺少地块编号"
            "（需要 properties.parcel_id 或数字类型的 Feature id）。",
        )
    if isinstance(raw, bool) or not isinstance(raw, int):
        if isinstance(raw, float) and math.isfinite(raw) and raw == int(raw):
            raw = int(raw)
        else:
            raise BoundaryError(
                "invalid_parcel_id",
                f"第 {index + 1} 个 Feature 的地块编号必须是正整数，实际为 {raw!r}。",
            )
    if raw <= 0:
        raise BoundaryError(
            "invalid_parcel_id",
            f"第 {index + 1} 个 Feature 的地块编号必须是正整数，实际为 {raw}。",
        )
    return int(raw)


def _check_geojson_crs(payload: dict) -> None:
    if "crs" not in payload:
        return  # RFC 7946 default: WGS84
    crs_obj = payload.get("crs") or {}
    name = ""
    if isinstance(crs_obj, dict):
        name = str((crs_obj.get("properties") or {}).get("name") or "")
    normalized = name.lower().replace("epsg::", "epsg:")
    accepted = {
        "epsg:4326",
        "urn:ogc:def:crs:ogc:1.3:crs84",
        "urn:ogc:def:crs:epsg::4326",
        "wgs84",
        "wgs 84",
    }
    if normalized not in accepted:
        raise BoundaryError(
            "invalid_geojson",
            f"GeoJSON 坐标系必须是 WGS84（EPSG:4326），实际声明为 {name or '未知'}。",
        )


def _densify_ring(ring: Sequence[Sequence[float]], max_segment_deg: float) -> list[list[float]]:
    """Insert vertices so no ring segment spans more than ``max_segment_deg``."""

    out: list[list[float]] = []
    for start, end in zip(ring, ring[1:]):
        out.append([start[0], start[1]])
        lon_span = abs(end[0] - start[0])
        lat_span = abs(end[1] - start[1])
        steps = max(
            int(math.ceil(lon_span / max_segment_deg)),
            int(math.ceil(lat_span / max_segment_deg)),
        )
        if steps > 1:
            for step in range(1, steps):
                ratio = step / steps
                out.append(
                    [
                        start[0] + (end[0] - start[0]) * ratio,
                        start[1] + (end[1] - start[1]) * ratio,
                    ]
                )
    if ring:
        out.append([ring[-1][0], ring[-1][1]])
    return out


def _densify_geometry(geometry: dict, max_segment_deg: float) -> dict:
    geometry_type = geometry.get("type")
    if geometry_type == "Polygon":
        return {
            "type": "Polygon",
            "coordinates": [
                _densify_ring(ring, max_segment_deg)
                for ring in geometry["coordinates"]
            ],
        }
    if geometry_type == "MultiPolygon":
        return {
            "type": "MultiPolygon",
            "coordinates": [
                [_densify_ring(ring, max_segment_deg) for ring in polygon]
                for polygon in geometry["coordinates"]
            ],
        }
    return geometry


def _reproject_polygon(feature: dict, parcel_id: int, target_crs: CRS) -> object:
    geometry = feature.get("geometry")
    if not isinstance(geometry, dict):
        raise BoundaryError(
            "invalid_geometry",
            f"地块 {parcel_id} 的 geometry 缺失或不是 GeoJSON 几何对象。",
            parcel_id=parcel_id,
        )
    if geometry.get("type") != "Polygon":
        raise BoundaryError(
            "invalid_geometry",
            f"地块 {parcel_id} 的几何必须是 Polygon，"
            f"实际为 {geometry.get('type')!r}。",
            parcel_id=parcel_id,
        )

    # Densify in geographic space first so long edges do not cut corners after
    # the nonlinear reprojection.
    dense = _densify_geometry(geometry, DENSIFY_MAX_SEGMENT_DEG)
    try:
        projected = transform_geom(_GEOJSON_CRS, target_crs, dense)
    except Exception as exc:  # driver-specific transformation failures
        raise BoundaryError(
            "reprojection_failed",
            f"地块 {parcel_id} 的边界无法重投影到底图坐标系：{exc}",
            parcel_id=parcel_id,
        ) from exc

    try:
        polygon = shapely_shape(projected)
    except Exception as exc:
        raise BoundaryError(
            "invalid_geometry",
            f"地块 {parcel_id} 重投影后的几何无效：{exc}",
            parcel_id=parcel_id,
        ) from exc
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    if polygon.is_empty or polygon.geom_type != "Polygon":
        raise BoundaryError(
            "invalid_geometry",
            f"地块 {parcel_id} 的几何重投影后不是有效的单个 Polygon。",
            parcel_id=parcel_id,
        )
    return polygon


def parse_geojson(raw_bytes: bytes, meta: RasterMeta) -> list[ParcelPolygon]:
    """Parse and reproject every feature; validate the batch structurally."""

    try:
        payload = json.loads(raw_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BoundaryError("invalid_geojson", f"GeoJSON 无法解析：{exc}") from exc

    if not isinstance(payload, dict) or payload.get("type") != "FeatureCollection":
        raise BoundaryError(
            "invalid_geojson",
            "GeoJSON 顶层必须是 FeatureCollection。",
        )
    features = payload.get("features")
    if not isinstance(features, list) or not features:
        raise BoundaryError(
            "invalid_geojson",
            "FeatureCollection 至少要包含一个 Feature。",
        )
    if len(features) > MAX_PARCELS:
        raise BoundaryError(
            "too_many_parcels",
            f"一次最多更正 {MAX_PARCELS} 块地，本次提交了 {len(features)} 个。",
        )
    _check_geojson_crs(payload)

    parcels: list[ParcelPolygon] = []
    seen: set[int] = set()
    for index, feature in enumerate(features):
        if not isinstance(feature, dict) or feature.get("type") != "Feature":
            raise BoundaryError(
                "invalid_geojson",
                f"第 {index + 1} 个元素不是 Feature 对象。",
            )
        parcel_id = _feature_id(feature, index)
        if parcel_id in seen:
            raise BoundaryError(
                "duplicate_parcel_id",
                f"地块编号 {parcel_id} 在 GeoJSON 中出现了多次，每个编号只能有一个 Polygon。",
                parcel_id=parcel_id,
            )
        seen.add(parcel_id)
        polygon = _reproject_polygon(feature, parcel_id, meta.crs)
        parcels.append(
            ParcelPolygon(
                parcel_id=parcel_id,
                geometry=polygon,
                prepared=shapely_prep(polygon),
            )
        )
    parcels.sort(key=lambda p: p.parcel_id)
    return parcels


def _pixel_box(
    polygon: object, transform: Affine, width: int, height: int
) -> tuple[int, int, int, int] | None:
    """Return (col_off, row_off, width, height) of candidate pixel centres.

    The box spans every pixel whose centre lies within the geometry's bounds;
    point-in-polygon testing later decides actual membership.  Centres exactly
    on a bound are included because boundary coverage counts as a hit.
    """

    minx, miny, maxx, maxy = polygon.bounds
    # Centre of pixel (r, c): x = c0 + a*(c + 0.5), y = r0 + e*(r + 0.5).
    first_col = math.ceil((minx - transform.c) / transform.a - 0.5)
    last_col = math.floor((maxx - transform.c) / transform.a - 0.5)
    # transform.e is negative for north-up images, so the y inequality flips.
    first_row = math.ceil((maxy - transform.f) / transform.e - 0.5)
    last_row = math.floor((miny - transform.f) / transform.e - 0.5)

    col_off = max(0, first_col)
    row_off = max(0, first_row)
    col_end = min(width - 1, last_col)
    row_end = min(height - 1, last_row)
    if col_end < col_off or row_end < row_off:
        return None
    return col_off, row_off, col_end - col_off + 1, row_end - row_off + 1


def build_plan(
    meta: RasterMeta,
    parcels: Sequence[ParcelPolygon],
    block_size: int = DEFAULT_BLOCK_SIZE,
) -> CorrectionPlan:
    plan = CorrectionPlan(
        raster_meta=meta,
        parcels=list(parcels),
        block_size=max(1, block_size),
    )
    for parcel in parcels:
        plan.boxes[parcel.parcel_id] = _pixel_box(
            parcel.geometry, meta.transform, meta.width, meta.height
        )
    return plan


def _centre_points(
    window: Window, transform: Affine
) -> tuple[np.ndarray, np.ndarray]:
    """World coordinates of every pixel centre in a window."""

    col_idx, row_idx = np.meshgrid(
        np.arange(int(window.col_off), int(window.col_off) + int(window.width)),
        np.arange(int(window.row_off), int(window.row_off) + int(window.height)),
    )
    xs = transform.c + transform.a * (col_idx + 0.5)
    ys = transform.f + transform.e * (row_idx + 0.5)
    return xs, ys


def _box_intersects(
    box: tuple[int, int, int, int], window: Window
) -> tuple[int, int, int, int] | None:
    """Intersect a pixel box with a window, returned in window-local coords."""

    col_off, row_off, width, height = box
    local_col = max(0, col_off - int(window.col_off))
    local_row = max(0, row_off - int(window.row_off))
    end_col = min(int(window.width), col_off + width - int(window.col_off))
    end_row = min(int(window.height), row_off + height - int(window.row_off))
    if end_col <= local_col or end_row <= local_row:
        return None
    return local_row, local_col, end_row, end_col


def apply_correction(plan: CorrectionPlan, destination: str) -> None:
    """Validate the whole batch and write the corrected GeoTIFF.

    The raster is read and written in bounded tiles.  On any batch error the
    destination file is left for the caller to delete and no image is served.
    """

    meta = plan.raster_meta
    corrected_ids = {p.parcel_id for p in plan.parcels}
    dtype_max = int(np.iinfo(np.dtype(meta.dtype)).max)
    errors: list[BoundaryError] = []

    # IDs that cannot be represented in the output dtype must fail up front,
    # before a single tile is written.
    for parcel in plan.parcels:
        if parcel.parcel_id > dtype_max:
            errors.append(
                BoundaryError(
                    "id_out_of_range",
                    f"地块编号 {parcel.parcel_id} 超出 {meta.dtype} 可表示的"
                    f"最大值 {dtype_max}。",
                    parcel_id=parcel.parcel_id,
                )
            )
    if errors:
        raise BoundaryError(
            "invalid_batch",
            "批次中存在无法写入的地块编号，未执行更正。",
            related=errors,
        )

    hit_counts: dict[int, int] = {p.parcel_id: 0 for p in plan.parcels}
    present_ids: set[int] = set()
    # One representative conflicting pixel per parcel pair keeps the report
    # bounded (at most 32*31/2 pairs) while still naming every clashing pair.
    new_new_pairs: set[tuple[int, int]] = set()
    new_old_pairs: set[tuple[int, int]] = set()
    corrected_id_array = np.array(sorted(corrected_ids), dtype=np.int64)

    with rasterio.open(meta.upload.path) as source:
        profile = source.profile.copy()
        # Preserve the original dimensions, grid, CRS and data type verbatim.
        profile.update(
            driver="GTiff",
            count=1,
            width=meta.width,
            height=meta.height,
            dtype=meta.dtype,
            crs=meta.crs,
            transform=meta.transform,
            nodata=0,
        )
        with rasterio.open(destination, "w", **profile) as dst:
            for _idx, block in source.block_windows(1):
                for window in _split_window(block, max(1, plan.block_size)):
                    original = source.read(1, window=window)
                    corrected = original.copy()
                    xs, ys = _centre_points(window, meta.transform)

                    members: list[tuple[ParcelPolygon, np.ndarray]] = []
                    coverage_count = np.zeros(original.shape, dtype=np.uint8)
                    combined = np.zeros(original.shape, dtype=bool)
                    for parcel in plan.parcels:
                        box = plan.box(parcel)
                        if box is None:
                            continue
                        intersection = _box_intersects(box, window)
                        if intersection is None:
                            continue
                        local_row, local_col, end_row, end_col = intersection
                        # Test only the candidate points inside this parcel's
                        # bounding box, then scatter the result back.
                        sub_points = shapely.points(
                            xs[local_row:end_row, local_col:end_col].ravel(),
                            ys[local_row:end_row, local_col:end_col].ravel(),
                        )
                        # ``covers`` (not ``contains``) makes a centre lying
                        # exactly on the polygon edge count as a hit.  The
                        # prepared geometry's predicate is vectorized over the
                        # candidate point array.
                        sub_hits = np.asarray(
                            parcel.prepared.covers(sub_points), dtype=bool
                        ).reshape((end_row - local_row, end_col - local_col))
                        if not sub_hits.any():
                            continue
                        hits = np.zeros(original.shape, dtype=bool)
                        hits[local_row:end_row, local_col:end_col] = sub_hits
                        members.append((parcel, hits))
                        coverage_count += hits.astype(np.uint8)
                        combined |= hits
                        hit_counts[parcel.parcel_id] += int(hits.sum())

                    # Two replacement polygons claiming the same pixel centre.
                    if (coverage_count > 1).any():
                        clash_locations = np.argwhere(coverage_count > 1)
                        reported_here: set[tuple[int, int]] = set()
                        for local in clash_locations:
                            owners = tuple(
                                parcel.parcel_id
                                for parcel, mask in members
                                if mask[local[0], local[1]]
                            )
                            pair = tuple(sorted(owners))
                            if pair in new_new_pairs or pair in reported_here:
                                continue
                            reported_here.add(pair)
                            row = int(window.row_off) + int(local[0])
                            col = int(window.col_off) + int(local[1])
                            errors.append(
                                BoundaryError(
                                    "new_new_conflict",
                                    f"地块 {pair[0]} 与 {pair[1]} 的新边界同时命中"
                                    f"像元（行 {row}，列 {col}）。",
                                    parcel_ids=(pair[0], pair[1]),
                                    row=row,
                                    column=col,
                                )
                            )
                        new_new_pairs.update(reported_here)

                    # Old cells of corrected parcels are released to this
                    # batch: any not re-hit by a new boundary become NoData 0.
                    old_corrected = np.isin(
                        original.astype(np.int64), corrected_id_array
                    )
                    corrected[old_corrected & ~combined] = 0
                    for parcel, mask in members:
                        corrected[mask] = parcel.parcel_id

                    # A new boundary taking over a pixel that still belongs to
                    # a parcel outside this correction batch.
                    if combined.any():
                        foreign_map = original.astype(np.int64).copy()
                        foreign_map[old_corrected] = 0
                        foreign_hits = combined & (foreign_map > 0)
                        if foreign_hits.any():
                            reported_here: set[tuple[int, int]] = set()
                            for local in np.argwhere(foreign_hits):
                                foreign_id = int(foreign_map[local[0], local[1]])
                                new_owner = next(
                                    parcel.parcel_id
                                    for parcel, mask in members
                                    if mask[local[0], local[1]]
                                )
                                pair = tuple(sorted((new_owner, foreign_id)))
                                if pair in new_old_pairs or pair in reported_here:
                                    continue
                                reported_here.add(pair)
                                row = int(window.row_off) + int(local[0])
                                col = int(window.col_off) + int(local[1])
                                errors.append(
                                    BoundaryError(
                                        "new_old_conflict",
                                        f"地块 {new_owner} 的新边界占用了未更正地块 "
                                        f"{foreign_id} 的像元（行 {row}，列 {col}）。",
                                        parcel_ids=(new_owner, foreign_id),
                                        row=row,
                                        column=col,
                                    )
                                )
                            new_old_pairs.update(reported_here)

                    for value in np.unique(original):
                        value = int(value)
                        if value > 0:
                            present_ids.add(value)
                    dst.write(corrected, 1, window=window)

    for parcel_id in sorted(corrected_ids - present_ids):
        errors.append(
            BoundaryError(
                "unknown_parcel_id",
                f"地块编号 {parcel_id} 在原地块图中不存在。",
                parcel_id=parcel_id,
            )
        )
    for parcel in plan.parcels:
        if hit_counts[parcel.parcel_id] == 0:
            errors.append(
                BoundaryError(
                    "no_pixels_hit",
                    f"地块 {parcel.parcel_id} 的新边界在地块图范围内未命中任何像元。",
                    parcel_id=parcel.parcel_id,
                )
            )

    if errors:
        raise BoundaryError(
            "invalid_batch",
            "边界更正批次核对未通过，未输出影像。",
            related=errors,
        )