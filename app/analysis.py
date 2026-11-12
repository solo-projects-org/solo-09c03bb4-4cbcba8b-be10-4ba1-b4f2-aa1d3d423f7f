"""Windowed GeoTIFF validation and parcel cloud-free coverage analysis."""

from __future__ import annotations

import math
from contextlib import ExitStack
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.windows import Window

# Read the parcel map in bounded tiles.  A GeoTIFF can legally have one very
# large strip; splitting natural blocks keeps memory use independent of the
# raster size.  This can be overridden for unusually small containers.
_DEFAULT_BLOCK_SIZE = 1024
_GRID_EPS = 1e-6  # pixels
_RES_EPS_REL = 1e-9
_RES_EPS_ABS = 1e-9


@dataclass(frozen=True)
class UploadedRaster:
    """An uploaded file that has been copied to a private temporary directory."""

    path: str
    filename: str
    role: str  # "parcel" or "mask"


@dataclass(frozen=True)
class RasterMeta:
    upload: UploadedRaster
    driver: str
    width: int
    height: int
    band_count: int
    crs: CRS | None
    transform: Affine
    nodata: float | int | None
    dtype: str

    @property
    def filename(self) -> str:
        return self.upload.filename


@dataclass(frozen=True)
class IncompatibleFile:
    filename: str
    role: str
    reason: str

    def as_dict(self) -> dict[str, str]:
        return {
            "filename": self.filename,
            "role": self.role,
            "reason": self.reason,
        }


@dataclass
class ParcelResult:
    parcel_id: int
    total_pixels: int
    covered_pixels: int

    @property
    def coverage_ratio(self) -> float:
        if self.total_pixels == 0:
            return 0.0
        return self.covered_pixels / self.total_pixels

    def as_dict(self, target_coverage: float) -> dict[str, int | float]:
        return {
            "parcel_id": self.parcel_id,
            "total_pixels": self.total_pixels,
            "covered_pixels": self.covered_pixels,
            "coverage_ratio": self.coverage_ratio,
            "below_target": self.coverage_ratio < target_coverage,
        }


def _fmt_num(value: float) -> str:
    return f"{value:g}" if math.isfinite(value) else str(value)


def _display_crs(crs: CRS | None) -> str:
    if crs is None:
        return "无"
    return crs.to_string()


def _inspect(upload: UploadedRaster) -> tuple[RasterMeta | None, IncompatibleFile | None]:
    try:
        with rasterio.open(upload.path) as ds:
            meta = RasterMeta(
                upload=upload,
                driver=ds.driver,
                width=ds.width,
                height=ds.height,
                band_count=ds.count,
                crs=ds.crs,
                transform=ds.transform,
                nodata=ds.nodata,
                dtype=ds.dtypes[0] if ds.count else "unknown",
            )
    except Exception as exc:  # rasterio raises several driver-specific errors
        message = str(exc).replace(upload.path, upload.filename)
        return None, IncompatibleFile(
            filename=upload.filename,
            role=upload.role,
            reason=f"无法作为栅格影像打开：{message}",
        )

    reason: str | None = None
    if meta.driver.lower() != "gtiff":
        reason = f"格式不是 GeoTIFF（driver={meta.driver}）"
    elif meta.band_count != 1:
        reason = f"必须为单波段，实际有 {meta.band_count} 个波段"
    elif meta.width <= 0 or meta.height <= 0:
        reason = "影像尺寸无效"
    elif meta.crs is None:
        reason = "缺少坐标系（CRS）"
    else:
        t = meta.transform
        if abs(t.b) > 1e-12 or abs(t.d) > 1e-12:
            reason = "影像不是北向上（仿射变换含旋转/剪切）"
        elif t.a <= 0 or t.e >= 0:
            reason = "影像不是北向上（像元宽应为正、高应为负）"
        elif not math.isfinite(t.a) or not math.isfinite(t.e):
            reason = "像元分辨率无效"

    if reason is not None:
        return meta, IncompatibleFile(
            filename=upload.filename, role=upload.role, reason=reason
        )
    return meta, None


def _pixel_offset(delta: float, resolution: float) -> int | None:
    quotient = delta / resolution
    nearest = round(quotient)
    if abs(quotient - nearest) > _GRID_EPS:
        return None
    return int(nearest)


def _same_resolution(a: float, b: float) -> bool:
    return math.isclose(
        abs(a), abs(b), rel_tol=_RES_EPS_REL, abs_tol=_RES_EPS_ABS
    )


def validate_uploads(
    parcel: UploadedRaster, masks: Sequence[UploadedRaster]
) -> tuple[RasterMeta | None, list[RasterMeta], list[IncompatibleFile]]:
    """Validate file structure and that every mask lies on the parcel grid.

    Returns the parcel metadata, mask metadata, and all incompatibilities.  A
    non-empty error list means the request must be rejected: returning partial
    coverage would silently understate coverage for the rejected scene(s).
    """

    errors: list[IncompatibleFile] = []
    parcel_meta, parcel_error = _inspect(parcel)
    if parcel_error is not None:
        errors.append(parcel_error)

    # The parcel map is the reference grid.  If it is structurally unusable we
    # can still report every mask's own structural problems.
    reference = parcel_meta if parcel_error is None else None

    mask_metas: list[RasterMeta] = []
    for mask in masks:
        meta, own_error = _inspect(mask)
        if own_error is not None:
            errors.append(own_error)
        if meta is None:
            continue

        if reference is not None:
            reasons: list[str] = []
            if meta.crs is None or reference.crs is None:
                reasons.append("坐标系缺失")
            elif meta.crs != reference.crs:
                reasons.append(
                    "坐标系与地块图不一致"
                    f"（{_display_crs(meta.crs)} vs 地块图 {_display_crs(reference.crs)}）"
                )

            if not _same_resolution(meta.transform.a, reference.transform.a) or not (
                _same_resolution(meta.transform.e, reference.transform.e)
            ):
                reasons.append(
                    "分辨率与地块图不一致"
                    f"（{_fmt_num(abs(meta.transform.a))}×{_fmt_num(abs(meta.transform.e))}"
                    f" vs {_fmt_num(abs(reference.transform.a))}"
                    f"×{_fmt_num(abs(reference.transform.e))}）"
                )
            else:
                # Origins may differ, but only by an integral number of pixels
                # in each direction so the pixel grids stay coincident.
                ox = _pixel_offset(
                    reference.transform.c - meta.transform.c, reference.transform.a
                )
                oy = _pixel_offset(
                    reference.transform.f - meta.transform.f, reference.transform.e
                )
                if ox is None or oy is None:
                    qx = (
                        reference.transform.c - meta.transform.c
                    ) / reference.transform.a
                    qy = (
                        reference.transform.f - meta.transform.f
                    ) / reference.transform.e
                    reasons.append(
                        "原点与地块图未对齐到整数像元"
                        f"（偏移 {qx:.6f}, {qy:.6f} 像元）"
                    )

            # A structural problem in the mask itself is reported even if
            # comparisons with the reference can still be attempted.
            if own_error is None and reasons:
                errors.append(
                    IncompatibleFile(
                        filename=mask.filename,
                        role=mask.role,
                        reason="；".join(reasons),
                    )
                )

        mask_metas.append(meta)

    if errors:
        return None, [], errors
    return parcel_meta, mask_metas, []


def _split_window(window: Window, max_size: int) -> Iterable[Window]:
    width = int(window.width)
    height = int(window.height)
    col_off = int(window.col_off)
    row_off = int(window.row_off)
    for row in range(0, height, max_size):
        tile_height = min(max_size, height - row)
        for col in range(0, width, max_size):
            tile_width = min(max_size, width - col)
            yield Window(
                col_off=col_off + col,
                row_off=row_off + row,
                width=tile_width,
                height=tile_height,
            )


def _valid_data(values: np.ndarray, nodata: float | int | None) -> np.ndarray:
    """Boolean valid-data mask, treating the dataset nodata value as invalid."""

    if nodata is None:
        return np.ones(values.shape, dtype=bool)
    if np.issubdtype(values.dtype, np.floating) and math.isnan(float(nodata)):
        return ~np.isnan(values)
    return values != nodata


def _parcel_ids(
    values: np.ndarray, nodata: float | int | None
) -> tuple[np.ndarray, np.ndarray]:
    """Return int64 parcel IDs and a validity mask for positive integers only."""

    valid = _valid_data(values, nodata)
    if np.issubdtype(values.dtype, np.floating):
        valid &= np.isfinite(values)
        valid &= values > 0
        valid &= values == np.floor(values)
        # Avoid casting NaN/infs: replace invalid entries before the cast below.
        safe = np.where(valid, values, 0)
        ids = safe.astype(np.int64, copy=False)
        return ids, valid

    if np.issubdtype(values.dtype, np.unsignedinteger):
        valid &= values > 0
    else:
        valid &= values > 0
    ids = values.astype(np.int64, copy=False)
    return ids, valid


def _clip_to_scene(
    tile: Window, offset_x: int, offset_y: int, scene_width: int, scene_height: int
) -> tuple[int, int, int, int, int, int] | None:
    """Map a parcel tile into a scene and return intersection geometry.

    Result is (scene_col, scene_row, parcel_dx, parcel_dy, width, height) in
    local tile coordinates.  ``None`` means no intersection.
    """

    tile_col = int(tile.col_off)
    tile_row = int(tile.row_off)
    tile_width = int(tile.width)
    tile_height = int(tile.height)

    scene_col = tile_col + offset_x
    scene_row = tile_row + offset_y
    dest_col = max(0, -scene_col)
    dest_row = max(0, -scene_row)
    src_col = max(0, scene_col)
    src_row = max(0, scene_row)

    read_width = min(tile_width - dest_col, scene_width - src_col)
    read_height = min(tile_height - dest_row, scene_height - src_row)
    if read_width <= 0 or read_height <= 0:
        return None
    return src_col, src_row, dest_col, dest_row, read_width, read_height


def analyze_coverage(
    parcel_meta: RasterMeta,
    mask_metas: Sequence[RasterMeta],
    block_size: int = _DEFAULT_BLOCK_SIZE,
) -> list[ParcelResult]:
    """Count parcel pixels and cloud-free coverage without loading whole rasters."""

    totals: dict[int, int] = {}
    covered: dict[int, int] = {}

    with ExitStack() as stack:
        parcel_ds = stack.enter_context(rasterio.open(parcel_meta.upload.path))
        mask_datasets = [
            stack.enter_context(rasterio.open(m.upload.path)) for m in mask_metas
        ]
        offsets = [
            (
                int(
                    round(
                        (parcel_ds.transform.c - ds.transform.c)
                        / parcel_ds.transform.a
                    )
                ),
                int(
                    round(
                        (parcel_ds.transform.f - ds.transform.f)
                        / parcel_ds.transform.e
                    )
                ),
            )
            for ds in mask_datasets
        ]

        # Iterate over the parcel map's own blocks, splitting exceptionally
        # large strips/tiles so no more than block_size² pixels are resident.
        for _block_index, block_window in parcel_ds.block_windows(1):
            for tile in _split_window(block_window, max(1, block_size)):
                parcel_values = parcel_ds.read(1, window=tile)
                ids, valid = _parcel_ids(parcel_values, parcel_meta.nodata)
                if not valid.any():
                    continue
                coverage = np.zeros(parcel_values.shape, dtype=bool)

                for mask_ds, mask_meta, (offset_x, offset_y) in zip(
                    mask_datasets, mask_metas, offsets
                ):
                    geom = _clip_to_scene(
                        tile,
                        offset_x,
                        offset_y,
                        mask_ds.width,
                        mask_ds.height,
                    )
                    if geom is None:
                        continue
                    (
                        src_col,
                        src_row,
                        dest_col,
                        dest_row,
                        read_width,
                        read_height,
                    ) = geom
                    values = mask_ds.read(
                        1,
                        window=Window(
                            src_col, src_row, read_width, read_height
                        ),
                    )
                    mask_valid = _valid_data(values, mask_meta.nodata)
                    # A pixel is covered when any scene's valid mask equals 1.
                    # The boolean OR makes overlapping clear scenes count once.
                    coverage[
                        dest_row : dest_row + read_height,
                        dest_col : dest_col + read_width,
                    ] |= (values == 1) & mask_valid

                flat_ids = ids[valid]
                if flat_ids.size == 0:
                    continue

                parcel_ids, counts = np.unique(flat_ids, return_counts=True)
                for parcel_id, count in zip(parcel_ids.tolist(), counts.tolist()):
                    totals[int(parcel_id)] = totals.get(int(parcel_id), 0) + int(
                        count
                    )

                if coverage.any():
                    covered_ids = flat_ids[coverage[valid]]
                    parcel_ids, counts = np.unique(
                        covered_ids, return_counts=True
                    )
                    for parcel_id, count in zip(
                        parcel_ids.tolist(), counts.tolist()
                    ):
                        covered[int(parcel_id)] = covered.get(
                            int(parcel_id), 0
                        ) + int(count)

    results: list[ParcelResult] = []
    for parcel_id in sorted(totals):
        total_pixels = totals[parcel_id]
        covered_pixels = covered.get(parcel_id, 0)
        results.append(
            ParcelResult(
                parcel_id=parcel_id,
                total_pixels=total_pixels,
                covered_pixels=covered_pixels,
            )
        )
    return results
