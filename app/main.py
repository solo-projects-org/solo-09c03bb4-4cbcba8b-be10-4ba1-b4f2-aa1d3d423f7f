"""FastAPI application for finding parcels with insufficient cloud-free coverage."""

from __future__ import annotations

import math
import shutil
import tempfile
from pathlib import Path
from typing import Annotated

from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from app.analysis import (
    IncompatibleFile,
    UploadedRaster,
    analyze_coverage,
    validate_uploads,
)
from app.boundaries import (
    BoundaryError,
    apply_boundary_updates,
    open_parcel_raster,
    parse_boundary_geojson,
)

app = FastAPI(
    title="农情遥感补采分析 API",
    version="1.0.0",
    description=(
        "上传一张地块编号图和任一数量的单波段无云掩膜，统计每个地块的"
        "总像元数与至少一景无云覆盖的像元数，并返回覆盖率低于目标值的地块。"
    ),
)


class ParcelCoverage(BaseModel):
    parcel_id: int = Field(description="地块编号（地块图中的正整数）")
    total_pixels: int = Field(description="该地块在地块图中的总有效像元数")
    covered_pixels: int = Field(
        description="至少一景无云掩膜值为 1 的像元数；多景重叠只计一次"
    )
    coverage_ratio: float = Field(description="covered_pixels / total_pixels")
    below_target: bool = Field(description="覆盖率是否低于请求中的目标覆盖率")


class AnalysisResponse(BaseModel):
    target_coverage: float
    total_parcels_analyzed: int
    below_target_count: int
    parcels: list[ParcelCoverage] = Field(
        description="低于目标覆盖率的地块，按地块编号升序排列"
    )


def _safe_upload_name(upload: UploadFile, index: int) -> str:
    raw = upload.filename or f"上传文件_{index}"
    # Keep a readable label for error messages but never use it as a path.
    name = raw.replace("\\", "/").rsplit("/", 1)[-1].strip()
    return name or f"上传文件_{index}"


def _save_upload(upload: UploadFile, destination: Path) -> None:
    upload.file.seek(0)
    with destination.open("wb") as out:
        shutil.copyfileobj(upload.file, out, length=1024 * 1024)


def _incompatible_response(errors: list[IncompatibleFile]) -> HTTPException:
    return HTTPException(
        status_code=400,
        detail={
            "message": "输入影像不兼容，未执行分析。",
            "incompatible_files": [e.as_dict() for e in errors],
        },
    )


@app.get("/health", tags=["系统"], summary="健康检查")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/", include_in_schema=False)
def root() -> dict[str, str]:
    return {
        "name": "农情遥感补采分析 API",
        "docs": "/docs",
        "openapi": "/openapi.json",
    }


@app.post(
    "/analyze",
    response_model=AnalysisResponse,
    tags=["分析"],
    summary="找出无云观测覆盖不足的地块",
    responses={
        400: {
            "description": "存在不兼容的 GeoTIFF，响应中列出具体文件与原因",
        }
    },
)
def analyze(
    parcel: Annotated[
        UploadFile,
        File(
            description=(
                "地块编号图：单波段 GeoTIFF，北向上；有效地块编号为正整数"
            )
        ),
    ],
    masks: Annotated[
        list[UploadFile],
        File(
            description=(
                "各景无云掩膜：单波段 GeoTIFF，北向上，坐标系和分辨率与地块图"
                "相同；掩膜值 1 表示无云。可一次上传多景，景幅可错开重叠。"
            )
        ),
    ],
    target_coverage: Annotated[
        float,
        Form(
            description="目标覆盖率，范围 0 到 1；覆盖率严格低于该值的地块被返回"
        ),
    ],
) -> AnalysisResponse:
    if not math.isfinite(target_coverage) or not (0.0 <= target_coverage <= 1.0):
        raise HTTPException(
            status_code=400,
            detail="target_coverage 必须是 0 到 1 之间的有限数值。",
        )
    if not masks:
        raise HTTPException(status_code=400, detail="至少上传一景无云掩膜。")

    with tempfile.TemporaryDirectory(prefix="agri-coverage-") as tmp:
        tmp_dir = Path(tmp)
        parcels_dir = tmp_dir / "parcels"
        masks_dir = tmp_dir / "masks"
        parcels_dir.mkdir()
        masks_dir.mkdir()

        parcel_name = _safe_upload_name(parcel, 0)
        parcel_dest = parcels_dir / "parcel_0000.tif"
        _save_upload(parcel, parcel_dest)
        parcel_upload = UploadedRaster(
            path=str(parcel_dest), filename=parcel_name, role="parcel"
        )

        mask_uploads: list[UploadedRaster] = []
        for index, mask in enumerate(masks, start=1):
            mask_name = _safe_upload_name(mask, index)
            mask_dest = masks_dir / f"mask_{index:04d}.tif"
            _save_upload(mask, mask_dest)
            mask_uploads.append(
                UploadedRaster(
                    path=str(mask_dest), filename=mask_name, role="mask"
                )
            )

        parcel_meta, mask_metas, errors = validate_uploads(
            parcel_upload, mask_uploads
        )
        if errors:
            raise _incompatible_response(errors)

        assert parcel_meta is not None
        assert len(mask_metas) == len(mask_uploads)
        results = analyze_coverage(parcel_meta, mask_metas)

    below = [r for r in results if r.coverage_ratio < target_coverage]
    return AnalysisResponse(
        target_coverage=target_coverage,
        total_parcels_analyzed=len(results),
        below_target_count=len(below),
        parcels=[
            ParcelCoverage(
                parcel_id=r.parcel_id,
                total_pixels=r.total_pixels,
                covered_pixels=r.covered_pixels,
                coverage_ratio=r.coverage_ratio,
                below_target=True,
            )
            for r in below
        ],
    )


@app.post(
    "/update-parcel-boundaries",
    tags=["边界更正"],
    summary="按复测 GeoJSON 批量更正地块编号图边界",
    response_class=FileResponse,
    responses={
        200: {"content": {"image/tiff": {}}, "description": "更正后的地块编号图 GeoTIFF"},
        400: {"description": "整批核对未通过，不生成影像；detail 中列出编号、冲突像元行列号等"},
    },
)
def update_parcel_boundaries(
    background_tasks: BackgroundTasks,
    parcel: Annotated[
        UploadFile,
        File(
            description=(
                "底图地块编号图：北向上的单波段 UInt16/UInt32 GeoTIFF，NoData 为 0"
            )
        ),
    ],
    boundaries: Annotated[
        UploadFile,
        File(
            description=(
                "WGS84 GeoJSON（FeatureCollection）：每个原地块编号对应一个 Polygon，"
                "编号取 properties.parcel_id（或要素 id），最多 32 个要素"
            )
        ),
    ],
) -> FileResponse:
    temp_root = Path(tempfile.mkdtemp(prefix="agri-boundaries-"))
    background_tasks.add_task(shutil.rmtree, temp_root, ignore_errors=True)

    # Fixed names: an upload filename must never collide with another field
    # or be used as a path component.
    parcel_name = _safe_upload_name(parcel, 0)
    parcel_dest = temp_root / "parcel_map.tif"
    _save_upload(parcel, parcel_dest)

    geojson_dest = temp_root / "boundaries.geojson"
    _save_upload(boundaries, geojson_dest)

    try:
        spec = open_parcel_raster(str(parcel_dest))
        updates = parse_boundary_geojson(geojson_dest.read_bytes(), spec)
        output_path = temp_root / "parcel_map_updated.tif"
        apply_boundary_updates(spec, updates, str(output_path))
    except BoundaryError as exc:
        raise HTTPException(status_code=400, detail=exc.as_detail()) from exc

    return FileResponse(
        output_path,
        media_type="image/tiff",
        filename=f"{Path(parcel_name).stem}_updated.tif",
    )
