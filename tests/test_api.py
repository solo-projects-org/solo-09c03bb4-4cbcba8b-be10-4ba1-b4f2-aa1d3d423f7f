"""End-to-end API tests.  Run with ``pytest -q``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
from fastapi.testclient import TestClient
from rasterio.transform import from_origin

from app.main import app

SAMPLES = Path(__file__).resolve().parent.parent / "samples"


@pytest.fixture(scope="module")
def client() -> TestClient:
    return TestClient(app)


def _post(client: TestClient, parcel: Path, masks: list[Path], target: float = 0.9):
    files = [("parcel", (parcel.name, parcel.open("rb"), "image/tiff"))]
    files += [("masks", (m.name, m.open("rb"), "image/tiff")) for m in masks]
    try:
        return client.post(
            "/analyze", files=files, data={"target_coverage": str(target)}
        )
    finally:
        for _, (_, handle, _) in files:
            handle.close()


def _write(path: Path, values: np.ndarray, crs: str, transform) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=values.shape[1],
        height=values.shape[0],
        count=1,
        dtype=values.dtype,
        crs=crs,
        transform=transform,
        nodata=0,
    ) as dst:
        dst.write(values, 1)


def test_sample_request(client: TestClient):
    response = _post(
        client,
        SAMPLES / "parcel_map.tif",
        [
            SAMPLES / "scene_01_cloud_free.tif",
            SAMPLES / "scene_02_cloud_free.tif",
        ],
    )
    assert response.status_code == 200
    body = response.json()
    assert body["total_parcels_analyzed"] == 4
    assert body["below_target_count"] == 2
    assert [p["parcel_id"] for p in body["parcels"]] == [102, 103]
    assert body["parcels"][0]["total_pixels"] == 2200
    assert body["parcels"][0]["covered_pixels"] == 600
    assert body["parcels"][1]["total_pixels"] == 2750
    assert body["parcels"][1]["covered_pixels"] == 500


def test_overlapping_masks_count_once(client: TestClient, tmp_path: Path):
    values = np.ones((10, 10), dtype=np.uint16)
    _write(
        tmp_path / "parcel.tif",
        values,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )

    mask_a = np.zeros((10, 10), dtype=np.uint8)
    mask_a[:, :7] = 1
    _write(
        tmp_path / "a.tif",
        mask_a,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )

    mask_b = np.zeros((10, 10), dtype=np.uint8)
    mask_b[:, 3:] = 1
    _write(
        tmp_path / "b.tif",
        mask_b,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )

    # 7 + 7 columns overlap in 4, so the union is all 10 columns.
    response = _post(client, tmp_path / "parcel.tif", [tmp_path / "a.tif", tmp_path / "b.tif"])
    assert response.status_code == 200
    assert response.json()["below_target_count"] == 0


def test_partial_coverage_is_below_target(client: TestClient, tmp_path: Path):
    parcel = np.full((10, 10), 7, dtype=np.uint16)
    _write(
        tmp_path / "parcel.tif",
        parcel,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )
    mask = np.zeros((10, 10), dtype=np.uint8)
    mask[:5, :] = 1
    _write(
        tmp_path / "mask.tif",
        mask,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )

    response = _post(client, tmp_path / "parcel.tif", [tmp_path / "mask.tif"])
    assert response.status_code == 200
    row = response.json()["parcels"][0]
    assert row == {
        "parcel_id": 7,
        "total_pixels": 100,
        "covered_pixels": 50,
        "coverage_ratio": 0.5,
        "below_target": True,
    }


def test_incompatible_crs_reports_file(client: TestClient, tmp_path: Path):
    parcel = np.ones((10, 10), dtype=np.uint16)
    _write(
        tmp_path / "parcel.tif",
        parcel,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )
    mask = np.ones((10, 10), dtype=np.uint8)
    _write(
        tmp_path / "wrong_crs.tif",
        mask,
        "EPSG:4326",
        from_origin(110, 30, 0.001, 0.001),
    )

    response = _post(client, tmp_path / "parcel.tif", [tmp_path / "wrong_crs.tif"])
    assert response.status_code == 400
    errors = response.json()["detail"]["incompatible_files"]
    assert errors[0]["filename"] == "wrong_crs.tif"
    assert "坐标系" in errors[0]["reason"]


def test_non_integer_origin_reports_file(client: TestClient, tmp_path: Path):
    parcel = np.ones((10, 10), dtype=np.uint16)
    _write(
        tmp_path / "parcel.tif",
        parcel,
        "EPSG:32650",
        from_origin(500000, 3000000, 10, 10),
    )
    mask = np.ones((10, 10), dtype=np.uint8)
    _write(
        tmp_path / "shifted.tif",
        mask,
        "EPSG:32650",
        from_origin(500005, 3000000, 10, 10),
    )

    response = _post(client, tmp_path / "parcel.tif", [tmp_path / "shifted.tif"])
    assert response.status_code == 400
    errors = response.json()["detail"]["incompatible_files"]
    assert errors[0]["filename"] == "shifted.tif"
    assert "整数像元" in errors[0]["reason"]
