"""End-to-end tests for POST /update-parcel-boundaries.  Run with ``pytest -q``."""

from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pytest
import rasterio
from fastapi.testclient import TestClient
from rasterio.transform import from_origin
from rasterio.warp import transform_geom

from app.main import app

SAMPLES = Path(__file__).resolve().parent.parent / "samples"


@pytest.fixture(scope="module")
def client() -> TestClient:
    return TestClient(app)


def _write_parcel(path: Path, values: np.ndarray, crs: str = "EPSG:32650",
                  origin=(500_000.0, 3_000_000.0), pixel: float = 10.0,
                  dtype: str | None = None) -> None:
    if dtype is not None:
        values = values.astype(dtype)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=values.shape[1],
        height=values.shape[0],
        count=1,
        dtype=values.dtype,
        crs=crs,
        transform=from_origin(origin[0], origin[1], pixel, pixel),
        nodata=0,
        tiled=True,
        blockxsize=16,
        blockysize=16,
        compress="deflate",
    ) as dst:
        dst.write(values, 1)


def _polygon_feature(parcel_id: int, utm_ring: list[tuple[float, float]],
                     crs: str = "EPSG:32650") -> dict:
    geom = {
        "type": "Polygon",
        "coordinates": [[[x, y] for x, y in utm_ring]],
    }
    projected = transform_geom(crs, "EPSG:4326", geom)
    ring = [[float(lon), float(lat)] for lon, lat in projected["coordinates"][0]]
    ring.append(ring[0])
    return {
        "type": "Feature",
        "properties": {"parcel_id": parcel_id},
        "geometry": {"type": "Polygon", "coordinates": [ring]},
    }


def _post(client: TestClient, parcel: Path, geojson: dict | bytes):
    if isinstance(geojson, dict):
        payload = json.dumps(geojson).encode("utf-8")
    else:
        payload = geojson
    with parcel.open("rb") as handle:
        return client.post(
            "/update-parcel-boundaries",
            files={
                "parcel": (parcel.name, handle, "image/tiff"),
                "boundaries": ("boundaries.geojson", io.BytesIO(payload),
                               "application/geo+json"),
            },
        )


def _read_response_tiff(response) -> tuple[np.ndarray, dict]:
    with rasterio.open(io.BytesIO(response.content)) as ds:
        return ds.read(1), {
            "dtype": ds.dtypes[0],
            "crs": ds.crs,
            "transform": ds.transform,
            "width": ds.width,
            "height": ds.height,
            "count": ds.count,
            "nodata": ds.nodata,
        }


def _square_ring(x0: float, y0: float, x1: float, y1: float) -> list[tuple[float, float]]:
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def test_sample_boundary_swap(client: TestClient, tmp_path: Path):
    geojson = json.loads((SAMPLES / "boundary_update_sample.geojson").read_text())
    response = _post(client, SAMPLES / "parcel_map.tif", geojson)
    assert response.status_code == 200, response.text
    assert response.headers["content-type"] == "image/tiff"

    data, meta = _read_response_tiff(response)
    with rasterio.open(SAMPLES / "parcel_map.tif") as src:
        assert (meta["width"], meta["height"]) == (src.width, src.height)
        assert meta["crs"] == src.crs
        assert meta["transform"] == src.transform
        assert meta["dtype"] == src.dtypes[0]
        assert meta["nodata"] == 0 and meta["count"] == 1

    # Parcel 101: rows 0-19 cols 0-59 (1200 px) + rows 20-39 cols 0-44 (900 px).
    assert data[0:20, 0:60].min() == 101 and np.unique(data[0:20, 0:60]).tolist() == [101]
    assert data[20:40, 0:45].min() == 101 and np.unique(data[20:40, 0:45]).tolist() == [101]
    # Parcel 102: rows 0-19 cols 60-109 and rows 20-39 cols 45-109.
    assert np.unique(data[0:20, 60:110]).tolist() == [102]
    assert np.unique(data[20:40, 45:110]).tolist() == [102]
    # Each parcel occupies part of the other's former territory.
    with rasterio.open(SAMPLES / "parcel_map.tif") as original:
        base = original.read(1)
    assert ((data == 101) & (base == 102)).sum() == 100
    assert ((data == 102) & (base == 101)).sum() == 100
    # Old pixels of the corrected parcels are fully reallocated; none left
    # dangling under the old ID.
    assert ((base == 101) & (data == 0)).sum() == 0
    assert ((base == 102) & (data == 0)).sum() == 0
    # Untouched parcels remain exactly as they were; the inter-row gap and
    # the blank far-right columns are still 0.
    assert np.array_equal(data[45:100, :], base[45:100, :])
    assert np.array_equal(data[40:45, :], base[40:45, :])
    assert (data[:, 110:120] == 0).all()
    assert (data == 101).sum() == 2100
    assert (data == 102).sum() == 2300


def test_untouched_parcels_preserved(client: TestClient, tmp_path: Path):
    values = np.zeros((12, 12), dtype=np.uint16)
    values[0:6, 0:6] = 1
    values[0:6, 6:12] = 2
    values[6:12, 0:6] = 3
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    # Shrink parcel 1 to the top-left 3x3 block.
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, _square_ring(500_000.0, 3_000_000.0,
                                             500_030.0, 2_999_970.0))
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 200, response.text
    data, _ = _read_response_tiff(response)
    assert np.unique(data[0:3, 0:3]).tolist() == [1]
    assert (data == 1).sum() == 9
    # The rest of parcel 1's old pixels are released to 0.
    assert data[3:6, 0:6].sum() == 0
    # Parcels 2 and 3 are untouched.
    assert np.unique(data[0:6, 6:12]).tolist() == [2]
    assert np.unique(data[6:12, 0:6]).tolist() == [3]


def test_new_new_conflict_reports_ids_and_pixel(client: TestClient, tmp_path: Path):
    values = np.zeros((10, 10), dtype=np.uint16)
    values[:, :5] = 1
    values[:, 5:] = 2
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    geojson = {
        "type": "FeatureCollection",
        "features": [
            # Both new polygons cover columns 2-4 -> overlapping pixels.
            _polygon_feature(1, _square_ring(500_000.0, 3_000_000.0,
                                             500_050.0, 2_999_900.0)),
            _polygon_feature(2, _square_ring(500_020.0, 3_000_000.0,
                                             500_100.0, 2_999_900.0)),
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    detail = response.json()["detail"]
    assert len(detail["conflicts"]) == 1
    conflict = detail["conflicts"][0]
    assert conflict["type"] == "new_new"
    assert conflict["parcel_ids"] == [1, 2]
    # First overlapping pixel center is column 2, row 0.
    assert conflict["column"] == 2 and conflict["row"] == 0
    assert detail["unknown_parcel_ids"] == []
    assert detail["no_hit_parcel_ids"] == []


def test_new_old_conflict_with_untouched_parcel(client: TestClient, tmp_path: Path):
    values = np.zeros((10, 10), dtype=np.uint16)
    values[:5, :] = 1
    values[5:, :] = 2
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    # New boundary of parcel 1 reaches two rows into parcel 2, which is not
    # part of the batch.
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, _square_ring(500_000.0, 3_000_000.0,
                                             500_100.0, 2_999_930.0))
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    conflict = response.json()["detail"]["conflicts"][0]
    assert conflict["type"] == "new_old"
    assert conflict["parcel_ids"] == [1, 2]
    assert conflict["row"] == 5 and conflict["column"] == 0


def test_reusing_batch_parcel_old_pixels_is_allowed(client: TestClient, tmp_path: Path):
    values = np.zeros((10, 10), dtype=np.uint16)
    values[:5, :5] = 1
    values[:5, 5:] = 2
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    # Both parcels are corrected as simple rectangles whose pixel sets are
    # disjoint but each reaches into the other's former territory.  Parcel 1
    # keeps rows 0-4 across its old columns plus one column of parcel 2's old
    # land (col 5); parcel 2 keeps the remaining old-2 columns (6-9) and drops
    # one row band into parcel 1's old rows below.  Because both IDs are in the
    # batch, reusing each other's old pixels is allowed and there is no
    # new-new conflict.
    ring_1 = [
        (500_000.0, 3_000_000.0),
        (500_060.0, 3_000_000.0),
        (500_060.0, 2_999_950.0),
        (500_000.0, 2_999_950.0),
    ]
    ring_2 = [
        (500_060.0, 2_999_950.0),
        (500_100.0, 2_999_950.0),
        (500_100.0, 2_999_900.0),
        (500_060.0, 2_999_900.0),
    ]
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, ring_1),
            _polygon_feature(2, ring_2),
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 200, response.text
    data, _ = _read_response_tiff(response)
    # Parcel 1: rows 0-4, cols 0-5 (includes old parcel-2 column 5).
    assert np.unique(data[0:5, 0:6]).tolist() == [1]
    assert (data == 1).sum() == 30
    # Parcel 2: rows 5-9, cols 6-9 (its old rows 0-4 are released).
    assert np.unique(data[5:10, 6:10]).tolist() == [2]
    assert (data == 2).sum() == 20
    # Everything else is released to 0.
    assert data.sum() == 30 * 1 + 20 * 2
    assert (data[0:5, 6:10] == 0).all()
    assert (data[5:10, 0:6] == 0).all()


def test_unknown_parcel_id_is_rejected(client: TestClient, tmp_path: Path):
    values = np.ones((8, 8), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(99, _square_ring(500_000.0, 3_000_000.0,
                                              500_040.0, 2_999_960.0))
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    detail = response.json()["detail"]
    assert detail["unknown_parcel_ids"] == [99]


def test_duplicate_parcel_ids_rejected(client: TestClient, tmp_path: Path):
    values = np.ones((8, 8), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    ring = _square_ring(500_000.0, 3_000_000.0, 500_040.0, 2_999_960.0)
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, ring),
            _polygon_feature(1, ring),
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["geojson_errors"]
    assert any("重复" in e for e in errors)


def test_geometry_type_and_duplicate_ids_rejected(client: TestClient, tmp_path: Path):
    values = np.ones((8, 8), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"parcel_id": 1},
                "geometry": {"type": "MultiPolygon", "coordinates": []},
            },
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["geojson_errors"]
    assert any("Polygon" in e for e in errors)


def test_polygon_without_hits_is_rejected(client: TestClient, tmp_path: Path):
    values = np.ones((8, 8), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    # Polygon sits entirely beyond the raster extent (raster spans 100 m).
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, _square_ring(500_500.0, 3_000_000.0,
                                             500_540.0, 2_999_960.0))
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    assert response.json()["detail"]["no_hit_parcel_ids"] == [1]


def test_more_than_32_features_rejected(client: TestClient, tmp_path: Path):
    values = np.ones((40, 40), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    features = [
        _polygon_feature(
            i + 1,
            _square_ring(500_000.0, 3_000_000.0, 500_010.0, 2_999_990.0),
        )
        for i in range(33)
    ]
    # Make IDs unique even though only id 1 exists on the map.
    for i, feature in enumerate(features):
        feature["properties"]["parcel_id"] = i + 1
    response = _post(client, parcel, {"type": "FeatureCollection", "features": features})
    assert response.status_code == 400
    assert any("32" in e for e in response.json()["detail"]["geojson_errors"])


def test_wrong_dtype_and_nodata_rejected(client: TestClient, tmp_path: Path):
    int32_path = tmp_path / "int32.tif"
    _write_parcel(int32_path, np.ones((6, 6), dtype=np.int32), dtype="int32")
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, _square_ring(500_000.0, 3_000_000.0,
                                             500_030.0, 2_999_970.0))
        ],
    }
    response = _post(client, int32_path, geojson)
    assert response.status_code == 400
    assert "UInt16" in response.text or "UInt32" in response.text

    no_nodata_path = tmp_path / "no_nodata.tif"
    with rasterio.open(
        no_nodata_path, "w", driver="GTiff", width=6, height=6, count=1,
        dtype="uint16", crs="EPSG:32650",
        transform=from_origin(500_000, 3_000_000, 10, 10),
    ) as dst:
        dst.write(np.ones((6, 6), dtype=np.uint16), 1)
    response = _post(client, no_nodata_path, geojson)
    assert response.status_code == 400
    assert "NoData" in response.text


def test_pixel_center_on_edge_counts_as_hit(client: TestClient, tmp_path: Path):
    # A 10 m raster; polygon edges pass exactly through pixel centers
    # (half-pixel offsets).  Centers on the edge must count as hits.
    values = np.zeros((8, 8), dtype=np.uint16)
    values[1:6, 1:6] = 7
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    # Edges run through the centers of the boundary pixel rows/columns:
    # pixel (r, c) center maps to (500000 + 10c + 5, 3000000 - 10r - 5).
    # The box [1.5, 5.5) in pixel space covers exactly rows/cols 2..5 by
    # center inclusion, with every boundary center sitting on an edge.
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(7, _square_ring(500_015.0, 2_999_985.0,
                                             500_055.0, 2_999_945.0))
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 200, response.text
    data, _ = _read_response_tiff(response)
    # Pixels whose centers are inside or on the polygon: rows/cols 2-5.
    assert np.unique(data[2:6, 2:6]).tolist() == [7]
    assert (data == 7).sum() == 16
    # The neighboring row/column centers are strictly outside.
    assert data[1, 2:6].sum() == 0 and data[6, 2:6].sum() == 0
    assert data[2:6, 1].sum() == 0 and data[2:6, 6].sum() == 0


def test_wgs84_geojson_reprojected_onto_utm_grid(client: TestClient, tmp_path: Path):
    # Build the polygon directly in lon/lat around the UTM raster area.
    values = np.zeros((10, 10), dtype=np.uint32)
    values[:, :] = 5
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    utm_geom = {
        "type": "Polygon",
        "coordinates": [[
            [500_000.0, 3_000_000.0],
            [500_040.0, 3_000_000.0],
            [500_040.0, 2_999_960.0],
            [500_000.0, 2_999_960.0],
        ]],
    }
    projected = transform_geom("EPSG:32650", "EPSG:4326", utm_geom)
    ring = [[float(lon), float(lat)] for lon, lat in projected["coordinates"][0]]
    ring.append(ring[0])
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"parcel_id": 5},
                "geometry": {"type": "Polygon", "coordinates": [ring]},
            }
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 200, response.text
    data, meta = _read_response_tiff(response)
    assert meta["dtype"] == "uint32"
    assert np.unique(data[0:4, 0:4]).tolist() == [5]
    assert (data == 5).sum() == 16
    assert data[4:, :].sum() == 0


def test_conflict_pixel_coordinates_use_full_raster_origin(
    client: TestClient, tmp_path: Path
):
    # The raster is large enough that the conflict falls inside a block whose
    # origin is not (0, 0); the reported row/column must still be absolute.
    values = np.zeros((40, 40), dtype=np.uint16)
    values[20:30, 0:20] = 1
    values[20:30, 20:40] = 2
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    def ring(x0, y0, x1, y1):
        return _square_ring(
            500_000.0 + x0 * 10.0,
            3_000_000.0 - y0 * 10.0,
            500_000.0 + x1 * 10.0,
            3_000_000.0 - y1 * 10.0,
        )

    geojson = {
        "type": "FeatureCollection",
        "features": [
            # Both cover pixel (row 25, col 18).
            _polygon_feature(1, ring(10, 20, 25, 30)),
            _polygon_feature(2, ring(18, 20, 30, 30)),
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    conflict = response.json()["detail"]["conflicts"][0]
    assert conflict["type"] == "new_new"
    assert conflict["parcel_ids"] == [1, 2]
    assert conflict["row"] == 20 and conflict["column"] == 18


def test_polygon_holes_are_excluded(client: TestClient, tmp_path: Path):
    values = np.ones((10, 10), dtype=np.uint16)
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)

    outer = _square_ring(500_000.0, 3_000_000.0, 500_100.0, 2_999_900.0)
    hole = _square_ring(500_030.0, 2_999_970.0, 500_070.0, 2_999_930.0)
    projected_outer = transform_geom(
        "EPSG:32650",
        "EPSG:4326",
        {"type": "Polygon", "coordinates": [outer + [outer[0]]]},
    )
    projected_hole = transform_geom(
        "EPSG:32650",
        "EPSG:4326",
        {"type": "Polygon", "coordinates": [hole + [hole[0]]]},
    )
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"parcel_id": 1},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        projected_outer["coordinates"][0],
                        projected_hole["coordinates"][0],
                    ],
                },
            }
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 200, response.text
    data, _ = _read_response_tiff(response)
    # 10x10 exterior minus the 4x4 hole (pixel rows/cols 3-6) -> 84 pixels.
    assert data[4, 4] == 0
    assert data[1, 1] == 1
    assert (data == 1).sum() == 84


def test_invalid_geojson_payload_rejected(client: TestClient, tmp_path: Path):
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, np.ones((4, 4), dtype=np.uint16))
    response = _post(client, parcel, b"{not valid json")
    assert response.status_code == 400
    assert "GeoJSON" in response.text


def test_no_output_file_when_batch_fails(client: TestClient, tmp_path: Path):
    # Conflict case must not produce a raster body.
    values = np.zeros((6, 6), dtype=np.uint16)
    values[:3, :] = 1
    values[3:, :] = 2
    parcel = tmp_path / "parcel.tif"
    _write_parcel(parcel, values)
    geojson = {
        "type": "FeatureCollection",
        "features": [
            _polygon_feature(1, _square_ring(500_000.0, 3_000_000.0,
                                             500_060.0, 2_999_940.0)),
            _polygon_feature(2, _square_ring(500_020.0, 3_000_000.0,
                                             500_060.0, 2_999_940.0)),
        ],
    }
    response = _post(client, parcel, geojson)
    assert response.status_code == 400
    assert response.headers["content-type"] != "image/tiff"
    with pytest.raises(Exception):
        rasterio.open(io.BytesIO(response.content))
