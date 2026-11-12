"""End-to-end tests for the /correct-boundaries endpoint.  Run with ``pytest -q``."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import rasterio
from fastapi.testclient import TestClient
from rasterio.transform import from_origin
from rasterio.warp import transform as warp_transform

from app.analysis import UploadedRaster
from app.boundaries import apply_correction, build_plan, inspect_parcel_raster, parse_geojson
from app.main import app

client = TestClient(app)

ORIGIN_X = 500_000.0
ORIGIN_Y = 3_000_000.0
PIXEL = 10.0
CRS = "EPSG:32650"


def _write_raster(path: Path, values: np.ndarray, *, nodata=0, crs=CRS) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=values.shape[1],
        height=values.shape[0],
        count=1,
        dtype=values.dtype,
        crs=crs,
        transform=from_origin(ORIGIN_X, ORIGIN_Y, PIXEL, PIXEL),
        nodata=nodata,
        tiled=True,
        blockxsize=64,
        blockysize=64,
    ) as dst:
        dst.write(values, 1)


def _rect_geojson(parcel_id: int, col0: int, row0: int, col1: int, row1: int) -> dict:
    """A WGS84 rectangle covering UTM pixel columns [col0, col1), rows [row0, row1)."""

    x0 = ORIGIN_X + col0 * PIXEL
    x1 = ORIGIN_X + col1 * PIXEL
    y_top = ORIGIN_Y - row0 * PIXEL
    y_bottom = ORIGIN_Y - row1 * PIXEL
    lon, lat = warp_transform(
        CRS,
        "EPSG:4326",
        [x0, x1, x1, x0, x0],
        [y_top, y_top, y_bottom, y_bottom, y_top],
    )
    return {
        "type": "Feature",
        "properties": {"parcel_id": parcel_id},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[round(x, 10), round(y, 10)] for x, y in zip(lon, lat)]],
        },
    }


def _collection(*features: dict) -> bytes:
    return json.dumps(
        {"type": "FeatureCollection", "features": list(features)},
        ensure_ascii=False,
    ).encode("utf-8")


def _post(parcel: Path, geojson: bytes):
    return client.post(
        "/correct-boundaries",
        files=[
            ("parcel", (parcel.name, parcel.open("rb"), "image/tiff")),
            ("boundaries", ("bounds.geojson", geojson, "application/geo+json")),
        ],
    )


def _base_raster(tmp_path: Path, dtype=np.uint16) -> Path:
    # Two side-by-side parcels on top, parcel 3 untouched at the bottom.
    values = np.zeros((10, 12), dtype=dtype)
    values[0:6, 0:6] = 1
    values[0:6, 6:12] = 2
    values[7:10, 0:12] = 3
    path = tmp_path / "parcel.tif"
    _write_raster(path, values)
    return path


def test_swap_reuses_old_cells_and_preserves_grid(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    # New 1 takes columns 3-8 (part of old 1 and part of old 2);
    # new 2 takes columns 0-2 (part of old 1).  No overlap between them.
    geojson = _collection(
        _rect_geojson(1, 3, 0, 9, 7),
        _rect_geojson(2, 0, 0, 3, 7),
    )

    response = _post(parcel, geojson)
    assert response.status_code == 200, response.text
    assert response.headers["content-type"] == "image/tiff"

    out = tmp_path / "out.tif"
    out.write_bytes(response.content)
    with rasterio.open(out) as ds:
        assert ds.width == 12 and ds.height == 10
        assert ds.dtypes[0] == "uint16"
        assert ds.nodata == 0
        assert ds.crs.to_epsg() == 32650
        assert ds.transform == from_origin(ORIGIN_X, ORIGIN_Y, PIXEL, PIXEL)
        result = ds.read(1)

    expected = np.zeros((10, 12), dtype=np.uint16)
    expected[0:7, 0:3] = 2
    expected[0:7, 3:9] = 1
    # Old cells of 1/2 not re-hit become NoData; parcel 3 is untouched.
    expected[7:10, 0:12] = 3
    np.testing.assert_array_equal(result, expected)


def test_new_new_conflict_reports_two_ids_and_pixel(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    geojson = _collection(
        _rect_geojson(1, 0, 0, 8, 6),
        _rect_geojson(2, 4, 0, 12, 6),
    )
    response = _post(parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["errors"]
    conflict = next(e for e in errors if e["kind"] == "new_new_conflict")
    assert conflict["parcel_ids"] == [1, 2]
    assert isinstance(conflict["row"], int)
    assert isinstance(conflict["column"], int)
    # The reported pixel really is inside both rectangles.
    assert 0 <= conflict["row"] < 6
    assert 4 <= conflict["column"] < 8


def test_new_old_conflict_reports_uncorrected_parcel(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    # New 1 reaches into parcel 3's rows; parcel 3 is not in the batch.
    geojson = _collection(_rect_geojson(1, 0, 0, 6, 9))
    response = _post(parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["errors"]
    conflict = next(e for e in errors if e["kind"] == "new_old_conflict")
    assert conflict["parcel_ids"] == [1, 3]
    assert conflict["row"] in (7, 8)
    assert 0 <= conflict["column"] < 6


def test_unknown_parcel_id_is_reported(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    geojson = _collection(_rect_geojson(99, 0, 0, 6, 6))
    response = _post(parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["errors"]
    err = next(e for e in errors if e["kind"] == "unknown_parcel_id")
    assert err["parcel_id"] == 99


def test_polygon_outside_map_hits_no_pixels(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    geojson = _collection(_rect_geojson(1, 40, 40, 46, 46))
    response = _post(parcel, geojson)
    assert response.status_code == 400
    errors = response.json()["detail"]["errors"]
    assert any(e["kind"] == "no_pixels_hit" and e["parcel_id"] == 1 for e in errors)


def test_invalid_dtype_rejected(tmp_path: Path):
    values = np.ones((5, 5), dtype=np.uint8)
    path = tmp_path / "u8.tif"
    _write_raster(path, values)
    response = _post(path, _collection(_rect_geojson(1, 0, 0, 5, 5)))
    assert response.status_code == 400
    assert response.json()["detail"]["errors"][0]["kind"] == "invalid_raster"


def test_nodata_must_be_zero(tmp_path: Path):
    values = np.ones((5, 5), dtype=np.uint16)
    path = tmp_path / "nodata.tif"
    _write_raster(path, values, nodata=65535)
    response = _post(path, _collection(_rect_geojson(1, 0, 0, 5, 5)))
    assert response.status_code == 400
    assert "NoData" in response.json()["detail"]["errors"][0]["message"]


def test_more_than_32_features_rejected(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    features = [_rect_geojson(i + 1, 0, 0, 1, 1) for i in range(33)]
    response = _post(parcel, _collection(*features))
    assert response.status_code == 400
    assert response.json()["detail"]["errors"][0]["kind"] == "too_many_parcels"


def test_pixel_centre_on_edge_counts_as_hit(tmp_path: Path):
    # Parcel 1 fills the whole 5x5 map.
    values = np.ones((5, 5), dtype=np.uint32)
    path = tmp_path / "edge.tif"
    _write_raster(path, values)

    # The rectangle's left edge runs exactly along the grid line between
    # columns 1 and 2 (an integral pixel boundary, robust to WGS84 round
    # trips).  Columns 2-4 are strictly inside; column 1 is outside.  This
    # pins down which side of the edge counts and the boundary points
    # themselves are covered via ``covers`` (see the centre-point probe
    # below).
    x_edge = ORIGIN_X + 2 * PIXEL
    x_right = ORIGIN_X + 5 * PIXEL
    y_top = ORIGIN_Y
    y_bottom = ORIGIN_Y - 5 * PIXEL
    lon, lat = warp_transform(
        CRS,
        "EPSG:4326",
        [x_edge, x_right, x_right, x_edge, x_edge],
        [y_top, y_top, y_bottom, y_bottom, y_top],
    )
    feature = {
        "type": "Feature",
        "properties": {"parcel_id": 1},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[x, y] for x, y in zip(lon, lat)]],
        },
    }
    response = _post(path, _collection(feature))
    assert response.status_code == 200, response.text
    out = tmp_path / "edge_out.tif"
    out.write_bytes(response.content)
    with rasterio.open(out) as ds:
        assert ds.dtypes[0] == "uint32"
        result = ds.read(1)
    np.testing.assert_array_equal(result[:, 0:2], 0)  # outside, old cells zeroed
    np.testing.assert_array_equal(result[:, 2:5], 1)  # inside columns kept

    # Direct probe: a point lying exactly on the polygon boundary must be
    # covered (``covers``), whereas ``contains`` would exclude it.
    from shapely.geometry import Point, box
    from shapely.prepared import prep

    polygon = prep(box(x_edge, y_bottom, x_right, y_top))
    edge_point = Point(x_edge, y_top - 2.5 * PIXEL)
    assert polygon.covers(edge_point)
    assert not polygon.contains(edge_point)


def test_windowed_writing_with_small_blocks(tmp_path: Path):
    # Exercise the tiled read/write path directly with a tiny block size.
    values = np.zeros((10, 10), dtype=np.uint16)
    values[:, :] = 7
    path = tmp_path / "small_blocks.tif"
    _write_raster(path, values)

    upload = UploadedRaster(path=str(path), filename="small_blocks.tif", role="parcel")
    meta = inspect_parcel_raster(upload)
    parcels = parse_geojson(_collection(_rect_geojson(7, 2, 2, 9, 9)), meta)
    plan = build_plan(meta, parcels, block_size=3)
    out = tmp_path / "small_blocks_out.tif"
    apply_correction(plan, str(out))

    with rasterio.open(out) as ds:
        result = ds.read(1)
        assert ds.block_shapes[0]  # tiled output still readable
    expected = np.zeros((10, 10), dtype=np.uint16)
    expected[2:9, 2:9] = 7
    np.testing.assert_array_equal(result, expected)


def test_duplicate_parcel_id_rejected(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    geojson = _collection(
        _rect_geojson(1, 0, 0, 2, 2),
        _rect_geojson(1, 3, 3, 5, 5),
    )
    response = _post(parcel, geojson)
    assert response.status_code == 400
    assert response.json()["detail"]["errors"][0]["kind"] == "duplicate_parcel_id"


def test_uint32_round_trip_keeps_dtype_and_large_id(tmp_path: Path):
    values = np.full((8, 8), 100_000, dtype=np.uint32)
    path = tmp_path / "u32.tif"
    _write_raster(path, values)
    response = _post(path, _collection(_rect_geojson(100_000, 0, 0, 8, 8)))
    assert response.status_code == 200, response.text
    out = tmp_path / "u32_out.tif"
    out.write_bytes(response.content)
    with rasterio.open(out) as ds:
        assert ds.dtypes[0] == "uint32"
        assert int(ds.read(1).flat[0]) == 100_000


def test_rotated_raster_rejected(tmp_path: Path):
    from affine import Affine

    path = tmp_path / "rotated.tif"
    transform = Affine(10, 0.1, ORIGIN_X, 0, -10, ORIGIN_Y)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=4,
        height=4,
        count=1,
        dtype=np.uint16,
        crs=CRS,
        transform=transform,
        nodata=0,
    ) as dst:
        dst.write(np.ones((4, 4), dtype=np.uint16), 1)

    response = _post(path, _collection(_rect_geojson(1, 0, 0, 4, 4)))
    assert response.status_code == 400
    assert response.json()["detail"]["errors"][0]["kind"] == "invalid_raster"


def test_non_polygon_geometry_rejected(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    payload = json.dumps(
        {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "properties": {"parcel_id": 1},
                    "geometry": {"type": "MultiPolygon", "coordinates": []},
                }
            ],
        }
    ).encode("utf-8")
    response = _post(parcel, payload)
    assert response.status_code == 400
    assert response.json()["detail"]["errors"][0]["kind"] == "invalid_geometry"


def test_non_wgs84_geojson_crs_rejected(tmp_path: Path):
    parcel = _base_raster(tmp_path)
    payload = json.dumps(
        {
            "type": "FeatureCollection",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::32650"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"parcel_id": 1},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [116.0, 30.0],
                                [116.01, 30.0],
                                [116.01, 30.01],
                                [116.0, 30.01],
                                [116.0, 30.0],
                            ]
                        ],
                    },
                }
            ],
        }
    ).encode("utf-8")
    response = _post(parcel, payload)
    assert response.status_code == 400
    assert "WGS84" in response.json()["detail"]["errors"][0]["message"]
