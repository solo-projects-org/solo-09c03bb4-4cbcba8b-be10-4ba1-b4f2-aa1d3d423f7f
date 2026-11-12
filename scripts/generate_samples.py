#!/usr/bin/env python3
"""Generate the demonstration GeoTIFFs and GeoJSON in samples/.

Run from the repository root:

    python scripts/generate_samples.py

The two mask scenes have different, non-zero tile origins.  Their grids are
offset by integral pixels from the parcel map and overlap each other, which is
exactly the condition the service must handle without resampling.

The boundary GeoJSON redefines parcels 101 and 102 so that each takes over
part of the other's old footprint; it exercises the /correct-boundaries
endpoint, including old-cell reuse and NoData zeroing.
"""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import transform as warp_transform

OUT = Path(__file__).resolve().parent.parent / "samples"

PIXEL_SIZE = 10.0
PARCEL_ORIGIN_X = 500_000.0
PARCEL_ORIGIN_Y = 3_000_000.0
WIDTH = 120
HEIGHT = 100

# (scene origin in parcels' pixel space, width, height)
SCENE_A = ((-7, -5), 100, 90)
SCENE_B = ((15, 10), 95, 85)


def write_raster(path: Path, data: np.ndarray, origin, nodata=0, crs="EPSG:32650"):
    height, width = data.shape
    transform = from_origin(
        origin[0], origin[1], xsize=PIXEL_SIZE, ysize=PIXEL_SIZE
    )
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=width,
        height=height,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
        tiled=True,
        blockxsize=64,
        blockysize=64,
        compress="deflate",
    ) as dst:
        dst.write(data, 1)


def global_grid(offset, width, height):
    """Coordinates in the parcel map's pixel index space."""
    rows = np.arange(offset[1], offset[1] + height)[:, None]
    cols = np.arange(offset[0], offset[0] + width)[None, :]
    return cols, rows


def write_boundary_sample():
    """Write the WGS84 boundary-correction GeoJSON.

    Parcels 101 and 102 swap a strip of each other's old ground:

    * new 101 = old 101 rows 0-29 + old 102 rows 0-29 (columns 55-109)
    * new 102 = old 102 rows 30-39 + old 101 rows 30-39 (columns 0-49)

    The two new polygons do not overlap, so the batch validates; old cells of
    101/102 that neither new boundary covers become NoData 0, while parcels
    103 and 104 stay untouched.
    """

    import json

    def rect_utm(col0, row0, col1, row1):
        x0 = PARCEL_ORIGIN_X + col0 * PIXEL_SIZE
        x1 = PARCEL_ORIGIN_X + col1 * PIXEL_SIZE
        y0 = PARCEL_ORIGIN_Y - row0 * PIXEL_SIZE
        y1 = PARCEL_ORIGIN_Y - row1 * PIXEL_SIZE
        return (x0, y0), (x1, y1)

    def to_wgs84_polygon(col0, row0, col1, row1):
        (x0, y_top), (x1, y_bottom) = rect_utm(col0, row0, col1, row1)
        xs = [x0, x1, x1, x0, x0]
        ys = [y_top, y_top, y_bottom, y_bottom, y_top]
        lon, lat = warp_transform(
            "EPSG:32650", "EPSG:4326", xs, ys
        )
        return [[[round(x, 10), round(y, 10)] for x, y in zip(lon, lat)]]

    features = []
    for parcel_id, (col0, row0, col1, row1) in [
        (101, (0, 0, 110, 30)),
        (102, (0, 30, 110, 40)),
    ]:
        features.append(
            {
                "type": "Feature",
                "id": parcel_id,
                "properties": {"parcel_id": parcel_id},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": to_wgs84_polygon(col0, row0, col1, row1),
                },
            }
        )

    collection = {
        "type": "FeatureCollection",
        "name": "parcel_boundary_corrections",
        "crs": {
            "type": "name",
            "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"},
        },
        "features": features,
    }
    path = OUT / "boundary_corrections_101_102.geojson"
    path.write_text(
        json.dumps(collection, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    parcels = np.zeros((HEIGHT, WIDTH), dtype=np.uint16)
    parcels[0:40, 0:50] = 101
    parcels[0:40, 55:110] = 102
    parcels[45:100, 0:50] = 103
    parcels[45:100, 55:110] = 104
    write_raster(
        OUT / "parcel_map.tif",
        parcels,
        (PARCEL_ORIGIN_X, PARCEL_ORIGIN_Y),
        nodata=0,
    )

    offset, width, height = SCENE_A
    cols, rows = global_grid(offset, width, height)
    mask_a = ((cols < 70) & (rows < 55)).astype(np.uint8)
    write_raster(
        OUT / "scene_01_cloud_free.tif",
        mask_a,
        (
            PARCEL_ORIGIN_X + offset[0] * PIXEL_SIZE,
            PARCEL_ORIGIN_Y - offset[1] * PIXEL_SIZE,
        ),
        nodata=0,
    )

    offset, width, height = SCENE_B
    cols, rows = global_grid(offset, width, height)
    mask_b = ((cols >= 50) & (rows >= 40)).astype(np.uint8)
    write_raster(
        OUT / "scene_02_cloud_free.tif",
        mask_b,
        (
            PARCEL_ORIGIN_X + offset[0] * PIXEL_SIZE,
            PARCEL_ORIGIN_Y - offset[1] * PIXEL_SIZE,
        ),
        nodata=0,
    )

    boundary_path = write_boundary_sample()

    print("Generated:")
    for path in sorted(OUT.glob("*.tif")):
        print(f"  {path.relative_to(OUT.parent)} ({path.stat().st_size} bytes)")
    print(
        f"  {boundary_path.relative_to(OUT.parent)} "
        f"({boundary_path.stat().st_size} bytes)"
    )


if __name__ == "__main__":
    main()
