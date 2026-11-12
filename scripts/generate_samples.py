#!/usr/bin/env python3
"""Generate the demonstration GeoTIFFs in samples/.

Run from the repository root:

    python scripts/generate_samples.py

The two mask scenes have different, non-zero tile origins.  Their grids are
offset by integral pixels from the parcel map and overlap each other, which is
exactly the condition the service must handle without resampling.
"""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

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

    print("Generated:")
    for path in sorted(OUT.glob("*.tif")):
        print(f"  {path.relative_to(OUT.parent)} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
