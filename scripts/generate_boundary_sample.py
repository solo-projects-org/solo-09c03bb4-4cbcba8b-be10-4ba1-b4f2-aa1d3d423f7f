#!/usr/bin/env python3
"""Generate the WGS84 boundary-update sample for /update-parcel-boundaries.

Run from the repository root:

    python scripts/generate_boundary_sample.py

The sample corrects parcels 101 and 102 of samples/parcel_map.tif.  Both new
boundaries are L-shaped (a single Polygon each) so that:

* parcel 101 extends east into the former western strip of parcel 102
  (rows 0-19, columns 55-59), and
* parcel 102 extends west into the former eastern strip of parcel 101
  (rows 20-39, columns 45-49).

The two new polygons share no pixel, so the batch passes the conflict checks
while exercising the "old pixels of batch parcels may be reused" rule.
Vertices are defined in the raster CRS (EPSG:32650) on pixel corners and
reprojected to WGS84 for the GeoJSON output.
"""

import json
from pathlib import Path

from rasterio.warp import transform_geom

OUT = Path(__file__).resolve().parent.parent / "samples"

PIXEL_SIZE = 10.0
ORIGIN_X = 500_000.0
ORIGIN_Y = 3_000_000.0

# UTM (x, y) corners of the two L-shaped polygons.  Each tuple is one
# polygon vertex in EPSG:32650; rings are closed below.
PARCEL_101_UTM = [
    (500_000.0, 3_000_000.0),  # top-left
    (500_600.0, 3_000_000.0),  # top-right (columns 0-59, rows 0-19)
    (500_600.0, 2_999_800.0),
    (500_450.0, 2_999_800.0),  # step west (columns 0-44, rows 20-39)
    (500_450.0, 2_999_600.0),
    (500_000.0, 2_999_600.0),  # bottom-left
]

PARCEL_102_UTM = [
    (500_600.0, 3_000_000.0),  # top-left (columns 60-109, rows 0-19)
    (501_100.0, 3_000_000.0),  # top-right
    (501_100.0, 2_999_600.0),  # bottom-right
    (500_450.0, 2_999_600.0),  # bottom-left (columns 45-109, rows 20-39)
    (500_450.0, 2_999_800.0),
    (500_600.0, 2_999_800.0),
]


def _to_wgs84(ring_utm: list[tuple[float, float]]) -> list[list[float]]:
    geom = {
        "type": "Polygon",
        "coordinates": [[[x, y] for x, y in ring_utm]],
    }
    projected = transform_geom("EPSG:32650", "EPSG:4326", geom)
    ring = [[float(lon), float(lat)] for lon, lat in projected["coordinates"][0]]
    ring.append(ring[0])  # GeoJSON rings are closed
    return ring


def main() -> None:
    features = []
    for parcel_id, ring in (
        (101, PARCEL_101_UTM),
        (102, PARCEL_102_UTM),
    ):
        features.append(
            {
                "type": "Feature",
                "id": parcel_id,
                "properties": {"parcel_id": parcel_id},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [_to_wgs84(ring)],
                },
            }
        )

    collection = {
        "type": "FeatureCollection",
        "name": "boundary_update_sample",
        "features": features,
    }
    out_path = OUT / "boundary_update_sample.geojson"
    out_path.write_text(
        json.dumps(collection, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Generated: {out_path.relative_to(OUT.parent)}")


if __name__ == "__main__":
    main()
