"""
Tappa 7 -- convert the fauna .npy raster outputs to GeoTIFF so QGIS (or any
standard GIS) can actually open them. QGIS has no .npy reader: a raw NumPy
array carries no georeferencing (no CRS, no transform, no band metadata) --
it's just a binary blob of numbers to any tool outside this project's own
Python scripts. GeoTIFF is the fix, not a workaround: same numbers, same
grid, now self-describing.

Grid geometry -- same as every other Tappa 5/6/7 raster in this project,
confirmed against data/processed/biomes/biome_id.hdr and .prj:
  shape (1334, 1084), 120m cells, origin (xmin=-65000.0, ymax=80000.0),
  CRS: project LCC proj4 (see CRS_PROJ4 below).

Converts every 2D .npy file in data/processed/fauna/ whose shape matches
that grid. Boolean arrays (outpost_search_mask, povo_livre_zone) are written
as uint8 (0/1) -- GeoTIFF has no native bool dtype. int16 label rasters
(nacre_massif_labels: 0=none/1=main spine/2=South Branch) keep their
integer dtype. float32 suitability/exposure rasters keep float32, with 0.0
land-background values left as-is (they are real "not suitable here" zeros
from the composite math, not nodata -- there is no separate ocean mask
applied to these composites beyond each species' own exclusao, so no nodata
value is set here; ocean cells read as the same low/zero values as
unsuitable land, a known limitation already flagged for the source rasters
in 07_tappa7_regional_scenario.md, not introduced by this conversion).

Output -> data/processed/fauna/*.tif (one per source .npy, same basename).
"""
import glob
import os

import numpy as np
import rasterio

BASE = 'data/processed'
FAUNA_DIR = 'data/processed/fauna'
CRS_PROJ4 = ("+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
             "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
EXPECTED_SHAPE = (1334, 1084)
RES_X, RES_Y = 120.0, 120.0
XMIN, YMAX = -65000.0, 80000.0


def main():
    transform = rasterio.transform.from_origin(XMIN, YMAX, RES_X, RES_Y)
    npy_files = sorted(glob.glob(f'{FAUNA_DIR}/*.npy'))
    if not npy_files:
        print(f'No .npy files found under {FAUNA_DIR}')
        return

    written = []
    skipped = []
    for f in npy_files:
        arr = np.load(f)
        if arr.shape != EXPECTED_SHAPE:
            skipped.append((f, arr.shape))
            continue

        if arr.dtype == np.bool_:
            out_arr = arr.astype(np.uint8)
            dtype = 'uint8'
        elif arr.dtype in (np.float32, np.float64):
            out_arr = arr.astype(np.float32)
            dtype = 'float32'
        else:
            out_arr = arr
            dtype = str(arr.dtype)

        out_path = os.path.splitext(f)[0] + '.tif'
        with rasterio.open(
            out_path, 'w', driver='GTiff',
            height=out_arr.shape[0], width=out_arr.shape[1],
            count=1, dtype=dtype,
            crs=CRS_PROJ4, transform=transform,
            compress='deflate',
        ) as dst:
            dst.write(out_arr, 1)
        written.append(out_path)
        print(f'{f} ({arr.dtype}, {arr.shape}) -> {out_path}')

    if skipped:
        print('\nSkipped (shape does not match the standard 1334x1084 grid -- inspect manually):')
        for f, shp in skipped:
            print(f'  {f}: {shp}')

    print(f'\n{len(written)} GeoTIFF(s) written to {FAUNA_DIR}/.')
    print('In QGIS: Layer -> Add Layer -> Add Raster Layer, point at the .tif files directly -- '
          'CRS should auto-detect from the embedded proj4 string this time (unlike the GeoJSON '
          'legacy-crs-member issue, GeoTIFF CRS embedding is standard and GDAL reads it natively).')


if __name__ == '__main__':
    main()
