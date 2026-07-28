"""
Raster export: ESRI ASCII Grid (.asc + .prj), and a more compact ENVI raw
binary (.bin + .hdr) for the full-resolution grid. Single-band or multi-band
(the twelve monthly climate layers of Tappa 2).

Why not GeoTIFF via rasterio: rasterio (and its GDAL dependency) isn't
installable in this sandbox (see noise.py's docstring for the network
constraint).

ESRI ASCII Grid is plain text, needs zero dependencies to write, and QGIS
opens it natively with the .prj sidecar assigning the CRS -- but at
23M cells and "%.2f" per value it comes out to roughly 8 bytes/cell
(~185 MB for this DEM's full 5334x4334 grid), which is unwieldy to hand
over as a chat attachment.

ENVI's raw-binary format (a flat array of numbers with NO header baked
into the same file, described instead by a separate plain-text .hdr) is
the standard compact alternative GDAL/QGIS already understand with no
extra plugins -- QGIS opens the .hdr directly. Two element types are
offered:
  - float32 ("f4"): full precision, 4 bytes/cell (~92 MB here).
  - int16   ("i2"): elevation rounded to the nearest metre, 2 bytes/cell
    (~46 MB here). For a fictional-world DEM whose whole point is visual
    terrain shape (not, say, surveying), 1 m vertical quantization is far
    below what's perceptible at 30 m horizontal resolution -- the same
    order as the noise/erosion process's own numerical noise. Rounding
    error introduced: uniform in [-0.5, 0.5] m, i.e. under 0.02% of this
    DEM's ~4000 m relief.
"""

import numpy as np


def write_esri_ascii_grid(path: str, array: np.ndarray, xmin: float, ymin: float, cellsize: float, nodata: float = -9999.0):
    """`array` must be oriented with row 0 = NORTH (top), matching how it
    was built from a np.meshgrid over an ascending y axis and then
    flipped -- see generate.py, which flips before calling this."""
    ny, nx = array.shape
    header = (
        f"ncols {nx}\n"
        f"nrows {ny}\n"
        f"xllcorner {xmin}\n"
        f"yllcorner {ymin}\n"
        f"cellsize {cellsize}\n"
        f"NODATA_value {nodata}\n"
    )
    out = np.where(np.isnan(array), nodata, array)
    with open(path, "w") as f:
        f.write(header)
        np.savetxt(f, out, fmt="%.2f")


def write_prj(path: str, proj4_string: str):
    with open(path, "w") as f:
        f.write(proj4_string.strip() + "\n")


_ENVI_DTYPE = {
    "f4": (4, "<f4"),   # (ENVI data type code, numpy dtype string) -- float32
    "i2": (2, "<i2"),   # int16
}


def write_envi_raw(
    path_stem: str,
    array: np.ndarray,
    xmin: float,
    ymin: float,
    cellsize: float,
    description: str,
    dtype: str = "i2",
    nodata: float = -9999.0,
    band_names: list[str] | None = None,
):
    """Write `path_stem`.bin (raw binary, row-major, BSQ) + `path_stem`.hdr
    (ENVI plain-text header). `array` must have row 0 = NORTH, same
    convention as write_esri_ascii_grid.

    Accepts either a single (ny, nx) band or a (bands, ny, nx) stack — the
    latter added in Tappa 2 so the twelve monthly climate layers ship as one
    multi-band raster QGIS opens as a single layer with a band selector,
    instead of twenty-four separate files. 2-D input behaves exactly as
    before.

    `dtype`: "f4" (float32, full precision) or "i2" (int16, metres rounded
    to the nearest integer -- see module docstring for why that's an
    acceptable tradeoff here). CRS is NOT embedded (ENVI's "coordinate
    system string" needs a full WKT this sandbox has no way to generate
    without GDAL) -- `map info` gives the affine georeferencing (origin +
    pixel size) so the raster lands in the right place/scale, but QGIS
    will import it with an unknown/no CRS. Assign the project's existing
    "Fictional World LCC" custom CRS to the layer manually after loading
    (Layer Properties -> Source -> Assigned CRS), the same one already
    defined in config/parameters.yml.
    """
    if dtype not in _ENVI_DTYPE:
        raise ValueError(f"dtype must be one of {list(_ENVI_DTYPE)}")
    code, np_dtype = _ENVI_DTYPE[dtype]
    if array.ndim == 2:
        array = array[None, ...]
    elif array.ndim != 3:
        raise ValueError("array must be (ny, nx) or (bands, ny, nx)")
    nbands, ny, nx = array.shape

    out = np.where(np.isnan(array), nodata, array)
    if dtype == "i2":
        out = np.clip(np.round(out), -32768, 32767)
    out.astype(np_dtype).tofile(path_stem + ".bin")   # BSQ = band-sequential,
    # which for a C-ordered (bands, ny, nx) array is exactly its memory layout

    # map info: pixel (1,1)'s (upper-left corner) maps to (xmin, y_top),
    # where y_top is the top edge of the grid AS ACTUALLY BUILT (ymin +
    # ny*cellsize) -- not necessarily the original xmax/ymax passed to
    # generate_dem, since np.arange(ymin, ymax, cellsize) doesn't always
    # land exactly on ymax. Deriving from the array's own shape keeps this
    # consistent with write_esri_ascii_grid's xllcorner/yllcorner, which
    # is anchored the same way (at xmin/ymin, corner convention).
    y_top = ymin + ny * cellsize
    hdr = (
        "ENVI\n"
        f"description = {{{description}}}\n"
        f"samples = {nx}\n"
        f"lines = {ny}\n"
        f"bands = {nbands}\n"
        "header offset = 0\n"
        "file type = ENVI Standard\n"
        f"data type = {code}\n"
        "interleave = bsq\n"
        "byte order = 0\n"
        f"map info = {{Arbitrary, 1, 1, {xmin}, {y_top}, {cellsize}, {cellsize}, units=Meters}}\n"
        f"data ignore value = {nodata}\n"
        + (f"band names = {{{', '.join(band_names)}}}\n" if band_names else "")
    )
    with open(path_stem + ".hdr", "w") as f:
        f.write(hdr)
