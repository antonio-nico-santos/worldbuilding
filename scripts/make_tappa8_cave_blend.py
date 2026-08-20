"""
Tappa 8 -- blend all five cave eligibility masks into ONE raster.

Nico asked whether the five separate cave rasters (lava_tube, talus_pseudokarst,
glacier_moulin, sea_cave, karst) could be one file with the pixel value carrying the
type info, instead of five files.

Checked before picking an encoding: the five masks are NOT mutually exclusive. Real
overlap exists -- talus_pseudokarst alone shares 111.0 km2 with glacier_moulin, 142.6
km2 with sea_cave, and 55.2 km2 with karst (all three make geological sense: steep
terrain near a stream can simultaneously sit in the alpine snow zone, the coastal
buffer, or a soluble-rock zone). 342,673 cells (308.4 km2) carry 2+ types
simultaneously, and a handful carry 3. A single mutually-exclusive category code would
have to silently pick a winner and discard the others -- real information loss.

So this is a BITMASK raster, not a category-code raster: each bit is one cave type,
independently 0/1, and a cell's value is the OR of every type present. Value range is
0-31 (5 bits), fits uint8 (ENVI "u1") with no precision loss. Any GIS user recovers a
single type with a bitwise AND against its bit value (documented in the .hdr
description and the meta JSON below), or filters "any cave" with `value > 0`.

Run after run_tappa8_caves_v4.py (needs cave_karst.npy) and with the four existing
cave_*.npy already on disk.
"""
import sys, time, json
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUTPUT_DIR = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

# bit -> (name, source file, value)
BITS = [
    ("lava_tube", "cave_lava_tube_v3.npy", 1),
    ("talus_pseudokarst", "cave_talus_pseudokarst.npy", 2),
    ("glacier_moulin", "cave_glacier_moulin_v2.npy", 4),
    ("sea_cave", "cave_sea_cave.npy", 8),
    ("karst", "cave_karst.npy", 16),
]

log("loading all five cave eligibility masks...")
masks = {}
for name, fname, bit in BITS:
    masks[name] = np.load(f"{OUTPUT_DIR}/{fname}").astype(bool)

ny, nx = next(iter(masks.values())).shape
blend = np.zeros((ny, nx), dtype=np.uint8)
for name, fname, bit in BITS:
    blend |= (masks[name].astype(np.uint8) * bit)

log("computing overlap diagnostics (documented, not just asserted)...")
names = [b[0] for b in BITS]
pairwise_km2 = {}
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        n_cells = int((masks[names[i]] & masks[names[j]]).sum())
        if n_cells:
            pairwise_km2[f"{names[i]}+{names[j]}"] = round(n_cells * cell_km2, 3)

n_types_per_cell = np.stack(list(masks.values()), axis=0).sum(axis=0)
multi_type_km2 = float((n_types_per_cell >= 2).sum() * cell_km2)
max_simultaneous = int(n_types_per_cell.max())
area_km2 = {name: float(masks[name].sum() * cell_km2) for name, _, _ in BITS}
any_cave_km2 = float((blend > 0).sum() * cell_km2)

log(f"per-type area km2: {area_km2}")
log(f"pairwise overlap km2 (only nonzero pairs shown): {pairwise_km2}")
log(f"cells with 2+ simultaneous types: {multi_type_km2:.2f} km2, max simultaneous on one cell: {max_simultaneous}")
log(f"any-cave-type union: {any_cave_km2:.2f} km2")

log("exporting blended bitmask raster...")
np.save(f"{OUTPUT_DIR}/cave_blend.npy", blend)
description = (
    "Tappa 8 cave eligibility BITMASK (not mutually exclusive -- OR of all five types): "
    + ", ".join(f"{bit}={name}" for name, _, bit in BITS)
    + ". Decode a single type with (value & bit); value==0 means no cave type eligible; "
      "value>0 means at least one type eligible."
)
write_envi_raw(f"{OUTPUT_DIR}/cave_blend", blend.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=description, dtype="u1")
write_prj(f"{OUTPUT_DIR}/cave_blend.prj", CRS_PROJ4)

meta = {
    "encoding": "bitmask, not category code -- the five types genuinely overlap (see overlap_km2 below), "
                "so a single mutually-exclusive code would silently discard information",
    "bits": {name: bit for name, _, bit in BITS},
    "area_km2": area_km2,
    "overlap_km2": pairwise_km2,
    "cells_with_2plus_types_km2": multi_type_km2,
    "max_simultaneous_types_on_one_cell": max_simultaneous,
    "any_cave_type_union_km2": any_cave_km2,
    "how_to_decode": "single_type_mask = (raster_value & bit_value) > 0; any_cave_mask = raster_value > 0",
}
with open(f"{OUTPUT_DIR}/tappa8_cave_blend_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
