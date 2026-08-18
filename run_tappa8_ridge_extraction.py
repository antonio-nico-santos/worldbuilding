import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from hydrology.flow import priority_flood_d8, accumulate_flow

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- real ridge-crest extraction (inverted-DEM flow accumulation) ===")
log("Same technique as Tappa 4's stream extraction, run on -DEM: 'streams' of the")
log("inverted surface are the ridges of the real one. Validated first on a 400x400")
log("crop around the Spine's midpoint (visual match against hillshade, see chat) --")
log("this is the full-domain commit run.")

OUTPUT_DIR = "data/processed/geomorphology"
os.makedirs(OUTPUT_DIR, exist_ok=True)

dem = np.load("data/processed/dem_v3_final_30m_eroded.npy").astype(np.float64)
ny, nx = dem.shape
land = dem > 0.0
log(f"DEM shape={dem.shape} land_fraction={land.mean():.4f}")

inv = dem.max() - dem

# Seed = domain boundary only, NOT "ocean" (dem<=0). Ocean cells are the
# lowest points of the ORIGINAL surface, i.e. the HIGHEST points of the
# inverted one -- not a valid drainage exit for the inverted problem the way
# they are for the real one. The domain boundary alone is a generic, always-
# valid seed set for the fill algorithm regardless of which surface it's
# fed.
seed = np.zeros((ny, nx), dtype=bool)
seed[0, :] = seed[-1, :] = seed[:, 0] = seed[:, -1] = True

log("running priority-flood fill + D8 direction on inverted DEM (long step)...")
filled_inv, receiver, pop_order = priority_flood_d8(inv, seed, epsilon=1e-4)
log("fill + direction done")

log("accumulating (plain cell-count, land only -- topology, not a real quantity)...")
w = land.astype(np.float64)
ridge_accum_cells = accumulate_flow(receiver, pop_order, w).reshape(ny, nx)
log(f"done. max={ridge_accum_cells.max():.0f} p99={np.percentile(ridge_accum_cells, 99):.1f}")

np.save(f"{OUTPUT_DIR}/ridge_accum_cells.npy", ridge_accum_cells.astype(np.float32))

meta = {
    "method": "inverted-DEM (dem.max()-dem) priority-flood + D8 + plain accumulation, "
              "boundary-only seed (not ocean -- see script docstring)",
    "runtime_s": time.time() - t_start,
    "max_accum_cells": float(ridge_accum_cells.max()),
    "p99_accum_cells": float(np.percentile(ridge_accum_cells, 99)),
}
with open(f"{OUTPUT_DIR}/tappa8_ridge_extraction_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
