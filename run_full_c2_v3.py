import sys, time, json
sys.path.insert(0, 'src')
import numpy as np
from params import load_params
from terrain.generate import generate_dem
from terrain.erosion import erode
from terrain.raster_io import write_envi_raw, write_prj

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 1 full regeneration: v3 (C2 real-detail + warp-fold fix) ===")

# All tuning parameters below come from config/parameters.yml -- see
# docs/decisions/01_tappa1_terrain.md S11c for the rationale behind each
# value. Previously these were hardcoded here, which had drifted out of
# sync with parameters.yml at least once already; loading them avoids
# relying on this script (or memory) staying in lockstep with the config
# file, which the file's own header declares to be the single source of
# truth.
params = load_params("config/parameters.yml")
domain = params["domain"]
terrain = params["terrain"]
erosion_cfg = terrain["erosion"]
real_detail_cfg = terrain["real_detail"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]

DOMAIN = dict(
    xmin=domain["xmin"], xmax=domain["xmax"],
    ymin=domain["ymin"], ymax=domain["ymax"],
    resolution_m=domain["resolution_m"],
)

# Ridge/zone skeleton paths and the two per-feature override dicts are NOT
# in parameters.yml: they're keyed by the feature names authored in
# data/input/terrain_ridges.geojson and terrain_zones.geojson, so they
# belong with the skeleton, not the generic tuning-parameter config. See
# decision doc S3b/S8 for what each override does.
RIDGES_PATH = "data/input/terrain_ridges.geojson"
ZONES_PATH = "data/input/terrain_zones.geojson"
SHELF_MULTIPLIERS = {
    "Spine": 1.6,
    "North branch (Big Brother)": 1.3,
    "West Branch (Little Brother)": 1.3,
    "South Branch": 1.3,
}
ZONE_BASE_LIFT_M = {"North plains": 300.0}

log("generating DEM (parameters loaded from config/parameters.yml)...")
dem = generate_dem(
    **DOMAIN,
    ridges_path=RIDGES_PATH, zones_path=ZONES_PATH,
    seed=terrain["seed"],
    noise_octaves=terrain["noise_octaves"],
    noise_base_wavelength_m=terrain["noise_base_wavelength_m"],
    noise_amplitude_m=terrain["noise_amplitude_m"],
    sea_level_offset_m=terrain["sea_level_offset_m"],
    lacunarity=terrain["lacunarity"],
    persistence=terrain["persistence"],
    noise_plain_octaves=terrain["noise_plain_octaves"],
    noise_persistence_fine=terrain["noise_persistence_fine"],
    noise_crossover_wavelength_m=terrain["noise_crossover_wavelength_m"],
    noise_crossover_width_factor=terrain["noise_crossover_width_factor"],
    shelf_multipliers=SHELF_MULTIPLIERS,
    zone_base_lift_m=ZONE_BASE_LIFT_M,
    warp_wavelength_m=terrain["warp_wavelength_m"],
    warp_amplitude_m=terrain["warp_amplitude_m"],
    noise_warp_wavelength_m=terrain["noise_warp_wavelength_m"],
    noise_warp_amplitude_m=terrain["noise_warp_amplitude_m"],
    warp_octaves=terrain["warp_octaves"],
    warp_lacunarity=terrain["warp_lacunarity"],
    warp_persistence=terrain["warp_persistence"],
    noise_warp_octaves=terrain["noise_warp_octaves"],
    noise_warp_lacunarity=terrain["noise_warp_lacunarity"],
    noise_warp_persistence=terrain["noise_warp_persistence"],
    real_detail_path=real_detail_cfg["path"],
    real_detail_xmin=real_detail_cfg["xmin"],
    real_detail_ymax=real_detail_cfg["ymax"],
    real_detail_cellsize_m=real_detail_cfg["cellsize_m"],
    real_detail_fine_supplement_weight=real_detail_cfg["fine_supplement_weight"],
    real_detail_fine_min_wavelength_m=real_detail_cfg["fine_min_wavelength_m"],
    row_chunk_size=250, verbose=True,
)
log(f"DEM generated: shape={dem.shape} range={dem.min():.0f}..{dem.max():.0f}m land={100*(dem>0).mean():.1f}%")
assert not np.isnan(dem).any() and not np.isinf(dem).any(), "NaN/Inf in raw DEM!"
np.save("dem_v3_final_30m_raw.npy", dem)

land_cells = int((dem > 0).sum())
n_droplets = int(land_cells * erosion_cfg["n_droplets_per_land_cell"])
log(f"land cells={land_cells:,} -> n_droplets={n_droplets:,} ({erosion_cfg['n_droplets_per_land_cell']}/land-cell)")

log("running erosion (this is the long step)...")
eroded = erode(
    dem, cell_size_m=domain["resolution_m"], n_droplets=n_droplets,
    seed=erosion_cfg["seed"], max_steps=erosion_cfg["max_steps"],
    erode_speed=erosion_cfg["erode_speed"], deposit_speed=erosion_cfg["deposit_speed"],
    verbose=True,
)
log(f"erosion done: range={eroded.min():.0f}..{eroded.max():.0f}m")
assert not np.isnan(eroded).any() and not np.isinf(eroded).any(), "NaN/Inf after erosion!"
np.save("dem_v3_final_30m_eroded.npy", eroded)

diff = eroded - dem
log(f"erosion diff: min={diff.min():.2f} max={diff.max():.2f} mean={diff.mean():.4f}")

land_frac = (eroded > 0).mean()
log(f"FINAL sanity: shape={eroded.shape} range={eroded.min():.1f}..{eroded.max():.1f}m land={100*land_frac:.1f}% nan=False inf=False")

log("exporting ENVI raw binary (int16) + hdr + prj...")
write_envi_raw(
    "dem_v3_final_30m_eroded", eroded, xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 1 final DEM v3, 30m, eroded, C2 real-detail-textured + warp-fold fix (Fictional World LCC domain)",
    dtype="i2",
)
write_prj("dem_v3_final_30m_eroded.prj", CRS_PROJ4)
log("export done")

meta = {
    "shape": list(eroded.shape), "elev_min": float(eroded.min()), "elev_max": float(eroded.max()),
    "land_fraction": float(land_frac), "n_droplets": n_droplets, "land_cells": land_cells,
}
with open("dem_v3_run_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== ALL DONE ===")
