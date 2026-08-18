import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_VOLCANIC
from geomorphology.caves import lava_tube_candidates_v2

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- cave candidates v3 (lava tube ONLY, re-run against lithology v5 "
    "[DEM-native, island now split volcanic/basin_fill]; talus/sea/glacier UNCHANGED, "
    "none of them consume lithology at all) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]

OUTPUT_DIR = "data/processed/geomorphology"

log("loading DEM + lithology v5 (DEM-native, island volcanic/basin_fill split) + reused v1/v2 fields...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
lithology_v5 = np.load(f"{OUTPUT_DIR}/lithology_v5.npy")
slope_pct = np.load(f"{OUTPUT_DIR}/slope_pct_30m.npy")
vent_weight = np.load(f"{OUTPUT_DIR}/lava_tube_vent_weight.npy")  # vent-location field, independent of lithology -- unchanged

log("--- lava tube v3: SAME test as v2 (vent proximity AND gentle slope, both lithology-independent), "
    "just re-masked against lithology_v5's CLASS_VOLCANIC (island shrank from 100% to 75% volcanic; "
    "the other 25%, the flat isthmus/apron margin, is now basin_fill and confirmed by Nico as the "
    "correct reading) ---")
volcanic_mask_v5 = (lithology_v5 == CLASS_VOLCANIC) & land_mask
volcanic_slope_p25 = float(np.percentile(slope_pct[volcanic_mask_v5], 25.0))
VENT_WEIGHT_THRESHOLD = 0.1  # unchanged placeholder from v2, not independently calibrated
log(f"  volcanic-zone-local gentle threshold (p25, recomputed within the SMALLER v5 volcanic zone) = "
    f"{volcanic_slope_p25:.2f}%, vent_weight threshold = {VENT_WEIGHT_THRESHOLD} (unchanged placeholder)")
lava_tube_v3 = lava_tube_candidates_v2(
    lithology_v5, land_mask, vent_weight, VENT_WEIGHT_THRESHOLD, slope_pct, volcanic_slope_p25,
)

cell_km2 = (RES_M / 1000.0) ** 2
lava_tube_v2 = np.load(f"{OUTPUT_DIR}/cave_lava_tube_v2.npy")
areas = {
    "lava_tube_v2_km2_lithology_v2_full_island_volcanic": float(lava_tube_v2.sum() * cell_km2),
    "lava_tube_v3_km2_lithology_v5_island_75pct_volcanic": float(lava_tube_v3.sum() * cell_km2),
    "volcanic_zone_km2_v2": float(((np.load(f'{OUTPUT_DIR}/lithology_v2.npy') == CLASS_VOLCANIC) & land_mask).sum() * cell_km2),
    "volcanic_zone_km2_v5": float(volcanic_mask_v5.sum() * cell_km2),
}
log(f"areas: {areas}")

log("checking whether any v2 lava-tube candidate cells fell OUTSIDE the new (smaller) v5 volcanic zone "
    "-- i.e. cells that used to be eligible and no longer are...")
lost_to_reclassification = lava_tube_v2 & ~volcanic_mask_v5
gained = lava_tube_v3 & ~lava_tube_v2
log(f"  lost (was candidate under v2, now outside volcanic zone under v5): "
    f"{int(lost_to_reclassification.sum())} cells ({lost_to_reclassification.sum()*cell_km2:.3f} km2)")
log(f"  gained (newly a candidate under v3, wasn't under v2): "
    f"{int(gained.sum())} cells ({gained.sum()*cell_km2:.3f} km2)")

log("exporting raster...")
np.save(f"{OUTPUT_DIR}/cave_lava_tube_v3.npy", lava_tube_v3)
write_envi_raw(
    f"{OUTPUT_DIR}/cave_lava_tube_v3", lava_tube_v3.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description="Tappa 8 cave candidate eligibility mask v3: lava tube, re-masked against lithology_v5",
    dtype="i2",
)
write_prj(f"{OUTPUT_DIR}/cave_lava_tube_v3.prj", CRS_PROJ4)

meta = {
    "change_reason": "lithology v5 (DEM-native elevation+relief classification, extended to split the "
                      "SW Island into volcanic core (75%) / basin_fill flat isthmus-apron margin (25%), "
                      "confirmed correct by Nico) replaces lithology v2 as the lava-tube eligibility input. "
                      "talus/pseudokarst, sea caves, and glacier/moulin are UNCHANGED and NOT re-run here: "
                      "none of their candidate functions take lithology as an input at all.",
    "volcanic_zone_local_slope_p25_pct": volcanic_slope_p25,
    "vent_weight_threshold": VENT_WEIGHT_THRESHOLD,
    "areas_km2": areas,
    "lost_to_reclassification_km2": float(lost_to_reclassification.sum() * cell_km2),
    "gained_km2": float(gained.sum() * cell_km2),
}
with open(f"{OUTPUT_DIR}/tappa8_caves_v3_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
