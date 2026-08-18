import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_VOLCANIC, identify_landmasses
from geomorphology.caves import (
    resample_nearest, lava_tube_candidates_v2, vent_proximity_weight,
    talus_pseudokarst_candidates, glacier_moulin_candidates_v2,
    sea_cave_candidates, distance_to_ocean_km, compute_slope,
)

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- cave candidates v2 (narrow glacier/moulin + lava tube; "
    "talus + sea caves UNCHANGED, already approved) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]

OUTPUT_DIR = "data/processed/geomorphology"
os.makedirs(OUTPUT_DIR, exist_ok=True)

log("loading DEM + lithology v2 + stream_mask + v1 cached fields...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
lithology = np.load(f"{OUTPUT_DIR}/lithology_v2.npy")   # v2 lithology (real crest), not v1
stream_mask = np.load("data/processed/hydrology/stream_mask.npy")
ny, nx = dem.shape

# reuse fields already computed in the v1 caves run -- unchanged inputs
slope_pct = np.load(f"{OUTPUT_DIR}/slope_pct_30m.npy")
dist_to_stream_km = np.load(f"{OUTPUT_DIR}/dist_to_stream_km_30m.npy")
dist_to_ocean_km = np.load(f"{OUTPUT_DIR}/dist_to_ocean_km_30m.npy")
vent_weight = np.load(f"{OUTPUT_DIR}/lava_tube_vent_weight.npy")

log("identifying mainland (excludes SW Island + every other islet)...")
labeled, mainland_label, sw_island_label, areas_km2 = identify_landmasses(land_mask, cellsize_m=RES_M)
mainland_mask = labeled == mainland_label
steep_threshold_pct = float(np.percentile(slope_pct[mainland_mask], 75.0))
log(f"steep threshold (p75 of mainland land slope, SAME as v1 talus/sea) = {steep_threshold_pct:.2f}%")

log("resampling permanent_snow_mask (120m climate grid -> native 30m)...")
with open("data/processed/climate/permanent_snow_mask.hdr") as f:
    hdr_txt = f.read()
def hdr_val(key):
    for line in hdr_txt.splitlines():
        if line.strip().startswith(key):
            return line.split("=", 1)[1].strip()
    raise KeyError(key)
psm_samples = int(hdr_val("samples"))
psm_lines = int(hdr_val("lines"))
psm_mapinfo = hdr_val("map info").strip("{}").split(",")
psm_xmin = float(psm_mapinfo[3])
psm_ytop = float(psm_mapinfo[4])
psm_cellsize = float(psm_mapinfo[5])
psm = np.fromfile("data/processed/climate/permanent_snow_mask.bin", dtype="<i2").reshape(psm_lines, psm_samples)
psm_native = resample_nearest(
    psm, coarse_xmin=psm_xmin, coarse_ymax=psm_ytop, coarse_cellsize=psm_cellsize,
    ny=ny, nx=nx, xmin=domain["xmin"], ymax=domain["ymax"], cellsize=RES_M,
)
snow_mask = psm_native.astype(bool) & land_mask
log(f"snow mask (unchanged from v1): {snow_mask.sum() * (RES_M/1000)**2:.2f} km2")

log("--- glacier/moulin v2: slope (p75 within snow mask) AND margin proximity ---")
snow_slope_p75 = float(np.percentile(slope_pct[snow_mask], 75.0))
margin_depth_km = ndimage.distance_transform_edt(snow_mask, sampling=(RES_M, RES_M)) / 1000.0
margin_depth_p50 = float(np.percentile(margin_depth_km[snow_mask], 50.0))
log(f"  snow-mask-local steep threshold (p75) = {snow_slope_p75:.2f}%, "
    f"margin-depth threshold (median) = {margin_depth_p50:.2f} km")
glacier_moulin_v2 = glacier_moulin_candidates_v2(
    psm_native, land_mask, slope_pct, snow_slope_p75, margin_depth_km, margin_depth_p50,
)

log("--- lava tube v2: vent proximity (weight>=0.1) AND gentle slope (p25 within volcanic zone) ---")
volcanic_mask_v2 = (lithology == CLASS_VOLCANIC) & land_mask
volcanic_slope_p25 = float(np.percentile(slope_pct[volcanic_mask_v2], 25.0))
VENT_WEIGHT_THRESHOLD = 0.1  # placeholder -- not independently calibrated, see write-up
log(f"  volcanic-zone-local gentle threshold (p25) = {volcanic_slope_p25:.2f}%, "
    f"vent_weight threshold = {VENT_WEIGHT_THRESHOLD} (placeholder)")
lava_tube_v2 = lava_tube_candidates_v2(
    lithology, land_mask, vent_weight, VENT_WEIGHT_THRESHOLD, slope_pct, volcanic_slope_p25,
)

log("--- talus/pseudokarst + sea caves: UNCHANGED from v1 (approved as-is) ---")
talus = talus_pseudokarst_candidates(slope_pct, dist_to_stream_km, mainland_mask, steep_threshold_pct, stream_buffer_km=0.5)
sea_caves = sea_cave_candidates(land_mask, dist_to_ocean_km, slope_pct, steep_threshold_pct, coastal_buffer_km=0.5)

cell_km2 = (RES_M / 1000.0) ** 2
areas_v1 = {
    "lava_tube_v1_km2": 794.2850999999999,
    "glacier_moulin_v1_km2": 960.1632,
}
areas_v2 = {
    "lava_tube_v2_km2": float(lava_tube_v2.sum() * cell_km2),
    "talus_pseudokarst_km2": float(talus.sum() * cell_km2),
    "glacier_moulin_v2_km2": float(glacier_moulin_v2.sum() * cell_km2),
    "sea_cave_km2": float(sea_caves.sum() * cell_km2),
}
log(f"v1 areas (for comparison): {areas_v1}")
log(f"v2 areas: {areas_v2}")

log("exporting rasters...")
for name, arr in [
    ("cave_lava_tube_v2", lava_tube_v2),
    ("cave_glacier_moulin_v2", glacier_moulin_v2),
]:
    np.save(f"{OUTPUT_DIR}/{name}.npy", arr)
    write_envi_raw(
        f"{OUTPUT_DIR}/{name}", arr.astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 cave candidate eligibility mask v2: {name}",
        dtype="i2",
    )
    write_prj(f"{OUTPUT_DIR}/{name}.prj", CRS_PROJ4)

meta = {
    "resolution_m": RES_M,
    "steep_threshold_pct_p75_mainland_unchanged": steep_threshold_pct,
    "glacier_moulin_v2": {
        "snow_mask_local_slope_p75_pct": snow_slope_p75,
        "margin_depth_median_km": margin_depth_p50,
        "method": "snow_mask & land & slope>=p75(within snow mask) & margin_depth_km<=median(within snow mask)",
    },
    "lava_tube_v2": {
        "volcanic_zone_local_slope_p25_pct": volcanic_slope_p25,
        "vent_weight_threshold": VENT_WEIGHT_THRESHOLD,
        "method": "volcanic(v2 lithology) & land & vent_weight>=0.1 & slope<=p25(within volcanic zone) "
                  "-- p25 mirrors talus/sea's p75-steep test, inverted",
    },
    "areas_v1_km2": areas_v1,
    "areas_v2_km2": areas_v2,
    "talus_sea_unchanged": "talus_pseudokarst and sea_cave use IDENTICAL logic/thresholds to v1, "
                            "just recomputed here for a consistent single meta file -- approved as-is",
}
with open(f"{OUTPUT_DIR}/tappa8_caves_v2_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
