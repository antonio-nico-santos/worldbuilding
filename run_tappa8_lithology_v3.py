import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    load_ridges, load_zones, classify, jade_eligible_mask, place_jade_pods,
    LITHOLOGY_CLASSES, CLASS_SCHIST,
)
from geomorphology.real_crest import extract_real_crest_cross_section, height_normalized_falloff_m
from geomorphology.lithology import _grid_xy
from terrain.skeleton import load_geojson

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology v3 (Option A: cross-section snapping + height-normalized falloff) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]

OUTPUT_DIR = "data/processed/geomorphology"

log("loading DEM, skeleton...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape
land_mask = dem > 0
ridges = load_ridges("data/input/terrain_ridges.geojson")
zones = load_zones("data/input/terrain_zones.geojson")

log("building full grid coordinates...")
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

# PLACEHOLDER, first pass -- not independently calibrated:
CROSS_SECTION_HALFWIDTH_M = 1000.0   # tightest of the 3 values tested (1/2/3 km); kept
                                       # jump between consecutive stations to <=2.1km on
                                       # every ridge, with zero empty bins.
FALLOFF_LOW_MULT = 0.6   # Nico's explicit values: lowest crest point -> 0.6x base falloff
FALLOFF_HIGH_MULT = 1.0  #                          highest crest point -> 1.0x base falloff

log(f"extracting real crest per ridge (Option A: cross-section snapping, "
    f"halfwidth={CROSS_SECTION_HALFWIDTH_M}m, DEM elevation argmax per station)...")
raw_ridge_features = load_geojson("data/input/terrain_ridges.geojson")
ridge_coords_by_name = {f["properties"]["name"]: np.array(f["geometry"]["coordinates"], dtype=np.float64) for f in raw_ridge_features}

real_crest_trees = {}
real_crest_falloff_m = {}
crest_info = {}
for ridge in ridges:
    tree, elev, info = extract_real_crest_cross_section(
        ridge, ridge_coords_by_name[ridge.name], dem, land_mask,
        xmin=domain["xmin"], ymax=domain["ymax"], cellsize_m=domain["resolution_m"],
        bin_spacing_m=1000.0, cross_section_halfwidth_m=CROSS_SECTION_HALFWIDTH_M,
    )
    crest_info[ridge.name] = info
    log(f"  {ridge.name}: {info}")
    if tree is not None:
        real_crest_trees[ridge.name] = tree
        falloff_arr = height_normalized_falloff_m(elev, ridge.falloff_km, FALLOFF_LOW_MULT, FALLOFF_HIGH_MULT)
        real_crest_falloff_m[ridge.name] = falloff_arr
        log(f"    elev range [{elev.min():.0f},{elev.max():.0f}]m -> "
            f"falloff range [{falloff_arr.min():.0f},{falloff_arr.max():.0f}]m "
            f"(base {ridge.falloff_km*1000:.0f}m)")

log("classifying lithology against REAL crest (height-normalized falloff)...")
result = classify(
    dem, ridges, zones,
    xmin=domain["xmin"], xmax=domain["xmax"], ymin=domain["ymin"], ymax=domain["ymax"],
    row_chunk_size=500, real_crest_trees=real_crest_trees, real_crest_falloff_m=real_crest_falloff_m,
)
lithology = result["lithology"]
basin_fill_grounded = result["basin_fill_grounded"]
schist_grade = result["schist_grade"]

cell_km2 = (domain["resolution_m"] / 1000.0) ** 2
class_areas_km2 = {name: float((lithology == code).sum() * cell_km2) for code, name in LITHOLOGY_CLASSES.items() if code != 0}
log(f"class areas km2: {class_areas_km2}")

basin_fill_total_km2 = class_areas_km2["sedimentary_basin_fill"]
basin_fill_grounded_km2 = float(basin_fill_grounded.sum() * cell_km2)
log(f"basin fill grounded fraction: {100*basin_fill_grounded_km2/basin_fill_total_km2:.1f}%")

log("placing jade/pounamu pods (stochastic, weighted by grade, seed=13)...")
pod_mask, suitable_mask, pod_centers_xy, pod_radii_m = place_jade_pods(
    lithology, schist_grade, xx, yy, domain["resolution_m"],
    n_pods=10, min_separation_km=5.0, radius_range_m=(300.0, 800.0),
    grade_percentile=80.0, seed=13,
)
jade_pod_area_km2 = float(pod_mask.sum() * cell_km2)
suitable_area_km2 = float(suitable_mask.sum() * cell_km2)
log(f"jade: {len(pod_centers_xy)} pods placed, pod area={jade_pod_area_km2:.2f} km2, "
    f"suitability zone={suitable_area_km2:.1f} km2")

log("exporting rasters...")
np.save(f"{OUTPUT_DIR}/lithology_v3.npy", lithology)
np.save(f"{OUTPUT_DIR}/basin_fill_grounded_v3.npy", basin_fill_grounded)
np.save(f"{OUTPUT_DIR}/schist_grade_v3.npy", schist_grade)
np.save(f"{OUTPUT_DIR}/jade_suitable_v3.npy", suitable_mask)
np.save(f"{OUTPUT_DIR}/jade_pods_v3.npy", pod_mask)

write_envi_raw(f"{OUTPUT_DIR}/lithology_v3", lithology.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 lithology v3 (Option A cross-section crest, height-normalized falloff): 0=ocean 1=basin_fill 2=greywacke 3=schist 4=volcanic",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/lithology_v3.prj", CRS_PROJ4)

write_envi_raw(f"{OUTPUT_DIR}/jade_pods_v3", pod_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 jade/pounamu pods v3: discrete stochastic pod footprints, NOT the full suitability zone",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/jade_pods_v3.prj", CRS_PROJ4)

meta = {
    "method": "Option A: cross-section snapping (DEM elevation argmax per perpendicular station) "
              "+ height-normalized falloff (low=0.6x base at lowest crest point, high=1.0x base at highest)",
    "cross_section_halfwidth_m": CROSS_SECTION_HALFWIDTH_M,
    "falloff_low_multiplier": FALLOFF_LOW_MULT,
    "falloff_high_multiplier": FALLOFF_HIGH_MULT,
    "crest_extraction": crest_info,
    "class_areas_km2": class_areas_km2,
    "basin_fill_grounded_km2": basin_fill_grounded_km2,
    "basin_fill_grounded_fraction": basin_fill_grounded_km2 / basin_fill_total_km2,
    "jade": {
        "n_pods": len(pod_centers_xy),
        "pod_centers_xy": [list(map(float, c)) for c in pod_centers_xy],
        "pod_radii_m": [float(r) for r in pod_radii_m],
        "pod_area_km2": jade_pod_area_km2,
        "suitability_zone_km2": suitable_area_km2,
        "grade_percentile_threshold": 80.0,
        "min_separation_km": 5.0,
        "seed": 13,
    },
}
with open(f"{OUTPUT_DIR}/tappa8_lithology_v3_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
