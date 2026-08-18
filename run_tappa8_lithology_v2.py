import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    load_ridges, load_zones, classify, jade_eligible_mask, place_jade_pods,
    LITHOLOGY_CLASSES, CLASS_SCHIST,
)
from geomorphology.real_crest import extract_real_crest
from geomorphology.lithology import _grid_xy
from terrain.skeleton import load_geojson

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology v2 (real DEM-grounded crest, no elevation-falloff modulation) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]

OUTPUT_DIR = "data/processed/geomorphology"

log("loading DEM, skeleton, ridge_accum...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape
land_mask = dem > 0
ridges = load_ridges("data/input/terrain_ridges.geojson")
zones = load_zones("data/input/terrain_zones.geojson")
ridge_accum_cells = np.load(f"{OUTPUT_DIR}/ridge_accum_cells.npy")

log("building full grid coordinates...")
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

log("extracting real crest per ridge (arclength-binned, 1km bins, local-window "
    "argmax of inverted-flow-accum within each bin's own falloff_km radius)...")
raw_ridge_features = load_geojson("data/input/terrain_ridges.geojson")
ridge_coords_by_name = {f["properties"]["name"]: np.array(f["geometry"]["coordinates"], dtype=np.float64) for f in raw_ridge_features}

# PLACEHOLDER, first pass: search_radius_km=2.5 fixed for every ridge, NOT
# ridge.falloff_km. Diagnosed after the first full run (radius=falloff_km,
# 15-24km) produced occasional huge jumps between consecutive bin points
# (Spine max step 15.4km, South Branch max step 20.7km) -- a window that
# wide can and does snap onto a real terrain high point belonging to a
# totally different spur, not the named massif's own crest. The union-of-
# circles buffer between two such widely separated points draws a straight
# tangent chord, which is *visually worse* than v1's smooth curve (sharp
# angular facets instead of a rounded artificial curve -- still artificial,
# differently). Shrinking the search window to a few km keeps the "real
# DEM texture" idea while forcing the point sequence to stay coherent along
# the authored line. 2.5km is not calibrated against anything -- it is the
# first value tried that killed the >5km jumps (max step dropped to ~4.3km
# across all 5 ridges) without creating empty bins. Flagging for review.
CREST_SEARCH_RADIUS_KM = 2.5
real_crest_trees = {}
crest_info = {}
for ridge in ridges:
    tree, elev, info = extract_real_crest(
        ridge, ridge_coords_by_name[ridge.name], dem, ridge_accum_cells, land_mask,
        xmin=domain["xmin"], ymax=domain["ymax"], cellsize_m=domain["resolution_m"],
        bin_spacing_m=1000.0, search_radius_km=CREST_SEARCH_RADIUS_KM,
    )
    crest_info[ridge.name] = info
    log(f"  {ridge.name}: {info}")
    if tree is not None:
        real_crest_trees[ridge.name] = tree

log("classifying lithology against REAL crest (unmodulated falloff/shelf)...")
result = classify(
    dem, ridges, zones,
    xmin=domain["xmin"], xmax=domain["xmax"], ymin=domain["ymin"], ymax=domain["ymax"],
    row_chunk_size=500, real_crest_trees=real_crest_trees,
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
np.save(f"{OUTPUT_DIR}/lithology_v2.npy", lithology)
np.save(f"{OUTPUT_DIR}/basin_fill_grounded_v2.npy", basin_fill_grounded)
np.save(f"{OUTPUT_DIR}/schist_grade_v2.npy", schist_grade)
np.save(f"{OUTPUT_DIR}/jade_suitable_v2.npy", suitable_mask)
np.save(f"{OUTPUT_DIR}/jade_pods_v2.npy", pod_mask)

write_envi_raw(f"{OUTPUT_DIR}/lithology_v2", lithology.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 lithology v2 (real DEM-grounded crest): 0=ocean 1=basin_fill 2=greywacke 3=schist 4=volcanic",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/lithology_v2.prj", CRS_PROJ4)

write_envi_raw(f"{OUTPUT_DIR}/jade_pods_v2", pod_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 jade/pounamu pods v2: discrete stochastic pod footprints, NOT the full suitability zone",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/jade_pods_v2.prj", CRS_PROJ4)

meta = {
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
with open(f"{OUTPUT_DIR}/tappa8_lithology_v2_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
