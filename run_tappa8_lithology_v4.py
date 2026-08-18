import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    load_ridges, load_zones, classify, jade_eligible_mask, place_jade_pods,
    LITHOLOGY_CLASSES, CLASS_SCHIST,
)
from geomorphology.real_crest import extract_real_crest_network_local, height_normalized_falloff_m
from geomorphology.boundary_noise import boundary_warp_field
from geomorphology.lithology import _grid_xy
from terrain.skeleton import load_geojson

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology v4 (DEM ridge NETWORK as source geometry + height-normalized falloff + noise warp) ===")

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

# PLACEHOLDERS, first pass -- not independently calibrated:
CORRIDOR_HALFWIDTH_KM = 2.0        # scale that already avoided cross-ridge jump artifacts in v2/v3
WINDOW_SPACING_M = 1000.0          # along-line station spacing
WINDOW_HALFLENGTH_M = 750.0        # along-line window half-length (some overlap between adjacent stations)
ACCUM_PERCENTILE = 97.0            # LOCAL percentile within each station's window, NOT global over the whole
                                     # corridor -- see real_crest.py's extract_real_crest_network_local
                                     # docstring: a GLOBAL corridor-wide percentile (tried first) left 84% of
                                     # the Spine's 64km length with zero qualifying cells (a real bug, caught
                                     # by checking arclength coverage before this was shown to Nico), because
                                     # inverted-flow accumulation is an upstream-area proxy, not locally
                                     # uniform along a ridge trunk. Windowed+local fixes that; 97th percentile
                                     # (vs. an initial 90th) was chosen after checking crest_cells/corridor_land
                                     # ratio -- 90th let a 3x-overlapping window union swallow 47-88% of the
                                     # entire corridor as "crest", which would have made the schist band as
                                     # blocky as the corridor itself. 97th keeps it a much thinner real network
                                     # (~13% of corridor land for the Spine) while keeping zero coverage gaps.
FALLOFF_LOW_MULT = 0.6             # Nico's explicit values
FALLOFF_HIGH_MULT = 1.0
WARP_AMPLITUDE_M = 300.0           # noise displacement std, applied to ridge-distance query points only
WARP_BASE_WAVELENGTH_M = 800.0
WARP_OCTAVES = 3

log(f"extracting real crest NETWORK per ridge (DEM ridge line as source, windowed/local percentile, "
    f"corridor={CORRIDOR_HALFWIDTH_KM}km, percentile={ACCUM_PERCENTILE})...")
raw_ridge_features = load_geojson("data/input/terrain_ridges.geojson")
ridge_coords_by_name = {f["properties"]["name"]: np.array(f["geometry"]["coordinates"], dtype=np.float64) for f in raw_ridge_features}

real_crest_trees = {}
real_crest_falloff_m = {}
crest_info = {}
for ridge in ridges:
    tree, elev, info = extract_real_crest_network_local(
        ridge, ridge_coords_by_name[ridge.name], dem, ridge_accum_cells, land_mask,
        xmin=domain["xmin"], ymax=domain["ymax"], cellsize_m=domain["resolution_m"],
        corridor_halfwidth_km=CORRIDOR_HALFWIDTH_KM, window_spacing_m=WINDOW_SPACING_M,
        window_halflength_m=WINDOW_HALFLENGTH_M, accum_percentile_in_window=ACCUM_PERCENTILE,
    )
    crest_info[ridge.name] = info
    log(f"  {ridge.name}: n_crest_cells={info.get('n_crest_cells')} "
        f"n_stations={info.get('n_stations')} n_empty_windows={info.get('n_empty_windows')} "
        f"line_length_km={info.get('total_line_length_km'):.1f} status={info['status']}")
    if tree is not None:
        real_crest_trees[ridge.name] = tree
        falloff_arr = height_normalized_falloff_m(elev, ridge.falloff_km, FALLOFF_LOW_MULT, FALLOFF_HIGH_MULT)
        real_crest_falloff_m[ridge.name] = falloff_arr
        log(f"    elev range [{elev.min():.0f},{elev.max():.0f}]m -> "
            f"falloff range [{falloff_arr.min():.0f},{falloff_arr.max():.0f}]m "
            f"(base {ridge.falloff_km*1000:.0f}m)")

log(f"generating boundary noise (domain-warp of query points only, amplitude={WARP_AMPLITUDE_M}m, "
    f"wavelength={WARP_BASE_WAVELENGTH_M}m, {WARP_OCTAVES} octaves) -- see boundary_noise.py, "
    f"NOT a reuse of Tappa 1's own domain_warp (that code wasn't staged into this session)...")
warp_dx, warp_dy = boundary_warp_field(
    ny, nx, domain["resolution_m"], amplitude_m=WARP_AMPLITUDE_M,
    base_wavelength_m=WARP_BASE_WAVELENGTH_M, octaves=WARP_OCTAVES,
)

log("classifying lithology against REAL crest NETWORK (height-normalized falloff + noise warp)...")
result = classify(
    dem, ridges, zones,
    xmin=domain["xmin"], xmax=domain["xmax"], ymin=domain["ymin"], ymax=domain["ymax"],
    row_chunk_size=500, real_crest_trees=real_crest_trees, real_crest_falloff_m=real_crest_falloff_m,
    warp_dx=warp_dx, warp_dy=warp_dy,
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
np.save(f"{OUTPUT_DIR}/lithology_v4.npy", lithology)
np.save(f"{OUTPUT_DIR}/basin_fill_grounded_v4.npy", basin_fill_grounded)
np.save(f"{OUTPUT_DIR}/schist_grade_v4.npy", schist_grade)
np.save(f"{OUTPUT_DIR}/jade_suitable_v4.npy", suitable_mask)
np.save(f"{OUTPUT_DIR}/jade_pods_v4.npy", pod_mask)

write_envi_raw(f"{OUTPUT_DIR}/lithology_v4", lithology.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 lithology v4 (DEM ridge network source, height-normalized falloff, noise warp): 0=ocean 1=basin_fill 2=greywacke 3=schist 4=volcanic",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/lithology_v4.prj", CRS_PROJ4)

write_envi_raw(f"{OUTPUT_DIR}/jade_pods_v4", pod_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 jade/pounamu pods v4: discrete stochastic pod footprints, NOT the full suitability zone",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/jade_pods_v4.prj", CRS_PROJ4)

meta = {
    "method": "DEM ridge NETWORK as source geometry (extract_real_crest_network: percentile-thresholded "
              "inverted-flow-accum cells inside a narrow corridor around the authored line) "
              "+ height-normalized falloff (0.6x-1.0x) + noise domain-warp on query points",
    "corridor_halfwidth_km": CORRIDOR_HALFWIDTH_KM,
    "accum_percentile_in_corridor": ACCUM_PERCENTILE,
    "falloff_low_multiplier": FALLOFF_LOW_MULT,
    "falloff_high_multiplier": FALLOFF_HIGH_MULT,
    "warp_amplitude_m": WARP_AMPLITUDE_M,
    "warp_base_wavelength_m": WARP_BASE_WAVELENGTH_M,
    "warp_octaves": WARP_OCTAVES,
    "warp_note": "independent value-noise implementation for this stage, NOT Tappa 1's actual domain_warp code "
                 "(terrain.generate/terrain.erosion were never staged into this session)",
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
with open(f"{OUTPUT_DIR}/tappa8_lithology_v4_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
