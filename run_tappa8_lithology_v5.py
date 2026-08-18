import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    _grid_xy, identify_landmasses, place_jade_pods, LITHOLOGY_CLASSES,
)
from geomorphology.terrain_relief import compute_local_relief, classify_from_terrain
from terrain.skeleton import load_geojson, build_zone_fields

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology v5 (DEM-NATIVE: elevation + local relief, NO ridge/zone geometry in the decision) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]

OUTPUT_DIR = "data/processed/geomorphology"

log("loading DEM...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape
land_mask = dem > 0

log("identifying landmasses (mainland vs SW Island -- still landmass id, not ridge geometry)...")
labeled, mainland_label, sw_island_label, landmass_areas_km2 = identify_landmasses(land_mask, cellsize_m=domain["resolution_m"])
mainland_mask = labeled == mainland_label
volcanic_mask = labeled == sw_island_label
log(f"mainland={landmass_areas_km2[mainland_label]:.1f} km2, sw_island(volcanic)={landmass_areas_km2[sw_island_label]:.1f} km2")

# PLACEHOLDERS, first pass -- not independently calibrated:
RELIEF_WINDOW_M = 2000.0     # windowed max-min elevation, 2km side
ELEV_PERCENTILE = 60.0       # mainland-land elevation percentile floor for schist eligibility (AND branch)
RELIEF_SCHIST_PERCENTILE = 75.0    # mainland-land relief percentile floor for schist (AND'd with elevation)
RELIEF_GREYWACKE_PERCENTILE = 50.0 # mainland-land relief percentile floor for greywacke (flank, not schist)
HIGH_ELEV_PERCENTILE = 90.0  # OR branch: elevation alone above this qualifies for schist regardless of
                               # relief -- fixes rounded-summit crest points that fail the relief test
                               # despite being the actual highest terrain (Nico's catch: 33-47% of the
                               # authored crest lines came out basin_fill without this). Checked it doesn't
                               # reopen the SE-plains fix: mainland p90 elevation (~2612m) sits above every
                               # authored plateau's own max (Central plateau tops out at 2463m).
ISLAND_RELIEF_PERCENTILE = 25.0  # SW Island (volcanic landmass) was unconditionally 100% volcanic through
                               # v5 -- a Tappa 7 S1 scope lock, not something this relief redesign touched
                               # until Nico's explicit request to extend the same logic to the island.
                               # Calibrated against the ISLAND's OWN elevation/relief population (island
                               # elevation range 0-737m sits below every mainland threshold, so reusing
                               # mainland thresholds would make the whole island basin_fill, wrong the other
                               # direction). Island land below this relief percentile (of island land only)
                               # -> basin_fill (flat coastal apron/isthmus reading, cf. Banks Peninsula's real
                               # flat margins); at/above -> stays volcanic. No greywacke/schist tier on the
                               # island (young monogenetic shield, not a metamorphic-flank story).
                               # UNREVIEWED first-pass placeholder -- not independently calibrated against
                               # any authored island zone composition target.

log(f"computing local relief (windowed max-min, window={RELIEF_WINDOW_M}m)...")
relief = compute_local_relief(dem, domain["resolution_m"], RELIEF_WINDOW_M)

log(f"classifying from terrain (elev>=p{ELEV_PERCENTILE} AND relief>=p{RELIEF_SCHIST_PERCENTILE} -> schist; "
    f"OR elev>=p{HIGH_ELEV_PERCENTILE} regardless of relief; "
    f"relief>=p{RELIEF_GREYWACKE_PERCENTILE} -> greywacke; else basin_fill)...")
result = classify_from_terrain(
    dem, land_mask, mainland_mask, volcanic_mask, relief,
    elev_percentile=ELEV_PERCENTILE,
    relief_schist_percentile=RELIEF_SCHIST_PERCENTILE,
    relief_greywacke_percentile=RELIEF_GREYWACKE_PERCENTILE,
    high_elev_percentile=HIGH_ELEV_PERCENTILE,
    island_relief_percentile=ISLAND_RELIEF_PERCENTILE,
)
lithology = result["lithology"]
schist_grade = result["schist_grade"]
log(f"resolved thresholds: elev>={result['elev_threshold_m']:.0f}m, "
    f"relief_schist>={result['relief_schist_threshold_m']:.0f}m, "
    f"relief_greywacke>={result['relief_greywacke_threshold_m']:.0f}m, "
    f"high_elev_or>={result['high_elev_threshold_m']:.0f}m "
    f"(recovered {result['n_schist_from_high_elev_only']} cells via the OR branch alone)")

cell_km2 = (domain["resolution_m"] / 1000.0) ** 2
class_areas_km2 = {name: float((lithology == code).sum() * cell_km2) for code, name in LITHOLOGY_CLASSES.items() if code != 0}
log(f"class areas km2: {class_areas_km2}")

# island diagnostics -- note the classify_from_terrain return uses "_km2_cells" naming but returns raw
# CELL COUNTS, not km2 (misleading name in the current terrain_relief.py; converting here)
island_volcanic_km2 = result["island_volcanic_km2_cells"] * cell_km2
island_basin_km2 = result["island_basin_km2_cells"] * cell_km2
island_total_km2 = island_volcanic_km2 + island_basin_km2
island_relief_threshold_m = result["island_relief_threshold_m"]
if island_relief_threshold_m is not None:
    log(f"island (SW Island, volcanic landmass): relief_threshold>={island_relief_threshold_m:.0f}m -> "
        f"volcanic={island_volcanic_km2:.1f} km2 ({100*island_volcanic_km2/island_total_km2:.1f}%), "
        f"basin_fill={island_basin_km2:.1f} km2 ({100*island_basin_km2/island_total_km2:.1f}%)")
else:
    log("island: island_relief_percentile=None, island left 100% volcanic (pre-extension behavior)")

log("checking basin_fill_grounded against authored plateau/plains zones (QA only, NOT a classification input)...")
zones = build_zone_fields(load_geojson("data/input/terrain_zones.geojson"))
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)
xy_flat = np.column_stack([xx.ravel(), yy.ravel()])
zone_any = np.zeros((ny, nx), dtype=bool)
zone_composition = {}
for zone in zones:
    if zone.feature_type != "plateau":
        continue
    inside = zone.path.contains_points(xy_flat).reshape(ny, nx) & land_mask
    zone_any |= inside
    n = int(inside.sum())
    if n == 0:
        continue
    comp = {name: float((lithology[inside] == code).sum() / n * 100.0) for code, name in LITHOLOGY_CLASSES.items()}
    zone_composition[zone.name] = comp
    log(f"  {zone.name}: {comp}")
from geomorphology.lithology import CLASS_BASIN_FILL
basin_fill_grounded = zone_any & (lithology == CLASS_BASIN_FILL)
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
np.save(f"{OUTPUT_DIR}/lithology_v5.npy", lithology)
np.save(f"{OUTPUT_DIR}/schist_grade_v5.npy", schist_grade)
np.save(f"{OUTPUT_DIR}/basin_fill_grounded_v5.npy", basin_fill_grounded)
np.save(f"{OUTPUT_DIR}/jade_suitable_v5.npy", suitable_mask)
np.save(f"{OUTPUT_DIR}/jade_pods_v5.npy", pod_mask)
np.save(f"{OUTPUT_DIR}/relief_2km.npy", relief)

write_envi_raw(f"{OUTPUT_DIR}/lithology_v5", lithology.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 lithology v5 (DEM-native: elevation + local relief, no ridge/zone geometry): 0=ocean 1=basin_fill 2=greywacke 3=schist 4=volcanic",
    dtype="i2")
write_prj(f"{OUTPUT_DIR}/lithology_v5.prj", CRS_PROJ4)

write_envi_raw(f"{OUTPUT_DIR}/jade_pods_v5", pod_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 jade/pounamu pods v5", dtype="i2")
write_prj(f"{OUTPUT_DIR}/jade_pods_v5.prj", CRS_PROJ4)

meta = {
    "method": "DEM-native terrain classification (elevation + local relief), NO ridge or zone authored geometry "
              "used in the classification decision -- see terrain_relief.py",
    "relief_window_m": RELIEF_WINDOW_M,
    "elev_percentile": ELEV_PERCENTILE,
    "relief_schist_percentile": RELIEF_SCHIST_PERCENTILE,
    "relief_greywacke_percentile": RELIEF_GREYWACKE_PERCENTILE,
    "high_elev_percentile": HIGH_ELEV_PERCENTILE,
    "resolved_thresholds": {
        "elev_threshold_m": result["elev_threshold_m"],
        "relief_schist_threshold_m": result["relief_schist_threshold_m"],
        "relief_greywacke_threshold_m": result["relief_greywacke_threshold_m"],
        "high_elev_threshold_m": result["high_elev_threshold_m"],
        "island_relief_threshold_m": island_relief_threshold_m,
    },
    "n_schist_from_high_elev_only": result["n_schist_from_high_elev_only"],
    "island_relief_percentile": ISLAND_RELIEF_PERCENTILE,
    "island_volcanic_km2": island_volcanic_km2,
    "island_basin_fill_km2": island_basin_km2,
    "island_volcanic_fraction": island_volcanic_km2 / island_total_km2 if island_total_km2 > 0 else None,
    "class_areas_km2": class_areas_km2,
    "basin_fill_grounded_km2": basin_fill_grounded_km2,
    "basin_fill_grounded_fraction": basin_fill_grounded_km2 / basin_fill_total_km2,
    "plateau_zone_composition_pct": zone_composition,
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
with open(f"{OUTPUT_DIR}/tappa8_lithology_v5_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
