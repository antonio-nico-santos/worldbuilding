import sys, time, json, os
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    load_ridges, load_zones, classify, jade_eligible_mask,
    LITHOLOGY_CLASSES,
)

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]

OUTPUT_DIR = "data/processed/geomorphology"
os.makedirs(OUTPUT_DIR, exist_ok=True)

log("loading DEM (native 30m, v3 eroded)...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape
log(f"DEM shape={dem.shape} land_fraction={100*(dem>0).mean():.2f}%")

log("loading skeleton (ridges + plateau zones)...")
ridges = load_ridges("data/input/terrain_ridges.geojson")
zones = load_zones("data/input/terrain_zones.geojson")
log(f"ridges: {[r.name for r in ridges]}")
log(f"zones (plateau only counted): {[z.name for z in zones if z.feature_type == 'plateau']}")

log("classifying lithology (chunked, this is the long step)...")
result = classify(
    dem, ridges, zones,
    xmin=domain["xmin"], xmax=domain["xmax"], ymin=domain["ymin"], ymax=domain["ymax"],
    row_chunk_size=500,
)
lithology = result["lithology"]
basin_fill_grounded = result["basin_fill_grounded"]
schist_grade = result["schist_grade"]

log(f"landmass areas (km2, top 5): {sorted(result['landmass_areas_km2'], reverse=True)[:5]}")
log(f"mainland_label={result['mainland_label']} sw_island_label={result['sw_island_label']}")
sw_island_area_km2 = float(result["landmass_areas_km2"][result["sw_island_label"]])
mainland_area_km2 = float(result["landmass_areas_km2"][result["mainland_label"]])

assert not np.isnan(lithology).any()
cell_km2 = (domain["resolution_m"] / 1000.0) ** 2
land_mask = dem > 0
land_km2 = float(land_mask.sum() * cell_km2)

class_areas_km2 = {}
for code, name in LITHOLOGY_CLASSES.items():
    if code == 0:
        continue
    area = float((lithology == code).sum() * cell_km2)
    class_areas_km2[name] = area

basin_fill_total_km2 = class_areas_km2["sedimentary_basin_fill"]
basin_fill_grounded_km2 = float(basin_fill_grounded.sum() * cell_km2)
basin_fill_leftover_km2 = basin_fill_total_km2 - basin_fill_grounded_km2

log(f"class areas km2: {class_areas_km2}")
log(f"basin fill: grounded={basin_fill_grounded_km2:.1f} leftover={basin_fill_leftover_km2:.1f} "
    f"({100*basin_fill_grounded_km2/basin_fill_total_km2:.1f}% grounded)")

log("computing jade/pounamu eligible mask (top 20% schist grade)...")
jade_mask = jade_eligible_mask(lithology, schist_grade, grade_percentile=80.0)
jade_area_km2 = float(jade_mask.sum() * cell_km2)
log(f"jade-eligible area: {jade_area_km2:.2f} km2 "
    f"({100*jade_area_km2/class_areas_km2['schist']:.1f}% of schist)")

log("exporting rasters (npy + ENVI int16/int8 + prj)...")
np.save(f"{OUTPUT_DIR}/lithology.npy", lithology)
np.save(f"{OUTPUT_DIR}/basin_fill_grounded.npy", basin_fill_grounded)
np.save(f"{OUTPUT_DIR}/schist_grade.npy", schist_grade)
np.save(f"{OUTPUT_DIR}/jade_eligible.npy", jade_mask)

write_envi_raw(
    f"{OUTPUT_DIR}/lithology", lithology.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 lithology: 0=ocean 1=sedimentary_basin_fill 2=greywacke_argillite "
                 "3=schist 4=volcanic (Fictional World LCC domain)",
    dtype="i2",
)
write_prj(f"{OUTPUT_DIR}/lithology.prj", CRS_PROJ4)

write_envi_raw(
    f"{OUTPUT_DIR}/jade_eligible", jade_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=domain["resolution_m"],
    description="Tappa 8 jade/pounamu eligible mask (top 20th percentile schist grade)",
    dtype="i2",
)
write_prj(f"{OUTPUT_DIR}/jade_eligible.prj", CRS_PROJ4)

meta = {
    "resolution_m": domain["resolution_m"],
    "shape": list(lithology.shape),
    "land_km2": land_km2,
    "mainland_area_km2": mainland_area_km2,
    "sw_island_area_km2": sw_island_area_km2,
    "class_areas_km2": class_areas_km2,
    "basin_fill_grounded_km2": basin_fill_grounded_km2,
    "basin_fill_leftover_km2": basin_fill_leftover_km2,
    "basin_fill_grounded_fraction": basin_fill_grounded_km2 / basin_fill_total_km2,
    "jade_eligible_km2": jade_area_km2,
    "jade_grade_percentile_threshold": 80.0,
    "shelf_multipliers_used": {r.name: r.shelf_multiplier for r in ridges},
    "ridges": [{"name": r.name, "falloff_km": r.falloff_km, "shelf_multiplier": r.shelf_multiplier} for r in ridges],
    "plateau_zones": [z.name for z in zones if z.feature_type == "plateau"],
}
with open(f"{OUTPUT_DIR}/tappa8_lithology_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
