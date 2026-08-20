import sys, time, json
sys.path.insert(0, 'src')
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import LITHOLOGY_CLASSES, load_authoral_zones, apply_authoral_zones

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- lithology v6 (v5 DEM-native base + hand-authored marble/"
    "sedimentary_limestone/granite zones, priority-tiered composite) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUTPUT_DIR = "data/processed/geomorphology"

log("loading v5 lithology (base) + authored zones...")
lithology_v5 = np.load(f"{OUTPUT_DIR}/lithology_v5.npy")
zones = load_authoral_zones("data/input/lithology_authoral.geojson")
log(f"  {len(zones)} authored zones: " +
    ", ".join(f"{z['name']!r}({z['feature_type']}, rank{z['priority_rank']})" for z in zones))

log("compositing (weakest rank painted first, rank 1 last)...")
lithology_v6, per_zone_stats = apply_authoral_zones(
    lithology_v5, zones, domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], RES_M,
)

for s in per_zone_stats:
    log(f"  {s['name']!r} ({s['feature_type']}, rank {s['priority_rank']}, "
        f"grounded_claimed={s['grounded_claimed']}): nominal={s['nominal_polygon_area_km2']:.2f} km2, "
        f"painted={s['painted_area_km2']:.2f} km2, pre-authoral composition={s['pre_authoral_composition_pct']}")

cell_km2 = (RES_M / 1000.0) ** 2
class_areas_km2 = {name: float((lithology_v6 == code).sum() * cell_km2)
                    for code, name in LITHOLOGY_CLASSES.items() if code != 0}
log(f"class areas km2 (v6): {class_areas_km2}")

log("exporting rasters...")
np.save(f"{OUTPUT_DIR}/lithology_v6.npy", lithology_v6)
write_envi_raw(f"{OUTPUT_DIR}/lithology_v6", lithology_v6.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=("Tappa 8 lithology v6 (v5 DEM-native base + hand-authored marble/sedimentary_limestone/"
                  "granite): 0=ocean 1=basin_fill 2=greywacke 3=schist 4=volcanic 5=marble "
                  "6=sedimentary_limestone 7=granite"),
    dtype="u1")
write_prj(f"{OUTPUT_DIR}/lithology_v6.prj", CRS_PROJ4)

meta = {
    "method": ("v5 DEM-native base (elevation+relief, unchanged) + hand-authored zone overlay from "
               "data/input/lithology_authoral.geojson, composited by apply_authoral_zones() using each "
               "zone's priority_rank tier (see lithology.py PRIORITY_RANK_BEATS)."),
    "authored_zones": per_zone_stats,
    "class_areas_km2": class_areas_km2,
}
with open(f"{OUTPUT_DIR}/tappa8_lithology_v6_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
