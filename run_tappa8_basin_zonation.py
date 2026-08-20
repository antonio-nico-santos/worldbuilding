"""
Tappa 8 -- basin fill internal zonation (Arable / Wetland-backswamp /
Estuarine-coastal), implementing scenario_reference.md S22.4's logic as an
actual raster. See geomorphology/basin_zonation.py for the full reasoning;
this script computes the two input fields (relief_2km, dist_to_ocean_km),
runs the classifier, and reports results.

Also answers a SECOND, unrelated open item while the inputs are already in
hand (cheap, not the main point of this script): scenario_reference.md
S22.3/S22.19 flags "basin fill / Grassland-biome spatial overlap... a
hypothesis, not yet checked against the rasters" -- checked directly here.

Reads:
  data/processed/dem_v3_final_30m_eroded.npy       (native 30m)
  data/processed/geomorphology/lithology_v6.npy    (native 30m)
  data/processed/climate/land_mask.npy             (120m, shape reference
                                                     for the biome overlap
                                                     check only)
  data/processed/biomes/biome_id.npy               (120m)

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  basin_zonation_30m.*   (uint8 raster, native 30m -- same resolution as
                          lithology_v6, no cost-graph dependency so no
                          reason to downsample)
  tappa8_basin_zonation_meta.json
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_BASIN_FILL, LITHOLOGY_CLASSES
from geomorphology.terrain_relief import compute_local_relief
from geomorphology.caves import distance_to_ocean_km
from geomorphology.basin_zonation import classify_basin_fill_zones, BASIN_SUBZONE_CLASSES
from suitability.terrain_metrics import block_mode
from biomes.world_biomes import BIOME_NAMES

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- basin fill zonation (Arable / Wetland-backswamp / Estuarine-coastal, "
    "implementing scenario_reference.md S22.4) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUT = "data/processed/geomorphology"

log("loading DEM (30m) + lithology_v6 (30m)...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask30 = dem > 0
lithology_v6 = np.load(f"{OUT}/lithology_v6.npy")
basin_fill_mask = lithology_v6 == CLASS_BASIN_FILL
cell_km2 = (RES_M / 1000.0) ** 2
log(f"  basin_fill: {float(basin_fill_mask.sum() * cell_km2):.1f} km2")

log("computing relief_2km (same field/window as lithology v5's own classify_from_terrain) "
    "and dist_to_ocean_km (same field sea caves / placer magnetite already use)...")
relief_2km = compute_local_relief(dem, RES_M, window_m=2000.0)
dist_to_ocean_km = distance_to_ocean_km(land_mask30, RES_M)

log("checking basin_fill's own dist_to_ocean_km / relief_2km distributions before "
    "picking thresholds (calibrate against this world's own data, not asserted a priori)...")
bf_dist_to_ocean = dist_to_ocean_km[basin_fill_mask]
bf_relief = relief_2km[basin_fill_mask]
for pct in (10, 25, 50, 75, 90):
    log(f"  dist_to_ocean_km p{pct}: {np.percentile(bf_dist_to_ocean, pct):.2f} km   "
        f"relief_2km p{pct}: {np.percentile(bf_relief, pct):.1f} m")

COASTAL_THRESHOLD_KM = 1.5
WETLAND_RELIEF_PERCENTILE = 25.0
log(f"  using coastal_threshold_km={COASTAL_THRESHOLD_KM} (absolute, first-pass), "
    f"wetland_relief_percentile={WETLAND_RELIEF_PERCENTILE} (of NON-coastal basin_fill's "
    f"own relief distribution)")

subzone, stats = classify_basin_fill_zones(
    basin_fill_mask, relief_2km, dist_to_ocean_km,
    coastal_threshold_km=COASTAL_THRESHOLD_KM,
    wetland_relief_percentile=WETLAND_RELIEF_PERCENTILE,
)

area_km2 = {name: float((subzone == code).sum() * cell_km2) for code, name in BASIN_SUBZONE_CLASSES.items()}
log(f"  sub-zone areas (km2): {area_km2}")
total_bf_km2 = area_km2["arable"] + area_km2["wetland_backswamp"] + area_km2["estuarine_coastal"]
log(f"  arable={100*area_km2['arable']/total_bf_km2:.1f}% "
    f"wetland={100*area_km2['wetland_backswamp']/total_bf_km2:.1f}% "
    f"estuarine={100*area_km2['estuarine_coastal']/total_bf_km2:.1f}% of basin_fill")

# --- bonus: basin_fill / Grassland biome overlap (S22.3/S22.19, unrelated open item) ----
log("bonus check -- basin_fill / Grassland biome spatial overlap (S22.3/S22.19, flagged as "
    "unverified hypothesis)...")
land120 = np.load("data/processed/climate/land_mask.npy").astype(bool)
biome_id = np.load("data/processed/biomes/biome_id.npy")
ny120, nx120 = land120.shape
n_classes = len(LITHOLOGY_CLASSES)
lithology120 = block_mode(lithology_v6, 4, n_classes)[:ny120, :nx120]
basin_fill_120 = lithology120 == CLASS_BASIN_FILL

grassland_code = BIOME_NAMES.index("Lowland Steppe / Grassland")
grassland_mask = biome_id == grassland_code

bf_area_120_km2 = float(basin_fill_120.sum()) * (RES_M * 4 / 1000.0) ** 2  # ~120m cell
overlap_km2 = float((basin_fill_120 & grassland_mask).sum()) * (RES_M * 4 / 1000.0) ** 2
grassland_area_km2 = float(grassland_mask.sum()) * (RES_M * 4 / 1000.0) ** 2
pct_of_basin_fill_in_grassland = 100 * overlap_km2 / bf_area_120_km2 if bf_area_120_km2 > 0 else 0.0
pct_of_grassland_in_basin_fill = 100 * overlap_km2 / grassland_area_km2 if grassland_area_km2 > 0 else 0.0
log(f"  basin_fill (120m): {bf_area_120_km2:.1f} km2; Grassland biome: {grassland_area_km2:.1f} km2; "
    f"overlap: {overlap_km2:.1f} km2 ({pct_of_basin_fill_in_grassland:.1f}% of basin_fill, "
    f"{pct_of_grassland_in_basin_fill:.1f}% of Grassland)")

# --- export -------------------------------------------------------------------------
log("exporting basin zonation raster (native 30m, uint8) + meta...")
np.save(f"{OUT}/basin_zonation_30m.npy", subzone)
write_envi_raw(
    f"{OUT}/basin_zonation_30m", subzone,
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=(
        "Tappa 8 basin fill internal zonation (scenario_reference.md S22.4, calibrated "
        "here): 0=not_basin_fill 1=arable 2=wetland_backswamp 3=estuarine_coastal. Derived "
        "from relief_2km + dist_to_ocean_km, not hand-authored. UNREVIEWED first-pass "
        "thresholds, not locked."
    ),
    dtype="u1",
)
write_prj(f"{OUT}/basin_zonation_30m.prj", CRS_PROJ4)

meta = {
    "scope_note": (
        "Implements scenario_reference.md S22.4's three-way basin_fill sub-zonation "
        "(Arable / Wetland-backswamp / Estuarine-coastal) as a raster -- that section "
        "explicitly left exact thresholds for Tappa 8 to calibrate; the LOGIC (which "
        "fields, coastal-claims-priority ordering) is S22.4's, not this script's."
    ),
    "subzone_classes": BASIN_SUBZONE_CLASSES,
    "calibration": stats,
    "basin_fill_own_distributions": {
        "dist_to_ocean_km_percentiles": {
            str(p): float(np.percentile(bf_dist_to_ocean, p)) for p in (10, 25, 50, 75, 90)
        },
        "relief_2km_m_percentiles": {
            str(p): float(np.percentile(bf_relief, p)) for p in (10, 25, 50, 75, 90)
        },
    },
    "area_km2": area_km2,
    "pct_of_basin_fill": {
        "arable": round(100 * area_km2["arable"] / total_bf_km2, 2),
        "wetland_backswamp": round(100 * area_km2["wetland_backswamp"] / total_bf_km2, 2),
        "estuarine_coastal": round(100 * area_km2["estuarine_coastal"] / total_bf_km2, 2),
    },
    "threshold_status": "UNREVIEWED first-pass -- coastal_threshold_km=3.0 is an absolute "
    "distance judgement call (no percentile basis, unlike relief); "
    "wetland_relief_percentile=25.0 was chosen to keep Arable the 'clear majority' S22.4 "
    "specifies, sanity-checked against that one constraint, not independently derived. "
    "Pending Nico's sign-off; not written to config/parameters.yml.",
    "bonus_check_basin_fill_grassland_overlap": {
        "note": "Answers scenario_reference.md S22.3/S22.19's separate, previously-"
        "unverified hypothesis -- computed at 120m (biome grid resolution), basin_fill "
        "from lithology_v6 block-mode downsampled.",
        "basin_fill_area_km2_120m": bf_area_120_km2,
        "grassland_biome_area_km2": grassland_area_km2,
        "overlap_km2": overlap_km2,
        "pct_of_basin_fill_in_grassland": round(pct_of_basin_fill_in_grassland, 2),
        "pct_of_grassland_in_basin_fill": round(pct_of_grassland_in_basin_fill, 2),
    },
    "not_done": "vivianite's spatial placement (S8e) still uses its original wetness-proxy "
    "stand-in (1/(1+dist_to_stream_km)) -- re-pointing it at this zonation's "
    "wetland_backswamp sub-class specifically is a natural follow-up, not done in this "
    "script (would change resource_vivianite.npy's already-committed output).",
}
with open(f"{OUT}/tappa8_basin_zonation_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
