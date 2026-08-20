"""
Tappa 8 -- spatial pod placement for the four Vertice materials that were only ever a
non-spatial per-class lookup (resources.py's VERTICE_MATERIALS), generalizing jade's
own place_jade_pods() method (Nico's request: "map other resources possible location
the same way we did with jade").

Deliberately NOT all six materials get a new pod raster here -- two are excluded on
purpose, not by oversight (see the log lines and meta.json for the reasoning):
- schist's gold-bearing quartz veins / muscovite mica: resources.py already says
  these co-locate with the jade/pounamu high-grade band ("spatial_note"). Placing
  them independently would contradict that note. They reuse jade_pods_v5.npy
  directly -- no new raster.
- volcanic's magnetite (primary): titanomagnetite in basalt is a disseminated bulk
  mineral (a few percent through the rock), not a vein/joint/bog-localized
  phenomenon the way the other four are. Treating it as discrete pods would be LESS
  geologically honest than what resources.py already does (a class-level fact) --
  so it stays a class-level fact, no new raster.

The four that DO get pods, and what each is weighted by:
- laumontite (greywacke, primary): UNIFORM. No citation supports a within-greywacke
  gradient (unlike schist_grade's real metamorphic-grade rank story) -- inventing
  one just to have a weight field would be worse than admitting none exists.
- vivianite (basin_fill, primary): weighted by a wetness proxy (inverse distance to
  stream). Citation is "low-oxygen, organic-rich floodplain/bog sediment" --
  genuinely a wetland-adjacency fact, not uniform across all 4437 km2 of basin_fill.
  This is a LIGHTWEIGHT stand-in for the fuller basin-fill zonation (Arable/Wetland/
  Estuarine, scenario_reference.md S22.4) discussed but not yet built -- flagged as
  provisional, to be revisited once that zonation actually exists as a raster.
- reworked placer magnetite (basin_fill, secondary_weak): weighted by coastal
  proximity x proximity to the volcanic landmass. CORRECTS a geological
  inconsistency caught while doing this spatial pass: resources.py's existing text
  says "reconcentrated by rivers", but mainland basin_fill and the volcanic zone are
  on two different landmasses separated by 16.44 km of open water (Tappa 7 S1) --
  there is no shared river catchment for fluvial reworking to happen through. Real
  NZ ironsand beaches (Taranaki, Westport) actually form by COASTAL/longshore
  redistribution of eroded volcanic material, not river transport -- a mechanism
  that DOES work across the water gap. resources.py's citation text should be
  corrected to this (not done automatically here -- flagging for Nico's sign-off,
  see the decision doc).
- native silver / native copper (volcanic, secondary): weighted by vent_weight (the
  SAME Gaussian vent-proximity field lava tubes already use, from
  geothermal.geojson) -- epithermal veins are a real, direct hydrothermal/volcanic-
  vent association, so reusing this existing field is a clean fit, not a stretch.
  Each placed pod is independently rolled silver (rarer) vs copper (common) --
  ratio is an UNREVIEWED placeholder (25/75), no citation, flagged the same way as
  every other first-pass numeric choice in this project.
"""
import sys, time, json
sys.path.insert(0, 'src')
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    CLASS_GREYWACKE, CLASS_BASIN_FILL, CLASS_VOLCANIC, _grid_xy, place_material_pods,
)
from geomorphology.caves import vent_proximity_weight, distance_to_ocean_km
from terrain.skeleton import load_geojson

t_start = time.time()
def log(msg):
    print(f"[{time.time()-t_start:7.1f}s] {msg}", flush=True)

log("=== Tappa 8 -- resource pod placement (laumontite, vivianite, placer magnetite, "
    "silver/copper), generalizing place_jade_pods() to the other four Vertice materials ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUTPUT_DIR = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

log("loading DEM + lithology v6 + stream_mask + geothermal vents...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
ny, nx = dem.shape
lithology_v6 = np.load(f"{OUTPUT_DIR}/lithology_v6.npy")
stream_mask = np.load("data/processed/hydrology/stream_mask.npy").astype(bool)
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

log("computing shared fields (dist_to_stream, dist_to_ocean, dist_to_volcanic, vent_weight)...")
dist_to_stream_km = ndimage.distance_transform_edt(~stream_mask, sampling=(RES_M, RES_M)) / 1000.0
dist_to_ocean_km = distance_to_ocean_km(land_mask, RES_M)
volcanic_mask = lithology_v6 == CLASS_VOLCANIC
dist_to_volcanic_km = ndimage.distance_transform_edt(~volcanic_mask, sampling=(RES_M, RES_M)) / 1000.0

vents_geojson = load_geojson("data/input/geothermal.geojson")
vents = [
    (feat["geometry"]["coordinates"][0], feat["geometry"]["coordinates"][1], feat["properties"]["falloff_km"])
    for feat in vents_geojson
]
vent_weight = vent_proximity_weight(xx, yy, vents, RES_M)

# shared placeholder pod parameters -- same status/rationale as jade's own (S5 of the decision
# doc): first-pass, not independently calibrated, reusing jade's 5 km min-separation scale.
N_PODS = 8
MIN_SEPARATION_KM = 5.0
RADIUS_RANGE_M = (300.0, 800.0)
SEED_BASE = 130  # distinct from jade's seed=13, still in the same "documented, not random" spirit

results = {}

log("--- laumontite (greywacke, UNIFORM -- no citable within-class gradient) ---")
greywacke_mask = lithology_v6 == CLASS_GREYWACKE
pod_mask, centers, radii = place_material_pods(
    greywacke_mask, None, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M, seed=SEED_BASE + 1,
)
results["laumontite"] = {"pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii,
                          "eligible_km2": float(greywacke_mask.sum() * cell_km2), "weight": "uniform"}
log(f"  {len(centers)} pods placed, pod area={float(pod_mask.sum()*cell_km2):.2f} km2 "
    f"within {results['laumontite']['eligible_km2']:.1f} km2 eligible greywacke")

log("--- vivianite (basin_fill primary, weighted by wetness proxy = 1/(1+dist_to_stream_km)) ---")
basin_fill_mask = lithology_v6 == CLASS_BASIN_FILL
wetness = 1.0 / (1.0 + dist_to_stream_km)
pod_mask, centers, radii = place_material_pods(
    basin_fill_mask, wetness, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M, seed=SEED_BASE + 2,
)
results["vivianite"] = {"pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii,
                         "eligible_km2": float(basin_fill_mask.sum() * cell_km2), "weight": "wetness_proxy"}
log(f"  {len(centers)} pods placed, pod area={float(pod_mask.sum()*cell_km2):.2f} km2")

log("--- placer magnetite (basin_fill secondary_weak, weighted by coastal x volcanic proximity "
    "-- CORRECTS resources.py's river-transport citation, see module docstring) ---")
coastal_island_weight = (1.0 / (1.0 + dist_to_ocean_km)) * (1.0 / (1.0 + dist_to_volcanic_km))
pod_mask, centers, radii = place_material_pods(
    basin_fill_mask, coastal_island_weight, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M, seed=SEED_BASE + 3,
)
results["placer_magnetite"] = {"pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii,
                                "eligible_km2": results["vivianite"]["eligible_km2"], "weight": "coastal_x_volcanic_proximity"}
log(f"  {len(centers)} pods placed, pod area={float(pod_mask.sum()*cell_km2):.2f} km2")

log("--- silver/copper (volcanic secondary, weighted by vent_weight -- same field lava tubes use) ---")
pod_mask, centers, radii = place_material_pods(
    volcanic_mask, vent_weight, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M, seed=SEED_BASE + 4,
)
rng = np.random.default_rng(SEED_BASE + 40)
SILVER_FRACTION = 0.25  # unreviewed placeholder, no citation for this specific ratio
pod_kinds = ["silver" if rng.random() < SILVER_FRACTION else "copper" for _ in centers]
results["silver_copper"] = {"pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii,
                             "eligible_km2": float(volcanic_mask.sum() * cell_km2), "weight": "vent_weight",
                             "pod_kinds": pod_kinds}
log(f"  {len(centers)} pods placed ({pod_kinds.count('silver')} silver, {pod_kinds.count('copper')} copper), "
    f"pod area={float(pod_mask.sum()*cell_km2):.2f} km2")

log("gold-bearing quartz veins / muscovite mica (schist): NOT independently placed -- "
    "resources.py's own spatial_note says these co-locate with jade_pods_v5.npy. No new raster.")
log("magnetite primary (volcanic): NOT independently placed -- disseminated bulk basalt mineral, "
    "not vein/joint-hosted like the other four. Stays a class-level fact in resources.py. No new raster.")

log("exporting rasters...")
NAME_MAP = {"laumontite": "resource_laumontite", "vivianite": "resource_vivianite",
            "placer_magnetite": "resource_placer_magnetite", "silver_copper": "resource_silver_copper"}
for key, out_name in NAME_MAP.items():
    r = results[key]
    np.save(f"{OUTPUT_DIR}/{out_name}.npy", r["pod_mask"])
    write_envi_raw(f"{OUTPUT_DIR}/{out_name}", r["pod_mask"].astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 resource pods: {key} (weighted by {r['weight']})",
        dtype="u1")
    write_prj(f"{OUTPUT_DIR}/{out_name}.prj", CRS_PROJ4)

meta = {
    "method": "place_material_pods() (lithology.py) -- generalizes place_jade_pods() to an "
              "arbitrary eligibility mask + optional weight field. n_pods, min_separation_km, "
              "radius_range_m are first-pass placeholders, same status as jade's own.",
    "n_pods_target": N_PODS,
    "min_separation_km": MIN_SEPARATION_KM,
    "radius_range_m": list(RADIUS_RANGE_M),
    "excluded_from_this_pass": {
        "gold_quartz_mica_schist": "co-locates with jade_pods_v5.npy per resources.py's own spatial_note -- no new raster",
        "magnetite_primary_volcanic": "disseminated bulk basalt mineral, not vein/joint-hosted -- stays class-level, no new raster",
    },
    "citation_correction_flagged_not_yet_applied_to_resources_py": {
        "placer_magnetite": "resources.py currently says 'reconcentrated by rivers', but mainland "
            "basin_fill and the volcanic zone are separate landmasses (16.44 km water gap, Tappa 7 S1) "
            "with no shared river catchment. Real NZ ironsand beaches form by coastal/longshore "
            "redistribution, not river transport -- this pass weighted placement by coastal x "
            "volcanic-landmass proximity accordingly, but resources.py's own citation text has NOT "
            "been edited yet, pending Nico's sign-off.",
    },
    "materials": {
        key: {
            "weight_basis": r["weight"],
            "eligible_km2": r["eligible_km2"],
            "n_pods_placed": len(r["centers_xy"]),
            "pod_area_km2": float(r["pod_mask"].sum() * cell_km2),
            "pod_centers_xy": [list(map(float, c)) for c in r["centers_xy"]],
            "pod_radii_m": [float(x) for x in r["radii_m"]],
            **({"pod_kinds": r["pod_kinds"], "silver_fraction_placeholder": SILVER_FRACTION} if "pod_kinds" in r else {}),
        }
        for key, r in results.items()
    },
}
with open(f"{OUTPUT_DIR}/tappa8_resource_pods_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
