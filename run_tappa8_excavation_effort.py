"""
Tappa 8 -- excavation effort (Nico's follow-up to the transport friction
multiplier, S8f): a SEPARATE relative index for "how much harder is it to
build on/dig into this ground", after Nico clarified his original
"construction" framing wasn't asking for a full engineering-cost model
(excavation crew-hours, blasting budgets, haulage) -- just a relative,
per-class number, "how easier it is to work on basins than schist or
granite."

This is deliberately a MUCH simpler script than
run_tappa8_transport_friction.py: excavation effort is a pure per-cell
function of lithology class (no Dijkstra graph, no direction, no DEM/lake/
boat-speed inputs at all) -- see geomorphology/lithology.py's
EXCAVATION_EFFORT_MULTIPLIER / excavation_effort_multiplier() for the
actual values, real-world grounding (relative mineral hardness + historical
stoneworking accounts), and an important caveat repeated here because it
matters: this does NOT capture karst's foundation-STABILITY risk (sinkholes/
voids under marble and sedimentary_limestone) -- that's a separate hazard
from raw workability. Taken alone, this index would make marble/limestone
look like easy building ground, which is only true for the quarrying half
of the question.

Computed at the SAME 120m grid as S8f's transport friction (block-mode
downsample of lithology_v6, reusing terrain_metrics.block_mode) purely for
consistency/easy overlay in QGIS -- not because anything here needs Tappa
6's grid specifically; unlike friction, this has no cost-graph dependency
and no reason it couldn't be computed at lithology_v6's native 30m instead.

Reads:
  data/processed/climate/land_mask.npy          (shape reference only, for
                                                 the same [:ny,:nx] cropping
                                                 convention S8f uses)
  data/processed/geomorphology/lithology_v6.npy  (native 30m, categorical)

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  excavation_effort_multiplier_120m.*   (float32 raster)
  tappa8_excavation_effort_meta.json
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    EXCAVATION_EFFORT_MULTIPLIER,
    excavation_effort_multiplier,
    LITHOLOGY_CLASSES,
)
from suitability.terrain_metrics import block_mode

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- excavation effort (construction-difficulty follow-up to S8f) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

OUT = "data/processed/geomorphology"

log("loading lithology_v6 (30m) + land_mask (120m shape reference)...")
land = np.load("data/processed/climate/land_mask.npy").astype(bool)
lithology30 = np.load(f"{OUT}/lithology_v6.npy")

ny, nx = land.shape
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny

log("downsampling lithology_v6 to 120m via block-mode majority vote (same as S8f)...")
n_classes = len(LITHOLOGY_CLASSES)
lithology120 = block_mode(lithology30, 4, n_classes)[:ny, :nx]

effort = excavation_effort_multiplier(lithology120)

class_areas_km2 = {
    name: float((lithology120 == code).sum() * cs_x * cs_y / 1e6)
    for code, name in LITHOLOGY_CLASSES.items()
}
log(f"  120m lithology class areas (km2): {class_areas_km2}")

land_areas = {name: a for name, a in class_areas_km2.items() if name != "ocean"}
total_land_km2 = sum(land_areas.values())
area_weighted_mean_effort = sum(
    EXCAVATION_EFFORT_MULTIPLIER.get(code, 1.0) * class_areas_km2[LITHOLOGY_CLASSES[code]]
    for code in LITHOLOGY_CLASSES
    if code != 0
) / total_land_km2
log(f"  area-weighted mean excavation effort across all land: {area_weighted_mean_effort:.3f}")
log(f"  (dominated by basin_fill+greywacke+schist, ~92% of land, all <=2.3 -- "
    f"granite/volcanic's harsher values cover a small footprint)")

log("exporting excavation effort raster (120m, float32) + meta...")
np.save(f"{OUT}/excavation_effort_multiplier_120m.npy", effort.astype(np.float32))
write_envi_raw(
    f"{OUT}/excavation_effort_multiplier_120m", effort.astype(np.float32),
    xmin=XMIN, ymin=YMIN, cellsize=(cs_x + cs_y) / 2,
    description=(
        "Tappa 8 excavation effort: per-cell relative index of how hard the ground is "
        "to dig/cut, derived from lithology_v6 (block-mode downsampled to 120m). "
        "1.0 = basin_fill baseline (easiest), higher = more effort. Does NOT model "
        "foundation stability (karst voids), haulage, or labor -- see "
        "EXCAVATION_EFFORT_MULTIPLIER in geomorphology/lithology.py for values, "
        "citations, and caveats. UNREVIEWED first-pass estimates, not locked."
    ),
    dtype="f4",
)
write_prj(f"{OUT}/excavation_effort_multiplier_120m.prj", CRS_PROJ4)

meta = {
    "scope_note": (
        "A relative, per-lithology-class excavation-EFFORT index (how hard the raw rock "
        "is to break/cut), NOT a construction-cost model -- no foundation bearing "
        "capacity, haulage distance, or labor/skill availability data exists anywhere in "
        "this project to build that from. A distinct, complementary axis to S8f's travel "
        "friction: the two orderings deliberately diverge (see below) since they measure "
        "different physical properties."
    ),
    "excavation_effort_table": {
        LITHOLOGY_CLASSES[code]: mult for code, mult in EXCAVATION_EFFORT_MULTIPLIER.items()
    },
    "effort_status": "UNREVIEWED first-pass estimates -- grounded in relative mineral "
    "hardness (calcite ~Mohs 3 for limestone/marble vs. quartz/feldspar ~Mohs 6-7 for "
    "greywacke/schist/granite/basalt) and well-documented historical stoneworking "
    "accounts, not a cited numeric engineering source. Pending Nico's sign-off; not "
    "written to config/parameters.yml.",
    "important_caveat": "Does NOT capture karst foundation-stability risk (sinkholes/"
    "voids under marble and sedimentary_limestone) -- that's a separate hazard from raw "
    "workability. Taken alone this index makes marble/limestone look like easy building "
    "ground, which is only true for the quarrying half of the question.",
    "divergence_from_travel_friction": (
        "Granite: mild travel penalty (0.90, good footing on slab) but the HARDEST "
        "excavation value (3.0, no foliation weakness) -- easy to walk, brutal to quarry. "
        "Marble/limestone: the WORST travel penalty (0.60, karst hazard) but among the "
        "EASIEST excavation values (1.3-1.6, soft calcite) -- hard to walk, comparatively "
        "easy to cut. These are opposite orderings on purpose: traversability and raw-rock "
        "workability are different physical properties, not two views of the same fact."
    ),
    "resolution_m": [cs_x, cs_y],
    "class_areas_km2_120m": class_areas_km2,
    "area_weighted_mean_effort_all_land": area_weighted_mean_effort,
}
with open(f"{OUT}/tappa8_excavation_effort_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
