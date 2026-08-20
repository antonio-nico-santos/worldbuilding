"""
Tappa 8 -- separating quartz (Onda), mica (Energia), and gold (mundane) into their own
independent resource-pod rasters, per feedback relayed from the Scenario chat: these
three currently share ONE bundled "co-locates with jade" note in `resources.py`
(`spatial_note`), reusing `jade_pods_v5.npy` directly with no material-specific raster
of their own -- S8e's own docstring explicitly named this as one of only two deliberate
exclusions from the pod-placement generalization, on the reasoning that giving them
independent pods "would contradict" the co-location note. That reasoning holds for
"are these minerals found in the schist high-grade zone" (yes, unchanged here) but not
for "must their DISCOVERABLE footprints be pixel-identical to jade's" -- a real vein
system is heterogeneous at the pod scale (the same logic jade's own docstring already
uses against a smooth grade-distance function: "if [deposits] were [smooth], prospecting
would be trivial"), so three independent stochastic draws within the SAME plausible zone
is the more honest representation, not a contradiction of the co-location fact.

Eligible mask and weight field are UNCHANGED from jade's own method
(`jade_eligible_mask` + `schist_grade`, both real, already-locked machinery) -- this is
not inventing a new criterion, only generalizing `place_jade_pods`'s already-established
"real veins in this vein system aren't uniformly distributed" logic to gold/mica/quartz,
the same way S8e generalized it to laumontite/vivianite/etc.

ONE deliberate deviation from jade's own existing placement, flagged explicitly: this
recomputes `jade_eligible_mask` against `lithology_v6` (current, authoritative) rather
than `lithology_v5` (what `jade_pods_v5.npy` itself still uses, unchanged since S5/S8a).
Checked directly: 9,013 cells (8.11 km^2, ~0.46% of v5's schist area) that were schist in
v5 became marble in v6's priority-tier compositing (marble's priority_rank=1 outranks
schist) -- using v6 avoids seeding a new mineral pod on ground that is, per the CURRENT
lithology, no longer schist at all. `jade_pods_v5.npy` itself carries this same staleness
(unchanged, out of scope here) -- worth flagging as a real inconsistency between jade's
raster and the three new ones this script produces, not something this script silently
reconciles.

`schist_grade` was NOT in this session's data commit (S6's "continuous float32 fields...
dropped, cheap to regenerate" list) -- regenerated here by re-running
`run_tappa8_lithology_v5.py` in full (fully deterministic: connected-component landmass
ID, windowed max-min relief, percentile-threshold classification, no RNG until jade's own
seeded pod placement at the very end) and verified byte-for-byte identical to the already-
committed `lithology_v5.npy` and `jade_pods_v5.npy` before trusting the regenerated
`schist_grade_v5.npy` for anything.

Reads:
  data/processed/geomorphology/lithology_v6.npy
  data/processed/geomorphology/schist_grade_v5.npy   (regenerated this run, see above)

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  resource_quartz.npy, resource_mica.npy, resource_gold.npy   (all new)
  tappa8_schist_vein_materials_meta.json
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import _grid_xy, jade_eligible_mask, place_material_pods

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- quartz/mica/gold, separated into their own resource-pod rasters ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUT = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

log("loading lithology_v6 (current, authoritative) + schist_grade_v5 (regenerated, "
    "verified against the committed lithology_v5/jade_pods_v5 this same session)...")
lithology_v6 = np.load(f"{OUT}/lithology_v6.npy")
schist_grade = np.load(f"{OUT}/schist_grade_v5.npy")
ny, nx = lithology_v6.shape
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

log("recomputing jade_eligible_mask against v6 (SAME 80th-percentile schist_grade test "
    "jade itself uses) -- flagging the v5-vs-v6 divergence directly, not assuming none...")
eligible_v6 = jade_eligible_mask(lithology_v6, schist_grade, grade_percentile=80.0)
log(f"  eligible (high-grade schist) ground: {float(eligible_v6.sum()*cell_km2):.2f} km2")

# same placeholder pod parameters as the S8e/S10 material family (n_pods=8, not jade's own
# 10 -- these are independent materials generalized the same way laumontite/vivianite/etc
# were, not treated as jade's own special-cased raster). New seed offsets (SEED_BASE=130,
# +10/+11/+12) chosen to avoid colliding with S8e/S10/S12's existing +1..+6 offsets.
N_PODS = 8
MIN_SEPARATION_KM = 5.0
RADIUS_RANGE_M = (300.0, 800.0)
SEED_BASE = 130

MATERIALS = [
    ("quartz", "Onda -- birefringence", SEED_BASE + 10),
    ("mica", "Energia -- piezoelectric", SEED_BASE + 11),
    ("gold", "mundane -- no domain-gating property", SEED_BASE + 12),
]

results = {}
for name, domain_note, seed in MATERIALS:
    log(f"--- {name} ({domain_note}), weighted by schist_grade (same field jade uses), "
        f"seed={seed} ---")
    pod_mask, centers, radii = place_material_pods(
        eligible_v6, schist_grade, xx, yy, RES_M,
        n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M,
        seed=seed,
    )
    results[name] = {"pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii}
    log(f"  {len(centers)}/{N_PODS} pods placed, pod area={float(pod_mask.sum()*cell_km2):.2f} km2")

log("checking pairwise overlap among the three -- checked directly, not assumed. All "
    "three draw from the SAME eligible ground/weight field with only the seed differing, "
    "so nonzero overlap would be unremarkable if found (same status as any other "
    "independently-seeded pair in this project's resource layer)...")
names = list(results.keys())
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        n_cells = int((results[names[i]]["pod_mask"] & results[names[j]]["pod_mask"]).sum())
        log(f"  {names[i]}+{names[j]} overlap: {n_cells * cell_km2:.2f} km2")

log("exporting rasters...")
for name, r in results.items():
    np.save(f"{OUT}/resource_{name}.npy", r["pod_mask"])
    write_envi_raw(
        f"{OUT}/resource_{name}", r["pod_mask"].astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 resource pods: {name} (schist high-grade zone, weighted by schist_grade)",
        dtype="u1",
    )
    write_prj(f"{OUT}/resource_{name}.prj", CRS_PROJ4)

meta = {
    "scope_note": "Separates quartz/mica/gold from S8e's bundled 'co-locates with jade' "
    "note into their own independent pod rasters, per Scenario-chat feedback relayed by "
    "Nico. Eligible mask + weight field are jade's own established machinery "
    "(jade_eligible_mask + schist_grade), not a new criterion.",
    "v5_vs_v6_divergence_flag": "jade_eligible_mask recomputed here against lithology_v6 "
    "(current), not v5 (what jade_pods_v5.npy itself still uses) -- 9013 cells (8.11 km2) "
    "that were schist in v5 are marble in v6. jade_pods_v5.npy carries this staleness "
    "unchanged; these three new rasters do not.",
    "n_pods_target": N_PODS,
    "seed_offsets": {name: seed for name, _, seed in MATERIALS},
    "materials": {
        name: {
            "domain": domain_note,
            "n_pods_placed": len(r["centers_xy"]),
            "pod_area_km2": float(r["pod_mask"].sum() * cell_km2),
        }
        for (name, domain_note, _), r in zip(MATERIALS, results.values())
    },
    "eligible_km2": float(eligible_v6.sum() * cell_km2),
}
with open(f"{OUT}/tappa8_schist_vein_materials_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
