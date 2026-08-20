"""
Tappa 8 -- pod-level resource rasters for the two v6 rock classes that so far only had a
domain/tier ASSIGNMENT (locked in `scenario_reference.md` S21, relayed verbatim by Nico from
the Scenario chat) but no spatial raster: sedimentary_limestone (calcite, Onda secondary)
and granite (mica, Energia secondary; quartz, Onda tertiary). `marble` is DELIBERATELY
excluded -- the Scenario chat's own S21 assignment is an explicit null (no Vertice domain at
all, not a gap), recorded as a plain-text entry in resources.py, not a raster.

Three materials, same `place_material_pods()` machinery already used for the other ten:

- calcite (sedimentary_limestone, Onda secondary): eligible = lithology_v6 ==
  CLASS_SEDIMENTARY_LIMESTONE (87.72 km2 across "North Coast Limestone" + "Sedimentary Bay",
  the two authored zones compositing to this one class -- no need to treat them separately,
  same as schist/greywacke/basin_fill never distinguish which authored sub-zone a cell came
  from). UNIFORM weight: the handoff's "clear, intact calcite crystals... are much rarer to
  find in limestone cavities" is a statement about real-world SCARCITY, not a citable
  within-class spatial gradient (no analog to schist_grade exists for limestone) -- inventing
  one to have "something to weight by" would be worse than admitting none exists, same
  reasoning laumontite's placement already used.

- mica (granite, Energia secondary) + quartz (granite, Onda tertiary): eligible =
  lithology_v6 == CLASS_GRANITE (13.84 km2, "Granite South" alone). Both UNIFORM for the same
  reason as calcite -- granite is a compositionally uniform felsic rock in this model, no
  metamorphic-grade field like schist_grade exists for it. Handoff's own framing is explicit
  that this is "literally the same minerals as the schist assignment, just diversified
  sourcing" -- no new physical property claim, so no new weighting logic either.

DELIBERATELY NOT rescaling n_pods/radius for calcite's "rarity" framing: schist's own mica
("rarer than the quartz itself, same vein system") got the IDENTICAL n_pods/radius as quartz
and gold in that class (S13) -- the primary/secondary/tertiary TIER already encodes relative
scarcity in-fiction; pod count/footprint is a separate axis this project has not conflated
with tier before. Same standard placeholder parameters as every other generalized material
(n_pods=8, min_separation_km=5.0, radius_range_m=(300,800)) are used here for consistency,
NOT reduced -- if Nico/Scenario chat wants calcite's rarity reflected in footprint too, that
is a follow-up decision, not invented here.

Granite's actual result is worth flagging BEFORE it's run, not just reported after: 13.84 km2
is roughly 43x smaller than volcanic's 595.7 km2 (which itself only fit 6/8 silver_copper
pods under this same 5 km separation rule). A zone this small, compact, single-blob shape may
only geometrically admit 1-2 pods regardless of target n_pods=8 -- checked directly below,
not assumed either way.

Reads:
  data/processed/geomorphology/lithology_v6.npy

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  resource_calcite.npy, resource_mica_granite.npy, resource_quartz_granite.npy   (all new)
  tappa8_limestone_granite_materials_meta.json

Naming note: mica/quartz already have rasters from schist (S13, resource_mica.npy /
resource_quartz.npy) -- these are DIFFERENT, independent rasters for a DIFFERENT host rock,
suffixed `_granite` to avoid overwriting or confusing the two. The resource blend (S14) packs
both schist and granite sources as separate bits; a downstream consumer wanting "all mica
regardless of host rock" ORs the two rasters together, not something this script does for it.
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import (
    _grid_xy, place_material_pods, CLASS_SEDIMENTARY_LIMESTONE, CLASS_GRANITE,
)

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- calcite (sedimentary_limestone) + mica/quartz (granite) resource pods ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUT = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

log("loading lithology_v6...")
lithology_v6 = np.load(f"{OUT}/lithology_v6.npy")
ny, nx = lithology_v6.shape
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

# same standard placeholder parameters as every other generalized material (S8e/S10/S13) --
# NOT rescaled for calcite's rarity framing, see module docstring for why.
N_PODS = 8
MIN_SEPARATION_KM = 5.0
RADIUS_RANGE_M = (300.0, 800.0)
SEED_BASE = 130  # same family as S8e/S10/S13; offsets +20/+21/+22 are new, don't collide with
                  # +1..+4 (S8e), +3 reused (S12 placer magnetite re-placement), +10..+12 (S13)

MATERIALS = [
    ("calcite", CLASS_SEDIMENTARY_LIMESTONE, "sedimentary_limestone", "Onda -- secondary (optical-grade Iceland spar)", SEED_BASE + 20),
    ("mica_granite", CLASS_GRANITE, "granite", "Energia -- secondary (same piezoelectric mineral as schist's mica)", SEED_BASE + 21),
    ("quartz_granite", CLASS_GRANITE, "granite", "Onda -- tertiary (same birefringent mineral as schist's quartz)", SEED_BASE + 22),
]

results = {}
for name, class_code, host_name, domain_note, seed in MATERIALS:
    eligible = lithology_v6 == class_code
    eligible_km2 = float(eligible.sum() * cell_km2)
    log(f"--- {name} (host={host_name}, {domain_note}), UNIFORM weight (no citable "
        f"within-class gradient), eligible={eligible_km2:.2f} km2, seed={seed} ---")
    pod_mask, centers, radii = place_material_pods(
        eligible, None, xx, yy, RES_M,
        n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M,
        seed=seed,
    )
    results[name] = {
        "pod_mask": pod_mask, "centers_xy": centers, "radii_m": radii,
        "host": host_name, "domain_note": domain_note, "eligible_km2": eligible_km2,
    }
    n_placed = len(centers)
    pod_km2 = float(pod_mask.sum() * cell_km2)
    log(f"  {n_placed}/{N_PODS} pods placed, pod area={pod_km2:.2f} km2"
        + ("  <-- FEWER THAN TARGET, checked directly (see module docstring's granite-size flag)"
           if n_placed < N_PODS else ""))

log("checking mica_granite vs quartz_granite overlap directly -- SAME tiny eligible zone "
    "(13.84 km2), only the seed differs, so nonzero/high overlap would be unsurprising given "
    "how little room the 5 km separation rule leaves in a zone this small...")
overlap_cells = int((results["mica_granite"]["pod_mask"] & results["quartz_granite"]["pod_mask"]).sum())
overlap_km2 = overlap_cells * cell_km2
log(f"  mica_granite+quartz_granite overlap: {overlap_km2:.2f} km2 "
    f"({'zero' if overlap_cells == 0 else 'nonzero -- see meta for interpretation'})")

log("exporting rasters...")
NAME_MAP = {"calcite": "resource_calcite", "mica_granite": "resource_mica_granite",
            "quartz_granite": "resource_quartz_granite"}
for key, out_name in NAME_MAP.items():
    r = results[key]
    np.save(f"{OUT}/{out_name}.npy", r["pod_mask"])
    write_envi_raw(
        f"{OUT}/{out_name}", r["pod_mask"].astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 resource pods: {key} ({r['host']}, {r['domain_note']}, UNIFORM weight)",
        dtype="u1",
    )
    write_prj(f"{OUT}/{out_name}.prj", CRS_PROJ4)

meta = {
    "scope_note": "Pod rasters for the two S21-assigned rock classes (sedimentary_limestone "
    "-> calcite/Onda-secondary; granite -> mica/Energia-secondary + quartz/Onda-tertiary) "
    "per the Vertice material assignments Nico relayed verbatim from the Scenario chat's "
    "scenario_reference.md S21. marble is explicitly excluded -- S21 assigns it NO Vertice "
    "domain (recorded as a plain-text null entry in resources.py, not a raster here).",
    "weight_basis": "uniform for all three -- no citable within-class spatial gradient for "
    "either host rock (no metamorphic-grade analog to schist_grade exists for limestone or "
    "granite in this model). Not rescaled for calcite's real-world rarity framing -- that "
    "distinction is already carried by its 'secondary' tier, same treatment schist's own "
    "'rarer' mica got in S13 (identical n_pods/radius to quartz/gold there).",
    "n_pods_target": N_PODS,
    "min_separation_km": MIN_SEPARATION_KM,
    "radius_range_m": list(RADIUS_RANGE_M),
    "seed_offsets": {name: seed for name, _, _, _, seed in MATERIALS},
    "mica_granite_vs_quartz_granite_overlap_km2": overlap_km2,
    "materials": {
        name: {
            "host": r["host"],
            "domain_note": r["domain_note"],
            "eligible_km2": r["eligible_km2"],
            "n_pods_placed": len(r["centers_xy"]),
            "pod_area_km2": float(r["pod_mask"].sum() * cell_km2),
        }
        for name, r in results.items()
    },
    "marble_excluded_deliberately": "scenario_reference.md S21: marble carries no Vertice "
    "domain -- recrystallization destroys calcite's optical clarity (marble is cloudy/"
    "light-scattering, not transparent), and calcite's centrosymmetric crystal structure "
    "means it was never piezoelectric either, so it doesn't qualify for Energia by that "
    "route. Stays purely mundane (construction/monumental stone). No raster produced.",
}
with open(f"{OUT}/tappa8_limestone_granite_materials_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
