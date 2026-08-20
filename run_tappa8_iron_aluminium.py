"""
Tappa 8 -- iron (bog iron) and aluminium (bauxite) resource additions, resolving the
two items decision doc S8e left explicitly open ("Iron -- currently confined to a
small footprint..." / "Aluminium -- genuinely unresolved...").

Both are Nico's explicit calls (not this script's invention):
- Iron: bog iron added as sedimentary_basin_fill's mundane_only material, CO-LOCATED
  with vivianite -- scenario_reference.md S22.4 frames the two as "co-located uses of
  one zone... a single wetland patch plausibly yields both, worked by different
  specialists." Taken literally: bog iron's pod footprint is IDENTICAL to vivianite's,
  not an independent stochastic placement. Doing this "right" also means vivianite's
  own pods move from S8e's original whole-basin_fill placement to S9's
  wetland_backswamp sub-class specifically (S9 didn't exist yet when S8e ran) --
  this REPLACES resource_vivianite.npy's already-committed output, done here on
  Nico's explicit go-ahead, not silently.
- Aluminium: Nico's call was "add the resource; electrolysis is in development with
  Vertice assistance" -- i.e. treat NZ's real (subeconomic-in-reality) Northland
  bauxite as an ACTIVE resource here, with the historical Hall-Heroult technology
  gate narratively closed by Vertice-assisted electrolysis rather than left as an
  unresolved "no ore / can't extract it" ambiguity. The Vertice-assistance MECHANIC
  itself is a Scenario-chat-level claim -- this script only records it as citation
  context for why the resource is viable, no domain-gating property is assigned to
  the ore itself (see resources.py's mundane_only entry).

Bauxite is placed on the volcanic class (matching the real Northland citation:
Otoroa/Matauri Bay bauxite is a lateritic weathering product of the Tangihua Complex's
basaltic rock, not of granite or any other class this world has), weighted by
FLATNESS (1/(1+slope_pct)) -- real lateritic bauxite caps form and survive on
low-relief surfaces and get stripped off steep slopes by erosion, the same
slope-based logic this codebase already uses for lava tubes (inverted) and
talus/sea caves (direct). Placed as a deliberately SMALL footprint (2 pods, tighter
radius range than the other five materials) to match the real citation's own framing
("small, ~20 Mt largest, historically subeconomic") -- narratively now workable, but
not suddenly abundant.

Reads:
  data/processed/dem_v3_final_30m_eroded.npy
  data/processed/geomorphology/lithology_v6.npy
  data/processed/geomorphology/basin_zonation_30m.npy
  data/processed/hydrology/stream_mask.npy

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  resource_vivianite.npy       (OVERWRITES S8e's version -- now wetland_backswamp-only)
  resource_bog_iron.npy        (new -- identical pod geometry to resource_vivianite.npy)
  resource_bauxite.npy         (new)
  tappa8_iron_aluminium_meta.json
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_BASIN_FILL, CLASS_VOLCANIC, _grid_xy, place_material_pods
from geomorphology.basin_zonation import CLASS_WETLAND_BACKSWAMP
from suitability.terrain_metrics import compute_slope_pct

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- iron (bog iron) + aluminium (bauxite), resolving S8e's two open items ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUT = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

log("loading DEM + lithology_v6 + basin_zonation + stream_mask...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
ny, nx = dem.shape
lithology_v6 = np.load(f"{OUT}/lithology_v6.npy")
basin_zonation = np.load(f"{OUT}/basin_zonation_30m.npy")
stream_mask = np.load("data/processed/hydrology/stream_mask.npy").astype(bool)
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

# same placeholder pod parameters as run_tappa8_resource_pods.py, same seed base/spirit
N_PODS = 8
MIN_SEPARATION_KM = 5.0
RADIUS_RANGE_M = (300.0, 800.0)
SEED_BASE = 130

# ---------------------------------------------------------------------------------
# Iron -- re-place vivianite within wetland_backswamp (S9 didn't exist during S8e),
# then co-locate bog iron on the IDENTICAL pod geometry.
# ---------------------------------------------------------------------------------
log("--- vivianite RE-PLACEMENT: basin_fill (S8e) -> wetland_backswamp (S9), same seed/weight ---")
wetland_backswamp_mask = basin_zonation == CLASS_WETLAND_BACKSWAMP
old_basin_fill_mask = lithology_v6 == CLASS_BASIN_FILL
dist_to_stream_km = ndimage.distance_transform_edt(~stream_mask, sampling=(RES_M, RES_M)) / 1000.0
wetness = 1.0 / (1.0 + dist_to_stream_km)

log(f"  eligible ground shrinks from {float(old_basin_fill_mask.sum()*cell_km2):.1f} km2 (whole "
    f"basin_fill, S8e) to {float(wetland_backswamp_mask.sum()*cell_km2):.1f} km2 "
    f"(wetland_backswamp only, S9)")

vivianite_pod_mask, vivianite_centers, vivianite_radii = place_material_pods(
    wetland_backswamp_mask, wetness, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M,
    seed=SEED_BASE + 2,  # SAME seed as S8e's original vivianite placement
)
log(f"  {len(vivianite_centers)}/{N_PODS} pods placed (was 8/8 under the old, larger mask), "
    f"pod area={float(vivianite_pod_mask.sum()*cell_km2):.2f} km2")

log("--- bog iron: CO-LOCATED with vivianite -- identical pod geometry, not an independent "
    "placement (scenario_reference.md S22.4: 'a single wetland patch plausibly yields both') ---")
bog_iron_pod_mask = vivianite_pod_mask.copy()
log(f"  bog iron pod area={float(bog_iron_pod_mask.sum()*cell_km2):.2f} km2 (== vivianite's, by design)")

# ---------------------------------------------------------------------------------
# Aluminium -- bauxite on the volcanic class (real Northland citation is lateritic
# weathering of basaltic rock), weighted by flatness (real laterite caps form/survive
# on low-relief ground, stripped off steep slopes by erosion).
# ---------------------------------------------------------------------------------
log("--- bauxite (volcanic, weighted by flatness = 1/(1+slope_pct) -- real lateritic caps "
    "form/survive on low-relief ground, same slope-based logic caves.py already uses) ---")
volcanic_mask = lithology_v6 == CLASS_VOLCANIC
slope_pct = compute_slope_pct(dem, RES_M)
flatness = 1.0 / (1.0 + slope_pct)

N_PODS_BAUXITE = 2  # deliberately smaller than the other five materials' 8 -- matches the real
# citation's own framing (Northland bauxite: small, ~20 Mt largest, historically subeconomic).
# Narratively workable now (Vertice-assisted electrolysis), not suddenly abundant.
RADIUS_RANGE_M_BAUXITE = (200.0, 500.0)  # tighter than the 300-800m default, same "small" reasoning

bauxite_pod_mask, bauxite_centers, bauxite_radii = place_material_pods(
    volcanic_mask, flatness, xx, yy, RES_M,
    n_pods=N_PODS_BAUXITE, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M_BAUXITE,
    seed=SEED_BASE + 6,
)
log(f"  {len(bauxite_centers)}/{N_PODS_BAUXITE} pods placed, "
    f"pod area={float(bauxite_pod_mask.sum()*cell_km2):.2f} km2 within "
    f"{float(volcanic_mask.sum()*cell_km2):.1f} km2 eligible volcanic")

# --- export -------------------------------------------------------------------------
log("exporting rasters (resource_vivianite.npy OVERWRITTEN, two new: bog_iron, bauxite)...")
EXPORTS = {
    "resource_vivianite": vivianite_pod_mask,
    "resource_bog_iron": bog_iron_pod_mask,
    "resource_bauxite": bauxite_pod_mask,
}
for name, mask in EXPORTS.items():
    np.save(f"{OUT}/{name}.npy", mask)
    write_envi_raw(
        f"{OUT}/{name}", mask.astype(np.int16),
        xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
        description=f"Tappa 8 resource pods: {name}",
        dtype="u1",
    )
    write_prj(f"{OUT}/{name}.prj", CRS_PROJ4)

meta = {
    "scope_note": (
        "Resolves decision doc S8e's two explicitly-open items (iron, aluminium), both on "
        "Nico's explicit direction, not unilateral implementation choices."
    ),
    "iron": {
        "change": "vivianite pods RE-PLACED from S8e's whole-basin_fill mask to S9's "
        "wetland_backswamp sub-class (S9 postdates S8e) -- same seed/weight, only the "
        "eligible mask changed. OVERWRITES resource_vivianite.npy.",
        "vivianite_pods_placed": len(vivianite_centers),
        "vivianite_pod_area_km2": float(vivianite_pod_mask.sum() * cell_km2),
        "bog_iron": "co-located with vivianite -- IDENTICAL pod geometry, not independently "
        "placed, per scenario_reference.md S22.4's 'single wetland patch yields both' framing.",
        "bog_iron_domain": "mundane_only -- no citable Vertice-domain-gating physical property "
        "(no analog to titanomagnetite's real magnetism), same category as schist's gold/jade.",
        "citation_honesty_flag": "bog iron (goethite/limonite, anoxic floodplain precipitation) "
        "is a well-documented GENERAL pre-industrial process (Iron Age Scandinavia through "
        "colonial North America), but NOT NZ-specific -- Te Ara's NZ iron history is exclusively "
        "the ironsand/titanomagnetite story. First exception in this project's resource layer "
        "to the NZ-specific-citation norm.",
    },
    "aluminium": {
        "decision_source": "Nico, explicit: 'Adicionar recurso. Eletrolise esta em "
        "desenvolvimento com auxilio de Vertices.'",
        "material": "bauxite (gibbsite/boehmite lateritic weathering crust)",
        "class": "volcanic -- matches the real citation (Northland/Otoroa/Matauri Bay bauxite "
        "is a lateritic weathering product of the Tangihua Complex's BASALTIC rock, not granite "
        "or any other class this world has).",
        "spatial_weighting": "flatness = 1/(1+slope_pct) within the volcanic class -- real "
        "lateritic caps form and are preserved on low-relief surfaces and get stripped off "
        "steep slopes by erosion; reuses the existing slope field, same logic direction "
        "caves.py already uses for lava tubes (gentle-slope preference).",
        "n_pods": N_PODS_BAUXITE,
        "pod_area_km2": float(bauxite_pod_mask.sum() * cell_km2),
        "deliberately_small_footprint": "2 pods / tighter 200-500m radius range vs the other "
        "five materials' 8 pods / 300-800m -- matches the real citation's own framing (small, "
        "~20 Mt largest, historically subeconomic). Vertice-assisted electrolysis (Nico's "
        "framing) closes the historical Hall-Heroult TECHNOLOGY gate; it does not make the "
        "ore itself abundant -- these stay two separate constraints.",
        "domain": "mundane_only -- the ore itself has no citable domain-gating physical "
        "property; the Vertice involvement is at the EXTRACTION/PROCESS level (assisted "
        "electrolysis), not a property of the mineral. That process-level mechanic is a "
        "Scenario-chat-level claim -- recorded here as citation context only, not elaborated.",
        "not_resolved_here": "whether the historically-real 'subeconomic, never mined' framing "
        "still applies narratively (i.e. how RARE/hard-won this bauxite should feel in play) is "
        "a tone question, not a geology question -- left for Scenario chat if it matters "
        "mechanically.",
    },
    "unchanged_flagged_for_consistency": {
        "placer_magnetite": "still placed on the whole basin_fill mask (S8e), weighted by "
        "coastal x volcanic proximity -- NOT restricted to S9's estuarine_coastal sub-class, "
        "even though that sub-class now exists and its own definition is coastal-proximity-"
        "driven. Consistent with this script's actual scope (iron, aluminium only) -- flagged "
        "as a natural follow-up, not done here to avoid unrequested scope creep.",
    },
}
with open(f"{OUT}/tappa8_iron_aluminium_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
