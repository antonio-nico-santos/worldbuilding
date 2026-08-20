"""
Tappa 8 -- placer magnetite re-placement: basin_fill (S8e) -> estuarine_coastal (S9),
closing the last of S10/S11's flagged consistency items ("placer magnetite still uses
the whole basin_fill mask rather than S9's estuarine_coastal sub-class, even though that
sub-class's own definition is coastal-proximity-driven and would be a natural fit" --
flagged in S10, flagged again in S11, not done until now on Nico's explicit go-ahead).

Same treatment S10 already gave vivianite: SAME seed, SAME weight field, only the
eligible mask changes -- from S8e's whole basin_fill (4334.5 km2) to S9's
estuarine_coastal sub-class (640.8 km2). Placer magnetite's own weight field (coastal x
volcanic-landmass proximity) was ALREADY a coastal-proximity signal before this -- S9's
estuarine_coastal sub-class is *itself* defined by the same kind of coastal-proximity
threshold (S9's coastal_threshold_km test), so restricting the eligible ground to that
sub-class is tightening an already-coastal-biased placement to the zone that's
authoritatively coastal, not introducing a new criterion.

Expected structural consequence, checked directly below rather than assumed: S9's three
basin_fill sub-zones (arable / wetland_backswamp / estuarine_coastal) are mutually
exclusive by construction (classify_basin_fill_zones assigns exactly one code per basin_
fill cell). Vivianite/bog_iron are restricted to wetland_backswamp (S10); placer
magnetite is now restricted to estuarine_coastal here -- these two eligible masks CANNOT
overlap, so placer_magnetite's pods and vivianite/bog_iron's pods can no longer overlap
either, even incidentally. S11's resource blend found zero incidental overlap between
these pairs anyway (with the old whole-basin_fill mask) -- this closes that possibility
structurally, not just empirically.

Reads:
  data/processed/dem_v3_final_30m_eroded.npy
  data/processed/geomorphology/lithology_v6.npy
  data/processed/geomorphology/basin_zonation_30m.npy

Writes to data/processed/geomorphology/ (gitignored, regenerate locally):
  resource_placer_magnetite.npy   (OVERWRITES S8e's version -- now estuarine_coastal-only)
  tappa8_placer_magnetite_restrict_meta.json
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np
from scipy import ndimage

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj
from geomorphology.lithology import CLASS_BASIN_FILL, CLASS_VOLCANIC, _grid_xy, place_material_pods
from geomorphology.basin_zonation import CLASS_ESTUARINE_COASTAL
from geomorphology.caves import distance_to_ocean_km

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- placer magnetite RE-PLACEMENT: basin_fill (S8e) -> estuarine_coastal (S9) ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUT = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

log("loading DEM + lithology_v6 + basin_zonation...")
dem = np.load("data/processed/dem_v3_final_30m_eroded.npy")
land_mask = dem > 0
ny, nx = dem.shape
lithology_v6 = np.load(f"{OUT}/lithology_v6.npy")
basin_zonation = np.load(f"{OUT}/basin_zonation_30m.npy")
xx, yy = _grid_xy(domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"], ny, nx)

# same placeholder pod parameters and seed as run_tappa8_resource_pods.py's original
# placer magnetite placement (SEED_BASE=130, +3) -- only the eligible mask changes here.
N_PODS = 8
MIN_SEPARATION_KM = 5.0
RADIUS_RANGE_M = (300.0, 800.0)
SEED_BASE = 130

log("recomputing the SAME weight field as S8e (coastal x volcanic-landmass proximity) -- "
    "unchanged, only the eligible mask is new...")
old_basin_fill_mask = lithology_v6 == CLASS_BASIN_FILL
estuarine_coastal_mask = basin_zonation == CLASS_ESTUARINE_COASTAL
dist_to_ocean_km = distance_to_ocean_km(land_mask, RES_M)
volcanic_mask = lithology_v6 == CLASS_VOLCANIC
dist_to_volcanic_km = ndimage.distance_transform_edt(~volcanic_mask, sampling=(RES_M, RES_M)) / 1000.0
coastal_island_weight = (1.0 / (1.0 + dist_to_ocean_km)) * (1.0 / (1.0 + dist_to_volcanic_km))

log(f"  eligible ground shrinks from {float(old_basin_fill_mask.sum()*cell_km2):.1f} km2 (whole "
    f"basin_fill, S8e) to {float(estuarine_coastal_mask.sum()*cell_km2):.1f} km2 "
    f"(estuarine_coastal only, S9)")

pod_mask, centers, radii = place_material_pods(
    estuarine_coastal_mask, coastal_island_weight, xx, yy, RES_M,
    n_pods=N_PODS, min_separation_km=MIN_SEPARATION_KM, radius_range_m=RADIUS_RANGE_M,
    seed=SEED_BASE + 3,  # SAME seed as S8e's original placer magnetite placement
)
log(f"  {len(centers)}/{N_PODS} pods placed, pod area={float(pod_mask.sum()*cell_km2):.2f} km2")

log("checking the structural non-overlap claim directly against vivianite/bog_iron "
    "(both restricted to wetland_backswamp, S10) -- should be zero, not just small...")
vivianite_mask = np.load(f"{OUT}/resource_vivianite.npy").astype(bool)
overlap_cells = int((pod_mask & vivianite_mask).sum())
log(f"  placer_magnetite vs vivianite/bog_iron pod overlap: {overlap_cells} cells "
    f"({'confirmed zero, as expected from S9 sub-zone exclusivity' if overlap_cells == 0 else 'UNEXPECTED -- investigate'})")
assert overlap_cells == 0, (
    "Expected zero overlap between placer_magnetite (estuarine_coastal) and vivianite/"
    "bog_iron (wetland_backswamp) pods -- S9's sub-zones are supposed to be mutually "
    "exclusive by construction. A nonzero count here means that assumption is wrong."
)

log("exporting resource_placer_magnetite.npy (OVERWRITTEN)...")
np.save(f"{OUT}/resource_placer_magnetite.npy", pod_mask)
write_envi_raw(
    f"{OUT}/resource_placer_magnetite", pod_mask.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=(
        "Tappa 8 resource pods: placer_magnetite, RESTRICTED to S9's estuarine_coastal "
        "sub-class (re-placed from S8e's whole-basin_fill mask), weighted by coastal x "
        "volcanic-landmass proximity."
    ),
    dtype="u1",
)
write_prj(f"{OUT}/resource_placer_magnetite.prj", CRS_PROJ4)

meta = {
    "scope_note": "Closes the placer-magnetite/estuarine_coastal consistency item flagged "
    "in decision doc S10 and S11 -- Nico's explicit go-ahead, same treatment as S10's "
    "vivianite re-placement (same seed, same weight field, only the eligible mask changed).",
    "change": "eligible mask restricted from whole basin_fill (S8e, 4334.5 km2) to S9's "
    "estuarine_coastal sub-class (640.8 km2). Weight field (coastal x volcanic-landmass "
    "proximity) and seed (133) unchanged.",
    "n_pods_placed": len(centers),
    "pod_area_km2": float(pod_mask.sum() * cell_km2),
    "structural_non_overlap_check": {
        "against": "resource_vivianite.npy (S10, restricted to wetland_backswamp)",
        "overlap_cells_found": overlap_cells,
        "expected": 0,
        "reasoning": "S9's basin_fill sub-zones (arable/wetland_backswamp/estuarine_coastal) "
        "are mutually exclusive by construction -- confirmed here by assertion, not assumed.",
    },
}
with open(f"{OUT}/tappa8_placer_magnetite_restrict_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
