"""
Tappa 8 -- blend all seven Vertice/mundane resource-pod eligibility masks into ONE
raster, mirroring `make_tappa8_cave_blend.py`'s exact bitmask convention (Nico's
request there: "one file with the pixel value carrying type info" -- same ask, applied
to the resource-pod layer instead of the cave layer).

The seven materials that HAVE a pod raster (four from S8e, jade from S5, two more from
S10 -- see docs/decisions/08_tappa8_geomorphology.md):
  jade              (schist -- ALSO carries schist's gold-bearing quartz veins/mica,
                      which resources.py's own spatial_note says co-locate with jade,
                      no separate raster)
  laumontite        (greywacke, primary)
  vivianite         (basin_fill, primary)
  bog_iron          (basin_fill, mundane_only -- co-located with vivianite, S10)
  placer_magnetite  (basin_fill, secondary_weak)
  silver_copper     (volcanic, secondary -- silver+copper together, one reservoir,
                      differentiated by kind in tappa8_resource_pods_meta.json's
                      pod_kinds list, not by separate spatial footprints)
  bauxite           (volcanic, mundane_only, S10)

Excluded, same reasoning as S8e/S10: volcanic's magnetite (primary) is a disseminated
bulk mineral through the whole basalt class, not vein/pod-localized -- inventing a pod
footprint for it would be LESS geologically honest than what resources.py already does
(a class-level fact). No raster, no bit.

UNLIKE the five cave types (independent phenomena that can genuinely co-occur on the
same ground), overlap here is structurally CONSTRAINED, not open -- worth checking
directly rather than assuming either "none" or "same as caves":
- jade/laumontite/{vivianite,bog_iron,placer_magnetite}/{silver_copper,bauxite} each
  draw their eligible_mask from a DIFFERENT lithology class (schist / greywacke /
  basin_fill / volcanic), and lithology classes are mutually exclusive by construction
  (one code per cell) -- so a jade pod and a laumontite pod CANNOT overlap, structurally,
  not just empirically. Checked below to confirm the data agrees with the geometry.
- WITHIN basin_fill, vivianite/bog_iron/placer_magnetite are three INDEPENDENT
  stochastic draws over the same eligible ground (bog_iron is the one exception --
  identical to vivianite BY DESIGN, S10 -- so its bit is fully redundant with
  vivianite's, kept separate only for per-material catalogue completeness, not because
  it carries new spatial information). vivianite and placer_magnetite use different
  weight fields (wetness vs. coastal-x-volcanic proximity) and different seeds, so any
  overlap between THEM is incidental, not structural -- checked below, not assumed.
- WITHIN volcanic, silver_copper and bauxite are two independent stochastic draws
  (vent-proximity vs. flatness weighting) -- same "incidental, not structural" status.

Run after run_tappa8_iron_aluminium.py (needs resource_bog_iron.npy/resource_bauxite.npy)
and with jade_pods_v5.npy staged (Tappa 8 S5, closed stage, committed separately from
this session's other resource work).
"""
import sys, time, json

sys.path.insert(0, "src")
import numpy as np

from params import load_params
from terrain.raster_io import write_envi_raw, write_prj

t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.1f}s] {msg}", flush=True)


log("=== Tappa 8 -- blend all 7 resource-pod eligibility masks into one bitmask raster ===")

params = load_params("config/parameters.yml")
domain = params["domain"]
CRS_PROJ4 = params["crs"]["PROJ string parameter"]
RES_M = domain["resolution_m"]
OUTPUT_DIR = "data/processed/geomorphology"
cell_km2 = (RES_M / 1000.0) ** 2

# bit -> (name, source file, lithology class it's drawn from, value)
BITS = [
    ("jade", "jade_pods_v5.npy", "schist", 1),
    ("laumontite", "resource_laumontite.npy", "greywacke", 2),
    ("vivianite", "resource_vivianite.npy", "basin_fill", 4),
    ("bog_iron", "resource_bog_iron.npy", "basin_fill", 8),
    ("placer_magnetite", "resource_placer_magnetite.npy", "basin_fill", 16),
    ("silver_copper", "resource_silver_copper.npy", "volcanic", 32),
    ("bauxite", "resource_bauxite.npy", "volcanic", 64),
]

log("loading all seven resource-pod masks...")
masks = {}
for name, fname, lith_class, bit in BITS:
    masks[name] = np.load(f"{OUTPUT_DIR}/{fname}").astype(bool)

ny, nx = next(iter(masks.values())).shape
blend = np.zeros((ny, nx), dtype=np.uint8)
for name, fname, lith_class, bit in BITS:
    blend |= (masks[name].astype(np.uint8) * bit)

log("computing overlap diagnostics -- checked directly, not assumed, since the "
    "structural-vs-incidental distinction above needs real numbers to back it up...")
names = [b[0] for b in BITS]
lith_of = {b[0]: b[2] for b in BITS}
pairwise_km2 = {}
cross_class_overlap_found = False
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        n_cells = int((masks[names[i]] & masks[names[j]]).sum())
        if n_cells:
            pairwise_km2[f"{names[i]}+{names[j]}"] = round(n_cells * cell_km2, 3)
            if lith_of[names[i]] != lith_of[names[j]]:
                cross_class_overlap_found = True

n_types_per_cell = np.stack(list(masks.values()), axis=0).sum(axis=0)
multi_type_km2 = float((n_types_per_cell >= 2).sum() * cell_km2)
max_simultaneous = int(n_types_per_cell.max())
area_km2 = {name: float(masks[name].sum() * cell_km2) for name, _, _, _ in BITS}
any_resource_km2 = float((blend > 0).sum() * cell_km2)

log(f"per-material area km2: {area_km2}")
log(f"pairwise overlap km2 (only nonzero pairs shown): {pairwise_km2}")
log(f"cells with 2+ simultaneous materials: {multi_type_km2:.2f} km2, "
    f"max simultaneous on one cell: {max_simultaneous}")
log(f"any-resource-material union: {any_resource_km2:.2f} km2")
log(f"cross-lithology-class overlap found: {cross_class_overlap_found} "
    "(should be False -- lithology classes are mutually exclusive by construction, "
    "so this would indicate a real bug if True)")
assert not cross_class_overlap_found, (
    "Cross-lithology-class overlap detected between resource pods -- this should be "
    "structurally impossible given lithology_v6's mutual-exclusivity; investigate "
    "before shipping, do not silently proceed."
)
bog_iron_vivianite_km2 = float((masks["vivianite"] & masks["bog_iron"]).sum() * cell_km2)
log(f"vivianite+bog_iron overlap: {bog_iron_vivianite_km2:.4f} km2 -- expected to equal "
    f"BOTH materials' own area exactly (co-located by design, S10), i.e. bog_iron's bit "
    f"carries zero new spatial information beyond vivianite's -- checking that now...")
assert np.array_equal(masks["vivianite"], masks["bog_iron"]), (
    "bog_iron was supposed to be IDENTICAL to vivianite's pod geometry (S10) -- the raw "
    "masks don't match cell-for-cell, something diverged since that script ran."
)
log("  confirmed: bog_iron is a byte-for-byte structural subset (in fact exact match) "
    "of vivianite -- flagged in the meta JSON below, not hidden.")

log("round-trip verification: decoding all 7 bits back out of the blend must reproduce "
    "every source mask exactly...")
for name, fname, lith_class, bit in BITS:
    decoded = (blend & bit) > 0
    assert np.array_equal(decoded, masks[name]), f"round-trip mismatch for {name}!"
log("  all 7 materials round-trip exactly.")

log("exporting blended bitmask raster...")
np.save(f"{OUTPUT_DIR}/resource_blend.npy", blend)
description = (
    "Tappa 8 resource-pod eligibility BITMASK (OR of all seven materials that have a "
    "pod raster): " + ", ".join(f"{bit}={name}" for name, _, _, bit in BITS)
    + ". Decode a single material with (value & bit); value==0 means no pod material "
      "present; value>0 means at least one. NOTE: bit 8 (bog_iron) is by construction "
      "IDENTICAL to bit 4 (vivianite) -- co-located deposit (S10), not independent "
      "information. Volcanic's magnetite (primary) has no bit -- disseminated bulk "
      "mineral, not pod-localized, stays a class-level fact in resources.py."
)
write_envi_raw(
    f"{OUTPUT_DIR}/resource_blend", blend.astype(np.int16),
    xmin=domain["xmin"], ymin=domain["ymin"], cellsize=RES_M,
    description=description, dtype="u1",
)
write_prj(f"{OUTPUT_DIR}/resource_blend.prj", CRS_PROJ4)

meta = {
    "encoding": "bitmask, not category code -- mirrors cave_blend.npy's exact convention "
    "(S8d). Unlike the five cave types, overlap here is structurally CONSTRAINED (each "
    "material's eligible_mask is drawn from one lithology class, and lithology classes "
    "are mutually exclusive by construction) rather than open -- verified directly, see "
    "cross_lithology_class_overlap_found below, not assumed.",
    "bits": {name: bit for name, _, _, bit in BITS},
    "lithology_class_per_material": {name: lith_class for name, _, lith_class, _ in BITS},
    "area_km2": area_km2,
    "overlap_km2": pairwise_km2,
    "cells_with_2plus_materials_km2": multi_type_km2,
    "max_simultaneous_materials_on_one_cell": max_simultaneous,
    "any_resource_material_union_km2": any_resource_km2,
    "cross_lithology_class_overlap_found": cross_class_overlap_found,
    "bog_iron_flag": "bit 8 (bog_iron) is an EXACT structural subset of bit 4 (vivianite) "
    "-- co-located by design (S10's 'single wetland patch yields both' framing), not an "
    "independent placement. Kept as its own bit for per-material catalogue completeness "
    "(consistent with resources.py listing it as its own mundane_only material), not "
    "because it adds spatial information beyond vivianite's own bit.",
    "excluded_no_raster": {
        "schist_gold_quartz_mica": "co-locates with jade (bit 1) per resources.py's own "
        "spatial_note -- no separate bit, would be double-counting the same footprint.",
        "volcanic_magnetite_primary": "disseminated bulk mineral in basalt, not vein/"
        "joint/pod-localized like the other seven -- stays a class-level fact in "
        "resources.py, same reasoning as S8e/S10.",
    },
    "how_to_decode": "single_material_mask = (raster_value & bit_value) > 0; "
    "any_resource_mask = raster_value > 0",
    "round_trip_verified": True,
}
with open(f"{OUTPUT_DIR}/tappa8_resource_blend_meta.json", "w") as f:
    json.dump(meta, f, indent=2)

log("=== DONE ===")
