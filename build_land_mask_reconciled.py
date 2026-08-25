"""
Reconciles data/processed/climate/land_mask.npy (1334x1084, ~120m, Tappa-1-era)
against a true-ocean mask derived from the newer, more carefully built
lithology_v6.npy (5334x4334, ~30m native), at Nico's explicit request
(2026-08-23, following up on the Coastal_Village_03/11 coordinate-nudge
investigation in docs/decisions/10_tappa10_auxiliary_settlements.md).

WHY A SEPARATE FILE, NOT AN IN-PLACE OVERWRITE:
land_mask.npy is a shared input across the whole pipeline -- biomes,
lithology-120m, site suitability (Tappa 6/8), and now transport (Tappa 9/10)
were all computed against it. Overwriting it in place risks silently
invalidating every one of those already-locked outputs. This script instead
writes a NEW, transport-scoped mask to data/processed/transport/, which only
run_tappa9_road_network.py / run_tappa10_network_connections.py /
run_tappa9_sensitivity_analysis.py are switched to load.

METHOD:
1. true_ocean (native ~30m) = border-connected 8-connectivity components of
   (lithology_v6 == 0), via scipy.ndimage.label -- identical method already
   used and documented for the Tappa 10 third-pass coastal-village fix.
   Only components touching the raster's edge count as real/connected ocean;
   enclosed depressions (lakes etc.) are NOT true ocean.
2. true_land_native = ~true_ocean.
3. Aggregate to the 120m routing grid with block_any (NOT a >=0.5 majority
   mean) -- see FINDING below for why.
4. reconciled = land_mask OR true_land_120_any -- MONOTONIC / ADDITIVE ONLY.
   Never removes a cell land_mask already calls land, only adds cells the
   true-ocean analysis confirms are real land that land_mask currently
   misses. This protects the already-verified backbone property ("0 edges
   touching ocean or lake", checked cell-by-cell in Tappa 9's third pass):
   that check can only get MORE true after an additive-only change, never
   less -- any edge that was previously fully on land_mask=True cells stays
   exactly as valid.

FINDING THAT DROVE THE any-vs-majority CHOICE (2026-08-23):
The initial design (matching the general-reconciliation diagnostic) used a
>=0.5 majority-mean threshold, which recovers 10,468 false-negative cells
(0.724% of the grid) with 99.242% overall agreement with land_mask. But
checking it DIRECTLY against the two cells Nico actually asked about showed
it does NOT fix either one:
    Coastal_Village_03: true_land_120_mean=0.3125 (12/16 native subcells sea)
    Coastal_Village_11: true_land_120_mean=0.1875 (13/16 native subcells sea)
Both villages sit, at native 30m resolution, on real land (true_land_native
is True at their exact coordinate in both cases) -- but they're placed right
at the shoreline, which is exactly where a 120m cell is most likely to be
majority-sea even though it genuinely touches land. A majority vote
systematically fails shoreline points by design, not as an edge-case bug.
block_any instead asks "does this 120m cell contain ANY real land pixel" --
true for both villages -- which is the geometrically honest question for
"can a settlement plausibly sit on shore here." The cost: block_any recovers
26,062 cells (1.802% of the grid) instead of 10,468 (0.724%) -- 2.5x the
footprint, still 0% false-positive risk to already-verified land_mask=True
cells (monotonic OR never removes), but a larger area where a follow-up
Tappa 9 backbone re-run COULD find a cheaper route than before (never a more
expensive or newly-impossible one, since this is purely additive). That
re-run + a direct topology diff against the currently-locked backbone is a
separate, explicit verification step -- see run_tappa9_road_network.py's
re-run report.
"""
import json
import numpy as np
from scipy import ndimage

from params import load_params
from suitability.terrain_metrics import block_any

OUT_MASK = "data/processed/transport/land_mask_reconciled_v1.npy"
OUT_REPORT = "data/processed/transport/land_mask_reconciled_v1_report.json"

params = load_params("config/parameters.yml")
domain = params["domain"]
XMIN, XMAX, YMIN, YMAX = domain["xmin"], domain["xmax"], domain["ymin"], domain["ymax"]

land = np.load("data/processed/climate/land_mask.npy").astype(bool)
ny, nx = land.shape

lithology_v6 = np.load("data/processed/geomorphology/lithology_v6.npy")
ocean_mask_raw = lithology_v6 == 0
labeled, num = ndimage.label(ocean_mask_raw, structure=np.ones((3, 3)))
border_labels = set(labeled[0, :]) | set(labeled[-1, :]) | set(labeled[:, 0]) | set(labeled[:, -1])
border_labels.discard(0)
true_ocean_native = np.isin(labeled, list(border_labels))
true_land_native = ~true_ocean_native

true_land_120_any = block_any(true_land_native, 4)[:ny, :nx]

reconciled = land | true_land_120_any

added = reconciled & ~land
n_added = int(added.sum())
total = int(land.size)

report = {
    "source_land_mask": "data/processed/climate/land_mask.npy",
    "source_true_ocean_basis": "data/processed/geomorphology/lithology_v6.npy",
    "aggregation_rule": "block_any (see docstring for why not majority-mean >=0.5)",
    "reconciliation_rule": "monotonic OR -- only adds cells, never removes",
    "grid_shape": [ny, nx],
    "total_cells": total,
    "cells_added": n_added,
    "cells_added_pct": round(100 * n_added / total, 4),
    "cells_removed": 0,
    "agreement_with_land_mask_pct": round(100 * (land == reconciled).sum() / total, 4)
    if False
    else None,  # not meaningful post-hoc (reconciled is a superset by construction); see cells_added_pct instead
}

np.save(OUT_MASK, reconciled)
with open(OUT_REPORT, "w") as f:
    json.dump(report, f, indent=2)

print(f"reconciled mask saved: {OUT_MASK}  shape={reconciled.shape}  cells_added={n_added} ({report['cells_added_pct']}%)")
print(f"report saved: {OUT_REPORT}")

# quick sanity check on the two villages that motivated this
with open("data/processed/suitability/auxiliary_settlements_tappa10_v2.geojson") as f:
    aux = json.load(f)
by_name = {ft["properties"]["name"]: ft for ft in aux["features"]}
cs_x = (XMAX - XMIN) / nx
cs_y = (YMAX - YMIN) / ny


def xy_to_rc(x, y):
    col = int(round((x - XMIN) / cs_x - 0.5))
    row = int(round((YMAX - y) / cs_y - 0.5))
    return max(0, min(ny - 1, row)), max(0, min(nx - 1, col))


for name in ["Coastal_Village_03", "Coastal_Village_11"]:
    x, y = by_name[name]["geometry"]["coordinates"]
    r, c = xy_to_rc(x, y)
    print(f"  {name}: row,col=({r},{c})  land_mask={land[r, c]}  reconciled={reconciled[r, c]}")
