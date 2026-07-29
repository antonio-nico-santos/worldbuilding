"""
Tappa 6 (step 3) -- biome suitability lookup, at the 120 m working grid.
See src/suitability/biome_lookup.py's module docstring for the full method,
the economic-base and ice/rock decisions, and honest limitations.

Reads:
  data/processed/biomes/biome_id.npy                (Tappa 5, 120 m)
  data/processed/climate/land_mask.npy               (Tappa 2, 120 m)
  data/processed/suitability/slope_suitability_120m.npy  (step 1, for the
                                                          ice/rock cross-check
                                                          logged in meta only)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  biome_id_smoothed_120m       (majority-filtered biome_id, int8 -- kept
                                 alongside the raw biome_id so the smoothing
                                 step is itself inspectable, same "keep
                                 component layers visible" precedent as
                                 Tappa 5's biotemperature_c/pet_ratio)
  biome_suitability_120m       (0-1, float32)
  tappa6_biome_meta.json
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

from src.biomes.world_biomes import BIOME_NAMES
from src.suitability.biome_lookup import (
    BIOME_NOTES,
    BIOME_SUITABILITY,
    biome_suitability_from_id,
    majority_filter_biome_id,
)
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)
MAJORITY_WINDOW = 3


def _area_km2_by_class(biome_id, cs_x, cs_y, n_classes):
    cell_km2 = (cs_x / 1000.0) * (cs_y / 1000.0)
    counts = np.bincount(biome_id.ravel().astype(np.int64), minlength=n_classes)
    return {BIOME_NAMES[k]: float(counts[k] * cell_km2) for k in range(1, n_classes)}


def main():
    t0 = time.time()
    biome_id = np.load("data/processed/biomes/biome_id.npy")
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)
    slope_suit = np.load("data/processed/suitability/slope_suitability_120m.npy")

    ny, nx = land.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny
    n_classes = len(BIOME_NAMES)

    biome_smoothed, frac_changed = majority_filter_biome_id(biome_id, land, MAJORITY_WINDOW)
    suit = biome_suitability_from_id(biome_smoothed, land).astype(np.float32)

    out = "data/processed/suitability"
    os.makedirs(out, exist_ok=True)

    biome_smoothed_i8 = biome_smoothed.astype(np.int8)
    np.save(f"{out}/biome_id_smoothed_120m.npy", biome_smoothed_i8)
    write_envi_raw(
        f"{out}/biome_id_smoothed_120m", biome_smoothed_i8, XMIN, YMIN, cs_x,
        "Tappa6 biome_id_smoothed_120m", dtype="i2",
    )
    write_prj(f"{out}/biome_id_smoothed_120m.prj", PROJ4)

    np.save(f"{out}/biome_suitability_120m.npy", suit)
    write_envi_raw(
        f"{out}/biome_suitability_120m", suit, XMIN, YMIN, cs_x,
        "Tappa6 biome_suitability_120m", dtype="f4",
    )
    write_prj(f"{out}/biome_suitability_120m.prj", PROJ4)

    # ice/rock cross-check, logged for the record (this is what motivated
    # keeping ice/rock as nucleo-only vs. also-exclusao -- see BIOME_NOTES)
    ice_rock_check = {}
    for name, bid in [("Permanent Snow & Ice", 1), ("Alpine Fellfield", 2)]:
        mask = biome_id == bid
        ice_rock_check[name] = {
            "n_cells": int(mask.sum()),
            "frac_slope_suitability_gt_0.3": float((slope_suit[mask] > 0.3).mean()),
            "slope_suitability_p50_p90": [
                float(v) for v in np.percentile(slope_suit[mask], [50, 90])
            ],
        }

    area_before = _area_km2_by_class(biome_id, cs_x, cs_y, n_classes)
    area_after = _area_km2_by_class(biome_smoothed, cs_x, cs_y, n_classes)

    meta = {
        "resolution_m": [cs_x, cs_y],
        "biome_suitability_lut": {
            BIOME_NAMES[k]: BIOME_SUITABILITY[k] for k in range(1, n_classes)
        },
        "decisions": BIOME_NOTES,
        "majority_filter": {
            "window": MAJORITY_WINDOW,
            "fraction_land_cells_changed": frac_changed,
            "area_km2_before": area_before,
            "area_km2_after": area_after,
        },
        "ice_rock_slope_crosscheck": ice_rock_check,
        "biome_suitability_land_stats": {
            "mean": float(suit[land].mean()),
            "p5_p50_p95": [float(v) for v in np.percentile(suit[land], [5, 50, 95])],
        },
    }
    with open(f"{out}/tappa6_biome_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
