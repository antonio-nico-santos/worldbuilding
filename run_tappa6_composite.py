"""
Tappa 6 (final step) -- weighted composite suitability, one per population
(Circulo / Povo Livre), at the 120 m working grid. See
src/suitability/composite.py's module docstring for the weight tables and
the reasoning behind them (including the honest limitations already flagged
in the layers this combines -- water_suitability's near-flatness, solar's
annual-vs-seasonal mismatch for a nomadic population, etc.).

Reads (all from data/processed/, 120 m, already delivered by the earlier
Tappa 6 runners):
  suitability/slope_suitability_120m.npy
  suitability/water_suitability_120m.npy
  suitability/agriculture_suitability_120m.npy
  suitability/solar_suitability_annual_120m.npy
  suitability/biome_suitability_120m.npy
  suitability/biome_suitability_povo_livre_120m.npy
  suitability/povo_silencioso_exclusion_120m.npy
  climate/land_mask.npy

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  suitability_circulo_120m       (0-1, float32)
  suitability_povo_livre_120m    (0-1, float32)
  tappa6_composite_meta.json
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

from src.suitability.composite import WEIGHTS_CIRCULO, WEIGHTS_POVO_LIVRE, weighted_composite
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)


def main():
    t0 = time.time()
    src = "data/processed/suitability"
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)

    layers = {
        "slope": np.load(f"{src}/slope_suitability_120m.npy").astype(np.float64),
        "water": np.load(f"{src}/water_suitability_120m.npy").astype(np.float64),
        "agriculture": np.load(f"{src}/agriculture_suitability_120m.npy").astype(np.float64),
        "solar": np.load(f"{src}/solar_suitability_annual_120m.npy").astype(np.float64),
    }
    biome_circulo = np.load(f"{src}/biome_suitability_120m.npy").astype(np.float64)
    biome_povo_livre = np.load(f"{src}/biome_suitability_povo_livre_120m.npy").astype(np.float64)
    exclusion = np.load(f"{src}/povo_silencioso_exclusion_120m.npy").astype(np.float64)

    ny, nx = land.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny

    suit_circulo = weighted_composite(
        {**layers, "biome": biome_circulo}, WEIGHTS_CIRCULO, land, exclusion
    ).astype(np.float32)
    suit_povo_livre = weighted_composite(
        {**layers, "biome": biome_povo_livre}, WEIGHTS_POVO_LIVRE, land, exclusion
    ).astype(np.float32)

    out = src
    os.makedirs(out, exist_ok=True)
    for name, arr in [
        ("suitability_circulo_120m", suit_circulo),
        ("suitability_povo_livre_120m", suit_povo_livre),
    ]:
        np.save(f"{out}/{name}.npy", arr)
        write_envi_raw(f"{out}/{name}", arr, XMIN, YMIN, cs_x, f"Tappa6 {name}", dtype="f4")
        write_prj(f"{out}/{name}.prj", PROJ4)

    # spatial comparison, same style check as the biome variant comparison --
    # confirm the two composites actually favour different regions, not just
    # different scores in the same place
    xc = XMIN + (np.arange(nx) + 0.5) * cs_x
    X = np.broadcast_to(xc, (ny, nx))
    top_circulo = land & (suit_circulo >= np.nanpercentile(suit_circulo[land], 90))
    top_povo_livre = land & (suit_povo_livre >= np.nanpercentile(suit_povo_livre[land], 90))
    overlap = (top_circulo & top_povo_livre).sum()

    meta = {
        "resolution_m": [cs_x, cs_y],
        "weights_circulo": WEIGHTS_CIRCULO,
        "weights_povo_livre": WEIGHTS_POVO_LIVRE,
        "povo_silencioso_exclusion_applied_to_both": True,
        "suitability_circulo_land_stats": {
            "mean": float(suit_circulo[land].mean()),
            "p5_p50_p95": [float(v) for v in np.percentile(suit_circulo[land], [5, 50, 95])],
        },
        "suitability_povo_livre_land_stats": {
            "mean": float(suit_povo_livre[land].mean()),
            "p5_p50_p95": [float(v) for v in np.percentile(suit_povo_livre[land], [5, 50, 95])],
        },
        "top_decile_comparison": {
            "circulo_top10pct_mean_x_m": float(X[top_circulo].mean()),
            "povo_livre_top10pct_mean_x_m": float(X[top_povo_livre].mean()),
            "overlap_cells": int(overlap),
            "overlap_pct_of_circulo_top": float(100 * overlap / top_circulo.sum()),
            "overlap_pct_of_povo_livre_top": float(100 * overlap / top_povo_livre.sum()),
        },
    }
    with open(f"{out}/tappa6_composite_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
