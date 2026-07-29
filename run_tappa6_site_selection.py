"""
Tappa 6 (final deliverable) -- Circulo candidate site selection, greedy by
population, on top of suitability_circulo_120m. See
src/suitability/site_selection.py's module docstring for the full method,
the density assumption, and honest limitations (coarse resolution for the
smallest Circulos, greedy-not-jointly-optimal placement, square-window
approximation of a circular footprint pending Tappa 7's actual layout).

Reads:
  data/processed/suitability/suitability_circulo_120m.npy
  data/processed/climate/land_mask.npy
  data/processed/biomes/biome_id.npy          (informational only -- which
                                               biome each site lands in,
                                               a preview of the Tappa 7
                                               architecture-style lookup,
                                               NOT a formal Tappa 6 output)

Writes to data/processed/suitability/ (gitignored, regenerate locally):
  circulo_candidate_sites.geojson   (Point features, one per Circulo)
  circulo_claimed_footprint_120m    (bool raster, which cells got reserved
                                     by ANY Circulo's site+buffer -- a
                                     sanity-check visual, not a formal layer)
  tappa6_site_selection_meta.json
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

from src.biomes.world_biomes import BIOME_NAMES
from src.suitability.site_selection import (
    BUFFER_FACTOR,
    DENSITY_PPL_KM2,
    MIN_LAND_FRACTION,
    place_circulos,
)
from src.terrain.raster_io import write_envi_raw, write_prj

XMIN, XMAX = -65000.0, 65000.0
YMIN, YMAX = -80000.0, 80000.0
PROJ4 = (
    "+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 "
    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

# Populations from the Tappa 6 planning chat. The 8 smallest only had a
# combined total (5,000) specified, not individual sizes -- split evenly
# here (625 each), a simplifying assumption, see site_selection.py docstring.
CIRCULOS = [
    ("Circulo_A_40k", 40000),
    ("Circulo_B_35k", 35000),
    ("Circulo_C_25k", 25000),
    ("Circulo_D_20k", 20000),
] + [(f"Circulo_E{i+1}_2k", 2000) for i in range(5)] + [
    (f"Circulo_F{i+1}_small", 625) for i in range(8)
]


def main():
    t0 = time.time()
    suit = np.load("data/processed/suitability/suitability_circulo_120m.npy").astype(np.float64)
    land = np.load("data/processed/climate/land_mask.npy").astype(bool)
    biome_id = np.load("data/processed/biomes/biome_id.npy")

    ny, nx = land.shape
    cs_x = (XMAX - XMIN) / nx
    cs_y = (YMAX - YMIN) / ny
    cellsize_km = (cs_x + cs_y) / 2 / 1000.0

    results = place_circulos(suit, land, cellsize_km, CIRCULOS)

    features = []
    claimed_viz = np.zeros((ny, nx), dtype=bool)
    for r in results:
        if not r["placed"]:
            continue
        row, col = r["row"], r["col"]
        x = XMIN + (col + 0.5) * cs_x
        y = YMAX - (row + 0.5) * cs_y
        biome_name = BIOME_NAMES[int(biome_id[row, col])]
        features.append({
            "type": "Feature",
            "properties": {
                "name": r["name"],
                "population": r["population"],
                "radius_km": round(r["radius_km"], 3),
                "window_cells": r["window_cells"],
                "mean_suitability": round(r["mean_suitability"], 4),
                "land_fraction": round(r["land_fraction"], 4),
                "biome_at_site": biome_name,
            },
            "geometry": {"type": "Point", "coordinates": [x, y]},
        })
        rc = int(np.ceil(r["radius_km"] / cellsize_km * BUFFER_FACTOR))
        r0, r1 = max(0, row - rc), min(ny, row + rc + 1)
        c0, c1 = max(0, col - rc), min(nx, col + rc + 1)
        yy, xx = np.ogrid[r0:r1, c0:c1]
        dist2 = (yy - row) ** 2 + (xx - col) ** 2
        claimed_viz[r0:r1, c0:c1] |= dist2 <= (r["radius_km"] / cellsize_km * BUFFER_FACTOR) ** 2

    geojson = {
        "type": "FeatureCollection",
        "name": "circulo_candidate_sites",
        "crs": {"type": "proj4", "properties": {"proj4": PROJ4}},
        "features": features,
    }

    out = "data/processed/suitability"
    os.makedirs(out, exist_ok=True)
    with open(f"{out}/circulo_candidate_sites.geojson", "w") as f:
        json.dump(geojson, f, indent=2)

    claimed_i2 = claimed_viz.astype(np.int16)
    np.save(f"{out}/circulo_claimed_footprint_120m.npy", claimed_i2)
    write_envi_raw(
        f"{out}/circulo_claimed_footprint_120m", claimed_i2, XMIN, YMIN, cs_x,
        "Tappa6 circulo_claimed_footprint_120m", dtype="i2",
    )
    write_prj(f"{out}/circulo_claimed_footprint_120m.prj", PROJ4)

    n_placed = sum(1 for r in results if r["placed"])
    n_failed = len(results) - n_placed
    meta = {
        "assumptions": {
            "density_ppl_km2": DENSITY_PPL_KM2,
            "buffer_factor": BUFFER_FACTOR,
            "min_land_fraction": MIN_LAND_FRACTION,
            "note": "8 smallest Circulos' individual populations were not specified "
            "(only their 5,000 combined total) -- split evenly (625 each) here.",
        },
        "resolution_m": [cs_x, cs_y],
        "n_circulos": len(CIRCULOS),
        "n_placed": n_placed,
        "n_failed_to_place": n_failed,
        "sites": [
            {k: v for k, v in r.items() if k not in ("row", "col")} for r in results
        ],
    }
    with open(f"{out}/tappa6_site_selection_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Done in {time.time() - t0:.1f}s -- {n_placed}/{len(CIRCULOS)} placed, {n_failed} failed")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
