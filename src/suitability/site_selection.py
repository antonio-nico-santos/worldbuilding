"""
Tappa 6 (final deliverable) -- greedy, size-aware Circulo site selection on
top of suitability_circulo_120m (or any other composite 0-1 surface).

This is the "settlement-size downstream filter" left open since the Tappa 6
architecture was first locked: rather than a separate suitability SURFACE
per settlement size, each Circulo's population is converted into a required
footprint AREA (via an assumed density), and site candidates are screened
by the MEAN suitability over a square window of that footprint size, not
just the single best point -- a large Circulo needs a large contiguous
patch of good land, not one lucky pixel surrounded by mediocre ground.

SCOPE, on purpose: this picks a SITE per Circulo (a point + an indicative
radius), not a footprint shape. The circular/organic layout each Circulo
actually has is explicitly Tappa 7's job ("urban zoom", deferred until
Tappa 6 picks a location) -- this module's square-window approximation of a
circular footprint (window side = diameter, window area ~1.27x the true
disk area) is an acceptable coarseness for THIS stage, not a shortcut being
hidden.

HONEST LIMITATIONS:
- DENSITY_PPL_KM2 (population -> required area) is an assumed constant, not
  derived from anything in this project or the scenario -- a "dense but
  green, walkable, with internal agriculture/parks" solarpunk town, denser
  than car-dependent suburbia (~10-30 ppl/ha) but far short of a packed
  historic core (Barcelona's Eixample ~360 ppl/ha). Changing it rescales
  every radius by 1/sqrt(density); nothing else in this module depends on
  its specific value.
- The 8 smallest Circulos' individual populations were never specified
  (only their 5,000-person TOTAL) -- split evenly here (625 each), a
  simplifying assumption, not a scenario fact.
- For those same smallest Circulos, the required window is only ~4x4 120 m
  cells (16 cells, ~0.23 km2) against a ~0.156 km2 target footprint -- the
  120 m grid is genuinely coarse relative to their size. Treat their
  reported site as indicative/approximate, not precise, more so than for
  the larger Circulos.
- Placement is GREEDY (largest population first, since large Circulos are
  the most area-constrained) and irrevocable -- once a site is chosen it is
  never revisited even if a later, smaller Circulo's placement would have
  produced a better global arrangement. A jointly-optimal placement (e.g.
  simulated annealing over all sites at once) was not attempted; greedy is
  standard practice for this class of problem and is transparent about
  which choice happened first.
- BUFFER_FACTOR controls spacing between Circulos: after a site is chosen,
  a disk of BUFFER_FACTOR x its own radius is marked claimed (off-limits to
  every later Circulo), so two settlements can't end up adjacent or
  overlapping. This is a placeholder "keep some countryside between
  settlements" choice, not a scenario-specified minimum distance.
- A candidate window must be >= MIN_LAND_FRACTION land (default 0.98) to be
  considered at all -- keeps sites off the coast where the window would
  otherwise average in ocean cells (which don't carry a suitability value).
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter

__all__ = [
    "DENSITY_PPL_KM2",
    "BUFFER_FACTOR",
    "MIN_LAND_FRACTION",
    "circulo_radius_km",
    "window_stats",
    "place_circulos",
]

DENSITY_PPL_KM2 = 4000.0
BUFFER_FACTOR = 2.0
MIN_LAND_FRACTION = 0.98


def circulo_radius_km(population: float, density_ppl_km2: float = DENSITY_PPL_KM2) -> float:
    """Equivalent-circle radius (km) of the land area a settlement of this
    population needs, at the given density."""
    area_km2 = population / density_ppl_km2
    return float(np.sqrt(area_km2 / np.pi))


def window_stats(
    suitability: np.ndarray, land_mask: np.ndarray, window_cells: int
) -> tuple[np.ndarray, np.ndarray]:
    """For a square window of `window_cells` side, return (mean_suitability
    over LAND cells only, land_fraction) at every pixel, via two box-filter
    passes (fast, exact, separable -- no need for a true circular kernel at
    this stage, see module docstring).
    """
    land_f = land_mask.astype(np.float64)
    land_frac = uniform_filter(land_f, size=window_cells, mode="nearest")
    suit_masked = np.where(land_mask, np.nan_to_num(suitability, nan=0.0), 0.0)
    suit_box_mean = uniform_filter(suit_masked, size=window_cells, mode="nearest")
    mean_suit_land = np.divide(
        suit_box_mean, land_frac, out=np.zeros_like(suit_box_mean), where=land_frac > 1e-9
    )
    return mean_suit_land, land_frac


def _mark_disk_claimed(claimed: np.ndarray, row: int, col: int, radius_cells: float) -> None:
    ny, nx = claimed.shape
    r = int(np.ceil(radius_cells))
    r0, r1 = max(0, row - r), min(ny, row + r + 1)
    c0, c1 = max(0, col - r), min(nx, col + r + 1)
    yy, xx = np.ogrid[r0:r1, c0:c1]
    dist2 = (yy - row) ** 2 + (xx - col) ** 2
    claimed[r0:r1, c0:c1] |= dist2 <= radius_cells**2


def place_circulos(
    suitability: np.ndarray,
    land_mask: np.ndarray,
    cellsize_km: float,
    circulos: list[tuple[str, float]],
    density_ppl_km2: float = DENSITY_PPL_KM2,
    buffer_factor: float = BUFFER_FACTOR,
    min_land_fraction: float = MIN_LAND_FRACTION,
) -> list[dict]:
    """Greedy site selection, largest population first. `circulos` is a
    list of (name, population) pairs -- order in the input does not matter,
    it is re-sorted internally by descending population.

    Returns a list of dicts (same order as input `circulos`, NOT placement
    order): name, population, radius_km, window_cells, row, col,
    mean_suitability, land_fraction, placed (False if no valid candidate
    remained -- e.g. every unclaimed site failed min_land_fraction).
    """
    ny, nx = land_mask.shape
    claimed = np.zeros((ny, nx), dtype=bool)
    order = sorted(range(len(circulos)), key=lambda i: -circulos[i][1])
    results: list[dict | None] = [None] * len(circulos)

    for i in order:
        name, population = circulos[i]
        radius_km = circulo_radius_km(population, density_ppl_km2)
        radius_cells = radius_km / cellsize_km
        window_cells = max(3, int(round(2 * radius_cells)))
        if window_cells % 2 == 0:
            window_cells += 1  # odd window -> well-defined center cell

        mean_suit, land_frac = window_stats(suitability, land_mask, window_cells)
        valid = (land_frac >= min_land_fraction) & ~claimed & land_mask
        # also require the window itself doesn't touch a claimed cell: a
        # claimed pixel anywhere in the window is disallowed, not just at
        # the center -- checked via the same box-filter trick (claimed
        # fraction in window must be exactly 0)
        claimed_frac, _ = window_stats(claimed.astype(np.float64), land_mask, window_cells)
        valid &= claimed_frac <= 1e-9

        if not valid.any():
            results[i] = {
                "name": name, "population": population, "radius_km": radius_km,
                "window_cells": window_cells, "placed": False,
            }
            continue

        score = np.where(valid, mean_suit, -1.0)
        flat_idx = int(np.argmax(score))
        row, col = np.unravel_index(flat_idx, score.shape)

        results[i] = {
            "name": name,
            "population": population,
            "radius_km": radius_km,
            "window_cells": window_cells,
            "row": int(row),
            "col": int(col),
            "mean_suitability": float(mean_suit[row, col]),
            "land_fraction": float(land_frac[row, col]),
            "placed": True,
        }
        _mark_disk_claimed(claimed, row, col, radius_cells * buffer_factor)

    return results
