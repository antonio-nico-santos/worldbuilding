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
- BUFFER_FACTOR controls the DEFAULT spacing every Circulo gets, tier rule
  or not: each site enforces at least radius_km * BUFFER_FACTOR straight-
  line separation from every other already-placed site (whichever site's
  own buffer is larger wins), so nothing ever overlaps or lands right on
  top of something else even without an explicit tier_min_distance rule.
  This is a placeholder "keep some countryside between settlements"
  choice, not a scenario-specified minimum distance.
- Per-tier min_distance (added after Nico's "greek city-state" note, later
  UPGRADED from straight-line km to cost-distance HOURS -- see below): each
  Circulo can carry a `tier` label ("large"/"medium"/"small"/None). Two
  SEPARATE minimum-distance mechanisms now coexist, on purpose:

  1. `buffer_factor` / own footprint non-overlap (straight-line, km,
     ALWAYS enforced regardless of tier): a purely physical constraint --
     two settlements' built-up areas cannot literally overlap in space, no
     matter how hard the terrain between them is. Small and cheap, a few
     km at most.
  2. `tier_min_hours` (cost-distance, HOURS, only between sites that BOTH
     have a tier AND whose tier-pair is a key in the dict): a
     socio-cultural/travel-time separation intent -- "how independent/
     self-sufficient does this settlement need to feel from its peers" --
     which is what "greek city-states" and "not practically integrated"
     are actually about, not raw map distance. See cost_distance.py for
     the underlying Tobler's-hiking-function + boat-crossing cost graph.
     A cross-tier pair (e.g. large-medium), or either site untiered, is
     NOT covered by an explicit rule and only gets mechanism 1 above.

  This REPLACES the earlier straight-line-km tier_min_distance (60 km
  large-large / 25 km medium-medium) with cost-distance hour thresholds,
  after Nico asked to test whether that changes the same-biome clustering
  outcome, and to also give the small Circulos (previously untiered) their
  own minimum -- "2-3 horas de caminhada, senao e praticamente integrado" --
  which is naturally an HOURS quantity, not a km one. See
  run_tappa6_site_selection.py for the specific tested/chosen hour values
  per tier and the feasibility sweep that picked them.

  Earlier history, kept for context: the first straight-line tier_min_distance
  attempt (checking each site's OWN minimum against every already-placed
  site regardless of tier) collapsed placement from 17/17 to 5/17 -- a
  60 km-radius disk (~11,310 km2) is LARGER than this world's entire land
  area (~9,904 km2). Fixed by scoping the rule to same-tier pairs only.
  That km-based, same-tier-pair version is what got 17/17 placed with 0
  violations before this cost-distance upgrade -- see git history /
  the earlier delivered version of this docstring for those numbers.
- Distance for mechanism 1 (footprint non-overlap) is straight-line
  (Euclidean) in the LCC plane -- a genuinely physical, terrain-independent
  constraint (a footprint's built area is the same size regardless of what
  is between it and its neighbour). Distance for mechanism 2 (tier rules)
  is cost-distance (hours), via a graph spanning the WHOLE domain (land
  walking at Tobler's-function speed, any edge touching non-land at a flat
  boat speed -- see cost_distance.py's docstring for the exact model and
  its own placeholder assumptions, esp. BOAT_SPEED_KMH). Islands are NOT
  excluded from candidacy either way: land_mask includes them like any
  other land, and checked directly, the SW Island (label 135, connected-
  component analysis) has a HIGHER mean suitability_circulo (0.593) than
  the mainland (0.564) -- a live candidate, not a token inclusion.
- A candidate window must be >= MIN_LAND_FRACTION land (default 0.98) to be
  considered at all -- keeps sites off the coast where the window would
  otherwise average in ocean cells (which don't carry a suitability value).
- Tier-rule ENFORCEMENT (the actual placement constraint) is a per-pair
  lookup and was verified directly to be correct (0 violations at the
  hour thresholds this project actually uses). Do not confuse that with
  VISUALIZING a tier's reach: a single large Circulo's isochrone (all
  cells reachable within its own tier's hour threshold) already covers
  roughly a quarter of this world's land on its own, and the union of all
  17 sites' isochrones covers effectively the entire grid (99.998%,
  checked directly) -- not a useful "what area is still open" picture.
  run_tappa6_site_selection.py's claimed_footprint output therefore only
  shades each site's small straight-line footprint-buffer circle
  (mechanism 1), not the tier isochrones, for exactly this reason.
"""

from __future__ import annotations

from typing import Callable

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


def place_circulos(
    suitability: np.ndarray,
    land_mask: np.ndarray,
    cellsize_km: float,
    circulos: list[tuple[str, float]] | list[tuple[str, float, str | None]],
    tier_min_hours: dict[frozenset, float] | None = None,
    cost_graph=None,
    cost_distance_fn: Callable | None = None,
    density_ppl_km2: float = DENSITY_PPL_KM2,
    buffer_factor: float = BUFFER_FACTOR,
    min_land_fraction: float = MIN_LAND_FRACTION,
    xmin_km: float = 0.0,
    ymax_km: float = 0.0,
) -> list[dict]:
    """Greedy site selection, largest population first. `circulos` is a
    list of (name, population) or (name, population, tier) triples -- order
    in the input does not matter, it is re-sorted internally by descending
    population. `tier` is any hashable label (e.g. "large"/"medium"/"small")
    or None for "no tier rule, just avoid crowding".

    Two independent minimum-distance mechanisms (see module docstring for
    the full rationale):
    1. Footprint non-overlap, ALWAYS enforced, straight-line km: each site
       requires >= max(own radius_km, other's radius_km) * buffer_factor
       separation from every already-placed site, tier or no tier.
    2. `tier_min_hours` maps frozenset({tier_a, tier_b}) -> required
       COST-DISTANCE minimum separation, in HOURS, between a site of
       tier_a and one of tier_b (pass frozenset({"large"}) i.e. a
       single-element set for a same-tier rule). Only checked between
       sites that BOTH have a tier AND whose tier-pair is a key in this
       dict -- a cross-tier pair or an untiered site only gets mechanism 1.
       Requires `cost_graph` (a scipy.sparse cost graph, see
       cost_distance.py's build_cost_graph) and `cost_distance_fn` (a
       callable (row, col) -> full-grid hours array from that source,
       e.g. functools.partial(cost_distance_from_source, cost_graph,
       shape=land_mask.shape)) to be given; if either is None,
       tier_min_hours is ignored (mechanism 1 alone still applies).

    `xmin_km`/`ymax_km`: real-world offset (km) of the grid's (row=0,col=0)
    corner, so returned/compared coordinates are in real map km, not grid
    indices. Pass 0.0/0.0 (default) to work purely in grid-index space.

    Returns a list of dicts (same order as input `circulos`, NOT placement
    order): name, population, tier, radius_km, window_cells, row, col,
    x_km, y_km, mean_suitability, land_fraction, placed (False if no valid
    candidate remained).
    """
    ny, nx = land_mask.shape
    tier_min_hours = tier_min_hours or {}
    use_cost = cost_graph is not None and cost_distance_fn is not None
    entries = [(c[0], c[1], c[2] if len(c) > 2 else None) for c in circulos]
    order = sorted(range(len(entries)), key=lambda i: -entries[i][1])
    results: list[dict | None] = [None] * len(entries)
    # (x_km, y_km, tier, footprint_min_km, cost_hours_raster_or_None) per
    # already-placed site -- the cost raster is this SITE's own full-grid
    # hours-from-here array, computed once when it was placed, reused by
    # every later (smaller) Circulo's candidate check.
    placed_sites: list[tuple[float, float, object, float, np.ndarray | None]] = []

    xc_km = xmin_km + (np.arange(nx) + 0.5) * cellsize_km
    yc_km = ymax_km - (np.arange(ny) + 0.5) * cellsize_km
    Xg, Yg = np.meshgrid(xc_km, yc_km)

    for i in order:
        name, population, tier = entries[i]
        radius_km = circulo_radius_km(population, density_ppl_km2)
        radius_cells = radius_km / cellsize_km
        window_cells = max(3, int(round(2 * radius_cells)))
        if window_cells % 2 == 0:
            window_cells += 1  # odd window -> well-defined center cell
        own_footprint_min_km = radius_km * buffer_factor

        mean_suit, land_frac = window_stats(suitability, land_mask, window_cells)
        valid = (land_frac >= min_land_fraction) & land_mask

        for (px, py, other_tier, other_footprint_min_km, other_cost_raster) in placed_sites:
            required = max(own_footprint_min_km, other_footprint_min_km)
            dist = np.hypot(Xg - px, Yg - py)
            valid &= dist >= required

            if use_cost and tier is not None and other_tier is not None:
                hr_req = tier_min_hours.get(frozenset({tier, other_tier}))
                if hr_req is not None and other_cost_raster is not None:
                    valid &= other_cost_raster >= hr_req

        if not valid.any():
            results[i] = {
                "name": name, "population": population, "tier": tier,
                "radius_km": radius_km, "window_cells": window_cells, "placed": False,
            }
            continue

        score = np.where(valid, mean_suit, -1.0)
        flat_idx = int(np.argmax(score))
        row, col = np.unravel_index(flat_idx, score.shape)
        x_km, y_km = float(Xg[row, col]), float(Yg[row, col])

        results[i] = {
            "name": name,
            "population": population,
            "tier": tier,
            "radius_km": radius_km,
            "window_cells": window_cells,
            "row": int(row),
            "col": int(col),
            "x_km": x_km,
            "y_km": y_km,
            "mean_suitability": float(mean_suit[row, col]),
            "land_fraction": float(land_frac[row, col]),
            "placed": True,
        }
        # only bother computing the (cheap, but not free) cost-distance
        # raster for this site if it's tiered and cost-distance is in use
        # at all -- an untiered site never participates in a tier rule
        own_cost_raster = (
            cost_distance_fn(row, col) if (use_cost and tier is not None) else None
        )
        placed_sites.append((x_km, y_km, tier, own_footprint_min_km, own_cost_raster))

    return results
