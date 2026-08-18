"""
Tappa 8 -- lithology v5: DEM-NATIVE classification, no authored ridge/zone
geometry involved in the decision at all (Nico's catch, with the attached
real NZ geology reference map as the clue: the Otago/Haast Schist body is
an irregular BLOB tied to actual relief/uplift, not a uniform-width band
that hugs a drawn ridge line -- and it does NOT respect "plateau breaks the
spine logic": a flat high plateau next to a mountain range reads as
sedimentary basin fill in the real map, regardless of how close it sits to
the range).

v2-v4 all shared one structural flaw despite fixing the crest-GEOMETRY
progressively: schist/greywacke were always defined as "within some
distance of a line" (authored, or DEM-grounded, or a whole real ridge
NETWORK) -- so a flat authored plateau sitting geometrically close to a
ridge's footprint (the SE plains zone, concretely) still got painted
schist/greywacke by the falloff rule, even though it's flat terrain that
should read as basin fill. Distance-to-line can never fix that, no matter
how the line itself is derived, because flatness isn't a positional fact.

v5 instead classifies directly off two DEM-derived fields, no line/corridor
of any kind:
- elevation (dem itself)
- local relief (windowed max-min elevation, `compute_local_relief` below) --
  a real ruggedness signal: a flat high plateau (Central plateau, 1392m
  mean elevation) reads LOW on this even though it's high, exactly the
  distinction "distance to a ridge" can never make.

schist = high elevation AND high local relief (real exhumed mountainous
core); greywacke = moderate relief flank, not already schist; basin fill =
everything else (flat, regardless of elevation -- this is what fixes the
SE plains case). Thresholds are percentiles of MAINLAND land's own
elevation/relief distributions (same "calibrate against this world's own
data" convention already used for jade's grade percentile, the biome
moisture terciles, etc.), not absolute values.

Bonus, not the reason this was built but worth stating: because the class
boundary is now a threshold on two raw DEM fields instead of a distance
field, it inherits ALL of the real terrain's own irregularity for free --
no synthetic noise warp needed (v4's boundary_noise.py is not used by this
version). This is arguably a more direct fix for the original "too smooth,
no noise" complaint than any of v2-v4's geometry changes were.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def compute_local_relief(dem: np.ndarray, cellsize_m: float, window_m: float = 2000.0):
    """Windowed local relief = max(dem) - min(dem) in a square window of
    side `window_m`, via separable max/min filters (fast: O(n) per pixel
    regardless of window size, scipy's van Herk/Gil-Werman implementation).
    NOT restricted to land -- ocean cells (elevation <= 0) will show up as
    high "relief" near the coastline; callers mask to land afterward."""
    win_cells = int(round(window_m / cellsize_m))
    if win_cells % 2 == 0:
        win_cells += 1
    win_cells = max(3, win_cells)
    mx = ndimage.maximum_filter(dem, size=win_cells, mode="nearest")
    mn = ndimage.minimum_filter(dem, size=win_cells, mode="nearest")
    return (mx - mn).astype(np.float32)


def classify_from_terrain(
    dem: np.ndarray,
    land_mask: np.ndarray,
    mainland_mask: np.ndarray,
    volcanic_mask: np.ndarray,
    relief: np.ndarray,
    elev_percentile: float = 60.0,
    relief_schist_percentile: float = 75.0,
    relief_greywacke_percentile: float = 50.0,
    high_elev_percentile: float | None = 90.0,
    island_relief_percentile: float | None = 25.0,
):
    """Percentile thresholds computed over MAINLAND LAND ONLY for the
    schist/greywacke/basin split.

    `island_relief_percentile`: Nico's second catch -- the SW Island was
    unconditionally 100% volcanic in every version through v5, a Tappa 7
    scope lock (S1), not something this relief redesign touched. Real
    Banks Peninsula (the cited analog) is mostly-but-not-entirely volcanic
    bedrock -- it has a real flat isthmus/harbour-margin apron. Applying
    the SAME kind of relief test used on the mainland, but calibrated
    against the ISLAND's OWN elevation/relief population (its whole
    elevation range, 0-737m, sits below every mainland threshold, so
    reusing mainland thresholds would make the entire island basin_fill,
    which is wrong the other direction): island land with relief below
    this percentile (of ISLAND land only) becomes basin_fill (the flat
    coastal apron/isthmus reading); at or above stays volcanic. No
    greywacke/schist tier on the island -- a young monogenetic shield
    doesn't have that metamorphic-flank story, just "volcanic edifice" vs
    "flat sedimentary apron at its margin." Set to None to keep the
    island 100% volcanic (pre-this-request behavior).

    `high_elev_percentile`: Nico's catch -- sampling lithology_v5 directly
    along the AUTHORED ridge crest lines showed 33-47% of the actual named
    crest classified as basin_fill, at elevations up to 3350m. Root cause,
    verified numerically: this world's ridges (Tappa 1's Gaussian falloff
    decay) generate SMOOTHLY ROUNDED summits, not knife edges -- a fixed
    2km window centered exactly on a rounded high point often shows LESS
    relief than a window centered on the steep flank just below it (the
    window can't "see" the full up-and-down swing from a true extremum the
    way it can from a mid-slope position). Sampled basin_fill-classified
    crest points had elevation 2700-3350m but only 159-265m of local
    relief -- below even the greywacke threshold. The elevation+relief AND
    test alone can never recover these, no matter how the relief window is
    tuned, because the terrain really is locally smooth exactly there.

    Fix: an OR escape hatch -- elevation alone, above a MUCH higher bar
    than the AND-branch's elev_percentile, also qualifies for schist
    regardless of relief. Checked this doesn't reopen the SE-plains/plateau
    problem before adding it: mainland p90 elevation is ~2612m, while the
    highest authored plateau (Central plateau) tops out at 2463m and every
    other plateau zone is far lower -- so this branch cannot be triggered
    by any authored flat zone, only by genuinely high, DEM-real terrain.
    Set `high_elev_percentile=None` to disable (pre-fix v5 behavior).
    """
    # thresholds are CALIBRATED against mainland-only population (avoids
    # contaminating the percentiles with the volcanic island's own
    # elevation/relief character), but then APPLIED to every non-volcanic
    # land cell, including the ~200+ minor islets -- restricting both
    # calibration AND application to mainland_mask left every minor islet
    # cell unclassified (silently defaulted to the ocean class code, a real
    # bug caught before this was shown: 268,101 land cells / 241 km2 with
    # no lithology at all).
    mainland_land = land_mask & mainland_mask
    elev_pop = dem[mainland_land]
    relief_pop = relief[mainland_land]

    elev_threshold = float(np.percentile(elev_pop, elev_percentile))
    relief_thresh_schist = float(np.percentile(relief_pop, relief_schist_percentile))
    relief_thresh_grey = float(np.percentile(relief_pop, relief_greywacke_percentile))

    classifiable_land = land_mask & ~volcanic_mask
    schist_and = (dem >= elev_threshold) & (relief >= relief_thresh_schist)
    if high_elev_percentile is not None:
        high_elev_threshold = float(np.percentile(elev_pop, high_elev_percentile))
        schist_or_high = dem >= high_elev_threshold
    else:
        high_elev_threshold = None
        schist_or_high = np.zeros_like(schist_and)
    schist_mask = classifiable_land & (schist_and | schist_or_high)
    greywacke_mask = classifiable_land & ~schist_mask & (relief >= relief_thresh_grey)
    basin_mask = classifiable_land & ~schist_mask & ~greywacke_mask

    from .lithology import CLASS_OCEAN, CLASS_BASIN_FILL, CLASS_GREYWACKE, CLASS_SCHIST, CLASS_VOLCANIC

    # island: same relief-based logic, own calibration population
    island_land = land_mask & volcanic_mask
    island_relief_threshold = None
    island_volcanic_mask = island_land
    island_basin_mask = np.zeros_like(island_land)
    if island_relief_percentile is not None and island_land.any():
        island_relief_pop = relief[island_land]
        island_relief_threshold = float(np.percentile(island_relief_pop, island_relief_percentile))
        island_basin_mask = island_land & (relief < island_relief_threshold)
        island_volcanic_mask = island_land & (relief >= island_relief_threshold)

    lithology = np.full(dem.shape, CLASS_OCEAN, dtype=np.uint8)
    lithology[basin_mask] = CLASS_BASIN_FILL
    lithology[greywacke_mask] = CLASS_GREYWACKE
    lithology[schist_mask] = CLASS_SCHIST
    lithology[island_basin_mask] = CLASS_BASIN_FILL
    lithology[island_volcanic_mask] = CLASS_VOLCANIC

    # rank-based grade within the schist population: 0.5*elevation_rank +
    # 0.5*relief_rank, both ranked ONLY among schist cells themselves (not
    # all mainland land) so the grade field spans the full 0-1 range across
    # whatever the schist body turns out to contain
    schist_grade = np.full(dem.shape, np.nan, dtype=np.float32)
    if schist_mask.any():
        elev_s = dem[schist_mask]
        relief_s = relief[schist_mask]
        elev_rank = np.argsort(np.argsort(elev_s)) / max(1, len(elev_s) - 1)
        relief_rank = np.argsort(np.argsort(relief_s)) / max(1, len(relief_s) - 1)
        grade_s = (0.5 * elev_rank + 0.5 * relief_rank).astype(np.float32)
        schist_grade[schist_mask] = grade_s

    return {
        "lithology": lithology,
        "schist_grade": schist_grade,
        "elev_threshold_m": elev_threshold,
        "relief_schist_threshold_m": relief_thresh_schist,
        "relief_greywacke_threshold_m": relief_thresh_grey,
        "high_elev_threshold_m": high_elev_threshold,
        "n_schist_from_high_elev_only": int((classifiable_land & schist_or_high & ~schist_and).sum()),
        "island_relief_threshold_m": island_relief_threshold,
        "island_volcanic_km2_cells": int(island_volcanic_mask.sum()),
        "island_basin_km2_cells": int(island_basin_mask.sum()),
    }
