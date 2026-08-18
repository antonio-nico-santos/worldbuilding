"""
Tappa 8 -- Geomorphology: cave candidates.

Four types confirmed for this pass (docs/decisions/07_tappa7_regional_scenario.md
S2); schist fracture caves stay explicitly NOT built (parked pending visual
review of these four -- see that doc). All four are candidate ELIGIBILITY
masks (boolean, "could plausibly host this cave type here"), not a curated
cave catalogue -- same honesty framing Tappa 4 already used for
`lake_mask` (S4: "a floor against tiny numerical fill noise, not a curated
lake catalogue").

- lava tubes: volcanic lithology class only (S2 -- strong citation, Banks
  Peninsula). Optionally weighted by proximity to a known geothermal vent
  (data/input/geothermal.geojson, the 3 points authored for Tappa 6's
  geothermal hazard layer) -- lava tubes form along actual flow conduits
  from a vent, so distance-to-vent is a real plausibility signal, not
  present in S2's text. Flagged as this module's own addition.
- talus/pseudokarst: steep MAINLAND relief (excludes the SW Island and
  every other islet -- S2 says "mainland") intersected with a stream
  cutting beneath it. "Steep" = top quartile of mainland land slope (same
  method already used for Nacre's cave candidates in Tappa 7 S7's v6
  addendum -- "top quartile within the alpine band"); "near a stream"
  reuses Tappa 6's own 0.5 km `water_gentle_km` threshold directly (S1 of
  that doc: "0.5 km ~ where trivially close stops covering nearly
  everyone") rather than inventing a new distance.
- glacier/moulin ice: Tappa 3's `permanent_snow_mask` directly, zero new
  geometry (S2) -- resampled from its native 120 m climate grid onto this
  stage's native-30m grid by nearest-neighbor coordinate lookup, the same
  cross-grid method already used repeatedly in Tappa 7 (biome_id vs
  stream_mask, dist_to_subalpine, etc.) since the two grids' cell counts
  don't divide evenly.
- sea caves: land within a coastal buffer of the 0 m contour, reusing the
  SAME steep-slope threshold as talus caves (this stage's own locked
  decision: reuse the existing slope field on a coastal buffer rather than
  inventing a new local-relief statistic). Buffer width and the shared
  slope threshold are both placeholders, not independently calibrated --
  see the decision doc.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from suitability.terrain_metrics import compute_slope_pct


def resample_nearest(coarse: np.ndarray, coarse_xmin, coarse_ymax, coarse_cellsize, ny, nx, xmin, ymax, cellsize):
    """Nearest-neighbor lookup from a coarser grid onto this module's
    native grid, by real-world coordinate (not an index ratio) -- the two
    grids' cell counts don't divide evenly (5334x4334 vs 1334x1084, a
    padding difference between the DEM/hydrology pipeline and the climate
    pipeline, not a bug -- same note already on record in Tappa 7 for the
    identical situation)."""
    y = ymax - (np.arange(ny) + 0.5) * cellsize
    x = xmin + (np.arange(nx) + 0.5) * cellsize
    row_idx = np.clip(((coarse_ymax - y) / coarse_cellsize).astype(int), 0, coarse.shape[0] - 1)
    col_idx = np.clip(((x - coarse_xmin) / coarse_cellsize).astype(int), 0, coarse.shape[1] - 1)
    return coarse[np.ix_(row_idx, col_idx)]


def lava_tube_candidates(lithology: np.ndarray, land_mask: np.ndarray):
    from .lithology import CLASS_VOLCANIC

    return (lithology == CLASS_VOLCANIC) & land_mask


def lava_tube_candidates_v2(
    lithology: np.ndarray,
    land_mask: np.ndarray,
    vent_weight: np.ndarray,
    vent_weight_threshold: float,
    slope_pct: np.ndarray,
    gentle_threshold_pct: float,
):
    """v2 -- narrows v1's "entire volcanic class" candidate (794.29 km2,
    literally the whole island) using two real lava-tube plausibility
    signals instead of the bare lithology test:

    - `vent_weight >= vent_weight_threshold`: real lava tubes are conduits
      FROM an actual eruptive vent, not a property of "being on volcanic
      rock" anywhere on the island. `vent_weight` was already computed in
      v1 (Gaussian falloff from the 3 geothermal.geojson vents) but never
      actually used as a filter -- this is that fix.
    - `slope_pct <= gentle_threshold_pct`: the INVERSE of talus/sea caves'
      steep-preference test. Talus/sea want steep relief (mass-wasting,
      wave-cut cliffs); lava tubes are the opposite -- they form in
      channelized surface flows down gentle-to-moderate flanks, and don't
      survive on cliff faces or crater rims. Mirrors talus/sea's own
      "top quartile" convention exactly, just inverted to "bottom
      quartile" of the volcanic zone's own slope population.

    Both thresholds are first-pass placeholders (p25 slope within the
    volcanic zone / vent_weight>=0.1 in the driver), not independently
    calibrated -- see the checkpoint write-up."""
    from .lithology import CLASS_VOLCANIC

    return (
        (lithology == CLASS_VOLCANIC)
        & land_mask
        & (vent_weight >= vent_weight_threshold)
        & (slope_pct <= gentle_threshold_pct)
    )


def vent_proximity_weight(xx, yy, vents, cellsize_m):
    """Optional plausibility weight in [0, 1] for lava tube candidates:
    max over vents of a Gaussian falloff using each vent's own
    `falloff_km` (geothermal.geojson properties), the exact same
    half-max-distance convention Tappa 1's ridge decay uses. This is an
    ADDITION beyond S2's literal text (which only names the volcanic
    class), not itself the eligibility test."""
    ln2 = np.log(2.0)
    weight = np.zeros(xx.shape, dtype=np.float32)
    for vx, vy, falloff_km in vents:
        d = np.hypot(xx - vx, yy - vy)
        w = np.exp(-ln2 * (d / (falloff_km * 1000.0)) ** 2)
        weight = np.maximum(weight, w.astype(np.float32))
    return weight


def talus_pseudokarst_candidates(
    slope_pct: np.ndarray,
    dist_to_stream_km: np.ndarray,
    mainland_mask: np.ndarray,
    steep_threshold_pct: float,
    stream_buffer_km: float = 0.5,
):
    return mainland_mask & (slope_pct >= steep_threshold_pct) & (dist_to_stream_km <= stream_buffer_km)


def glacier_moulin_candidates(permanent_snow_mask_native: np.ndarray, land_mask: np.ndarray):
    return permanent_snow_mask_native.astype(bool) & land_mask


def glacier_moulin_candidates_v2(
    permanent_snow_mask_native: np.ndarray,
    land_mask: np.ndarray,
    slope_pct: np.ndarray,
    steep_threshold_pct: float,
    margin_depth_km: np.ndarray,
    margin_depth_threshold_km: float,
):
    """v2 -- narrows v1's "whole permanent-snow-mask" candidate (960.16 km2,
    Nico's "everything included, nothing included" complaint) using two
    real glaciological signals instead of the raw mask:

    - `slope_pct >= steep_threshold_pct`: moulins form where meltwater
      finds crevasses to descend into, and crevassing concentrates on
      steeper ice, not flat snowfield/ice-cap interiors. SAME method as
      talus/sea caves (top quartile of the relevant land population's
      slope) -- here that population is the snow mask itself, not all
      mainland land.
    - `margin_depth_km <= margin_depth_threshold_km`: moulins also
      concentrate in the ablation zone near the ice margin where melt is
      actively flowing, not deep in the accumulation-zone interior.
      `margin_depth_km` is each snow cell's own distance-to-nearest-non-
      snow-cell (computed by the driver via distance_transform_edt on the
      snow mask itself) -- "how far into the interior this cell sits."

    Both thresholds are first-pass placeholders (p75 slope / median depth
    in the driver), not independently calibrated -- see the checkpoint
    write-up."""
    return (
        permanent_snow_mask_native.astype(bool)
        & land_mask
        & (slope_pct >= steep_threshold_pct)
        & (margin_depth_km <= margin_depth_threshold_km)
    )


def sea_cave_candidates(
    land_mask: np.ndarray,
    dist_to_ocean_km: np.ndarray,
    slope_pct: np.ndarray,
    steep_threshold_pct: float,
    coastal_buffer_km: float = 0.5,
):
    return land_mask & (dist_to_ocean_km <= coastal_buffer_km) & (slope_pct >= steep_threshold_pct)


def distance_to_ocean_km(land_mask: np.ndarray, cellsize_m: float):
    """Euclidean distance from each LAND cell to the nearest ocean cell --
    the mirror of every other distance_transform_edt usage in this
    project (Tappa 2 continentality, Tappa 6 stream/Povo Silencioso
    distance), just computed on the ocean mask instead of a feature mask,
    and at this stage's native 30 m rather than the 120 m climate grid
    those used."""
    return ndimage.distance_transform_edt(land_mask, sampling=(cellsize_m, cellsize_m)) / 1000.0


def compute_slope(dem: np.ndarray, cellsize_m: float):
    return compute_slope_pct(dem, cellsize_m)
