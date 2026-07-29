"""
Tappa 6 -- core (non-exclusion) suitability criteria derived from terrain and
hydrology: slope, distance-to-stream, and a soil-free agriculture proxy
(Topographic Wetness Index). Also provides the Povo Silencioso exclusion
distance field, since it reuses the same block-reduce/EDT machinery.

Working resolution: 120 m, matching the climate/biome grid (Tappa 2 S2,
Tappa 5) rather than the DEM/hydrology native 30 m (Tappa 4) -- Tappa 6 is a
siting decision, not a per-building layout, so 120 m is the right level for
the FINAL suitability index. Sources that only exist at 30 m (the DEM,
contributing_area_km2, stream_mask) are downsampled here.

Downsampling convention: reuses the edge-replicate-to-a-multiple-of-4 padding
from src/climate/grid.py's block_mean, but slope and stream presence use
block_max / block_any instead of block_mean. Averaging first would let a
locally steep sub-cell or a thin stream disappear inside an otherwise gentle/
dry 120 m cell -- exactly the wrong direction of error for criteria meant to
flag hazards or resources. This is the same reasoning Tappa 2 S2 used for
"clip bathymetry before averaging, not after" -- averaging can hide the thing
you were trying to measure if applied before the operation that cares about
extremes, not after.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

__all__ = [
    "block_max",
    "block_any",
    "block_mean",
    "compute_slope_pct",
    "slope_suitability",
    "distance_to_stream_km",
    "water_suitability",
    "topographic_wetness_index",
    "agriculture_suitability",
    "povo_silencioso_distance_km",
    "povo_silencioso_exclusion_factor",
]


def _smootherstep(x: np.ndarray) -> np.ndarray:
    """Ken Perlin's smootherstep, x already clipped to [0, 1]."""
    return x * x * x * (x * (x * 6 - 15) + 10)


def _pad_to_factor(a: np.ndarray, factor: int) -> np.ndarray:
    ny, nx = a.shape
    py = (-ny) % factor
    px = (-nx) % factor
    if py or px:
        a = np.pad(a, ((0, py), (0, px)), mode="edge")
    return a


def block_max(a: np.ndarray, factor: int) -> np.ndarray:
    """Block-maximum downsample -- see module docstring for why this, not
    block_mean, is used for slope and stream presence."""
    if factor == 1:
        return a.astype(np.float64, copy=True)
    a = _pad_to_factor(a, factor)
    ny, nx = a.shape
    return a.reshape(ny // factor, factor, nx // factor, factor).max(axis=(1, 3))


def block_any(a: np.ndarray, factor: int) -> np.ndarray:
    """Block-OR downsample for boolean presence masks (e.g. stream_mask)."""
    if factor == 1:
        return a.astype(bool, copy=True)
    a = _pad_to_factor(a.astype(bool), factor)
    ny, nx = a.shape
    return a.reshape(ny // factor, factor, nx // factor, factor).any(axis=(1, 3))


def block_mean(a: np.ndarray, factor: int) -> np.ndarray:
    """Same as src/climate/grid.py's block_mean -- reproduced here so this
    module has no import-order dependency on src/climate. Used for
    contributing_area_km2, where a total-catchment-area quantity is
    reasonably approximated by simple resampling."""
    if factor == 1:
        return a.astype(np.float64, copy=True)
    a = _pad_to_factor(a, factor)
    ny, nx = a.shape
    return a.reshape(ny // factor, factor, nx // factor, factor).mean(axis=(1, 3))


def compute_slope_pct(dem: np.ndarray, cellsize_m: float) -> np.ndarray:
    """Slope magnitude as percent (rise/run * 100) via central differences."""
    gy, gx = np.gradient(dem, cellsize_m, cellsize_m)
    return np.sqrt(gx**2 + gy**2) * 100.0


def slope_suitability(
    slope_pct: np.ndarray, gentle_pct: float = 5.0, hard_limit_pct: float = 30.0
) -> np.ndarray:
    """0-1 score, 1.0 at/below `gentle_pct`, smootherstep decay to 0.0 at
    `hard_limit_pct`.

    HONEST CAVEAT: gentle_pct/hard_limit_pct are placeholder generic land-use
    -planning buildability figures, not calibrated against this project's own
    data the way e.g. Tappa 3's t50_snow_rain_c is. They are also almost
    certainly too permissive for the largest Circulos (dezenas de milhares)
    and too strict for the smallest (which the scenario explicitly places in
    steep Terrapedra-style terrain) -- the downstream settlement-size filter
    (not yet built) is meant to tighten/loosen the *effective* requirement per
    size class; this function only supplies the generic core layer.
    """
    x = np.clip((hard_limit_pct - slope_pct) / (hard_limit_pct - gentle_pct), 0.0, 1.0)
    return _smootherstep(x)


def distance_to_stream_km(
    stream_mask: np.ndarray, cellsize_m: tuple[float, float], factor: int = 1
) -> np.ndarray:
    """Euclidean distance (km) to the nearest stream cell. `stream_mask` may
    already be at the target resolution (factor=1) or native 30 m needing a
    block_any downsample first (factor=4)."""
    mask = block_any(stream_mask, factor) if factor > 1 else stream_mask.astype(bool)
    cs_y, cs_x = cellsize_m
    dist = ndimage.distance_transform_edt(~mask, sampling=(cs_y, cs_x))
    return dist / 1000.0


def water_suitability(
    dist_to_stream_km: np.ndarray, gentle_km: float = 0.5, hard_limit_km: float = 5.0
) -> np.ndarray:
    """0-1 score, 1.0 at/below `gentle_km`, smootherstep decay to 0.0 at
    `hard_limit_km` -- same shape family as slope_suitability, mirrored
    (closer is better here, instead of gentler).

    HONEST CAVEAT, and an important one for how much this layer actually
    matters: checked directly against this world's own hydrology, the
    stream network is dense enough that distance-to-stream barely
    discriminates across most of the land. 82.5% of land is already within
    0.5 km of a stream, 97.7% within 1 km, 98.9% within 2 km -- only the
    most remote 0.4% of land sits beyond 5 km (max 13.46 km). With
    gentle_km/hard_limit_km set from those percentiles (0.5 km ~ the point
    where "trivially close" stops being nearly everyone; 5 km ~ p99.5, so
    only the genuine outlier tail gets meaningfully penalised), this layer
    will read as ~flat (near 1.0) over the large majority of land and only
    pull scores down at the margin -- for the rare gentle/sunny/otherwise-
    good site that happens to sit unusually far from any stream. That is
    the correct behaviour for THIS world's hydrology (dense drainage, not a
    scarce-water setting), not a sign the function is miscalibrated; do not
    expect this layer to be a strong discriminator in the final composite.
    """
    x = np.clip((hard_limit_km - dist_to_stream_km) / (hard_limit_km - gentle_km), 0.0, 1.0)
    return _smootherstep(x)


def topographic_wetness_index(
    contributing_area_km2: np.ndarray,
    slope_pct: np.ndarray,
    cellsize_m: float,
    eps_slope_pct: float = 0.1,
) -> np.ndarray:
    """TWI = ln(specific catchment area / tan(slope)).

    `contributing_area_km2` is TOTAL catchment area (D8 accumulation), not
    the per-unit-contour-width specific catchment area TWI formally wants --
    dividing by `cellsize_m` approximates that, standard practice when only
    total accumulation is available (no multi-flow-direction routing with an
    explicit contour width in this pipeline, see Tappa 4). `slope_pct` is
    floored at `eps_slope_pct` to avoid a divide-by-zero / -inf on perfectly
    flat cells -- HONEST CAVEAT: this is a geomorphic proxy standing in for
    agricultural suitability with NO real soil/pedology data anywhere in this
    project (no texture, depth, drainage class, or pH); it ranks wet, low-
    relief, well-drained-looking land higher, nothing more.
    """
    slope_rad = np.arctan(np.maximum(slope_pct, eps_slope_pct) / 100.0)
    specific_area_m = (contributing_area_km2 * 1e6) / cellsize_m
    with np.errstate(divide="ignore"):
        return np.log(specific_area_m / np.tan(slope_rad))


def agriculture_suitability(twi: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """Normalize TWI to 0-1 via the 2nd/98th land percentiles (robust to the
    handful of extreme values right at channel heads)."""
    land_finite = land_mask & np.isfinite(twi)
    p2, p98 = np.percentile(twi[land_finite], [2, 98])
    suit = np.clip((twi - p2) / (p98 - p2), 0.0, 1.0)
    return np.where(land_mask, np.nan_to_num(suit, nan=0.0), np.nan)


def povo_silencioso_distance_km(
    land_mask: np.ndarray, archipelago_labels: list[int], cellsize_m: tuple[float, float]
) -> np.ndarray:
    """Distance (km) to the nearest cell of the Povo Silencioso's NE
    archipelago. `archipelago_labels` are the specific connected-component
    label ids for that cluster (identified once, by hand, against
    `land_mask` with `scipy.ndimage.label(..., structure=np.ones((3,3)))`) --
    NOT re-derived by a geometric heuristic each run, since a bounding-box
    rule over-collects unrelated tiny islets elsewhere in the domain.
    """
    labeled, _ = ndimage.label(land_mask, structure=np.ones((3, 3)))
    target = np.isin(labeled, archipelago_labels)
    cs_y, cs_x = cellsize_m
    dist = ndimage.distance_transform_edt(~target, sampling=(cs_y, cs_x))
    return dist / 1000.0


def povo_silencioso_exclusion_factor(
    dist_to_povo_silencioso_km: np.ndarray,
    hard_buffer_km: float = 5.0,
    soft_buffer_km: float = 15.0,
) -> np.ndarray:
    """0-1 MULTIPLIER (not an additive nucleo score): 0.0 at/within
    `hard_buffer_km` of the Povo Silencioso archipelago, smootherstep ramp
    up to 1.0 (no penalty) by `soft_buffer_km` -- the mirror shape of
    water_suitability (near is bad here, instead of good).

    HONEST CAVEAT: this only affects a small, localised area by
    construction -- checked directly, only 1.02% of non-archipelago land
    (~100 km2) falls within 10 km of the archipelago at all (land as close
    as 0.24 km exists, so the buffer is not vacuous, but it is inherently a
    corner-of-the-map effect, not a domain-wide one). hard_buffer_km/
    soft_buffer_km are a placeholder "respect the territory" buffer size,
    not derived from any treaty/lore distance specified so far -- revisit
    if the scenario ever specifies one.
    """
    x = np.clip(
        (dist_to_povo_silencioso_km - hard_buffer_km) / (soft_buffer_km - hard_buffer_km),
        0.0,
        1.0,
    )
    return _smootherstep(x)
