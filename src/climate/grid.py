"""
Grid utilities for Tappa 2 — coarsening the 30 m DEM to the climate working
resolution, and the inverse Lambert Conformal Conic needed for the latitude
baseline of the temperature model.

Why a separate working resolution (Tappa 2 decision doc, S2): twelve monthly
temperature + precipitation layers at the DEM's native 30 m would be ~2.2 GB.
Climate fields are intrinsically smooth — the Smith & Barstad model's own
advection length (U * tau, ~15 km) is three orders of magnitude coarser than
30 m — so nothing is lost by solving at 120 m. Temperature is a closed-form
function of the DEM and can be re-evaluated at 30 m for hero renders at any
time (see temperature.temperature_month).
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "block_mean",
    "coarsen",
    "lcc_inverse_lat",
    "latitude_grid",
]


def block_mean(a: np.ndarray, factor: int) -> np.ndarray:
    """Area-average `a` by `factor` in both axes.

    The array is edge-replicated up to a multiple of `factor` first, so no
    cells are silently dropped. For the project's 5334 x 4334 grid at
    factor 4 this replicates 2 rows and 2 columns (60 m) at the S and E
    edges — negligible, and documented rather than trimmed so the coarse
    grid still starts exactly at the ROI's (xmin, ymax) origin.
    """
    if factor == 1:
        return a.astype(np.float64, copy=True)
    ny, nx = a.shape
    py = (-ny) % factor
    px = (-nx) % factor
    if py or px:
        a = np.pad(a, ((0, py), (0, px)), mode="edge")
    ny, nx = a.shape
    return a.reshape(ny // factor, factor, nx // factor, factor).mean(axis=(1, 3))


def coarsen(dem30: np.ndarray, factor: int = 4):
    """Coarsen the eroded 30 m DEM to the climate working grid.

    Returns ``(surface, uplift_h, land_fraction)``:

    * ``surface``   — area-mean elevation including bathymetry. This is the
      land surface the temperature model applies its lapse rate to.
    * ``uplift_h``  — area-mean of ``max(dem, 0)``. This is what the
      orographic model sees: air flows over the sea *surface*, not over the
      sea floor, so bathymetry must be clipped BEFORE averaging. Averaging
      first and clipping after would let a deep trench next to a headland
      cancel the headland's uplift.
    * ``land_fraction`` — fraction of each coarse cell above 0 m, so callers
      can pick their own land/sea threshold rather than inheriting one.
    """
    surface = block_mean(dem30, factor)
    uplift_h = block_mean(np.maximum(dem30, 0.0), factor)
    land_fraction = block_mean((dem30 > 0).astype(np.float32), factor)
    return surface, uplift_h, land_fraction


# --- inverse Lambert Conformal Conic (Snyder 1987, eqs. 14-x / 15-x) --------
#
# The latitude term in the temperature model is worth only ~1 C across the
# whole 160 km domain, so a spherical approximation would be defensible. The
# ellipsoidal inverse is implemented anyway because this is a GIS portfolio
# piece and getting the projection maths right is part of the point; it costs
# about twenty lines and runs once per grid.

_A_WGS84 = 6378137.0
_F_WGS84 = 1.0 / 298.257223563


def _m(phi: float, e: float) -> float:
    return math.cos(phi) / math.sqrt(1.0 - e * e * math.sin(phi) ** 2)


def _t(phi: float, e: float) -> float:
    s = math.sin(phi)
    return math.tan(math.pi / 4.0 - phi / 2.0) / ((1.0 - e * s) / (1.0 + e * s)) ** (e / 2.0)


def lcc_inverse_lat(
    y: np.ndarray,
    lat_1: float,
    lat_2: float,
    lat_0: float,
    a: float = _A_WGS84,
    f: float = _F_WGS84,
    x: np.ndarray | float = 0.0,
) -> np.ndarray:
    """Latitude (degrees) of northing `y` (metres) in a two-parallel LCC.

    `x` is accepted because the parallels of an LCC are arcs, not straight
    lines: at the domain's east/west edges (65 km off the central meridian)
    a given `y` is very slightly further from the pole than on-axis. The
    effect here is under 0.01 deg, but taking `x` makes that explicit rather
    than hidden.
    """
    e = math.sqrt(2 * f - f * f)
    p1, p2, p0 = math.radians(lat_1), math.radians(lat_2), math.radians(lat_0)
    m1, m2 = _m(p1, e), _m(p2, e)
    t1, t2, t0 = _t(p1, e), _t(p2, e), _t(p0, e)
    n = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
    F = m1 / (n * t1**n)
    rho0 = a * F * t0**n

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = np.hypot(x, rho0 - y) * np.sign(n)
    t = (rho / (a * F)) ** (1.0 / n)

    # Snyder eq. 3-5: iterate to invert t -> phi (converges in ~4 passes)
    phi = math.pi / 2.0 - 2.0 * np.arctan(t)
    for _ in range(8):
        s = np.sin(phi)
        phi = math.pi / 2.0 - 2.0 * np.arctan(t * ((1.0 - e * s) / (1.0 + e * s)) ** (e / 2.0))
    return np.degrees(phi)


def latitude_grid(ny: int, nx: int, ymax: float, xmin: float, res: float, **lcc) -> np.ndarray:
    """Per-cell latitude (degrees) for a north-up grid, cell-centre sampled."""
    yc = ymax - (np.arange(ny) + 0.5) * res
    xc = xmin + (np.arange(nx) + 0.5) * res
    return lcc_inverse_lat(yc[:, None], x=xc[None, :], **lcc)
