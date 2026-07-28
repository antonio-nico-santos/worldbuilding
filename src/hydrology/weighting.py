"""
Upsample Tappa 2's 120 m annual precipitation field onto the 30 m DEM grid,
for precipitation-weighted flow accumulation (see docs/decisions/
04_tappa4_hydrology.md).

Bilinear (order=1), not nearest/block-repeat: Tappa 2 S2 already established
precipitation has no signal below ~1 km (the Smith & Barstad advection length
is ~24 km), so nothing physically meaningful is added by interpolating rather
than repeating -- but repeating would stamp visible 120 m block-boundary
steps onto the accumulated-flow field, a purely numerical artifact bilinear
interpolation avoids for free.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import zoom


def upsample_precip_to_dem(precip_annual_mm: np.ndarray, dem_shape: tuple[int, int]) -> np.ndarray:
    dem_ny, dem_nx = dem_shape
    p_ny, p_nx = precip_annual_mm.shape
    zy, zx = dem_ny / p_ny, dem_nx / p_nx
    up = zoom(precip_annual_mm.astype(np.float64), (zy, zx), order=1, mode="nearest")
    # zoom's output size can be off by a cell or two from a non-integer
    # ratio; pad (edge-replicate, same convention as grid.block_mean) or
    # crop to match the DEM exactly.
    out = np.empty(dem_shape, dtype=np.float64)
    uy, ux = up.shape
    y = min(uy, dem_ny)
    x = min(ux, dem_nx)
    out[:y, :x] = up[:y, :x]
    if uy < dem_ny:
        out[uy:, :x] = up[-1:, :x]
    if ux < dem_nx:
        out[:y, ux:] = up[:y, -1:]
    if uy < dem_ny and ux < dem_nx:
        out[uy:, ux:] = up[-1, -1]
    return out
