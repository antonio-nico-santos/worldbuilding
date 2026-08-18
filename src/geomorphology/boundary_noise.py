"""
Tappa 8 -- boundary domain-warp noise ("apply noise at the end", Nico's
third request on lithology geometry).

IMPORTANT HONESTY NOTE: Tappa 1 already has its own domain-warp noise
(`noise_warp_wavelength_m` / `noise_warp_amplitude_m` / `noise_warp_octaves`
/ `noise_warp_lacunarity` / `noise_warp_persistence` in parameters.yml,
implemented in `terrain.generate`/`terrain.erosion`). Those modules were
NEVER staged into this session -- only the specific files this stage needed
were uploaded (see the uploads listing). This module is an INDEPENDENT,
comparable value-noise implementation built for this stage's own need, not
a reuse of Tappa 1's actual code. If reusing the literal Tappa 1 noise
function matters, that has to happen on the actual device repo, not here.

Method: band-limited value noise -- generate a coarse grid of iid standard
normal values (grid spacing = one noise wavelength), bicubic-upsample it to
the full raster resolution, and sum multiple octaves at halving wavelength /
decaying amplitude (persistence) for a fractal-ish look. This is NOT
Perlin/simplex noise (no gradient-vector construction, no lattice
gradients) -- it's the simpler "smoothed random field" family, which is
sufficient for a boundary-roughening displacement field where the exact
noise flavor is not the object being verified against a citation, unlike
(say) Tappa 1's actual terrain amplitude.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import zoom


def value_noise_field(ny: int, nx: int, cellsize_m: float, base_wavelength_m: float,
                       octaves: int = 3, persistence: float = 0.5, lacunarity: float = 2.0,
                       seed: int = 0) -> np.ndarray:
    """Returns an (ny, nx) float32 field, zero mean, roughly unit std
    (exact std depends on octave count/persistence, NOT re-normalized after
    summing -- callers scale by their own target amplitude)."""
    rng = np.random.default_rng(seed)
    field = np.zeros((ny, nx), dtype=np.float64)
    wavelength_m = base_wavelength_m
    amp = 1.0
    for _o in range(octaves):
        coarse_ny = max(2, int(np.ceil((ny * cellsize_m) / wavelength_m)) + 2)
        coarse_nx = max(2, int(np.ceil((nx * cellsize_m) / wavelength_m)) + 2)
        coarse = rng.standard_normal(size=(coarse_ny, coarse_nx))
        zy, zx = ny / coarse.shape[0], nx / coarse.shape[1]
        up = zoom(coarse, (zy, zx), order=3, mode="nearest")[:ny, :nx]
        field += amp * up
        amp *= persistence
        wavelength_m /= lacunarity
    return field.astype(np.float32)


def boundary_warp_field(ny: int, nx: int, cellsize_m: float, amplitude_m: float,
                         base_wavelength_m: float, octaves: int = 3,
                         persistence: float = 0.5, lacunarity: float = 2.0,
                         seed_x: int = 101, seed_y: int = 202):
    """Two independent noise fields (dx, dy), each scaled to
    `amplitude_m` std, for use as `classify()`'s warp_dx/warp_dy -- a
    displacement applied to the query points before measuring distance to
    the real crest geometry, roughening the resulting class boundary
    without needing any changes to the crest-extraction method itself."""
    fx = value_noise_field(ny, nx, cellsize_m, base_wavelength_m, octaves, persistence, lacunarity, seed=seed_x)
    fy = value_noise_field(ny, nx, cellsize_m, base_wavelength_m, octaves, persistence, lacunarity, seed=seed_y)
    fx = fx / (fx.std() + 1e-9) * amplitude_m
    fy = fy / (fy.std() + 1e-9) * amplitude_m
    return fx, fy
