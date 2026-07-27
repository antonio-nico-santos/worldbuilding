"""
Self-contained 2D Simplex noise + ridged fractal Brownian motion (fBm).

Why self-contained rather than the `opensimplex` package listed in
requirements.txt: this generation session runs in a sandbox whose network
egress does not reach PyPI (pip install opensimplex/rasterio/geopandas all
fail with host_not_allowed). Rather than deliver untested code, this module
reimplements classic 2D Simplex noise (Gustavson's algorithm) directly in
numpy, fully vectorized (no per-pixel Python loops), so the whole pipeline
can run and be verified end-to-end in this sandbox. It has no external
dependency beyond numpy. If you'd rather depend on the `opensimplex`
package when running locally (where PyPI is reachable), the two should be
visually similar (same noise family) but will NOT produce bit-identical
output for the same seed -- pick one and stick with it for reproducibility.
"""

import numpy as np

_F2 = 0.5 * (np.sqrt(3.0) - 1.0)
_G2 = (3.0 - np.sqrt(3.0)) / 6.0

# 8-direction gradient set (unit vectors at 45-degree increments). Classic
# simplex implementations often reuse the 12 edge-midpoints of a cube
# (grad3) projected to 2D; an 8-direction set is simpler and produces
# equivalent-quality 2D noise.
_GRADIENTS = np.array(
    [
        [1, 0], [-1, 0], [0, 1], [0, -1],
        [1, 1], [-1, 1], [1, -1], [-1, -1],
    ],
    dtype=np.float64,
)
_GRADIENTS[4:] /= np.sqrt(2.0)


class Simplex2D:
    """Vectorized 2D Simplex noise, seeded independently of global numpy state."""

    def __init__(self, seed: int = 0):
        rng = np.random.RandomState(seed)
        perm = rng.permutation(256).astype(np.int64)
        self.perm = np.concatenate([perm, perm])  # length 512, avoids wraparound
        self.perm_mod8 = self.perm % 8

    def noise2(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Evaluate noise at arbitrary (x, y) float coordinates (any shape).
        Output is approximately in [-1, 1]."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)

        s = (x + y) * _F2
        i = np.floor(x + s)
        j = np.floor(y + s)
        t = (i + j) * _G2
        X0 = i - t
        Y0 = j - t
        x0 = x - X0
        y0 = y - Y0

        i1 = (x0 > y0).astype(np.int64)
        j1 = 1 - i1

        x1 = x0 - i1 + _G2
        y1 = y0 - j1 + _G2
        x2 = x0 - 1.0 + 2.0 * _G2
        y2 = y0 - 1.0 + 2.0 * _G2

        ii = i.astype(np.int64) & 255
        jj = j.astype(np.int64) & 255

        perm = self.perm
        pmod8 = self.perm_mod8

        gi0 = pmod8[ii + perm[jj]]
        gi1 = pmod8[ii + i1 + perm[jj + j1]]
        gi2 = pmod8[ii + 1 + perm[jj + 1]]

        def corner(gi, cx, cy):
            t = 0.5 - cx * cx - cy * cy
            t = np.maximum(t, 0.0)
            t4 = t * t * t * t
            grad = _GRADIENTS[gi]
            dot = grad[..., 0] * cx + grad[..., 1] * cy
            return t4 * dot

        n0 = corner(gi0, x0, y0)
        n1 = corner(gi1, x1, y1)
        n2 = corner(gi2, x2, y2)

        return 70.0 * (n0 + n1 + n2)


def ridged_fbm(
    noise: Simplex2D,
    x: np.ndarray,
    y: np.ndarray,
    octaves: int,
    base_freq: float,
    lacunarity: float = 2.0,
    persistence: float = 0.5,
    plain_octaves: int = 0,
    persistence_fine: float = None,
    crossover_wavelength_m: float = None,
    crossover_width_factor: float = 2.5,
) -> np.ndarray:
    """Signed ridged/turbulence fBm.

    `persistence_fine` / `crossover_wavelength_m` (optional, both None by
    default = old single-persistence behavior): real terrain's power
    spectrum is NOT one clean power law -- measured directly from a real
    SRTM tile (Piemonte/western Alps, 90m, two independent clean 1024x1024
    patches), beta (power ~ freq^-beta) is ~2.5 at wavelengths above ~5km
    (rougher -- tectonic/regional-scale relief) but climbs to ~4.5 below
    ~1.2km (much smoother -- almost certainly hillslope diffusion, a real
    physical process that erodes/smooths short-wavelength terrain, acting
    on top of the rougher large-scale structure). A single `persistence`
    implies one Hurst exponent everywhere and can only match one of these
    two regimes at a time -- tuned for the large-scale roughness (like the
    2026-07-27 fix), the fine-scale octaves come out rougher/grainier than
    real terrain actually is. `persistence_fine` lets the fine octaves
    decay faster (steeper, smoother) than the dominant/large-scale ones,
    blended smoothly (smootherstep in log-wavelength space, width
    `crossover_width_factor` as a multiplicative factor around
    `crossover_wavelength_m`) rather than an abrupt per-octave switch,
    which would itself add a spectral discontinuity of the same kind this
    project has already had to fix twice (see domain_warp/plain_octaves
    docstrings below).

    Per octave: fold the raw simplex value with `1 - 2*abs(n)`, which is
    +1 exactly where the underlying noise crosses zero (producing sharp
    ridge creases along a network of zero-contours) and -1 at the noise's
    extremes (smooth, low valley floors). This is the signed variant of
    the standard abs()-fold ridge technique -- the plain (unsigned)
    `1 - abs(n)` fold is always >= 0 and has a strongly positive mean,
    which would bias every far-field cell upward and make it impossible
    to get clean ocean; folding to [-1, 1] keeps the result roughly
    zero-mean so it can push terrain both above and below the sea-level
    offset, while keeping the ridged/turbulence character (sharp peaks,
    smooth troughs) that was chosen for this project's background noise.

    `plain_octaves`: the first N octaves (lowest frequency, largest
    amplitude -- the ones that decide WHERE land exists at all) skip the
    ridge fold and use the raw noise value instead. This matters because
    zero-crossings of a smooth field are, by construction, a CONNECTED
    network (like river deltas or cracked mud), never isolated blobs --
    folding at every octave, including the dominant one, is what made
    "scattered islands" come out as one continuous winding filament
    instead. Reserving the fold for higher (finer-detail, smaller-
    amplitude) octaves keeps the jagged ridged texture as fine roughness
    while letting the large-scale island/no-island placement come from
    plain noise, whose excursion sets ARE compact, separate blobs.

    Returns an array in approximately [-1, 1], NOT further scaled to
    metres -- multiply by a chosen amplitude afterwards.
    """
    total = np.zeros_like(x, dtype=np.float64)
    freq = base_freq
    amp = 1.0
    norm = 0.0
    for octave_i in range(octaves):
        n = noise.noise2(x * freq, y * freq)
        n = np.clip(n, -1.0, 1.0)
        if octave_i < plain_octaves:
            contribution = n  # already ~zero-mean, no folding
        else:
            contribution = 1.0 - 2.0 * np.abs(n)
        total += contribution * amp
        norm += amp
        if persistence_fine is None:
            octave_persistence = persistence
        else:
            wavelength_m = 1.0 / freq
            half_width_log = np.log(crossover_width_factor) / 2.0
            t = (np.log(crossover_wavelength_m) - np.log(wavelength_m) + half_width_log) / (2.0 * half_width_log)
            t = np.clip(t, 0.0, 1.0)
            t = t * t * t * (t * (t * 6 - 15) + 10)  # smootherstep -- t=0 large-scale, t=1 fine-scale
            octave_persistence = persistence * (1 - t) + persistence_fine * t
        freq *= lacunarity
        amp *= octave_persistence
    return total / norm


def domain_warp(x: np.ndarray, y: np.ndarray, seed: int, wavelength_m: float, amplitude_m: float, octaves: int = 3, lacunarity: float = 2.2, persistence: float = 0.5):
    """Displace (x, y) by a smooth pseudo-random vector field before it's
    used for anything else (ridge distance, zone membership, the terrain
    noise itself).

    Why this exists: distance-to-a-curve (what drives the ridge decay and
    zone boundaries) produces perfectly smooth, parallel offset contours --
    mathematically a "buffer" of the original line/polygon. Real
    coastlines are never that regular (rivers, headlands, embayments,
    differential erosion all roughen them at multiple scales). The
    standard fix (used throughout procedural terrain generation, e.g.
    Inigo Quilez's fBm domain warping) is to NOT query the distance field
    at the true (x, y), but at (x, y) shifted by another noise field --
    so straight buffer contours become wiggly, organic ones. Two
    independently-seeded noise channels drive the x and y displacement so
    the warp isn't just a uniform rotation/skew.

    octaves > 1 sums the displacement at several frequencies (like fBm),
    NOT a single noise2 call. A single frequency has one characteristic
    bump spacing, which is exactly what produced the too-regular
    "scalloped" coastline (every lobe roughly the same size, roughly
    evenly spaced) -- the same underlying issue as the round islands
    before the secondary noise warp. Multiple octaves mix several
    spacings so peninsula/bay size varies rather than repeating.
    """
    warp_noise_x = Simplex2D(seed=seed + 9001)
    warp_noise_y = Simplex2D(seed=seed + 9002)
    freq = 1.0 / wavelength_m
    dx = np.zeros_like(x, dtype=np.float64)
    dy = np.zeros_like(y, dtype=np.float64)
    amp = 1.0
    norm = 0.0
    for _ in range(octaves):
        dx += warp_noise_x.noise2(x * freq, y * freq) * amp
        dy += warp_noise_y.noise2(x * freq, y * freq) * amp
        norm += amp
        freq *= lacunarity
        amp *= persistence
    dx = dx / norm * amplitude_m
    dy = dy / norm * amplitude_m
    return x + dx, y + dy
