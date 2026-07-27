"""
Tappa 1: assemble the DEM from the hand-authored skeleton (ridges + zones)
plus background noise, calibrated so a coastline emerges as the 0 m
contour. See docs/decisions/01_tappa1_terrain.md for the full rationale
behind each formula/constant here.

Processes the grid in ROW CHUNKS rather than one 23-million-point array.
Discovered the hard way: at full 30m resolution the whole-domain version
of this got OOM-killed partway through (each noise2() call -- and there
are several per octave, times several octaves, times two warp layers --
allocates a double handful of 23M-float64 temporaries; they don't all
coexist for long, but enough do at once to exceed this sandbox's ~7GB).
Chunking bounds peak memory to one chunk's worth regardless of total
domain size, which is the standard approach for large-raster geoprocessing.
"""

import time

import numpy as np
from scipy.ndimage import map_coordinates

from .noise import Simplex2D, ridged_fbm, domain_warp
from .skeleton import load_geojson, build_ridge_fields, build_zone_fields


def _estimate_noise_normalization(seed, noise_octaves, noise_base_wavelength_m, lacunarity, persistence, noise_plain_octaves, xmin, xmax, ymin, ymax,
                                   persistence_fine=None, crossover_wavelength_m=None, crossover_width_factor=2.5, sample_n=400):
    """Ridged fBm's fold is NOT zero-mean by construction (see noise.py) --
    every chunk needs to subtract the SAME mean/std to normalize
    consistently, or each chunk would center on its own local statistics
    and produce a visible seam at every chunk boundary. Estimate once from
    a modest, domain-spanning sample (independent of chunk boundaries) and
    reuse everywhere."""
    rng = np.random.RandomState(seed + 777)
    sx = rng.uniform(xmin, xmax, size=sample_n * sample_n)
    sy = rng.uniform(ymin, ymax, size=sample_n * sample_n)
    noise = Simplex2D(seed=seed)
    base_freq = 1.0 / noise_base_wavelength_m
    raw = ridged_fbm(noise, sx, sy, noise_octaves, base_freq, lacunarity, persistence, plain_octaves=noise_plain_octaves,
                      persistence_fine=persistence_fine, crossover_wavelength_m=crossover_wavelength_m, crossover_width_factor=crossover_width_factor)
    return raw.mean(), raw.std()


def generate_dem(
    xmin, xmax, ymin, ymax, resolution_m,
    ridges_path, zones_path,
    seed, noise_octaves, noise_base_wavelength_m,
    noise_amplitude_m, sea_level_offset_m,
    shelf_multipliers=None, default_shelf_multiplier=3.0,
    zone_base_lift_m=None, default_zone_base_lift_m=0.0,
    warp_wavelength_m=None, warp_amplitude_m=0.0,
    noise_warp_wavelength_m=None, noise_warp_amplitude_m=0.0,
    noise_plain_octaves=0,
    lacunarity=2.0, persistence=0.5,
    noise_persistence_fine=None, noise_crossover_wavelength_m=None, noise_crossover_width_factor=2.5,
    warp_octaves=3, warp_lacunarity=2.2, warp_persistence=0.5,
    noise_warp_octaves=3, noise_warp_lacunarity=2.2, noise_warp_persistence=0.5,
    real_detail_path=None, real_detail_xmin=None, real_detail_ymax=None, real_detail_cellsize_m=None,
    real_detail_fine_supplement_weight=0.0, real_detail_fine_min_wavelength_m=60.0,
    row_chunk_size=250,
    verbose=True,
):
    def log(msg):
        if verbose:
            print(f"[{time.time()-t_start:6.1f}s] {msg}")

    t_start = time.time()

    x = np.arange(xmin, xmax, resolution_m)
    y = np.arange(ymin, ymax, resolution_m)
    nx, ny = len(x), len(y)
    log(f"grid ({ny}, {nx}) ({nx*ny:,} cells) -- processing in row chunks of {row_chunk_size}")

    ridge_features = load_geojson(ridges_path)
    zone_features = load_geojson(zones_path)
    ridge_fields = build_ridge_fields(ridge_features, shelf_multipliers=shelf_multipliers, default_shelf_multiplier=default_shelf_multiplier)
    zone_fields = build_zone_fields(zone_features, base_lift_m=zone_base_lift_m, default_base_lift_m=default_zone_base_lift_m)
    log(f"loaded {len(ridge_fields)} ridges, {len(zone_fields)} zones")

    use_real_detail = real_detail_path is not None
    if use_real_detail:
        # C2: grafted real-terrain detail instead of synthetic noise. The array is
        # already detrended (its own large-scale "which massif is where" placement
        # removed -- see data/real_detail/piemonte_detail_meta.json) so only its
        # anisotropic ridge/valley TEXTURE remains; our own ridges/zones still own
        # 100% of where land/mountains/coastline actually are. Its exact mean/std
        # are used directly (no sampling estimate needed -- we HAVE the full array,
        # unlike the synthetic case where a sample was the only practical option).
        real_detail = np.load(real_detail_path).astype(np.float64)
        real_mean, real_std = real_detail.mean(), real_detail.std()
        real_ny, real_nx = real_detail.shape
        log(f"real-detail texture loaded: {real_detail.shape}, mean={real_mean:.2f} std={real_std:.2f}")

        if real_detail_fine_supplement_weight:
            # same seam-avoidance reasoning as the synthetic path: this MUST be
            # estimated once globally, not per chunk, or each chunk normalizes to
            # its own local statistics and a visible seam appears at every
            # chunk boundary (bit us once already during Tappa 1 -- see
            # _estimate_noise_normalization's docstring).
            # starts at real data's own Nyquist (coarsest gap it can't fill) and
            # gets finer each octave via lacunarity, down to real_detail_fine_min_wavelength_m
            fine_start_wavelength_m = real_detail_cellsize_m * 2
            fine_octaves = max(1, int(np.ceil(np.log(fine_start_wavelength_m / real_detail_fine_min_wavelength_m) / np.log(lacunarity))))
            fine_mean, fine_std = _estimate_noise_normalization(
                seed, fine_octaves, fine_start_wavelength_m, lacunarity, 0.5, 0,
                xmin, xmax, ymin, ymax,
            )
            log(f"fine-detail supplement: {fine_octaves} octaves down to {real_detail_fine_min_wavelength_m}m, mean={fine_mean:.3f} std={fine_std:.3f}")
    else:
        raw_mean, raw_std = _estimate_noise_normalization(
            seed, noise_octaves, noise_base_wavelength_m, lacunarity, persistence, noise_plain_octaves,
            xmin, xmax, ymin, ymax,
            persistence_fine=noise_persistence_fine, crossover_wavelength_m=noise_crossover_wavelength_m,
            crossover_width_factor=noise_crossover_width_factor,
        )
        log(f"noise normalization estimated once: mean={raw_mean:.3f} std={raw_std:.3f} (reused for every chunk)")

    elevation = np.empty((ny, nx), dtype=np.float32)
    noise = Simplex2D(seed=seed)
    base_freq = 1.0 / noise_base_wavelength_m

    n_chunks = int(np.ceil(ny / row_chunk_size))
    for chunk_i in range(n_chunks):
        r0 = chunk_i * row_chunk_size
        r1 = min(r0 + row_chunk_size, ny)
        y_chunk = y[r0:r1]
        X, Y = np.meshgrid(x, y_chunk)  # row 0 of this chunk = its first y value
        chunk_shape = X.shape
        xy = np.column_stack([X.ravel(), Y.ravel()])

        if warp_amplitude_m:
            wx, wy = domain_warp(xy[:, 0], xy[:, 1], seed, warp_wavelength_m, warp_amplitude_m,
                                  octaves=warp_octaves, lacunarity=warp_lacunarity, persistence=warp_persistence)
            xy_query = np.column_stack([wx, wy])
        else:
            xy_query = xy

        structure = np.zeros(xy.shape[0], dtype=np.float64)
        for rf in ridge_fields:
            structure = np.maximum(structure, rf.contribution(xy_query))

        amplitude_multiplier = np.ones(xy.shape[0], dtype=np.float64)
        for zf in zone_fields:
            w = zf.blend_weight(xy_query)
            if zf.feature_type == "plateau":
                structure = structure * (1 - w) + zf.target_elevation_m * w
            elif zf.feature_type == "amplitude_zone" and zf.base_lift_m:
                structure = structure + w * zf.base_lift_m
            amplitude_multiplier = amplitude_multiplier * (1 - w) + zf.amplitude_scale * w

        xy_noise = xy_query
        if noise_warp_amplitude_m:
            nwx, nwy = domain_warp(xy_query[:, 0], xy_query[:, 1], seed + 500, noise_warp_wavelength_m, noise_warp_amplitude_m,
                                    octaves=noise_warp_octaves, lacunarity=noise_warp_lacunarity, persistence=noise_warp_persistence)
            xy_noise = np.column_stack([nwx, nwy])

        if use_real_detail:
            # map world (x,y) -> fractional (row,col) into the real-detail array.
            # row 0 = north (real_detail_ymax), col 0 = west (real_detail_xmin) --
            # same convention as this function's own final flipped output.
            real_col = (xy_noise[:, 0] - real_detail_xmin) / real_detail_cellsize_m
            real_row = (real_detail_ymax - xy_noise[:, 1]) / real_detail_cellsize_m
            real_col = np.clip(real_col, 0, real_nx - 1)
            real_row = np.clip(real_row, 0, real_ny - 1)
            sampled = map_coordinates(real_detail, [real_row, real_col], order=1, mode="nearest")
            noise_z = (sampled - real_mean) / real_std
            if real_detail_fine_supplement_weight:
                # real data has no information below its own native resolution --
                # this DEM's target resolution is finer, so upsampling alone would
                # leave the last mile (native_cellsize down to ~2 cells) smoother
                # than genuine terrain roughness at that scale. A SMALL amount of
                # synthetic fine noise, confined to wavelengths the real data can't
                # supply, closes that specific gap without displacing the real
                # signal anywhere it actually has coverage. fine_octaves/fine_mean/
                # fine_std were estimated ONCE above (same seam-avoidance reasoning
                # as everywhere else in this file) -- not recomputed per chunk.
                fine_raw = ridged_fbm(noise, xy_noise[:, 0], xy_noise[:, 1], fine_octaves, 1.0 / fine_start_wavelength_m, lacunarity, 0.5)
                fine_z = (fine_raw - fine_mean) / fine_std
                noise_z = noise_z + fine_z * real_detail_fine_supplement_weight
        else:
            raw = ridged_fbm(noise, xy_noise[:, 0], xy_noise[:, 1], noise_octaves, base_freq, lacunarity, persistence, plain_octaves=noise_plain_octaves,
                              persistence_fine=noise_persistence_fine, crossover_wavelength_m=noise_crossover_wavelength_m,
                              crossover_width_factor=noise_crossover_width_factor)
            noise_z = (raw - raw_mean) / raw_std

        chunk_elevation = structure + noise_z * noise_amplitude_m * amplitude_multiplier - sea_level_offset_m
        elevation[r0:r1, :] = chunk_elevation.reshape(chunk_shape).astype(np.float32)

        if verbose and (chunk_i % max(1, n_chunks // 10) == 0 or chunk_i == n_chunks - 1):
            log(f"  chunk {chunk_i+1}/{n_chunks} (rows {r0}-{r1}) done")

    land_frac = (elevation > 0).mean()
    log(f"elevation range {elevation.min():.0f}..{elevation.max():.0f} m, land fraction {land_frac:.1%}")

    # Flip so row 0 = north, matching ESRI ASCII Grid's top-to-bottom
    # convention (row 0 as built above = ymin = SOUTH).
    elevation = np.flipud(elevation)

    return elevation
