# Tappa 1 — Procedural terrain (DEM) generation

Status: **DEM finalized, v3**, at 30 m resolution, hydraulic erosion
applied. §§1–9 below document the original (v1) run; §10 documents a
second pass (v2) that replaced v1's background noise with a texture
grafted from a real DEM; §11 documents a third pass (v3) that fixed a
domain-warp folding bug found in v2 by closer inspection in QGIS.
**§11's parameters and outputs are the ones actually in use now** — §7's
and §10's `.bin`/`.npy` files are superseded and kept only as the
historical record of how this stage evolved. This document is the
decision record for this stage, per the project's established per-stage
workflow (see `README.md`). It is written to serve either as a log of
what was decided here, or as a reusable recipe for applying the same
skeleton+noise+erosion approach to a different domain.

> **Superseded in part.** The skeleton GeoJSONs were edited after this
> document was written — every `amplitude_zone` became a `plateau` (the fix
> §13 lists as planned for the North plains), the Spine peak was lowered,
> the South Branch falloff shortened, and an `Island` ridge added — and the
> pipeline was re-run with them. §1's skeleton table and §11c's output
> figures are therefore stale. See **`02_tappa2_climate.md` §0** for the
> skeleton and run figures actually on disk. The `base_lift_m` follow-up in
> §13 is now moot, since no `amplitude_zone` features remain.

## 1. Inputs from Tappa 0

- Domain: 130 km (E–W) x 160 km (N–S), `xmin=-65000, xmax=65000,
  ymin=-80000, ymax=80000` (metres, in the project's custom CRS).
- CRS: "Fictional World LCC" — Lambert Conformal Conic,
  `+proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 +x_0=0 +y_0=0
  +datum=WGS84 +units=m +no_defs` (from `config/parameters.yml`).
- Hand-authored terrain skeleton (`data/input/terrain_ridges.geojson`,
  `terrain_zones.geojson`):

  | ridge | peak_elevation_m | falloff_km |
  |---|---|---|
  | Spine | 3794 | 18 |
  | North branch (Big Brother) | 1650 | 13 |
  | West Branch (Little Brother) | 930 | 10 |
  | South Branch | 2250 | 23 |

  | zone | type | target_elevation_m | amplitude_scale | edge_transition_km |
  |---|---|---|---|---|
  | NW Plateau | plateau | 1200 | 0.2 | 5 |
  | Central plateau | plateau | 1900 | 0.3 | 10 |
  | North plains | amplitude_zone | — | 0.3 | 12 |

## 2. Why no external GIS/noise libraries

This session's sandbox has no PyPI egress (`pip install` fails with
`host_not_allowed` for `pypi.org`/`files.pythonhosted.org`). `opensimplex`,
`rasterio`, `geopandas`, `shapely` were all unavailable. Everything needed
was reimplemented on top of numpy/scipy/matplotlib only:

- 2D Simplex noise (Gustavson's algorithm), vectorized — `src/terrain/noise.py`.
- Ridge/zone distance fields via `scipy.spatial.cKDTree` (nearest-vertex
  distance on a densified polyline, not exact point-to-segment — error
  bound documented in `skeleton.py`, negligible at the spacing used).
- Point-in-polygon via `matplotlib.path.Path` instead of shapely.
- Raster export as ESRI ASCII Grid or ENVI raw binary instead of GeoTIFF
  (no rasterio/GDAL Python bindings available — see §7).

**In a normal Python environment with PyPI access**, `opensimplex`/
`shapely`/`rasterio` could be swapped in directly — the noise family is
visually equivalent but will NOT reproduce this exact DEM bit-for-bit for
the same seed. Whichever implementation is used should stay fixed for
reproducibility.

## 3. Elevation model

`elevation(x, y) = structure(x, y) + noise(x, y) * noise_amplitude_m *
amplitude_multiplier(x, y) - sea_level_offset_m`

### 3a. Ridge contribution (`structure`)

Each ridge is a Gaussian-style decay from its polyline, `falloff_km` set
so that it's literally the **half-max distance**:

```
base(d) = peak_elevation_m * exp(-ln2 * (d / falloff_m)^2)
```

Windowed to exactly zero beyond `shelf_km = falloff_km * shelf_multiplier`
via a smootherstep taper (else the Gaussian tail never truly reaches zero,
and a single global sea-level offset would be the only thing separating
land from ocean everywhere). `structure` at any point is the `max` across
all ridges (not a sum), so overlapping ridges don't stack additively.

Locked-in `shelf_multiplier` per ridge (tuned visually, see §8):

| ridge | shelf_multiplier | shelf reach |
|---|---|---|
| Spine | 1.6 | 28.8 km |
| North branch (Big Brother) | 1.3 | 16.9 km |
| West Branch (Little Brother) | 1.3 | 13.0 km |
| South Branch | 1.3 | 29.9 km |

### 3b. Zones

- **`plateau`**: hard elevation override, blended in across
  `edge_transition_km` via smootherstep: `structure = structure*(1-w) +
  target_elevation_m*w`. Deep inside the polygon, elevation IS
  `target_elevation_m`, regardless of what ridge decay was doing there.
- **`amplitude_zone`**: never touches elevation directly — only scales
  `amplitude_multiplier` (noise roughness). The underlying ridge/noise
  trend still has to reach that spot for it to read as land. See
  `docs/terrain_skeleton_attributes_ADDENDUM.md` for the full writeup of
  this distinction and the `base_lift_m` mechanism that patches the gap
  (below).
- **`base_lift_m`** (North plains: 300 m): additive-only elevation nudge
  for `amplitude_zone` features that sit beyond every ridge's shelf reach
  (North plains is 36–52 km from the nearest ridge — well past any shelf
  above). Blended by the same zone weight `w`:
  `structure = structure + w * base_lift_m`. **Not yet a real attribute**
  in `terrain_skeleton.gpkg`/GeoJSON — currently a Python-side dict in
  the generation call, keyed by zone name (`zone_base_lift_m={"North
  plains": 300.0}`). Recommended follow-up (not done here): add a real
  `base_lift_m` numeric column to the zones layer so it is visible/
  editable from QGIS like `amplitude_scale` is.

### 3c. Background noise

Ridged/turbulence Simplex fBm (`ridged_fbm` in `noise.py`): per octave,
fold the raw noise with `1 - 2*abs(n)` (signed fold, ranges [-1, 1]) —
NOT the naive `1 - abs(n)` fold, which is positively biased (empirically
mean ≈ 0.34, 93.5% of area > 0) and made clean ocean nearly impossible.
Even the signed fold isn't exactly zero-mean in practice, so the
generator additionally centers/scales using an empirically measured
mean/std (see §6 on why this has to be measured ONCE globally, not per
chunk).

`plain_octaves=2`: the first 2 (lowest-frequency, dominant) octaves skip
the fold entirely and use raw noise. Zero-crossings of a smooth noise
field are mathematically a **connected network** (like river deltas or
cracked mud), never isolated blobs — folding at every octave, including
the dominant one, produced "islands" that were actually one continuous
winding filament. Reserving the fold for finer octaves keeps the jagged
ridged texture as fine detail while letting large-scale island placement
come from plain noise, whose excursion sets are compact and separate.

### 3d. Domain warping

Distance-to-a-line/polygon (driving ridge decay and zone membership)
produces perfectly smooth, parallel "buffer" contours — real coastlines
are never that regular. Fix: query all three fields (ridge distance, zone
membership, background noise) at `(x, y)` displaced by a smooth
pseudo-random vector field (`domain_warp`, Inigo Quilez-style), not at
the true `(x, y)`.

Two independent warp layers, applied in sequence:

| layer | wavelength | amplitude | octaves | applied to |
|---|---|---|---|---|
| coastline warp | 15 km | 5000 m | 3 | ridge distance, zone membership, and as the base for the noise warp below |
| noise-detail warp | 4 km | 1500 m | 3 | background noise sampling coordinates only |

Both are multi-octave (fBm-style summed displacement), not single-
frequency — a single frequency has one characteristic bump spacing, which
produced a too-regular "scalloped" coastline (every lobe about the same
size, evenly spaced). Multiple octaves mix several spacings so
peninsula/bay size varies.

### 3e. Locked-in generation parameters

```
seed = 42
noise_octaves = 6
noise_base_wavelength_m = 25000
noise_amplitude_m = 200
sea_level_offset_m = 300
noise_plain_octaves = 2
lacunarity = 2.0, persistence = 0.5   (defaults, not overridden)
warp_wavelength_m = 15000, warp_amplitude_m = 5000
noise_warp_wavelength_m = 4000, noise_warp_amplitude_m = 1500
```

`sea_level_offset_m` was the main knob for island density/size — the
mainland is barely affected by it since its elevation is dominated by
ridge peaks (930–3794 m), but it strongly controls how much of the
background noise pokes above zero out in the open ocean.

## 4. Hydraulic erosion

Vectorized droplet-based erosion (Beyer/Lague algorithm — rainfall runoff
carving channels, picking up/depositing sediment). This is **fluvial
erosion only**: rain, gravity, flow, sediment capacity. It does **not**
model wave/tidal/coastal erosion — the unnaturally smooth peninsulas and
over-connected islands the DEM had before this stage were fixed upstream,
via the domain-warp and `plain_octaves` changes in §3, not by simulating
a separate coastal process.

Locked-in parameters for the final run:

```
n_droplets = 300000
seed = 1
max_steps = 427
erode_speed = 0.45
deposit_speed = 0.35
# left at library defaults:
inertia=0.05, sediment_capacity_factor=4.0, min_sediment_capacity=0.01,
evaporate_speed=0.02, gravity=4.0, initial_water=1.0, initial_speed=1.0,
spawn_above_m=0.0 (droplets spawn on land only), max_erosion_per_step_m=8.0
```

`max_steps=427` is not arbitrary: the reference algorithm's tunings were
validated at 200 m/pixel with `max_steps=64`. A droplet's real-world reach
is `max_steps * cell_size_m`, so preserving the same physical flow-path
length at 30 m/pixel requires `max_steps ≈ 64 * (200/30) ≈ 427`.

**Two numerical-stability bugs found and fixed during this stage**
(both in `src/terrain/erosion.py`, see its docstring for the full
mechanism):

1. **Overflow to `inf`**: the reference constants (gravity=4,
   sediment_capacity_factor=4, etc.) assume heightmaps normalized to
   roughly [0, 1]. Feeding this DEM's raw metre-scale elevations (up to
   ~3800 m) directly into the same constants caused a genuine runaway
   feedback loop (erosion deepens a pit → next droplet visiting it sees a
   bigger drop → more erosion → unbounded). Fixed by dividing elevation
   by `height_scale_m` before simulating and multiplying back after
   (`height_scale_m` defaults to this DEM's own elevation range), plus an
   independent `max_erosion_per_step_m` hard cap as defense in depth.
2. **NaN poisoning from dead droplets**: droplets that go out-of-bounds
   still had `delta_height` computed from their (clipped) position before
   being zeroed by the `alive` mask — and an unmasked value could reach
   `inf`, and `inf * 0.0 == NaN`, which then got scattered into the shared
   grid via `np.add.at` and corrupted every neighboring bilinear sample
   from then on. Fixed by masking `delta_height` to zero for dead
   droplets immediately after computing it, before it feeds capacity/
   speed calculations.

## 5. Resolution-independence testing before scaling up

Before running at the full 30 m grid, the ridge shelf widths, warp
parameters, island density, and erosion strength were all tuned at
200 m/64-step previews (~360x faster to iterate on), validated via
side-by-side comparison renders, then the erosion step count scaled per
§4's formula and re-validated with a full-resolution zoom crop before
committing to the full-domain run.

## 6. Chunked generation (memory)

At full 30 m resolution the grid is 5334 x 4334 = ~23.1M cells. The first
two attempts to generate it as one array were OOM-killed (`exit 137`)
around 5 minutes in — each `noise2()` call (several per octave, several
octaves, two warp layers) allocates several 23M-float64 temporaries; not
all of them coexist for long, but enough do at once to exceed this
sandbox's ~7 GB limit.

Fix: `generate_dem` processes the grid in row chunks (`row_chunk_size=250`
rows, ~1.08M points/chunk), bounding peak memory regardless of total
domain size. The one thing that has to happen OUTSIDE the chunk loop:
the ridged-fbm noise's mean/std normalization (§3c) is estimated ONCE from
a 400x400 domain-spanning random sample, then reused for every chunk —
normalizing per-chunk instead would center each chunk on its own local
statistics and produce a visible seam at every chunk boundary.

## 7. Final run and outputs

```
generate_dem(
    xmin=-65000, xmax=65000, ymin=-80000, ymax=80000, resolution_m=30,
    ridges_path="data/input/terrain_ridges.geojson",
    zones_path="data/input/terrain_zones.geojson",
    seed=42, noise_octaves=6, noise_base_wavelength_m=25000,
    noise_amplitude_m=200, sea_level_offset_m=300,
    shelf_multipliers={"Spine": 1.6, "North branch (Big Brother)": 1.3,
                        "West Branch (Little Brother)": 1.3, "South Branch": 1.3},
    zone_base_lift_m={"North plains": 300.0},
    warp_wavelength_m=15000, warp_amplitude_m=5000,
    noise_warp_wavelength_m=4000, noise_warp_amplitude_m=1500,
    noise_plain_octaves=2, row_chunk_size=250,
)
# -> (5334, 4334), 214.4s, land fraction 46.8%, no NaN/Inf

erode(dem, cell_size_m=30, n_droplets=300000, seed=1, max_steps=427,
      erode_speed=0.45, deposit_speed=0.35)
# -> 122.7s, no NaN/Inf, diff min=-149.03 max=+131.72 mean=-0.42 m
```

Final elevation range: **-853.3 m to +3843.7 m**; land fraction 46.7%.
(The peak sits ~50 m above the Spine's authored 3794 m — plausible: the
noise term still contributes up to ±~100 m even where the ridge structure
term dominates.)

### Export format

The natural full-precision export (ESRI ASCII Grid, already had a working
writer from earlier prototyping) comes out to ~185 MB for this grid
(`%.2f` text, ~8 bytes/cell) — too large to hand over comfortably as a
single deliverable. Exported instead as **ENVI raw binary** (`.bin` +
plain-text `.hdr`, no GDAL/rasterio needed to write, and QGIS reads it
natively):

- `dem_final_30m_eroded.bin` / `.hdr` — **int16**, elevation rounded to
  the nearest metre. 46 MB. Rounding error: uniform in [-0.5, 0.5] m,
  verified by round-tripping the file (max abs error 0.5 m, mean 0.25 m)
  — about 0.01% of this DEM's ~4700 m total relief, well below anything
  visible at 30 m horizontal resolution or introduced by the noise/
  erosion process itself.
- A float32 (full-precision, ~92 MB) export is available on request for
  downstream processing that needs to skip even that quantization.

**CRS caveat**: the `.hdr`'s `map info` field carries the affine
georeferencing (origin + 30 m pixel size) so the raster lands at the
right place/scale, but ENVI's header format has no simple way to embed a
full CRS from this sandbox without GDAL. A `.prj` sidecar with the
project's PROJ4 string is included alongside it — GDAL's ENVI driver
picks up a same-named `.prj` in many cases, but if QGIS loads the layer
with no/unknown CRS, assign "Fictional World LCC" manually via **Layer
Properties → Source → Assigned CRS**.

## 8. Notable tuning decisions along the way

- Spine's shelf tightened to 1.6x (from a wider first attempt) to narrow
  the main landmass's extent — chosen from a 3-way side-by-side
  comparison.
- Ridge endpoints for the North branch and West Branch were shortened in
  QGIS (not by editing generation parameters) after their falloff
  Gaussians were found to reach the domain edge — a genuine authored-
  geometry issue, not something noise/shelf tuning could fix.
- `sea_level_offset_m=300` chosen from an island-density sweep (260–320
  tested) to get "one large central landmass, several smaller islands"
  rather than either a single blob or an overly scattered archipelago.
- Erosion strength: "strong" preset (`erode_speed=0.45,
  deposit_speed=0.35`) chosen over default/light after a zoomed side-by-
  side comparison of drainage channel visibility.

## 10. v2 revision — closing the naturalism gap vs. a real reference

v1 (§§1–9) looked reasonable at a glance but was rejected on closer
inspection ("still some unnatural curvy shapes", visible once viewed at
proper zoom on a full monitor rather than a phone screen). This section
documents the diagnosis, the two approaches tried, and why the second
one was adopted for the DEM now in use.

### 10a. Diagnosis: it's the spectrum, not the erosion

The initial hypothesis (insufficient ocean-floor erosion coverage) was
checked and ruled out — land-vs-ocean erosion coverage was already
comparable. The real cause was found via 2D FFT radial power-spectral-
density analysis: terrain's PSD as a function of spatial frequency
follows a power law, `power ~ frequency^-β`, and β (the fractal
roughness/Hurst-related exponent) for this DEM's noise came out to
β≈3.72–3.74 — far steeper (i.e. smoother/more "curvy", less rough) than
any real terrain, at every wavelength, because `ridged_fbm` used one
fixed `persistence` for all octaves, implying one exponent everywhere.

Measuring the *actual* target required real data: two independent clean
1024×1024 patches from a real SRTM tile (see §10c for how that tile was
obtained) gave β≈2.50 above ~5 km wavelength (rougher — tectonic/
regional-scale relief) climbing to β≈4.54 below ~1.2 km (much smoother —
attributed to hillslope diffusion, a real physical smoothing process
acting at short range on top of the rougher large-scale structure). A
single-persistence fBm structurally cannot match both regimes.

### 10b. First fix tried: two-band persistence (kept in the code, not used in the final run)

`ridged_fbm` (`src/terrain/noise.py`) was extended with optional
`persistence_fine` / `crossover_wavelength_m` / `crossover_width_factor`
params: per-octave persistence is blended via a smootherstep function of
log-wavelength between `persistence` (large-scale) and `persistence_fine`
(fine-scale), centered on `crossover_wavelength_m` over a
`crossover_width_factor`-wide multiplicative range — smooth blending
rather than an abrupt per-octave switch, for the same reason `domain_warp`
uses multiple octaves rather than one (see §3d): a hard discontinuity is
itself a spectral artifact.

Calibrated against the real measurements above: `persistence=0.87`
(large-scale), `persistence_fine=0.50` (fine-scale), crossover at 4 km.
Validated on a 15 km test crop: measured β = 2.38 (large, vs. 2.50
target) and 4.03 (fine, vs. 4.54 target) — a good match, and a visible
improvement over v1 in side-by-side renders (`compare_twoband_final.png`,
`compare_real_vs_ours_final.png`).

**This alone was still not good enough.** Direct visual comparison
against the real SRTM tile (not just matching spectral numbers) showed
real terrain has a consistent directional "grain" — parallel, similarly-
oriented ridge/valley structures — that isotropic noise cannot produce,
*no matter how well its radial (direction-blind, by construction) power
spectrum is calibrated*. This is a structural ceiling of any noise field
built the way `ridged_fbm` is, not a tuning problem — confirmed by
checking the skeleton's own ridge bearings (Spine 149.6°, North branch
90.5°, West Branch 157.2°, South Branch 76.9°) against the real texture's
angular spectrum (broadly spread ~70–170°, not one sharp direction) and
finding no clean rotation would fix it either way.

The two-band persistence code is left in `noise.py` (harmless, opt-in via
`None` defaults) since it's a genuine, validated improvement over a flat
persistence and may be useful background noise for some other terrain
that doesn't have real reference data available.

### 10c. Adopted fix: grafting real DEM detail as texture

Rather than continue tuning synthetic noise, the fine-scale **detail**
term was replaced with real elevation data, while the **placement** of
land/mountains/coastline stayed 100% controlled by the hand-authored
ridge/zone skeleton — same principle as v1's `structure(x,y)` term,
just a different source for the noise/detail term added on top of it.

**Source data**: a real SRTM GL3 (3 arc-second, ~90 m) GeoTIFF tile
covering the Piemonte/western Alps, supplied from
`C:\projects\worldbuilding\data\raw\SRTMGL3.tiff` (12.7 MB, LZW-
compressed, EPSG:4326). Processing chain (all intermediate files stayed
in the generation sandbox, not committed to the repo):

1. Resampled to isotropic ~92.767 m/px (SRTM GL3 pixels are NOT square in
   real-world metres away from the equator — this tile's cell is
   ≈65.60 m E–W × 92.77 m N–S at ~45°N).
2. Voids (NODATA) filled via nearest-valid-cell (`scipy.ndimage.
   distance_transform_edt`) — one ~253k-cell void blob found, confined to
   one tile corner and avoided entirely by the window used below.
3. **Detrended**: subtracted a 15 km-sigma Gaussian blur of itself, which
   removes "where the big massifs are" (that's this fictional world's
   job, via the skeleton) while keeping 100% of the anisotropic ridge/
   valley texture from 15 km down to the source's native resolution.
4. A clean 1725×1401 window (160.02 km N–S × 129.97 km E–W, zero void
   cells) — coincidentally large enough to cover this project's entire
   130×160 km domain in one piece, so **no tiling or seam-blending was
   needed** — was saved as `data/real_detail/piemonte_detail_15km_
   detrend.npy` (float32, mean=16.00, std=518.74), with metadata in
   `piemonte_detail_meta.json`.

**Coherence/erosion/3D-compatibility** (raised as a direct question
before committing to this approach): preserved by construction, for
three independent reasons —
- The final elevation is still a single well-behaved sum,
  `structure(x,y) + real_detail(x,y) * amplitude`, mathematically
  indistinguishable in shape from v1's `structure + synthetic_noise`.
  Real and structure fields are added, not stitched or masked.
- No seams: the real tile covers the whole domain in one contiguous
  window (point 4 above), so there's no tile-boundary blending logic to
  get wrong.
- The erosion algorithm (`erosion.py`) only ever reads local elevation
  differences between grid neighbors — it has no notion of "synthetic"
  vs. "real" and cannot behave differently based on where a value came
  from. Same for any 3D-model/raster consumer downstream: they see one
  standard single-band elevation grid, full stop.

Validated first on the same 15 km test crop used for the two-band
approach — `compare_c2_realtexture.png` shows a dramatic, unmistakable
improvement in matching the real reference's parallel-ridge grain versus
both v1 and the two-band synthetic attempt.

### 10d. Wiring into `generate_dem()`

`generate.py` gained: `real_detail_path`, `real_detail_xmin`,
`real_detail_ymax`, `real_detail_cellsize_m` (world→array coordinate
mapping — row 0 = north/`real_detail_ymax`, col 0 = west/
`real_detail_xmin`, same convention as this function's own output),
and an optional `real_detail_fine_supplement_weight` /
`real_detail_fine_min_wavelength_m` pair.

**Fine-detail supplement**: real data has no information below its own
native resolution (~185 m Nyquist here). Since this DEM's target
resolution (30 m) is finer than that, a *small* amount of synthetic
ridged-fBm noise is added, confined strictly to wavelengths the real data
can't supply (from real data's own Nyquist down to
`real_detail_fine_min_wavelength_m=60 m`), so the last mile isn't
artificially smoother than genuine terrain roughness at that scale
without displacing real signal anywhere it actually has coverage. Its
octave count is derived (`ceil(log(fine_start/fine_min) / log(lacunarity))`),
and — same seam-avoidance reasoning as the noise-normalization step in
§6 — its mean/std are estimated ONCE globally, not per chunk.

**Zone `amplitude_scale` still applies for free**: `chunk_elevation =
structure + noise_z * noise_amplitude_m * amplitude_multiplier -
sea_level_offset_m` is *exactly the same line* as v1 — `noise_z` just
comes from a different source depending on whether `real_detail_path` is
set. No special-casing needed downstream of that assignment.

One side-effect worth flagging: the existing "noise-detail" domain warp
(§3d) is applied to the *query coordinates* before they're used to sample
the real-detail array too — meaning the real texture also gets a small
organic displacement where it's read from, not just synthetic noise as
in v1. This was left as-is (not a bug) since it adds the same kind of
irregularity domain warping was already providing.

### 10e. Final v2 run and outputs

```
generate_dem(
    xmin=-65000, xmax=65000, ymin=-80000, ymax=80000, resolution_m=30,
    ridges_path="data/input/terrain_ridges.geojson",
    zones_path="data/input/terrain_zones.geojson",
    seed=42, noise_amplitude_m=200, sea_level_offset_m=300,
    shelf_multipliers={"Spine": 1.6, "North branch (Big Brother)": 1.3,
                        "West Branch (Little Brother)": 1.3, "South Branch": 1.3},
    zone_base_lift_m={"North plains": 300.0},
    warp_wavelength_m=15000, warp_amplitude_m=5000,
    noise_warp_wavelength_m=4000, noise_warp_amplitude_m=1500,
    warp_octaves=6, warp_lacunarity=1.75, warp_persistence=0.7,
    noise_warp_octaves=5, noise_warp_lacunarity=1.75, noise_warp_persistence=0.7,
    real_detail_path="data/real_detail/piemonte_detail_15km_detrend.npy",
    real_detail_xmin=-65000, real_detail_ymax=80000, real_detail_cellsize_m=92.767,
    real_detail_fine_supplement_weight=0.15, real_detail_fine_min_wavelength_m=60.0,
    lacunarity=1.75, row_chunk_size=250,
)
# -> (5334, 4334), 237.1s, elevation -953..4116 m, land fraction 46.0%, no NaN/Inf
```

**Note on unused-looking arguments**: `noise_octaves`, `noise_base_
wavelength_m`, `persistence`, `noise_persistence_fine`, `noise_
crossover_wavelength_m`/`_width_factor`, `noise_plain_octaves` only take
effect when `real_detail_path` is `None` (the v1/synthetic-only path,
§§1–9) — when real-detail texture is supplied they're bypassed entirely,
*except* `lacunarity`, which is still used by the fine-detail supplement
above. `run_full_c2.py` (the script actually run) still passes v1-era
values for the now-bypassed ones; harmless, but worth cleaning up if this
script is reused, to avoid implying they still matter.

Erosion was re-run at a **denser droplet count** than v1: `0.3` droplets
per land cell (3,189,482 droplets total, vs. v1's fixed 300,000 ≈
0.028/land-cell) — a deliberate choice for this run, not a default,
chosen to make sure the much sharper real-detail texture gets properly
carved rather than under-eroded relative to its own resolution:

```
erode(dem, cell_size_m=30, n_droplets=3189482, seed=1, max_steps=427,
      erode_speed=0.45, deposit_speed=0.35)
# -> ~1275s, no NaN/Inf, diff min=-332.42 max=+362.76 mean=-5.51 m
```

The larger mean erosion diff versus v1 (-5.51 m vs. -0.42 m) is an
expected consequence of the ~10x-denser droplet count doing more net
carving, not a new instability — sanity-checked the same way as v1 (no
NaN/Inf, land fraction stable).

**Final v2 elevation range: -951.7 m to +4043.6 m; land fraction 45.3%**
(total runtime 1512s: 237s generation + 1275s erosion). Exported the same
way as v1 (§7): `dem_c2_final_30m_eroded.bin`/`.hdr`/`.prj`, int16, same
rounding-error profile.

### 10f. Methodological note: this DEM's fine texture is now tied to a real place

Worth being explicit about for anyone extending this pipeline later:
below the ~15 km scale the skeleton controls, this fictional world's
terrain texture *is*, literally, the detrended shape of a real place
(Piemonte/western Alps) — not an abstraction of "alpine-like terrain in
general". That's a reasonable trade for visual quality, but it means:
reusing `piemonte_detail_15km_detrend.npy` for a *different* region of
this same fictional world would reproduce the same real drainage pattern
twice (a repeated "fingerprint"), and anyone doing close, feature-level
scrutiny of the terrain (rather than an overview) could in principle
recognize the source landscape. Neither concern applies to normal use
(overview maps, gameplay, most 3D renders), but both are worth keeping in
mind before, say, publishing a very-high-resolution crop side-by-side
with real Piemonte imagery, or reusing this same detail asset for a
second, unrelated part of the map.

## 11. v3 revision — domain-warp folding bug found after v2

Direct visual inspection of v2 in QGIS (rather than the render crops alone)
surfaced a problem described at the time as "a lot of curvy shapes that
seem more like a natural pattern (like some coral or a biological
pattern), than a terrain surface." That description turned out to be a
precise, technically accurate diagnosis of a real bug, not a vague
quality complaint.

### 11a. Diagnosis

`domain_warp` (§3d) displaces sampling coordinates by a smooth
pseudo-random vector field before they're used for anything else —
including, since §10, sampling the real-detail texture. A domain warp
like this is only well-behaved (a smooth, locally-invertible coordinate
change) if its own spatial gradient stays below a threshold; past that
threshold the map **folds onto itself** — some region gets sampled from
two different directions at once, the mathematical signature of which is
a Jacobian determinant that goes negative.

Measured directly (finite-difference Jacobian of the actual warp field
used in production): **~35% of the entire domain had a folded warp map**
(measured at 100 m sampling, whole-domain). Root cause, in order of
contribution:

1. `noise_warp_persistence (0.7) * noise_warp_lacunarity (1.75) = 1.225 >
   1` for both warp layers — meaning each successive octave contributes
   *more* gradient than the last, the opposite of how fBm-style noise is
   supposed to behave (higher octaves are meant to add fine detail with
   *diminishing* weight, not dominate).
2. Amplitude too large relative to wavelength on top of that (coastline
   layer: 5000 m amplitude at 15000 m wavelength = 33%; noise-detail
   layer: 1500 m at 4000 m = 37.5%).

**Why this was invisible in v1 and even in early v2 renders**: folding a
*synthetic, directionless* simplex noise field is optically undetectable
— a folded region just looks like more noise, indistinguishable from an
unfolded one. Once real, recognizable DEM texture (§10c) is what's being
sampled through that folded map, the same distinctive real ridge/valley
shapes get sampled from multiple folded directions at once, producing
exactly the kind of closed rings, blobs, and maze-like patterns described
as "coral" or "biological" above — this is a known visual signature of
coordinate-map folding (sometimes called caustics), not a new/different
kind of artifact. Confirmed three independent ways before touching any code:
- Direct visual correlation between a fold map (cells with negative
  Jacobian) and the "coral" appearance in the actual v2 DEM, in a
  low-relief crop (`diag_fold_vs_coral_plains.png`) — the two patterns
  match essentially cell-for-cell.
- A whole-domain fold map at 100 m resolution
  (`fold_map_full_domain.png`): the folded-cell pattern for the
  production warp, viewed on its own with no elevation data at all, is
  already a dense, uniform, domain-covering maze/coral pattern.
- Closed circular ring/blob artifacts visible directly in eroded hillshade
  crops (`compare_warpfix_eroded.png`, `compare_coastline_warpfix.png`) —
  real drainage never forms a fully closed ring (water always has to
  reach an edge); a closed ring is a structural tell of folding, not of
  any real erosional or tectonic process.

### 11b. Fix

Scaled both warp layers' amplitude to 40% of their v2 values, leaving
wavelength/octaves/lacunarity/persistence untouched:

```
warp_amplitude_m: 5000 -> 2000        (wavelength unchanged, 15000 m)
noise_warp_amplitude_m: 1500 -> 600   (wavelength unchanged, 4000 m)
```

This was chosen empirically (a parameter sweep over amplitude scale and
over lacunarity/persistence combinations), not derived in closed form —
closed-form fold-free bounds exist for single-octave warps but get messy
fast for multi-octave summed noise. Validated at native 30 m resolution
before committing to a full run, across five separate regions chosen to
cover different parts of the domain (a plains/low-relief zone, the
highest peak, a coastline stretch, the largest offshore island, and the
West Branch ridge area) plus a whole-domain fold map:

| region | fold %, v2 (production) | fold %, fixed |
|---|---|---|
| whole domain (100 m sample) | 35.1% | 1.3% |
| north plains | 40.3% | 2.2% |
| peak | 39.2% | 2.1% |
| coastline | 39.1% | 1.8% |

Residual folding (~1.3-2.2%) is not mathematically zero, but is ~18-20x
lower than production, evenly spread across the domain (1.2-1.4% in every
quadrant, no hotspots) rather than concentrated, and produced no visible
ring/blob artifacts in any of the five test crops, eroded or not. A more
conservative amplitude (25-30% of v2, not yet tested) would very likely
reduce it further, at the cost of somewhat less organic coastline
variation — not pursued here since 40% already tested clean, but worth
knowing if a future close inspection still finds something.

Land fraction was essentially unaffected in every test crop (e.g. island:
27.9% -> 28.0%; West Branch: 21.4% -> 20.1%) — this fix changes *where
exactly* warped coordinates land, not the overall balance of the
elevation model, so the previously-approved coastline/island layout
should stay recognizable, just without the fold artifacts. Exact fine
contours did shift somewhat (this is a genuinely different coordinate
sampling, not a texture-only patch), which is expected and unavoidable.

### 11c. Final v3 run and outputs

Identical to §10e's `run_full_c2.py` in every parameter except the two
warp amplitudes above (see `run_tappa1_terrain.py`, renamed from
`run_full_c2_v3.py` once Tappa 1 was locked in):

```
generate_dem(..., warp_amplitude_m=2000, noise_warp_amplitude_m=600, ...)
# -> (5334, 4334), 232.0s, elevation -953..4145 m [pre-erosion], land ~45.9%, no NaN/Inf

erode(dem, cell_size_m=30, n_droplets=3183498, seed=1, max_steps=427,
      erode_speed=0.45, deposit_speed=0.35)
# -> ~1324s, no NaN/Inf, diff min=-272.92 max=+234.20 mean=-5.56 m
```

**Final v3 elevation range: -945.6 m to +4020.2 m; land fraction 45.3%**
— both numbers essentially unchanged from v2 (-951.7..4043.6 m, 45.3%),
consistent with the "layout preserved, fold artifacts removed" expectation
above. Exported the same way as v1/v2: `dem_v3_final_30m_eroded.bin`/
`.hdr`/`.prj`, int16, same rounding-error profile as before.

**This revision (v3) is the current, active DEM.** §10's v2 files are
superseded in turn, kept only as the historical record of the C2
real-detail-texture decision (which itself is still correct and unchanged
in v3 — only the domain-warp amplitude changed).

## 12. NZ real-detail texture investigated and rejected

This project's own `validation_reference` in `config/parameters.yml` names
South Island, NZ — not Piemonte — so once v3 shipped, the obvious question
was raised: should the real-detail texture (§10c) come from an actual
South Island DEM instead, for methodological consistency? Investigated
properly rather than assumed either way.

**Data**: a real SRTM GL3 GeoTIFF (`SRTM_NZ.tif`, uncompressed, 4°×3° box,
169.0–173.0°E / 45.0–42.0°S, covering the Aoraki/Mount Cook massif and
southward) was sourced via the device bridge, same workflow as
`SRTMGL3.tiff`. Processed identically to the Piemonte asset (isotropic
resample to 92.767 m/px, void-fill, 15 km-sigma detrend).

**Problem found — not a bad window, a regional characteristic**: this
tile's NODATA fraction is far higher than Piemonte's (39.7% vs. Piemonte's
near-zero), and unlike Piemonte, no 130×160 km sub-window in the tile is
anywhere near void-free — a systematic scan found exactly one distinct
low-void region (0.23% void, near Lake Hawea/Wanaka, ~169.5–171.1°E /
43.55–45.0°S). Grafting it onto the fictional skeleton produced an
obvious, singular, dead-straight valley/lake gash cutting across the
terrain (`compare_nz_graft.png`) — visually nothing like Piemonte's more
varied, less singular texture. A second, geographically distinct
candidate near the Mackenzie Basin (6.1% void, requiring more fill) showed
the *same* problem — at least two more large, similarly dead-straight
lake-troughs (`nz_mackenzie_preview.png`), plus a nearest-neighbor void-
fill radiating-artifact near a larger void region. Conclusion: this isn't
an unlucky window — the Southern Alps corridor near Aoraki/Mount Cook is
genuinely dominated by a handful of very large, fault-controlled, dead-
straight glacial troughs (Pukaki, Tekapo, Hawea, Wanaka, and others), at a
scale and singularity Piemonte's sampled terrain didn't happen to have.
Any large clean-enough window from this specific tile is likely to catch
one.

Same-method spectral comparison (a quick reimplementation, not the exact
script behind §10a's published 2.50/4.54 figures — re-measuring Piemonte
with it gave 2.87/4.91, so treat these as *relative*, not consistent with
the historical absolute numbers): NZ came out smoother in this window
(β=3.12 large-scale, 5.09 fine-scale, residual std 349 m) than Piemonte
(β=2.87, 4.91, std 519 m) — consistent with the visual impression of a
texture dominated by one huge trough plus otherwise gentler surroundings,
rather than Piemonte's more uniformly-distributed roughness.

**Mitigation tried**: an asymmetric soft-clip (`tanh`-based) on the
detrended residual, compressing negative (valley) excursions harder
(limit 1.2σ) than positive (ridge) ones (limit 2.5σ) — targeting the
valleys specifically since they were the visually dominant problem, not
the peaks. This roughly halved the trough's depth (residual min went from
-889 m to -407 m) and removed the stark white "gash" contrast in a plain
hillshade render (`compare_nz_valleysuppress_texture.png`). **Only a
partial fix**: amplitude compression changes how deep the trough reads,
not its shape — the exact same dead-straight, single, unusually long
line is still there, just shallower (`compare_nz_3way.png`). A real fix
would need to break the trough's *geometry* (e.g., a local warp targeted
specifically at the most linear/elongated low-lying features), which is
materially more engineering with an uncertain payoff — the same
diminishing-returns territory as §13's (below) deferred coastal
micro-warp, just applied to an interior feature instead of the coast.

**Decision: kept Piemonte.** The methodological inconsistency (this
world's validation reference is South Island, NZ, but its mountain
texture is Piemonte's) stays, documented and open-eyed rather than
silently accepted — see `piemonte_detail_meta.json` /  §10f for what that
means in practice. Assets from this investigation
(`data/real_detail/nz_detail_15km_detrend.npy`,
`nz_detail_wanaka_valleysuppressed.npy`, `nz_detail_meta.json`) are kept
in the sandbox as a record but are **not** wired into `generate_dem()` —
`piemonte_detail_15km_detrend.npy` remains the one actually used by
`run_tappa1_terrain.py` and `config/parameters.yml`.

## 13. Open follow-ups (not done in this stage, deliberately left open)

- **Coastline still reads as somewhat "blobby"/regular** compared to v2's
  (pre-fold-fix) coastline — a direct, known side-effect of §11's fix
  (safe warp amplitude = less organic coastline variation than the
  unsafe, folding amplitude gave). Investigated one alternative during
  this stage (a vertical elevation-band perturbation near the 0 m
  contour) and rejected it after testing: masking by elevation value
  alone seeds spurious islands out in open ocean (any cell that happens
  to sit near 0 m by noise chance gets perturbed, not just true coastal
  cells); masking by true geometric distance-to-coast fixes that but
  still reads as a chain of small volcanic-looking "pimples" right at
  the shore, because a *vertical* bump is the wrong mechanism — real
  coastline irregularity (bays, headlands) is fundamentally a
  *horizontal* displacement effect, not a vertical one. A **horizontal,
  distance-to-coast-masked micro-warp** (small, safe amplitude, applied
  only within a band near the existing coastline, decoupled from the
  main domain-warp so it can't reintroduce the real-detail-texture fold
  bug in §11) tested visually clean in a quick preview
  (`compare_coastal_microwarp.png` — via cheap post-hoc resampling of the
  already-built v3 DEM, NOT yet a real `generate_dem()` implementation).
  Deliberately not pursued further/wired into the pipeline this stage —
  v3's coastline was judged safe-but-plain, good enough to close on. If
  picked up later: implement as a proper third warp layer (own
  wavelength/amplitude/reach params) applied to `xy_query` before ridge/
  zone/real-detail sampling, masked by `scipy.ndimage.distance_transform_edt`
  on the land/sea boundary, and validate across the same 5 regions used
  for §11's fix before committing to a full run.
- **North plains zone reads as mountainous again**, same root cause as
  before this revision: `amplitude_zone` only scales noise amplitude —
  it never flattens the ridge-distance `structure` term, and this zone
  sits within a nearby ridge's shelf reach, so the underlying structure
  still dominates. Not a code bug — the planned fix is a GeoJSON
  geometry/type edit (move the zone further from ridge reach, switch it
  to `plateau` with a low `target_elevation_m`, or both), rather than
  patching this with more Python-side special-casing.
- Promote `base_lift_m` to a real attribute column on the zones layer
  (see §3b).
- `config/parameters.yml`'s `terrain:` block still has `spine_path: null`
  — that field is now superseded by the GeoJSON skeleton
  (`terrain_ridges.geojson`/`terrain_zones.geojson`) and should probably
  be removed or flagged as deprecated. `seed`, `noise_octaves`,
  `noise_frequency` (→ `noise_base_wavelength_m` here),
  `erosion_iterations` (→ `max_steps`/`n_droplets` here) should be filled
  in with this stage's locked-in values so the config file stays the
  single source of truth the README says it is. Since v2 is now the
  active DEM, also consider adding fields for `real_detail_path` and
  `real_detail_fine_supplement_weight` so the "which real texture asset
  is this world using" fact lives in the same config file, not just in
  this doc and `piemonte_detail_meta.json`.
- Verify against the Tappa 0 validation reference (South Island, NZ) for
  a numeric plausibility check beyond visual comparison — not attempted
  here.
- If a *second* region of this world later also needs grafted real
  detail, source a different real tile for it (see §10f) rather than
  reusing the Piemonte asset, to avoid a repeated "fingerprint".
