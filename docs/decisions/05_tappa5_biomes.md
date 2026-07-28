# Tappa 5 — Biome classification

Status: **closed**. A single biome layer at the 120 m climate working
resolution, derived from Tappa 2's monthly temperature/precipitation stacks
via the Holdridge Life Zone system (Holdridge 1967), fine-tuned to this
world's own data, with permanent snow/ice overridden from Tappa 3's
mass-balance mask. Two-axis classification as scoped in
`00_pre_project_planning.md` (temperature + precipitation), correcting the
original (May 2025) project's temperature-only limitation.

This document is the decision record for the stage, per the project's
per-stage workflow. As with the four prior documents it is written to work
either as a log of what was decided, or as a reusable recipe.

## 0. Why Holdridge, not the more commonly tutorial-cited Whittaker diagram

Both take temperature and precipitation as input. The difference is whether
the result is reproducible:

* **Whittaker's biome diagram** has no closed-form boundaries. Its regions
  were drawn by inspection on a 1970s textbook figure; every digital
  redrawing in circulation disagrees slightly with the others because there
  is no authoritative numeric definition to trace back to. Adopting it would
  mean digitizing polygon vertices off somebody's redrawn image -- the same
  eyeballing problem this project already declined once, for a comparable
  reason (`04_tappa4_hydrology.md` S2, the rejected area x slope^2
  channel-initiation threshold).
* **Holdridge Life Zones** reduce to two numbers -- biotemperature and
  precipitation -- placed on a log-spaced grid, with potential
  evapotranspiration (PET) computed *from* biotemperature by a fixed
  constant, not measured independently. Fully reproducible: the same inputs
  always give the same zone, no image involved.

Holdridge was implemented as two separate modules, deliberately:

* `src/biomes/holdridge.py` -- a general-purpose, world-agnostic
  implementation. Nothing in it is specific to this project; it would work
  unchanged on a different domain.
* `src/biomes/world_biomes.py` -- everything specific to THIS world, built by
  looking at what this world's own data actually produces rather than
  importing a scheme wholesale.

## 1. The general Holdridge module

Formulas and constants, and what was actually verified this session (no
primary access to Holdridge's 1967 monograph):

* **Biotemperature**: annual mean of monthly means, each month clipped to
  [0, 30] deg C first (most plants are dormant below freezing and above
  ~30 deg C). Confirmed against two independent secondary sources: the
  Wikipedia summary of Holdridge (1967), and Aguirre et al. (2025,
  *Scientific Data*, a global climate-classification uncertainty paper)'s
  methods section, verbatim: "those months with a mean temperature above
  30.0 C and below 0.0 C are considered as 30.0 C and 0.0 C, respectively."
* **PET (mm/yr) = 58.93 x biotemperature**. The 58.93 constant is Holdridge's
  own empirical calibration, confirmed verbatim in the same 2025 paper.
  **Known weak spot, not silently ignored**: this constant was fit mostly
  against tropical/subtropical stations; the resulting PET ratio is known to
  get numerically fragile at low biotemperature (below ~3 deg C), because
  PET itself collapses toward zero there. This world's land-mean
  biotemperature (7.39 deg C) sits in a range where this is a live concern,
  not a hypothetical one -- see S3.
* **Belt** (temperature): 7 named bands (Polar, Subpolar, Boreal, Cool
  Temperate, Warm Temperate, Subtropical, Tropical), boundaries at 1.5 / 3 /
  6 / 12 / 18 / 24 deg C biotemperature. The 1.5-3-6-12-24 sequence (a clean
  doubling each step) is directly confirmed via an EPA/ORNL life-zones
  report; that report's own summary table collapses Warm Temperate and
  Subtropical into one row, so the interior 18 deg C split is reproduced
  from secondary-literature consensus, not independently re-derived.
* **Humidity province** (moisture): 8 named bands by PET-ratio
  (Superhumid..Superarid), boundaries at powers of 2 from 0.5 to 32. Log2
  spacing and the 0.125-32 axis range corroborated by multiple secondary
  sources; Holdridge's original monograph not itself accessed.
* **No official compound zone names** (e.g. "Boreal Wet Forest"): Holdridge's
  diagram assigns ~30-38 traditionally-named zones to specific belt x
  humidity cells, read directly off the printed diagram. Without pixel-exact
  access to it, inventing those names risked getting them wrong. The module
  returns belt and humidity-province as two independently-verified axes;
  compound official-style names are not fabricated.
* **Latitudinal vs. altitudinal naming**: Holdridge's diagram uses different
  NAMES for the same biotemperature bands depending on whether the cold
  comes from latitude (Polar/Subpolar/Boreal/...) or elevation on a
  warmer-based mountain (Nival/Alpine/Subalpine/...) -- numeric boundaries
  identical either way, and which label set is "correct" depends on the
  biotemperature of the region's own base. Reproducing that correctly needs
  the primary diagram; `holdridge.py` always returns the latitudinal name
  set, and the renaming appropriate to this specific world happens in
  `world_biomes.py` (S4) rather than being guessed generically.

## 2. What running the unmodified module on this world's data found

Applying the general module directly (before any world-specific tuning) to
the Tappa 2 monthly stack (land mean biotemperature 7.39 deg C, range
0.00-12.00; land mean PET 435.6 mm, range 0.0-707.3):

| belt | area |
|---|---|
| Cool Temperate | 6,487 km² (65%) |
| Boreal | 1,445 km² (15%) |
| Subpolar | 730 km² (7%) |
| Polar | 1,243 km² (13%) |
| Warm Temperate | ~0.06 km² (4 cells, right at the 12°C belt edge) |

The belt axis is fine -- well distributed across elevation, no fix needed.

**The moisture axis is not.** Land PET ratio ranges 0.00-2.18 (mean 0.56).
Every land cell falls in Superhumid/Perhumid/Humid/Subhumid --
**nothing on this island clears Holdridge's Semiarid threshold (ratio >=
4)**, including the driest cell on record (324 mm/yr precipitation --
the same cell Tappa 2-4 validated against real Alexandra's ~300 mm). That
cell's biotemperature is 11.72 deg C, giving PET = 58.93 x 11.72 ~= 690 mm,
ratio ~= 2.13 -- Subhumid, not Semiarid.

**Diagnosis, not a bug**: Holdridge's PET scales linearly with
biotemperature, and this world's biotemperature never exceeds 12 deg C
(elevation-controlled, per `02_tappa2_climate.md`), capping PET near
700 mm/yr even at the driest, lowest-elevation cell. Holdridge's global
calibration implicitly assumes dry places are usually also hot (deserts); a
cool maritime island can be genuinely dry in absolute precipitation terms
without ever generating enough PET to register as arid on Holdridge's own
scale. (For context, not as a verified citation: applying Köppen's BSk
aridity formula, which scales more gently with temperature, to the same
real Alexandra numbers -- ~10-11 deg C, ~300 mm, spread fairly evenly through
the year -- gives a threshold near 490 mm and a semi-arid classification;
Köppen and Holdridge appear to disagree on whether a cool, modestly-dry
place like this counts as "arid" at all. Derived, not independently
cross-checked against a named source this session.)

Taken at face value this would erase the entire windward/leeward moisture
contrast Tappa 2-4 spent three stages establishing and validating (5.6:1
precipitation ratio `02_tappa2_climate.md` S5, 1136 m ELA differential
`03_tappa3_snow.md` S5, 5.9:1 discharge ratio `04_tappa4_hydrology.md` S3) --
the biome map would show zero moisture-driven variation, only elevation.

Separately: Holdridge's own Polar belt (biotemperature < 1.5 deg C) is
1,243 km², **~30% larger** than Tappa 3's physically-modelled permanent-snow
area (960 km², 68% cell overlap) -- the identical naive-temperature-threshold
error `02_tappa2_climate.md` / `03_tappa3_snow.md` already caught once (naive
permanent snow 108 km² vs. mass-balance 960 km²), now reappearing if Polar
were used directly as the ice/nival class.

## 3. World-specific fix: rebin the moisture axis by this world's own data

The PET ratio is kept as the underlying quantity -- it still correctly
*ranks* wet vs. dry cells -- but binned by this world's own land-cell
percentiles rather than Holdridge's desert-calibrated absolute scale. This is
the same move `03_tappa3_snow.md` / `04_tappa4_hydrology.md` already made for
their windward/leeward split (annual-precipitation tercile), applied here to
the PET ratio instead of raw precipitation because the ratio also carries
the temperature-driven part of "how available is moisture," not just the
precipitation total.

Permanent snow/ice is overridden from Tappa 3's mass-balance mask rather
than read off Holdridge's Polar belt, for the reason quantified in S2.

## 4. Two rejected drafts before the locked scheme

**Draft 1**: quartiles (Superwet/Wet/Moist/Dry) x belt, Polar/Subpolar
unsplit, no Warm Temperate fold -- 10 land classes + permanent snow. Rendered
as a hillshade + biome overlay and reviewed visually. Two problems surfaced:

* A "Temperate Rainforest" class (Cool Temperate x Superwet) came out at
  96 km² -- a sliver that was, in the delivered render, visually
  indistinguishable from the adjacent Subalpine Cloud Forest class.
* Running the palette through the project's dataviz color-validation
  tooling (`scripts/validate_palette.js`, all-pairs mode -- the correct mode
  for a choropleth, where any two categories can end up spatially adjacent)
  failed on multiple axes even before the rainforest problem: worst
  normal-vision pair (Subalpine Forest vs. Temperate Forest) measured
  Delta-E 4.2, below the "hard to tell apart even with full color vision"
  floor of 15.

**Draft 2**: folded the Superwet tier into Wet (so the rainforest sliver
merges into "Temperate Forest" instead of vanishing as noise), 4 tiers ->
per the user's direction ("2 moisture is fine ... 3 if you feel it
simplifies too much; merging [the rainforest/cloud-forest pair] is also fine
by me"). Evaluated 2 vs. 3 tiers on the actual data: 2 tiers (median split)
would have collapsed the Cool Temperate belt's "transitional woodland" and
"open grassland" distinction into one class, discarding legible structure in
the single belt carrying 65% of this world's land area and the clearest
windward/leeward narrative signal established over three prior stages.
3 tiers (terciles), applied uniformly to Boreal and Cool Temperate and left
unsplit for Polar/Subpolar (S5), was adopted instead -- same method
Tappa 3/4 already used for windward/leeward, extended here to a third band.
Final palette iteration (S6) also separated Alpine Tundra (recolored to a
heather/purple hue) and Subalpine Dry Scrub (recolored to rust) out of the
yellow-green cluster that was crowding the draft-1 palette.

## 5. Locked scheme

| id | biome | belt | moisture tercile | area (km²) | % of land |
|---|---|---|---|---|---|
| 1 | Permanent Snow & Ice | (mass-balance override) | -- | 960.2 | 9.7% |
| 2 | Alpine Fellfield | Polar | unsplit | 396.0 | 4.0% |
| 3 | Alpine Tundra | Subpolar | unsplit | 616.5 | 6.2% |
| 4 | Subalpine Wet Forest | Boreal | Wet | 942.4 | 9.5% |
| 5 | Subalpine Woodland | Boreal | Moist | 167.2 | 1.7% |
| 6 | Subalpine Dry Scrub | Boreal | Dry | 335.1 | 3.4% |
| 7 | Temperate Forest | Cool Temperate | Wet | 689.8 | 7.0% |
| 8 | Woodland / Shrubland | Cool Temperate | Moist | 2,834.8 | 28.6% |
| 9 | Lowland Steppe / Grassland | Cool Temperate | Dry | 2,962.2 | 29.9% |

Moisture tercile edges (PET ratio, land cells): Wet < 0.1184, Moist
0.1184-0.5360, Dry > 0.5360. Warm Temperate (4 cells, ~0.06 km²) folded into
Cool Temperate. Polar and Subpolar are not split by moisture: checked
directly, 0 km² of Polar land falls in the driest tercile, and the areas
involved are modest enough that splitting would add categories without
adding legible distinction.

**Fragmentation** (connected-component check, 4-connectivity, before
shipping): every class's area sits overwhelmingly in one or two large
patches -- under 0.1% of any class's area falls in patches smaller than
1 km². The map reads as coherent zonation, not classification noise.

## 6. Palette: an accepted, documented trade-off

Colors (`src/biomes/world_biomes.py::BIOME_COLORS_HEX`) were chosen for a
natural, continuous wet-green -> dry-tan cartographic read, per the
pre-project plan's "qualitative palette, nominal data" brief for this layer.
Run through the project's dataviz color-validation tooling in all-pairs mode
(the correct mode for a choropleth): **fails** strict colorblind-safety
checks -- worst pair (Temperate Forest vs. Subalpine Wet Forest) measures
Delta-E 5.6 normal-vision, below the 15 floor. This is structural, not a
fixable hex tweak: several categories are deliberately shades of green
because they represent a forest-density gradient, and a natural biome map
reads better with a green-to-tan continuum than with unrelated hues spread
around the color wheel for CVD-safety's sake alone. Two iterations did pull
the worst offenders (Alpine Tundra, Subalpine Dry Scrub) out of the
yellow-green cluster into distinct hue families (heather, rust) once they
were found to be genuinely unreadable against their neighbors on the actual
map -- this is documented as a partial, evidence-driven fix, not a full
pass. **Mitigation**: the pre-project plan's own click-to-select detail
panel and always-visible legend (`00_pre_project_planning.md`) mean no
reading of the interactive map depends on color alone. Stated here rather
than silently accepted.

## 7. Validation against South Island, NZ

Same wind-direction decision as Tappa 2-4 (`02_tappa2_climate.md` S1, 250°
WSW locked): SW is windward/wet, NE is leeward/dry. Cells split by annual
precipitation tercile (same method as Tappa 3/4):

| | windward (wet tercile) | leeward (dry tercile) |
|---|---|---|
| Permanent Snow & Ice | 21.9% | 0.3% |
| Subalpine Wet Forest | 26.5% | 0.0% |
| Temperate Forest | 20.9% | 0.0% |
| Woodland / Shrubland | 20.4% | 0.3% |
| Subalpine Dry Scrub | 0.0% | 10.2% |
| Lowland Steppe / Grassland | 0.0% | 73.4% |

Subalpine Dry Scrub is **335.1 km² on the leeward side and 0.0 km² on the
windward side** -- a clean, complete separation, not a mild skew. Lowland
Steppe/Grassland dominates the leeward tercile (73.4%) and is entirely
absent from the windward tercile. This is the strongest confirmation so far
across any stage that the wind-direction decision (`02_tappa2_climate.md`
S1) propagates all the way through to a visually legible map feature: the
rust-colored dry-scrub band is visibly wider on the NE lee side than the SW
windward side in the overview render, matching the physical asymmetry every
upstream stage already quantified numerically.

Qualitative match against the real Southern Alps / Westland / Central Otago
pattern (`00_pre_project_planning.md` validation reference): windward wet
forest -> subalpine scrub/tundra -> permanent snow on the wet side, and
leeward tussock grassland/scrub on the dry side, both recognizable in the
locked scheme. As with every prior stage's validation, this island is
smaller and steeper than South Island, so the contrast is sharper here than
the real reference (same structural point `02_tappa2_climate.md` S5 and
`03_tappa3_snow.md` S5 already made about precipitation ratio and ELA
differential) -- restated rather than tuned away.

## 8. Locked-in parameters

```
pet_constant_mm_per_c        = 58.93   # Holdridge (1967)
biotemperature_clip_c        = [0.0, 30.0]
belt_boundaries_c            = [1.5, 3.0, 6.0, 12.0, 18.0, 24.0]
moisture_tiers                = 3       # world-specific tercile of PET ratio, NOT Holdridge's own bands
belts_split_by_moisture       = ["Boreal", "Cool Temperate"]
belts_not_split                = ["Polar", "Subpolar"]
warm_temperate_folded_into     = "Cool Temperate"
permanent_snow_source          = "Tappa 3 mass-balance mask (balance_mm >= 0)"
```

## 9. Outputs

`run_tappa5_biomes.py` reads `data/processed/climate/` (Tappa 2's outputs,
recomputing Tappa 3's permanent-snow mask directly rather than depending on
a possibly-stale saved copy -- cheap, ~1 s) and writes to
`data/processed/biomes/` (gitignored, regenerates in ~15 s):

| file | contents |
|---|---|
| `biome_id.npy` / `.bin+.hdr+.prj` | int16, 0=ocean, 1-9=land biomes per S5's table |
| `biotemperature_c.npy` / `.bin+.hdr+.prj` | float32, Holdridge biotemperature, °C |
| `pet_ratio.npy` / `.bin+.hdr+.prj` | float32, PET / annual precipitation |
| `moisture_idx.npy` | int8, 0=Wet/1=Moist/2=Dry tercile index |
| `tappa5_biomes_meta.json` | parameters, areas, fragmentation stats, NZ validation numbers |

Same CRS caveat as Tappa 1-4: `.hdr`'s `map info` carries the affine
georeferencing, a `.prj` sidecar carries the PROJ4 string, and QGIS needs
"Fictional World LCC" assigned manually if it loads with an unknown CRS.

Summary of the locked run: land mean biotemperature 7.39°C (range 0.00-12.00),
land mean PET ratio 0.56 (range 0.00-2.18), moisture tercile edges (PET
ratio) 0.1184 / 0.5360, permanent snow 960.2 km² (9.7% of land), 9 land
biomes ranging from 167.2 km² (Subalpine Woodland) to 2,962.2 km² (Lowland
Steppe/Grassland).

## 10. Consequences for the site's visual plan

`00_pre_project_planning.md` calls for a "qualitative palette, nominal data"
biome layer in the interactive map, plus a narrative sequence (terrain ->
climate -> hydrology -> biomes). The overview render
(`assets/05_tappa5_biomes/05_overview.png`) is a candidate for that
narrative sequence's final panel. S6's documented palette limitation should
carry into the site implementation: the click-to-select detail panel is not
optional polish for this layer, it is the accessibility mitigation for a
deliberately natural (not fully CVD-spread) color scheme.

## 11. Open follow-ups (not done in this stage, deliberately left open)

* **No compound official-style Holdridge zone names** (S1) -- if a primary
  copy of Holdridge's diagram is obtained later, the belt x humidity-province
  cells used here could be cross-checked against his ~30-38 traditional
  names; not done this session to avoid fabricating a citation.
* **Subalpine Woodland is the smallest land biome** (167.2 km², 1.7% of
  land) -- worth a visual QA pass in QGIS to confirm it reads as a
  legible transitional band rather than a thin, easy-to-miss line, the same
  kind of check `04_tappa4_hydrology.md` S10 flagged for `lake_mask`.
* **Palette is not CVD-safe** (S6) -- documented and mitigated via the
  planned click-to-select panel, not fixed at the color level. Revisit if
  the final site implementation drops that interaction pattern for this
  layer.
* **No lakes/rivers overlay interaction considered** -- `00_pre_project_planning.md`
  treats rivers as an always-on context layer separate from biome, and this
  stage did not check whether Tappa 4's `lake_mask` cells fall sensibly
  within a single biome class or straddle boundaries oddly; not attempted
  here, out of scope for a biome-classification stage.
* **Tappa 6 (suitability) can now reference `biome_id`** directly as one of
  its inputs alongside slope, water proximity, and solar exposure -- criteria
  still deliberately deferred by the user (`00_pre_project_planning.md`).
