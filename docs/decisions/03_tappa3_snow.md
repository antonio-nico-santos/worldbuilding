# Tappa 3 — Derived climate metrics (snow, seasonality, permanent snow line)

Status: **closed**. Three derived layers computed from Tappa 2's monthly
temperature and precipitation stacks at the existing 120 m working grid:
months-with-snow (precipitation-aware), a seasonality/escursion index, and
a permanent-snow / equilibrium-line-altitude (ELA) proxy built on a minimal
snow mass-balance model, replacing the temperature-only threshold Tappa 2
used provisionally and flagged as wrong.

This document is the decision record for the stage, per the project's
per-stage workflow. As with the two prior documents it is written to work
either as a log of what was decided, or as a reusable recipe.

## 0. What this stage inherits, and what it had to fix

`00_pre_project_planning.md` scoped Tappa 3 as three things: months-with-snow
(count), a seasonality index, and a permanent snow line "useful as a
validation hook against real Southern Alps glaciers." `02_tappa2_climate.md`
S5/S8 then flagged, in its own validation table, that its provisional
permanent-snow number (1.1%, 108 km², "all twelve monthly means below 0 °C")
sits **above** the real glaciological equilibrium line — real Southern Alps
glaciers persist below the elevation where the warmest month's mean crosses
0 °C, because snowfall accumulation there is large enough to outlast the melt
season. A pure temperature threshold has no way to represent that; it has no
notion of how much snow actually falls. That is the one methodological
question this stage had to resolve before the count and index metrics
(which don't have that problem) were worth computing at all.

## 1. Rain/snow partition

Monthly precipitation is split into snow and rain by monthly-mean
temperature, using a **Normal-CDF transition** rather than an all-or-nothing
0 °C cutoff:

```
snow_fraction(T) = Phi((t50_snow_rain_c - T) / sigma_day_c)
```

`t50_snow_rain_c = 1.0 °C` (dry-bulb) is the 50% rain/snow threshold from
Jennings et al. (2018), *Nature Communications* 9:1148 — a ~17,000-station
Northern Hemisphere analysis building on Dai (2008), reporting 95% of
stations' local thresholds between −0.4 and 2.4 °C. 1.0 °C is their
central estimate using screen-level (dry-bulb) air temperature, which is the
variable this project's temperature model actually produces; wet-bulb-based
schemes (closer to the physical mechanism, since evaporative cooling of a
falling snowflake matters more than the ambient air temperature) run
systematically colder, but wiring in a wet-bulb estimate would need a
humidity field this model doesn't have. Using dry-bulb air temperature is a
documented, one-directional simplification, not an oversight.

`sigma_day_c = 3.0 °C` stands in for information a monthly-mean field cannot
supply on its own: within any given month, some days are colder than the
mean and some warmer, and precipitation on the cold days is more likely to
be snow than precipitation on the warm ones. This is modelled by treating
daily temperature within a month as Normal(monthly_mean, sigma_day_c²) —
`sigma_day_c` is *this stage's* one calibrated-by-assumption parameter, in
the same spirit as Tappa 2's `orographic_duty_cycle` (S4d there): a single
number standing in for sub-resolution behaviour the model doesn't otherwise
represent. **Honest limitation**: unlike `orographic_duty_cycle`, this one is
not calibrated against a named external number. A NIWA-sourced
day-to-day standard deviation for West Coast South Island stations
(Hokitika, Franz Josef) was searched for and not found; the value used
bridges a global generic daily-temperature standard deviation (~2 °C,
reported in passing by a degree-day-modelling paper) and a common generic
ice-sheet-model default (~5 °C). Treat 3.0 °C as a documented placeholder
open to recalibration if a better regional figure turns up, not as a
verified constant.

`snow_fraction` is monotonic in T by construction (colder always gives a
fraction at least as high as warmer) and is driven by `sigma_day_c` rather
than a separately-fit curve width — one parameter, one physical
interpretation (how variable are the days within a month), rather than two
that could drift out of sync.

A month is counted for **months-with-snow** if its modelled snow water
equivalent clears `min_snow_month_mm = 5.0` — a floor against counting
months where the Normal-CDF tail is technically nonzero but the actual
snow amount is negligible (the same kind of "effectively zero" cutoff as
Tappa 1's shelf taper, `01_tappa1_terrain.md` S3a).

## 2. Seasonality / escursion index

`seasonality_index = max(monthly T) - min(monthly T)`, per cell. No
methodological choice needed — a direct read of the stack Tappa 2 already
produced. Range on this domain: **7.60 – 10.72 °C**, land mean 8.76 °C,
consistent with Tappa 2's own coastal (8.05 °C) to inland (10.74 °C)
figures (`02_tappa2_climate.md` S3b) — the same physical quantity, computed
independently here from the full monthly stack rather than the two-point
sample Tappa 2 quoted, and landing in the same range. That agreement is a
useful internal consistency check: nothing about the pipeline between
these two stages introduced drift.

## 3. Permanent snow / ELA — a minimal mass-balance model

**Accumulation**: annual sum of monthly modelled snow (§1).

**Ablation**: a positive-degree-day (PDD) model,
`melt_mm = degree_day_factor_mm_per_c_day * sum(monthly PDD)`, with
`degree_day_factor_mm_per_c_day = 4.5` — Hock (2003), *J. Hydrol.* 282,
Table 1: observed snow degree-day factors cluster 3.2–7.6 mm w.e. per °C per
day across the reviewed studies (full reported range 2.5–16.9, generally
lower than ice's 5.5–20.0, since snow's higher albedo means it absorbs less
of the incoming radiation driving the melt). 4.5 sits centrally in the
clustered range.

Monthly PDD is **not** `max(monthly_mean_T, 0) * days_in_month` — that would
silently return zero melt for every month whose mean sits just below
freezing, which is exactly the elevation band the permanent-snow boundary
sits at, i.e. the place this choice matters most. Instead, monthly PDD uses
the same day-to-day Normal(mean, `sigma_day_c`²) assumption as §1, giving a
closed-form expectation:

```
E[max(T_daily, 0)] = mu * Phi(mu / sigma) + sigma * phi(mu / sigma)
```

(the standard rectified-Gaussian / half-truncated-Normal mean — the same
expression used to price a financial call option under a Normal, applied
here to temperature instead of an asset price). A month averaging exactly
0 °C returns ≈1.2 degree-days rather than 0, because roughly half its days
were above freezing even though the mean was not.

**Permanent snow** is redefined as `annual_balance = accumulation - melt ≥ 0`
— accumulation keeping pace with melt over a full repeating annual cycle —
rather than any temperature threshold. This is a **one-year steady-state
proxy** for an equilibrium line altitude under a stationary climatology, not
a real multi-year glacier mass-balance model: no ice dynamics or flow, no
multi-year snowpack carry-over beyond "did this year's accumulation exceed
this year's melt", no meltwater refreezing, no rain-on-snow melt
enhancement, no wind/avalanche redistribution. That is a deliberate scope
cut, not an oversight — see §6.

## 4. Result and the sensitivity check Tappa 2 asked for

| metric | naive (T-only, Tappa 2) | mass-balance (this stage) |
|---|---|---|
| permanent snow area | 108.0 km² | **960.2 km²** |
| fraction of land | 1.09% | **9.69%** |

The mass-balance model gives a permanent-snow area **~9x larger** than the
naive threshold — the correction goes in exactly the direction Tappa 2's own
S5 predicted (the temperature-only definition sits above the real
equilibrium line, i.e. underestimates area / overestimates the boundary's
elevation), now with a mechanism and a number attached rather than just a
caveat.

`02_tappa2_climate.md` S8 flagged `lapse_seasonal_amplitude_c_per_km`
(locked at 0.3) as needing a re-check here, since it "moves summit
seasonality and therefore the permanent-snow area." Re-run at 0.0:

| metric | amplitude = 0.3 (locked) | amplitude = 0.0 (sensitivity) | change |
|---|---|---|---|
| naive permanent snow area | 108.0 km² | 12.5 km² | **−88%** |
| mass-balance permanent snow area | 960.2 km² | 885.8 km² | **−8%** |

The mechanism, checked directly (summit-cell monthly temperatures, both
settings): `lapse_seasonal_amplitude` steepens the lapse rate in February
(SH summer, `lapse_peak_month=2.0`) and shallows it symmetrically six months
later. The naive metric is governed entirely by each cell's **warmest**
month (`all 12 months < 0` reduces to `max(months) < 0`), and the warmest
month is Feb/Jan almost everywhere on this domain — so the naive metric is
maximally exposed to exactly the parameter that steepens Feb specifically.
Turning the amplitude off warms the modelled summer aloft, which is why the
naive area collapses by an order of magnitude. The mass-balance metric
integrates accumulation and melt across all twelve months, so summer's shift
and winter's opposite shift partially cancel in the annual sum — an ~8%
change instead of an ~88% one.

This is worth stating plainly: switching to a mass-balance definition didn't
just fix the naive metric's known elevation bias, it also made the
permanent-snow estimate an order of magnitude less sensitive to a parameter
this project has already admitted (`02_tappa2_climate.md` S3a) it cannot
pin down from first principles. That is a genuine, unplanned side benefit of
this stage's choice of recommended method over the simpler one.

## 5. Validation against Southern Alps, NZ

`00_pre_project_planning.md`'s stated validation reference is South Island,
NZ; `02_tappa2_climate.md` §1 locked the wind at 250° (WSW), making the
south-west-facing flank windward/wet and the north-east side leeward/dry.
Cells are split into windward/leeward proxies by **annual precipitation
tercile** (top third = wet/windward, bottom third = dry/leeward) — a
data-driven proxy using a field the model already produces, rather than a
new geometric wind-facing computation. Within each group, cells are binned
by 100 m elevation band and the ELA is the elevation where the binned mean
annual balance crosses zero (`assets/03_tappa3_snow/03_ela_bands.png`).

| metric | this world | real Southern Alps |
|---|---|---|
| ELA, windward/wet side | 2151 m | ~1600 m (west of the Main Divide) |
| ELA, leeward/dry side | 3288 m | ~2000–2200 m (east) |
| windward:leeward ELA differential | **1136 m** | **400–600 m** |

Source for the real figures: Chinn (1995) and Lamont et al. (1999), via a
2016 Frontiers in Earth Science synthesis — modern Southern Alps ELAs rise
from ~1600 m west of the Main Divide to 2000–2200 m on eastern glaciers, a
400–600 m precipitation-driven differential (rising further, to ~2500 m, at
the drier northern end of the range).

Honest reading: the **direction** is right (wet side lower ELA, dry side
higher, by several hundred metres to over a kilometre) but the
**differential is roughly 2–3x too large**, and both absolute ELAs sit
somewhat higher than their real counterparts. This is the same structural
story Tappa 2's own validation table told (`02_tappa2_climate.md` S5): this
island is smaller and steeper than South Island, so its rain shadow is
sharper and its windward:leeward *precipitation* contrast is already
overstated (5.6:1 modelled vs. 4.0:1 real, per that table) — an
over-steep precipitation gradient feeding directly into an over-steep ELA
gradient is the expected, not surprising, consequence. It is geography
propagating through two stages of the same model, not a new error
introduced here. The leeward ELA (3288 m) sitting only ~307 m below this
DEM's own summit (3595.8 m, `dem_v3_run_meta.json`) is consistent with that
same picture: the dry side of this island barely has enough elevation to
sustain permanent snow at all, mirroring how the real Southern Alps' drier
eastern glaciers are smaller and more marginal than the western ones.

## 6. What this model deliberately does not do

* No ice dynamics or flow — a real glacier's terminus sits well below its
  ELA because ice flows downhill before it melts; this model has no notion
  of flow, so "permanent snow" here means the accumulation zone only, not
  where ice would actually be visible at the surface. Expect this model's
  permanent-snow polygon to sit higher/smaller than a real glacier outline
  at the same ELA.
* No multi-year snowpack memory beyond the annual balance sign — a cell
  with balance exactly at zero is modelled identically to one that just
  barely tips positive, with no notion of firn compaction or multi-year
  accumulation building toward glacier ice.
* No meltwater refreezing (which would reduce effective melt in cold
  snowpacks) and no rain-on-snow melt enhancement (which would increase
  it) — the two omissions push in opposite directions and are not assumed
  to cancel, just both left out for the same reason: neither has the data
  (snowpack temperature profile, rain-during-snow-season timing) this model
  carries.
* No wind redistribution or avalanching, both of which measurably relocate
  snow in real alpine terrain independent of the local energy balance.

None of these are needed for this project's stated purposes (worldbuilding
overview maps, a TTRPG scenario, a portfolio validation hook) — they would
matter for anyone trying to predict this fictional world's actual glacier
outlines at a scale where ice flow is visible.

## 7. Locked-in parameters

```
t50_snow_rain_c                 = 1.0   # Jennings et al. 2018, dry-bulb, ~50th pct threshold
sigma_day_c                     = 3.0   # documented estimate, not independently verified (S1)
degree_day_factor_mm_per_c_day  = 4.5   # Hock 2003 Table 1, snow, clustered range 3.2-7.6
min_snow_month_mm               = 5.0   # "effectively zero" floor for months-with-snow
```

## 8. Outputs

`run_tappa3_snow.py` reads `data/processed/climate/` (Tappa 2's outputs) and
writes into the same directory (gitignored, regenerates in ~18 s):

| file | contents |
|---|---|
| `months_with_snow.npy` / `.bin/.hdr/.prj` | int16, 0-12, count of months clearing `min_snow_month_mm` |
| `seasonality_index_c.npy` / `.bin/.hdr/.prj` | float32, warmest − coldest month mean, °C |
| `snow_accum_mm.npy`, `snow_melt_mm.npy` | float32, annual sums, mm w.e. |
| `mass_balance_mm.npy` | float32, annual accumulation − melt, mm w.e. |
| `permanent_snow_mask.npy` / `.bin/.hdr/.prj` | int16 (0/1), `mass_balance_mm >= 0` |
| `tappa3_snow_meta.json` | parameters, summary stats, sensitivity run, elevation-balance band tables |

Same CRS caveat as Tappa 1 §7 / Tappa 2 §6: `.hdr`'s `map info` carries the
affine georeferencing, a `.prj` sidecar carries the PROJ4 string, and QGIS
needs "Fictional World LCC" assigned manually if it loads with an unknown
CRS.

Summary of the locked run: months-with-snow land mean 4.6 (median 3),
seasonality land range 7.60–10.72 °C (mean 8.76), permanent snow 960 km²
(9.69% of land), windward ELA 2151 m, leeward ELA 3288 m.

## 9. Consequences for the site's visual plan

`00_pre_project_planning.md` already confirmed both "months with snow"
(graduated) and "permanent snow line" (sharp boundary) as switcher options —
both are now backed by a defensible model rather than a placeholder. The
naive-vs-mass-balance comparison map
(`assets/03_tappa3_snow/03_overview.png`, right panel) is worth keeping as a
methods figure for the technical write-up: it makes the ~9x area correction
visually obvious in a way the summary table alone doesn't.

Seasonality's narrow range (7.6–10.7 °C) confirms Tappa 2 §7's point about
temperature choropleths being redundant with a hypsometric map — the same
applies here; an isoline treatment stays the better choice for this
variable too.

## 10. Open follow-ups

* `sigma_day_c = 3.0` is the one parameter in this stage without a named
  external calibration source (§1). If a real NZ West Coast daily-variability
  figure surfaces later, re-run the sensitivity here the same way §4 did for
  `lapse_seasonal_amplitude`.
* The ELA windward:leeward differential (1136 m) is roughly 2-3x the real
  Southern Alps figure (400-600 m) — stated and attributed to this island's
  smaller/steeper geometry (§5), not tuned away, consistent with how Tappa 2
  handled its own oversized windward:leeward precipitation ratio.
* No ice-flow term (§6) — if Tappa 7's urban zoom or any future high-detail
  render needs a plausible glacier *outline* rather than just an
  accumulation-zone mask, that is a materially different (and larger) model
  to build, not an extension of this one.
* `config/parameters.yml`'s `climate:` block should get a `snow:` subsection
  with §7's four parameters, matching how `temperature:`/`precipitation:`
  are already the single source of truth for Tappa 2.
