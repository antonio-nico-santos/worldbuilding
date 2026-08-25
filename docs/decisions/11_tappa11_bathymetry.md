# Tappa 11 — Bathymetry

Status: **v4 shipped**, supersedes v3's synthetic-noise background and
hardcoded Povo rectangle (SS1-6 are v1's original design record, SS11 is
the v2 revision, SS13 is v3, SS14 is v4 — left intact as history, not
rewritten). v3 introduced the authored-shape workflow (Nico hand-draws
ridges/zones in QGIS, a script turns them into the raster); v4 keeps that
workflow but replaces the *background* texture source — Nico's own DEM
detrended and reused instead of synthetic noise — and replaces Povo
Silencioso's hardcoded bounding-box extent with a derived convex hull, see
SS14. One open flag for Nico carried over unchanged since v1 (ferry
corridor shelter, see SS5 — still no authored shape for the ferry corridor
itself as of v4). Opened from a locked handoff
(`tappa11_bathymetry_prompt.md`, from the Scenario chat) combining a base
spec from Tappa 9 (the provisional B_35k↔C_25k rail bridge corridor) with
two Scenario-chat lore additions (the mainland↔SW-island ferry corridor,
`scenario_reference.md` §18; the Povo Silencioso/NE-archipelago "guardian"
isolation, §19). Placement/raster authorship was this chat's own job, not
sourced from the handoff.

## 0. Inputs consulted before designing anything

- `01_tappa1_terrain.md` for the DEM's exact grid/CRS/export conventions —
  domain `xmin=-65000, xmax=65000, ymin=-80000, ymax=80000`, 30 m
  resolution, shape `(5334, 4334)`, row 0 = north / col 0 = west (same
  convention `write_envi_raw` uses for every other raster in this project).
- `data/processed/dem_v3_final_30m_eroded.npy` (the actual land DEM, not
  just its meta) and `data/processed/geomorphology/lithology_v6.npy` (for
  the volcanic-zone cross-check).
- `data/processed/transport/tappa9_road_network_meta.json` and
  `data/processed/suitability/circulo_candidate_sites.geojson`, to anchor
  the ferry corridor at the project's own already-established
  `Circulo_D_20k<->Circulo_E3_2k` cheapest boat link (10.8757 h / 57.7 km,
  06 row of `index.md`) rather than inventing a new mainland anchor point.

## 1. A finding that changed the scope: there is no existing bathymetry to build on

`01_tappa1_terrain.md` documents `structure(x,y) + noise(x,y)*amplitude -
sea_level_offset_m` for the whole domain, land and ocean alike, plus fluvial
(land-only) erosion — nothing about the seafloor specifically. Checked
directly rather than assumed: binned the existing DEM's ocean depth by
distance-to-coast (90 m sample) and got a **non-monotonic** profile —
mean depth is *shallower* 30-40 km from shore (-357 m) than 5-10 km from
shore (-449 m), the opposite of any real continental margin. That's the
signature of the raw terrain-generation function extended below zero with
no shelf/slope logic, not a designed bathymetry that happens to need three
zones added.

**Consequence**: this stage doesn't carve three zones into an existing
bathymetry — it replaces every true-ocean cell's elevation with a designed
bathymetry (distance-to-coast shelf → slope → basin profile, textured with
noise, same "structure + noise" shape Tappa 1 used for land), then
composites the three lore-driven zones on top the same way Tappa 8
composited its authored lithology zones onto a DEM-native base. Land cells
and **lake** cells are untouched — verified byte-for-byte, see SS6.

## 2. True ocean vs. lake

Reused the project's own established distinction (Tappa 9/10, `04` row of
`index.md`): label the `elevation <= 0` mask's connected components,
border-touching components = true ocean, everything else = a lake. True
ocean = 51.70% of the domain, lakes = 0.66% (matches the existing
`dem<=0`/`lithology_v6==0` mask exactly — checked, 0 mismatched cells).
Only true-ocean cells get a new value; lakes keep their Tappa-1 elevation
exactly.

## 3. Background bathymetry (everywhere not in one of the three zones)

`depth(x,y) = shelf_slope_basin_profile(dist_to_coast_km) + noise(x,y) *
roughness(dist_to_coast_km)`, both terms blended with the project's own
smootherstep function (same one `01_tappa1_terrain.md` §3b uses for zone
blending):

- **Shelf**, 0–18 km from coast: ramps 0 → -140 m.
- **Slope**, 18–55 km: ramps -140 m → -880 m.
- **Basin**, beyond 55 km: holds ~-880 m (domain's max distance-to-coast is
  only 32.8 km, so in practice nothing in this world reaches the pure basin
  plateau — the whole open sea is shelf/slope by this domain's own
  geometry, which is worth knowing for any future open-water design).
- **Roughness**: 8–28 m near the coast, peaking ~70 m mid-slope (canyons/
  scarps concentrate there in real margins), settling to ~35 m in deep
  water.
- **Noise**: a 5-octave value-noise fBm (seed 1104) — not Tappa 1's
  ridged-fBm/domain-warp machinery (that exists to place ridges/coastlines,
  which are already fixed by the DEM; this only needs plausible seafloor
  texture) — normalized to unit std, reused (at different amplitudes) for
  all three zones below rather than drawing a fresh field each time.
- **Hard floor**: `background = min(background, -0.5)` — right at the coast
  the profile approaches 0 and noise alone pushed a few cells to +35.7 m in
  the first run (checked, not caught by construction). Same "independent
  hard cap, defense in depth" pattern Tappa 1's erosion step uses. A second,
  final clamp on the fully-composited output re-asserts this domain-wide.

## 4. Zone 1 — Bridge corridor (B_35k↔C_25k rail crossing, provisional)

Handoff bbox `x∈[-4.7,11.2] km, y∈[42.0,42.5] km`, design footprint buffered
~1-1.5 km, 2 km smootherstep taper back to background.

**What's actually there** (checked directly): the west shore ends at
x=-4.10 km, then a clean ~10.5 km open-water gap, then a small stepping-
stone islet/island cluster (x 6.4–8.35 km, up to 45 m elevation), then a
further ~1.6 km gap, then the east shore begins at x=9.94 km — a real
strait with a natural mid-channel island, close to the handoff's "~11.2 km
primary + ~4.0 km secondary crossing" framing (the islet naturally
partitions the crossing into two hops; I didn't force a precise 11.2+4.0
split — see SS8). **These islets were left untouched** (they're land, not
ocean) — a bridge landing on a real natural island partway across is exactly
how large real sea bridges are built (Chesapeake Bay, Confederation
Bridge), so this is a bonus, not a compromise.

- **Depth**: target -18 m, heavily flattened relief (noise amplitude cut to
  4 m), clamped to [-35, -3] m. Actual result inside the exact handoff bbox:
  **-22.2 to -13.2 m, mean -17.3 m** — comfortably inside "commonly under
  30-50 m," no deep trench anywhere in the zone.
- **Shelter**: a local land-fraction-within-12km index (0=open ocean,
  reference open-ocean points score 0.00-0.09) reads **0.43-0.62** across
  the corridor — a genuinely enclosed bay, not open coast.
- **Substrate**: no offshore geotechnical layer exists in this project, so
  this is inferred, not measured — flagged as such. The immediate shores on
  both sides are firm rock, not soft basin-fill: west shore is greywacke,
  east shore is sedimentary_limestone (checked directly against
  `lithology_v6`). Consolidated bedrock at both landings is a real, if
  indirect, point in favor of a firm nearshore substrate; it says nothing
  about the mid-channel seabed itself.
- **Tidal current**: no tidal model exists in this project (as the handoff
  itself flags) — not addressed here, same open item.
- **Seismic/volcanic**: checked directly — 0 volcanic-class
  (`lithology_v6==4`) cells within a 5 km buffer of the bridge bbox. The
  entire volcanic lithology class in this project is one landmass (the SW
  island, 595.7 km², see `08_tappa8_geomorphology.md`), 50+ km away. Clean
  by a wide margin, not a close call.

## 5. Zone 2 — Ferry corridor (mainland ↔ SW volcanic island)

Anchored at the project's own already-established route:
`Circulo_E3_2k` (7.02, -15.29 km) to `Circulo_D_20k` (-31.12, -56.18 km) —
confirmed the SW island (label 3065 in a connected-components pass) is
literally the entire `volcanic` lithology class (595.7135 km², exact match
to the Tappa 8 meta figure), so "SW volcanic island" = `Circulo_D_20k`'s
island, unambiguously. The straight line between the two Círculos crosses a
single clean **41.74 km** ocean stretch, from (0.82,-21.94) to
(-27.65,-52.46) km, with land only near each settlement itself — used as
the corridor centerline, buffered ±4 km, 3 km taper.

- **Depth**: target -70 m, clamped [-35,-150] m, roughness at half the
  background's amplitude. **Along the actual centerline** (not just the
  buffered swath): -83.4 to -59.0 m, mean -70.2 m, including right up to
  both shore ends (-66 to -74 m there too) — moderate and consistent for
  the full 41.7 km, no near-zero grounding risk anywhere on the sailing
  line itself.
- **Shelter — honest finding, not fully resolved**: a land-fraction-
  within-12km shelter index sampled along the centerline gave **0.53 near
  the mainland end, 0.45 near the island end, but 0.00-0.02 across the
  middle third** (t≈0.29-0.71, roughly x∈[-7,-23], y∈[-31,-48] km) — that
  middle stretch reads as open as this project's own open-ocean reference
  points (0.00-0.09). **This is a real tension with §18's "sheltered...
  not exposed to full open-sea swell" framing, not something bathymetry
  design can fix** — depth and wave/current exposure are different
  physical properties, and no amount of seabed shaping manufactures
  shelter from open coastline geometry. A chain of small islands does
  exist further west (labels 2140/2839, 10.7/42.1 km², roughly 20-45 km
  off the direct line) that a real ferry captain might route through for
  genuine lee shelter, at the cost of a materially longer, more indirect
  crossing — **not built here**, since re-routing the corridor is a bigger
  decision than this stage's own scope and wasn't requested. Recorded as
  an open item for Nico/the Scenario chat (SS8) rather than silently
  asserted as resolved.
- **Consistency**: kept the corridor's depth in one moderate, non-extreme
  band end-to-end (nothing near the existing -698 m trough that sat on this
  exact line in the unmodified DEM) — satisfies "nothing here should read
  as more hazardous than the general sea it's an exception to" as far as
  depth goes, independent of the shelter finding above.

## 6. Zone 3 — Povo Silencioso (NE archipelago)

**Checked directly rather than hand-picking a few islands**: the "NE
archipelago" is not a handful of named islands — a connected-components
pass over the NE bounding box (`x∈[27,65] km, y∈[58,80] km`) found **300+
separate land components**, from a few 15-30 km² islands down to
single-cell skerries, none carrying volcanic lithology (all `basin_fill`/
`greywacke`). Footprint = true-ocean cells within 20 km of any land cell in
that box (211,885 land cells / 190.7 km² of islets total), full target
reached within 4 km of any islet (steep near-shore drop-off), tapering back
to background by 20 km out.

- **Depth**: target -1200 m, clamped [-150,-1900] m, roughness amplified to
  160 m (vs. the background's 8-70 m). Result within 4 km of any islet:
  **-1755.9 to ~0 m** (the ~0 end is the immediate shoreline of the
  islets themselves, expected), **mean -1197 m**. This is a deliberate new
  domain extreme — the input DEM's deepest point anywhere was -945.6 m;
  this zone alone reaches -1755.9 m, clamped no deeper than -1900 m.
  **Flagged explicitly, not hidden**: nowhere else in this world is
  authored this deep.
- **Steepness**: median local seafloor slope within 4 km of an islet is
  **0.194** (≈11°) vs. **0.0055** for both the bridge and ferry zones —
  roughly **35x steeper**, a real, measured drop-off signature, not just a
  deeper number.
- **Design intent, not a resolution of the guardian question**: this makes
  the water itself a sufficient, mundane explanation for the archipelago's
  isolation (steep, deep, genuinely dangerous to navigate) without asserting
  or denying anything about "the guardian" — exactly the ambiguity §19
  asks to be preserved. A boat could still physically attempt this crossing;
  it would just be a bad idea on the merits of the seafloor alone.

## 7. Verification (Task-list style, all checked directly against the output array)

- No NaN/Inf anywhere in the output (assert, passed).
- Land cells byte-identical to `dem_v3_final_30m_eroded.npy` (max abs diff
  = 0.0, assert, passed).
- Lake cells byte-identical to the same input (assert, passed) — this stage
  only ever touches true-ocean cells.
- Every true-ocean cell ≤ -0.5 m (assert, passed) — no cell reads as dry
  land or exactly at sea level.
- Bridge zone: 0 volcanic-lithology cells within a 5 km buffer (checked).
- Bridge zone depth -22.2..-13.2 m inside the exact handoff bbox (checked).
- Ferry centerline depth -83.4..-59.0 m end-to-end, including both shore
  approaches (checked) — no unexpected shallow points.
- Povo Silencioso median slope 35x the bridge/ferry zones' (checked).
- A full-domain preview render (`bathymetry_v1_overview.png`) visually
  confirms all three zones read distinctly from their surroundings — pale/
  shallow strait at the bridge, a visibly deep/rough patch specifically
  around the NE islet cluster, a lighter mid-depth band along the ferry
  line through otherwise darker water.

## 8. Explicitly not done here (open items)

- The bridge corridor's "primary ~11.2 km + secondary ~4.0 km" framing was
  **not** forced into two separately-routed polylines through the natural
  mid-channel islet — the whole buffered bbox was treated as one shallow-
  shelf footprint, which necessarily covers both crossings without
  requiring a speculative reconstruction of exactly how a not-yet-designed
  rail line would split around the islet. If an actual rail alignment is
  authored later, this raster already supports it either way.
- The ferry corridor's mid-channel shelter gap (SS5) is flagged, not fixed.
  Two honest paths forward: accept that this route's real-world safety
  comes from operational factors (scheduled sailings, larger vessels,
  experienced crews — a legitimate real-world pattern, not a cop-out) rather
  than pure geography; or re-route the corridor through the western island
  chain for genuine shelter at the cost of a longer crossing. Not decided
  here.
- No tidal/current model exists in this project at all (bridge's own
  requirement #4) — unaddressed, same status the handoff itself flags.
- No offshore geotechnical/substrate raster exists — the bridge's firm-
  substrate claim rests on adjacent-shore lithology only, disclosed as an
  inference, not a measurement.
- The background shelf/slope/basin profile is a reasonable default for the
  ~83% of true-ocean cells outside all three zones, not an authored,
  per-region bathymetry the way the three requested zones are — in scope
  for this handoff, but worth knowing before treating, say, the open sea
  southwest of the domain as narratively significant without a closer look.

## 9. Cross-reference

A condensed handoff summary was written to the claude.ai project as
`claude/tappa11_bathymetry_report.md` (same role as `tappa8_new_rocks_
vertice_handoff.md` played for Tappa 8's rock-property handoff). The
project's big `index.md` navigation doc was **deliberately not rewritten**
this session — it's a ~40K-word, carefully-maintained cross-chat history,
and full-content-replace-only editing (no in-place patch available from
this session) made a transcription-safe edit impractical here. Whoever next
touches `index.md` should add: a status-table row (`11 | Bathymetry | v1
shipped, one open flag (ferry corridor shelter) | ...`) and a short
per-document-detail paragraph pointing at this file and the project report
above.

## 11. v2 revision — Nico's review and the authored-shape redesign

Nico reviewed v1 in QGIS and reported three concrete, specific problems, plus
a direct question:

1. "Despite the area around Povo Silencioso, the water is too shallow
   everywhere and seems flat. Mix some deepen areas around the island in a
   natural way, specially in the areas facing the open sea."
2. "There is a clearly stripe on south, connecting the island to the
   Mainland. It should be a more organic and real shape."
3. "Even in the part that water should be more shallow, in the North area,
   they could be around -20 to -40m and it is almost like -10 or even more
   shallow." — i.e. the bridge corridor.
4. "Could authoral shapes as basis help to draw it, in a process more
   similar to Tappa 1 process?"

**Answer to (4), stated plainly rather than folded into the fixes below**:
yes, and it's not a stylistic preference — it's the actual root cause of (1)
and (2). v1 used pure distance-field logic (distance-to-coast,
distance-to-nearest-islet, distance-to-a-straight-line) for every zone. A
distance field to a straight line is, by construction, a straight-line
buffer — that's not "a stripe-like result," it *is* a stripe, exactly and
only. Likewise a single radial falloff from "nearest islet" with one target
depth has no mechanism to distinguish an open-sea-facing arc from a
sheltered one — nothing in the math *could* produce that asymmetry, no
matter how the noise was tuned. These weren't tuning mistakes; the chosen
technique structurally couldn't produce what was asked for. Tappa 1's actual
toolkit (`src/terrain/skeleton.py`'s `RidgeField`, `src/terrain/noise.py`'s
`domain_warp` + `ridged_fbm`) exists specifically because Tappa 1 hit the
same wall authoring coastlines and ridges, and the fix is the same fix: (a)
domain-warp the query coordinates so straight buffers become organic
channels, (b) use real hand-placed lines/points as the shape's skeleton
instead of "nearest feature of a whole class," (c) drive asymmetry off a
real directional quantity (here, distance-to-mainland) instead of a
symmetric radial one. Item (3) wasn't a technique problem — v1's bridge zone
used the same clamp-and-blend approach v2 does, just retargeted at the wrong
depth — but doing all three fixes as one pass made sense rather than three
separate small patches.

### 11.1 Fix for (3): bridge corridor retargeted

Target changed -18 m → -30 m, clamp changed `[-35,-3]` → `[-40,-20]` m.
Geometry (bbox, taper) is unchanged from SS4, so the substrate/volcanic/
shelter findings there still hold as-is — only depth moved. Verified by
recomputing the bbox's true-ocean depth distribution directly against the
new output (not reused from the generation script's own asserts):

- v1: min -24.5 m, max -11.1 m, mean -17.7 m — only 24.9% of cells actually
  fell inside [-20,-40] m.
- v2: min -36.5 m, max -23.1 m, mean -29.7 m — **100%** of cells inside
  [-20,-40] m.

### 11.2 Fix for (2): ferry corridor as an authored, warped polyline

Replaced the 2-point straight line + constant ±4 km buffer with:

- A **4-waypoint centerline** (mainland shore → two interior bends → island
  shore) instead of 2 points. The interior bends were placed a few km off
  the original straight line, in the direction v1's own output showed as
  the strait's genuinely deeper water (checked directly in v1: the deepest
  existing water along the straight line sat slightly west of its
  midpoint) — so the curve tracks a plausible real channel, not an
  arbitrary wiggle.
- **Domain-warped query coordinates** before the distance-to-centerline
  lookup (`domain_warp`, seed 4402, 6 km wavelength, 1.8 km amplitude, 3
  octaves) — this is the actual mechanism that turns the buffer's edges
  from two parallel curves into an organic, textured coastline-like
  boundary. Same function, same technique Tappa 1 used for its coastlines.
- A **breathing halfwidth** (2.85–4.65 km via a `Simplex2D` field along the
  route's arc length) instead of a constant 4 km.

Verified visually (`bathymetry_v2_ferry_zoom.png`, v1 vs v2 side by side at
native resolution) — v1 is a clean parallel-edged band; v2 has a genuinely
irregular, branching boundary with internal texture, not just a wider or
narrower version of the same shape. Numerically, sampling perpendicular
transects at each of the 3 route segments (width measured at the -35 m
threshold) gave 10.6–13.3 km, a ~25% variation along the route that a
constant-width buffer cannot produce by construction.

### 11.3 Fix for (1): Povo Silencioso — asymmetry + hand-placed trenches

Three additive pieces, all reusing Tappa 1's real math rather than ad hoc
approximations of it:

- **Mainland-distance asymmetry** (the actual fix for "especially in the
  areas facing the open sea"): target depth is no longer one uniform
  number. `povo_target = -320 m + (-1400 - -320) * smootherstep(dist_to_
  mainland_km, 2, 15)` — shallower near the mainland-facing side of the
  archipelago, deep on the far side. `dist_to_mainland` is computed against
  the single largest landmass component specifically (not "any land"),
  which is the correct referent for "facing the open sea" vs. "facing the
  mainland."
- **Two hand-placed trench lines** — `RidgeField.contribution`'s exact
  Gaussian-decay-from-a-polyline math (`src/terrain/skeleton.py`), reused
  directly with a negative peak instead of reimplemented: Trench A
  ("east approach," peak -750 m, hugs the domain's open-ocean edge) and
  Trench B ("inter-island channel," peak -550 m, threads the open water
  between the two largest northern islands). These are real, structured
  deep features at chosen locations, not a bigger noise number.
- **A coarse `ridged_fbm` band** (9 km wavelength, seed 7301, Tappa 1's own
  noise family) layered under the existing fine-texture band, for genuine
  seamount/basin-scale structure rather than only high-frequency grain.

Verified numerically (recomputed from scratch against the delivered
raster, splitting the Povo bbox at the median distance-to-mainland found
inside it, 6.24 km):

- v1: sheltered-half mean -1194.5 m, open-half mean -1194.2 m — **no
  asymmetry at all** (a 0.3 m difference, i.e. noise, not signal). Whole-box
  std 173.4 m.
- v2: sheltered-half mean -146.6 m, open-half mean -585.9 m — a **439 m**
  real difference. Whole-box std 495.6 m (2.9x v1's), min -1900.0 m (the
  clamp floor, actually reached, vs. v1's -1755.9 m natural max) down to
  -0.5 m right at an islet's shore (the fringe taper — expected, islands
  need beaches).

Verified visually (`bathymetry_v2_povo_zoom.png`) — v1 reads as a uniform
teal wash; v2 shows a clear shallow shelf hugging the mainland-facing/SW
side and dark, structured deep water (down to the -1900 m clamp) on the
NE/open-sea-facing side, with visible trench-like features rather than
smooth radial fade.

### 11.4 What did *not* change

- The guardian-legend ambiguity is preserved exactly as before — the deep
  water is still a sufficient mundane explanation for isolation without
  asserting or denying anything about it. If anything, v2 makes this
  *stronger*: real trenches and a genuine open-sea/sheltered-side contrast
  read as more deliberately hazardous than v1's uniform depth did, without
  adding a single narrative claim.
- The background shelf/slope/basin profile (SS3) — Nico didn't flag it, and
  it wasn't touched.
- The ferry corridor's mid-channel shelter gap (SS5's open item) — **not
  re-verified for v2's shifted centerline.** The waypoint bends are only a
  few km off the original straight line, well inside the stretch SS5
  measured as open (t≈0.29–0.71 along a 41.7 km line), so the finding almost
  certainly still holds, but this was not re-run against the new geometry
  and shouldn't be treated as re-confirmed — flagging the gap rather than
  quietly assuming it's fine.
- Land and lake cells: still byte-identical to the input DEM (re-asserted
  against the v2 output specifically, not just inherited from v1's checks).

### 11.5 A generation-time fix worth recording (not a design decision, but a real constraint hit)

The first v2 run was killed by the sandbox's own OOM killer partway through
the coarse `ridged_fbm` call. Cause, confirmed via `dmesg`: `Simplex2D.
noise2` (`src/terrain/noise.py`) force-casts to float64 and holds roughly
20 full-grid temporary arrays alive at once inside a single call; at this
domain's size (5,334×4,334 = 23.1M cells, ~185 MB per float64 array) that
peaks well past this sandbox's ~6 GB memory ceiling when called on the
whole grid at once. Fixed by adding a `chunked_apply` helper in the
generation script that calls `ridged_fbm` / `domain_warp` / raw `.noise2`
on ~250-row horizontal bands instead of the full grid — identical math
(each cell is still evaluated at its true (x, y)), bounded transient
memory — plus explicit `del` of intermediate zone arrays once each zone's
blend is folded into `result`, since Python doesn't free module-level names
until reassigned. Mentioned here because it's a real constraint of this
generation environment, not this project's domain, and the same pattern
(row-band `chunked_apply`) will apply to any future stage that calls
`ridged_fbm`/`domain_warp` across this full domain size.

## 13. v3 — authored-shape workflow (Nico's hand-drawn ridges/zones)

Nico's request: draw the shapes themselves in QGIS instead of me hand-coding
zone geometry, mirroring Tappa 1's own authored-skeleton pipeline
(`src/terrain/skeleton.py`'s `RidgeField`/`ZoneField`, `src/terrain/
generate.py`'s compositing formula). Inputs: `data/input/
bathymetry_ridges.geojson` (4 lines) and `bathymetry_zones.geojson` (2
polygons), authored in the device repo's QGIS project.

**Loader deviations from the Tappa 1 reference (`src/terrain/skeleton.py`),
by design, both in `run_tappa11_bathymetry_v3.py` directly rather than
touching the shared module:**
- `build_ridge_fields()` requires the literal `feature_type=='ridge'`; the
  custom loader here also accepts `'creek'`, `'valley'`, `'trench'` as
  equivalent labels (a negative `peak_elevation_m` already encodes the
  semantic difference), since Nico's "Guardian creek" trench was authored
  with a more descriptive label than the code's literal discriminator.
- `build_ridge_fields()` never reads `shelf_multiplier` from GeoJSON
  properties at all (only via an external `{name: value}` dict at call
  time). The custom loader reads the authored `shelf_muliplier` property
  (typo, as drawn) directly per feature.

**Background profile changes (from a chat discussion on realism, before
generation):**
- `SHELF_BREAK_KM` 18→12 km (Nico's explicit request).
- Coast→shelf-break segment switched from `smootherstep` to an ease-out
  curve `1-(1-t)^2`. Smootherstep has zero 1st *and* 2nd derivative at
  t=0, so it was giving unrealistic -0.7..-7.2 m depths 1-2.3 km offshore
  — checked directly against Nico's own QGIS observation of the same
  numbers. The ease-out curve gives ~-22 m at 1 km, ~-43 m at 2 km instead.

**Compositing order** (`structure` = deterministic shape, noise added once
at the very end): background profile → ferry corridor (still v2's
hardcoded domain-warped waypoint corridor — no authored equivalent yet) →
Povo Silencioso mainland/open-sea asymmetry (still hardcoded — no authored
equivalent yet) → authored ridges, **summed** (not `np.maximum`, unlike
Tappa 1's mountain-building ridges — max() is for unidirectional
peak-competition and makes a negative-peak trench a no-op against a
less-negative background; v2 already established sum as the bathymetry
convention for this reason) → authored zones, exact Tappa 1 plateau formula
(`structure = structure*(1-w) + target_elevation_m*w`), applied
largest-area-first so "Underground Lake" composites after, and wins inside,
"Bridge Base" (checked: 169/169 of its vertices sit inside Bridge Base —
genuinely nested, the "undersea lake" technique discussed in chat) → noise
(fine + gated coarse, scaled by distance-based roughness × each zone's own
`amplitude_scale`) → defense-in-depth clamp (true_ocean cells never exceed
-0.5 m, same as v1/v2 — a ridge that crests above sea level gets flattened
to -0.5 m rather than becoming new land).

Guardian creek **replaces** v2's hand-placed `trench_a`/`trench_b` (both
lived in the same NE-archipelago region) rather than adding to them — the
whole point of authoring was for Nico to place this by hand.

**Two real bugs found and fixed while generating (not design surprises —
both confirmed by sampling the actual composited grid, not just the
pointwise validator):**

1. *Coastal "bathtub ring."* First run: 368,507 true_ocean cells (331.6 km²,
   661 separate patches strung along nearly every coastline in the domain)
   hit the -0.5 m clamp — nowhere near the four authored ridges. Cause:
   roughness had a flat 8 m floor right at the coastline, oversized against
   a background that's only -2..-20 m deep in the first km (v2 likely had
   the same issue, silently absorbed by its own early clamp with no
   reporting). Fixed in two steps: (a) ramp roughness's near-shore term
   0→8 m over the first 1.5 km instead of a flat floor (331.6→35.8 km²);
   (b) gate the coarse (9 km wavelength, basin-scale) noise band out of the
   surf zone entirely, 0→full weight over 0.5-2.5 km offshore — a slow
   regional "high" in that band had no business breaching sea level at the
   coast (35.8→25.6 km²). Final clamp footprint: 34.9 km² (post all fixes,
   see below), 95.7% of it within 500 m of shore, mean distance 94 m — a
   thin, physically-plausible waterline fringe, not a broad artifact.
2. *Povo box-membership leak.* Sampling the real raster along "North
   camling ridge 1" found its crest at ~-287 m instead of the intended
   ~-20..-30 m reef depth. Cause: Povo's falloff (`dist_to_ne_land`)
   measures distance to the nearest land cell *inside* the NE-archipelago
   box, with no check on whether the query point itself is near the box —
   this ridge sits ~20 km outside the box but ~7 km from a land pixel just
   inside its edge, so `w_povo` reached ~0.92 there, overriding the plain
   shelf background with Povo's -320 m sheltered target almost entirely.
   This is a latent v2 bug too (identical logic, just never exercised at a
   point anyone checked). Fixed by gating on distance to the box itself
   (0→full weight killed beyond 3 km past the box edge), leaving in-box
   behavior untouched. Povo's own core-effect footprint dropped from a
   (leaky) 375.1 km² to a (correct) 188.4 km² as a result.

**Final per-ridge crest depth, sampled from the actual generated raster
along each authored line** (not the pointwise approximation from the
pre-generation validator):

| ridge | `peak_elevation_m` | actual crest depth (min / max / mean) | clamp hits |
|---|---|---|---|
| Guardian creek | -300 | -1536.1 / -322.1 / -598.7 m | 0/252 |
| North camling ridge 1 | 20 | -27.6 / -0.5 / -14.5 m | 1/25 |
| North Caming ridge 2 | 23 | -37.2 / -0.5 / -16.8 m | 6/43 |
| South Barrier | 20 | -92.4 / -7.1 / -51.0 m | 0/36 |

Guardian creek's much deeper range than its own -300 m peak alone would
suggest is expected and correct once the Povo leak is fixed properly (not
removed): it sits genuinely *inside* the NE box, so it inherits Povo's own
-320..-1400 m mainland/open-sea asymmetric target underneath it, plus its
own -300 m on top. North camling ridge 1/2 both still touch the clamp at
their shallowest points (1/25 and 6/43 vertices) — a small residual, not
zeroed, since Nico's own peak-elevation edits already cut per-vertex
surfacing risk from 44%/28% to these single-digit counts; further
lowering `peak_elevation_m` on either would remove it entirely if wanted.

Zone centroids sampled from the final raster: Bridge Base -36.7 m (target
-27 m, difference is roughness noise at `amplitude_scale=0.8`), Underground
Lake -70.6 m (target -75 m) — both in the expected range around their
authored targets.

**Verification:** same numeric checks as v1/v2 (no NaN/inf, land/lake
pixels bit-identical to the DEM, true_ocean max < 0), plus the per-ridge
sampling and clamp-footprint analysis above (independent of the generation
script's own asserts — done by reloading the delivered raster and
resampling, not trusting in-script logging alone). Visual: full-domain
overview plus 4 zoomed panels (Bridge Base/Underground Lake,
North camling ridges 1/2 — now visibly two distinct crests rather than
invisible at this ridge's ~250-350 m footprint, Guardian creek — visibly a
continuous deep channel distinct from the surrounding shelf/basin once
plotted on a -400..0 m scale, South Barrier — a subtle but visible shallow
band).

Files: `run_tappa11_bathymetry_v3.py` (script, committed to repo root);
`data/processed/bathymetry/bathymetry_v3_30m.{bin.gz,hdr,prj}` +
`bathymetry_v3_meta.json` (committed, same convention as v1/v2 — gunzip the
`.bin.gz` in place before opening in QGIS); `bathymetry_v3_30m.npy` (full
float32 precision, sandbox/chat-download only, same as v1/v2's pattern);
`docs/decisions/assets/11_tappa11_bathymetry/bathymetry_v3_{overview,
zooms}.png` (committed).

## 14. v4 — DEM's own texture, and a derived Povo hull

Nico's QGIS review of v3 flagged two things: the synthetic background read
as shallow and flat, and a visible square artifact around Povo Silencioso.

**The square artifact**, checked directly (rendered a crop straddling
`NE_X0=27000`): a hard vertical seam exactly at the box's edge, made worse
by v3's own Povo-leak fix, which added a hard-ish 3km falloff past the
rectangle — smoothing a background-to-Povo-target jump that can exceed
1000m over just 3km is still effectively a cliff. Fixed by dropping the
rectangle as the thing that gates the raster (kept only as a seed filter
for "which land pixels count as Povo") and deriving a convex hull from the
archipelago's actual land pixels instead — 11 vertices, 655.3 km², no
straight edges — blended via the same edge_transition mechanism an
authored zone uses (5km band straddling the hull boundary). Confirmed by
re-rendering the same crop: seam gone, Guardian creek's trench now follows
the real channel between islands with no discontinuity.

**The flat/shallow background.** Bin the DEM's own raw ocean depth by
distance-to-coast and it deepens fast, peaks around 5-8km (-452m mean),
then gets *shallower* again out to 35km (-348m) — non-monotonic, the
opposite of a real shelf/slope/basin. That's Tappa 1's land-generation
noise field extended below sea level with no oceanic constraint, not
designed bathymetry — confirmed not safe to build authored shapes directly
on top of (Nico's own proposed "option A"). But its *local* texture (the
part correlated with nearby land, not the bad radial trend) is exactly
what v3's synthetic noise structurally couldn't provide, being a pure
function of distance-to-coast with no spatial correlation to anything.

Fix (Nico's own proposed "option B", implemented): detrend the DEM's ocean
signal (fine-binned median by distance-to-coast, Gaussian-smoothed,
removed regardless of its shape) to isolate the residual texture; replace
the radial component with the same shelf/slope/basin curve v3 used;
recombine. First attempt applied the residual at literal full strength
(Nico's initial choice) — broke immediately: 31% of all true_ocean
(3396 km², out to 20km offshore) clamped flat, because the residual is
land-terrain-scale noise (std 66-140m near shore, p95 up to +224m within
3-8km) and no uniform multiplier fits that inside a ~20-140m-deep shelf.
Swept a headroom-based scale (full strength once the deterministic
background/zone/ridge structure reaches some threshold depth, ramped down
in shallower water) against the actual composited grid rather than
guessing: 150m threshold → 427 km² still clamped; 250m → 14.5 km² (already
better than v3's own 34.9 km²) at ~17% mean residual strength; 350m → zero
clamp but only ~9% mean strength, texture barely survives. Shipped at
250m — smallest clamp footprint of any version so far, while keeping
meaningfully more DEM character than v3's tuned synthetic bands ever gave.
Actual delivered clamp: 23.5 km² (slightly above the isolated sweep,
since ridges/zones/ferry/povo locally reduce headroom in a few places).

**Per-ridge crest depth, re-sampled from the v4 raster** (all four still
land-free, no Povo leak, no clamp hits on any ridge):

| ridge | `peak_elevation_m` | actual crest depth (min / max / mean) |
|---|---|---|
| Guardian creek | -300 | -1100.5 / -295.7 / -666.1 m |
| North camling ridge 1 | 20 | -26.0 / -2.3 / -17.3 m |
| North Caming ridge 2 | 23 | -46.7 / -1.2 / -24.8 m |
| South Barrier | 20 | -105.2 / -12.1 / -58.5 m |

Zone centroids: Bridge Base -39.5m (target -27m), Underground Lake -73.3m
(target -75m) — both shifted somewhat from v3's numbers since the
background under them is now DEM-textured rather than synthetic-noise, an
expected consequence of the texture-source swap, not a regression.

Ferry corridor and Povo's mainland/open-sea asymmetric target formula are
both still unchanged/hardcoded (same open item as v3, see SS5 and the
Outputs section).

Files: `run_tappa11_bathymetry_v4.py` (script);
`data/processed/bathymetry/bathymetry_v4_30m.{bin.gz,hdr,prj}` +
`bathymetry_v4_meta.json`; `bathymetry_v4_30m.npy` (full precision,
sandbox/chat-download only, same pattern as v1-v3);
`docs/decisions/assets/11_tappa11_bathymetry/bathymetry_v4_{overview,
zooms,box_check}.png`.

## 12. Outputs

**v2 is the current deliverable.** v1's files are kept in place as the
historical record (same "supersede, don't delete" convention
`01_tappa1_terrain.md` uses across its own v1→v2→v3) — anything reading
`bathymetry_v1_*` going forward should treat it as superseded, not current.

**Committed to this repo** (device-bridge size limits meant only these made
it back automatically this session — see note below):
- `run_tappa11_bathymetry.py` (v1) and `run_tappa11_bathymetry_v2.py` (v2,
  current) — both kept; v2 imports and reuses Tappa 1's actual
  `src/terrain/noise.py` / `src/terrain/skeleton.py` machinery rather than
  v1's own simplified distance-field logic (see SS11). Both fully
  deterministic (fixed seeds) — re-running either reproduces its output
  bit-for-bit, including files not committed below.
- `docs/decisions/11_tappa11_bathymetry.md` — this file.
- `docs/decisions/assets/11_tappa11_bathymetry/bathymetry_v1_overview.png`
  (superseded) and `bathymetry_v2_overview.png` (current) — full-domain
  preview renders with the three zones annotated.
- `docs/decisions/assets/11_tappa11_bathymetry/bathymetry_v2_ferry_zoom.png`,
  `bathymetry_v2_povo_zoom.png`, `bathymetry_v2_bridge_hist.png` — v1-vs-v2
  comparison renders backing SS11's verification claims directly.
- `data/processed/bathymetry/bathymetry_v1_meta.json` and
  `bathymetry_v2_meta.json` — grid, method, all zone parameters,
  verification numbers, and (v2 only) an explicit `changes_from_v1` list.
- `data/processed/bathymetry/bathymetry_v{1,2}_30m.hdr` + `.prj` — ENVI
  header + CRS sidecar, same convention as `dem_v3_final_30m_eroded`.
- `data/processed/bathymetry/bathymetry_v{1,2}_30m.bin.gz` — the actual
  int16 raster data, gzip-compressed (46 MB → ~17 MB) to clear the device
  bridge's per-file transfer cap. **Run `gunzip` on it in place** (so it
  becomes e.g. `bathymetry_v2_30m.bin` next to the `.hdr`/`.prj` above)
  before opening in QGIS — the `.hdr`'s filename convention expects the
  uncompressed name.

**Not committed** — the device bridge used this session caps file transfers
well under this data's actual size (a 20-30 MB limit vs. this domain's
23.1M-cell grid; int16 compresses under that limit, float32 does not, even
gzipped). These exist in the generation sandbox and were delivered as chat
downloads instead, or can be regenerated in one command:
- `data/processed/bathymetry/bathymetry_v{1,2}_30m.npy` (float32, full
  precision, ~92 MB each) — same "full precision available on request,
  int16 is the standard deliverable" pattern `01_tappa1_terrain.md` §7
  already established for the DEM itself, not a new gap this stage
  introduced.
- `data/processed/bathymetry/bathymetry_depth_only_30m.npy` (v1 only,
  true-ocean cells only, NaN elsewhere — a convenience layer for anyone who
  wants depths without re-deriving the ocean mask; regenerate via either
  script, or derive it from a committed raster with the same
  border-connected-component method §2 describes).
- `verify_v2.py`, `render_v2_preview.py` — the independent verification and
  rendering scripts behind SS11's numbers and images (independent in the
  sense of recomputing `true_ocean`/`lake` and all zone boxes from scratch
  against the delivered rasters, rather than trusting the generation
  script's own internal asserts).
