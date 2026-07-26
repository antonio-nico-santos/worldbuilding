# Pre-project planning — World-building GIS case study

This file supersedes the earlier mid-conversation summary
(`5_pianificazione_pipeline_worldbuilding.md`) — several numbers changed
after it was written (domain size, CRS parameters). Treat this as the
current, authoritative handoff for starting Tappa 1 in a new chat.

## Project identity
- Second case study for a GIS freelance portfolio, dual-purposed as a TTRPG
  worldbuilding scenario (solarpunk setting)
- Rebuild of a prior failed attempt (May 2025): no version control, notes
  scattered in Notion, oversized files from real-world bathymetric data,
  single-axis (temperature-only) biome classification
- Repository: `worldbuilding-gis`, separate from the `gis-portfolio` site
  repo — only final lightweight exports get copied into the site's
  `public/data/`
- Language: English, for this project (code, docs, comments)
- Workflow: this chat = pre-project planning; each pipeline stage gets its
  own chat, closed out with a decision summary in `docs/decisions/`

## Domain and CRS — final, confirmed values

- **Domain**: 160 km (N-S / height) × 130 km (E-W / width), centered on the
  reference point. Area ≈ 20,800 km². (History: started at 150×150,
  briefly 210×100 for wind-alignment reasons, settled here for a better
  on-screen proportion while keeping the area budget roughly constant.)
- **Reference point**: 44°S, 42°E (lon chosen as "the answer to the
  ultimate question of life, the universe, and everything") — confirmed as
  the exact center of the domain, not a corner.
- **CRS**: custom Lambert Conformal Conic, standard parallels recalculated
  for the final 160 km height via the 1/6 rule:
  ```
  +proj=lcc +lat_1=-44.48 +lat_2=-43.52 +lat_0=-44 +lon_0=42 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
  ```
  Defined in QGIS under Settings → Custom Projections.
- **ROI extent** (in the custom CRS, meters):
  ```
  xmin = -65000   xmax = 65000   ymin = -80000   ymax = 80000
  ```
  Created via Processing → "Extent to Layer" (not "Create Grid" — that
  tool's spacing parameter validates against the extent and threw an error
  with the huge-spacing trick originally suggested; avoid it).
- **Validation reference region**: South Island, New Zealand (Southern
  Alps) — chosen for young/active orogeny ("dobramentos modernos"), true
  island status, four distinct seasons, strong biome diversity, and one of
  the best real-world orographic rain-shadow examples on Earth. Data
  sources for later validation: NIWA, or WorldClim/CHELSA.
  **Still open**: exact transect to use (full coast-to-coast vs. partial).

## Pipeline stages (Tappe 0-7)

- **Tappa 0** (parameters) — closed, values above.
- **Tappa 1** (procedural terrain) — next up. Multi-octave Simplex/Perlin
  noise (`opensimplex`/`noise` in Python) + optional hydraulic erosion pass.
  Macro shape (mountain spine, possibly branching into a Y, plateaus,
  regional plains) comes from hand-authored control vectors, not from noise
  alone — see "Terrain skeleton" section below. Coastline is NOT authored;
  it emerges as the zero-elevation contour of the generated DEM. Extra small
  islands can either emerge naturally from background noise amplitude, or
  be authored the same way as the main skeleton if a specific location
  matters narratively.
  **Open**: seed, noise octaves/frequency, erosion iterations — to be
  decided during the Tappa 1 chat itself. **Open**: ridge orientation — N-S
  vs. a diagonal (NW-SE) axis was under discussion; if diagonal, the wind
  vector (Tappa 2) must be set perpendicular to whatever exact angle is
  drawn. Note: total "effect area" of a ridge scales with its drawn
  *length*, not its orientation — a diagonal ridge is only bigger in
  footprint if it's also drawn longer (e.g. stretched corner-to-corner);
  keeping the same length as an N-S version and just rotating it keeps the
  footprint roughly the same.
- **Tappa 2** (climate, monthly) — temperature: latitude baseline → lapse
  rate → continentality correction (`r.grow.distance` + decay function).
  Precipitation: same baseline + orographic/wind effect (windward wetter,
  leeward rain shadow). **Open**: wind direction vector (depends on final
  ridge orientation, see above).
- **Tappa 3** (derived climate metrics) — months-with-snow (count),
  seasonality/escursion index (hottest − coldest month), permanent snow
  line (months where even the warmest stays below 0°C, useful validation
  hook against real Southern Alps glaciers).
- **Tappa 4** (hydrology) — `r.watershed` + `r.stream.extract` on the DEM.
- **Tappa 5** (biomes) — two-axis classification (temperature + precipitation,
  correcting the original project's temperature-only limitation).
- **Tappa 6** (suitability, solarpunk) — deliberately still open; user wants
  to define criteria with their own narrative/scenario context rather than
  the proposed defaults (slope, water proximity, solar exposure via `r.sun`,
  biome).
- **Tappa 7** (urban zoom) — deferred, location depends on Tappa 6. Technique
  already agreed: same noise seed/function as Tappa 1, sampled at finer
  resolution (1-5 m) with additional high-frequency octaves never computed
  at 30 m, plus a hydrological conditioning pass (e.g. `r.fill.dir`) after.

## Visual products (planned, not yet built)

- Hero: hillshade of the full domain — candidate to eventually replace the
  site's placeholder `ContourBackground.astro` with real generated contours.
- Narrative sequence (terrain → climate → hydrology → biomes), mirroring
  Torino's Overview structure.
- Climate display: isolines proposed over choropleth (continuous data,
  consistent with the isochrone pattern already used for Torino) —
  **still not explicitly confirmed by the user**.
- Seasonal snapshots (2-4 representative months) instead of a full 12-month
  slider — reuses the existing radio layer-switcher pattern, no new UI needed.
- Snow: both "months with snow" (graduated) and "permanent snow line"
  (sharp boundary) confirmed as switcher options.
- Interactive map: extension of `InteractiveMap.astro` — biome (qualitative
  palette, nominal data), climate variables (sequential palettes, new hue
  families distinct from Torino's), snow, rivers as always-on context layer.
- 3D model: confirmed, but deliberately sequenced *after* Tappa 6 — draping
  suitability zones on real 3D terrain (via MapLibre's native `raster-dem`
  terrain support) is where 3D adds analytical value, not just aesthetics.
  Requires converting the DEM to a Terrain-RGB-style encoding (e.g.
  `rio-rgbify`).
- Design system: keep constant across both case studies (typography,
  `.site-sheet`, breakout pattern, HTML-only legends, no external tiles,
  radio switcher, click-to-select detail panel). Vary the accent color and
  decorative motif per case study — not yet assigned specific hex values.

## Data authoring workflow (established this session)

- Hand-authored control vectors live in `data/input/`, distinct from
  `data/raw/` (external reference data) and `data/processed/` (generated,
  gitignored).
- `terrain_skeleton.gpkg` (single GeoPackage, since it natively supports
  multiple layers) contains `ridges` (LineString) and `zones` (Polygon)
  tables — kept together because both drive Tappa 1 as one coherent
  "terrain skeleton" concern. `roi` kept as its own file (different role:
  reference frame, not a shape-driving control feature).
- Workflow: create as an editable GeoPackage (safe from crashes, unlike a
  temporary scratch layer) → digitize/edit → export via "Save Features
  As..." to GeoJSON, which is what the actual pipeline script consumes.
  GeoJSON stays split by geometry type (`terrain_ridges.geojson`,
  `terrain_zones.geojson`) — GeoJSON doesn't support multi-layer files
  cleanly the way GeoPackage does.
- **`terrain_ridges` and `terrain_zones` GeoJSON files are now finished**
  (as of this message).

### Attribute schema and calibration (full detail in `docs/terrain_skeleton_attributes.md`)

`terrain_ridges`: `feature_type` ("ridge", controlled vocabulary),
`peak_elevation_m`, `falloff_km`, `name` (human label only).

`terrain_zones`: `feature_type` ("plateau" or "amplitude_zone"),
`target_elevation_m` (plateau only, blank for amplitude zones),
`amplitude_scale` (0-1, low = flat), `edge_transition_km` (small = sharp
edge, large = soft blend), `name`.

Worked example used as a calibration reference: main spine 2800 m peak /
20 km falloff; eastern branch 1600 m / 15 km; coastal foothill 700 m /
12 km; central plateau 1100 m target / 0.15 amplitude / 5 km edge; southern
plains (amplitude_zone) 0.3 amplitude / 12 km edge.

## Folder structure (as built)

```
worldbuilding-gis/
├── README.md
├── .gitignore
├── requirements.txt
├── config/parameters.yml       — needs updating: height_km=160, width_km=130,
│                                  lon=42.0, lat_1=-44.48, lat_2=-43.52
├── scripts/
│   └── setup_structure.py       — one-off bootstrap, moved out of src/
├── src/{terrain,climate,hydrology,biomes,suitability}/
├── notebooks/
├── data/
│   ├── input/                   — terrain_skeleton.gpkg, roi.gpkg,
│   │                               *.geojson exports (this session's work)
│   ├── raw/                     — external reference data (NZ normals)
│   ├── processed/                — generated, gitignored
│   └── exports/                  — final lightweight web-ready assets
├── qgis/                         — .qgz project, .qml styles
└── docs/
    ├── decisions/                — this file lives here
    └── terrain_skeleton_attributes.md
```

## README note (intentional, not an oversight)

The README states both purposes (GIS portfolio *and* TTRPG scenario) with
equal weight, deliberately — aimed partly at a game-design audience that
might find the repo through GitHub. This differs intentionally from the
site's own case-study copy, which should keep the GIS angle primary and
the game use as secondary/complementary context.

## Documentation approach (agreed)

Capture raw parameters and decisions continuously (this file,
`config/parameters.yml`) — cheap, and avoids reconstructing exact values
from memory later, the same lesson the Torino technical appendix already
reflects. Write the final polished site prose (`technical.astro`, etc.)
only once the pipeline and interactive map are stable — writing narrative
against a pipeline still likely to change means rewriting it repeatedly.

## Open items carried into Tappa 1 and beyond

1. Ridge orientation: N-S or diagonal — not finalized
2. Wind direction vector — depends on the above
3. Noise generation parameters (seed, octaves, frequency, erosion) — Tappa 1's job
4. Exact NZ transect for climate validation
5. Isolines vs. choropleth for climate display — proposed, not confirmed
6. Seasonal snapshot months — not chosen
7. Palette hex values — only families proposed so far
8. Tappa 6 suitability criteria — deliberately deferred by the user
9. Tappa 7 exact zoom location — depends on Tappa 6

## Working preferences (for continuity)

Critical, non-validating feedback preferred over agreement. Explain the
*why* behind technical choices, not just the answer. Precision matters —
several corrections in this session were about imprecise or incorrect
claims (e.g. "known CRS" reasoning, "official topo map standard" claim),
and catching those was valued over polite agreement. Language corrections
(Italian and English) requested at the end of messages, ongoing —
recurring patterns worth remembering: capitalizing "I", "didn't + base
verb" (not past participle), question-word-order inversion, singular
nouns after "every".
