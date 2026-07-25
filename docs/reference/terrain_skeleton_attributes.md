# Terrain skeleton — attribute reference

Reference for filling in `terrain_ridges` and `terrain_zones` (the hand-authored
control geometry that drives Tappa 1 procedural terrain generation), stored
together in `data/input/terrain_skeleton.gpkg`, exported separately as
`terrain_ridges.geojson` and `terrain_zones.geojson`.

These are starting-point calibrations, not final values — real tuning happens
once the actual generation script is running and you can see rendered output.

## `terrain_ridges` (LineString layer)

| field | meaning |
|---|---|
| `feature_type` | fixed value `"ridge"` — controlled vocabulary, exact match, no free text |
| `peak_elevation_m` | target elevation at the ridge crest |
| `falloff_km` | distance over which elevation decays from the crest toward baseline |
| `name` | human-readable label only — never parsed by the script |

**`peak_elevation_m`** — anchor against the real Southern Alps validation
reference: Aoraki/Mt Cook is 3,724 m, with most of the main crest in the
2,000–3,500 m range. Main spine should land in that band; secondary branches
(remember the range can fork into a Y) should sit noticeably lower — both for
visual distinction and because real secondary ranges are usually lower than
the main uplift.

**`falloff_km`** — matters relative to the domain's half-width (65 km, since
domain width is 130 km). The real Southern Alps are exceptionally steep —
sea level to 3,000+ m within roughly 20–30 km horizontal distance, one of the
fastest elevation gradients on Earth. Setting `falloff_km` close to the full
65 km half-width would let the mountains dominate the entire cross-section
with no room for coastal lowland. Keep it well under that: 15–25 km for a
dramatic, young-uplift feel (consistent with the "dobramentos modernos"
criterion from the original region choice); 30–40 km for something gentler
and older-looking.

## `terrain_zones` (Polygon layer)

| field | meaning |
|---|---|
| `feature_type` | `"plateau"` or `"amplitude_zone"` — controlled vocabulary |
| `target_elevation_m` | flat target elevation — **plateau only**, leave blank for amplitude zones |
| `amplitude_scale` | 0–1 multiplier on fine-noise roughness within the zone |
| `edge_transition_km` | how sharp (small) or soft (large) the boundary blend is |
| `name` | human-readable label only |

**`amplitude_scale`** (0 = completely flat/no texture, 1 = full unmodified
roughness, same as background terrain) — plateaus want this low (0.1–0.2) to
read as clearly flat-topped. Plains (`amplitude_zone`) can sit a bit higher
(0.2–0.35) for a "rolling," less billiard-table-flat feel. Above roughly 0.5
the effect is barely distinguishable from doing nothing, so rarely useful for
either type.

**`edge_transition_km`** — small (2–5 km) gives a sharp, almost
escarpment-like boundary; large (10–20 km) blends gradually into surrounding
terrain. Keep this noticeably smaller than the `falloff_km` of any nearby
ridge, or the zone dissolves into the ridge's own gradient instead of reading
as a distinct feature.

**Plateau vs. amplitude_zone**: a true plateau is flat-topped *and* elevated
— it needs both a low `amplitude_scale` and a fixed `target_elevation_m`.
Plains are the amplitude-only case: low roughness, but no elevation clamp —
the broader trend (e.g. a ridge's falloff) continues underneath, just
smoothed. That's why `target_elevation_m` is left blank for
`amplitude_zone` rows.

## Worked example

**`terrain_ridges`**

| name | feature_type | peak_elevation_m | falloff_km |
|---|---|---|---|
| Main spine | ridge | 2800 | 20 |
| Eastern branch | ridge | 1600 | 15 |
| Coastal foothill | ridge | 700 | 12 |

**`terrain_zones`**

| name | feature_type | target_elevation_m | amplitude_scale | edge_transition_km |
|---|---|---|---|---|
| Central plateau | plateau | 1100 | 0.15 | 5 |
| Southern plains | amplitude_zone | *(blank)* | 0.3 | 12 |
