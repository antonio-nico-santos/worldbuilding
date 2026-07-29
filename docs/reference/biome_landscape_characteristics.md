# Biome landscape characteristics — visualization reference

Reference for what each of Tappa 5's 9 land biomes (plus Permanent Snow & Ice)
should actually *look like* on the ground: vegetation structure, typical plant
and tree sizes, common species, and the temperature/precipitation character
that produces them. Intended to help with concept art, prose description, or
any other creative pass that needs more than a biome name and a hex color.

**How to read this doc.** Each section has two parts. The first is this
world's own numbers, computed directly from the Tappa 5 output
(`data/processed/biomes/biome_id.npy`, `biotemperature_c.npy`) and Tappa 2's
climate stack (`data/processed/climate/`) — these are exact, not estimated.
The second is a real-world analog, almost entirely South Island New Zealand,
chosen because it's the same analog the whole project has validated against
since Tappa 1 (a young, steep, maritime mountain range at cool-temperate
latitude with a strong windward/leeward precipitation contrast). **The
species and figures below are real-world facts about New Zealand ecology,
not in-world botanical canon.** Use them the way you'd use a reference
photo — for scale, structure, and plausibility — not as a claim that this
world's flora is literally identical to Earth's. If your world has different
plant lineages, keep the sizes/structure/density and swap the species names.

Belt and moisture-tercile terminology (Polar, Subpolar, Boreal, Cool
Temperate; Wet/Moist/Dry) follows `src/biomes/holdridge.py` and
`src/biomes/world_biomes.py`; see `docs/decisions/05_tappa5_biomes.md` for
the classification derivation. "Biotemperature" is Holdridge's own quantity
(monthly mean temperature, clipped to 0–30°C, then averaged) — it runs
noticeably cooler than a plain annual mean because it zeroes out any
below-freezing months rather than letting them pull the average down.

---

## Permanent Snow & Ice

| | this world |
|---|---|
| Elevation | 1,981–3,562 m (mean 2,782 m) |
| Biotemperature | 0.0–2.4°C (mean 0.6°C) |
| Warmest month mean | 2.1°C |
| Coldest month mean | −6.8°C |
| Annual precipitation | 410–12,320 mm (mean 5,912 mm) |

No vegetation cover — permanently glaciated or snow-covered ground, per
Tappa 3's mass-balance model (accumulation ≥ ablation year-round), not a
temperature threshold. Landscape is bare ice, firn, crevassed glacier
tongues, and exposed rock ribs (nunataks) above the snowline. The very wide
precipitation range (410–12,320 mm) reflects that this class spans both the
wet windward icefields and drier leeward summits — visually, expect the
windward glaciers to be broader and more heavily crevassed (more
accumulation feeding them), leeward ones smaller and more wind-scoured, the
same asymmetry validated in Tappa 3.

**Real-world analog**: the permanent snow and ice fields of the central
Southern Alps/Kā Tiritiri o te Moana — Aoraki/Mt Cook (3,724 m) and the
Fox/Franz Josef/Tasman glacier systems. Tasman Glacier is New Zealand's
longest at roughly 23 km.

---

## Alpine Fellfield

| | this world |
|---|---|
| Elevation | 2,214–3,331 m (mean 2,689 m) |
| Biotemperature | 0.0–1.5°C (mean 0.67°C) |
| Warmest month mean | 2.5°C |
| Coldest month mean | −6.2°C |
| Annual precipitation | 324–4,540 mm (mean 725 mm) |

Holdridge's "Polar" belt, unsplit by moisture (checked directly: essentially
none of this belt falls in the driest tercile, so a separate dry-fellfield
class wasn't warranted). This is the coldest vegetated ground in the domain
— true closed forest and even continuous tundra sod give out here. Expect
open, patchy ground: bare rock, scree, and frost-shattered talus dominate,
with vegetation confined to sheltered pockets between stones rather than
forming a continuous cover. Plants are prostrate, cushion-forming, or
rosette-shaped — growth forms that shed wind and hug the surface for
warmth rather than growing upright.

**Real-world analog**: New Zealand alpine fellfield above roughly
1,800–2,000 m. Vegetation is sparse and low — cushion plants, mat-forming
herbs, and scattered rosettes rather than any woody cover. Characteristic
genera: cushion *Raoulia* ("vegetable sheep," dense grey-green mounds),
prostrate *Hebe* and *Dracophyllum* (mountain heaths reduced to
ground-hugging shrublets here, versus 2–4 m shrubs lower down), rosette
*Celmisia* (mountain daisies), and scattered *Ranunculus lyallii* ("Mount
Cook lily," actually the largest buttercup in the world, with leaves up to
30 cm across) in the wetter, better-watered fellfield hollows. Bare rock and
scree cover a large fraction of total ground area — this is not a "green"
biome even where vegetated.

---

## Alpine Tundra

| | this world |
|---|---|
| Elevation | 1,727–2,376 m (mean 2,008 m) |
| Biotemperature | 1.5–3.0°C (mean 2.3°C) |
| Warmest month mean | 6.0°C |
| Coldest month mean | −3.1°C |
| Annual precipitation | 324–10,625 mm (mean 3,497 mm) |

Holdridge's "Subpolar" belt, also unsplit by moisture. Warmer and lower than
Fellfield, cold enough that trees still can't form a closed canopy, but mild
enough to support a near-continuous low vegetation mat rather than fellfield's
patchy cover between bare stones. Precipitation here spans an enormous range
(324–10,625 mm) — visually this should read as the belt where you can most
dramatically see windward-vs-leeward contrast within a single biome class: on
the wet side, dense, tall tussock and boggy herbfield; on the dry side, a
much shorter, sparser sward closer to fellfield in character.

**Real-world analog**: New Zealand's alpine tussock grassland and herbfield
zone (roughly 1,200/1,500 m up to treeline, treeline itself varying with
local exposure and moisture — Southern Alps treeline is highly variable and
not a single fixed elevation). Dominant plants are the snow tussocks
(*Chionochloa* spp.), forming dense tufted clumps 0.5–1.5 m tall on the
wetter sites; drier or more exposed tundra thins to shorter tussock and
herbfield with speargrass (*Aciphylla* spp., spike-leaved rosettes, some
species over 1 m including the flower spike), mountain daisies (*Celmisia*),
and low mat shrubs (*Dracophyllum*, *Coprosma*). Snow tōtara
(*Podocarpus nivalis*), a genuinely dwarfed conifer under 1 m tall, and
alpine *Carex* sedges also appear here, especially in damper hollows.

---

## Subalpine Wet Forest

| | this world |
|---|---|
| Elevation | 1,035–1,850 m (mean 1,447 m) |
| Biotemperature | 3.0–6.0°C (mean 4.35°C) |
| Warmest month mean | 8.9°C |
| Coldest month mean | −0.4°C |
| Annual precipitation | 1,508–9,857 mm (mean 5,594 mm) |

Boreal belt, wettest tercile. This is the first belt down where the
biotemperature (3–6°C) is high enough to support closed-canopy forest again,
and the precipitation (1,508–9,857 mm, i.e. never actually dry) keeps it
dense and wet. Expect a genuinely closed canopy, heavy epiphyte and moss
load, and a damp, low-light understory — closer to a cloud forest than a
typical temperate woodland.

**Real-world analog**: high-elevation Nothofagus (southern beech) forest,
specifically the wetter mountain-beech and silver-beech stands that persist
right up to treeline on the wet, western side of the Southern Alps. Mountain
beech (*Fuscospora cliffortioides*) grows 12–15 m at these higher, harsher
elevations (it reaches much taller — up to 25 m — lower down, but is
consistently stunted near its own treeline). Silver beech (*Lophozonia
menziesii*) tolerates the widest elevation range of any NZ beech, sea level
to treeline, and reaches 20–25 m even in wetter, higher sites. Understory is
thick with ferns, mosses, and epiphytes given the near-continuous moisture.

---

## Subalpine Woodland

| | this world |
|---|---|
| Elevation | 1,057–1,814 m (mean 1,381 m) |
| Biotemperature | 3.0–6.0°C (mean 4.86°C) |
| Warmest month mean | 9.4°C |
| Coldest month mean | 0.3°C |
| Annual precipitation | 334–2,954 mm (mean 1,301 mm) |

Boreal belt, middle tercile — same temperature band as Subalpine Wet Forest
but roughly a quarter the precipitation on average (1,301 mm vs. 5,594 mm).
Expect a more open canopy than the wet-forest class: trees still dominant,
but with visible gaps, a drier and more sparse understory, and less moss/
epiphyte load. This is the transitional "woodland" reading between closed
wet forest and open dry scrub, not a distinct third structure.

**Real-world analog**: drier montane Nothofagus stands and beech/shrub
mosaics on the same general elevation band, but away from the wettest
windward faces — mountain beech again, but in a more open, lower-density
stand than its wet-forest counterpart, often with a tussock or shrub
understory rather than a closed fern layer. Black beech (*Fuscospora
solandri*, 20–25 m where well-grown) tends toward drier sites lower in this
range than silver beech does.

---

## Subalpine Dry Scrub

| | this world |
|---|---|
| Elevation | 1,067–1,839 m (mean 1,401 m) |
| Biotemperature | 3.0–6.0°C (mean 4.59°C) |
| Warmest month mean | 8.8°C |
| Coldest month mean | 0.3°C |
| Annual precipitation | 324–659 mm (mean 332 mm) |

Boreal belt, driest tercile — the same cool temperature band as the two
classes above, but genuinely arid (332 mm mean, barely a fifth of Subalpine
Wet Forest's). This is the class that specifically validated the
windward/leeward wind-direction decision from Tappa 2: it comes out **100%
leeward** in the current run (335.1 km² on the dry/leeward side, 0.0 km² on
the wet/windward side of the domain) — there is, physically, no wet-side
equivalent. Expect this to look distinctly different from the other two
Boreal classes: no closed canopy at all, low sparse woody scrub, exposed
ground and rock between clumps, a rain-shadow landscape rather than a forest
that's simply been thinned.

**Real-world analog**: the driest, most rain-shadowed subalpine slopes of
inland/eastern South Island (the same rain-shadow effect that produces
Central Otago's semi-arid climate, just at higher elevation here). Expect
low, hardy shrubland — matagouri (*Discaria toumatou*, a thorny, densely
branched shrub to ~3–4 m but usually shorter and wind-shorn at this
elevation), dwarfed mountain beech where it persists at all, *Dracophyllum*
and *Hebe* scrub, and bare stony ground between plants. This is the driest
forest-belt class in the whole scheme and should read visually as scrubby
and open, not wooded.

---

## Temperate Forest

| | this world |
|---|---|
| Elevation | 194–1,154 m (mean 815 m) |
| Biotemperature | 6.0–10.65°C (mean 7.51°C) |
| Warmest month mean | 12.0°C |
| Coldest month mean | 3.0°C |
| Annual precipitation | 3,014–7,323 mm (mean 4,936 mm) |

Cool Temperate belt, wettest tercile. This is also the class that absorbed
the "Temperate Rainforest" sliver from the rejected quartile draft (the old
Superwet quartile), so the wettest, most rainforest-like ground in the
domain sits at the wet tail of this class rather than forming its own
category — expect the wettest sites within this biome (upper end of the
3,014–7,323 mm range, lower elevation) to look noticeably lusher and taller
than the class's drier/higher sites.

**Real-world analog**: South Island podocarp-broadleaf temperate rainforest,
the classic wet-west-coast forest type. This is genuinely tall forest by
world standards: rimu (*Dacrydium cupressinum*) reaches 50–60 m, kahikatea
(*Dacrycarpus dacrydioides*), New Zealand's tallest native tree, has been
recorded up to 60 m and is claimed at 80 m in some older records. Beneath
the podocarp emergent layer sits a broadleaf canopy — kāmahi (*Weinmannia
racemosa*) and rātā (*Metrosideros* spp.) — with a dense fern, moss, and
epiphyte understory. At the cooler/higher end of this class's elevation
range, expect a transition toward beech-dominated forest instead (see
Subalpine Wet Forest above) as the true podocarp giants drop out.

---

## Woodland / Shrubland

| | this world |
|---|---|
| Elevation | 0–1,176 m (mean 293 m) |
| Biotemperature | 6.0–12.0°C (mean 10.01°C) |
| Warmest month mean | 14.3°C |
| Coldest month mean | 5.7°C |
| Annual precipitation | 664–5,291 mm (mean 2,385 mm) |

Cool Temperate belt, middle tercile. Warmest of the land biomes along with
Lowland Steppe/Grassland (both reach the 12°C belt ceiling), but wetter on
average (2,385 mm vs. 531 mm) — this is the mid-moisture lowland class,
sitting between full closed forest and open grassland. Expect a genuine
woodland/shrubland mosaic: patches of lower-stature forest and scrub
interspersed with more open ground, not a uniform canopy.

**Real-world analog**: lowland-to-montane mixed shrubland and open forest
mosaics of moderate-rainfall eastern/central South Island valleys — broadleaf
shrub species (*Coprosma*, *Olearia*, *Hebe*), scattered beech or kānuka
groves, tussock openings. Structurally intermediate between the tall wet
forest above and the open grassland below; expect visible patchiness at map
scale rather than a single uniform texture.

---

## Lowland Steppe / Grassland

| | this world |
|---|---|
| Elevation | 0–1,176 m (mean 230 m) |
| Biotemperature | 6.0–12.0°C (mean 10.45°C) |
| Warmest month mean | 14.7°C |
| Coldest month mean | 6.2°C |
| Annual precipitation | 324–1,310 mm (mean 531 mm) |

Cool Temperate belt, driest tercile — the warmest and driest land biome in
the scheme (mean 531 mm, warmest-month mean 14.7°C, the highest of any
class). This is the direct lowland counterpart of Subalpine Dry Scrub: same
rain-shadow logic, much lower elevation and warmer. Expect open, largely
treeless grassland/steppe — golden-tan rather than green for most of the
year, the driest and most open-looking land biome on the map.

**Real-world analog**: Central Otago and Mackenzie Basin dry tussock
grassland, New Zealand's most rain-shadowed lowland environment. Dominant
species are hard/fescue tussock (*Festuca novae-zelandiae*) and silver
tussock (*Poa cita*), both forming dense golden-tan clumps roughly
0.3–0.6 m tall (taller in flower) rather than a continuous turf — expect a
tufted, semi-open ground cover with bare soil visible between tussocks,
scattered low shrubs (*Discaria toumatou*/matagouri, *Coprosma*), and a
strongly seasonal, sun-bleached color that shifts from green in spring to
gold/tan through summer and autumn.

---

## A note on treeline

The elevation ranges above show real overlap between forest classes
(Subalpine Wet Forest/Woodland/Dry Scrub) and the open alpine classes above
them (Alpine Tundra tops out at 2,376 m; Subalpine Wet Forest's own range
starts at 1,035 m) — this is expected, not an error. Treeline on this world,
as in the real Southern Alps, isn't a fixed elevation contour; it shifts
with local moisture and exposure the same way the belt boundaries in
`holdridge.py` are driven by biotemperature rather than elevation directly.
A real, specific treeline-elevation figure for the Southern Alps was sought
during this research pass but not obtained — the sources found
(Wikipedia's Southern Alps overview) describe the zonation qualitatively
(subalpine scrub and tussock give way to alpine herbfield and fellfield with
increasing elevation) without citing an exact meters figure, so none is
asserted here. Treat the boundary in any given spot as a visual judgment
call informed by the surrounding biotemperature and moisture class, not a
single number to draw on a cross-section.

---

## Sources

Real-world reference facts above are drawn from:

- [Fergus Murray Sculpture — New Zealand's Temperate Rain Forests](https://www.fergusmurraysculpture.com/new-zealand/temperate-rain-forests/) — podocarp species, rimu/kahikatea heights, canopy structure
- [Te Ara Encyclopedia of New Zealand — Southern Beech Forest](https://teara.govt.nz/en/southern-beech-forest/print) — beech species (mountain, silver, black, red), heights, elevation/moisture preferences
- [nznatureguy.com — Kahikatea, tallest native tree](https://www.nznatureguy.com/2019/08/21/8-kahikatea-tallest-native-tree/) — kahikatea height records
- [nznatureguy.com — Rimu tree facts](https://www.nznatureguy.com/2019/09/08/7-rimu-tree-facts/) — rimu height
- [One Earth — New Zealand South Island Montane Grasslands ecoregion](https://www.oneearth.org/ecoregions/new-zealand-south-island-montane-grasslands/) — alpine tussock, herbfield, and fellfield species
- [One Earth — Canterbury-Otago Tussock Grasslands ecoregion](https://www.oneearth.org/ecoregions/canterbury-otago-tussock-grasslands/) — dry lowland tussock species, rain-shadow context
- [Wikipedia — Southern Alps](https://en.wikipedia.org/wiki/Southern_Alps) — general zonation, Aoraki/Mt Cook and glacier figures

This world's own elevation/biotemperature/precipitation figures are computed
directly from `data/processed/climate/` and `data/processed/biomes/` — not
from any external source.
