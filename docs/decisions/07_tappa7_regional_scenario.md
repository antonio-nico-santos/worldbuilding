# Tappa 7 — Regional scenario deepening (planning and decisions)

Status: **RESTRUCTURED this session.** Tappa 1-6 each paired a set of decisions with a delivered
cartographic output. Tappa 7 was scoped far wider than any predecessor — geology/lithology, caves,
Vértice materials, dangerous creatures, dangerous seas, transportation, and the fauna build-out all
in one stage — and most of that content still can't honestly produce a map, because the layers it
would map (lithology classification, road/rail cost surfaces, the hazard layer) don't exist yet.
Rather than force a premature map, this stage is now explicitly framed like **Tappa 0**: a
decisions-and-planning pass, not a closed implementation-plus-map stage. See "Scope and structure"
below for what that means section by section. See `docs/reference/scenario_reference.md` for the
non-technical scenario content this stage draws on (society, fauna, culture) and
`06_tappa6_suitability.md` for the machinery being reused.

## Scope and structure — revised this session

Sections 1-6 below (geology/lithology, caves, jade/Vértice materials, dangerous creatures/conflict
zones, dangerous seas, transportation/infrastructure) stay **decision-only**, the same way Tappa 0
was decision-only. None of them can produce a real map yet — each depends on a layer that hasn't
been computed (lithology classification, road/rail cost functions, the wind-shadow mask, the
hazard layer). Their eventual cartographic output will each become its own future Tappa once its
blocking layer actually gets built. Deliberately **not assigning specific future Tappa numbers to
these yet** — the point of treating this stage as decisions-first is to keep that ordering open
until work resumes, rather than locking in a sequence now. (Tappa 8, urban zoom, already has a
committed scope from Tappa 0/6 and is unaffected by this change — it keeps its position regardless
of how the geology/transport/seas content above eventually gets split out.)

Section 7 (fauna) is the exception: it's the one piece of Tappa 7 content that's actually
**unblocked** right now. Species range and habitat placement only need layers that already exist —
Tappa 5's biome polygons, Tappa 1's terrain — so, unlike the other five domains, it can honestly
produce a cartographic deliverable without waiting on anything else in this document. That map
(fauna distribution/range) is scoped at the end of §7 below as this stage's one closing
deliverable, the same way Tappa 1-6 each closed with their own map. Worth being precise about what
kind of map it is: unlike Tappa 1-6's algorithmically-derived outputs (noise functions, watershed
extraction, suitability composites), a fauna range map is **authored from existing lore**, not
computed — each species' biome and range were already decided in prose in `scenario_reference.md`;
building the map means encoding those existing decisions as polygons/attributes and rendering them,
not writing new geoprocessing code.

## 1. Geology / lithology layer

Four rock classes, reusing existing Tappa 1 skeleton geometry rather than authoring anything new:

- **Schist** (metamorphic core) — within each ridge's own `falloff_km` of its axis. Real,
  citable mechanism: metamorphic grade increases toward the Alpine Fault in the Haast
  Schist/Southern Alps analog (this project's own validation reference). No real number exists
  for the schist belt's actual width, so reusing each ridge's own already-authored `falloff_km`
  as its schist radius is the zero-new-parameter solution — physically motivated (bigger uplift
  → plausibly wider exhumed zone), not invented.
- **Greywacke/argillite** (flank) — beyond the schist radius, out to each ridge's shelf reach.
  Same citation: unmetamorphosed Torlesse/Caples-terrane protolith.
- **Sedimentary basin fill — now cited.** The existing plains/plateau zone polygons (NW Plateau,
  North plains, Central plateau, South plains, SE plains, Island plateau). Real analog: the
  Canterbury Plains (NZ) — Quaternary glacial outwash gravel, transported from the Southern Alps
  during Pleistocene glacial cycles (roughly 3 million to 10,000 years ago) and reworked into
  coalescing shingle fans by braided alpine rivers (Waimakariri, Rakaia, Rangitātā). Mechanism:
  glacially-fed rivers draining mountains build broad gravel plains downstream — a real, common
  process, not specific to Canterbury, and it costs nothing new in this project's own citation
  set since Canterbury is already the validation reference for the schist/volcanic classes above.
  **Caveat that still stands**: this supports the *mechanism*, not specific numbers — the "3
  million years / gravel depth" figures describe Canterbury's real depositional history, not
  something this project's static DEM models or should claim to match.

  **Hydrology check — RESOLVED, checked directly against Tappa 4 output, not just asserted.**
  Rasterized all six plains/plateau zone polygons (`data/input/terrain_zones.geojson`) against the
  actual `stream_mask.npy` (377,281 stream cells, area ≥0.3 km² threshold, per
  `04_tappa4_hydrology.md`). Every zone is genuinely threaded with stream network, not just
  adjacent to it: mean distance to the nearest mapped stream is under 0.5 km for five of six zones
  (Island plateau highest at 0.72 km, still close), and 84.5-99.7% of each zone's area sits within
  1 km of a stream:

  | zone | area (km²) | stream km inside | mean dist to stream (km) | % of area within 1 km |
  |---|---|---|---|---|
  | NW Plateau | 572.8 | 749.1 | 0.27 | 99.7% |
  | North plains | 1385.1 | 1397.9 | 0.46 | 90.6% |
  | Central plateau | 778.9 | 1049.9 | 0.25 | 99.7% |
  | South plains | 939.2 | 1021.0 | 0.35 | 96.2% |
  | SE plains | 720.8 | 900.2 | 0.28 | 99.3% |
  | Island plateau | 246.5 | 185.6 | 0.72 | 84.5% |

  This is exactly the geomorphic pattern the Canterbury Plains citation requires — braided,
  ridge-fed rivers actually threading through the lowland zones, not an isolated dry plain — so
  the sedimentary basin fill class moves from "plausible" to "grounded and checked" against this
  project's own computed data, not just against the real-world analog.
- **Volcanic** — the Caldária/SW Island zone. Strongly cited via Banks Peninsula (Canterbury,
  NZ): real intraplate volcanism geologically unrelated to the Alpine Fault orogeny, erupted
  through the same Torlesse-terrane basement, with a cause still debated in reality too.

**Precedence rule — RESOLVED.** Tappa 1's `max` compositing doesn't translate directly, because
this isn't two continuous signals of the same kind (like two ridges' elevation contributions) — it's
a categorical classification where three of the four classes are positive evidence (schist,
greywacke, volcanic each require a specific geometric condition to be true) and the fourth
(sedimentary basin fill) is not: it's the default label for anything the other three don't claim.
So the resolution is a straightforward priority order, not a numeric max: **volcanic > schist >
greywacke > sedimentary basin fill**, where basin fill is assigned only to cells left over after
the other three are painted, not competed against them as an equal polygon. Wherever a ridge's
schist radius geometrically overlaps a plains/basin polygon, schist simply wins — no blending, no
sum, nothing new to implement beyond the paint order. Not yet coded, but no longer an open
decision.

## 2. Caves

Four types grounded, one deliberately deferred:

- **Lava tubes** — volcanic class only. Strong citation (Banks Peninsula).
- **Talus/pseudokarst caves** — steep mainland relief (slope) intersected with a stream cutting
  beneath it (`stream_mask`). Rock-type agnostic by citation (forms in "granites, gneiss, or any
  fractured rock that forms large angular boulders") — works on both schist and greywacke without
  waiting for the lithology layer to be built.
- **Glacier/moulin ice caves** — reuses Tappa 3's `permanent_snow_mask` directly, zero new
  geometry. Real mechanism (meltwater exploiting weaknesses in glacier ice) is independent of
  bedrock, so the mainland's non-carbonate lithology doesn't disqualify it.
- **Sea caves — new, added this session.** Mechanically distinct from all three above: hydraulic
  wave pressure exploiting zones of structural weakness (joints, fractures) in a coastal cliff,
  enlarging the cavity over repeated wave impact. Real citation, and notably rock-type agnostic in
  the sources found — occurs on "almost every cliffed headland or coast where the waves break
  directly on a rock cliff," not tied to a specific lithology, so it doesn't need to wait on the
  lithology layer either. Motivated by Povo Silencioso's NE archipelago dwelling in coastal
  caves/simple stone construction (see `scenario_reference.md`). Geometry needed: coastline +
  cliff steepness (coastal relief), not yet built. Possible cross-reference worth checking later:
  if the archipelago's assigned lithology turns out to be greywacke, the same joints/shatter zones
  that host laumontite (section 3) are plausibly the same weaknesses waves would exploit — not
  confirmed, just a coherent connection if it comes up.
- **Schist fracture caves** — explicitly **not adopted**. The real mechanism (tectonic/fracture
  caves) is documented specifically for sedimentary rock with contrasting layer properties;
  nothing confirms it for schist. Parked pending Nico's visual review of the other three
  (now four) distributions.

## 3. Jade / pounamu resource

Haast Schist pounamu occurs as isolated pods specifically in the highest metamorphic grades near
the Alpine Fault — maps directly onto the highest-grade schist band (closest to each ridge axis),
i.e. the steepest, most dangerous terrain. Real citation, zero new geometry. Now doubles as the
in-world source material for Vértice interpreter crystals (see `scenario_reference.md`) — ties
the resource's economic/narrative weight to a location the pipeline already produces.

**Vértice materials — RESOLVED.** Full system (domains, verbs, tiers) lives in
`scenario_reference.md`; this is the geology-facing summary of which rock hosts which material.

- **Schist — richest of the four.** Gold-bearing quartz veins (orogenic, fault-zone; real
  citation: Otago Schist's Invincible Vein — quartz, gold, pyrite, arsenopyrite, muscovite,
  chlorite, calcite, albite). The quartz itself is the Onda-domain Vértice crystal (birefringence);
  mica (muscovite/biotite) from the same veins is the rarer Energia-domain crystal (piezoelectric).
  Gold and jade/pounamu (see above) are mundane precious materials from the same terrane, not
  Vértice crystals themselves — neither has an electrical, magnetic, or optical property that
  gates a domain.
- **Greywacke/argillite — laumontite (a zeolite).** Fills veinlets, joints, and shatter zones
  throughout NZ's Mesozoic greywacke ranges — genuinely specific to this rock class. This is the
  Matéria-domain crystal (real zeolite molecular-sieve structure, a direct fit for
  solid/liquid/gas).
- **Sedimentary basin fill — vivianite (primary) + reworked placer magnetite (weak secondary).**
  Vivianite is a real authigenic iron-phosphate mineral that forms in low-oxygen, organic-rich
  floodplain/bog sediment — the Bios-domain crystal, tied directly to organic decay in exactly
  this depositional environment. The weak Campo (magnetism) tier here is reworked titanomagnetite,
  eroded from the volcanic zone upstream and reconcentrated by rivers (real NZ ironsand-beach
  provenance mechanism) — lower-grade than the volcanic primary below, not a distinct source.
- **Volcanic (Caldária/Banks Peninsula) — magnetite (primary Campo crystal), replacing the old
  olivine claim.** Basalt's titanomagnetite content is the real basis of paleomagnetism — a much
  stronger citation than "olivine is present, gem quality unconfirmed," which is dropped. Volcanic
  epithermal veins (real NZ analog: Hauraki Goldfield) also carry native silver as a rare
  accessory — the higher-quality half of the Mente-domain "tuning" pair (native copper, common,
  from volcanic amygdules — real Keweenaw Peninsula analog — being the lower/common half). Mente's
  crystal isn't purely mineral: silver or copper alone gives a weak connection, refined by pairing
  with a black fungus found in Wet Forest (combination method still in-world "under study" as of
  present day; the fungal refinement was only discovered ~3-4 years before the TTRPG's present).

## 4. Dangerous creatures / conflict zones

Reclassified in Tappa 6 as an **informational conflict layer**, not a suitability penalty — still
holds. Status per seed:

| seed | status |
|---|---|
| Alpine terrestrial apex predator | **RESOLVED as Reaper.** Real dromaeosaurid-anchored (Utahraptor-scale), territorial, alpine/tundra band, confined there by real anthropogenic range-contraction logic rather than preference. Seasonal snowmelt-triggered lowland descents, resource-scarcity-driven — see `scenario_reference.md` §16 |
| Aerial "sky" danger | Resolved as a **narrative device**, not an ecology model — see `scenario_reference.md`; do not build as ordinary wildlife territory |
| Volcanic/geothermal hazard (Caldária) | Buildable — ties directly to the volcanic lithology class above |
| Seasonal toxic bloom | **RESOLVED: Wet Forest.** External hazard for Povo Livre rather than internal Círculo-economy tension — see `scenario_reference.md` for the narrative tie-in (Povo Livre's own knowledge of bloom timing/location becomes a point of local expertise). Should still anchor to an existing Tappa 2/3 seasonal month, not an invented one — not yet picked |
| Windward storm corridor | Buildable — physical-only hazard, consistent with the "independent" half of the sea-danger model below |
| Círculo–Povo Livre frontier tension | **Direction set, not yet built.** Map it as a conflict-zone layer: candidate cells are those with high Povo Livre suitability (>80%, reusing Tappa 6's forestry-foraging composite) that also sit close to Círculo territory and infrastructure connections (roads once built). To be added once regional-scale maps are actually being produced, not before |

## 5. Dangerous seas

Decided structure: **creature-driven hazard layer + physical/independent hazard layer, summed**
— not two competing systems. Carries forward a known caveat from Tappa 6:
`BOAT_SPEED_KMH=6.0` is ≈ or faster than Tobler's flat walking pace, so sea isn't a real barrier
in the current cost model.

**Navigability — RESOLVED.** Rather than treat the whole sea as passable-with-a-speed-penalty (the
Tappa 6 approach), a small, explicitly-designated set of navigable corridors will be authored
directly (e.g. the mainland↔SW-island ferry route) instead of derived from the hazard layer. This
doesn't touch the Tappa 6 suitability map at all — it's scoped entirely to the danger-seas layer
once that's built, which is also where the `BOAT_SPEED_KMH` caveat above finally gets addressed
(a designated-corridors model sidesteps it rather than fixing the underlying speed number, which
is an acceptable resolution since open-ocean routing was never the intent).

## 6. Transportation / infrastructure

- **Rail** — solar-electric, in scope. Needs a genuinely new cost function (hard grade ceiling
  ~2-4%), not a relabeling of Tobler friction. Not built.
- **Road** — reuses `cost_distance.py`'s Tobler-friction + Dijkstra graph directly; needs
  predecessor extraction for actual path reconstruction (current code only returns costs). Not
  built.
- **Biome-differentiated network cost — new idea, not yet built.** Wet Forest should cost more to
  build/maintain a road or rail segment through than Grassland — real-world justification is
  straightforward (denser vegetation clearance, wetter/softer ground, more maintenance against
  regrowth) and gives the transportation-network layer a biome-sensitivity it currently lacks
  entirely (both rail and road cost functions above are terrain-slope-only so far). Applies as a
  per-biome cost multiplier on top of whichever grade/friction function each mode already uses.
  To be implemented when the road/rail code is actually written, not before.
- **Trams** — urban-scale, large Círculos only. Deferred to Tappa 8 (urban zoom) — Tappa 6's
  site-selection windows are too coarse for tram-layout detail below the four large Círculos.
- **Electric utility vehicles, bikes** — ride the road network once it exists; no new modeling.
- **Domesticated animals** — yak (alpine/Spine niche), reindeer (Círculo
  Grassland ecotype + Povo Livre/nomad forest ecotype), sheep, dogs, chickens — each diverged into
  biome/ecotype variants. Full grounding in `scenario_reference.md` §7.
- **Boats** — slow, scheduled ferries, running only on the small set of explicitly-designated
  navigable corridors decided above (e.g. mainland↔SW volcanic island), not the open sea generally.
  Should respect the existing Povo Silencioso 5km/15km exclusion buffer (already built into Tappa 6
  suitability) if any route runs near the archipelago. **Motorization — RESOLVED**: boats may be
  motorized, but only within the designated navigable corridors — motorization doesn't extend
  their range into the danger-seas zones; the corridor boundary, not engine power, is what limits
  where a boat can go.
- **Kite transport** (land + water) — real technology (kite buggy/landboarding), not invented.
  Water routes (lakes, coastline) are ready now — open water has no obstruction. Land routes need
  one new piece: a directional wind-shadow mask (ridges mechanically shelter their lee, given the
  fixed 250° wind bearing — same cost class as Tappa 1's shelf falloff) plus a favorable/
  requires-tacking label per route segment. **Scope deliberately kept simple, per Nico**: no need
  to resolve wind *speed* at ground level if that's expensive to derive — the map should show
  where and in which direction kite travel is most viable, with an approximate/qualitative speed
  label, not a fully resolved near-surface wind field. Not built.
- **Aircraft** — explicitly excluded. Gated narratively, not modeled as a routing option; see
  `scenario_reference.md` for the in-world explanation.

## 7. Native and migratory fauna — niche build-out

New multi-session initiative, run niche by niche starting at the coast and working inland toward
the tundra, ending with the still-open alpine apex predator from §4's dangerous-creatures table.
Species counts are set per niche before any naming/trait work happens. Two different design
methods apply depending on category: native/resident species use the same method already used for
Povo Silencioso and the aerial sky-danger — a real trait from a different Earth lineage,
transplanted into a combination Earth itself never produced. Migratory/seasonal-visitor species are
handled differently on purpose, kept closer to their real Earth templates, since their narrative
job is the opposite of the residents': proving the world extends past this island, not marking it
as exotic. Migratory species are explicitly not counted as island fauna.

**Status:**

- **Coastal/marine niche — DONE.** Two resident species (Shadowless Sardine, Moonfur) plus six
  migratory visitors (Green Mackerel, Sea Terror, Bronzeshell, Stained Penguin, Skydrifter, and 
  SummerVisitant/Wetducks). Full species detail in `scenario_reference.md` §8.
- **Grassland large-grazer niche — DONE.** One species, Grassmothers — bison-anchored matriarchal
  bovid with a windiest-region eye-fringe ecotype (real citation: Highland cattle forelock). Full
  detail in `scenario_reference.md` §9.
- **Mid-size browser niche — DONE.** Two species, one per biome: Blacknose (Temperate Forest,
  tapir-anchored) and Tailstand (Woodland-Shrubland, gerenuk-anchored). Full detail in
  `scenario_reference.md` §10.
- **Reptile/amphibian niche — DONE.** Three species (expanded from the originally scoped two):
  Mudlizard (riverbank/wetland, Nile-monitor-anchored, real platypus electroreception), Flashfrog
  (Wet Forest, poison-dart-frog-anchored, toxin sequestered from the toxic bloom's insect
  consumers, cephalopod-chromatophore flash display), and Clicksnake (Grassland,
  prairie-rattlesnake-anchored, venomous, human-audible echolocation clicks doubling as an
  unintentional warning, aposematic red-eyed black/yellow coloring echoing male Grassmothers). Full
  detail in `scenario_reference.md` §11.
- **Small mammal niche — DONE.** Three species (expanded from the originally scoped two):
  Scattermouse (Woodland-Shrubland, jerboa-anchored, meerkat-style mobbing plus confusion-effect
  scatter defense), Quillhog (wide-ranging, hedgehog-anchored, barbed porcupine-style quills
  anointed with the Wet Forest toxic bloom where present, milder local toxins elsewhere), and
  Snaketail (Temperate-Forest/Grassland ecotone, tree-hyrax-convergent climbing rabbit-analog,
  tail mimics Clicksnake's aposematic coloring — deliberately placed at the biome edge so the
  mimicry's predator-overlap requirement actually holds). Full detail in `scenario_reference.md`
  §12.
- **Mesopredator niche — DONE.** Three species (expanded from the originally scoped two),
  deliberately diverse per mesopredator release given the apex predator's alpine restriction:
  Treefox (red-fox-anchored generalist, raccoon-style dexterous paws, elephant-derived
  ground-vibration sense, near-180° ankle rotation for headfirst tree descent, real threat to
  Snaketail), Furypack (pine-marten-anchored, wolf-style coordinated pack hunting grafted onto a
  normally solitary mustelid), and Twinshadows (fossa-anchored cougar-convergent predator scaled
  to ~50-70kg, melanistic countershaded coloring, silver eyes, hunts in bonded pairs of exactly
  two). Full detail in `scenario_reference.md` §13.
- **Bird niche — DONE.** Two resident species plus the previously-deferred sixth coastal migrant,
  now resolved: Trinketbird/Wingedthief (bowerbird-anchored object-collector, cryptic plumage per
  the real "extended phenotype" citation, magnetoreception repurposed for bower-site memory, a
  genuine nuisance targeting Vértice crystals and wire/metal fittings) and Rustowl (pygmy-owl-
  anchored, small and non-intimidating, reddish-bronze/rufous plumage, real pit-viper infrared
  heat-sensing grafted onto the facial disc, Beau Geste-style multi-perch calling to sound like a
  chorus, real false eye-markings, known publicly as the "owl choir"). The sixth coastal migrant is
  now named and closed out: SummerVisitant/Wetducks, a steamer-duck-anchored seasonal visitor to the
  island's interior lakes, publicly documented only as an unresolved migration-origin mystery. Full
  species detail in `scenario_reference.md` §14.
- **Alpine prey base niche — DONE.** Two species, deliberately paired around complementary
  predator-detection strategies: Cryburrow (groundhog-anchored, scaled to ~1.3x real size — a
  slower-digging trade-off that forces more above-ground exposure — grizzled pale-gray/white fur
  extending the real hoary marmot's own coloring, prairie-dog-style referential alarm calls
  encoding predator type rather than a generic warning) and Deergoat (alpine-ibex-anchored, kept at
  real scale and standard coloring, real mule-deer-proportion cupped/funnel ears for long-range
  passive sound detection, horse-style ear-pinning alert response). Together they give whatever
  hunts them two distinct warning channels — near-range social/acoustic-specific and long-range
  solitary/acoustic-directional — a deliberate constraint carried forward into the final niche.
  Full detail in `scenario_reference.md` §15.
- **Tundra/alpine apex predator niche — DONE. This closes the fauna build-out initiative.**
  Reaper: a real dromaeosaurid-anchored (Utahraptor-scale, ~7m/~500kg) apex predator, the sole
  fully-fictional-lineage (extinct-Earth-anchored, not modern-Earth) resident species in this
  fauna set, and the only one given a single proper-noun name rather than a trait-compound name —
  deliberately marking it as this arc's capstone. Feathered/insulated, mesothermic (real, debated
  citation: Grady et al. 2014), countershaded gold-brown coloring tuned to Grassland/foothill
  terrain (its prey there hunts by sound, not sight, so alpine camouflage would have been
  functionally wasted), small structurally-iridescent blue-green accent patches (a real color
  range unavailable to any mammal in this world). Confined to the alpine/tundra band by real
  anthropogenic range-contraction, not preference. Four known family packs of 4-6, informally
  tracked/named by direction ("the South family") — real wolf-style dispersal-and-challenge social
  structure, not patricide. Baseline diet is yaks (direct territorial overlap) plus opportunistic
  alpine fauna; rare seasonal snowmelt-triggered descents into the Grassland lowlands, driven by
  real resource-scarcity prey-switching, are multi-day lingering events, not single-night raids —
  the actual "hurricane"-scale danger, distinct from the chronic low-level yak-predation risk
  herders already manage. Current hook: the southern pack's recent leadership change is driving
  unusually frequent incursions, breaking the normal seasonal rarity pattern. Full detail in
  `scenario_reference.md` §16.

**All seven planned fauna niches are now complete**, from the coast to the tundra apex predator —
the multi-session initiative opened at the start of this section is closed.

**Post-closure addition — decomposers/scavengers, DONE.** A genuine structural gap identified only
after the arc above had already closed: nothing built accounted for what happens to the biomass the
arc's own large-bodied residents produce (Reaper's multi-day kills, Grassmothers mortality,
livestock losses). Two species, deliberately wide-ranging across biomes rather than tied to one:
Farsmell (condor-scale vulture-anchored aerial scavenger, dark plumage, real turkey-vulture
olfaction extended past any real vulture's range to detect carrion scent under forest canopy) and
Meatcleaner (Tasmanian-devil-anchored ground scavenger, kept at real devil scale, total-carcass
consumption as real disease control, communal screech display extended into deliberate rival
intimidation, real snake vomeronasal chemoreception grafted on for ground-trail tracking — a
distinct detection mechanism from Farsmell's airborne search). Full detail in
`scenario_reference.md` §17.

**Naming convention — locked.** Species names follow a strict register split. Domesticated animals
(yak, reindeer, sheep, dogs, chickens — §7 of `scenario_reference.md`) keep their plain, unmodified
real Earth names. Every wild species instead uses an invented trait-compound name (Cryburrow,
Deergoat, Clicksnake, Grassmothers, and so on). A third, higher register is reserved specifically
for entities whose danger sets them apart from ordinary wildlife: a single proper noun rather than
a compound descriptor. Reaper is the only fauna example; the sky-watcher (§4 above) uses the same
register. This is a naming rule, not a taxonomy change — Reaper is still explicitly ordinary
wildlife, not myth or Vértice-tied; the naming choice alone marks the danger tier.

**Cartographic deliverable — this stage's one map product, SPEC IN PROGRESS, NOT YET BUILT.** Fauna
distribution/range maps, authored from the biome/range placements already decided above and in
`scenario_reference.md` §8-§17, rendered against the existing Tappa 5 biome layer and Tappa 1
terrain. This is the only map Tappa 7 will close with — see "Scope and structure" at the top of
this document for why the other five domains stay decision-only for now.

**Map structure — four content buckets, RESOLVED.** Species split by which existing machinery their
distribution can honestly be derived from, not by taxonomy:

- **Regular/resident land fauna** — a species × biome percent table (independent per-cell "chance
  to find the species in this biome," not row-normalized — a species absent from a biome is 0% in
  that cell, and there's no requirement that a row sum to 100%), joined against the real Tappa 5
  `biome_id` polygons. Covers most residents: Grassmothers, Blacknose, Mudlizard, Flashfrog,
  Scattermouse, Quillhog, Snaketail, Treefox, Furypack, Trinketbird, Rustowl (mainland ecotype),
  Cryburrow, Deergoat, Farsmell, Meatcleaner, **and Clicksnake, reclassified into this bucket this
  session** — see below.
- **High-threat species** — Reaper and Twinshadows only (Clicksnake reclassified out, see below):
  Tappa-6-style suitability composites (weighted layer stacks), not the percent table, since each has
  a hard biome/habitat requirement plus graded preferences within it that a flat percent table can't
  express (Reaper: alpine/tundra + anthropogenic range-contraction + prey-base proximity; Twinshadows:
  dense low-light forest + prey-base proximity). Composite weights LOCKED this session, see below.
- **Domesticated animals** — three primitives, LOCKED this session (method only — not yet computed,
  see below): point-per-Círculo for Círculo-kept ecotypes (straight lookup against the existing
  17-Círculo layer, `biome_at_site == Grassland` — reindeer's Círculo ecotype, sheep's lowland
  ecotype; dogs' herding type and chickens' dooryard type are ubiquitous, no filter, all 17); a zone
  against the existing `suitability_povo_livre_120m` surface (Tappa 6 built it but never used it —
  reindeer's forest ecotype, dogs' companion/tracking type, chickens' free-ranging landrace,
  threshold >80%, reusing the same cutoff already proposed for the frontier-tension conflict layer
  in §4); and a new third primitive — **alpine resource outposts** — see below.
- **Coastal/marine and migratory species** — authorial shapes, not computed, the same "hand-author
  when there's no real derivation path" method as Tappa 1's ridge/zone control vectors. Covers every
  migratory visitor (Green Mackerel, Sea Terror, Bronzeshell, Stained Penguin, Skydrifter,
  SummerVisitant) *and* the coastal/marine residents (Moonfur, Shadowless Sardine), since there's no
  coastal-biome polygon to join percentages against either way. Shape sub-types: seasonal coastal
  zones (Green Mackerel, Bronzeshell, Stained Penguin), a transit corridor line (Skydrifter, since it
  overflies rather than settling), a point/polygon over specific interior lakes (SummerVisitant), and
  deliberately uncertain-styled (dashed/lower-opacity) sighting zones for irregular, non-annual
  arrivals (Sea Terror — overlapping Green Mackerel's zone, since the text already has it arriving
  *after* Green Mackerel, opportunistically). Moonfur gets two shapes: a steady SW-coast range plus a
  secondary relocation zone that only activates when Sea Terror is present, so the map can actually
  show the folk-warning behavior already written for it.

**Windward/leeward asymmetry — CONFIRMED against real computed data, not assumed.** Checked
`05_tappa5_biomes.md`'s own validation table rather than guessing: the 250° WSW wind bearing locked
in Tappa 2 already produces a near-total biome split, SW (windward/wet) vs. NE (leeward/dry) — not
north/south. Lowland Steppe/Grassland is 73.4% of the leeward tercile and 0.0% of the windward
tercile; Temperate Forest and Subalpine Wet Forest run the reverse (20.9%/26.5% windward, 0.0%
leeward). This means the percent-join bucket inherits this asymmetry automatically, for free, as
long as it joins against the real biome polygons rather than an abstract "this biome exists
somewhere" assumption — no species-level north/south logic needed.

**Biome fragmentation — CHECKED via connected-component analysis on the actual `biome_id` raster**,
not assumed from the branching terrain skeleton alone. Most biomes are one contiguous mass (Grassland
83.8% single-patch, Woodland/Shrubland 78.4%, Subalpine Wet Forest 96.9%, Subalpine Dry Scrub 98.1%)
— no fragmentation concern, the percent join is safe as designed. Three are genuinely, substantially
split across separate massifs on different arms of the branching spine (main spine / North branch /
West branch / South branch, per `terrain_ridges.geojson`): Alpine Tundra (480 km² on the main spine,
136.5 km² ~40 km away on the South Branch), Permanent Snow & Ice (684.5 km² main spine, 263.2 km²
further SE), and Temperate Forest (456.2 km² near the main spine, 196.6 km² ~35 km away near the West
Branch). This directly informs the still-open Reaper suitability composite: the "4 known family
packs" lore already fits a main-spine-massif / South-Branch-massif split that wasn't authored with
this data in mind — worth anchoring specific packs to specific massifs rather than painting one
undifferentiated alpine range. Same logic applies to Cryburrow/Deergoat (same biomes) and to Blacknose
on the Temperate Forest side. Exact pack-to-massif assignment: still open.

**South/SW island — CONFIRMED as a real, separate landmass**, not just named in prose. Landmass
connectivity check on `biome_id != Ocean`: mainland is 8,879 km², a second landmass is 794.6 km² with
a **16.44 km minimum water gap** to the mainland — matches the "SW volcanic island"/ferry-route lore
in §6 above and the Caldária framing in §3. Land cover: 69.1% Woodland/Shrubland, 30.9% Grassland, no
forest/alpine/snow (island tops out at 580m, per the "Island ridge" feature in
`terrain_ridges.geojson`) — and it does carry a real stream network (~569 km), confirmed by checking
`stream_mask.npy` against the island's footprint, so wetland/riverbank habitat is real there too, not
assumed.

**Mainland vs. island — separate joins, RESOLVED as the model; species-by-species calls IN
PROGRESS.** A naive join would mechanically paint mainland species onto the island wherever the
biome type matches, but 16 km of open water is a real dispersal barrier real island biogeography
(MacArthur & Wilson) says most terrestrial fauna doesn't cross without help. Default: no species gets
island presence via the ordinary join unless its real-world anchor supports the crossing (strong
flight/swimming) or there's an authored in-world reason (introduction). Calls made so far, all
grounded in each species' existing real-world anchor rather than a blanket rule:

- *Excluded on habitat alone* (biome doesn't exist on the island): Blacknose, Snaketail, Furypack
  (need Temperate Forest/its ecotone), Cryburrow, Deergoat, Twinshadows (need Alpine/Tundra or dense
  forest).
- *Excluded on dispersal despite matching habitat*: Flashfrog (real amphibians essentially never
  cross open water naturally — also Wet-Forest-locked anyway), Grassmothers both ecotypes (no
  real-world precedent for megafauna-scale open-water crossing).
- *Strong presence*: Mudlizard (real monitor lizards are well-documented open-water swimmers —
  Komodo dragons/relatives are the standard citation — and the island's real stream network gives it
  habitat too), Farsmell (condor-scale soaring flight makes the gap trivial; scavenger, doesn't need a
  resident prey base, just occasional carrion).
- *Moderate*: Scattermouse (small rodents are the classic real storm-rafting colonizer).
- *Low-moderate*: Treefox (real foxes are documented swimmers over several km, plus the
  ecological-flexibility citation already in its own design), Quillhog (small-bodied, rafting-
  plausible, plus its own "wide-ranging" design already implies habitat flexibility).
- *Excluded on body-plan grounds specifically*, not genericially: Tailstand — see below, moved out of
  the natural-dispersal model entirely.
- *Excluded, CONFIRMED*: Meatcleaner (ground-bound scavenger, no swim/flight precedent in its real
  anchor).

**High-threat bucket needs the same per-species island review, CONFIRMED after Nico caught the gap** —
Reaper and Twinshadows happen to already be excluded from the island (no alpine, no forest), but that
was incidental to their own habitat requirements, not because the bucket had actually been reviewed.
**Clicksnake — LOCKED as present on the island**, same High threat level as the mainland population.
Stronger case than most of the regular-bucket "low-moderate" calls above: reptiles are disproportionately
successful natural island colonizers compared to similarly-sized mammals, specifically because low
metabolic rate lets them survive a debris-rafting crossing without feeding — a standard real
herpetological-biogeography citation — and unlike Reaper/Twinshadows, Clicksnake's habitat requirement
(near-uniform Grassland) is actually present on the island (30.9% of its area). One non-blocking detail:
Clicksnake's coloring was written as "a distant echo of male Grassmothers' coloring, since the two
species share the same Grassland habitat" (§11/13 of `scenario_reference.md`) — Grassmothers is
excluded from the island, so that specific justification is mainland-only flavor, not a functional
requirement; doesn't affect whether Clicksnake belongs there. **Reclassified this session** into the
regular/resident percent-table bucket (see "Map structure" above) — its mainland/island values now
live in `tappa7_fauna_biome_percent.xlsx` as an ordinary Specialist row (Grassland 100%, Woodland/
Shrubland 15% mainland; same values with no colonization discount on the island, per the parity call
above) rather than a separate high-threat composite.

**Rustowl — island-origin twist, LOCKED.** Two ecotypes rather than a uniform species, same pattern
as Grassmothers' windy ecotype and the domesticated-animal ecotypes: **Island Rustowl** is the
ancestral, original population (and the *stronger* one on the map — core range, not a minor
secondary population), **Mainland Rustowl** descends from a rare wind-driven founder event reaching
the mainland (real, documented mechanism in ornithology — storm-displaced vagrants establishing new
breeding populations, distinct from ordinary migration), and is the smaller, more recently established
population. Chosen over Trinketbird specifically: Trinketbird's defining trait (bower-site memory,
stealing Vértice crystals/metal fittings) only makes sense already living alongside Círculo
infrastructure, which rules out an island-first origin; Rustowl's traits (heat-sensing, vocal
chorus-deception hunting) have no such dependency, and real small owls are actually weak, sedentary
dispersers — which is what makes the "rare accidental crossing" framing earn its weight rather than
being redundant with ordinary flight capability. Trinketbird stays mainland-only, no island presence.

**Ecotype distinction — LOCKED, functional not visual.** Deliberately no external/physical difference
between Island and Mainland Rustowl. Two internal distinctions instead:

- **Heat-sense recalibration.** The island's geothermal activity (Caldária-linked volcanism) raises
  ambient thermal background and adds noise, which matters specifically because pit-organ-style
  infrared sensing detects a *contrast* between target and background temperature, not an absolute
  value — real thermal-sensing biology, not an invented mechanic. Island Rustowl's heat sense is
  recalibrated to filter that geothermal noise, keeping it functional near volcanic terrain, but this
  is a genuine trade-off, not a strict upgrade: it's *less* sensitive to subtle prey heat signatures in
  the island's cooler, non-volcanic ground, unlike Mainland Rustowl, which stays tuned to an ordinary
  cool background throughout its range.
- **The choir, genuinely real on the island.** Refined from Nico's original framing to stay consistent
  with the established solitary-hunter design (real pygmy owls, and Mainland Rustowl as written, don't
  cooperatively hunt) rather than overriding it: Island Rustowl's higher population density (consistent
  with being the ancestral, core population) means individual territories overlap enough that several
  real birds are independently running their own solo perch-switching mimicry trick at the same time,
  in range of each other — not a new cooperative behavior, just real numbers colliding with the
  existing trick, compounding into something louder and more chaotic than any one performer intended.
  This produces a folklore inversion worth keeping: on the mainland, the specialist-knowledge fact is
  "the owl choir is usually just one bird." On the island, that same specialist knowledge would be
  *wrong* — there, it really is several.

**Tailstand — human introduction, LOCKED, moved out of the natural-dispersal model entirely.**
Recently introduced by Círculo settlers via the existing ferry route (aesthetic motive — found them
beautiful), not natural colonization; population growing because its two mainland checks, Twinshadows
and Furypack, are both already habitat-excluded from the island for the same Temperate-Forest/dense-
forest reasons listed above — the predator-release mechanism was already sitting in the existing
design, just needed the right species named (an earlier draft of this idea used Rustowl as the
missing predator, which doesn't hold up: Rustowl is ~15-17cm, built for small prey, with no plausible
path to threatening an 80kg animal — corrected before locking). Framed as a live, modest, actively-
debated community concern, not a crisis — and the Círculo blind spot behind it is deliberately narrow:
not general ecological naivety (which would sit badly against everything else established about
Círculo society — science-specialized councils, solarpunk framing, natural-cycles culture), but the
specific fact that cross-water species introduction has never been a live problem for a civilization
that has only ever existed on one connected landmass. Island presence for Tailstand is authored
(point-source-plus-growth), the same category as the migratory/marine shapes above, not the percent
join.

**Mainland/island percentages — LOCKED and CLOSED this session**, converting every qualitative
tier above into a number (`tappa7_fauna_biome_percent.xlsx`, "Mainland → Island %" sheet, sent to
Nico, not committed — working deliverable like the rest of this workbook). Only Woodland/Shrubland
and Grassland apply — the island's only two biomes. Method: reuse each species' own mainland
per-biome percent-table value as a baseline, scaled by a colonization-strength multiplier per
tier — **Strong = 90%, Moderate = 50%, Low-moderate = 25%, Excluded = 0%** — e.g. Treefox
(Low-moderate) goes from 100%/100% mainland to 25%/25% island; Scattermouse (Moderate) goes from
100%/15% to 50%/8%.

Two species don't fit that formula and get a direct value instead, flagged explicitly rather than
forced through the multiplier:

- **Mudlizard** — its mainland baseline (15%/0%) is itself a weak proxy, a known consequence of
  being a hydrology-driven species modeled on a biome-distance grid (same issue already flagged
  for the mainland table). Checked the island's own stream network directly rather than scaling a
  bad number: 99.1% of Woodland/Shrubland and 88.9% of Grassland sit within 1 km of a mapped
  stream (mean distance 357 m vs. 540 m) — both genuinely well-watered, Woodland/Shrubland modestly
  more so. Locked at Woodland/Shrubland 80% / Grassland 65%, a direct value grounded in that check,
  not 90%-of-15%.
- **Rustowl — Island ecotype** — a reversed-origin case, not a colonized fringe: the island is
  this ecotype's point of origin. Its two applicable named anchors (Woodland/Shrubland, Grassland
  — Temperate Forest doesn't exist on the island) carry over at full value, 100%/100%, no
  colonization discount. Mainland Rustowl is the smaller, later-founded population and keeps its
  own values in the mainland-only regular-bucket table, unaffected.

Clicksnake is now a row in this sheet directly (reclassified out of the high-threat bucket, see
below) — Grassland 100%/Woodland-Shrubland 15%, same as its mainland values, no colonization
discount, per its existing parity lock above.

**High-threat suitability composites — LOCKED this session, Reaper and Twinshadows only.**
Clicksnake reclassified out of this bucket entirely: its stated design constraint ("near-uniform
across Grassland") is exactly what the regular-bucket Specialist archetype already expresses — a
weighted composite wasn't adding anything a flat percent-table row couldn't already say, so it
moved to the regular-bucket table above instead of getting invented weights it doesn't need.

Both composites use a **hard biome mask**, not a weighted layer, for each species' absolute
habitat requirement — the source text says "confined to" and similar categorical language, not a
soft preference, so blending biome into the weighted sum the way Tappa 6 did for human-settlement
suitability would understate a constraint that's actually binary. This directly reuses Tappa 6's
own **exclusão** concept (`06_tappa6_suitability.md` §0): a 0/1 multiplier on the composite, not
an additive núcleo term — the same architectural role the Povo Silencioso buffer already plays
there, just gated by biome membership instead of distance-to-archipelago.

- **Reaper** — exclusão: Alpine Fellfield ∪ Alpine Tundra. Within that, núcleo layers: distance
  from Círculo/settlement infrastructure, farther = more suitable (45% — its single most
  load-bearing citation, "confined... by real anthropogenic range contraction," not preference);
  prey-base proximity, reusing the regular-bucket percent table's Cryburrow + Deergoat values as an
  input layer (35% — stated baseline diet); ruggedness/remoteness (20% — real large-carnivore
  precedent for marginal, rugged terrain once anthropogenic pressure pushes a species out).
- **Twinshadows** — exclusão: Temperate Forest (the "dense, low-light forest" habitat named for
  its melanism citation). Within that, núcleo layers: low-light proxy (40% — a loose stand-in for
  canopy density/"low-light," flagged honestly as terrain insolation, not vegetation cover; this
  pipeline has no real canopy-density layer); prey-base proximity, reusing Blacknose + Snaketail's
  percent-table values (35%); ruggedness/remoteness (25% — ambush cover for a stalk-and-pounce
  hunter).

**Verified against `06_tappa6_suitability.md`'s actual layer inventory** (checked once the device
reconnected — two corrections against the original draft of this section):

- **No `r.sun`/GRASS layer exists in this pipeline** — Tappa 6's solar module
  (`src/suitability/solar.py`) is a custom FAO-56 elevation-corrected clear-sky computation, built
  from scratch because GRASS wasn't available in this sandbox. Twinshadows' low-light proxy reuses
  `solar_suitability_annual_120m` (or the raw `annual_insolation_MJm2_120m`), inverted, not `r.sun`
  — same underlying caveat chain Tappa 6 already documented for that layer (clear-sky only, no
  cloud climatology, 120 m/16-direction horizon shading is coarse), on top of the new
  insolation-≠-canopy-cover caveat already flagged above.
- **`slope_suitability_120m` can't be reused directly for "ruggedness/remoteness."** That layer is
  tuned for human habitability — 1.0 at gentle grade, decaying to 0.0 by 30% — the opposite sense
  from what a predator favoring remote, rugged terrain needs. Both composites reuse the underlying
  raw `slope_pct_120m` field, but need a new suitability curve applied to it (rising with slope,
  not falling) rather than reusing `slope_suitability_120m`'s existing 0-1 output.

One incidental cross-check worth noting: Tappa 6's own water layer (`water_suitability_120m`)
came out "nearly flat in practice" (land mean 0.9907) because this world's stream network is dense
enough that distance-to-stream barely discriminates anywhere (82.5% of land within 0.5 km of a
stream) — independently corroborating this session's own Mudlizard finding (93.9-100% of every
biome within 1 km of a stream on the mainland, similarly dense on the island). Two separate checks,
same conclusion: this world's hydrology doesn't vary enough to use as a differentiating layer on
its own, for settlements or for fauna.

One dependency still open: "distance from Círculo/settlement" doesn't exist in the pipeline yet —
a new raster to compute (distance-to-nearest-Círculo-point), not something Tappa 6 already
produced.

**Alpine resource outposts — domesticated-animal bucket's third primitive, LOCKED this session
(method only, coordinates NOT YET COMPUTED — see below).** Grew out of a real gap: yak (alpine/
Spine niche) had no Círculo to attach to under the original point-per-Círculo model, since none of
the real 17 sites sit in an alpine biome at all (Tappa 6's site-selection algorithm structurally
never favors that terrain). Resolved by dropping yak from the original-17 lookup entirely — its
lore-established "herders already living in the range" (§16, Reaper) was never actually tied to
the 17 mapped Círculos in the first place — and giving it a new settlement type instead.

**One settlement type, not two.** Small, closed outposts sited near the Spine for resource
extraction (the same schist-hosted gold-quartz-mica veins already locked in the geology/Vértice
sections above) and Vértice study alike — extraction-vs-protection is a live social/political
debate *among the people who live there*, staying entirely off the map; a single outpost can house
both purposes, or people with either purpose can share the same site. No spatial split needed.

**Suitability composite, three núcleo layers, no exclusão beyond the standard land/lake hygiene
already inherited from Tappa 6:**

- Resource proximity (~40%) — continuous distance decay from each ridge's own axis, reusing the
  exact geometry the geology section already locked (`terrain_ridges.geojson`, the same
  `falloff_km`-based schist radius used for jade/Vértice-material siting). Deliberately soft, not a
  hard cutoff — a site just outside peak-grade terrain but much safer to build on is a realistic
  real-world mining-camp pattern, so this trades off against the other two layers rather than
  gating them.
- Slope (~25%) — reuses `slope_suitability_120m` directly, in its native human-buildability sense
  (gentle = suitable) — unlike the Reaper/Twinshadows composites above, no inverted curve needed
  here.
- Climate mildness (~35%) — reuses `biotemperature_c` (Tappa 5) as a survivability proxy within the
  alpine band: warmer within that band = more livable, avoiding the coldest ground nearest
  permanent snow.

**Site generation**: local suitability maxima within each confirmed alpine/subalpine massif (same
connected-component segmentation already run for the biome-fragmentation check), not a dense
placement algorithm — a handful of specific candidates, matched in scale to Reaper's own "4 known
family packs," not saturating the range.

**Status — active / abandoned / temporary refuge — is a separate, later authorial decision, NOT
derived from the composite.** The composite only answers where a viable site could be; which
candidates are currently inhabited, historically abandoned, or seasonal/temporary refuges is
assigned by hand afterward, per site — mirrors Tappa 6's own relationship between its suitability
surface (where) and its post-hoc architecture-style lookup (what kind). Flagged risk, not a flaw:
two candidates with near-identical suitability scores could end up with very different fates (one
active, one abandoned) since history doesn't purely follow terrain quality — worth remembering when
status actually gets assigned so it doesn't read as arbitrary.

**Feral fauna ties to abandoned instances specifically.** Yak, reindeer, and sheep can all go feral
at an abandoned outpost (real citation: feral cattle/horses/goats/sheep are all well-documented —
Camargue horses, Chillingham cattle — feral dogs/chickens sustaining independently of human refuse
is a much weaker real pattern, so those two are excluded). Reindeer gets no dedicated ecotype for
this — unlike yak (pack animal) and sheep (already has an established Terrapedra mountain ecotype),
reindeer has no existing alpine-adjacent reason to be there; its feral presence is simply
opportunistic historical keeping at a since-abandoned site, not a distinct locked variant.

**Narrative payoff, LOCKED**: abandonment ties directly to Reaper's own anthropogenic-range-
contraction citation — sustained predation pressure on an outpost is the mechanism, reinforcing
existing lore (mutual human/Reaper retreat from the same marginal frontier) rather than adding an
unrelated cause. Made concrete via the status-derivation method below rather than left abstract.

**Outpost status (active / temporary refuge / abandoned) — LOCKED as derived from a Reaper-
exposure overlay, not freehand.** Not sampled as Reaper's raw suitability *value* at each outpost
— outposts need buildable slope and mild climate to exist at all, which puts most candidates just
below the harshest terrain, outside Reaper's hard exclusão mask (Alpine Fellfield ∪ Alpine
Tundra) entirely; sampling suitability there would mostly read 0 and give no differentiation.
Instead: **distance from each outpost to the nearest Reaper pack territory** (see massif
assignment below) is the exposure signal. Three tiers, each with a real mechanism, not just three
bins on one gradient: **active** (low exposure, occupied year-round); **temporary refuge**
(moderate exposure, occupied most of the year but evacuated specifically during Reaper's real
snowmelt-triggered descent window, late winter/early spring — the same real transhumance pattern
of seasonally dodging a known danger window, not an invented mechanic); **abandoned** (exposure
became untenable entirely). Kept as a strong prior, not a hard threshold — real abandonment
history has causes beyond predation exposure (resource depletion, economic shifts), and forcing
every status purely off a distance number would trade one kind of arbitrariness for another;
exposure informs the authorial call, doesn't automate it.

**Reaper pack-to-massif assignment — LOCKED, closing this open item.** Checked against the real
fragmentation numbers rather than assumed: Reaper's actual habitat mask (Alpine Fellfield ∪ Alpine
Tundra, not Permanent Snow & Ice) has only **two** real massifs, not four. Alpine Fellfield is
essentially unfragmented (395.8 of 396.0 km² in one patch); Alpine Tundra is the one that splits —
480 km² on the main spine, 136.5 km² on the South Branch, ~40 km apart. Main spine ≈876 km²
combined, South Branch 136.5 km² — a real ~6.4:1 area ratio. Packs split unevenly to match:
**South Branch = 1 pack ("the South family," the only pack already named in the text — a small,
isolated, area-constrained massif is exactly where one unstable pack's behavior would be most
locally disruptive, matching its already-established leadership-change hook non-coincidentally);
main spine = the other 3 packs**, sub-territories deliberately left without further GIS precision
— the text already frames pack ranges as "semi-mapped... informally tracked," so manufacturing
exact boundaries there would overspecify what the lore never wanted precise.

**South Branch carries extra exposure weight right now, not equal weight.** The leadership-change
hook means an outpost equidistant from both massifs isn't equally exposed at this moment in the
story — South Branch is currently the sharper danger. Outpost-status exposure = distance to
nearest massif, with a South-Branch-specific multiplier on top reflecting the current instability.

Not yet computed: both composites (outposts and Reaper) need real data staged and run —
`slope_suitability_120m`, `biotemperature_c`, and Reaper's own suitability layers aren't staged
yet (only `terrain_ridges.geojson` is). Method is fully locked for both; coordinates are not.

**Still open, not yet decided:**

- Domesticated-animal Círculo-point/Povo-Livre-zone/alpine-outpost rendering — method LOCKED for
  all three primitives, none actually computed yet (needs `slope_suitability_120m` and
  `biotemperature_c` staged for the outpost composite).
- Alpine outposts: exact number/location of candidates (needs the composite actually run), and the
  active/temporary-refuge/abandoned status per candidate — exposure-derivation method LOCKED
  (distance to nearest Reaper pack territory, South-Branch-weighted), still needs both composites
  actually computed to produce real numbers.

**Regular/resident bucket — species × biome percent table, LOCKED and delivered
(`tappa7_fauna_biome_percent.xlsx`, sent to Nico, not committed — a working spreadsheet
deliverable, same status as `tappa7_fauna.xlsx`).** Closes the "still open, not yet decided"
percentages item above for this one bucket.

Method: every biome sits at a real (belt, moisture-tercile) coordinate already locked in
`05_tappa5_biomes.md` §5 — biome-distance = belt-step difference + moisture-tercile-step
difference (moisture only applies within Boreal/Cool Temperate). Permanent Snow & Ice = 0% for
every species in this bucket; no resident regular fauna is placed there. One archetype
classification per species (not per species-per-biome), read from each species' own written
description: **Specialist** (one or two named core biomes = 100%, distance-1 neighbors = 15%,
else 0% — Grassmothers, Blacknose, Tailstand, Flashfrog, Scattermouse, Snaketail + Furypack
[dual: Temperate Forest/Grassland ecotone], Cryburrow + Deergoat [dual: Alpine Fellfield/Alpine
Tundra]); **Named multi-biome** (text explicitly lists 2-3 biomes, each = 100%, same decay —
Treefox: Grassland/Woodland-Shrubland/Temperate Forest); **Explicit generalist** (text's own
words say "wide-ranging"/"not confined to a single biome" — flat baseline across all land
biomes, 0% Permanent Snow & Ice: Quillhog 80%, Farsmell 80%, Meatcleaner 70%).

Three gaps resolved, none invented unilaterally:

- **Tailstand** was missing from the bucket-1 species list entirely — added as a standard
  Specialist (Woodland/Shrubland). Its island presence stays separate and authorial (human
  introduction, already locked above), not this join.
- **Mudlizard** ("riverbank/wetland" in the text — a hydrology dependency, not a biome one).
  Checked directly against `stream_mask.npy` (distance-to-nearest-stream per biome) rather than
  assumed: streams turned out near-uniformly dense across every biome (93.9-100% of area within
  1 km of a stream everywhere), so stream proximity doesn't actually differentiate biomes here —
  a real finding that overturned the originally planned method, not a confirmation of it. Pivoted
  to Named multi-biome using the two warm/wet biomes (Subalpine Wet Forest, Temperate Forest) as
  the ecologically honest anchor instead.
- **Trinketbird** has no stated biome — its ecology is Círculo-infrastructure-driven (targets
  Vértice crystals/metal fittings), not habitat-driven. Nico's resolution: weighted by the real
  distribution of the 17 placed Círculos (13 Grassland/2 Temperate Forest/2 Woodland-Shrubland) —
  Grassland 100%, Temperate Forest 15%, Woodland/Shrubland 15%, 0% elsewhere. A proxy, not the
  literal mechanism — a true point-distance-to-Círculo model belongs with bucket 3's machinery,
  not this join.
- **Rustowl (mainland ecotype)** has no stated biome or prey base. Nico's call: a flat generalist
  baseline would undercut the "owl choir" rarity/mystery framing (§7 above), so anchored instead
  to the small-mammal niche's own biomes as an implied prey base — Named multi-biome, same three
  biomes as Treefox (Grassland, Woodland/Shrubland, Temperate Forest).

Out of scope for this table: the high-threat bucket, domesticated animals, and coastal/marine +
migratory species (each its own bucket's machinery, per "Map structure" above), and every
species' mainland/island split (separate model, already locked species-by-species above).

## 8. Open follow-ups

Direction resolved on all of these; none are implemented in code yet. Per the restructuring above,
items from sections 1-6 are decision-only for now — each becomes its own future Tappa, number not
yet assigned, once its blocking layer is built. The fauna item is different: its map is this
stage's actual open deliverable, not deferred to a future Tappa.

- Lithology layer (now including its overlap precedence rule), cave layers, and the geothermal
  hazard: direction fully set, none implemented in code yet. Cartographic output deferred to a
  future Tappa (number TBD) once the lithology classification is actually computed.
- Sedimentary basin-fill lithology class: RESOLVED AND VERIFIED. Canterbury Plains glacial-outwash
  analog checked directly against Tappa 4's `stream_mask.npy` — all six plains/plateau zones are
  genuinely threaded with stream network (mean distance to nearest stream 0.25-0.72 km, 84.5-99.7%
  of each zone's area within 1 km of a mapped stream). See geology section above for the table.
- Círculo–Povo Livre frontier tension: mapping approach set (Povo Livre suitability >80% vs.
  Círculo territory/infrastructure proximity) — deferred to a future Tappa alongside the rest of
  the conflict-zone/dangerous-creatures layer.
- Dangerous-seas layer: navigability model set (designated corridors, not open-sea routing) —
  deferred to a future Tappa (number TBD) once the hazard layer exists.
- Road path reconstruction (predecessor extraction) and the rail grade-ceiling cost function: not
  built. Both need a biome-differentiated cost multiplier (Wet Forest costs more than Grassland) —
  not built. Deferred to a future Tappa (number TBD) alongside the rest of the
  transportation/infrastructure layer.
- Kite-transport wind-shadow mask and directional route labeling: scope simplified (qualitative
  speed, not a resolved wind field) — deferred with the rest of transportation/infrastructure.
- **Fauna distribution/range maps — this stage's actual open deliverable, not deferred.** Scoped at
  the end of §7 above; exact map content still to be decided.
- Terracota architecture style has no assigned biome despite covering the majority of actually
  placed Círculos (13/17 sit in Grassland) — deliberately deferred by Nico to Tappa 8, where
  architecture will be revisited in full. Tappa 8's scope/position is unaffected by this session's
  restructuring.
- Roadmap sequencing for the next scenario-depth passes, in the order Nico wants them: (1) close
  out geology (this document) — DONE, (2) develop Povo Silencioso further — DONE, (3) Vértice
  types/mechanics — DONE, (4) other Vértice-adjacent materials beyond jade/pounamu — DONE, see
  geology section above and `scenario_reference.md` for the full domain/verb/tier system.
- Alpine terrestrial apex predator: **RESOLVED as Reaper** — the last open item from the original
  dangerous-creatures seed list, closing the entire fauna build-out initiative. Ordinary wildlife
  (not a Visitante), full design in `scenario_reference.md` §16. The aerial "sky" danger side of
  that same seed list is resolved separately — see `scenario_reference.md`'s public framing; the
  full design (which is Visitante-tied) lives only in the GM-only secret file, never this repo.
