# New rock types (marble, sedimentary limestone, granite) — geological facts for the Vértice-nature re-evaluation

Prepared by the Tappa 8 chat, 2026-08-20, for the Scenario chat. Purpose: Nico wants to
revisit the "which real mineral/device property best embodies each Vértice domain's
nature" discussion from Tappa 7, now that lithology v6 adds three rock types that didn't
exist when that discussion happened. This document is a **factual handoff, not a
decision** — geology only. Which domain (if any) a fact below should actually be assigned
to is a Scenario-chat call; the original tuner-selection reasoning from Tappa 7 isn't
reproduced here because this chat doesn't have full access to it, only to its final
output (`resources.py`'s `VERTICE_MATERIALS` table).

## Current domain → material assignments (recap, for context)

| domain | current material/rock | real property cited | status |
|---|---|---|---|
| Onda | schist quartz veins | birefringence | locked |
| Energia | schist mica (same veins) | piezoelectric | locked |
| Matéria | greywacke laumontite (zeolite) | molecular-sieve structure | locked |
| Bios | basin_fill vivianite | tied to organic decay (anoxic floodplain) | locked |
| Campo | volcanic magnetite | paleomagnetism | locked |
| Mente | volcanic silver/copper pair | a "tuning" pair — refinement method still "under study" in-world | locked, unrefined |
| Tempo | — none | — | **no known access method at all** |

## Marble

**Composition**: recrystallized calcite (CaCO₃), sometimes with dolomite admixture — same
base mineral as sedimentary limestone below, but metamorphosed into an interlocking
crystalline fabric (Ordovician Takaka-terrane analog, already cited in the decision doc).
Mohs hardness ~3-4.

**Optical**: marble's recrystallized fabric scatters light a few millimeters into the
stone (real, well-documented subsurface scattering — the property classical sculptors
exploited for the "living" look of Carrara marble). This is a translucency/diffusion
effect, NOT the same thing as birefringence (see calcite note below) — worth keeping the
two properties distinct if either comes up.

**Chemical**: reacts readily with dilute acid (the standard field test for carbonate
rock) — this reactivity is also literally what drives karst dissolution (S8c).

**Piezoelectric — an important negative fact**: calcite is NOT piezoelectric. Its crystal
structure has a center of symmetry, which piezoelectricity requires the ABSENCE of
(quartz lacks that center of symmetry, which is why it works for Energia). If a future
discussion considers a calcite-family material for Energia specifically, that citation
would not hold up the way it does for quartz — flagging this now so it doesn't get
assumed by default.

## Sedimentary limestone

**Composition**: same calcite base as marble, but NOT recrystallized — retains original
depositional texture, often bioclastic (contains fossils/shell/skeletal fragments —
literally compressed ancient marine organisms). Real citation already in this project:
Oamaru stone, which is exactly this kind of bioclastic limestone.

**A genuinely strong candidate fact for Onda specifically**: clear calcite crystals
("Iceland spar") are one of the most famous naturally birefringent materials in
science — their double refraction (Δn ≈ 0.17) is roughly 15-20× stronger than quartz's
(Δn ≈ 0.009), which is why clear calcite, not quartz, is the textbook material for
demonstrating birefringence (Nicol prisms; the Viking "sunstone" navigation legend). If
Onda's nature is meant to track birefringence as a real physical property rather than
being tied to schist/quartz specifically, calcite is a stronger real-world citation than
the current one.

**Important caveat on that, so it isn't overclaimed**: this would need the SAME framing
already used for schist — it's not that "limestone the bulk rock is birefringent," it's
that clear, optical-grade calcite crystals occur in cavities/veins WITHIN or associated
with carbonate deposits (a real, well-documented occurrence, same structural logic as
"quartz veins within schist," not a new kind of claim). And per the marble note above,
this source would give birefringence but NOT piezoelectricity — if Onda and Energia were
ever to be re-sourced together from one material the way they currently are (both from
schist's veins), calcite doesn't reproduce that pairing.

**Porosity**: more porous/permeable than marble (retains original grain structure) —
consistent with why sedimentary_limestone and marble are both karst-forming but the
project treats them as textural variants of one dissolution story (S8c), not two.

## Granite

**Composition**: coarse-grained igneous rock — quartz + feldspar (orthoclase/plagioclase)
+ mica (biotite/muscovite). Real citation already in this project: Coromandel granite,
dimension/monumental stone. Mohs ~6-7.

**Shares minerals with schist, worth noting explicitly**: granite contains real quartz and
real mica, the same two minerals schist's veins are cited for (Onda/Energia). This is a
genuine overlap, not a coincidence to ignore — if Onda/Energia's "nature" is really about
those MINERALS rather than about schist as a rock unit specifically, granite is a
legitimate second source, not a competing claim.

**A genuinely strong candidate fact for Tempo specifically — the one domain with no
material at all**: granite is well known to carry the highest natural background
radioactivity of any common rock type, from trace uranium, thorium, and potassium-40
concentrated in accessory minerals (zircon, monazite, potassium feldspar) — the same real
phenomenon behind granite countertop radon-gas discussions. Radioactive decay is a literal
physical clock (a fixed, measurable half-life) — if Tempo's nature is meant to track an
intrinsic, unstoppable, measurable rate of change, this is about as direct a real-world
analogy as exists in common rock. Flagging this as the single most interesting new-rock
possibility from this batch, precisely because Tempo currently has nothing.

**Magnetism — a negative check, for completeness**: granite has low magnetic
susceptibility relative to volcanic basalt's titanomagnetite content. Doesn't compete with
Campo's existing volcanic assignment; checked so this document isn't just cherry-picking
positive matches.

## Supplementary reference (already built, not about "nature" but may be useful context)

Two relative indices already exist per rock class, from this chat's other work this
session — not mechanically about Vértice tuning, but potentially useful texture/flavor
reference if a domain's "how it feels to work with" matters narratively:

| class | travel friction (§8f, ≤1.0=slower) | excavation effort (§8g, ≥1.0=harder) |
|---|---|---|
| basin_fill | 1.00 | 1.0 |
| sedimentary_limestone | 0.60 (worst — karst hazard) | 1.3 (soft, easy to cut) |
| marble | 0.60 (worst — karst hazard) | 1.6 (harder than limestone, still soft) |
| greywacke / schist | 0.85 | 2.3 |
| volcanic | 0.70 | 2.6 |
| granite | 0.90 (mild — good footing) | 3.0 (hardest — no foliation weakness) |

Note the granite/marble inversion again here: granite is easy underfoot but the hardest
to quarry; marble/limestone are hazardous underfoot but comparatively easy to cut once
you're there.

## What this document is not

No conclusions about which domain should end up with which material — that's the
Scenario chat's call, informed by the fuller Tappa 7 "nature of each Vértice" reasoning
this chat doesn't have. This is geology and citation-strength notes only.
