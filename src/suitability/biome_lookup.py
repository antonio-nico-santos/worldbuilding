"""
Tappa 6 -- biome suitability lookup, the 6th nucleo criterion.

This is the one nucleo layer that is NOT derived from a continuous physical
field (unlike slope/distance-to-stream/TWI/solar, which are all measurements
or physically-modelled quantities): it is a hand-authored 0-1 score per
biome class, an explicit value judgement about what a Circulo settlement
economy needs, not something that can be "computed" from the terrain. See
docs/decisions/06_tappa6_suitability.md (or wherever this stage's decision
doc lands) for the discussion; the short version is recorded in BIOME_NOTES
below so the rationale travels with the code.

What this layer measures, and why it isn't redundant with the other nucleo
layers: TWI (agriculture_suitability_120m) is a LOCAL topographic wetness
proxy (soil moisture from upslope drainage); the biome's own moisture axis
(Holdridge PET ratio, see src/biomes/world_biomes.py) is a REGIONAL climate
quantity (precip vs. potential evapotranspiration) -- a different physical
mechanism, not a duplicate measurement. What only biome_id carries is (a)
the regional thermal regime / growing-season length, and (b) vegetation
STRUCTURE (dense forest vs. scrub vs. bare rock vs. ice), which stands in
for timber/biomass availability. Slope and solar exposure are already their
own layers, so biome scores below try not to re-encode "steep" or "shaded"
via biome class.

HONEST LIMITATIONS:
- The relative order of the three warm/dominant classes (Temperate Forest,
  Woodland/Shrubland, Lowland Steppe/Grassland -- together 58.5% of all
  land) was resolved by an explicit economic-base decision (open-field
  agriculture as the Circulos' primary subsistence mode), not derived from
  any data in this project. This is NOT a cosmetic choice: checked directly
  against Tappa 2's own climate output, Temperate Forest's land cells
  average 4936 mm/yr precipitation vs. Lowland Steppe/Grassland's 531 mm/yr
  (9.3x) and sit on opposite sides of the domain (mean X -17.5 km vs.
  +19.4 km, in a 130 km-wide domain) -- the two orderings point at two
  geographically disjoint "best" regions (0% cell overlap between the two
  classes, by construction), not a subtle reweighting of the same place.
  See BIOME_SUITABILITY_POVO_LIVRE below for the mirrored ordering used for
  the Povo Livre (forestry/foraging-centred, not open-field-agriculture-
  centred) rather than treating this as a single fixed truth.
- Permanent Snow & Ice and Alpine Fellfield are scored near-zero (0.00/0.05)
  but participate ONLY in the nucleo weighted sum, not as an exclusao
  multiplier -- a deliberate decision made *after* checking directly that
  slope does NOT reliably zero these cells out on its own (41% of Permanent
  Snow & Ice cells and 21% of Alpine Fellfield cells have
  slope_suitability_120m > 0.3, i.e. they sit on plateaus/ice-cap terrain,
  not just steep faces). This means a flat, well-watered, sun-exposed patch
  of permanent ice CAN still pull a moderate composite score once the other
  nucleo layers are combined -- known and accepted, not an oversight.
- All 9 scores are hand-set on a 0-1 scale with no independent calibration
  data (no real economy to check against), same status as slope_suitability's
  gentle/hard-limit knees in terrain_metrics.py.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from src.biomes.world_biomes import BIOME_NAMES

__all__ = [
    "BIOME_SUITABILITY",
    "BIOME_SUITABILITY_POVO_LIVRE",
    "BIOME_NOTES",
    "majority_filter_biome_id",
    "biome_suitability_from_id",
]

# index matches BIOME_NAMES / biome_id (0 = Ocean, unused -- caller's land_mask
# decides which cells get a value at all).
BIOME_SUITABILITY = [
    float("nan"),  # 0 Ocean -- not scored, masked out by land_mask downstream
    0.00,  # 1 Permanent Snow & Ice
    0.05,  # 2 Alpine Fellfield
    0.20,  # 3 Alpine Tundra
    0.55,  # 4 Subalpine Wet Forest
    0.45,  # 5 Subalpine Woodland
    0.30,  # 6 Subalpine Dry Scrub
    0.85,  # 7 Temperate Forest
    0.90,  # 8 Woodland / Shrubland
    1.00,  # 9 Lowland Steppe / Grassland
]

# Alternate scoring for the Povo Livre (isolated-by-choice, less structured
# than a Circulo -- see the scenario chat, not committed to this repo): only
# the two endpoints are mirrored (Temperate Forest <-> Lowland Steppe/
# Grassland), Woodland/Shrubland is left at 0.90 in the middle in both
# versions, i.e. a straight reflection, not an independently-reasoned
# 9-class table. The other 6 classes (ice/rock/tundra/subalpine belt) are
# UNEXAMINED for this variant -- kept identical to BIOME_SUITABILITY because
# nothing in the Forest-vs-Grassland discussion bears on them, not because
# they were checked and found to already fit a foraging/forestry economy.
BIOME_SUITABILITY_POVO_LIVRE = [
    float("nan"),  # 0 Ocean
    0.00,  # 1 Permanent Snow & Ice        (unchanged)
    0.05,  # 2 Alpine Fellfield            (unchanged)
    0.20,  # 3 Alpine Tundra               (unchanged)
    0.55,  # 4 Subalpine Wet Forest        (unchanged)
    0.45,  # 5 Subalpine Woodland          (unchanged)
    0.30,  # 6 Subalpine Dry Scrub         (unchanged)
    1.00,  # 7 Temperate Forest            <- was 0.85
    0.90,  # 8 Woodland / Shrubland        (unchanged, stays the "middle" class)
    0.85,  # 9 Lowland Steppe / Grassland  <- was 1.00
]

BIOME_NOTES = {
    "economic_base_decision": (
        "open-field agriculture assumed as the Circulos' primary subsistence "
        "mode -> Lowland Steppe/Grassland (1.00) > Woodland/Shrubland (0.90) "
        "> Temperate Forest (0.85). These 3 classes cover 58.5% of all land, "
        "so this single ordering choice dominates the layer's effect more "
        "than any of the other 6 classes combined."
    ),
    "povo_livre_variant_decision": (
        "BIOME_SUITABILITY_POVO_LIVRE mirrors the Forest/Grassland endpoints "
        "for a forestry/foraging-centred economy instead of open-field "
        "agriculture. Verified this is a real geographic fork, not a cosmetic "
        "tweak: Temperate Forest averages 4936 mm/yr precip vs. Grassland's "
        "531 mm/yr (9.3x), and their land cells sit on opposite sides of the "
        "domain (mean X -17.5 km vs. +19.4 km) with 0% overlap between the "
        "two classes. The score swing itself is bounded (+-0.15, affecting "
        "36.9% of land) since Woodland/Shrubland stays fixed at 0.90 as the "
        "middle class in both versions."
    ),
    "ice_rock_decision": (
        "Permanent Snow & Ice / Alpine Fellfield kept as nucleo-only (not "
        "also an exclusao multiplier), a decision made AFTER verifying that "
        "slope_suitability_120m > 0.3 on 41% / 21% of their cells "
        "respectively -- slope does not reliably zero these areas out on "
        "its own. Accepted trade-off, not an oversight."
    ),
}

assert len(BIOME_SUITABILITY) == len(BIOME_NAMES)


def majority_filter_biome_id(
    biome_id: np.ndarray, land_mask: np.ndarray, window: int = 3
) -> tuple[np.ndarray, float]:
    """Smooth single/few-cell biome_id speckle (many of this world's biome
    patches are single-digit-km2 fragments at the 120 m grid -- see
    tappa5_biomes_meta.json's fragmentation block, e.g. Lowland Steppe/
    Grassland has 209 patches with a 0.058 km2 MEDIAN patch size) before
    turning it into a siting-relevant score. Settlement suitability is a
    coarse decision; pixel-level "sale-e-pepe" noise in the class map has no
    business driving it.

    Majority vote is restricted to LAND neighbours only (ocean cells are
    never counted, so coastal land cells aren't pulled toward "ocean" and
    ocean itself is left untouched) within a `window` x `window` block
    (default 3x3, edge-replicate-equivalent boundary via 'nearest' mode --
    same "don't invent data past the edge" spirit as this project's other
    block-reduce helpers). Ties (including the original class itself tying
    the new majority) keep the ORIGINAL class rather than picking an
    arbitrary winner.

    Returns (smoothed_biome_id, fraction_of_land_cells_changed).
    """
    n_classes = len(BIOME_NAMES)  # 0..9
    kernel = np.ones((window, window))
    land_f = land_mask.astype(np.float64)
    land_neighbor_count = ndimage.convolve(land_f, kernel, mode="nearest")

    counts = np.zeros((n_classes,) + biome_id.shape, dtype=np.float64)
    for k in range(n_classes):
        mask_k = ((biome_id == k) & land_mask).astype(np.float64)
        counts[k] = ndimage.convolve(mask_k, kernel, mode="nearest")

    majority = np.argmax(counts, axis=0).astype(biome_id.dtype)
    max_count = counts.max(axis=0)
    orig_count = np.take_along_axis(
        counts, biome_id.astype(np.int64)[None, ...], axis=0
    )[0]
    tie_with_original = orig_count >= max_count  # keep original on any tie

    smoothed = np.where(land_mask & ~tie_with_original, majority, biome_id)
    changed = land_mask & (smoothed != biome_id)
    fraction_changed = float(changed.sum()) / float(land_mask.sum())
    return smoothed, fraction_changed


def biome_suitability_from_id(
    biome_id: np.ndarray, land_mask: np.ndarray, lut: list[float] | None = None
) -> np.ndarray:
    """Map (smoothed) biome_id -> a suitability score, NaN on ocean.

    `lut` defaults to BIOME_SUITABILITY (the Circulo/open-field-agriculture
    scoring); pass BIOME_SUITABILITY_POVO_LIVRE (or any other same-length
    list) to score a different population's economy against the same
    smoothed biome map."""
    lut_arr = np.array(BIOME_SUITABILITY if lut is None else lut, dtype=np.float64)
    suit = lut_arr[biome_id.astype(np.int64)]
    return np.where(land_mask, suit, np.nan)
