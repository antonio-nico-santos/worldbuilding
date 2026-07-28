"""
World-specific biome layer for Tappa 5 (fine-tuning on top of `holdridge.py`).

Where `holdridge.py` is deliberately world-agnostic (it would work unchanged
on any other project's climate stack), everything in this module is specific
to THIS domain, and encodes decisions made by looking at what this world's
own data actually produces -- not an off-the-shelf scheme imported wholesale.
See docs/decisions/05_tappa5_biomes.md for the full record, including the
two rejected drafts. Summary of what's fine-tuned here and why:

1. **Moisture axis is rebinned, not read off Holdridge's stock scale.**
   Applying Holdridge's own global PET-ratio bands to this world's climate
   stack finds that literally no land cell clears the Semiarid threshold
   (ratio >= 4) -- the driest cell on record (324 mm/yr, the same cell
   Tappa 2-4 validated against real Alexandra's ~300 mm) still only reaches
   ratio ~2.1, "Subhumid". This is structural, not a bug: Holdridge's PET is
   linear in biotemperature, and this world's biotemperature never exceeds
   12 deg C, capping PET near 700 mm/yr even at the driest, lowest-elevation
   spot. Used at face value, Holdridge would erase the entire windward/
   leeward moisture contrast Tappa 2-4 spent three stages validating (5.6:1
   precipitation ratio, 1136 m ELA differential, 5.9:1 discharge ratio).
   Fix: the PET ratio is kept as the underlying quantity (it still correctly
   *ranks* wet vs. dry cells) but binned by this world's own land-cell
   percentiles -- the same move Tappa 3/4 already made for their windward/
   leeward split -- rather than Holdridge's desert-calibrated absolute bins.
2. **3 moisture tiers (terciles), not the quartiles first tried.** A first
   draft used quartiles (Superwet/Wet/Moist/Dry) and surfaced two problems:
   a 96 km^2 "Temperate Rainforest" sliver (Cool Temperate x Superwet) that
   was visually unreadable at map scale, and a 4-way "which shade of green"
   palette collision flagged by the project's color-validation tooling
   (worst-case normal-vision Delta-E 4.2, i.e. hard to tell apart even with
   full color vision). Folding to 3 tiers removes both: the old Superwet
   category merges into the Wet tier (so the rainforest sliver folds into
   "Temperate Forest" rather than disappearing as noise), and the total land
   biome count drops from 10 to 9.
3. **Polar and Subpolar belts are NOT split by moisture.** Both are cold
   enough that a "dry polar" class barely exists on this domain (checked
   directly: 0 km^2 of Polar land falls in the driest tercile) and the areas
   involved are modest -- splitting them would only add categories without
   adding legible distinctions. Warm Temperate (four cells, ~0.06 km^2,
   right at the 12 deg C belt edge) is folded into Cool Temperate for the
   same reason.
4. **Permanent snow is overridden from Tappa 3's mass-balance mask, not
   read off Holdridge's Polar belt.** Checked directly: Holdridge's own
   Polar belt (biotemperature < 1.5 deg C) is 1243 km^2, ~30% larger than
   Tappa 3's physically-modelled permanent-snow area (960 km^2, 68% cell
   overlap) -- the identical naive-temperature-threshold error Tappa 2/3
   already caught once (naive permanent-snow: 108 km^2 vs. mass-balance:
   960 km^2), now showing up again if Polar were used directly as the
   ice/nival class.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .holdridge import HoldridgeResult, classify

__all__ = [
    "BIOME_NAMES",
    "BIOME_COLORS_HEX",
    "MOISTURE_NAMES",
    "WorldBiomeResult",
    "classify_world_biomes",
]

# index -> name. 0 is reserved for ocean/non-land, filled in by the caller's
# land mask rather than by this module (this module only ever writes >=1
# into biome_id for land cells).
BIOME_NAMES = [
    "Ocean",                       # 0 (not assigned here -- caller's land_mask)
    "Permanent Snow & Ice",        # 1
    "Alpine Fellfield",            # 2  Polar, unsplit by moisture
    "Alpine Tundra",               # 3  Subpolar, unsplit by moisture
    "Subalpine Wet Forest",        # 4  Boreal x Wet tercile
    "Subalpine Woodland",          # 5  Boreal x Moist tercile
    "Subalpine Dry Scrub",         # 6  Boreal x Dry tercile
    "Temperate Forest",            # 7  Cool Temperate x Wet tercile (absorbs the old Superwet/"rainforest" sliver)
    "Woodland / Shrubland",        # 8  Cool Temperate x Moist tercile
    "Lowland Steppe / Grassland",  # 9  Cool Temperate x Dry tercile
]

# Colors chosen for a natural, continuous wet-green -> dry-tan cartographic
# read (matching the pre-project plan's "qualitative palette, nominal data"
# brief for the biome layer) rather than a fully hue-spread categorical set.
# Run through the project's dataviz color validator (see decision doc S6):
# FAILS strict all-pairs CVD safety (worst pair, Temperate Forest vs.
# Subalpine Wet Forest, Delta-E 5.6 normal-vision -- both are, deliberately,
# shades of green representing a forest-density gradient). Mitigation is the
# pre-project plan's own click-to-select detail panel + always-visible
# legend (00_pre_project_planning.md), so no reading of the map ever depends
# on color alone -- documented here rather than silently accepted.
BIOME_COLORS_HEX = [
    "#bcdcee",  # 0 Ocean
    "#f5f7f8",  # 1 Permanent Snow & Ice
    "#8f8579",  # 2 Alpine Fellfield (gray-brown rock)
    "#a89bb0",  # 3 Alpine Tundra (heather -- deliberately NOT in the yellow-green family)
    "#1f6f54",  # 4 Subalpine Wet Forest
    "#4f8a5c",  # 5 Subalpine Woodland
    "#b8622e",  # 6 Subalpine Dry Scrub (rust -- deliberately NOT green)
    "#2e7d46",  # 7 Temperate Forest
    "#7a9c4a",  # 8 Woodland / Shrubland
    "#e0a83f",  # 9 Lowland Steppe / Grassland
]

MOISTURE_NAMES = ("Wet", "Moist", "Dry")  # 0, 1, 2 -- ratio terciles, ascending (0=wettest)


@dataclass
class WorldBiomeResult:
    biome_id: np.ndarray                  # int8, 0=ocean/non-land
    holdridge: HoldridgeResult
    moisture_idx: np.ndarray              # 0..2 into MOISTURE_NAMES
    moisture_tercile_edges: tuple[float, float]
    permanent_snow_mask: np.ndarray
    notes: dict = field(default_factory=dict)


def classify_world_biomes(
    temp_c_monthly: np.ndarray,
    precip_mm_monthly: np.ndarray,
    land: np.ndarray,
    permanent_snow_mask: np.ndarray,
) -> WorldBiomeResult:
    """Full Tappa 5 classification: general Holdridge fields -> this world's
    rebinned moisture tercile -> the 9-class scheme above -> permanent-snow
    override. All four inputs are Tappa 2/3's own working-resolution (120 m)
    arrays; shapes must match."""
    annual_precip = precip_mm_monthly.sum(axis=0)
    hr = classify(temp_c_monthly, annual_precip)

    ratio_land = hr.pet_ratio[land]
    edges = np.percentile(ratio_land, [33.333, 66.667])
    moisture_idx = np.searchsorted(edges, hr.pet_ratio, side="right").astype(np.int8)

    belt = hr.belt_idx.copy()
    belt[belt == 4] = 3  # Warm Temperate (negligible, ~0.06 km2) folds into Cool Temperate

    biome_id = np.zeros(temp_c_monthly.shape[1:], dtype=np.int8)
    biome_id[land & (belt == 0)] = 2                                     # Polar -> fellfield
    biome_id[land & (belt == 1)] = 3                                     # Subpolar -> tundra
    biome_id[land & (belt == 2) & (moisture_idx == 0)] = 4               # Boreal Wet
    biome_id[land & (belt == 2) & (moisture_idx == 1)] = 5               # Boreal Moist
    biome_id[land & (belt == 2) & (moisture_idx == 2)] = 6               # Boreal Dry
    biome_id[land & (belt == 3) & (moisture_idx == 0)] = 7               # Cool Temperate Wet
    biome_id[land & (belt == 3) & (moisture_idx == 1)] = 8               # Cool Temperate Moist
    biome_id[land & (belt == 3) & (moisture_idx == 2)] = 9               # Cool Temperate Dry
    biome_id[permanent_snow_mask & land] = 1                             # override, applied last

    return WorldBiomeResult(
        biome_id=biome_id,
        holdridge=hr,
        moisture_idx=moisture_idx,
        moisture_tercile_edges=(float(edges[0]), float(edges[1])),
        permanent_snow_mask=permanent_snow_mask & land,
    )
