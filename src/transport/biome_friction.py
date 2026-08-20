"""
Tappa 9 -- biome-differentiated travel-friction multiplier for
cost_distance.py's Tobler-hiking-function graph.

Same architectural role as geomorphology/lithology.py's LAND_TRAVEL_FRICTION
(Tappa 8 S8f): a per-cell 0-1-ish multiplier on LAND-LAND edge speed, passed
to build_cost_graph's `friction_multiplier` parameter -- but keyed on
biome_id (vegetation/ground cover) instead of lithology (rock type). The two
are independent physical properties (what's growing on the ground vs. what
the ground is made of underneath) and are meant to be COMBINED
multiplicatively by the caller, not to replace one another -- see
run_tappa9_road_network.py for how they're actually combined into one
friction field.

Direction was already decided before this module existed, not invented here:
`07_tappa7_regional_scenario.md` S6 and `scenario_reference.md` S18 both
state, as a locked narrative fact, that Wet Forest costs more to build/
maintain a road through than Grassland -- neither ever attached a number to
that. This module supplies numbers, for every land biome, not just the one
named Wet-Forest-vs-Grassland contrast.

Values are grounded in this project's own
`docs/reference/biome_landscape_characteristics.md` (vegetation structure
per biome, already computed from this world's real Tappa 5 data) plus real
off-track-travel accounts for the closest matching real vegetation
structure, not a cited per-biome hiking-speed dataset (none exists).
**UNREVIEWED first-pass estimates, same status as LAND_TRAVEL_FRICTION --
pending Nico's sign-off, not written to config/parameters.yml.**
"""
from __future__ import annotations

import numpy as np

__all__ = ["BIOME_TRAVEL_FRICTION", "biome_friction_multiplier"]

# Keyed on biome_id (src/biomes/world_biomes.py's BIOME_NAMES index).
# 1.00 = Lowland Steppe/Grassland baseline (open, treeless, minimal
# clearance) -- same "cheapest class is the anchor" convention
# LAND_TRAVEL_FRICTION uses for sedimentary_basin_fill. Biome 0 (Ocean) is
# deliberately absent -- sea edges never consult this table at all (see
# biome_friction_multiplier's docstring).
BIOME_TRAVEL_FRICTION: dict[int, float] = {
    1: 0.35,  # Permanent Snow & Ice -- crevasse/avalanche hazard; real
              # glacier travel is dramatically slower than ordinary hiking,
              # often needing roping and crevasse-probing even on gentle
              # slope, independent of whatever the DEM slope term already
              # charges for elevation
    2: 0.90,  # Alpine Fellfield -- open, no canopy to clear, but loose
              # scree/talus underfoot slows travel; a footing penalty only,
              # the mildest of the three alpine/subalpine classes
    3: 0.75,  # Alpine Tundra -- open, no canopy, but real NZ snow-tussock
              # (Chionochloa, 0.5-1.5 m, per biome_landscape_
              # characteristics.md) is well documented as slow, tiring
              # walking despite reading as "open" from a distance
    4: 0.55,  # Subalpine Wet Forest -- closed canopy, heavy epiphyte/moss
              # load, damp low-light understory (biome_landscape_
              # characteristics.md); real NZ off-track native-bush travel
              # is notoriously slow (commonly well under 1-1.5 km/h)
    5: 0.75,  # Subalpine Woodland -- same temperature band as Wet Forest,
              # but a visibly more open canopy and drier/sparser understory
              # (biome_landscape_characteristics.md) -- transitional value
    6: 0.85,  # Subalpine Dry Scrub -- no closed canopy, low sparse woody
              # scrub with exposed ground/rock between clumps; real
              # matagouri scrub is thorny in patches but not canopy-dense,
              # so only a mild penalty
    7: 0.55,  # Temperate Forest -- the tallest, wettest, most rainforest-
              # like class (real podocarp-broadleaf structure, rimu/
              # kahikatea heights); SAME value as Subalpine Wet Forest on
              # purpose -- both are closed-canopy, heavy-understory forest,
              # and nothing in this project's sources differentiates their
              # off-track difficulty specifically (mirrors
              # LAND_TRAVEL_FRICTION's own schist=greywacke tie for the
              # identical "nothing distinguishes them" reason)
    8: 0.90,  # Woodland/Shrubland -- open mosaic, patches of lower-stature
              # forest/scrub interspersed with open ground, not a uniform
              # closed canopy (biome_landscape_characteristics.md)
    9: 1.00,  # Lowland Steppe/Grassland -- baseline/anchor, open treeless
              # tussock grassland, minimal clearance needed
}


def biome_friction_multiplier(
    biome_id: np.ndarray,
    multipliers: dict[int, float] = BIOME_TRAVEL_FRICTION,
    default: float = 1.0,
) -> np.ndarray:
    """Map each biome_id cell to a travel-friction multiplier. Cells whose
    biome code isn't a key in `multipliers` -- biome_id == 0 (Ocean/
    non-land) above all -- fall back to `default` (neutral). This is
    harmless by construction, not just by luck: `build_cost_graph` only
    ever applies `friction_multiplier` to LAND-LAND edges (see its own
    docstring), and every land cell in this project's biome_id carries a
    real 1-9 code (0 is reserved for the caller's land_mask, per
    world_biomes.py) -- so a non-1.0 default would never even be read for
    an ocean cell; the fallback exists purely so this function never
    raises on an out-of-table code, same convention as
    `travel_friction_multiplier` in geomorphology/lithology.py.
    """
    out = np.full(biome_id.shape, default, dtype=np.float32)
    for code, mult in multipliers.items():
        out[biome_id == code] = mult
    return out
