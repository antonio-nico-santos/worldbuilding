"""
Tappa 6 -- final weighted composite, combining the 5 nucleo layers into one
0-1 suitability index per population (Circulo / Povo Livre), then applying
the Povo Silencioso exclusion multiplier. See docs/decisions/ (Tappa 6) for
the full discussion; the short version is recorded here so the weights
travel with the code that uses them.

WEIGHTS_CIRCULO and WEIGHTS_POVO_LIVRE are both explicit VALUE JUDGEMENTS
about what each population's economy needs, exactly like BIOME_SUITABILITY
in biome_lookup.py -- there is no data-derived "correct" weighting, only a
documented, revisable choice.

Rationale (from the Tappa 6 planning chat):
- slope: reduced for Povo Livre (0.20 -> 0.10), NOT zeroed -- a temporary
  camp tolerates far more relief than a permanent Circulo foundation does,
  but rugged terrain still matters for mobility/safety even for a nomadic
  group, so this isn't dropped to zero either.
- water: low weight (0.05) for both, and this barely matters in practice --
  water_suitability_120m is empirically near-flat (land mean 0.991) because
  this world's stream network is dense; its weight mostly controls how hard
  the rare far-from-water outlier (<1% of land) gets penalised, not how the
  bulk of land ranks.
- agriculture (TWI proxy): high for Circulo (0.25, open-field-agriculture
  economy) vs. low for Povo Livre (0.10, foraging/hunting economy -- wetter
  land still means more biomass/game, so not zero, but far from central).
- solar: reduced for Povo Livre (0.20 -> 0.15), NOT zeroed. Two independent
  reasons solar matters less for a nomadic population: portable panels can
  be re-angled by hand, so the tilted-surface aspect/slope optimisation in
  solar_suitability_annual_120m matters less to them; and the layer is an
  ANNUAL average, which implicitly assumes occupying one site year-round --
  the wrong temporal frame for a population that moves seasonally (a known,
  documented limitation, not fixed here -- fixing it would need an actual
  migration-route/season model that doesn't exist yet). Solar is NOT zeroed
  because horizon shading and latitude still apply to any panel regardless
  of portability, and this is nominally a solarpunk setting for both
  populations, not just the Circulos.
- biome: dominant for Povo Livre (0.60) -- for a forest/foraging economy,
  vegetation type is close to the whole story. Moderate-high for Circulo
  (0.30) -- important, but shares the "economic base" signal with
  agriculture (0.25) rather than carrying it alone.

Povo Silencioso exclusion is applied to BOTH populations' composites, not
just Circulo's -- a default choice (the archipelago is framed as the Povo
Silencioso's own territory to respect, not specifically a Circulo-only
concern), not something Nico was asked to confirm. Revisit if that turns
out wrong.
"""

from __future__ import annotations

import numpy as np

__all__ = ["WEIGHTS_CIRCULO", "WEIGHTS_POVO_LIVRE", "weighted_composite"]

WEIGHTS_CIRCULO = {
    "slope": 0.20,
    "water": 0.05,
    "agriculture": 0.25,
    "solar": 0.20,
    "biome": 0.30,
}

WEIGHTS_POVO_LIVRE = {
    "slope": 0.10,
    "water": 0.05,
    "agriculture": 0.10,
    "solar": 0.15,
    "biome": 0.60,
}

for _name, _w in [("WEIGHTS_CIRCULO", WEIGHTS_CIRCULO), ("WEIGHTS_POVO_LIVRE", WEIGHTS_POVO_LIVRE)]:
    _total = sum(_w.values())
    assert abs(_total - 1.0) < 1e-9, f"{_name} must sum to 1.0, got {_total}"


def weighted_composite(
    layers: dict[str, np.ndarray],
    weights: dict[str, float],
    land_mask: np.ndarray,
    exclusion: np.ndarray | None = None,
) -> np.ndarray:
    """Weighted sum of `layers` (each already 0-1, NaN off-land is fine --
    treated as 0 contribution there since land_mask masks the output anyway)
    by `weights` (same keys, must match WEIGHTS_CIRCULO/WEIGHTS_POVO_LIVRE's
    keys and sum to 1.0), then multiplied by `exclusion` (0-1 multiplier,
    e.g. povo_silencioso_exclusion_120m) if given. Returns float64, NaN off
    land_mask, 0-1 on it (weights sum to 1 and both inputs are 0-1, so no
    explicit clipping is needed -- left unclipped deliberately so a values
    bug upstream would show up as an out-of-range number instead of being
    silently hidden).
    """
    total_w = sum(weights.values())
    assert abs(total_w - 1.0) < 1e-9, f"weights must sum to 1.0, got {total_w}"
    assert set(layers) == set(weights), f"layers/weights key mismatch: {set(layers)} vs {set(weights)}"

    acc = np.zeros_like(land_mask, dtype=np.float64)
    for key, w in weights.items():
        acc = acc + w * np.nan_to_num(layers[key], nan=0.0)

    if exclusion is not None:
        acc = acc * exclusion

    return np.where(land_mask, acc, np.nan)
