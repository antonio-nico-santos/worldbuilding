"""
Tappa 8 -- Geomorphology: jade/pounamu + Vertice materials.

Per-lithology-class material attributes -- non-spatial (a lookup table,
not a raster field) except where noted. Full domain/verb/tier system for
Vertice mechanics lives in docs/reference/scenario_reference.md; this is
the geology-facing summary already locked in
docs/decisions/07_tappa7_regional_scenario.md S3, reproduced here as
structured data so it ships alongside the lithology/jade rasters instead
of only existing in prose.
"""

VERTICE_MATERIALS = {
    "schist": {
        "richness": "richest of the four classes",
        "primary": {
            "material": "gold-bearing quartz veins (orogenic, fault-zone)",
            "domain": "Onda (birefringence)",
            "citation": "Otago Schist Invincible Vein: quartz, gold, pyrite, "
            "arsenopyrite, muscovite, chlorite, calcite, albite",
        },
        "secondary": {
            "material": "muscovite/biotite mica, same veins",
            "domain": "Energia (piezoelectric)",
            "citation": "rarer than the quartz itself, same vein system",
        },
        "mundane_only": [
            "gold (no electrical/magnetic/optical gating property, not a Vertice crystal)",
            "jade/pounamu (see jade_eligible_mask spatial subset -- same rock, not a "
            "Vertice crystal itself)",
        ],
        "spatial_note": "co-locates with the jade/pounamu high-grade band, see lithology.jade_eligible_mask",
    },
    "greywacke_argillite": {
        "primary": {
            "material": "laumontite (a zeolite, fills veinlets/joints/shatter zones)",
            "domain": "Materia (molecular-sieve structure -- solid/liquid/gas)",
            "citation": "genuinely specific to NZ's Mesozoic greywacke ranges",
        },
    },
    "sedimentary_basin_fill": {
        "primary": {
            "material": "vivianite (authigenic iron-phosphate)",
            "domain": "Bios (tied to organic decay)",
            "citation": "forms in low-oxygen, organic-rich floodplain/bog sediment -- "
            "matches the Canterbury-Plains-analog depositional environment directly",
        },
        "secondary_weak": {
            "material": "reworked placer magnetite (titanomagnetite eroded from the "
            "volcanic zone, reconcentrated by rivers)",
            "domain": "Campo (weak tier -- lower-grade than the volcanic primary)",
            "citation": "real NZ ironsand-beach provenance mechanism",
        },
    },
    "volcanic": {
        "primary": {
            "material": "magnetite (basalt titanomagnetite content)",
            "domain": "Campo (paleomagnetism)",
            "citation": "Banks Peninsula basalt; replaces an earlier, weaker olivine claim",
        },
        "secondary": {
            "material": "native silver (rare, epithermal veins) / native copper "
            "(common, volcanic amygdules)",
            "domain": "Mente (a 'tuning' pair -- silver higher-quality, copper the common half)",
            "citation": "silver: real NZ analog Hauraki Goldfield; copper: real Keweenaw "
            "Peninsula analog",
            "refinement": "the mineral alone gives a weak connection; paired with a black "
            "fungus found in Wet Forest for the real Mente crystal -- combination method "
            "still in-world 'under study', discovered ~3-4 years before the TTRPG present",
        },
    },
}
