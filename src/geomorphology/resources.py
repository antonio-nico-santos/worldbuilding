"""
Tappa 8 -- Geomorphology: jade/pounamu + Vertice materials.

Per-lithology-class material attributes -- non-spatial (a lookup table,
not a raster field) except where noted. Full domain/verb/tier system for
Vertice mechanics lives in docs/reference/scenario_reference.md; this is
the geology-facing summary already locked in
docs/decisions/07_tappa7_regional_scenario.md S3, reproduced here as
structured data so it ships alongside the lithology/jade rasters instead
of only existing in prose.

`spatial` fields added this session (run_tappa8_resource_pods.py, decision doc S8e) --
four of the six materials now have an actual raster, generalizing jade's own
place_jade_pods() via the new place_material_pods(). Schist's gold/mica and volcanic's
magnetite (primary) deliberately do NOT get a new raster -- see each entry's `spatial`
note for why, and S8e for the full reasoning.
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
            "spatial": "resource_laumontite.npy -- 8 pods, UNIFORM placement (no citable "
            "within-greywacke gradient, unlike schist_grade's real metamorphic-rank story; "
            "see run_tappa8_resource_pods.py rather than inventing one)",
        },
    },
    "sedimentary_basin_fill": {
        "primary": {
            "material": "vivianite (authigenic iron-phosphate)",
            "domain": "Bios (tied to organic decay)",
            "citation": "forms in low-oxygen, organic-rich floodplain/bog sediment -- "
            "matches the Canterbury-Plains-analog depositional environment directly",
            "spatial": "resource_vivianite.npy -- 8 pods, weighted by a wetness proxy "
            "(inverse distance to stream), RESTRICTED to S9's wetland_backswamp sub-class "
            "(re-placed from S8e's original whole-basin_fill mask now that S9's Arable/"
            "Wetland-backswamp/Estuarine-coastal zonation actually exists as a raster -- "
            "see decision doc S8h).",
        },
        "mundane_only": [
            "bog iron (goethite/limonite, precipitated by iron-oxidizing processes in the "
            "SAME anoxic, organic-rich floodplain/bog setting vivianite occupies -- no "
            "electrical/magnetic/optical gating property, not a Vertice crystal. Co-located "
            "with vivianite's exact pod footprint (resource_bog_iron.npy), per "
            "scenario_reference.md S22.4's framing that a single wetland patch plausibly "
            "yields both, worked by different specialists. CITATION FLAG: well-documented as "
            "a GENERAL pre-industrial process (Iron Age Scandinavia through colonial North "
            "America), but NOT NZ-specific -- Te Ara's NZ iron history is exclusively the "
            "ironsand/titanomagnetite story. First exception in this catalogue to the "
            "NZ-specific-citation norm.",
        ],
        "secondary_weak": {
            "material": "reworked placer magnetite (titanomagnetite eroded from the "
            "volcanic zone, reconcentrated by rivers)",
            "domain": "Campo (weak tier -- lower-grade than the volcanic primary)",
            "citation": "real NZ ironsand-beach provenance mechanism",
            "citation_flagged_for_review": "'reconcentrated by rivers' doesn't hold up: "
            "mainland basin_fill and the volcanic zone are separate landmasses, 16.44 km of "
            "open water apart (Tappa 7 S1), no shared river catchment. Real NZ ironsand "
            "beaches actually form by COASTAL/longshore redistribution, which DOES work "
            "across a water gap -- the spatial placement below already uses that mechanism; "
            "this citation text itself hasn't been edited to match, pending Nico's sign-off.",
            "spatial": "resource_placer_magnetite.npy -- 8 pods, weighted by coastal "
            "proximity x proximity to the volcanic landmass (not by any river field), "
            "consistent with the coastal/longshore reframing above. RESTRICTED to S9's "
            "estuarine_coastal sub-class (re-placed from S8e's original whole-basin_fill "
            "mask, S12) -- a natural fit, since that sub-class is itself defined by "
            "coastal proximity. Structurally cannot overlap vivianite/bog_iron's pods "
            "(both restricted to wetland_backswamp, S10) -- S9's sub-zones are mutually "
            "exclusive by construction, confirmed by assertion when this ran.",
        },
    },
    "volcanic": {
        "primary": {
            "material": "magnetite (basalt titanomagnetite content)",
            "domain": "Campo (paleomagnetism)",
            "citation": "Banks Peninsula basalt; replaces an earlier, weaker olivine claim",
            "spatial": "no pod raster -- titanomagnetite is a disseminated bulk mineral in "
            "basalt (a few percent through the rock), not vein/joint-localized like the "
            "other five materials. Stays a class-level fact (the whole volcanic zone "
            "qualifies), which is MORE geologically honest than inventing discrete deposits.",
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
            "spatial": "resource_silver_copper.npy -- 6 pods (fewer than the target 8; the "
            "volcanic zone is small, 595.7 km2, and the 5 km minimum-separation constraint "
            "limits how many fit), weighted by vent_weight (the same geothermal-vent-proximity "
            "field lava tubes already use -- epithermal veins are a direct real hydrothermal/"
            "volcanic association, not a stretch to reuse it here). Each pod independently "
            "rolled silver vs copper, 25%/75% -- an UNREVIEWED placeholder ratio, no citation.",
        },
        "mundane_only": [
            "bauxite (gibbsite/boehmite lateritic weathering crust) -- no citable "
            "Vertice-domain-gating property, not a Vertice crystal. Real NZ citation: "
            "Northland (Otoroa/Matauri Bay) relict Pliocene-Pleistocene lateritic bauxite, "
            "a weathering product of the SAME basaltic rock this class represents -- small "
            "(largest real deposit ~20 Mt) and historically subeconomic/never mined in the "
            "real world. Nico's explicit call: treated here as an ACTIVE resource, with the "
            "real-world Hall-Heroult (1886) electrolysis technology gate closed narratively "
            "by Vertice-assisted electrolysis -- a process-level mechanic (not a property of "
            "the ore itself), recorded here as citation context only; further mechanical "
            "detail on HOW Vertices assist electrolysis is a Scenario-chat-level question, "
            "not resolved here. Spatial: resource_bauxite.npy -- deliberately small footprint "
            "(2 pods, 200-500 m radius vs the other materials' 8 pods/300-800 m), weighted by "
            "flatness (1/(1+slope_pct)) since real lateritic caps form/survive on low-relief "
            "ground and are stripped off steep slopes by erosion. Vertice-assisted extraction "
            "closes the TECHNOLOGY gate; it does not make the ore itself abundant -- these "
            "stay two separate constraints.",
        ],
    },
}
