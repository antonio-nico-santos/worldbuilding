"""
Tappa 8 -- Geomorphology: jade/pounamu + Vertice materials.

Per-lithology-class material attributes -- non-spatial (a lookup table,
not a raster field) except where noted. Full domain/verb/tier system for
Vertice mechanics lives in docs/reference/scenario_reference.md; this is
the geology-facing summary already locked in
docs/decisions/07_tappa7_regional_scenario.md S3, reproduced here as
structured data so it ships alongside the lithology/jade rasters instead
of only existing in prose.

`spatial` fields added across this session (run_tappa8_resource_pods.py S8e; S10 iron/
aluminium; S12 placer-magnetite re-placement; S13 quartz/mica/gold, generalizing jade's
own place_jade_pods() via the shared place_material_pods() helper; S14 sedimentary_
limestone/granite). Thirteen materials now have an actual pod raster. Two deliberately do
NOT: volcanic's magnetite (primary, a disseminated bulk mineral, not vein-localized -- see
its own `spatial` note) and marble (S21, NO Vertice domain at all, an explicit null -- see
its own entry below, not a placeholder).
"""

VERTICE_MATERIALS = {
    "schist": {
        "richness": "richest of the four classes",
        "primary": {
            "material": "quartz, gold-bearing (orogenic, fault-zone veins)",
            "domain": "Onda (birefringence)",
            "citation": "Otago Schist Invincible Vein: quartz, gold, pyrite, "
            "arsenopyrite, muscovite, chlorite, calcite, albite",
            "spatial": "resource_quartz.npy -- 8 pods, weighted by schist_grade (same "
            "field jade uses), INDEPENDENT of jade's own pods (S13, generalizing "
            "place_jade_pods()'s own reasoning: a real vein system is heterogeneous at "
            "the pod scale, so quartz/mica/gold each get their own stochastic draw within "
            "the same high-grade eligible zone, not one shared footprint).",
        },
        "secondary": {
            "material": "muscovite/biotite mica, same veins",
            "domain": "Energia (piezoelectric)",
            "citation": "rarer than the quartz itself, same vein system",
            "spatial": "resource_mica.npy -- 8 pods, same eligible zone/weight field as "
            "quartz, independent seed (S13).",
        },
        "mundane_only": [
            "gold (no electrical/magnetic/optical gating property, not a Vertice "
            "crystal) -- resource_gold.npy, 8 pods, same eligible zone/weight field as "
            "quartz/mica, independent seed (S13). Checked directly: zero pairwise "
            "overlap among quartz/mica/gold this run -- incidental (all three draw from "
            "the same ~354.6 km2 zone), not by design.",
            "jade/pounamu (see jade_eligible_mask spatial subset -- same rock, not a "
            "Vertice crystal itself) -- resource_jade still uses jade_pods_v5.npy "
            "(S5, unchanged, its own locked raster).",
        ],
        "spatial_note": "quartz/mica/gold (S13) and jade (S5) all draw eligible ground "
        "from the SAME high-grade schist test (jade_eligible_mask + schist_grade, "
        "80th percentile), but each is an INDEPENDENT pod raster now, not one shared "
        "footprint -- this note previously said they all reuse jade_pods_v5.npy "
        "directly; that's no longer true for quartz/mica/gold. FLAG: quartz/mica/gold "
        "use lithology_v6's current schist extent; jade_pods_v5.npy still uses v5's "
        "(9013 cells / 8.11 km2 now marble in v6) -- a pre-existing staleness in jade's "
        "own raster, not reconciled here, out of this fix's scope.",
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
    "sedimentary_limestone": {
        "secondary": {
            "material": "optical-grade calcite (\"Iceland spar\" -- clear calcite found "
            "in cavities within the limestone, not the bulk rock itself)",
            "domain": "Onda (birefringence) -- rare superior-grade SECONDARY crystal "
            "alongside quartz. Does NOT replace quartz/schist as Onda's primary material; "
            "quartz stays primary, calcite is a rarer upgrade.",
            "citation": "calcite's birefringence (delta-n ~= 0.17) is roughly 20x "
            "stronger than quartz's (delta-n ~= 0.009) -- the real physical basis of the "
            "Iceland-spar/Viking-sunstone effect. Clear, intact calcite crystals large "
            "enough to use are much rarer to find in limestone cavities than vein quartz "
            "is in schist, which is why this is modeled as a prized upgrade rather than a "
            "bulk replacement -- that rarity is carried by the primary/secondary TIER, not "
            "by a reduced pod count/footprint (S14 deliberately used the SAME n_pods/"
            "radius as every other generalized material; see run_tappa8_limestone_"
            "granite_materials.py for why).",
            "citation_source": "scenario_reference.md S21 (Scenario chat, relayed "
            "verbatim by Nico, 2026-08-20).",
            "spatial": "resource_calcite.npy -- 7/8 target pods placed (S14), UNIFORM "
            "weight (no citable within-class spatial gradient -- no analog to "
            "schist_grade exists for limestone), eligible ground = whole "
            "sedimentary_limestone class (87.72 km2, spanning the authored \"North Coast "
            "Limestone\" + \"Sedimentary Bay\" zones -- the raster does not distinguish "
            "which authored zone a cell came from, same convention as every other class).",
        },
    },
    "granite": {
        "secondary": {
            "material": "muscovite/biotite mica -- SAME mineral as schist's Energia "
            "source, DIFFERENT host rock. Granite is a felsic igneous rock naturally "
            "composed of quartz + feldspar + mica.",
            "domain": "Energia (piezoelectric) -- same property, same mechanic as "
            "schist's mica; this is a second, independent PLACE to find it, not a new "
            "gating rule.",
            "citation": "no new physical-property claim needed -- literally the same "
            "mineral as the schist assignment, diversified sourcing so a Circulo doesn't "
            "have to sit on schist specifically to reach Energia.",
            "citation_source": "scenario_reference.md S21 (Scenario chat, relayed "
            "verbatim by Nico, 2026-08-20).",
            "spatial": "resource_mica_granite.npy -- ONLY 3/8 target pods placed (S14), "
            "UNIFORM weight (no citable within-class gradient, granite is compositionally "
            "uniform in this model). NAMED SEPARATELY from resource_mica.npy (schist's "
            "mica raster, S13) -- same mineral, different host rock, different footprint, "
            "not a duplicate. FLAG: granite's eligible ground (13.84 km2, \"Granite "
            "South\" alone) is far smaller than any other material's, so the project's "
            "standard 5 km min_separation_km left room for only 3 of the usual 8 pods -- "
            "checked directly, not forced to 8. Worth deciding whether granite warrants a "
            "smaller separation constant if 8 comparable pods matters for play.",
        },
        "tertiary": {
            "material": "quartz -- SAME mineral as schist's Onda primary source, "
            "DIFFERENT host rock.",
            "domain": "Onda (birefringence) -- minor/tertiary source, behind schist "
            "(primary) and limestone-calcite (rare secondary).",
            "citation": "same reasoning as granite's mica entry -- no new physical "
            "property claim, diversified sourcing only.",
            "citation_source": "scenario_reference.md S21 (Scenario chat, relayed "
            "verbatim by Nico, 2026-08-20).",
            "spatial": "resource_quartz_granite.npy -- ONLY 2/8 target pods placed "
            "(S14), UNIFORM weight, same 13.84 km2 eligible ground as mica_granite "
            "(independent seed -- checked directly: 0.0 km2 overlap between the two "
            "granite materials this run). NAMED SEPARATELY from resource_quartz.npy "
            "(schist's quartz raster, S13) for the same reason as mica_granite above.",
        },
    },
    "marble": {
        "no_vertice_domain": True,
        "material": "marble itself (metamorphosed/recrystallized limestone) -- stays "
        "purely mundane: construction and monumental stone.",
        "citation": "marble's recrystallization destroys the very thing that makes "
        "limestone's calcite useful to Onda -- the parallel-lattice clarity birefringence "
        "depends on. Marble is optically cloudy/light-scattering, not transparent. "
        "Separately, calcite's crystal structure has a center of symmetry, so it was "
        "never piezoelectric either -- it wouldn't have qualified for Energia by that "
        "route either.",
        "citation_source": "scenario_reference.md S21 (Scenario chat, relayed verbatim "
        "by Nico, 2026-08-20). EXPLICIT null, not a placeholder or an oversight -- do "
        "NOT add a spatial pod raster for marble without a new decision superseding this "
        "one.",
        "spatial": "none -- no pod raster, no bit in resource_blend.npy (S14). marble's "
        "class-level facts (77.8 km2, three authored zones -- Marble North/Forest "
        "Marble/Marble Wall) live in lithology_v6/resources for excavation-effort (S8g) "
        "and transport-friction (S8f) purposes only, unrelated to this entry.",
    },
}
