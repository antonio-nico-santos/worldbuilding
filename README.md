# World-building GIS — Portfolio Case Study 2
 
Procedurally generated fictional terrain, climate, hydrology and biomes,
built as a from-scratch GIS pipeline (not a pre-built world generator) to
demonstrate technical GIS skills for a freelance portfolio and a TTRPG scenario.

This is the second time I've applied GIS methods to a worldbuilding pipeline. The first attempt (May 2025) had no version control, working notes scattered across Notion instead of alongside the code, and relied on real-world bathymetric data at a resolution/extent that produced oversized, undocumented files. The biome classification also used temperature alone, without precipitation as a second axis — a methodological gap this rebuild corrects. This time, every parameter and decision is tracked under version control from the start, and the pipeline is extended with hydrology and settlement-suitability analysis not present in the original.
 
Real-world reference used for climate model validation: South Island, New
Zealand (Southern Alps).
 
See `docs/decisions/` for the per-stage planning summaries and
`config/parameters.yml` for the current pipeline parameters.
 
This repository is intentionally kept separate from the portfolio site
repository (`gis-portfolio`) — only final lightweight exports (simplified
GeoJSON, pre-rendered images) are copied into the site's `public/data/`.
