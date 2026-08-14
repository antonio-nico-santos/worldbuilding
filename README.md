# World-building GIS — Portfolio Case Study 2
 
Procedurally generated fictional terrain, climate, hydrology, biomes, and
settlement suitability, built as a from-scratch GIS pipeline (not a
pre-built world generator) to demonstrate technical GIS skills for a
freelance portfolio and a TTRPG scenario.

This is the second time I've applied GIS methods to a worldbuilding pipeline. The first attempt (May 2025) had no version control, working notes scattered across Notion instead of alongside the code, and relied on real-world bathymetric data at a resolution/extent that produced oversized, undocumented files. Biome classification already used the Holdridge system (temperature and precipitation as the two axes) in that first attempt. This time, every parameter and decision is tracked under version control from the start, and the pipeline is extended with hydrology and settlement-suitability analysis not present in the original.
 
Real-world reference used for climate model validation: South Island, New
Zealand (Southern Alps).
 
See `docs/decisions/` for the per-stage planning summaries and
`config/parameters.yml` for the current pipeline parameters.

## Pipeline status

Tappe 0-6 are closed: procedural terrain, monthly climate, snow/ELA
metrics, hydrology, biomes, and solarpunk settlement suitability (a
weighted multi-criteria composite plus greedy, cost-distance-aware site
selection for 17 candidate settlements). This is considered complete,
sufficient content for the portfolio case study as it stands.

Two further stages are planned but currently paused, and are scoped as
depth for the TTRPG scenario side of this project rather than additional
portfolio material: **Tappa 7** (regional scenario deepening, still at
macro/domain-wide scale — dangerous-creature conflict zones, road/rail
infrastructure, dangerous seas) and **Tappa 8** (urban zoom, per-settlement
layout at the locations Tappa 6 already chose).
 
This repository is intentionally kept separate from the portfolio site
repository (`gis-portfolio`) — only final lightweight exports (simplified
GeoJSON, pre-rendered images) are copied into the site's `public/data/`.

## License

Code in this repository is licensed under MIT (see `LICENSE.md`). Creative
and narrative content — world, species, and scenario material under
`docs/decisions/`, `docs/reference/`, and elsewhere — is all rights
reserved (see `CONTENT-LICENSE.md`) and is not licensed for reuse.
