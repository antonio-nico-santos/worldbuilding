"""
Tappa 7 -- fauna suitability composites (high-threat predators, alpine
outposts, domesticated-animal point/zone lookups).

Reuses layers already produced by Tappa 5/6 (biome_id, biotemperature_c,
slope_suitability_120m, slope_pct_120m, solar_suitability_annual_120m,
suitability_povo_livre_120m) and Tappa 4 (stream_mask.npy, 30m -- resampled
to this grid) plus the real 17-Circulo site table and terrain_ridges.geojson.
Method locked in docs/decisions/07_tappa7_regional_scenario.md; this script
is the first real computation of it -- run once this session, after several
search-zone/weight corrections documented inline (see NACRE/TWINSHADOWS/
OUTPOST sections below).

v6 (Nico's redesign): Nacre's whole architecture was rebuilt this pass --
caves -> coves (dens) -> a wild-yak prey layer -> one graduated suitability
field, no hard exclusao. This RETIRES nacre_exposure_120m.npy and
nacre_massif_labels_120m.npy from earlier in the session -- delete those
two files (and their .tif conversions) locally, they no longer reflect the
current model and this script no longer regenerates them.

Outputs -> data/processed/fauna/:
  nacre_suitability_120m.npy (v6 -- graduated, land-masked only, no exclusao)
  nacre_threat_band_120m.npy (0=safe <0.15, 1=graduated, 2=danger >=0.9)
  nacre_pack_territory_120m.npy (4-zone, nearest-cove assignment)
  nacre_coves.geojson (4 den points: 3 main-spine, 1 South Branch)
  yak_suitability_120m.npy (wild/feral yak distribution, biome x slope)
  twinshadows_suitability_120m.npy (unchanged this pass)
  outpost_composite_120m.npy
  outpost_search_mask_120m.npy
  outpost_candidates.geojson (18 sites, status_prior now sampled directly from
    nacre_suitability_120m.npy against Nico's 0.9/0.15 thresholds)
  povo_livre_zone_120m.npy (suitability_povo_livre_120m > 0.80)
  nacre_grassmother_incursion_120m.npy (v8 addition -- seasonal, slope-weighted
    cost-distance from each cove into the Grassland biome, see v8 note below)
  tappa7_fauna_composite_meta.json

v8 addition (Nico's call): a seasonal Nacre incursion layer specific to hunting
Grassmothers in Lowland Steppe/Grassland, modeled as central-place foraging from
each cove -- real citation: wolves during pup-rearing/denning season are
range-restricted around a fixed den, typically ~10-20km foraging radius, with
rare longer forays; used here as the decay anchor. Kept as a SEPARATE raster,
not blended into nacre_suitability_120m.npy -- the 0.9/0.15 threat bands and
outpost status_prior are already locked and this doesn't touch them. Distance
is a real slope-weighted cost-distance (skimage.graph.MCP_Geometric), not a
straight Euclidean decay like the rest of this pipeline -- Nico asked
specifically for slope to shape how far the incursion penetrates, which a
plain distance-to-cove reuse would not have captured.

v10/v11 note: a computed Skydrifter flight-corridor raster (skydrifter_route_
area_120m.npy) was built in v10, then RETIRED in v11 at Nico's call -- a species
that overflies the whole crossing "mostly instantaneously" and never stops has
no behavioral reason to detour around terrain, so the computed route's
terrain-hugging curve overstated what the lore actually supports. No GIS layer
for this species anymore; see scenario_reference.md for the plain descriptive
lore instead (wide, roughly straight band over the North plains). DELETE the
stale skydrifter_route_area_120m.npy/.tif files locally -- this script no
longer regenerates them.
"""
import numpy as np, json, math
from scipy import ndimage
from skimage.graph import MCP_Geometric

BASE = 'data/processed'
INPUT = 'data/input'
OUT = 'data/processed/fauna'

meta = json.load(open(f'{BASE}/biomes/tappa5_biomes_meta.json'))
g = meta['grid']
ny, nx = g['ny'], g['nx']
xmin, ymax = g['xmin'], g['ymax']
res_x, res_y = 119.92619926199262, 119.9400299850075

BIOME_NAMES = {
    0: 'Ocean', 1: 'Permanent Snow & Ice', 2: 'Alpine Fellfield', 3: 'Alpine Tundra',
    4: 'Subalpine Wet Forest', 5: 'Subalpine Woodland', 6: 'Subalpine Dry Scrub',
    7: 'Temperate Forest', 8: 'Woodland / Shrubland', 9: 'Lowland Steppe / Grassland',
}


def xy_to_rc(x, y):
    return (ymax - y) / res_y, (x - xmin) / res_x


def rc_to_xy(r, c):
    return xmin + c * res_x, ymax - r * res_y


def lut_apply(biome_id, lut):
    out = np.zeros(biome_id.shape, dtype=np.float32)
    for k, v in lut.items():
        out[biome_id == k] = v
    return out


def pick_maxima(mask, composite, n_target, min_sep_km, res_x, res_y, taken=None):
    """Greedy top-value picker with a minimum-separation constraint. `taken` is an
    optional list of (r, c) already-chosen points from OTHER calls (e.g. the other
    massif's candidates) that also count against the separation constraint -- fixes
    the v3 bug where main-spine and South-Branch picks were only separation-checked
    against their own zone, so both zones could independently converge on the same
    ridge-overlap hot spot."""
    if taken is None:
        taken = []
    ys, xs = np.where(mask)
    vals = composite[ys, xs]
    order = np.argsort(vals)[::-1]
    chosen = []
    for idx in order:
        r, c = ys[idx], xs[idx]
        if all(math.hypot((r - cr) * res_y, (c - cc) * res_x) / 1000.0 >= min_sep_km
               for cr, cc in list(taken) + [(cr2, cc2) for cr2, cc2, _ in chosen]):
            chosen.append((r, c, vals[idx]))
        if len(chosen) >= n_target:
            break
    return chosen


def rasterize_line_with_arclength(coords, ny, nx, xy_to_rc):
    """Rasterize a ridge polyline, and for every grid cell return the arc-length (km)
    of its NEAREST ridge sample point -- used to bin candidates evenly along the
    ridge's length rather than letting them cluster wherever the composite peaks."""
    pts_rc = []
    arclens = []
    cum = 0.0
    for i in range(len(coords) - 1):
        x0, y0 = coords[i]
        x1, y1 = coords[i + 1]
        seg_len_m = math.hypot(x1 - x0, y1 - y0)
        n_steps = max(1, int(seg_len_m / 60))
        for j, t in enumerate(np.linspace(0, 1, n_steps + 1)):
            x = x0 + (x1 - x0) * t
            y = y0 + (y1 - y0) * t
            r, c = xy_to_rc(x, y)
            pts_rc.append((int(round(r)), int(round(c))))
            arclens.append((cum + seg_len_m * t) / 1000.0)
        cum += seg_len_m

    mask = np.zeros((ny, nx), dtype=bool)
    arclen_at_seed = {}
    for (r, c), a in zip(pts_rc, arclens):
        if 0 <= r < ny and 0 <= c < nx:
            mask[r, c] = True
            arclen_at_seed[(r, c)] = a  # last write wins, fine for a dense polyline sample

    seed_rs = np.array([rc[0] for rc in arclen_at_seed])
    seed_cs = np.array([rc[1] for rc in arclen_at_seed])
    seed_vals = np.array([arclen_at_seed[rc] for rc in arclen_at_seed])

    _, (ind_r, ind_c) = ndimage.distance_transform_edt(~mask, return_indices=True)
    # nearest-seed lookup: for every pixel, find which seed index it was assigned to
    seed_lookup = np.full((ny, nx), -1, dtype=np.int64)
    seed_lookup[seed_rs, seed_cs] = np.arange(len(seed_vals))
    nearest_seed_idx = seed_lookup[ind_r, ind_c]
    arclen_raster = seed_vals[nearest_seed_idx]
    total_km = arclens[-1]
    return arclen_raster, total_km


def pick_along_ridge(zone_mask, composite, arclen_raster, total_km, n_target, min_sep_km,
                      res_x, res_y, taken=None):
    """Bin the zone into n_target equal-arclength segments along the ridge and pick the
    single best composite cell per bin -- guarantees spread along the ridge's full
    length instead of relying on min-separation alone to break up a cluster."""
    if taken is None:
        taken = []
    chosen = []
    bin_edges = np.linspace(0, total_km, n_target + 1)
    for i in range(n_target):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        bin_mask = zone_mask & (arclen_raster >= lo) & (arclen_raster < hi)
        cand = pick_maxima(bin_mask, composite, 1, min_sep_km, res_x, res_y,
                            taken=list(taken) + [(r, c) for r, c, _ in chosen])
        if cand:
            chosen.append(cand[0])
    return chosen


def main():
    biome_id = np.load(f'{BASE}/biomes/biome_id.npy')
    land = biome_id != 0
    biotemp = np.load(f'{BASE}/biomes/biotemperature_c.npy')
    slope_suit = np.load(f'{BASE}/suitability/slope_suitability_120m.npy')
    slope_pct = np.load(f'{BASE}/suitability/slope_pct_120m.npy')
    solar_suit = np.load(f'{BASE}/suitability/solar_suitability_annual_120m.npy')
    povo_livre = np.load(f'{BASE}/suitability/suitability_povo_livre_120m.npy')
    circulo_sites = json.load(open(f'{BASE}/suitability/circulo_candidate_sites.geojson'))
    ridges = json.load(open(f'{INPUT}/terrain_ridges.geojson'))
    ridge_by_name = {f['properties']['name']: f for f in ridges['features']}

    # Ruggedness: raw slope_pct with an inverted curve (NOT slope_suitability_120m,
    # which is tuned for human buildability, opposite sense -- corrected against
    # the real Tappa-6 layer inventory earlier this session).
    p95 = np.nanpercentile(slope_pct[land], 95)
    ruggedness = np.nan_to_num(np.clip(slope_pct / p95, 0, 1))

    # Circulo point lookup (kept for the domesticated-animal section below -- Nacre's
    # nucleo no longer uses distance-to-Circulo, see v6 NACRE block).
    circulo_mask = np.zeros((ny, nx), dtype=bool)
    circulo_rc = {}
    for f in circulo_sites['features']:
        x, y = f['geometry']['coordinates']
        r, c = (int(round(v)) for v in xy_to_rc(x, y))
        circulo_mask[r, c] = True
        circulo_rc[f['properties']['name']] = (r, c, f['properties']['biome_at_site'])

    # Prey-base proximity (reusing regular-bucket percent-table values, locked earlier
    # this session in tappa7_fauna_biome_percent.xlsx). Twinshadows only, as of v6 --
    # Nacre's own prey-base term is now the new wild-yak layer (see below).
    LUT_BLACKNOSE_SNAKETAIL = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0.15, 5: 0.0, 6: 0.075, 7: 1.0, 8: 0.15, 9: 0.5}
    preybase_twin = lut_apply(biome_id, LUT_BLACKNOSE_SNAKETAIL)

    DECAY_KM, SOUTH_BRANCH_WEIGHT = 20.0, 1.3  # shared decay constant, reused throughout this
    # script wherever a "distance from a territorial anchor" term is needed (coves, alpine
    # affinity, outpost ridge-exposure) -- one real mechanism (predator/resource influence
    # fading with distance), not independently tuned per use.

    # --- Alpine massifs (real terrain fact, independent of any species model) ---
    # ids 2/3 = Alpine Fellfield + Alpine Tundra. Used for: cave-candidate search zone,
    # cove siting, the outpost search zone (v5, already locked), and the main-spine/
    # South-Branch geographic split for outposts. NOT used as a hard suitability
    # exclusao for Nacre anymore -- see v6 below.
    alpine_mask = np.isin(biome_id, [2, 3])
    labeled, n = ndimage.label(alpine_mask, structure=np.ones((3, 3), dtype=int))
    sizes = ndimage.sum(alpine_mask, labeled, index=range(1, n + 1))
    order = np.argsort(sizes)[::-1]
    main_spine_label, south_branch_label = order[0] + 1, order[1] + 1
    main_spine_mask = labeled == main_spine_label
    south_branch_mask = labeled == south_branch_label
    dist_ms_km = ndimage.distance_transform_edt(~main_spine_mask) * ((res_x + res_y) / 2) / 1000.0
    dist_sb_km = ndimage.distance_transform_edt(~south_branch_mask) * ((res_x + res_y) / 2) / 1000.0
    zone_south = dist_sb_km < dist_ms_km

    # --- NACRE v6 (Nico's redesign, replacing the exclusao/nucleo Tier-1/Tier-2
    # architecture entirely): caves -> coves (dens) -> wild yak prey base -> a single
    # graduated suitability field, no hard biome mask. ---
    #
    # 1. Cave candidates: reuses the "talus/pseudokarst cave" mechanism already
    # documented (decision-only) in 07_tappa7_regional_scenario.md sec.2 -- steep
    # relief intersected with a nearby stream, rock-type agnostic. stream_mask.npy is
    # a Tappa 4 layer at 30m resolution (this grid is 120m) -- resampled by nearest-
    # neighbor coordinate lookup, not a block-reduce, since the two grids' cell counts
    # don't divide evenly (5334/4334 vs 1334/1084 -- off by a few cells, a small
    # padding difference between the two pipelines, not a bug in either).
    stream_mask_30m = np.load(f'{BASE}/hydrology/stream_mask.npy')
    dist_stream_km_30m = ndimage.distance_transform_edt(~stream_mask_30m) * 30.0 / 1000.0
    ny_s, nx_s = stream_mask_30m.shape
    rows_120, cols_120 = np.arange(ny), np.arange(nx)
    y_centers = ymax - (rows_120 + 0.5) * res_y
    x_centers = xmin + (cols_120 + 0.5) * res_x
    r_s = np.clip(((ymax - y_centers) / 30.0).astype(int), 0, ny_s - 1)
    c_s = np.clip(((x_centers - xmin) / 30.0).astype(int), 0, nx_s - 1)
    dist_stream_km = dist_stream_km_30m[np.ix_(r_s, c_s)]
    # Honest finding, consistent with every other stream-density check already in this
    # doc (Mudlizard, water_suitability_120m): 99.9% of the alpine band is already
    # within 1km of a mapped stream, so "near stream" barely discriminates candidates
    # here -- it's kept as a real criterion (matches the documented cave mechanism)
    # but in practice the cave score below is driven almost entirely by slope.
    CAVE_SLOPE_PCTILE = 75  # steep relative to the alpine band itself, not the whole map
    cave_slope_threshold = np.nanpercentile(slope_pct[alpine_mask], CAVE_SLOPE_PCTILE)
    CAVE_STREAM_KM = 1.5
    cave_candidate_mask = alpine_mask & (slope_pct >= cave_slope_threshold) & (dist_stream_km <= CAVE_STREAM_KM)

    # v8 (Nico's call, lore consistency catch): checked the v6 cove placements directly --
    # 3 of the 4 sat just 0.12km from the Subalpine boundary (real max possible depth into
    # the alpine band, via distance-to-Subalpine, is 8.09km; median across the whole band is
    # 1.98km), because "steep + near-stream" alone doesn't push siting away from the biome
    # edge -- the same clustering-at-the-boundary dynamic already flagged for outposts (v5)
    # showed up again here for the same underlying reason (transition terrain is often
    # locally steep). A den that close to Subalpine Wet Forest/Woodland/Dry Scrub would
    # plausibly hunt Blacknose/Snaketail/Furypack territory as routine behavior, undercutting
    # the locked lore ("baseline diet is yaks... rare seasonal descents into the lowlands are
    # the actual danger event, not routine") -- Nico's read: coves should sit deeper in the
    # true alpine interior, up to and including the Permanent Snow & Ice border (unlike
    # Subalpine proximity, PSI proximity isn't a lore problem -- it's the harshest, most
    # "confined by anthropogenic range contraction" end of the range, consistent with the
    # existing citation, not an inconsistency).
    #
    # Fix: added an interior-depth term (distance to the nearest SUBALPINE cell specifically,
    # not to any non-alpine cell -- proximity to Permanent Snow & Ice is deliberately NOT
    # penalized) to the cove-picking score, alongside ruggedness. cave_candidate_mask (steep +
    # near-stream) still gates real cave eligibility; interior depth now shapes WHICH of those
    # candidates wins.
    dist_to_subalpine_km = ndimage.distance_transform_edt(~np.isin(biome_id, [4, 5, 6])) * \
        ((res_x + res_y) / 2) / 1000.0
    COVE_INTERIOR_NORM_KM = 5.0  # ~p90 of dist_to_subalpine_km within the alpine band --
    # candidates at or beyond this depth get full interior credit, not a literal cap on how
    # deep a cove can be.
    interior_depth_score = np.clip(dist_to_subalpine_km / COVE_INTERIOR_NORM_KM, 0, 1)
    COVE_W_RUGGED, COVE_W_INTERIOR = 0.40, 0.60
    cove_score = (COVE_W_RUGGED * ruggedness + COVE_W_INTERIOR * interior_depth_score).astype(np.float32)

    # 2. Coves (dens): 3 on the main spine, 1 on South Branch -- matches the already-
    # locked pack-to-massif split (South Branch = 1 pack, "the South family"; main
    # spine = the other 3), with a minimum separation so the 3 main-spine dens read as
    # distinct territories, not one cluster. COVE_MIN_SEP_KM reuses FALLOFF_KM's value
    # (15km, already locked for schist/jade/outpost siting) rather than inventing a new
    # distance parameter.
    COVE_MIN_SEP_KM = 15.0
    main_spine_coves = pick_maxima(cave_candidate_mask & main_spine_mask, cove_score, 3,
                                    COVE_MIN_SEP_KM, res_x, res_y)
    south_branch_coves = pick_maxima(cave_candidate_mask & south_branch_mask, cove_score, 1,
                                      COVE_MIN_SEP_KM, res_x, res_y)
    all_coves = [(r, c, v, 'main_spine', i) for i, (r, c, v) in enumerate(main_spine_coves)] + \
                [(r, c, v, 'south_branch', i) for i, (r, c, v) in enumerate(south_branch_coves)]
    cove_mask = np.zeros((ny, nx), dtype=bool)
    for r, c, v, massif, i in all_coves:
        cove_mask[r, c] = True

    # v8 (Nico's follow-up catch): dist_to_cove_km was pure Euclidean
    # (ndimage.distance_transform_edt on the point mask) -- a straight-line ruler
    # distance that can't tell a short walk across flat ground from an equally
    # short line drawn straight through a cliff or ravine. Nico's own example:
    # an outpost that reads "closest" by straight-line distance could actually be
    # much harder for a pack to reach if steep terrain sits between it and the
    # cove -- which a Euclidean distance structurally cannot express. Replaced
    # with the same real slope-weighted cost-distance already built for the
    # Grassmother incursion layer below (skimage.graph.MCP_Geometric): friction
    # per cell = 1.0 + SLOPE_FRICTION_WEIGHT * ruggedness, accumulated along the
    # actual path, not a straight line. Computed once here and reused for the
    # Grassmother layer too (same cost surface, same source points -- no reason
    # to run the cost-distance transform twice). This changes dist_to_cove_suit,
    # which is 35% of nacre_suitability's weight and therefore the exact signal
    # outpost status_prior samples -- a real recalibration, not cosmetic; see the
    # before/after comparison in the delivery notes.
    SLOPE_FRICTION_WEIGHT = 2.0  # modeling choice (not a citation), reused from
    # the ruggedness field already computed above -- a cell at ~p95 slope or
    # steeper costs SLOPE_FRICTION_WEIGHT+1x a flat cell to cross.
    slope_friction = (1.0 + SLOPE_FRICTION_WEIGHT * ruggedness).astype(np.float64)
    cove_mcp = MCP_Geometric(slope_friction, fully_connected=True)
    cove_starts = [(r, c) for r, c, v, massif, i in all_coves]
    cove_cumulative_cost, _ = cove_mcp.find_costs(cove_starts)
    # cumulative cost is in "cost units" (1 unit = 1 flat cell-width at
    # friction=1); converted to an effective-km figure via the grid's real cell
    # size, same convention as every other distance term in this script, so the
    # existing DECAY_KM (already locked, unrelated to this change) still means
    # a real distance.
    dist_to_cove_km = cove_cumulative_cost * ((res_x + res_y) / 2) / 1000.0
    dist_to_cove_suit = np.exp(-dist_to_cove_km / DECAY_KM).astype(np.float32)

    # v8c (Nico's follow-up, extending the dist_to_cove fix for consistency):
    # pack-territory assignment was still Euclidean-nearest-cove -- "whose
    # family's turf is this cell" decided by ruler distance, the same blind
    # spot as the old dist_to_cove_suit had. A cell that's Euclidean-closest to
    # one cove could be genuinely easier to reach FROM a different cove if a
    # ravine/cliff sits on the straight-line path to the "nearest" one -- real
    # parallel in carnivore ecology, territory/home-range delineation is
    # standardly done with cost-distance or resistance surfaces, not straight-
    # line Voronoi, for exactly this reason. The multi-source cost-distance
    # computed above (cove_cumulative_cost) only gives the cost to the NEAREST
    # cove, not which one -- doesn't carry a label -- so per-cove assignment
    # needs one cost-distance run per cove (4 calls, same friction surface,
    # each with a single source point) and a per-cell argmin across the 4
    # resulting cost rasters, not the multi-source shortcut used above.
    cove_cost_per_source = np.stack([
        MCP_Geometric(slope_friction, fully_connected=True).find_costs([(r, c)])[0]
        for r, c, v, massif, i in all_coves
    ], axis=0)
    pack_territory = (np.argmin(cove_cost_per_source, axis=0) + 1).astype(np.int16)

    # --- v8 addition: Grassmother incursion layer (Nico's call) ---
    # Grassmothers are locked at 100% Lowland Steppe/Grassland (id 9), with only
    # marginal (15%) presence in Subalpine Dry Scrub/Woodland-Shrubland -- those
    # two are edge/transition habitat for the species, not its actual hunting
    # ground, so this incursion layer targets Grassland only, per Nico's choice.
    # Reuses the exact same cost-distance computed above for dist_to_cove_km
    # (same friction surface, same 4 cove source points) -- no need to recompute.
    gm_effective_km = dist_to_cove_km
    # Decay anchor: real central-place-foraging citation, not an arbitrary
    # number or a reused pipeline constant (DECAY_KM=20 and FALLOFF_KM=15 were
    # both derived for unrelated things -- ridge falloff/schist radius -- not
    # predator foraging range). Wolves during denning/pup-rearing season are
    # commonly reported restricted to routine foraging trips within ~10-20km of
    # the den, with rarer, longer forays beyond that -- 15km (mid-range) used as
    # the exponential decay constant, so effective_km=15 reads ~37% suitability,
    # matching a "still plausible but declining" reading, and by ~45km (3x decay)
    # it's under 5% -- rare event range, not routine, consistent with the locked
    # "rare seasonal descents" lore.
    GM_DECAY_KM = 15.0
    grassland_mask = biome_id == 9
    nacre_grassmother_incursion = np.where(
        grassland_mask, np.exp(-gm_effective_km / GM_DECAY_KM), 0.0
    ).astype(np.float32)

    # 3. Wild yak distribution -- feral populations, NOT tied to Circulo/outposts (a
    # deliberate departure from the "yak goes feral at an abandoned outpost" point-
    # source narrative already locked elsewhere in the doc; that narrative hook still
    # stands as extra local color, this is the ambient baseline population). "Biome x
    # slope", literally: a habitat mask (Alpine Fellfield/Tundra + the Subalpine band,
    # ids 2-6 -- real yak habitat spans alpine meadow and subalpine shrub margins, the
    # Tibetan-plateau analog) multiplied by a slope-favorability curve MUCH more
    # permissive than the human slope_suitability_120m curve: real yaks are sure-
    # footed high-altitude grazers, tolerant of much steeper ground than people build
    # on. Full suitability up to 15% grade, exponential decay with a 60%-grade
    # characteristic scale beyond that (never hits exactly 0 -- consistent with the
    # graduated, non-hard-cutoff philosophy Nico set for Nacre below).
    yak_biome_mask = np.isin(biome_id, [2, 3, 4, 5, 6])
    yak_slope_suit = np.exp(-np.maximum(0, slope_pct - 15.0) / 60.0).astype(np.float32)
    wild_yak_suit = np.where(yak_biome_mask, yak_slope_suit, 0.0).astype(np.float32)

    # 4. Alpine-biome affinity -- soft, not the old hard exclusao: 1.0 inside Alpine
    # Fellfield/Tundra, decaying with distance outside it (same DECAY_KM mechanism as
    # dist_to_cove_suit). This is what gives Nacre's suitability a graduated presence
    # in the lowlands too, folding the old separate Tier-2 "seasonal incursion" raster
    # into this single field instead of maintaining it as a second layer.
    dist_to_alpine_km = ndimage.distance_transform_edt(~alpine_mask) * ((res_x + res_y) / 2) / 1000.0
    alpine_affinity = np.exp(-dist_to_alpine_km / DECAY_KM).astype(np.float32)

    # 5. Final suitability -- graduated, land-masked only (no exclusao). Weights:
    # den proximity 35% (highest -- real territorial predators concentrate activity
    # around natal dens), ruggedness/remoteness 20%, wild-yak prey availability 25%,
    # alpine-biome affinity 20%. Proposed defaults, not independently re-confirmed with
    # Nico line by line -- flagged as adjustable.
    NACRE_W_COVE, NACRE_W_RUGGED, NACRE_W_PREY, NACRE_W_ALPINE = 0.35, 0.20, 0.25, 0.20
    nacre_nucleo_raw = (NACRE_W_COVE * dist_to_cove_suit + NACRE_W_RUGGED * ruggedness +
                          NACRE_W_PREY * wild_yak_suit + NACRE_W_ALPINE * alpine_affinity)
    # Calibration step, needed for Nico's 0.9/0.15 thresholds to mean anything: a convex
    # combination of 4 independently-peaking terms essentially never hits 1.0 (all 4 would
    # have to peak at the exact same cell) -- checked directly, the raw composite's land-wide
    # max is 0.82, so >=0.9 would be permanently unreachable and every candidate would read
    # 'temporary refuge' by default, collapsing status differentiation right back to the
    # problem this whole redesign was meant to fix. Rescaled by the land-wide 99th percentile
    # (same normalization pattern already used elsewhere in this pipeline for distance-based
    # terms, e.g. the old dist_to_circulo_suit) so the genuinely most favorable ~1% of land
    # approaches 1.0 -- preserves relative ordering, makes the absolute thresholds real.
    nacre_p99 = np.nanpercentile(nacre_nucleo_raw[land], 99)
    nacre_nucleo = np.clip(nacre_nucleo_raw / nacre_p99, 0, 1)
    nacre_suitability = np.where(land, nacre_nucleo, 0.0).astype(np.float32)

    # Graduated threat bands, per Nico's spec -- for reading/visualization only, the
    # continuous nacre_suitability raster stays the primary output:
    # >=0.9 = high danger (an outpost or Circulo there reads as untenable/"abandoned"-
    # tier risk); <0.15 = comparatively safe; in between = graduated threat, the
    # "constant threat, dangerous and safer areas" Nico asked for, in place of a binary
    # in/out territory.
    NACRE_DANGER_THRESHOLD, NACRE_SAFE_THRESHOLD = 0.9, 0.15
    nacre_threat_band = np.select(
        [nacre_suitability >= NACRE_DANGER_THRESHOLD, nacre_suitability < NACRE_SAFE_THRESHOLD],
        [2, 0], default=1,
    ).astype(np.int16)  # 0=safe, 1=graduated, 2=danger

    # --- v10 Skydrifter flight corridor: RETIRED, v11 (Nico's call) ---
    # v10 computed a mountain-avoiding least-cost path (see git history / decisions doc for the
    # full method) and buffered it into skydrifter_route_area_120m.npy. Nico's real-world read,
    # once shown the result: a species that overflies the whole crossing "mostly instantaneously"
    # and never stops has no behavioral reason to detour around terrain the way a ground-bound or
    # even a soaring species would -- the curved, terrain-hugging route implied more deliberate
    # routing than the locked lore supports. Retired the computed raster entirely; replaced with a
    # plain descriptive lore line in scenario_reference.md (wide, roughly straight band over the
    # North plains region, position varying flight to flight) -- no GIS layer needed for this one.
    # skydrifter_route_area_120m.npy/.tif are now STALE -- delete them locally alongside the
    # retired reaper_* files (see decisions doc "still open"/cleanup notes).

    # --- TWINSHADOWS ---
    # v2 fix (Nico's call): a real cougar's territory runs 150-1000 km2, but Temperate
    # Forest total area is only 689 km2 -- the original strict-forest exclusao could fit
    # under 1 to ~4.6 individuals in its ENTIRE range, which is why the layer read as an
    # insignificant patch rather than a territorial predator. Cougars use dense forest for
    # cover/ambush but routinely range and hunt into adjacent edge habitat -- extended
    # exclusao to Temperate Forest + Woodland/Shrubland within a buffer of the forest edge
    # (mirrors the ecotone concept already used for Snaketail/Furypack). Buffer distance:
    # 10km, the middle of a ~7-18km characteristic territory radius implied by 150-1000km2
    # (radius = sqrt(area/pi)) -- a documented estimate, not a precise figure.
    low_light = 1.0 - solar_suit
    temperate_forest_mask = biome_id == 7
    TWIN_EDGE_BUFFER_KM = 10.0
    dist_to_temp_forest_km = ndimage.distance_transform_edt(~temperate_forest_mask) * ((res_x + res_y) / 2) / 1000.0
    twin_exclusao = temperate_forest_mask | ((biome_id == 8) & (dist_to_temp_forest_km <= TWIN_EDGE_BUFFER_KM))
    twin_nucleo = 0.40 * low_light + 0.35 * preybase_twin + 0.25 * ruggedness
    twin_suitability = np.where(twin_exclusao, twin_nucleo, 0.0).astype(np.float32)

    # --- Alpine resource outposts ---
    # v1 (40% resource / 25% slope / 35% climate, 1-falloff-radius search) always
    # maximized toward valley-floor Grassland/Woodland -- climate+slope outvoted
    # resource proximity, so the composite reproduced the Circulo suitability
    # surface instead of siting anything alpine-specific. Fixed by (a) restricting
    # the search zone to the real Subalpine biome band (ids 4/5/6) instead of a
    # flat km radius, and (b) reweighting so resource proximity dominates
    # (60/15/25). Also fixed along the way: Permanent Snow & Ice sits outside
    # Nacre's exclusao mask (which only covers Alpine Fellfield/Alpine Tundra)
    # but is obviously unbuildable -- excluded explicitly rather than relying on
    # the exclusao mask alone.
    #
    # v5 (Nico's call, image review): v2-v4's Subalpine-band search zone sited every
    # candidate at the BASE of the alpine region, never inside it -- "some at the
    # peak is fine, but not all should be at the base." Fixed by moving the search
    # zone itself into Alpine Fellfield/Alpine Tundra (ids 2/3). Real consequence
    # surfaced immediately: ids 2/3 ARE, cell-for-cell, Nacre's own exclusao mask
    # (confirmed: alpine_mask.sum() == nacre_exclusao.sum(), both 1012.5-1012.6km2)
    # -- every outpost now sits INSIDE Nacre's core territory, so the original
    # status-exposure signal (distance to the whole massif MASK) saturates to
    # ~1.0-1.3 for 100% of alpine cells, killing all status differentiation.
    # Nico's chosen fix: derive outpost-specific exposure from distance to the
    # ridge CREST LINE instead (dist_spine_km/dist_sbridge_km, already computed
    # below for resource_suit -- zero new geometry), not from massif-mask
    # membership. Real range within the alpine band: 0-14.7km on the main spine,
    # 0-8.8km on South Branch -- plenty of spread. This is deliberately kept
    # SEPARATE from nacre_exposure_120m.npy (Tier 2, the seasonal-incursion
    # raster) -- that raster's own semantics (decay from the whole territory
    # outward, for LOWLAND incursion probability) are still correct and untouched;
    # only the per-outpost status-differentiation signal changes. Narrative
    # consequence, not a modeling accident: outposts closest to the ridge crest
    # (the same places resource_suit favors most) now also read as the most
    # exposed -- richest sites are the most dangerous, a real trade-off.
    def rasterize_line(coords):
        pts = []
        for i in range(len(coords) - 1):
            x0, y0 = coords[i]
            x1, y1 = coords[i + 1]
            n_steps = max(1, int(math.hypot(x1 - x0, y1 - y0) / 60))
            for t in np.linspace(0, 1, n_steps + 1):
                r, c = xy_to_rc(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t)
                pts.append((int(round(r)), int(round(c))))
        m = np.zeros((ny, nx), dtype=bool)
        for r, c in pts:
            if 0 <= r < ny and 0 <= c < nx:
                m[r, c] = True
        return m

    # Search-zone radius v2: widened from 8km to the full falloff_km=15 radius after
    # Nico compared candidate density to real alpine hut networks (Aoraki/Mt Cook NP:
    # 15 huts / ~715km2 =~ 1 per 48km2). The 8km zone (182.9km2) already implied ~3.8
    # outposts at that density -- matching the original 4 almost exactly -- so "more
    # outposts" required widening the zone, not just loosening min-separation. 15km
    # reuses the SAME falloff_km already locked for schist/jade siting elsewhere in
    # the doc, not a new arbitrary number; it raises the zone to 852.1km2, implying
    # ~17.9 at Aoraki density -- rounded to a target of 18.
    FALLOFF_KM = 15.0  # both Spine and South Branch share falloff_km=15 (geology section)
    dist_spine_km = ndimage.distance_transform_edt(
        ~rasterize_line(ridge_by_name['Spine']['geometry']['coordinates'])
    ) * ((res_x + res_y) / 2) / 1000.0
    dist_sbridge_km = ndimage.distance_transform_edt(
        ~rasterize_line(ridge_by_name['South Branch']['geometry']['coordinates'])
    ) * ((res_x + res_y) / 2) / 1000.0
    resource_suit = np.maximum(np.exp(-dist_spine_km / FALLOFF_KM), np.exp(-dist_sbridge_km / FALLOFF_KM)).astype(np.float32)
    climate_mildness = np.clip(biotemp / 12.0, 0, 1).astype(np.float32)

    RESOURCE_W, SLOPE_W, CLIMATE_W = 0.60, 0.15, 0.25
    outpost_composite = (RESOURCE_W * resource_suit + SLOPE_W * slope_suit + CLIMATE_W * climate_mildness).astype(np.float32)

    # v5 (Nico's call): search zone moved from the Subalpine band (ids 4-6, "base of the
    # alpine region") into Alpine Fellfield/Alpine Tundra itself (ids 2-3, "within the alpine
    # region") -- see the long comment above this section for why that also required a new
    # exposure signal. SEARCH_RADIUS_KM kept as a belt-and-suspenders bound (it's near-
    # non-binding now: max dist_spine_km within the alpine band is 14.7km, already under 15).
    #
    # v7 (Nico's call): search zone widened again to include Permanent Snow & Ice (id 1) --
    # "distribute the huts between the two alpine biomes AND the Permanent Ice biome." This
    # explicitly REVERSES the v1-era correction that excluded PSI as "obviously unbuildable"
    # -- checked directly before implementing rather than assuming: 41% of PSI area has
    # slope_suit > 0.3 (comparable to the 26.9% already found suitable within Alpine
    # Fellfield/Tundra back when v5 was being scoped), so real buildable ground exists there
    # too, just colder (mean biotemp 0.6C vs ~1.7C in ids 2-3, climate_mildness maxes out at
    # only 0.196 in PSI vs higher elsewhere) -- resource proximity (60% weight) still
    # dominates the composite, so PSI candidates are feasible, just colder-scoring on
    # average, not automatically excluded. Deliberately a SEPARATE mask from Nacre's own
    # alpine_mask (ids 2-3 only, unchanged) -- Permanent Snow & Ice has no vegetation/prey
    # base, so it stays out of the cave/cove/wild-yak model, this widening is outposts-only.
    OUTPOST_BIOME_IDS = [1, 2, 3]
    outpost_biome_mask = np.isin(biome_id, OUTPOST_BIOME_IDS)
    SEARCH_RADIUS_KM = 15.0
    search_mask = land & outpost_biome_mask & ((dist_spine_km <= SEARCH_RADIUS_KM) | (dist_sbridge_km <= SEARCH_RADIUS_KM))

    # Candidate count split proportional to each massif's own search-zone area
    # (not the fixed 3:1 pack-count mirror used at n=4) -- more defensible once the
    # zone is wide enough that area, not pack count, should drive the split.
    TARGET_TOTAL = 18
    area_ms = (search_mask & ~zone_south).sum() * res_x * res_y / 1e6
    area_sb = (search_mask & zone_south).sum() * res_x * res_y / 1e6
    n_sb = max(1, round(TARGET_TOTAL * area_sb / (area_ms + area_sb)))
    n_ms = TARGET_TOTAL - n_sb

    # v4 fix (Nico's call): v3's pick_maxima ran main-spine and South-Branch picks as two
    # INDEPENDENT greedy searches, each only separation-checked against its own zone's
    # picks -- both zones converged on the same hot spot where the Spine and South Branch
    # ridges' 15km falloff radii overlap (resource_suit is boosted there by BOTH ridges at
    # once), so most candidates piled into one strip instead of spreading along the ~50km
    # Spine. Fixed two ways: (1) bin each zone into n_target equal-arclength segments along
    # its OWN ridge and pick one candidate per segment, guaranteeing spread along the full
    # ridge length; (2) separation is now checked GLOBALLY across both zones, not just
    # within each one.
    MIN_SEP_KM = 3.0
    spine_arclen, spine_total_km = rasterize_line_with_arclength(
        ridge_by_name['Spine']['geometry']['coordinates'], ny, nx, xy_to_rc)
    sbridge_arclen, sbridge_total_km = rasterize_line_with_arclength(
        ridge_by_name['South Branch']['geometry']['coordinates'], ny, nx, xy_to_rc)

    # v7 (Nico's call): "distribute the huts between the two alpine biomes and the
    # Permanent Ice biome" -- checked what plain area-proportional + greedy-per-arclength-
    # bin picking would do first, and it skewed hard: 13/18 candidates landed in Permanent
    # Snow & Ice (72%), because PSI cells sit directly on parts of the ridge crest
    # (maximizing resource_suit, 60% of the composite) and climate's 25% weight barely
    # penalizes it relative to Alpine Fellfield/Tundra. That's a real result, not a bug, but
    # it reads as "PSI-dominated" rather than "distributed between" the three biomes -- so
    # each massif's target count is now split as evenly as possible across whichever of the
    # 3 biomes actually has area in that massif's search zone (base n//k + remainder to the
    # first biomes, in Permanent Snow & Ice / Alpine Fellfield / Alpine Tundra order), with
    # each biome's picks still spread along the ridge via pick_along_ridge and separation
    # checked globally across the whole outpost set.
    def split_targets_by_biome(n_target, zone_mask):
        available = [b for b in OUTPOST_BIOME_IDS if (zone_mask & (biome_id == b)).any()]
        if not available or n_target <= 0:
            return {}
        base, rem = divmod(n_target, len(available))
        return {b: base + (1 if i < rem else 0) for i, b in enumerate(available)}

    def pick_balanced_by_biome(zone_mask, arclen_raster, total_km, n_target, taken):
        chosen = []
        for b, n_b in split_targets_by_biome(n_target, zone_mask).items():
            if n_b <= 0:
                continue
            picked = pick_along_ridge(
                zone_mask & (biome_id == b), outpost_composite, arclen_raster, total_km,
                n_b, MIN_SEP_KM, res_x, res_y,
                taken=list(taken) + [(r, c) for r, c, _ in chosen])
            chosen.extend(picked)
        return chosen

    main_spine_candidates = pick_balanced_by_biome(
        search_mask & (~zone_south), spine_arclen, spine_total_km, n_ms, taken=[])
    south_branch_candidates = pick_balanced_by_biome(
        search_mask & zone_south, sbridge_arclen, sbridge_total_km, n_sb,
        taken=[(r, c) for r, c, _ in main_spine_candidates])

    # v6 (supersedes v5's ad-hoc ridge-distance signal): now that Nacre's own
    # suitability is a single graduated field defined everywhere (not just inside a
    # hard mask), sample it directly at each outpost and reuse Nico's own thresholds
    # (>=0.9 danger, <0.15 safe) instead of a separate relative min-max banding. One
    # real threat signal instead of two parallel ones.
    def band_status(e):
        if e >= NACRE_DANGER_THRESHOLD:
            return 'abandoned'
        elif e < NACRE_SAFE_THRESHOLD:
            return 'active'
        return 'temporary refuge'

    all_candidates = [(r, c, v, 'main_spine', i) for i, (r, c, v) in enumerate(main_spine_candidates)] + \
                      [(r, c, v, 'south_branch', i) for i, (r, c, v) in enumerate(south_branch_candidates)]
    exposure_vals = [float(nacre_suitability[r, c]) for r, c, v, massif, i in all_candidates]

    features = []
    for (r, c, v, massif, i), ev in zip(all_candidates, exposure_vals):
        x, y = rc_to_xy(r, c)
        name = f"Outpost_{'MainSpine' if massif=='main_spine' else 'SouthBranch'}_{i+1}"
        features.append({
            'type': 'Feature',
            'properties': {
                'name': name, 'massif': massif,
                'composite_suitability': round(float(v), 3),
                'resource_suit': round(float(resource_suit[r, c]), 3),
                'slope_suit': round(float(slope_suit[r, c]), 3),
                'climate_mildness': round(float(climate_mildness[r, c]), 3),
                'biome': BIOME_NAMES[int(biome_id[r, c])],
                'nacre_suitability_at_site': round(ev, 3),
                'status_prior': band_status(ev),
                'status_is_authorial_final': False,
            },
            'geometry': {'type': 'Point', 'coordinates': [x, y]},
        })
    outpost_geojson = {'type': 'FeatureCollection', 'name': 'tappa7_alpine_outposts', 'features': features}

    # --- Domesticated points/zones ---
    povo_livre_zone = povo_livre > 0.80

    # --- Coves (Nacre dens) geojson ---
    cove_features = []
    for idx, (r, c, v, massif, i) in enumerate(all_coves):
        x, y = rc_to_xy(r, c)
        cove_features.append({
            'type': 'Feature',
            'properties': {
                'name': f"Cove_{'MainSpine' if massif=='main_spine' else 'SouthBranch'}_{i+1}",
                'pack_id': idx + 1, 'massif': massif,
                'cove_score': round(float(v), 3),
                'slope_pct': round(float(slope_pct[r, c]), 1),
                'dist_to_stream_km': round(float(dist_stream_km[r, c]), 3),
                'dist_to_subalpine_km': round(float(dist_to_subalpine_km[r, c]), 3),
                'biome': BIOME_NAMES[int(biome_id[r, c])],
            },
            'geometry': {'type': 'Point', 'coordinates': [x, y]},
        })
    coves_geojson = {'type': 'FeatureCollection', 'name': 'tappa7_nacre_coves', 'features': cove_features}

    import os
    os.makedirs(OUT, exist_ok=True)
    np.save(f'{OUT}/nacre_suitability_120m.npy', nacre_suitability)
    np.save(f'{OUT}/twinshadows_suitability_120m.npy', twin_suitability)
    np.save(f'{OUT}/nacre_threat_band_120m.npy', nacre_threat_band)
    np.save(f'{OUT}/nacre_pack_territory_120m.npy', pack_territory)
    np.save(f'{OUT}/yak_suitability_120m.npy', wild_yak_suit)
    np.save(f'{OUT}/outpost_composite_120m.npy', outpost_composite)
    np.save(f'{OUT}/outpost_search_mask_120m.npy', search_mask)
    np.save(f'{OUT}/povo_livre_zone_120m.npy', povo_livre_zone)
    np.save(f'{OUT}/nacre_grassmother_incursion_120m.npy', nacre_grassmother_incursion)
    json.dump(outpost_geojson, open(f'{OUT}/outpost_candidates.geojson', 'w'), indent=2)
    json.dump(coves_geojson, open(f'{OUT}/nacre_coves.geojson', 'w'), indent=2)

    n_grassland = sum(1 for r, c, b in circulo_rc.values() if b == 'Lowland Steppe / Grassland')
    run_meta = {
        'nacre': {
            'version': 'v6 -- Nico\'s redesign, replaces the exclusao/Tier-1-Tier-2 architecture entirely.',
            'method': 'Caves (steep + near-stream within the alpine band) -> 4 coves/dens (3 main-spine, '
                      '1 South Branch, min-sep 15km) -> wild yak prey base (biome x slope, independent of '
                      'Circulo/outposts) -> a single graduated suitability field, land-masked only, no hard '
                      'biome exclusao.',
            'weights': {'dist_to_cove': NACRE_W_COVE, 'ruggedness': NACRE_W_RUGGED,
                        'wild_yak_prey': NACRE_W_PREY, 'alpine_affinity': NACRE_W_ALPINE},
            'weights_note': 'Proposed defaults, not independently re-confirmed line by line with Nico -- '
                             'adjustable.',
            'dist_to_cove_method': 'v8 fix (Nico\'s catch): dist_to_cove was pure Euclidean straight-line '
                                     'distance through v7 -- couldn\'t distinguish an easy walk from a line '
                                     'drawn straight through a ravine/cliff. Replaced with real slope-weighted '
                                     'cost-distance (skimage.graph.MCP_Geometric, friction = 1.0 + '
                                     f'{SLOPE_FRICTION_WEIGHT}*ruggedness per cell), same surface reused for the '
                                     'Grassmother incursion layer below. This is 35% of nacre_suitability\'s '
                                     'weight, so it materially changed which outpost sites read as close to a '
                                     'cove -- see outpost status counts for the before/after.',
            'calibration': 'Raw weighted composite land-wide max was 0.82 -- >=0.9 was unreachable before '
                            'rescaling, which would have made every outpost read as the same status band. '
                            f'Rescaled by land-wide p99 (raw={round(float(nacre_p99), 3)}) so the thresholds '
                            'are real; relative ordering preserved.',
            'pack_territory_method': 'v8c fix (Nico\'s follow-up: "a change in the pack would also be more '
                                       'coherent, right?"): territory assignment was Euclidean-nearest-cove '
                                       '(straight-line Voronoi) through v8b, even after dist_to_cove itself was '
                                       'fixed -- same blind spot, a cell could read as closest to one cove by '
                                       'ruler distance while actually being cheaper to reach from a different '
                                       'cove once a ravine/cliff on the straight-line path is accounted for. '
                                       'Real parallel: carnivore territory/home-range delineation is standardly '
                                       'done with cost-distance/resistance surfaces in wildlife ecology, not '
                                       'straight-line Voronoi. Fixed by running the same slope-weighted cost-'
                                       'distance separately from each of the 4 coves (one MCP_Geometric call '
                                       'per cove, since the multi-source version used for dist_to_cove only '
                                       'returns the cost to the NEAREST cove, not which one) and assigning each '
                                       'cell to its argmin. Real effect: 6.14% of land cells reassigned to a '
                                       'different pack\'s territory, including a real coastal-strip flip near '
                                       'the SW island from South Branch to the westernmost main-spine pack. '
                                       'Honest trade-off, not hidden: cost-distance boundaries are noisier than '
                                       'the old clean Euclidean lines -- each territory keeps 81-99.5% of its '
                                       'area as one contiguous block, with minor speckling (small patches under '
                                       '10 cells, ~1.2km2 or less, invisible at map scale) at contested edges '
                                       'where local terrain noise creates locally-cheaper detached patches; not '
                                       'smoothed/cleaned, left as the model\'s honest output.',
            'decay_km': DECAY_KM, 'south_branch_weight': SOUTH_BRANCH_WEIGHT,
            'threat_bands': {'danger_threshold': NACRE_DANGER_THRESHOLD, 'safe_threshold': NACRE_SAFE_THRESHOLD,
                              'description': '>=danger_threshold reads as untenable/"abandoned"-tier risk for '
                                              'human presence; <safe_threshold reads as comparatively safe; '
                                              'in between is graduated threat -- for reading/status derivation, '
                                              'not a new hard cutoff.'},
            'land_area_mean_suitability': round(float(nacre_suitability[land].mean()), 3),
            'land_area_pct_in_danger_band': round(float(100 * (nacre_threat_band[land] == 2).mean()), 2),
            'land_area_pct_in_safe_band': round(float(100 * (nacre_threat_band[land] == 0).mean()), 2),
            'cave_candidates': {
                'slope_percentile_within_alpine': CAVE_SLOPE_PCTILE,
                'stream_proximity_km': CAVE_STREAM_KM,
                'honest_finding': '99.9% of the alpine band is already within 1km of a mapped stream (same '
                                   'dense-hydrology pattern already documented elsewhere in this doc), so the '
                                   'stream-proximity criterion barely discriminates candidates here -- the cave '
                                   'score is driven almost entirely by slope in practice.',
                'n_candidate_cells': int(cave_candidate_mask.sum()),
            },
            'coves': {
                'n_main_spine': len(main_spine_coves), 'n_south_branch': len(south_branch_coves),
                'min_sep_km': COVE_MIN_SEP_KM,
                'scoring': {'ruggedness_weight': COVE_W_RUGGED, 'interior_depth_weight': COVE_W_INTERIOR,
                            'interior_depth_norm_km': COVE_INTERIOR_NORM_KM},
                'v8_lore_fix': 'Nico caught this: v6 coves sat on the Subalpine edge (3 of 4 within 0.12km), '
                               'which would make Blacknose/Snaketail/Furypack routine prey rather than the '
                               'locked "baseline diet is yaks, rare seasonal descents" lore. Added an interior-'
                               'depth term (distance to nearest Subalpine cell specifically -- NOT penalizing '
                               'proximity to Permanent Snow & Ice, which is lore-consistent) to the cove score.',
                'file': 'nacre_coves.geojson',
            },
            'wild_yak': {
                'file': 'yak_suitability_120m.npy',
                'biome_habitat': 'Alpine Fellfield/Tundra + Subalpine Wet Forest/Woodland/Dry Scrub (ids 2-6)',
                'slope_curve': 'full suitability to 15% grade, exp decay with 60%-grade characteristic scale '
                                'beyond that -- much more permissive than the human slope_suitability_120m '
                                'curve (real yaks are sure-footed high-altitude grazers)',
                'not_tied_to': 'Circulo sites or outposts -- an ambient baseline population, independent of the '
                                'existing "yak goes feral at an abandoned outpost" point-source narrative, '
                                'which still stands as additional local color, not replaced by this.',
                'mean_suitability_within_habitat': round(float(wild_yak_suit[np.isin(biome_id, [2,3,4,5,6])].mean()), 3),
            },
            'grassmother_incursion': {
                'file': 'nacre_grassmother_incursion_120m.npy',
                'purpose': 'Nico\'s ask: model how far Nacre packs would plausibly penetrate into Lowland '
                           'Steppe/Grassland from their coves to hunt Grassmothers, considering slope, as a '
                           'decreasing-value seasonal zone (early spring) rather than a hard boundary.',
                'target_biome': 'Lowland Steppe / Grassland (id 9) only -- Grassmothers are locked at 100% '
                                 'there, with only marginal (15%) presence in Subalpine Dry Scrub/Woodland-'
                                 'Shrubland, which read as edge habitat, not real hunting ground, so excluded '
                                 'from the target mask (Nico\'s choice).',
                'method': 'Real slope-weighted cost-distance from each cove (skimage.graph.MCP_Geometric), not '
                          'a masked Euclidean decay like the rest of this pipeline -- Nico asked specifically '
                          'for slope to shape penetration depth. Reuses the exact same cost surface as '
                          'dist_to_cove_km above (friction = 1.0 + SLOPE_FRICTION_WEIGHT*ruggedness per cell, '
                          'SLOPE_FRICTION_WEIGHT=2.0) -- computed once, not twice.',
                'decay_anchor': 'GM_DECAY_KM=15.0, from real wolf denning/pup-rearing-season central-place-'
                                 'foraging literature (routine trips commonly ~10-20km from den, rarer longer '
                                 'forays beyond) -- NOT a reused pipeline constant (DECAY_KM=20/FALLOFF_KM=15 '
                                 'were both derived for unrelated geometry, ridge falloff and schist radius).',
                'kept_separate': 'Deliberately NOT blended into nacre_suitability_120m.npy -- the 0.9/0.15 '
                                  'threat bands and outpost status_prior are already locked by Nico and this '
                                  'doesn\'t touch them; this is an additional, seasonal layer.',
                'mean_suitability_within_grassland': round(float(nacre_grassmother_incursion[grassland_mask].mean()), 4),
                'pct_grassland_above_0.1': round(float(100 * (nacre_grassmother_incursion[grassland_mask] > 0.1).mean()), 2),
            },
            'superseded': [
                'Tier 1 (hard exclusao = Alpine Fellfield/Tundra) and Tier 2 (nacre_exposure_120m.npy, '
                'distance-decay from the massif mask) are both retired -- replaced by the single graduated '
                'field above. nacre_massif_labels_120m.npy (2-zone) is retired too, replaced by '
                'nacre_pack_territory_120m.npy (4-zone, nearest-cove assignment, one per real den).',
                'Outpost status_prior now samples nacre_suitability directly at each site (see outposts '
                'section) instead of the v5 ad-hoc ridge-crest-distance signal -- one real threat signal '
                'instead of two parallel ones.',
            ],
        },
        'skydrifter': {
            'status': 'v10 computed-corridor approach RETIRED (v11, Nico\'s call) -- no raster/GIS layer '
                      'for this species. A species that overflies the whole crossing "mostly '
                      'instantaneously" and never stops has no behavioral reason to detour around terrain '
                      '(unlike Nacre\'s deliberate cost-distance movement) -- the v10 route\'s '
                      'terrain-hugging curve implied more routing intent than the locked lore supports. '
                      'Replaced with a plain descriptive line in scenario_reference.md (wide, roughly '
                      'straight band over the North plains, varies flight to flight) -- see that doc, '
                      'not this meta block, for the current lore.',
        },
        'twinshadows': {
            'exclusao_area_km2': round(float(twin_exclusao.sum() * res_x * res_y / 1e6), 1),
            'exclusao': 'Temperate Forest UNION Woodland/Shrubland within 10km of the forest edge '
                        '(v2 fix -- was Temperate Forest only, 689 km2, which could fit under 1 to ~4.6 '
                        'real cougar-scale territories (150-1000 km2 each) in its ENTIRE range).',
            'edge_buffer_km': TWIN_EDGE_BUFFER_KM,
            'weights': {'low_light': 0.40, 'prey_base': 0.35, 'ruggedness': 0.25},
            'suitability_mean_within_mask': round(float(twin_nucleo[twin_exclusao].mean()), 3),
        },
        'pack_to_massif': {
            'main_spine_km2': round(float(main_spine_mask.sum() * res_x * res_y / 1e6), 1),
            'south_branch_km2': round(float(south_branch_mask.sum() * res_x * res_y / 1e6), 1),
            'ratio': round(float(main_spine_mask.sum() / south_branch_mask.sum()), 2),
            'n_components_found': int(n),
            'decay_km': DECAY_KM, 'south_branch_weight': SOUTH_BRANCH_WEIGHT,
        },
        'outposts': {
            'weights_final': {'resource': RESOURCE_W, 'slope': SLOPE_W, 'climate': CLIMATE_W},
            'search_zone': f'Permanent Snow & Ice + Alpine Fellfield/Alpine Tundra (ids 1-3), within {SEARCH_RADIUS_KM:.0f}km of Spine or South Branch ridge axis (v7 -- see corrections_applied)',
            'search_zone_area_km2': {'main_spine': round(float(area_ms), 1), 'south_branch': round(float(area_sb), 1)},
            'falloff_km': FALLOFF_KM, 'min_sep_km': MIN_SEP_KM,
            'n_candidates': len(features), 'n_target_main_spine': n_ms, 'n_target_south_branch': n_sb,
            'n_filled_main_spine': len(main_spine_candidates), 'n_filled_south_branch': len(south_branch_candidates),
            'status_exposure_signal': 'v6: nacre_suitability_120m.npy sampled directly at each site, against '
                                       'Nico\'s 0.9/0.15 thresholds -- see the nacre meta block above.',
            'unfilled_bins_note': (
                f'{n_ms - len(main_spine_candidates)} main-spine and {n_sb - len(south_branch_candidates)} '
                'arclength bins had no valid cell within the search zone (Permanent Snow & Ice + Alpine '
                'Fellfield/Tundra) + global separation constraint and were left empty rather than silently '
                'backfilled -- real gaps in the along-ridge distribution, not a bug.'
            ) if (n_ms - len(main_spine_candidates) or n_sb - len(south_branch_candidates)) else None,
            'real_world_benchmark': {
                'source': 'Aoraki/Mount Cook National Park, NZ (707-722 km2, 15 huts)',
                'km2_per_hut': round(715 / 15, 1),
                'v1_zone_182.9km2_implied_outposts': round(182.9 / (715 / 15), 1),
                'v2_zone_implied_outposts': round((area_ms + area_sb) / (715 / 15), 1),
                'caveat': 'Real hut density reflects tourist day-hike spacing, not resource-extraction/'
                          'research-station siting (which is real-world sparser) -- used as an order-of-'
                          'magnitude anchor, not a literal target, per Nico\'s explicit steer to widen the zone.',
            },
            'corrections_applied': [
                'v1 weights (40/25/35) + 15km flat radius always favored valley-floor Grassland/Woodland '
                '(climate+slope outvoted resource proximity) -- not usable, model reproduced Circulo suitability instead.',
                'Permanent Snow & Ice explicitly excluded from the search mask (sits outside the Alpine Fellfield/'
                'Alpine Tundra exclusao mask but is obviously unbuildable).',
                'v2 (Nico\'s call): search zone restricted to real Subalpine biome band + reweighted to 60/15/25.',
                'v3 (Nico\'s call): search radius widened 8km->15km (reuses the already-locked falloff_km, not '
                'a new number) after comparing candidate density to Aoraki/Mt Cook NP\'s real hut density -- '
                'raised target count 4->18, split by massif search-zone area instead of a fixed 3:1 pack-count '
                'mirror, min-separation tightened 5km->3km so 18 slots are actually fillable.',
                'v4 BUG FIX (Nico caught it): v3 ran main-spine and South-Branch picks as two independent '
                'greedy searches, each only checking separation against its own zone -- both converged on the '
                'same ridge-overlap hot spot (Spine and South Branch fall-off radii overlap there, boosting '
                'resource_suit for both), clustering most candidates into one strip instead of along the whole '
                'Spine. Fixed: bin each zone into equal-arclength segments along its own ridge and pick one '
                'candidate per segment (guarantees spread), plus a GLOBAL separation check across both zones.',
                'v5 (Nico\'s call, image review): v2-v4 sited every candidate at the BASE of the alpine region '
                '(Subalpine band, ids 4-6) -- "some at the peak is fine, but not all should be at the base." '
                'Search zone moved into Alpine Fellfield/Alpine Tundra itself (ids 2-3). Real consequence found '
                'immediately: ids 2-3 are, cell for cell, Nacre\'s own exclusao mask (both ~1012.5km2) -- every '
                'outpost now sits inside Nacre\'s core territory, so the original status signal (distance to '
                'the whole massif mask, nacre_exposure_120m.npy) saturates to ~1.0-1.3 for 100% of alpine cells, '
                'killing all status differentiation. Nico\'s chosen fix, from 3 options presented: derive a '
                'SEPARATE outpost-only exposure signal from distance to the ridge crest line (dist_spine_km/'
                'dist_sbridge_km, already computed for resource_suit -- zero new geometry), real range 0-14.7km '
                'on the main spine. nacre_exposure_120m.npy (Tier 2, lowland seasonal-incursion probability) is '
                'untouched -- this only changes how outpost status_prior is derived. Narrative consequence, not '
                'an accident: the ridge-closest sites (best resource_suit) now also read as most exposed -- '
                'richest sites are the most dangerous.',
                'v7 (Nico\'s call): search zone widened again to include Permanent Snow & Ice (id 1), '
                'reversing the original v1-era exclusion of PSI as "obviously unbuildable" -- checked before '
                'implementing: 41% of PSI has slope_suit > 0.3 (comparable to the 26.9% already found in ids '
                '2-3), so real buildable ground exists there too, just colder on average (climate_mildness '
                'maxes at 0.196 in PSI). Kept as a separate outpost-only mask (OUTPOST_BIOME_IDS) -- Nacre\'s '
                'own alpine_mask (ids 2-3) is unchanged, since PSI has no vegetation/prey base and stays out '
                'of the cave/cove/wild-yak model.',
            ],
        },
        'domesticated_points_zones': {
            'n_circulo_grassland': n_grassland, 'n_circulo_total': 17,
            'povo_livre_zone_km2': round(float(povo_livre_zone.sum() * res_x * res_y / 1e6), 1),
            'povo_livre_zone_pct_land': round(float(100 * povo_livre_zone.sum() / land.sum()), 1),
        },
    }
    json.dump(run_meta, open(f'{OUT}/tappa7_fauna_composite_meta.json', 'w'), indent=2)
    print(json.dumps(run_meta, indent=2))


if __name__ == '__main__':
    main()
