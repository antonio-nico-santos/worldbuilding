"""
Tappa 9 -- road network topology and route geometry over the 17 already-
placed, LOCKED Circulo sites (Tappa 6; permanently not re-placed under any
downstream friction layer -- see `08_tappa8_geomorphology.md` S12).

REVISED after Nico's direct review of the first pass (2026-08-20). Two
problems in that pass, both fixed here:

1. **The graph allowed open-water crossings.** The first pass ran the MST
   over `cost_distance.py`'s general-purpose graph, which -- correctly,
   for its original Tappa 6 use case (isochrone/tier-distance siting) --
   treats any non-land cell (ocean OR lake) as boat-traversable at a flat
   speed. For an actual ROAD network that's wrong: 6 of 16 edges ran
   straight across open sea, and roughly half of the rest crossed a lake,
   because a road can't just be "built" over open water the way a distance
   check can cross it. `run_tappa9_road_network.py` now builds this
   module's graph with `sea_mode="impassable"` (see `cost_distance.py`) --
   land-only, no boat fallback at all.
2. **That means some site pairs may be genuinely unreachable by road.**
   Cutting sea/lake edges can split the 17 sites into more than one
   connected component (an island Circulo, or one separated from the
   mainland cluster by a lake with no land route around it). The original
   `build_mst_edges` assumed one connected graph and would silently break
   (or, worse, look like it worked while quietly building on top of
   `np.inf` costs) if that assumption stopped holding. This module now
   builds a MINIMUM SPANNING FOREST (`build_mst_forest`): one MST per
   connected component, plus an explicit list of which sites ended up in
   which component -- any component with more than one member is a
   genuine "needs a ferry to reach the rest of the network" finding for
   the caller to surface, not a bug to hide.

TOPOLOGY DECISION: within each connected component, connect sites with a
MINIMUM SPANNING TREE over pairwise cost-distance (hours) on the combined
(lithology x biome x river-crossing) friction-adjusted, land-only graph --
the minimum-total-travel-time way to connect every site in the component.
On top of that skeleton, `add_redundant_edges` adds EXTRA edges -- this is
what fixes Nico's "always linear, even where several Circulos cluster
together" complaint from the first review.

**REVISED AGAIN after Nico's SECOND review (2026-08-20, same day):** the
redundancy pass's first version (relative to each site's own cheapest
neighbour, `redundancy_factor=1.4` ALONE) added 8 edges, 3 of which Nico
flagged directly by GeoJSON feature ID as excessive -- two cost almost
EXACTLY what the existing MST tree path between the same two sites
already cost (pure clutter, a second line on the map that gets you
nowhere faster). `add_redundant_edges` now ALSO requires a genuine
tree-path-shortcut improvement (`min_shortcut_improvement`, default
0.20) on top of the original candidate scoping -- see that function's own
docstring for the two intermediate attempts (an O(n^2) candidate search,
and a straight swap instead of an AND) that were tested and rejected
before landing on requiring both conditions together. This is what a
"redundant connection actually helps" check should measure; the original
rule never checked whether the edge shortened anything over the tree at
all.

Still a heuristic, not a formal network-design optimization (Steiner tree,
facility-location, etc.) -- picked for being cheaply explainable and
directly responsive to Nico's specific complaints across two review
rounds, not because it's the only or the "correct" way to add redundancy.

DIRECTED-GRAPH NOTE: cost_distance.py's graph is directed (uphill != downhill
Tobler cost over the same edge), so the raw pairwise-hours matrix is not
symmetric. Both the MST and the redundancy pass need one scalar weight per
UNORDERED pair, so this module uses the mean of the two directions
throughout -- the natural choice for a road that gets walked/ridden both
ways, not just once.
"""
from __future__ import annotations

import numpy as np

from suitability.cost_distance import (
    cost_distance_from_source_with_predecessors,
    reconstruct_path,
)


def compute_pairwise_cost_distance(sites: list[dict], graph, shape: tuple[int, int]):
    """All-pairs cost-distance (hours) among `sites` (each a dict with at
    least "row", "col"), by running one single-source Dijkstra PER SITE
    (n runs, not C(n,2) -- every other site's distance is read off each
    run's full-grid result for free) and keeping every run's predecessor
    array so paths can be reconstructed later without re-running Dijkstra.

    Returns (hours, predecessors_by_source) where `hours` is an (n, n)
    array (hours[i, j] = cost FROM site i TO site j -- NOT guaranteed
    symmetric, see module docstring) and `predecessors_by_source` is a
    list of n raw predecessor arrays, index-aligned with `sites`. On an
    "impassable"-sea graph, an unreachable pair reads back as np.inf --
    expected, not an error (see module docstring's point 2).
    """
    n = len(sites)
    hours = np.full((n, n), np.inf)
    predecessors_by_source = []
    for i, s in enumerate(sites):
        dist, pred = cost_distance_from_source_with_predecessors(
            graph, s["row"], s["col"], shape
        )
        predecessors_by_source.append(pred)
        for j, t in enumerate(sites):
            hours[i, j] = dist[t["row"], t["col"]]
    return hours, predecessors_by_source


def connected_components(hours: np.ndarray) -> list[list[int]]:
    """Group site indices 0..n-1 into connected components using the
    SYMMETRIZED cost matrix -- two sites are in the same component iff a
    finite-cost path connects them (in EITHER direction; a one-way-only
    connection would be a strange result for a walking/road graph and
    doesn't occur in practice here, but treating either direction as
    sufficient for "connected" is the conservative, inclusive choice for
    grouping purposes). Returns components as lists of indices, ordered by
    each component's smallest member index (component 0 contains site 0).
    """
    n = hours.shape[0]
    reachable = np.isfinite(hours) | np.isfinite(hours.T)
    seen = [False] * n
    components = []
    for start in range(n):
        if seen[start]:
            continue
        stack, comp = [start], []
        seen[start] = True
        while stack:
            node = stack.pop()
            comp.append(node)
            for nxt in range(n):
                if not seen[nxt] and reachable[node, nxt]:
                    seen[nxt] = True
                    stack.append(nxt)
        components.append(sorted(comp))
    return components


def build_mst_forest(hours: np.ndarray) -> tuple[list[tuple[int, int, float]], list[list[int]]]:
    """Prim's algorithm MST, run independently within EACH connected
    component of the symmetrized cost matrix (mean of hours[i,j] and
    hours[j,i]) -- a minimum spanning FOREST, not a single tree, since
    cutting sea/lake edges (see module docstring) can leave more than one
    component among the 17 sites. Returns (edges, components): `edges` is
    the concatenation of every component's own MST edges
    (i, j, weight_hours), i < j; `components` is `connected_components`'s
    own output, so the caller can tell which sites ended up isolated
    together (a component of size 1 = a site the land-only graph couldn't
    connect to ANY other site) without re-deriving it.
    """
    n = hours.shape[0]
    sym = 0.5 * (hours + hours.T)
    components = connected_components(hours)
    edges: list[tuple[int, int, float]] = []
    for comp in components:
        if len(comp) < 2:
            continue  # isolated site -- no edge to draw, flagged by the caller
        comp_set = set(comp)
        in_tree = {comp[0]}
        remaining = set(comp[1:])
        while remaining:
            best = None
            for i in in_tree:
                for j in remaining:
                    w = sym[i, j]
                    if best is None or w < best[2]:
                        best = (i, j, w)
            i, j, w = best
            edges.append((min(i, j), max(i, j), float(w)))
            in_tree.add(j)
            remaining.discard(j)
    return edges, components


def _tree_path_cost(mst_edges: list[tuple[int, int, float]], i: int, j: int) -> float | None:
    """Sum of edge weights along the (unique, since it's a tree) path from
    `i` to `j` through `mst_edges`. Returns None if `i`/`j` aren't
    connected through these edges at all (different component, or one of
    them isn't in the tree) -- the caller treats that as "not this
    function's job," not an error."""
    adj: dict[int, list[tuple[int, float]]] = {}
    for a, b, w in mst_edges:
        adj.setdefault(a, []).append((b, w))
        adj.setdefault(b, []).append((a, w))
    if i not in adj and i != j:
        return None
    visited = {i}
    stack = [(i, 0.0)]
    while stack:
        node, cost = stack.pop()
        if node == j:
            return cost
        for nxt, w in adj.get(node, []):
            if nxt not in visited:
                visited.add(nxt)
                stack.append((nxt, cost + w))
    return None


def add_redundant_edges(
    hours: np.ndarray,
    mst_edges: list[tuple[int, int, float]],
    components: list[list[int]],
    max_extra_per_site: int = 2,
    redundancy_factor: float = 1.4,
    min_shortcut_improvement: float = 0.20,
) -> list[tuple[int, int, float]]:
    """REVISED THREE TIMES, 2026-08-20, after two rounds of Nico's review.

    v1 (first review) added a redundant edge whenever a site had a second
    connection within `redundancy_factor` of its OWN cheapest one --
    purely relative to that one site's cheapest neighbour, with no check
    on whether the edge actually shortened anything. That produced three
    edges Nico flagged directly as excessive by GeoJSON feature ID (15,
    18, 20 in that run's export): two (`Circulo_A_40k`<->`Circulo_E4_2k`,
    `Circulo_E3_2k`<->`Circulo_F8_small`) turned out to cost almost
    EXACTLY the same as the path already available through the MST tree
    (0.1% and 0.003% "savings" -- pure visual clutter, a second line on
    the map that gets you nowhere faster); the third
    (`Circulo_C_25k`<->`Circulo_E2_2k`) saved a real but modest 16.2%.

    v2 (never shipped, caught in testing before delivery) REPLACED v1's
    acceptance rule with a tree-path-shortcut check but ALSO widened the
    candidate search to literally every non-MST pair in the graph
    (O(n^2)). That overshot badly the other way -- 26 edges for 16 sites,
    because a direct edge between two sites on opposite sides of the
    network will almost always beat their multi-hop tree-path sum by a
    wide margin (more hops accumulate more "detour"), a property of the
    ARITHMETIC, not evidence the edge is a locally useful redundant
    connection -- Nico's original complaint was about Circulos that
    cluster closely getting loops, not arbitrary long-distance chords
    across the whole map. Caught by actually running it and counting
    before it ever shipped.

    v3 (also tested, also not what shipped) kept v1's local candidate
    scoping (each site's `max_extra_per_site` next-cheapest neighbours
    only) but REPLACED the acceptance rule outright with the tree-path-
    shortcut check. Tested directly rather than assumed: this did NOT
    just prune v1's 8 candidates down to a subset -- it dropped the 3 bad
    ones but ALSO newly accepted 3 different candidates v1 had rejected
    (cheap enough in tree-shortcut terms, but not "close to that site's
    own cheapest neighbour" in v1's sense), landing at a different
    8-edge set nobody had reviewed yet. Not shipped specifically BECAUSE
    it introduced unreviewed churn rather than just removing the 3
    reported edges.

    **v4 (current, shipped): require BOTH v1's and v2's conditions
    together**, not one replacing the other -- `redundancy_factor` still
    scopes candidates to a site's genuinely near neighbours, and
    `min_shortcut_improvement` still requires the edge to be a real
    shortcut over the MST tree path, not just a nearby-cost restatement
    of it. Checked directly: this combination keeps exactly the 5 edges
    from v1's original 8 that Nico did NOT flag (21.0%-43.4% tree-path
    improvement) and drops exactly the 3 he did (0.003%-16.2%), while
    introducing NO new candidates beyond that original reviewed set --
    the conservative fix, not a redesign.

    Returns only the qualifying NEW edges (not already in `mst_edges`),
    each as (i, j, weight_hours) with i < j.
    """
    n = hours.shape[0]
    sym = 0.5 * (hours + hours.T)
    mst_edge_set = {(i, j) for i, j, _ in mst_edges}
    comp_of = {}
    for comp in components:
        for idx in comp:
            comp_of[idx] = frozenset(comp)

    extra: dict[tuple[int, int], float] = {}
    for i in range(n):
        same_comp = comp_of.get(i, frozenset())
        candidates = sorted(
            (sym[i, j], j) for j in same_comp if j != i and np.isfinite(sym[i, j])
        )
        if not candidates:
            continue
        cheapest_cost = candidates[0][0]
        for direct_cost, j in candidates[1 : 1 + max_extra_per_site]:
            edge = (min(i, j), max(i, j))
            if edge in mst_edge_set or edge in extra:
                continue
            if direct_cost > redundancy_factor * cheapest_cost:
                continue
            tree_cost = _tree_path_cost(mst_edges, i, j)
            if tree_cost is None or tree_cost <= 0:
                continue
            improvement = (tree_cost - direct_cost) / tree_cost
            if improvement >= min_shortcut_improvement:
                extra[edge] = direct_cost

    return [(i, j, w) for (i, j), w in sorted(extra.items())]


def edge_path_cells(
    edges: list[tuple[int, int, float]],
    sites: list[dict],
    predecessors_by_source: list[np.ndarray],
    shape: tuple[int, int],
) -> list[list[tuple[int, int]]]:
    """For each edge (i, j, hours), reconstruct the actual route as a list
    of (row, col) cells from site i to site j, reusing site i's already-
    computed predecessor array (no new Dijkstra runs needed -- the all-pairs
    step already paid for every shortest-path tree this network will ever
    need). Direction (i->j vs j->i) follows the edge's own (i, j) order,
    i.e. the DIRECTED path actually taken FROM i, which may differ in exact
    cell sequence (though not much in practice, at 120 m) from the reverse
    walk -- consistent with the directed-graph note in the module
    docstring: this is one concrete realization of the edge, not a claim
    that the reverse trip follows the identical cells. Works for MST edges,
    redundant edges, or any other (i, j, hours) list built from the same
    `sites`/`predecessors_by_source`/`shape`.
    """
    paths = []
    for i, j, _ in edges:
        si, sj = sites[i], sites[j]
        path = reconstruct_path(
            predecessors_by_source[i], shape, si["row"], si["col"], sj["row"], sj["col"]
        )
        if path is None:
            raise RuntimeError(
                f"Edge {si['name']}<->{sj['name']} came from a finite-cost "
                f"pairwise distance but reconstruct_path couldn't walk it back -- "
                f"this should not happen and indicates a real bug, not an "
                f"unreachable-site edge case (those are excluded before edge "
                f"construction, see run_tappa9_road_network.py)."
            )
        paths.append(path)
    return paths
