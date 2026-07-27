"""
Load the hand-authored terrain skeleton (ridges + zones) and turn it into
raster fields: a ridge elevation-influence field, a zone elevation-override
field (plateaus), and a noise-amplitude multiplier field (zones).

Distance calculations use scipy.spatial.cKDTree against densely resampled
vertices of each line/polygon boundary, rather than exact point-to-segment
distance. This is an approximation, not exact geometry -- see
`_densify_polyline`'s docstring for the error bound and why it's negligible
here. geopandas/shapely are not used because they aren't installable in
this sandbox (see noise.py's module docstring for the same network
constraint); point-in-polygon uses matplotlib.path.Path instead, which is
already a project dependency (matplotlib) and vectorizes over the whole
grid in one call.
"""

import json
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree
from matplotlib.path import Path


def load_geojson(path):
    with open(path) as f:
        return json.load(f)["features"]


def _densify_polyline(coords: np.ndarray, max_spacing_m: float) -> np.ndarray:
    """Resample a polyline (N,2) so consecutive points are at most
    `max_spacing_m` apart, by linear interpolation along each segment.

    Why this is an acceptable stand-in for true point-to-segment distance:
    nearest-vertex distance from a resampled polyline overstates the true
    distance to the nearest point on the original segment by at most
    (max_spacing_m / 2) in the worst case (a query point exactly opposite
    the midpoint of a long straight segment). With max_spacing_m set well
    below the smallest falloff_km in play (falloffs here are 10-23 km;
    using ~250 m spacing), the worst-case error is under 1% of any
    falloff_km distance -- negligible next to the visual effect of the
    decay curve itself.
    """
    out = [coords[0]]
    for i in range(len(coords) - 1):
        p0, p1 = coords[i], coords[i + 1]
        seg_len = np.hypot(*(p1 - p0))
        if seg_len <= max_spacing_m:
            out.append(p1)
            continue
        n_steps = int(np.ceil(seg_len / max_spacing_m))
        for s in range(1, n_steps + 1):
            out.append(p0 + (p1 - p0) * (s / n_steps))
    return np.array(out)


@dataclass
class RidgeField:
    peak_elevation_m: float
    falloff_km: float
    name: str
    tree: cKDTree
    shelf_multiplier: float = 3.0  # shelf_km = falloff_km * shelf_multiplier

    def contribution(self, xy: np.ndarray) -> np.ndarray:
        """Gaussian-style decay near the crest -- falloff_km is the
        HALF-MAX distance, i.e. contribution == peak/2 exactly at
        d == falloff_km*1000 metres (k = ln(2) makes that literal) --
        windowed to exactly zero beyond shelf_km = falloff_km *
        shelf_multiplier. Without this window the Gaussian tail never
        truly reaches zero, so a single global sea-level offset would be
        the ONLY thing separating land from ocean everywhere, and every
        ridge would implicitly claim the same *relative* shelf regardless
        of how far its narrative role should reach. shelf_multiplier is
        the per-ridge knob for that: a wide multiplier gives a ridge a
        generous coastal shelf (its land keeps going well past the
        falloff-defined crest zone); a tight multiplier makes its
        coastline hug close to the crest itself.
        """
        dist_m, _ = self.tree.query(xy, k=1)
        k = np.log(2.0)
        falloff_m = self.falloff_km * 1000.0
        shelf_m = falloff_m * self.shelf_multiplier
        base = self.peak_elevation_m * np.exp(-k * (dist_m / falloff_m) ** 2)
        t = np.clip((shelf_m - dist_m) / (shelf_m - falloff_m), 0.0, 1.0)
        taper = t * t * t * (t * (t * 6 - 15) + 10)  # smootherstep
        taper = np.where(dist_m <= falloff_m, 1.0, taper)
        return base * taper


@dataclass
class ZoneField:
    feature_type: str  # "plateau" or "amplitude_zone"
    target_elevation_m: float  # None for amplitude_zone
    amplitude_scale: float
    edge_transition_km: float
    name: str
    path: Path
    boundary_tree: cKDTree
    base_lift_m: float = 0.0  # amplitude_zone only, see generate.py

    def blend_weight(self, xy: np.ndarray) -> np.ndarray:
        """Smootherstep weight in [0, 1]: 1.0 deep inside the polygon,
        0.0 far outside, transitioning across a band of width
        2*edge_transition_km straddling the drawn boundary (i.e.
        edge_transition_km on the inside AND edge_transition_km on the
        outside -- documented explicitly since the source attribute
        schema doesn't say which)."""
        inside = self.path.contains_points(xy)
        dist_m, _ = self.boundary_tree.query(xy, k=1)
        signed = np.where(inside, dist_m, -dist_m)
        band_m = self.edge_transition_km * 1000.0
        t = np.clip((signed + band_m) / (2 * band_m), 0.0, 1.0)
        # smootherstep (Ken Perlin's improved smoothstep: zero 1st & 2nd
        # derivative at both ends, avoids a visible seam at the band edge)
        return t * t * t * (t * (t * 6 - 15) + 10)


def build_ridge_fields(ridge_features, densify_spacing_m: float = 250.0, shelf_multipliers: dict = None, default_shelf_multiplier: float = 3.0):
    """`shelf_multipliers`: optional {ridge name: multiplier} overrides.
    Any ridge not listed falls back to `default_shelf_multiplier`."""
    shelf_multipliers = shelf_multipliers or {}
    fields = []
    for feat in ridge_features:
        props = feat["properties"]
        if props.get("feature_type") != "ridge":
            raise ValueError(f"unexpected feature_type in ridges layer: {props}")
        coords = np.array(feat["geometry"]["coordinates"], dtype=np.float64)
        dense = _densify_polyline(coords, densify_spacing_m)
        name = props.get("name", "unnamed ridge")
        fields.append(
            RidgeField(
                peak_elevation_m=props["peak_elevation_m"],
                falloff_km=props["falloff_km"],
                name=name,
                tree=cKDTree(dense),
                shelf_multiplier=shelf_multipliers.get(name, default_shelf_multiplier),
            )
        )
    return fields


def build_zone_fields(zone_features, densify_spacing_m: float = 250.0, base_lift_m: dict = None, default_base_lift_m: float = 0.0):
    """`base_lift_m`: optional {zone name: metres} overrides for
    amplitude_zone features (ignored for plateau, which already gets its
    own target_elevation_m). See generate.py's docstring for why this
    exists -- an amplitude_zone drawn far from every ridge has nothing to
    "smooth", since the spec's own definition (low roughness, broader
    trend continues underneath) assumes a trend is actually present."""
    base_lift_m = base_lift_m or {}
    fields = []
    for feat in zone_features:
        props = feat["properties"]
        ftype = props.get("feature_type")
        if ftype not in ("plateau", "amplitude_zone"):
            raise ValueError(f"unexpected feature_type in zones layer: {props}")
        ring = np.array(feat["geometry"]["coordinates"][0], dtype=np.float64)
        dense = _densify_polyline(ring, densify_spacing_m)
        name = props.get("name", "unnamed zone")
        fields.append(
            ZoneField(
                feature_type=ftype,
                target_elevation_m=props.get("target_elevation_m"),
                amplitude_scale=props["amplitude_scale"],
                edge_transition_km=props["edge_transition_km"],
                name=name,
                path=Path(ring),
                boundary_tree=cKDTree(dense),
                base_lift_m=base_lift_m.get(name, default_base_lift_m) if ftype == "amplitude_zone" else 0.0,
            )
        )
    return fields
