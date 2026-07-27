"""
Vectorized droplet-based hydraulic erosion (rainfall/river carving), after
Hans Theodor Beyer's "Implementation of a Fast Method for Hydraulic Erosion
Simulation Using CUDA" (2015) -- the same droplet algorithm widely used in
procedural terrain tools (e.g. Sebastian Lague's implementation, which this
follows closely for parameter names/defaults).

This is FLUVIAL erosion only: rain falls, flows downhill under gravity,
picks up sediment where it can carry more, deposits where it can't. It has
nothing to do with wave/tidal/coastal erosion (a different physical
process this project doesn't model -- see the Tappa 1 decision doc for
that discussion).

Vectorization strategy: a true droplet simulation is inherently sequential
STEP BY STEP for a single droplet (each step depends on the last), but
droplets don't depend on each other within a step. So instead of a Python
loop over droplets (slow: hundreds of thousands of iterations), this runs
a Python loop over STEPS only (~30-70 iterations) with every droplet's
state held as one numpy array and advanced together. Erosion/deposition
writes to the shared height grid via np.add.at, which correctly
accumulates when multiple droplets touch the same cell in the same step
(not a race condition on the CPU/numpy side -- add.at is defined to
accumulate for repeated indices).
"""

import numpy as np


def _bilinear_gather(grid, gx, gy):
    """Sample `grid` (ny, nx) at fractional coordinates (gx, gy) -- gx along
    columns, gy along rows -- via bilinear interpolation. Returns the
    sampled values and the interpolation weights/corner indices, since the
    caller (erosion) needs the same corners to scatter deposits/erosion
    back with matching weights."""
    ny, nx = grid.shape
    x0 = np.floor(gx).astype(np.int64)
    y0 = np.floor(gy).astype(np.int64)
    x0 = np.clip(x0, 0, nx - 2)
    y0 = np.clip(y0, 0, ny - 2)
    x1 = x0 + 1
    y1 = y0 + 1
    fx = gx - x0
    fy = gy - y0

    w00 = (1 - fx) * (1 - fy)
    w10 = fx * (1 - fy)
    w01 = (1 - fx) * fy
    w11 = fx * fy

    values = (
        grid[y0, x0] * w00 + grid[y0, x1] * w10 +
        grid[y1, x0] * w01 + grid[y1, x1] * w11
    )
    corners = (y0, x0, y1, x1)
    weights = (w00, w10, w01, w11)
    return values, corners, weights


def _scatter_add(grid, corners, weights, amount):
    y0, x0, y1, x1 = corners
    w00, w10, w01, w11 = weights
    np.add.at(grid, (y0, x0), amount * w00)
    np.add.at(grid, (y0, x1), amount * w10)
    np.add.at(grid, (y1, x0), amount * w01)
    np.add.at(grid, (y1, x1), amount * w11)


def erode(
    elevation: np.ndarray,
    cell_size_m: float,
    n_droplets: int,
    seed: int,
    max_steps: int = 64,
    inertia: float = 0.05,
    sediment_capacity_factor: float = 4.0,
    min_sediment_capacity: float = 0.01,
    erode_speed: float = 0.3,
    deposit_speed: float = 0.3,
    evaporate_speed: float = 0.02,
    gravity: float = 4.0,
    initial_water: float = 1.0,
    initial_speed: float = 1.0,
    spawn_above_m: float = 0.0,
    height_scale_m: float = None,
    max_erosion_per_step_m: float = 8.0,
    verbose: bool = True,
):
    """Returns a NEW eroded elevation array (float32); `elevation` is not
    modified in place. Droplets spawn only on land (elevation > spawn_above_m,
    default sea level) -- rain doesn't erode a landmass from droplets that
    spawned over open ocean, and it wastes compute given roughly half this
    domain is underwater. A random land cell is picked per droplet (with
    sub-cell jitter), rather than uniform-over-the-whole-domain sampling
    then rejecting.

    height_scale_m: the reference Beyer/Lague algorithm (and the default
    constants above, taken from it) was tuned against heightmaps normalized
    to roughly [0, 1]. This terrain's elevation is in raw metres up to
    ~4000 -- feeding that directly into the same constants (gravity=4,
    sediment_capacity_factor=4, etc) inflates delta_height/capacity/erosion
    amounts by orders of magnitude relative to what those constants assume,
    and produces a genuine runaway feedback loop (confirmed by reproducing
    it: erosion deepens a pit, which increases next-step delta_height for
    whichever droplet visits it next, which increases erosion further,
    without bound, until the grid overflows to inf). Dividing by
    height_scale_m before simulating and multiplying back after fixes the
    units mismatch; defaults to this DEM's own elevation range if not given.
    max_erosion_per_step_m is a second, independent safety cap (defense in
    depth, not a substitute for the units fix) on how much any single
    droplet can remove/deposit in one step, in real metres.
    """
    rng = np.random.RandomState(seed)
    ny, nx = elevation.shape
    if height_scale_m is None:
        height_scale_m = max(float(np.nanmax(elevation) - np.nanmin(elevation)), 1.0)
    grid = elevation.astype(np.float64).copy() / height_scale_m
    max_erosion_per_step = max_erosion_per_step_m / height_scale_m

    land_rows, land_cols = np.where(elevation[1:ny-1, 1:nx-1] > spawn_above_m)
    land_rows = land_rows + 1  # undo the border trim in the where() call
    land_cols = land_cols + 1
    if len(land_rows) == 0:
        raise ValueError("no land cells above spawn_above_m -- nothing to erode")
    chosen = rng.randint(0, len(land_rows), size=n_droplets)
    pos_x = land_cols[chosen].astype(np.float64) + rng.uniform(-0.5, 0.5, size=n_droplets)
    pos_y = land_rows[chosen].astype(np.float64) + rng.uniform(-0.5, 0.5, size=n_droplets)
    pos_x = np.clip(pos_x, 1.0, nx - 2.0)
    pos_y = np.clip(pos_y, 1.0, ny - 2.0)
    dir_x = np.zeros(n_droplets)
    dir_y = np.zeros(n_droplets)
    speed = np.full(n_droplets, initial_speed, dtype=np.float64)
    water = np.full(n_droplets, initial_water, dtype=np.float64)
    sediment = np.zeros(n_droplets)
    alive = np.ones(n_droplets, dtype=bool)

    for step in range(max_steps):
        if not alive.any():
            break

        height, corners, weights = _bilinear_gather(grid, pos_x, pos_y)

        # gradient via central-ish finite difference from the same corners
        # (cheap re-use: gradient x ~ (right samples - left samples))
        y0, x0, y1, x1 = corners
        w00, w10, w01, w11 = weights
        grad_x = (grid[y0, x1] - grid[y0, x0]) * (1 - (pos_y - y0)) + (grid[y1, x1] - grid[y1, x0]) * (pos_y - y0)
        grad_y = (grid[y1, x0] - grid[y0, x0]) * (1 - (pos_x - x0)) + (grid[y1, x1] - grid[y0, x1]) * (pos_x - x0)

        dir_x = dir_x * inertia - grad_x * (1 - inertia)
        dir_y = dir_y * inertia - grad_y * (1 - inertia)
        norm = np.sqrt(dir_x ** 2 + dir_y ** 2)
        norm = np.where(norm < 1e-8, 1.0, norm)
        dir_x /= norm
        dir_y /= norm

        new_x = pos_x + dir_x
        new_y = pos_y + dir_y

        out_of_bounds = (new_x < 1.0) | (new_x > nx - 2.0) | (new_y < 1.0) | (new_y > ny - 2.0)
        alive &= ~out_of_bounds
        new_x = np.clip(new_x, 1.0, nx - 2.0)
        new_y = np.clip(new_y, 1.0, ny - 2.0)

        new_height, new_corners, new_weights = _bilinear_gather(grid, new_x, new_y)
        delta_height = new_height - height
        # Zero out dead droplets' delta_height HERE, before it feeds capacity/
        # speed -- not just at the final scatter. A droplet clipped onto the
        # boundary can sit at a cell with a large height difference every
        # remaining step; left unmasked, its capacity/speed can grow to inf
        # over many steps, and inf * active(=0) is NaN, not zero -- that NaN
        # then gets written into the grid and poisons every neighboring
        # bilinear sample from then on. Confirmed this by reproducing the
        # blow-up and tracing it to exactly this path.
        delta_height = np.where(alive, delta_height, 0.0)

        capacity = np.maximum(-delta_height * speed * water * sediment_capacity_factor, min_sediment_capacity)

        depositing = (delta_height > 0) | (sediment > capacity)
        deposit_amount = np.where(
            delta_height > 0,
            np.minimum(delta_height, sediment),
            (sediment - capacity) * deposit_speed,
        )
        deposit_amount = np.where(depositing, np.maximum(deposit_amount, 0.0), 0.0)

        erode_amount = np.where(~depositing, np.minimum((capacity - sediment) * erode_speed, -delta_height), 0.0)
        erode_amount = np.maximum(erode_amount, 0.0)

        # Hard safety cap, independent of the height_scale_m fix above --
        # see the erode() docstring for why both exist.
        deposit_amount = np.minimum(deposit_amount, max_erosion_per_step)
        erode_amount = np.minimum(erode_amount, max_erosion_per_step)

        active = alive.astype(np.float64)
        _scatter_add(grid, corners, weights, deposit_amount * active)
        _scatter_add(grid, corners, weights, -erode_amount * active)

        sediment = sediment + (erode_amount - deposit_amount) * active
        speed = np.sqrt(np.maximum(speed ** 2 + (-delta_height) * gravity, 0.0))
        speed = np.clip(speed, 0.0, 100.0)  # defense in depth against any other runaway path
        water = water * (1 - evaporate_speed)

        pos_x, pos_y = new_x, new_y
        alive &= water > 1e-3

        if verbose and (step % 10 == 0 or step == max_steps - 1):
            print(f"  erosion step {step+1}/{max_steps}: {alive.sum():,} droplets still active")

    return (grid * height_scale_m).astype(np.float32)
