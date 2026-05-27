"""Scene: translates SceneConfig into C++ Grid + BoundaryManager objects."""
from __future__ import annotations
import math
from typing import TYPE_CHECKING

from .config import SceneConfig, FanConfig, ObstacleConfig, OpeningConfig

if TYPE_CHECKING:
    pass


def _load_core():
    """Lazy import of compiled C++ module."""
    try:
        import _ffroom_core as _c
        return _c
    except ImportError:
        import sys
        import pathlib
        # The .so lives next to this file when installed via build.sh
        pkg_dir = str(pathlib.Path(__file__).parent)
        if pkg_dir not in sys.path:
            sys.path.insert(0, pkg_dir)
        try:
            import _ffroom_core as _c
            return _c
        except ImportError as e:
            raise ImportError(
                "C++ extension _ffroom_core not found. "
                "Run: bash build.sh"
            ) from e


class Scene:
    """Assembled simulation scene backed by C++ Grid and BoundaryManager."""

    def __init__(self, config: SceneConfig):
        self.config = config
        self._c = _load_core()
        self.grid = None
        self.bm   = None
        self._build()

    def _build(self) -> None:
        cfg = self.config
        r   = cfg.room

        self.grid = self._c.Grid(
            r.grid[0], r.grid[1], r.grid[2],
            r.size[0], r.size[1], r.size[2],
        )
        self.bm = self._c.BoundaryManager()

        for obs in cfg.obstacles:
            self._add_obstacle(obs)

        for fan in cfg.fans:
            self._add_fan(fan)

        for op in cfg.openings:
            self._add_opening(op)

        # Auto outflow: mark the full wall opposite the first fan as OUTFLOW.
        # This provides Dirichlet p=0 for the Poisson solver.
        # Skipped when user-defined openings already provide pressure BCs,
        # because auto-outflow cells would inject air at T_initial (room temp)
        # instead of T_outside, preventing correct thermal cooling simulation.
        if cfg.fans and not cfg.openings:
            self._add_auto_outflow(cfg.fans[0])

    def _add_auto_outflow(self, fan: FanConfig) -> None:
        """Mark the entire opposite wall as OUTFLOW.

        Physically: open face on the wall the fan blows toward (like a window or vent).
        Using the full wall gives enough Dirichlet constraints for the Poisson solver
        to be well-conditioned.
        """
        r  = self.config.room
        ax = fan.snapped_axis()
        sg = fan.snapped_sign()

        obc = self._c.OutflowBC()
        obc.axis = ax
        obc.side = 1 if sg > 0 else 0  # wall that fan blows toward

        if ax == 0:
            obc.a_min, obc.a_max = 0, r.grid[1] - 1
            obc.b_min, obc.b_max = 0, r.grid[2] - 1
        elif ax == 1:
            obc.a_min, obc.a_max = 0, r.grid[0] - 1
            obc.b_min, obc.b_max = 0, r.grid[2] - 1
        else:
            obc.a_min, obc.a_max = 0, r.grid[0] - 1
            obc.b_min, obc.b_max = 0, r.grid[1] - 1

        self.bm.add_outflow(self.grid, obc)

    def _add_obstacle(self, obs: ObstacleConfig) -> None:
        g = self.grid
        r = self.config.room
        def to_cell(x, L, N): return int(x / L * N)
        i0 = to_cell(obs.bbox_min[0], r.size[0], r.grid[0])
        i1 = to_cell(obs.bbox_max[0], r.size[0], r.grid[0]) - 1
        j0 = to_cell(obs.bbox_min[1], r.size[1], r.grid[1])
        j1 = to_cell(obs.bbox_max[1], r.size[1], r.grid[1]) - 1
        k0 = to_cell(obs.bbox_min[2], r.size[2], r.grid[2])
        k1 = to_cell(obs.bbox_max[2], r.size[2], r.grid[2]) - 1
        self.bm.add_solid_box(self.grid, i0, i1, j0, j1, k0, k1)

    def _add_fan(self, fan: FanConfig) -> None:
        g  = self.grid
        r  = self.config.room
        ax = fan.snapped_axis()
        sg = fan.snapped_sign()

        def to_cell(x, L, N):   return int(x / L * N)
        def to_cell_r(x, r, L, N):  # range from center ± radius
            lo = max(0, int((x - r) / L * N))
            hi = min(N-1, int((x + r) / L * N))
            return lo, hi

        cx, cy, cz = fan.position
        rad = fan.radius

        bc = self._c.FanBC()

        if ax == 0:  # blowing in x
            ix = to_cell(cx, r.size[0], r.grid[0])
            ix = max(0, min(r.grid[0]-1, ix))
            jlo, jhi = to_cell_r(cy, rad, r.size[1], r.grid[1])
            klo, khi = to_cell_r(cz, rad, r.size[2], r.grid[2])
            bc.i_min, bc.i_max = ix, ix
            bc.j_min, bc.j_max = jlo, jhi
            bc.k_min, bc.k_max = klo, khi
            bc.vel = [sg * fan.velocity, 0.0, 0.0]
        elif ax == 1:  # blowing in y
            ilo, ihi = to_cell_r(cx, rad, r.size[0], r.grid[0])
            jy = to_cell(cy, r.size[1], r.grid[1])
            jy = max(0, min(r.grid[1]-1, jy))
            klo, khi = to_cell_r(cz, rad, r.size[2], r.grid[2])
            bc.i_min, bc.i_max = ilo, ihi
            bc.j_min, bc.j_max = jy, jy
            bc.k_min, bc.k_max = klo, khi
            bc.vel = [0.0, sg * fan.velocity, 0.0]
        else:  # blowing in z
            ilo, ihi = to_cell_r(cx, rad, r.size[0], r.grid[0])
            jlo, jhi = to_cell_r(cy, rad, r.size[1], r.grid[1])
            kz = to_cell(cz, r.size[2], r.grid[2])
            kz = max(0, min(r.grid[2]-1, kz))
            bc.i_min, bc.i_max = ilo, ihi
            bc.j_min, bc.j_max = jlo, jhi
            bc.k_min, bc.k_max = kz, kz
            bc.vel = [0.0, 0.0, sg * fan.velocity]

        self.bm.add_fan(self.grid, bc)

    def _add_opening(self, op: OpeningConfig) -> None:
        r   = self.config.room
        sc  = self.config.solver

        # Wall name → axis, side, and which room dimensions are transverse
        wall_map = {
            "west":    (0, 0, r.size[1], r.grid[1], r.size[2], r.grid[2]),
            "east":    (0, 1, r.size[1], r.grid[1], r.size[2], r.grid[2]),
            "south":   (1, 0, r.size[0], r.grid[0], r.size[2], r.grid[2]),
            "north":   (1, 1, r.size[0], r.grid[0], r.size[2], r.grid[2]),
            "floor":   (2, 0, r.size[0], r.grid[0], r.size[1], r.grid[1]),
            "ceiling": (2, 1, r.size[0], r.grid[0], r.size[1], r.grid[1]),
        }
        ax, sd, La, Na, Lb, Nb = wall_map[op.wall.lower()]

        ca, cb = op.center
        wa, wb = op.size

        def cell_range(c, w, L, N):
            lo = max(0, int((c - w/2) / L * N))
            hi = min(N - 1, int((c + w/2) / L * N))
            return lo, hi

        a_lo, a_hi = cell_range(ca, wa, La, Na)
        b_lo, b_hi = cell_range(cb, wb, Lb, Nb)

        T_out = op.T_outside if op.T_outside is not None else sc.T_outside

        obc = self._c.OpeningBC()
        obc.axis      = ax
        obc.side      = sd
        obc.a_min     = a_lo
        obc.a_max     = a_hi
        obc.b_min     = b_lo
        obc.b_max     = b_hi
        obc.T_outside = T_out

        self.bm.add_opening(self.grid, obc)

    def rebuild(self) -> None:
        """Re-create grid and BCs from config (after config change)."""
        self._build()
