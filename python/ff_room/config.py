"""Configuration dataclasses for ff-room scenes.

All physical quantities in SI units: meters, m/s, kg/m³, Pa.
Angles in degrees (converted to unit vectors internally).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Tuple, List
import math
import json
import yaml


@dataclass
class RoomConfig:
    """Axis-aligned rectangular room."""
    size:    Tuple[float, float, float] = (5.0, 4.0, 2.5)  # Lx, Ly, Lz [m]
    grid:    Tuple[int, int, int]       = (25, 20, 12)      # Nx, Ny, Nz cells
    stretch: float                      = 1.0  # grid stretch ratio near walls
    # stretch=1.0: uniform (default). >1.0: geometric refinement near walls.
    # e.g. stretch=1.3 gives ~3x finer cells at walls vs. room center.


@dataclass
class FanConfig:
    """Single fan modeled as a rectangular inflow patch.

    position: center of fan face [m] in room coordinates
    direction: blowing direction as (azimuth_deg, elevation_deg)
               azimuth: 0=+x, 90=+y, 180=-x, 270=-y; elevation: 0=horizontal, 90=up
               Full 3-component oblique directions supported.
    velocity:  bulk inflow speed [m/s]
    radius:    half-width of fan face [m] (square patch: side = 2*radius)
    """
    position:  Tuple[float, float, float] = (0.5, 2.0, 1.0)
    direction: Tuple[float, float] = (0.0, 0.0)  # (azimuth_deg, elevation_deg)
    velocity:  float = 3.0
    radius:    float = 0.15
    name:      str   = "fan0"

    def direction_vector(self) -> Tuple[float, float, float]:
        """Return unit vector from azimuth/elevation angles."""
        az  = math.radians(self.direction[0])
        el  = math.radians(self.direction[1])
        x   = math.cos(el) * math.cos(az)
        y   = math.cos(el) * math.sin(az)
        z   = math.sin(el)
        return (x, y, z)

    def snapped_axis(self) -> int:
        """Return dominant axis index (0=x, 1=y, 2=z) for MVP axis-aligned BC."""
        v = self.direction_vector()
        return int(max(range(3), key=lambda i: abs(v[i])))

    def snapped_sign(self) -> float:
        """Return +1 or -1 for dominant axis direction."""
        v = self.direction_vector()
        ax = self.snapped_axis()
        return 1.0 if v[ax] >= 0 else -1.0


@dataclass
class ObstacleConfig:
    """Axis-aligned solid box obstacle (furniture, etc.)."""
    bbox_min: Tuple[float, float, float] = (2.0, 1.5, 0.0)
    bbox_max: Tuple[float, float, float] = (2.5, 2.5, 0.8)
    name:     str = "obstacle0"


@dataclass
class OpeningConfig:
    """Window or door: pressure-outflow BC + bidirectional temperature exchange.

    wall:   which room wall the opening is on
            "west"/"east"   (x=0 / x=Lx)
            "south"/"north" (y=0 / y=Ly)
            "floor"/"ceiling" (z=0 / z=Lz)
    center: (coord1, coord2) center of the opening on that wall [m]
            east/west walls:    (y_center, z_center)
            south/north walls:  (x_center, z_center)
            floor/ceiling:      (x_center, y_center)
    size:   (width, height) of the opening [m] on the wall plane
    T_outside: outside temperature [°C]; None = use solver.T_outside
    """
    wall:      str                      = "east"
    center:    Tuple[float, float]      = (2.0, 1.2)
    size:      Tuple[float, float]      = (0.8, 1.0)
    T_outside: float                    = None   # type: ignore[assignment]
    name:      str                      = "opening0"


@dataclass
class SolverConfig:
    """Numerical solver parameters.

    Each option is an independent toggle — combine freely:

      method:         "projection" (default) | "lbm"
                      projection = Chorin's projection method (MAC grid, CG/multigrid Poisson)
                      lbm        = D3Q19 BGK Lattice Boltzmann (isothermal, bounce-back BCs)
      tau:            LBM relaxation time; 0.5 < τ < 2.0 (method="lbm" only, default 0.8)

      advection:      "upwind" (default) | "quick" | "lax_wendroff"
                      upwind       = 1st-order upwind (fast, diffusive)
                      quick        = 3rd-order QUICK (accurate, same cost)
                      lax_wendroff = 2nd-order Lax-Wendroff (accurate, conditionally stable)
                      (projection method only; ignored for lbm)

      parallel:       False (default) | True
                      True = OpenMP multi-thread on advection/diffusion loops
                      (projection method only; silently ignored for lbm)

      poisson_method: "cg" (default) | "multigrid"
                      cg        = Conjugate Gradient O(N^1.5) — robust
                      multigrid = V-cycle Gauss-Seidel O(N)   — ~5x faster on large grids
                      (projection method only; ignored for lbm)

    Grid refinement (in RoomConfig):
      room.stretch:   1.0 (default, uniform) | >1.0 geometric wall refinement
                      e.g. 1.1 = mild refinement, 1.3 = strong refinement near all walls
    """
    dt:               float = 0.01    # time step [s]
    max_steps:        int   = 1000    # max iterations
    convergence_tol:  float = 1.0     # max |u_new - u_old| / dt
    rho:              float = 1.2     # air density [kg/m3]
    nu:               float = 1.5e-5  # kinematic viscosity [m2/s]
    poisson_max_iter: int   = 1000
    poisson_tol:      float = 1e-6
    # --- Solver selection ---
    method:         str   = "projection"  # "projection" | "lbm"
    tau:            float = 0.8           # LBM relaxation time (lbm only)
    # --- Projection-method sub-options ---
    advection:      str   = "upwind"      # "upwind" | "quick" | "lax_wendroff"
    parallel:       bool  = False         # True = OpenMP on
    poisson_method: str   = "cg"          # "cg" | "multigrid"
    # --- Thermal simulation ---
    T_initial:    float = 20.0        # initial indoor temperature [°C]
    T_outside:    float = 20.0        # outdoor temperature [°C]
    T_target:     float = None        # type: ignore[assignment]  stop at this mean T [°C]
    buoyancy:     bool  = False       # Boussinesq buoyancy (hot air rises)
    k_thermal:    float = 0.026       # thermal conductivity of air [W/m/K]
    cp:           float = 1005.0      # specific heat of air [J/kg/K]


@dataclass
class SceneConfig:
    """Complete reproducible scene specification.

    Serializable to/from JSON and YAML for full reproducibility.
    """
    room:      RoomConfig    = field(default_factory=RoomConfig)
    fans:      List[FanConfig]      = field(default_factory=list)
    obstacles: List[ObstacleConfig] = field(default_factory=list)
    openings:  List[OpeningConfig]  = field(default_factory=list)
    solver:    SolverConfig  = field(default_factory=SolverConfig)
    metadata:  dict          = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Serialization helpers
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        import dataclasses
        def _cvt(obj):
            if dataclasses.is_dataclass(obj):
                return {k: _cvt(v) for k, v in dataclasses.asdict(obj).items()}
            if isinstance(obj, (list, tuple)):
                return [_cvt(x) for x in obj]
            return obj
        return _cvt(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SceneConfig":
        def _cvt(v):
            return tuple(v) if isinstance(v, list) else v

        room = RoomConfig(**{k: _cvt(v) for k, v in d.get("room", {}).items()})
        fans = [FanConfig(**{k: _cvt(v) for k, v in f.items()})
                for f in d.get("fans", [])]
        obstacles = [ObstacleConfig(**{k: _cvt(v) for k, v in o.items()})
                     for o in d.get("obstacles", [])]
        openings  = [OpeningConfig(**{k: _cvt(v) for k, v in op.items()})
                     for op in d.get("openings", [])]
        solver = SolverConfig(**d.get("solver", {}))
        return cls(room=room, fans=fans, obstacles=obstacles, openings=openings,
                   solver=solver, metadata=d.get("metadata", {}))

    def save_yaml(self, path: str) -> None:
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    def save_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load_yaml(cls, path: str) -> "SceneConfig":
        with open(path) as f:
            return cls.from_dict(yaml.safe_load(f))

    @classmethod
    def load_json(cls, path: str) -> "SceneConfig":
        with open(path) as f:
            return cls.from_dict(json.load(f))
