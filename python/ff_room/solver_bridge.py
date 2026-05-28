"""Python wrapper around C++ FluidSolver."""
from __future__ import annotations
import time
from dataclasses import dataclass, field
from typing import List, Optional, Callable

import numpy as np

from .config import SceneConfig
from .scene import Scene, _load_core


@dataclass
class FieldSnapshot:
    """Full field snapshot captured during a simulation run.

    step:    solver step index
    time_s:  physical time = step * dt [s]
    T_mean:  mean indoor temperature [°C]
    speed:   (Nx, Ny, Nz) velocity magnitude [m/s]
    T_field: (Nx, Ny, Nz) temperature [°C]; zeros if thermal disabled
    vel:     (Nx, Ny, Nz, 3) cell-centered velocity vectors [m/s]
    """
    step:    int
    time_s:  float
    T_mean:  float
    speed:   np.ndarray
    T_field: np.ndarray
    vel:     np.ndarray


@dataclass
class SimResult:
    """Velocity, pressure, and (optionally) temperature fields after simulation.

    Shaped arrays:
      u_field:   (Nx+1, Ny, Nz)
      v_field:   (Nx,   Ny+1, Nz)
      w_field:   (Nx,   Ny,   Nz+1)
      p_field:   (Nx,   Ny,   Nz)
      T_field:   (Nx,   Ny,   Nz)  temperature [°C]; all T_initial if thermal disabled
      cell_type: (Nx,   Ny,   Nz)  uint8
    """
    u_field:   np.ndarray
    v_field:   np.ndarray
    w_field:   np.ndarray
    p_field:   np.ndarray
    T_field:   np.ndarray
    cell_type: np.ndarray
    config:    SceneConfig
    # Diagnostics
    steps:           int   = 0
    converged:       bool  = False
    velocity_change: float = 0.0
    divergence_max:  float = 0.0
    T_mean:          float = 0.0
    wall_time_s:     float = 0.0

    def velocity_magnitude(self) -> np.ndarray:
        """Cell-centered velocity magnitude (Nx, Ny, Nz)."""
        uc = 0.5 * (self.u_field[:-1, :, :] + self.u_field[1:, :, :])
        vc = 0.5 * (self.v_field[:, :-1, :] + self.v_field[:, 1:, :])
        wc = 0.5 * (self.w_field[:, :, :-1] + self.w_field[:, :, 1:])
        return np.sqrt(uc**2 + vc**2 + wc**2)

    def velocity_cell_centered(self) -> np.ndarray:
        """Returns (Nx, Ny, Nz, 3) cell-centered velocity vector field."""
        uc = 0.5 * (self.u_field[:-1, :, :] + self.u_field[1:, :, :])
        vc = 0.5 * (self.v_field[:, :-1, :] + self.v_field[:, 1:, :])
        wc = 0.5 * (self.w_field[:, :, :-1] + self.w_field[:, :, 1:])
        return np.stack([uc, vc, wc], axis=-1)


def _apply_solver_flags(sc, params, c) -> None:
    """Map SolverConfig advection/parallel onto C++ FluidSolverParams."""
    adv = (sc.advection or 'upwind').lower()
    if adv == 'quick':
        params.advection = c.AdvectionScheme.QUICK
    elif adv == 'lax_wendroff':
        params.advection = c.AdvectionScheme.LAX_WENDROFF
    else:
        params.advection = c.AdvectionScheme.UPWIND1
    params.use_openmp = bool(sc.parallel)


class SolverBridge:
    """Run C++ FluidSolver from Python, return SimResult."""

    def __init__(self, config: SceneConfig):
        self.config = config
        self._c = _load_core()

    def run(
        self,
        progress_callback: Optional[Callable[[int, float, float], None]] = None,
        print_interval: int = 50,
    ) -> SimResult:
        """Build scene, run solver to convergence, return SimResult.

        progress_callback(step, vel_change, div_max, T_mean) called each step if set.
        print_interval: print diagnostics every N steps (0=silent).
        """
        scene = Scene(self.config)
        g     = scene.grid
        bm    = scene.bm
        sc    = self.config.solver

        params = self._c.FluidSolverParams()
        params.rho             = sc.rho
        params.nu              = sc.nu
        params.dt              = sc.dt
        params.max_steps       = sc.max_steps
        params.convergence_tol = sc.convergence_tol
        params.poisson.max_iter = sc.poisson_max_iter
        params.poisson.tol      = sc.poisson_tol
        # Thermal
        has_openings = bool(self.config.openings)
        params.thermal    = has_openings or sc.buoyancy or (sc.T_target is not None)
        params.T_initial  = sc.T_initial
        params.T_target   = sc.T_target if sc.T_target is not None else 1e30
        params.buoyancy   = sc.buoyancy
        params.k_thermal  = sc.k_thermal
        params.cp         = sc.cp
        params.beta       = 1.0 / (sc.T_initial + 273.15)

        _apply_solver_flags(sc, params, self._c)
        params.use_multigrid = (getattr(sc, 'poisson_method', 'cg') == 'multigrid')

        if sc.method == 'lbm':
            solver = self._c.LBMSolver(g, bm, sc.tau, sc.max_steps, sc.convergence_tol, sc.dt)
        else:  # default: projection
            solver = self._c.FluidSolver(g, bm, params)

        diagnostics = []
        t0 = time.perf_counter()

        def _cb(res):
            diagnostics.append(res)
            if progress_callback:
                progress_callback(res.step, res.velocity_change, res.divergence_max,
                                  res.T_mean)
            if print_interval > 0 and res.step % print_interval == 0:
                line = (f"  step={res.step:4d}  vel_change={res.velocity_change:.3e}"
                        f"  p_res={res.pressure_residual:.3e}"
                        f"  div_max={res.divergence_max:.3e}")
                if sc.method != 'lbm' and params.thermal:
                    line += f"  T_mean={res.T_mean:.2f}°C"
                line += "  [CONVERGED]" if res.converged else ""
                print(line)

        final = solver.run(_cb)
        wall_time = time.perf_counter() - t0

        # Copy fields out (numpy copies so grid can be freed)
        result = SimResult(
            u_field   = np.array(g.u_shaped()),
            v_field   = np.array(g.v_shaped()),
            w_field   = np.array(g.w_shaped()),
            p_field   = np.array(g.p_shaped()),
            T_field   = np.array(g.T_shaped()),
            cell_type = g.cell_type_shaped(),
            config    = self.config,
            steps           = final.step,
            converged       = final.converged,
            velocity_change = final.velocity_change,
            divergence_max  = final.divergence_max,
            T_mean          = final.T_mean,
            wall_time_s     = wall_time,
        )
        return result

    def run_animated(
        self,
        save_every: int = 50,
        progress_callback: Optional[Callable[[int, float, float], None]] = None,
        print_interval: int = 50,
    ):
        """Run solver step-by-step, capturing field snapshots every `save_every` steps.

        Returns (SimResult, list[FieldSnapshot]).

        save_every:        capture a full snapshot every N steps (step 0 is always captured)
        progress_callback: same signature as run() — called every step
        print_interval:    print diagnostics every N steps (0=silent)
        """
        scene = Scene(self.config)
        g     = scene.grid
        bm    = scene.bm
        sc    = self.config.solver

        params = self._c.FluidSolverParams()
        params.rho             = sc.rho
        params.nu              = sc.nu
        params.dt              = sc.dt
        params.max_steps       = sc.max_steps
        params.convergence_tol = sc.convergence_tol
        params.poisson.max_iter = sc.poisson_max_iter
        params.poisson.tol      = sc.poisson_tol
        # Thermal
        has_openings = bool(self.config.openings)
        params.thermal    = has_openings or sc.buoyancy or (sc.T_target is not None)
        params.T_initial  = sc.T_initial
        params.T_target   = sc.T_target if sc.T_target is not None else 1e30
        params.buoyancy   = sc.buoyancy
        params.k_thermal  = sc.k_thermal
        params.cp         = sc.cp
        params.beta       = 1.0 / (sc.T_initial + 273.15)

        _apply_solver_flags(sc, params, self._c)
        params.use_multigrid = (getattr(sc, 'poisson_method', 'cg') == 'multigrid')

        if sc.method == 'lbm':
            solver = self._c.LBMSolver(g, bm, sc.tau, sc.max_steps, sc.convergence_tol, sc.dt)
        else:  # default: projection
            solver = self._c.FluidSolver(g, bm, params)

        diagnostics: list = []
        snapshots: List[FieldSnapshot] = []
        t0 = time.perf_counter()
        dt = sc.dt

        def _make_snapshot(step: int, T_mean: float) -> FieldSnapshot:
            u = np.array(g.u_shaped())
            v = np.array(g.v_shaped())
            w = np.array(g.w_shaped())
            T = np.array(g.T_shaped()) if (sc.method != 'lbm' and params.thermal) else np.zeros(
                (u.shape[0]-1, u.shape[1], u.shape[2]))
            # Cell-center averages
            uc = 0.5 * (u[:-1, :, :] + u[1:, :, :])
            vc = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
            wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
            speed = np.sqrt(uc**2 + vc**2 + wc**2)
            vel   = np.stack([uc, vc, wc], axis=-1)
            return FieldSnapshot(
                step    = step,
                time_s  = step * dt,
                T_mean  = T_mean,
                speed   = speed,
                T_field = T,
                vel     = vel,
            )

        def _cb(res):
            diagnostics.append(res)
            if res.step % save_every == 0:
                snapshots.append(_make_snapshot(res.step, res.T_mean))
            if progress_callback:
                progress_callback(res.step, res.velocity_change, res.divergence_max,
                                  res.T_mean)
            if print_interval > 0 and res.step % print_interval == 0:
                line = (f"  step={res.step:4d}  vel_change={res.velocity_change:.3e}"
                        f"  p_res={res.pressure_residual:.3e}"
                        f"  div_max={res.divergence_max:.3e}")
                if sc.method != 'lbm' and params.thermal:
                    line += f"  T_mean={res.T_mean:.2f}°C"
                line += "  [CONVERGED]" if res.converged else ""
                print(line)

        final = solver.run(_cb)
        wall_time = time.perf_counter() - t0

        # Always capture a final snapshot if not already at the last step
        if not snapshots or snapshots[-1].step != final.step:
            snapshots.append(_make_snapshot(final.step, final.T_mean))

        result = SimResult(
            u_field   = np.array(g.u_shaped()),
            v_field   = np.array(g.v_shaped()),
            w_field   = np.array(g.w_shaped()),
            p_field   = np.array(g.p_shaped()),
            T_field   = np.array(g.T_shaped()),
            cell_type = g.cell_type_shaped(),
            config    = self.config,
            steps           = final.step,
            converged       = final.converged,
            velocity_change = final.velocity_change,
            divergence_max  = final.divergence_max,
            T_mean          = final.T_mean,
            wall_time_s     = wall_time,
        )
        return result, snapshots
