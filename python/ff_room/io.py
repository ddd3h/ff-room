"""Result I/O: save and load simulation results.

Formats:
  .npz  — NumPy compressed (velocity, pressure, cell_type, metadata)
  .vts  — VTK structured grid via PyVista (for paraview / re-visualization)
  .jsonl — experiment log (append-only)
"""
from __future__ import annotations
import json
import time
import hashlib
from pathlib import Path

import numpy as np

from .solver_bridge import SimResult
from .config import SceneConfig


class ResultStore:
    """Save and load SimResult objects."""

    @staticmethod
    def save_npz(result: SimResult, path: str) -> None:
        """Save fields + config as compressed NumPy archive."""
        cfg_dict = result.config.to_dict()
        np.savez_compressed(
            path,
            u_field   = result.u_field,
            v_field   = result.v_field,
            w_field   = result.w_field,
            p_field   = result.p_field,
            T_field   = result.T_field,
            cell_type = result.cell_type,
            _steps          = result.steps,
            _converged      = result.converged,
            _velocity_change = result.velocity_change,
            _divergence_max  = result.divergence_max,
            _T_mean          = result.T_mean,
            _wall_time_s     = result.wall_time_s,
            _config_json     = np.array(json.dumps(cfg_dict)),
        )

    @staticmethod
    def load_npz(path: str) -> SimResult:
        """Load SimResult from .npz archive."""
        data = np.load(path, allow_pickle=True)
        cfg  = SceneConfig.from_dict(json.loads(str(data["_config_json"])))
        T_field = data["T_field"] if "T_field" in data else np.zeros_like(data["p_field"])
        return SimResult(
            u_field   = data["u_field"],
            v_field   = data["v_field"],
            w_field   = data["w_field"],
            p_field   = data["p_field"],
            T_field   = T_field,
            cell_type = data["cell_type"],
            config    = cfg,
            steps           = int(data["_steps"]),
            converged       = bool(data["_converged"]),
            velocity_change = float(data["_velocity_change"]),
            divergence_max  = float(data["_divergence_max"]),
            T_mean          = float(data["_T_mean"]) if "_T_mean" in data else 0.0,
            wall_time_s     = float(data["_wall_time_s"]),
        )

    @staticmethod
    def save_vtk(result: SimResult, path: str) -> None:
        """Save as VTK StructuredGrid (.vts) via PyVista."""
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista required: pip install pyvista")

        from .visualization import Visualizer
        viz  = Visualizer(result)
        grid = viz._make_structured_grid()
        grid.save(path)

    @staticmethod
    def save_velocity_csv(result: SimResult, path: str,
                          component: str = "all") -> None:
        """Export cell-centered velocity to CSV (x,y,z,u,v,w,speed).

        Useful for post-processing in other tools.
        component: 'all' | 'u' | 'v' | 'w' | 'speed'
        """
        r   = result.config.room
        Nx, Ny, Nz = r.grid
        Lx, Ly, Lz = r.size
        vel = result.velocity_cell_centered()  # (Nx, Ny, Nz, 3)
        spd = result.velocity_magnitude()      # (Nx, Ny, Nz)

        xs = np.linspace(Lx/(2*Nx), Lx*(1 - 1/(2*Nx)), Nx)
        ys = np.linspace(Ly/(2*Ny), Ly*(1 - 1/(2*Ny)), Ny)
        zs = np.linspace(Lz/(2*Nz), Lz*(1 - 1/(2*Nz)), Nz)
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

        rows = np.column_stack([
            X.ravel(), Y.ravel(), Z.ravel(),
            vel[:,:,:,0].ravel(), vel[:,:,:,1].ravel(), vel[:,:,:,2].ravel(),
            spd.ravel()
        ])
        header = "x_m,y_m,z_m,u_ms,v_ms,w_ms,speed_ms"
        np.savetxt(path, rows, delimiter=",", header=header, comments="")


class ExperimentLog:
    """Append-only JSON Lines experiment log.

    Each entry records config hash + key metrics for reproducibility.
    """

    def __init__(self, path: str):
        self.path = Path(path)

    def append(self, result: SimResult, notes: str = "") -> str:
        """Append log entry. Returns config hash (8 chars)."""
        cfg_str   = json.dumps(result.config.to_dict(), sort_keys=True)
        cfg_hash  = hashlib.sha256(cfg_str.encode()).hexdigest()[:8]
        entry = {
            "timestamp":      time.strftime("%Y-%m-%dT%H:%M:%S"),
            "config_hash":    cfg_hash,
            "steps":          result.steps,
            "converged":      result.converged,
            "velocity_change": result.velocity_change,
            "divergence_max": result.divergence_max,
            "T_mean":         result.T_mean,
            "wall_time_s":    result.wall_time_s,
            "notes":          notes,
            "config":         result.config.to_dict(),
        }
        with open(self.path, "a") as f:
            f.write(json.dumps(entry) + "\n")
        return cfg_hash

    def load_all(self) -> list:
        if not self.path.exists():
            return []
        with open(self.path) as f:
            return [json.loads(line) for line in f if line.strip()]
