"""Visualization for SimResult.

Two backends:
  - matplotlib (always available, 2D slices/quivers, saves PNG)
  - pyvista     (3D interactive, requires OpenGL; falls back gracefully)
"""
from __future__ import annotations
from typing import Optional, Tuple
from pathlib import Path

import numpy as np

from .solver_bridge import SimResult


# ---------------------------------------------------------------------------
# Matplotlib backend (headless-safe)
# ---------------------------------------------------------------------------

class MatplotlibVisualizer:
    """2D slice plots using matplotlib. Works on WSL / headless servers."""

    def __init__(self, result: SimResult):
        self.result = result
        self.config = result.config
        r = self.config.room
        self.Lx, self.Ly, self.Lz = r.size
        self.Nx, self.Ny, self.Nz = r.grid

    def _speed(self) -> np.ndarray:
        return self.result.velocity_magnitude()

    def _vel(self) -> np.ndarray:
        return self.result.velocity_cell_centered()

    def _has_temperature(self) -> bool:
        T = self.result.T_field
        return float(T.max() - T.min()) > 0.05

    def slice_xy(self, z_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True, quiver_stride: int = 2):
        """Horizontal slice at z = z_frac * Lz."""
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        kz   = int(z_frac * self.Nz)
        kz   = max(0, min(self.Nz-1, kz))
        spd  = self._speed()[:, :, kz].T      # (Ny, Nx)
        vel  = self._vel()[:, :, kz, :]       # (Nx, Ny, 3)

        xs = np.linspace(0, self.Lx, self.Nx)
        ys = np.linspace(0, self.Ly, self.Ny)

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.pcolormesh(xs, ys, spd, cmap='viridis', shading='auto')
        plt.colorbar(im, ax=ax, label='speed [m/s]')

        # Quiver
        s = quiver_stride
        Xs, Ys = np.meshgrid(xs[::s], ys[::s])
        U = vel[::s, ::s, 0].T
        V = vel[::s, ::s, 1].T
        mag = np.sqrt(U**2 + V**2) + 1e-12
        ax.quiver(Xs, Ys, U/mag, V/mag, alpha=0.5, color='white',
                  scale=25, width=0.003)

        ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
        ax.set_title(f'Speed — horizontal slice z={kz*self.Lz/self.Nz:.2f} m')
        ax.set_aspect('equal')
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=150)
            print(f"Saved: {save}")
        if show:
            plt.show()
        return fig

    def slice_xz(self, y_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True, quiver_stride: int = 2):
        """Vertical slice at y = y_frac * Ly."""
        import matplotlib.pyplot as plt

        jy  = int(y_frac * self.Ny)
        jy  = max(0, min(self.Ny-1, jy))
        spd = self._speed()[:, jy, :].T     # (Nz, Nx)
        vel = self._vel()[:, jy, :, :]      # (Nx, Nz, 3)

        xs = np.linspace(0, self.Lx, self.Nx)
        zs = np.linspace(0, self.Lz, self.Nz)

        fig, ax = plt.subplots(figsize=(8, 5))
        im = ax.pcolormesh(xs, zs, spd, cmap='viridis', shading='auto')
        plt.colorbar(im, ax=ax, label='speed [m/s]')

        s = quiver_stride
        Xs, Zs = np.meshgrid(xs[::s], zs[::s])
        U = vel[::s, ::s, 0].T
        W = vel[::s, ::s, 2].T
        mag = np.sqrt(U**2 + W**2) + 1e-12
        ax.quiver(Xs, Zs, U/mag, W/mag, alpha=0.5, color='white',
                  scale=25, width=0.003)

        ax.set_xlabel('x [m]'); ax.set_ylabel('z [m]')
        ax.set_title(f'Speed — vertical slice y={jy*self.Ly/self.Ny:.2f} m')
        ax.set_aspect('equal')
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=150)
            print(f"Saved: {save}")
        if show:
            plt.show()
        return fig

    def slice_yz(self, x_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True, quiver_stride: int = 2):
        """Cross-section at x = x_frac * Lx."""
        import matplotlib.pyplot as plt

        ix  = int(x_frac * self.Nx)
        ix  = max(0, min(self.Nx-1, ix))
        spd = self._speed()[ix, :, :].T     # (Nz, Ny)
        vel = self._vel()[ix, :, :, :]      # (Ny, Nz, 3)

        ys = np.linspace(0, self.Ly, self.Ny)
        zs = np.linspace(0, self.Lz, self.Nz)

        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.pcolormesh(ys, zs, spd, cmap='viridis', shading='auto')
        plt.colorbar(im, ax=ax, label='speed [m/s]')

        s = quiver_stride
        Ys, Zs = np.meshgrid(ys[::s], zs[::s])
        V = vel[::s, ::s, 1].T
        W = vel[::s, ::s, 2].T
        mag = np.sqrt(V**2 + W**2) + 1e-12
        ax.quiver(Ys, Zs, V/mag, W/mag, alpha=0.5, color='white',
                  scale=25, width=0.003)

        ax.set_xlabel('y [m]'); ax.set_ylabel('z [m]')
        ax.set_title(f'Speed — cross-section x={ix*self.Lx/self.Nx:.2f} m')
        ax.set_aspect('equal')
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=150)
            print(f"Saved: {save}")
        if show:
            plt.show()
        return fig

    def multi_panel(self, save: Optional[str] = None, show: bool = True):
        """2×2 panel: 3 slices + speed profile along center line."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        fig.suptitle('ff-room airflow overview', fontsize=13)

        spd = self._speed()
        vel = self._vel()

        def _slice(ax, data_xy, xs, ys, vel_u, vel_v, xlabel, ylabel, title, s=2):
            im = ax.pcolormesh(xs, ys, data_xy.T, cmap='viridis', shading='auto')
            plt.colorbar(im, ax=ax, label='m/s')
            Xs, Ys = np.meshgrid(xs[::s], ys[::s])
            U = vel_u[::s, ::s].T
            V = vel_v[::s, ::s].T
            mag = np.sqrt(U**2 + V**2) + 1e-12
            ax.quiver(Xs, Ys, U/mag, V/mag, alpha=0.45, color='white',
                      scale=28, width=0.003)
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.set_title(title); ax.set_aspect('equal')

        xs = np.linspace(0, self.Lx, self.Nx)
        ys = np.linspace(0, self.Ly, self.Ny)
        zs = np.linspace(0, self.Lz, self.Nz)
        kz = self.Nz // 2
        jy = self.Ny // 2
        ix = self.Nx // 2

        _slice(axes[0,0], spd[:,:,kz], xs, ys,
               vel[:,:,kz,0], vel[:,:,kz,1],
               'x [m]', 'y [m]', f'z={kz*self.Lz/self.Nz:.2f}m (top view)')
        _slice(axes[0,1], spd[:,jy,:], xs, zs,
               vel[:,jy,:,0], vel[:,jy,:,2],
               'x [m]', 'z [m]', f'y={jy*self.Ly/self.Ny:.2f}m (side view)')
        _slice(axes[1,0], spd[ix,:,:], ys, zs,
               vel[ix,:,:,1], vel[ix,:,:,2],
               'y [m]', 'z [m]', f'x={ix*self.Lx/self.Nx:.2f}m (front view)')

        # Bottom-right: temperature slice (if thermal) or centerline speed profile
        ax = axes[1,1]
        if self._has_temperature():
            T = self.result.T_field
            T_slice = T[:, jy, :].T  # (Nz, Nx)
            im2 = ax.pcolormesh(xs, zs, T_slice, cmap='RdBu_r', shading='auto')
            plt.colorbar(im2, ax=ax, label='T [°C]')
            ax.set_xlabel('x [m]'); ax.set_ylabel('z [m]')
            ax.set_title(f'Temperature y={jy*self.Ly/self.Ny:.2f}m (T_mean={self.result.T_mean:.1f}°C)')
            ax.set_aspect('equal')
        else:
            profile = spd[:, self.Ny//2, self.Nz//2]
            ax.plot(xs, profile, 'b-o', markersize=3)
            ax.set_xlabel('x [m]'); ax.set_ylabel('speed [m/s]')
            ax.set_title('Centerline speed profile')
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=150)
            print(f"Saved: {save}")
        if show:
            plt.show()
        return fig


# ---------------------------------------------------------------------------
# PyVista backend (3D interactive)
# ---------------------------------------------------------------------------

def _try_pyvista():
    """Return pyvista module only if it can actually render. Subprocess-tests to avoid SIGSEGV."""
    import os, sys, subprocess

    # Suppress VTK stderr chatter (libGL errors etc.) during import
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    old2 = os.dup(2)
    os.dup2(devnull_fd, 2)
    os.close(devnull_fd)
    try:
        import pyvista as pv
        imported = True
    except ImportError:
        imported = False
    finally:
        os.dup2(old2, 2)
        os.close(old2)

    if not imported:
        return None

    # Validate that an off-screen renderer actually works (catches broken OpenGL).
    # Use subprocess so a SIGSEGV in VTK doesn't kill the main process.
    try:
        r = subprocess.run(
            [sys.executable, '-c',
             'import pyvista as pv; p = pv.Plotter(off_screen=True); p.close()'],
            capture_output=True, timeout=12,
        )
        return pv if r.returncode == 0 else None
    except Exception:
        return None


class PyVistaVisualizer:
    """3D interactive visualization using PyVista."""

    def __init__(self, result: SimResult):
        self.result = result
        self.config = result.config
        r = self.config.room
        self.Lx, self.Ly, self.Lz = r.size
        self.Nx, self.Ny, self.Nz = r.grid

    def _make_structured_grid(self):
        pv = _try_pyvista()
        r  = self.result
        x = np.linspace(self.Lx/self.Nx*0.5, self.Lx*(1-0.5/self.Nx), self.Nx)
        y = np.linspace(self.Ly/self.Ny*0.5, self.Ly*(1-0.5/self.Ny), self.Ny)
        z = np.linspace(self.Lz/self.Nz*0.5, self.Lz*(1-0.5/self.Nz), self.Nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        grid = pv.StructuredGrid(X, Y, Z)
        vel = r.velocity_cell_centered()
        grid.point_data["velocity"] = vel.reshape(-1, 3, order='C')
        grid.point_data["speed"]    = r.velocity_magnitude().ravel(order='C')
        grid.point_data["pressure"] = r.p_field.ravel(order='C')
        return grid

    def slice_xz(self, y_frac: float = 0.5, show: bool = True):
        pv = _try_pyvista()
        grid = self._make_structured_grid()
        slc  = grid.slice(normal="y", origin=(0, y_frac*self.Ly, 0))
        pl = pv.Plotter(title=f"y={y_frac*self.Ly:.2f}m slice")
        pl.add_mesh(slc, scalars="speed", cmap="viridis", show_scalar_bar=True)
        arrows = slc.glyph(orient="velocity", scale="speed", factor=0.3, tolerance=0.05)
        pl.add_mesh(arrows, color="white", opacity=0.6)
        pl.add_axes()
        if show: pl.show()
        return pl

    def streamlines(self, source_center=None, max_time: float = 50.0, show: bool = True):
        pv = _try_pyvista()
        grid = self._make_structured_grid()
        if source_center is None and self.config.fans:
            source_center = self.config.fans[0].position
        if source_center is None:
            source_center = (self.Lx/2, self.Ly/2, self.Lz/2)
        radius = min(self.Lx, self.Ly, self.Lz) * 0.1
        sphere = pv.Sphere(radius=radius, center=source_center,
                           theta_resolution=6, phi_resolution=6)
        streams = grid.streamlines_from_source(sphere, vectors="velocity",
                                               max_time=max_time,
                                               integration_direction="forward")
        pl = pv.Plotter(title="Streamlines")
        pl.add_mesh(streams.tube(radius=0.02), scalars="speed", cmap="plasma")
        box = pv.Box(bounds=(0, self.Lx, 0, self.Ly, 0, self.Lz))
        pl.add_mesh(box, style="wireframe", color="gray", opacity=0.3)
        pl.add_axes()
        if show: pl.show()
        return pl

    def multi_panel(self, show: bool = True):
        pv = _try_pyvista()
        grid = self._make_structured_grid()
        pl = pv.Plotter(shape=(2, 2), title="ff-room airflow overview")

        for (row, col, normal, pos, title) in [
            (0, 0, "z", (0, 0, self.Lz/2),  f"z={self.Lz/2:.1f}m"),
            (0, 1, "y", (0, self.Ly/2, 0),   f"y={self.Ly/2:.1f}m"),
            (1, 0, "x", (self.Lx/2, 0, 0),   f"x={self.Lx/2:.1f}m"),
        ]:
            pl.subplot(row, col)
            slc = grid.slice(normal=normal, origin=pos)
            pl.add_mesh(slc, scalars="speed", cmap="viridis", show_scalar_bar=True)
            arrows = slc.glyph(orient="velocity", scale="speed",
                               factor=0.3, tolerance=0.05)
            pl.add_mesh(arrows, color="white", opacity=0.5)
            pl.add_text(title, font_size=10)
            pl.add_axes()

        pl.subplot(1, 1)
        src = self.config.fans[0].position if self.config.fans \
              else (self.Lx/2, self.Ly/2, self.Lz/2)
        sphere = pv.Sphere(radius=min(self.Lx,self.Ly,self.Lz)*0.08, center=src,
                           theta_resolution=5, phi_resolution=5)
        try:
            streams = grid.streamlines_from_source(sphere, vectors="velocity",
                                                   max_time=30.0,
                                                   integration_direction="forward")
            pl.add_mesh(streams.tube(radius=0.015), scalars="speed", cmap="plasma")
        except Exception:
            pass
        pl.add_text("Streamlines", font_size=10)
        pl.add_axes()
        if show: pl.show()
        return pl


# ---------------------------------------------------------------------------
# Unified Visualizer: auto-selects backend
# ---------------------------------------------------------------------------

class Visualizer:
    """Auto-selects PyVista (3D) or matplotlib (2D) based on availability.

    Force a backend:
        Visualizer(result, backend='matplotlib')
        Visualizer(result, backend='pyvista')
    """

    def __init__(self, result: SimResult, backend: str = 'auto'):
        self.result = result
        self._mpl = MatplotlibVisualizer(result)
        self._pv_viz = None
        self._use_pyvista = False

        if backend == 'matplotlib':
            pass  # stay matplotlib
        elif backend == 'pyvista':
            self._use_pyvista = True
            self._pv_viz = PyVistaVisualizer(result)
        else:  # auto: only use pyvista when a display is likely available
            import os
            has_display = bool(os.environ.get('DISPLAY') or
                               os.environ.get('WAYLAND_DISPLAY'))
            if has_display and _try_pyvista() is not None:
                self._use_pyvista = True
                self._pv_viz = PyVistaVisualizer(result)

        backend_name = 'pyvista' if self._use_pyvista else 'matplotlib'
        print(f"[Visualizer] backend: {backend_name}")

    # --- delegate to active backend ---

    def multi_panel(self, save: Optional[str] = None, show: bool = True):
        if self._use_pyvista:
            return self._pv_viz.multi_panel(show=show)
        return self._mpl.multi_panel(save=save, show=show)

    def slice_xy(self, z_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True):
        if self._use_pyvista:
            pv = _try_pyvista()
            grid = self._pv_viz._make_structured_grid()
            slc  = grid.slice(normal="z",
                              origin=(0, 0, z_frac*self._mpl.Lz))
            pl = pv.Plotter()
            pl.add_mesh(slc, scalars="speed", cmap="viridis")
            if show: pl.show()
            return pl
        return self._mpl.slice_xy(z_frac=z_frac, save=save, show=show)

    def slice_xz(self, y_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True):
        if self._use_pyvista:
            return self._pv_viz.slice_xz(y_frac=y_frac, show=show)
        return self._mpl.slice_xz(y_frac=y_frac, save=save, show=show)

    def slice_yz(self, x_frac: float = 0.5, save: Optional[str] = None,
                 show: bool = True):
        return self._mpl.slice_yz(x_frac=x_frac, save=save, show=show)

    def streamlines(self, source_center=None, show: bool = True):
        if self._use_pyvista:
            return self._pv_viz.streamlines(source_center=source_center, show=show)
        print("[Visualizer] streamlines require pyvista; showing slice_xz instead")
        return self._mpl.slice_xz(show=show)

    # Always-available matplotlib methods
    def plot(self, save: Optional[str] = None, show: bool = True):
        """Alias for multi_panel."""
        return self.multi_panel(save=save, show=show)

    def save_panels(self, path: str = "results/overview.png"):
        """Save multi-panel PNG (no display needed)."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        return self._mpl.multi_panel(save=path, show=False)
