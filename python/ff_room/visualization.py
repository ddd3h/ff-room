"""Visualization for SimResult and SceneConfig.

Two backends for flow results:
  - matplotlib (always available, 2D slices/quivers, saves PNG)
  - pyvista     (3D interactive, requires OpenGL; falls back gracefully)

ScenePlotter: room setup diagram (no simulation needed).
"""
from __future__ import annotations
from typing import Optional, TYPE_CHECKING
from pathlib import Path

import numpy as np

from .solver_bridge import SimResult

if TYPE_CHECKING:
    from .config import SceneConfig


# ---------------------------------------------------------------------------
# Scene setup diagram (SceneConfig only — no simulation needed)
# ---------------------------------------------------------------------------

class ScenePlotter:
    """Room layout diagram from SceneConfig: floor plan + two cross-sections.

    Usage::

        from ff_room import SceneConfig
        from ff_room.visualization import ScenePlotter

        cfg = SceneConfig.load_yaml("examples/summer_cooling.yaml")
        ScenePlotter(cfg).plot(save="results/scene.png", show=False)
    """

    # Visual constants
    _C_OBS  = '#8B7355'   # obstacle: brown
    _C_FAN  = '#1E88E5'   # fan marker: blue
    _C_ARR  = '#0D47A1'   # fan arrow: dark blue
    _C_OPEN = '#43A047'   # opening: green

    def __init__(self, config: "SceneConfig"):
        from .config import SceneConfig  # avoid circular at module level
        self.cfg = config
        r = config.room
        self.Lx, self.Ly, self.Lz = r.size

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def plot(self, save: Optional[str] = None, show: bool = True):
        """3-panel setup diagram: top view + YZ section + XZ section."""
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D

        fig, (ax_top, ax_yz, ax_xz) = plt.subplots(1, 3, figsize=(16, 5))

        Lx, Ly, Lz = self.Lx, self.Ly, self.Lz
        cfg = self.cfg
        arrow_len = max(Lx, Ly, Lz) * 0.12

        # Cross-section planes pass through first fan (or room center)
        if cfg.fans:
            sx, sy = cfg.fans[0].position[0], cfg.fans[0].position[1]
        else:
            sx, sy = Lx / 2, Ly / 2

        # --- room outlines ---
        self._room_box(ax_top, 0, Lx, 0, Ly, 'x [m]', 'y [m]', 'Top view (XY — floor plan)',
                       labels={'bottom': 'South (y=0)', 'top': 'North (y=Ly)',
                               'left': 'West', 'right': 'East'})
        self._room_box(ax_yz, 0, Ly, 0, Lz, 'y [m]', 'z [m]',
                       f'Section YZ  (x = {sx:.2f} m)',
                       labels={'bottom': 'South', 'top': 'North'})
        self._room_box(ax_xz, 0, Lx, 0, Lz, 'x [m]', 'z [m]',
                       f'Section XZ  (y = {sy:.2f} m)',
                       labels={'left': 'West', 'right': 'East'})

        # --- obstacles ---
        for obs in cfg.obstacles:
            bmi, bma = obs.bbox_min, obs.bbox_max
            self._rect(ax_top, bmi[0], bmi[1], bma[0]-bmi[0], bma[1]-bmi[1],
                       obs.name, self._C_OBS)
            if bmi[0] <= sx <= bma[0]:
                self._rect(ax_yz, bmi[1], bmi[2], bma[1]-bmi[1], bma[2]-bmi[2],
                           obs.name, self._C_OBS)
            if bmi[1] <= sy <= bma[1]:
                self._rect(ax_xz, bmi[0], bmi[2], bma[0]-bmi[0], bma[2]-bmi[2],
                           obs.name, self._C_OBS)

        # --- openings ---
        for op in cfg.openings:
            self._draw_opening(ax_top, ax_yz, ax_xz, op, Lx, Ly, Lz)

        # --- fans ---
        for fan in cfg.fans:
            fx, fy, fz = fan.position
            dv = fan.direction_vector()  # (dx, dy, dz) unit vector

            self._fan_marker(ax_top, fx, fy, dv[0], dv[1], fan.name, arrow_len)
            self._fan_marker(ax_yz,  fy, fz, dv[1], dv[2], fan.name, arrow_len)
            self._fan_marker(ax_xz,  fx, fz, dv[0], dv[2], fan.name, arrow_len)

        # --- figure title ---
        sc = cfg.solver
        name = cfg.metadata.get('name', cfg.metadata.get('description', 'room'))
        if isinstance(name, str) and len(name) > 40:
            name = name[:37] + '...'
        thermal = ''
        if sc.T_initial != sc.T_outside or sc.T_target is not None:
            thermal = f'  |  T_in={sc.T_initial}°C  T_out={sc.T_outside}°C'
            if sc.T_target is not None:
                thermal += f'  →  target={sc.T_target}°C'
        title = (f'{name}  |  {Lx:.1f}×{Ly:.1f}×{Lz:.1f} m'
                 f'  |  {len(cfg.fans)} fan(s)  {len(cfg.openings)} opening(s){thermal}')
        fig.suptitle(title, fontsize=10)

        # --- legend ---
        legend_elems = [
            mpatches.Patch(facecolor=self._C_OBS, edgecolor='#5D4037',
                           alpha=0.7, label='Obstacle'),
            Line2D([0], [0], color=self._C_OPEN, lw=5, label='Opening (window/door)'),
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=self._C_FAN, markersize=9, label='Fan'),
            Line2D([0], [0], color=self._C_ARR, lw=2.5, marker='>',
                   markersize=8, label='Fan direction'),
        ]
        fig.legend(handles=legend_elems, loc='lower center', ncol=4,
                   bbox_to_anchor=(0.5, -0.06), fontsize=8)

        plt.tight_layout()
        if save:
            Path(save).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save, dpi=150, bbox_inches='tight')
            print(f"Saved scene diagram: {save}")
        if show:
            plt.show()
        return fig

    # ------------------------------------------------------------------
    # Drawing helpers
    # ------------------------------------------------------------------

    def _room_box(self, ax, x0, x1, y0, y1,
                  xlabel, ylabel, title, labels=None):
        import matplotlib.patches as mpatches
        ax.add_patch(mpatches.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            fill=False, edgecolor='black', lw=2.0, zorder=2))
        ax.set_xlim(x0 - 0.4, x1 + 0.4)
        ax.set_ylim(y0 - 0.3, y1 + 0.4)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=9)
        ax.set_aspect('equal')
        if labels:
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            kw = dict(fontsize=6, color='#888888')
            if 'bottom' in labels:
                ax.text(cx, y0 - 0.12, labels['bottom'], ha='center', va='top', **kw)
            if 'top' in labels:
                ax.text(cx, y1 + 0.12, labels['top'], ha='center', va='bottom', **kw)
            if 'left' in labels:
                ax.text(x0 - 0.12, cy, labels['left'], ha='right', va='center',
                        rotation=90, **kw)
            if 'right' in labels:
                ax.text(x1 + 0.12, cy, labels['right'], ha='left', va='center',
                        rotation=90, **kw)

    def _rect(self, ax, x, y, w, h, label, color):
        import matplotlib.patches as mpatches
        ax.add_patch(mpatches.Rectangle(
            (x, y), w, h,
            facecolor=color, edgecolor='#5D4037', alpha=0.65, lw=1.2, zorder=3))
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=6, color='white', zorder=4)

    def _fan_marker(self, ax, px, pz, dvx, dvz, name, arrow_len):
        ax.plot(px, pz, 'o', color=self._C_FAN, ms=9, zorder=6)
        mag = (dvx**2 + dvz**2) ** 0.5
        if mag > 0.05:
            ax.annotate(
                '', xy=(px + dvx / mag * arrow_len, pz + dvz / mag * arrow_len),
                xytext=(px, pz),
                arrowprops=dict(arrowstyle='->', color=self._C_ARR, lw=2.2),
                zorder=6)
        ax.text(px, pz + arrow_len * 0.15 + 0.08, name,
                ha='center', va='bottom', fontsize=6.5,
                color=self._C_FAN, fontweight='bold', zorder=7)

    def _draw_opening(self, ax_top, ax_yz, ax_xz, op, Lx, Ly, Lz):
        """Draw opening on all relevant panels."""
        wall = op.wall.lower()
        lw = 5.5
        c  = self._C_OPEN

        def _label(ax, x, y, name, ha='center', va='bottom'):
            ax.text(x, y, name, ha=ha, va=va, fontsize=6, color=c, zorder=7)

        if wall == 'south':    # y=0  center=(x_c, z_c)  size=(w_x, h_z)
            x_c, z_c = op.center;  w_x, h_z = op.size
            ax_top.plot([x_c - w_x/2, x_c + w_x/2], [0, 0], '-', color=c, lw=lw, zorder=5)
            _label(ax_top, x_c, 0.08, op.name)
            ax_yz.plot([0, 0], [z_c - h_z/2, z_c + h_z/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_yz, 0.08, z_c, op.name, ha='left', va='center')

        elif wall == 'north':  # y=Ly center=(x_c, z_c)  size=(w_x, h_z)
            x_c, z_c = op.center;  w_x, h_z = op.size
            ax_top.plot([x_c - w_x/2, x_c + w_x/2], [Ly, Ly], '-', color=c, lw=lw, zorder=5)
            _label(ax_top, x_c, Ly - 0.08, op.name, va='top')
            ax_yz.plot([Ly, Ly], [z_c - h_z/2, z_c + h_z/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_yz, Ly - 0.08, z_c, op.name, ha='right', va='center')

        elif wall == 'west':   # x=0  center=(y_c, z_c)  size=(w_y, h_z)
            y_c, z_c = op.center;  w_y, h_z = op.size
            ax_top.plot([0, 0], [y_c - w_y/2, y_c + w_y/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_top, 0.08, y_c, op.name, ha='left', va='center')
            ax_xz.plot([0, 0], [z_c - h_z/2, z_c + h_z/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_xz, 0.08, z_c, op.name, ha='left', va='center')

        elif wall == 'east':   # x=Lx center=(y_c, z_c)  size=(w_y, h_z)
            y_c, z_c = op.center;  w_y, h_z = op.size
            ax_top.plot([Lx, Lx], [y_c - w_y/2, y_c + w_y/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_top, Lx - 0.08, y_c, op.name, ha='right', va='center')
            ax_xz.plot([Lx, Lx], [z_c - h_z/2, z_c + h_z/2], '-', color=c, lw=lw, zorder=5)
            _label(ax_xz, Lx - 0.08, z_c, op.name, ha='right', va='center')

        elif wall == 'floor':   # z=0  center=(x_c, y_c)  size=(w_x, w_y)
            x_c, y_c = op.center;  w_x, w_y = op.size
            import matplotlib.patches as mpatches
            ax_top.add_patch(mpatches.Rectangle(
                (x_c - w_x/2, y_c - w_y/2), w_x, w_y,
                facecolor=c, alpha=0.35, edgecolor=c, lw=1.5, zorder=5))
            ax_xz.plot([x_c - w_x/2, x_c + w_x/2], [0, 0], '-', color=c, lw=lw, zorder=5)
            _label(ax_xz, x_c, 0.06, op.name)

        elif wall == 'ceiling':  # z=Lz center=(x_c, y_c)  size=(w_x, w_y)
            x_c, y_c = op.center;  w_x, w_y = op.size
            import matplotlib.patches as mpatches
            ax_top.add_patch(mpatches.Rectangle(
                (x_c - w_x/2, y_c - w_y/2), w_x, w_y,
                facecolor=c, alpha=0.35, edgecolor=c, lw=1.5, zorder=5))
            ax_xz.plot([x_c - w_x/2, x_c + w_x/2], [Lz, Lz], '-', color=c, lw=lw, zorder=5)
            _label(ax_xz, x_c, Lz - 0.06, op.name, va='top')


def plot_scene(config: "SceneConfig", save: Optional[str] = None,
               show: bool = True):
    """Convenience wrapper: draw and optionally save room setup diagram."""
    return ScenePlotter(config).plot(save=save, show=show)


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

    def save_scene(self, path: str = "results/scene.png"):
        """Save room setup diagram (fan/obstacles/openings, no flow data needed)."""
        return ScenePlotter(self.result.config).plot(save=path, show=False)

    def save_animation(self, snapshots, path: str = "results/animation.gif",
                       fps: int = 10):
        """Save animated GIF/MP4 from snapshots captured with run_animated().

        Parameters
        ----------
        snapshots : list[FieldSnapshot]
            Snapshots returned by SolverBridge.run_animated().
        path : str
            Output file path.  Use .gif for GIF, .mp4 for MP4.
        fps : int
            Frames per second.
        """
        from .animation import Animator
        Animator(snapshots, self.result.config).save_gif(path, fps=fps)
