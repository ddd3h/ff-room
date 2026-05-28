"""Animation / video export for ff-room simulation results.

Animator takes a list of FieldSnapshot objects (captured by
SolverBridge.run_animated) and renders each as a 2-panel matplotlib figure:

  Left panel : top view  (XY slice at z = Nz//2) — speed colormap + quiver
  Right panel: side view (YZ slice at x = Nx//2) — speed colormap, or
               temperature colormap (RdBu_r) when T varies > 0.1 °C

Output:
  save_gif(path)  — animated GIF via Pillow (always available)
  save_mp4(path)  — MP4 via FFMpegWriter; falls back to GIF if ffmpeg absent
"""
from __future__ import annotations

from pathlib import Path
from typing import List, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .config import SceneConfig
    from .solver_bridge import FieldSnapshot


class Animator:
    """Render a sequence of FieldSnapshot objects as an animated GIF or MP4.

    Parameters
    ----------
    snapshots : list[FieldSnapshot]
        Snapshots captured by SolverBridge.run_animated().
    config : SceneConfig
        Scene configuration (used for physical dimensions).
    """

    def __init__(self, snapshots: "List[FieldSnapshot]", config: "SceneConfig"):
        self.snapshots = snapshots
        self.config    = config
        r = config.room
        self.Lx, self.Ly, self.Lz = r.size
        self.Nx, self.Ny, self.Nz = r.grid

        # Decide once whether this is a thermal run
        T_range = max(s.T_field.max() - s.T_field.min() for s in snapshots) \
                  if snapshots else 0.0
        self._is_thermal: bool = T_range > 0.1

        # Pre-compute global color limits for consistent colorbar across frames
        all_speeds = np.concatenate([s.speed.ravel() for s in snapshots])
        self._vmin_spd = float(all_speeds.min())
        self._vmax_spd = float(np.percentile(all_speeds, 99))  # clip outliers

        if self._is_thermal:
            all_T = np.concatenate([s.T_field.ravel() for s in snapshots])
            self._vmin_T = float(all_T.min())
            self._vmax_T = float(all_T.max())

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def save_gif(self, path: str, fps: int = 10, dpi: int = 100) -> None:
        """Save animation as an animated GIF (requires Pillow)."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.animation as manim

        fig, axes = self._make_figure()
        anim = manim.FuncAnimation(
            fig, self._update_frame, frames=len(self.snapshots),
            fargs=(axes,), interval=1000 // fps, blit=False,
        )

        Path(path).parent.mkdir(parents=True, exist_ok=True)
        writer = manim.PillowWriter(fps=fps)
        anim.save(path, writer=writer, dpi=dpi)
        plt.close(fig)
        print(f"GIF saved: {path}  ({len(self.snapshots)} frames, {fps} fps)")

    def save_mp4(self, path: str, fps: int = 15, dpi: int = 120) -> None:
        """Save animation as MP4 (requires ffmpeg).

        Falls back to GIF (same path with .gif extension) if ffmpeg is not
        available and prints a warning.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.animation as manim

        if not manim.FFMpegWriter.isAvailable():
            gif_path = str(Path(path).with_suffix('.gif'))
            print(f"[Animator] WARNING: ffmpeg not available — saving GIF instead: {gif_path}")
            self.save_gif(gif_path, fps=fps, dpi=dpi)
            return

        import matplotlib.pyplot as plt

        fig, axes = self._make_figure()
        anim = manim.FuncAnimation(
            fig, self._update_frame, frames=len(self.snapshots),
            fargs=(axes,), interval=1000 // fps, blit=False,
        )

        Path(path).parent.mkdir(parents=True, exist_ok=True)
        writer = manim.FFMpegWriter(fps=fps, bitrate=1800)
        anim.save(path, writer=writer, dpi=dpi)
        plt.close(fig)
        print(f"MP4 saved: {path}  ({len(self.snapshots)} frames, {fps} fps)")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_figure(self):
        """Create the 2-panel figure and return (fig, (ax_left, ax_right))."""
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        fig.subplots_adjust(wspace=0.35)
        return fig, axes

    def _update_frame(self, frame_idx: int, axes) -> None:
        """Draw one frame into *axes* (called by FuncAnimation)."""
        import matplotlib.pyplot as plt

        snap = self.snapshots[frame_idx]
        ax_l, ax_r = axes
        ax_l.cla()
        ax_r.cla()

        xs = np.linspace(0, self.Lx, self.Nx)
        ys = np.linspace(0, self.Ly, self.Ny)
        ys2 = np.linspace(0, self.Ly, self.Ny)
        zs = np.linspace(0, self.Lz, self.Nz)

        kz = self.Nz // 2   # top-view z index
        ix = self.Nx // 2   # side-view x index

        # ---- Left: top view (XY at z=kz) ----
        spd_xy = snap.speed[:, :, kz].T    # (Ny, Nx)
        vel_xy = snap.vel[:, :, kz, :]     # (Nx, Ny, 3)

        im_l = ax_l.pcolormesh(
            xs, ys, spd_xy,
            cmap='viridis', shading='auto',
            vmin=self._vmin_spd, vmax=self._vmax_spd,
        )
        plt.colorbar(im_l, ax=ax_l, label='speed [m/s]', fraction=0.046, pad=0.04)

        s = max(1, min(self.Nx, self.Ny) // 15)   # quiver stride (2–3 cells)
        s = max(2, s)
        Xs, Ys = np.meshgrid(xs[::s], ys[::s])
        U = vel_xy[::s, ::s, 0].T
        V = vel_xy[::s, ::s, 1].T
        mag = np.sqrt(U**2 + V**2) + 1e-12
        ax_l.quiver(Xs, Ys, U / mag, V / mag,
                    alpha=0.55, color='white', scale=25, width=0.003)

        ax_l.set_xlabel('x [m]')
        ax_l.set_ylabel('y [m]')
        ax_l.set_title(f'Top view  z={kz * self.Lz / self.Nz:.2f} m', fontsize=9)
        ax_l.set_aspect('equal')

        # ---- Right: side view (YZ at x=ix) ----
        if self._is_thermal:
            T_yz = snap.T_field[ix, :, :].T    # (Nz, Ny)
            im_r = ax_r.pcolormesh(
                ys2, zs, T_yz,
                cmap='RdBu_r', shading='auto',
                vmin=self._vmin_T, vmax=self._vmax_T,
            )
            plt.colorbar(im_r, ax=ax_r, label='T [°C]', fraction=0.046, pad=0.04)
            ax_r.set_title(f'Temp side view  x={ix * self.Lx / self.Nx:.2f} m', fontsize=9)
        else:
            spd_yz = snap.speed[ix, :, :].T    # (Nz, Ny)
            im_r = ax_r.pcolormesh(
                ys2, zs, spd_yz,
                cmap='viridis', shading='auto',
                vmin=self._vmin_spd, vmax=self._vmax_spd,
            )
            plt.colorbar(im_r, ax=ax_r, label='speed [m/s]', fraction=0.046, pad=0.04)
            ax_r.set_title(f'Side view  x={ix * self.Lx / self.Nx:.2f} m', fontsize=9)

        ax_r.set_xlabel('y [m]')
        ax_r.set_ylabel('z [m]')
        ax_r.set_aspect('equal')

        # ---- Super-title with time / temperature ----
        if self._is_thermal:
            suptitle = (f't = {snap.time_s:.1f} s    '
                        f'T_mean = {snap.T_mean:.2f} °C    '
                        f'step {snap.step}')
        else:
            suptitle = f't = {snap.time_s:.1f} s    step {snap.step}'
        ax_l.get_figure().suptitle(suptitle, fontsize=10, y=1.01)
