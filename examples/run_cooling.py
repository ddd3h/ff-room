#!/usr/bin/env python3
"""
Summer ventilation example: hot room (32°C) + cool outside (25°C) + window + fan.

Usage:
  python examples/run_cooling.py [config.yaml]

Default config: examples/summer_cooling.yaml
To compare fan positions, copy the YAML, change fan.position or fan.direction,
run again, and compare T_mean convergence curves and PNG overviews.
"""
import matplotlib
matplotlib.use('Agg')

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from ff_room import SceneConfig, SolverBridge, Visualizer, ResultStore, ScenePlotter
from ff_room.io import ExperimentLog

config_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
              Path(__file__).parent / "summer_cooling.yaml"
out_dir = Path("results")
out_dir.mkdir(exist_ok=True)

print(f"Config: {config_path}")
config = SceneConfig.load_yaml(str(config_path))

sc = config.solver
print(f"Room:    {config.room.size} m   Grid: {config.room.grid}")
print(f"Fan:     {[f.name for f in config.fans]}")
print(f"Window:  {[o.name for o in config.openings]}")
print(f"Thermal: T_initial={sc.T_initial}°C  T_outside={sc.T_outside}°C  "
      f"T_target={sc.T_target}°C  buoyancy={sc.buoyancy}")
print()

stem = config_path.stem
scene_png = str(out_dir / f"{stem}_scene.png")
ScenePlotter(config).plot(save=scene_png, show=False)
print(f"Scene diagram: {scene_png}")
print()

bridge = SolverBridge(config)
result = bridge.run(print_interval=100)

sim_time_s = result.steps * sc.dt
print(f"\nDone: steps={result.steps}  sim_time={sim_time_s:.1f}s  "
      f"wall_time={result.wall_time_s:.1f}s  converged={result.converged}")
print(f"  T_mean={result.T_mean:.2f}°C  (target={sc.T_target}°C)  "
      f"div_max={result.divergence_max:.3e}")

npz_path = out_dir / f"{stem}.npz"
ResultStore.save_npz(result, str(npz_path))
print(f"Saved: {npz_path}")

log = ExperimentLog(str(out_dir / "cooling_log.jsonl"))
cfg_hash = log.append(result, notes=f"config={stem}")
print(f"Log: hash={cfg_hash}")

print("\nGenerating visualization...")
viz = Visualizer(result)
png_path = str(out_dir / f"{stem}_overview.png")
viz.save_panels(png_path)
print(f"Overview: {png_path}")
print("To compare fan positions: copy YAML, change fan.position, re-run, compare T_mean.")
