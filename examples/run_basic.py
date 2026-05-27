#!/usr/bin/env python3
"""Run basic_room example and visualize results."""
# Force non-interactive backend before any import that might pull in matplotlib.
# Without this, TkAgg on WSL / headless servers crashes with SIGSEGV.
import matplotlib
matplotlib.use('Agg')

import sys
from pathlib import Path

# Allow running from project root without install
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from ff_room import SceneConfig, SolverBridge, Visualizer, ResultStore
from ff_room.io import ExperimentLog

config_path = Path(__file__).parent / "basic_room.yaml"
out_dir     = Path("results")
out_dir.mkdir(exist_ok=True)

print(f"Loading config: {config_path}")
config = SceneConfig.load_yaml(str(config_path))

print(f"Room: {config.room.size} m  Grid: {config.room.grid}")
print(f"Fans: {[f.name for f in config.fans]}")
print(f"Obstacles: {[o.name for o in config.obstacles]}")
print(f"Max steps: {config.solver.max_steps}")
print()

bridge = SolverBridge(config)
result = bridge.run(print_interval=50)

print(f"\nDone: steps={result.steps}, converged={result.converged}, "
      f"wall_time={result.wall_time_s:.1f}s")
print(f"  max divergence: {result.divergence_max:.3e}")

# Save results
npz_path = out_dir / "basic_room.npz"
ResultStore.save_npz(result, str(npz_path))
print(f"Saved: {npz_path}")

csv_path = out_dir / "basic_room_velocity.csv"
ResultStore.save_velocity_csv(result, str(csv_path))
print(f"Saved: {csv_path}")

# Experiment log
log = ExperimentLog(str(out_dir / "experiment_log.jsonl"))
cfg_hash = log.append(result, notes="basic single-fan run")
print(f"Log entry: hash={cfg_hash}")

# Visualize — always saves PNG (Agg backend set at top of file)
print("\nGenerating visualization...")
viz = Visualizer(result)

png_path = str(out_dir / "basic_room_overview.png")
viz.save_panels(png_path)
print(f"Overview saved: {png_path}")
print("To view interactively: open the PNG or run from a Jupyter notebook.")
