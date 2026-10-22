"""Exercise the installed extension, NumPy views, and both solver bindings."""
import numpy as np
import pytest

from ff_room import _ffroom_core as core


def test_grid_numpy_views():
    grid = core.Grid(4, 5, 6, 1.0, 1.0, 1.0)
    velocity = grid.u_shaped()
    assert velocity.shape == (5, 5, 6)
    assert velocity.dtype == np.float64
    velocity[1, 2, 3] = 2.5
    assert np.shares_memory(velocity, grid.u)
    assert grid.u_shaped()[1, 2, 3] == 2.5
    grid.zero_velocity()
    assert np.count_nonzero(velocity) == 0
    assert grid.cell_type_shaped().dtype == np.uint8


@pytest.mark.parametrize("method", ["projection", "lbm"])
def test_solver_callback(method):
    grid = core.Grid(4, 4, 4, 1.0, 1.0, 1.0)
    boundary = core.BoundaryManager()
    if method == "projection":
        params = core.FluidSolverParams()
        params.max_steps = 2
        solver = core.FluidSolver(grid, boundary, params)
    else:
        solver = core.LBMSolver(grid, boundary, 0.8, 2, 1e-6)
    steps = []
    result = solver.run(lambda step: steps.append(step.step))
    assert steps
    assert result.step == steps[-1]
    assert 1 <= result.step <= 2
    assert np.isfinite(result.divergence_max)
    assert np.isfinite(grid.u_shaped()).all()
    if method == "projection":
        assert result.converged
        np.testing.assert_allclose(grid.u_shaped(), 0.0, atol=1e-12)


def test_public_api_and_matplotlib(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ff_room import RoomConfig, SceneConfig, SolverBridge, Visualizer

    config = SceneConfig(room=RoomConfig(size=(1.0, 1.0, 1.0), grid=(4, 4, 4)))
    config.solver.max_steps = 2
    result = SolverBridge(config).run(print_interval=0)
    assert result.converged
    assert result.velocity_magnitude().shape == (4, 4, 4)
    path = tmp_path / "slice.png"
    figure = Visualizer(result, backend="matplotlib").slice_xz(
        save=str(path), show=False
    )
    assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    plt.close(figure)
