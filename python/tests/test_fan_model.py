"""Physical units, outlet geometry, integrated flow, and solver compatibility."""
import math

import numpy as np
import pytest

from ff_room import FanConfig, RoomConfig, Scene, SceneConfig, SolverBridge
from ff_room import _ffroom_core as core


@pytest.mark.parametrize("kwargs", [
    {"radius": 0}, {"velocity": -1}, {"hub_ratio": 1},
    {"swirl_ratio": float("nan")}, {"flow_rate": -1},
    {"profile": "unknown"}, {"model": "unknown"},
    {"position": (0, float("inf"), 0)}, {"direction": (float("nan"), 0)},
])
def test_invalid_fan_config(kwargs):
    with pytest.raises(ValueError):
        FanConfig(**kwargs)


def test_flow_rate_units_and_serialization():
    fan = FanConfig(radius=0.25, hub_ratio=0.2, flow_rate=0.3, velocity=99,
                    direction=(35, 20), swirl_ratio=-0.15)
    assert fan.bulk_velocity() == pytest.approx(0.3 / (math.pi*0.25**2*0.96))
    config = SceneConfig(fans=[fan])
    restored = SceneConfig.from_dict(config.to_dict()).fans[0]
    assert restored == fan
    assert restored.bulk_velocity() == fan.bulk_velocity()
    patch = FanConfig(model="patch", radius=0.25, flow_rate=0.3)
    assert patch.bulk_velocity() == pytest.approx(0.3 / 0.25)


def disk():
    bc = core.DiskFanBC()
    bc.center = [0.5, 0.5, 0.5]
    bc.radius = 0.25
    bc.velocity = 2.0
    return bc


@pytest.mark.parametrize("uniform", [False, True])
def test_disk_flux_and_circular_mask(uniform):
    grid = core.Grid(4, 100, 100, 1.0, 1.0, 1.0)
    bm = core.BoundaryManager()
    bc = disk()
    bc.uniform_profile = uniform
    bm.add_disk_fan(grid, bc)
    bm.apply_inflow(grid)
    u = grid.u_shaped()[2]
    flux = u.sum() * 0.01**2
    target = bc.velocity * math.pi*bc.radius**2 * (1-bc.hub_ratio**2)
    assert flux == pytest.approx(target, rel=0.015)
    assert u[50, 50] == 0.0  # motor hub
    assert u[72, 72] == 0.0  # corner of the old square patch
    assert u.max() > 0
    np.testing.assert_allclose(grid.v_shaped(), 0.0)
    np.testing.assert_allclose(grid.w_shaped(), 0.0)


def test_swirl_reverses_without_changing_axial_flow():
    bc = disk()
    point = [0.5, 0.65, 0.5]
    bc.swirl_ratio = 0.25
    positive = np.array(bc.velocity_at(point))
    bc.swirl_ratio = -0.25
    negative = np.array(bc.velocity_at(point))
    assert positive[0] == negative[0]
    assert positive[2] > 0
    assert negative[2] == -positive[2]
    assert bc.velocity_at([0.5, 0.5, 0.5]) == [0.0, 0.0, 0.0]
    assert bc.velocity_at([0.5, 0.75, 0.5]) == [0.0, 0.0, 0.0]


def test_oblique_disk_on_stretched_grid():
    grid = core.Grid(20, 20, 20, 1.0, 1.0, 1.0)
    coords = np.linspace(0, 1, 21)**1.2
    grid.set_face_coords(coords, coords, coords)
    bc = disk()
    bc.direction = [1.0, 1.0, 1.0]
    bm = core.BoundaryManager()
    bm.add_disk_fan(grid, bc)
    indices = np.argwhere(grid.cell_type_shaped() == int(core.CellType.INFLOW))
    assert len(indices) > 0
    assert np.ptp(indices[:, 0]) > 1  # the plane is not snapped onto x
    centers = 0.5*(coords[:-1]+coords[1:])
    points = centers[indices]
    normal = np.ones(3)/np.sqrt(3)
    distances = (points-0.5) @ normal
    thickness = 0.5*(np.diff(coords)[indices] @ normal)
    assert np.all(np.abs(distances) <= thickness+1e-12)
    radial = points-0.5-distances[:, None]*normal
    assert np.all(np.linalg.norm(radial, axis=1) < bc.radius)


def test_fan_respects_solids_and_room_walls():
    grid = core.Grid(10, 40, 40, 1.0, 1.0, 1.0)
    bm = core.BoundaryManager()
    bm.add_solid_box(grid, 4, 4, 15, 24, 15, 24)
    bc = disk()
    bm.add_disk_fan(grid, bc)
    bm.apply_inflow(grid)
    assert np.all(grid.cell_type_shaped()[4, 15:25, 15:25] == int(core.CellType.SOLID))
    assert np.count_nonzero(grid.u_shaped()[4:6, 15:25, 15:25]) == 0
    bc.center = [0.0, 0.5, 0.5]
    bm.add_disk_fan(grid, bc)
    bm.apply_inflow(grid)
    assert np.count_nonzero(grid.u_shaped()[0]) == 0


def test_unresolved_disk_is_reported_and_off_fan_is_inactive():
    grid = core.Grid(4, 4, 4, 1.0, 1.0, 1.0)
    bc = disk()
    bc.radius = 0.01
    bm = core.BoundaryManager()
    with pytest.raises(ValueError, match="refine the grid"):
        bm.add_disk_fan(grid, bc)
    bc.velocity = 0.0
    bm.add_disk_fan(grid, bc)
    assert np.all(grid.cell_type_shaped() == int(core.CellType.FLUID))


@pytest.mark.parametrize("method", ["projection", "lbm"])
@pytest.mark.parametrize("model", ["disk", "patch"])
def test_fan_runs_in_both_solvers(method, model):
    config = SceneConfig(
        room=RoomConfig(size=(1.0, 1.0, 1.0), grid=(12, 12, 12)),
        fans=[FanConfig(position=(0.5, 0.5, 0.5), radius=0.25, velocity=0.1,
                        direction=(30, 15), model=model, swirl_ratio=0.1)],
    )
    config.solver.method = method
    config.solver.max_steps = 2
    result = SolverBridge(config).run(print_interval=0)
    assert result.steps > 0
    assert np.isfinite(result.velocity_cell_centered()).all()
    assert np.isfinite(result.p_field).all()


def test_scene_warns_when_fan_is_underresolved():
    config = SceneConfig(room=RoomConfig(size=(1, 1, 1), grid=(10, 10, 10)),
                         fans=[FanConfig(position=(0.5, 0.5, 0.5), radius=0.15)])
    with pytest.warns(UserWarning, match="fewer than four cells"):
        Scene(config)


@pytest.mark.parametrize("second_model", ["disk", "patch"])
def test_overlapping_fans_are_reported_without_mutating_the_grid(second_model):
    grid = core.Grid(10, 40, 40, 1.0, 1.0, 1.0)
    bm = core.BoundaryManager()
    bc = disk()
    bm.add_disk_fan(grid, bc)
    before = grid.cell_type_shaped().copy()
    if second_model == "patch":
        bc = core.FanBC()
        bc.i_min, bc.i_max = 4, 5
        bc.j_min, bc.j_max = 15, 24
        bc.k_min, bc.k_max = 15, 24
        bc.vel = [2.0, 0.0, 0.0]
    with pytest.raises(ValueError, match="outlets overlap"):
        if second_model == "disk":
            bm.add_disk_fan(grid, bc)
        else:
            bm.add_fan(grid, bc)
    np.testing.assert_array_equal(grid.cell_type_shaped(), before)
