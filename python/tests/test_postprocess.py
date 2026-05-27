"""Tests for SimResult post-processing (no C++ required)."""
import numpy as np
import pytest
from ff_room.config import RoomConfig, SceneConfig, FanConfig
from ff_room.solver_bridge import SimResult


def _mock_result(Nx=4, Ny=4, Nz=4) -> SimResult:
    """Create a SimResult with synthetic velocity data for testing."""
    config = SceneConfig(
        room=RoomConfig(size=(1.0, 1.0, 1.0), grid=(Nx, Ny, Nz)),
        fans=[FanConfig()],
    )
    # Uniform flow u=1.0, v=0, w=0
    u = np.ones((Nx+1, Ny, Nz))
    v = np.zeros((Nx, Ny+1, Nz))
    w = np.zeros((Nx, Ny, Nz+1))
    p = np.zeros((Nx, Ny, Nz))
    ct = np.zeros((Nx, Ny, Nz), dtype=np.uint8)
    return SimResult(u_field=u, v_field=v, w_field=w, p_field=p,
                     cell_type=ct, config=config)


class TestSimResultHelpers:
    def test_velocity_magnitude_uniform(self):
        r = _mock_result()
        speed = r.velocity_magnitude()
        assert speed.shape == (4, 4, 4)
        assert np.allclose(speed, 1.0)

    def test_velocity_cell_centered_shape(self):
        r = _mock_result(Nx=5, Ny=6, Nz=7)
        vel = r.velocity_cell_centered()
        assert vel.shape == (5, 6, 7, 3)

    def test_velocity_cell_centered_values(self):
        r = _mock_result()
        vel = r.velocity_cell_centered()
        assert np.allclose(vel[:,:,:,0], 1.0)
        assert np.allclose(vel[:,:,:,1], 0.0)
        assert np.allclose(vel[:,:,:,2], 0.0)


class TestResultStoreNpz:
    def test_save_load_roundtrip(self, tmp_path):
        from ff_room.io import ResultStore
        r = _mock_result()
        path = str(tmp_path / "test.npz")
        ResultStore.save_npz(r, path)
        r2 = ResultStore.load_npz(path)
        assert np.allclose(r2.u_field, r.u_field)
        assert np.allclose(r2.p_field, r.p_field)
        assert r2.config.room.size == r.config.room.size

    def test_save_velocity_csv(self, tmp_path):
        from ff_room.io import ResultStore
        r = _mock_result()
        path = str(tmp_path / "vel.csv")
        ResultStore.save_velocity_csv(r, path)
        import csv
        with open(path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 4*4*4
        assert "speed_ms" in rows[0]
        # All speeds should be ~1.0
        for row in rows:
            assert abs(float(row["speed_ms"]) - 1.0) < 1e-6
