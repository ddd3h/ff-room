"""Tests for config serialization and FanConfig helpers."""
import math
import json
import tempfile
from pathlib import Path

import pytest
from ff_room.config import (
    RoomConfig, FanConfig, ObstacleConfig, SolverConfig, SceneConfig
)


def make_scene() -> SceneConfig:
    return SceneConfig(
        room=RoomConfig(size=(5.0, 4.0, 2.5), grid=(25, 20, 12)),
        fans=[
            FanConfig(position=(0.5, 2.0, 1.2), direction=(0.0, 0.0),
                      velocity=3.0, radius=0.15, name="fan0"),
            FanConfig(position=(4.5, 2.0, 1.2), direction=(180.0, 0.0),
                      velocity=2.0, radius=0.15, name="fan1"),
        ],
        obstacles=[ObstacleConfig(bbox_min=(2.0, 0.5, 0.0),
                                  bbox_max=(3.0, 1.8, 0.9))],
    )


class TestFanConfig:
    def test_direction_vector_plus_x(self):
        fan = FanConfig(direction=(0.0, 0.0))
        v = fan.direction_vector()
        assert abs(v[0] - 1.0) < 1e-9
        assert abs(v[1]) < 1e-9
        assert abs(v[2]) < 1e-9

    def test_direction_vector_minus_x(self):
        fan = FanConfig(direction=(180.0, 0.0))
        v = fan.direction_vector()
        assert abs(v[0] + 1.0) < 1e-6
        assert fan.snapped_axis() == 0
        assert fan.snapped_sign() == -1.0

    def test_direction_vector_plus_y(self):
        fan = FanConfig(direction=(90.0, 0.0))
        v = fan.direction_vector()
        assert abs(v[1] - 1.0) < 1e-6
        assert fan.snapped_axis() == 1

    def test_direction_vector_up(self):
        fan = FanConfig(direction=(0.0, 90.0))
        v = fan.direction_vector()
        assert abs(v[2] - 1.0) < 1e-6
        assert fan.snapped_axis() == 2

    def test_unit_vector_length(self):
        for az in [0, 30, 90, 135, 270]:
            for el in [-30, 0, 45, 90]:
                fan = FanConfig(direction=(az, el))
                v = fan.direction_vector()
                length = math.sqrt(sum(x**2 for x in v))
                assert abs(length - 1.0) < 1e-9, f"az={az} el={el}: |v|={length}"


class TestSceneConfigRoundTrip:
    def test_to_from_dict(self):
        s = make_scene()
        d = s.to_dict()
        s2 = SceneConfig.from_dict(d)
        assert s2.room.size == s.room.size
        assert len(s2.fans) == 2
        assert s2.fans[0].name == "fan0"
        assert s2.fans[1].velocity == 2.0
        assert len(s2.obstacles) == 1

    def test_yaml_roundtrip(self):
        s = make_scene()
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            path = f.name
        try:
            s.save_yaml(path)
            s2 = SceneConfig.load_yaml(path)
            assert s2.room.grid == s.room.grid
            assert s2.solver.nu == pytest.approx(s.solver.nu, rel=1e-9)
        finally:
            Path(path).unlink(missing_ok=True)

    def test_json_roundtrip(self):
        s = make_scene()
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            s.save_json(path)
            s2 = SceneConfig.load_json(path)
            assert s2.fans[0].position == s.fans[0].position
        finally:
            Path(path).unlink(missing_ok=True)
