"""ff-room: 3D indoor airflow simulation toolkit."""
from .config import RoomConfig, FanConfig, ObstacleConfig, SceneConfig
from .scene import Scene
from .solver_bridge import SolverBridge
from .visualization import Visualizer, ScenePlotter, plot_scene
from .io import ResultStore

__all__ = [
    "RoomConfig", "FanConfig", "ObstacleConfig", "SceneConfig",
    "Scene", "SolverBridge", "Visualizer", "ScenePlotter", "plot_scene",
    "ResultStore",
]
