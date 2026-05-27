# ff-room: 3D室内気流シミュレーション

部屋に扇風機を置いたときの気流を3Dで可視化するCFDツール。
扇風機の位置・向き・風量を変えながら流れの違いを比較できる。研究用途向けに再現性を重視。

**数値手法**: 非圧縮性Navier-Stokes方程式 + Chorin射影法 on MAC格子  
**言語**: C++17 (ソルバー) + Python 3.9+ (UI・可視化)  
**連携**: pybind11

---

## 目次

- [インストール](#インストール)
- [クイックスタート](#クイックスタート)
- [設定ファイル](#設定ファイル)
- [可視化](#可視化)
- [結果の保存と再利用](#結果の保存と再利用)
- [実験ログ](#実験ログ)
- [ディレクトリ構成](#ディレクトリ構成)
- [テスト実行](#テスト実行)
- [ドキュメント一覧](#ドキュメント一覧)

---

## インストール

### 前提条件

| ツール | バージョン |
|--------|-----------|
| Python | >= 3.9 |
| C++ コンパイラ | GCC 9+ / Clang 10+ (C++17対応) |
| CMake | >= 3.15 |
| pybind11 | >= 2.11 (自動DL可) |

### 手順

```bash
# 1. Python依存パッケージ
pip install pybind11 pyyaml numpy pyvista

# 2. C++拡張モジュールをビルド
bash build.sh

# または pip install (scikit-build-core経由)
pip install -e .
```

`build.sh` は `cpp/build/` でcmakeを走らせ、生成した `_ffroom_core*.so` を
`python/ff_room/` にインストールする。

---

## クイックスタート

```bash
# サンプル実行 (ビルド後)
python examples/run_basic.py
```

コードから使う場合:

```python
from ff_room import SceneConfig, SolverBridge, Visualizer

# 設定ファイル読み込み
config = SceneConfig.load_yaml("examples/basic_room.yaml")

# シミュレーション実行
result = SolverBridge(config).run(print_interval=50)
print(f"収束: {result.converged}, ステップ数: {result.steps}")

# 可視化
viz = Visualizer(result)
viz.multi_panel()          # 2×2パネル (断面×3 + 流線)
viz.streamlines()          # 3D流線
viz.slice_xz(y_frac=0.5)  # 中央垂直断面
```

---

## 設定ファイル

YAML形式 (JSONも可)。設定ファイルから**完全に**再現できることを保証する。

```yaml
metadata:
  name: my_experiment
  description: 説明を書く
  author: ""
  date: ""

room:
  size: [5.0, 4.0, 2.5]   # 部屋サイズ Lx, Ly, Lz [m]
  grid: [25, 20, 12]       # 格子分割数 Nx, Ny, Nz

fans:
  - name: fan0
    position: [0.5, 2.0, 1.2]   # 扇風機中心座標 [m]
    direction: [0.0, 0.0]        # [方位角_deg, 仰角_deg]
    velocity: 3.0                 # 吹き出し風速 [m/s]
    radius: 0.15                  # 扇風機半径 [m] (正方形近似: 辺=2×radius)

obstacles:
  - name: sofa
    bbox_min: [2.0, 0.5, 0.0]   # 直方体下端 [m]
    bbox_max: [3.0, 1.8, 0.9]   # 直方体上端 [m]

solver:
  dt: 0.01                # タイムステップ [s]
  max_steps: 500          # 最大反復回数
  convergence_tol: 1.0e-4 # 収束判定 (max |Δu| / dt)
  rho: 1.2                # 空気密度 [kg/m³]
  nu: 1.5e-5              # 動粘度 [m²/s]
  poisson_max_iter: 1000  # Poisson CG最大反復
  poisson_tol: 1.0e-6     # Poisson収束閾値
```

### 扇風機の向き

| direction | 意味 |
|-----------|------|
| `[0, 0]`   | +x 方向 (東) |
| `[90, 0]`  | +y 方向 (北) |
| `[180, 0]` | -x 方向 (西) |
| `[270, 0]` | -y 方向 (南) |
| `[0, 90]`  | 上向き (+z) |
| `[0, -90]` | 下向き (-z) |

MVP: 方向は最近軸にスナップされる (斜め向きは非対応)。

### 座標系

```
z (高さ・上)
│
│   y (奥行き・北)
│  /
│ /
└──── x (横・東)

原点: 部屋の南西床角 (0, 0, 0)
```

---

## 可視化

```python
from ff_room import Visualizer

viz = Visualizer(result)

# --- 断面スライス ---
viz.slice_xy(z_frac=0.5)       # 水平断面 z=Lz/2
viz.slice_xz(y_frac=0.5)       # 垂直断面 y=Ly/2
viz.slice_yz(x_frac=0.5)       # 横断面   x=Lx/2

# z_frac=0.0〜1.0 で断面位置を変える
viz.slice_xy(z_frac=0.48)      # ほぼ床面近く

# --- 流線 ---
viz.streamlines(
    source_center=(0.5, 2.0, 1.2),  # 種点の中心 (省略時: 最初の扇風機)
    n_points=50,
    max_time=50.0,
)

# --- ボリュームレンダリング ---
viz.volume_speed()

# --- 2×2まとめパネル ---
viz.multi_panel()
```

PyVistaウィンドウ: マウスドラッグで回転、スクロールでズーム。

---

## 結果の保存と再利用

```python
from ff_room.io import ResultStore

# 保存
ResultStore.save_npz(result, "results/run1.npz")         # NumPy圧縮
ResultStore.save_vtk(result, "results/run1.vts")          # ParaView用VTK
ResultStore.save_velocity_csv(result, "results/run1.csv") # CSV (x,y,z,u,v,w,speed)

# 読み込んで再可視化
r2 = ResultStore.load_npz("results/run1.npz")
Visualizer(r2).multi_panel()
```

### NumPy配列として取得

```python
import numpy as np

# 速度場 (各軸のface-centered)
u = result.u_field   # shape: (Nx+1, Ny, Nz)
v = result.v_field   # shape: (Nx, Ny+1, Nz)
w = result.w_field   # shape: (Nx, Ny, Nz+1)
p = result.p_field   # shape: (Nx, Ny, Nz)

# セル中心ベクトル場 (Nx, Ny, Nz, 3)
vel = result.velocity_cell_centered()
speed = result.velocity_magnitude()   # (Nx, Ny, Nz)

# 任意断面を自分で切る
mid_z = speed.shape[2] // 2
import matplotlib.pyplot as plt
plt.imshow(speed[:, :, mid_z].T, origin='lower')
plt.colorbar(label='speed [m/s]')
plt.show()
```

---

## 実験ログ

配置を変えながら複数回実行するとき、ログで比較できる。

```python
from ff_room.io import ExperimentLog

log = ExperimentLog("results/log.jsonl")

# 実行のたびにappend
hash1 = log.append(result1, notes="扇風機: 南壁側")
hash2 = log.append(result2, notes="扇風機: 中央")

# 全ログ読み込み
entries = log.load_all()
for e in entries:
    print(e["config_hash"], e["converged"], e["wall_time_s"], e["notes"])
```

ログはJSON Lines形式 (`.jsonl`)。設定の完全なスナップショットを含む。

---

## ディレクトリ構成

```
ff-room/
├── build.sh                        ← C++ビルドスクリプト
├── pyproject.toml                  ← pip installエントリポイント
│
├── cpp/
│   ├── CMakeLists.txt
│   ├── src/
│   │   ├── core/
│   │   │   ├── Grid.hpp/.cpp           ← MACスタガード格子
│   │   │   ├── PoissonSolver.hpp/.cpp  ← 共役勾配法
│   │   │   ├── BoundaryManager.hpp/.cpp← 境界条件
│   │   │   └── FluidSolver.hpp/.cpp    ← 射影法メインループ
│   │   └── bindings/
│   │       └── module.cpp              ← pybind11バインディング
│   └── tests/                          ← Catch2ユニットテスト
│
├── python/
│   ├── ff_room/
│   │   ├── config.py        ← RoomConfig, FanConfig, SceneConfig (YAML/JSON)
│   │   ├── scene.py         ← SceneConfig → C++オブジェクト変換
│   │   ├── solver_bridge.py ← SolverBridge → SimResult
│   │   ├── visualization.py ← Visualizer (PyVista)
│   │   └── io.py            ← ResultStore, ExperimentLog
│   └── tests/               ← pytestテスト
│
├── examples/
│   ├── basic_room.yaml      ← 単一扇風機サンプル
│   ├── two_fans.yaml        ← 向き合う扇風機サンプル
│   └── run_basic.py         ← サンプル実行スクリプト
│
└── docs/
    ├── model_assumptions.md ← 数値モデルの仮定・境界条件・参考文献
    ├── user_guide.md        ← ユーザーガイド (詳細)
    └── developer_guide.md   ← 開発者ガイド (拡張方法)
```

---

## テスト実行

### Python テスト (C++不要)

```bash
cd python
python -m pytest tests/ -v
```

### C++ テスト (Catch2が必要)

```bash
cd cpp/build
cmake .. -DBUILD_TESTS=ON
make && ./ffroom_tests
```

---

## ドキュメント一覧

| ファイル | 内容 |
|---------|------|
| `docs/model_assumptions.md` | 数値モデルの仮定・境界条件・単位系・参考文献 |
| `docs/user_guide.md`        | 使い方詳細・パラメータ選び方・トラブルシューティング |
| `docs/developer_guide.md`   | アーキテクチャ・拡張方法・「○○したい場合はここを直す」 |
