# 開発者ガイド

## アーキテクチャ概要

```
┌────────────────────────────────────────────────────┐
│ Python レイヤー                                     │
│                                                    │
│  SceneConfig (config.py)                          │
│       ↓ to_dict / from_dict                       │
│  Scene (scene.py)          ResultStore (io.py)    │
│       ↓ builds              ↑ saves                │
│  SolverBridge (solver_bridge.py)                  │
│       ↓ calls               ↑ returns SimResult   │
├────────────┬───────────────────────────────────────┤
│ pybind11   │  cpp/src/bindings/module.cpp          │
├────────────┴───────────────────────────────────────┤
│ C++ レイヤー                                        │
│                                                    │
│  Grid               FluidSolver    LBMSolver      │
│  (MAC格子, 非均一)   (射影法)        (D3Q19 BGK)    │
│       ↕                   ↕                       │
│  BoundaryManager    PoissonSolver                 │
│  (境界条件適用)      (CG + Multigrid V-cycle)      │
└────────────────────────────────────────────────────┘
```

**C++の責務**: 格子データ管理・数値計算・境界条件の適用  
**Pythonの責務**: 設定解析・C++呼び出し・可視化・I/O・実験管理

データの流れ:
1. `SceneConfig` (YAML) → `Scene` が C++ `Grid` + `BoundaryManager` を構築
2. `SolverBridge` が C++ `FluidSolver` を呼び出し、収束まで実行
3. 結果を `numpy` 配列にコピーして `SimResult` として返す
4. `Visualizer` / `ResultStore` が Python 側で処理

---

## C++コア詳解

### Grid (MAC格子)

```
u[i][j][k]: x方向速度 - x面中心 (i=0..Nx, j=0..Ny-1, k=0..Nz-1)
v[i][j][k]: y方向速度 - y面中心 (i=0..Nx-1, j=0..Ny, k=0..Nz-1)
w[i][j][k]: z方向速度 - z面中心 (i=0..Nx-1, j=0..Ny-1, k=0..Nz)
p[i][j][k]: 圧力      - セル中心 (i=0..Nx-1, j=0..Ny-1, k=0..Nz-1)
```

インデックスはFlat配列。変換マクロ: `g.u_idx(i,j,k)`, `g.c_idx(i,j,k)` など。

### FluidSolver (射影法)

1ステップの流れ:
```
compute_intermediate()          ← advection(UPWIND1/QUICK/LAX_WENDROFF) + diffusion(中心差分)
apply_noslip()                  ← 壁・障害物: u=0
apply_inflow()                  ← 扇風機: u=U_fan
[apply_buoyancy()]              ← (thermal) w* += g*β*(T_face - T_ref)*dt
build_poisson_rhs()             ← (ρ/dt) * ∇·u*
poisson_.solve()                ← CG法で ∇²p = rhs
apply_correction()              ← u -= (dt/ρ) ∇p
apply_noslip()                  ← 再適用
apply_inflow()                  ← 再適用
[advect_temperature()]          ← (thermal) UPWIND1/QUICK/LAX_WENDROFF + 拡散; T ← T_tmp
[apply_opening_temperature()]   ← (thermal) 開口部: 流入=T_outside, 流出=ゼロ勾配
```

`[]` は `thermal=true` のときのみ実行。  
移流スキームは `FluidSolverParams.advection` (デフォルト UPWIND1) で切り替え。  
OpenMP 並列化は `FluidSolverParams.use_openmp` (デフォルト false) で制御。

### PoissonSolver (CG法)

標準的な共役勾配法。初期推定=0。  
境界条件は `apply_laplacian()` 内でゴーストセル法により実装:
- SOLID隣接: ∂p/∂n = 0 (Neumann) → ゴーストセル = 内側と同値
- OUTFLOW隣接: p = 0 (Dirichlet) → ゴーストセル = -内側

---

## 「○○したい場合はここを直す」

### A. 扇風機を斜め方向に向けたい — 実装済み

**実装済み**: `scene.py` の `_add_fan()` が `direction_vector()` から 3成分速度を計算し、
`FanBC.vel = [vx*V, vy*V, vz*V]` として C++ `apply_inflow()` に渡す。
`apply_inflow()` は `FanBC.vel` の3成分を対応するface速度 (u/v/w) すべてに設定する。

使い方:

```yaml
fans:
  - direction: [45, 0]    # 斜め45° (X-Y平面)
  - direction: [0, 30]    # 水平+30°上向き
  - direction: [0, -90]   # 真下
```

`snapped_axis()` はBC面の決定 (扇風機パッチをどの壁に置くか) にのみ使用。
速度ベクトル自体はスナップしない。

---

### B. 温度・浮力 (Boussinesq近似) — 実装済み

温度シミュレーションは既に実装されている。YAML で有効にする:

```yaml
solver:
  T_initial: 32.0    # 初期室内温度 [°C]
  T_outside: 25.0    # 外気温 [°C]
  T_target:  30.0    # 目標温度 (nullで速度収束を使用)
  buoyancy:  true    # 熱対流を有効化
```

**実装されている方程式**:
```
∂T/∂t + (u·∇)T = α∇²T            (温度移流拡散)
f_z = g β (T_face - T_ref)        (Boussinesq浮力, +z = 上方向)
```

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/Grid.hpp/.cpp` | `T[]`, `T_tmp[]` フィールド追加 |
| `cpp/src/core/FluidSolver.cpp` | `apply_buoyancy()`, `advect_temperature()`, `compute_T_mean()` |
| `cpp/src/core/BoundaryManager.hpp/.cpp` | `OpeningBC` (T_outside付き), `apply_opening_temperature()` |
| `cpp/src/bindings/module.cpp` | `T_shaped()`, `OpeningBC`, thermal params のバインディング |
| `python/ff_room/config.py` | `OpeningConfig`, `SolverConfig.T_initial/T_outside/T_target/buoyancy` |
| `python/ff_room/scene.py` | `_add_opening()`, 開口部がある場合の自動アウトフロー無効化 |
| `python/ff_room/solver_bridge.py` | `FluidSolverParams` への thermal params 設定, `SimResult.T_field/T_mean` |

**温度境界条件の詳細** (`apply_opening_temperature()`):

```
面の法線方向速度 u_f を確認:
  u_f が外→内の方向 (流入) → T[セル] = T_outside  (外気温を適用)
  u_f が内→外の方向 (流出) → T[セル] = T[内側セル] (ゼロ勾配)
```

起動時 (u=0) はすべて流入扱い → 開口部セルが T_outside に初期化される。

**拡張: 扇風機に吹き出し温度を設定したい場合**

現状は扇風機の吹き出し空気温度は `T_initial` (室温) と同じ。
加熱/冷却付き扇風機 (エアコン) を模擬したい場合:

1. `BoundaryManager` の `FanBC` に `T_inflow` フィールドを追加
2. `apply_opening_temperature()` と同様のロジックで INFLOW セルに `T_inflow` を適用
3. `config.py` の `FanConfig` に `temperature` フィールドを追加

**難易度**: 易 (1〜2時間)

---

### C-0. アニメーション出力 (GIF/MP4) — 実装済み

**実装済み** (TASK-01): `SolverBridge.run_animated()` でフィールドスナップショットを収集し、
`Animator` クラスで GIF/MP4 に書き出せる。

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `python/ff_room/solver_bridge.py` | `FieldSnapshot` データクラス, `SolverBridge.run_animated(save_every, progress_callback, print_interval)` |
| `python/ff_room/animation.py` | `Animator(snapshots, config)` クラス — `save_gif(path, fps, dpi)`, `save_mp4(path, fps, dpi)` |
| `python/ff_room/visualization.py` | `Visualizer.save_animation(snapshots, path, fps)` ラッパー |
| `python/ff_room/__init__.py` | `FieldSnapshot`, `Animator` をエクスポート |

**`FieldSnapshot` フィールド**:

```python
@dataclass
class FieldSnapshot:
    step:    int         # ソルバーステップ番号
    time_s:  float       # 物理時間 = step * dt [s]
    T_mean:  float       # 室内平均温度 [°C]
    speed:   np.ndarray  # (Nx, Ny, Nz) 速度スカラー [m/s]
    T_field: np.ndarray  # (Nx, Ny, Nz) 温度場 [°C]; 非熱計算時はゼロ
    vel:     np.ndarray  # (Nx, Ny, Nz, 3) セル中心速度ベクトル [m/s]
```

**`Animator` のレンダリング仕様**:
- 左パネル: 上面図 (XY スライス, z=Nz//2) — 速度カラーマップ + 方向矢印
- 右パネル: 側面図 (YZ スライス, x=Nx//2) — 非熱計算: 速度, 熱計算: 温度 (RdBu_r)
- カラーバー範囲: 全フレーム共通 (99パーセンタイルでクリップ)
- `save_mp4()` は ffmpeg がない環境では GIF に自動フォールバック

---

### C. 移流スキームを高精度化したい — 実装済み

**実装済み**: `compute_intermediate()` と `advect_temperature()` で QUICK/Lax-Wendroff が選択可能。
設定: `solver.advection` フィールドで直接指定。

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/FluidSolver.hpp` | `enum class AdvectionScheme { UPWIND1, QUICK, LAX_WENDROFF }`, `FluidSolverParams.advection` フィールド |
| `cpp/src/core/FluidSolver.cpp` | `quick_adv()`, `lax_wendroff_adv()` ヘルパー, `compute_intermediate()` / `advect_temperature()` 内のスキームディスパッチ |
| `python/ff_room/config.py` | `SolverConfig.advection: str = "upwind"` |
| `python/ff_room/solver_bridge.py` | `_apply_solver_flags()` が `sc.advection` → `params.advection` に変換 |

**設定** (YAML または Python):

```yaml
solver:
  advection: quick          # upwind | quick | lax_wendroff
```

```python
config.solver.advection = "quick"
result = SolverBridge(config).run()
```

**選択肢** (実装済み):

#### QUICK (3次精度 upwind-biased)
```cpp
static inline double quick_adv(double phi_mm, double phi_m, double phi_c,
                                double phi_p, double uc) {
    if (uc > 0)
        return (3*phi_c + 6*phi_m - phi_mm) / 8.0;  // upwind-biased
    else
        return (3*phi_m + 6*phi_c - phi_p) / 8.0;
}
```

#### Lax-Wendroff (2次精度, 条件付き安定)
```cpp
static inline double lax_wendroff_adv(double phi_m, double phi_c, double phi_p,
                                       double uc, double h, double dt) {
    double flux_r = 0.5*uc*(phi_p+phi_c) - 0.5*uc*uc*dt/h*(phi_p-phi_c);
    double flux_l = 0.5*uc*(phi_c+phi_m) - 0.5*uc*uc*dt/h*(phi_c-phi_m);
    return (flux_r - flux_l) / h;
}
```

境界セル (i=1, i=Nx-1 など) では1次にフォールバック。

---

### D. 圧力ソルバーを高速化したい — マルチグリッド実装済み

**実装済み**: `poisson_method: multigrid` で V-cycle Gauss-Seidel ソルバーに切り替えられる。
25×20×12格子で **CG比約5倍高速** (2.25s → 0.44s)。

**設定**:

```yaml
solver:
  poisson_method: multigrid   # cg (デフォルト) | multigrid
```

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/PoissonSolver.hpp` | `MultigridPoissonSolver` クラス宣言 |
| `cpp/src/core/PoissonSolver.cpp` | V-cycle: `smooth()` (Gauss-Seidel), `restrict_residual()` (injection), `prolongate()` (nearest-neighbor), `compute_residual()` |
| `cpp/src/core/FluidSolver.hpp` | `FluidSolverParams.use_multigrid`, `multigrid_poisson_` メンバ |
| `cpp/src/core/FluidSolver.cpp` | step()内: `use_multigrid ? multigrid_poisson_.solve() : poisson_.solve()` |
| `python/ff_room/config.py` | `SolverConfig.poisson_method: str = "cg"` |
| `python/ff_room/solver_bridge.py` | `params.use_multigrid = (sc.poisson_method == 'multigrid')` |

**残る選択肢** (将来): FFTベース O(N log N) — FFTW3が必要、全境界Dirichlet専用。

---

### E. 並列化したい (OpenMP) — 実装済み

**実装済み**: `solver.parallel: true` で移流・拡散ループが OpenMP 並列実行される。

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/FluidSolver.cpp` | `#pragma omp parallel for collapse(3) schedule(static) if(params_.use_openmp)` を各主要ループに追加 |
| `cpp/CMakeLists.txt` | `find_package(OpenMP)` (グレースフルフォールバック付き), `FFROOM_HAS_OPENMP` define, `OpenMP::OpenMP_CXX` リンク |
| `python/ff_room/config.py` | `SolverConfig.parallel: bool = False` |
| `python/ff_room/solver_bridge.py` | `_apply_solver_flags()`: `params.use_openmp = bool(sc.parallel)` |

**設定**:

```yaml
solver:
  parallel: true   # false がデフォルト
```

OpenMP が利用できない環境では自動的にシングルスレッドにフォールバックする。  
Poissonソルバー (CG/マルチグリッド) はデータ依存があるため並列化対象外。

CMake ビルドログで確認:
```
OpenMP found — enabling parallelization   # 有効
OpenMP not found — building without parallelization  # フォールバック
```

---

### F. GUIを追加したい

`SolverBridge` はコールバックを受け取れるので、GUI側で進捗バーを更新できる。

#### tkinterの場合

```python
import tkinter as tk
from tkinter import ttk
import threading
from ff_room import SceneConfig, SolverBridge, Visualizer

class SimGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.progress = ttk.Progressbar(self.root, maximum=500)
        self.progress.pack()
        self.btn = tk.Button(self.root, text="Run", command=self.run)
        self.btn.pack()

    def run(self):
        config = SceneConfig.load_yaml("examples/basic_room.yaml")
        def cb(step, vel, div):
            self.progress['value'] = step
            self.root.update()
        def _run():
            result = SolverBridge(config).run(progress_callback=cb)
            Visualizer(result).multi_panel()
        threading.Thread(target=_run, daemon=True).start()
        self.root.mainloop()
```

#### Jupyter Notebookの場合

```python
from ipywidgets import interact, FloatSlider, Button, Output
import ipywidgets as widgets

out = Output()
run_btn = widgets.Button(description="Run Simulation")
fan_vel = FloatSlider(min=1.0, max=6.0, step=0.5, value=3.0, description="Fan m/s")
fan_x   = FloatSlider(min=0.1, max=4.9, step=0.1, value=0.5, description="Fan x [m]")

def on_run(b):
    with out:
        out.clear_output()
        cfg = SceneConfig.load_yaml("examples/basic_room.yaml")
        cfg.fans[0].velocity  = fan_vel.value
        cfg.fans[0].position  = (fan_x.value, 2.0, 1.2)
        r = SolverBridge(cfg).run(print_interval=0)
        # Jupyter内でインライン可視化
        viz = Visualizer(r)
        pl = viz.slice_xz(show=False)
        pl.show(jupyter_backend='trame')

run_btn.on_click(on_run)
display(widgets.VBox([fan_vel, fan_x, run_btn, out]))
```

---

### G. 非均一格子 (壁面付近を細かくしたい) — 実装済み

**実装済み**: `room.stretch > 1.0` で幾何学的格子伸張が有効になる。
Nx/Ny/Nz は変わらないため計算量は増えない。

**設定**:

```yaml
room:
  stretch: 1.1    # 壁面付近を軽く細分化 (推奨: 1.05〜1.2)
  # stretch: 1.3  # 強い細分化
```

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/Grid.hpp/.cpp` | `dxv[Nx]`, `dyv[Ny]`, `dzv[Nz]` (可変セルサイズ), `xs[Nx+1]` 等の面座標, `set_face_coords()` |
| `cpp/src/core/FluidSolver.cpp` | 移流・拡散・Poisson RHS・補正・発散計算すべてに `dxv[i]`/`dyv[j]`/`dzv[k]` を使用 |
| `cpp/src/core/PoissonSolver.cpp` | 非均一Laplacianステンシル: `scale = 1/(h_face * dx_cell)` |
| `python/ff_room/config.py` | `RoomConfig.stretch: float = 1.0` |
| `python/ff_room/scene.py` | `_stretched_faces(N, L, r)` — 対称幾何学伸張, `Scene._build()` で `set_face_coords()` 呼び出し |

**注意**: stretch格子では圧力ソルバーの条件数が悪化するため、
`poisson_method: multigrid` との組み合わせを推奨。

---

### H. 出力フォーマットを追加したい

`python/ff_room/io.py` に新メソッドを追加するだけ。

例: HDF5形式

```python
# io.py に追加
@staticmethod
def save_hdf5(result: SimResult, path: str) -> None:
    import h5py
    with h5py.File(path, 'w') as f:
        f.create_dataset('u', data=result.u_field, compression='gzip')
        f.create_dataset('v', data=result.v_field, compression='gzip')
        f.create_dataset('w', data=result.w_field, compression='gzip')
        f.create_dataset('p', data=result.p_field, compression='gzip')
        f.attrs['config'] = json.dumps(result.config.to_dict())
```

---

### I. LBM (Lattice Boltzmann Method) に切り替えたい — 実装済み

**実装済み**: `solver.method: lbm` で D3Q19 BGK LBM ソルバーに切り替えられる。

**射影法との比較**:

| 項目 | 射影法 | LBM D3Q19 |
|------|--------|-----------|
| 並列化 | OpenMP | 非常に容易 (セル独立) |
| 圧縮性 | 非圧縮 | 弱圧縮性 (div_max高め) |
| 温度計算 | 対応 | 非対応 (等温のみ) |
| 境界条件 | 自然 | bounce-back (簡潔) |
| メモリ | 低 | 高 (19分布関数×格子) |
| 収束速度 | ~529 steps | ~119 steps (basic_room) |

**設定**:

```yaml
solver:
  method: lbm
  tau: 0.8        # 緩和時間 (0.5 < τ < 2.0 が安定)
  max_steps: 500  # LBMは射影法より少ないstepで収束することが多い
```

**実装の場所**:

| ファイル | 役割 |
|---------|------|
| `cpp/src/core/LBMSolver.hpp` | クラス宣言, D3Q19格子定数, 物理↔格子単位変換 |
| `cpp/src/core/LBMSolver.cpp` | BGK衝突, streaming (double-buffer), bounce-back BC, 流入eq初期化, 流出ゼロ勾配BC, マクロ変数抽出 |
| `cpp/CMakeLists.txt` | `LBMSolver.cpp` を `CORE_SOURCES` に追加 |
| `cpp/src/bindings/module.cpp` | `LBMSolver` pybind11バインディング |
| `python/ff_room/config.py` | `SolverConfig.method = "lbm"`, `SolverConfig.tau = 0.8` |
| `python/ff_room/solver_bridge.py` | `sc.method == 'lbm'` → `LBMSolver(g, bm, sc.tau, ...)` |

**注意**: LBMは `div_max` が射影法より高くなる (弱圧縮性のため)。
`T_field` は zeros、`buoyancy` は無効。温度シミュレーションには射影法を使うこと。

---

## テスト戦略

### Python テスト

```bash
cd python && python -m pytest tests/ -v
```

- `test_config.py`: 設定シリアライズ、FanConfig角度計算、YAML/JSON往復
- `test_postprocess.py`: SimResultのNumPy操作、CSV書き出し、npz往復

C++不要。任意の環境で走る。

### C++ テスト (Catch2)

```bash
cd cpp/build && cmake .. -DBUILD_TESTS=ON && make && ./ffroom_tests
```

- `test_grid.cpp`: Grid構築、インデックス一意性、ゼロ化
- `test_poisson.cpp`: 製造解でCGの収束と精度を検証
- `test_couette.cpp`: Couette流の定常プロファイルと解析解を比較

### ベンチマーク追加方法

`cpp/tests/test_lid_driven_cavity.cpp` を新規作成し、Ghia (1982) の値と比較:

```cpp
// Re=100の蓋駆動空洞流
// 解析: Ghia et al. (1982), Table 1
// y=0.5 断面の u 速度プロファイルを比較
TEST_CASE("Lid-driven cavity Re=100", "[benchmark][cavity]") {
    // 2D-likeドメイン: Ny=2の薄いスライス
    Grid g(64, 2, 64, 1.0, 1.0/32, 1.0);
    // ...
    REQUIRE(max_err_u < 0.05);  // Ghia値との最大誤差 5%以内
}
```

---

## 新機能を追加するときの手順

1. **C++変更が必要か判断**
   - 数値計算・格子・境界条件 → C++
   - 設定・可視化・I/O → Python のみで可

2. **C++変更の場合**
   - `cpp/src/core/` で実装
   - `cpp/tests/` にテスト追加
   - `cpp/src/bindings/module.cpp` に pybind11バインディング追加
   - `bash build.sh` でビルド確認

3. **Python側の変更**
   - `config.py`: パラメータ追加
   - `scene.py`: C++オブジェクトへの変換
   - `solver_bridge.py`: 呼び出し追加
   - `python/tests/` にテスト追加

4. **ドキュメント更新**
   - `docs/model_assumptions.md`: 仮定・BCの変更を記載
   - `docs/user_guide.md`: 新パラメータの使い方
   - `examples/` にサンプル設定追加

---

## よくある落とし穴

### MAC格子のインデックス範囲
u は `(Nx+1)*Ny*Nz` 要素。`u_idx(Nx, j, k)` が有効 (東壁面)。
`p_idx(Nx, j, k)` はアウトオブバウンズ (`Nx-1` まで)。

### 境界条件の適用順序
`apply_noslip()` の後に `apply_inflow()` を呼ぶこと。  
逆にすると扇風機面が no-slip で上書きされる。

### pybind11のゼロコピービュー
`g.u_shaped()` が返す numpy 配列は C++ の `Grid.u` のメモリを直接参照する。
`SimResult` に格納するときは `np.array(g.u_shaped())` でコピーを取ること
(`solver_bridge.py` 参照)。Grid が解放された後も使えるようにするため。

### 収束判定は最大値ノルム
`convergence_tol` は `max |u^{n+1} - u^n| / dt`。  
RMS (二乗平均平方根) ではないので、局所的な振動があると収束しにくい。
その場合は `dt` を小さくするか、局所的に不安定なBCがないか確認する。
