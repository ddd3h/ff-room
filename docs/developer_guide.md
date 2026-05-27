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
│  Grid               FluidSolver                   │
│  (MAC格子)           (射影法メインループ)             │
│       ↕                   ↕                       │
│  BoundaryManager    PoissonSolver                 │
│  (境界条件適用)      (CG法)                         │
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
compute_intermediate()          ← advection(upwind1) + diffusion(中心差分)
apply_noslip()                  ← 壁・障害物: u=0
apply_inflow()                  ← 扇風機: u=U_fan
[apply_buoyancy()]              ← (thermal) w* += g*β*(T_face - T_ref)*dt
build_poisson_rhs()             ← (ρ/dt) * ∇·u*
poisson_.solve()                ← CG法で ∇²p = rhs
apply_correction()              ← u -= (dt/ρ) ∇p
apply_noslip()                  ← 再適用
apply_inflow()                  ← 再適用
[advect_temperature()]          ← (thermal) upwind + 拡散; T ← T_tmp
[apply_opening_temperature()]   ← (thermal) 開口部: 流入=T_outside, 流出=ゼロ勾配
```

`[]` は `thermal=true` のときのみ実行。

### PoissonSolver (CG法)

標準的な共役勾配法。初期推定=0。  
境界条件は `apply_laplacian()` 内でゴーストセル法により実装:
- SOLID隣接: ∂p/∂n = 0 (Neumann) → ゴーストセル = 内側と同値
- OUTFLOW隣接: p = 0 (Dirichlet) → ゴーストセル = -内側

---

## 「○○したい場合はここを直す」

### A. 扇風機を斜め方向に向けたい

**現状**: `scene.py` の `_add_fan()` が `snapped_axis()` で最近軸にスナップ。

**必要な変更**:

1. `BoundaryManager` に新しいBC型を追加 (`cpp/src/core/BoundaryManager.hpp`):
   ```cpp
   struct ObliqueInflowBC {
       // セルインデックスリスト + 速度ベクトル
       std::vector<std::array<int,3>> cells;
       std::array<double,3> vel;
   };
   ```

2. `FluidSolver::compute_intermediate()` でINFLOW近傍の速度を補間処理。  
   斜め流入は速度の分解 (u_x, u_y, u_z) を各faceに分配する。
   具体的には、INFLOWセルの隣接面すべてに `vel[ax] * face_normal[ax]` を設定。

3. `scene.py` の `_add_fan()` で `snapped_axis()` 分岐を削除し、
   3成分すべて設定するパスに変更。

4. バインディングに `ObliqueInflowBC` を追加 (`module.cpp`)。

**難易度**: 中 (2〜3時間)

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

### C. 移流スキームを高精度化したい (2次精度)

**現状**: `compute_intermediate()` で `upwind1()` (1次upwind)。数値拡散が強い。

**選択肢**:

#### Lax-Wendroff (2次精度, 条件付き安定)
```cpp
// φの1次元移流
static inline double lax_wendroff(double phi_m, double phi_c, double phi_p,
                                   double uc, double h, double dt) {
    double flux_r = 0.5*uc*(phi_p+phi_c) - 0.5*uc*uc*dt/h*(phi_p-phi_c);
    double flux_l = 0.5*uc*(phi_c+phi_m) - 0.5*uc*uc*dt/h*(phi_c-phi_m);
    return (flux_r - flux_l) / h;
}
```

#### QUICK (3次精度 upwind-biased)
```cpp
static inline double quick(double phi_mm, double phi_m, double phi_c,
                            double phi_p, double uc) {
    if (uc > 0)
        return (3*phi_c + 6*phi_m - phi_mm) / 8.0;  // upwind-biased
    else
        return (3*phi_m + 6*phi_c - phi_p) / 8.0;
}
```

`compute_intermediate()` 内の `upwind1(...)` を置き換えるだけ。
境界セル (i=1, i=Nx-1 など) では1次にフォールバックが必要。

**`FluidSolverParams` に追加**:
```cpp
enum class AdvectionScheme { UPWIND1, LAX_WENDROFF, QUICK };
AdvectionScheme advection = AdvectionScheme::UPWIND1;
```

**難易度**: 易〜中 (1〜3時間)

---

### D. 圧力ソルバーを高速化したい

**現状**: CG法 O(N^{1.5}) 程度。N=25×20×12≈6000なら十分速い。  
格子を大きくする (100×80×50≈400,000セル) と数十秒かかる。

**選択肢**:

#### FFTベース (直交格子+均一BC専用)
全境界がDirichletなら `FFTW3` ライブラリを使い O(N log N) で解ける。
`PoissonSolver` を差し替えるだけ。

```cpp
// PoissonSolver の派生クラス or 別実装
class FFTPoissonSolver {
    fftw_plan plan_r2c, plan_c2r;
public:
    double solve(Grid& grid, const std::vector<double>& rhs);
};
```

`CMakeLists.txt` に `find_package(FFTW3)` を追加。

#### マルチグリッド (境界条件非依存)
収束が O(N) 。実装は重いが長期的に推奨。
既存ライブラリ: Hypre, AMGcl, もしくは Gauss-Seidel + V-cycle を手書き。

**難易度**: FFT=中, マルチグリッド=難

---

### E. 並列化したい (OpenMP)

`FluidSolver.cpp` の内側ループにOpenMP追加:

```cpp
// FluidSolver.cpp - compute_intermediate() の u* 計算
#pragma omp parallel for collapse(3) schedule(static)
for (int i = 1; i < grid_.Nx; i++)
for (int j = 0; j < grid_.Ny; j++)
for (int k = 0; k < grid_.Nz; k++) {
    ...
}
```

`CMakeLists.txt` に追加:
```cmake
find_package(OpenMP REQUIRED)
target_link_libraries(_ffroom_core PRIVATE OpenMP::OpenMP_CXX)
```

Poissonのコードはデータ依存があるためCGのdot product部分のみ並列化可能。

**難易度**: 易 (1時間)

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

### G. 非均一格子 (壁面付近を細かくしたい)

現状は均一格子。壁面境界層を解像したい場合は格子を引き伸ばす。

**必要な変更** (`Grid.hpp`):

```cpp
// 均一: dx = Lx/Nx
// 非均一: x[i] を任意の座標列に変更
std::vector<double> x_face;  // x-面の座標列 (長さ Nx+1)
std::vector<double> y_face;
std::vector<double> z_face;
std::vector<double> dx_cell; // dx[i] = x_face[i+1] - x_face[i]
```

差分演算子すべてを `dx` → `dx_cell[i]` に変更する (大規模な変更)。

**代替案**: 格子を細かくしてセル数を増やす (実装変更なしで可)。

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

### I. LBM (Lattice Boltzmann Method) に切り替えたい

**現状の射影法との比較**:

| 項目 | 射影法 (現状) | LBM D3Q19 |
|------|------------|-----------|
| 実装難易度 | 中 | 中 |
| 並列化 | OpenMP容易 | 非常に容易 |
| 圧縮性 | 非圧縮 | 弱圧縮性 |
| 境界条件 | 自然 | bounce-back (簡潔) |
| 複雑形状 | やや複雑 | 容易 |
| メモリ | 低 | 高 (19分布関数) |

LBMに切り替える場合: `cpp/src/core/` に `LBMSolver.hpp/.cpp` を新規作成し、
`FluidSolver` と同じインターフェース (`step()`, `run()`) を持たせる。
`SolverBridge` 側でソルバー種別を切り替えるだけで Python 側の変更なし。

```yaml
solver:
  method: lbm   # "projection" (default) | "lbm"
  tau: 0.8      # LBM緩和時間 (0.5 < τ < 2.0 が安定)
```

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
