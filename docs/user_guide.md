# ユーザーガイド

## 1. インストールと初回実行

### 必要なもの

```
Python >= 3.9
pip install pybind11 pyyaml numpy pyvista
GCC/Clang (C++17対応), CMake >= 3.15
```

### ビルド

```bash
bash build.sh
```

`python/ff_room/_ffroom_core*.so` が生成されれば成功。

### 動作確認

```bash
python examples/run_basic.py
```

収束ログが流れ、`results/basic_room_overview.png` が生成されれば OK。

---

## 2. 設定ファイルの作り方

### 2-1. 部屋サイズと格子解像度

```yaml
room:
  size: [5.0, 4.0, 2.5]   # Lx, Ly, Lz [m]
  grid: [25, 20, 12]       # Nx, Ny, Nz
```

**格子解像度の選び方**

| 格子 | セルサイズ (5×4×2.5m部屋) | 計算時間目安 | 用途 |
|------|--------------------------|-------------|------|
| 15×12×8 | ~0.2m | 数秒 | 動作確認 |
| 25×20×12 | ~0.2m | 数十秒 | 標準 |
| 50×40×25 | ~0.1m | 数分 | 精度重視 |
| 100×80×50 | ~0.05m | 数十分 | 研究用詳細 |

セルサイズはできれば各軸で均等にする (dx ≈ dy ≈ dz)。
格子数を増やすとメモリは O(N³)、計算は O(N³) × 反復数 で増える。

**CFL安定条件**  
`dt` が大きすぎると発散する。目安:

```
dt < dx / U_max     (移流安定性、CFL≦1)
```

U_max = 扇風機風速 [m/s]、dx = 最小セルサイズ [m]。  
例: U_max=3 m/s, dx=0.2m → dt < 0.067 s。dt=0.01 s で十分余裕あり。

---

### 2-2. 扇風機の設定

```yaml
fans:
  - name: fan0
    position: [0.5, 2.0, 1.2]   # 扇風機フェイス中心 [m]
    direction: [0.0, 0.0]        # [方位角, 仰角] deg
    velocity: 3.0                 # 吹き出し速度 [m/s]
    radius: 0.15                  # 半径 [m]
```

**directionの指定方法**

方位角: 0=+x, 90=+y, 180=-x, 270=-y  
仰角: 0=水平, 90=真上, -90=真下

```
direction: [0, 0]    → 真東に吹く (扇風機が西壁に設置)
direction: [90, 0]   → 真北に吹く (扇風機が南壁に設置)
direction: [180, 0]  → 真西に吹く (扇風機が東壁に設置)
direction: [0, 30]   → 東に向けて少し上に吹く
direction: [0, -90]  → 床に向けて真下に吹く
```

> **注意 (MVP)**: 方向は最近接の軸にスナップされる。例えば `[45, 0]` は
> `[0, 0]` (x軸) と同等になる。斜め方向対応は developer_guide.md 参照。

**positionの注意点**  
- 壁際に置く場合: position の x 座標を `radius` より大きく設定する
- 壁面に接する場合: 例えば西壁 (x=0) に設置するなら `position[0] = radius + dx/2`
- 格子の外に出ないこと (0〜Lx, 0〜Ly, 0〜Lz の範囲内)

**velocityの目安**

| 扇風機強さ | velocity [m/s] |
|-----------|---------------|
| 弱風 | 1.0〜2.0 |
| 中風 | 2.5〜4.0 |
| 強風 | 4.0〜6.0 |
| 業務用 | 6.0〜 |

---

### 2-3. 障害物 (家具など)

```yaml
obstacles:
  - name: sofa
    bbox_min: [2.0, 0.5, 0.0]   # 直方体の最小隅 [m]
    bbox_max: [3.0, 1.8, 0.9]   # 直方体の最大隅 [m]
```

床から浮かせた障害物も可:
```yaml
  - name: table
    bbox_min: [1.0, 1.0, 0.7]   # テーブル面の高さ0.7mから
    bbox_max: [2.0, 2.5, 0.72]  # 板厚2cm
```

障害物を床に接触させる場合は `bbox_min[2] = 0.0`。

---

### 2-4. 窓・ドアの設定

窓やドアを開けると「開口部」になる。開口部には2つの効果がある:
1. **気圧境界条件**: 圧力が外気と等しくなり空気が出入りできる
2. **温度交換**: 外気が流入する面では外気温を適用、流出する面では室内温度をそのまま使用

```yaml
openings:
  - name: window_north
    wall: north            # 壁の名前 (west/east/south/north/floor/ceiling)
    center: [3.0, 1.2]    # 壁面上の中心座標 [m]
                           #   east/west → (y_center, z_center)
                           #   south/north → (x_center, z_center)
                           #   floor/ceiling → (x_center, y_center)
    size: [0.9, 1.1]      # 幅 × 高さ [m]
    T_outside: 25.0        # この開口部の外気温 (省略時はsolver.T_outsideを使用)
  - name: door_south
    wall: south
    center: [1.5, 1.0]
    size: [0.8, 2.0]
```

**壁名の対応 (x=東西, y=南北, z=上下)**

| 壁名 | 面 | 典型的な用途 |
|------|------|------------|
| `west` | x=0 | 西壁の窓 |
| `east` | x=Lx | 東壁の窓 |
| `south` | y=0 | 南向き入口ドア |
| `north` | y=Ly | 北壁の窓 |
| `floor` | z=0 | (通常は使わない) |
| `ceiling` | z=Lz | (通常は使わない) |

**換気のための2開口部設定 (クロスベンチレーション)**

1開口部だけでは外気の取り込みが難しい。典型的な一人暮らしワークフロー:

```yaml
openings:
  - name: window          # 窓: 涼しい外気が入る側
    wall: north
    center: [cx, 1.2]
    size: [0.9, 1.1]
  - name: door_crack      # ドア少し開け: 熱い室内空気が出る側
    wall: south
    center: [1.5, 1.0]
    size: [0.8, 2.0]      # ドア全開なら 0.8×2.0m
```

> **注意**: 開口部を設定すると内部の自動アウトフロー(圧力調整用)は無効になる。
> 最低2つ以上の開口部を設定すると空気の流入・流出が成立する。

---

### 2-5. 温度シミュレーション設定

夏場の「室内が暑い、外が涼しい」ような温度差シナリオを扱える。

```yaml
solver:
  dt: 0.05              # 温度計算では少し大きなdtでも安定 (CFL ≦ 1)
  max_steps: 5000
  # ---熱流体設定---
  T_initial: 32.0       # 初期室内温度 [°C]
  T_outside: 25.0       # 外気温 [°C] (全開口部のデフォルト)
  T_target:  28.0       # この温度に達したら計算終了 [°C]
                        # 省略時 (null) は速度収束を使用
  buoyancy:  true       # 熱対流 (暖かい空気が上に集まる現象)
```

**T_targetの方向判定**

- `T_initial > T_target` → 冷却シナリオ。T_mean ≦ T_targetで終了
- `T_initial < T_target` → 加熱シナリオ。T_mean ≧ T_targetで終了

**buoyancyの効果**

`buoyancy: true` にすると暖かい空気(密度小)が上昇する自然対流が計算に含まれる。  
窓の上半分から熱い空気が出て、下半分から冷たい外気が入る「スタック効果」が再現される。

**dtの選び方 (温度計算あり)**

温度計算は流速計算と同じdtを使う。CFL ≦ 1 が安定条件:

```
dt < dx / U_max   (移流安定性)
```

温度計算自体は熱拡散が遅いため、より大きなdtも安定。  
例: U_max=3 m/s, dx=0.2m → dt < 0.067s → `dt: 0.05` で十分。

**冷却に必要な計算時間の目安**

時定数 τ ≈ 部屋体積 / 開口部流量 [s]

| 設定 | τ目安 | T目標達成に必要なmax_steps (dt=0.05) |
|------|-------|--------------------------------------|
| 小窓1つ (0.5×0.8m) | ~2000s | 難しい (流量不足) |
| 普通の窓+ドア | ~500s | 10000〜20000 |
| 大きい窓+ドア全開 | ~200s | 4000〜8000 |

> **参考**: 6×4.5×2.5m の部屋 (67.5m³) で 窓0.9m×1.1m + ドア0.8m×2.0m の場合、
> 実効時定数は約1000秒。250秒(max_steps=5000, dt=0.05)で2°C程度の冷却。

---

### 2-6. ソルバーパラメータ

```yaml
solver:
  dt: 0.01              # タイムステップ [s]
  max_steps: 1000       # 最大反復回数
  convergence_tol: 1.0  # 収束判定閾値 (速度変化量/dt)
  rho: 1.2              # 空気密度 [kg/m³]
  nu: 1.5e-5            # 動粘度 [m²/s]
```

**max_stepsの目安**

| シナリオ | 推奨max_steps | 推奨dt |
|---------|--------------|--------|
| 速度場のみ (等温) | 1000 | 0.01 |
| 温度あり、短時間 | 3000〜5000 | 0.05 |
| 温度あり、長時間 | 10000〜20000 | 0.05 |

**convergence_tolの目安**

`vel_change = max|Δu|/dt` [m/s²]。1.0 ≈ 0.05 m/s/step の変化。

- `1.0`: 標準 (流れが発達すれば十分)
- `0.1`: 厳密

**ρ, ν の変更**

| 条件 | ρ [kg/m³] | ν [m²/s] |
|------|-----------|---------|
| 20°C, 1atm (標準) | 1.204 | 1.516×10⁻⁵ |
| 10°C, 1atm | 1.247 | 1.418×10⁻⁵ |
| 30°C, 1atm | 1.165 | 1.608×10⁻⁵ |

---

## 2-7. ソルバーオプション設定

各オプションは独立したフィールドで制御できる。組み合わせも自由。

### ソルバーアルゴリズム (`method`)

```yaml
solver:
  method: projection   # デフォルト: 射影法 (MAC格子, CG/マルチグリッド圧力ソルバー)
  # method: lbm        # D3Q19 BGK Lattice Boltzmann法 (等温, bounce-back境界)
```

| method | 特徴 | 用途 |
|--------|------|------|
| `projection` | 非圧縮, 温度計算対応, Poissonソルバー | 標準・推奨 |
| `lbm` | 弱圧縮性, 等温のみ, 並列化容易 | 複雑形状・高速プロトタイプ |

LBM使用時は `tau` (緩和時間) も設定可能 (0.5 < τ < 2.0):
```yaml
solver:
  method: lbm
  tau: 0.8   # デフォルト。τ=1/√3 + 0.5 ≈ 1.08 で安定
  max_steps: 500
```

---

### 移流スキーム (`advection`)

```yaml
solver:
  advection: upwind        # デフォルト: 1次精度, 安定, 数値拡散あり
  # advection: quick       # 3次精度 QUICK: 数値拡散を低減 (推奨)
  # advection: lax_wendroff # 2次精度: 境界付近で1次にフォールバック
```

| advection | 精度 | 計算コスト | 備考 |
|-----------|------|-----------|------|
| `upwind` | 1次 | 基準 | 安定, 数値拡散大 |
| `quick` | 3次 | ≒同等 | 精度/コスト最良 |
| `lax_wendroff` | 2次 | ≒同等 | 条件付き安定 |

(`method: lbm` 時は無視)

---

### OpenMP並列化 (`parallel`)

```yaml
solver:
  parallel: false   # デフォルト: シングルスレッド
  # parallel: true  # マルチコアで高速化 (移流・拡散ループ)
```

- `true` にするとマルチコアCPUで移流・拡散ループが並列実行される
- ビルド時に OpenMP が見つからない場合はシングルスレッドにフォールバック
- Poisson ソルバー (CG/マルチグリッド) はデータ依存のため対象外
- (LBM時は無視)

---

### 圧力ソルバー (`poisson_method`)

```yaml
solver:
  poisson_method: cg          # デフォルト: 共役勾配法 O(N^1.5)
  # poisson_method: multigrid  # マルチグリッドV-cycle O(N) — 約5倍高速
```

| poisson_method | アルゴリズム | 計算量 | 用途 |
|----------------|------------|--------|------|
| `cg` | 共役勾配法 | O(N^1.5) | 標準・安定 |
| `multigrid` | V-cycle Gauss-Seidel | O(N) | 高解像度格子で有効 |

(LBM時は無視 — LBMに圧力Poissonソルバーなし)

---

### 格子細分化 (`room.stretch`)

```yaml
room:
  stretch: 1.0    # デフォルト: 均一格子
  # stretch: 1.1  # 壁面付近を軽く細分化
  # stretch: 1.3  # 壁面付近を強く細分化 (境界層解像向け)
```

壁面付近のセルが細かくなり、室内中央のセルは相対的に粗くなる幾何学的格子伸張。
Nx, Ny, Nz は変わらないため計算量は増えない。

---

### 組み合わせ例

```yaml
# 精度重視 + 並列化 (推奨バランス)
solver:
  method: projection
  advection: quick
  parallel: true
  poisson_method: multigrid

# 高速プロトタイプ (LBM)
solver:
  method: lbm
  tau: 0.8
  max_steps: 300

# 壁面境界層を解像したい
solver:
  advection: quick
  poisson_method: multigrid   # stretch格子で圧力収束を助ける
room:
  stretch: 1.15
```

Pythonから直接設定:

```python
from ff_room import SceneConfig, SolverBridge

config = SceneConfig.load_yaml("my_room.yaml")

# QUICK + OpenMP + multigrid
config.solver.advection      = "quick"
config.solver.parallel       = True
config.solver.poisson_method = "multigrid"

# または LBM
config.solver.method = "lbm"
config.solver.tau    = 0.8

result = SolverBridge(config).run()
```

---

## 3. シミュレーション実行

### Pythonスクリプトから

```python
from ff_room import SceneConfig, SolverBridge

config = SceneConfig.load_yaml("my_room.yaml")

# print_interval=0 で進捗非表示
result = SolverBridge(config).run(print_interval=100)

print(f"収束: {result.converged}")
print(f"ステップ数: {result.steps}")
print(f"計算時間: {result.wall_time_s:.1f} s")
print(f"最大発散残差: {result.divergence_max:.2e}")
```

### コールバックで進捗を取得

温度シミュレーションの場合、コールバックに `T_mean` も渡される:

```python
history = []

def cb(step, vel_change, div_max, T_mean):
    history.append((step, vel_change, div_max, T_mean))

result = SolverBridge(config).run(progress_callback=cb, print_interval=0)

# 温度冷却カーブをプロット
import matplotlib.pyplot as plt
steps, changes, divs, temps = zip(*history)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.semilogy(steps, changes, label='vel_change')
ax1.semilogy(steps, divs, label='div_max')
ax1.set_xlabel('step'); ax1.legend()
ax2.plot(steps, temps, label='T_mean [°C]')
ax2.axhline(config.solver.T_target, ls='--', label='T_target')
ax2.set_xlabel('step'); ax2.legend()
plt.tight_layout(); plt.savefig('cooling_curve.png')
```

温度シミュレーションの標準出力例:
```
step=  100  vel_change=1.234e+00  div_max=3.21e-04  T_mean=31.80°C
step=  200  vel_change=8.543e-01  div_max=1.87e-04  T_mean=31.51°C
step=  500  vel_change=4.210e-01  div_max=9.32e-05  T_mean=30.91°C
...
Done: steps=5000  sim_time=250.0s  wall_time=38.2s  converged=False
  T_mean=30.41°C  (target=30.0°C)  div_max=6.12e-05
```

### コードから設定を変える

```python
from ff_room.config import SceneConfig, FanConfig, RoomConfig

config = SceneConfig.load_yaml("basic_room.yaml")

# 扇風機を動かして再実行
config.fans[0].position = (1.0, 2.0, 1.2)
config.fans[0].velocity = 4.0
result2 = SolverBridge(config).run()
```

---

## 3-2. アニメーション / 動画書き出し

シミュレーション中のフィールド変化を GIF または MP4 動画として保存できる。

### 基本的な使い方 (GIF)

```python
from ff_room import SceneConfig, SolverBridge, Animator

config = SceneConfig.load_yaml("my_room.yaml")
config.solver.max_steps = 500

# run_animated() でスナップショットを収集しながら実行
result, snapshots = SolverBridge(config).run_animated(
    save_every=25,        # 25 ステップごとにスナップショットを保存
    print_interval=50,    # 50 ステップごとに進捗を表示
)

print(f"スナップショット数: {len(snapshots)}")
print(f"最終ステップ: {result.steps}, 収束: {result.converged}")

# Animator で GIF に書き出す (Pillow が必要)
anim = Animator(snapshots, config)
anim.save_gif("results/flow_animation.gif", fps=10, dpi=100)
```

### MP4 動画の書き出し

```python
# MP4 (ffmpeg が必要; 未インストールの場合は GIF に自動フォールバック)
anim.save_mp4("results/flow_animation.mp4", fps=15, dpi=120)
```

### 温度シミュレーションの場合

```python
config = SceneConfig.load_yaml("examples/summer_cooling.yaml")
result, snapshots = SolverBridge(config).run_animated(save_every=100)

# 温度変化が 0.1°C 以上あると右パネルが温度カラーマップ (RdBu_r) になる
Animator(snapshots, config).save_gif("results/cooling.gif", fps=5)
```

### `FieldSnapshot` の中身

各スナップショットには以下のフィールドが含まれる:

```python
snap = snapshots[0]
snap.step     # int: ソルバーステップ番号
snap.time_s   # float: 物理時間 = step * dt [s]
snap.T_mean   # float: 室内平均温度 [°C]
snap.speed    # ndarray (Nx, Ny, Nz): 速度スカラー [m/s]
snap.T_field  # ndarray (Nx, Ny, Nz): 温度場 [°C]
snap.vel      # ndarray (Nx, Ny, Nz, 3): 速度ベクトル [m/s]
```

### 動画内容

| パネル | 内容 |
|--------|------|
| 左 (Top view) | XY 断面 (z=中央高さ) の速度カラーマップ + 方向矢印 |
| 右 (Side view) | YZ 断面 (x=中央) の速度 or 温度カラーマップ |

**パッケージ要件**:
- GIF: `pip install Pillow`
- MP4: システムに `ffmpeg` が必要 (`apt install ffmpeg` など)

---

## 4. 可視化

### 部屋セットアップ図 (シミュレーション前に確認)

YAML 設定を読み込んだあと、シミュレーションを実行する前に部屋の構成を図で確認できる。
扇風機の位置・向き、障害物、窓・ドアの配置を3方向から表示する。

```python
from ff_room import SceneConfig, ScenePlotter

config = SceneConfig.load_yaml("examples/summer_cooling.yaml")
ScenePlotter(config).plot(save="results/scene.png", show=False)
```

または:

```python
from ff_room import SceneConfig, plot_scene

plot_scene(config, save="results/scene.png", show=False)
```

生成される図:

| パネル | 内容 |
|--------|------|
| **Top view (XY)** | 床面からの俯瞰。ベッド・デスクの位置、扇風機の向き、南北壁の窓・ドア位置 |
| **Section YZ** | 扇風機の x 座標を通る縦断面。南北壁の窓・ドアの高さ範囲を確認 |
| **Section XZ** | 扇風機の y 座標を通る側断面。東西壁の開口部高さを確認 |

図タイトルに部屋サイズ、扇風機数、開口部数、温度設定 (T_in/T_out/target) が表示される。

`run_cooling.py` はシミュレーション前に自動でセットアップ図を保存する
(`results/<stem>_scene.png`)。

---

### 流れ場の可視化

PyVistaウィンドウの操作:
- マウス左ドラッグ: 回転
- マウス右ドラッグ / スクロール: ズーム
- マウス中ボタンドラッグ: 平行移動
- `q` または閉じるボタン: 終了

### 断面スライス

```python
viz = Visualizer(result)

# z方向 (床から fraction の高さ)
viz.slice_xy(z_frac=0.5)   # 部屋の中間高さ
viz.slice_xy(z_frac=0.2)   # 床近く
viz.slice_xy(z_frac=0.9)   # 天井近く

# y方向 (部屋の奥行き)
viz.slice_xz(y_frac=0.5)   # 中央縦断面

# x方向 (横断面)
viz.slice_yz(x_frac=0.3)   # 扇風機前
viz.slice_yz(x_frac=0.7)   # 扇風機後

# scalar: "speed", "pressure", "cell_type"
viz.slice_xz(y_frac=0.5, scalar="pressure")
```

### 複数断面を比較したい場合

```python
import pyvista as pv

pl = pv.Plotter(shape=(1, 3))
for col, z_frac in enumerate([0.2, 0.5, 0.8]):
    pl.subplot(0, col)
    # Visualizerの内部メソッドを活用
    grid = viz._make_structured_grid()
    slc = grid.slice(normal="z", origin=(0, 0, z_frac * viz.Lz))
    pl.add_mesh(slc, scalars="speed", cmap="viridis")
    pl.add_text(f"z={z_frac*viz.Lz:.1f}m", font_size=9)
pl.show()
```

---

## 5. 結果の保存

```python
from ff_room.io import ResultStore

ResultStore.save_npz(result, "results/run_fan_west.npz")
ResultStore.save_vtk(result, "results/run_fan_west.vts")    # ParaView
ResultStore.save_velocity_csv(result, "results/run_fan_west.csv")
```

VTKファイルは ParaView や VisIt で開ける。CSVはExcel/Pythonで集計可。

---

## 6. 配置比較ワークフロー

### 等温 (速度場のみ) 比較

```python
from ff_room import SceneConfig, SolverBridge, Visualizer
from ff_room.io import ResultStore, ExperimentLog

base = SceneConfig.load_yaml("examples/basic_room.yaml")
log  = ExperimentLog("results/comparison.jsonl")

configs = {
    "west_wall":   [0.5, 2.0, 1.2],
    "center":      [2.5, 2.0, 1.2],
    "corner_sw":   [0.5, 0.5, 1.2],
}

results = {}
for name, pos in configs.items():
    cfg = SceneConfig.from_dict(base.to_dict())  # deep copy
    cfg.fans[0].position = tuple(pos)
    cfg.metadata["name"] = name

    r = SolverBridge(cfg).run(print_interval=100)
    ResultStore.save_npz(r, f"results/{name}.npz")
    log.append(r, notes=name)
    results[name] = r
    print(f"{name}: converged={r.converged}, steps={r.steps}")

# 各配置の最大風速を比較
for name, r in results.items():
    print(f"{name}: max_speed={r.velocity_magnitude().max():.3f} m/s")
```

### 夏の換気: 扇風機位置ごとの冷却速度比較

`summer_cooling.yaml` を雛形にして `fan.direction` を変えて比較する例:

```python
from ff_room import SceneConfig, SolverBridge, Visualizer
from ff_room.io import ResultStore, ExperimentLog
import copy, matplotlib.pyplot as plt

base = SceneConfig.load_yaml("examples/summer_cooling.yaml")
log  = ExperimentLog("results/cooling_log.jsonl")

# 比較する扇風機設定
# direction: 270=窓から室内に吹き込む(-y), 90=窓に向かって吹く(+y)
variants = {
    "fan_blowing_in":  {"direction": [270.0, 0.0], "position": [3.0, 3.8, 0.9]},
    "fan_blowing_out": {"direction": [90.0, 0.0],  "position": [3.0, 3.8, 0.9]},
    "fan_at_center":   {"direction": [270.0, 0.0], "position": [3.0, 2.25, 0.9]},
}

results = {}
for name, fan_cfg in variants.items():
    cfg = SceneConfig.from_dict(base.to_dict())
    cfg.fans[0].direction = tuple(fan_cfg["direction"])
    cfg.fans[0].position  = tuple(fan_cfg["position"])
    cfg.metadata["name"]  = name

    print(f"\n=== {name} ===")
    r = SolverBridge(cfg).run(print_interval=500)
    ResultStore.save_npz(r, f"results/{name}.npz")
    log.append(r, notes=name)
    results[name] = r
    sim_t = r.steps * cfg.solver.dt
    print(f"T_mean={r.T_mean:.2f}°C  sim_time={sim_t:.0f}s  converged={r.converged}")

# T_mean を比較 (最終値)
print("\n=== 冷却効果比較 ===")
for name, r in results.items():
    delta = base.solver.T_initial - r.T_mean
    print(f"{name:20s}: T_mean={r.T_mean:.2f}°C  ΔT={delta:.2f}°C")
```

**結果の読み方**:
- `T_mean` が低いほど冷却効果が高い
- `converged=True` は T_target に達したことを意味する
- `sim_time` (物理時間) が短いほど早く冷える

---

## 7. トラブルシューティング

### シミュレーションが収束しない

1. `max_steps` を増やす (500 → 2000)
2. `dt` を小さくする (0.01 → 0.005)
3. 格子解像度を下げる (粗くして動作確認)
4. 扇風機風速が大きすぎないか確認 (6 m/s超は要注意)
5. `convergence_tol` を緩める (1e-4 → 1e-3)

### 計算が発散する (速度が無限大になる)

CFL条件違反。`dt` を減らす。目安: `dt = 0.3 * dx / U_max`

### Poissonソルバーが収束しない

`poisson_max_iter` を増やす (1000 → 5000)、または `poisson_tol` を緩める。

### `_ffroom_core` が見つからない

```
ImportError: C++ extension _ffroom_core not found
```
→ `bash build.sh` を再実行。`python/ff_room/` に `.so` ファイルが存在するか確認。

### PyVistaウィンドウが開かない

ヘッドレス環境 (サーバーなど) の場合:

```python
import pyvista as pv
pv.start_xvfb()   # Xvfbが必要: apt install xvfb
```

または off-screen レンダリングでPNGに書き出す:

```python
pl = viz.slice_xz(y_frac=0.5, show=False)
pl.screenshot("slice.png")
```

---

### 温度が下がらない (T_mean が初期値のまま)

**原因1**: 開口部が1つしかない → 空気の流入・流出ができない  
→ `openings` に最低2つ設定する (窓 + ドア少し開け)

**原因2**: 扇風機の向きが逆 → 室内の空気を外に向かって押している  
→ 窓から室内に吹き込む方向に設定する  
例: 北壁の窓なら `direction: [270, 0]` (-y 方向、室内に向かう)

**原因3**: `max_steps` が少ない  
→ 6×4.5×2.5m 部屋で 2°C 冷やすには約 5000 ステップ (dt=0.05) 以上が目安

**原因4**: `T_outside >= T_initial` になっている → 外が室内より暑い設定  
→ `T_initial` と `T_outside` の値を確認する

### T_target に到達しない

物理的な制限。部屋の体積に対して開口部の流量が不足している。

```
τ (時定数) ≈ 部屋体積 [m³] / 有効換気流量 [m³/s]
```

有効換気流量は `ファン風速 × ファン面積` の概算だが、開口部面積に制限される。  
→ 扇風機の `velocity` を上げる、または開口部の `size` を大きくする  
→ それでも難しければ `max_steps` を増やし長時間シミュレーション

### 温度場が可視化されない

PNG の温度パネルが表示されない場合:
- 温度変化が 0.05°C 未満だと温度パネルは非表示になる
- `T_initial`, `T_outside` が同じ値になっていないか確認
- `buoyancy: true` でも流れが起きないと温度変化しない (開口部の確認)
