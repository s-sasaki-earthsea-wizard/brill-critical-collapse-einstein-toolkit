# Brill Wave Vacuum Critical Collapse (Einstein Toolkit)

## プロジェクト概要

Einstein Toolkit (ET) で**純粋な真空時空における Brill 波の臨界崩壊**を計算する。
振幅 A の二分探索で臨界振幅 A* を追い込み、超臨界ランの見かけの地平線質量
M_AH から `M_AH ∝ (A − A*)^γ` の臨界指数 γ (目標 0.3〜0.4 程度) を測る。

### 成果物

1. `templates/brill_collapse.par.template` — `__AMPLITUDE__` / `__FINAL_TIME__` を置換するパラメータファイル
2. `scripts/bisection_runner.py` — 二分探索の制御 (実行・判定・`bisection_history.json`)
3. `scripts/analyze_critical_exponent.py` — γ の線形回帰と `critical_scaling.png`
4. `README.md` — 前提条件、手順、ディレクトリ構成

加えて **全ランでメモリ監視 (OOM 対策) と ETA 見積もりを必須**とする
(`scripts/run_monitor.py`, `make status` / `make monitor`)。

### 環境

- 兄弟プロジェクト [../gw230529-einstein-toolkit](../gw230529-einstein-toolkit) の
  イメージ `gw230529-et:local` (ET_2025_05 Kruskal, MPICH 単独, フル thornlist) に
  **薄いオーバーレイ**を被せる。必要な thorn は全てコンパイル済みであることを確認済み
  (`make docker-check` でアサート)。ソースからの再ビルド (60〜120 分, 17 GB) を避けるため
- コンテナ: `brill-collapse-et`、Jupyter は 8890 (8888: GW150914, 8889: GW230529)
- 出力は `SIM_OUTPUT_DIR` (既定 `~/brill-collapse-output`、ローカル SSD) →
  コンテナ内 `/home/etuser/simulations`
- コンテナは `memswap_limit = mem_limit` でスワップ無効。暴走ランはホストを
  スワップ地獄に落とさずコンテナ内で OOM kill される (runner のメモリガードが先に止める設計)

## 仕様の読み替え (ET の実在パラメータへの対応)

依頼仕様の名前の一部は ET に存在しない。ET_2025_05 の param.ccl で確認した対応:

| 仕様の記述 | ET での実体 |
| --- | --- |
| `BrillWave` ソーン | `EinsteinInitialData/IDBrillData`: `initial_data = "brilldata"`, `q_function = "gundlach"`, `gundlach_a = A`, `gundlach_b = gundlach_c = 2`, `gundlach_r0 = 0`, `gundlach_sigmarho = gundlach_sigmar = 1` で Holz 型 `q = A ρ² e^{−r²}` |
| `lapse_advection_coeff = 1.0` | `ML_BSSN::advectLapse = 1` (`LapseAdvectionCoeff` は OUTDATED 扱い) |
| `Shift::Spatial_Gauge = "Gamma-driver"` | `ML_BSSN::shiftGammaCoeff = 0.75`, `advectShift = 1`, `betaDriver` |
| 1+log スライシング | `ML_BSSN::harmonicN = 1`, `harmonicF = 2.0` |
| AHFD の「クォータニオン/質量ファイル」 | 形状 `h.t<it>.ah1.gp`、質量は `BH_diagnostics.ah1.gp` の `m_irreducible` 列 |
| 八分儀対称 + AHFD | `patch_system_type = "+xyz octant (mirrored)"` を**明示必須**。AHFD の自動判定は `CartGrid3D::domain = "octant"` しか見ず、ReflectionSymmetry 構成では full sphere を選んでしまう |

## 調査・検証で判明した事実 (2026-10-15)

1. **IDBrillData は Carpet の FMR 上で使えない** (最大の未決事項)。
   ハミルトニアン拘束を PUGH 時代の EllBase/EllSOR で解く実装で、EllSOR は
   単一ローカルパッチをループし、local mode のルーチン内から大域リダクションと
   sync を呼ぶ。Carpet 上では各リファインメントレベルを**そのレベル自身の外縁に
   Robin 条件を置いて独立に解く**ことになり、物理的に誤り。テストスイートも
   PUGH unigrid のみ。dx=0.025 の unigrid は八分儀でも ~640³ 点で非現実的。
   `CT_Analytic` に Brill モードは無い。→ 代替方式の決定が Phase 1
2. **ベースイメージに `OMP_NUM_THREADS` と Cactus の PATH が焼かれていない**
   (GW230529 の知見 7 の修正より前のビルド)。オーバーレイの Dockerfile で
   `ENV OMP_NUM_THREADS=1` と PATH を設定、runner も全 rank に `-genv` で渡す
3. **テンプレートはパラメータレベルで有効**。初期データ未設定 (= ADMBase 既定の
   平坦時空) で `t_final = 0` 起動し exit 0、パラメータ警告ゼロを確認
4. **TerminationTrigger のファイルは `IO::out_dir` 下** (`output/TERMINATE`) に作られる
5. **メモリ**: Carpet 申告 10.9 GB (207 GF)。レベル別の点数は
   L0 126³ = 2.00M, **L1 133³ = 2.35M**, L2 93³ = 0.80M, L3 93³ = 0.80M。
   半径 12 の L1 が最大の消費源。np=4 で 1 rank あたり maxrss 2.8 GB、コンテナ合計 11 GiB
6. **暫定スループット (平坦時空, np=4 × OMP=4, cpuset 0-15)**:
   Carpet の `physical_time_per_hour` ≈ **2.07** (64 iteration の間ずっと安定)、
   粗いステップ (8 iteration, Δt=0.05) あたり約 87 秒。
   **t_max = 30 で 1 ランあたり約 14.5 時間、7 回の二分探索で 4 日強**。
   CPU 使用率は 848% / 1600% で、並列効率が悪い。np/OMP/cpuset の最適化と
   グリッド (L1 半径) の見直しの余地が大きい。AHFD の挙動が本番と異なるので暫定値
7. **ホスト CPU は異種コア構成** (Core Ultra 9 285H): CPU 0-5 が P コア 5.4 GHz、
   6-13 が E コア 4.5 GHz、14-15 が LP E コア 2.5 GHz。バルク同期の MPI は最遅 rank に
   律速されるため、cpuset の選び方がそのまま ETA に効く

## 科学的な注意点

- **A\* の文献値**: 中心 Holz 型 (σ=1) で A* ∈ [4.6966875, 4.696703125]
  (擬スペクトル法、Hilditch らのグループ)。本設定 (FD, dx=0.025) の A* は解像度依存でずれる
- **moving-puncture ゲージは臨界近傍の Brill データで破綻することが報告されている**
  ([arXiv:1309.5008](https://arxiv.org/abs/1309.5008))。A* 近傍の crash は想定内。
  runner は crash を勝手にどちらかの側に分類せず、キャンペーンを止めて人間に判断を仰ぐ
- **臨界近傍では地平線が原点ではなく z 軸上の対称な 2 点に別々にできる**
  ([arXiv:1706.01829](https://arxiv.org/abs/1706.01829))。原点中心・半径 0.2 の
  AHFD 1 本では見逃しうる。軸上の地平線 (`"+xy quadrant (mirrored)"` を (0,0,z0) 中心) の追加が必要。
  その場合の「M_AH」の定義 (2 つの和か、共通地平線か) も決める必要がある
- **普遍性は無いとする結果が出ている**
  ([Baumgarte et al. 2023, PRL 131, 181401](https://doi.org/10.1103/PhysRevLett.131.181401))。
  測定した γ は「この族・この設定の値」であり、0.3〜0.4 は合否基準ではなく比較対象
  (歴史的な値は Abrahams & Evans 1993 の γ ≈ 0.37)
- **質量レンジは解像度で頭打ちになる**。γ = 0.37 なら M_AH の 1 桁は ΔA の約 2.7 桁に相当する。
  最細 dx = 0.025 で分解できる最小の地平線を考えると、M_AH は 1 桁程度しか振れない見込み
- **7 回の二分探索では最終幅は 2/2⁷ ≈ 0.016**。A* に最も近い超臨界点の ΔA は
  この幅と同程度で、ln ΔA の誤差が支配的になる。γ を真面目に測るなら
  (a) 二分探索をもっと深く回し、(b) A* + 対数等間隔の ΔA で専用の超臨界スキャンを行い、
  (c) A* を自由パラメータにした非線形フィットで比較する

## Phase 計画

| Phase | 内容 | 状態 |
| --- | --- | --- |
| 0 | スケルトン、テンプレートのパラメータ検証、暫定ベンチ、init commit | ✅ 2026-10-15 |
| 1 | 初期データ方式の決定と実装 (IDBrillData の代替) + 検証 (ADM 質量の文献値比較、拘束違反) | 未着手 |
| 2 | runner の単発実行 + メモリ監視 + ETA を実装、np/OMP/cpuset とグリッドのベンチ | 未着手 |
| 3 | 単発 smoke: A=4.0 (分散) / A=6.0 (崩壊)、判定基準と AHFD 設定 (軸上地平線) の確定 | 未着手 |
| 4 | 判定・二分探索ロジックを実装してキャンペーン実行 | 未着手 |
| 5 | 超臨界スキャン + γ 解析 | 未着手 |

## 開発ルール

- スクリプトは**コンテナの Python 3.8** で動くこと (`from __future__ import annotations`、
  `match` 文や組み込みジェネリクスを実行時に使わない)。`make check` / `make test` はコンテナ内で実行する
- コード内のコメント、ログ、docstring (Google style) は英語
- **全ランでメモリ監視と ETA を必須とする**。監視なしで cactus_sim を直接起動しない
- cactus_sim の二重起動禁止 (`_precheck-idle`)
- テンプレートに新しいパラメータを書くときは、ベースイメージの param.ccl で実在を確認する
- セッション記録は `.claude-notes/` (git 管理外)
