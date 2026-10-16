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

初期データは `scripts/brill_initial_data.py` (2D ソルバ、`make id-table`) と
自作 thorn `thorns/BrillCollapse/BrillPsiTable` (表を Cactus の格子に補間) の組。

加えて **全ランでメモリ監視 (OOM 対策) と ETA 見積もりを必須**とする
(`scripts/run_monitor.py`, `make status` / `make monitor`)。

### 環境

- 兄弟プロジェクト [../gw230529-einstein-toolkit](../gw230529-einstein-toolkit) の
  イメージ `gw230529-et:local` (ET_2025_05 Kruskal, MPICH 単独, フル thornlist) に
  **薄いオーバーレイ**を被せる。必要な thorn は全てコンパイル済みであることを確認済み
  (`make docker-check` でアサート)。ソースからの再ビルド (60〜120 分, 17 GB) を避けるため。
  自作 thorn (`thorns/`) はオーバーレイでインクリメンタルビルド (1 分弱、レイヤー +2 GB、イメージ 19 GB)
- ベースイメージ `gw230529-et:local` は ECR タグ
  `575907694783.dkr.ecr.us-west-2.amazonaws.com/gw230529/einstein-toolkit:latest` と同一。
  ローカルタグが消えたら `docker tag` で戻す (ECR タグまで消すと全再ビルド)
- コンテナは `init: true` (tini が PID 1)。中断した MPI ランのゾンビを回収させるため
- コンテナの TZ は `Asia/Tokyo` (`.env` の `CONTAINER_TZ`)。status.json やログの時刻は JST
- コンテナ: `brill-collapse-et`、Jupyter は 8890 (8888: GW150914, 8889: GW230529)
- 出力は `SIM_OUTPUT_DIR` (既定 `~/brill-collapse-output`、ローカル SSD) →
  コンテナ内 `/home/etuser/simulations`
- コンテナは `memswap_limit = mem_limit` でスワップ無効。暴走ランはホストを
  スワップ地獄に落とさずコンテナ内で OOM kill される (runner のメモリガードが先に止める設計)

## 仕様の読み替え (ET の実在パラメータへの対応)

依頼仕様の名前の一部は ET に存在しない。ET_2025_05 の param.ccl で確認した対応:

| 仕様の記述 | ET での実体 |
| --- | --- |
| `BrillWave` ソーン | 自作 `BrillCollapse/BrillPsiTable` (`initial_data = "brill_psi_table"`) + `scripts/brill_initial_data.py`。IDBrillData (`q_function = "gundlach"`, `gundlach_b = gundlach_c = 2` 等で同じ Holz 型) は Carpet 上で使えない (下記 1) |
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
   `CT_Analytic` に Brill モードは無い。→ 方式 (a) に決定 (下記「決定事項」)
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

## Phase 1a: Cartoon2D feasibility の結果 (2026-10-15)

`par/feasibility/` の Brill-Lindquist (BL) テスト。グリッド、時間積分、ゲージ、散逸は Brill テンプレートと同じ
(4 レベル、最細 dx=0.025、RK4、Courant 0.25、1+log + Γ-driver η=1、KO ε=0.1)。共形因子は W 法。
`make run-par PAR=... CAMPAIGN=feasibility`、np=4 × OMP=4。

| parfile | 内容 | 結果 |
| --- | --- | --- |
| `cartoon_bl_single.par` | m=0.5 を原点に、t=5 | 完走 708 s、**25 /h**、0.9 GiB。M_irr=0.50000 が全区間一定 |
| `octant_bl_single.par` | 同じ物理を 3D 八分儀で、t=1 | 1578 s、**2.34 /h**、11.6 GiB。t=1 で lapse の中心値と最小値は cartoon と 5e-5 以内、地平線の M_irr と座標半径も一致 |
| `cartoon_bl_pair.par` | m=0.25 を z=±0.6125 に、t=6 | 完走 908 s、23 /h。軸上 finder (ah2) が M_irr=0.2755 (解析見積もり m(1+m'/2d)=0.2755) で t≈2.9 まで追跡。原点 finder (ah1) は t≈3 に共通地平線を捕捉、M_irr は 0.467→0.4991 |
| `cartoon_bl_small_l4.par` | m=0.1 (地平線の半径 ≈ 2 セル)、4 レベル | t≤0.6 は未検出、検出後も M_irr 0.095〜0.098 (誤差 2〜5%)、形が数値的に 20% 歪む |
| `cartoon_bl_small_l6.par` | 同、最細を 2 レベル追加 (dx=0.00625、8〜17 セル) | 完走 1065 s、**6.7 /h**、1.15 GiB。M_irr=0.10000 が t=0 から一定、rmin=rmax |

**結論: ML_BSSN + Carpet FMR + time refinement + AHFD (原点の octant finder と軸上の quadrant finder) は
Cartoon2D 上で動く。** コストは 3D 八分儀の約 1/10、メモリは 1/13。最細レベルを 2 つ足しても 3D の 4 レベル版より 3 倍速い。

Cartoon2D で必要な設定と注意点:

- `CoordBase::ymin = ymax = 0`、`boundary_shiftout_{x_lower,y_lower,y_upper,z_lower} = 1`、反射対称は z のみ
- **`CarpetLib::poison_value = 0`**。Cartoon2D は y=0 面の x<0 ゴーストを Lagrange 補間で埋めるが、そのステンシルに
  ゴースト自身を重み 0 で含む。t=0 の `ML_BSSN::Xt` のゴーストは未設定なので、NaN poison だと 0×NaN=NaN で停止する。
  代償として、未初期化の値を読んでも NaNChecker では検出できなくなる
- 軸上 finder は `origin_z[i]` とは別に **`initial_guess__coord_sphere__z_center[i]`** が要る (既定の 0 だと abort)
- CarpetReduce は Cartoon2D 格子で「領域体積 0」と警告し、norm2 系の値は 3D と比べられない (BL の |H|₂ が 14 倍違った)。
  x 軸上の H は地平線の外で cartoon/3D 比 0.83〜1.12 と一致。パンクチャ近傍 (急峻な場の回転補間) だけ 3〜5 倍大きい
- 地平線の分解能の床: 半径 2〜3 セルでは不正確、8 セル以上なら M_irr は 5 桁正確 → R_min ≈ 6〜8 セルが目安

## Phase 1b: 初期データの結果 (2026-10-16)

**ソルバ** (`scripts/brill_initial_data.py`): Δψ + (k/8)(q,ρρ + q,zz)ψ = 0 を極座標 (r,θ) で解く。
r = L tan(πs/2) (L=2) でコンパクト化、s と θ はセル中心、6 次差分、ゴーストは解の偶奇性
(原点・軸・赤道で偶、s=1 で奇: 外部解 Σ a_l P_l / r^{l+1} の 1/r は (1−s) の奇関数)。
SuperLU の直接解法で 512×64 が 0.2 秒。

- M_ADM の 3 つの推定量 (単極子の漸近係数、ソース積分 (1/2π)∫Vψ dV、Brill の (1/2π)∫|∇lnψ|² dV) が
  512×64 で 1e-10 以内で一致し、どれも 6 次で収束
- 解析解のテスト (erf(r)/r + r²e^{−r²}P₂) でラプラシアンが 6 次収束
- **文献値 (正振幅、中心 Holz 型、σ=1、e^{2q}) は全て報告された桁で一致**:
  Hilditch et al. 2017 Table 1 (A=4.6 → 0.597、4.69 → 0.619、4.697 → 0.621、4.7 → 0.622)、
  同 Fig. 3 (A=6.073 → 1.02)、Khirnov & Ledvinka 2019 (A=5 → 0.70、5.5 → 0.84)。
  本ソルバでは M(A*≈4.6967) = 0.62115、M(5) = 0.69912108、M(6) = 0.99170
- 規約: k=1 の A は k=2 の A/2 と完全に同じデータ (テスト済み)。文献は全て k=2
- 小振幅で M ≈ 0.0372 A²。負振幅は A≈−6.5 付近で演算子が特異になり、それ以上では ψ<0 (初期データが存在しない)
- **未解決の食い違い**: Hilditch et al. 2013 §IV.3 は中心・負振幅 A=−1, −2.5, −5 の質量を 0.61, 1.37, 3.15 と書くが、
  本ソルバは 0.042, 0.336, 3.21 (逆算すると 0.61, 1.37, 3.15 は A≈−3.13, −4.07, −4.98)。
  小振幅展開と整合するのは本ソルバの値で、同論文のオフセンター族の質量 (0.61, 1.36, 3.15) とほぼ同じ数字なので
  転記の誤りを疑っている。本番は正振幅なので保留

**thorn** (`thorns/BrillCollapse/BrillPsiTable`): PARAMCHECK で表を読み、A・σ・k・格子と偶奇性を parfile と照合
(不一致なら起動前に停止)。全格子点 (ゴースト含む) で 8 点の Lagrange 補間、計量は ψ と解析的な q から組み立て、
軸上の (e^{kq}−1)/ρ² は `expm1` で厳密。

- `par/initial_data/brill_id_check{,_hires}.par` (Cartoon2D、A=5、iteration 0):
  thorn の ψ と Python の `evaluate_table` の差は最大 9.9e-15。
  **t=0 のハミルトニアン拘束は解像度 2 倍で 15〜16.7 倍減少** (ML_ADMConstraints の 4 次差分どおり)。
  原点・軸上・全レベルで同じなので、表の補間誤差は効いていない
- HDF5 の固定長文字列は、読み出し側のメモリ型を「長さ + 1」にしないと最後の 1 文字が NUL で潰れる

## Phase 2: ゲージ feasibility の途中結果 (2026-10-16)

`templates/brill_cartoon.par.template` (Cartoon2D、4 レベル、最細 dx=0.025、W 法、KO ε=0.1、初期データは runner が
ランごとに生成、WeylScal4 の不変量、AHFD) のゲージ節だけを差し替えた `par/gauge/*.par.template` を
`make single-run AMP=... CAMPAIGN=gauge_<v> EXTRA="--template ..."` で回した。

| 変種 | ゲージ | A=5 の結果 |
| --- | --- | --- |
| S0 | 1+log + Γ-driver η=1 (仕様) | t=6.2 で crash、AH なし。曲率 |Re I| のピークは赤道面 x≈0.4〜0.9 |
| S0z | 1+log + シフト 0 | t=5.9 で crash、AH なし (KL 2019 の「1+log は AH 前に破綻」を再現) |
| S1 / S1z | harmonic + Γ-driver / シフト 0 | t=1.55 / 3.75 で crash (特異点回避が弱い) |
| **S2** | **shock-avoiding κ=1 + Γ-driver η=1** | **t=8.8 に AH (原点、初期推定 0.5)、M_irr 0.567→0.587**。扁平 (座標半径 0.3〜0.7)。形成後 t=10.05 に地平線内部で crash |
| **S2z** | **shock-avoiding κ=1 + シフト 0** | 初期推定 0.5 では未検出。初期推定 1.5 で **t=13.3 に AH、M_irr 0.5755→0.5775**。t=25 まで crash なし。座標半径が 1.1→2.2 と膨張して最細レベル (半径 2) からはみ出し、質量が 0.5705 まで下がる |

A=4.0 (S2) は **dispersion と判定、t=5.15 で早期終了** (曲率 463→9.6、lapse の最小値 0.14→0.59)。

**crash を遅らせる条件の探索** (A=5、早期終了なし。「安定」は変化率基準で後から判定した時刻):

| 変種 | AH 出現 / 確定 | crash | 形成後の生存 | 安定と判定 (M_irr) |
| --- | --- | --- | --- | --- |
| κ=1, η=1 | 8.8 / 9.0 | 10.05 | 1.25 | t=10.0 (0.5850) |
| κ=1, η=1 + 地平線内の散逸 ε→0.3 | 8.8 / 9.0 | 10.0 | 1.2 | 効果なし |
| κ=1, η=11.4 | なし | 9.05 | — | — |
| **κ=2/3, η=1** | 8.5 / 8.7 | 10.9 | 2.4 | **t=9.5 (0.5959)** |
| **κ=2/3, η=1, 6 レベル** | 8.5 / 8.7 | 11.85 | 3.35 | **t=9.5 (0.5954)** |
| **κ=1/2, η=1** | 7.4 / 7.6 | 11.25 | 3.85 | **t=8.6 (0.5967)** |
| κ=2/3, η=0 | 6.8 / 7.0 | なし (t=15) | 8.2+ | t=8.1 (0.6404)。その後 0.692 まで増加 (|H|∞ 20〜75、数値的と判断) |
| κ=1, シフト 0 | 13.3 / 13.5 | なし (t=25) | 11.7+ | t=14.3 (0.5756) |

- κ を下げると形成後の生存が延びる (1 → 1.25、2/3 → 2.4、1/2 → 3.85)。η を上げると悪化、地平線内の散逸は無効、
  解像度を上げると ~1 延びる
- **第一候補: shock-avoiding κ=2/3 (または 1/2) + Γ-driver η=1 + 変化率基準の早期終了**。crash 前に安定した質量が取れ、
  4 レベルと 6 レベルで 0.08% 一致
- ゲージ間で「安定した」M_irr は 0.576〜0.640 とばらつく (見かけの地平線はスライス依存、降着でゆっくり増える)。
  γ には族の中での一貫性が効くので G2 には許容
- 原点 finder の初期推定 0.5 は、κ ≤ 2/3 で内側の MOTS (質量が減る) を拾う。C2 (質量非減少) で棄却される

**臨界近傍 (κ=2/3, η=1)**: A=4.0 は dispersion (t=5.1 で早期終了)。**A=4.6 は t=6.35、A=4.8 は t=7.5 で判定前に crash**
(曲率のピークは赤道面のリング x≈1.0→0.68 が内側へ縮みつつ増大)。**A=4.8 は 6 レベルでも t=7.5 で crash → 解像度に依らず、
ゲージの病理**。このままでは二分探索が ΔA≈0.1 で止まる

- **地平線の質量は KL 2019 (quasi-maximal、A=5: t≈10.25 に AH、M 0.56→0.58) と一致**。M_irr はゲージ不変なので、
  初期データ→発展→AH 検出のパイプラインの定量検証になった
- 1+log と harmonic は AH 前に破綻。shock-avoiding は AH まで届くが、Γ-driver では形成後すぐ crash、
  シフト 0 では地平線が座標上で膨張する
- どちらも判定器の「質量が落ち着いた」(窓 max(10M,1) で相対 1e-3) を満たさない。Brill の BH は残りの波の降着で
  ゆっくり育つ (KL でも t=10→20 で 0.56→0.58)。M_AH の定義と閾値は Phase 2.5 / 3 で決める
- 原点 finder の初期推定は 1 本では足りない (シフト 0 では座標半径 > 1)
- shock-avoiding には ML_BSSN へのパッチ (`patches/ML_BSSN-shock-avoiding-slicing.patch`、`ML_BSSN::shockAvoidingKappa`)、
  WeylScal4 を Cartoon2D で使うにはテンソル型のパッチ (`patches/WeylScal4-cartoon-tensortypes.patch`) が要る。
  どちらもイメージに焼き込み済み (325e1689)。κ=0 で未パッチと、イメージ版とコンテナ内ビルドでビット一致を確認

## 判定・二分探索・解析の実装 (2026-10-16)

- **`scripts/outcome_classifier.py`**: 三値判定 (基準はモジュールの docstring)。崩壊は C1 連続検出 (既定 3 回)、
  C2 質量が減らない (相対 1e-3)、C3 地平線の最小座標半径 (既定 0 = 無効)。分散は D1 曲率不変量の大域最大値が
  初期値の 1/10 以下を窓 2.0 のあいだ維持、D2 一度も検出なし、D3 t ≥ 5、D4 その後に再成長なし、D5 最小 lapse ≥ 0.5 で
  下がっていない。分散条件の後に地平線 → undecided。OOM > 崩壊 (後の crash より優先) > crash > 分散 > timeout > undecided。
  **閾値は全て暫定値**で、Phase 2.5 のパイロットで較正する (`--classifier-config` の JSON で上書き、フィンガープリントに入る)
- **曲率の入力**: `weylscal4-curvir_group.norm_inf.asc` (WeylScal4 の Re I の全域 |max|、真空で Kretschmann = 16 Re I)。
  **テンプレートはまだこれを出力しないので、今のテンプレートでは分散と判定されない** (常に undecided)。
  ファイル名は WeylScal4 を有効にした最初のランで確認すること
- **早期終了**: 実行中 60 秒ごとに同じ判定を回し、崩壊で M_irr が安定した (直近 1.0 の最小二乗の傾きで
  |dM/dt|/M ≤ 2e-3 /時間) か分散が確定したら TerminationTrigger で止める (`--no-early-stop` で無効)。
  当初は「窓 max(10M,1) での相対変化 ≤ 1e-3」だったが、Brill の BH は残りの波の降着でゆっくり育つため
  (A=5 で 6 の時間に 0.3%) 変化率の基準に変更 (Syota さんの提案、2026-10-16)。さかのぼって当てると
  S2 は t=10.0 (M=0.5850)、S2z wide は t=14.3 (M=0.5756) で安定と判定される
- **`bisection_history.json`**: 原子的に保存。数値的な族のフィンガープリント (テンプレートの sha256、判定の閾値、
  t_final、cactus_sim のサイズと mtime、初期データの格子) を持ち、`--resume` で一致しなければ拒否。
  cactus_sim を再ビルドすると mtime が変わるので、キャンペーン中は `make docker-build` / `make thorn-build` をしない
- 非決定的な結果・単調性の破れ・端点の不一致で停止 (exit 2)。人が判断したら
  `make bisection-resume EXTRA="--set-outcome ITER=collapse"` (元の判定は `human_override` に残る) か `EXTRA="--retry"`
- `--verify-endpoints`: A_MIN が分散、A_MAX が崩壊であることを先に確認
- テンプレートに `__TABLE_FILE__` があれば、ランごとに初期データを解いて `psi_table.h5` を置く (`--id-n-s` 等)
- **`scripts/analyze_critical_exponent.py`**: 超臨界 ln M_AH = γ ln(A−A*) + c と亜臨界 ln K_max = −4γ ln(A*−A) + c' の
  両側でフィット (`make analyze SIDE=both MASS=settled`)。A* は lower / mid / value / fit (非線形、片側 4 点以上)
- テストは合成出力 (`tests/fake_outputs.py`) で 103 件。A*=4.7 の偽 Cactus で二分探索が端から端まで収束することも確認

## 決定事項 (2026-10-15)

| 項目 | 決定 | 理由・備考 |
| --- | --- | --- |
| 初期データ | **(a)** Python で (ρ,z) の線形楕円型方程式を解き、表を小さな thorn で各 Carpet レベルへ点ごとに補間 | 試行錯誤の速さ優先。ID がグリッド非依存、Cartoon2D と両立、M_ADM が副産物。ET 標準だけで完結させる意義は薄いと判断 |
| スライシング | 仕様の 1+log で破綻したら **harmonic / shock-avoiding への逸脱を許容** | KL 2019 で A=5 の 1+log 破綻が報告済み。ダメなら別方策でやり直す |
| 粗視化 | **パイロットキャンペーンのみ**。本番の Courant 数・L1 半径は据え置き | 手戻りリスク。早期終了と np/OMP/cpuset は数値解を変えないので対象外 |
| 軸対称化 | Cartoon2D を第一候補。Phase 1a の feasibility で判断 | 点数 1/15〜1/20 → 浮いたコストを最細レベルの追加へ |
| 合格基準 | **G2** で確定 (2026-10-16、Phase 1a で Cartoon2D が成立したため)。G3 はオプション | G1: 閾値のブラケット / G2: Type II らしさ (M_AH→0、亜臨界の曲率 max 増大) / G3: 超臨界・亜臨界の γ 推定が一致 |
| 判定 | **三値判定** (collapse / dispersion / undecided)。分散も積極的な証拠で示す | 見逃しも早すぎる分散判定も、どちらも A* を過大にする |
| 数値的な族 | 族 N (ID 設定、グリッド、Courant、FD 次数、ゲージ、ε、AHFD、判定閾値、イメージ) を本番前に凍結し、フィンガープリントで混在を拒否 | 途中で変えるなら新しい族として二分探索をやり直し、旧データは頑健性チェックに回す |

## 科学的な注意点

- **A\* の文献値**: 中心 Holz 型 (σ=1) で A* ∈ [4.6966875, 4.696703125]
  (擬スペクトル法、Hilditch らのグループ)。本設定 (FD, dx=0.025) の A* は解像度依存でずれる
- **1+log は臨界近傍の Brill データで地平線ができる前に破綻する**。
  Khirnov & Ledvinka 2019 ([arXiv:1908.06034](https://arxiv.org/abs/1908.06034), ET + Carpet + McLachlan、
  独自の analytic cartoon、8 次 FD): A=5 (M_ADM = 0.70) で 1+log は t≲4.9 に slice が非滑らかになり crash、
  **解像度に依らない**。Γ-driver も t≈10 で衝撃波状のシフトを作るので A=5 はシフト 0、
  検証ランは η = 11.4/σ。quasi-maximal slicing なら A=5 で t≈10.25 に AH (M_irr 0.56→0.58)、
  A=5.5 で t≈7.5 (0.63→0.73)。Hilditch et al. 2013 ([arXiv:1309.5008](https://arxiv.org/abs/1309.5008)) も同様の破綻を報告。
  代替: harmonic (`harmonicN=2, harmonicF=1`、パラメータのみ) / shock-avoiding f=1+κ/α²
  ([arXiv:2207.06376](https://arxiv.org/abs/2207.06376)、ML_BSSN の lapse 右辺へのパッチが要る) /
  quasi-maximal (楕円型 solve、重い)。BaikalVacuum, LeanBSSNMoL, ML_CCZ4 にも該当パラメータは無い
- **臨界近傍の曲率の大域ピークは z 軸上で原点から離れた位置**
  ([arXiv:1706.01829](https://arxiv.org/abs/1706.01829), [arXiv:2205.04379](https://arxiv.org/abs/2205.04379))。
  地平線は (0,0,±z0) に別々にできうるので、原点中心の AHFD 1 本では見逃しうる。
  軸上 finder (`"+xy quadrant (mirrored)"` を (0,0,z0) 中心) を追加する。
  スライシング次第でピークが軸から外れる報告もある (1706.01829)。トーラス型の地平線は AHFD では原理的に見つからない
- **M_AH の定義**: 個別地平線の落ち着いた M_irr (J=0 なので Christodoulou 質量)。
  フィットには 1 種類の地平線だけを使い、中心型と軸上ペア型が切り替わる ΔA は窓の端とする。
  全 finder の BH_diagnostics を残せば定義は後からでも変えられる
- **文献でも誤判定が起きている**: 1706.01829 は、他で亜臨界とされた振幅の一部が AH を作ると報告
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
| 1a | 監視付き単発実行 (メモリガード + ETA) の実装。**Cartoon2D feasibility**: Brill-Lindquist BH で ML_BSSN + FMR + time refinement + AHFD (軸上 finder 含む) を 3D 八分儀と比較 (完走不要) | ✅ 2026-10-15 成立 (上記) |
| 1b | 初期データ (a): Python ソルバ + 補間 thorn。M_ADM(A=5) = 0.70 で検証 (1a と並行可) | ✅ 2026-10-16 (下記) |
| 2 | ゲージ feasibility: A=5 で 1+log の破綻を再現 → harmonic → shock-avoiding。np/OMP/cpuset のベンチ | 進行中: 1+log の破綻を再現、shock-avoiding κ=2/3 + η=1 + 早期終了を第一候補に (上記)。臨界近傍 (A=4.8, 4.6) の確認とベンチは未 |
| 2.5 | パイロットキャンペーン (粗視化): z0、判定閾値の較正、A* の概算、crash の出方 | 未着手 |
| 3 | 数値的な族の凍結 (η, ε, FD 次数, finder 配置, 判定閾値)、三値判定の実装、既知解テスト | 一部: 三値判定・二分探索・両側の γ 解析を実装 (2026-10-16、合成データで検証)。閾値の較正、族の凍結、テンプレートの Cartoon2D 化は未 |
| 4 | 本番の二分探索キャンペーン | 未着手 |
| 5 | 超臨界の穴埋めスキャン + γ 解析 (超臨界 M_AH と亜臨界 K_max の両側) | 未着手 |

## 開発ルール

- スクリプトは**コンテナの Python 3.8** で動くこと (`from __future__ import annotations`、
  `match` 文や組み込みジェネリクスを実行時に使わない)。`make check` / `make test` はコンテナ内で実行する
- コード内のコメント、ログ、docstring (Google style) は英語
- **全ランでメモリ監視と ETA を必須とする**。監視なしで cactus_sim を直接起動しない。
  テンプレート以外の parfile も `make run-par PAR=... RUN_NAME=... CAMPAIGN=...` で監視下で回す
- cactus_sim の二重起動禁止 (`_precheck-idle`、ゾンビは数えない)
- 標準 thorn への変更は `patches/*.patch` に置き、`scripts/apply_cactus_patches.sh` が冪等に当てる
  (Dockerfile と `make thorn-build` の両方)。パッチの頭に理由を書く
- 自作 thorn の開発ループ: `thorns/` を編集 → `make thorn-build` (稼働中コンテナで再コンパイル、約 15 秒)
  → ラン。確定したら `make docker-build` でイメージに焼き込む (コンテナを作り直すとビルドは消える)。
  古いイメージは `docker image prune` で掃除する
- 初期データの表は `make id-table A=...` で `SIM_OUTPUT_DIR/tables/holz_A<A>_k2_<n_s>x<n_theta>.h5` に作る
- テンプレートに新しいパラメータを書くときは、ベースイメージの param.ccl で実在を確認する
- セッション記録は `.claude-notes/` (git 管理外)
