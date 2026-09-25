# hfauto の設計(現状)

- 本書は現在のコード(`hfauto/`)の説明である。計画時の設計書は `docs/reviews/2026-09-25_refactor_design.md`(凍結)にあり、コードの docstring にある「design §x.y」はその節番号を指す。
- 原則: stage は薄いアダプタで、化学の判定は `hfauto/chemistry/` の純関数だけが行う。外部ジョブの観測事実は型付きの `Evidence` で受け渡し、合否はゲート関数 1 か所で決める。バックエンドは能力 Protocol の実装で、subprocess を起動するのは `hfauto/execution/process.py` だけである。fail-closed で、互換層・ダミーエンジン・内部フォールバックは持たない。エンジン名の文字列は `hfauto/backends/` と `configs/` の外に書かない。

## 1. 層構造と import 契約

依存の向きは `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core` とする。複数の層で使う型は、使う層のうち最も下の層に置く。

| 層 | パッケージ(行数) | 役割 |
|---|---|---|
| cli | `hfauto/cli/`(186) | 5 コマンド: `doctor`、`run`、`status`、`report`、`case` |
| pipeline | `hfauto/pipeline/`(828) | 設定の読み込み、preflight、`RunLayout`、`Runtime`(`StageRuntime` の実装)、`run_pipeline` |
| stages | `hfauto/stages/`(1,362) | 8 つの stage。入力を得て driver か純関数を呼び、artifact を返すだけ |
| drivers | `hfauto/drivers/`(1,166) | `MinimumDriver`(`hfauto/drivers/minimum.py`)と `ReactionCaseDriver`(`hfauto/drivers/reaction_case/`) |
| reporting | `hfauto/reporting/`(563) | 順位・被覆率・手法パネルの表(`hfauto/reporting/summary.py`)と静的 HTML(`hfauto/reporting/html.py`) |
| backends | `hfauto/backends/`(2,443) | 能力 Protocol(`hfauto/backends/protocols.py`)、registry(`hfauto/backends/engines.py`)、6 つのアダプタ |
| execution | `hfauto/execution/`(744) | `run_command`、`JobStore`、`JobRunner`、`SiteLock`、worker |
| chemistry | `hfauto/chemistry/`(2,659) | ゲートと化学カーネル(純関数。ファイル IO は `hfauto/chemistry/xyz.py` と `hfauto/chemistry/xyz_trajectory.py` だけ) |
| core | `hfauto/core/`(755) | `Evidence`、payload、manifest、`MethodSpec`、system のモデル。依存は pydantic、numpy、PyYAML だけ |

import-linter の契約(`pyproject.toml`、`lint-imports` で検査。5 契約とも KEPT):
1. `hfauto-layers`: 上の層構造(exhaustive)。
2. `adapters-private`: cli・pipeline・stages・drivers・reporting はアダプタ(nwchem、xtb、pysis、crest、readuct、goodvibes)を import しない。
3. `protocols-only`: stages・drivers・reporting は `hfauto.backends.engines` と `hfauto.execution` を直接 import しない。エンジンは `StageRuntime.engine()` から Protocol 型で受け取る。
4. `adapters-independent`: 6 つのアダプタは互いに import しない。
5. `stages-independent`: 8 つの stage モジュールは互いに import しない。

ruff は `C901`(max-complexity 15)と `TID251`(`subprocess` は `hfauto/execution/process.py` だけ)を有効にしている。型検査は `pyrefly check`(0 errors)。

## 2. Evidence とゲート

`Evidence`(`hfauto/core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すること自体が次を保証する: (1) 正常終了と SCF 収束、(2) opt / saddle なら構造最適化の収束、(3) freq なら 3N − n_external 本の振動数、(4) 原子順序の保持、(5) 入力フレームの保持(出力に出た初期構造が入力と 1e-4 Å 以内)、(6) 出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること(`level_mismatches`)。満たさなければアダプタは `Failure` を返す。振動数と正準形式の Hessian(`.npy`、Eh/bohr²、入力フレーム)は、どのエンジンでも `hfauto/chemistry/vibrations.py` の補空間射影(並進・回転空間を SVD で決め、主同位体の質量を使う)で求める。経路は `PathProfile`、失敗は `Failure(kind: FailureKind, reason)` で表す。FailureKind は executable_missing、input_invalid、timeout、stagnated、nonzero_exit、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、budget_exhausted、gate_rejected の 11 種である。

ゲートは `hfauto/chemistry/gates.py` にだけ置く。閾値は `Policy` にまとめ、pipeline YAML の `gates:` だけが上書きできる。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=True)` | program・version・method・basis・dispersion・solvation・charge・multiplicity・電子温度が一致。`numerics=True` なら grid と scf_tol も。停留点・freq・QRC・熱化学は True、DFT string(種にだけ使う)は False |
| `spin_ok` | `s2` が None か、\|<S²> − S(S+1)\| ≤ 0.1。不合格は注記 `spin_contaminated`(rankable の blocker) |
| `imaginary_tier` | ν < −50 は saddle、−50 ≤ ν < −10 は soft、−10 ≤ ν < 0 は noise |
| `is_minimum(freq, opt=)` | task が freq / opt、`freq.start.fingerprint == opt.final.fingerprint`、numerics を含めて同一 PES、本数が 3N − n_external、tier が saddle でない(soft / noise は注記) |
| `is_first_order_saddle(freq, saddle=)` | 指紋と同一 PES、本数、最低モード < −50(ねじれは −20)、ほかに < −50 がない(あれば `higher_order`)、−50〜−10 は注記 `soft_secondary_mode`、E_TS − max(端点) ≥ max(2e-5 Eh, 20 × scf_tol)。モードの重なりは合否に使わない |
| `barrier_verdict` | 基準は高い方の端点。陰性証拠があれば `negative_evidence`、DFT 3 点未満は `unavailable`、DFT か低レベルの内部極大 ≥ 1.0 kcal/mol なら `proceed`、どちらも未満なら `barrierless`(ノード間隔を解像度として記録、接線モードの ½hν 未満なら `below_zpe`) |
| `connection` | 両側で drop = max(1e-5 Eh, 20 × scf_tol): 変位直後に E_TS − drop 以下、軌跡の最大も以下、最終値が先頭 − drop 以下、同一 PES。割付けの集合が期待どおりなら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつ両側の構造が別なら degenerate |
| `thermo_consistent` | GoodVibes の実振動の本数、ZPE(許容 1e-5 Eh)、E(許容 1e-6 Eh)が自前の振動数と一致 |
| `rankable` | outcome が elementary / degenerate / reassigned、ΔG‡ がある、blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated、method_sign_disagreement)がない、\|ΔE0‡ − ΔE‡\| ≤ 1 + 4 × max(1, 移動する H の数) kcal/mol |
| `reaction_tier` | connected(ConnectionClaim)> saddle(SaddleClaim)> minima(driver が BLOCKED 以外で終えた)> screening |

claim を作るのは 1 か所だけである。`MinimumRecord` は minima stage と、reaction-paths で新しい basin を見つけたときの `Registry` が作る。`SaddleClaim` と `ConnectionClaim` は reaction-paths の driver だけが作る。thermo と report は claim の有無だけを見て再検証しない(例外は thermo の `same_pes` と `thermo_consistent`)。

## 3. payload(artifact の型 8 種)

各 stage の `<run>/<stage_id>/manifest.json` には、その stage の出力だけを入れる。artifact は `status: success | failed` で、payload は `kind` で判別する frozen なタグ付き共用体(`hfauto/core/records.py`、`extra="forbid"`)である。

| 型 | レコード | 出す stage | 主な中身 |
|---|---|---|---|
| species | `SpeciesRecord` | structures、conformers、minima、explore、reaction-paths | 組成キー、電荷、多重度、`Geometry`、source、状態ラベル |
| calculation | `Evidence` | minima、reaction-paths、sp | 1 ジョブの観測事実。id は `calc_` + job_key の先頭 16 桁 |
| minimum | `MinimumRecord` | minima、reaction-paths | basin、tier(screen / dft)、level_key、opt / freq の calc id、members、notes |
| discovery | `DiscoveryRecord` | explore、minima(mode_follow) | 機構、trial、product / negative / failed と理由、低レベル TS |
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、source、`BarrierVerdict`、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome` |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None)、Boltzmann 重み、注記 |
| reaction_thermo | `ReactionThermo` | thermo | ΔE‡、ΔG‡、ΔG_rxn、ΔG_assoc、単量体基準の ΔG‡、感度の幅、blockers |
| report | `ReportRecord` | report | 順位の行と表(ranking.csv、coverage.csv、method_panel.csv、report.html) |

id の規則: species は `species_<species_id>`、minimum は `min_<species>_<level_key[:8]>`、失敗した極小は `min_<species>_<stage_id>`、分割した子反応は `<parent>_split<n>`、report は `report_<stage_id>`。

## 4. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス(`hfauto/backends/engines.py`) | 結果 |
|---|---|---|---|
| qm | `energy`、`optimize(tight, init_hessian)`、`frequencies` | `nwchem` → `NWChemEngine`(DFT と MP2 / CCSD(T) の SP)、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path, refine_ts)` | `nwchem_string` → `NWChemString`(ZTS)、`pysis_gs` → `PysisGrowingString`(xTB ネイティブの GS と TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode_index)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2 / AFIR、worker で実行) | `DiscoveryResult` |
| thermo | `thermo(freq, settings)` | `goodvibes` → `GoodVibesEngine`(4.3.0 の API、worker で実行) | `ThermoResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作り、`requirements()` は classmethod とする。`engines.create` は結果が能力の Protocol を満たさなければ `TypeError`、表にないキーは `KeyError` にする。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult` / `odft`、xTB と CREST の `--chrg` / `--uhf`、pysis の calc、ReaDuct の計算器設定)。
- 初期 Hessian は、対象と同じ構造の freq `Evidence` として渡す(Level は問わない)。指紋が違えば `INPUT_INVALID("hessian_geometry_mismatch")`。
- SCINE、pysisyphus、GoodVibes は `python -m hfauto.execution.worker <module:function> <job.json>` の子プロセスの中だけで import する。

エンジンを追加する手順:
1. アダプタ(`prepare` / `parse` / `monitor` / `continuation` と `result_type`)とエンジンクラスを `hfauto/backends/<engine>/` に書く。Level は出力から観測し、`level_mismatches` があれば `METHOD_MISMATCH`、振動数は `projected_frequencies` で求める。
2. `hfauto/backends/engines.py` の `_TABLE` に 1 行足し、`pyproject.toml` の `adapters-private` と `adapters-independent` にモジュールを加える。
3. **golden(必須)**: 実出力を `tests/golden/excerpt.py` の規則で抜粋して `tests/golden/data/<engine>/` に置き、`tests/golden/SOURCES.json` に出所と sha を記録し、パーサのテストを `tests/golden/` に書く。
4. **smoke(必須)**: `tests/smoke/` に `pytest.mark.real` のテストを書き、`real_engine` fixture で WSL の実エンジンを 1 回ずつ動かす(`docs/environment.md`)。
5. `configs/sites/wsl_local.yaml` に `EngineSite`(版数の pin、実行ファイル、`ExecutionSpec`)を加え、`hfauto doctor` で確認する。

## 5. MinimumDriver(`hfauto/drivers/minimum.py`)

`relax_to_minimum(mol, method, qm, *, known, init_hessian, policy, deadline, load_xyz)` の手順:
1. opt(`init_hessian` があれば使う)。`known` の Registry で最終構造が既存の極小に一意に割り付けば、freq を省いて `known`。
2. 別ジョブの freq → `is_minimum` と `spin_ok`。
3. saddle 級の虚振動なら最低モードに沿って ±0.1 Å 変位し、両側を opt → freq(最大 2 サイクル)。両側が同じ極小なら置き換え、別々の極小なら `ts_candidate`(minima stage が `mode_follow` の discovery にする)、片側だけならその側を採る。
4. soft なら 1 回だけ押し、残れば `soft_minimum`。noise はそのまま `minimum`。

`Registry` は組成 × level_key ごとの basin の集合で、`identity.compare_minima`(same: RMSD ≤ 0.02 Å かつ \|ΔE\| ≤ 1e-5 Eh、distinct: > 0.05 Å か > 5e-5 Eh)と `identity.assign`(≤ 0.05 Å かつ ≤ 5e-5 Eh で次点と十分離れていること)で判定する。minima stage は全構造を relax した後、species id の順に登録し、`ambiguous` は `tight` で再最適化して 1 回だけ再判定する。未解決事項: freq を省けるのは stage の開始時に登録済みの極小だけである(`docs/validation.md` §6)。

## 6. ReactionCaseDriver(`hfauto/drivers/reaction_case/`)

`decide(case, state, policy)` は `hfauto/drivers/reaction_case/state.py` の 17 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。`<run>/<stage_id>/cases/<reaction_id>/log.jsonl` は書くだけ(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 両端が同じ level_key の dft 極小でない | BLOCKED(`endpoints_not_on_one_pes`) |
| 2 | 両端が同じ basin で縮退でない | SAME_BASIN(`same_basin`) |
| 3 | ΔE_rxn > 40 kcal/mol | OUT_OF_WINDOW(`out_of_window`) |
| 4 | 予算(1 反応 6 h)切れ | UNRESOLVED(`walltime`) |
| 5 | 接続を判定済み | elementary / degenerate / reassigned で完了(`connection:<label>`)。failed で 2 回未満なら CONNECT(振幅 × 2、`connection_retry`)、それ以外は UNRESOLVED(`connection_failed`) |
| 6 | SaddleClaim がある | CONNECT(`ts_validated`) |
| 7 | saddle が収束し未検証 | VALIDATE_TS(`saddle_converged`) |
| 8 | TS 検証で崩壊(n_imag = 0) | VALIDATE_INTERMEDIATE(`saddle_collapsed`) |
| 9 | 中間体が両端と別の basin | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(深さ 2 まで) |
| 10 | 高次の saddle で試行 < 2 | REFINE_SADDLE(`higher_order_retry`、2 本目のモードに沿って 0.1 Å 下った種) |
| 11 | SCREEN が有効で未実行 | SCREEN(`screen`) |
| 12 | SCREEN が barrierless / negative_evidence | BARRIERLESS / NO_PRODUCT(`screen:<verdict>`) |
| 13 | 直近の DFT 経路が multi_max で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_multi_max`) |
| 14 | 直近の DFT 経路が monotonic | 2 回続けば BARRIERLESS(`dft_path_monotonic`)、そうでなければ 15 beads で FIND_PATH(`confirm_monotonic`) |
| 15 | 未使用の種があり saddle の試行 < 2 | REFINE_SADDLE(`seed:<source>`) |
| 16 | DFT string が未実行 | FIND_PATH(`no_dft_path`) |
| 17 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

- SCREEN: 低レベル TS(discovery / mode-follow)があれば xTB freq で確認し、その構造と両端の DFT SP の 3 点で判定する。なければ DFT の両端を xTB で再最適化し、崩壊しなければ `pysis_gs`(11 ノード、CI、同じ入力で TSOpt)の全ノードと TS で DFT SP、崩壊したら DFT の IDPP 11 点で判定する。種は xTB の TS(`screen_ts`)か DFT の HEI(`screen_hei`)。低レベルの失敗は `unavailable` で、FIND_PATH へ進む。
- REFINE_SADDLE: 1 回目は種の xTB Hessian(負の固有値が 1 本で接線との重なり ≥ 0.3 のとき)、それ以外は DFT Hessian を渡す。`mode_index` は接線(または宣言座標の勾配)と最も重なる負モード。
- FIND_PATH: 初期経路は低レベル GS の経路(写像付きで整列して再標本化)→ SCREEN の IDPP → 新しい IDPP の順。HEI を種に加えるのは `single_max` のときだけ。
- CONNECT(QRC): 振幅はエネルギー目標 max(3 × drop, 3e-4 Eh) と曲率から決めて 0.05〜0.4 Å に収める。± 変位を opt し、`Registry.find`(ねじれ反応は周期を考慮した最近傍)、未知なら `relax_to_minimum` で新しい basin にしてから `connection` で判定する。

## 7. 化学プロトコルの既定値(設計書 §8 を実装に合わせて更新)

| 段 | 既定値と規則 |
|---|---|
| 共通 | 虚振動 −50 / −10 cm⁻¹(ねじれ TS は −20)。DFT は PBE0-D3BJ/def2-SVPD、grid fine。**NWChem の DFT デッキはすべて `convergence energy` を書く**(既定 1e-7。opt・freq・SP を同じ数値層にそろえる)。geometry は `units angstrom nocenter noautosym`、autoz が失敗したら Cartesian で再投入。予算は 1 反応 6 h、ジョブの timeout 4 h |
| structures | xyz はそのまま、SMILES は RDKit ETKDG で 1 構造(`'.'` を含むものは不可)。`check_electronic_state` で電子数と多重度を検査し、宣言反応の両端は元素列が一致しなければ `INPUT_INVALID` |
| conformers | 単量体: `crest --gfn2 --quick -T 4 --ewin 6 --chrg --uhf`。重原子 3 個以下で回転可能結合 0 本なら省く。トポロジー停止は `crest_topology` の species として残し、`--noreftopo` で 1 回だけ再実行。**CREST は attempt ディレクトリで `--scratch` なしに実行する**(3.0.2 は scratch との往復で取り込んだ stdout を上書きする)。錯体: 受容原子の lone-pair 円錐(結合から 110〜120°、方位 0/120/240°、±30° 傾け)に極性 H を置く 6 seed。**結合が複数ある受容原子では、すべての結合から円錐角以上離れた方向(sp3 では lone-pair 軸)を使う**。donor / 受容原子がなければ剛体のランダム配置(vdW 半径の和 + 0.5 Å)。衝突は H···H ≥ 1.2 Å、重原子間 ≥ 2.2 Å。最良の seed から `--nci`、酸塩基の組成は labile H と受容原子に `--notopo`。状態ラベルごとに 4 kcal/mol 以内の 6 構造 |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル、`thread_map` で並列。状態ラベルが変わった seed は members に記録し、DFT には流さない |
| explore | 出発点は組成 × 状態ラベルごとの screen 最低 2 極小。trial は極性 H → C–H の 1,2- / 1,3- 移動 → 重原子の結合組み換え → 会合の順で出発点あたり 10 件まで。直線分子は 10° 曲げ 0.05 Å 乱してから。NT2 → AFIR(γ 125 → 300 kJ/mol)。SCC 失敗は 1000 K で再試行し 300 K の SP で評価し直す。窓(ΔE‡ > 150 kJ/mol、ΔE_rxn > 100 kJ/mol は `out_of_window`)は stage が適用する。**ReaDuct の `spin_mode` は常に `restricted_open_shell`**(scine-xtb-wrapper 3.0.2 が受け付けるのは any と restricted_open_shell だけ) |
| minima(dft) | 選択 `window`: 組成 × 状態ラベルごとに 6 kcal/mol 以内の 3 構造 + explore の生成物 + 入力単量体の最低構造。`rerank_sp` は組成 × 状態ラベルごとの上位 8 構造を DFT SP で並べ直す。`all` は入力 species 全部。開始構造は、screen 極小があればその最適化構造。**opt は `trust 0.1`**(既定 0.3 では QRC や mode-follow の変位から TS より上へ戻ることがある)、maxiter 100。2 断片以上は xTB の Hessian を初期 Hessian にする(W7 の A/B で点数 −33%) |
| 仮説(`hypotheses.select`) | 優先順: 宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補 → conformer ペア(結合変化なし、二面角差 ≥ 30°)。宣言以外は別 basin、ΔE_rxn ≤ 40 kcal/mol、変化 ≥ 0.2 Å か ≥ 30°、組成あたり 6 件まで。端点は各 basin の members から恒等写像 RMSD が最小の組(`pick_endpoints`)。`mapped_equivalent` なら縮退反応。同じ組成・結合変化の negative を陰性証拠に集める。**宣言反応の端点が screen の basin に崩壊した seed なら、その basin の DFT 極小を使う**(`_Pool.endpoint_basin`、W7 の修正 F1) |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`。最低以外のモードは `noautoz` と `moddir`。ZTS は nbeads 9(確認は 15)、stepsize 0.05、interpol 3、impose、freeze1 / freezeN、`tol 1e-5`、maxiter 20 × 3 チャンク、初期経路を必ず渡す。採否は gmax ≤ 1e-3 Eh/bohr で直近 3 反復が悪化なし。停滞(10 反復で改善 5% 未満、Xmax = 0 が 5 反復)は `STAGNATED`。`pysis_gs` は xTB 専用で、直線(全原子が主軸から 0.3 Å 以内)か 2 断片以上なら cart、それ以外は dlc。**厳密に直線の端点は主軸と垂直に ±0.01 Å ずらす**。TSOpt は rsirfo(2 断片以上は tric) |
| sp | `methods:` の各 LOT で、反応の極小と TS(`targets: all_minima` なら全 dft 極小も)の SP |
| thermo | GoodVibes 4.3.0 の API に自前の振動数を渡す。qs grimme、cutoff 100 cm⁻¹、symm。スケール因子は Truhlar DB の完全一致 → 同じ汎関数(PBE0 は PBE0/MG3S の 0.989 / 0.975、注記 `scale_factor_from:`)→ 1.0(`scale_factor_unverified`)。1 atm で計算して 1 bar / 1 M に換算。composite G = E_SP + (G_GV − E_GV)。会合量、感度の幅(qs × cutoff 50/100/150)、Boltzmann 重み。温度は pipeline の `conditions`。失敗は G = None(`thermo_unavailable`)、LOT の不一致は `mixed_level_of_theory` |
| report | rankable な反応だけを ΔG‡ で並べ、感度の幅が重なれば同順位。被覆率(機構別の試行・生成物・陰性理由・FailureKind)。手法パネルで ΔE_rxn の符号が割れたら `method_sign_disagreement` |

## 8. JobStore と attempt ladder

- JobStore(`hfauto/execution/jobstore.py`)は run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合し、`adapter.result_type` で復元する(検証に失敗したら miss)。
- 終端的な失敗(input_invalid、method_mismatch、executable_missing、incomplete_output)も記録し、`--retry-failed` で指定した種類だけを再実行する。ジョブ単位の排他は `.lock` を O_EXCL で作る。
- ladder(`hfauto/execution/jobs.py` の `LADDER`。これ以外の自動再試行はしない): timeout と geometry_maxiter は `continuation` で 2 回まで継続(NWChem は最新構造・movecs・drv.hess、xTB は最終構造)、input_invalid(autoz)は Cartesian で 1 回、scf_not_converged は前回の vectors に damping と level shift をかけて 1 回。attempt の timeout は min(設定値, 残り時間) で、残り 60 s 未満なら `BUDGET_EXHAUSTED`。
- 並行実行: `JobRunner` のセマフォで ranks × threads の合計を site の cores 以下に保つ。DFT と反応ケースは直列で、xTB と CREST だけを `thread_map` で並列にする。`SiteLock` は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。
- run ディレクトリ: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/diagnostics.json`(conformers と minima が書く。読むコードはない)、`<run>/run_state.json`(実行順の stage 状態と JobStore の統計)、`<run>/resolved_config.yaml`。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に union したもの。done で入力と設定の sha が変わらない stage は飛ばし(再開)、`--from X` は X 以降を stale にする。

## 9. 設定の 4 分割

層どうしはマージせず、優先順位もない。絶対パスは site ファイルにだけ書く。

| 層 | 置き場所 | 中身 |
|---|---|---|
| site | `configs/sites/wsl_local.yaml` | scratch_root、cores、memory_mb、エンジン別の `EngineSite`(版数の pin、実行ファイル、`ExecutionSpec`、worker の python、scratch_dir) |
| method | `configs/methods/` | `MethodSpec`(JobStore のキーに入る): `gfn2`、`pbe0-d3bj_def2-svpd`、`pbe0-d3bj_def2-svp`、`pbe0-d3bj_def2-tzvp`、`wb97x-d3_def2-tzvpd`、`ccsd-t_def2-tzvp` |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | species(xyz / SMILES、電荷、多重度、役割)、compositions、宣言反応: hcn、hono、nh3_inversion、formaldehyde、tma_hf2、amine_hf_panel、water_same_basin |
| pipeline | `configs/pipelines/` | stage の並びと設定、`conditions`、`gates:`(`Policy` の上書き): discover、known_endpoints、method_panel |

## 10. stage 表(生成)

次のコマンドの出力をそのまま貼った(`hfauto` を入れた Python で実行する)。

```bash
python - <<'EOF'
from hfauto.stages import catalog
print("| stage | クラス | consumes | produces | 設定キー(* は必須) |\n|---|---|---|---|---|")
for name, target in catalog.STAGES.items():
    spec = catalog.get(name).spec
    keys = ", ".join(k + "*" * f.is_required() for k, f in spec.config.model_fields.items())
    io = [", ".join(t.value for t in ts) or "—" for ts in (spec.consumes, spec.produces)]
    print(f"| `{name}` | `{target}` | {io[0]} | {io[1]} | {keys or '—'} |")
EOF
```

| stage | クラス | consumes | produces | 設定キー(* は必須) |
|---|---|---|---|---|
| `structures` | `hfauto.stages.structures:StructuresStage` | — | species | — |
| `conformers` | `hfauto.stages.conformer_search:ConformersStage` | species | species | engine*, method*, settings, seeds_per_composition, keep_per_state, window_kcal |
| `minima` | `hfauto.stages.minima:MinimaStage` | species | calculation, minimum, species, discovery | level*, engine*, method*, select, init_hessian, mode_follow |
| `explore` | `hfauto.stages.explore:ExploreStage` | minimum, species | discovery, species | engine*, method*, sources_per_state, max_trials_per_source, settings, window |
| `reaction-paths` | `hfauto.stages.reaction_paths:ReactionPathsStage` | minimum, species | reaction, calculation, minimum, species | method*, engines*, screen, policy, reaction_ids |
| `sp` | `hfauto.stages.single_point:SinglePointStage` | minimum | calculation | engine*, methods*, targets |
| `thermo` | `hfauto.stages.thermochemistry:ThermoStage` | minimum, calculation | species_thermo, reaction_thermo | engine*, settings, energy_method |
| `report` | `hfauto.stages.report:ReportStage` | — | report | T_K, standard_state, metric |

## 11. 設計書から変えた点

| 項目 | 設計書 | 現状(理由) |
|---|---|---|
| 削除したもの | core の units モジュール、`Level.surface_key`、`vibrations.curvature_along` と `frame_residual`、NWChem の最終振動ブロックの読み取り、`SpeciesRecord.formula`、`MinimumRecord.n_fragments`、`ReactionRecord.parent_id` | 使われていなかった。経路の PES 比較は `same_pes(numerics=False)` で行う。子反応の親は id(`<parent>_split<n>`)で分かる |
| MinimumDriver | `relax_to_minimum` の引数と `Registry(minima)` | `load_xyz=` を追加し、`Registry` は (MinimumRecord, Geometry) の組で作る。`MinimumOutcome.notes` を `MinimumRecord.notes` に写す |
| `CaseState` | 17 行の判断に必要な状態 | `minima`(両端の MinimumRecord)を追加した。行 1〜3 は tier・level_key・basin・エネルギーを要する |
| 電子状態 | `resolve_electronic_state` | `check_electronic_state(symbols, charge, multiplicity)` は検査だけを行う |
| minima stage | species と discovery を消費 | species だけを消費する(explore の生成物は species として渡る) |
| DFT 極小の freq の省略 | 既存の極小と same なら省く | stage の開始時に登録済みの極小だけが対象(未解決。§5) |
| pyproject | production extra、dev に mypy | production extra に `setuptools>=84` を加えた(scine_utilities が宣言せずに使う)。型検査は pyrefly |

## 12. 指標(Wave 8 の検証、設計書 §1.3)

| 指標 | 開始時(`5d76201`) | W8 前(`761c736`) | W8 後(`37af80a`) | 目標 |
|---|---:|---:|---:|---:|
| 本体の Python 行数 | 46,784(add-on と scripts を含む) | 11,050(81 ファイル) | 10,707(80 ファイル) | ≤ 12,500 |
| テスト行数 | 10,533 | 4,803 | 4,819 | ≤ 5,000 |
| golden データ | — | 566,624 B | 566,624 B | ≤ 1.2 MB |
| `configs/` 行数 | 1,714 | 307 | 307 | ≤ 600 |
| 文書(README を含み、レビューと設計書を除く) | 約 1,230 行・18 ファイル | 1,042 行・11 ファイル | 1,042 行・11 ファイル | ≤ 600 行・4 ファイル |
| stage の型 / artifact の型 / `qc` のキー | 41 / 65 / 287 | 8 / 8 / 0 | 8 / 8 / 0 | 8 / 8 / 0 |
| import-linter の契約 | 11(違反 7) | 5(違反 0) | 5(違反 0) | 5(違反 0) |
| 循環的複雑度の最大(15 超の関数) | 186 | 19(5) | 15(0) | ≤ 15 |
| 1 ファイルの最大行数 | 1,296 | 494 | 499 | ≤ 500 |
| 既定のテスト実行 | 270 passed / 2 skipped、2.1 s | 273 passed、10.7 s | 274 passed、10.8〜11.2 s | 15 s 未満 |
| pyrefly | — | 0 errors / 17 warnings | 0 errors / 10 warnings | — |

パッケージ別の行数(W8 前 → 後): core 804 → 755、chemistry 2,950 → 2,659、backends 2,452 → 2,443、execution 749 → 744、drivers 1,162 → 1,166、stages 1,352 → 1,362、pipeline 831 → 828、reporting 563、cli 186。
文書は W9 で 4 ファイル・475 行になった(README 110、design 222、environment 70、validation 73)。
