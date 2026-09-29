# hfauto の設計

本書は現在のコード(`hfauto/`)の説明である。対象範囲と使い方は README、環境は [environment.md](environment.md)、実計算の記録は [validation.md](validation.md) にある。コードの docstring にある「design §x.y」は、計画時の設計書 `docs/reviews/2026-09-25_refactor_design.md`(凍結)の節番号を指す。

## 0. 原則

- 化学の判定は `hfauto/chemistry/` の純関数だけが行う。stage は入力を集めて driver か純関数を呼び、artifact を返す薄いアダプタである。
- 外部ジョブの観測事実は型付きの `Evidence` で受け渡し、合否はゲート関数(`hfauto/chemistry/gates.py`)だけが決める。
- バックエンドは能力 Protocol の実装で、subprocess を起動するのは `hfauto/execution/process.py` だけである。エンジン名の文字列は `hfauto/backends/` と `configs/` の外に書かない。
- fail-closed: 互換層、ダミーエンジン、内部フォールバックは持たない。失敗は `FailureKind` 付きで最小の単位に閉じ込め、run は続ける(§8)。
- 停留点(opt・saddle・string・freq)は大域混成か GGA に D3 を足した RKS / UKS で計算する(NWChem の解析 Hessian が使える範囲。RSH・meta-GGA では数値 Hessian になり 3〜4 倍かかる)。RSH と CCSD(T) は SP のエネルギー層だけで使う。気相専用で、溶媒和は持たない。

## 1. 層構造と import 契約

依存の向きは `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core`。複数の層で使う型は、使う層のうち最も下の層に置く。

| 層 | 役割 |
|---|---|
| cli | 5 コマンド(`doctor`、`run`、`status`、`report`、`case`) |
| pipeline | 設定の読み込み、preflight、run ディレクトリ(`RunLayout`)、`StageRuntime` の実装、`run_pipeline` |
| stages | 8 つの stage |
| drivers | `MinimumDriver`(`drivers/minimum.py`)と `ReactionCaseDriver`(`drivers/reaction_case/`) |
| reporting | 順位・被覆率・手法パネルの表と静的 HTML |
| backends | 能力 Protocol(`backends/protocols.py`)、registry(`backends/engines.py`)、5 つのアダプタ |
| execution | `run_command`、`JobStore`、`JobRunner`、`SiteLock`、worker |
| chemistry | ゲートと化学カーネル(純関数。ファイル IO は xyz の読み書きだけ) |
| core | `Evidence`、payload、manifest、`MethodSpec`、system のモデル |

import-linter の契約(`pyproject.toml`、`lint-imports` で検査):
1. `hfauto-layers`: 上の層構造(exhaustive)。
2. `adapters-private`: cli・pipeline・stages・drivers・reporting はアダプタ(nwchem、xtb、pysis、crest、readuct)を import しない。
3. `protocols-only`: stages・drivers・reporting は `hfauto.backends.engines` と `hfauto.execution` を直接 import しない。エンジンは `StageRuntime.engine()` から Protocol 型で受け取る。
4. `adapters-independent`: アダプタどうしは import しない。
5. `stages-independent`: stage モジュールどうしは import しない。

ruff は `C901`(max-complexity 15)と `TID251`(`subprocess` の禁止)を有効にし、型検査は `pyrefly check` で行う。

## 2. 処理区分の責務

レビュー `2026-09-27_platform_review.md` §1.1 の区分(U0〜U9)とコードの対応。

| 区分 | 責務 | 実装 | 出力 |
|---|---|---|---|
| U0 プロトコル | 段ごとの理論レベルと順位の面を決め、LOT をエンジンに渡して出力の観測で照合する(Z>36 の def2-ECP は自動) | `configs/methods/`、`core/method.py`、アダプタの入力生成 | `Level`、デッキ |
| U1 入力・電子状態 | 化学種と組成を、対応元素・電荷・多重度・状態ラベルを持つ構造にする。範囲外は入口で 1 回だけ拒否する | structures、`chemistry/electronic_state.py`・`elements.py`・`topology.py` | `SpeciesRecord` |
| U2 配座・錯体配置 | screen に渡す候補構造を作る(同一性の確定と窓の選抜はしない) | conformers、`chemistry/placement.py`、CREST | 状態ラベルごとの候補 |
| U3 極小の確定 | 同じ LOT で確かめた極小を basin に登録する。鞍点に落ちた入力は mode-follow で解消する | minima、`MinimumDriver`、`chemistry/identity.py`・`selection.py` | `MinimumRecord` |
| U4 反応探索 | 生成物を仮定しない片端探索で、低レベルの固有な生成物と TS を出す | explore、`chemistry/trials.py`、ReaDuct | `DiscoveryRecord` |
| U5 仮説・経路 | 極小対から仮説を作り、DFT 極小を固定端とする経路を分類して種を渡す | `chemistry/hypotheses.py`・`profile.py`、SCREEN・FIND_PATH | `ReactionRecord`、`BarrierVerdict` |
| U6 鞍点・接続 | 種を 1 次の鞍点に精密化し、別ジョブの freq で TS を確かめ、QRC で両側の極小につなぐ | REFINE_SADDLE・VALIDATE_TS・CONNECT・VALIDATE_INTERMEDIATE | `SaddleClaim`、`ConnectionClaim`、outcome |
| U7 熱化学 | 停留点の G(qRRHO、キラリティ、エネルギー層)と順位の量 δG_eff を出す | thermo、`chemistry/thermo.py`(GoodVibes) | `SpeciesThermo`、`ReactionThermo` |
| U8 一点計算・報告 | エネルギー層を順位に使う点だけで計算し、δG_eff で並べ、手法の幅を列で示す | sp、report、`reporting/` | ranking.csv、method_panel.csv、report.html |
| U9 実行基盤 | stage を動かし、ジョブを再利用・継続し、失敗を閉じ込め、コアと予算を管理する | `pipeline/`、`execution/`、`cli/` | manifest、JobStore、終了コード |

## 3. Evidence とゲート

`Evidence`(`core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すれば次が成り立つ: 正常終了と SCF 収束、opt / saddle の収束、freq なら 3N − n_external 本の振動数(単原子は 0 本)、原子順序と入力フレームの保持、出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること。満たさなければアダプタは `Failure(kind, reason)` を返す(maxiter で止まった saddle は最終フレームも持つ)。振動数と正準形式の Hessian(Eh/bohr²)は、どのエンジンでも `chemistry/vibrations.py` の補空間射影(並進・回転を除き、主同位体の質量)で求めるので、NWChem の表示値とは低振動モードで数 cm⁻¹ ずれる(G で約 0.02 kcal/mol)。

FailureKind は executable_missing、input_invalid、timeout、nonzero_exit、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、budget_exhausted、gate_rejected の 10 種。

化学の閾値は `Policy` の 5 つ(§9)で、数値の許容幅 `qrc_drop` = max(1e-5 Eh, 20 × scf_tol)(同じ幾何での SCF の雑音)は gates のモジュール定数である。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=, state=)` | program・version・method・basis・dispersion・電子温度・電荷・多重度が一致。`numerics=True` なら grid と scf_tol も(停留点・freq・QRC・熱化学)。`state=False` は会合量の錯体と単量体の比較だけ |
| `spin_ok` | ⟨S²⟩ がないか、\|⟨S²⟩ − S(S+1)\| ≤ spin_tol。freq と、thermo が使うエネルギー層の SP にかける。不合格は blocker `spin_contaminated` |
| `imaginary_tier` | 最低の ν < −saddle_cm1 は saddle、< −noise_cm1 は soft、< 0 は noise |
| `is_minimum(freq, opt=)` | freq が opt の最終構造で同一 PES(numerics を含む)、\|E_freq − E_opt\| ≤ qrc_drop(外れたら `state_mismatch`)、tier が saddle でない(soft・noise は注記) |
| `is_first_order_saddle(freq, saddle=)` | 同じ連結条件で、最低モード < −noise_cm1(なければ `no_imaginary_mode`)、2 本目 < −saddle_cm1 なら `higher_order`(2 本目が −saddle_cm1〜−noise_cm1 なら注記 `soft_secondary_mode`)。モードの大きさ・重なり・端点とのエネルギー比較は見ない |
| `barrier_verdict` | 両端が DFT 極小のエネルギーであるプロファイルを分類する。山と井戸を同じ関数で resolution_kcal の深さから数え、井戸があれば intermediate、井戸がなく山があれば single、どちらもなければ barrierless、3 点未満は unavailable |
| `connection` | QRC の両側が TS と同一 PES で E_TS − qrc_drop より下の極小に割り付くこと(途中の軌跡は見ない)。期待の組なら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつラベル付きで別構造なら degenerate(結合が変わる縮退反応では両側の結合グラフの組が {bonds(R), bonds(P)} であること。違えば `bond_change_missing`)、それ以外で両側が同じ basin なら failed(`sides_same_basin`) |
| `discovery_verdict` | explore の採否: NT2 は検証済みの TS と IRC が要る。ΔE_rxn ≤ reaction_window_kcal(評価できなければ窓の外)、低レベルの障壁 ≤ 50 kcal/mol |
| `rankable` | outcome が elementary / degenerate / reassigned / barrierless_at_resolution で、ReactionThermo に blocker がない |
| `reaction_tier` | connected(ConnectionClaim)> saddle(SaddleClaim)> minima > screening |

claim を作るのは 1 か所だけである。`MinimumRecord` は minima stage と reaction-paths の `Registry`、`SaddleClaim` と `ConnectionClaim` は `ReactionCaseDriver` が作る。blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated)は thermo stage だけが作る。thermo と report は claim を再検証しない(thermo の `same_pes` は例外)。

## 4. payload

各 stage の `<run>/<stage_id>/manifest.json` にはその stage の出力だけを入れる。artifact は `status: success | failed` を持ち、payload は `kind` で判別する frozen なタグ付き共用体(`core/records.py`、`extra="forbid"`)である。

| 型 | レコード | 出す stage | 主な中身 |
|---|---|---|---|
| species | `SpeciesRecord` | structures、conformers、minima、explore、reaction-paths | 組成キー、電荷、多重度、構造、source、状態ラベル |
| calculation | `Evidence` | minima、reaction-paths、sp | 1 ジョブの観測事実(id は `calc_` + job_key の先頭 16 桁) |
| minimum | `MinimumRecord` | minima、reaction-paths | basin、tier(screen / dft)、level_key、opt / freq の calc、members、chiral(m = 2)、注記 |
| discovery | `DiscoveryRecord` | explore、minima(mode_follow) | 機構、trial、product / negative / failed と理由、低レベル TS、ΔE‡・ΔE_rxn。DFT 段の mode-follow は鞍点の opt(`ts_calc`) |
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、source、検証済みの DFT 鞍点(`ts_calc`)、`BarrierVerdict`、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome` |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None) |
| reaction_thermo | `ReactionThermo` | thermo | δG_eff と感度の幅、ΔE‡、ΔE_rxn、ΔZPE‡、ΔG‡、ΔG_rxn、ΔG_assoc、分離反応物基準の ΔG‡、blockers、`energy_level`、注記 |
| report | `ReportRecord` | report | 順位の行と表のファイル |

CaseOutcome は elementary_step、degenerate_rearrangement、reassigned_step、multi_step、barrierless_at_resolution、same_basin、out_of_window、unresolved_within_budget、blocked_upstream の 9 種。id は species が `species_<id>`、minimum が `min_<species>_<level_key[:8]>`、分割した子反応が `<parent>_split<n>`。

## 5. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス | 結果 |
|---|---|---|---|
| qm | `energy`、`optimize(init_hessian)`、`frequencies` | `nwchem` → `NWChemEngine`(DFT と CCSD(T))、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path)` | `nwchem_string` → `NWChemString`(ZTS)、`pysis_neb` → `PysisNEB`(xTB の CI-NEB と TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2 / AFIR) | `DiscoveryResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作る。`engines.create` は Protocol を満たさなければ `TypeError`、表にない名前は `KeyError`。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult`、xTB と CREST の `--chrg` / `--uhf`、pysisyphus と ReaDuct の計算器設定)。
- 初期 Hessian は freq の `Evidence` で渡す。saddle は種と同じ構造、optimize は同じ原子順序で各原子が 0.5 Å 以内の構造のものだけを受け付け、外れれば `INPUT_INVALID(hessian_geometry_mismatch)`。xTB は初期 Hessian を使わない。
- SCINE と pysisyphus は `python -m hfauto.execution.worker <module:function> <job.json>` の子プロセスの中だけで import する。
- 熱化学はエンジンではない。qRRHO は閉じた式なので、`chemistry/thermo.py` の `species_thermo` が GoodVibes 4.3.0 を同じプロセスで呼ぶ(ジョブにしない)。版数は preflight が照合し、結果は golden で保証する。

エンジンを足す手順:
1. アダプタ(`prepare` / `parse` / `continuation` と `result_type`)とエンジンクラスを `hfauto/backends/<engine>/` に書く。Level は出力から観測し、振動数は `projected_frequencies` で求める。
2. `backends/engines.py` の表に 1 行足し、`pyproject.toml` の `adapters-private` と `adapters-independent` にモジュールを加える。
3. golden: 実出力を `tests/golden/excerpt.py` の規則で抜粋して `tests/golden/data/<engine>/` に置き、出所と sha を `tests/golden/SOURCES.json` に書き、パーサのテストを書く。
4. smoke: `tests/smoke/` に `pytest.mark.real` のテストを書く(`real_engine` fixture)。
5. `configs/sites/wsl_local.yaml` に `EngineSite` を加え、`hfauto doctor` で確かめる。

## 6. driver

### 6.1 MinimumDriver(`drivers/minimum.py`)

`relax_to_minimum` の手順:
1. opt(初期 Hessian があれば使う)。最終構造が `Registry` の既存の basin(同じ組成・電荷・多重度・level_key)に一意に割り付けば、freq を省いて `known`。
2. 別ジョブの freq → `is_minimum` と `spin_ok`。
3. saddle 級なら ν < −saddle_cm1 の虚モードに沿って `modes.amplitude`(½κs² = 3e-4 Eh、0.05〜0.4 Å)だけ変位し、変位元の freq を初期 Hessian にして opt → freq(最大 `mode_follow` サイクル)。虚モード 1 本なら ± 両側、2 本以上なら全モードの和の方向へ片側だけ押す。両側がラベル付きで別の構造なら `ts_candidate`(縮退転位を含む。minima stage が `mode_follow` の discovery にし、元の化学種は入力構造に近い側の basin に注記 `endpoint_was_saddle` 付きで加わる)、同じなら置き換え、片側だけならその側を採る。
4. soft なら ν < −noise_cm1 の虚モードの和の方向へ 1 回だけ押し、残れば `soft_minimum`。noise はそのまま `minimum`。

`Registry` は組成 × level_key ごとの basin の集合で、同一性の基準は `identity.assign` の 1 つだけである(\|ΔE\| ≤ 5e-5 Eh の候補のうち、回転と鏡映を許す置換不変 RMSD ≤ 0.05 Å で、次点と十分離れていること)。鏡像は同じ basin で、鏡像と重ならない構造の basin は `chiral`(m = 2)を持つ。ラベル付きの比較(`mapped_rmsd`)は真の回転だけなので、NH3 の反転のような縮退転位は恒等と区別される。minima stage はジョブを species id の順に直列で relax し、直後に登録するので、同じ stage の中でも登録済みの basin に落ちたジョブは freq を省く。mode-follow の経過は `<run>/<stage_id>/diagnostics.json` に書く(記録用で、読むコードはない)。

### 6.2 ReactionCaseDriver(`drivers/reaction_case/`)

`decide(case, state, rules)` は `state.py` の 17 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。判断は `<run>/<stage_id>/cases/<reaction_id>/log.jsonl` に書く(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 両端が同じ level_key の DFT 極小でない | BLOCKED(`endpoints_not_on_one_pes`) |
| 2 | 両端が同じ basin で縮退反応でない | SAME_BASIN(`same_basin`) |
| 3 | ΔE_rxn > reaction_window_kcal | OUT_OF_WINDOW(`out_of_window`) |
| 4 | 接続が elementary / degenerate / reassigned | 完了(`connection:<label>`) |
| 5 | 中間体が両端と別の basin | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(深さ `max_split_depth` まで) |
| 6 | 最新の DFT プロファイル(SCREEN か string)が barrierless | BARRIERLESS(`screen:barrierless` / `string:barrierless`) |
| 7 | 予算 `walltime_h` を使い切った(残りが 60 s 未満で、どのジョブも始められない) | UNRESOLVED(`walltime`) |
| 8 | 接続を判定したが完了していない | 両側が同じ basin で 2 回未満なら CONNECT(振幅 × 2、`connection_retry`)。\|ν\| < saddle_cm1 の TS なら行 11 へ、それ以外は UNRESOLVED(`connection_failed`) |
| 9 | 未接続の SaddleClaim がある | CONNECT(`ts_validated`) |
| 10 | saddle が収束し未検証(`ts_calc` の case は最初から) | VALIDATE_TS(`saddle_converged`) |
| 11 | TS 検証で崩壊した(最低モード ≥ −noise_cm1)か、柔らかい TS の QRC が失敗した | VALIDATE_INTERMEDIATE(`saddle_collapsed`) |
| 12 | SCREEN が有効で低レベル経路が未実行、かつプロファイルがないか種が尽きた(近道の種が失敗したら string の前に 1 回だけ NEB) | SCREEN(`screen`) |
| 13 | 最新の DFT プロファイルが intermediate で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_intermediate`、最も低い井戸) |
| 14 | 未使用の種があり、saddle の試行 < `max_saddle_attempts`(先頭が `saddle_restart` なら上限でも可) | REFINE_SADDLE(`seed:<source>`) |
| 15 | DFT string が未実行で、saddle の試行が残っている | FIND_PATH(`no_dft_path`) |
| 16 | 最新の string が single か intermediate で、種がなく、string < 3 チャンク、saddle の試行が残っている | FIND_PATH(`next_chunk`) |
| 17 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

**予算**: 行 4〜6 は手元の証拠だけで決まる完了なので行 7 より前に置き、締め切りの後に確定した接続・中間体・障壁なしを捨てない。行 8 以降は新しい計算を始めるか打ち切る。saddle の試行は case ごと(子反応は 0 から)。高次 saddle を押した種は数え、停滞の再開(1 回)は数えない。種を精密化できない string は走らせない。仮説と子反応が共有するのは `walltime_h` の締め切りだけなので、reaction-paths の最悪の所要時間は仮説数 × `walltime_h` である。

action:
- **SCREEN**: 低レベル TS(explore か mode-follow)があれば xTB freq で確かめ、その構造の DFT SP と両端の DFT 極小の 3 点が single なら種にする(`discovery_ts`)。なければ DFT 極小の間の IDPP 11 点を初期経路に `pysis_neb`(xTB の CI-NEB、端は固定、未収束でも経路として使い、CI から TSOpt)で緩和し、内部 9 点の DFT SP と両端の DFT 極小のエネルギーを `barrier_verdict` で分ける。single の種は NEB の TS(`screen_ts`)か、山を放物線補間した構造(`screen_hei`)。NEB が失敗したら IDPP のまま分ける。経路の両端は常に DFT 極小なので、低レベルの PES に端点の極小は要らない。barrierless は最も高い内部点の両隣の区間の中点 2 点の DFT SP を加えて分類し直してから受け入れる(SP が失敗したら unavailable、`midpoint_single_point`)。FIND_PATH も同じ。
- **REFINE_SADDLE**: 反応方向は 1 つの規則で決める。低レベル TS の種はその虚モード、高次 saddle を押した種はその反応モード、それ以外は経路の接線のうち反応中心(結合が変わる原子とその隣接原子)の成分。結合が変わらなければ宣言座標の勾配、なければ最も変わる二面角の勾配、それもなければ全原子の接線。初期 Hessian は、種が 0.5 Å 以内の検証済み TS freq を持てばそれ、なければ種の xTB freq で負モードのどれかと方向の重なりが 0.3 以上ならそれ、そうでなければ種の DFT freq。adapter はその Hessian を、方向と最も重なる固有ベクトルだけが負になるように整えて渡す(`vibrations.shape_hessian`)。saddle が maxiter で止まったら、最終フレームを種(`saddle_restart`)にして新しい Hessian で 1 回だけやり直す。
- **VALIDATE_TS**: 別ジョブの DFT freq → `is_first_order_saddle`。`ts_calc`(mode-follow の鞍点か親が検証した TS)はまずそれを検証し(freq は JobStore の再利用)、ゲートを通らなければ SCREEN へ。higher_order なら、方向と最も重なる負モード以外で最も負のモードに沿って片側に押した構造を種(`higher_order_retry`)にし、検証済みの freq をその Hessian にする。
- **FIND_PATH**: ZTS を 1 チャンク(9 beads、maxiter 20)だけ走らせ、収束を問わず `barrier_verdict` で分ける。初期経路は最新の DFT プロファイルの経路、なければ IDPP。端点は DFT 極小にし、画像を隣へ逐次整列する。種は single のときだけ山の放物線補間(`path_hei`)。種が尽きれば行 16 がその経路から次のチャンクを続ける。
- **CONNECT(QRC)**: TS の虚モードに沿って ± に変位し(振幅はエネルギー目標 max(3 × qrc_drop, 3e-4 Eh) と ν・質量から 0.05〜0.4 Å。再試行は 2 倍で 0.4 Å で切る)、TS の freq Hessian で opt し、`Registry` に割り付け(未知なら `relax_to_minimum` で新しい basin)、`connection` で判定する。reassigned のうち片側だけが端点の basin で他方が別の DFT basin なら行 5 で分割し、その 2 つを結ぶ子反応に TS を渡す(`ts_calc`)。どちらの端点も含まない TS は reassigned のまま。
- **VALIDATE_INTERMEDIATE**: 崩壊した saddle か最新のプロファイルの最も低い井戸を `relax_to_minimum` にかける。両端と別の basin なら行 5 で分割する。緩和が失敗したら結果とせず(端点扱いも barrierless もしない)、最も高い山を種にする。端点に落ちたら、内部の最大 − 高い方の端点 < resolution_kcal なら barrierless、そうでなければ最も高い山を種にする。

## 7. 化学プロトコル

### 7.1 stage ごとの既定値

| 段 | 既定値と規則 |
|---|---|
| 共通 | DFT は PBE0-D3BJ/def2-SVPD、grid fine、`convergence energy` 1e-7(opt・freq・SP で同じ)。反復上限は DFT `iterations 100`、WFT の SCF `maxiter 100`。SCF 未収束の救済は閉殻 DFT が `cgmin`、WFT が前回の vectors からの再開で、1 回だけ(smear・fon は PES を変えるので使わない)。開殻 DFT は cgmin が ⟨S²⟩ を出さず Evidence にならないので救済しない。geometry は `units angstrom nocenter noautosym`、autoz が失敗したら Cartesian で再投入。Z>36 の元素には def2 系の基底のときだけ `<元素> library def2-ecp` を元素ごとに書く |
| structures | 化学種は xyz と SMILES のちょうど一方。SMILES は RDKit ETKDG の 1 配座(同位体と `'.'` は不可)。宣言反応の端点は xyz(SMILES では原子の対応と配座が決まらない)で、両端の組成・電荷・多重度と原子順序が一致すること。元素・電荷・多重度はここで 1 回だけ検査する(範囲外の元素は `unsupported_element`、d ブロック元素を含む化学種かスピン結合で多重度が 1 つに決まらない組成で未宣言なら `declare_multiplicity`、電子数とのパリティ違いは `INPUT_INVALID`)。組成の電荷は成分の和、多重度は宣言値かスピン結合で 1 つに決まる値。開殻の成分から多重度 1 を作る宣言は `open_shell_singlet_unsupported` |
| conformers | 単量体: `crest --gfn2 --quick -T <n> --ewin 6 --chrg --uhf`(重原子 3 個以下で回転可能結合のない分子は省く)。トポロジー変化で止まったら停止構造を残し、そこから 1 回だけ再実行する。組成: 受容原子の lone-pair 円錐に極性 H を置く seed(donor / 受容原子がなければ vdW 接触 + 0.5 Å の剛体配置)を `seeds_per_composition` 個作り、先頭の seed から `--nci --quick --notopo <全原子> --noopt` で探索する(状態は hfauto の状態ラベルが決め、CREST はサンプリングだけを行う。`--noopt` は CREST 3.0.2 の初期トポロジー検査が `--notopo` を見ないため)。CREST が失敗した組成(開殻など)は seed をそのまま出す。選抜は状態ラベルごとに GFN2 エネルギーの低い `keep_per_state` 個。CREST は attempt ディレクトリで `--scratch` なしに実行する |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル。状態ラベルが変わった seed は members に記録し、DFT には流さない |
| explore | 出発点は組成 × 状態ラベルごとの screen 最低 `sources_per_state` 個。trial は元素に依らない列挙器 1 つで作る: 形成の候補対は結合しておらず r ≤ Σr_vdW でグラフ上 3 結合以上離れた対。T1 移動・置換(形成 1 + 切断 1)、T2 リレー(H 移動 2 つの連鎖)、T3 形成だけ、T4 切断(開殻か電荷系だけ)を順に巡回し、出発点あたり `max_trials_per_source` 件。変化後の結合数が元素の最大配位数を超える drive は作らず、原子クラスと距離が同じ drive は 1 件にする。各 trial で NT2、極大がなければ同じ drive で AFIR(γ 125 → 300 kJ/mol)。単位は並列に走らせ、入力順に記録する。採否は `discovery_verdict`。残した生成物は同じ組成の screen 極小と既出の生成物に状態ラベルと置換不変 RMSD で 1 回だけ同定し、新しい basin だけを species にする(縮退の生成物は添字付きの結合の一致も要る)。ReaDuct の `spin_mode` は `restricted_open_shell` |
| minima(dft) | 選択 `window`: 組成 × 状態ラベルごとに `window_kcal` 以内の `per_state` 構造 + explore の生成物 + 入力単量体の最低構造(`rerank_sp` なら上位 `rerank_top` 構造を DFT SP で並べ直してから同じ窓で切る)。`all` は入力の化学種全部。opt は初期 Hessian なしなら `trust 0.1`、ありなら `trust 0.3`、maxiter 100。2 断片以上は `init_hessian` の xTB Hessian を使う |
| 仮説 | 優先順は宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補。宣言反応は必ず評価する。それ以外は同じレベルの DFT 極小の対(または 1 つの basin)で、ΔE_rxn ≤ reaction_window_kcal、端点の間で結合が変わるものだけ(組成あたり 6 件まで)。未宣言のねじれ・配座変化・鏡像化は仮説にしない(Curtin–Hammett)。両端が同じ basin の発見は、結合の相手が入れ替わるときだけ縮退反応にする。DFT 段の mode-follow の鞍点は `ts_calc` として同じ basin の対(宣言反応を含む)に貸す。縮退・結合変化・ねじれは basin の最適化構造を端点の原子順と掌性に並べた構造で判定する。explore の陰性結果は仮説を棄却しない |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`、`moddir 1`(整えた Hessian の負のモード)。maxiter は継続せず、最終フレームを持つ `Failure` を返す。ZTS は `nbeads 9`、`stepsize 0.05`、`interpol 3`、`freeze1` / `freezeN`、maxiter 20 のチャンクで、NWChem の収束判定は使わない。`pysis_neb` は IDPP を初期経路とする Cartesian の CI-NEB(`opt: lbfgs`、max_cycles 100)と rsprfo の TSOpt。QRC の振幅の上限 0.4 Å は初期 Hessian を受け付ける距離 0.5 Å 以下にする |
| sp | `methods` の各 LOT で、順位に使う点だけを計算する: 順位を付けられる outcome の反応の TS、その反応物・生成物と同じ状態の DFT 極小すべて、組成の単量体の状態の DFT 極小。freq の最終構造で計算し、parents に対象(minimum_id か TS の freq calc)を書く。thermo と手法パネルはこの parents だけで SP と停留点を対応付ける。CCSD(T) は閉殻が RHF の ccsd モジュール、開殻が ROHF 参照の TCE(`2eorb 2emet 13`)。`freeze atomic` は ECP 原子の軌道を凍結しない |
| thermo | GoodVibes 4.3.0 に自前の振動数を渡す。H は RRHO、S は Grimme の qRRHO(`qs`、`cutoff_cm1` 100)。対称数は pymsym(libmsym)の点群から、等価判定のしきい値を既定の 5e-4 から 2e-3 に緩めて求める(最適化の終点のわずかなずれで C3v が Cs になるのを防ぐ)。スケール因子は `vib_scale` の 1 つ(既定 1.0、振動数と ZPE の両方)。負モードは固定の規則: 極小は全モードを \|ν\| に、TS は最低モードを除いて \|ν\| にする。キラルな極小と TS の G に −RT ln 2 を加える(鏡像対は 1 つの basin で m = 2。σ の比は対称数で入る)。種の値は 1 atm で、反応は `standard_states` ごとに換算する。状態の G は、その状態の極小のうち錯体(または基準)と同じ LOT で spin_contaminated でないものの最小の G で、反応物・生成物と会合量の単量体で同じ定義(LOT 判定は `same_pes(state=False)`)。感度の幅は qs × cutoff 50/100/150。G がなければ `thermo_unavailable`、LOT の不一致は `mixed_level_of_theory`。トンネル補正はない |
| report | rankable な反応を δG_eff で並べ、感度の幅が重なれば同順位。(T, 標準状態) は report の `T_K`・`standard_state`、なければ thermo の最初の組。ranking.csv の列は rank、reaction_id、outcome、tier、rankable、T_K、standard_state、energy_level、dG_eff_kcal、band_low_kcal、band_high_kcal、dG_act_kcal、dG_rxn_kcal、dG_act_vs_separated_kcal、torsional(結合変化のない段)、blockers、notes。エネルギーが 2 つ以上の LOT にあれば `dE_act_panel_min_kcal`・`dE_act_panel_max_kcal` を足す(手法の幅は列で示し、順位は止めない)。coverage.csv は機構別の試行・生成物・陰性理由・FailureKind |

結合と状態(`chemistry/topology.py`): 結合は r < r_cov,i + r_cov,j + 0.4 Å(Cordero の共有結合半径)で 1 つの構造だけから決め、`bond_changes` は結合集合の差で向きに対称である。状態ラベルは断片の組成式と WL ハッシュ。FHF⁻ と I3⁻ は 1 断片、ハロゲン結合(I···N 2.8 Å)は非結合。イオン–双極子錯体は許容値の内側なら 1 断片になり、GFN2 の強い H 結合錯体は閾値の近くでラベルが分かれうる(DFT 極小には影響しない)。

### 7.2 順位の量 δG_eff

δG_eff = max(G_TS, G_R, G_P) − G_R(`thermo.effective_barrier`)。G_R・G_P は反応物・生成物と同じ状態(組成と状態ラベル)で、freq とエネルギーの LOT が同じ DFT 極小のうち最小の G である(速い配座平衡の下で反応物の状態は 1 つ。両端が同じ状態の配座変化は、どちらの向きも状態の最低配座から測る)。TS が反応物か生成物より下にあれば律速ではない(微視的可逆性)。順方向か逆方向の ΔE0‡ = ΔE‡ + ΔZPE‡ が 0 以下なら鞍点の TST は意味を持たないので TS を除いて max(ΔG_rxn, 0) とし、注記 `submerged_barrier` を付ける(blocker ではない)。barrierless_at_resolution も max(ΔG_rxn, 0) で、解像度未満の障壁であることは outcome 列で分かる。`dG_act_kcal` はその反応の反応物極小から見た TST の値で、表示だけに使う。分離反応物の G は基準に混ぜず、`dG_act_vs_separated_kcal` の列に出す。反応はすべて Δn = 0 なので δG_eff は標準状態に依らない。

### 7.3 エネルギー層と手法パネル

既定の順位は停留点レベル(PBE0-D3BJ/def2-SVPD)の G で付ける。エネルギー層(composite G = E_SP + (G_GV − E_GV))は、既存の run に `method_panel` を `--run-dir` で追記して使う: panel_sp → panel_thermo(`energy_method` = 参照手法)→ panel_report。thermo の結果の id は stage id を含まないので、追記した層が元の thermo の値を上書きし、ranking.csv の `energy_level` 列に並べた層が出る。元の値を残したいときは run ディレクトリを複製してから追記する(`cp -a`)。SP の対象は単量体を含むので ΔG_assoc も同じ層で出る。エネルギー層の SP がない対象は `energy_layer_missing`(G = None)で、その SP にも `spin_ok` をかける。

- ωB97X-D3/def2-TZVPD は閉殻 4 系中 3 系で CCSD(T) に近い(平均絶対誤差 2.1 → 1.3 kcal/mol)が、17 原子の SP が 1 点 29 分かかるので既定の pipeline には入れない。
- CCSD(T)/def2-TZVPD は重原子約 4 個までの分子で `methods` に足して `energy_method` にする(5 個では (T) が 1,200 MB/rank で配列を確保できない)。多参照性のゲートはない。

### 7.4 結果の読み方

- 停留点レベルを def2-SVPD にしたのは、FHF⁻ の De が拡散関数のない SVP では 21.7 kcal/mol 過大になるが SVPD では実験値の誤差内に入るからである。
- PBE0 の障壁の系統誤差は反応クラスで違う(平均符号付き誤差: H 移動 −4.2、重原子移動 −6.6、求核置換 −1.9 kcal/mol。Zhao & Truhlar, J. Phys. Chem. A 109, 2012 (2005))。CCSD(T)/def2-TZVPD(SVPD 構造での SP)との差の実測は HCN −1.2、HONO +2.0、SN2 −3.0 kcal/mol。順位は δG_eff の序数として読み、反応クラスが混ざるときは数 kcal/mol 未満の差の順序は入れ替わり得る。
- SCREEN・FIND_PATH の barrierless は停留点レベルの PES 上の判定で、エネルギー層では確かめない。
- 電荷系では PBE0 の非局在化誤差で電荷の広がった TS・錯体が低く出やすく、D3 は電荷に依らない。電荷系の障壁と ΔG_assoc はエネルギー層で確かめる。
- ΔG_assoc は BSSE を補正しない値で、結合を過大に見積もることがある。qRRHO の扱いによる幅は 298 K で約 1.2 kcal/mol(TMA·(HF)₂)。
- 電子分配関数は g = 2S+1 だけで、開殻原子と ²Π ラジカルの軌道縮重とスピン軌道補正は含めない。原子の SMILES は高スピン仮定なので、基底状態と違うことがある([C] は五重項になる)。

## 8. 実行基盤

- **JobStore**(`execution/jobstore.py`): run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、結果を変える入力(method、構造の指紋、パラメータ、入力ファイルの sha)をすべて含み、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合する。終端的な失敗も記録し、`--retry-failed` で指定した種類だけを再実行する。キーにコードの版数は入らないので、パーサを直した後は `--retry-failed incomplete_output` で取り直す。
- **ladder**(`execution/jobs.py` の `LADDER`。これ以外の自動再試行はない): timeout と geometry_maxiter は最新構造からの継続を 2 回まで(saddle の maxiter は継続しない)、input_invalid(autoz)は Cartesian で 1 回、scf_not_converged は閉殻 DFT の `cgmin` か WFT の再開で 1 回(開殻 DFT はなし)。attempt の timeout は min(設定値, 予算の残り) で、残りが 60 s 未満なら `budget_exhausted`。
- **失敗の閉じ込め**: `adapter.parse` の例外はそのジョブの `Failure(incomplete_output, 'parse:<型>: <1 行目>')` になる。反応ケースの例外はそのケースを unresolved(`error:<型>`)にし、log.jsonl に 1 行書いてほかのケースを続ける。消費する型の入力がない stage は 0 artifact で done になり、上流の失敗を warning に出して後続の stage(report を含む)へ進む。`HFAUTO_STRICT=1`(テスト専用、`tests/conftest.py`)は前 2 つを再送出にする。
- **並行実行**: `JobRunner` のセマフォで実行中のジョブの ranks × threads の合計を site の `cores` 以下に保つ。同時に走るジョブの数はこれだけで決まる。極小(xTB と DFT)と反応ケースは直列(Registry の決定性のため)、CREST と explore の単位は `thread_map` で並列(結果は入力順)。`SiteLock`(`<scratch_root>/.hfauto_site.lock`)は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。
- **シグナル**: `run_command` は外部プログラムを新しいセッションで起動して登録する。`hfauto run` は SIGINT・SIGTERM・SIGHUP を受けると(POSIX)登録済みのプロセスグループをすべて止め、実行中の stage を failed と記録して 128 + signum で終わる。止めたジョブは JobStore に残らないので、再開すると取り直す。
- **run ディレクトリ**: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/cases/`、`<run>/run_state.json`(実行順の stage 状態とジョブの統計)、`<run>/resolved_config.yaml`(設定と code_version = git の sha、差分があれば -dirty)。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に合わせたもの(別の pipeline を同じ run に追記できる)。done で入力と設定の sha が変わらない stage は飛ばし(再開)、`--from X` は X 以降を取り直す。

## 9. 設定

4 つの層は中身が重ならないので、マージも優先順位もない。結果を変える設定は method・system・pipeline に置き(ジョブの鍵か stage の設定 sha に入る)、site には実行設定だけを置く。system と pipeline の YAML は未知のキーを実行前に拒否する。

| 層 | 置き場所 | キー |
|---|---|---|
| site | `configs/sites/` | `site`、`scratch_root`、`cores`、`engines.<registry 名>`: `version`(版数の pin)、`executables`、`execution`(`ranks`、`threads`、`memory_mb_per_rank`、`timeout_s`、`env`)、`python`(worker)、`scratch_dir`。絶対パスはこの層だけに書く |
| method | `configs/methods/` | `id`(ファイル名と同じ)、`kind`(dft / xtb / wft)、`functional`、`wft_method`(ccsd(t))、`basis`、`dispersion`(d3zero / d3bj)、`gfn`、`electronic_temperature_K`、`grid`、`scf_energy_tol`。gfn2(低レベル)、pbe0-d3bj_def2-svpd(停留点)、pbe0-d3bj_def2-tzvpd と wb97x-d3_def2-tzvpd(手法パネル)、ccsd-t_def2-tzvpd(opt-in) |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | `system_id`、`species`(`id`、`xyz` か `smiles`、`charge`、`multiplicity`、`role`: monomer / endpoint)、`compositions`(`id`、`components`、`multiplicity`)、`reactions`(`id`、`reactant`、`product`、`coordinate`、`torsional`) |
| pipeline | `configs/pipelines/` | `pipeline_id`、`gates`(knob)、`stages`(`id`、`stage` と下の stage 表の設定キー) |

**knob は 8 つ**だけである(検証は既定値で行っている)。

| knob | 置き場所 | 既定値 | 意味 |
|---|---|---|---|
| `noise_cm1` | `gates:` | 10 | これより小さい虚振動(cm⁻¹)は雑音 |
| `saddle_cm1` | `gates:` | 50 | これより大きい虚振動(cm⁻¹)は鞍点級 |
| `resolution_kcal` | `gates:` | 1.0 | 経路の山と井戸を数える深さ(kcal/mol) |
| `reaction_window_kcal` | `gates:` | 40 | 発見・仮説・判断表(行 3)の ΔE_rxn の上限(kcal/mol) |
| `spin_tol` | `gates:` | 0.1 | \|⟨S²⟩ − S(S+1)\| の許容幅 |
| `walltime_h` | reaction-paths の `policy:` | 6 | 1 仮説(分割した子を含む)の予算 |
| `max_saddle_attempts` | reaction-paths の `policy:` | 2 | 1 case の saddle の試行数(子反応は 0 から、再開は数えない) |
| `max_split_depth` | reaction-paths の `policy:` | 2 | 多段反応の分割の深さ |

- 停留点の method は minima(dft)と reaction-paths の 2 か所に書き、食い違えば読み込み時に拒否する。手法を変えるときは別の method ファイルを使う。
- method は pipeline ファイルの隣の `methods/<id>.yaml`、なければ作業ディレクトリの `configs/methods/<id>.yaml` を読む。
- SCREEN の有無は reaction-paths の `screen:` の有無で決まる。string の bead 数とチャンク数、NEB の画像数、QRC の振幅はモジュール定数である。
- 温度と標準状態は thermo の `temperatures_K`(既定 [298.15])と `standard_states`(1atm / 1bar / 1M、既定 [1atm])だけで指定する。

stage 表(次のコマンドの出力。`hfauto` を入れた Python で実行する):

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
| `conformers` | `hfauto.stages.conformer_search:ConformersStage` | species | species | engine*, method*, quick, ewin_kcal, seeds_per_composition, keep_per_state |
| `minima` | `hfauto.stages.minima:MinimaStage` | species | calculation, minimum, species, discovery | level*, engine*, method*, select, init_hessian, mode_follow |
| `explore` | `hfauto.stages.explore:ExploreStage` | minimum, species | discovery, species | engine*, method*, sources_per_state, max_trials_per_source, settings |
| `reaction-paths` | `hfauto.stages.reaction_paths:ReactionPathsStage` | minimum, species | reaction, calculation, minimum, species | method*, engines*, screen, policy, reaction_ids |
| `sp` | `hfauto.stages.single_point:SinglePointStage` | minimum, reaction | calculation | engine*, methods* |
| `thermo` | `hfauto.stages.thermochemistry:ThermoStage` | minimum, calculation | species_thermo, reaction_thermo | settings, energy_method, temperatures_K, standard_states |
| `report` | `hfauto.stages.report:ReportStage` | — | report | T_K, standard_state |

入れ子のキー: minima の `select`(`include`: window / all、`per_state` 3、`window_kcal` 6、`rerank_sp`、`rerank_top` 8)と `init_hessian`(`engine`、`method`)、explore の `settings`(ReaDuct の `max_scf_iterations`、`electronic_temperature_K`、`scc_retry_temperature_K`、`afir_gamma_kj_mol`、`afir_gamma_retry_kj_mol`、`imag_cutoff_cm1`、`timeout_s`)、reaction-paths の `engines`(`qm`、`saddle`、`path`)・`screen`(`method`、`qm`、`path`)・`policy`(上の knob)、thermo の `settings`(`qs`、`cutoff_cm1`、`vib_scale`、`sensitivity`)。
