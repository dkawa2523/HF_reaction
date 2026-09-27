# hfauto の設計(現状)

- 本書は現在のコード(`hfauto/`)の説明である。計画時の設計書は `docs/reviews/2026-09-25_refactor_design.md`(凍結)にあり、コードの docstring にある「design §x.y」はその節番号を指す。
- 原則: stage は薄いアダプタで、化学の判定は `hfauto/chemistry/` の純関数だけが行う。外部ジョブの観測事実は型付きの `Evidence` で受け渡し、合否はゲート関数 1 か所で決める。バックエンドは能力 Protocol の実装で、subprocess を起動するのは `hfauto/execution/process.py` だけである。fail-closed で、互換層・ダミーエンジン・内部フォールバックは持たない。エンジン名の文字列は `hfauto/backends/` と `configs/` の外に書かない。

## 1. 層構造と import 契約

依存の向きは `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core` とする。複数の層で使う型は、使う層のうち最も下の層に置く。

| 層 | パッケージ(行数) | 役割 |
|---|---|---|
| cli | `hfauto/cli/`(186) | 5 コマンド: `doctor`、`run`、`status`、`report`、`case` |
| pipeline | `hfauto/pipeline/`(834) | 設定の読み込み、preflight、`RunLayout`、`Runtime`(`StageRuntime` の実装)、`run_pipeline` |
| stages | `hfauto/stages/`(1,361) | 8 つの stage。入力を得て driver か純関数を呼び、artifact を返すだけ |
| drivers | `hfauto/drivers/`(1,138) | `MinimumDriver`(`hfauto/drivers/minimum.py`)と `ReactionCaseDriver`(`hfauto/drivers/reaction_case/`) |
| reporting | `hfauto/reporting/`(520) | 順位・被覆率・手法パネルの表(`hfauto/reporting/summary.py`)と静的 HTML(`hfauto/reporting/html.py`) |
| backends | `hfauto/backends/`(2,460) | 能力 Protocol(`hfauto/backends/protocols.py`)、registry(`hfauto/backends/engines.py`)、6 つのアダプタ |
| execution | `hfauto/execution/`(744) | `run_command`、`JobStore`、`JobRunner`、`SiteLock`、worker |
| chemistry | `hfauto/chemistry/`(2,582) | ゲートと化学カーネル(純関数。ファイル IO は `hfauto/chemistry/xyz.py` と `hfauto/chemistry/xyz_trajectory.py` だけ) |
| core | `hfauto/core/`(779) | `Evidence`、payload、manifest、`MethodSpec`、system のモデル。依存は pydantic、numpy、PyYAML だけ |

import-linter の契約(`pyproject.toml`、`lint-imports` で検査。5 契約とも KEPT):
1. `hfauto-layers`: 上の層構造(exhaustive)。
2. `adapters-private`: cli・pipeline・stages・drivers・reporting はアダプタ(nwchem、xtb、pysis、crest、readuct、goodvibes)を import しない。
3. `protocols-only`: stages・drivers・reporting は `hfauto.backends.engines` と `hfauto.execution` を直接 import しない。エンジンは `StageRuntime.engine()` から Protocol 型で受け取る。
4. `adapters-independent`: 6 つのアダプタは互いに import しない。
5. `stages-independent`: 8 つの stage モジュールは互いに import しない。

ruff は `C901`(max-complexity 15)と `TID251`(`subprocess` は `hfauto/execution/process.py` だけ)を有効にしている。型検査は `pyrefly check`(0 errors)。

## 2. Evidence とゲート

`Evidence`(`hfauto/core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すること自体が次を保証する: (1) 正常終了と SCF 収束、(2) opt / saddle なら構造最適化の収束、(3) freq なら 3N − n_external 本の振動数(原子は n_external 3 で振動なし)、(4) 原子順序の保持、(5) 入力フレームの保持(出力に出た初期構造が入力と 1e-4 Å 以内)、(6) 出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること(`level_mismatches`)。満たさなければアダプタは `Failure` を返す。振動数と正準形式の Hessian(`.npy`、Eh/bohr²、入力フレーム)は、どのエンジンでも `hfauto/chemistry/vibrations.py` の補空間射影(並進・回転空間を SVD で決め、主同位体の質量を使う)で求める。経路は `PathProfile`、失敗は `Failure(kind: FailureKind, reason)` で表す(maxiter で止まった saddle は最終フレーム `final` も持つ)。FailureKind は executable_missing、input_invalid、timeout、nonzero_exit、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、budget_exhausted、gate_rejected の 10 種である。

ゲートは `hfauto/chemistry/gates.py` にだけ置く。閾値は `Policy` にまとめ、pipeline YAML の `gates:` だけが上書きできる。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=True, state=True)` | program・version・method・basis・dispersion・solvation・charge・multiplicity・電子温度が一致。`numerics=True` なら grid と scf_tol も。停留点・freq・QRC・熱化学は True、DFT string(種にだけ使う)は False。`state=False` は電荷と多重度を比べず、会合量の錯体と単量体の LOT 判定だけで使う |
| `spin_ok` | `s2` が None か、\|<S²> − S(S+1)\| ≤ 0.1。freq と、thermo が使うエネルギー層の SP(UKS)にかける(ROHF-CCSD(T) には <S²> がない)。不合格は blocker `spin_contaminated` |
| `imaginary_tier` | ν < −50 は saddle、−50 ≤ ν < −10 は soft、−10 ≤ ν < 0 は noise |
| `is_minimum(freq, opt=)` | task が freq / opt、`freq.start.fingerprint == opt.final.fingerprint`、numerics を含めて同一 PES、本数が 3N − n_external、tier が saddle でない(soft / noise は注記)。同じ幾何で \|E_freq − E_親\| ≤ qrc_drop(max(1e-5, 20×scf_tol))、外れたら `state_mismatch` |
| `is_first_order_saddle(freq, saddle=)` | TS の定義(負の固有値がちょうど 1 本、大きさは問わない)だけを見る: 指紋と同一 PES、本数、最低モード < −10(noise_cm1。なければ `no_imaginary_mode`)、2 本目 < −50 なら `higher_order`、2 本目が −50〜−10 なら注記 `soft_secondary_mode`。宣言端点とのエネルギー比較、ねじれ用の閾値、モードの重なりは合否に使わない。同じ幾何で \|E_freq − E_親\| ≤ qrc_drop、外れたら `state_mismatch` |
| `barrier_verdict` | 両端が 2 つの DFT 極小のエネルギーであるプロファイルを `profile.classify` で分ける(解像度 1.0 kcal/mol)。山(内部極大)と井戸(内部極小)を同じ関数 `interior_maxima` で数え、井戸があれば `intermediate`(山が 2 つあれば間に必ず井戸がある)、井戸がなく山があれば `single`、どちらもなければ `barrierless`。3 点未満は `unavailable`。記録は verdict、source(screen / string)、内部の最大 − 高い方の端点、ノード間隔(barrierless の解像度)、理由だけ |
| `connection` | 接続の定義(両側が別々の極小へ下る)だけを見る: 各側の最終エネルギー < E_TS − drop(drop = `qrc_drop(level)` = max(1e-5 Eh, 20 × scf_tol)。CONNECT の振幅もこの値から決める)で同一 PES。途中の軌跡は見ない(直線変位の出発点が TS より上でもよい)。割付けの集合が期待どおりなら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつ両側の構造がラベル付きで別(`mapped_equivalent`。鏡像化を含む)なら degenerate、それ以外で両側が同じ basin なら failed(`sides_same_basin`)。結合が変わる縮退反応では QRC 両側の原子番号付き結合グラフの組が {bonds(R), bonds(P)} と一致しなければ failed(`bond_change_missing`) |
| `thermo_consistent` | GoodVibes の実振動の本数、ZPE(許容 1e-5 Eh)、E(許容 1e-6 Eh)が、渡した振動数(`thermo_frequencies`、`vib_scale` 倍)と freq の E に一致 |
| `rankable` | outcome が elementary / degenerate / reassigned / barrierless_at_resolution で、ReactionThermo に blocker がない。blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated)を作るのは thermo stage だけで、rankable が足すのは `outcome:<o>` と、記録そのものがないときの `thermo_record_missing` だけ。順位の量でないもの(ΔZPE‡ の大きさ、手法パネルの ΔE_rxn の符号や幅)では順位を止めない |
| `reaction_tier` | connected(ConnectionClaim)> saddle(SaddleClaim)> minima(driver が BLOCKED 以外で終えた)> screening |

claim を作るのは 1 か所だけである。`MinimumRecord` は minima stage と、reaction-paths で新しい basin を見つけたときの `Registry` が作る。`SaddleClaim` と `ConnectionClaim` は reaction-paths の driver だけが作る。thermo と report は claim の有無だけを見て再検証しない(例外は thermo の `same_pes` と `thermo_consistent`)。

## 3. payload(artifact の型 8 種)

各 stage の `<run>/<stage_id>/manifest.json` には、その stage の出力だけを入れる。artifact は `status: success | failed` で、payload は `kind` で判別する frozen なタグ付き共用体(`hfauto/core/records.py`、`extra="forbid"`)である。

| 型 | レコード | 出す stage | 主な中身 |
|---|---|---|---|
| species | `SpeciesRecord` | structures、conformers、minima、explore、reaction-paths | 組成キー、電荷、多重度、`Geometry`、source、状態ラベル |
| calculation | `Evidence` | minima、reaction-paths、sp | 1 ジョブの観測事実。id は `calc_` + job_key の先頭 16 桁 |
| minimum | `MinimumRecord` | minima、reaction-paths | basin、tier(screen / dft)、level_key、opt / freq の calc id、members、chiral(鏡像が別の構造、m = 2)、notes |
| discovery | `DiscoveryRecord` | explore、minima(mode_follow) | 機構、trial、product / negative / failed と理由、低レベル TS |
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、source、`BarrierVerdict`(barrierless / single / intermediate / unavailable)、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome`(elementary_step、degenerate_rearrangement、reassigned_step、multi_step、barrierless_at_resolution、same_basin、out_of_window、unresolved_within_budget、blocked_upstream) |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None)、注記 |
| reaction_thermo | `ReactionThermo` | thermo | 順位の量 δG_eff とその感度の幅、ΔE‡、ΔE_rxn、ΔZPE‡、ΔG‡、ΔG_rxn、ΔG_assoc、単量体基準の ΔG‡、blockers、エネルギー層のラベル(`energy_level`、例 `wb97x-d3/def2-tzvpd`)、注記(`submerged_barrier`) |
| report | `ReportRecord` | report | 順位の行と表(ranking.csv、coverage.csv、method_panel.csv、report.html) |

id の規則: species は `species_<species_id>`、minimum は `min_<species>_<level_key[:8]>`、失敗した極小は `min_<species>_<stage_id>`、分割した子反応は `<parent>_split<n>`、report は `report_<stage_id>`。

## 4. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス(`hfauto/backends/engines.py`) | 結果 |
|---|---|---|---|
| qm | `energy`、`optimize(init_hessian)`、`frequencies` | `nwchem` → `NWChemEngine`(DFT と、MP2 / CCSD(T) の SP。開殻の CCSD(T) は ROHF 参照の TCE)、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path)`(初期経路は必須) | `nwchem_string` → `NWChemString`(ZTS 1 チャンク)、`pysis_neb` → `PysisNEB`(xTB ネイティブの固定端 CI-NEB と CI からの TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2 / AFIR、worker で実行) | `DiscoveryResult` |
| thermo | `thermo(freq, settings, temperatures_K, saddle)` | `goodvibes` → `GoodVibesEngine`(4.3.0 の API、worker で実行) | `ThermoResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作り、`requirements()` は classmethod とする。`engines.create` は結果が能力の Protocol を満たさなければ `TypeError`、表にないキーは `KeyError` にする。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult` / `odft`、xTB と CREST の `--chrg` / `--uhf`、pysis の calc、ReaDuct の計算器設定)。
- 初期 Hessian は freq `Evidence` として渡す(Level は問わない)。saddle は種と同じ構造(指紋が一致)、optimize は原子順序とフレームが同じで各原子が 0.5 Å 以内の構造の Hessian を受け付け(QRC の両側は TS の Hessian を使う)、外れれば `INPUT_INVALID("hessian_geometry_mismatch")`。xTB は初期 Hessian を使わない。
- SCINE、pysisyphus、GoodVibes は `python -m hfauto.execution.worker <module:function> <job.json>` の子プロセスの中だけで import する。

エンジンを追加する手順:
1. アダプタ(`prepare` / `parse` / `continuation` と `result_type`)とエンジンクラスを `hfauto/backends/<engine>/` に書く。Level は出力から観測し、`level_mismatches` があれば `METHOD_MISMATCH`、振動数は `projected_frequencies` で求める。
2. `hfauto/backends/engines.py` の `_TABLE` に 1 行足し、`pyproject.toml` の `adapters-private` と `adapters-independent` にモジュールを加える。
3. **golden(必須)**: 実出力を `tests/golden/excerpt.py` の規則で抜粋して `tests/golden/data/<engine>/` に置き、`tests/golden/SOURCES.json` に出所と sha を記録し、パーサのテストを `tests/golden/` に書く。
4. **smoke(必須)**: `tests/smoke/` に `pytest.mark.real` のテストを書き、`real_engine` fixture で WSL の実エンジンを 1 回ずつ動かす(`docs/environment.md`)。
5. `configs/sites/wsl_local.yaml` に `EngineSite`(版数の pin、実行ファイル、`ExecutionSpec`)を加え、`hfauto doctor` で確認する。

## 5. MinimumDriver(`hfauto/drivers/minimum.py`)

`relax_to_minimum(mol, method, qm, *, known, init_hessian, max_mode_follow, gates, deadline, load_xyz)` の手順:
1. opt(`init_hessian` があれば使う)。`known` の Registry で最終構造が既存の極小(同じ組成・電荷・多重度・level_key の basin だけが候補)に一意に割り付けば、freq を省いて `known`。
2. 別ジョブの freq → `is_minimum` と `spin_ok`。
3. saddle 級なら ν < −saddle_cm1 の虚モードに沿って `modes.amplitude`(½κs² = 3e-4 Eh、κ = ω²·Σm\|u\|²、0.05〜0.4 Å)だけ変位し、変位元の freq の Hessian をそのまま初期 Hessian(trust 0.3)にして opt → freq(最大 2 サイクル)。虚モード 1 本なら ± 両側、2 本以上なら全モードの和の方向へ片側 1 回だけで `ts_candidate` は作らない。両側がラベル付きで同じ構造(`same_basin` かつ `mapped_equivalent` でない)なら置き換え、別々の構造なら `ts_candidate`(NH3 の反転・鏡像化などの縮退転位を含む。minima stage が `mode_follow` の discovery にし、親の species は入力構造に置換不変 RMSD で近い側(0.05 Å 以内なら低い方)の basin に注記 `endpoint_was_saddle` 付きで加わる)、片側だけならその側を採る。
4. soft なら ν < −noise_cm1 の虚モードの和の方向へ同じ振幅規則で 1 回だけ押し、残れば `soft_minimum`。noise はそのまま `minimum`。

`Registry` は組成 × level_key ごとの basin の集合で、同一性の基準は `identity.assign` の 1 つだけである(\|ΔE\| ≤ 5e-5 Eh の候補に先に絞り、真の回転と鏡映(det = ±1)を許す置換不変 RMSD ≤ 0.05 Å で、次点と十分離れていること)。鏡像は同じ basin で、鏡像と重ならない構造(鏡像との置換不変 RMSD > 0.05 Å、`identity.is_chiral`)の basin は `MinimumRecord.chiral` を持つ。ラベル付きの比較(`mapped_rmsd`)は真の回転だけなので、NH3 の反転のような縮退転位は恒等と区別される。`find` は Evidence の Level から組成(電荷・多重度)と level_key を作って必ず照合し、`add` は known の basin か `find` が割り付く basin に加え、なければ新しい basin を作る。minima stage はジョブを species id の順に直列で relax し、relax の直後に登録する(mode-follow の両側は親の直後)。このため同じ stage の中でも、登録済みの basin に落ちたジョブは `known` になり、freq は basin あたり 1 本で済む。

振動数は質量加重の Eckart 射影(`hfauto/chemistry/vibrations.py`)で求めるので、NWChem が表示する Projected Frequencies とは低振動モードで数 cm⁻¹ ずれる(G では約 0.02 kcal/mol)。

## 6. ReactionCaseDriver(`hfauto/drivers/reaction_case/`)

`decide(case, state, policy)` は `hfauto/drivers/reaction_case/state.py` の 16 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。`<run>/<stage_id>/cases/<reaction_id>/log.jsonl` は書くだけ(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 両端が同じ level_key の dft 極小でない | BLOCKED(`endpoints_not_on_one_pes`) |
| 2 | 両端が同じ basin で縮退でない | SAME_BASIN(`same_basin`) |
| 3 | ΔE_rxn > 40 kcal/mol | OUT_OF_WINDOW(`out_of_window`) |
| 4 | 予算(1 反応 6 h)切れ | UNRESOLVED(`walltime`) |
| 5 | 接続を判定済み | elementary / degenerate / reassigned で完了(`connection:<label>`)。両側が同じ basin で 2 回未満なら CONNECT(振幅 × 2、上限 0.4 Å、`connection_retry`)。\|ν\| < 50 の TS なら行 8 へ、それ以外は UNRESOLVED(`connection_failed`) |
| 6 | 未接続の SaddleClaim がある | CONNECT(`ts_validated`) |
| 7 | saddle が収束し未検証 | VALIDATE_TS(`saddle_converged`) |
| 8 | TS 検証で崩壊(最低モード ≥ −10)か、\|ν\| < 50 の TS の QRC が失敗 | VALIDATE_INTERMEDIATE(`saddle_collapsed`) |
| 9 | 中間体が両端と別の basin | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(深さ 2 まで。子の torsional は子の両端の結合グラフで決め直す) |
| 10 | SCREEN が有効で未実行 | SCREEN(`screen`) |
| 11 | 最新の DFT プロファイル(SCREEN か string)が barrierless | BARRIERLESS(`screen:barrierless` / `string:barrierless`) |
| 12 | 最新の DFT プロファイルが intermediate で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_intermediate`、最も低い井戸のノード) |
| 13 | 未使用の種があり saddle の試行 < 2(高次 saddle を押した種と、maxiter で止まった saddle の最終フレームも同じ予算) | REFINE_SADDLE(`seed:<source>`) |
| 14 | DFT string が未実行 | FIND_PATH(`no_dft_path`) |
| 15 | 最新の string が single か intermediate で、種が残っておらず、string < 3 チャンク | FIND_PATH(`next_chunk`、その string の経路から) |
| 16 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

- SCREEN: 低レベル TS(discovery / mode-follow)があれば xTB freq で確認し、その構造の DFT SP と両端の DFT 極小の 3 点が single なら種にする。なければ DFT 極小の間の IDPP 11 点(sin(πt) の対称性破り)を初期経路にして `pysis_neb`(xTB の CI-NEB、端は固定、未収束でも経路として使う。同じ入力で CI から TSOpt)で緩和し(隣の画像へ逐次整列する)、内部ノード 9 点の DFT SP と両端の DFT 極小のエネルギーを `barrier_verdict` で分ける。経路の両端は常に DFT 極小なので、低レベルの PES に端点の極小は要らず(H2O2 gauche、H + H2)、barrierless はそのまま上界として case を閉じる。single の種は xTB freq で確認した NEB の TS(`screen_ts`)か、検出した山を放物線補間した構造(`screen_hei`)、intermediate は行 12 で最も低い井戸を検証する。NEB が失敗したら IDPP のまま分ける(note `screen_neb:<kind>`)。IDPP か DFT SP の失敗は `unavailable` で、FIND_PATH へ進む。below_zpe(高い方の端点の DFT freq で接線モードの ½hν と比べる)は読む側がなかったので削除した。ZPE で消える障壁は熱化学の量として扱う。
- VALIDATE_INTERMEDIATE: 崩壊した saddle か、最新のプロファイルの最も低い井戸を `relax_to_minimum` にかける。きっかけ(崩壊した saddle の ts_check・last_saddle、QRC が失敗した柔らかい TS の claim・connection)はこの action が消費するので、行 8 は次の saddle でだけ発火する。両端と別の basin なら MULTI_STEP。井戸が端点に落ちたら、内部の最大 − 高い方の端点 < 解像度なら barrierless、そうでなければ最も高い山を種にする(REFINE_SADDLE は判定を解かないので、その山が失敗しても同じ井戸と種は繰り返さない)。新しい DFT プロファイル(SCREEN か string)を記録すると、前の saddle の判定(last_saddle・ts_check・intermediate)は解かれ、そのプロファイルの井戸を改めて検証する。解像度は `gates: resolution_kcal`(既定 1.0)。
- REFINE_SADDLE: 反応方向は 1 つの規則で決める(`Ctx.direction`)。低レベル TS(`screen_ts`・`discovery_ts`)の種は、SCREEN で確認した xTB freq の唯一の虚モード、高次 saddle を押した種(`higher_order_retry`)はその反応モード。それ以外は経路の接線(経路外の種は、写像付きで整列した端点の差)のうち、反応中心(`topology.bond_changes` で結合が変わる原子とその隣接原子)の成分だけを残して正規化したもの。結合が変わらなければ(ねじれ・反転)宣言座標の勾配、なければ両端で最も変わる二面角(結合の鎖 i–j–k–l)の勾配、二面角もなければ全原子の接線。傍観者の自由度が接線を支配しないためで、acac の PT では全原子の差と PT モードの重なり 0.22 が、反応中心に限ると 0.96 になる。初期 Hessian は、種が近傍(0.5 Å 以内)の freq を持てばそれ(検証済みの TS freq、`saddle_hessian:ts_freq`)、なければ種で xTB freq をとり、ν < −10 の負モードのどれかと方向の重なりが 0.3 以上ならそれを、そうでなければ種で DFT freq を使う(負モードの本数は問わない)。saddle には方向を `mode` として渡し、adapter が方向だけを負にした Hessian で探索する(§7)。log には `saddle_hessian:<xtb|dft>:overlap:<値>` を書く。saddle が maxiter で止まったら(更新 Hessian の小さい固有値で停滞する)drv.hess で継続せず、最終フレームを種(`saddle_restart`、方向は同じ)として先頭に積み、新しい Hessian で 1 回だけやり直す。VALIDATE_TS が higher_order なら、方向と最も重なる負モードを反応モードとし、それ以外で最も負のモードに沿って片側に 1 回だけ `modes.amplitude` の振幅で押した構造を種にし、検証済みの TS freq をその Hessian にする(DFT Hessian を計算し直さない)。
- FIND_PATH: 1 回の呼び出しで ZTS を 1 チャンク(maxiter 20)だけ走らせ、収束を問わずすぐ `barrier_verdict` で分ける(連続な DFT 経路の最大値は saddle の上界になり、山は種にすぎない。種の質は鞍点探索の成否で確かめる)。種を加えるのは single のときだけで、検出した山を放物線補間する。初期経路は最新の DFT プロファイルの経路(SCREEN の NEB か前のチャンク)を再標本化したもの、なければ新しい IDPP。種が尽きて最新の string が single か intermediate なら、行 15 がその経路から次のチャンクを続ける(case あたり最大 3 チャンク)。NWChem の `freezeN` は凍結した bead の座標を保たない(画像の間に剛体の跳びがあると凍結端点も回され、ゆがむ)ので、初期経路と各チャンクの経路では端点を DFT 極小にし、画像を隣へ逐次整列する。NWChem の converged と gmax は使わない。
- CONNECT(QRC): 振幅はエネルギー目標 max(3 × drop, 3e-4 Eh) と TS の ν・質量の曲率(`modes.amplitude`)から決めて 0.05〜0.4 Å に収め、再試行(両側が同じ basin のときだけ)では 2 倍にしてから 0.4 Å で切る。± 変位を TS の freq Hessian(`inhess 2`、trust 0.3)で opt し、`Registry.find`(ねじれ反応は周期を考慮した最近傍)、未知なら `relax_to_minimum` で新しい basin にしてから `connection` で判定する。

## 7. 化学プロトコルの既定値(設計書 §8 を実装に合わせて更新)

| 段 | 既定値と規則 |
|---|---|
| 共通 | 虚振動 −50 / −10 cm⁻¹。DFT は PBE0-D3BJ/def2-SVPD、grid fine。**NWChem の DFT デッキはすべて `convergence energy` を書く**(既定 1e-7。opt・freq・SP を同じ数値層にそろえる)。SCF の反復上限は DFT が `iterations 100`、WFT の scf が `maxiter 100`(NWChem 7.2.3 の既定は DFT 50、SCF 20 で、開殻・陰イオン・TS 近傍で打ち切られる)。SCF 未収束の救済は §8 の ladder の 1 つだけで、smear・fon は使わない(エントロピー項で PES が変わる)。geometry は `units angstrom nocenter noautosym`、autoz が失敗したら Cartesian で再投入。予算は 1 反応 6 h、ジョブの timeout 4 h |
| structures | 読み込み時の検査: 化学種は xyz と SMILES のちょうど一方、組成の個数 ≥ 1、化学種と組成の id は合わせて一意、反応 id も一意、宣言座標は種類ごとに 2 / 3 / 4 個の相異なる添字 ≥ 0、**宣言反応の端点は xyz**(SMILES では原子の対応と配座を決められない)。xyz はそのまま、SMILES は RDKit ETKDG で 1 配座(未指定の立体は任意に選ばれる。同位体を指定した SMILES は `INPUT_INVALID`(xyz の `D` も未対応の元素として拒否)。`'.'` を含むものは不可)。元素・電荷・多重度はここで 1 回だけ `check_electronic_state` で検査する(下流では検査しない): 元素表にない元素は `unsupported_element:<元素>`、d ブロック元素を含む化学種で多重度が未宣言なら `declare_multiplicity`、電子数と多重度のパリティが合わなければ `INPUT_INVALID`。未宣言の多重度は、SMILES なら RDKit の不対電子数 + 1、xyz なら 1。宣言反応の両端は composition id(組成式・電荷・多重度)と原子順序が一致し、座標の添字が原子数の範囲内でなければ `INPUT_INVALID`。組成(conformers で適用)の電荷は成分の和、多重度は宣言値で、未宣言ならスピン結合で許される値が 1 つのときだけその値(複数なら `declare_multiplicity: candidates (…)`)。宣言値が結合で許される集合になければ `INPUT_INVALID`、開殻の成分から多重度 1 を作る宣言は `open_shell_singlet_unsupported`。多重度 1 は閉殻 RKS |
| conformers | 単量体: `crest --gfn2 --quick -T <n> --ewin 6 --chrg --uhf`。重原子 3 個以下で回転可能結合 0 本なら省く。トポロジー停止は `crest_topology` の species として残し、停止構造から同じ設定で 1 回だけ再実行(2 回目も停止したら失敗とし、停止構造と入力を出す。単量体の最低構造には未緩和の入力より停止構造を使う)。**CREST は attempt ディレクトリで `--scratch` なしに実行する**(3.0.2 は scratch との往復で取り込んだ stdout を上書きする)。錯体: 受容原子の lone-pair 円錐(結合から 110〜120°、方位 0/120/240°、±30° 傾け)に極性 H を置く 6 seed。**結合が複数ある受容原子では、すべての結合から円錐角以上離れた方向(sp3 では lone-pair 軸)を使う**。donor / 受容原子がなければ剛体のランダム配置(vdW 半径の和 + 0.5 Å)。衝突は H を含む分子間の対すべてで ≥ 1.2 Å、重原子どうし ≥ 2.2 Å。CREST に渡すのは seed00(生成順の先頭で、最良の seed ではない)だけで、残りの seed は CREST が失敗したときの出力にだけ使う。argv は `--nci --quick … --noopt`、酸塩基の組成は labile H と受容原子に `--notopo`。**`--noopt` を付けるのは、CREST 3.0.2 の初期最適化後のトポロジー検査が `--notopo` を参照せず(setuptest.f90)、酸塩基や陰イオンの錯体がそこで止まるため**。開殻の組成は `--noopt` で trial MTD が失敗することがあり、そのときは placement の seed が出る。`--nci --quick` の順なので実効設定は quick(MTD の長さ × 0.5 で 2.5〜3 ps × 6 本。壁は NCI のまま。`--ewin 6` が quick の 5 を上書きする)。quick にするのは、同じ最安状態が 1.5〜2.8 倍速く得られるからである。選抜は状態ラベルごとに GFN2 エネルギーの低い順の 6 構造で、エネルギーの窓はない(エネルギーのない候補は、そのラベルに数値の候補がないときだけ残す)。スレッド数(`-T` と OMP)は site の `engines.crest.execution.threads` だけから決まる |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル、species id の順に直列(§5)。状態ラベルが変わった seed は members に記録し、DFT には流さない |
| explore | 出発点は組成 × 状態ラベルごとの screen 最低 2 極小で、その screen の最適化構造(`opt_calc` の最終構造)から始める。ReaDuct は source との一致を、source の原子順で比べた結合の組(添字付き、1.15〜1.45 Σr_cov の中間帯は結合)で判定する。配座・立体・E/Z だけ違う IRC の端は same_as_source として捨て、状態ラベルが同じでも原子ごとの結合が違う端(縮退転位。NH3·HF の二重 H 交換、恒等 SN2 など)は生成物とし、原子順を保って残す。trial は元素に依らない結合変化の列挙器 1 つで作る。形成の候補対 (a, b) は、結合しておらず r ≤ Σr_vdW で、グラフ上 3 結合以上離れた対(断片の内外を問わない)。テンプレートは 4 つ: T1 移動・置換(a–b を形成し b–c を 1 本切る。H 移動、1,2- / 1,3- 移動、SN2、H 引き抜き。2 結合しか離れていない対は c が両者をつなぐ 1,2 移動だけ。c が複数なら角 a–b–c が直線に近い順)、T2 リレー(H の T1 を 2 つ連鎖。2 つ目の H が 1 つ目の受容原子から供与原子へ移る交換(NH3·HF の二重 H 交換)は、1 つ目の移動が両原子を近づけるので 2 つ目の接触を問わない)、T3 形成だけ(付加・会合・環化)、T4 切断(多重度 > 1 か電荷 ≠ 0 の source だけ。閉殻・中性の均一開裂は RKS で記述できない)。**原子価の規則**: 結合が増える原子の変化後の結合数が、元素表の最大配位数(H 1、B/C/N 4、O 3、ハロゲン 1、第 3 周期以降のほかの p ブロック 6、金属 9、He/Ne 0。超原子価ヨウ素は対象外)を超える drive は作らない。T1 → T2 → T3 → T4 を巡回して出発点あたり 10 件まで。各テンプレート内は r_ab/Σr_cov の小さい順(リレーは遠い方の接触、切断は最も伸びた結合から)。1,2 移動は空間を隔てた接触の後に並べる(1,2 移動の距離は結合角で決まり、出発構造の配置を表さない。マロンアルデヒドでは骨格の 1,2 移動(約 1.6)が O–H···O のプロトン移動(1.7)より前に来て、上限 10 件から外れていた)。形成対と切断対の WL 原子クラス(`topology.wl_classes`)と形成対の距離(0.1 Å に丸める)が同じ drive は、同じ source の中で 1 件にする。直線分子は 10° 曲げ 0.05 Å 乱してから。NT2 → AFIR(γ 125 → 300 kJ/mol)。SCC 失敗は 1000 K で再試行し 300 K の SP で評価し直す。窓(ΔE‡ > 150 kJ/mol、ΔE_rxn > 100 kJ/mol は `out_of_window`)は stage が適用する。**ReaDuct の `spin_mode` は常に `restricted_open_shell`**(scine-xtb-wrapper 3.0.2 が受け付けるのは any と restricted_open_shell だけ) |
| minima(dft) | 選択 `window`: 組成 × 状態ラベルごとに 6 kcal/mol 以内の 3 構造 + explore の生成物 + 入力単量体の最低構造。`rerank_sp` は組成 × 状態ラベルごとの上位 8 構造を DFT SP で並べ直した後、同じ窓(window_kcal)と per_state で切り直す(always の構造は残す)。`all` は入力 species 全部。開始構造は、screen 極小があればその最適化構造。**opt は初期 Hessian なしなら `trust 0.1`**(対角推定の Hessian では mode-follow の変位から TS より上へ戻ることがある)、**初期 Hessian を渡すなら `trust 0.3`**(NWChem の既定)、maxiter 100。2 断片以上は xTB の Hessian を初期 Hessian にする(trust 0.3。W7 の A/B で点数 −33%) |
| 仮説(`hypotheses.select`) | 優先順: 宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補。同じ状態の配座対は仮説にしない(Curtin–Hammett。ねじれは宣言したときだけ)。宣言以外は別 basin(explore・mode-follow の発見で両端が同じ DFT basin に入った場合は、発見の source と生成物を端点にし、basin 構造で端点の写像が恒等でなければ(`mapped_equivalent`)縮退反応として宣言反応と同じ経路で評価する。恒等なら仮説にしない。同じ basin の宣言反応があれば、そちらが低レベル TS を借りる)、ΔE_rxn ≤ 40 kcal/mol、変化 ≥ 0.2 Å、組成あたり 6 件まで。端点は各 basin の members から恒等写像 RMSD が最小の組(`pick_endpoints`)。縮退・結合変化・ねじれは、入力座標ではなく、basin の最適化構造を各端点の原子順と掌性に並べ替えた構造(`identity.basin_coords`。driver の経路の端点と同じ)で判定する。両端が同じ basin で 2 つの構造がラベル付きで異なれば(`mapped_equivalent`)縮退反応で、宣言された鏡像化(R→S、H2O2 の gauche 鏡像化)もここに入る。explore の negative は仮説を棄却せず、被覆率にだけ使う。**宣言反応の端点が screen の basin に崩壊した seed なら、その basin の DFT 極小を使う**(`_Pool.endpoint_basin`、W7 の修正 F1) |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`、常に autoz で `moddir 1`。adapter は初期 Hessian を `vibrations.shape_hessian` で整える: 質量加重しない P·H·P(並進・回転を射影)を対角化し、反応方向と \|cos\| 最大の固有ベクトルだけを −max(\|λ\|, 1e-3)、ほかを max(\|λ\|, 1e-3) Eh/bohr² にする(合同変換で慣性は保たれるので、`inhess 2` の内部座標でも負は 1 本で、`moddir 1` がそのモードになる)。方向の sha は JobStore のキーに入る。autoz が失敗したら同じデッキを Cartesian(`noautoz`)で再投入し、TIMEOUT の継続は drv.hess を引き継いで `moddir 1` のまま。maxiter は継続せず、最終フレームを持つ `Failure` を返す(§6 の REFINE_SADDLE が新しい Hessian で 1 回やり直す)。初期 Hessian は種から原子あたり 0.5 Å 以内(`HESSIAN_NEAR_A`)の freq なら受け付ける。最小化(QRC の両側、mode-follow)には生の Hessian を渡す(NWChem が負の固有値を \|e\| として扱う)。NWChem 7.2.3 の saddle は \|固有値\| < 1e-4 で勾配が gmax を超えるモードに最急降下の歩を取る(opt_drv.F の smalleig)。床 1e-3 で最初の数歩はこれを避けるが、更新でできた小さい固有値は防げないので、柔らかい反応モード(\|ν\| が数十 cm⁻¹ の TS)は収束しにくい。ZTS は nbeads 9、stepsize 0.05、interpol 3、impose、freeze1 / freezeN、`tol 1e-5`、maxiter 20 のチャンクを FIND_PATH 1 回に 1 つ(case あたり最大 3 チャンク)、初期経路を必ず渡す。impose は xyz_path を渡すと NWChem が使わない(string.F の newchain 分岐の中だけ)。`tol` は NWChem のソースでは変位の閾値で、文書の「最大勾配」とは違う(hfauto は NWChem の収束判定を使わない)。`pysis_neb` は xTB 専用で、hfauto の IDPP を `initial.trj` として Cartesian の CI-NEB(`cos: neb, climb`、端は COS の既定で固定、`opt: lbfgs, max_cycles 100`)で緩和し、同じ入力で CI から rsprfo(xTB Hessian)の TSOpt を 1 回だけ行う。例外で終わっても NEB の画像(`final_geometries.trj`、エラーで止まったら最後のサイクルの `current_geometries.trj`)があれば TS なしの経路として返す。座標系の分岐(DLC・tric・直線の端点のずらし)はない。分割した子反応(§6 の行 9)はそれぞれ新しい 6 h の予算で実行する(深さ 2 までなので、1 反応から最悪 7 case)。`qrc_bounds_A` の上限は `HESSIAN_NEAR_A`(0.5 Å)以下にする(振幅は最大の原子変位で、QRC の opt は TS の Hessian を使うので、超えると `hessian_geometry_mismatch` になる) |
| sp | `methods:` の各 LOT で、順位に使う点だけを計算する(1 つの規則): 順位を付けられる outcome(elementary / degenerate / reassigned / barrierless_at_resolution)の反応の TS と、その反応物・生成物と同じ状態の DFT 極小すべて(δG_eff の最小 G)、組成の単量体の状態の DFT 極小(分離基準)。multi_step・same_basin・blocked・out_of_window・unresolved・outcome なしの反応の点は計算しない。各点は freq の最終構造で計算し、parents に subject(minimum_id か TS の SaddleClaim.freq_calc。構造が同じなら複数)を書く。thermo と手法パネルは、この parents だけで SP と停留点を対応付ける(構造の指紋は使わない)。CCSD(T) は閉殻を RHF の ccsd モジュール、開殻を ROHF 参照の TCE(`scf; rohf; nopen <m−1>; maxiter 100; end`、`tce; 2eorb; 2emet 13; ccsd(t); freeze atomic; end`、`task tce energy`)で計算する(`2eorb 2emet 13` がないと CH3O• の v2 積分の GA 確保が 1,200 MB/rank で失敗する)。T1 診断などの多参照性のゲートはない。CCSD のブロックには `maxiter 50`(既定 20)を書く。`freeze atomic` は ECP 原子の軌道を凍結しない(I の 4s4p4d は相関に入る。golden G27) |
| thermo | GoodVibes 4.3.0 の API に自前の振動数を渡す。QH=False(H は RRHO、S は Grimme の qRRHO)、qs grimme、cutoff 100 cm⁻¹、symm。スケール因子は `vib_scale` の 1 つ(既定 1.0)で、振動数と ZPE の両方に使う。1.0 は未スケールの調和振動数で、PBE0/def2-SVPD の文献の ZPE 因子 0.985(Kesharwani 2015)にしても ΔG‡ は 0.06 以下、ΔG_assoc は 0.09 kcal/mol 以下しか変わらない。負モードは固定の規則(`thermo_frequencies`)で扱う: 極小は全モードを \|ν\| にし、TS は最低モード(反応座標)を除いて残りを \|ν\| にする(GoodVibes の反転は使わない)。経路縮重度 L = σ_R m‡ / (σ‡ m_R)(Fernández-Ramos 2007)のうち、σ の比は種の σ を通して、m は chiral を通して入る: キラルな極小(`MinimumRecord.chiral`)と TS(構造から `identity.is_chiral`)の G に −RT ln 2(298.15 K で −0.41 kcal/mol)を加える。鏡像対は 1 つの basin なので、アンサンブルでも 1 回だけ m = 2 で数える。トンネル補正はない。1 atm で計算して 1 bar / 1 M に換算。温度と標準状態は thermo の `temperatures_K`(既定 [298.15])と `standard_states`(既定 [1atm])だけから決まる。会合量の単量体は、system の組成の各成分と同じ状態(組成と状態ラベル)の DFT 極小のアンサンブルで、異性体を混ぜない。錯体と単量体の LOT 判定は電荷と多重度を除く(`same_pes(state=False)`)。G のない単量体があれば ΔG_assoc は None。反応を持たない錯体では ΔG_assoc は出ない。δG_eff の感度の幅(qs × cutoff 50/100/150)。ΔG_assoc には幅を出さないが、qRRHO の扱いによる幅は 298 K で約 1.2 kcal/mol(TMA·(HF)₂)。失敗は G = None(`thermo_unavailable`)、LOT の不一致は `mixed_level_of_theory`。composite G = E_SP + (G_GV − E_GV) は `energy_method` を指定したときだけで、SP は parents がその subject を指すもの(view で後のものが勝つ)。ない subject は `energy_layer_missing`(G = None)で、その subject を含む ΔE も None。エネルギー層の SP にも `spin_ok` をかけ、不合格は `spin_contaminated` に合流する。blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated)を作るのは thermo stage だけである。`energy_level` は参加する停留点に共通するエネルギー層の method/basis(層が欠けるか混ざれば空)。エネルギー層で ΔE0‡ ≤ 0 になる TS(PBE0 の PES 上の停留点なので起こりうる)は `submerged_barrier` の規則で扱い、ゲートは足さない。JobStore のキーは freq の job_key、Hessian の sha、渡す振動数、温度、settings |
| report | rankable な反応(障壁なしを含む)を δG_eff で並べ、感度の幅が重なれば同順位。(T, 標準状態) は ReportConfig の指定、なければ view の最初の ReactionThermo のもの。ranking.csv の列は rank、reaction_id、outcome、tier、rankable、T_K、standard_state、energy_level(並べたエネルギー層)、dG_eff_kcal、band_low/high_kcal(δG_eff の幅)、dG_act_kcal、dG_rxn_kcal、dG_act_vs_separated_kcal、torsional(結合変化のないねじれの段。mode-follow 由来のねじれを化学反応と見分ける。宣言したねじれ反応も順位に残す)、blockers、notes。エネルギーが 2 つ以上の LOT にある(手法パネルがある)ときだけ `dE_act_panel_min_kcal`・`dE_act_panel_max_kcal` を足す(手法の不確かさは列で示し、順位は止めない。method_panel.csv は各 LOT の値と min/max)。被覆率(機構別の試行・生成物・陰性理由・FailureKind)。method_panel を追記すると、panel_thermo(`energy_method` はパネルの参照手法)の値が同じ id の thermo の値を上書きし(後から追記した層が勝つ)、panel_report はその層で並べる |

- 対応範囲(入口で宣言する範囲): 元素は `hfauto/chemistry/elements.py` の表 1 つ(Z=1〜57、72〜86。Ce〜Lu と Z>86 は範囲外)。これは計算できる範囲で、遷移金属は検証していない。Z>36 の元素は def2 系の基底だけで扱い、NWChem のデッキ(DFT と WFT)に `<元素> library def2-ecp` を元素ごとに自動で書く(`*` は使わない。def2 以外の基底では `INPUT_INVALID`(`no_ecp_for_basis`))。単原子の化学種も同じ経路(opt → freq → thermo)を通り、freq の Evidence は n_external 3・振動なし、熱化学は並進と電子(g = 2S+1。L>0 の縮重とスピン軌道は含めない)だけになる。explore の出発点にはならない。スピン状態は利用者が宣言し、閉殻と高スピン(UKS)を対象とする。開殻一重項(ビラジカル、ラジカル対)と MECP・スピン交差は範囲外で、BS-UKS は使わない。ラジカル再結合の生成物は単量体として宣言する。原子の SMILES は高スピン仮定なので、実際の基底状態と違うことがある([C] は RDKit で五重項、実際は ³P)。
- 順位の量: δG_eff = max(G_TS, G_R, G_P) − G_R(`thermo.effective_barrier`)。G_R・G_P は、反応物・生成物と同じ状態(組成と状態ラベル)で freq と energy の LOT が同じ DFT 極小のうち最小の G(速い配座平衡の下では反応物は状態として 1 つ。窓で DFT に流す配座は数件なので、アンサンブルの G は採らない。両端が同じ状態の配座変化(HONO の trans → cis など)は、どちらの向きも状態の最低配座から測る)。TS が反応物か生成物より下にあれば律速ではない(微視的可逆性)。順方向か逆方向の ΔE0‡ = ΔE‡ + ΔZPE‡(ΔE は energy 層)が 0 以下なら、鞍点の TST は意味を持たないので TS を除いて max(ΔG_rxn, 0) とし、注記 `submerged_barrier` を付ける(blocker ではない)。barrierless_at_resolution も max(ΔG_rxn, 0) で、解像度未満の障壁の下限値であることは outcome 列で分かる。ΔG‡(`dG_act_kcal`)はその反応の反応物極小から見た TST の値で、表示だけに使う。分離反応物の G は基準に混ぜない(`dG_act_vs_separated_kcal` は別の列。反応はすべて Δn = 0 なので δG_eff は標準状態に依らない)。トンネル補正と σ の自前判定は入れない。
- 結果の読み方: 停留点レベルを PBE0-D3BJ/def2-SVPD にしたのは、FHF⁻ の De が diffuse 関数のない SVP では 21.7 kcal/mol 過大になるが SVPD では実験値の誤差内に入るからである。既定の順位はこのレベルの δG_eff で、PBE0 の障壁の系統誤差は反応クラスで違う(平均符号付き誤差: H 移動 −4.2、重原子移動 −6.6、求核置換 −1.9 kcal/mol。Zhao & Truhlar, J. Phys. Chem. A 109, 2012 (2005))。CCSD(T)/def2-TZVPD(SVPD 構造での SP)との差は HCN −1.2、HONO +2.0、SN2 −3.0 kcal/mol(validation.md の S-A)。順位は δG_eff の序数として読み、反応クラスが混ざるときは、差が数 kcal/mol 未満の順序は手法誤差で入れ替わり得る。SCREEN・FIND_PATH の barrierless は停留点レベル(PBE0)の PES 上の判定で、エネルギー層では確かめない。ΔG_assoc は BSSE を補正しない SVPD の値で、結合を過大に見積もることがある。amine·(HF)n の結論(プロトンの位置、同一 basin)は PBE0 の PES 上のものである。
- エネルギー層(composite G)は 1 つの手順で使う: 既存の run に `method_panel` を `--run-dir` で追記する(panel_sp → panel_thermo(`energy_method` = 参照手法)→ panel_report)。thermo の結果の id は stage id を含まないので、後から追記した層が元の thermo の値を上書きし(view で後のものが勝つ)、ranking.csv の `energy_level` 列に並べた層が出る。元の値を残したいときだけ run-dir を複製する(`cp -a`。JobStore と相対パスの FileRef ごと移る)。SP の対象は順位に使う点(sp の行)で単量体を含むので、ΔG_assoc も同じ層で出る。ωB97X-D3/def2-TZVPD の層は閉殻 4 系中 3 系で CCSD(T) に近い(平均絶対誤差 2.1 → 1.3 kcal/mol)が、17 原子(TMA·(HF)₂)の SP が 1 点 29 分(PBE0/SVPD の SP 35 秒、freq 約 15 分)かかるので、既定の pipeline には入れていない(validation.md の S-A)。

## 8. JobStore と attempt ladder

- JobStore(`hfauto/execution/jobstore.py`)は run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合し、`adapter.result_type` で復元する(検証に失敗したら miss)。
- 終端的な失敗(input_invalid、method_mismatch、executable_missing、incomplete_output)も記録し、`--retry-failed` で指定した種類だけを再実行する。ジョブ単位の排他は `.lock` を O_EXCL で作る。
- ladder(`hfauto/execution/jobs.py` の `LADDER`。これ以外の自動再試行はしない): timeout と geometry_maxiter は `continuation` で 2 回まで継続(NWChem は最新構造・movecs・drv.hess、xTB は最終構造。saddle の geometry_maxiter だけは継続しない、§6)、input_invalid(autoz)は Cartesian で 1 回、scf_not_converged は前回の vectors から `cgmin`(二次収束)で 1 回(P3c の H3 bead で damping・rabuck は収束せず、cgmin だけが収束した)。cgmin は ⟨S²⟩ を出さないので、⟨S²⟩ のない開殻 DFT の Evidence は incomplete_output(`s2_not_reported`)にする。初期 Hessian ありの opt の継続は、drv.hess を INHESS 0 の restart で引き継ぎ、trust 0.1 で続ける。attempt の timeout は min(設定値, 残り時間) で、残り 60 s 未満なら `BUDGET_EXHAUSTED`。
- 失敗の閉じ込め: `adapter.parse` の例外は `Failure(incomplete_output, 'parse:<型>: <1行目>')` になる(終端的で JobStore に残り、ladder には入らない。パーサを直したら `--retry-failed incomplete_output`)。反応ケースの例外はそのケースを `unresolved`(`error:<型>`)にし、log.jsonl に error の行を 1 行書いて、ほかのケースを続ける。消費する型の入力がない stage は 0 artifact で done になり、上流の失敗を最大 3 件 warning に出す。後続の stage(report を含む)は同じコマンドで実行される。`HFAUTO_STRICT=1`(tests/conftest.py が設定)は前 2 つを再送出にする(テスト専用)。
- 並行実行: `JobRunner` のセマフォで ranks × threads の合計を site の cores 以下に保つ。極小(xTB と DFT)と反応ケースは直列で、CREST と GoodVibes だけを `thread_map` で並列にする。`SiteLock` は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。
- run ディレクトリ: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/diagnostics.json`(conformers と minima が書く。読むコードはない)、`<run>/run_state.json`(実行順の stage 状態と JobStore の統計)、`<run>/resolved_config.yaml`。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に union したもの。done で入力と設定の sha が変わらない stage は飛ばし(再開)、`--from X` は X 以降を stale にする。resolved_config.yaml の code_version は git の sha(追跡中のファイルに差分があれば -dirty、改行コードだけの差は数えない、git がなければ unknown)で、pipeline_id ごとに最後に実行したときの値に上書きされる。done で飛ばした stage を計算した版は残らないので、版をまたいだ結果が要るときは新しい run dir か `--from` で取り直す。

## 9. 設定の 4 分割

層どうしはマージせず、優先順位もない。絶対パスは site ファイルにだけ書く。

| 層 | 置き場所 | 中身 |
|---|---|---|
| site | `configs/sites/wsl_local.yaml` | scratch_root、cores、エンジン別の `EngineSite`(版数の pin、実行ファイル、`ExecutionSpec`、worker の python、scratch_dir) |
| method | `configs/methods/` | `MethodSpec`(JobStore のキーに入る): `gfn2`、`pbe0-d3bj_def2-svpd`(停留点)、`pbe0-d3bj_def2-tzvpd` と `wb97x-d3_def2-tzvpd`(手法パネル)、`ccsd-t_def2-tzvpd`(小さい分子で opt-in。開殻は ROHF-CCSD(T))。停留点の method は minima(dft)と reaction-paths に書き、食い違えば読み込み時に拒否する |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | species(xyz / SMILES、電荷、多重度、役割)、compositions(多重度を宣言できる)、宣言反応: hcn、hono、nh3_inversion、formaldehyde、tma_hf2、amine_hf_panel、water_same_basin と、レビュー §6 の検証セット(sn2_cl ほか 17 系。対応は docs/validation.md の r6) |
| pipeline | `configs/pipelines/` | stage の並びと設定、`gates:`(`Policy` の上書き): discover、known_endpoints、method_panel。温度と標準状態は thermo stage の設定。method は pipeline の隣の `methods/` にあればそれを、なければ作業ディレクトリの `configs/methods/` を使う(CLI の名前解決と同じ) |

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
| `conformers` | `hfauto.stages.conformer_search:ConformersStage` | species | species | engine*, method*, quick, ewin_kcal, seeds_per_composition, keep_per_state |
| `minima` | `hfauto.stages.minima:MinimaStage` | species | calculation, minimum, species, discovery | level*, engine*, method*, select, init_hessian, mode_follow |
| `explore` | `hfauto.stages.explore:ExploreStage` | minimum, species | discovery, species | engine*, method*, sources_per_state, max_trials_per_source, settings, window |
| `reaction-paths` | `hfauto.stages.reaction_paths:ReactionPathsStage` | minimum, species | reaction, calculation, minimum, species | method*, engines*, screen, policy, reaction_ids |
| `sp` | `hfauto.stages.single_point:SinglePointStage` | minimum, reaction | calculation | engine*, methods* |
| `thermo` | `hfauto.stages.thermochemistry:ThermoStage` | minimum, calculation | species_thermo, reaction_thermo | engine*, settings, energy_method, temperatures_K, standard_states |
| `report` | `hfauto.stages.report:ReportStage` | — | report | T_K, standard_state |

## 11. 設計書から変えた点

| 項目 | 設計書 | 現状(理由) |
|---|---|---|
| 削除したもの | core の units モジュール、`Level.surface_key`、`vibrations.curvature_along` と `frame_residual`、NWChem の最終振動ブロックの読み取り、`SpeciesRecord.formula`、`MinimumRecord.n_fragments`、`ReactionRecord.parent_id` | 使われていなかった。経路の PES 比較は `same_pes(numerics=False)` で行う。子反応の親は id(`<parent>_split<n>`)で分かる |
| MinimumDriver | `relax_to_minimum` の引数と `Registry(minima)` | `load_xyz=` を追加し、`Registry` は (MinimumRecord, Geometry) の組で作る。`MinimumOutcome.notes` を `MinimumRecord.notes` に写す |
| `CaseState` | 16 行の判断に必要な状態 | `minima`(両端の MinimumRecord)を追加した。行 1〜3 は tier・level_key・basin・エネルギーを要する |
| 電子状態 | `resolve_electronic_state` | `check_electronic_state(symbols, charge, multiplicity, declared=)` は検査だけを行う。組成の多重度は `composition_multiplicity` が宣言とスピン結合の集合から決める |
| minima stage | species と discovery を消費 | species だけを消費する(explore の生成物は species として渡る) |
| DFT 極小の freq の省略 | 既存の極小と same なら省く | relax の直後に登録するので、同じ stage の中でも登録済みの basin に落ちたジョブは freq を省く(§5) |
| 計算仕様レビュー(2026-09-26)の改良 | §7.3 の判断表と §8 の既定値 | 低レベルの陰性結果による棄却、barrierless の多解像度確認、同一性の 2 段判定と再最適化、スケール因子の表引き、conformers のエネルギー窓を削除し、読み込み時の入力検査を加えた(§6、§7。根拠は `docs/reviews/2026-09-26_calculation_spec_review.md` §4) |
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
