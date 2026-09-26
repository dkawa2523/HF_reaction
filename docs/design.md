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
| reporting | `hfauto/reporting/`(566) | 順位・被覆率・手法パネルの表(`hfauto/reporting/summary.py`)と静的 HTML(`hfauto/reporting/html.py`) |
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

`Evidence`(`hfauto/core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すること自体が次を保証する: (1) 正常終了と SCF 収束、(2) opt / saddle なら構造最適化の収束、(3) freq なら 3N − n_external 本の振動数、(4) 原子順序の保持、(5) 入力フレームの保持(出力に出た初期構造が入力と 1e-4 Å 以内)、(6) 出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること(`level_mismatches`)。満たさなければアダプタは `Failure` を返す。振動数と正準形式の Hessian(`.npy`、Eh/bohr²、入力フレーム)は、どのエンジンでも `hfauto/chemistry/vibrations.py` の補空間射影(並進・回転空間を SVD で決め、主同位体の質量を使う)で求める。経路は `PathProfile`、失敗は `Failure(kind: FailureKind, reason)` で表す。FailureKind は executable_missing、input_invalid、timeout、stagnated、nonzero_exit、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、budget_exhausted、gate_rejected の 11 種である。

ゲートは `hfauto/chemistry/gates.py` にだけ置く。閾値は `Policy` にまとめ、pipeline YAML の `gates:` だけが上書きできる。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=True, state=True)` | program・version・method・basis・dispersion・solvation・charge・multiplicity・電子温度が一致。`numerics=True` なら grid と scf_tol も。停留点・freq・QRC・熱化学は True、DFT string(種にだけ使う)は False。`state=False` は電荷と多重度を比べず、会合量の錯体と単量体の LOT 判定だけで使う |
| `spin_ok` | `s2` が None か、\|<S²> − S(S+1)\| ≤ 0.1。不合格は注記 `spin_contaminated`(rankable の blocker) |
| `imaginary_tier` | ν < −50 は saddle、−50 ≤ ν < −10 は soft、−10 ≤ ν < 0 は noise |
| `is_minimum(freq, opt=)` | task が freq / opt、`freq.start.fingerprint == opt.final.fingerprint`、numerics を含めて同一 PES、本数が 3N − n_external、tier が saddle でない(soft / noise は注記) |
| `is_first_order_saddle(freq, saddle=)` | 指紋と同一 PES、本数、最低モード < −50(ねじれは −20)、ほかに < −50 がない(あれば `higher_order`)、−50〜−10 は注記 `soft_secondary_mode`、E_TS − max(端点) ≥ max(2e-5 Eh, 20 × scf_tol)。モードの重なりは合否に使わない |
| `barrier_verdict` | 基準は高い方の端点。DFT 3 点未満は `unavailable`、DFT か低レベルの内部極大 ≥ 1.0 kcal/mol なら `proceed`、どちらも未満なら `barrierless`(ノード間隔を解像度として記録、接線モードの ½hν 未満なら `below_zpe`) |
| `connection` | 両側で drop = `qrc_drop(level)` = max(1e-5 Eh, 20 × scf_tol)(CONNECT の振幅もこの値から決める): 変位直後に E_TS − drop 以下、軌跡の最大も以下、最終値が先頭 − drop 以下、同一 PES。割付けの集合が期待どおりなら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつ両側の構造が別なら degenerate |
| `thermo_consistent` | GoodVibes の実振動の本数、ZPE(許容 1e-5 Eh)、E(許容 1e-6 Eh)が、渡した振動数(`thermo_frequencies`、`vib_scale` 倍)と freq の E に一致 |
| `rankable` | outcome が elementary / degenerate / reassigned、ΔG‡ がある、blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated、method_sign_disagreement)がない、\|ΔE0‡ − ΔE‡\| ≤ 1 + 4 × max(1, 移動する H の数) kcal/mol。順位を付けられない outcome は `outcome:<o>` と本当の blocker だけを返し(ΔG‡ と ΔE0‡ − ΔE‡ は検査しない)、`dzpe_out_of_tolerance` は ΔE0‡ − ΔE‡ に値があるときだけ付ける |
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
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、source、`BarrierVerdict`(proceed / barrierless / unavailable)、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome`(elementary_step、degenerate_rearrangement、reassigned_step、multi_step、barrierless_at_resolution、same_basin、out_of_window、unresolved_within_budget、blocked_upstream) |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None)、Boltzmann 重み、注記 |
| reaction_thermo | `ReactionThermo` | thermo | ΔE‡、ΔG‡、ΔG_rxn、ΔG_assoc、単量体基準の ΔG‡、感度の幅、blockers |
| report | `ReportRecord` | report | 順位の行と表(ranking.csv、coverage.csv、method_panel.csv、report.html) |

id の規則: species は `species_<species_id>`、minimum は `min_<species>_<level_key[:8]>`、失敗した極小は `min_<species>_<stage_id>`、分割した子反応は `<parent>_split<n>`、report は `report_<stage_id>`。

## 4. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス(`hfauto/backends/engines.py`) | 結果 |
|---|---|---|---|
| qm | `energy`、`optimize(init_hessian)`、`frequencies` | `nwchem` → `NWChemEngine`(DFT と、閉殻の MP2 / CCSD(T) の SP)、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path, refine_ts)` | `nwchem_string` → `NWChemString`(ZTS)、`pysis_gs` → `PysisGrowingString`(xTB ネイティブの GS と TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode_index)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2 / AFIR、worker で実行) | `DiscoveryResult` |
| thermo | `thermo(freq, settings, temperatures_K, saddle)` | `goodvibes` → `GoodVibesEngine`(4.3.0 の API、worker で実行) | `ThermoResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作り、`requirements()` は classmethod とする。`engines.create` は結果が能力の Protocol を満たさなければ `TypeError`、表にないキーは `KeyError` にする。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult` / `odft`、xTB と CREST の `--chrg` / `--uhf`、pysis の calc、ReaDuct の計算器設定)。
- 初期 Hessian は freq `Evidence` として渡す(Level は問わない)。saddle は種と同じ構造(指紋が一致)、optimize は原子順序とフレームが同じで各原子が 0.5 Å 以内の構造の Hessian を受け付け(QRC の両側は TS の Hessian を使う)、外れれば `INPUT_INVALID("hessian_geometry_mismatch")`。xTB は初期 Hessian を使わない。
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
3. saddle 級の虚振動なら最低モードに沿って ±0.1 Å 変位し、両側を opt → freq(最大 2 サイクル)。両側が同じ極小(`identity.same_minimum`)なら置き換え、別々の極小なら `ts_candidate`(minima stage が `mode_follow` の discovery にする)、片側だけならその側を採る。
4. soft なら 1 回だけ押し、残れば `soft_minimum`。noise はそのまま `minimum`。

`Registry` は組成 × level_key ごとの basin の集合で、同一性の基準は `identity.assign` の 1 つだけである(置換不変 RMSD ≤ 0.05 Å かつ \|ΔE\| ≤ 5e-5 Eh で、次点と十分離れていること)。`add` は known の basin か、同じ composition_id・level_key・元素列の basin のうち assign が一意に割り付くものに加え、なければ新しい basin を作る。minima stage はジョブを species id の順に直列で relax し、relax の直後に登録する(mode-follow の両側は親の直後)。このため同じ stage の中でも、登録済みの basin に落ちたジョブは `known` になり、freq は basin あたり 1 本で済む。

振動数は質量加重の Eckart 射影(`hfauto/chemistry/vibrations.py`)で求めるので、NWChem が表示する Projected Frequencies とは低振動モードで数 cm⁻¹ ずれる(G では約 0.02 kcal/mol)。

## 6. ReactionCaseDriver(`hfauto/drivers/reaction_case/`)

`decide(case, state, policy)` は `hfauto/drivers/reaction_case/state.py` の 17 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。`<run>/<stage_id>/cases/<reaction_id>/log.jsonl` は書くだけ(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 両端が同じ level_key の dft 極小でない | BLOCKED(`endpoints_not_on_one_pes`) |
| 2 | 両端が同じ basin で縮退でない | SAME_BASIN(`same_basin`) |
| 3 | ΔE_rxn > 40 kcal/mol | OUT_OF_WINDOW(`out_of_window`) |
| 4 | 予算(1 反応 6 h)切れ | UNRESOLVED(`walltime`) |
| 5 | 接続を判定済み | elementary / degenerate / reassigned で完了(`connection:<label>`)。failed で 2 回未満なら CONNECT(振幅 × 2、上限 0.4 Å、`connection_retry`)、それ以外は UNRESOLVED(`connection_failed`) |
| 6 | SaddleClaim がある | CONNECT(`ts_validated`) |
| 7 | saddle が収束し未検証 | VALIDATE_TS(`saddle_converged`) |
| 8 | TS 検証で崩壊(n_imag = 0) | VALIDATE_INTERMEDIATE(`saddle_collapsed`) |
| 9 | 中間体が両端と別の basin | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(深さ 2 まで) |
| 10 | 高次の saddle で試行 < 2 | REFINE_SADDLE(`higher_order_retry`、2 本目のモードに沿って 0.1 Å 下った種) |
| 11 | SCREEN が有効で未実行 | SCREEN(`screen`) |
| 12 | SCREEN が barrierless | BARRIERLESS(`screen:barrierless`) |
| 13 | 直近の DFT 経路が multi_max で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_multi_max`) |
| 14 | 直近の DFT 経路が monotonic | BARRIERLESS(`dft_path_monotonic`) |
| 15 | 未使用の種があり saddle の試行 < 2 | REFINE_SADDLE(`seed:<source>`) |
| 16 | DFT string が未実行 | FIND_PATH(`no_dft_path`) |
| 17 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

- SCREEN: 低レベル TS(discovery / mode-follow)があれば xTB freq で確認し、その構造と両端の DFT SP の 3 点で判定する。なければ DFT の両端を xTB で再最適化し、崩壊しなければ `pysis_gs`(11 ノード、CI、同じ入力で TSOpt)の全ノード(隣の画像へ逐次整列する)と TS で DFT SP、崩壊したら DFT の IDPP 11 点で判定する。種は xTB の TS(`screen_ts`)か DFT の HEI(`screen_hei`)。低レベルの失敗は `unavailable` で、FIND_PATH へ進む。
- REFINE_SADDLE: 1 回目は種の xTB Hessian(負の固有値が 1 本で接線との重なり ≥ 0.3 のとき)、それ以外は DFT Hessian を渡す。`mode_index` は接線(または宣言座標の勾配)と最も重なる負モード。
- FIND_PATH: 初期経路は低レベル GS の経路(写像付きで整列して再標本化)→ SCREEN の IDPP → 新しい IDPP の順。ZTS のチャンク(maxiter 20)は、bead エネルギーの変化が直近 3 反復で 0.1 × 解像度(0.1 kcal/mol)未満になるまで、最大 3 回続ける。NWChem の `freezeN` は凍結した bead の座標を保たない(画像の間に剛体の跳びがあると凍結端点も回され、ゆがむ)ので、初期経路、チャンクの継ぎ目(次のチャンクの初期経路)と最終経路では端点を DFT 極小にし、画像を隣へ逐次整列する。1 チャンク内の端点のずれ(逐次整列の後は 0.004 Å 以下)は note `string<n>:c<k>:end_drift:<写像 RMSD Å>` で監視する(ゲートにはしない)。NWChem の string は gmax を記録しない(pysis_gs だけ)。NWChem の converged は使わない。経路の形は収束を問わずそのまま記録する(連続な DFT 経路の最大値は saddle の上界になる)。HEI を種に加えるのは `single_max` のときだけ。
- CONNECT(QRC): 振幅はエネルギー目標 max(3 × drop, 3e-4 Eh) と曲率から決めて 0.05〜0.4 Å に収め、再試行では 2 倍にしてから 0.4 Å で切る。± 変位を TS の freq Hessian(`inhess 2`、trust 0.3)で opt し、`Registry.find`(ねじれ反応は周期を考慮した最近傍)、未知なら `relax_to_minimum` で新しい basin にしてから `connection` で判定する。

## 7. 化学プロトコルの既定値(設計書 §8 を実装に合わせて更新)

| 段 | 既定値と規則 |
|---|---|
| 共通 | 虚振動 −50 / −10 cm⁻¹(ねじれ TS は −20)。DFT は PBE0-D3BJ/def2-SVPD、grid fine。**NWChem の DFT デッキはすべて `convergence energy` を書く**(既定 1e-7。opt・freq・SP を同じ数値層にそろえる)。geometry は `units angstrom nocenter noautosym`、autoz が失敗したら Cartesian で再投入。予算は 1 反応 6 h、ジョブの timeout 4 h |
| structures | 読み込み時の検査: 化学種は xyz と SMILES のちょうど一方、組成の個数 ≥ 1、化学種と組成の id は合わせて一意、反応 id も一意、宣言座標は種類ごとに 2 / 3 / 4 個の相異なる添字 ≥ 0、**宣言反応の端点は xyz**(SMILES では原子の対応と配座を決められない)。xyz はそのまま、SMILES は RDKit ETKDG で 1 配座(未指定の立体は任意に選ばれる。同位体を指定した SMILES は `INPUT_INVALID`(xyz の `D` も未知の元素として拒否)。`'.'` を含むものは不可)。`check_electronic_state` で電子数と多重度を検査し(未知の元素も `INPUT_INVALID`)、宣言反応の両端は composition id(組成式・電荷・多重度)と原子順序が一致し、座標の添字が原子数の範囲内でなければ `INPUT_INVALID`。組成の電荷は和、多重度は高スピン。多重度 1 は閉殻 RKS で、開殻一重項(ジラジカル)は扱わず、検出もしない |
| conformers | 単量体: `crest --gfn2 --quick -T <n> --ewin 6 --chrg --uhf`。重原子 3 個以下で回転可能結合 0 本なら省く。トポロジー停止は `crest_topology` の species として残し、`--noreftopo` で 1 回だけ再実行。**CREST は attempt ディレクトリで `--scratch` なしに実行する**(3.0.2 は scratch との往復で取り込んだ stdout を上書きする)。錯体: 受容原子の lone-pair 円錐(結合から 110〜120°、方位 0/120/240°、±30° 傾け)に極性 H を置く 6 seed。**結合が複数ある受容原子では、すべての結合から円錐角以上離れた方向(sp3 では lone-pair 軸)を使う**。donor / 受容原子がなければ剛体のランダム配置(vdW 半径の和 + 0.5 Å)。衝突は H を含む分子間の対すべてで ≥ 1.2 Å、重原子どうし ≥ 2.2 Å。CREST に渡すのは seed00(生成順の先頭で、最良の seed ではない)だけで、残りの seed は CREST が失敗したときの出力にだけ使う。argv は `--nci --quick … --noopt`、酸塩基の組成は labile H と受容原子に `--notopo`。**`--noopt` を付けるのは、CREST 3.0.2 の初期最適化後のトポロジー検査が `--notopo` を参照せず(setuptest.f90)、酸塩基や陰イオンの錯体がそこで止まるため**。開殻の組成は `--noopt` で trial MTD が失敗することがあり、そのときは placement の seed が出る。`--nci --quick` の順なので実効設定は quick(MTD の長さ × 0.5 で 2.5〜3 ps × 6 本。壁は NCI のまま。`--ewin 6` が quick の 5 を上書きする)。quick にするのは、同じ最安状態が 1.5〜2.8 倍速く得られるからである。選抜は状態ラベルごとに GFN2 エネルギーの低い順の 6 構造で、エネルギーの窓はない(エネルギーのない候補は、そのラベルに数値の候補がないときだけ残す)。スレッド数(`-T` と OMP)は site の `engines.crest.execution.threads` だけから決まる |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル、species id の順に直列(§5)。状態ラベルが変わった seed は members に記録し、DFT には流さない |
| explore | 出発点は組成 × 状態ラベルごとの screen 最低 2 極小で、その screen の最適化構造(`opt_calc` の最終構造)から始める。ReaDuct は source との一致を状態ラベルで判定する。状態ラベルは結合グラフ(1-WL、立体と E/Z を区別しない、1.15〜1.45 Σr_cov の中間帯は結合)なので、配座・立体・E/Z だけ違う IRC の端は same_as_source として捨てる。trial は極性 H → C–H の 1,2- / 1,3- 移動 → 重原子の結合組み換え → 会合の順で出発点あたり 10 件まで。直線分子は 10° 曲げ 0.05 Å 乱してから。NT2 → AFIR(γ 125 → 300 kJ/mol)。SCC 失敗は 1000 K で再試行し 300 K の SP で評価し直す。窓(ΔE‡ > 150 kJ/mol、ΔE_rxn > 100 kJ/mol は `out_of_window`)は stage が適用する。**ReaDuct の `spin_mode` は常に `restricted_open_shell`**(scine-xtb-wrapper 3.0.2 が受け付けるのは any と restricted_open_shell だけ) |
| minima(dft) | 選択 `window`: 組成 × 状態ラベルごとに 6 kcal/mol 以内の 3 構造 + explore の生成物 + 入力単量体の最低構造。`rerank_sp` は組成 × 状態ラベルごとの上位 8 構造を DFT SP で並べ直す。`all` は入力 species 全部。開始構造は、screen 極小があればその最適化構造。**opt は初期 Hessian なしなら `trust 0.1`**(対角推定の Hessian では mode-follow の変位から TS より上へ戻ることがある)、**初期 Hessian を渡すなら `trust 0.3`**(NWChem の既定)、maxiter 100。2 断片以上は xTB の Hessian を初期 Hessian にする(trust 0.3。W7 の A/B で点数 −33%) |
| 仮説(`hypotheses.select`) | 優先順: 宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補 → conformer ペア(結合変化なし、二面角差 ≥ 30°)。宣言以外は別 basin、ΔE_rxn ≤ 40 kcal/mol、変化 ≥ 0.2 Å か ≥ 30°、組成あたり 6 件まで。端点は各 basin の members から恒等写像 RMSD が最小の組(`pick_endpoints`)。`mapped_equivalent` なら縮退反応。explore の negative は仮説を棄却せず、被覆率にだけ使う。**宣言反応の端点が screen の basin に崩壊した seed なら、その basin の DFT 極小を使う**(`_Pool.endpoint_basin`、W7 の修正 F1) |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`。追うモードは常に明示する: 負モードが 1 本なら autoz のまま `moddir 1`、2 本以上なら `noautoz` にし、`moddir` の番号は P·H·P(質量加重しない Cartesian で、並進・回転を射影)の負の固有ベクトルの順で決める。autoz / Cartesian の分岐は ν < 0 のモード(noise 以下を含む)の本数で決める。NWChem 7.2.3 の saddle は \|固有値\| < 1e-4 かつ勾配 < gmax_tol のモードをゼロとして扱う(opt_drv.F の smalleig)ので、柔らかい負モードは追えないことがある。TIMEOUT / MAXITER の継続では `moddir` を 1 に抑える。`moddir` と座標系は JobStore のキーに入る。ZTS は nbeads 9、stepsize 0.05、interpol 3、impose、freeze1 / freezeN、`tol 1e-5`、maxiter 20 × 最大 3 チャンク(打ち切りは §6 の FIND_PATH)、初期経路を必ず渡す。impose は xyz_path を渡すと NWChem が使わない(string.F の newchain 分岐の中だけ)。`tol` は NWChem のソースでは変位の閾値で、文書の「最大勾配」とは違う(hfauto は NWChem の収束判定を使わない)。`pysis_gs` は xTB 専用で、直線(全原子が主軸から 0.3 Å 以内)か 2 断片以上なら cart、それ以外は dlc。**厳密に直線の端点は主軸と垂直に ±0.01 Å ずらす**。TSOpt は rsirfo(2 断片以上は tric)で、例外で終わったら GS だけをやり直し、TS なし(`tsopt_error`)の経路を返す。分割した子反応(§6 の行 9)はそれぞれ新しい 6 h の予算で実行する(深さ 2 までなので、1 反応から最悪 7 case)。`qrc_bounds_A` の上限は `HESSIAN_NEAR_A`(0.5 Å)以下にする(振幅は最大の原子変位で、QRC の opt は TS の Hessian を使うので、超えると `hessian_geometry_mismatch` になる) |
| sp | `methods:` の各 LOT で、反応の極小と TS(`targets: all_minima` なら全 dft 極小も)の SP。MP2 と CCSD(T) は閉殻(RHF)専用で、開殻は `INPUT_INVALID`(`wft_closed_shell_only`)。CCSD のブロックには `maxiter 50`(既定 20)を書く |
| thermo | GoodVibes 4.3.0 の API に自前の振動数を渡す。QH=False(H は RRHO、S は Grimme の qRRHO)、qs grimme、cutoff 100 cm⁻¹、symm。スケール因子は `vib_scale` の 1 つ(既定 1.0)で、振動数と ZPE の両方に使う。1.0 は未スケールの調和振動数で、PBE0/def2-SVPD の文献の ZPE 因子 0.985(Kesharwani 2015)にしても ΔG‡ は 0.06 以下、ΔG_assoc は 0.09 kcal/mol 以下しか変わらない。負モードは固定の規則(`thermo_frequencies`)で扱う: 極小は全モードを \|ν\| にし、TS は最低モード(反応座標)を除いて残りを \|ν\| にする(GoodVibes の反転は使わない)。σ は外部回転の対称数だけで、経路縮重度 L(キラリティ、対称反応の因子)は含まない。トンネル補正はない。1 atm で計算して 1 bar / 1 M に換算。温度は pipeline の `conditions` だけから決まる。会合量の単量体のキーは (Hill 式, 電荷, 多重度) で、錯体と単量体の LOT 判定は電荷と多重度を除く(`same_pes(state=False)`)。G のない単量体があれば ΔG_assoc は None。単原子の単量体を含む組成と、反応を持たない錯体では ΔG_assoc は出ない。単量体のアンサンブルには同じ (Hill 式, 電荷, 多重度) の異性体も入る。感度の幅(qs × cutoff 50/100/150)、Boltzmann 重み。ΔG_assoc には幅を出さないが、qRRHO の扱いによる幅は 298 K で約 1.2 kcal/mol(TMA·(HF)₂)。失敗は G = None(`thermo_unavailable`)、LOT の不一致は `mixed_level_of_theory`。composite G = E_SP + (G_GV − E_GV) は `energy_method` を指定したときだけで、同じ構造の SP がない subject は `energy_layer_missing`(G = None)。JobStore のキーは freq の job_key、Hessian の sha、渡す振動数、温度、settings |
| report | rankable な反応だけを ΔG‡ で並べ、感度の幅が重なれば同順位。被覆率(機構別の試行・生成物・陰性理由・FailureKind)。手法パネルで、ΔE_rxn が +1 kcal/mol を超える手法と −1 kcal/mol を下回る手法が両方あれば `method_sign_disagreement`(±1 以内の値は符号を問わない)。method_panel を追記した後は、panel_report の順位が最終である(method_panel.yaml の `conditions` を元の run とそろえた場合。違えば、順位を付けられる反応は thermo_unavailable になる。same_basin などは outcome だけ) |

- 結果の読み方: 停留点レベルを PBE0-D3BJ/def2-SVPD にしたのは、FHF⁻ の De が diffuse 関数のない SVP では 21.7 kcal/mol 過大になるが SVPD では実験値の誤差内に入るからである。障壁を CCSD(T)/def2-TZVPD(SVPD 構造での SP)と比べると、差(SVPD − CCSD(T))は HCN −1.2、HONO +2.0 kcal/mol で(validation.md §3 のプローブ)、誤差の符号は反応で異なる(2 系での見積もり)。順位は ΔG‡ の序数として読み、2 反応の ΔG‡ の差が約 3 kcal/mol 未満なら、順序は手法誤差で入れ替わり得る。ΔG_assoc は BSSE を補正しない SVPD の値で、結合を過大に見積もることがある。amine·(HF)n の結論(プロトンの位置、同一 basin)は PBE0 の PES 上のものである。
- composite G を使うには、元の run-dir を複製し(`cp -a`)、複製に `--run-dir` で sp(`targets: all_minima`)→ thermo(`energy_method`)→ report を追記する。stage id は元と別にする(例: `sp_tz` / `thermo_tz` / `report_tz`)。JobStore も複製されるので、計算済みの SP と GoodVibes の結果は再利用される。FileRef は run からの相対パスなので、複製しても参照は壊れない。thermo の結果の id は stage id を含まないので、複製した run では thermo_tz の値が元の thermo の値を上書きする(複製で元の report を再実行すると composite の値が出る)。sp の既定の targets では単量体が `energy_layer_missing` になり、ΔG_assoc は None になる。SP 層で ΔE‡ ≤ 0 になった反応の ΔG‡ は順位として読まない(これを外す blocker はない)。陰イオンや H 結合の系では TZVPD の層(`pbe0-d3bj_def2-tzvpd` など)を使う。

## 8. JobStore と attempt ladder

- JobStore(`hfauto/execution/jobstore.py`)は run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合し、`adapter.result_type` で復元する(検証に失敗したら miss)。
- 終端的な失敗(input_invalid、method_mismatch、executable_missing、incomplete_output)も記録し、`--retry-failed` で指定した種類だけを再実行する。ジョブ単位の排他は `.lock` を O_EXCL で作る。
- ladder(`hfauto/execution/jobs.py` の `LADDER`。これ以外の自動再試行はしない): timeout と geometry_maxiter は `continuation` で 2 回まで継続(NWChem は最新構造・movecs・drv.hess、xTB は最終構造)、input_invalid(autoz)は Cartesian で 1 回、scf_not_converged は前回の vectors に damping と level shift をかけて 1 回。初期 Hessian ありの opt の継続は、drv.hess を INHESS 0 の restart で引き継ぎ、trust 0.1 で続ける。attempt の timeout は min(設定値, 残り時間) で、残り 60 s 未満なら `BUDGET_EXHAUSTED`。
- 並行実行: `JobRunner` のセマフォで ranks × threads の合計を site の cores 以下に保つ。極小(xTB と DFT)と反応ケースは直列で、CREST と GoodVibes だけを `thread_map` で並列にする。`SiteLock` は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。
- run ディレクトリ: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/diagnostics.json`(conformers と minima が書く。読むコードはない)、`<run>/run_state.json`(実行順の stage 状態と JobStore の統計)、`<run>/resolved_config.yaml`。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に union したもの。done で入力と設定の sha が変わらない stage は飛ばし(再開)、`--from X` は X 以降を stale にする。resolved_config.yaml の code_version は git の sha(追跡中のファイルに差分があれば -dirty、改行コードだけの差は数えない、git がなければ unknown)で、pipeline_id ごとに最後に実行したときの値に上書きされる。done で飛ばした stage を計算した版は残らないので、版をまたいだ結果が要るときは新しい run dir か `--from` で取り直す。

## 9. 設定の 4 分割

層どうしはマージせず、優先順位もない。絶対パスは site ファイルにだけ書く。

| 層 | 置き場所 | 中身 |
|---|---|---|
| site | `configs/sites/wsl_local.yaml` | scratch_root、cores、memory_mb、エンジン別の `EngineSite`(版数の pin、実行ファイル、`ExecutionSpec`、worker の python、scratch_dir) |
| method | `configs/methods/` | `MethodSpec`(JobStore のキーに入る): `gfn2`、`pbe0-d3bj_def2-svpd`(停留点)、`pbe0-d3bj_def2-tzvpd` と `wb97x-d3_def2-tzvpd`(手法パネル)、`ccsd-t_def2-tzvpd`(小さい閉殻分子で opt-in)。停留点の method は minima(dft)と reaction-paths に書き、食い違えば読み込み時に拒否する |
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
| `conformers` | `hfauto.stages.conformer_search:ConformersStage` | species | species | engine*, method*, quick, ewin_kcal, seeds_per_composition, keep_per_state |
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
