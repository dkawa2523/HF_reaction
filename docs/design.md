# hfauto の設計

本書は現在の作業ツリーの実装(`hfauto/`)の説明であり、全体検証の合格を示すものではない。今後の設計は [improvement-plan.md](improvement-plan.md)、実装・統合・検証の進捗は [roadmap.md](roadmap.md) に分ける。対象範囲と使い方はREADME、既知の制限は§10、環境は [environment.md](environment.md)、実計算の証拠は [validation.md](validation.md) にある。古いdocstringの「design §x.y」は凍結した初期設計書の番号であり、今後の仕様を指さない。

## 0. 原則

- 化学の判定は `hfauto/chemistry/` の純関数だけが行う。stage は入力を集めて driver か純関数を呼び、artifact を返す薄いアダプタである。
- 外部ジョブの観測事実は型付きの `Evidence` で受け渡す。化学的な判定は`chemistry/`にまとめ、driverが判定に必要な観測と計算手順を組み合わせる。
- バックエンドは能力 Protocol の実装で、subprocess を起動するのは `hfauto/execution/process.py` だけである。エンジン名の文字列は `hfauto/backends/` と `configs/` の外に書かない。
- 必要な観測が欠けた結果を合格や別理論の値で埋めない。明示した計算救済は来歴を残し、失敗は`FailureKind`付きで項目に閉じ込める(§8)。契約・診断は入力、科学的判断、再利用の整合に必要なものに絞る。
- 停留点(opt・saddle・string・freq)は大域混成か GGA に D3 を足した RKS / UKS で計算する(NWChem の解析 Hessian が使える範囲。RSH・meta-GGA では数値 Hessian になり 3〜4 倍かかる)。meta-GGA 混成(M06-2X)、RSH と CCSD(T) は SP だけで使う(既定の順位のエネルギー層は M06-2X、§7.3)。気相専用で、溶媒和は持たない。

## 1. 層構造と import 契約

依存の向きは `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core`。複数の層で使う型は、使う層のうち最も下の層に置く。

| 層 | 役割 |
|---|---|
| cli | 5 コマンド(`doctor`、`run`、`status`、`report`、`case`) |
| pipeline | 設定の読み込み、preflight、run ディレクトリ(`RunLayout`)、`StageRuntime` の実装、`run_pipeline` |
| stages | 8 つの stage |
| drivers | `MinimumDriver`(`hfauto/drivers/minimum.py`)と `ReactionCaseDriver`(`hfauto/drivers/reaction_case/`) |
| reporting | 順位・被覆率・手法パネルの表と静的 HTML、検証の比較(`reporting/validation.py`) |
| backends | 能力 Protocol(`hfauto/backends/protocols.py`)、registry(`hfauto/backends/engines.py`)、5 つのアダプタ |
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

レビュー `docs/reviews/2026-09-27_platform_review.md` §1.1 の区分(U0〜U9)とコードの対応。

| 区分 | 責務 | 実装 | 出力 |
|---|---|---|---|
| U0 プロトコル | 段ごとの理論レベルと順位の面を決め、LOT をエンジンに渡して出力の観測で照合する(Z>36 の def2-ECP は自動) | `configs/methods/`、`hfauto/core/method.py`、アダプタの入力生成 | `Level`、デッキ |
| U1 入力・電子状態 | 化学種と組成を、対応元素・電荷・多重度・状態ラベルを持つ構造にする。範囲外は入口で 1 回だけ拒否する | structures、`hfauto/chemistry/` の electronic_state・elements・topology | `SpeciesRecord` |
| U2 配座・錯体配置 | screen に渡す候補構造を作る(同一性の確定と窓の選抜はしない) | conformers、`hfauto/chemistry/placement.py`、CREST | 状態ラベルごとの候補 |
| U3 極小の確定 | 同じ LOT で確かめた極小を basin に登録する。鞍点に落ちた入力は mode-follow で解消する | minima、`MinimumDriver`、`hfauto/chemistry/` の identity・selection | `MinimumRecord` |
| U4 反応探索 | 生成物を仮定しない片端探索で、状態の間の低レベルの辺(TS あり・なし)を世代を重ねて出す | explore、`hfauto/chemistry/trials.py`、ReaDuct | `DiscoveryRecord`(辺) |
| U5 仮説・経路 | 極小対から仮説を作り、DFT 極小を固定端とする経路(会合は分離した単量体からの緩和スキャン)を分類して種を渡す | `hfauto/chemistry/` の hypotheses・profile、SCREEN・FIND_PATH | `ReactionRecord`、`BarrierVerdict` |
| U6 鞍点・接続 | 種を 1 次の鞍点に精密化し、別ジョブの freq で TS を確かめ、QRC で両側の極小につなぐ | REFINE_SADDLE・VALIDATE_AND_CONNECT・VALIDATE_INTERMEDIATE | `SaddleClaim`、`ConnectionClaim`、outcome |
| U7 熱化学 | 停留点の G(qRRHO の S と QH の H、電子準位、点群、キラリティ、エネルギー層)と順位の量 δG_eff を出す | thermo、`hfauto/chemistry/thermo.py`(GoodVibes) | `SpeciesThermo`、`ReactionThermo` |
| U8 一点計算・報告 | エネルギー層を反応が読む点(`thermo.reaction_points`)だけで計算し、δG_eff で並べ、手法の幅を列で示す | sp、report、`hfauto/reporting/` | ranking.csv、method_panel.csv、report.html |
| U9 実行基盤 | stage を動かし、ジョブを再利用・継続し、失敗を閉じ込め、コアと予算を管理する | `hfauto/pipeline/`、`hfauto/execution/`、`hfauto/cli/` | manifest、JobStore、終了コード |

## 3. Evidence とゲート

`Evidence`(`hfauto/core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すれば次が成り立つ: 正常終了と SCF 収束、opt / saddle の収束、freq なら 3N − n_external 本の振動数(単原子は 0 本)とその Hessian、原子順序と入力フレームの保持、出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること。満たさなければアダプタは `Failure(kind, reason)` を返す(maxiter で止まった saddle は最終フレームと最後の歩のエネルギー `energy_hartree` も持つ)。opt と saddle の Evidence は最終構造の DFT の勾配 `gradient`(Eh/bohr、3N、入力フレーム)も持ち、最後の勾配のブロックが最終フレームにない出力は `incomplete_output`(`gradient_not_at_final`)になる。極小の停留性はこの勾配で判定する(§6.1)。振動数と正準形式の Hessian(Eh/bohr²)は、どのエンジンでも `hfauto/chemistry/vibrations.py` の補空間射影(並進・回転を除き、主同位体の質量)で求めるので、NWChem の表示値とは低振動モードで数 cm⁻¹ ずれる(G で約 0.02 kcal/mol)。thermo は同じ Hessian を freq の構造そのもので射影し(直線性は点群から)、慣性モーメントだけを点群の対称化した構造からとる(§7.1)。

FailureKind は executable_missing、input_invalid、timeout、nonzero_exit、out_of_memory、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、gate_rejected、error の 11 種。out_of_memory は資源の確保の失敗(SIGKILL の rc −9/137 か、出力の末尾に確保の失敗の文言がある nonzero_exit)、error は `StageRuntime.contain` が閉じ込めた項目の例外(§8)。

化学の閾値は `Policy` の 5 つ(§9)で、数値の許容幅 `qrc_drop` = max(1e-5 Eh, 20 × scf_tol)(同じ幾何での SCF の雑音)と反応モード性の下限 `REACTION_MODE_MIN` = 0.3(質量加重の χ の分布の空白 (0.065, 0.635) の中。[validation.md](validation.md) §11.2)は gates のモジュール定数で、YAML では変えない。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=, state=)` | program・version・method・basis・dispersion・電子温度・電荷・多重度が一致。`numerics=True` なら grid と scf_tol も(停留点・freq・QRC・熱化学)。`state=False` は電荷と多重度の違う点どうし(鎖の井戸と分離した単量体、blocker の LOT の判定、手法パネルの LOT のキー) |
| `spin_ok` | ⟨S²⟩ がないか、\|⟨S²⟩ − S(S+1)\| ≤ spin_tol。freq と、thermo が使うエネルギー層の SP にかける。不合格は blocker `spin_contaminated`。手法パネルでは不合格のエネルギーを使う値を最小・最大から外す(§7.3) |
| `imaginary_tier` | 最低の ν < −saddle_cm1 は saddle、< −noise_cm1 は soft、< 0 は noise(停留点の振動数の段。停留性は §6.1 の信頼領域の認証で別に判定する) |
| `is_minimum(freq, opt=)` | freq が opt の最終構造で同一 PES(numerics を含む。NWChem の freq は opt・saddle の収束 vectors から SCF を始めるので同じ電子状態に留まる)、\|E_freq − E_opt\| ≤ qrc_drop(外れたら `state_mismatch`)、tier が saddle でない(soft・noise は注記) |
| `is_first_order_saddle(freq, saddle=)` | 同じ連結条件で、最低モード < −noise_cm1(なければ `no_imaginary_mode`)、2 本目 < −saddle_cm1 なら `higher_order`(2 本目が −saddle_cm1〜−noise_cm1 なら注記 `soft_secondary_mode`)。モードの大きさ・重なり・端点とのエネルギー比較は見ない。最低モードが −saddle_cm1〜−noise_cm1 の停留点も TS として受理し、この case の TS かどうかは χ と QRC の接続で決める(下限を saddle_cm1 に上げると、\|ν\| < 50 cm⁻¹ の重い回転子の TS を受理できない) |
| `reaction_mode_character(symbols, mode, coords, bonds, gradient)` | 一次の鞍点が仮説の TS であることの必要条件。χ = ‖QᵀL̂‖ ≥ `REACTION_MODE_MIN`(`reaction_mode_chi`)。χ はモードと IRC が定義された質量加重の計量で測る: L̂ は Cartesian の虚モード q(M^-½L を正規化したもの)から作る単位ベクトル M^½q(重みを 2 度かけない)、Q は変わる結合の質量加重の Wilson 伸縮ベクトル M^-½∂r/∂x が張る空間の正規直交基底(SVD)。結合が変わらない仮説は宣言座標の質量加重の勾配との \|cos\|、どちらもなければかけない。結合は case のラベル付きの両端の変化 `Ctx.change()` から取り、QRC の側からは取らない(縮退反応では側の結合変化が消えて見える)。不合格は `not_reaction_mode:<χ>`(向き替えや回転子の鞍点)。限界: 無関係なモードでも χ はおよそ √(k/(3N−6))(k は変わる結合の数)なので、原子が少なく変わる結合が多い反応では χ は弱く、QRC が決める |
| `same_spin_state(ev, ref)`、`energy_spin_ok(energy, freq)` | 同じ構造の 2 つのジョブ(LOT は問わない)の ⟨S²⟩ の差 ≤ spin_tol(どちらかが ⟨S²⟩ を持たなければ判定しない。不合格は `spin_state_mismatch`)。`energy_spin_ok` は freq の `spin_ok` とこの照合の両方で、層の SP と手法パネルのエネルギーが使えるかの唯一の定義(不合格の対象は `spin_contaminated`) |
| `profile.judge` | DFT のプロファイル(`profile.Profile`: 節点の構造・エネルギー・出所・⟨S²⟩)の型と判定はそれぞれ 1 つで、SCREEN・string・会合のスキャン・VALIDATE_INTERMEDIATE の部分プロファイルが共有する。両端は DFT 極小のエネルギー(会合のスキャンでは分離した単量体の和と付加体)。順に: 3 点未満は unavailable(`too_few_points`)、⟨S²⟩ をもつプロファイルで内部の最大点の SCF の枝が跳んでいれば unavailable(`scf_branch_jump`、§7.1)、山と井戸を resolution_kcal の深さから数えて、井戸があれば intermediate、なく山があれば single、どちらもなければ barrierless。barrierless なら最も高い内部節点の両隣の区間の中点に、そのプロファイルのエネルギー関数で節点を足して 1 回だけ分類し直す(長さ 0 の区間には足さない。SP が失敗したら unavailable、`midpoint_single_point`) |
| `connection` | QRC の両側が TS と同一 PES で E_TS − qrc_drop より下の極小に割り付くこと(途中の軌跡は見ない)。組は case の判定の粒度のキー(§6.2)で比べる。期待の組なら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつラベル付きで別構造なら degenerate(結合が変わる縮退反応では、QRC の両側の間の帯付きの結合変化が case の変化 `Ctx.change()` かその逆と等しいこと。違えば `bond_change_missing`)、それ以外で両側が同じ basin なら failed(`sides_same_basin`) |
| `rankable` | outcome が elementary / degenerate / reassigned で、ReactionThermo に blocker がない。barrierless_at_resolution は TST の障壁を持たないので序数の順位に入れず、ΔG_rxn とともに別の表(捕獲律速)に出す |
| `reaction_tier` | connected(ConnectionClaim)> saddle(SaddleClaim)> minima > screening |

claim を作るのは 1 か所だけである。`MinimumRecord` は minima stage と reaction-paths の `Registry`、`SaddleClaim` と `ConnectionClaim` は `ReactionCaseDriver` が作る。blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated、not_stationary)は thermo stage だけが作る。thermo と report は claim を再検証しない(thermo の `same_pes` は例外)。

## 4. payload

各 stage の `<run>/<stage_id>/manifest.json` にはその stage の出力だけを入れる。artifact は `status: success | failed` を持ち、payload は `kind` で判別する frozen なタグ付き共用体(`hfauto/core/records.py`、`extra="forbid"`)である。

| 型 | レコード | 出す stage | 主な中身 |
|---|---|---|---|
| species | `SpeciesRecord` | structures、conformers、minima、explore、reaction-paths | 組成キー、電荷、多重度、構造、source、状態ラベル |
| calculation | `Evidence` | minima、reaction-paths、sp | 1 ジョブの観測事実(id は `calc_` + job_key の先頭 16 桁) |
| minimum | `MinimumRecord` | minima、reaction-paths | basin、tier(screen / dft)、level_key、opt / freq の calc、members、注記 |
| discovery | `DiscoveryRecord` | explore、minima | 状態の間の辺: 機構(nt2 / relaxation / mode_follow)、trial(編集の類)、outcome(product = 辺、unconnected = 出発の状態から届かない辺、negative、failed、not_attempted = 予算で切った類か残った状態)と理由、両端の species(`source_species`、`product_species`、同じ原子順)、低レベル TS(TS のない辺は None)、出発側の端からの ΔE‡・ΔE_rxn、`generation`。DFT 段の mode-follow は鞍点の opt(`ts_calc`)。minima(dft) の入口で落ちた辺は negative として書き直す。`source_minimum`・`ts_imag_cm1`・`electronic_temperature_K` は古い manifest を読むためだけの欄(W7-F3 で削除) |
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、会合なら分離した単量体の極小(`monomers`)、source、低レベル TS(`low_level_ts`、優先順)、検証済みの DFT 鞍点(`ts_calc`)、`BarrierVerdict`(出所 screen / string / scan と、プロファイルの節点の calc id `points`。単量体の和と string の bead は空)、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome`、多段の親の子の id `steps`(W6-1 で埋める) |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None)と、来歴の列の点群 `point_group`・σ・m |
| reaction_thermo | `ReactionThermo` | thermo | δG_eff と感度の幅、そのゼロの点(`reference`: complex / separated)、ΔE‡、ΔE_rxn、ΔG‡、ΔG_rxn、ΔG_assoc、分離反応物基準の ΔG‡、blockers、`energy_level`、注記、捕獲律速の層の判定 `layer_verdict`・`dE_act_path_kcal`(W6-3 で埋める) |
| report | `ReportRecord` | report | 順位の行と表のファイル |

CaseOutcome は elementary_step、degenerate_rearrangement、reassigned_step、multi_step、barrierless_at_resolution、same_basin、out_of_window、unresolved_within_budget、blocked_upstream の 9 種。id は species が `species_<id>`、minimum が `min_<species>_<level_key[:8]>`、分割した子反応が `<parent>_split<n>`。

## 5. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス | 結果 |
|---|---|---|---|
| qm | `energy(scf_guess)`、`optimize(init_hessian, hessian_model, fixed_bond, scf_guess)`、`frequencies(scf_guess)` | `nwchem` → `NWChemEngine`(DFT と CCSD(T))、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path)` | `nwchem_string` → `NWChemString`(ZTS)、`pysis_neb` → `PysisNEB`(xTB の CI-NEB と TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode, scf_guess)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2) | `DiscoveryResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作る。`engines.create` は Protocol を満たさなければ `TypeError`、表にない名前は `KeyError`。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult`、xTB と CREST の `--chrg` / `--uhf`、pysisyphus と ReaDuct の計算器設定)。
- 初期 Hessian は freq の `Evidence` で渡す。saddle は種と同じ構造、optimize は同じ原子順序で各原子が 0.5 Å 以内の構造のものだけを受け付け、外れれば `INPUT_INVALID(hessian_geometry_mismatch)`。xTB は初期 Hessian を使わない。
- optimize の `fixed_bond` (i, j, r) は原子対の距離を r Å に固定する(会合のスキャンの点)。開始構造はその距離を持っていなければならない(`fixed_bond_mismatch`)。NWChem は `zcoord` の `bond i j r constant`、Cartesian の継続(autoz の後)では `constraints` の `spring bond`(k = 20 Eh/bohr²。距離は力 F の向きに F/2k だけずれ、CH3 + O2 のスキャンの最大の力で 7e-4 Å。束縛の項は全エネルギーに入らない)。xTB は `constraint_unsupported` で拒否する。`scf_guess`(qm の 3 つと saddle)は親のジョブ(同じ原子と電子状態。LOT は問わない)の収束した vectors から DFT の SCF を始める。基底が違えば親の基底を deck に `parent` と宣言して射影する(`vectors input project parent guess.movecs output job.movecs`。出力を名指ししないと NWChem は vectors を guess.movecs に上書きし、そのジョブから派生するジョブが guess を失う)。親は、freq なら opt・saddle、層とパネルの SP なら対象の freq、スキャンの点なら前の点、開殻のプロファイルの SP なら端の DFT 極小の opt、高次の鞍点を押した探索ならその TS の freq である(§7.1)。CCSD(T) の参照は atomic guess から始め、xTB は受け取って無視する。`fixed_bond` と `scf_guess` は渡したときだけジョブ鍵に入る(scf_guess は DFT の SCF を始めるときだけ)。
- optimize の `init_hessian` は、呼び出し側が宣言するモデル `hessian_model`(`HessianModel`: positive | as_is、既定 positive)で書く。positive は剛体運動を除いて固有値を max(\|λ\|, 1e-3) にした正定値のモデル H₊(`vibrations.shape_hessian(h, x)`、ジョブ鍵に `hessian_model: positive`)、as_is はそのまま。アダプタは虚モードの数から推測しない。MinimumDriver は高次の鞍点からの側だけ as_is を宣言し、それ以外(一次の鞍点からの QRC と mode-follow の側、信頼領域の緩和、低レベルの錯体と R6 の seed の xTB Hessian)は positive である。NWChem は負の固有値に沿った最小化の歩幅を trust の 0.03〜0.3 倍に切り、BFGS 更新は符号を保つので、負のまま渡すと下りが遅い(acac の QRC 34/27 歩 → 15/14 歩)。高次の鞍点からの側は別の鞍点方向へ出られる必要がある(DME C2v はそのままで 29 歩、モデルでは 207 歩で未収束)。宣言が以前の推測と同じジョブの鍵は変わらない。xTB は受け取って無視する。
- SCINE と pysisyphus は `python -m hfauto.execution.worker <module:function> <job.json>` の子プロセスの中だけで import する。
- 熱化学はエンジンではない。qRRHO は閉じた式なので、`hfauto/chemistry/thermo.py` の `species_thermo` が GoodVibes 4.3.0 を同じプロセスで呼ぶ(ジョブにしない)。版数は preflight が照合し、結果は golden で保証する。

エンジンを足す手順:
1. アダプタ(`prepare` / `parse` / `continuation` と `result_type`)とエンジンクラスを `hfauto/backends/<engine>/` に書く。Level は出力から観測し、振動数は `projected_frequencies` で求める。
2. `hfauto/backends/engines.py` の表に 1 行足し、`pyproject.toml` の `adapters-private` と `adapters-independent` にモジュールを加える。
3. golden: 実出力を `tests/golden/excerpt.py` の規則で抜粋して `tests/golden/data/<engine>/` に置き、出所と sha を `tests/golden/SOURCES.json` に書き、パーサのテストを書く。
4. smoke: `tests/smoke/` に `pytest.mark.real` のテストを書く(`real_engine` fixture)。
5. `configs/sites/wsl_local.yaml` に `EngineSite` を加え、`hfauto doctor` で確かめる。

## 6. driver

### 6.1 MinimumDriver(`hfauto/drivers/minimum.py`)

`relax_to_minimum` の手順:
1. opt(初期 Hessian があれば使う。収束済みの opt を渡せば省く: QRC の側、`opt:reused`)。最終構造が `Registry` の既存の basin(同じ組成・電荷・多重度・level_key)に一意に割り付けば、freq を省いて `known`。
2. 別ジョブの freq → `is_minimum` と `spin_ok`。
3. saddle 級なら ν < −saddle_cm1 の虚モードに沿って `modes.off_saddle` で変位し(各モードは QRC と同じ `modes.qrc_step`: 目標 max(3 × qrc_drop, 3e-4 Eh)、0.05〜0.4 Å。和は 0.4 Å で切る)、下りる(最大 `mode_follow` サイクル。履歴 `follow<n>:<判定>`、側は `opt:follow<n>`)。虚モード 1 本なら ± 両側、2 本以上なら全モードの和の方向へ片側だけ押す。下りは QRC と共有する `Relaxer.descend` 1 つである: 側の opt(変位元の freq を初期 Hessian に)を同時に走らせ、側ごとに、既知の basin なら freq なし、先の側と同じ basin ならその settle を共有(鏡像の側)、それ以外は freq と認証(下の 4)。判定は、両側がラベル付きで別の構造なら `ts_candidate`(縮退転位を含む。両側の結果 opt・freq・注記を返す)、両側が同じ構造なら `replace`、片側だけ動けば `one_side`(その側を採る)。minima stage は `ts_candidate` の両側をジョブなしで同一性だけで登録し(鏡像の側は同じ basin)、側 1 を `source_species`、側 2 を `product_species` とする `mode_follow` の discovery にする(DFT 段では鞍点の opt が `ts_calc`)。元の化学種は入力構造に近い側の basin に注記 `endpoint_was_saddle` 付きで加わる。
4. 停留性の認証: saddle 級でない点では、opt の最終構造の DFT 勾配 g と freq の Hessian から、信頼領域の部分問題 min gᵀs + ½sᵀ\|H\|s(‖s‖ ≤ r_B。Moré–Sorensen、`vibrations.trust_region_step`。剛体運動を除いた Cartesian の全内部モード、\|H\| に床なし)の降下量 ΔE_TR を求め、ΔE_TR ≤ 5e-5 Eh なら停留点と認める(`minimum.not_stationary`)。半径 r_B = BASIN_A·√N(0.05 Å の RMSD)と幅 5e-5 Eh(`BASIN_DE_HARTREE`)は、2 つの極小を 1 つの basin とみなす `identity` の基準そのもので、新しい数はない。半径があるので平らなモードの寄与は傾き × 半径で止まる(以前の無制限の二次モデル Σg²/2\|λ\| は acac の TS で 8.6e-3 Eh、ΔE_TR は 2.1e-4)。B1 の DFT 極小 73 個はすべて ≤ 1.2e-5 Eh、受理した TS 31 本では acac の TS だけが境界の外([validation.md](validation.md) §11.5)。勾配のない Evidence(勾配を記録する前の JobStore の結果)は判定しない。noise_cm1 は停留点での数値の雑音の幅で、認証した点の soft・noise の虚モードは注記(`soft_imaginary_mode`)だけで押さない。
5. 認証に落ちた点は、x + s から freq を正定値のモデルにして 1 回だけ opt し(`tr_relax`)、もう一度 settle する。なお落ちれば注記 `not_stationary` を足し(ほかの注記は残す)、thermo はそれを blocker にする(§7.1)。結果はどちらも `minimum` である。TS の認証は §6.2 の VALIDATE_AND_CONNECT が同じ関数で行う。

`Registry` は組成 × level_key ごとの basin の集合で、同一性の基準は `identity.assign` の 1 つだけである(\|ΔE\| ≤ 5e-5 Eh の候補のうち回転と鏡映を許す置換不変 RMSD が最小のものが ≤ 0.05 Å なら、それに入る。同じ距離なら id の小さい方。重心からの距離を元素ごとに並べた列の RMS の差は RMSD の下限なので、それが 0.05 Å を超える候補は RMSD を計算せずに除く。結果は変わらず、explore の TS の重複判定(S8 で数千の TS)が約 100 倍速くなる)。鏡像は同じ basin である(m は thermo が点群から決める。§7.1)。reaction-paths では `Registry` が DFT 極小の唯一のストアで、`Registry.minima` で minimum_id から記録と代表の最適化構造(FileRef)を引く。basin の構造をメンバーの原子順に写すのは `identity.member_coords` だけである(代表自身なら入力の座標をビット単位でそのまま返す)。ラベル付きの比較(`identity.same_as_labelled`: \|ΔE\| ≤ 5e-5 Eh かつ `mapped_rmsd` ≤ 0.05 Å)は真の回転だけなので、NH3 の反転のような縮退転位は恒等と区別される。minima stage は、screen 極小が持つ species を先に、それぞれ species id の順にジョブを並べ、opt の群を一度に走らせ、既知の basin にも先の opt の basin にもない opt の freq と mode-follow を一度に走らせてから、species id の順に直列で登録する(登録簿は 1 つずつ緩和して登録するのと同じ。偽の 3 つの PES、cores 4 と 1 で確かめた)。同じ stage の中でも登録済みの basin に落ちたジョブは freq を省く。開始構造が先のジョブの厳密な置換・鏡像の像(`identity.is_image`: `carry` の RMSD ≤ `IMAGE_A` 0.005 Å。QRC の − 側の省略と同じ判定。対称な像は 0.0009 Å 以下、それ以外の組は 0.025 Å 以上で、basin の 0.05 Å とは別の基準)であるジョブは計算しない。PES は同種核の置換と反転で不変なので、代表がそのまま極小か既知の basin に落ちたときだけその basin に入る(`image_of:<sid>`)。代表が鞍点・mode-follow の歩・失敗なら自分で緩和する。mode-follow の経過は `<run>/<stage_id>/diagnostics.json` に書く(記録用で、読むコードはない)。

### 6.2 ReactionCaseDriver(`hfauto/drivers/reaction_case/`)

静的な判定は判断表の前に `hypotheses.select` が行い、ジョブも log もない閉じた記録を出す: 端点に DFT 極小がない(宣言した端点が鞍点に落ちたなど)なら BLOCKED(`endpoint_without_dft_minimum`)、両端が 1 つの basin で縮退反応でも会合でもなければ SAME_BASIN(`same_basin`)、反応物の漸近(会合は分離した単量体の和、それ以外は反応物の極小)から生成物までの ΔE が reaction_window_kcal を超えれば OUT_OF_WINDOW(`out_of_window`。`hypotheses.uphill` 1 つ)。窓の外の未宣言の候補も閉じた記録として coverage に出る。両端が 1 つの PES にあることは、設定の読み込みで `PipelineConfig._check_stages` が保証する(minima(dft) と reaction-paths の method が 1 つで、stage は DFT 極小だけを読む)。

`decide(case, state, rules)` は `hfauto/drivers/reaction_case/state.py` の 9 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。判断は `<run>/<stage_id>/cases/<reaction_id>/log.jsonl` に書く(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 接続が elementary / degenerate / reassigned | 完了(`connection:<label>`) |
| 2 | 中間体が両端と別(判定の粒度で) | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(下の「分割の子」) |
| 3 | 最新の DFT プロファイル(SCREEN、会合のスキャンか string)が barrierless | BARRIERLESS(`screen:barrierless` / `scan:barrierless` / `string:barrierless`) |
| 4 | 収束した saddle が未検証(`saddle_pending`。`ts_calc` の case は最初から) | VALIDATE_AND_CONNECT(`saddle_converged`) |
| 5 | SCREEN が有効(会合は常に)で、低レベル経路(会合ではスキャン)が未実行、かつ待っている種がない(近道の種がすべて失敗したら string の前に 1 回だけ NEB) | SCREEN(`screen`) |
| 6 | 最新の DFT プロファイルが intermediate で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_intermediate`、最も低い井戸) |
| 7 | 未使用の種があり、失敗の印 < `max_saddle_attempts` | REFINE_SADDLE(`seed:<source>`) |
| 8 | 会合でなく、印が残り、走らせた string の数 ≤ 印の数 | FIND_PATH(`dft_path`) |
| 9 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

**止める規則**: 件数だけで、壁時計の締め切りは持たない。行 1〜3 は手元の証拠だけで決まる完了で、行 4 以降は新しい計算を始めるか打ち切る。答えを出さなかった探索(接続の失敗、TS の棄却、未収束、経路のエンジンの失敗、scf_unavailable)はどれも失敗の印 1 つ(`CaseState.attempts`)で、claim と接続を消し、`max_saddle_attempts` が残る間は次の種か経路へ進む(以前の行 7 の 3 値の扱いと 2 回目の QRC はない)。印は case ごと(子反応はそれぞれ 0 から)。種 1 つの saddle の探索が 1 つで、高次の鞍点と認証に落ちた TS の継続も次の種として数える。maxiter で止まった探索の継続は同じ印の中で行う(数えると 1 回の停滞ごとに 1 つ使い、W4 の S6 の d403cb0271 が既定の予算 2 で unresolved になる)。FIND_PATH は、失敗か unavailable の string だけが印 1 つで、barrierless の string は case を閉じ(行 3)、intermediate は行 6、single の種は精密化で数える。string は、それまでの string の種をすべて試した後(string の数 ≤ 印の数)だけ次を走らせるので、種を出さなかった string の後は打ち切る(行 9)。分割の深さは `max_split_depth` で切る。1 つの attempt の長さは site の `timeout_s`(ハングの検出)だけが決め、run 全体を止めるのは外からのシグナルである(§8)。

**継続**: 継続は 1 つの仕組み(`actions.continuation`)で、新しい種(SCREEN・FIND_PATH の種と低レベル TS)から 1 回だけ行う(`Seed.depth` 1。継続はもう継続しない)。その点の DFT freq を種の Hessian と SCF の guess にし、case の結合変化(なければ反応方向)に対する χ が最大の虚モードに沿って整える(`shape_hessian`)。始点はその点 x で、その勾配が P-RFO を導く(\|H\| の信頼領域の歩は反応モードに沿って下ってしまうので、x + s は極小にだけ使う)。停留した高次の鞍点だけは、それ以外の ν < −saddle_cm1 のモードの和の方向に 1 回押す(`modes.off_saddle`)。継続を起こすのは 2 つ: maxiter で止まった探索(最後のフレームで DFT freq をとり、同じ試行の中で続ける)と、高次の鞍点か認証に落ちた一次の鞍点(次の種として先頭に入れ、次の印に数える)。継続でもなお認証に落ちる一次の鞍点は受理し、注記 `not_stationary`(thermo の blocker)を付ける。継続の高次の鞍点は棄却する。

**継続の上限**: maxiter で止まった探索の最後の歩のエネルギー(`Failure.energy_hartree`)が、最新の DFT プロファイル(SCREEN、string か会合のスキャン。低スピン結合のスキャンは AP のエネルギー)の最大 + resolution_kcal を超えるなら、継続しない(注記 `saddle:above_path_bound:<screen|string|scan>`)。両端を結ぶ連続した DFT の経路の最大は、探している鞍点のエネルギーの上限なので、それより上へ登った探索は、その鞍点を通り過ぎている(VAL7 の S5 は、プロファイルの最大より 27〜32 kcal/mol 上のフレームから続けていた。[validation.md](validation.md) §8)。エネルギーのない Failure(エネルギーを記録する前の JobStore の結果)とプロファイルのない case には上限をかけない。

**分割の子**: 行 2 の子反応は、親の直後に、それぞれ自分の試行の予算で駆動する。次の子は駆動せず、ジョブも log もない記録にする(`classification.undriven`)。
- 同じ case キー(`hypotheses.pair_key`: 両端の状態が違えば 2 つの状態、同じなら 2 つの極小)の case が駆動済みか待ち行列にある子: その case の駆動を待ち、結論が出ていれば outcome と claim を写す(`same_as:<reaction_id>`)。同じ問いを 2 度探索しない。その case が unresolved で終わったら、結論のなさは共有せず、子を自分の予算で駆動する(待った子は、待ち行列が空になってから)。
- 親の深さが `max_split_depth` に達している子(仮説が深さ 0): unresolved_within_budget(`split_depth`)。

**原子の対応**: 反応の両端の添字は TS を挟んで連続していなければならない。basin は置換について不変な同値類で、代表の species の添字は任意なので、派生する端点には代表を使わず、実際に到達した添字付きの構造を使う。
- 中間体(QRC の reassigned と VALIDATE_INTERMEDIATE): 到達した構造(QRC の側か、緩和した井戸の opt の最終構造)をケースの原子順のまま子反応の端点にする。case が新しい basin を登録したときは、その basin の代表になった case 自身の species(opt・freq を持つ)を使う。既知の basin に落ちたときだけ、その構造をジョブなしの `SpeciesRecord`(`spc_<reaction_id>_intermediate`、source `intermediate`)にして basin のメンバーに加える(`Registry.join`)。この処理は `drivers/reaction_case/connection.py` にある。
- 発見: 仮説の両端は辺(outcome product)自身の `source_species` と `product_species` で、どちらも `_Pool.basin` で DFT の basin に引く(nt2、relaxation、mode_follow で同じ)。TS のない辺は、その source が DFT でも自分の状態を保つときだけ仮説にする。unconnected の辺は仮説にしない。
- case の端点は、basin の最適化構造を端点の species の添字と掌性に並べたもの(`identity.basin_coords`)である。species の構造が basin と同じ結合グラフ(状態ラベル)を持つときは、原子を同じ正準クラス(RDKit の `CanonicalRankAtoms(breakTies=False)`。等分割で、軌道より細かくならない)の原子にだけ対応させるので、species の添字付きの結合が保たれる。basin から遠い構造(DFT の basin に入った xTB の生成物など)でも置換を誤らない。結合グラフが違う(basin への緩和で状態が変わった)ときは元素だけで対応させる。そのときの結合変化は化学的な結果で、添字の誤りではない。したがって R → I → P の添字は連続し、分割の子の結合変化に添字の置換は入らない([validation.md](validation.md) §8)。

**判定の粒度**: case は、両端が異なる粒度で判定する(`drivers/reaction_case/connection.py` の `_key`)。
- 両端の化学状態 (composition_id, state_label) が異なる case(結合が変わる case)は状態で判定する。thermo は状態の G を最小の G とする(§7.2、Curtin–Hammett)ので、それと同じ粒度である。端点と同じ状態の井戸は、basin が違っても中間体にしない。
- 同じ状態の case(宣言したねじれ、配座変化、縮退転位)は basin で判定する。ただし端点と同じ状態の別の basin で、その端との \|ΔE\| < resolution_kcal、かつそれを運んだ経路(プロファイルの井戸なら端までのプロファイル、QRC の側なら TS)に resolution_kcal 以上の山がないものは、分解能では同じ端点とみなす。
- 実際に結ばれた極小は `ConnectionClaim.minima` に残す。elementary の `ReactionRecord.minima` は仮説の両端のままで書き換えない(重複除去と reaction_id を保つ)。

action:
- **SCREEN**: 仮説の低レベル TS(`low_level_ts`。explore か mode-follow の、同じキーの相異なる TS すべて、§7.1 の仮説)を順に xTB freq で確かめ、DFT SP が両端の DFT 極小のどちらよりも resolution_kcal 以上高いものを、その順に種にする(`discovery_ts`。近道は 1 回だけで(`CaseState.shortcut_done`)、3 点はプロファイルではないので判定も ⟨S²⟩ の検査も記録しない。DFT SP の順に並べ替えるのは W6-1)。それが尽きたら DFT 極小の間の IDPP 11 点(自前の実装。pysisyphus の `interpolate.IDPP` は平面の HONO の経路を面内の H–O–N 177.8° に通し、作業ディレクトリにログを書くので採らなかった。[validation.md](validation.md) §11.3)を初期経路に `pysis_neb`(xTB の CI-NEB、端は固定、未収束でも経路として使い、CI から TSOpt)で緩和し、内部 9 点の DFT SP と両端の DFT 極小のエネルギーを `profile.judge` で分ける。single の種は NEB の TS(`screen_ts`)か、山を放物線補間した構造(`screen_hei`)。NEB が失敗したら IDPP のまま分ける。経路の両端は常に DFT 極小なので、低レベルの PES に端点の極小は要らない。barrierless は最も高い内部点の両隣の区間の中点 2 点の DFT SP を加えて分類し直してから受け入れる(`profile.judge`。FIND_PATH・会合のスキャン・VALIDATE_INTERMEDIATE の部分も同じ)。根拠(分析 G3-P7): 単峰のプロファイルなら最大は最も高い節点の両隣の区間にあり、節点が見落とす高さは、山を高さ V・底の幅 w の放物線で近似すると、節点の間隔 h に対して (h/w)²V 以下である。中点でその区間の間隔が半分になり、上限は 1/4 になる(実測の節点の間隔は 0.04〜0.30 Å)。多峰のプロファイルで節点の間に隠れた別の山は保証しない。開殻の case では、これらの SP を端の SCF の解から始め、枝の跳んだプロファイルは unavailable にする(§7.1 の「開殻のプロファイルの SCF の連続性」)。会合(`ReactionRecord.monomers` を持つ case)では、SCREEN は低レベルの経路の代わりに、分離した単量体から付加体までの緩和スキャンを走らせる(低レベルのエンジンは要らない。§7.1 の「会合と低スピン結合の組成」)。
- **REFINE_SADDLE**: 反応方向は 1 つの規則で決める(`Ctx.direction`)。
  - TS の種(低レベル TS の `screen_ts` と `discovery_ts`、継続 `continuation`)は、自分の虚モード(`Seed.mode`)。
  - それ以外の種(`screen_hei`、`path_hei`)は、ラベル付きの両端の結合変化(`Ctx.change()`。ρ・χ・高次の鞍点の反応モードの選択・縮退の QRC の判定が共有する 1 つの変化)から作る ρ = ∇(Σ_切れる r − Σ_できる r)を種の構造で評価したもの。結合が変わらなければ宣言座標の勾配、なければ最も変わる二面角の勾配、それもなければ種に重ねた両端の差(NH3 の反転など)。
  - 経路の接線は使わない。接線には回転子が混ざる(VAL7 の S6 の screen_hei で ρ との cos 0.69)。
  - 初期 Hessian は、種が自分の DFT freq(継続の点の freq)を持てばそれ、なければ種の xTB freq(あれば常に)、なければ種の DFT freq。
  - adapter はそれを H₀ = P·H₊·P − κd̂d̂ᵀ にして渡す(`vibrations.shape_hessian`)。H₊ は剛体運動を除いて固有値を max(\|λ\|, 1e-3) にしたもの、d̂ は方向の内部運動の成分を正規化したもの、P = I − d̂d̂ᵀ、κ = max(d̂ᵀH₊d̂, 0.05 Eh/bohr²)。負の曲率は d̂ の 1 本だけで、ほぼ縮退した固有ベクトルのどれを選ぶかに依らない。低レベルの Hessian の負のモードは \|λ\| として正の部分に入るだけなので、xTB Hessian の採否の検査はない。
  - saddle のジョブ鍵はモデル名(`hessian_model: negative_along_mode`)を含むので、以前のモデルの saddle は再利用しない。注記は `saddle_hessian:<seed_freq|xtb|dft>:<mode|rho|coordinate|chord>`。
  - saddle が maxiter で止まったら(柔らかい反応モードでは、2 本目の負の固有値が十数歩残って停滞しやすい)、同じ試行の中で最終フレームの DFT freq から 1 回だけ継続する(上の継続と上限の範囲で)。
- **VALIDATE_AND_CONNECT**: TS の検証と QRC を 1 つの action で行う(`connection.validate_and_connect`)。
  - 検証: 別ジョブの DFT freq → `is_first_order_saddle` → `reaction_mode_character` → 停留性の認証(§6.1 と同じ `minimum.not_stationary`。\|H\| の全モードでかける: 反応モードを除くと尾根から外れた TS が通る)。χ で拒否した鞍点(`ts_rejected:not_reaction_mode:<χ>`)と虚モードのない点(`ts_rejected:no_imaginary_mode`)は、失敗の印 1 つで QRC にかけず、探索を続ける。仮説の正味の結合変化を動かさない第 1 段の TS(後の段で消える一時的な結合だけを動かすもの)もここで拒否され、その中間体は行 6 の `path_intermediate`(プロファイルの井戸)で拾う。中間体を出すのは、プロファイルの井戸(行 6)と QRC の側(行 2)だけである。`ts_calc`(mode-follow の鞍点か親が検証した TS)はまずそれを検証する(freq は JobStore の再利用)。拒否されるのはこの case の両端に対する χ か freq の失敗のときだけで、そのときはほかの拒否と同じく探索の行に進む。高次の鞍点と、新しい種の認証だけに落ちた一次の鞍点は、上の継続を次の種にする(`ts_rejected:higher_order` / `ts_rejected:not_stationary`)。
  - QRC: 受理した TS の虚モードに沿って ± に変位し(振幅は `modes.qrc_step`: エネルギー目標 max(3 × qrc_drop, 3e-4 Eh) と ν・質量から 0.05〜0.4 Å)、§6.1 の `Relaxer.descend` で下りる(TS の freq を正定値のモデルにした opt を同時に走らせ、側ごとに既知の basin なら freq なし、先の側と同じ basin ならその freq を共有、それ以外は freq と認証)。対称な TS では − の開始構造が + の厳密な像(`IMAGE_A` 以内)なので + 側だけを opt し、− 側はその像を運んだ構造で割り付ける(注記 `qrc:minus_is_image`、`side_calcs` は (plus, plus)、ゲートは同じ)。鞍点に止まった側は、その鞍点の虚モードに沿って TS から離れる向きに 1 回だけ下る(`opt:qrc_past`。W3 の S5 split2 の `unassigned_side`)。側を `Registry` に割り付け(未知なら新しい basin)、`connection` で判定する。両側のキーが端点のキーの組なら elementary、片側だけが端点のキーで他方が別のキーの DFT 極小なら、その側の構造を中間体(上の原子の対応)にして行 2 で 2 つの子反応に分割し(注記 `qrc1:end<i>_to_new_basin`)、その 2 つを結ぶ子に TS を渡す(`ts_calc`。子は検証と QRC を JobStore から再生する)。結合変化の case の中間体は、端点と別の化学状態に限る。どちらも端点のキーでない TS は reassigned のまま。結合の変わる縮退反応は、QRC の両側の間の結合変化(帯付き)が case の変化かその逆と等しいときだけ degenerate にする(メチルの回転の TS は移る H を動かさない)。答えのない接続(側の失敗、両側が同じ basin・同じキー)は失敗の印で、claim と接続を消して探索を続ける。振幅を変えた 2 回目の QRC はない(B1・W3・W4 で 2 回目が答えた case は 0)。多段の親は順位を持たず(outcome `multi_step`)、子反応がそれぞれ順位に載る。
- **FIND_PATH**: ZTS を 1 チャンク(9 beads、maxiter 20)だけ走らせ、収束を問わず `profile.judge` で分ける。初期経路は最新の DFT プロファイルの経路、なければ IDPP。端点は DFT 極小にし、画像を隣へ逐次整列する。種は single のときだけ山の放物線補間(`path_hei`)。失敗か unavailable の string は失敗の印 1 つ。その string の種を試した後も印が残れば、行 8 がその経路から次のチャンクを続ける。
- **VALIDATE_INTERMEDIATE**: 最新のプロファイルの最も低い井戸を `relax_to_minimum` にかける。端点と別のキー(上の分解能の規則で端点とみなす井戸を除く)なら、緩和した構造を中間体にして行 2 で分割する。緩和が失敗したら結果とせず(端点扱いも barrierless もしない)、最も高い山を種にする。端点に落ちたら、もう一方の端からその井戸までの部分プロファイルを最新のエネルギー関数で `profile.judge` にかけ、barrierless なら case を閉じ、そうでなければその最も高い山を種にする(log に `<中間体>:end<k>:<判定>`)。

## 7. 化学プロトコル

### 7.1 stage ごとの既定値

| 段 | 既定値と規則 |
|---|---|
| 共通 | DFT は PBE0-D3BJ/def2-SVPD、grid fine、`convergence energy` 1e-7(opt・freq・SP で同じ)。反復上限は DFT `iterations 100`、WFT の SCF `maxiter 100`。DFT の SCF は親の解から始める(§5 の `scf_guess`。基底が違えば射影)。SCF の段は上限付きの 3 段だけ(P0a): (1) 親の(なければ atomic の)guess + DIIS、失敗したら 1 回だけ (2) atomic guess から `cgmin`(`iterations 40`。NWChem が自分で 50 に上げる)、(3) `cgmin` を外した通常の SCF をその vectors から回して ⟨S²⟩ とエネルギーを読む(cgmin は ⟨S²⟩ を出さない)。成否は (3) だけで決め、cgmin の収束では決めない: (3) が cgmin のエネルギーから qrc_drop より離れれば `scf_unavailable:scf_noise`、sp・freq(親の構造にいるジョブ)で親との `same_spin_state` に落ちれば `scf_unavailable:spin_state`、(2)(3) が収束しなければ `scf_unavailable:scf`(いずれも `scf_not_converged`)。駆動系のジョブ(QRC の側、スキャンの点、継続の saddle)は ⟨S²⟩ が正当に変わる(CH3 + O2 で 1.71 → 0.75)ので scf_noise だけを見る。発散した attempt の vectors からは再開しない。WFT は atomic guess から始め、再試行しない。smear・fon は PES を変えるので使わない。geometry は `units angstrom nocenter noautosym`。autoz の致命的な失敗(`insufficient internal variables` など)は Cartesian で続ける(§8)。Cartesian に切り替えたという無害な注記(`AUTOZ failed to generate good internal coordinates`)は失敗にしない。Z>36 の元素には def2 系の基底のときだけ `<元素> library def2-ecp` を元素ごとに書く |
| structures | 化学種は xyz と SMILES のちょうど一方。SMILES は RDKit ETKDG の 1 配座(同位体と `'.'` は不可)。宣言反応の端点は xyz(SMILES では原子の対応と配座が決まらない)で、両端の電荷と原子順序が一致すること。多重度だけが違う両端はスピン交差で、`spin_crossing_reaction_unsupported`(MECP は探さない。それぞれの面は別の反応として宣言すれば評価できる)。化学種の多重度は必ず宣言する(省略は system の読み込みで拒否。SMILES のラジカル電子数も xyz も多重度を決めない)。元素と多重度のパリティはここで 1 回だけ検査する(範囲外の元素は `unsupported_element`、電子数とのパリティ違いは `INPUT_INVALID`)。組成は conformers の入口で決める: 電荷は成分の和、多重度は宣言値かスピン結合で 1 つに決まる値(決まらなければ `declare_multiplicity`)。低スピン結合の組成は受け付け、その一重項だけを `low_spin_singlet_unsupported` で止める(下の「会合と低スピン結合の組成」) |
| conformers | 単量体: `crest --gfn2 --quick -T <n> --ewin 6 --chrg --uhf`(重原子 3 個以下で回転可能結合のない分子は省く)。トポロジー変化で止まったら停止構造を残し、そこから 1 回だけ再実行する。組成: ゲストを乱数(種は固定)の向きと方向の剛体として vdW 接触 + 0.5 Å に置く seed を `seeds_per_composition` 個作り、先頭の seed から `--nci --quick --notopo <全原子> --noopt` で探索する(状態は hfauto の状態ラベルが決め、CREST はサンプリングだけを行う。`--noopt` は CREST 3.0.2 の初期トポロジー検査が `--notopo` を見ないため)。CREST の候補に入力(単量体)か seed(組成)の状態ラベルを持つものが 1 つもなければ(CREST の失敗を含む)、それらの構造もそのまま出し、理由は項目の Failure(CREST の Failure か `GATE_REJECTED input_state_lost:<ラベル>`)に書く。GFN2 の screen で最適化を通るのは別の結合状態に崩れる seed だけで(S5 は 6 個中 1 個が CH3O2、S6 は 2 個が CH3 + H2O)、接触の seed 1 つではどちらも失敗する(W3-3 のプローブ)。CREST 3.0.2 は CH3·O2(rc −11、SIGSEGV)と OH·CH4(MTD が収束せず rc 1)で失敗した。CH3·O2 は低スピン結合の対で、スピン分極のない GFN2 では対の二重項の SCC が収束しない(四重項なら `--nci` で rc 0)。SIGSEGV がこの SCC の不安定さから来るという因果は推論である(出力に SCC 未収束の記録はない)。OH·CH4 は二重項と一重項の対なので、この機構には当たらない。開殻そのものは原因ではない(六重項の FeCl3·CH4 は rc 0)。選抜は状態ラベルごとに GFN2 エネルギーの低い `keep_per_state` 個。CREST は attempt ディレクトリで `--scratch` なしに実行する |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル。状態ラベルが変わった seed は落ちた basin の members に記録する(失われた状態は explore の TS のない辺になる) |
| explore | 状態は (組成, 状態ラベル)、出発の状態は screen 極小と入力の化学種の状態。第 1 世代は screen の各状態の、最低から `SCREEN_WINDOW_KCAL` 6 以内の極小で走らせ、以後の世代は届いた未実行の状態を、その IRC の端の極小で走らせる(CREST なし。順は screen の状態からの GFN2 のボトルネックが低いものから。順だけで、閾値ではない)。編集の類(`chemistry/trials.py`)は状態ごとに 1 回列挙する: 形成 ≤ 2、切断 ≤ 2、形成 ≥ 1、変化した原子が 6 員環以下の環か和のグラフの鎖の上、形成の対はグラフ距離 ≥ 2、正味の配位数が `fits()` に収まり、生成物のグラフに Lewis 構造がある(形式電荷 ±1、不対電子 m − 1、P {3,5}、S {2,4,6}、周期 3 以上の 15・16 族は電荷で価数をずらす。電荷分離は制限しない: P0f の規則は P0c の 1,5-H 移動の生成物を失ったので棄却した)。出発のグラフには求めない。類は編集のラベル(保つ・できる・切れる)を付けた出発のグラフの正準 SMILES で、グラフの自己同型の軌道そのもの(上限も距離のキーもない)。反応の型は列挙の結果で、機構ごとのコードはない。類は状態の構造のうち drive Σ(r/Σr_cov) が最小のものに作る(1 つの状態に原子順の違う構造があれば(同じ原子の組成を 2 つ宣言した S6)、原子順ごとに列挙して類でまとめる)。形成する分子間の対が接触していなければ、小さい方の断片を、反応する原子が相手の外向きを向くよう回して 1.5 × Σr_cov に剛体で置く(`PLACE_RATIO`、調整しない定数)。3 原子以上の直線分子は決まった 10° の曲げ。予算は件数 1 つ `max_trials`(全世代の合計、既定 3000)で、(ボトルネック、結合変化の少ない順、drive、類)の順に試し、1 世代の状態はボトルネックの低い順に、予算を満たすまでだけ列挙する。切った類と、列挙しないまま残った状態(状態ごとに 1 件)は `not_attempted`。trial_id は組成と類だけから作る。各 trial で NT2 を 1 回。worker は、TS の射影した虚モードが 1 本で IRC の両端が Hessian で確かめた極小のときだけ辺を返す(両端が同じ結合なら `no_bond_change`)。辺の端は、出発の構造と同じ結合ならその構造、さもなければその端を持つ既知の極小(`same_bonding` と `identity.assign`)、なければ新しい species。TS のない辺は、試行の開始構造が別の状態の極小に緩んだもの(開始構造を自分の species にする)と、screen で失われた状態(同じ組成の screen 極小にそのラベルがなく、seed と落ちた basin の間に帯付きの結合変化がある)の最低の seed から落ちた basin への `relaxation`。TS が既出の辺の TS と 1 つの basin なら辺を足さない(`same_edge:<id>`)。辺は出発の状態に近い端から向け、どちらの端も届かなければ `unconnected`。`generation` は出発の状態からの辺の数 + 1。単位は並列に走らせ、入力順に記録する。ReaDuct の `spin_mode` は `restricted_open_shell`、版数の pin は scine-readuct・scine-utilities・scine-xtb-wrapper の 3 つ。範囲外: 立体だけが違う経路は 1 つの類、d ブロックの原子が結合した出発点は類を持たない |
| minima(dft) | 入口(`chemistry/selection.py`。stage は開始構造とジョブだけを持つ): 低レベルの辺(explore の product)ごとに、出発の端・TS・生成物の端の低レベルの構造で DFT SP をとり(構造 × 電荷 × 多重度ごとに 1 回)、`selection.admit` が 1 回だけ決める。SP は、出発の状態から高さの低い順にたどって窓の内側で届いた状態を出る辺にだけとる(高さは道に沿って増えるだけなので、窓の外でしか届かない状態の先の辺は SP なしで `out_of_window`。S10 は予算の後に辺が 850 本ある)。辺の高さは、出発の状態(第 1 世代の辺の source)から辺の連鎖をたどった最も低い道の上の最高点(節点は状態で、状態のエネルギーはその端の構造の SP の最小。配座は辺より低く入れ替わるとみなす)(TS、なければ高い方の端。生成物も入るので ΔE_rxn も同じ数で抑える)で、出発の状態のうち最も低くなるものから測る。高さの順(同じ高さは同じ道の上のエネルギーの幅(通った最も低い状態からの最大の上り)、次に道の辺の数の順。S5 では未緩和の seed が高く、その下の辺はすべて高さ 0 になる)に、`reaction_window_kcal` を超えれば `out_of_window`、組成あたり `max_edges`(既定 6、端の組の数)を超えれば `over_cap`、SP が失敗すれば `not_evaluated:<失敗の種類>`(届かなければ `not_evaluated:unreached`)で、いずれも negative として書き直す(窓の外に黙って落とさない)。判定と高さは diagnostics.json に書く。GFN2 のエネルギーはこの入口で何も決めない。精密化する構造: 反応する組成(宣言反応の端点、入った辺の端、それらの錯体の宣言した単量体 = `thermo.declared_monomers`)だけを、組成 × 状態ラベルごとに `window_kcal` 以内の `per_state` 構造(宣言した端点と入った辺の端は必ず含む)。候補が `per_state` を超える組は上位 `rerank_top` を DFT SP で並べ直してから同じ窓で切る。ほかの組成は diagnostics に `not_reacting`。screen 極小が持たない辺の端は自分の構造から、持つものの後に opt する。`all` は入力の化学種全部。opt は初期 Hessian なしなら `trust 0.1`、ありなら `trust 0.3`、maxiter 100。2 断片以上は `init_hessian` の xTB Hessian を使う(§5 の規則で書く) |
| 仮説 | 優先順は宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補。宣言反応は必ず評価する。それ以外は同じレベルの DFT 極小の対(または 1 つの basin)で、ΔE_rxn ≤ reaction_window_kcal、端点の間で結合が変わるものだけ(件数の上限は minima(dft) の入口の `max_edges` 1 つ)。未宣言のねじれ・配座変化・鏡像化は仮説にしない(Curtin–Hammett)。発見の仮説の両端はその発見の `source_species` と `product_species`(§6.2 の原子の対応)で、両端が同じ basin の発見は、結合の相手が入れ替わるときだけ縮退反応にする。仮説の単位は case キー(`pair_key`。§6.2 の判定の粒度と同じで、両端の状態が違えば 2 つの状態の組、同じなら 2 つの極小の組)で、宣言反応はそれぞれ残し、すでにある仮説のキーの発見は新しい仮説にしない。どの仮説も、そのキーの全候補の低レベル TS を優先順に `low_level_ts` に持ち(持っている TS と構造で 1 つの basin に入るものは除く)、DFT 段の mode-follow の鞍点(検証済み)の最初のものを `ts_calc` にする(宣言反応にも貸す)。縮退・結合変化・ねじれは basin の最適化構造を端点の原子順と掌性に並べた構造で判定する。explore の陰性結果は仮説を棄却しない(GFN2 で障壁がなく発見にならない反応も、宣言すれば評価する) |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`、`moddir 1`(整えた Hessian の負のモード)。maxiter は継続せず、最終フレームと最後の歩のエネルギーを持つ `Failure` を返す(継続は case が §6.2 の上限の内側で 1 回だけ決める)。ZTS は `nbeads 9`、`stepsize 0.05`、`interpol 3`、`freeze1` / `freezeN`、maxiter 20 のチャンクで、NWChem の収束判定は使わない。`pysis_neb` は IDPP を初期経路とする Cartesian の CI-NEB(`opt: lbfgs`、max_cycles 100)と rsprfo の TSOpt。QRC の振幅の上限 0.4 Å は初期 Hessian を受け付ける距離 0.5 Å 以下にする |
| sp | 既定の pipeline(discover、known_endpoints)は paths の後に M06-2X-D3(0)/def2-TZVPD の 1 つ(§7.3)。`methods` の各 LOT で、反応が読む点(`thermo.reaction_points`)だけを計算する: connected か barrierless の outcome の反応の TS と、その自分の点と鎖の井戸の状態(分離した単量体を含む)の DFT 極小すべて。ΔG_assoc だけが読む点(R_sep と barrierless の前駆錯体)は補助で、その SP の失敗は artifact にしない(値が空になるだけで、stage と rc を止めない)。freq の最終構造で計算し、parents に対象(minimum_id か TS の freq calc)を書く。thermo と手法パネルはこの parents だけで SP と停留点を対応付ける。DFT の SP は対象の freq の解から始める(`scf_guess`、基底が違えば射影。同じ freq を共有する対象は SP も 1 つ)。CCSD(T) は閉殻が RHF の ccsd モジュール、開殻が UHF 参照の TCE(P0d: NWChem の TCE ROHF-(T) は ROHF の正準軌道と固有値で f_ov 項を持たない定義で、標準の Watts–Gauss–Bartlett と OH で 0.108、CH3O で 0.195 kcal/mol 違う。TCE の UHF-CCSD(T) は PySCF と 2.1e-7 Eh 以内で一致する)。開殻の CCSD(T) のジョブ鍵は `reference: uhf` を持つ(ROHF の旧 deck の結果を再利用しない)。CCSD(T) は DFT の freq が `spin_ok` を満たす点だけにかける(BS の点には単一参照の状態がない。「低スピン結合でない」という条件は置かない: BS の点は必ず `spin_ok` に落ち、spin_ok を通る低スピン結合の点は CH3OO• のような結合したラジカルで、UHF-CCSD(T) が正しく書く)。凍結する芯は `freeze <n>` で明示する(n は原子ごとの NWChem の `freeze atomic` の芯の軌道数から def2-ECP が置き換えた分を引いた和。I は 4s4p の 4 で、全電子の Br の 3d と同じく 4d は相関させる。I⁻ + CH3I は 9、CH3O は 2)。Kr より重い原子を含む CCSD(T) のジョブ鍵は `frozen_core` を持つ(`freeze atomic` の旧 deck の結果を再利用しない)。WFT の deck は `memory_mb_per_rank` を GA 寄りに分ける(heap 5%、stack 25%、global 70%。NWChem は大きさごとに単位を要する。DFT は `memory total`) |
| thermo | GoodVibes 4.3.0 に自前の振動数を渡す。S は Grimme の qRRHO(`qs`、`cutoff_cm1` 100)、H は QH(Head-Gordon、h_freq_cutoff = cutoff_cm1。truhlar の変種では S が RRQHO なので H と S は別のモデル)。σ・m・直線性は、下の「対称性」の点群 1 つから決める。振動数は freq の構造で射影した Hessian のもの(点群が直線なら 3N − 5 本)、回転定数は点群の対称化した構造から。GoodVibes には多重度 1 を渡し、電子項は `thermo.electronic`(G_el = −RT ln q_el)と `thermo.spin_orbit`(E_SO)だけで入れる。準位は、宣言した化学種の状態の極小なら `electronic_levels`(Λ > 0 の直線ラジカル。OH は oh_ch4.yaml)、単離原子なら NIST ASD の基底項(`electronic_state.atom_levels`、(Z, 電荷) の表: B C O F Al Si S Cl Br I。表の 2S+1 と宣言が違えば structures で INPUT_INVALID)、それ以外はスピン多重項だけ。E_SO は E に足す(ΔE・G・submerged の判定)。錯体と TS ではスピン軌道の分裂は消光とみなす。BH76 の比較(`validation/bh76`)は層の SP の電子エネルギーを読み、E_SO を含めない。スケール因子は `vib_scale` の 1 つ(既定 1.0、振動数と ZPE の両方)。負モードは固定の規則: 極小は全モードを \|ν\| に、TS は最低モードを除いて \|ν\| にする。m = 2 の極小と TS の G に −RT ln 2 を加える(鏡像は同じ basin なので、鏡像対を 1 つとして数える)。σ の比は対称数で入る。種の値は 1 atm で、反応は `standard_states` ごとに換算する。層の SP は `thermo.single_points`(対象と Level ごとに 1 つ。2 つあれば誤り)で対応付ける。状態の G は、その状態の極小のうち反応の端と同じ LOT(`same_pes(state=False)`)で spin_contaminated でないものの最小の G で、鎖のどの井戸でも同じ定義。ただし、その状態の freq の LOT の極小(spin_contaminated でないもの)が 1 つでも層の SP か G を欠けば、状態の G は None(`energy_layer_missing`)で、その反応は `thermo_unavailable` になる(欠けた配座を黙って飛ばさない)。感度の幅は qs × cutoff 50/100/150。行の値(δG_eff、なければ ΔG_rxn)がなければ `thermo_unavailable`、自分の点の LOT の不一致は `mixed_level_of_theory`、自分の点の注記 `not_stationary`(§6.1)は同名の blocker。トンネル補正はない(§10) |
| report | rankable な反応を δG_eff で並べ、感度の幅が重なれば同順位。(T, 標準状態) は report の `T_K`・`standard_state`、なければ thermo の最初の組。ranking.csv の列は rank、reaction_id、outcome、tier、rankable、T_K、standard_state、energy_level、dG_eff_kcal、reference、band_low_kcal、band_high_kcal、dG_act_kcal、dG_rxn_kcal、dG_act_vs_separated_kcal、torsional(結合変化のない段)、blockers、notes。barrierless の行は順位を持たず(blocker `outcome:barrierless_at_resolution`)、report.html の「Barrierless steps」の表(ΔE_rxn、ΔG_rxn、ΔG_assoc。会合は捕獲律速)に出る。エネルギーが 2 つ以上の LOT にあれば `dE_act_panel_min_kcal`・`dE_act_panel_max_kcal` を足す(既定の pipeline では停留点レベルとエネルギー層の 2 つ。手法の幅は列で示し、順位は止めない)。coverage.csv は機構別の試行・生成物・陰性理由・FailureKind |

結合と状態(`hfauto/chemistry/topology.py`): 結合は r < r_thr = r_cov,i + r_cov,j + 0.4 Å(Cordero の共有結合半径)で 1 つの構造だけから決める。状態ラベルは、断片の Hill 式を '+' でつないだものと、接続の分子(全結合を単結合、sanitize なし)の RDKit 正準 SMILES の sha256 の先頭 16 桁。原子クラス(`atom_classes`)は `CanonicalRankAtoms(breakTies=False)`。RDKit の版はラベルを決めるので環境の一部として固定する([environment.md](environment.md))。FHF⁻ と I3⁻ は 1 断片、ハロゲン結合(I···N 2.8 Å)は非結合。イオン–双極子錯体は許容値の内側なら 1 断片になり、GFN2 の強い H 結合錯体は閾値の近くでラベルが分かれうる(DFT 極小には影響しない)。2 つの構造の間の結合変化(`bond_changes`)は、r − r_thr が一方で ≥ +`RESOLVED_A`、他方で ≤ −`RESOLVED_A`(0.1 Å、モジュール定数)の対だけで、向きに対称である。1 つの basin の中をしきい値が横切るだけの違いは変化としない(根拠は S19 の N···H で、GFN2 と PBE0 の差が 0.15 Å)。結合変化の定義はこの 1 つで、結合が同じかは `same_bonding`(= 結合変化なし。推移的ではない)1 つで判定する。仮説と分割の子のねじれの判定(`torsional`。未宣言のねじれは仮説にしない)、会合の形、ρ と χ の結合、explore の辺の端の同定と失われた状態がこれを使う。結合グラフと状態ラベルは帯を持たないので、結合変化のない 2 つの構造が別の状態ラベルを持つことがある。ラベルが割れても別の状態になるだけで、別の状態を誤って 1 つにまとめることはない。立体(U1-P4): 断片ごとに正準の原子順で `DetermineBondOrders`(電荷 0)を走らせ、Lewis 構造があれば、電荷のない原子の間の二重結合の E/Z と四面体中心を 3D から読む。Lewis 構造がなければ(ラジカル・イオン)四面体だけ。鏡像は小さい方のラベル(min(L(x), L(−x)))にまとめる。四面体中心が 1 つだけなら立体なし(鏡像が同じ状態)、TB/OH/SP の印は捨て、電荷のある原子への二重結合(1 つの共鳴構造の結合)は E/Z を持たない。立体があれば正準 SMILES に立体の文字列を足して hash する。立体のない構造のラベルは W4 と同じ。

対称性(`hfauto/chemistry/symmetry.py` の `analyze`): 熱化学の点群は、freq の構造と Hessian からここで 1 つだけ決める。
- 候補: 同種原子の置換(原子間距離の変化 ≤ `CANDIDATE_A` 0.15 Å の後戻り探索。探索を打ち切るだけで受理には使わない)× 真回転と回映の写像。直線の候補は自前で対称化する(重心を通る主軸への射影が C∞v、それと反転の平均が D∞h)。libmsym(pymsym 0.3.5)は受理した群の名前を付けるだけ(VAL9 の S6 の TS には、どのしきい値でも候補を出さなかった)。

- 受理: 群の対称化(像の平均)y への変位 d を質量加重し、外部運動を射影で除いたとき、½dᵀ\|H\|d ≤ ½ħω_d(ω_d² = dᵀ\|H\|d / dᵀMd)。対称化のエネルギーがその方向の零点振幅より小さければ、その方向の歪みは零点運動の中にあり、平均の構造は高い群を持つ。例外: 開殻が縮退表現を持つ群(位数 3 以上の操作、直線群)に入るのは、`identity` の雑音の基準(5e-5 Eh、0.05 Å)を満たすときだけ(Jahn–Teller の歪みを零点で平均しない)。

- 選択: 受理した巡回群を上昇の小さい順に、閉包が受理される間だけ合わせる(直線群が先)。構造そのもの(C1)は常に受理する。熱化学は E と振動数を x から、σ・m・慣性を群からとる。

- 規約: σ は、受理した点群の対称操作のうち恒等操作と真回転(C_n)の数である(直線は C∞v 1、D∞h 2。Fernández-Ramos et al. 2007、Gilson & Irikura 2010。名前の表は使わない)。m は、点群に回映操作(鏡映、反転、S_n)がなければ 2、あれば 1(Fernández-Ramos et al., Theor. Chem. Acc. 118, 813 (2007))。直線性と外部自由度の数も同じ点群から決める。basin は鏡像を含むので、m = 2 の basin は鏡像対の両方を表す。m は鏡像の basin の数え方と同じ規則である。
- 点群のない場面(エンジンが自分の freq を解析するとき、試行の方向、整えた Hessian)だけは、直線性を外部自由度の数値の階数(`vibrations._EXTERNAL_RANK_TOL`)で決める。

会合と低スピン結合の組成(分析 §3 X4。`electronic_state.low_spin_coupled`、`hypotheses`、`drivers/reaction_case/paths.py`):
- **低スピン結合の組成**: 開殻の成分が 2 つ以上あり、宣言した多重度が高スピンの結合 Σ(m_i − 1) + 1 より小さい組成(CH3• + O2 の二重項。四重項は高スピン)。拒否しない。結合した異性体(CH3OO• など)は UKS で正しく書けるので、その間の反応は通常の case である。成分が離れた領域だけは UKS が broken-symmetry(BS)解になり、⟨S²⟩ は (S_A − S_B)(S_A − S_B + 1) + 2S_B(CH3·O2 で 1.75。二重項 2/3 と四重項 1/3 の混合)になる。これは誤りではなく、この電子構造から必ず出る値である。そのため非結合の錯体の極小は `spin_ok` で落ち、熱化学には使えない。
- **一重項**: 一重項は RKS で計算するので BS 解がない。このクラスの一重項は入口で止める(`low_spin_singlet_unsupported`。OH• + OH•、O2 + O2 の一重項)。再結合の生成物は単量体として宣言すれば計算できる。
- **会合の仮説**: 反応物から生成物への結合変化が、反応物の 2 つの断片を結ぶ形成 1 本だけで切断がなく、反応物の断片が宣言した組成の単量体と一致し(`thermo.separated_states`)、各単量体の状態の DFT 極小が錯体と同じ LOT にあれば、その仮説を会合にする。起動条件は仮説の形だけで決め、出所(宣言、explore、R6 の relaxation)には依らない。反応物側は分離した単量体である(`ReactionRecord.monomers` に各状態の最低の DFT 極小を個数だけ並べ、`reactants` もその組成にする)。錯体の極小は `minima[0]`(仮説の同一性、重複除去、状態の参照)に残るが、エネルギーと熱化学の端点にはしない。障壁のない会合では反応物は漸近であって極小ではない。低スピン結合の組成では、錯体の極小は BS の汚染による人工物でもある(S5 の CH3···O2、C···O 2.84 Å、⟨S²⟩ 1.71)。会合の形は錯体の DFT 極小の構造で判定する。錯体が自分の DFT 極小を持たず付加体に落ちた宣言反応(BH3 + NH3)は、錯体の入力構造で判定する。このとき `minima[0]` は付加体の極小で、case の反応物端はその入力構造である(`hypotheses.select` の same_basin は会合には効かない)。
- **スキャン**(行 5 の SCREEN): 形成する原子対 (i, j) の距離 r_P + k·0.214 Å(`SCAN_STEP_A`、r_P は付加体の値)の点を、まず k = 7 から 1 まで内側へたどる(最初の 7 点の距離と順は以前と同じ)。最初の点は、付加体の j を含む断片(錯体での断片)を i→j の方向に剛体で引き離した構造で、次の点は前の点の最適構造を同じように押し込んだ構造である。次に最も外の点から 1 点ずつ外へ伸ばし(その点の構造と SCF から)、外端の点が単量体のエネルギーの和から resolution_kcal 以内に入ったら止める(前駆錯体を通り過ぎ、スキャンの外に分解能の山や井戸を残さない。`SCAN_MAX_POINTS` 24 点までに届かなければ unavailable、`scan_reach`)。各点は r_ij を固定した opt(§5 の `fixed_bond`)で、SCF は前の点の vectors から始める(BS の対が連続した枝に留まる)。点の名前 `scan<k>` の k は付加体からの歩数。山はその前の井戸から局所に数え、障壁の値は反応物の漸近(単量体の和)から測る(`profile.judge`)。外側から始めた以前の固定の到達 r_P + 1.5 Å(`SCAN_REACH_A`)は削除した。
- **判定**: プロファイル [Σ E(単量体)、スキャンの点、E(付加体)] を `profile.judge` で分ける。単調なら barrierless_at_resolution(行 3、`scan:barrierless`)。山があれば最大点を放物線補間した構造(`scan_hei`、種の出所は `path_hei`)を種にして REFINE_SADDLE に進み、その後は通常の case と同じである。string は走らせない(string は 2 つの極小を結ぶもので、単量体の端は極小ではない)。点が 1 つでも失敗したら unavailable(`scan_point`)で、点を落として判定することはしない。中点の密化は他と同じく 1 回: 新しい節点は隣の 2 点の線形補間の中点での SP で、外側の点の SCF の解からその点の座標系で始め、低スピン結合の対では AP で射影する(注記 `scan_mid:ap:<多重度>`)。単量体の和は最初のスキャン点と同じ構造なので、その区間には足さない。
- **AP**(分析 G2-P3): 低スピン結合の組成(単量体の多重度と case の多重度で `low_spin_coupled`)のスキャンだけ、各点の最適構造で高スピン(多重度 Σ(m_i − 1) + 1)の SP を 1 本とり、Yamaguchi の近似スピン射影 E_AP = E_BS + α(E_BS − E_HS)、α = (⟨S²⟩_BS − S(S+1)) / (⟨S²⟩_HS − ⟨S²⟩_BS)(S は宣言した多重度のスピン。Yamaguchi et al., Chem. Phys. Lett. 149, 537 (1988))で高スピンの混入を除き、そのエネルギーで判定する。高スピンの SP は各点で並列に atomic guess から解き(注記 `scan:ap:<多重度>`)、失敗したか `spin_ok` に落ちた点は BS のエネルギーのまま(注記 `scan<k>:ap_skipped:<理由>`)にする。付加体に近い点は ⟨S²⟩_BS ≈ S(S+1) で α ≈ 0 なので、AP の曲線は付加体へ連続につながる(CH3·O2 のプローブで、C–O 2.84 Å の BS −1.06 → AP −1.94、1.45 Å の −36.79 → −37.03 kcal/mol)。AP は判定のエネルギーだけで、`spin_ok`、停留点と熱化学の blocker、エネルギー層、手法パネルには入れない(BS の錯体の極小は `spin_contaminated` のまま)。高スピンの宣言(四重項の CH3·O2)、開殻の成分が 1 つの組成(H + C2H4)、閉殻(BH3 + NH3)にはかけない。
- **範囲外**(分析 §0.4): 会合の速度定数(VRC-TST、圧力依存)。多参照法による定量的なエネルギー。低スピン結合の対で答えるのは、AP-UKS による定性的な障壁の有無と ΔG_assoc(分離した単量体と付加体の G、どちらも純粋なスピン状態なので AP に依らない)までで、BS の点(錯体の極小、スキャンの点)は CCSD(T) では校正しない(freq が `spin_ok` に落ちるので、手法パネルも CCSD(T) の SP をとらない)。低スピン結合の一重項。形成する結合が 2 本以上の付加(環化付加。1 本の結合のスキャンではたどれないので、通常の仮説のまま)。多重度の違う端点のスピン交差(MECP)。

開殻のプロファイルの SCF の連続性(分析 G2-P1、`Ctx.sps`、`profile.judge`):
- **端の解からの SP**: 開殻(多重度 > 1)の case のプロファイルの SP(近道の低レベル TS、SCREEN の内部点、barrierless の中点)は atomic guess から独立に解かず、構造が近い端(`mapped_rmsd`)の DFT 極小の opt の vectors から始める(`energy(scf_guess)`)。両端の ⟨S²⟩ が spin_tol より離れている(スピンの結合が違う)ときだけ、両端の guess で 2 回解いて低い方をとる。どの点も端から直接始めるので、点の SP は並列のまま走る(前の点からの逐次の継続はしない)。NWChem の `vectors input` は MO 係数を回転しないので、guess は計算したときの原子順と向きで使う: 各点をその端に重ね、端を guess の opt の最終構造に運ぶ置換と回転で運んでから計算する(`identity.carry`。エネルギーと ⟨S²⟩ は変わらず、ジョブ鍵だけが変わる)。閉殻の case の SP は変わらない。
- **理由**: atomic guess では、点ごとに別の SCF の枝に落ちて偽の山ができる。S5 の CH3···O2 → CH3OO• の IDPP では、端の BS(⟨S²⟩ 1.71)に近い点がほぼ純粋な二重項の枝(0.76〜0.81)に落ち、+28.9 kcal/mol の山になった。端の BS の vectors を直接渡すと、⟨S²⟩ は 1.71 → 0.75 と単調で、山は +0.80 kcal/mol(解像度未満)だった([validation.md](validation.md) §8)。
- **枝の跳び**: ⟨S²⟩ をもつ開殻のプロファイルのすべて(SCREEN、中点を足したもの、会合のスキャン(単量体の和は ⟨S²⟩ なし、点は BS の値、付加体は opt の値))で、内部の最大点の ⟨S²⟩ が両隣の値の間になく、かつどちらかの隣と spin_tol より離れていれば、分類せずに unavailable(`scf_branch_jump`)にし、種を作らない(`profile.branch_jump`、`profile.judge`)。端の ⟨S²⟩ は DFT 極小の opt のもので、両端のスピンの結合は問わない。跳んだ点を落として判定し直すことはしない。理由: 本物のラジカルの TS では UKS のスピン汚染が滑らかな極大をとるので、最大点の ⟨S²⟩ はしばしば両隣の外にあるが、隣との差は小さい(記録のプロファイルで ≤ 0.022。H + H2 の TS で 0.767、両隣 0.757)。SCF の枝の跳びでは差が 0.77〜0.95 になる(VAL7 の S5 で 0.760、隣 1.711)。spin_tol(0.1)はその間にあり、新しいしきい値は置かない([validation.md](validation.md) §8)。Coulson–Fischer 領域の急な変化(S5 のスキャンの 2.3 → 1.9 Å で 1.46 → 0.77)は単調なので、最大点は両隣の間にあり跳びにならない。string(bead は NWChem の string が解き、⟨S²⟩ を持たないので guess も渡さない)には、この判定をかけない。近道の 3 点はプロファイルではないので判定しない。

### 7.2 順位の量 δG_eff

反応は共通のゼロを持つ点の鎖 [R_sep?, R, TS, P, P_sep?] で評価する(`thermo.reaction_points`。sp・thermo・手法パネルが同じ定義を読む)。井戸の点は状態の組(分子ごとに 1 つ)で、その G は各状態の G(反応の端と同じ LOT の、状態の最小の G。速い配座平衡の下で状態は 1 つ)の和に、点の分子数 n_k の標準状態の換算を足したものである(TS は 1 分子)。

- δG_eff = max_k(G_k − min_{井戸 j ≤ k} G_j)(`thermo.chain_barrier`)。[R, TS, P] では max(G_TS, G_R, G_P) − G_R に等しく、単分子の素過程の値は変わらない。TS が自分より前の最も低い井戸より下にあれば律速ではない(微視的可逆性)。順方向か逆方向の ΔE0‡ = ΔE‡ + ΔZPE‡ が自分の点で 0 以下なら鞍点の TST は意味を持たないので、TS を鎖から外し、注記 `submerged_barrier` を付ける(blocker ではない)。
- 分離した点 R_sep・P_sep は、端の断片の状態ラベルが宣言した組成の単量体と一致するときだけ持つ(`thermo.separated_states`)。断片の極小を自動では作らず、宣言した単量体の G が欠ければ、錯体に戻さずに値を空にする。会合の反応物端は分離した単量体(`ReactionRecord.monomers`)で、潰れた錯体は R にしない。分割の子の端は親の鎖の中にあるので、分離した点を持たない。
- ゼロは G_ref = min(G_錯体, ΣG_単量体 + 換算) で、`reference`(complex / separated)に出す。SN2 は錯体、OH + CH4 は分離した単量体がゼロになる。
- `dG_act_kcal`・`dG_rxn_kcal`・ΔE は反応の自分の点(`thermo.participants`)の値、ΔG_assoc = G(R) − G(R_sep)、`dG_act_vs_separated_kcal` = G_TS − G(R_sep) で、表示だけに使う。
- barrierless_at_resolution は TST の障壁を持たないので δG_eff はなく(序数の順位の外。§3 の `rankable`)、ΔG_rxn(前駆錯体があれば ΔG_assoc も)を出す。会合の速度は捕獲律速で、範囲外である(§10)。

### 7.3 エネルギー層と手法パネル

既定の順位は M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD である。discover と known_endpoints は paths の後の sp stage で、順位に使う点(§7.1 の sp)を M06-2X-D3(0)/def2-TZVPD(grid fine、`convergence energy` 1e-7)で計算し、thermo の `energy_method` でその SP をエネルギー層にする: composite G = E_SP + (G_GV − E_GV)。構造と振動・回転の項は PBE0-D3BJ/def2-SVPD の freq のままである。層の定義は thermo の `energy_method` の 1 つだけで、ranking.csv の `energy_level` 列に `m06-2x-d3zero/def2-tzvpd` と出る。SP の対象は鎖の井戸の単量体を含むので、分離した点の値も同じ層で出る。

- **正直な停止**: 層の SP がない対象(SP の失敗を含む)は `energy_layer_missing`(G = None)で、その反応は `thermo_unavailable` で順位から外れる。PBE0 のエネルギーで代用しない(層の混ざった δG_eff は序数の意味を失う)。層の SP は対象の freq の解から始め(射影)、`gates.energy_spin_ok`(freq の `spin_ok` と、SP が freq と同じスピン状態にあること `same_spin_state`)を満たさない対象は `spin_contaminated` である。以前は atomic guess から解いていたので、S5 の BS 錯体の層の SP は ⟨S²⟩ 0.759 の別の解に落ち、ΔE_rxn が −83.15 kcal/mol になっていた。射影した guess では BS 解(⟨S²⟩ 1.735)に留まり −36.0 になる(その極小は freq が `spin_ok` に落ちるので `spin_contaminated` のまま)。
- **層の選び方**(分析 [2026-09-30_remaining_issues_analysis.md](reviews/2026-09-30_remaining_issues_analysis.md) §2.6、§3 X7): PBE0 は非局在化誤差で電子の広がった TS を過安定化する(障壁の平均符号付き誤差は H 移動 −4.2、重原子移動 −6.6、求核置換 −1.9 kcal/mol。Zhao & Truhlar, J. Phys. Chem. A 109, 2012 (2005))。誤差の主因は汎関数で、SVPD → TZVPD の変化は ≤ 1.1 kcal/mol だった。4 反応の ΔE‡ の CCSD(T)/def2-TZVPD との平均絶対誤差は PBE0/SVPD 2.68、PBE0/TZVPD 2.36、ωB97X-D3/TZVPD 1.56、M06-2X/TZVPD 1.03 kcal/mol(§7.4 の表。ROHF-CCSD(T) で測った値で、CH3O を UHF にすると各汎関数の差は約 0.1 動く)。TMA·(HF)₂(17 原子、327 基底関数、4 rank)の 1 点は M06-2X 255 s、PBE0/TZVPD 312 s、ωB97X-D3 1,689 s。M06-2X の grid は fine と xfine で障壁の差が ≤ 0.01 kcal/mol。BH76 の部分集合(P0e)での選び直しは [validation/bh76/table.md](../validation/bh76/table.md)(`table.py` が `measured.json` から生成)にある。選定の行の類別 MUE(PBE0 構造での電子障壁、kcal/mol)は M06-2X HTBH38 1.20・NHTBH38 0.68、revM06 1.91・0.60、ωB97X-D3 3.60・1.44 で、g は全行 ≤ 0.27(geometry-limited の行はない)。事前の規則(選定の類ごとの MUE が現行以下、開殻の弱結合錯体 H···C2H4・CH3···O2・OH···CH4 で SCF の失敗 0、費用 2 倍以内)を満たす候補はなく、M06-2X を残す。兼ね合い: revM06 は原子 guess から 3 つの錯体をすべて 1 段目で解き(M06-2X は H···C2H4 で両段とも失敗)、費用 1.12 倍だが H 移動の MUE が 0.7 悪い。ωB97X-D3 は費用 2.32 倍(TMA·(HF)₂ では 6.6 倍)。
- **手法パネル**(`method_panel`): 既存の run に `--run-dir` で追記し、panel_sp(PBE0/def2-TZVPD、ωB97X-D3/def2-TZVPD)→ panel_report。表示だけで順位は変えない: method_panel.csv に LOT ごとの ΔE_rxn・ΔE‡、ranking.csv に ΔE‡ の最小・最大が出る。`energy_spin_ok` に落ちるエネルギーを使う値は行に残し、notes に `spin_contaminated:dE_rxn`・`spin_contaminated:dE_act` と書いて最小・最大から外す(順位と同じ定義)。パネルの DFT の SP も対象の freq の解から始める。パネルは反応の自分の点(`thermo.participants`。会合は分離した単量体から)を読み、LOT は電荷と多重度を除いて対応させる(列の `level_key` は生成物の Level)。TS は connected の outcome だけが持つ。層の定義を 2 つにしないため、パネルは thermo を持たない。
- CCSD(T)/def2-TZVPD は校正用で、パネルの `methods` に足す。5 重原子で 1 点約 23 分(4 ranks × 2,000 MB、マロンアルデヒド 227 基底関数の実測)。多参照性のゲートはない。

### 7.4 結果の読み方

- **順位の主張**: M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD の qRRHO δG_eff の序数として読む。CCSD(T)/def2-TZVPD(同じ PBE0/SVPD の停留点での SP)との実測の差は、障壁で平均 0.8、最大 1.8 kcal/mol(5 反応、下表)。ΔE_rxn は層で良くならない(平均 1.2、最大 2.3。PBE0/SVPD は 1.0、2.3)ので、障壁のない段(序数の順位の外で、ΔG_rxn を示す)の ΔG_rxn は停留点レベルと同程度の精度である。δG_eff は CCSD(T) の composite と比べ SN2 +0.46、HCN −1.79、CH3O• −0.82(UHF。ROHF では −0.72)、HONO +0.09(validation.md §3.3)。参照そのものも CBS から約 0.5 ずれる。HTBH38/NHTBH38 の文献の MUE は約 1.2(Zhao & Truhlar, Theor. Chem. Acc. 120, 215 (2008))。1 kcal/mol 未満の差は序数を保証しない。遷移金属・溶媒・多参照性の強い系は未検証である。

  | kcal/mol(括弧は CCSD(T) との差。B1 の値。開殻は UHF-CCSD(T)、OH + CH4 だけは B1 にパネルがなく ROHF-CCSD(T) の旧値) | ΔE‡ CCSD(T) | ΔE‡ PBE0/SVPD | ΔE‡ M06-2X | ΔE_rxn CCSD(T) | ΔE_rxn M06-2X |
  |---|---|---|---|---|---|
  | HCN → HNC | 47.81 | 46.62(−1.19) | 46.02(−1.79) | 15.22 | 12.92(−2.30) |
  | CH3O• → CH2OH• | 33.47 | 33.01(−0.46) | 32.66(−0.81) | −7.81 | −8.41(−0.60) |
  | HONO trans → cis | 11.65 | 13.67(+2.02) | 11.64(−0.01) | 0.47 | −0.04(−0.50) |
  | Cl⁻ + CH3Cl(錯体基準。ΔE_rxn は分離 → 錯体) | 13.44 | 10.45(−2.99) | 13.90(+0.46) | −10.60 | −11.27(−0.67) |
  | OH + CH4(分離基準) | 6.96 | 0.79(−6.17) | 5.82(−1.14) | −12.77 | −14.54(−1.77) |

- **PES の形は PBE0 のまま**: 鞍点の有無と位置、PES の平らさ、SCREEN・FIND_PATH・会合のスキャン・行 3 の barrierless の判定は PBE0-D3BJ/def2-SVPD の PES で決まり、層では確かめない。PBE0 に鞍点がない(山が解像度未満の)反応の障壁は層では取り戻せない(H + C2H4: PBE0 の山 +0.60、参照 1.72 kcal/mol。validation.md §4)。PBE0 の非局在化誤差は障壁の値だけでなく PES の形も変え、H 引き抜きの TS を早く平らにしうる(推論。S6 の平らな PES、G1)。停留点の汎関数は変えない: M06-2X と RSH では NWChem が解析 Hessian を使わず freq が高くつき、SP//PBE0 の構造の誤差は S6 で +0.31 kcal/mol(PBE0 の QRC の軌跡上の CCSD(T) の最大と PBE0 の TS の差)と、汎関数の誤差より 1 桁小さい。
- **発見の段の窓は層より前のエネルギーで切る**: minima(dft) の `window_kcal` 6.0 は GFN2(混んだ状態は PBE0/SVPD の SP)、`reaction_window_kcal` 40 は minima(dft) の入口では DFT//xTB の高さ、仮説では PBE0/SVPD の ΔE(反応物の漸近から)で決まる(explore は GFN2 で何も決めない)。窓で落ちた構造と仮説は層に届かず、層で順位に戻らないので、窓は停留点レベルの誤差より広くとる。PBE0/SVPD の CCSD(T)/def2-TZVPD との差は障壁で最大 −6.2(OH + CH4)、S6 の分離基準の ΔE_rxn で −2.2 kcal/mol(PBE0 −15.02、W3 の freq のエネルギー)で、40 はその数倍ある。6.0 は同じ状態(同じ結合)の配座の差に効くので、結合の組み替えに伴う非局在化誤差は主に入らない(推論。配座のエネルギー差の誤差は GFN2 でも PBE0 でも測っていない)。
- 低スピン結合の組成のスキャンは AP-UKS のエネルギーで判定するので、障壁の有無の定性的な判定として読む。AP は BS 解への高スピンの混入だけを除く 2 状態の近似で、多参照性(静的相関)は直さない。
- 電荷系では PBE0 の非局在化誤差で電荷の広がった TS・錯体が低く出やすく(SN2 の ΔE‡ −2.99)、D3 は電荷に依らない。層はこの値を直すが(+0.46)、PES の形は PBE0 のままである。
- 停留点レベルを def2-SVPD にしたのは、FHF⁻ の De が拡散関数のない SVP では 21.7 kcal/mol 過大になるが SVPD では実験値の誤差内に入るからである。
- ΔG_assoc は BSSE を補正しない値で、結合を過大に見積もることがある。qRRHO の扱いによる幅は 298 K で約 1.2 kcal/mol(TMA·(HF)₂)。
- 電子分配関数は、宣言した準位と単離原子の NIST の基底項のほかはスピン多重項 2S+1 だけである。錯体と TS のスピン軌道の分裂は消光とみなし、軌道縮重を持つ未宣言のラジカルは 2S+1 で数える。

## 8. 実行基盤

- **JobStore**(`hfauto/execution/jobstore.py`): run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、結果を変える入力(method、構造の指紋、パラメータ、入力ファイルの sha)をすべて含み、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合する。Evidence の来歴は `job_key` と `FileRef`(run ディレクトリからの相対パスと sha256。作るのは `JobStore.file_ref` だけ)の 2 つで持つが、ジョブの出力と Hessian の FileRef は常にそのジョブの `jobs/<k[:2]>/<key>/attempt_NN/` を指すので、2 つは食い違わない。JobStore は run の中の記憶で、ladder の最後の失敗も(saddle の最後の構造ごと)記録して同じ鍵に同じ結果を返し、`--retry-failed` で指定した種類だけを再実行する。環境に依る失敗 `VOLATILE`(timeout、executable_missing、out_of_memory)は記録せず、以前に記録したものも読まない(site を直せば `--retry-failed` なしで再実行される)。パーサの例外はジョブの失敗ではなく項目の error(下の閉じ込め)で、記録しない(パーサを直せば同じジョブを読み直す)。run の間の決定性は保証しない(並列の NWChem と CREST はビット単位では再現しない)。
- **ladder**(`hfauto/execution/jobs.py` の `LADDER`。これ以外の自動再試行はない): 駆動系ジョブの継続はすべて最後の frame から: timeout(opt と saddle)と opt の geometry_maxiter は 2 回まで、opt の autoz(input_invalid)は Cartesian で 1 回(vectors だけを引き継ぎ、driver の Hessian は捨てる。固定結合は spring bond にする)。frame がまだない autoz と saddle の autoz は同じ始点から Cartesian で 1 回。saddle の maxiter は継続せず、最後の frame とその歩のエネルギーを持つ Failure を返す(その frame から新しい Hessian でやり直すかは case が決める: 経路のエネルギーの上限の内側で、種の続きの段数が残るときだけ。§6.2)。scf_not_converged は DFT だけ 1 回: 同じ始点から atomic guess + `cgmin` + 通常の SCF(§7.1。駆動系のジョブも attempt の始点からで、vectors だけを捨て継続の driver の Hessian は保つ)。発散した vectors からは再開せず、WFT は再試行しない。結果を残さなかった前の attempt(止めた・殺したジョブ)は、timeout の後と同じく 1 回だけ最後の frame から継続する(`adapter.continuation(..., Failure(TIMEOUT))`。ladder の回数に数えない。読めない(途中で切れた)movecs と drv.hess は使わない)。attempt の timeout は site の `timeout_s`(既定 86,400 s)で、ハングの検出だけに使う(止まった opt・saddle は再実行で最後の frame から続く)。
- **失敗の閉じ込め**: 項目の例外を閉じ込めるのは `StageRuntime.contain` の 1 か所だけで、`thread_map` の項目にも直列のループにも同じ呼び出しを使う(minima の化学種ごとの opt・緩和・SP、sp の (method, freq) ごとのジョブ、反応ケース、explore の試行、conformers の探索)。例外は traceback 付きで log に出し、`StageState.n_errors` に数え、`Failure(error, 'error:<型>: <1 行目>')` を返す。`HFAUTO_STRICT=1`(テスト専用、`tests/conftest.py`)はここで再送出する。例外を出した反応ケースは failed の reaction artifact になり、そのケースが共有の `Registry` と species に登録したものは取り除く(後のケースが出力されない basin に入らないように)。ジョブは JobStore に残る。minima の登録(共有の Registry)は閉じ込めない。消費する型の入力がない stage は 0 artifact で done になり、上流の失敗を warning に出して後続の stage(report を含む)へ進む。
- **並行実行**: `JobRunner` のセマフォで実行中のジョブの ranks × threads の合計を site の `cores` 以下に保つ。CREST と explore の単位は `thread_map` で並列(結果は入力順)。反応ケースは直列だが、ケースの中の独立な NWChem ジョブ(SCREEN の内部 9 点と barrierless の中点 2 点の SP(開殻でも各点は端の vectors から直接始めるので並列のまま)、非対称な QRC の両側)は `CaseRuntime.map`(`thread_map`)で同時に走らせる。項目は `cores` 個ずつの波で走り、各項目は cores // (その波の項目数) の MPI rank を使い(4 コアで 9 点の SP は 4×1、4×1、1×4)、timeout は同じ倍率で延ばす。NWChem は常に mpirun に `--bind-to none` を付けて走らせる(Open MPI 5 は np ≤ 2 の mpirun をすべてコア 0 から束ねるので、rank を減らしたジョブや、ranks < cores の site で同時に走る全 rank のジョブがコアを取り合う)。rank は `ExecutionSpec` なので鍵に入らず、同じ鍵のエネルギーは直列と 1e-8 Eh 以内で一致する。minima は opt の群を一度に、次に新しい basin(既知の basin にも先の opt にもない opt)の freq と mode-follow を一度に走らせ、登録は species id の順に直列で行う(1 つずつと同じ登録簿)。sp は (method, freq) のジョブをすべて一度に走らせる。`thread_map` は 1 つの項目が例外を出すと、まだ始まっていない項目を取り消す。`SiteLock`(`<scratch_root>/.hfauto_site.lock`)は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。ロックは保持者のホストを記録し、別のホストの保持者は生きているとみなす(死んだホストのロックは手で消す)。
- **シグナル**: `run_command` は外部プログラムを新しいセッションで起動して登録する。`hfauto run` は SIGINT・SIGTERM・SIGHUP を受けると(POSIX)登録済みのプロセスグループをすべて止め、実行中の stage を `incomplete` と記録して rc 1 で終わる(Windows の Ctrl-C も同じ)。止めたジョブは JobStore に残らず、再実行すると done でない stage(failed・incomplete)を JobStore から取り直し、止めた opt・saddle はその最後の frame から続ける(acac の QRC の側の opt で、attempt_01 の最初の frame が attempt_00 の最後の frame と一致: max \|dx\| 0.0 Å)。
- **終了コード**: 0 = すべての stage が done(化学の未解決、CREST・SCC・scf_unavailable の失敗、補助の負の結果は coverage に出す)、1 = 止まった stage(incomplete か running)、error の項目、`VOLATILE` のジョブの失敗のどれか(再実行に仕事が残る。その stage は fresh でない)、2 = stage の失敗(入力の拒否を含む)。`run` と `status` は同じ規則(`cli.exit_code`)である。
- **run ディレクトリ**: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/cases/`、`<run>/run_state.json`(実行順の stage 状態とジョブの統計)、`<run>/resolved_config.yaml`(設定と code_version = git の sha、差分があれば -dirty)。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に合わせたもの(別の pipeline を同じ run に追記できる)。done で入力と設定の sha が変わらない stage は飛ばし(再開。pending・failed・incomplete は走らせる)、`--from X` は X 以降を取り直す。記録の型は未知のフィールドを拒否する(`extra="forbid"`)ので、型が変わる前の run は途中の stage からは再開できない。round 7 より前の run は新しい run dir で流し直す。それより後の run は、記録の欄が変わっていても、最初の stage からの `--from` なら JobStore を使って取り直せる(再生ゲートはこの形)。
- **検証**(`validation/`): 期待する結論は `validation/cases.yaml`(pipeline の連鎖、反応ごとの両端の断片の組成式の式・outcome・δG_eff と許容、境界の系の `twice`、置き換えた case の `superseded`)に置く。BH76 の参照は `validation/bh76/subset.yaml`(P0e の宣言ファイルそのまま。GMTKN55 の構造は `bh76/xyz/`、出所と sha は `bh76/SOURCES.json`)で、case は `bh76: <反応>.<向き>` で行を指す。比較は純関数 `compare`(結論と式)と `deviations`(参照と矛盾する outcome。合格に数えない)の 2 つで、`validation/check.py RUNS_DIR` が run dir の集まりにかける。BH76 の表は `validation/bh76/table.py`(行ごとの電子障壁を PBE0 構造と GMTKN55 構造で比べ、g・類別 MUE・SCF・費用と事前の規則による層の判定を table.md に書く)と `geometry_sp.py`(hfauto の NWChem エンジン、原子 guess と SCF の梯子、QM のロック。単点を measured.json に書く)で作る。`validation/replay.py OUT --src DIR --strict` は記録済みの run を JobStore から再生し(新しいジョブ 0 件と同じ記録)、`--src` なしなら新しい run を流す(`twice` の case は `<name>.2` にも)。OUT は `/home/user/hfauto_r10` の下で、QM のロックを取る。

## 9. 設定

4 つの層は中身が重ならないので、マージも優先順位もない。結果を変える設定は method・system・pipeline に置き(ジョブの鍵か stage の設定 sha に入る)、site には実行設定だけを置く。system と pipeline の YAML は未知のキーを実行前に拒否する。

| 層 | 置き場所 | キー |
|---|---|---|
| site | `configs/sites/` | `site`、`scratch_root`、`cores`、`engines.<registry 名>`: `version`(版数の pin)、`executables`、`execution`(`ranks`、`threads`、`memory_mb_per_rank`、`timeout_s`、`env`)、`python`(worker)、`scratch_dir`。絶対パスはこの層だけに書く(`scratch_root` と `scratch_dir` の `$VAR` と `~` は読み込みで 1 回だけ展開する) |
| method | `configs/methods/` | `id`(ファイル名と同じ)、`kind`(dft / xtb / wft)、`functional`、`wft_method`(ccsd(t))、`basis`、`dispersion`(d3zero / d3bj)、`gfn`、`electronic_temperature_K`、`grid`、`scf_energy_tol`。gfn2(低レベル)、pbe0-d3bj_def2-svpd(停留点)、m06-2x-d3_def2-tzvpd(順位のエネルギー層)、pbe0-d3bj_def2-tzvpd と wb97x-d3_def2-tzvpd(手法パネル)、ccsd-t_def2-tzvpd(opt-in) |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | `system_id`、`species`(`id`、`xyz` か `smiles`、`charge`、`multiplicity`、`role`: monomer / endpoint)、`compositions`(`id`、`components`、`multiplicity`)、`reactions`(`id`、`reactant`、`product`、`coordinate`、`torsional`) |
| pipeline | `configs/pipelines/` | `pipeline_id`、`gates`(knob)、`stages`(`id`、`stage` と下の stage 表の設定キー) |

**knob は 7 つ**だけである(検証は既定値で行っている)。

| knob | 置き場所 | 既定値 | 意味 |
|---|---|---|---|
| `noise_cm1` | `gates:` | 10 | 停留点で、これより小さい虚振動(cm⁻¹)は数値の雑音(停留性は §6.1 で別に判定する) |
| `saddle_cm1` | `gates:` | 50 | これより大きい虚振動(cm⁻¹)は鞍点級 |
| `resolution_kcal` | `gates:` | 1.0 | 経路の山と井戸を数える深さ(kcal/mol) |
| `reaction_window_kcal` | `gates:` | 40 | DFT の入口の辺の高さ(DFT//xTB)と、仮説の反応物の漸近から生成物までの ΔE(`hypotheses.select`)の上限(kcal/mol) |
| `spin_tol` | `gates:` | 0.1 | \|⟨S²⟩ − S(S+1)\| の許容幅(開殻のプロファイルの両端のスピンの結合が違うかと、最大点の ⟨S²⟩ が隣から離れているか(枝の跳び)の判定にも使う。§7.1) |
| `max_saddle_attempts` | reaction-paths の `policy:` | 2 | 1 case の失敗の印の数(答えのない探索・接続・string が 1 つずつ。子反応は 0 から。maxiter からの継続は同じ印の中) |
| `max_split_depth` | reaction-paths の `policy:` | 2 | 多段反応の分割の深さ(それより深い子は駆動せず `split_depth` で記録) |

- 停留点の method は minima(dft)と reaction-paths の 2 か所に書き、食い違えば読み込み時に拒否する。手法を変えるときは別の method ファイルを使う。
- method は pipeline ファイルの隣の `methods/<id>.yaml`、なければ作業ディレクトリの `configs/methods/<id>.yaml` を読む。
- SCREEN の有無は reaction-paths の `screen:` の有無で決まる。string の bead 数と 1 チャンクの反復数、NEB の画像数、QRC の振幅、会合のスキャンの幅と点数はモジュール定数である。
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
| `explore` | `hfauto.stages.explore:ExploreStage` | minimum, species | discovery, species | engine*, method*, max_trials, settings |
| `reaction-paths` | `hfauto.stages.reaction_paths:ReactionPathsStage` | minimum, species | reaction, calculation, minimum, species | method*, engines*, screen, policy, reaction_ids |
| `sp` | `hfauto.stages.single_point:SinglePointStage` | minimum, reaction | calculation | engine*, methods* |
| `thermo` | `hfauto.stages.thermochemistry:ThermoStage` | minimum, calculation | species_thermo, reaction_thermo | settings, energy_method, temperatures_K, standard_states |
| `report` | `hfauto.stages.report:ReportStage` | — | report | T_K, standard_state |

入れ子のキー: minima の `select`(`include`: window / all、`per_state` 3、`window_kcal` 6、`rerank_top` 8、`max_edges` 6)と `init_hessian`(`engine`、`method`)、explore の `settings`(ReaDuct の `max_scf_iterations`、`electronic_temperature_K`、`scc_retry_temperature_K`、`imag_cutoff_cm1`、`timeout_s`)、reaction-paths の `engines`(`qm`、`saddle`、`path`)・`screen`(`method`、`qm`、`path`)・`policy`(上の knob)、thermo の `settings`(`qs`、`cutoff_cm1`、`vib_scale`、`sensitivity`)。

## 10. 範囲外と既知の制限

範囲外(分析 [2026-09-30_remaining_issues_analysis.md](reviews/2026-09-30_remaining_issues_analysis.md) §0.4)。入口で止められるものは理由を付けて止め、止められないものは結果をそのまま出す。

| 項目 | 扱い |
|---|---|
| 多参照法の定量的なエネルギー(CASPT2、MRCI)、会合の速度定数(VRC-TST、圧力依存) | 範囲外。低スピン結合の対で答えるのは AP-UKS による障壁の有無と ΔG_assoc まで(§7.1)。多参照性のゲートはない |
| スピン交差(MECP) | 多重度だけが違う端点は `spin_crossing_reaction_unsupported` で止める |
| 開殻一重項 | 三重項の TDDFT/TDA による安定性の診断がないので範囲外。低スピン結合の組成の一重項は `low_spin_singlet_unsupported` で止める。単一の分子の一重項ビラジカルは RKS で計算され、不安定性は検出しない |
| トンネル | 順位に入れない。H 移動の δG_eff はその分だけ過大になりうる。1 次元 Eckart の見積もりで、298 K の κ は HCN 5.4、NH3 1.8、HONO 1.7(RT ln κ で 1.0、0.35、0.31 kcal/mol)。交差温度 T_c = hc\|ν̃‡\|/(2πk_B) ≈ 0.229 K·cm × \|ν̃‡\| より低い温度では TST そのものが成り立たない |
| 溶媒 | 気相だけ。NWChem の COSMO の freq は勾配の有限差分になる |
| 気相で束縛されない小さな多価陰イオン | 範囲外。電荷だけでは拒否しない(大きな多価陰イオンは束縛される) |
| 陽イオン | コードの制限はないが、検証セットにない |
| 遷移金属の DFT | 計算はできるが未検証。PBE0 は高スピンを過度に安定化する(Reiher et al., Theor. Chem. Acc. 107, 48 (2001)) |
| ビット単位の再現性 | 範囲外。並列の NWChem(NXTVAL の和の順序)と CREST(OS の乱数の種)は run ごとに軌跡が分かれる。再現性は結論の水準(生成物の状態の集合、各段の区分、感度の幅の中の δG_eff)で定義し、id ではなく構造の対応で比べる。NT2 の負の理由の件数は観測量ではない |
| 極小の stage の並列 | opt の群と新しい basin の freq を同時に走らせる(§8)。重い化学種 1 つが支配する組では速くならない: s19 の minima(dft) は直列の見積もりより約 14% 長く(W5、[validation.md](validation.md) §11.5)、round 7 の TMA は 105 → 206 s。core 秒は半分になる。VAL10 の時間で直列と比べて決める |
| 壁時計 | 共有ホストの壁時計は負荷に依る。比較は歩数と core 秒で行う |
| Windows のアプリケーション制御 | 品質ゲートは WSL の prod venv で判定する([environment.md](environment.md) §4) |
| 物理行数の目標 | 指標は SLOC(radon raw)と、1 つの概念を実装している場所の数 |

既知の制限。R9時点の結果報告は過去の証拠であり、現在の進捗は [roadmap.md](roadmap.md)、W5の部分結果は [validation.md](validation.md) §11.5を参照する。

- 停留性の認証の半径と幅は basin の基準(0.05 Å·√N、5e-5 Eh)で、NWChem の saddle の既定の収束(gmax 4.5e-4 Eh/bohr)より厳しい。大振幅の柔らかいモード(40 cm⁻¹ 未満)に傾きを残して収束する TS は、継続の後も認証に落ちて `not_stationary`(順位の外)になる(acac の TS: 継続の前後とも ΔE_TR 2.05e-4 Eh。[validation.md](validation.md) §11.5)。
- 状態ラベルは、イオンの E/Z と電荷分離した二重結合(イリド、HCOH の cis/trans)の E/Z を区別しない(立体の知覚は断片ごとに電荷 0。§7.1)。explore の類は立体を区別しないので、立体だけが違う経路は 1 つの類になる。
- explore の網羅は件数の予算で限られる。ibuprofen は 10,708 類で、3000 件の予算では f2b2 の 7,708 類が `not_attempted` になる(report に budget-limited と出る)。d ブロックの原子が結合した出発点は Lewis の価数を持たないので類を持たない。
- DFT の入口の高さは、未緩和の seed の SP を出発の点にすると甘くなりうる(S5 の二重項の seed の SP は CH3 + O2 の漸近より約 19 kcal/mol 高い。[validation.md](validation.md) §11.4)。案(未実装): 宣言した単量体の seed は単量体の SP の和から測る(`thermo.separated_states`)。
- ⟨S²⟩ による状態の照合(§7.3)は、同じ多重度の別の空間状態(CH3OO• の X̃ と Ã など)を区別しない。射影の guess は親の解に留まるので、親の freq が励起の SCF 解にあれば層の SP もそこに留まる([validation.md](validation.md) §11.2)。
- 零点振幅の基準でも点群の境界は遠ざかるだけで、なくならない(VAL9 の S6 の TS は上昇 5.9e-5 Eh に対し ½ħω_d 7.7e-4 Eh で Cs)。
- string(ZTS)の bead は ⟨S²⟩ を持たず、端の SCF の解を guess にしない。SCF の救済の後の通常の SCF も string の deck には入らない(§7.1)。

## 11. 次の改良

今後の設計・実装順・受入条件は [improvement-plan.md](improvement-plan.md) に集約する。現在の統合・検証状況は [roadmap.md](roadmap.md)。旧R10の未統合案と変更履歴は [history/round10.md](history/round10.md) に保存した。
