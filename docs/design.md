# hfauto の設計

本書は現在のコード(`hfauto/`)の説明である。改良で目指す設計は §11 にあり、各 wave の統合で §0〜§10 に取り込む。対象範囲と使い方は README、範囲外と既知の制限は §10、環境は [environment.md](environment.md)、実計算の記録は [validation.md](validation.md) にある。コードの docstring にある「design §x.y」は、計画時の設計書 `docs/reviews/2026-09-25_refactor_design.md`(凍結)の節番号を指す。

## 0. 原則

- 化学の判定は `hfauto/chemistry/` の純関数だけが行う。stage は入力を集めて driver か純関数を呼び、artifact を返す薄いアダプタである。
- 外部ジョブの観測事実は型付きの `Evidence` で受け渡し、合否はゲート関数(`hfauto/chemistry/gates.py`)だけが決める。
- バックエンドは能力 Protocol の実装で、subprocess を起動するのは `hfauto/execution/process.py` だけである。エンジン名の文字列は `hfauto/backends/` と `configs/` の外に書かない。
- fail-closed: 互換層、ダミーエンジン、内部フォールバックは持たない。失敗は `FailureKind` 付きで最小の単位に閉じ込め、run は続ける(§8)。
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
| U4 反応探索 | 生成物を仮定しない片端探索で、低レベルの固有な生成物と TS を出す | explore、`hfauto/chemistry/trials.py`、ReaDuct | `DiscoveryRecord` |
| U5 仮説・経路 | 極小対から仮説を作り、DFT 極小を固定端とする経路(会合は分離した単量体からの緩和スキャン)を分類して種を渡す | `hfauto/chemistry/` の hypotheses・profile、SCREEN・FIND_PATH | `ReactionRecord`、`BarrierVerdict` |
| U6 鞍点・接続 | 種を 1 次の鞍点に精密化し、別ジョブの freq で TS を確かめ、QRC で両側の極小につなぐ | REFINE_SADDLE・VALIDATE_AND_CONNECT・VALIDATE_INTERMEDIATE | `SaddleClaim`、`ConnectionClaim`、outcome |
| U7 熱化学 | 停留点の G(qRRHO、キラリティ、エネルギー層)と順位の量 δG_eff を出す | thermo、`hfauto/chemistry/thermo.py`(GoodVibes) | `SpeciesThermo`、`ReactionThermo` |
| U8 一点計算・報告 | エネルギー層を反応が読む点(`thermo.reaction_points`)だけで計算し、δG_eff で並べ、手法の幅を列で示す | sp、report、`hfauto/reporting/` | ranking.csv、method_panel.csv、report.html |
| U9 実行基盤 | stage を動かし、ジョブを再利用・継続し、失敗を閉じ込め、コアと予算を管理する | `hfauto/pipeline/`、`hfauto/execution/`、`hfauto/cli/` | manifest、JobStore、終了コード |

## 3. Evidence とゲート

`Evidence`(`hfauto/core/evidence.py`)は正常終了した外部ジョブ 1 本の観測事実で、存在すれば次が成り立つ: 正常終了と SCF 収束、opt / saddle の収束、freq なら 3N − n_external 本の振動数(単原子は 0 本)とその Hessian、原子順序と入力フレームの保持、出力から観測した `Level` が要求した `MethodSpec` と site の版数 pin に一致すること。満たさなければアダプタは `Failure(kind, reason)` を返す(maxiter で止まった saddle は最終フレームと最後の歩のエネルギー `energy_hartree` も持つ)。opt と saddle の Evidence は最終構造の DFT の勾配 `gradient`(Eh/bohr、3N、入力フレーム)も持ち、最後の勾配のブロックが最終フレームにない出力は `incomplete_output`(`gradient_not_at_final`)になる。極小の停留性はこの勾配で判定する(§6.1)。振動数と正準形式の Hessian(Eh/bohr²)は、どのエンジンでも `hfauto/chemistry/vibrations.py` の補空間射影(並進・回転を除き、主同位体の質量)で求めるので、NWChem の表示値とは低振動モードで数 cm⁻¹ ずれる(G で約 0.02 kcal/mol)。thermo は同じ Hessian を、点群の対称化した構造で射影し直して使う(§7.1)。

FailureKind は executable_missing、input_invalid、timeout、nonzero_exit、scf_not_converged、geometry_maxiter、incomplete_output、method_mismatch、gate_rejected の 9 種。

化学の閾値は `Policy` の 5 つ(§9)で、数値の許容幅 `qrc_drop` = max(1e-5 Eh, 20 × scf_tol)(同じ幾何での SCF の雑音)と反応モード性の下限 `REACTION_MODE_MIN` = 0.3(質量加重の χ の分布の空白 (0.065, 0.635) の中。[validation.md](validation.md) §11.2)は gates のモジュール定数で、YAML では変えない。

| ゲート | 判定 |
|---|---|
| `same_pes(*levels, numerics=, state=)` | program・version・method・basis・dispersion・電子温度・電荷・多重度が一致。`numerics=True` なら grid と scf_tol も(停留点・freq・QRC・熱化学)。`state=False` は電荷と多重度の違う点どうし(鎖の井戸と分離した単量体、blocker の LOT の判定、手法パネルの LOT のキー) |
| `spin_ok` | ⟨S²⟩ がないか、\|⟨S²⟩ − S(S+1)\| ≤ spin_tol。freq と、thermo が使うエネルギー層の SP にかける。不合格は blocker `spin_contaminated`。手法パネルでは不合格のエネルギーを使う値を最小・最大から外す(§7.3) |
| `imaginary_tier` | 最低の ν < −saddle_cm1 は saddle、< −noise_cm1 は soft、< 0 は noise(停留点の振動数の段。停留していない極小の候補は §6.1 で soft にする) |
| `is_minimum(freq, opt=)` | freq が opt の最終構造で同一 PES(numerics を含む。NWChem の freq は opt・saddle の収束 vectors から SCF を始めるので同じ電子状態に留まる)、\|E_freq − E_opt\| ≤ qrc_drop(外れたら `state_mismatch`)、tier が saddle でない(soft・noise は注記) |
| `is_first_order_saddle(freq, saddle=)` | 同じ連結条件で、最低モード < −noise_cm1(なければ `no_imaginary_mode`)、2 本目 < −saddle_cm1 なら `higher_order`(2 本目が −saddle_cm1〜−noise_cm1 なら注記 `soft_secondary_mode`)。モードの大きさ・重なり・端点とのエネルギー比較は見ない。最低モードが −saddle_cm1〜−noise_cm1 の停留点も TS として受理し、この case の TS かどうかは χ と QRC の接続で決める(下限を saddle_cm1 に上げると、\|ν\| < 50 cm⁻¹ の重い回転子の TS を受理できない) |
| `reaction_mode_character(symbols, mode, coords, bonds, gradient)` | 一次の鞍点が仮説の TS であることの必要条件。χ = ‖QᵀL̂‖ ≥ `REACTION_MODE_MIN`(`reaction_mode_chi`)。χ はモードと IRC が定義された質量加重の計量で測る: L̂ は Cartesian の虚モード q(M^-½L を正規化したもの)から作る単位ベクトル M^½q(重みを 2 度かけない)、Q は変わる結合の質量加重の Wilson 伸縮ベクトル M^-½∂r/∂x が張る空間の正規直交基底(SVD)。結合が変わらない仮説は宣言座標の質量加重の勾配との \|cos\|、どちらもなければかけない。結合は case のラベル付きの両端の変化 `Ctx.change()` から取り、QRC の側からは取らない(縮退反応では側の結合変化が消えて見える)。不合格は `not_reaction_mode:<χ>`(向き替えや回転子の鞍点)。限界: 無関係なモードでも χ はおよそ √(k/(3N−6))(k は変わる結合の数)なので、原子が少なく変わる結合が多い反応では χ は弱く、QRC が決める |
| `same_spin_state(ev, ref)`、`energy_spin_ok(energy, freq)` | 同じ構造の 2 つのジョブ(LOT は問わない)の ⟨S²⟩ の差 ≤ spin_tol(どちらかが ⟨S²⟩ を持たなければ判定しない。不合格は `spin_state_mismatch`)。`energy_spin_ok` は freq の `spin_ok` とこの照合の両方で、層の SP と手法パネルのエネルギーが使えるかの唯一の定義(不合格の対象は `spin_contaminated`) |
| `profile.judge` | DFT のプロファイル(`profile.Profile`: 節点の構造・エネルギー・出所・⟨S²⟩)の型と判定はそれぞれ 1 つで、SCREEN・string・会合のスキャン・VALIDATE_INTERMEDIATE の部分プロファイルが共有する。両端は DFT 極小のエネルギー(会合のスキャンでは分離した単量体の和と付加体)。順に: 3 点未満は unavailable(`too_few_points`)、⟨S²⟩ をもつプロファイルで内部の最大点の SCF の枝が跳んでいれば unavailable(`scf_branch_jump`、§7.1)、山と井戸を resolution_kcal の深さから数えて、井戸があれば intermediate、なく山があれば single、どちらもなければ barrierless。barrierless なら最も高い内部節点の両隣の区間の中点に、そのプロファイルのエネルギー関数で節点を足して 1 回だけ分類し直す(長さ 0 の区間には足さない。SP が失敗したら unavailable、`midpoint_single_point`) |
| `connection` | QRC の両側が TS と同一 PES で E_TS − qrc_drop より下の極小に割り付くこと(途中の軌跡は見ない)。組は case の判定の粒度のキー(§6.2)で比べる。期待の組なら elementary、別の登録極小の組なら reassigned、縮退反応で両側が同じ basin かつラベル付きで別構造なら degenerate(結合が変わる縮退反応では両側の結合グラフの組が {bonds(R), bonds(P)} であること。違えば `bond_change_missing`)、それ以外で両側が同じ basin なら failed(`sides_same_basin`) |
| `discovery_verdict` | explore の採否: NT2 は検証済みの TS と IRC が要る。ΔE_rxn ≤ reaction_window_kcal(評価できなければ窓の外)、低レベルの障壁 ≤ 50 kcal/mol |
| `rankable` | outcome が elementary / degenerate / reassigned で、ReactionThermo に blocker がない。barrierless_at_resolution は TST の障壁を持たないので序数の順位に入れず、ΔG_rxn とともに別の表(捕獲律速)に出す |
| `reaction_tier` | connected(ConnectionClaim)> saddle(SaddleClaim)> minima > screening |

claim を作るのは 1 か所だけである。`MinimumRecord` は minima stage と reaction-paths の `Registry`、`SaddleClaim` と `ConnectionClaim` は `ReactionCaseDriver` が作る。blocker(thermo_unavailable、mixed_level_of_theory、spin_contaminated)は thermo stage だけが作る。thermo と report は claim を再検証しない(thermo の `same_pes` は例外)。

## 4. payload

各 stage の `<run>/<stage_id>/manifest.json` にはその stage の出力だけを入れる。artifact は `status: success | failed` を持ち、payload は `kind` で判別する frozen なタグ付き共用体(`hfauto/core/records.py`、`extra="forbid"`)である。

| 型 | レコード | 出す stage | 主な中身 |
|---|---|---|---|
| species | `SpeciesRecord` | structures、conformers、minima、explore、reaction-paths | 組成キー、電荷、多重度、構造、source、状態ラベル |
| calculation | `Evidence` | minima、reaction-paths、sp | 1 ジョブの観測事実(id は `calc_` + job_key の先頭 16 桁) |
| minimum | `MinimumRecord` | minima、reaction-paths | basin、tier(screen / dft)、level_key、opt / freq の calc、members、注記 |
| discovery | `DiscoveryRecord` | explore、minima(mode_follow) | 機構、trial、product / negative / failed と理由、両端の species(`source_species`、`product_species`)、低レベル TS、ΔE‡・ΔE_rxn。DFT 段の mode-follow は鞍点の opt(`ts_calc`) |
| reaction | `ReactionRecord` | reaction-paths | 両端の極小と構造、会合なら分離した単量体の極小(`monomers`)、source、低レベル TS(`low_level_ts`、優先順)、検証済みの DFT 鞍点(`ts_calc`)、`BarrierVerdict`(出所 screen / string / scan)、`SaddleClaim`、`ConnectionClaim`、`CaseOutcome` |
| species_thermo | `SpeciesThermo` | thermo | T ごとの G / H / ZPE(失敗は None) |
| reaction_thermo | `ReactionThermo` | thermo | δG_eff と感度の幅、そのゼロの点(`reference`: complex / separated)、ΔE‡、ΔE_rxn、ΔG‡、ΔG_rxn、ΔG_assoc、分離反応物基準の ΔG‡、blockers、`energy_level`、注記 |
| report | `ReportRecord` | report | 順位の行と表のファイル |

CaseOutcome は elementary_step、degenerate_rearrangement、reassigned_step、multi_step、barrierless_at_resolution、same_basin、out_of_window、unresolved_within_budget、blocked_upstream の 9 種。id は species が `species_<id>`、minimum が `min_<species>_<level_key[:8]>`、分割した子反応が `<parent>_split<n>`。

## 5. 能力 Protocol とエンジン

| 能力 | Protocol のメソッド | registry 名 → クラス | 結果 |
|---|---|---|---|
| qm | `energy(scf_guess)`、`optimize(init_hessian, fixed_bond, scf_guess)`、`frequencies(scf_guess)` | `nwchem` → `NWChemEngine`(DFT と CCSD(T))、`xtb` → `XTBEngine` | `Evidence` |
| path | `find_path(images, initial_path)` | `nwchem_string` → `NWChemString`(ZTS)、`pysis_neb` → `PysisNEB`(xTB の CI-NEB と TSOpt) | `PathProfile` |
| saddle | `refine(seed, hessian, mode, scf_guess)` | `nwchem_saddle` → `NWChemSaddle` | `Evidence` |
| conformers | `search(settings)` | `crest` → `CRESTEngine` | `ConformerEnsemble` |
| discovery | `explore(trial, settings)` | `readuct` → `ReaDuctEngine`(NT2) | `DiscoveryResult` |

- エンジンは `cls(jobs=JobRunner, site=EngineSite)` で作る。`engines.create` は Protocol を満たさなければ `TypeError`、表にない名前は `KeyError`。registry を import してよいのは pipeline だけである。
- すべてのアダプタに電荷と多重度を渡す(NWChem の `charge` / `mult`、xTB と CREST の `--chrg` / `--uhf`、pysisyphus と ReaDuct の計算器設定)。
- 初期 Hessian は freq の `Evidence` で渡す。saddle は種と同じ構造、optimize は同じ原子順序で各原子が 0.5 Å 以内の構造のものだけを受け付け、外れれば `INPUT_INVALID(hessian_geometry_mismatch)`。xTB は初期 Hessian を使わない。
- optimize の `fixed_bond` (i, j, r) は原子対の距離を r Å に固定する(会合のスキャンの点)。開始構造はその距離を持っていなければならない(`fixed_bond_mismatch`)。NWChem は `zcoord` の `bond i j r constant`、Cartesian の継続(autoz の後)では `constraints` の `spring bond`(k = 20 Eh/bohr²。距離は力 F の向きに F/2k だけずれ、CH3 + O2 のスキャンの最大の力で 7e-4 Å。束縛の項は全エネルギーに入らない)。xTB は `constraint_unsupported` で拒否する。`scf_guess`(qm の 3 つと saddle)は親のジョブ(同じ原子と電子状態。LOT は問わない)の収束した vectors から DFT の SCF を始める。基底が違えば親の基底を deck に `parent` と宣言して射影する(`vectors input project parent guess.movecs output job.movecs`。出力を名指ししないと NWChem は vectors を guess.movecs に上書きし、そのジョブから派生するジョブが guess を失う)。親は、freq なら opt・saddle、層とパネルの SP なら対象の freq、スキャンの点なら前の点、開殻のプロファイルの SP なら端の DFT 極小の opt、高次の鞍点を押した探索ならその TS の freq である(§7.1)。CCSD(T) の参照は atomic guess から始め、xTB は受け取って無視する。`fixed_bond` と `scf_guess` は渡したときだけジョブ鍵に入る(scf_guess は DFT の SCF を始めるときだけ)。
- optimize は渡された Hessian を、ν < −saddle_cm1(既定 50 cm⁻¹)のモードが 2 本以上(高次の鞍点)ならそのまま書き、それ以外はすべて固有値を max(\|λ\|, 1e-3) にした正定値のモデル H₊(`vibrations.shape_hessian(h, x)`、ジョブ鍵に `hessian_model: positive`)にして書く。規則はこの 1 行だけで、呼び出し側は選ばない。モデルになるのは、一次の鞍点(QRC と mode-follow の側)、soft の点、低レベルの錯体と R6 の seed の xTB Hessian である。NWChem は負の固有値に沿った最小化の歩幅を trust の 0.03〜0.3 倍に切り、BFGS 更新は符号を保つので、負のまま渡すと下りが遅い(acac の QRC 34/27 歩 → 15/14 歩)。xTB の停留点でない seed の負の固有値も \|λ\| になる。高次の鞍点からの側は別の鞍点方向へ出られる必要があるので、そのまま書く(DME C2v はそのままで 29 歩、モデルでは 207 歩で未収束)。
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
3. saddle 級なら ν < −saddle_cm1 の虚モードに沿って `modes.off_saddle` で変位し(各モードは QRC と同じ `modes.qrc_step`: 目標 max(3 × qrc_drop, 3e-4 Eh)、0.05〜0.4 Å。和は 0.4 Å で切る)、変位元の freq を初期 Hessian にして opt → freq(最大 `mode_follow` サイクル。履歴 `follow<n>:<判定>`、側は `opt:follow<n>`)。虚モード 1 本なら ± 両側、2 本以上なら全モードの和の方向へ片側だけ押す。判定は、両側がラベル付きで別の構造なら `ts_candidate`(縮退転位を含む。両側の結果 opt・freq・注記を返す)、両側が同じ構造なら `replace`、片側だけ動けば `one_side`(その側を採る)。minima stage は `ts_candidate` の両側をジョブなしで同一性だけで登録し(鏡像の側は同じ basin)、側 1 を `source_species`、側 2 を `product_species` とする `mode_follow` の discovery にする(DFT 段では鞍点の opt が `ts_calc`)。元の化学種は入力構造に近い側の basin に注記 `endpoint_was_saddle` 付きで加わる。
4. 停留性: saddle 級でない点では、opt の最終構造の DFT 勾配 g と freq の Hessian から二次モデルの降下量 ΔE_N = Σ g_i²/(2\|λ_i\|)(freq と同じ質量加重の内部モード、\|λ\| に床なし。`vibrations.stationarity_gap`)を求め、5e-5 Eh(`identity.BASIN_DE_HARTREE`、2 つの極小を 1 つの basin とみなすエネルギーの幅)を超えれば、振動数によらず soft(注記 `soft_imaginary_mode`)にする。平らな PES の opt は、振動数に何も出ない肩で止まりうる(W3 の S6 の CH3OH···H: −6.5i で noise、ΔE_N 2.7e-4 Eh)。したがって noise_cm1 は停留点での数値の雑音の幅である。勾配のない Evidence(勾配を記録する前の JobStore の結果)は判定せず、振動数だけで決める。鞍点では、2 本目の虚モード(`higher_order`)だけを同じ基準の停留点で判定する(§7.3 の検証)。一次の鞍点の受理には停留性をかけない(既知の制限: S6 の引き抜きの TS は、ほぼ直線の C–H···O をまたぐ autoz のねじれで内部座標の勾配が縮み、Cartesian の勾配 約 1e-3 Eh/bohr、ΔE_N 1e-4〜2e-4 Eh のまま収束と判定されうる。[validation.md](validation.md) §10)。
5. soft なら 1 回だけ押して緩和する。停留していなければ Newton の歩 −H₊⁺g(`vibrations.newton_step`、H₊ は §5 の正定値のモデルなので虚モードの方向にも下る)、停留していれば ν < −noise_cm1 の虚モードの和の方向で、どちらも最大の原子変位を 0.4 Å で切る。緩和した点が soft でなくなれば `soft:resolved` でその点を採り、なお soft なら元の点を採る(`soft:persisted`、注記 `soft_imaginary_mode`)。結果はどちらも `minimum` である。noise はそのまま `minimum`。

`Registry` は組成 × level_key ごとの basin の集合で、同一性の基準は `identity.assign` の 1 つだけである(\|ΔE\| ≤ 5e-5 Eh の候補のうち、回転と鏡映を許す置換不変 RMSD ≤ 0.05 Å で、次点と十分離れていること)。鏡像は同じ basin である(m は thermo が点群から決める。§7.1)。reaction-paths では `Registry` が DFT 極小の唯一のストアで、`Registry.minima` で minimum_id から記録と代表の最適化構造(FileRef)を引く。basin の構造をメンバーの原子順に写すのは `identity.member_coords` だけである(代表自身なら入力の座標をビット単位でそのまま返す)。ラベル付きの比較(`identity.same_as_labelled`: \|ΔE\| ≤ 5e-5 Eh かつ `mapped_rmsd` ≤ 0.05 Å)は真の回転だけなので、NH3 の反転のような縮退転位は恒等と区別される。minima stage はジョブを species id の順に直列で relax し(relaxation の seed は最後)、直後に登録するので、同じ stage の中でも登録済みの basin に落ちたジョブは freq を省く。開始構造が先のジョブの厳密な置換・鏡像の像(`identity.is_image`: `carry` の RMSD ≤ `IMAGE_A` 0.005 Å。QRC の − 側の省略と同じ判定。対称な像は 0.0009 Å 以下、それ以外の組は 0.025 Å 以上で、basin の 0.05 Å とは別の基準)であるジョブは計算しない。PES は同種核の置換と反転で不変なので、代表がそのまま極小か既知の basin に落ちたときだけその basin に入る(`image_of:<sid>`)。代表が鞍点・soft・失敗なら自分で緩和する。mode-follow の経過は `<run>/<stage_id>/diagnostics.json` に書く(記録用で、読むコードはない)。

### 6.2 ReactionCaseDriver(`hfauto/drivers/reaction_case/`)

`decide(case, state, rules)` は `hfauto/drivers/reaction_case/state.py` の 13 行を上から評価する純関数で、driver はその action を実行して `CaseState` を積み上げる。ジョブはすべて JobStore を通るので、再実行すると同じ判断を数秒でたどり直す。判断は `<run>/<stage_id>/cases/<reaction_id>/log.jsonl` に書く(`hfauto case` で表示する)。

| # | 条件 | 決定(reason) |
|---:|---|---|
| 1 | 端点に reaction-paths の DFT 極小がない(宣言した端点が鞍点に落ちたなど) | BLOCKED(`endpoint_without_dft_minimum`)。両端が 1 つの PES にあることは、設定の読み込みで `PipelineConfig._check_stages` が保証する(minima(dft) と reaction-paths の method が 1 つで、stage は DFT 極小だけを読む) |
| 2 | 両端が同じ basin で、縮退反応でも会合(反応物は分離した単量体)でもない | SAME_BASIN(`same_basin`) |
| 3 | ΔE_rxn > reaction_window_kcal | OUT_OF_WINDOW(`out_of_window`) |
| 4 | 接続が elementary / degenerate / reassigned | 完了(`connection:<label>`) |
| 5 | 中間体が両端と別(判定の粒度で) | MULTI_STEP(`intermediate_distinct`)。子反応 R→I、I→P に分ける(下の「分割の子」) |
| 6 | 最新の DFT プロファイル(SCREEN、会合のスキャンか string)が barrierless | BARRIERLESS(`screen:barrierless` / `scan:barrierless` / `string:barrierless`) |
| 7 | 接続を判定したが完了していない | 両側が同じ basin の TS(2 つの振幅の後)と、両側が同じキーの別 basin の TS(`same_state`)は、この case にとって別の過程の鞍点なので、saddle の試行が残れば探索を続ける(行 9〜12。次の saddle 探索がその TS の主張を捨てる)。柔らかい TS(\|ν\| < saddle_cm1)も同じ。それ以外は UNRESOLVED(`connection_failed`) |
| 8 | 収束した saddle が未検証(`ts_calc` の case は最初から) | VALIDATE_AND_CONNECT(`saddle_converged`) |
| 9 | SCREEN が有効(会合は常に)で、低レベル経路(会合ではスキャン)が未実行、かつ待っている種がない(近道の種がすべて失敗したら string の前に 1 回だけ NEB) | SCREEN(`screen`) |
| 10 | 最新の DFT プロファイルが intermediate で中間体が未判定 | VALIDATE_INTERMEDIATE(`path_intermediate`、最も低い井戸) |
| 11 | 未使用の種があり、saddle の試行 < `max_saddle_attempts` | REFINE_SADDLE(`seed:<source>`) |
| 12 | 会合でなく、saddle の試行が残り、走らせた string の数 ≤ saddle の試行数で、string が未実行か、最新のプロファイルが unavailable でない string | FIND_PATH(`dft_path`) |
| 13 | 上のどれでもない | UNRESOLVED(`attempts_exhausted`) |

**予算**: 件数だけで、壁時計の締め切りは持たない。行 4〜6 は手元の証拠だけで決まる完了で、行 7 以降は新しい計算を始めるか打ち切る。saddle の試行は case ごと(子反応はそれぞれ 0 から。兄弟の子が試行を使い切っても後の子は自分の試行を持つ)で、種 1 つが 1 試行である(高次の鞍点を押した種も数える)。maxiter で止まった探索の再開は、同じ試行の中で行う。string は、それまでの string の種をすべて試した後(string の数 ≤ saddle の試行数)だけ次を走らせるので、種を出さなかった string の後は打ち切る(行 13)。分割の深さは `max_split_depth` で切る。1 つの attempt の長さは site の `timeout_s` だけが決め、run 全体を止めるのは外からのシグナルである(§8)。

**続きの段数**: 種は、新しい種(SCREEN・FIND_PATH の種と低レベル TS)から数えた続きの段数(`Seed.depth`)を持つ。続きは、maxiter で止まった探索の最後のフレームからの再開と、高次の鞍点を押した種(`higher_order_retry`)の 2 つで、どちらも段数 + 1、2 段(`actions.MAX_DEPTH`)までである。押した種の探索が止まっても再開でき(W3 の S6 の split2_split2 の経路)、押す・再開の交互の連鎖は 2 段で終わる。

**再開の上限**: maxiter で止まった探索の最後の歩のエネルギー(`Failure.energy_hartree`)が、最新の DFT プロファイル(SCREEN、string か会合のスキャン。低スピン結合のスキャンは AP のエネルギー)の最大 + resolution_kcal を超えるなら、再開しない(注記 `saddle:above_path_bound:<screen|string|scan>`)。両端を結ぶ連続した DFT の経路の最大は、探している鞍点のエネルギーの上限なので、それより上へ登った探索は、その鞍点を通り過ぎている(VAL7 の S5 は、プロファイルの最大より 27〜32 kcal/mol 上のフレームから再開していた。[validation.md](validation.md) §8)。エネルギーのない Failure(エネルギーを記録する前の JobStore の結果)とプロファイルのない case には上限をかけない。

**分割の子**: 行 5 の子反応は、親の直後に、それぞれ自分の試行の予算で駆動する。次の子は駆動せず、ジョブも log もない記録にする(`classification.undriven`)。
- 同じ case キー(`hypotheses.pair_key`: 両端の状態が違えば 2 つの状態、同じなら 2 つの極小)の case が駆動済みか待ち行列にある子: その case の駆動を待ち、結論が出ていれば outcome と claim を写す(`same_as:<reaction_id>`)。同じ問いを 2 度探索しない。その case が unresolved で終わったら、結論のなさは共有せず、子を自分の予算で駆動する(待った子は、待ち行列が空になってから)。
- 親の深さが `max_split_depth` に達している子(仮説が深さ 0): unresolved_within_budget(`split_depth`)。

**原子の対応**: 反応の両端の添字は TS を挟んで連続していなければならない。basin は置換について不変な同値類で、代表の species の添字は任意なので、派生する端点には代表を使わず、実際に到達した添字付きの構造を使う。
- 中間体(QRC の reassigned と VALIDATE_INTERMEDIATE): 到達した構造(QRC の側か、緩和した井戸の opt の最終構造)をケースの原子順のまま子反応の端点にする。case が新しい basin を登録したときは、その basin の代表になった case 自身の species(opt・freq を持つ)を使う。既知の basin に落ちたときだけ、その構造をジョブなしの `SpeciesRecord`(`spc_<reaction_id>_intermediate`、source `intermediate`)にして basin のメンバーに加える(`Registry.join`)。この処理は `drivers/reaction_case/connection.py` にある。
- 発見: 仮説の両端は discovery 自身の `source_species` と `product_species` で、どちらかがなければ仮説にしない。source は、nt2 では trial を始めた screen の代表、relaxation では seed 自身の species(`hypotheses.seed_species_id`)、mode_follow では側 1 である。
- case の端点は、basin の最適化構造を端点の species の添字と掌性に並べたもの(`identity.basin_coords`)である。species の構造が basin と同じ結合グラフ(状態ラベル)を持つときは、原子を同じ WL クラスの原子にだけ対応させるので、species の添字付きの結合が保たれる。basin から遠い構造(DFT の basin に入った xTB の生成物など)でも置換を誤らない。結合グラフが違う(basin への緩和で状態が変わった)ときは元素だけで対応させる。そのときの結合変化は化学的な結果で、添字の誤りではない。したがって R → I → P の添字は連続し、分割の子の結合変化に添字の置換は入らない([validation.md](validation.md) §8)。

**判定の粒度**: case は、両端が異なる粒度で判定する(`drivers/reaction_case/connection.py` の `_key`)。
- 両端の化学状態 (composition_id, state_label) が異なる case(結合が変わる case)は状態で判定する。thermo は状態の G を最小の G とする(§7.2、Curtin–Hammett)ので、それと同じ粒度である。端点と同じ状態の井戸は、basin が違っても中間体にしない。
- 同じ状態の case(宣言したねじれ、配座変化、縮退転位)は basin で判定する。ただし端点と同じ状態の別の basin で、その端との \|ΔE\| < resolution_kcal、かつそれを運んだ経路(プロファイルの井戸なら端までのプロファイル、QRC の側なら TS)に resolution_kcal 以上の山がないものは、分解能では同じ端点とみなす。
- 実際に結ばれた極小は `ConnectionClaim.minima` に残す。elementary の `ReactionRecord.minima` は仮説の両端のままで書き換えない(重複除去と reaction_id を保つ)。
- 既知の制限: 状態ラベル(WL ハッシュ)は立体を区別しない。ジアステレオマー(SN2 の反転体と保持体など)は同じ状態になり、同じ端点とみなされうる(分析 G8-P6、未対応。今の検証系に不斉中心はない)。

action:
- **SCREEN**: 仮説の低レベル TS(`low_level_ts`。explore か mode-follow の、同じキーの相異なる TS すべて、§7.1 の仮説)を順に xTB freq で確かめ、DFT SP が両端の DFT 極小のどちらよりも resolution_kcal 以上高いものを、その順に種にする(`discovery_ts`。近道は 1 回だけで(`CaseState.shortcut_done`)、3 点はプロファイルではないので判定も ⟨S²⟩ の検査も記録しない。DFT SP の順に並べ替えるのは W6-1)。それが尽きたら DFT 極小の間の IDPP 11 点(自前の実装。pysisyphus の `interpolate.IDPP` は平面の HONO の経路を面内の H–O–N 177.8° に通し、作業ディレクトリにログを書くので採らなかった。[validation.md](validation.md) §11.3)を初期経路に `pysis_neb`(xTB の CI-NEB、端は固定、未収束でも経路として使い、CI から TSOpt)で緩和し、内部 9 点の DFT SP と両端の DFT 極小のエネルギーを `profile.judge` で分ける。single の種は NEB の TS(`screen_ts`)か、山を放物線補間した構造(`screen_hei`)。NEB が失敗したら IDPP のまま分ける。経路の両端は常に DFT 極小なので、低レベルの PES に端点の極小は要らない。barrierless は最も高い内部点の両隣の区間の中点 2 点の DFT SP を加えて分類し直してから受け入れる(`profile.judge`。FIND_PATH・会合のスキャン・VALIDATE_INTERMEDIATE の部分も同じ)。根拠(分析 G3-P7): 単峰のプロファイルなら最大は最も高い節点の両隣の区間にあり、節点が見落とす高さは、山を高さ V・底の幅 w の放物線で近似すると、節点の間隔 h に対して (h/w)²V 以下である。中点でその区間の間隔が半分になり、上限は 1/4 になる(実測の節点の間隔は 0.04〜0.30 Å)。多峰のプロファイルで節点の間に隠れた別の山は保証しない。開殻の case では、これらの SP を端の SCF の解から始め、枝の跳んだプロファイルは unavailable にする(§7.1 の「開殻のプロファイルの SCF の連続性」)。会合(`ReactionRecord.monomers` を持つ case)では、SCREEN は低レベルの経路の代わりに、分離した単量体から付加体までの緩和スキャンを走らせる(低レベルのエンジンは要らない。§7.1 の「会合と低スピン結合の組成」)。
- **REFINE_SADDLE**: 反応方向は 1 つの規則で決める(`Ctx.direction`)。
  - TS の種(低レベル TS の `screen_ts` と `discovery_ts`、高次 saddle を押した `higher_order_retry`)は、自分の虚モード(`Seed.mode`)。
  - それ以外の種(`screen_hei`、`path_hei`、再開した探索)は、ラベル付きの両端の結合変化(`Ctx.change()`。ρ・χ・高次の鞍点の反応モードの選択・縮退の QRC の判定が共有する 1 つの変化)から作る ρ = ∇(Σ_切れる r − Σ_できる r)を種の構造で評価したもの。結合が変わらなければ宣言座標の勾配、なければ最も変わる二面角の勾配、それもなければ種に重ねた両端の差(NH3 の反転など)。
  - 経路の接線は使わない。接線には回転子が混ざる(VAL7 の S6 の screen_hei で ρ との cos 0.69)。
  - 初期 Hessian は、種が 0.5 Å 以内の検証済み TS freq を持てばそれ、なければ種の xTB freq(あれば常に)、なければ種の DFT freq。
  - adapter はそれを H₀ = P·H₊·P − κd̂d̂ᵀ にして渡す(`vibrations.shape_hessian`)。H₊ は剛体運動を除いて固有値を max(\|λ\|, 1e-3) にしたもの、d̂ は方向の内部運動の成分を正規化したもの、P = I − d̂d̂ᵀ、κ = max(d̂ᵀH₊d̂, 0.05 Eh/bohr²)。負の曲率は d̂ の 1 本だけで、ほぼ縮退した固有ベクトルのどれを選ぶかに依らない。低レベルの Hessian の負のモードは \|λ\| として正の部分に入るだけなので、xTB Hessian の採否の検査はない。
  - saddle のジョブ鍵はモデル名(`hessian_model: negative_along_mode`)を含むので、以前のモデルの saddle は再利用しない。注記は `saddle_hessian:<ts_freq|xtb|dft>:<mode|rho|coordinate|chord>`。
  - saddle が maxiter で止まったら(柔らかい反応モードや押した種では、2 本目の負の固有値が十数歩残って停滞しやすい)、同じ試行の中で最終フレームから新しい Hessian で 1 回だけやり直す(上の続きの段数と再開の上限の範囲で)。
- **VALIDATE_AND_CONNECT**: TS の検証と QRC を 1 つの action で行う(`connection.validate_and_connect`)。
  - 検証: 別ジョブの DFT freq → `is_first_order_saddle` → `reaction_mode_character`。χ で拒否した鞍点(`ts_rejected:not_reaction_mode:<χ>`)と虚モードのない点(`ts_rejected:no_imaginary_mode`)は、試行に数えたまま QRC にかけず、探索を続ける。仮説の正味の結合変化を動かさない第 1 段の TS(後の段で消える一時的な結合だけを動かすもの)もここで拒否され、その中間体は行 10 の `path_intermediate`(プロファイルの井戸)で拾う。中間体を出すのは、プロファイルの井戸(行 10)と QRC の側(行 5)だけである。`ts_calc`(mode-follow の鞍点か親が検証した TS)はまずそれを検証する(freq は JobStore の再利用)。ts_calc は一次の DFT 鞍点か親が検証した TS なので、拒否されるのはこの case の両端に対する χ か freq の失敗のときだけで、そのときはほかの拒否と同じく探索の行(SCREEN など)に進む(専用の分岐はない)。higher_order なら、段数が 2 未満のとき、case の結合変化(なければ反応方向)に対する χ が最大の負モードを反応モードとし(TS の判定と同じ量)、それ以外の ν < −saddle_cm1 のモードの和の方向に片側に押した構造(`modes.off_saddle`。二次の鞍点では従来と同じ 1 本)を種(`higher_order_retry`、段数 + 1)にする。検証済みの freq をその Hessian と SCF の guess、反応モードをその `mode` にする。2 本目の虚モードは停留点でだけ数える: saddle の勾配と TS freq の Hessian で ΔE_N > 5e-5 Eh(§6.1 の極小と同じ基準)なら、符号付きの Newton 歩(二次モデルの停留点へ、\|λ\| の床 1e-3、最大 0.4 Å。`minimum.newton_push`)の点で freq を 1 本とり、そこで −saddle_cm1 より下が 1 本以下なら、押す代わりにその点を種にする(注記 `higher_order:not_stationary`。S6 の 790468f505 の ν2 −54.8i と −77.1i は、この 1 歩で実になった)。残れば従来どおり押す。勾配のない saddle(勾配を記録する前の JobStore の結果)も押す。
  - QRC: 受理した TS の虚モードに沿って ± に変位し(振幅は `modes.qrc_step`: エネルギー目標 max(3 × qrc_drop, 3e-4 Eh) と ν・質量から 0.05〜0.4 Å)、TS の freq Hessian(§5 の正定値のモデル)で opt し、`Registry` に割り付け(未知なら収束した側の opt から `relax_to_minimum` で新しい basin)、`connection` で判定する。対称な TS では − の開始構造が + の厳密な像(`IMAGE_A` 以内)なので + 側だけを opt し、− 側はその像を運んだ構造で割り付ける(注記 `qrc<n>:minus_is_image`、`side_calcs` は (plus, plus)、ゲートは同じ)。それ以外は両側を同時に走らせる(§8)。両側のキーが端点のキーの組なら elementary、片側だけが端点のキーで他方が別のキーの DFT 極小(freq で確かめた極小)なら、その側の構造を中間体(上の原子の対応)にして行 5 で 2 つの子反応に分割し(注記 `qrc<n>:end<i>_to_new_basin`)、その 2 つを結ぶ子に TS を渡す(`ts_calc`。子は検証と QRC を JobStore から再生する)。結合変化の case の中間体は、端点と別の化学状態に限る。どちらも端点のキーでない TS は reassigned のまま。両側が同じキーなら、同じ basin は `same_basin`(変位が basin を出るには小さすぎたかもしれないので、振幅を 2 倍(0.4 Å で切る)にしてもう 1 回だけ QRC を走らせる)、別 basin は `same_state`。結合の変わる縮退反応は、QRC の両側の間の結合変化(帯付き)が case の変化かその逆と等しいときだけ degenerate にする(メチルの回転の TS は移る H を動かさない)。残った判定は行 7 が扱う。多段の親は順位を持たず(outcome `multi_step`)、子反応がそれぞれ順位に載る。
- **FIND_PATH**: ZTS を 1 チャンク(9 beads、maxiter 20)だけ走らせ、収束を問わず `profile.judge` で分ける。初期経路は最新の DFT プロファイルの経路、なければ IDPP。端点は DFT 極小にし、画像を隣へ逐次整列する。種は single のときだけ山の放物線補間(`path_hei`)。その string の種を試した後も試行が残れば、行 12 がその経路から次のチャンクを続ける。
- **VALIDATE_INTERMEDIATE**: 最新のプロファイルの最も低い井戸を `relax_to_minimum` にかける。端点と別のキー(上の分解能の規則で端点とみなす井戸を除く)なら、緩和した構造を中間体にして行 5 で分割する。緩和が失敗したら結果とせず(端点扱いも barrierless もしない)、最も高い山を種にする。端点に落ちたら、もう一方の端からその井戸までの部分プロファイルを最新のエネルギー関数で `profile.judge` にかけ、barrierless なら case を閉じ、そうでなければその最も高い山を種にする(log に `<中間体>:end<k>:<判定>`)。

## 7. 化学プロトコル

### 7.1 stage ごとの既定値

| 段 | 既定値と規則 |
|---|---|
| 共通 | DFT は PBE0-D3BJ/def2-SVPD、grid fine、`convergence energy` 1e-7(opt・freq・SP で同じ)。反復上限は DFT `iterations 100`、WFT の SCF `maxiter 100`。DFT の SCF は親の解から始める(§5 の `scf_guess`。基底が違えば射影)。SCF の段は上限付きの 3 段だけ(P0a): (1) 親の(なければ atomic の)guess + DIIS、失敗したら 1 回だけ (2) atomic guess から `cgmin`(`iterations 40`。NWChem が自分で 50 に上げる)、(3) `cgmin` を外した通常の SCF をその vectors から回して ⟨S²⟩ とエネルギーを読む(cgmin は ⟨S²⟩ を出さない)。成否は (3) だけで決め、cgmin の収束では決めない: (3) が cgmin のエネルギーから qrc_drop より離れれば `scf_unavailable:scf_noise`、sp・freq(親の構造にいるジョブ)で親との `same_spin_state` に落ちれば `scf_unavailable:spin_state`、(2)(3) が収束しなければ `scf_unavailable:scf`(いずれも `scf_not_converged`)。駆動系のジョブ(QRC の側、スキャンの点、継続の saddle)は ⟨S²⟩ が正当に変わる(CH3 + O2 で 1.71 → 0.75)ので scf_noise だけを見る。発散した attempt の vectors からは再開しない。WFT は atomic guess から始め、再試行しない。smear・fon は PES を変えるので使わない。geometry は `units angstrom nocenter noautosym`。autoz の致命的な失敗(`insufficient internal variables` など)は Cartesian で続ける(§8)。Cartesian に切り替えたという無害な注記(`AUTOZ failed to generate good internal coordinates`)は失敗にしない。Z>36 の元素には def2 系の基底のときだけ `<元素> library def2-ecp` を元素ごとに書く |
| structures | 化学種は xyz と SMILES のちょうど一方。SMILES は RDKit ETKDG の 1 配座(同位体と `'.'` は不可)。宣言反応の端点は xyz(SMILES では原子の対応と配座が決まらない)で、両端の電荷と原子順序が一致すること。多重度だけが違う両端はスピン交差で、`spin_crossing_reaction_unsupported`(MECP は探さない。それぞれの面は別の反応として宣言すれば評価できる)。化学種の多重度は必ず宣言する(省略は system の読み込みで拒否。SMILES のラジカル電子数も xyz も多重度を決めない)。元素と多重度のパリティはここで 1 回だけ検査する(範囲外の元素は `unsupported_element`、電子数とのパリティ違いは `INPUT_INVALID`)。組成は conformers の入口で決める: 電荷は成分の和、多重度は宣言値かスピン結合で 1 つに決まる値(決まらなければ `declare_multiplicity`)。低スピン結合の組成は受け付け、その一重項だけを `low_spin_singlet_unsupported` で止める(下の「会合と低スピン結合の組成」) |
| conformers | 単量体: `crest --gfn2 --quick -T <n> --ewin 6 --chrg --uhf`(重原子 3 個以下で回転可能結合のない分子は省く)。トポロジー変化で止まったら停止構造を残し、そこから 1 回だけ再実行する。組成: ゲストを乱数(種は固定)の向きと方向の剛体として vdW 接触 + 0.5 Å に置く seed を `seeds_per_composition` 個作り、先頭の seed から `--nci --quick --notopo <全原子> --noopt` で探索する(状態は hfauto の状態ラベルが決め、CREST はサンプリングだけを行う。`--noopt` は CREST 3.0.2 の初期トポロジー検査が `--notopo` を見ないため)。CREST の候補に入力(単量体)か seed(組成)の状態ラベルを持つものが 1 つもなければ(CREST の失敗を含む)、それらの構造もそのまま出し、理由は項目の Failure(CREST の Failure か `GATE_REJECTED input_state_lost:<ラベル>`)に書く。GFN2 の screen で最適化を通るのは別の結合状態に崩れる seed だけで(S5 は 6 個中 1 個が CH3O2、S6 は 2 個が CH3 + H2O)、接触の seed 1 つではどちらも失敗する(W3-3 のプローブ)。CREST 3.0.2 は CH3·O2(rc −11、SIGSEGV)と OH·CH4(MTD が収束せず rc 1)で失敗した。CH3·O2 は低スピン結合の対で、スピン分極のない GFN2 では対の二重項の SCC が収束しない(四重項なら `--nci` で rc 0)。SIGSEGV がこの SCC の不安定さから来るという因果は推論である(出力に SCC 未収束の記録はない)。OH·CH4 は二重項と一重項の対なので、この機構には当たらない。開殻そのものは原因ではない(六重項の FeCl3·CH4 は rc 0)。選抜は状態ラベルごとに GFN2 エネルギーの低い `keep_per_state` 個。CREST は attempt ディレクトリで `--scratch` なしに実行する |
| minima(screen) | xTB `--opt vtight` → `--hess`、mode-follow 最大 2 サイクル。状態ラベルが変わった seed は落ちた basin の members に記録する(失われた状態は explore の `relaxation` が DFT に 1 回問う) |
| explore | 出発点は組成 × 状態ラベルごとの screen 最低 `sources_per_state` 個。trial は元素に依らない列挙器 1 つで作る: 形成の候補対は結合しておらず r ≤ Σr_vdW でグラフ上 3 結合以上離れた対。T1 移動・置換(形成 1 + 切断 1)、T2 リレー(H 移動 2 つの連鎖)、T3 形成だけ、T4 切断(開殻か電荷系だけ)を順に巡回し、出発点あたり `max_trials_per_source` 件。変化後の結合数が元素の最大配位数を超える drive は作らず、原子クラスと形成距離(0.1 Å)が同じ drive は 1 つの類、1 件の trial にする。類の代表は値(r/Σr_cov。移動は次に a–b–c 角、リレーは遠い接触→近い接触、切断は最も伸びた結合)が最小の drive で、同値は正準の順位(原子の WL クラスとクラスごとに並べた距離)で選ぶ。類は代表の値と類の記述の順に並べ、`max_trials_per_source` の境界で同値に分かれる類はまとめて落とす。trial_id は類の記述と source から作るので、trial は原子の番号付けに依らない([validation.md](validation.md) §8)。各 trial で NT2 を 1 回(AFIR は生成物を残さなかったので削除した。[validation.md](validation.md) §8)。IRC のどちらの端も source の添字付き結合を持たないときは、source と同じ状態ラベルで basin RMSD ≤ 0.05 Å の端を source の添字置換像とみなし、他端と TS を source の添字・座標系へ運ぶ。単位は並列に走らせ、入力順に記録する。採否は `discovery_verdict`。残した生成物は、同じ組成の screen 極小と既出の生成物のうち添字付きの結合が等しいものにだけ、置換不変 RMSD で 1 回だけ同定する。添字の違う写し(縮退の生成物を含む)は別の species にする。どの発見も出発点の species(`source_species`)を記録する。screen で失われた seed の状態は、組成 × 状態ラベルごとに 1 件の `relaxation` 生成物(結合変化をもつ seed のうち低レベルエネルギー最小(エネルギーのない seed は後ろで species id 順)、未緩和)にする。失われた状態とは、同じ組成の screen 極小にそのラベルがなく、seed と落ちた basin の構造(`identity.basin_coords` で seed の添字に並べたもの)の間に結合変化(下の「結合と状態」の帯で判定する)が 1 つ以上あるものである。minima はそれを `spc_<discovery_id>` として最後に未緩和の構造から 1 回 opt し(2 断片以上なら seed の xTB Hessian を初期 Hessian にする。seed は xTB の停留点ではなく負の固有値を持つが、−saddle_cm1 より下が 2 本以上でなければ §5 の規則で正定値のモデルになる)、別の状態に落ちたら diagnostics に `collapsed_at_dft_from_seed` と書く(barrierless とは呼ばない。障壁なしと言えるのは DFT のプロファイル(行 6)だけである)。状態を保てば seed の basin → 崩壊先の basin を仮説にする(結合変化・窓・`max_per_composition` は他と同じ、TS なし)。ReaDuct の `spin_mode` は `restricted_open_shell` |
| minima(dft) | 選択 `window`(`chemistry/selection.py`。stage は開始構造とジョブだけを持つ): 反応する組成(宣言反応の端点、発見の出発点と生成物、それらの錯体の宣言した単量体 = `thermo.declared_monomers`)だけを、組成 × 状態ラベルごとに `window_kcal` 以内の `per_state` 構造で DFT にかける(発見の出発点・生成物と relaxation の seed は必ず含む)。ほかの組成は diagnostics に `not_reacting`。その代償に、反応のない組成は DFT 極小も DFT の mode-follow の発見も持たない(検証セットでそうした発見は 0 件)。`rerank_sp` は候補が `per_state` を超える組(組成 × 状態)だけ、上位 `rerank_top` 構造を DFT SP で並べ直してから同じ窓で切る。`all` は入力の化学種全部。opt は初期 Hessian なしなら `trust 0.1`、ありなら `trust 0.3`、maxiter 100。2 断片以上は `init_hessian` の xTB Hessian を使う(relaxation の seed も同じ。§5 の規則で書く) |
| 仮説 | 優先順は宣言反応 → explore の生成物(低レベル TS を優先)→ mode-follow の TS 候補。宣言反応は必ず評価する。それ以外は同じレベルの DFT 極小の対(または 1 つの basin)で、ΔE_rxn ≤ reaction_window_kcal、端点の間で結合が変わるものだけ(組成あたり 6 件まで)。未宣言のねじれ・配座変化・鏡像化は仮説にしない(Curtin–Hammett)。発見の仮説の両端はその発見の `source_species` と `product_species`(§6.2 の原子の対応)で、両端が同じ basin の発見は、結合の相手が入れ替わるときだけ縮退反応にする。仮説の単位は case キー(`pair_key`。§6.2 の判定の粒度と同じで、両端の状態が違えば 2 つの状態の組、同じなら 2 つの極小の組)で、宣言反応はそれぞれ残し、すでにある仮説のキーの発見は新しい仮説にしない。どの仮説も、そのキーの全候補の低レベル TS を優先順に `low_level_ts` に持ち(持っている TS と構造で 1 つの basin に入るものは除く)、DFT 段の mode-follow の鞍点(検証済み)の最初のものを `ts_calc` にする(宣言反応にも貸す)。縮退・結合変化・ねじれは basin の最適化構造を端点の原子順と掌性に並べた構造で判定する。explore の陰性結果は仮説を棄却しない(GFN2 で障壁がなく発見にならない反応も、宣言すれば評価する) |
| reaction-paths | saddle は `trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2`、`moddir 1`(整えた Hessian の負のモード)。maxiter は継続せず、最終フレームと最後の歩のエネルギーを持つ `Failure` を返す(再開は case が §6.2 の上限と段数で決める)。ZTS は `nbeads 9`、`stepsize 0.05`、`interpol 3`、`freeze1` / `freezeN`、maxiter 20 のチャンクで、NWChem の収束判定は使わない。`pysis_neb` は IDPP を初期経路とする Cartesian の CI-NEB(`opt: lbfgs`、max_cycles 100)と rsprfo の TSOpt。QRC の振幅の上限 0.4 Å は初期 Hessian を受け付ける距離 0.5 Å 以下にする |
| sp | 既定の pipeline(discover、known_endpoints)は paths の後に M06-2X-D3(0)/def2-TZVPD の 1 つ(§7.3)。`methods` の各 LOT で、反応が読む点(`thermo.reaction_points`)だけを計算する: connected か barrierless の outcome の反応の TS と、その自分の点と鎖の井戸の状態(分離した単量体を含む)の DFT 極小すべて。ΔG_assoc だけが読む点(R_sep と barrierless の前駆錯体)は補助で、その SP の失敗は artifact にしない(値が空になるだけで、stage と rc を止めない)。freq の最終構造で計算し、parents に対象(minimum_id か TS の freq calc)を書く。thermo と手法パネルはこの parents だけで SP と停留点を対応付ける。DFT の SP は対象の freq の解から始める(`scf_guess`、基底が違えば射影。同じ freq を共有する対象は SP も 1 つ)。CCSD(T) は閉殻が RHF の ccsd モジュール、開殻が UHF 参照の TCE(P0d: NWChem の TCE ROHF-(T) は ROHF の正準軌道と固有値で f_ov 項を持たない定義で、標準の Watts–Gauss–Bartlett と OH で 0.108、CH3O で 0.195 kcal/mol 違う。TCE の UHF-CCSD(T) は PySCF と 2.1e-7 Eh 以内で一致する)。開殻の CCSD(T) のジョブ鍵は `reference: uhf` を持つ(ROHF の旧 deck の結果を再利用しない)。CCSD(T) は DFT の freq が `spin_ok` を満たす点だけにかける(BS の点には単一参照の状態がない。「低スピン結合でない」という条件は置かない: BS の点は必ず `spin_ok` に落ち、spin_ok を通る低スピン結合の点は CH3OO• のような結合したラジカルで、UHF-CCSD(T) が正しく書く)。凍結する芯は `freeze <n>` で明示する(n は原子ごとの NWChem の `freeze atomic` の芯の軌道数から def2-ECP が置き換えた分を引いた和。I は 4s4p の 4 で、全電子の Br の 3d と同じく 4d は相関させる。I⁻ + CH3I は 9、CH3O は 2)。Kr より重い原子を含む CCSD(T) のジョブ鍵は `frozen_core` を持つ(`freeze atomic` の旧 deck の結果を再利用しない)。WFT の deck は `memory_mb_per_rank` を GA 寄りに分ける(heap 5%、stack 25%、global 70%。NWChem は大きさごとに単位を要する。DFT は `memory total`) |
| thermo | GoodVibes 4.3.0 に自前の振動数を渡す。H は RRHO、S は Grimme の qRRHO(`qs`、`cutoff_cm1` 100)。σ・m・直線性は、下の「対称性」の点群 1 つから決める。振動数は freq の Hessian を点群の対称化した構造で射影し直したもの(直線なら 3N − 5 本)で、回転定数も同じ構造から求める。スケール因子は `vib_scale` の 1 つ(既定 1.0、振動数と ZPE の両方)。負モードは固定の規則: 極小は全モードを \|ν\| に、TS は最低モードを除いて \|ν\| にする。m = 2 の極小と TS の G に −RT ln 2 を加える(鏡像は同じ basin なので、鏡像対を 1 つとして数える)。σ の比は対称数で入る。種の値は 1 atm で、反応は `standard_states` ごとに換算する。層の SP は `thermo.single_points`(対象と Level ごとに 1 つ。2 つあれば誤り)で対応付ける。状態の G は、その状態の極小のうち反応の端と同じ LOT(`same_pes(state=False)`)で spin_contaminated でないものの最小の G で、鎖のどの井戸でも同じ定義。ただし、その状態の freq の LOT の極小(spin_contaminated でないもの)が 1 つでも層の SP か G を欠けば、状態の G は None(`energy_layer_missing`)で、その反応は `thermo_unavailable` になる(欠けた配座を黙って飛ばさない)。感度の幅は qs × cutoff 50/100/150。行の値(δG_eff、なければ ΔG_rxn)がなければ `thermo_unavailable`、自分の点の LOT の不一致は `mixed_level_of_theory`。トンネル補正はない(§10) |
| report | rankable な反応を δG_eff で並べ、感度の幅が重なれば同順位。(T, 標準状態) は report の `T_K`・`standard_state`、なければ thermo の最初の組。ranking.csv の列は rank、reaction_id、outcome、tier、rankable、T_K、standard_state、energy_level、dG_eff_kcal、reference、band_low_kcal、band_high_kcal、dG_act_kcal、dG_rxn_kcal、dG_act_vs_separated_kcal、torsional(結合変化のない段)、blockers、notes。barrierless の行は順位を持たず(blocker `outcome:barrierless_at_resolution`)、report.html の「Barrierless steps」の表(ΔE_rxn、ΔG_rxn、ΔG_assoc。会合は捕獲律速)に出る。エネルギーが 2 つ以上の LOT にあれば `dE_act_panel_min_kcal`・`dE_act_panel_max_kcal` を足す(既定の pipeline では停留点レベルとエネルギー層の 2 つ。手法の幅は列で示し、順位は止めない)。coverage.csv は機構別の試行・生成物・陰性理由・FailureKind |

結合と状態(`hfauto/chemistry/topology.py`): 結合は r < r_thr = r_cov,i + r_cov,j + 0.4 Å(Cordero の共有結合半径)で 1 つの構造だけから決める。状態ラベルは断片の組成式と WL ハッシュ。FHF⁻ と I3⁻ は 1 断片、ハロゲン結合(I···N 2.8 Å)は非結合。イオン–双極子錯体は許容値の内側なら 1 断片になり、GFN2 の強い H 結合錯体は閾値の近くでラベルが分かれうる(DFT 極小には影響しない)。2 つの構造の間の結合変化(`bond_changes`)は、r − r_thr が一方で ≥ +`RESOLVED_A`、他方で ≤ −`RESOLVED_A`(0.1 Å、モジュール定数)の対だけで、向きに対称である。1 つの basin の中をしきい値が横切るだけの違いは変化としない(根拠は S19 の N···H で、GFN2 と PBE0 の差が 0.15 Å)。結合変化の定義はこの 1 つで、仮説と分割の子のねじれの判定(`torsional`。未宣言のねじれは仮説にしない)、会合の形、ρ と χ の結合、explore の失われた状態がこれを使う。結合グラフと状態ラベルは帯を持たないので、結合変化のない 2 つの構造が別の状態ラベルを持つことがある。ラベルが割れても別の状態になるだけで、別の状態を誤って 1 つにまとめることはない。

対称性(`hfauto/chemistry/symmetry.py` の `analyze`): 熱化学の点群は、freq の構造と Hessian からここで 1 つだけ決める。
- 候補: libmsym(pymsym 0.3.5 の点群と対称操作。GoodVibes と同じ検出器)を 3 つの設定(equivalence 2e-3、7 つのしきい値すべて 1e-2、すべて 5e-2)で走らせて得た点群と、その対称化した構造。例外は候補なしとする。直線の候補は自前で対称化する(重心を通る主軸への射影が C∞v、それと反転の平均が D∞h)。
- 受理: 対称化した構造を、原子順を保ったまま回転だけで元の構造に重ね、その変位 δx が ½δxᵀ|H|δx ≤ 5e-5 Eh(|H| は Hessian の固有値を絶対値にしたもの)かつ RMSD ≤ 0.05 Å を満たすとき。2 つの極小を同じ basin とみなす `identity` の基準と同じで、新しい数はない。
- 選択: 受理した候補のうち位数が最大のもの(直線群を先に)。構造そのもの(C1)は常に受理する。
- 規約: σ は、受理した点群の対称操作のうち恒等操作と真回転(C_n)の数である(直線は C∞v 1、D∞h 2。Fernández-Ramos et al. 2007、Gilson & Irikura 2010。名前の表は使わない)。m は、点群に回映操作(鏡映、反転、S_n)がなければ 2、あれば 1(Fernández-Ramos et al., Theor. Chem. Acc. 118, 813 (2007))。直線性と外部自由度の数も同じ点群から決める。basin は鏡像を含むので、m = 2 の basin は鏡像対の両方を表す。受理が basin と同じ尺度なので、m と basin の数え方は食い違わない。
- 点群のない場面(エンジンが自分の freq を解析するとき、試行の方向、整えた Hessian)だけは、直線性を外部自由度の数値の階数(`vibrations._EXTERNAL_RANK_TOL`)で決める。

会合と低スピン結合の組成(分析 §3 X4。`electronic_state.low_spin_coupled`、`hypotheses`、`drivers/reaction_case/paths.py`):
- **低スピン結合の組成**: 開殻の成分が 2 つ以上あり、宣言した多重度が高スピンの結合 Σ(m_i − 1) + 1 より小さい組成(CH3• + O2 の二重項。四重項は高スピン)。拒否しない。結合した異性体(CH3OO• など)は UKS で正しく書けるので、その間の反応は通常の case である。成分が離れた領域だけは UKS が broken-symmetry(BS)解になり、⟨S²⟩ は (S_A − S_B)(S_A − S_B + 1) + 2S_B(CH3·O2 で 1.75。二重項 2/3 と四重項 1/3 の混合)になる。これは誤りではなく、この電子構造から必ず出る値である。そのため非結合の錯体の極小は `spin_ok` で落ち、熱化学には使えない。
- **一重項**: 一重項は RKS で計算するので BS 解がない。このクラスの一重項は入口で止める(`low_spin_singlet_unsupported`。OH• + OH•、O2 + O2 の一重項)。再結合の生成物は単量体として宣言すれば計算できる。
- **会合の仮説**: 反応物から生成物への結合変化が、反応物の 2 つの断片を結ぶ形成 1 本だけで切断がなく、反応物の断片が宣言した組成の単量体と一致し(`thermo.separated_states`)、各単量体の状態の DFT 極小が錯体と同じ LOT にあれば、その仮説を会合にする。起動条件は仮説の形だけで決め、出所(宣言、explore、R6 の relaxation)には依らない。反応物側は分離した単量体である(`ReactionRecord.monomers` に各状態の最低の DFT 極小を個数だけ並べ、`reactants` もその組成にする)。錯体の極小は `minima[0]`(仮説の同一性、重複除去、状態の参照)に残るが、エネルギーと熱化学の端点にはしない。障壁のない会合では反応物は漸近であって極小ではない。低スピン結合の組成では、錯体の極小は BS の汚染による人工物でもある(S5 の CH3···O2、C···O 2.84 Å、⟨S²⟩ 1.71)。会合の形は錯体の DFT 極小の構造で判定する。錯体が自分の DFT 極小を持たず付加体に落ちた宣言反応(BH3 + NH3)は、錯体の入力構造で判定する。このとき `minima[0]` は付加体の極小で、case の反応物端はその入力構造である(行 2 の same_basin は会合には効かない)。
- **スキャン**(行 9 の SCREEN): 形成する原子対 (i, j) の距離を、付加体の値 r_P + 1.5 Å(`SCAN_REACH_A`)から付加体まで、等間隔の 8 点(`SCAN_POINTS`。最後の点は付加体の極小)で内側へたどる。最初の点は、付加体の j を含む断片(錯体での断片)を i→j の方向に剛体で引き離した構造で、次の点は前の点の最適構造を同じように押し込んだ構造である。各点は r_ij を固定した opt(§5 の `fixed_bond`)で、SCF は前の点の vectors から始める(BS の対が連続した枝に留まる)。r_P + 1.5 Å は Coulson–Fischer 点(CH3·O2 の C–O で約 2.2 Å)より外側にある。内側から始めると単調性を作るだけになる。
- **判定**: プロファイル [Σ E(単量体)、スキャンの点、E(付加体)] を `profile.judge` で分ける。単調なら barrierless_at_resolution(行 6、`scan:barrierless`)。山があれば最大点を放物線補間した構造(`scan_hei`、種の出所は `path_hei`)を種にして REFINE_SADDLE に進み、その後は通常の case と同じである。string は走らせない(string は 2 つの極小を結ぶもので、単量体の端は極小ではない)。点が 1 つでも失敗したら unavailable(`scan_point`)で、点を落として判定することはしない。中点の密化は他と同じく 1 回: 新しい節点は隣の 2 点の線形補間の中点での SP で、外側の点の SCF の解からその点の座標系で始め、低スピン結合の対では AP で射影する(注記 `scan_mid:ap:<多重度>`)。単量体の和は最初のスキャン点と同じ構造なので、その区間には足さない。
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
- **層の選び方**(分析 [2026-09-30_remaining_issues_analysis.md](reviews/2026-09-30_remaining_issues_analysis.md) §2.6、§3 X7): PBE0 は非局在化誤差で電子の広がった TS を過安定化する(障壁の平均符号付き誤差は H 移動 −4.2、重原子移動 −6.6、求核置換 −1.9 kcal/mol。Zhao & Truhlar, J. Phys. Chem. A 109, 2012 (2005))。誤差の主因は汎関数で、SVPD → TZVPD の変化は ≤ 1.1 kcal/mol だった。4 反応の ΔE‡ の CCSD(T)/def2-TZVPD との平均絶対誤差は PBE0/SVPD 2.68、PBE0/TZVPD 2.36、ωB97X-D3/TZVPD 1.56、M06-2X/TZVPD 1.03 kcal/mol(§7.4 の表。ROHF-CCSD(T) で測った値で、CH3O を UHF にすると各汎関数の差は約 0.1 動く)。TMA·(HF)₂(17 原子、327 基底関数、4 rank)の 1 点は M06-2X 255 s、PBE0/TZVPD 312 s、ωB97X-D3 1,689 s。M06-2X の grid は fine と xfine で障壁の差が ≤ 0.01 kcal/mol。
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

- **PES の形は PBE0 のまま**: 鞍点の有無と位置、PES の平らさ、SCREEN・FIND_PATH・会合のスキャン・行 6 の barrierless の判定は PBE0-D3BJ/def2-SVPD の PES で決まり、層では確かめない。PBE0 に鞍点がない(山が解像度未満の)反応の障壁は層では取り戻せない(H + C2H4: PBE0 の山 +0.60、参照 1.72 kcal/mol。validation.md §4)。PBE0 の非局在化誤差は障壁の値だけでなく PES の形も変え、H 引き抜きの TS を早く平らにしうる(推論。S6 の平らな PES、G1)。停留点の汎関数は変えない: M06-2X と RSH では NWChem が解析 Hessian を使わず freq が高くつき、SP//PBE0 の構造の誤差は S6 で +0.31 kcal/mol(PBE0 の QRC の軌跡上の CCSD(T) の最大と PBE0 の TS の差)と、汎関数の誤差より 1 桁小さい。
- **発見の段の窓は層より前のエネルギーで切る**: minima(dft) の `window_kcal` 6.0 は GFN2(`rerank_sp` なら PBE0/SVPD の SP)、仮説・行 3 の `reaction_window_kcal` 40 は PBE0/SVPD の ΔE_rxn(explore は GFN2)で決まる。窓で落ちた構造と仮説は層に届かず、層で順位に戻らないので、窓は停留点レベルの誤差より広くとる。PBE0/SVPD の CCSD(T)/def2-TZVPD との差は障壁で最大 −6.2(OH + CH4)、S6 の分離基準の ΔE_rxn で −2.2 kcal/mol(PBE0 −15.02、W3 の freq のエネルギー)で、40 はその数倍ある。6.0 は同じ状態(同じ結合)の配座の差に効くので、結合の組み替えに伴う非局在化誤差は主に入らない(推論。配座のエネルギー差の誤差は GFN2 でも PBE0 でも測っていない)。
- 低スピン結合の組成のスキャンは AP-UKS のエネルギーで判定するので、障壁の有無の定性的な判定として読む。AP は BS 解への高スピンの混入だけを除く 2 状態の近似で、多参照性(静的相関)は直さない。
- 電荷系では PBE0 の非局在化誤差で電荷の広がった TS・錯体が低く出やすく(SN2 の ΔE‡ −2.99)、D3 は電荷に依らない。層はこの値を直すが(+0.46)、PES の形は PBE0 のままである。
- 停留点レベルを def2-SVPD にしたのは、FHF⁻ の De が拡散関数のない SVP では 21.7 kcal/mol 過大になるが SVPD では実験値の誤差内に入るからである。
- ΔG_assoc は BSSE を補正しない値で、結合を過大に見積もることがある。qRRHO の扱いによる幅は 298 K で約 1.2 kcal/mol(TMA·(HF)₂)。
- 電子分配関数は g = 2S+1 だけで、開殻原子と ²Π ラジカルの軌道縮重とスピン軌道補正は含めない。

## 8. 実行基盤

- **JobStore**(`hfauto/execution/jobstore.py`): run ごとの内容アドレス型キャッシュで、`<run>/jobs/<k[:2]>/<key>/` に job.json、result.json、attempt_NN/ を置く。キーは正規化 JSON {engine, version_pin, kind, key_payload} の sha256 で、結果を変える入力(method、構造の指紋、パラメータ、入力ファイルの sha)をすべて含み、`ExecutionSpec`(ranks、メモリ、timeout、パス)は含めない。再利用時はファイルの sha を照合する。Evidence の来歴は `job_key` と `FileRef`(run ディレクトリからの相対パスと sha256。作るのは `JobStore.file_ref` だけ)の 2 つで持つが、ジョブの出力と Hessian の FileRef は常にそのジョブの `jobs/<k[:2]>/<key>/attempt_NN/` を指すので、2 つは食い違わない。JobStore は run の中の記憶で、壁時計に依存する失敗(timeout)以外は ladder の最後の失敗も(saddle の最後の構造ごと)記録して同じ鍵に同じ結果を返し、`--retry-failed` で指定した種類だけを再実行する。run の間の決定性は保証しない(並列の NWChem と CREST はビット単位では再現しない)。キーにコードの版数は入らないので、パーサを直した後は `--retry-failed incomplete_output` で取り直す。
- **ladder**(`hfauto/execution/jobs.py` の `LADDER`。これ以外の自動再試行はない): 駆動系ジョブの継続はすべて最後の frame から: timeout(opt と saddle)と opt の geometry_maxiter は 2 回まで、opt の autoz(input_invalid)は Cartesian で 1 回(vectors だけを引き継ぎ、driver の Hessian は捨てる。固定結合は spring bond にする)。frame がまだない autoz と saddle の autoz は同じ始点から Cartesian で 1 回。saddle の maxiter は継続せず、最後の frame とその歩のエネルギーを持つ Failure を返す(その frame から新しい Hessian でやり直すかは case が決める: 経路のエネルギーの上限の内側で、種の続きの段数が残るときだけ。§6.2)。scf_not_converged は DFT だけ 1 回: 同じ始点から atomic guess + `cgmin` + 通常の SCF(§7.1。駆動系のジョブも attempt の始点からで、vectors だけを捨て継続の driver の Hessian は保つ)。発散した vectors からは再開せず、WFT は再試行しない。attempt の timeout は site の `timeout_s` である。
- **失敗の閉じ込め**: `adapter.parse` の例外はそのジョブの `Failure(incomplete_output, 'parse:<型>: <1 行目>')` になる。反応ケースの例外は `drive_case` が閉じ込め、そのケースを unresolved(`error:<型>`)にして log.jsonl に 1 行書き、それまでに登録した極小・species・計算とともに出して、ほかのケースを続ける(共有の `Registry` への登録は case のスレッドだけが行い、`CaseRuntime.map` はエンジンの呼び出しだけを並列にする)。消費する型の入力がない stage は 0 artifact で done になり、上流の失敗を warning に出して後続の stage(report を含む)へ進む。`HFAUTO_STRICT=1`(テスト専用、`tests/conftest.py`)は前 2 つを再送出にする。
- **並行実行**: `JobRunner` のセマフォで実行中のジョブの ranks × threads の合計を site の `cores` 以下に保つ。CREST と explore の単位は `thread_map` で並列(結果は入力順)。反応ケースは直列だが、ケースの中の独立な NWChem ジョブ(SCREEN の内部 9 点と barrierless の中点 2 点の SP(開殻でも各点は端の vectors から直接始めるので並列のまま)、非対称な QRC の両側)は `CaseRuntime.map`(`thread_map`)で同時に走らせる。項目は `cores` 個ずつの波で走り、各項目は cores // (その波の項目数) の MPI rank を使い(4 コアで 9 点の SP は 4×1、4×1、1×4)、timeout は同じ倍率で延ばす。NWChem は常に mpirun に `--bind-to none` を付けて走らせる(Open MPI 5 は np ≤ 2 の mpirun をすべてコア 0 から束ねるので、rank を減らしたジョブや、ranks < cores の site で同時に走る全 rank のジョブがコアを取り合う)。rank は `ExecutionSpec` なので鍵に入らず、同じ鍵のエネルギーは直列と 1e-8 Eh 以内で一致する。極小の stage(xTB と DFT)と mode-follow のような直列の連鎖は全 rank のまま直列にする(大きさの違う化学種に同じ rank を割ると遅くなる。validation.md §8)。`SiteLock`(`<scratch_root>/.hfauto_site.lock`)は `hfauto run` の間ずっと保持し、別の run の同時起動を拒否する。
- **シグナル**: `run_command` は外部プログラムを新しいセッションで起動して登録する。`hfauto run` は SIGINT・SIGTERM・SIGHUP を受けると(POSIX)登録済みのプロセスグループをすべて止め、実行中の stage を `incomplete` と記録して rc 1 で終わる(Windows の Ctrl-C も同じ)。止めたジョブは JobStore に残らないので、再実行すると done でない stage(failed・incomplete)を JobStore から取り直す。
- **run ディレクトリ**: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/cases/`、`<run>/run_state.json`(実行順の stage 状態とジョブの統計)、`<run>/resolved_config.yaml`(設定と code_version = git の sha、差分があれば -dirty)。stage の入力 view は、run_state でその stage より前にある done の manifest を実行順に合わせたもの(別の pipeline を同じ run に追記できる)。done で入力と設定の sha が変わらない stage は飛ばし(再開。pending・failed・incomplete は走らせる)、`--from X` は X 以降を取り直す。記録の型は未知のフィールドを拒否する(`extra="forbid"`)ので、型が変わる前の run は途中の stage からは再開できない。round 7 より前の run は新しい run dir で流し直す。それより後の run は、記録の欄が変わっていても、最初の stage からの `--from` なら JobStore を使って取り直せる(再生ゲートはこの形)。
- **検証**(`validation/`): 期待する結論は `validation/cases.yaml`(pipeline の連鎖、反応ごとの両端の断片の組成式の式・outcome・δG_eff と許容、境界の系の `twice`、置き換えた case の `superseded`)に置く。BH76 の参照は `validation/bh76/subset.yaml`(P0e の宣言ファイルそのまま。GMTKN55 の構造は `bh76/xyz/`、出所と sha は `bh76/SOURCES.json`)で、case は `bh76: <反応>.<向き>` で行を指す。比較は純関数 `compare`(結論と式)と `deviations`(参照と矛盾する outcome。合格に数えない)の 2 つで、`validation/check.py RUNS_DIR` が run dir の集まりにかける。`validation/replay.py OUT --src DIR --strict` は記録済みの run を JobStore から再生し(新しいジョブ 0 件と同じ記録)、`--src` なしなら新しい run を流す(`twice` の case は `<name>.2` にも)。OUT は `/home/user/hfauto_r10` の下で、QM のロックを取る。

## 9. 設定

4 つの層は中身が重ならないので、マージも優先順位もない。結果を変える設定は method・system・pipeline に置き(ジョブの鍵か stage の設定 sha に入る)、site には実行設定だけを置く。system と pipeline の YAML は未知のキーを実行前に拒否する。

| 層 | 置き場所 | キー |
|---|---|---|
| site | `configs/sites/` | `site`、`scratch_root`、`cores`、`engines.<registry 名>`: `version`(版数の pin)、`executables`、`execution`(`ranks`、`threads`、`memory_mb_per_rank`、`timeout_s`、`env`)、`python`(worker)、`scratch_dir`。絶対パスはこの層だけに書く |
| method | `configs/methods/` | `id`(ファイル名と同じ)、`kind`(dft / xtb / wft)、`functional`、`wft_method`(ccsd(t))、`basis`、`dispersion`(d3zero / d3bj)、`gfn`、`electronic_temperature_K`、`grid`、`scf_energy_tol`。gfn2(低レベル)、pbe0-d3bj_def2-svpd(停留点)、m06-2x-d3_def2-tzvpd(順位のエネルギー層)、pbe0-d3bj_def2-tzvpd と wb97x-d3_def2-tzvpd(手法パネル)、ccsd-t_def2-tzvpd(opt-in) |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | `system_id`、`species`(`id`、`xyz` か `smiles`、`charge`、`multiplicity`、`role`: monomer / endpoint)、`compositions`(`id`、`components`、`multiplicity`)、`reactions`(`id`、`reactant`、`product`、`coordinate`、`torsional`) |
| pipeline | `configs/pipelines/` | `pipeline_id`、`gates`(knob)、`stages`(`id`、`stage` と下の stage 表の設定キー) |

**knob は 7 つ**だけである(検証は既定値で行っている)。

| knob | 置き場所 | 既定値 | 意味 |
|---|---|---|---|
| `noise_cm1` | `gates:` | 10 | 停留点で、これより小さい虚振動(cm⁻¹)は数値の雑音(停留性は §6.1 で別に判定する) |
| `saddle_cm1` | `gates:` | 50 | これより大きい虚振動(cm⁻¹)は鞍点級 |
| `resolution_kcal` | `gates:` | 1.0 | 経路の山と井戸を数える深さ(kcal/mol) |
| `reaction_window_kcal` | `gates:` | 40 | 発見・仮説・判断表(行 3)の ΔE_rxn の上限(kcal/mol) |
| `spin_tol` | `gates:` | 0.1 | \|⟨S²⟩ − S(S+1)\| の許容幅(開殻のプロファイルの両端のスピンの結合が違うかと、最大点の ⟨S²⟩ が隣から離れているか(枝の跳び)の判定にも使う。§7.1) |
| `max_saddle_attempts` | reaction-paths の `policy:` | 2 | 1 case の saddle の試行数(種 1 つが 1 試行。子反応は 0 から。maxiter からの再開は同じ試行の中) |
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
| `explore` | `hfauto.stages.explore:ExploreStage` | minimum, species | discovery, species | engine*, method*, sources_per_state, max_trials_per_source, settings |
| `reaction-paths` | `hfauto.stages.reaction_paths:ReactionPathsStage` | minimum, species | reaction, calculation, minimum, species | method*, engines*, screen, policy, reaction_ids |
| `sp` | `hfauto.stages.single_point:SinglePointStage` | minimum, reaction | calculation | engine*, methods* |
| `thermo` | `hfauto.stages.thermochemistry:ThermoStage` | minimum, calculation | species_thermo, reaction_thermo | settings, energy_method, temperatures_K, standard_states |
| `report` | `hfauto.stages.report:ReportStage` | — | report | T_K, standard_state |

入れ子のキー: minima の `select`(`include`: window / all、`per_state` 3、`window_kcal` 6、`rerank_sp`、`rerank_top` 8)と `init_hessian`(`engine`、`method`)、explore の `settings`(ReaDuct の `max_scf_iterations`、`electronic_temperature_K`、`scc_retry_temperature_K`、`imag_cutoff_cm1`、`timeout_s`)、reaction-paths の `engines`(`qm`、`saddle`、`path`)・`screen`(`method`、`qm`、`path`)・`policy`(上の knob)、thermo の `settings`(`qs`、`cutoff_cm1`、`vib_scale`、`sensitivity`)。

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
| 4 コアでの極小の stage の同時実行 | 直列を正とする(§8) |
| 壁時計 | 共有ホストの壁時計は負荷に依る。比較は歩数と core 秒で行う |
| Windows のアプリケーション制御 | 品質ゲートは WSL の prod venv で判定する([environment.md](environment.md) §4) |
| 物理行数の目標 | 指標は SLOC(radon raw)と、1 つの概念を実装している場所の数 |

既知の制限(範囲の中で、結果に注記なしに効きうるもの。一覧と状態は最新の結果報告 [2026-10-02_round9_result.md](reviews/2026-10-02_round9_result.md) §5):

- 一次の鞍点の受理に停留性をかけない(§6.1)。40 cm⁻¹ 未満の大振幅のモードでは二次モデルが成り立たない(acac の TS で ΔE_N 8.6e-3 Eh)。
- 状態ラベルは立体を区別しない(§6.2)。
- ⟨S²⟩ による状態の照合(§7.3)は、同じ多重度の別の空間状態(CH3OO• の X̃ と Ã など)を区別しない。射影の guess は親の解に留まるので、親の freq が励起の SCF 解にあれば層の SP もそこに留まる([validation.md](validation.md) §11.2)。
- 点群の受理は basin の基準なので、ほぼ Cs の平らな TS では雑音で Cs と C1 が入れ替わり、δG_eff が RT ln 2 だけ動く(S6 の TS で 0.41 kcal/mol)。
- string(ZTS)の bead は ⟨S²⟩ を持たず、端の SCF の解を guess にしない。SCF の救済の後の通常の SCF も string の deck には入らない(§7.1)。

## 11. 改良の設計(段階 0 のプローブで確定)

レビュー [2026-10-04_platform_review_2.md](reviews/2026-10-04_platform_review_2.md) の計画のうち、段階 0 のプローブ(証拠は `/home/user/hfauto_r10/probe/P0a`〜`P0e`)で採否が決まった設計である。実装は wave の順に進め、統合のたびに該当する §0〜§10 を書き換えて本節から外す。

### 11.1 終了コード

- rc は 0 = すべての stage が done(化学の未解決や CREST・SCC の失敗は coverage に出す)、1 = incomplete か `error:*` の項目、2 = stage の失敗。項目の例外を閉じ込めるのは `StageRuntime` の 1 か所だけである。

### 11.2 電子状態の来歴(残り)

W2 で SCF の来歴、3 段の救済、層とパネルの状態の照合、開殻 CCSD(T) の UHF 参照を §5・§7.1・§7.3・§8 に取り込んだ。残りは次の 1 つ:
- QRC の側の opt も TS の freq の解から始める(`connection._sides` の `scf_guess=freq` の 1 行)。paths の全 QRC のジョブ鍵が変わるので、W7 の鍵の移行で入れる。

### 11.3 化学の判定

- 貸された TS(他の候補の低レベル TS、mode-follow の `ts_calc`)の χ は、貸し手のラベル付きの両端で測る。貸し手の両端は今の記録にないので、W6 の登録簿の辺に自分の両端を持たせ、そのとき種に結合変化を持たせる。それまでは case の両端の変化 `Ctx.change()` で測る(質量加重の χ は W2 で §3 に取り込んだ。等価な原子の Hungarian の付け替えは測って棄却した: [validation.md](validation.md) §11.2)。
- 点群の受理は、対称化の変位のエネルギーを零点振幅と比べて決める。電子分配関数とスピン軌道のエネルギーは NIST の準位から取る。H は QH にする。立体は状態に含める。
- 停留性は信頼領域のモデル 1 つで認証する。継続の仕組みと下りの仕組みはそれぞれ 1 つにし、止める規則も 1 つにする。

### 11.4 順位の量と計算点(残り)

W3 で反応が読む点、点の鎖と共通のゼロ、宣言した単量体の分離した点、捕獲律速の表、`profile.judge` を §3・§7 に取り込んだ。残り:
- 捕獲律速の表に、ランキングの LOT で経路の節点をたどった判定(`layer_verdict`、`dE_act_path_kcal`)を足す(X5-3)。
- 事実(状態の組の結論)は stage の中の登録簿 1 つに置き、種は 1 本の順序付きの流れで渡す(W6)。

### 11.5 配座と錯体(P0b)

- 組成の探索は GFN2 の CREST のまま続ける(`--nci --quick --notopo <全原子> --noopt`)。GFN-FF は採らなかった。s10 の 2 つの組成では、入力状態の候補が 1 つずつしか出なかった(GFN2 も同じ)。S5・S6 は rc 0 になったが、出た構造は断片がほぼ離れていた(最短接触の中央値 5.5 Å、GFN-FF のエネルギーの幅 0.7 kcal/mol 未満)。再検討するときは、ラベルではなく DFT の basin の集合で判定する。
- seed を組成ごとに 1 つ(重心を結ぶ線上の接触)にする案(U2-P4)は採らなかった。前提の U2-P1(GFN-FF)が落ち、接触の seed 1 つでは S5・S6 の GFN2 の screen がどちらも失敗して、両系の反応がすべて消える(W3-3)。錯体が screen の basin を経ずに explore・DFT に届く経路(U2-I1)ができたら再検討する。lone-pair 円錐と元素の表は W3 で削除した。

### 11.6 探索(P0c)

- 反応は、結合グラフの編集を列挙して作る。条件は次のとおり:
  - 形成 2 本以下、切断 2 本以下、形成は 1 本以上;
  - 正味の配位数が `fits()` に収まる;
  - 変化した原子が 6 員環以下の環か、和のグラフの鎖の上にある;
  - 形成する対はグラフ距離 2 以上。
  
  類は RDKit の正準化による原子クラスで定める。反応の型(移動、置換、[3,3]、[4+2]、1,2-脱離、SN2、リレー)は入力ではなく、列挙の結果として出る。
- 生成物のグラフには Lewis 構造があることを求める。形式電荷は ±1、不対電子は m − 1、価数の集合は S {2,4,6}、P {3,5} とし、周期 3 以上の 15・16 族は電荷に応じて価数をずらす。電荷分離をどこまで許すかは、P0f で事前に宣言した規則で決める。電荷分離を制限しないと、ibuprofen は 10,708 類になり(88.6% が電荷分離を要する)、全部を列挙すると DFT の freq 1 本の 2.5〜3 倍の費用がかかる。
- 予算は件数 1 つ(`max_trials`、全世代の合計、既定 3000)である。R2 の 1 試行の平均 6.35 core-s で約 19k core-s になり、ibuprofen 大の DFT freq 1 本(2.3 万〜2.8 万 core-s)より少ない。予算が尽きたら、結合変化の少ない類から順に、同数なら drive の値の順に試す。切った類は `not_attempted` として数える。薬物の大きさの分子では、網羅が予算で限られることを報告に出す。
- 発見は辺として持つ(低レベルの TS があるものもないものもある)。世代を重ねるのは低レベルの中だけである。DFT の入口では、DFT//xTB の一点計算で 1 回だけ窓をかける。
- 状態ラベルは Hill 式と、RDKit の正準 SMILES の hash で作る。生成物を比べるときは状態の水準で比べる(原子を対応させた編集には、傍観している H の入れ替えが入ることがある)。

### 11.7 エネルギー層と BH76(P0e)

- 参照は GMTKN55 の BH76 から取り、このうち 5 反応を使う。選定に使うのは OH + CH4、HCN ⇄ HNC、Cl⁻ + CH3Cl の 3 つで、H + C2H4 と H + H2 の 2 つは選定に使わず検証だけに使う。比べる量は層の電子エネルギーの障壁で、ZPE・熱・スピン軌道の項は含めない。障壁は pipeline の PBE0 の停留構造で求め、ゼロ(分離か錯体)は行ごとに決まっている。
- 許容は 1 行あたり 1.0 kcal/mol とする。このうち構造の寄与 g は 0.3 以下でなければならない。g は、BH76 の QCISD/MG3 構造での層の SP との差で測る。g > 0.3 の行は geometry-limited として報告し、合格にはしない。Minnesota/08 の値との差は横に並べて出す。
- エネルギー層の候補は M06-2X-D3(0)、revM06、ωB97X-D3(def2-TZVPD)である。採るのは次の 3 つをすべて満たす汎関数だけで、満たさなければ M06-2X を残し、その兼ね合いを書く:
  - 選定の行の類別 MUE が現行以下;
  - 開殻の弱結合錯体で SCF の失敗が 0;
  - 費用が 2 倍以内。
- ジョブ鍵を変える基盤の変更(deck の sha、版数の辞書、`canonical_json` の一本化、密度の収束、`levels`)は、最後に 1 回でまとめて移行する。
