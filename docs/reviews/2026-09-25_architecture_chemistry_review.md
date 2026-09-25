# hfauto アーキテクチャ・化学計算レビュー(2026-09-25)

- 対象: 作業ツリー(HEAD `423e77c` に対し未コミットの大規模リファクタを含む現状)。`hfauto`(約 34.6k 行)、`hfauto_ops`(約 3.4k 行)、`hfauto_viz`(約 2.1k 行)、`tests`(約 9.4k 行)、`configs/`、`runs/`(約 150 run)
- 方法: 静的解析 4 ツール(import-linter / pyrefly / radon / Ruff)を導入して計測した。そのうえで、アーキテクチャ 6 観点(A1〜A6)と化学計算 6 観点(C1〜C6)をレビューした。各指摘は別のエージェントが敵対的に検証した(コードの再読、ライブラリソースの確認、run 出力の再計算)。
- 検証結果: 127 件のうち confirmed 76 件、partially 51 件、refuted 0 件。partially の指摘は検証者が訂正した主張で記載し、重大度は検証後の値を使った。観点をまたぐ重複は統合し、本書では `AR-xx`(アーキテクチャ)と `CH-xx`(化学)の ID で扱う。元の ID は付録にある。
- 重大度: **致命**(科学的に誤った結果が公開または合格扱いになる)/ **高** / **中** / **低**
- 注意: import-linter の 11 契約は本レビューで新たに `pyproject.toml` へ追加したもので、HEAD 時点の契約ではない。違反数は「提案した層構造」に対する数値として読むこと。

---

## 0. エグゼクティブサマリー

### アーキテクチャ

1. **外枠は健全だが、実際の結合は型のない辞書で起きている。** `Stage.run(manifest, config, context) -> Manifest` の統一、41 stage の遅延登録、core が上位層を import しないこと、モジュール単位の循環 0 件は良い土台である。一方で stage 間の本当の契約は、文字列の `artifact_type` と `Artifact.qc/data`(`dict[str, Any]`)にある。qc キーは 260 種が書き込まれ、読まれるのは 47 種だけ。103 種(40%)はどこからも読まれない。`fallback_dummy` の判定式は 26 ファイルに複製されている。型付き Evidence モデルと 4 つの判定関数への集約が、最も効果の大きい構造改善である(AR-23, AR-24)。
2. **複雑度は少数の関数に集中している。** `ThermoStage.run` は CC 186・829 行で、実行テストは 0 件。ほかに `classify_reaction` CC 115、TS バックエンドの `search_ts` が CC 41〜55・350〜590 行ある。stages パッケージは関数平均 CC 15.0、F 評価 11 件で、backends/ts は 1 ファイル平均 570 SLOC である(AR-20, AR-12)。
3. **バックエンドは看板どおりには差し替えられない。** Protocol の参照は 0 件で、registry の戻り値に型がない。エンジン名は `workflow/reaction_state.py` と `core/environment.py` にハードコードされている。TS バックエンドは `NWChemEngine()` を 8 箇所で直接生成し、pysisyphus は NWChem 7.2.3/D3 に固定されている。chemistry/stages 層にも NWChem 固有のパーサと方針が入り込んでいる(AR-11〜AR-14)。
4. **連携計算は「判断」と「実行」の間が切れている。** `decide_reaction_case` は決定論的で良い状態機械だが、TS 探索以外の `next_action` を実行するディスパッチャがない。たとえば `run_fixed_geometry_frequency` は行き止まりになる。runs 137 件のうち runner 経由は 68 件だけで、m3 系列は `scripts/`(git 未追跡)と単段 `run-stage` の手動連鎖で回されていた。HPC の配列実行経路は、存在しない CLI コマンドを呼んでおり end-to-end で動かない(AR-05, AR-03)。
5. **長時間ジョブの保全がない。** opt+freq を 1 ジョブにまとめている、NWChem を常に `start` で起動する、`.movecs/.hess` を一度も再利用しない、再試行を初期構造からやり直す、の 4 点による。記録上の外部計算時間 124.1 h のうち 31%(38.0 h)が失敗またはタイムアウトしたジョブだった(AR-08)。
6. **契約・診断・テストが計算本体より重い。** 評価・監査系は 28 ファイル・6.9k SLOC あり、計算駆動系(17 ファイル・9.2k SLOC)の約 75% に当たる。IRC ステージにはゲートが 7 段あるのに、「分岐が TS から下降したか」「座標系が一致しているか」という決定的な検査は 0 件である。テストは assert の約 52% がフラグやラベルの形の確認で、実出力のゴールデンテストは 0 本、熱化学の数値コアを直接確かめるテストも 0 本である(AR-29, AR-30, CH-03)。
7. **削減余地は約 4,500 行ある**(ORCA を削除する場合はさらに約 1,300 行)。内訳は、到達不能モジュール 35 個(確認済み約 1,600 行)、どのパイプラインでも使われない stage 群(kinetics/calibrate/rank の legacy 部分/connector-audit)、hfauto.hpc と hfauto_ops の二重実装、legacy HF 熱化学分岐、未参照の configs 15 本、ソースが消えた孤立 `__pycache__` 75 モジュール分である(4.5 節)。

### 化学計算

1. **長時間化の最大の原因は化学そのもの、つまり探していた TS が存在しなかったことにある。** 気相の amine·(HF)n では、n=1 は中性の水素結合錯体、n=2 は共有プロトンの単一 basin、n=3 は障壁なしでイオン対になる。n=1,2 にプロトン移動 TS はない。TMA/amine 系 86 run・約 70〜75 h を費やして、検証済み TS/IRC は 0 件だった。xTB 緩和スキャン(0.1〜0.2 h)と NT2 の単調上昇から事前に分かっていたのに、それを前段ゲートとして使っていなかった(CH-06, CH-15)。
2. **【致命】pysisyphus+QCEngine+NWChem 経路は座標系が一致していない。** NWChem が分子を再配向(実測 82.9°)したフレームのまま、勾配と Hessian が返っている。root0 固有ベクトルの重なりは 0.126 しかない。加えて grid は medium、SCF は 1e-6 で、端点とは別の PES になっている。IRC は分岐の下降を検査しないため、`hcn_irc_v5` は両分岐とも TS より上で終わったのに `irc_validated=True` になった。現行の全 production パイプラインの IRC がこの経路を使っている(CH-01〜CH-03)。
3. **【致命】誤った値が公開されている。** HONO trans→cis の ΔG‡ = 24.11 kcal/mol は、GoodVibes が NWChem 出力内の 2 つの振動ブロックを連結して ZPE を二重に数えた値である。正しくは約 12.2 kcal/mol で、それにもかかわらず Q4・`production_rank_eligible=True` で公開されている。HCN の ΔG_rxn = 18.05 kcal/mol は、HNC の回転エントロピーが負(−5.0 cal/mol·K)になる NWChem RRHO フォールバックの値で、正しくは 12.74。これも `scientific_rank_eligible=True` になっている(CH-38, CH-39)。
4. **【致命】ゲートの偽陰性が有効な結果を捨てている。** ReaDuct の低レベル TS は、射影していない Cartesian Hessian の固有値で虚振動を数えているため、真の 1 次鞍点 8 件(ReaDuct 自身の解析ではすべて虚振動 1 本)を全件棄却し、低レベル IRC は一度も実行されていない。さらに `reaction-plan` が作る `basin_assessment` から必須キーが抜けた回帰があり、production の `reaction-classify` は全反応を `unresolved` にする(HONO で再現)(CH-26, CH-05)。
5. **数値設定と実行環境が過剰かつ不適切である。** xfine/SCF 1e-8 は Hessian を 3.6 倍にするが、結論は変わらない(HCN で fine との差は 6e-5 kcal/mol)。同じ構造の Hessian を 900〜3,400 s かけて 5〜6 回重複計算している。WSL は 4 vCPU/11 GB で、そこへ mpirun を同時起動して最初の SCF が 17 倍遅くなった。NWChem driver は `noautoz`(Cartesian)固定で、初期 Hessian は対角である(CH-11〜CH-14)。
6. **経路法の選び方と収束判定に欠陥がある。** NWChem NEB は impose・IDPP・climbing image を使っておらず(7.2.3 の CI 関数は呼ばれないコード)、実反応 9 件で収束は 0 件だった。NWChem ZTS は変位だけで収束を判定するので、gmax が 10 倍に悪化しても "converged" になる。単調だった経路は端点の線形補間ブラケットに回され、+51 kJ/mol の人工極大が繰り返し TS seed になった(CH-08, CH-09, CH-06)。
7. **有効だった計算もある。** HCN の string→NWChem saddle(849 s + 6 step。ΔE‡ 46.7 kcal/mol で CCSD(T) 44.6±1.0 とよく一致)、HONO の saddle(−680i)、NH3 反転、xTB 緩和スキャン、固定構造の method/grid パネル、CREST `--nci`、同一停留点の判定、minimum-mode-follow、GoodVibes qRRHO(ただしゲート付きで)。これらを中心にプロトコルを組み直すべきである(3.3 節の推奨プロトコル)。

---

## 1. 導入したツールと設定

### 1.1 導入と変更したファイル

- 既存の venv `.venv-win-qa`(Python 3.12.13。pip はないため uv を使用)に導入した:
  `uv pip install --python .venv-win-qa/Scripts/python.exe import-linter pyrefly radon ruff grimp`
  - 新規: import-linter 2.15、grimp 3.17、pyrefly 1.3.1、radon 6.0.1(依存: mando 0.7.1)、click 8.5.0
  - 既存: ruff 0.16.1、pytest 9.1.1、typer 0.27.1
- 変更したのは **`pyproject.toml` だけ**(スナップショットは `tools/pyproject_after.toml`)。`.py` ファイルは変更していない。
  1. `[project.optional-dependencies].dev` を `["pytest","ruff","mypy","types-PyYAML","types-requests","import-linter","pyrefly","radon"]` にした。**`uv lock` は未実行**なので、lockfile の更新が必要。
  2. `[tool.ruff.lint.mccabe] max-complexity = 15` を追加した。C901 は select していないため、現時点では効果はない(段階導入用)。
  3. `[tool.importlinter]` を追加した(root_packages: hfauto / hfauto_ops / hfauto_viz。契約 11 本、`ignore_imports` による緩和なし)。
  4. `[tool.pyrefly]` を追加した(project-includes は 3 パッケージ、tests と runs は除外、`python-interpreter-path = ".venv-win-qa/Scripts/python.exe"`)。このパスはこの Windows 機専用で、Linux/HPC では `--python-interpreter-path` の指定が必要。指定しないと site-packages を解決できず、エラー数が 168 から 250 に増える。
  5. `[tool.radon] exclude = "runs/*,scripts/*"`。radon 6.0.1 は `[tool]` テーブル全体を configparser で読むため、次の 2 点に注意する: 真偽値キー(`show_complexity` など)を書くと CLI フラグが反転する。値に裸の `%` があると壊れる。この注意は pyproject 内のコメントにも残した。
- キャッシュディレクトリは作っていない(`--no-cache`、`-p no:cacheprovider`)。

### 1.2 実行方法(リポジトリのルートで実行。`S=.venv-win-qa/Scripts`)

```bash
$S/lint-imports --no-cache                                   # 契約検査
$S/pyrefly check --summary=full                              # 型検査(Linux では --python-interpreter-path を付与)
$S/radon cc -s -a -n C hfauto hfauto_ops hfauto_viz          # 循環的複雑度(C 以上)
$S/radon mi -s hfauto hfauto_ops hfauto_viz                  # 保守性指標
$S/ruff check .                                              # 既定ルール(clean)
$S/ruff check hfauto hfauto_ops hfauto_viz tests --select C90,PLR0911,PLR0912,PLR0913,PLR0915,PLR2004,B,SIM,RET,ARG,PERF,ERA,TRY,FBT,PLW --statistics
$S/ruff format --check .
$S/python -m pytest --collect-only -q -p no:cacheprovider --continue-on-collection-errors
```

全コマンドは `tools/commands.txt` にある(付録 5.3 のパス)。

### 1.3 主要メトリクス

#### import-linter(258 ファイル、890 依存。4 kept / 7 broken)

| 契約 | 種別 | 結果 | 違反 | 主な内容 |
|---|---|---|---|---|
| hfauto-layers | layers(cli > workflow > stages > backends\|hpc\|reporting > chemistry > core) | **BROKEN** | 23 | stages→workflow 22(ドメインサービス側)、hpc.compare→reporting.html_report 1 |
| backend-categories-independent | independence | **BROKEN** | 8 | backends.ts→backends.qm(nwchem_neb/string/saddle/path_support/pysisyphus/pysisyphus_saddle→qm.nwchem、orca_nebts→qm.orca、ts.dummy→qm.dummy) |
| concrete-backends-via-registry | protected | **BROKEN** | 7 | build_complexes/relaxation_discovery→crest、kinetics→cantera/arkane、thermo→qm.nwchem、path_ensemble_assess/recover_path→ts.nwchem_path_support |
| stages-independent | independence(41 stage) | **BROKEN** | 4 | irc/ts_search/recover_path→reaction_plan.make_reaction_case_artifact、rank→reaction_rank |
| core-is-foundation | forbidden | KEPT | 0 | — |
| hfauto-without-addons | forbidden | **BROKEN** | 2 | stages.ops→hfauto_ops.bundle、stages.viz→hfauto_viz(lazy) |
| addons-consume-artifacts | forbidden | KEPT | 0 | — |
| addons-independent | independence | KEPT | 0 | — |
| viz-layers | layers | **BROKEN** | 1 | hfauto_viz.core.run_loader→data.loaders(run_loader の importer は 0) |
| ops-layers | layers | KEPT | 0 | — |
| acyclic-siblings | acyclic_siblings | **BROKEN** | 2 循環 | hfauto: stages↔workflow(runner の 2 import)、hfauto_viz: data↔core |

- 静的グラフから見えないもの: stages/registry.py が 41 stage を文字列から遅延 import するため、stage の読み込み経路は import-linter から見えない。
- モジュール単位の循環は 0。パッケージ単位の強連結成分は {stages, workflow, hfauto_ops.bundle, hfauto_ops.core, hfauto_ops.qcarchive} と {viz.core, viz.data}。
- 推移的依存: backends.registry は具象バックエンド 27 個を eager import するため、registry を使う stage は約 66 モジュールを引き込む。build_complexes は 20、minimum_mode_assess は 22 モジュール。

#### pyrefly(255 モジュール、43,927 行、0.44 s)

| 種別 | 件数 | | パッケージ | 件数 |
|---|---|---|---|---|
| bad-argument-type | 75 | | hfauto.stages | 45 |
| missing-attribute | 28 | | backends.ts | 27 |
| missing-import(環境起因: scine 10、plotly 9、h5py 1) | 20 | | chemistry | 23 |
| unsupported-operation | 10 | | viz.renderers | 16 |
| bad-function-definition | 7 | | backends.db | 14 |
| unbound-name | 7 | | backends.reaction_discovery | 10 |
| not-iterable / no-matching-overload | 5 / 5 | | backends.qm | 8 |
| bad-assignment / bad-index / bad-return / bad-override | 4 / 3 / 3 / 1 | | 合計(hfauto / ops / viz) | 139 / 10 / 19 |

- エラーは **168 件**(抑制 7 件)、警告は 313 件。警告のうち 312 件は unnecessary-type-conversion(`str()`/`float()` の防御的な二重変換。多いのは proton_transfer.py 32 件、stages/thermo.py 23 件)。
- 主なパターン: `Artifact.data/qc` が異種の値を入れる辞書なので、union 型の推論(`Literal[True] | dict | str`)と、`dict.get()` 由来の `X | None` が `int()/float()/read_xyz()` へ流れることによるエラーが多い。

**pyrefly とレビューで確認した実バグ**

| ID | 場所 | 内容 | 重大度 |
|---|---|---|---|
| BUG-01 | `hfauto/backends/db/comptox.py:51-53` | `HTTPResult` にない `resp.status/resp.value` を参照しており、CompTox の API 応答のたびに AttributeError になる | 中 |
| BUG-02 | `hfauto_ops/schedulers/slurm.py:19`、`lsf.py:19` | Slurm の `submit_command` が str を返し、`command[0]=='s'` として解決されて必ず失敗する。LSF の `["bsub","<",script]` はシェルを通さないので `<` がリダイレクトにならない | 中 |
| BUG-03 | `hfauto/stages/thermo.py:620`、`thermo_sensitivity.py:161` | `get_thermo_engine("arkane")` が `species_thermo` を持たない ArkaneEngine を返す(設定しだいで AttributeError) | 低(潜在) |
| BUG-04 | `hfauto/backends/ts/nwchem_path_support.py:85,87` | `kabsch_rmsd` が None を返したときに `None > float` で TypeError | 低 |
| BUG-05 | `hfauto/backends/registry.py:39` | name/provider/engine のない dict spec が `KeyError: Unknown ... None` になる | 低 |
| BUG-06 | `hfauto/core/thermo.py:198,229` | `wigner_tunneling_factor(arg1,arg2)` と `equilibrium_constant_from_delta_g(*args)` が引数の順序を値の大きさで推測する。(T=298.15, ΔG=150) で ΔG と T を取り違える(このモジュールは死蔵) | 中(削除で解消) |
| BUG-07 | `hfauto/backends/qm/xtb.py:47,284` | xtb 6.7.1 は `--cycles` 超過時も rc=0 と `normal termination` を出力し、xtbopt.xyz も書く(検証者が実行して確認)。ラッパーはこれを収束扱いにする | 高(CH-20) |
| BUG-08 | `hfauto/backends/conformer/crest.py:133-147` | エネルギーが欠けた配座を 0.0 kcal/mol(最安定)として扱う | 中(CH-19) |
| BUG-09 | `hfauto/backends/reaction_discovery/readuct.py:373` | Hessian のフォールバックキー `'ts'` が存在しない(短絡評価のため到達はまれ) | 低 |

#### radon

- CC: 1,148 関数・メソッドの評価分布は A 695 / B 224 / C 150 / D 32 / E 20 / F 27。関数とメソッドの平均 CC は 7.63。CC≥21 が 79 件、CC≥41 が 27 件。
- パッケージ別(関数数・平均 CC・最大 CC): stages 153・14.97・186(F 11 件) / chemistry 193・9.78・115(F 7 件) / backends.ts 93・8.10・55 / workflow 30・10.17・64 / backends.qm 平均 6.95・最大 47 / core 4.41・34 / ops.core 5.9・17 / viz 2.4〜6.7
- MI(保守性指標)で C 評価のファイル: `backends/qm/nwchem.py` 0.00(SLOC 1,152)、`stages/dft_minima.py` 0.00、`stages/thermo.py` 0.00、`chemistry/reaction_classification.py` 3.21、`stages/build_complexes.py` 6.12、`backends/ts/pysisyphus.py` 8.59
- 規模: LOC 43,940 / SLOC 38,790。コメントは 230 行、docstring は 768 行(約 1〜2%)

**ワースト 10(CC)**

| # | CC | 場所 | 関数 | 行数 | 備考 |
|---|---|---|---|---|---|
| 1 | 186 | `stages/thermo.py:339` | ThermoStage.run | 829 | 8 責務が同居、実行テスト 0(AR-20) |
| 2 | 115 | `chemistry/reaction_classification.py:249` | classify_reaction | 415 | 14 引数(AR-22) |
| 3 | 77 | `chemistry/minimum_recovery.py:20` | assess_minimum_mode_following_source | 204 | NWChem 固有(AR-14) |
| 4 | 73 | `stages/ts_search.py:94` | TSSearchStage.run | 397 | 状態機械と実行(AR-06) |
| 5 | 66 | `chemistry/nwchem_evidence.py:205` | assess_endpoint_nwchem_evidence | 159 | **呼び出し元 0(死蔵)** |
| 6 | 64 | `workflow/discovery_coverage.py:55` | build_discovery_coverage | 232 | — |
| 7 | 58 | `chemistry/path_convergence.py:98` | assess_path_ensemble | 157 | — |
| 8 | 58 | `stages/rank.py:51` | RankStage.run | 266 | legacy、どのパイプラインでも未使用 |
| 9 | 56 | `stages/thermo.py:150` | ThermoStage._endpoint_frequency_adapter | 187 | — |
| 10 | 55 | `backends/ts/nwchem_string.py:230` | NWChemStringEngine.search_ts | 521 | NEB と重複(AR-12) |

(11 位以下: NWChemNEBEngine.search_ts 54/589 行、GoodVibesEngine.species_thermo 52、ExploreReactionsStage.run 51、MinimumModeAssessStage.run 51、DFTMinimaStage.run 48、NWChemEngine._run 47 …)

#### Ruff

- 既定の実行(`ruff check .`、ruff 0.16 既定の 410 ルール、309 ファイル): **All checks passed**。
- レビュー用ルール群では 1,058 件(本体 851 件、tests 207 件):

| ルール | 件数 | 意味 | | ルール | 件数 |
|---|---|---|---|---|---|
| PLR2004 | 295 | マジックナンバー(1e-12 が 34 件、許容値が散在) | | C901(>15) | 24 |
| TRY003 | 283 | 例外メッセージの直書き | | B905 | 19(zip の strict なし。うち 12 件は幾何/XYZ) |
| PLR0913 | 66 | 引数が多すぎる(最大 14) | | PERF401 | 15 |
| FBT001/002/003 | 62/47/39 | 真偽値の位置引数 | | PLR0911 | 12 |
| PLR0912 | 43 | 分岐が多すぎる | | ARG001/002/005 | 33/23/27(未使用の config/manifest 引数) |
| PLR0915 | 37 | 文が多すぎる(最大 237) | | SIM105 | 11 |

- パッケージ別: chemistry 273、tests 207、stages 160、core 77、backends.ts 71。
- `ruff format --check`: 327 ファイル中 236 ファイルに差分がある(+7,526/−5,444 行)。フォーマットを適用するなら、リファクタとは別コミットで一括適用するのが望ましい。

#### pytest(収集のみ)

- 263 件を収集、収集エラー 1 件: `tests/test_amine_hf_panel.py` が Windows のアプリケーション制御による rdkit DLL のブロックで ImportError になった。環境要因だが、rdkit 依存テストに `importorskip` がないため、別の実行形態ではスイート全体が 0 件で中断する(AR-30)。
- `@pytest.fixture` は 0 件、`conftest.py` は 6 行、`monkeypatch` は 99 回。tests/`__pycache__` には、ソースが削除された 36 モジュール(`test_phase*` など)の pyc が残っている。

---

## 2. アーキテクチャレビュー

### 2.1 レイヤリングと依存

**評価**: 基盤の依存方向は保たれている。core は上位層を import せず、chemistry の依存先は core だけ(35 エッジ)、viz は core.io と core.schemas しか読まない。ただし、`workflow` パッケージ名が実際の責務と合っていないこと、add-on 境界が双方向になっていることの 2 つが構造上の課題である。

**AR-01【高】workflow が「オーケストレータ」と「stage 用ドメインサービス」を兼ね、stages↔workflow の循環を生んでいる**(元: A1-1 確認、A2-9/A3-8 一部訂正)
- 場所: `hfauto/workflow/runner.py:7-9`、`workflow/{reaction_state,reaction_evidence,reaction_inputs,reaction_artifacts,ts_execution,discovery_coverage}.py`、`stages/reaction_plan.py:23`(`make_reaction_case_artifact`)、`stages/irc.py:17`、`ts_search.py:24`、`recover_path.py:26`、`rank.py:20,79-82`
- 根拠: stages→workflow の 22 import はすべてドメインサービス向けで、逆向きは runner の 2 import だけ。workflow 1,256 行のうち runner は 94 行。stage 間の違反 4 件のうち 3 件は、共有ファクトリが stage モジュール内にあることが原因。`rank` は `ReactionRankStage(...).run()` を内部で呼ぶ「stage の中の stage」になっている。runner は `reporting.html_report.latest_manifest_path` にも依存している。
- 訂正: A2/A3 の検証では、既存 docs(`developer_architecture.md`)が workflow を「cross-stage services」と説明していることから、これは「層逆転」ではなく命名と配置の問題と判断している。重大度は A1 の確認値(高)を採った。
- 推奨: `runner.py` を `hfauto/pipeline/`(execute_stage、run_layout もここ)へ移す。残りは `hfauto/reaction_case/`(state/evidence/inputs/artifacts/ts_attempts/coverage)に改名する。`make_reaction_case_artifact` は reaction_case/artifacts.py へ移す。rank の委譲は削除する。層契約は `cli > pipeline > stages > reaction_case > backends|hpc|reporting > chemistry > core` とする。

**AR-02【中】add-on 境界が双方向で、配布物をまたぐ循環がある**(A1-7 確認)
- 場所: `stages/ops.py:8`、`stages/viz.py:26`、`hfauto_ops/core/array_jobs.py:31,271`、`hfauto_ops/retry/policy.py:9`、`pyproject.toml:72`(3 パッケージを 1 つの wheel に同梱)
- 根拠: hfauto→add-on が 2 import、add-on→hfauto.stages が 3 import。viz の import は lazy で、失敗時は FAILURE artifact を返すので結合は弱い。ops 方向の結合のほうが強い。
- 推奨: OpsStage/VizStage を各 add-on へ移し、`[project.entry-points."hfauto.stages"]` で登録する。ops は `hfauto.pipeline` の API だけを使う。`hfauto_viz/core/run_loader.py`(importer 0)を削除する。forbidden 契約「hfauto_ops → hfauto.stages / hfauto.workflow」を追加する。

**良い点**
- core-is-foundation / addons-consume-artifacts / addons-independent / ops-layers の 4 契約は守られている。モジュール単位の循環は 0。
- 重い任意依存(scine、rdkit、plotly)はバックエンドや描画器の関数内 import に閉じ込められていて、モジュールレベルで import されるものは 0 件。

### 2.2 計算単位の独立性

**評価**: 外側の Stage インターフェースは統一されていて、`hfauto run-stage --in manifest.json` で単独実行でき、22 クラスがテストから直接起動されている。一方で、計算時間と収束問題の大半を占める TS/IRC/thermo 系ほど結合が強い。そのため、収束改善のために NEB 1 回、saddle 1 回、IRC 1 回といった単位だけを差し替えたり単体で検証したりすることが難しい。

**stage 別評価表**(行数と最大 CC は radon による。「使用」は現行 `configs/pipelines/*.yaml` 13 本での使用有無)

| stage | 行数 | 最大CC | 判定 | 主な結合要因 | 使用 | 方針 |
|---|---:|---:|---|---|:-:|---|
| ingest | 216 | 41 | 独立 | — | ○ | run の分割(CC) |
| enrich | 94 | 18 | 独立 | DB provider(未テスト) | × | 要判断(公開 DB の位置づけ) |
| enumerate-states | 117 | 11 | 独立 | — | ○ | 維持 |
| detect-sites | 63 | 10 | 独立 | 出力を build-complexes が使っていない(元素順で選択) | × | 配線するか削除(CH-23) |
| conformers | 49 | 9 | 独立(理想形) | — | ○ | 手本 |
| build-complexes | 837 | 32 | 部分 | `CRESTConformerBackend()` を直接生成、配置アルゴリズムを内包 | ○ | registry 経由化、chemistry へ移設 |
| generate-reactions | 206 | 30 | 部分 | 独自の checkpoint | ○ | 共通 JobStore |
| explore-reactions | 635 | 51 | 部分 | 独自の checkpoint、端点選択アルゴリズム(約 170 行)を内包 | ○ | 同上 |
| discovery-audit | 135 | 16 | 部分 | workflow.discovery_coverage、他 run の manifest を読む | ○ | 維持 |
| minimum-registry | 217 | 16 | 独立 | chemistry の呼び出しのみ | ○ | 維持 |
| connect-minima | 231 | 29 | 独立 | 閾値 0.05 Å(CH-07) | ○ | 閾値の見直し |
| preopt | 405 | 45 | 部分 | 独自の再利用判定、`method.stage='preopt'` を他の stage が参照 | ○ | JobStore、--ohess(CH-20) |
| relaxation-discovery | 375 | 19 | 部分 | CREST を直接 import、失敗 artifact を再解析 | ○ | typed outcome |
| dft-minima | 877 | 48 | 部分 | NWChem 入力を stage 内で再パース(l.420-440)、stage 名を参照される | ○ | method_evidence 化(AR-14) |
| minimum-mode-follow | 357 | 20 | 部分 | `method.stage=='dft-minima'` で絞り込み | × | **配線する**(化学的に有効) |
| minimum-mode-assess | 471 | 51 | 部分 | 同上 | × | **配線する** |
| endpoint-seeds | 290 | 48 | 部分 | 入力 `endpoint_pair_selection` の生成元がない | × | 要判断(生成元を実装するか、入力を変更) |
| endpoint-seed-screen | 152 | 27 | 部分 | `method.stage=='preopt'` に依存 | × | 虚振動による除外を追加して配線 |
| reaction-plan | 435 | 20 | **結合** | 共有ファクトリを保持、`basin_assessment` のキー欠落(CH-05) | ○ | reaction_case へ移設 |
| recover-path | 276 | 15 | **結合** | reaction_plan、ts.nwchem_path_support、fail-open の既定値 | ○ | 修正(CH-06) |
| ts-search | 490 | 73 | **結合** | 状態機械、ハードコードされたエンジン表、workflow の 3 モジュール | ○ | 分割(AR-06) |
| path-ensemble | 364 | 39 | 部分 | 偽収束した経路を「独立経路」として数える | × | gmax 判定を直してから配線 |
| path-ensemble-assess | 316 | 24 | **結合** | `engine=='nwchem'` 分岐で入力を監査(エンジン間で非対称) | × | 同上 |
| path-intermediates | 146 | 30 | 部分 | workflow ヘルパー | × | 配線する |
| irc | 318 | 39 | **結合** | reaction_plan、backend の qc を AND で上書き | ○ | 分岐下降の判定(CH-03) |
| reaction-classify | 189 | 34 | 部分 | workflow ヘルパー | ○ | 維持 |
| reaction-segments | 264 | 38 | 部分 | stoichiometry のない reaction を出し、thermo の legacy 分岐へ落ちる | × | stoichiometry を必須にしてから配線 |
| sp | 101 | 8 | 独立(理想形) | — | ○ | 手本 |
| method-panel | 354 | 43 | 独立 | 鍵が method_id だけ(CH-25) | ○ | 完全な fingerprint を鍵に |
| basin-populations | 195 | 22 | 独立 | — | ○ | 維持 |
| thermo | 1167 | 186 | **結合** | NWChem backfill、stage 名フィルタ、legacy HF 分岐 | ○ | 4 分割(AR-20) |
| thermo-sensitivity | 273 | 39 | 部分 | goodvibes 固定、別の source 選択規則 | ○ | thermo_sources を共用 |
| descriptors | 60 | 13 | 独立 | hf_stretch の定義が誤り | × | 削除 |
| kinetics | 124 | 14 | **結合** | cantera/arkane を直接生成、エンジン名の読み替え | × | 削除(CH-43) |
| calibrate | 274 | 36 | 独立 | config 未使用、DB 支持スコア | × | 削除し、誤差評価で置き換え(CH-44) |
| connector-audit | 157 | 36 | 独立 | config 未使用 | × | 削除候補 |
| rank | 316 | 58 | **結合** | ReactionRankStage を内部で起動、legacy HF | × | 削除(reaction-rank の別名にする) |
| reaction-rank | 328 | 29 | 独立 | — | ○ | 維持 |
| viz | 91 | 8 | **結合** | hfauto_viz、`out_dir.parent` | × | add-on へ移設 |
| hpc-plan | 91 | 6 | **結合** | hfauto.hpc、`out_dir.parent` | × | 削除または ops へ統合 |
| ops | 39 | 9 | **結合** | hfauto_ops.bundle | × | add-on へ移設または削除 |

集計: 独立 14・部分 16・結合 11(A2-T の集計 13/16/12 を上表のとおり再分類)。現行パイプラインで使われているのは 23 stage、未使用は 18 stage。

**AR-21【高】科学ゲートが複数実装されていて、段によって判定基準が違う**(A4-6 確認)
- 場所: `stages/preopt.py:83`(ラベルに依存しない `hf_endpoint_metrics_from_xyz`)と `stages/dft_minima.py:148`(ラベル固定の `proton_transfer_metrics`、スペクテーター HF を無視)、`chemistry/geometry_qc.py:113-125`(絶対距離の閾値)。極小の受理判定は `minima.is_accepted_optimized_minimum`、`core/qc.minimum_promotion_gate`、`thermo._is_scientific_calc`、`dft_minima._reusable_minimum_calculation` の 4 実装。lineage は 6 箇所で再評価されている。
- 影響: (HF)n(n≥2)で等価な H が入れ替わると、preopt は受理し、dft-minima は棄却する。この食い違いは現行の m3 設定でも起こりうる。
- 推奨: `chemistry/gates/` に `proton_state()`、`minimum_evidence()`、`method_lineage()` を 1 つずつ置く。結果は gate_name・gate_version・入力ハッシュ付きで 1 回だけ記録し、下流はハッシュ照合だけにする。

**良い点**: conformers(49 行)と sp(101 行)は「registry からエンジンを得る → 入力を選ぶ → backend を呼ぶ → 記録する」だけの理想形。stage 内での subprocess 直接呼び出しは 0 件。失敗は `ArtifactStatus` として一級のデータになっている。

### 2.3 連携計算(オーケストレーション)の有用性

**評価**: 判断ロジック(`decide_reaction_case` の決定論的な純関数、同一戦略の再試行禁止、予算管理、`ts_attempt_key` による同一入力リトライの防止)は科学的に筋が良い。しかし、その判断を実行につなぐ仕組み(ディスパッチ、再開、並列化、HPC)が欠けているか壊れていて、再計画ループを人手で閉じている。

**AR-05【高】状態機械の `next_action` を実行するディスパッチャがなく、scripts/ による手動連鎖とゲートの迂回が起きている**(A3-1 一部訂正、C4 で見落とされた点を追加)
- 場所: `hfauto/workflow/reaction_state.py:52-257`、`stages/ts_search.py:185-495,453-465`、`scripts/run_nwchem_saddle.py:151-172`、`scripts/run_pysisyphus_irc.py`、`scripts/merge_manifests.py`
- 根拠: 10 種ある next_action のうち、TS 以外の `run_fixed_geometry_frequency / validate_intermediate_basins / recompute_endpoints_on_one_pes / revise_reaction_hypothesis / audit_basins` を処理する stage は 0。`run_fixed_geometry_frequency` になると、ts-search は failure も出さずに保留し、resume しても `reaction_case_not_ts_ready` で拒否される(行き止まり)。runs 137 件のうち runner 経由は 68 件、トップレベル manifest だけのものが 45 件、manifest なしが 24 件。m3_tma_* は 54 run。現行 registry にない stage 名が 19 種・42 manifest ある。`scripts/`(git 未追跡)は stage のゲートを通らない。`run_pysisyphus_irc.py` は `hcn_irc_v5` の偽陽性 `irc_validated=True` を作り、`run_nwchem_saddle.py` は `saddle_seed_evidence_validated=True` と `optimization_convergence='tight'` を固定値で入れる。
- 訂正: TS 探索の中(NEB→string→saddle_refinement→Hessian 更新による再開)は、ts-search 内の `while next_authorized_attempt` で自動的にディスパッチされている。「全部が人手」ではない。
- 推奨: `ReactionCaseDriver` を作り、next_action を WorkItem 種別に対応づける表を持たせる(freq、irc、dft-minima、segments)。そのうえで `while not complete and budget: dispatch → decide` のループを回す(AiiDA の WorkChain `while_` や jobflow の動的 Response と同じ構造)。scripts/run_* は stage か WorkItem に取り込んでから削除する。

**AR-06【中】状態機械が 4 stage に分散し、エンジン表が二重管理になっていて、語彙も一致しない**(A2-2 一部訂正、A3 で見落とされた点を追加)
- 場所: `workflow/reaction_state.py:7-28,136-160,213-216`、`workflow/ts_execution.py:41-52`、`core/schemas/path.py:15`、`backends/ts/nwchem_saddle.py:403-420`
- 根拠: `_ENGINES_FOR_STRATEGY` にエンジン名が直書きされ、別に `engine_order` があり、両者を事前に照合する仕組みがない。`endpoint_local_saddle_search`(pysisyphus)はどの diagnosis からも到達できない。orca_nebts と pysisyphus_saddle はどの engine_order にも含まれない。`final_frequency_missing` は先に return されるため死に項目になっている。`decide_reaction_case` の呼び出し元 4 つで、`max_path_attempts` の既定値が 4 と 3 に分かれている。PathDiagnosis の Literal は 6 値なのに、実際の語彙は 12 キーある。
- 訂正: `hono_isomerization_v1` の「同じ saddle_refinement の 2 回反復」は、現行コードより前の挙動だった。しかも 2 回とも seed 検証で即座に失敗しており、NWChem は実行されていない。現行コードで再現すると、`adaptive_double_ended_path` へ再計画される。
- 推奨: `Diagnosis`/`Strategy` を StrEnum にする。エンジンにはクラス属性 `strategies: frozenset[Strategy]` を宣言させ、registry が対応表を生成する。パイプライン開始時の preflight で、到達しうる全 Strategy に対応するエンジンがあることを検証する。予算チェックと非反復チェックは全分岐に先立って適用する。

**AR-03【高】stage の実行経路が 4 系統あり、HPC/配列実行の経路は end-to-end で動かない**(A1-2 確認、A3-2 一部訂正、BUG-02)
- 場所: `workflow/runner.py:72-88`、`cli/main.py:50,56-78`、`hfauto_ops/core/array_jobs.py:192-197,244-278`、`hfauto_ops/schedulers/production.py:80-88`、`hfauto/hpc/job_plan.py:45-51`、`hfauto/hpc/schedulers.py:80-143`
- 根拠: `StageContext(` を組み立てる箇所が 4 つある。ops の `run_array_task` は `global_config={}` で起動するため、production ゲートが外れる。生成されるスクリプトは `hfauto_ops.cli.main run-array-task` を呼ぶが、このコマンドは存在しない(実行して "No such command" を確認。`runs/foundation_smoke/16_ops/` に生成物が残っている)。Snakefile の 17 rule には依存辺がない。YAML に `stage_config` がないため、`hpc-plan` 経由のジョブは engine='dummy' で実行される。retry は未定義の `$HFAUTO_INPUT_MANIFEST` を参照する。ts-search のスライスは reaction_case を落とす。配列ヘッダは一律 1 CPU/4 GB/2 h。
- 訂正: `hfauto run-stage` は入力 manifest の `global_config` を継承するので、hpc-plan 経路では production ゲートは維持される(dummy の結果は失敗 artifact になる)。runs にこの経路の実行実績は 0 件で、「宣伝されている未使用機能が壊れている」という性質の問題である。
- 推奨: `hfauto/pipeline/execute.py` に `execute_stage(...)` を 1 つだけ置き、全経路がこれを通るようにする。global_config は「明示指定 → manifest.metadata → {}」の順で決める。配列経路は、Executor Protocol(Local=ProcessPool、Slurm=`--array`+`--dependency=afterok`)で作り直すか、当面は削除する。生成スクリプトに `bash -n` と CLI 解決の smoke テストを付ける。

**AR-04【中】hfauto.hpc と hfauto_ops が三重実装になっていて、資源表が矛盾し、死蔵モジュールも多い**(A1-9、A3-3 確認)
- 根拠: ts-search の資源既定値が 4 コア/4 GB/48 h(hpc)、16 コア/48 GB/2,880 分(ops core)、24 スレッド/64 GB/24 h(ops scheduler)と 3 通りある。retry 立案が 3 実装、artifact→stage の対応表が 4 つ(`species` の行き先が食い違う)、compare_runs が 3 実装、scheduler/ と schedulers/ の二重ディレクトリ。本番からの importer が 0 の ops モジュールは 8 本・約 1,010 行。retry/policy の 10 カテゴリのうち 7 つは hfauto 内で一度も emit されない。`configs/stages` 10 本、`configs/hpc` 3 本、`configs/ops`、`configs/viz` の計 15 ファイルは、どこからも参照されない。
- 推奨: ops の居場所を hfauto_ops に一本化する。`hfauto/hpc/*`、hpc-plan/ops stage、hfauto CLI の hpc 系コマンドは統合または削除する。資源表はサイトプロファイル 1 箇所にまとめる(4.5 節)。

**AR-07【高】部分再実行と manifest ポインタの意味が危うい**(A3-4 確認、A1-8 確認)
- 場所: `workflow/runner.py:12-14,47-61,92-93`、`reporting/html_report.py:12-21`、`hfauto_viz/core/paths.py:14-30`、`stages/ts_search.py:69-89,187`、`cli/main.py:253-273`
- 根拠:
  - `--from X` の入力は「前回の最終 stage」になり、carry_forward で古い下流の成果物(irc/thermo/rank)が残る。
  - ts-search の resume は reaction_id と qc しか照合しない。runner 経由の run の 53/55 が production(resume 既定で有効)なので、汎関数や基底を変えても旧 TS が「完了済み」としてスキップされうる。
  - `manifest.path` は全 stage が成功した後にしか書かれない。stage 番号は有効 stage の順序で決まるので、YAML を編集すると衝突する(`hcn_goodvibes_v4` に `12_` が 2 つある)。
  - パスは CWD 相対で、`latest_manifest_path` は reporting と viz で挙動が違う。リポジトリ外から呼ぶと、reporting 版は存在しないパスを返す(実測)。コード版の記録は 0。
- 推奨: `core/run_layout.py` に RunLayout を置く。stage ディレクトリは YAML の安定した `id:` で命名する。stage ごとの状態を `run_state.json` に逐次書く(入力 sha、resolved_config sha、code_version)。上流を再実行したら下流を stale にする。パスは run_dir 相対にする。resume キーは全 stage 共通の request fingerprint にする。

**AR-08【高】長時間ジョブの保全がない(opt+freq を単一ジョブにまとめ、restart がなく、初期構造から再実行する)**(A6-4 確認、A3-6 一部訂正、A3-5/A2-5/A4-7 確認)
- 場所: `backends/qm/nwchem.py:782,1133-1156`、`stages/dft_minima.py:75-146,783-797`、`stages/ts_search.py:43-66`、`stages/preopt.py:42-80`、`core/executables.py:162-170`
- 根拠:
  - command_result 428 件、124.1 h のうち、失敗・タイムアウト 31 件で 38.0 h(31%)。タイムアウト 11 件はすべて NWChem mpirun。
  - `spc_TMA_HF3_seed00_nci_0` は 38 step で最適化が収束した後、CPHF の途中(14,395 s)でタイムアウトし、ジョブ全体が failure になった。
  - NWChem は常に `start hfauto_job` で起動し、`.movecs/.hess/.gbw` の再利用は 0。
  - checkpoint と再開は 5 stage が個別に実装しており(ディレクトリ命名は 3 種、fingerprint も別々)、実際に発火したのは全 run で 2 回だけ。利用者は run_id を v3→v11 のように付け替えて回避している。
- 訂正: TS の saddle については `plan_saddle_recovery` が最終構造から trust を縮めて再開する上限付きの手順を持っていて、その範囲では再計画の段階が実装されている。失敗時間 33.1 h の約半分は、現行 registry にない旧 stage(proton-transfer-endpoints 系)が占めていた。
- 推奨: optimize と frequency を 2 ジョブに分ける(注: 現行 render は `task=optimize` にも freq ブロックを付けているので、先にこれを外す必要がある)。NWChem の restart(`restart` と `vectors input`、最終 final-NNN.xyz、`inhess 0`)を backend の責務にする。共通の JobStore(内容アドレス型キー = engine + 正規化 method + 入力ハッシュ)を作り、attempt ladder を FailureKind をキーに宣言的に持たせる。`run_command` に進捗監視コールバックを付け、SCF 失敗や NEB 停滞で早期に打ち切れるようにする。

**AR-09【中】並列化の手段がない**(A3-7 一部訂正)
- 根拠: concurrent.futures / multiprocessing / asyncio の使用は 0 件。ts-search は反応×戦略の直列ループ、dft-minima は化学種の直列ループで、15.7 h の dft-minima は 8 並列なら約 4 h で済む計算量だった。(訂正: `max_reactions: 1` は直列性の回避策ではなく、予算と resume の設計である)
- 推奨: Stage を `plan() → list[WorkItem]`、`execute(item)`、`reduce()` に分ける。Local は ProcessPool(`site.cores // item.ncores`)、HPC は 1 item を 1 array task にする。これは CH-13 のグローバル rank 予算と組み合わせることが前提になる。

**AR-10【中】設定の層分けがなく、再現性が低い**(A3-10 一部訂正、見落とし追加)
- 根拠: 13 パイプラインに `/home/user`・`/mnt/c` の絶対パスが 42 箇所ある。PBE0/def2-SVPD/D3/xfine の定義が 10〜11 ファイルにコピーされている。runs に記録された pipeline_config_path 21 種のうち 8 種は、もう存在しない。解決済み設定の本体も保存されない。`run-stage` は production の preflight を通らない。`resolve_executable` は YAML の明示パスを env より優先するので、`HFAUTO_*_EXECUTABLE` でサイト差を吸収できない。
- 推奨: 設定を pipeline(科学的なフロー)/ `configs/methods/*.yaml`(method profile)/ `configs/sites/<site>.yaml`(実行ファイル・MPI・scratch・資源・cache_root)の 3 層に分け、`--site` で切り替える。runner は `resolved_config.yaml` と code_version を書き出す。

**良い点**: `decide_reaction_case` は 24 テストで固定された純関数。`ts_attempt_key` は同一入力の無限リトライを構造的に防いでいる。予算超過は「反応が存在しない証拠ではない」`execution_budget_exhausted` として区別して記録される。dft-minima の再利用判定(入力ハッシュ、returncode、timed_out、fallback を照合)は、統一キャッシュを作るときの仕様としてそのまま使える。`Manifest.merge` は fan-out/gather の正しいプリミティブ。低レベル探索と DFT 追跡をパイプラインとして分けている運用も妥当である。

### 2.4 責務と複雑度

**評価**: 複雑度は「物理計算そのもの」ではなく、「証拠の辻褄合わせ、互換性の維持、ソース選択」に集中している。ホットスポットを「純関数カーネル + 薄いアダプタ + 薄い stage」に分解すれば、上位 10 関数の多くは CC 15 以下にできる。

**ホットスポットと分割案**

| 対象 | CC / 行 | 同居している責務 | 分割案 |
|---|---|---|---|
| `stages/thermo.py` ThermoStage.run | 186 / 829(MI 0) | 化学種の索引化、周波数系譜のアダプタ、NWChem スキーマ移行、ソース選択(E/熱補正)、分圧ロール、エンジン呼び出し、ΔG/K の 2 系統(化学量論/legacy HF)、科学ゲート、表出力 | `chemistry/thermo_sources.py`(検証済み (E, freq) の組を返す)/ `backends/thermo/goodvibes.py`(型付き設定)/ `chemistry/reaction_thermo.py`(ΔE, ΔE0, ΔH, ΔG, K)/ stage は I/O だけ。legacy は削除、移行は `hfauto migrate-manifest` へ |
| `chemistry/reaction_classification.py` classify_reaction | 115 / 415、14 引数 | 分類、上流検証の再実行(`validate_path_attempt_evidence`、lineage) | 型付きの `ReactionEvidence` を 1 引数で受け取り、`_classify_elementary/_stepwise/_unresolved` に分割(各 60 行以下) |
| `backends/ts/nwchem_{neb,string,saddle}.py`、`pysisyphus*.py` の search_ts | 41〜55 / 349〜589 | 入力生成、実行、収束判定、解析、経路 QC、推定 TS の選択、saddle 精密化、振動検証、`reaction_validated` の組み立て | `PathGenerator`(経路と HEI だけ)/ `SaddleRefiner`(QM を注入)/ `chemistry/ts_validation.py` / `ts/artifacts.py`(AR-12) |
| `backends/qm/nwchem.py` | ファイル MI 0、SLOC 1,152、C 以上 7 関数 | 入力生成、出力解析、実行、ダミー代替、スキーマ移行、HF 固有の後処理 | `backends/qm/nwchem/{input,output,runner,migrations,engine}.py`、HF 後処理は `chemistry/post_qc.py` |
| `chemistry/minimum_recovery.py` assess_minimum_mode_following_source | 77 / 204 | NWChem 入力パーサへの依存、`nwchem_real` 要求 | 正規化された method_evidence の比較に置き換える |
| `stages/dft_minima.py` DFTMinimaStage.run | 48 / 268(MI 0) | 候補選択 約 180 行、NWChem 入力監査、checkpoint | 選択は `chemistry/endpoint_selection.py`、再利用は JobStore |
| `stages/explore_reactions.py` | 51 / 387 | 端点のスコアリングと選択(約 170 行)、checkpoint、予算 | 同上 |
| `workflow/discovery_coverage.py` build_discovery_coverage | 64 / 232 | 被覆率の簿記 | 機構別の実行率・固有生成物数・陰性結果の分類に絞る(CH-29) |

**AR-20【高】ThermoStage.run の責務過多。legacy HF 分岐は誤ラベルの記録を生む潜在バグで、実行テストは 0**(A4-1 一部訂正、C6-9、A2-7、A5-3 確認)
- 場所: `stages/thermo.py:8,150,339,423-668,671-958,960-1157`、`stages/rank.py:51-316`、`tests/test_thermo_scientific_tasks.py`(定数 1 つを確認するだけ)
- 根拠: 化学量論版と legacy 版で食い違いがある(`ThermoRecord` の検証は片方だけ、`real_irc_executed` の判定式が別)。互換エイリアスが重複している(`delta_G_act_kcal_mol` と `delta_G_activation_standard_kcal_mol`、`main_values_are_fallback` と `main_values_are_dummy`)。
- 訂正: legacy 分岐は「到達不能」ではない。起動条件は `stoichiometry` が空であることだけで(l.672)、reaction-segments が出力する path_segment 反応(stoichiometry なし)を thermo に流すと、HF 固有名(`delta_G_ion`、`dg_act`)の低 tier 記録が生成される。現行 13 本の設定にはこの連結がないので、実質的には死蔵。
- 推奨: 上表の 4 分割。legacy 分岐と rank の旧スコアリングを削除し、stoichiometry のない reaction は `missing_stoichiometry` の failure にする。最小フィクスチャ(1 反応 + TS)によるゴールデンテストを追加する。

**AR-22【中】chemistry 層で純粋カーネルと証拠ゲート・I/O・プログラム固有処理が混在している**(A4-9 確認)
- 根拠: chemistry 34 モジュールのうち 17 がファイル I/O を行い、10 が Artifact を引数に取る。Ruff の指摘は chemistry が最多で 273 件(TRY003 が 115、PLR2004 が 98。許容値 1e-12 などが名前付き定数になっていない)。
- 推奨: `chemistry/kernels/`(XYZ と ndarray を受けて数値を返す)と `chemistry/evidence/` に分け、import-linter で「kernels は core.schemas.artifact / core.io を import しない」を強制する。許容値は `chemistry/tolerances.py` にまとめる。

**AR-18【中】QM バックエンドにプログラム I/O・解析・ダミー代替・HF 固有後処理・移行が同居している**(A4-5 一部訂正)
- 根拠: ダミー代替が 3 通りに分岐している(`nwchem._fallback`、`orca._dummy_fallback`、XTBEngine が DummyQMEngine を継承し代替ブロックを 2 回複製)。xtb だけが `subprocess.run` を直接呼んでいて(xtb.py:184)、`run_command` によるプロセスグループの kill と command_result の記録が揃わない。
- 訂正: `geometry_sane` の HF 閾値は HF を含む化学種にしか適用されない。xtb はスレッド並列なので、孤児プロセスのリスクは低い。
- 推奨: バックエンドを `render→run→parse→CalculationRecord` に限定する。geometry_qc と記述子は stage 側の後処理プラグインにする。xtb も `run_command` 経由にする。

**良い点**: 外部プログラムの起動は `core/executables.run_command`(逐次書き出し、プロセスツリーの kill、timeout の構造化、command_result.json)にほぼ集約されている。XYZ 入出力、Kabsch、共有結合の隣接判定、置換不変 RMSD、単位定数は hfauto 本体内で一元化されていて、重複はほとんどない。

### 2.5 契約・診断・テストの過剰性

**評価**: 「契約」は型ではなく辞書の約束事になっている。書かれる診断の量に比べて、機械が消費する情報は少ない。検証が層ごとに繰り返され、テストはその辞書を再現する作業が中心になっている。一方で、物理的に決定的な検査(IRC 分岐の下降、座標系、振動ブロックの整合)は欠けている。**量が多いのに要所が抜けている**、という状態である。

**AR-23【高】Artifact.qc/data が型のない辞書で、書かれるキーの約 4 割はどこからも読まれない。Pydantic スキーマも半分は使われていない**(A5-1、A4-2 確認)
- 根拠:
  - qc キーは 260 種が書かれ、103 種(40%)は書き込み箇所以外で一度も現れない。data キーは 248 種中 74 種が同様。読み取りは 47 キー・262 回で、上位 28 キーが 90.8% を占める。
  - 実 run の 1 manifest に qc キーが 169 種、data キーが 401 種ある。
  - `CalculationRecord/SpeciesRecord/SiteRecord/DescriptorRecord/ReactionRecord` はスキーマ外からの利用が 0。`ThermoRecord` は `extra="allow"` で、検証の直後に dict へ戻している。
  - 同じ事実を重複して持っている: `n_imag` は data と qc の両方にあり、`thermo.py:228` は両方が 0 であることを要求する(片方しか書かないバックエンドは黙って不採用になる)。xtb は geometry QC を 3 重に格納している。
- 推奨: `core/evidence.py` に frozen の `CalculationEvidence` を置く(engine, program_version, real_qm_executed, fallback_dummy, scf/geometry_converged, normal_termination, n_imag, imag_freq_cm1, frequency_count_complete, minimum_accepted, ts_validated_by_frequency, method_evidence_validated, dispersion_applied, irc_validated の約 15 フィールド)。qc は「Evidence の dump」と「非契約の `qc["diagnostics"]`」に分ける。未使用の 5 レコードクラスは削除する。

**AR-24【中】ゲート判定式の複製と欠損キーの扱いの揺れ**(A5-2 一部訂正、A2-3 一部訂正)
- 根拠: `fallback_dummy` の読み取りが 26 ファイルにあり、厳格(`is False`)が 15 箇所、寛容が 34 箇所。検証済み TS の判定式は `ts_search.py:78-89`、`ts_execution.py:130-144`、`reaction_classification.py:98-112` にコピーされている。IRC の判定は 3 系統(`method_evidence_validated` 必須 / `real_irc_executed` のみ / `real_orca_executed` も受理)。ops と viz は `real_orca_executed` を既定 True(または `is False`)で読むので、欠損は合格扱いになる。
- 訂正: core/qc と thermo は `real_qm_executed is True` を同時に要求し、NWChem は両キーを常に書く。したがって現行経路でダミーが本物扱いされた実例はない(潜在的な問題)。`real_orca_executed` は遺物ではなく、ORCA バックエンドが現在も書いている。
- 推奨: `is_real_calculation / is_validated_minimum / is_validated_ts / is_validated_irc` の 4 つの純関数を evidence.py にまとめ、欠損は常に不合格にする。qc キーを直接判定しているコードが残っていないかを CI の grep で監視する。

**AR-25【中】stage 契約(STAGE_CONTRACTS)が強制されず実装とずれていて、stage メタデータが 6 箇所に分散している**(A1-11、A2-4、A3-9、A4-8、A5-7 確認)
- 根拠: `core/artifact_types.py:51` は「intentionally not enforced at runtime」と明記している。STAGE_CONTRACTS には 4 stage が欠けていて、コード中の artifact_type のうち 22 種が定数表にない。唯一の利用者は死蔵の `hfauto_ops/scheduler/resources.py`。thermo の宣言(in: calculation, reaction_path_validated)と、実際に読むもの(species, species_preopt, species_optimized, irc, reaction_validated, reaction)がずれている。`STAGE_OUTPUT_HINTS` は中身が入力型になっている。`endpoint-seeds` の入力 `endpoint_pair_selection` には生成元がない(過去の run にはある。リファクタで消えた)。
- 推奨: `StageSpec(consumes, produces, config_model, external_requirements, default_resources)` を Stage のクラス変数にする。runner がまず warn モードで検証し、`stage_io_reference.md` は宣言から自動生成する。STAGE_CONTRACTS、STAGE_OUTPUT_HINTS、STAGE_TARGET_TYPES、job_plan の対応表は削除する。

**AR-26【中】失敗分類と状態が自由文字列で、retry が部分一致で誤った stage に振り分ける**(A5-5、A6-6 確認)
- 根拠: failure category は 40〜60 種あり、「実行されなかった」だけで 10 通りの綴りがある。NWChem の timeout、SCF 失敗、maxiter 到達はすべて `nwchem_failed` にまとめられている。`infer_retry_stage` は 60 カテゴリ中 44 を unknown にし、`xtb_opt_failed` を dft-minima へ、`dimer_saddle_frequency_rejected` を dft-minima へ振り分ける。`build_retry_plan` の呼び出し元は 0。`"partial"` を判定するコードも 0(amine pilot では 59 attempt 中 53 が partial)。
- 推奨: `status: Literal["success","failed","rejected","skipped"]` と、`FailureKind(StrEnum)` の約 12 種(DISABLED, EXECUTABLE_MISSING, INPUT_INVALID, TIMEOUT, SCF_NOT_CONVERGED, GEOMETRY_MAXITER, WRONG_STATIONARY_POINT, INCOMPLETE_HESSIAN, VERSION_MISMATCH, PARSE_ERROR, STAGNATED, BUDGET_EXHAUSTED)に集約する。失敗を作る側が `retry_stage` を明示し、部分一致による推論は削除する。

**AR-27【中】熱化学レコードに重複フラグと擬似的な confidence_score がある**(A5-4 確認、見落とし追加)
- 根拠: `confidence_score` は 5 引数のうち 3 つが定数で、値は 0.6575/0.7325/0.745/0.82 の 4 通りしかない。`activation_connectivity_validated` は `has_validated_irc` と常に同じ値、`main_values_are_fallback` は `main_values_are_dummy` の別名で、内部フォールバック熱化学でも False になる(CH-39)。IRC の `reactant/product_endpoint_match` も常に同じ値。信頼度は ops の `run_index.py:111` 経由で運用側に表示されている。
- 推奨: `evidence_tier: Literal[...]` と `blockers: list[str]` の 2 フィールドにまとめ、資格判定は 1 つの関数にする。

**AR-28【低】定数の免責フラグ 63 種、next_action 34 種(分岐に使うのは 2 種)、誰も読まない副次出力**(A5-6 確認)
- 根拠: `exhaustive_claim_allowed=False` や `transition_state_claimed=False` など、artifact 型から自明に決まる定数が個々の artifact に書かれ、テストでも 24 行で assert されている。`write_jsonl` が 89 回、`to_csv` が 23 回呼ばれるが、pipeline 内の読み手は 2 箇所だけ。
- 推奨: 型ごとの主張方針は 1 つの表に移す。next_action は分岐に使う値だけを Enum にする。副次出力は reporting 層が manifest から一括生成する。

**AR-29【低】下流 stage が上流の生ファイルを再ハッシュ・再パースし、method_evidence_* が 5 変種に分裂している**(A5-10 一部訂正)
- 根拠: `sha256_file` と `source_file_hashes` が 21 ファイル・約 40 箇所で呼ばれている。検証・監査系の関数は 47 個ある。(訂正: irc.py の上書きは AND 結合なので「格下げはできるが格上げはできない」。独立端点 QC は科学的に別の検査として意味がある)
- 推奨: 検証はバックエンド境界で 1 回だけ行い、`evidence_hash` と `method_signature` を不変値として発行する。下流は `same_surface(a,b)` で比較するだけにする。

**AR-30【中】テストが形状の確認に偏っていて、壊れやすく、要所が抜けている**(A5-8、A5-9 一部訂正、見落とし追加)
- 根拠:
  - assert 997 個のうち、フラグ・ラベル・状態の確認が約 52%、数値や幾何の確認は 12.6%。`test_reaction_classification.py` は数値 assert が 0。
  - `Artifact(` の手組みが 189 回、fixture は 0。類似のファクトリが 41 個あり、`_minimum` は 4 ファイルで別々に定義されている。
  - YAML 設定値を固定する assert が 58 個あり(nbeads==9 など)、収束対策のためにパラメータを調整するたびに赤くなる。
  - `parse_nwchem_output` の `scf_converged = energy is not None` は、合成スニペットでしか検証されていない。runs には実出力が 695 本あるが、ゴールデンテストは 0。(訂正: SCF 非収束の実例は runs にない)
  - Eyring/Wigner/qRRHO を直接確かめるテストは 0。
  - テスト済みの停滞検出器 `assess_optimizer_stagnation` は、本番コードから一度も呼ばれていない(本番は window=3 の別実装で、事後判定)。
  - `reaction_path_qc.py:306-329` は rdkit の ImportError/OSError を握りつぶし、「TS 未検証(overlap 0)」という科学的な診断に化ける。
  - rdkit の DLL がブロックされると `test_amine_hf_panel.py` の収集エラーでスイート全体が中断する。
- 推奨: `tests/factories.py` を作り、型付き Evidence から組み立てる。`tests/golden/nwchem/` に実出力 10〜20 本(opt+freq、SCF/geometry の最大反復到達、NEB 未収束、ZTS 失敗、kill されたもの)を置く。熱化学は教科書値で確かめる(Eyring: k(298.15 K, 20 kcal/mol) ≈ 1.35e-2 s⁻¹、Wigner: κ(−1000 cm⁻¹, 298.15 K) ≈ 1.97)。設定値を固定するテストは「メソッドが全段で一致する」などの不変条件テスト 2〜3 本にまとめる。原子質量は静的な表にして rdkit への依存をなくす。rdkit 依存テストは `importorskip` にする。

**AR-31【中】死蔵コードと重複実装**(A1-10、A6-10 確認・低、A4-3/C6-10 確認・中)
- 根拠: CLI の 3 エントリポイントと 41 stage を根とする到達可能性解析で、35 モジュールに到達できない。トップレベル定義のうち参照 0 が 112 個・1,890 行(テストからも参照 0 は 106 個・1,623 行)。重複は Wigner 3 実装、Eyring 5 実装、qRRHO proxy 4 種、`_reaction_id` 8 ファイル、`build_reuse_plan` 3 実装、`process_penalty_from_public_data` 2 実装(上限 1.0 と 1.5 で仕様も違う)。
- 推奨: 4.5 節の表のとおり削除する。CI に vulture(min-confidence 80)と到達可能性チェックを入れる。

**良い点**: 失敗を例外ではなく `Artifact.failure` として監査可能に残している。`Manifest.latest_artifacts/merge` の revision 規則は簡潔。虚振動カットオフの解決が `frequency_qc.resolve_imaginary_frequency_cutoff` の 1 箇所にまとまっている。`classify_reaction` の決定表テスト(19 件)は、振る舞いのテストとして有効。docs/current の分量(約 1,180 行)は過剰ではない。fail-closed の方針(欠損や曖昧さは不合格、AFIR のバイアス付きエネルギーは障壁に使わない)は正しく、維持すべきである。

### 2.6 バックエンド抽象化

**評価**: 外部実行の監査基盤(opt-in の subprocess、入力の常時書き出し、command_result)と、NWChemEngine の厳格な成功判定は質が高い。しかし能力単位の統一インターフェースがなく、「経路アルゴリズム × QM プログラム × 科学ゲート × Artifact の発行」が 1 つのメソッドに融合している。

**AR-11【高】バックエンドを差し替えられない(Protocol 未使用、registry に型がない、エンジン名のハードコード、能力が揃っていない)**(A1-3、A6-1 確認、A2-8/A4-10 確認・低)
- 場所: `backends/base.py:10-42`、`backends/registry.py:1-138`、`workflow/reaction_state.py:22-28,150-153,187`、`core/environment.py:63-170`(CC 34/31)
- 根拠:
  - QMBackend/TSBackend/ConformerBackend/DBProvider は定義以外の参照が 0。get_* 8 関数に戻り値型がない。
  - 能力が揃っていない: NWChem は 5 タスク、ORCA は 2、xTB は optimize_frequency だけで、しかも振動数は計算しない。TS 7 エンジンのうち 3 つは run_irc を持たず、nwchem_neb.run_irc は失敗スタブ。
  - thermo の `arkane`、kinetics の `cantera/arkane` は必要なメソッドを持たない。thermo の `dummy/internal/quasi_rrho` はすべて GoodVibes の別名。
  - registry は 27 の具象バックエンドを eager import する(0.40 s、pandas も読み込む)。
  - 訂正: kinetics.py:36-38 が tst に読み替えているため、kinetics.py:44 は実行時には失敗しない(pyrefly の偽陽性)。
- 推奨: `backends/contracts.py` に能力別の `@runtime_checkable` Protocol(SinglePoint / Optimize / Frequency / PathSearch / SaddleRefine / IRC / ThermoEngine / KineticsEngine / ReactionDiscovery)を置く。registry は `"module:Class"` 文字列による遅延ロードにし、型付き getter を用意する。エンジンは `capabilities` と `strategies` を宣言する。エンジン名は `required_executables(settings)` から導出し、core/environment の stage 名分岐は削除する。「backends/ の外に registry キーのリテラルがない」ことを確認するアーキテクチャテストを追加する。

**AR-12【高】TS バックエンドが経路探索・精密化・検証・成果物生成を丸抱えし、QM エンジンを直接生成して NWChem に固定されている**(A6-2 確認、A2-1/A4-4 一部訂正、A1 の見落とし)
- 場所: `backends/ts/nwchem_neb.py:214,694`、`nwchem_saddle.py:145,252,374,561`、`pysisyphus.py:189,239,366-373,404-406,661,832-833`、`pysisyphus_saddle.py:19-26,120-126,292,465`、`orca_nebts.py:533`、`ts/base.py:176`(620 行)
- 根拠: `NWChemEngine()` の直接生成が 8 箇所、`OrcaEngine()` が 1 箇所。pysisyphus は `program=nwchem`、`7.2.3`、`disp_vdw==3` 以外を ValueError にする。`PysisyphusEngine.search_ts` は既定で `NWChemNEBEngine` に委譲するので、エンジン名と実際の処理が一致しない。pysisyphus_saddle は pysisyphus の非公開関数を 5 つ import している。method dict のキーは 130〜144 種。timeout の既定値は 3,600〜172,800 s とばらつき、経路探索と内側の saddle で同じ値を共有している。
- 訂正: `method` 経由で ncores 等の設定は QM エンジンに渡っている(問題は NWChem への固定結合のほうである)。`validate_ts_frequency_calculation` は 4 エンジンで共有されていて、独自ゲートは orca と dummy だけ。重複は「数百行のコピペ」ではなく、「同じ骨格の並行実装。実質一致は NEB/String で約 170 行、Saddle 系で約 100 行」。
- 推奨: `PathGenerator`(render/run/parse → PathResult)/ `SaddleRefiner`(`get_ts_engine(spec, qm=get_qm_engine(...))` で QM を注入)/ `chemistry/ts_validation.py`(判定)/ `ts/artifacts.py` に分ける。`ts/settings.py` に pydantic の設定を置き(`extra='forbid'`、`path_timeout_s` と `saddle_timeout_s` を分離)、`_reaction_id` は core.ids に一本化する。backend-categories-independent 契約を ignore なしで通すことを受け入れ基準にする。

**AR-13【中】stage が具象バックエンドを直接 import している(7 件)。registry の名前と能力が一致しない**(A1-6 確認)
- 根拠: build_complexes.py:671 の `CRESTConformerBackend()` は Protocol のメソッドしか使っていないので置き換えられる。`recover_path` は NWChem 名のモジュールにある汎用関数 `validate_path_trajectory`(Kabsch と原子の並び)を import している。受け皿になるべき `chemistry/path_validation.py` は importer が 0 で、中身は HF 専用で矛盾する規則(RMSD≤1.25 Å で合格)になっている。
- 推奨: `get_conformer_engine(config.get('engine','crest'))` を使う。汎用の経路検証関数は chemistry へ移し、path_validation.py は置き換える。registry を能力別に分ける(`get_rate_engine` / `get_mechanism_exporter`)。

**AR-14【高】chemistry/stages に NWChem 固有のパーサ・既定値・ゲートが入り込んでいる**(A1-4 確認、A2-6/A6-8 確認・中)
- 場所: `chemistry/nwchem_evidence.py:18-24,71-160,205`、`chemistry/minimum_recovery.py:119-176`、`chemistry/method_lineage.py:127-152`、`chemistry/path_diagnostics.py:17-28`、`stages/dft_minima.py:9-14,421-444,515`、`stages/path_ensemble_assess.py:9-12,70-80`、`stages/thermo.py:8,295-310`
- 根拠: 文字列 "nwchem" の出現は chemistry に 56 回、stages に 16 回、core に 17 回。`REQUIRED_ENDPOINT_NWCHEM_METHOD`(7.2.3/PBE0/def2-SVPD/D3)は呼び出し元 0 の死蔵関数(CC 66)だけが使っていて、有効化すると B3LYP 系のパイプラインを誤って不合格にする。`developer_architecture.md:34` の「Stages do not parse external program output」に違反している。
- 推奨: バックエンドがパース時に正規化済みの `MethodEvidence(functional, basis, dispersion, grid, scf_tol, solvation, program_version, input/output sha256)` を 1 回だけ発行する。chemistry は `evidence_matches(MethodSpec, MethodEvidence)` で比較するだけにする。NWChem のパーサは `backends/qm/nwchem/` に移す。死蔵関数は削除する。chemistry/stages に `'nwchem'` のリテラルがないことを確認するテストを追加する。

**AR-15【中】TS の科学ゲートがエンジンごとに違い、ORCA 経路は production で検証を通れない**(A1-5 確認、A2 の見落とし)
- 根拠: orca_nebts は `validate_ts_frequency_calculation` を使わず、`mode_overlap_threshold`(0.7)を別フラグの判定にしか使っていない。キー名は `ts_mode_overlap_threshold` と不統一。orca_nebts の `reaction_validated.qc` には `method_evidence_validated` がないため、`has_validated_ts(require_real=True)` を必ず通れない。thermo は `real_orca_executed` を IRC 実行の証拠として扱っている。ORCA の IRC 経路にもダミー代替が埋め込まれている(`orca_nebts.py:626`)。
- 推奨: TS ゲートを `chemistry/ts_validation.py` に移し、ts-search で全エンジンに一律に適用する。バックエンドは生の証拠だけを返す。実行の来歴は `qc['execution']={engine, task, real_executed}` にまとめる。

**AR-16【中】ダミー実装が本番のクラス階層・既定値・registry に混入している**(A6-3 一部訂正)
- 根拠: `XTBEngine(DummyQMEngine)` は single_point をオーバーライドしていないので、`allow_subprocess=True` でも捏造エネルギーを `status=success, engine='xtb'` で返す(実行して確認)。`fallback_to_dummy` の既定値は xtb だけ True。stage の既定エンジンは 'dummy'(dft_minima.py:617、preopt.py:198、sp.py:58、ts_execution.py:26)。`_dummy_fallback` は 4 箇所にある。
- 訂正: `qc.engine_is_dummy=True` が付くので、core/qc と thermo のゲートで除外される。sp+xtb の組み合わせも現行設定にはない(潜在的な問題)。
- 推奨: ダミーは `hfauto/testing/fakes.py` に移して tests の conftest で登録する。本番 registry から 'dummy' を削除し、stage の engine を必須にする。production の preflight ではダミーとフォールバックを拒否する。

**AR-17【中】method dict が「科学的手法・資源・ワークフロー状態」を混ぜた約 136 キーの型のない契約になっている**(A6-5 一部訂正)
- 根拠: `method_for_reaction_case` が path_strategy や basin_assessment を method に注入している。別名が並存している(ncores/nprocs、memory_mb/maxcore_mb、nbeads/path_image_count など)。calc_id と checkpoint の再利用は method の全キーを含むので、timeout や実行ファイルパスを変えただけで再利用が無効になる。memory_mb の意味は経路によって rank あたり約 2 倍異なる(訂正: レビューでは 4 倍としていた)。
- 推奨: `MethodSpec`(hash 対象、allowlist 方式)/ `ExecutionSpec`(資源、hash 対象外、サイトプロファイルから注入)/ `SearchHints`(ワークフロー状態)の 3 つに分ける。

**AR-19【中】外部プログラム I/O の堅牢性がバックエンドによって不均一**(A6-7 一部訂正、見落とし追加)
- 根拠: ORCA は、振動数がないとき Gibbs に電子エネルギーを代入し、IRC 端点を mtime 順で選ぶ(TS 自身を端点にしうる)。電子状態は `resolve_electronic_state` を通していない。runs に実行実績は 0、テストもない。CREST は欠損エネルギーを **0.0 kcal/mol(最安定)** にする(訂正: +999 ではない。実害はさらに大きい)。ReaDuct は同一プロセスで実行され、timeout がない。
- 推奨: すべての外部プロセスを `run_command` 経由にする。パーサは型付きの ParsedOutput を返し、取れない値は None にする(代用しない)。未検証のバックエンドは golden テストを満たすまで experimental 扱いにする。

**良い点**: `run_command`(POSIX は killpg、Windows は taskkill /T)。NWChemEngine の成功判定(版数一致、出力上での D3 適用、3N 行の振動数、原子順序、final-NNN.xyz の検証)。`TSSearchResult`/`IRCResult`/`DiscoveryResult`(frozen dataclass で Artifact を返さない。TS/IRC 系の手本になる)。`resolve_electronic_state` による電荷・多重度・パリティの検証。DB provider がキャッシュ優先で、ネットワークは opt-in。

---

## 3. 化学計算レビュー

### 3.0 実行実績から見た非収束・長時間化の実態

**全体像**
- TMA/amine+HF/HF2 系: 86 run(2026-07-14〜08-16)、記録上の壁時計時間は約 70〜75 h。内訳は amine_hf_pilot_prod 約 19 h、m3_* 約 21 h、tma_hf2_* 約 10 h、amine_hf_pilot_dft_followup 約 23 h。**検証済みの TS/IRC は 0 件**。
- それ以外(HCN/HONO/NH3/SiO2 ほか): hf2 端点再評価 20.4 h、aniline topology ensemble 28.5 h(いずれも削除済み stage による)。
- NWChem(mpirun)の合計は 96.6 h で、非ゼロ終了またはタイムアウトが 14 件・33.1 h(34%)。command_result 428 件の集計では 124.1 h 中 38.0 h(31%)が失敗ジョブ。
- 化学的に妥当に完了したもの: HCN→HNC(ΔE‡ 46.67、ΔG‡ 42.28、ΔG_rxn 12.74 kcal/mol)、HONO trans→cis の TS(−680i、ΔE‡ 13.6 kcal/mol)、NH3 反転(−757i、4.25 kcal/mol)、same_basin_water、SiO2 v2 のスキャン+DFT 再重み付け(TS ではない)。

**パターン表**

| run(群) | 手法 | 失敗様式 | 推定原因 | 消費 |
|---|---|---|---|---|
| m3_real_qm_nh3_hf_001、amine_hf_pilot_prod_001 | xTB→NWChem PBE0/def2-TZVPD | イオン対 seed が中性錯体へ崩壊し、RC/IP の 5 対が同一構造に収束(RMSD ≤2.9e-4 Å) | 気相で HF が 1〜2 分子ではイオン対の極小が存在しない(Legon の分光結果と整合) | pilot 約 19 h |
| m3_real_qm_tma_hf2_001/002、m3_tma_*(54 run)、tma_hf2_path_v1、tma_hf2_recovered_path_v1 | NEB 3/5/11 bead、string 5/9/11/13 bead、bracket、NWChem saddle、RS-I-RFO、dimer | TS は 0 件。経路は単調下降に緩和し、最終判定は same_basin_relaxation | 両端点とも共有プロトン構造(N–H 1.317/1.271 Å、差 0.04 Å)。反応物は xfine で −27.09 cm⁻¹ の Cs 鞍点。gate が緩い(min_change 0.05 Å) | 約 19 h + 別パイプラインで約 10 h |
| m3_tma_p0_saddle_bracket_001、m3_tma_local_bracket_001 | 端点の直交座標線形補間ブラケット | f=0.5 に +51.3 kJ/mol の「極大」。seed Hessian は虚振動 4 本(−1324〜−77 cm⁻¹) | 線形補間による人工障壁。単調と分かった string 経路を捨てて線形補間をやり直した | Hessian 926 s + 再計算 |
| m3_tma_string_rsirfo_saddle_001〜005、m3_tma_p0_dimer_smoke | pysisyphus 1.0.0 + QCEngine + NWChem | 上り坂(+14.9 / +33.5 kJ/mol)、trust を絞ると極小へ崩壊、dimer は極小から始めて cycle 0 で「収束」 | **座標系の不一致(82.9°)**、grid medium/SCF 1e-6 の別 PES、cartesian 基底(−38.7 kJ/mol のずれ)、root 0 固定 | 各約 20 分 + Hessian 3 回 |
| m3_tma_mode3_saddle_001、string_saddle_001、mep_segment01_* | NWChem saddle(inhess 2) | n_imag=0 の中間極小(−374.7537 Eh)か反応物へ崩壊。trust 0.15→0.02 でも停滞するだけ | TS が存在しない。seed Hessian の二重計算(約 920 s/回)。prominence 閾値を SCF ノイズ以下(1e-8 Eh)まで下げた | 4,274 + 2,526 + 約 1,600 s |
| m3_tma_path_ensemble_001、string_continuation_001、refinement_009_003 | NWChem ZTS | 「converged」だが gmax は 9.1e-4→9.15e-3 と 10 倍に悪化 | **ZTS は xrms/xmax だけで収束を判定する**。Verlet 切替後の微小移動で偽収束 | 約 2 h + 1,127 s |
| m3_tma_endpoints_xfine_001 → endpoint_mode_minima / endpoint_seed_minima | xfine 再最適化、mode-follow | 反応物は Cs 鞍点で、±方向とも生成物極小へ。seed も同じ鞍点へ戻る | 対称な seed。fine grid では +6.83 cm⁻¹ に見えて問題が隠れていた。GFN2 preopt で対称性が戻る | 4.5 h + 3.2 h |
| tma_hf2_explore_v3〜v11 | SCINE ReaDuct 6.1.0 NT2/AFIR(GFN2) | `low_level_ts_validated` は全件 False。候補は 9→0 | **虚振動を射影なしの固有値で数えている**(常に 2〜4 本)。等価な置換 | 数分 × 9 run |
| amine_hf_pilot_hf1_hf3_v1/v2 | CREST 3.0.2 --nci | rc 1(トポロジー変化の検出) | 旧 seed 生成の不具合。現行 v3 では発生しない | 1 回 18〜1,983 s |
| amine_hf_pilot_hf1_hf3_v3_discovery | ReaDuct | TS guess なし 22 件(すべて単調上昇 +12〜247 kJ/mol)、虚振動 2〜4 本 7 件、AFIR 未収束 10 件 | 陰性結果(生成物 basin なし)を失敗として扱っている。数え方のバグ | 約 4.5 分 |
| amine_hf_pilot_dft_followup_v1(00_dft-minima) | NWChem xfine/1e-8、np 2 | 4 h タイムアウト 2 件。最初の SCF が 29 s→507 s | **WSL 4 vCPU 上で np4 と np2 を同時実行**、/mnt/c の scratch、Cartesian | 15.7 h |
| amine_hf_pilot_dft_followup_v1(04_ts-search) | NEB→saddle | NEB の Gmax が 1.085→1.066 で停滞、bead 間の幅 406 kcal/mol、seed の虚振動 6 本 | +225 kJ/mol の開環体をエネルギー上限なしで DFT へ昇格。線形補間で原子が衝突 | 7.3 h |
| hcn_ts_irc_validation(NEB) | NWChem NEB、impose なし | 150 反復で未収束、原子間最小距離 0.529 Å、中央 bead が +3.2 Eh | 生成物の向きを揃えていない線形補間 | 1,127 s(同じ系で string→saddle は 849+47 s で成功) |
| hcn_irc_v3/v4/v5、m3_real_qm_ammonia_inversion_003/004、hono_isomerization_v3(irc) | pysisyphus EulerPC + QCEngine | 遷移ベクトル −84.4 cm⁻¹(真値 −1131.6)、初期変位 2.5 bohr、端点最適化でエネルギー上昇。v5 は両分岐が TS より上で終わったのに合格 | **座標系の不一致**。IRC 合否が分岐の下降を見ていない | 107〜818 s/回 |
| hcn_alternative_pes_v1〜v4、hono_isomerization_v1〜v3 | NWChem NEB/string、直列、xfine、TZVPD | 手動 kill(rc −15)。HONO NEB は最後の 58 反復がステップ 0 で停滞 | 過剰な精度、ncores 1、停滞を検知しない、初期経路での原子衝突 | HONO NEB 2.5 h |
| hono_isomerization_v3(classify/thermo) | reaction-classify + GoodVibes | 正しい TS と IRC があるのに unresolved。ΔG‡ 24.11 kcal/mol を公開 | **`basin_assessment` のキー欠落(回帰)**、**ZPE の二重計上** | — |
| hcn_goodvibes_v1/v2、hcn_thermo_v5/v6_complete | GoodVibes / 内部フォールバック | 線形分子で GoodVibes が異常終了。フォールバックで ΔG_rxn 18.05(正しくは 12.74) | 回転定数 `**********`。NWChem RRHO で HNC の S_rot が負 | 秒単位 |
| hf2_endpoint_reassessment_001/002、aniline_hf2_topology_ensemble_001 | NWChem noautoz、PBE0/def2-SVPD | 50〜121 step、7,200 s タイムアウト 8 件、軟モードの虚振動 rescue に約 4 h | Cartesian 座標と平坦な H 結合 PES、同一計算の再実行(削除済み stage) | 20.4 h + 28.5 h |
| sio2_wet_etch_aniline_001/002 | xTB スキャン + DFT SP | 001 は無効化($fix が守られない)、002 は xTB と DFT で極大の位置が一致しない | xTB の拘束指定、GFN1 を etemp 1000 K で使用 | NWChem 0.48 h(破棄)+ 1.07 h |

**結論**: 長時間化の主因を順に並べると、(1) 存在しない TS を探索したこと(化学モデルの選択と前段ゲートの欠如)、(2) pysisyphus/QCEngine の座標系と PES の不一致、(3) 経路法の初期化と収束判定の欠陥、(4) 過剰な精度と資源競合、(5) Hessian の重複計算と再開機能の欠如、である。beads、stepsize、trust、grid、基底、座標系を試行ごとに変えても、(1)(2) は変わらない。

### 3.1 配座・初期構造(CREST / xTB preopt / build-complexes / RDKit)

**使用ライブラリと仕様上のポイント**
- **CREST 3.0.2**: CREGEN は既定で、入力と共有結合トポロジーが異なる構造を破棄する。`--notopo [atoms]` で無効化、`--noreftopo` で初期判定のみ無効化できる。`--ewin` は既定 6、`--quick` では 5 kcal/mol。`-T` はスレッド数。`--alpb <solvent>`、`--cinp`(拘束)、`--scratch` がある。トポロジー停止時の推奨は A) 事前に最適化する、B) `--noreftopo`、C) 拘束する([keywords](https://crest-lab.github.io/crest-docs/page/documentation/keywords.html)、[JCP 2024](https://pubs.aip.org/aip/jcp/article/160/11/114110/3278084/CREST-A-program-for-the-exploration-of-low-energy))。
- **xtb 6.7.1**: `--ohess [level]` は最適化してから Hessian を計算し、虚振動方向へ変位した `xtbhess.xyz` を書く。要約表示は虚振動カットオフ −20 cm⁻¹ で数える。`--cycles` を超えても rc=0 と `normal termination` を出力し、`FAILED TO CONVERGE` は stdout に出るだけ(検証者が実行して確認)([hessian](https://xtb-docs.readthedocs.io/en/latest/hessian.html)、[optimization](https://xtb-docs.readthedocs.io/en/latest/optimization.html))。GFN2 のプロトン移動ベンチマーク(PX13/WCPT18)の精度は PBE0 並みで、スクリーニングには妥当([JCTC 2019](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176))。
- **RDKit 2026.03.5**: ETKDGv3 の既定は `embedFragmentsSeparately=True`(断片が原点で重なる)。HF を含む系は MMFF94s のパラメータがなく、UFF になる([rdDistGeom](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html))。

**現状の設定**: `crest input.xyz --gfn2 --quick --nci` だけで、ewin・溶媒・トポロジーのオプションは渡らない。threads は OMP_NUM_THREADS 経由で効く(訂正済み)が、pilot は threads 1 に設定されていた。preopt は `xtb --opt <level>` だけで Hessian を計算しない。配置 seed は chain/star の交互、距離 1.7/2.1 Å の周期 2、twist 120°·seed。pilot は `max_nci_complexes 1`、`selected_max 1`。

**問題点**

- **CH-18【高】CREST の既定トポロジーフィルタがプロトン移動・イオン対の構造を捨てている。設定も伝わっていない**(C1-1 一部訂正、C1-6/C1-7 訂正で低)
  - 根拠: TMA_HF3 の CREGEN のトポロジー由来除去は累計 586 件(反復ごとの CREGEN 呼び出しの累積なので、586 個の固有構造という意味ではない)。残った 18 配座はすべて中性。最安の nci_0 を壁なしで緩和すると N–H が形成されて −4.64 kcal/mol になり、DFT でもイオン対が 29.9 kJ/mol 安定だった。
  - 訂正: イオン対は preopt/relaxation-discovery の経路で部分的に拾われている。したがって「探索集合に入らない」のではなく、「CREST のアンサンブルには入らず、偶然の緩和経路に依存している」。m3 の ewin 6.0 は CREST 内部の 5.0 で先に打ち切られる。v1/v2 のトポロジー停止は旧 seed によるもので、現行 v3 では発生していない。
- **CH-19【高】決定論的な配置 seed が対称で多様性を欠く**(C1-2 一部訂正、見落とし追加)
  - 根拠: 受容方向が「結合ベクトル和の逆向き」なので N···H–F···H–F が C3 軸上に一直線に並ぶ。twist は線形分子の HF に対しては恒等変換になり、seed00 と seed02 が完全に一致する。直線鎖は GFN2 で −14.6×2 cm⁻¹、DFT で −80×2 cm⁻¹ の鞍点で、TMA-(HF)2 に 5,518 s を使った。実験の (HF)2 は ∠F–F–H 約 115° に屈曲している。
  - 訂正: star 配置の衝突は旧コードで起きたもので、現行コードには 16 方向のフォールバックがある。ただしフォールバック方向は**実験室座標系に固定された方向**で、結果が入力 xyz の向きに依存する。
  - 見落とし: guest の donor に「最初の H 原子」を選ぶため、C–H が donor になりうる(`build_complexes.py:178-199`)。衝突閾値 0.55〜0.60 Å は原子の重なりに近い値で、閾値として小さすぎる。
- **CH-20【高】xTB preopt が Hessian を計算しないため、鞍点を極小として昇格させる。xtb の非収束も検出できない**(C1-3 確認、見落とし: BUG-07)
  - 根拠: tma_hf2_discovery_v3 の preopt 6 構造を GFN2 Hessian で調べると、**5 構造**が非極小(−14.5〜−22.2 cm⁻¹)。一方、`xtbhess.xyz` から `--ohess vtight` を 1 回かけるだけで、秒単位で真の極小(全振動数が実数)に到達した。`geometry_converged = parsed or opt_xyz.exists()` は実質常に真になる。
  - 注: C2 の検証者は「xtb 6.7.1 は非収束時に非ゼロ終了する」と述べたが、実行による確認はしていない。本書は C1 の検証者が実行した結果(rc=0)を採用する。
- **CH-07【高】共有結合半径 ×1.25 という単一閾値の結合判定が、水素結合・共有プロトン・FHF⁻ で反転する**(C1-4/C4-6 確認、C5-5 一部訂正で中)
  - 根拠: N–H のカットオフが 1.275 Å、H–F が 1.100 Å になる。FHF⁻(1.139 Å)は 3 つの断片に分割される。m3 の端点(差 0.04 Å)が「N–H 形成・H–F 切断」と判定される。RDKit の既定係数 1.3 なら両者は同じグラフになる。つまり「反応」があるかどうかが係数の 0.05 の差で決まっている。connect-minima は Σdelta ≥0.05 Å を受理し、24 件中 13 件が accepted(うち 8 件は 0.063〜0.084 Å の変化)。影響範囲は preopt、relaxation-discovery、explore の端点、basin 同一性、minimum 接続、path_initialization、reaction_trials、complexes の配置。
  - 訂正: m3 のペアそのものは、利用者が明示した仮説(min_change 0.05)から来たもので、グラフだけが原因ではない。
- **CH-22【中】アンサンブルの選別が GFN2 最安の 1 個に偏っている**(C1-5 一部訂正): 最初の seed が組成の上限を使い切る。preopt 後に幾何的な重複除去がない(nci_2 と nci_3 は 0.046 Å 差)。DFT での再ランキングがない。GFN2 最安の Cs 平面構造は DFT xfine では鞍点だった。(訂正: 現行設定の `required_source_data_keys` で配置 seed は下流へ行かない。`selected_max 1` は設定上の選択)
- **CH-21【中】geometry_qc の HF 判定が死んでいる**(C1-8 確認): ComponentRecord に `name` がないため HF 成分が見つからず、xtb.py の `setdefault` で `proton_transferred_unintentionally=False` と「問題なし」が記録される(TMA_HF3 では実際には N–H が形成されていた)。
- **CH-23【低】RDKit の複数断片は原点で重なる**(最小断片間距離 0.49〜0.69 Å で、0.55 Å の衝突判定も通過してしまう)。detect-sites は build-complexes に配線されていない(C1-9 一部訂正、C1-10 確認)。

**推奨設定値・プロトコル**
1. 手順を **配置 seed → `xtb --ohess vtight`(ν < −10 cm⁻¹ が残れば `xtbhess.xyz` から再最適化、最大 3 回。n_imag と最低振動数を記録)→ CREST `--nci` → `--ohess` → 重複除去** にする。収束判定は `FAILED TO CONVERGE` の文字列、NOT_CONVERGED ファイル、returncode で行う。
2. acid/base を含む組成では CREST に `--notopo <酸性 H と N/F/O>`(または `--noreftopo`)を付ける。出力は q=r(H–F)−r(B–H) で neutral/shared/ion_pair に分類し、状態ごとに上位 4〜6 個(ewin 4 kcal/mol)を保持する。CREGEN の除去数は qc に記録する。
3. CREST バックエンドにオプション対応表を持たせる: `energy_window_kcal_mol→--ewin`、`threads→-T`(既定 4)、`solvent→--alpb`、`topology_check:off→--notopo`、`scratch→--scratch`(WSL の ext4 上)。`crest=3.0.2`、`xtb=6.7.1` に固定し、`--version` の出力を provenance に記録する。
4. 配置 seed: 単結合の F/O 受容体は X–H 軸から 110〜120° の lone-pair 円錐上に置き、方位角を 0/120/240° で振る。sp3 の N は lone-pair 軸から 5〜10° 傾ける。seed 固定の RNG で ±30° の乱数回転を加える。衝突判定は非結合 H···H ≥1.2 Å、重原子間 ≥2.2 Å、約 50 回の rejection sampling。置換不変の距離 RMSD <0.1 Å で重複を除く。n=1 の直線 C3v/C2v は真の極小なので、現行のままでよい。
5. 重複除去は CREGEN 相当(ΔE <0.1 kcal/mol、置換不変 RMSD <0.05〜0.1 Å、回転定数差 1% 以内)とし、ewin 3〜4 kcal/mol 内の上位 5〜8 構造を安価な DFT 一点計算(r2SCAN-3c または PBE0-D3/def2-SVP の fine)で並べ直して、上位 2〜3 個を opt+freq に回す(CENSO 型、[JPCA 2021](https://pubs.acs.org/doi/10.1021/acs.jpca.1c00971))。
6. 結合判定: X–H は `r ≤ 1.15Σr_cov` なら結合、`r ≥ 1.45Σr_cov` なら非結合、中間は partial としてトポロジー変化に数えない。あるいは xTB の Wiberg 結合次数(BO ≥0.75 で結合、≤0.15 で非結合。NT2 と同じ基準)を使う。プロトン移動と認めるのは |q| ≥0.3 Å で符号が反転した場合だけ。FHF⁻/(HF)nF⁻ は 1 断片として扱う。
7. RDKit で `'.'` を含む SMILES を扱うときは failure にし、build-complexes へ回す。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| CREST `--nci`(GFN2、--quick) | **有効** | amine·(HF)2 の屈曲 H 結合ネットワークを 4 スレッドで約 33 s で発見した。nci_0 は DFT でも生成物極小に収束 |
| CREST(既定トポロジーフィルタのまま、酸塩基系) | 限定的 | 生成物側の basin を系統的に除外する |
| xTB preopt(`--opt` のみ) | 限定的 | 鞍点を昇格させる。`--ohess` にすれば有効 |
| GFN2 無バイアス緩和・relaxation-discovery | **有効** | TMA·(HF)3 のイオン対を数分で発見し、DFT(29.9 kJ/mol 安定)と整合 |
| 決定論的な配置 seed(現行) | 限定的 | 対称な鞍点と重複 seed を作る |
| RDKit(複合体) | 低価値 | 断片が重なる。単量体専用にすべき |
| detect-sites | 低価値(未配線) | 出力が使われていない |

### 3.2 電子状態計算設定(NWChem / xTB / ORCA)

**使用ライブラリと仕様上のポイント**
- **NWChem 7.2.3**: grid の目標精度は fine が約 1e-7、xfine が約 1e-8。SCF の収束を厳しくすると積分スクリーニング(itol2e、tol_rho)も自動的に厳しくなる。geometry の既定は冗長内部座標で、`noautoz` は Cartesian。4 原子以上の厳密に直線な鎖では autoz が失敗しうる(FAQ)。`inhess 2` は直前の frequency で得た Cartesian Hessian を使う。`scratch_dir`、`grid nodisk`、`vectors input` で再開できる。COSMO の振動数は有限差分で計算される([DFT](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html)、[Geometry-Optimization](https://nwchemgit.github.io/Geometry-Optimization.html)、[GEOMETRY keywords](https://nwchemgit.github.io/Keywords-for-the-GEOMETRY-directive.html)、[FAQ](https://nwchemgit.github.io/FAQ.html)、[Memory](https://nwchemgit.github.io/Memory.html)、[Scratch_Dir](https://nwchemgit.github.io/Scratch_Dir.html)、[COSMO](https://nwchemgit.github.io/COSMO-Solvation-Model.html))。
- **Open MPI 5.0.10**: 自分が oversubscribe されていると知らない rank は busy-poll し、性能が著しく落ちる。別々に起動した mpirun 同士は互いを知らない([scheduling](https://docs.open-mpi.org/en/main/launching-apps/scheduling.html))。WSL から /mnt/c にアクセスすると I/O が遅い([WSL filesystems](https://learn.microsoft.com/en-us/windows/wsl/filesystems))。
- **ORCA 5/6**: `Grid5/FinalGrid6` は廃止され、`DefGrid1-3` になった。hybrid では RIJCOSX が既定([数値積分](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/numericalintegration.html)、[ORCA 6.1 DFT](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html))。
- 汎関数の知見: PBE0(厳密交換 25%)は非局在化誤差で偽のプロトン移動や過大な H 結合を与えうる。系によって必要な交換割合は 8〜90%([PCCP 2024](https://pubs.rsc.org/en/content/articlehtml/2024/cp/d4cp00907j))。D3 は BJ damping が推奨([JCC 2011](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21759))。低振動モードの自由エネルギーは積分 grid と分子の向きに敏感([Bootsma & Wheeler](https://chemrxiv.org/engage/chemrxiv/article-details/60c7436f0f50dba3a3395ef2))。

**現状の設定**: 全本番パイプラインが PBE0 + `disp vdw 3`(D3 zero damping)/ def2-SVPD(一部 TZVPD)、grid xfine、SCF 1e-8、`nocenter noautoz noautosym`、`inhess 0`(対角)、`memory_mb 4000 × ncores 4`、scratch `.`(/mnt/c 上)。溶媒なし(SiO2 だけ `extra_blocks` で COSMO)。amine 系の driver は現行 default(訂正: tight 固定は HCN/HONO 系だけ)。

**問題点**

- **CH-13【高】WSL の実行資源の過剰割当**(C2-1 一部訂正、見落とし追加)
  - 根拠: `.wslconfig` は processors=4、memory=12 GB(実効 11 GB)。ホストは 20 コア / 63.6 GB。np4 と np2 の mpirun を同時に起動し、同一の 137 AO 入力でも最初の SCF が 29 s→507〜527 s に遅くなった(NH3_HF3 は 4 h でタイムアウト)。`memory 4000 mb/rank × np4 = 16 GB` は単独ジョブでも VM の上限を超える宣言になっている。`run_command` は `os.environ` をそのまま継承するので、OMP_NUM_THREADS が設定されていると rank ごとにスレッドが増える。
  - 訂正: 影響は非対称で、np4 の TMA 系 Hessian は約 5% しか遅くなっていない。recovered_path の SP 124 s は主に xfine/1e-8 の効果。
- **CH-12【高】xfine/SCF 1e-8 を全段に適用し、lineage が経路探索段にまで同じ設定を強制している**(C2-2 確認)
  - 根拠: HCN で fine と xfine の差は ΔE‡ 5.86e-5 kcal/mol。TMA·(HF)2 の解析 Hessian は 925 s→3,301/3,474 s(3.6 倍)。xfine の 11-bead NEB は 1 反復約 45 分で、maxiter 220 では約 7 日かかる見込みだった。`PATH_NUMERICAL_KEYS` は、path_ensemble 系で端点と経路の grid の完全一致を要求する。
  - 補足: 軟モードの符号は grid で変わる(TMA·(HF)2 の反応物は fine で +6.8、xfine で −27 cm⁻¹)。xfine は opt-in で残す価値がある。
- **CH-14【高】`noautoz` 固定、対角初期 Hessian、timeout 後の初期 seed からの再実行**(C2-3 一部訂正)
  - 根拠: 8 原子で 112〜122 step、aniline 系で 49〜72 step(7,200 s タイムアウト 7 件)、m3 の mode minima で 96/69 step。timeout 後は `final-NNN.xyz`・movecs・`.drv.hess` を引き継がずに seed からやり直す。TRIC は非共有結合系で既存の座標系より一貫して優れる([JCP 2016](https://pubs.aip.org/aip/jcp/article/144/21/214108/313176/Geometry-optimization-made-simple-with-translation))。
- **CH-15【高】化学モデル: 気相の PBE0-D3(0)/def2-SVPD で amine·(HF)n のプロトン移動状態を判定しており、溶媒の経路がない**(C2-6 確認、C1-11 一部訂正で中、C5-9 確認で中)
  - 根拠: 気相の H3N···HF と Me3N···HF は水素結合錯体で、プロトン移動の程度は F<Cl<Br<I の順に増え、気相でイオン対になるのは HI だけ([Legon 1993](https://pubs.rsc.org/en/content/articlelanding/1993/cs/cs9932200153))。NH3·HX のプロトン移動には水が必要([JPCA 1999](https://pubs.acs.org/doi/10.1021/jp991918j))。反応モード(1,000〜2,500 cm⁻¹)の ZPE は 1.4〜3.6 kcal/mol で、0.5 kcal/mol の古典障壁は低障壁水素結合の領域にある([Perrin & Nielson](https://www.annualreviews.org/doi/pdf/10.1146/annurev.physchem.48.1.511))。GFN2 のアミンのプロトン移動 MUE は 22.2 kJ/mol([PMC12288007](https://pmc.ncbi.nlm.nih.gov/articles/PMC12288007/))。組成の phase/environment は CREST・xtb・ReaDuct のどれにも渡らない。
  - 訂正: 「n=1,2 に TS がない」は現行レベルの PES での結論。n=2 の共有プロトン端点は PBE0 の非局在化誤差で過安定化している可能性があり、高レベル(ωB97X-D/TZ、DLPNO-CCSD(T))で確認するまでは「未確定」として扱うべき。ゲートを置く場所は build-complexes ではなく、generate-reactions(PT 仮説の生成)と TS 探索の前である。
- **CH-11【高】saddle の前に同じ構造の解析 Hessian を二重に計算している**(C2-5 確認・中、C3-6 確認・高): `seed_hessian/` で計算した後、`saddle_initial_hessian_only=True` の saddle 入力にも `task dft hessian` を出力する。外部から Hessian を供給しても同じ。1 回あたり 900〜1,000 s(xfine では 3,300 s)の重複が 5〜6 回、合計約 1.5 h 以上。pysisyphus も `hessian_init: calc` で再計算する。
- **CH-16【中】`solvation_model/dielectric` は method record に記録されるだけで、レンダラが使っていない**(C2-7 確認): COSMO は `extra_blocks` 経由でしか指定できず、`METHOD_KEYS` にも含まれないので、気相の端点と COSMO の TS の組み合わせが「同一 PES」として通過しうる。
- **CH-17【中】虚振動カットオフ −1 cm⁻¹ と自動変位の未配線**(C2-8 一部訂正、C4-8 確認、C6-6 一部訂正で中)
  - 根拠: 対称な seed が鞍点に収束する((HF)2 −427×2、TMA-(HF)2 −83×2、TMA·(HF)2 −27 cm⁻¹)。qRRHO での G への影響は 0.1〜0.4 kcal/mol で、汎関数の誤差(HNC−HCN で −2.0 kcal/mol)より小さいのに、同一条件の再計算と xfine rescue に数時間を使っていた。
  - 訂正: `minimum-mode-follow/assess`(± 変位 seed)は既に存在するが、どのパイプラインにも組み込まれていない。閾値は設定で変更できる。rescue ループの実装は現在の作業ツリーでは削除済み。
- **CH-24【中】ORCA バックエンドが古く、検証もされていない**(C2-10 一部訂正、C3-10 確認・低): `sp_orca.yaml` の `Grid5/FinalGrid6` は ORCA 5 以降で拒否される。Gibbs を電子エネルギーで代用している。`resolve_electronic_state` を通していない。成功判定は returncode と energy だけ。`environment.production.yml` に ORCA はなく、runs での実行実績も 0。`dft_minima_orca.yaml` は `fallback_to_dummy: true`。
- **CH-25【低】細部**(C2-11、C2-12 一部訂正): `hf_stretch_cm1 = max(振動数)` は錯体では C–H/N–H の伸縮を拾う(ただし descriptors stage は未使用)。method-panel の鍵が method_id だけで、同じ ID の別 grid の結果を上書きする。D3(0) の固定は QCEngine 経路だけ。xtb の `method_id` の既定値と env 未指定。

**推奨設定値・プロトコル**
1. **実行資源**: `.wslconfig` を processors=16、memory=48GB、swap=8GB にする。run_root と scratch は WSL の ext4 に置き(NWChem 入力に `scratch_dir /home/user/scratch/<job>`)、メモリに余裕があれば `grid nodisk` にする。200〜350 AO なら memory は 1,000〜1,500 MB/rank。起動時に「Σncores ≤ nproc かつ Σ(memory_mb×ncores) ≤ 0.8×MemTotal」を検査するグローバル予算(ファイルロック)を置く。NWChem には `OMP_NUM_THREADS=1`、xtb には `OMP_NUM_THREADS=<n>,1` と `OMP_STACKSIZE=4G` を明示する。
2. **数値精度**: 既定は grid fine、SCF 1e-7。lineage は 2 層に分け、「停留点層」(minima/TS saddle/freq/IRC/thermo SP)だけで完全一致を要求する。「経路探索層」(NEB/string/bracket)は fine(場合によっては medium)を許す。xfine は |ν_min| <50 cm⁻¹ の確認時だけ opt-in で使う。
3. **最適化**: `coordinates: auto|internal|cartesian` を設ける。auto は autoz で始め、`AUTOZ failed` が出たら noautoz で再投入する。178° を超える直線鎖は事前に 2〜3° 曲げる。初期 Hessian は GFN2 の `xtb --hess` を `.hess` に変換して `inhess 2` で読ませる。minima の driver は default(gmax 4.5e-4)を上限とし、maxiter 80〜100。極小かどうかは Hessian で判定する。timeout 時は `final-NNN.xyz`・`.movecs`・`.drv.hess` を新しい attempt にコピーして継続する。弱結合錯体には外部 TRIC(geomeTRIC または pysisyphus tric)を選べるようにする。
4. **Hessian 再利用**: saddle の前に `seed_hessian/hfauto_job.hess` と `.movecs` を `saddle_frequency/` にコピーして `inhess 2` だけにする(`task dft hessian` を出力しない)。再利用の可否は入力 xyz の sha256 と method fingerprint の一致で判定する。
5. **分散と溶媒**: `disp vdw 4`(D3BJ)を全経路で許可し、既定にする。`solvation: {model: none|cosmo|smd, solvent, dielectric}` を第一級の設定にし、NWChem(`cosmo; dielec 78.4` / `do_cosmo_smd true`)、ORCA(CPCM/SMD)、xtb(`--alpb water`)、CREST、ReaDuct に一貫して渡して、`METHOD_KEYS` に含める。
6. **化学モデルのゲート**: キャンペーンごとに「気相(分光と比較)」か「凝縮相(エッチング: COSMO/SMD + 明示的な H2O/HF 1〜3 分子)」かを決める。PT 状態の判定は ωB97X-D/def2-TZVPD の再最適化か、DLPNO-CCSD(T)/aug-cc-pVTZ の SP で確認する。エネルギーは DZ 構造上の def2-TZVPD SP で報告する。
7. **虚振動の方針(3 段)**: ν < −50 cm⁻¹ は鞍点として扱い、変位して再最適化する(最大 2 サイクル)。−50 ≤ ν < −10 は 1 回変位し、残れば「軟モード flag 付きの極小」として受理する。−10 ≤ ν < 0 はノイズとして flag を付けて受理する。これは minimum-mode-follow を標準フローに組み込むことで実現する。
8. **preopt で棄却された化学種**: dft-minima の既定を `dft_on_preopt_rejected: skip` にする(CH-06)。
9. **ORCA**: 削除するか、SP 専用の参照バックエンド(`! DLPNO-CCSD(T) aug-cc-pVTZ aug-cc-pVTZ/C TightSCF TightPNO`、または `! wB97X-D4 ma-def2-TZVP DefGrid3`、%maxcore は 0.75×RAM/nprocs、MOREAD)として作り直すかのどちらかにする。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| PBE0-D3/def2-SVPD(共有結合の異性化: HCN/HONO/NH3) | **有効** | HCN の ΔE‡ 46.7 / ΔE 13.2 kcal/mol が CCSD(T) の 44.6±1.0 / 14.4±1.0 と 1〜2 kcal/mol 以内で一致 |
| 同上(気相 amine·(HF)n の PT 状態判定) | 限定的 | 非局在化誤差、DZ 基底、ZPE が障壁を上回る。高レベルでの確認が必須 |
| grid xfine / SCF 1e-8(全段) | 低価値 | 結論は変わらず、コストは 2〜4 倍。軟モード確認時だけ有効 |
| 固定構造の method/grid パネル | **有効** | 0.05 h 未満で精度の要否を定量化した |
| xTB 緩和スキャン + DFT SP 数点 | **有効(最も費用対効果が高い)** | 0.1〜0.2 h で、約 49 h の DFT と同じ結論を出した |
| ORCA バックエンド(現行) | 低価値 | 未使用・未検証・古い文法 |
| SiO2 の B3LYP-D3/3-21G COSMO SP(xTB 構造上) | 限定的 | 極大の位置が DFT と一致しない。基底が小さすぎる |

### 3.3 TS 探索(NEB / String / Dimer / RS-I-RFO / NWChem saddle)

**使用ライブラリと仕様上のポイント**
- **NWChem 7.2.3 NEB**(`neb_drv.F`、`neb_utils.F`): `impose` は既定 false で、有効なときだけ endgeom を整列する。algorithm 3 は QN 固定点と damped Verlet の切り替え。Gmax/Grms/Xmax/Xrms の 4 条件ですべて収束を判定する。CI 風の `neb_gradient_get1` はソースにあるが **neb_drv から呼ばれず、入力から有効化できない**([neb_utils.F](https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/neb/neb_utils.F)、[neb_drv.F](https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/neb/neb_drv.F)、[マニュアル](https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html))。
- **NWChem ZTS**(`string.F`): l.23 のコメントは「RMS 座標変化による収束判定」、l.879 は `if ((xrms.lt.tol).and.(xmax.lt.tol))` で、**gmax は収束判定に使わない**(マニュアルの「勾配判定」という記述と実装が食い違う)。停滞すると damped Verlet に切り替わる([string.F](https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F))。
- **pysisyphus 1.0.0**: QCEngine calculator の `get_molecule` は `fix_com/fix_orientation` を指定しない。root の選択は `hessian_ref`(cart)か `rx_modes/prim_coord`(内部座標)でしかできず、それ以外は root 番号の最低固有値を追う。cart 座標では Hessian を射影しない。`trust_min` の既定は 0.1。TS 最適化は内部座標(redund/dlc)で行い、COS の後は HEI 接線との重なりで初期モードを選ぶことが推奨されている([tsoptimization](https://pysisyphus.readthedocs.io/en/latest/tsoptimization.html)、[chainofstates](https://pysisyphus.readthedocs.io/en/latest/chainofstates.html)、[QCEngine.py](https://github.com/eljost/pysisyphus/blob/master/pysisyphus/calculators/QCEngine.py))。
- **QCEngine 0.50.0 / QCElemental 0.50.4**: NWChem の harvester は `fix_com and fix_orientation` が両方 True のときだけ勾配と Hessian を入力フレームに戻す(harvester.py:1210-1226)。NWChem の geometry 文字列に `nocenter/noautosym` は付かない([harvester.py](https://github.com/MolSSI/QCEngine/blob/master/qcengine/programs/nwchem/harvester.py)、[Molecule model](https://molssi.github.io/QCElemental/model_molecule.html))。
- **ORCA NEB-TS**(参考): IDPP、四元数による整列、CI の自動開始(0.02 Eh/bohr)が既定([ORCA 6.1 NEB](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/neb.html))。
- 方法論: CI-NEB([Henkelman 2000](https://doi.org/10.1063/1.1329672))、IDPP([Smidstrup 2014](https://doi.org/10.1063/1.4878664))、geodesic 補間([Zhu 2019](https://doi.org/10.1063/1.5090303))、growing string([Zimmerman 2013](https://doi.org/10.1021/ct400319w))、pysisyphus([IJQC 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26390))。

**現状の設定**: `double_ended_path` は [nwchem_neb, orca_nebts] で、m3 の engine_order は NEB が先頭。NEB は nbeads/maxiter/stepsize/xyz_path だけで、impose なし。初期経路は二面角回転か結合移動の弧で、該当しなければ NWChem の Cartesian 線形補間。RS-I-RFO は cart、root 0、`hessian_init: calc`、`thresh: gau_tight`。NWChem saddle は trust 0.1〜0.15、sadstp 0.03〜0.05、inhess 2、tight(現行の m3/tma 設定は default)。

**問題点**

- **CH-01【致命】pysisyphus+QCEngine+NWChem で勾配と Hessian の座標系が入力と一致しない**(C3-1/C4-1/C2-4 確認)
  - 根拠: `rsirfo_003` では入力フレームと NWChem フレームが 82.9° 回転していた。`-forces[0]` は NWChem フレームの勾配と桁まで一致する。|H·R| は入力フレームで 0.05〜0.14、NWChem フレームで約 0.002。pysisyphus が追った root0 と正しいフレームの root0 の重なりは 0.126、Newton ステップの cos は 0.71。
  - IRC でも |H·R|/|R| が HCN [1.84, 0.20, 1.14]、HONO [0.10, 0.03, 0.68]、NH3 [0.001, 0.001, 0.36]。HCN の遷移ベクトルは −84.4 cm⁻¹(真値 −1131.6)。
  - NWChem の標準配向は構造ごとに決まるので、cycle 間でフレームが一定である保証もなく、Bofill 更新の勾配差分も不整合になる。
  - 訂正: 影響の大きさは系による。HONO の遷移ベクトルは −676 cm⁻¹ でほぼ正しく、勾配ゼロの点はフレームに依存しない。それでも軌跡は IRC として正当化できない。
- **CH-02【高】QCEngine 経由の PES が一致せず、provenance も誤記している**(C2-4、C3-4 一部訂正、C4-3 確認、見落とし追加)
  - 根拠: `_qcengine_keywords` は `dft__disp` と `basis__spherical` しか渡さないので、grid は medium、SCF は 1e-6、memory は 1 GB 固定になる。1 回の RS-I-RFO の中でも、初期の exact Hessian は fine/1e-7、勾配は medium/1e-6 と PES が混在している(rsirfo_002/003 で medium 33、fine 4)。IRC の成果物の method には `grid: xfine`、`scf_energy_tolerance: 1e-8` が 32 箇所記録されているのに、実際は medium/1e-6 で計算されていた。`_nwchem_qcengine_evidence` は grid、収束、フレームを検査しない。
- **CH-06【高】障壁の存在を事前に判定せず、単調経路を端点の線形補間ブラケットに回している**(C3-2 一部訂正、C4-7 確認、C5-6 一部訂正、C2-9 一部訂正)
  - 場所: `chemistry/path_diagnostics.py:247-257`、`workflow/reaction_state.py:7-8`、`backends/ts/nwchem_saddle.py:448-560`、`stages/recover_path.py:72-100,236-242`、`stages/dft_minima.py:177-191`
  - 根拠:
    - 「単調かつ distinct_basin」は `bracket_narrow_saddle` と診断され、収束した string 経路を捨てて Kabsch 整列した線形補間上の一点計算に回される。+51 kJ/mol の人工極大ができ、seed の虚振動は 4 本になる。
    - recover-path は未収束経路(NEB の iteration 1)を受理し、`basin_assessment` がないときの既定値が `distinct_basin/accepted`(fail-open)。
    - NT2 の「TS guess なし」22 件はすべて単調上昇(+12〜247 kJ/mol)で、生成物 basin がないことを示す陰性結果なのに、ゲートとして使われていない。
    - preopt で棄却された化学種(IP の崩壊)は、元の seed のまま dft-minima に流れる。
  - 訂正: 最終的な TS 判定(周波数、エネルギー順序)が偽 TS の採用を防いでいるので、問題の本質は「科学的に誤った結果」より「無効な計算とコスト」にある(重大度は高)。`validate_reaction_endpoint_pair` は宣言された bond_changes の方向と距離は検査している。欠けているのは共有結合グラフの変化と ΔE 窓。
- **CH-08【高】NWChem ZTS の偽収束**(C3-3 確認、見落とし追加)
  - 根拠: `m3_tma_path_ensemble_001` の images_11 では、gmax 9.14e-4 → 9.15e-3 と悪化した後に damped Verlet へ切り替わり、xrms 4.2e-7 で「converged」になった(images_13 も同様)。`parse_string_optimization_history` は `@zts` 行だけを読み、実際に出力されている `string: gmax,grms,xrms,xmax=` 行を無視している。このため、偽収束した 2 本を「独立な収束経路」と数えて `barrierless_at_resolution` を採択した。停滞判定は未収束のときしか働かないので、偽収束が停滞検知も無効にする。tol の既定値 4.5e-4 は勾配の閾値を流用したもので、変位の閾値としての意味がない。
- **CH-09【高】NWChem NEB が impose・IDPP・CI なしで先頭エンジンになっている**(C3-7 確認、A6-9 一部訂正)
  - 根拠: NEB の出力 12 本のうち「converged」は 3 本で、**すべて NH3 反転**。実反応 9 件の収束は 0。HCN では HNC が回転した向きのままで原子間距離 0.529 Å まで衝突し、150 反復後も中央 bead が +3.2 Eh。HONO v3 は 9,118 s のうち最後の 58 反復がステップ 0 で停滞し、CI なしの最高 bead(8.3 kcal/mol)は真の TS(13.6)を 5.3 kcal/mol 過小評価していた。停滞は事後にしか判定されない。
- **CH-10【高】RS-I-RFO が root 0 固定・cart・未射影 Hessian で動いている。endpoint dimer は価値が低い**(C3-5 確認)
  - 根拠: `rsirfo_003` で「Initial Hessian has 4 negative eigenvalue(s)」。NWChem の射影済み Hessian では虚振動 1 本(−110 cm⁻¹)。以後の重なりは 0.14→0.05→0.09。dimer は極小から 0.15 Å 変位した点から始め、smoke 1 は cycle 0 で「Converged!」、smoke 2 は曲率が正のまま。seed ゲートは射影済み Hessian(n_imag=1)で判定し、最適化は未射影 Hessian(負固有値 4 本)で行うという不整合がある。
- **CH-33【中】saddle が極小へ崩壊しても分類せず、ノイズ以下の prominence も許している**(C3-8 一部訂正): 崩壊(E < max(E_R, E_P)、または refresh Hessian で n_imag=0)を独立した診断として扱わず、「経路のやり直し」に回している。prominence の下限は method で 1e-8 Eh(SCF 収束より小さい)まで下げられる。(訂正: 自動再開の既定は 1 回だけで、trust を 4 段階縮めたのは主に手動 run)
- **CH-34【中】虚振動モードの重なり 0.5 を hard gate にしている**(C3-9 確認): HCN の正しい TS が 0.547 で一度棄却され、閾値だけを変えて同じ計算をやり直した(`hcn_saddle_v3`)。逆に TMA·(HF)2 では、HF 再配向の軟モード(−118 cm⁻¹、重なり 0.787)が「反応モード」として選ばれた。
- **CH-37【低】TS attempt ごとの timeout が stage の walltime を無視する**(pysisyphus 系の既定 172,800 s。recovered_path は上限 14,400 s に対して string だけで 12,089 s)。`PysisyphusEngine.search_ts` は NWChem NEB に委譲する(C3-10 確認)。

**推奨 TS プロトコル**(C3-11 を統合。数値は初期値)

- **P0 事前ゲート**(秒〜数十分):
  1. 端点は同一 PES(fine、SCF 1e-7)で最低実振動数 ≥30 cm⁻¹ とし、満たさなければ ±0.1 Å の mode-follow で再最適化を 1 回行う(minimum-mode-follow を利用)。
  2. 共有結合グラフの変化(ヒステリシス判定)を必須とし、配座変化だけなら TS 探索の対象外にする。
  3. プロトン移動は q=r(D–H)−r(A–H) の GFN2 緩和スキャン 15 点と DFT SP 5 点で、内部極大 ≥1 kcal/mol の場合だけ先へ進める。
  4. ΔE_rxn ≤ +40 kcal/mol。
  5. xTB 段の陰性結果(NT2 の単調上昇、preopt での IP 崩壊)があれば、明示的な override なしに DFT の TS 探索を起動しない。
  6. 凝縮相の問いには最初から ALPB/COSMO を同じ PES に入れる。
- **P1 低レベル経路**(座標系を修正した pysisyphus + GFN2、数分): `cos: {type: gs, max_nodes: 9, climb: True, climb_rms: 5e-3, reparam_every: 2}`、`geom: {type: dlc}`、`interpol: idpp`、`tsopt: {type: rsirfo, hessian_init: calc, thresh: gau}`。代替として ReaDuct B-spline(CH-30)。
- **P2 DFT 経路**(NWChem のみで行う場合): string を `nbeads 9, stepsize 0.05, interpol 3, impose, tol 2e-4`、`maxiter 20` × 2 回で区切り、**最終 gmax ≤1.0e-3 Eh/bohr かつ gmax が 3 反復連続で悪化しない**ことを収束の条件にする。基底は def2-SVP(中性系は diffuse なし)、grid fine。NEB を使う場合は必ず `impose` と IDPP/geodesic の `xyz_path`(最小原子間距離 ≥0.7 Å を QC で確認)を付け、maxiter 30 ずつ区切る。Xmax=0 が 5 反復続くか、10 反復で Gmax の改善が 5% 未満なら打ち切る。HEI は spline で決める。
- **P3 saddle**: NWChem driver は default 閾値、trust 0.1、sadstp 0.1、`inhess 2`(キャッシュした `.hess` をコピーし、再計算しない)、maxiter 50。pysisyphus RS-I-RFO は geom `redund`/`dlc`(2 断片以上は `tric`)、`hessian_init` に Eckart 射影したキャッシュ h5、root は経路接線との重なりで選び、trust 0.1/0.01/0.3、`thresh: gau`、max_cycles 50。崩壊を検出したら再試行せず、生成した極小を registry に登録して basin を再評価する。prominence の下限は max(2e-5 Eh, 20×SCF 許容値)、seed として採用する障壁は 0.5 kcal/mol 以上。
- **P4 判定**: fine で周波数を計算し、虚振動 1 本(|ν| >50 cm⁻¹、ねじれは >20 cm⁻¹)、エネルギー順序、両方向の IRC 接続で確定する。重なりは順位付けと記録にだけ使い、0.3 未満なら警告にする。
- **予算**: 経路 1 回 + saddle 2 回、1 反応 6 h まで。超えたら `unresolved_within_budget`。attempt の timeout は min(設定値, stage の残り walltime)。
- **注意**: C3-1 を直すまでは、pysisyphus 経路(P1/P3 の RS-I-RFO、P4 の IRC)は使わない。信頼できるのは「NWChem のみ(gmax 監視付き string → NWChem saddle)」の経路だけ。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| NWChem string(impose, interpol 3)→ NWChem saddle(inhess 2) | **有効** | HCN で 849+47 s、正しい TS(−1131.6 cm⁻¹) |
| 未収束 NEB の HEI → saddle | 限定的 | HONO では TS を得たが、NEB に 2.5 h かかり、障壁を 5.3 kcal/mol 過小評価 |
| NWChem NEB(現行設定) | 低価値 | 実反応 9 件で収束 0、原子衝突、停滞 |
| pysisyphus RS-I-RFO / dimer(QCEngine 経由) | 低価値(修正まで**無効**) | 座標系と PES が一致しない |
| 端点の線形補間ブラケット | 低価値 | 人工的な極大しか作らない |
| TS 探索全般(気相 amine·(HF)1,2、FHF⁻/F(HF)n⁻、HF/メチルの配向変化) | 低価値 | 障壁がないか ZPE 以下。スキャン・振動解析・配座アンサンブルで扱うべき |
| TS 探索(HCN/HONO の異性化、SiO2 エッチ) | **有効** | 結合の組み替えを伴う実在の障壁がある |

### 3.4 IRC・経路検証・端点接続

**使用ライブラリと仕様上のポイント**
- **pysisyphus EulerPC IRC**: 分岐ごとに `energy_increased`、`converged`、「Displaced geometry is higher in energy compared to TS!」を出力する。`energy_thresh` の既定は 1e-6 Eh。`hessian_init` がなければ TS の Hessian を再計算する([irc](https://pysisyphus.readthedocs.io/en/latest/irc.html))。
- **NWChem MEPGS**: `task dft mepgs`(Gonzalez–Schlegel)は同一 PES・同一フレームで動き、`inhess 2` で saddle_freq の `.hess` を再利用できる([mepgs](https://github.com/nwchemgit/nwchem-wiki/blob/master/mepgs.md))。
- IRC の終端は極小ではないので、制約なしで最適化して確認するのが標準的な実務である([Maeda らのレビュー](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24757)、[Hratchian & Schlegel](https://pubs.acs.org/doi/10.1021/ct0499783))。QRC(虚モード方向に ± 変位して最適化)は接続確認として安価([Goodman & Silva](https://www.sciencedirect.com/science/article/abs/pii/S0040403903021798))。
- 結合判定の既定係数は RDKit の `DetermineConnectivity` で covFactor 1.3([rdDetermineBonds](https://rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html))。

**現状の設定**: pysisyphus を QCEngine 経由で使い(現行の全 production パイプライン)、EulerPC(max_cycles 80〜150、step 0.08〜0.1、rms_grad 5e-4〜1e-3)、端点は NWChem native で再最適化する。端点は順序付き RMSD 0.35〜0.75 Å で basin に割り付け、両向きとも合格した場合は direct を採用する。`render_irc_input` は 6 キーだけを受け付け、`hessian_init`、`displ`、`energy_thresh` は設定できない。

**問題点**

- **CH-03【致命】IRC の合否が分岐の下降を見ておらず、端点最適化後の RMSD だけで決まる**(C4-2 確認、見落とし追加)
  - 根拠: `hcn_irc_v5` では、前進側が 54 step で TS 比 +190.68 kJ/mol、後退側は step 0 で +4.33 kJ/mol となり、両方とも「Energy increased!」で停止した。それでも `irc_validated=True`、`endpoint_rmsd_A=2.27e-5`。HCN の PES 上の極小は 2 つしかないので、どこから最適化しても合格する同語反復の検証になっている。NH3_004 の後退側も step 0 で −2.29 kJ/mol のまま合格した。`hcn_irc_v5` は `scripts/run_pysisyphus_irc.py` 経由で作られ、IRC stage の lineage 検査も迂回している。
- **CH-04【高】IRC 端点の basin 割付けが緩く、HONO の cis/trans を取り違えている**(C4-4 確認、見落とし追加)
  - 根拠: HONO v3 で 4 通りの割付けすべてが閾値 0.75 Å で合格し(0.639 Å)、`if direct_ok or not swapped_ok` によって誤った向きを採用した。distance spectrum 0.308 Å は判定に使われていない。NH3 反転の極小同士の RMSD は 0.326 Å で、0.35 Å の閾値でも両向き合格になる。さらに、宣言された反応座標(二面角の係数 +1、trans→cis)は符号付き Δφ で判定するため、正しい向き(Δφ≈−π)は必ず却下される。定義と判定式が矛盾しており、**取り違えた向きのときだけ合格する構成**になっている。
- **CH-05【致命】`reaction-plan` の `basin_assessment` から必須キーが抜け、production の `reaction-classify` が全件 unresolved になる(回帰)**(C4-5 確認)
  - 根拠: `reaction_plan.py:161-182` の `endpoint_evidence` は 3 キーしか持たないのに、`reaction_classification.py:282-291` は `real_qm_executed/fallback_dummy/method_evidence_validated` を要求する。HONO v3 は `unresolved`(`frequency_validated_distinct_endpoint_basins_missing`)。旧版の reaction-plan は 12 キーを出していた(`m3_tma_reaction_classify_seeds_001`)ので、今回のリファクタで生じた回帰である。テストの fixture にこの evidence がなく、`require_real_qm=True` のテストは unresolved を期待値にしているため、検出できなかった。
- **CH-35【中】置換等価な縮退素反応を same_basin に併合している**(C4-9 確認): NH3 反転の反応物と生成物は graph RMSD 1.8e-16、ΔE 0 で same_basin と判定される。一方で TS(−757i)は検証済みである。探索段階で置換等価を除外するのは仕様上の選択だが、basin 同一性の判定に暗黙に埋め込まれていて、契約として明示されていない。
- **CH-36【中】IRC 検証コードの重複・デッドフラグ・HF のハードコード**(C4-10 一部訂正): `allow_legacy_rmsd_only_irc_validation` は設定されるだけで、どこからも読まれない。`reaction_path_qc.py:600-667` は `ion_pair/neutral_complex` をハードコードしている。descriptors の `HF_fragment_charge/dipole_D` は常に None。(訂正: `render_irc_input` のテストは存在するが、grid/SCF を検査しない不完全な契約をそのまま固定している。独立端点 QC は nwchem_neb/dummy に対しては意味がある)

**推奨設定値・プロトコル**
1. **IRC の実装**: 修正するまでは pysisyphus+QCEngine を IRC に使わない。修正するなら、QCEngine calculator を継承して `fix_com=True, fix_orientation=True` を渡す driver(約 20 行)に替え、QCEngine のキーワードは method から生成する(`dft__grid`、`dft__convergence__energy`、`dft__convergence__density`、`dft__iterations`)。フレームガードとして、入力座標での |H·R_i| >1e-4 Eh/bohr² か、遷移ベクトルの波数が NWChem 値から 5% 以上ずれたら fail にする。より単純な代替は NWChem の `task dft mepgs`(`inhess 2`、stride 0.1、maxmep 60、opttol 3e-4)。接続の確認だけが目的なら、QRC(最大原子変位 0.05〜0.1 Å か ΔE 約 0.5 kcal/mol の ± 変位を NWChem で最適化)を既定にする。
2. **合否判定**(`chemistry/irc_validation.py` に純関数として 1 本化):
   - (a) 両分岐とも初期変位後の実際のエネルギー低下 >0(望ましくは ≥0.3×displ_energy)
   - (b) 各分岐が 5 step 以上進み、単調に下降している(許容 2e-6 Eh)。終了理由は `rms(grad) converged` か `Energy converged` に限り、5 step 未満の `Energy increased` は fail
   - (c) 分岐終点が E_TS − max(2 kcal/mol, 0.5×(E_TS − E_その側の極小)) 以下
   - (d) (a)〜(c) を満たした終点だけを最適化し、IRC による下降量と最適化による下降量を記録する
3. **端点の割付け**: 最適化済みの端点を、登録済み極小に「graph RMSD ≤0.05 Å かつ |ΔE| ≤5e-5 Eh」で一意に割り付ける(次点の RMSD が best の 3 倍以上、または 0.1 Å 以上離れていること)。forward と backward は別の basin でなければならず、曖昧なら `ambiguous_endpoint_assignment` で fail。ねじれ角は符号付き Δφ ではなく状態窓(cis: |φ| <30°、trans: |φ| >150°)で判定する。
4. **スキーマの一本化**: `basin_assessment` は pydantic の `BasinAssessment/EndpointEvidence`(必須キー付き)とし、reaction-plan は `basin_identity._minimum_evidence` を呼ぶ。HONO v3 の実データを golden として、reaction-plan → reaction-classify を `require_real_qm=True` でつなぐ統合テストを追加する。
5. **縮退反応**: 原子写像付きの basin 同一性と構造的同一性を分け、検証済み TS と IRC でつながる置換等価な極小の組は `degenerate_rearrangement`(ΔG_rxn=0、経路の縮重度を速度に反映)にする。
6. **golden テスト**: `hcn_irc_v5` は fail、`hono_v3` は orientation が曖昧なので fail、`same_basin_water_v1` は pass、NH3_004 の後退側は fail を期待値にする。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| pysisyphus IRC(QCEngine 経由、現行) | 低価値(**無効**) | 座標系と PES の不一致に加え、合否判定が同語反復 |
| IRC 端点を同一 PES で NWChem 再最適化 | **有効** | 方針は正しい。分岐判定を前段に加える必要がある |
| 同一停留点の判定(graph RMSD ≤0.02 Å かつ ΔE ≤1e-5 Eh) | **有効** | same_basin_water で RMSD 8e-6 Å、ΔE 2.6e-10 Eh |
| minimum-mode-follow / assess | **有効(最も決定的)** | TMA·(HF)2 の same_basin を確定させた。前段に移すべき |
| path-ensemble(多解像度) | 限定的 | 偽収束経路を数えている。gmax 判定を直せば有効 |
| 気相 H 結合錯体・HF2⁻ の IRC | 低価値 | TS が存在しないか人工物。端点の軟モード検証と下り坂接続テストを優先すべき |

### 3.5 反応探索(SCINE ReaDuct / NT2 / AFIR / 緩和探索)

**使用ライブラリと仕様上のポイント**
- **SCINE ReaDuct 6.1.0**: `HessianTask` は `NormalModeAnalysis::calculateNormalModes` で並進・回転を除去した質量加重の振動数を出力し、虚振動モードごとに vibmode ファイルを書く。`TsOptimizationTask` は `automatic_mode_selection` を与えたときだけ反応原子の寄与で追従モードを選び、与えなければ mode 0。出力は相対パス(CWD)へ書く。`run_bspline_task`(double-ended、`align_structures`、`extract_ts_guess`)がある([HessianTask.h](https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/HessianTask.h)、[TsOptimizationTask.h](https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/TsOptimizationTask.h)、[NtOptimization2Task.h](https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/NtOptimization2Task.h)、[manual v6.1.0](https://scine.ethz.ch/static/download/manuals/v6.1.0/readuct_manual.pdf))。
- **SCINE Utilities 10.1.0**: NT2 は Wiberg 結合次数で判定し(形成 >0.75、切断 <0.15)、Savitzky–Golay 平滑化の後に極大がなければ「No transition state guess」を投げる。AFIR の既定は γ=1000 kJ/mol、phase_in 30。`normal_modes.calculate(...).get_wave_numbers()` がある([NtOptimizer2.cpp](https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Utils/GeometryOptimization/NtOptimizer2.cpp)、[AfirOptimizerBase.h](https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Utils/GeometryOptimization/AfirOptimizerBase.h)、[NormalModesPython.cpp](https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Python/NormalModesPython.cpp))。
- **xtb_wrapper 3.0.2**: `solvation/solvent`(既定は気相)、`electronic_temperature`(300 K)、`max_scf_iterations`(100)を公開している([XtbSettings.cpp](https://github.com/qcscine/xtb_wrapper/blob/3.0.2/src/Xtb/Xtb/Wrapper/XtbSettings.cpp))。
- **参照実装 Puffin 2.1.0 / Chemoton 2.0**: 各ジョブを専用の work_dir に chdir し、`automatic_mode_selection` を NT の結合対から設定し、IRC 端点が出発構造と一致するかを検証する([scine_react_job.py](https://github.com/qcscine/puffin/blob/2.1.0/scine_puffin/jobs/templates/scine_react_job.py)、[Chemoton 2.0](https://arxiv.org/abs/2202.13011))。Chemoton 2.0 は、線形分子(HCN)の NT は力が回転軸上にかかるため難しいと明記している。
- **AFIR**: 熱反応の γ ≈ −RT ln(h/(t k_B T))(298 K、t=10 日で約 107 kJ/mol)([GRRM17](https://pmc.ncbi.nlm.nih.gov/articles/PMC5765425/))。autodE の `min_imag_freq` は −40 cm⁻¹([config](https://github.com/duartegroup/autodE/blob/master/autode/config.py))。

**現状の設定**: `hessian_negative_threshold=-1e-6` で未射影の固有値を数える。`nt_total_force_norm 0.1`、各 `max_iterations 200`、`afir_energy_allowance 300 kJ/mol`、TS 最適化は bofill→evf→dimer の順(モード選択なし)。trial は 1 化学種あたり 2 件(pilot の設定。コードの既定は 12)。xTB は気相。

**問題点**

- **CH-26【致命】低レベル TS の虚振動数を未射影の Cartesian Hessian 固有値で数え、真の 1 次鞍点を全件棄却している**(C5-1 確認)
  - 根拠: pilot の stdout にある「Rot. and trans. vib. removed」8 ブロックは、すべて負の振動数が 1 本(−1173.2、−506.7、−256.7、−1265.3、−601.9、−239.2、−246.8、−513.1 cm⁻¹)。一方 hfauto の記録は 2〜4 本。全 run 58 件で count=1 は 0 回で、IRC は一度も実行されていない。原因は並進・回転を射影していないことにある(質量加重は合同変換なので負固有値の数を変えない)。
- **CH-27【高】TS 最適化にモード選択がなく、IRC 端点が出発構造と一致するかも検証していない。低レベル端点の極小性も検証していない**(C5-2 確認、見落とし追加)
  - 根拠: TMA·(HF)2 の ts_bofill は −24.408736 Eh(極小 +2.6 kJ/mol)で、vibmode-00001 で大きく動くのはメチルの H(H5, H4, H12, H7)。つまりメチル回転の鞍点だった。この同じエネルギー点が tma_hf2_explore_v11 の生成物と m3 の seed 1ec1 にも現れ、DFT に流れる経路が開いている。C5-1 を直すとこの問題が表に出る。
- **CH-28【高】低レベル候補にエネルギー窓も化学的妥当性のゲートもない**(C5-3 確認、見落とし追加)
  - 根拠: PYR_HF1 の置換で環内の C–N を脱離基にした候補は、xTB で生成物 +183.4 kJ/mol、TS +196.2 kJ/mol だったのに、DFT(+224.7)まで昇格して 7.3 h を消費した。TS を棄却した後も、NT2 の強制到達構造を緩和したものが `success=True` で候補になる(pilot の NT2 候補 2 件はどちらも TS 検証なし)。AFIR の候補 3 件は `afir_converged=False` だった。(訂正: dft-minima は同一 PES 内で低レベルのエネルギー昇順に並べ、`candidate_endpoint_limit` を適用している。問題は絶対的な窓がないこと)
- **CH-29【中】trial 生成が role を無視し、協奏リレーの優先度が低い**(C5-4 一部訂正): 列挙 381 件のうち relay は 116 件で、優先度は 18−d(分子間移動は 20−d)。build-complexes が記録している role(base/acid/relay)を使っていない。coverage は列挙数・切り捨て数・機構別の内訳を記録していない(`trial_budget_saturated` の真偽値だけ)。(訂正: 2 件への切り捨ては pilot の設定値で、逆向き移動は 30 件中 2 件)
- **CH-30【中】NT2 は線形分子とねじれ異性化に向かない**(C5-7 一部訂正): HCN ではフレーム 0〜12 が軸上にとどまり(偏差 2e-8 Å)、C–N が 0.695 Å まで圧縮されて +1.53 Eh 上昇し、TS guess が破綻した。ReaDuct の `run_bspline_task` は使っていない。(訂正: HONO は NT2 ではなく明示的な二面角 trial で扱われている)
- **CH-31【中】ReaDuct/pysisyphus の出力が CWD(リポジトリ直下)に漏れ、一次証拠が失われる**(C5-8 確認): リポジトリ直下の `afir_biased/`、`nt2_guess/`、`product_unbiased/`、`ts_bofill/`、`ts_hessian/`、`internal_coords.log` は git 未追跡で、別の run(15 原子と 19 原子)の出力が混在して上書きされている。attempt ごとの timeout もない。
- **CH-32【低】AFIR の γ=300 kJ/mol と単一原子ペア**(C5-10 一部訂正。反復上限到達は 10 件)。
- 見落とし: `explore()` の包括的な `except Exception` が、コードの不具合(KeyError など)を化学的な「探索失敗」に変換してしまう。NT2 の抽出基準などの設定を config から渡せない。

**推奨設定値・プロトコル**
1. 虚振動: `modes = scine_utilities.normal_modes.calculate(hessian, atoms)`、`wn = modes.get_wave_numbers()` とし、`wn < −50 cm⁻¹`(設定可能)で数える。全振動数を attempt data に保存する。
2. TS 最適化: `automatic_mode_selection = sorted(associations ∪ dissociations の原子)` を渡す。収束後、虚モードでの反応原子の寄与が ≥0.5、または ± 変位で active bond が ≥0.3 Å 変わることを確認する。IRC の片端が出発構造とグラフ同型かつ置換不変 RMSD <0.1 Å で一致することを要求する(一致しなければ `irc_not_connected_to_source`)。
3. エネルギー窓: ΔE‡_xTB ≤150 kJ/mol かつ ΔE_rxn_xTB ≤100 kJ/mol を既定にする(設定可能)。置換の脱離基は、`fragments_after_bond_cut` が None でない非環結合に限り、芳香環の原子は除外する。TS が棄却された NT2 回収構造と未収束 AFIR の構造は、既定では DFT に昇格させない。
4. 陰性結果: NT2 の単調上昇は `monotonic_uphill_no_product_basin` として ΔE_stop を記録し、最終フレームも緩和する。source ごとの判定(`single_low_level_basin` / `exergonic_barrierless` / `candidates`)を DFT の前段ゲートにする。endpoint-seed-screen は、既知 basin との ΔE <0.05 kcal/mol かつ RMSD <0.05 Å の seed と、虚振動 < −20 cm⁻¹ の seed を除外する。
5. trial: donor は acid/relay、acceptor は base/relay に限り、n_HF ≥2 では relay の優先度を単純移動以上にする。`max_trials_per_species` は 8〜12(xTB 1 attempt 約 5 s)。
6. double-ended: 端点が既知の反応は `run_bspline_task`(GFN2、num_control_points 5〜8、num_integration_points 21〜41)を実行し、単調なら DFT NEB を起動しない。線形分子の NT2 には 0.05 Å のランダム変位か 10° 以上の折り曲げを加え、`nt_extraction_criterion='first_maximum'` を使う。
7. 計算器の設定: `load_system_into_calculator(..., max_scf_iterations=300, electronic_temperature=300, solvation='alpb', solvent='water')` を渡せるようにし、SCC が失敗したら 1000 K で 1 回だけ再試行する。溶媒は request fingerprint に含める。SiO2 は Q3 Si クラスタ + H2O 1〜2 個で、4 ペア同時駆動の relay trial を行う([Knotter 2000](https://pubs.acs.org/doi/abs/10.1021/ja993803z))。
8. AFIR: γ=100〜150 kJ/mol から始め、生成物が得られない場合だけ 300 に上げる。バイアスはフラグメント単位(MC-AFIR に近い形)にする。
9. 出力: `contextlib.chdir(attempt_dir)` で実行し、SCINE の Log を attempt に保存する。1 attempt を subprocess で動かして timeout と並列実行を可能にする。リポジトリ直下のディレクトリは削除し、`.gitignore` に追加する。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| ReaDuct NT2 低レベル探索 | **有効**(修正前提) | 59 attempt で約 4.5 分。正しい化学的結論(n=1,2 はイオン対 basin なし、n=3 で移動)を出していた |
| ReaDuct TS 最適化 + IRC | 限定的 → 修正後は有効 | 数え方のバグとモード選択の欠如で、現在は機能していない |
| AFIR | 限定的 | 単一ペアで γ が高く、未収束の昇格もある |
| relaxation-discovery | **有効** | TMA·(HF)3 のイオン対を自発的に生成 |
| NT2(線形分子・ねじれ) | 低価値 | 距離ドライブが原理的に合わない。B-spline やスキャンを使う |
| discovery-audit / coverage | 限定的 | 有限探索であることの明示は誠実だが、網羅性の指標になっていない |

### 3.6 熱化学・速度論・校正

**使用ライブラリと仕様上のポイント**
- **GoodVibes 4.3.0**(Linux)/ **3.2**(Windows): `-q` は Grimme qRRHO エントロピー + Head-Gordon qh-H(cutoff は `-f`)。`-v` は振動と ZPE の両方をスケールし、`--zpe-vscal` で ZPE だけ別に指定できる。ほかに `--symm`(pymsym)、`-c`/`--media`(標準状態)、`--temp`/`--ti`、`--invert`(既定 −50 cm⁻¹)、`--spc` がある。**NWChem パーサ `parse_nwchem_thermo` はファイル中の全 `P.Frequency` 行を連結し、ZPE とエネルギーは最後の値を使う**(io.py:1431-1441)。3.2 は `--csv` が引数を取らない([GoodVibes](https://github.com/patonlab/GoodVibes)、[F1000Research 2020](https://doi.org/10.12688/f1000research.22758.1)、[PyPI 3.2](https://pypi.org/project/goodvibes/3.2/))。
- qRRHO([Grimme 2012](https://doi.org/10.1002/chem.201200497))、qh-H([Li et al. 2015](https://doi.org/10.1021/jp509921r))、msRRHO([Grimme & Pracht 2021](https://doi.org/10.1039/D1SC00621E))。溶液の標準状態補正は 1 atm→1 M で Δn あたり 1.89 kcal/mol([Bryantsev 2008](https://doi.org/10.1021/jp802665d))。振動スケール因子は [Truhlar DB](https://comp.chem.umn.edu/freqscale/) と [CCCBDB](https://cccbdb.nist.gov/vibscalejust.asp) を参照。
- トンネル効果: Wigner は小補正の極限でしか有効でない([RMG Wigner](https://reactionmechanismgenerator.github.io/RMG-Py/reference/kinetics/wigner.html)、[OSTI 6485793](https://www.osti.gov/biblio/6485793))。非対称 Eckart([Johnston–Heicklen](https://doi.org/10.1021/j150606a003))。Arkane は NWChem のログを読めない([essfactory](https://reactionmechanismgenerator.github.io/RMG-Py/reference/arkane/essfactory.html))。
- 参照値: F⁻+HF→HF2⁻ 45.8±1.6 kcal/mol([Wenthold & Squires 1995](https://doi.org/10.1021/j100007a034))、アミンの PA(NIST: NH3 853.6、TMA 948.9、pyridine 930.0、aniline 882.5 kJ/mol、[NIST WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664393&Mask=8))、HNC−HCN ΔE0 = 14.83 kcal/mol(van Mourik らの CCSD(T))。

**現状の設定**: 8 パイプラインが `extra_args: -q -f 100 -v 0.985 --temp 298.15`。`settings.frequency_scale_factor` と `quasi_rrho_cutoff_cm1` は GoodVibes に渡らず、二重管理になっている。NWChem は `noautosym`(σ=1)、`temperature 1 T`。GoodVibes が失敗すると、NWChem RRHO に独自の低振動ペナルティを加えた内部フォールバックで続行する。

**問題点**

- **CH-38【致命】GoodVibes が NWChem 出力内の複数の振動ブロックを連結し、HONO の ΔG‡ を約 12 kcal/mol 過大評価したまま production として公開している**(C6-1 確認、見落とし追加)
  - 根拠: `saddle_initial_hessian=true` の saddle_freq 出力には振動ブロックが 2 つある(seed −333.54i で ZPE 0.019317、最終 TS −680.11i で ZPE 0.018747)。Goodvibes.csv の zpe は 0.0375115 で、0.985×両者の和と一致する。thermo_records は ΔE‡ 13.62 に対し ΔE0‡ 24.19、ΔG‡ 24.11 で、ΔZPE +10.6 kcal/mol は物理的にあり得ない。NWChem の最終ブロックから計算し直すと ΔG‡ ≈ 12.2 kcal/mol。この記録は `quality_tier Q4`、`production_rank_eligible=True` で reaction-rank まで流れている。
  - 範囲: runs で起きたのはこの 1 件だけだが、7 パイプラインが `saddle_initial_hessian: true` を使っていて、初回試行は 2 ブロック出力になるので今後も起こりうる。hfauto 自身のパーサは最終ブロックを正しく読んでいるのに、GoodVibes の結果との突き合わせをしていない(`external_zpe` をそのまま採用している)。
- **CH-39【高】内部フォールバック熱化学が誤った値を scientific ランクへ流す**(C6-2 確認、見落とし追加)
  - 根拠: 準線形の HNC が 'A= **********' で線形と判定されず、Rotational entropy = −5.047 cal/mol·K になった。その結果 ΔG_rxn は 18.05 kcal/mol(GoodVibes では 12.74)となり、`scientific_rank_eligible=True`、`confidence 0.82` が付いた。低振動ペナルティ 0.20·(1−ν/100)² は Grimme qRRHO より大幅に小さい(TMA·(HF)2 で 0.19 に対し qRRHO は 1.17 kcal/mol)。`main_values_are_fallback` は `main_values_are_dummy` の別名なので、フォールバックでも False になる(フラグが機能していない)。互換コピーは線形分子の回転定数を書き換えているのに `scientific_values_changed=False` を返し、`re.subn(count=1)` なので複数ブロックの出力では最終ブロックを修正しない。
- **CH-40【中】エネルギー源と熱補正源で理論レベルが混在しうる(潜在)**(C6-3 一部訂正): GoodVibes 経路では、`external_g` に freq 出力の scf_energy が含まれるため、SP のエネルギーは黙って捨てられる(`E_hartree` 列と `G` 列が別レベルになる)。反応内で method が一致しているかの検査もない。(訂正: 現行のパイプラインに sp+thermo の併用はない)
- **CH-41【中】回転対称数 σ=1 固定、1 M/溶媒の標準状態と配座アンサンブルに非対応**(C6-4 確認): NH3 反転は σ_R=3、σ_TS=6 なので RT ln2 = 0.41 kcal/mol 効く(3.42→約 3.83)。SiF4(σ=12)は 1.47 kcal/mol。溶液中の会合反応では Δn あたり 1.89 kcal/mol。
- **CH-42【中】温度依存を扱えず、温度判定が fail-open になっている**(C6-5 確認、見落とし追加): GoodVibes の CSV には温度列がないので、`external_temperature` に NWChem の温度が代入される。NWChem を 373 K にして `--temp 298.15` を残すと、298 K の G が「要求温度での外部 GoodVibes の値」として受理される。キャッシュキーは artifact_id だけ。
- **CH-43【中】kinetics 段は物理的に不正確で、どのパイプラインでも使われていない**(C6-7 確認): 'eckart' を指定すると κ=1 に置き換わる。クロスオーバー温度 T_c は 1131i で 259 K、1700i で 389 K で、HF プロトン移動では 298 K が深いトンネル領域に入るため、Wigner は大きく過小評価する。2 分子反応も 's⁻¹' で表示される。Cantera は h0=s0=cp0=0 の擬似機構、Arkane は雛形だけ(NWChem を読めない)。
- **CH-44【中】calibrate は手法の誤差を評価していない**(C6-8 確認): 任意の重みで DB 支持スコアを合算しているだけ。HNC−HCN の ΔE0 12.83 はベンチマーク 14.83 より 2.0 kcal/mol 低く、DFT 同士のパネル幅(0.98)の外にある。DFT だけのパネルでは真の誤差が見えない。
- **CH-45【低】GoodVibes のバージョンが OS で分かれ(3.2/4.3)、ZPE と振動にスケール因子を共用している。古い CSV が混入しうる**(C6-11 確認、見落とし追加): `parse_output_directory` は workdir の `*.csv` をすべて glob で拾うので、新しい CSV が書かれなかった場合に前回の値を採用しうる。
- **CH-46【低】感度解析が剛直な HCN だけで行われている**(C6-12 確認): cutoff を 50/100/150 cm⁻¹ と変えても差は 0.003 kcal/mol で、低振動モードがないので自明な結果。`--qs` と `--bav` は変えていない。thermo-sensitivity は thermo とは別のソース選択規則を持つ。

**推奨設定値・プロトコル**
1. **ZPE 連結の修正(即時)**: GoodVibes には最終振動ブロックだけを抜き出した thermo view ファイルを渡す(あるいは TS と極小で独立した `task dft frequencies` の出力を使う)。seed Hessian は `saddle_initial_hessian_only=true`(振動ブロックを出さない)を既定にする。整合性ゲートを必須にする: `n_real(GV) = 3N−6(線形は −5)−n_imag`、|ZPE_GV/scale − ZPE_NWChem(最終)| <2e-5 Ha、|E_GV − E_hfauto| <1e-6 Ha。反応レベルでは |ΔE0‡−ΔE‡| >5 kcal/mol を異常とする。`hono_isomerization_v3` の thermo と rank を再実行して、公開済みの 24.11 kcal/mol を撤回する。
2. **フォールバックの fail-closed 化**: GoodVibes が成功しなかった場合、G は None(`thermo_unavailable`)にし、`scientific_rank_eligible` に `production_thermo_ready` を必須にする。内部 proxy(4 種)は削除する。NWChem の熱化学を使う場合は、S_rot >0 であることと線形の扱いを検査する。
3. **GoodVibes の設定**: 全プラットフォームで `goodvibes==4.3.0` に固定する(Windows では熱化学を実行しない方針でもよい)。extra_args は廃止し、型付き設定(qs, cutoff, vscal, zpe_vscal, T, conc/media, symm, invert)から CLI を生成する。ZPE は `--zpe-vscal`(PBE0 系で約 0.975)、振動は `-v`(約 0.989)に分けて根拠を記録する。`--symm` を渡す。溶媒モデルを使うパイプラインでは `-c 1.0` か `--media H2O` を付ける。温度ごとに `--temp {T}`(多数なら `--ti`)を付け、キャッシュキーを (artifact_id, T, qs, cutoff, scale, conc, symm) にする。|ν_imag| ≤30 cm⁻¹ が 1 本以下の極小は `--invert 30` と qRRHO で受理し、`soft_imaginary_inverted` と ±0.5 kcal/mol の不確かさを付ける。古い CSV は mtime と command の対応で除外する。
4. **理論レベルの一致**: `energy_method_id` と `frequency_method_id` を必須にする。composite は G = E(SP) + [G_GV − E_GV](freq)(または `--spc`)とする。反応内の全化学種で method signature が一致しなければ `mixed_level_of_theory` の failure にする。
5. **配座アンサンブル**: 同一組成・同一 basin 系について G_ens = −RT ln Σ g_i exp(−G_i/RT) を会合 ΔG と ΔG‡ の基準にする。
6. **速度定数**(必要になった場合だけ): thermo の後処理として 1 関数だけ置く。k(T) = κ(T)·(k_BT/h)·(RT/p°)^{Δn‡}·exp(−ΔG‡°/RT)、κ は非対称 Eckart。T <1.2·T_c なら `deep_tunneling_needs_SCT` を記録する。
7. **校正**: 同じ method/thermo 設定で、アミンの PA/GB(NIST)、F⁻+HF→HF2⁻、(HF)2 と NH3·HF の De/D0、HNC−HCN ΔE0 の MAE/MSE を算出し、reaction-rank の uncertainty_policy に入力する。method panel には DLPNO-CCSD(T) の SP を最低 1 つ入れることを必須にする。
8. **感度解析**: PYR·HF、TMA·(HF)2、NH3·(HF)3 のイオン対で、`--qs grimme/truhlar × cutoff 50/100/150 × --bav global/conf × fine/xfine の Hessian` を比較する。ΔG_assoc の幅 <0.5 kcal/mol なら「fine + qRRHO」を既定として正式に採用する。

**この系に対する有効性判定**

| 計算 | 判定 | 理由 |
|---|---|---|
| GoodVibes 4.3(-q、100 cm⁻¹) | **有効**(整合ゲート付きで) | HCN の S=48.15(実測 48.2)、PYR·HF の qRRHO 補正 0.55(独立計算 0.54)を正しく再現 |
| 内部フォールバック熱化学 | 低価値(有害) | 準線形分子で S_rot が負、qRRHO を過小評価 |
| basin-populations(Curtin–Hammett の補正) | **有効**(条件付き) | 符号は正しく、非網羅であることも明示している |
| thermo-sensitivity(HCN) | 限定的 | 自明な系で行っている。対象錯体で実施すべき |
| kinetics(Wigner/Eyring、Cantera、Arkane) | 低価値 | 物理的に不正確で、対象化学には TST が適用できない |
| calibrate(DB 支持スコア) | 低価値 | 手法誤差を評価していない |

---

## 4. 推奨ロードマップ

工数の目安: S = 1 日以内、M = 2〜5 日、L = 1〜3 週。

### フェーズ A: 即効性の高い正しさ・収束の改善(1〜2 週)

| # | 項目 | 対応 ID | 期待効果 | 工数 |
|---|---|---|---|---|
| A1 | HONO の ΔG‡ 24.11 を撤回し、GoodVibes の最終ブロック抽出と整合ゲート、フォールバックの fail-closed 化 | CH-38, CH-39 | 公開値の誤り(12 kcal/mol、5 kcal/mol)を是正し、再発を防ぐ | S〜M |
| A2 | pysisyphus+QCEngine 経路を IRC/TS から外す(NWChem の mepgs/QRC、string→saddle)。使い続ける場合は `fix_com/fix_orientation`、grid/SCF の伝達、フレームガード | CH-01, CH-02 | IRC の偽陽性を根絶し、RS-I-RFO が機能するようになる | M |
| A3 | IRC を分岐の下降で合否判定し、端点を登録極小へ一意に割り付け、二面角は状態窓で判定する | CH-03, CH-04 | HCN/HONO/NH3 の偽陽性と取り違えを解消 | S〜M |
| A4 | `basin_assessment` の回帰を修正し、reaction-plan→classify の統合テストを追加 | CH-05 | production の分類が機能するようになる | S |
| A5 | ReaDuct: normal_modes による虚振動の計数、`automatic_mode_selection`、IRC の出発構造一致、chdir | CH-26, CH-27, CH-31 | 低レベル TS/IRC の証拠が得られるようになる(pilot で 8 TS が IRC に進む見込み) | S |
| A6 | 前段ゲート: 結合判定のヒステリシス、|q| ≥0.3 Å、min_change ≥0.2 Å、単調経路は端点安定性チェックへ、線形補間ブラケット禁止、recover-path の fail-closed、xTB 陰性結果のゲート、エネルギー窓、preopt 棄却種の DFT 除外 | CH-06, CH-07, CH-28 | m3 級(約 30 h)の空振りをスキャン数十分で打ち切れる。PYR 級の 7 h も防げる | M |
| A7 | ZTS に gmax 判定、NEB の打ち切り条件、string を先頭エンジンにし、impose と IDPP の初期経路を使う | CH-08, CH-09 | 偽収束と数時間の停滞を解消 | S〜M |
| A8 | 既定を fine/SCF 1e-7 にし、lineage を 2 層化、saddle の Hessian はコピー(再計算しない) | CH-11, CH-12 | Hessian 1 回 3.6 倍と TS 段で約 1.5 h/反応を削減 | S〜M |
| A9 | 実行資源: `.wslconfig`、ext4 の scratch、memory 1,000〜1,500 MB/rank、グローバル rank 予算、OMP の明示 | CH-13 | SCF の 17 倍遅延と 4 h タイムアウトを解消し、20 コアを活用 | S |
| A10 | preopt を `--ohess vtight` にして回復ループを付け、xtb の非収束を検出、CREST に `--notopo` / `-T` / `--ewin` | CH-18, CH-20 | 鞍点 seed の DFT 流入(1 構造 1.5〜4.5 h)を防ぐ | S〜M |
| A11 | minimum-mode-follow/assess を標準フローに配線し、虚振動の 3 段方針を導入 | CH-17 | 軟モード rescue の再計算(約 4 h/件)を削減 | S |
| A12 | opt と freq を分割し、NWChem restart(`final-NNN.xyz`、movecs、`.drv.hess`)に対応 | AR-08 | タイムアウトによる全損(失敗時間の 31%)を削減 | M |

### フェーズ B: 構造改善(1〜2 か月)

| # | 項目 | 対応 ID | 期待効果 | 工数 |
|---|---|---|---|---|
| B1 | 型付き `CalculationEvidence` と 4 つの判定関数、qc を diagnostics と分離 | AR-23, AR-24, AR-27 | ゲート修正が 1 箇所で済む。テストのフラグ辞書が不要になる | L |
| B2 | 能力別 Protocol、型付きの遅延 registry、エンジンの `strategies` 宣言、`required_executables` | AR-11, AR-06, AR-13 | バックエンドを追加するときに触るのが registry だけになる | M |
| B3 | TS 層の 3 分割(PathGenerator / SaddleRefiner / ts_validation)と QM の注入 | AR-12, AR-15 | NEB 1 回・saddle 1 回を単体で検証・差し替えできる | L |
| B4 | thermo の 4 分割と legacy の削除、NWChem backfill を `migrate-manifest` へ | AR-20, CH-40 | CC 186 → 15 以下。1 反応ゴールデンテストが可能になる | M |
| B5 | `hfauto/reaction_case/` と `ReactionCaseDriver`(next_action ディスパッチャ)、scripts の取り込み | AR-01, AR-05 | 夜間の無人再計画と監査可能な試行履歴 | L |
| B6 | `execute_stage` を単一 API にし、`StageSpec`(consumes/produces/config_model/resources)と RunLayout・run_state を導入 | AR-03, AR-07, AR-25 | 部分再実行が安全になり、契約を機械的に検査できる | M |
| B7 | 共通 JobStore、FailureKind、attempt ladder、run_command の監視フック | AR-08, AR-26 | 再開と再試行を統一し、停滞時に早期打ち切り | M |
| B8 | 設定を 3 層化(pipeline/method/site)し、`resolved_config.yaml` と code_version を記録 | AR-10 | WSL/HPC 間の移植性と過去 run の再現性 | M |
| B9 | WorkItem/Executor(Local ProcessPool、Slurm array + dependency) | AR-09, AR-03 | 反応・化学種単位の並列化(dft-minima 15.7 h → 約 4 h) | M |
| B10 | テスト: 実出力のゴールデン、熱化学の教科書値、factories、rdkit の importorskip、設定値固定の縮約 | AR-30 | 収束対策での調整を妨げない、振る舞いのテスト | M |
| B11 | CI: import-linter(契約を見直したうえで)、pyrefly(core/chemistry から段階導入)、ruff C901 を stages/thermo から、vulture | 全体 | 構造の後戻りを防ぐ | S |
| B12 | 化学モデルの方針決定(気相か凝縮相か)と高レベル参照(DLPNO-CCSD(T))、solvation を第一級の設定に | CH-15, CH-16, CH-44 | 対象の化学(エッチング)に合ったモデルで、誤差を定量化できる | M(研究判断を含む) |

### フェーズ C: 削減と簡素化(フェーズ A と並行可能。4.5 節の順)

| # | 項目 | 期待効果 | 工数 |
|---|---|---|---|
| C1 | 死蔵モジュールと関数の削除(4.5 節の「低リスク」行) | 約 2,200 行を削減。radon と pyrefly の指標から雑音を除く | S |
| C2 | 未使用 stage(kinetics/calibrate/connector-audit/rank/descriptors)と関連バックエンドの削除 | 約 1,700 行を削減 | S〜M |
| C3 | hfauto.hpc と hfauto_ops の一本化、未参照 configs 15 本の削除 | 約 1,000〜1,600 行を削減。資源表の矛盾を解消 | M |
| C4 | legacy HF 熱化学分岐、互換エイリアス、`real_orca_executed` の整理 | CC の大幅削減 | M |
| C5 | ORCA バックエンドを削除するか SP 専用に作り直すか決める | 約 1,300 行の削減、または参照計算能力の獲得 | S(決定)/ M(再構築) |
| C6 | 孤立した `__pycache__`、リポジトリ直下の ReaDuct 出力、古い docs の掃除 | 誤読の防止 | S |

---

## 4.5 削除・簡素化候補一覧

凡例: **種別** = 死蔵(参照 0)/ 未使用 stage / legacy・互換 / 重複 / 旧フェーズの名残 / 古い docs・スクリプト。**リスク** = 低(参照 0 で機能への影響なし)/ 中(置き換えや手当てが必要)/ 高。根拠の「importer 0」は、hfauto・hfauto_ops・hfauto_viz・tests・configs・pyproject への grep と grimp による到達可能性解析(CLI の 3 エントリポイント + 41 stage を根とする)で確認したもの。

### 4.5.1 死蔵モジュール(低リスク)

| 対象パス(関数名) | 行 | 種別 | 根拠 | リスク | 影響テスト |
|---|---:|---|---|---|---|
| `hfauto/core/thermochemistry.py` | 158 | 死蔵・重複 | importer 0。Wigner/Eyring/K の重複実装 | 低 | なし |
| `hfauto/core/thermo.py` | 273 | 死蔵・重複 | importer は死蔵の simple.py 2 本だけ。引数順を推測する関数(BUG-06)を含む | 低 | なし |
| `hfauto/backends/thermo/simple.py`、`backends/kinetics/simple.py` | 61 + 45 | 死蔵 | registry に未登録、importer 0 | 低 | なし |
| `hfauto/backends/db/fixtures.py`、`backends/db/reference_data.py` | 22 + 89 | 死蔵・重複 | importer 0(core/reference_data.py と重複) | 低 | なし |
| `hfauto/chemistry/path_validation.py` | 127 | 死蔵 | importer 0。RMSD ≤1.25 Å という矛盾した規則。汎用の経路検証関数の移設先として置き換える | 低 | なし |
| `hfauto/reporting/molecule_dossier.py` | 8 | 死蔵 | NotImplementedError だけのスケルトン(同名の実装が hfauto_viz にある) | 低 | なし |
| `hfauto_viz/renderers/reaction_path.py`、`structure/annotations.py`、`structure/reaction_paths.py`、`core/run_loader.py` | 17 + 15 + 31 + 11 | 死蔵 | importer 0。run_loader は viz-layers と acyclic 違反の原因 | 低 | なし |
| `hfauto_ops/cache/reuse.py`、`cache/registry.py`、`compare/runs.py`、`reuse/plan.py`、`core/artifact_plan.py`、`schedulers/factory.py` | 93 + 175 + 79 + 99 + 196 + 73 | 死蔵・重複 | 本番からの importer 0(registry は死蔵の reuse からだけ参照される) | 低 | なし |
| `hfauto_ops/retry/policy.py`(`build_retry_plan`、`DEFAULT_RETRY_POLICY`)、`hfauto_ops/scheduler/resources.py` | 88 + 207 | 死蔵(テストのみ)・重複 | 本番からの importer 0。ポリシー 10 キー中 7 キーは一度も emit されない。STAGE_CONTRACTS の唯一の利用者 | 低〜中 | `tests/test_ops_current_stages.py`(43 行)を削除または書き換え |

### 4.5.2 死蔵関数・定数・フラグ(低リスク)

| 対象パス(関数名) | 種別 | 根拠 | リスク | 影響テスト |
|---|---|---|---|---|
| `hfauto/chemistry/nwchem_evidence.py`(`assess_endpoint_nwchem_evidence` CC 66・約 161 行、`REQUIRED_ENDPOINT_NWCHEM_METHOD`、`nwchem_scientific_input_signature`) | 死蔵 | 定義以外の参照 0。有効化すると B3LYP 系を誤って不合格にする危険な残骸 | 低 | なし(`parse_nwchem_input_evidence` は現役なので残し、B2 で backends へ移す) |
| `hfauto/core/qc.py`(`quality_value`、`summarize_failures`、`process_penalty_from_public_data`) | 死蔵・重複 | 参照 0(rank は public_data 版を使用)。**`TIER_CONFIDENCE_CAP`、`minimum_qc_from_frequencies` は現役なので対象外** | 低 | なし |
| `hfauto/core/public_data.py`(`STATIC_*` 表、`first_static_match`、`identity_confidence_score`) | 死蔵 | 参照 0。**`gas_process_flags` は現役**(rank を削除するなら一緒に削除) | 低 | なし |
| `hfauto/chemistry/rdkit_utils.py`(`write_mol_conformer_xyz`、`has_3d_coordinates`、`smiles_to_3d_mol`、`molecule_identity`) | 死蔵 | 定義元以外の参照 0 | 低 | なし |
| `hfauto/backends/conformer/crest.py:99-109`(`parse_crest_energies` のファイル直接指定の互換分岐) | legacy・互換 | 内部からはディレクトリしか渡されない | 低 | なし |
| `hfauto/backends/db/common.py`(`CachedHTTPMixin`) | 死蔵・重複 | 定義のみ。`http.CachedHTTPClient` と二重 | 低 | なし |
| `hfauto/core/schemas/`(`CalculationRecord`、`SpeciesRecord`、`SiteRecord`、`DescriptorRecord`、`ReactionRecord`) | 死蔵 | スキーマ外からの利用 0(`__init__` の再エクスポートと viz の docstring だけ)。B1 の Evidence で置き換える | 低 | なし |
| `hfauto/chemistry/populations.py:10`(`R_KCAL_MOL_K`) | 重複 | core/constants.py と二重定義 | 低 | reaction_rank の import を付け替え |
| `hfauto/backends/thermo/goodvibes.py`(`species_correction`)、`backends/thermo/arkane.py`(`write_*_stub`)、`registry.get_arkane_engine` | 死蔵 | 呼び出し 0 | 低 | なし |
| `hfauto/stages/irc.py:123`(`allow_legacy_rmsd_only_irc_validation`) | 死蔵フラグ | 設定されるだけで読まれない | 低 | なし |
| `hfauto/backends/ts/nwchem_neb.py:805-829`(`run_irc` のスタブ) | 死蔵 | 常に失敗を返す。IRC registry から外す | 低 | なし |
| `hfauto/workflow/reaction_state.py`(`_STRATEGY_FOR_DIAGNOSIS['final_frequency_missing']`、`endpoint_local_saddle_search` のエンジン表) | 死蔵 | 到達不能(AR-06) | 低 | `tests/test_reaction_case_workflow.py` を確認 |

### 4.5.3 未使用 stage と関連バックエンド(現行 13 パイプラインで未使用、化学的価値も低い)

| 対象パス | 行 | 種別 | 根拠 | リスク | 影響テスト |
|---|---:|---|---|---|---|
| `hfauto/stages/kinetics.py`、`backends/kinetics/{cantera,tst}.py`、`backends/thermo/arkane.py`、`configs/stages/kinetics_*.yaml`、registry の該当エイリアス | 124 + 266 + 54 + 142 | 未使用 stage | パイプライン使用 0、物理的に不正確(CH-43)、import-linter 違反 3 件。tests はソース削除済み(pyc のみ) | 中(速度定数が必要なら thermo 後処理として 1 関数で再実装) | なし |
| `hfauto/stages/calibrate.py` | 274 | 未使用 stage | パイプライン使用 0、DB 支持スコアの集計にすぎない(CH-44)。誤差評価で置き換える | 中 | なし |
| `hfauto/stages/connector_audit.py` | 157 | 未使用 stage | パイプライン使用 0、config 引数も未使用 | 中(DB provider の扱いと合わせて判断) | なし |
| `hfauto/stages/rank.py` の legacy HF スコアリング(l.83-316)と `rank` stage | 316 | legacy・未使用 | パイプラインはすべて reaction-rank を使う。rank→reaction_rank の stage 間違反 | 中(`"rank"` を reaction-rank の別名にするか、登録から外す) | なし(`test_reaction_rank_uncertainty.py` は ReactionRankStage を使うので影響なし) |
| `hfauto/stages/descriptors.py`、`chemistry/descriptors.py` のプレースホルダ(`HF_fragment_charge`、`dipole_D`) | 60 + 一部 | 未使用 stage | パイプライン使用 0、`hf_stretch_cm1=max(ν)` の定義が誤り | 低 | なし |
| `hfauto/stages/hpc_plan.py`、`hfauto/hpc/*`(cache/compare/dashboard/job_plan/resources/schedulers)、hfauto CLI の hpc 系コマンド | 91 + 628 | 未使用・重複 | hfauto_ops と三重実装、Snakefile に依存辺なし、実績 0(AR-03/AR-04) | 中(ops へ一本化してから削除) | CLI の import(`hfauto/cli/main.py:12-17`) |
| `hfauto/stages/ops.py`、`stages/viz.py` | 39 + 91 | 未使用・境界違反 | hfauto-without-addons 違反。add-on 側へ移し、entry point で登録する | 中 | なし |

### 4.5.4 legacy・互換コード(中リスク)

| 対象パス(関数名) | 種別 | 根拠 | リスク | 影響テスト |
|---|---|---|---|---|
| `hfauto/stages/thermo.py:960-1157`(legacy HF の RC/IP/candidate/HF-cluster 分岐) | legacy | それを生成する stage(build_hf)は削除済み。stoichiometry がないと入り込み、誤ラベルの記録を作る | 中(stoichiometry がなければ `missing_stoichiometry` で failure にする) | なし(実行テスト 0) |
| `hfauto/core/schemas/thermo.py:28-38`、`species.py:18-20`、`reaction.py:21-26`(hf_n、candidate_species_id、hf_cluster_species_id、G_*_hartree、delta_G_assoc/ionpair/act) | legacy・互換 | 生成箇所 0(thermo/rank/schemas の外) | 中 | なし |
| `hfauto/stages/thermo.py:8,295-310`(`backfill_nwchem_thermochemistry_frequencies` の呼び出し) | legacy・移行 | 旧スキーマ移行を stage が実行している | 中(`hfauto migrate-manifest` へ移設) | なし |
| `hfauto/chemistry/reaction_path_qc.py:600-667`(legacy HF 分岐) | legacy | `ion_pair/neutral_complex` がハードコード | 中 | `tests/test_reaction_path_qc.py` を確認 |
| `hfauto/chemistry/geometry_qc.py:42-51`(`HF_` 成分名の分岐)、`backends/qm/xtb.py:79-81`(`setdefault False`) | legacy | ComponentRecord に `name` がない。評価していないのに「問題なし」と記録する | 低 | なし |
| `thermo_models.low_frequency_correction`、`scale_thermal_correction`、`process_adjusted_reaction_delta_g`、GoodVibes の内部フォールバック経路 | legacy・重複 | qRRHO proxy 4 種の 1 つ。フォールバックの fail-closed 化(CH-39)と同時に削除する | 中 | なし |
| `real_orca_executed` のフォールバック読み取り(thermo、dft_minima、ops、viz。計 27 参照) | legacy・互換 | ORCA を残すかどうかの判断と連動 | 中 | なし |
| `hfauto/backends/ts/pysisyphus.py`(`search_ts` の NWChem NEB 委譲、`_run_endpoint_dimer` 316 行) | 重複・低価値 | 名前と実際の処理が一致しない。dimer は極小付近で無効(CH-10) | 中 | `tests/test_nwchem_backend.py` の pysisyphus 関連を確認 |
| `hfauto/backends/qm/orca.py`(503)、`backends/ts/orca_nebts.py`(826)、`configs/stages/{sp_orca,dft_minima_orca}.yaml` | 要判断 | 実行実績 0、テスト 0、ORCA 4 の文法、Gibbs を E で代用、ダミー代替。**化学的には参照計算(DLPNO)の能力として有用なので、「削除」か「SP 専用に再構築」かを決める** | 中 | なし |

### 4.5.5 重複実装(統合)

| 対象 | 種別 | 根拠 | リスク | 影響テスト |
|---|---|---|---|---|
| `latest_manifest_path`(`reporting/html_report.py:12`、`hfauto_viz/core/paths.py:14`) | 重複 | 挙動が違う(CWD 依存)。15 モジュールが import している | 中(`core/run_layout.py` に統合) | ops/viz の関連テスト |
| `canonical_species_id` のインライン再実装(21〜25 ファイル、27〜36 箇所)、`hfauto_ops/core/array_jobs.py:52` の同名複製 | 重複 | `core.artifacts.canonical_species_id` がある | 低 | なし |
| `_reaction_id`(TS の 8 ファイル)、method マージの `_method`(6 箇所) | 重複 | 同一のヘルパー | 低 | なし |
| Wigner(3)/Eyring(5)/qRRHO proxy(4)/K(4) | 重複 | `core/thermo_models.py` 1 つにし、引数はキーワード専用にする | 低 | なし |
| `build_reuse_plan`(3)、compare_runs(3)、retry 立案(3)、資源表(3)、artifact→stage 対応表(4) | 重複 | AR-04 | 中 | `tests/test_ops_current_stages.py` |
| 停滞判定(`assess_optimizer_stagnation` と path_diagnostics の window=3 実装) | 重複 | テスト済みの方が本番で使われていない | 中(片方に統一) | `tests/test_reaction_case_workflow.py` |
| checkpoint/再利用判定(preopt、dft_minima、ts_search、build_complexes、explore_reactions の 5 実装) | 重複 | AR-08 | 中(JobStore へ) | 各 stage のテスト |

### 4.5.6 旧フェーズの名残・古い docs・スクリプト・生成物

| 対象 | 種別 | 根拠 | リスク | 影響テスト |
|---|---|---|---|---|
| 孤立した `__pycache__`: hfauto/hfauto_ops 側 27 モジュール分(`build_hf`、`proton_transfer_endpoints`、`proton_transfer_scan`、`path_plan`、`wet_etch*`、`hf_builder`、`reaction_generation`、`xtb_scan`、`scan_optts`、`ops_report` ほか)、tests 側 36 モジュール分(`test_phase*` ほか)、scripts 側 12 モジュール分(計 75) | 旧フェーズの名残 | 対応する .py が作業ツリーにない。レビュー時に「部位検出は使われている」「テストがある」と誤読させる | 低 | なし(`find . -name __pycache__ -prune -exec rm -rf {} +`) |
| `configs/stages/*.yaml`(10)、`configs/hpc/*.yaml`(3)、`configs/ops/default.yaml`、`configs/viz/default.yaml` | 旧設定 | コード・docs・テスト・他の設定から参照 0。`dft_minima_orca.yaml` は `fallback_to_dummy: true`、`preopt_xtb.yaml` は `allow_subprocess: false` | 低 | なし |
| `docs/architecture.md`、`docs/developer_guide.md`、`docs/stage_contracts.md`(各 3 行) | 古い docs | 内容が実態と矛盾する(「protocol を実装して registry に登録するだけ」など)。docs/current に統合 | 低 | なし |
| `docs/hf1_scan_pilot.md` | 古い docs | 参照している `scripts/run_hf1_scan_pilot.py` は削除済み。**スキャン機能自体は化学的に有用(前段ゲート)なので、stage として再実装し、そこに記述を移す** | 低 | なし |
| リポジトリ直下の `afir_biased/`、`nt2_guess/`、`product_unbiased/`、`ts_bofill/`、`ts_hessian/`、`internal_coords.log`、`runs/legacy_unscoped_raw_20260814/` | 生成物の漏れ | ReaDuct/pysisyphus が CWD に書いた出力。別の run の出力が混在している | 低(chdir 修正の後に削除し、`.gitignore` に追加) | なし |
| `scripts/run_nwchem_saddle.py`、`run_nwchem_string.py`、`run_nwchem_frequency.py`、`run_pysisyphus_irc.py`、`run_pysisyphus_saddle.py`、`merge_manifests.py`(git 未追跡) | ゲート迂回スクリプト | stage の lineage とゲートを通らない(`hcn_irc_v5` の偽陽性の原因、固定値の証拠フラグ) | 中(`fixed-geometry-frequency` stage と ReactionCaseDriver の WorkItem に取り込んでから削除) | なし |

### 4.5.7 削除せず「配線・修正」すべきもの(化学的に有用)

| 対象 | 理由 | 必要な作業 |
|---|---|---|
| `minimum-mode-follow`、`minimum-mode-assess` | 同一 basin を確定させる最も決定的な手法 | 標準フローに配線(CH-17) |
| `path-ensemble`、`path-ensemble-assess` | 多解像度で障壁なしを判定する | gmax 判定を修正し(CH-08)、engine=='nwchem' 分岐を除いてから配線 |
| `path-intermediates`、`reaction-segments` | 振動数で検証した第 3 の極小から子反応を作る | reaction-segments に stoichiometry を付与 |
| `endpoint-seeds`、`endpoint-seed-screen` | 端点の対称性破りと多様化 | 入力の生成元(`endpoint_pair_selection`)を実装するか、入力を変更。虚振動で除外 |
| `detect-sites`、`chemistry/site_detection.py` | 塩基性部位の選択(4-アミノピリジンなど N が複数ある分子) | build-complexes の anchor に配線する(配線しないなら削除) |
| `enrich` | 公開 DB 由来の同定情報 | DB provider の位置づけ(BUG-01 の修正を含む)を決める |
| xTB 緩和スキャン(旧 hf1/hf2 scan) | 0.1〜0.2 h で DFT 約 49 h と同じ結論を出した | 前段ゲート stage として再実装する |

**削減見込み**: 4.5.1〜4.5.3 と 4.5.4 の legacy 分岐で約 4,500 行(ORCA を削除する場合は +約 1,300 行)。

---

## 5. 付録

### 5.1 全指摘一覧

検証欄: 確認 = confirmed、一部訂正 = partially(本文は訂正後の主張)。refuted は 0 件。

**アーキテクチャ**

| ID | 重大度 | 検証 | 場所 | 一行要約 | 元 ID |
|---|---|---|---|---|---|
| AR-01 | 高 | 確認/一部訂正 | `hfauto/workflow/runner.py:7-9`、`stages/reaction_plan.py:23` | workflow がランナーとドメインサービスを兼ね、stages↔workflow が循環 | A1-1, A2-9, A3-8 |
| AR-02 | 中 | 確認 | `hfauto/stages/ops.py:8`、`hfauto_ops/core/array_jobs.py:31` | add-on 境界が双方向で、配布物をまたぐ循環 | A1-7 |
| AR-03 | 高 | 確認/一部訂正 | `hfauto_ops/core/array_jobs.py:256-278`、`hfauto/hpc/job_plan.py:45-51` | 実行経路が 4 系統あり、配列/HPC 経路は存在しない CLI を呼ぶ | A1-2, A3-2 |
| AR-04 | 中 | 確認 | `hfauto/hpc/resources.py:31-54`、`hfauto_ops/scheduler/resources.py` | hpc と ops の三重実装、資源表の矛盾、死蔵 ops | A1-9, A3-3 |
| AR-05 | 高 | 一部訂正 | `hfauto/workflow/reaction_state.py:52-257`、`scripts/run_*.py` | TS 以外の next_action にディスパッチャがなく、scripts による手動連鎖とゲート迂回 | A3-1 |
| AR-06 | 中 | 一部訂正 | `hfauto/workflow/reaction_state.py:7-28` | 状態機械の分散、エンジン表の二重管理、語彙の不一致、到達不能な戦略 | A2-2, A3 見落とし |
| AR-07 | 高 | 確認 | `hfauto/workflow/runner.py:47-61,92-93`、`stages/ts_search.py:69-89` | 部分再実行で古い下流が混入、resume がメソッド指紋を見ない、CWD 依存パス | A3-4, A1-8 |
| AR-08 | 高 | 確認/一部訂正 | `hfauto/backends/qm/nwchem.py:782,1133`、`stages/dft_minima.py:783-797` | opt+freq の単一ジョブ・restart なし・初期構造から再実行。失敗時間 31% | A6-4, A3-5, A3-6, A2-5, A4-7 |
| AR-09 | 中 | 一部訂正 | `hfauto/stages/ts_search.py:185-495` | 並列化の手段がない | A3-7 |
| AR-10 | 中 | 一部訂正 | `configs/pipelines/*.yaml`、`core/executables.py` | 絶対パス 42 箇所、メソッドのコピペ、解決済み設定の未保存、env で上書き不可 | A3-10 |
| AR-11 | 高 | 確認 | `hfauto/backends/base.py:10`、`registry.py`、`core/environment.py:63` | Protocol 未使用、registry に型がない、エンジン名のハードコード、能力が不揃い | A1-3, A6-1, A2-8, A4-10 |
| AR-12 | 高 | 確認/一部訂正 | `hfauto/backends/ts/nwchem_neb.py:214,694`、`pysisyphus.py:366-373` | TS バックエンドの丸抱え、QM の直接生成、NWChem への固定 | A6-2, A2-1, A4-4 |
| AR-13 | 中 | 確認 | `hfauto/stages/build_complexes.py:671`、`kinetics.py:7-9` | stage が具象バックエンドを直接 import(7 件)、registry の能力が不一致 | A1-6 |
| AR-14 | 高 | 確認 | `hfauto/chemistry/nwchem_evidence.py`、`stages/dft_minima.py:421` | chemistry/stages に NWChem 固有のパーサ・既定値・ゲート | A1-4, A2-6, A6-8 |
| AR-15 | 中 | 確認 | `hfauto/backends/ts/orca_nebts.py:545-616` | TS ゲートがエンジンごとに違い、ORCA は production で検証不能 | A1-5, A2 見落とし |
| AR-16 | 中 | 一部訂正 | `hfauto/backends/qm/xtb.py:85,103,118`、`stages/sp.py:58` | ダミー実装の本番混入(xTB single_point が捏造値を success で返す) | A6-3, A2-3 |
| AR-17 | 中 | 一部訂正 | `hfauto/workflow/ts_execution.py:96`、`chemistry/methods.py:9` | method dict が約 136 キーの型なし契約。資源が calc_id に混入 | A6-5 |
| AR-18 | 中 | 一部訂正 | `hfauto/backends/qm/nwchem.py`、`xtb.py:184` | QM バックエンドの責務混在。xtb が run_command を迂回 | A4-5 |
| AR-19 | 中 | 一部訂正 | `hfauto/backends/qm/orca.py:452`、`orca_nebts.py:662`、`conformer/crest.py:133` | 外部 I/O の堅牢性が不均一(Gibbs←E、mtime による端点選択、CREST の None→0.0) | A6-7 |
| AR-20 | 高 | 一部訂正/確認 | `hfauto/stages/thermo.py:339,960-1157` | ThermoStage.run CC 186。legacy 分岐は誤ラベルを生む潜在バグ。実行テスト 0 | A4-1, C6-9, A2-7 |
| AR-21 | 高 | 確認 | `hfauto/stages/preopt.py:83`、`dft_minima.py:148` | プロトン状態・極小受理・lineage のゲートが複数実装で基準が違う | A4-6 |
| AR-22 | 中 | 確認 | `hfauto/chemistry/reaction_classification.py:249` | chemistry でカーネルと証拠・I/O が混在。classify_reaction は 14 引数・CC 115 | A4-9 |
| AR-23 | 高 | 確認 | `hfauto/core/schemas/artifact.py:19-28` | qc/data が型なし辞書、キーの 4 割が未読、スキーマ 5 クラスが未使用 | A5-1, A4-2 |
| AR-24 | 中 | 一部訂正 | `hfauto/workflow/ts_execution.py:125-144` ほか 26 ファイル | ゲート判定式の複製と欠損キーの扱いの揺れ | A5-2, A2-3 |
| AR-25 | 中 | 確認 | `hfauto/core/artifact_types.py:51-225` | STAGE_CONTRACTS は強制されず実装とずれ、メタデータが 6 箇所に分散 | A1-11, A2-4, A3-9, A4-8, A5-7 |
| AR-26 | 中 | 確認 | `hfauto_ops/core/retry.py:30-52` | 失敗分類が自由文字列で、retry が部分一致で誤った stage に振り分ける | A5-5, A6-6 |
| AR-27 | 中 | 確認 | `hfauto/stages/thermo.py:713-731,926` | 重複フラグと擬似的な confidence_score(ops に表示される) | A5-4 |
| AR-28 | 低 | 確認 | `hfauto/stages/*`(定数フラグ 63 種) | 免責フラグ、未使用の next_action、副次出力の肥大 | A5-6 |
| AR-29 | 低 | 一部訂正 | `hfauto/chemistry/method_lineage.py` ほか 21 ファイル | 下流での再ハッシュ・再監査、method_evidence_* の 5 変種 | A5-10 |
| AR-30 | 中 | 一部訂正 | `tests/*`、`hfauto/chemistry/reaction_path_qc.py:306-329` | 形状確認に偏ったテスト、fixture 0、ゴールデン 0、rdkit による収集中断 | A5-8, A5-9 |
| AR-31 | 中 | 確認 | 4.5.1〜4.5.2 | 到達不能 35 モジュール、重複した物理関数(引数順の推測) | A1-10, A4-3, A6-10, C6-10 |
| BUG-01 | 中 | 確認 | `hfauto/backends/db/comptox.py:51-53` | CompTox の応答のたびに AttributeError | pyrefly, A6-10 |
| BUG-02 | 中 | 確認 | `hfauto_ops/schedulers/slurm.py:19`、`lsf.py:19` | Slurm の submit が str を返し必ず失敗。LSF の `<` が効かない | pyrefly |
| BUG-03 | 低 | 確認 | `hfauto/stages/thermo.py:620` | thermo の arkane 指定で AttributeError | pyrefly, A2-7 |
| BUG-04 | 低 | 確認 | `hfauto/backends/ts/nwchem_path_support.py:85` | `None > float` で TypeError | pyrefly |
| BUG-05 | 低 | 確認 | `hfauto/backends/registry.py:39` | 名前のない spec が KeyError(None) | pyrefly |
| BUG-06 | 中 | 確認 | `hfauto/core/thermo.py:198,229` | 引数順を値の大きさで推測(死蔵) | A4-3, C6 見落とし |

**化学計算**

| ID | 重大度 | 検証 | 場所 | 一行要約 | 元 ID |
|---|---|---|---|---|---|
| CH-01 | 致命 | 確認 | `hfauto/backends/ts/pysisyphus.py:846-870`、`pysisyphus_saddle.py:139-175` | QCEngine 経由で勾配と Hessian が再配向フレームのまま返る(IRC/RS-I-RFO/dimer が無効) | C3-1, C4-1, C2-4 |
| CH-02 | 高 | 確認/一部訂正 | `hfauto/backends/ts/pysisyphus.py:111-205` | QCEngine 経路の PES 不一致(medium/1e-6/1 GB)、provenance の誤記 | C2-4, C3-4, C4-3 |
| CH-03 | 致命 | 確認 | `hfauto/backends/ts/pysisyphus.py:995-1175`、`stages/irc.py:230-262` | IRC の合否が分岐の下降を見ず、hcn_irc_v5 が偽陽性 | C4-2 |
| CH-04 | 高 | 確認 | `hfauto/chemistry/reaction_path_qc.py:547-554`、`reactions.py:392-407` | 端点割付けの閾値が緩く、HONO の cis/trans を取り違え。二面角座標の定義に矛盾 | C4-4 |
| CH-05 | 致命 | 確認 | `hfauto/stages/reaction_plan.py:160-182`、`chemistry/reaction_classification.py:282-291` | basin_assessment のキー欠落(回帰)で production の分類が全件 unresolved | C4-5 |
| CH-06 | 高 | 一部訂正/確認 | `hfauto/chemistry/path_diagnostics.py:247-257`、`stages/recover_path.py:236-242` | 障壁の前段ゲートがなく、単調経路を線形補間ブラケットへ回す。fail-open | C3-2, C4-7, C5-6, C2-9 |
| CH-07 | 高 | 確認/一部訂正 | `hfauto/chemistry/connectivity.py:31-46`、`minimum_connections.py:176-205` | 単一閾値の結合グラフと min_change 0.05 Å で同一 basin を反応扱い | C1-4, C4-6, C5-5 |
| CH-08 | 高 | 確認 | `hfauto/chemistry/path_diagnostics.py:150-200`、`backends/ts/nwchem_string.py:59-73` | ZTS が変位だけで収束判定し、gmax の悪化を見逃す | C3-3 |
| CH-09 | 高 | 確認/一部訂正 | `hfauto/backends/ts/nwchem_neb.py:134-205`、`reaction_state.py:23` | NEB に impose/IDPP/CI がなく先頭エンジン。実反応 9 件で収束 0 | C3-7, A6-9 |
| CH-10 | 高 | 確認 | `hfauto/backends/ts/pysisyphus_saddle.py:127-172` | RS-I-RFO が root 0・cart・未射影 Hessian。dimer は低価値 | C3-5 |
| CH-11 | 高 | 確認 | `hfauto/backends/ts/nwchem_saddle.py:134-257`、`qm/nwchem.py:797-813` | 同じ構造の Hessian を二重計算(900〜3,400 s × 5〜6 回) | C2-5, C3-6 |
| CH-12 | 高 | 確認 | `configs/pipelines/*.yaml`、`chemistry/method_lineage.py:12-13` | xfine/1e-8 の過剰精度と、経路段への強制 | C2-2 |
| CH-13 | 高 | 一部訂正 | `C:/Users/user/.wslconfig`、`hfauto/backends/qm/nwchem.py:884-910` | WSL 4 vCPU 上の mpirun 同時起動、/mnt/c の scratch、メモリの過剰宣言 | C2-1 |
| CH-14 | 高 | 一部訂正 | `hfauto/backends/qm/nwchem.py:752-756`、`stages/dft_minima.py:135-148` | noautoz 固定、対角初期 Hessian、timeout 後は seed から再実行 | C2-3, A6-9 |
| CH-15 | 高 | 確認/一部訂正 | `configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml`、`backends/reaction_discovery/readuct.py:142-159` | 気相 PBE0-D3(0)/def2-SVPD で PT 状態を判定し、溶媒の経路がない | C2-6, C1-11, C5-9 |
| CH-16 | 中 | 確認 | `hfauto/chemistry/methods.py:37-67`、`method_lineage.py:12` | solvation_model は記録されるだけで、COSMO が lineage の外 | C2-7 |
| CH-17 | 中 | 一部訂正/確認 | `hfauto/core/frequency_qc.py:13`、`stages/minimum_mode_follow.py` | 虚振動 −1 cm⁻¹ と自動変位の未配線。熱化学感度と釣り合わない | C2-8, C4-8, C6-6 |
| CH-18 | 高 | 一部訂正 | `hfauto/backends/conformer/crest.py:235-260` | CREST の既定トポロジーフィルタが PT 構造を除去し、設定も伝わらない | C1-1, C1-6, C1-7 |
| CH-19 | 高 | 一部訂正 | `hfauto/chemistry/complexes.py:60-179`、`stages/build_complexes.py:174-279` | 配置 seed が対称で重複する。実験室座標系のフォールバック、donor 選択の誤り | C1-2, C1 見落とし |
| CH-20 | 高 | 確認 | `hfauto/backends/qm/xtb.py:47,150,284` | preopt に Hessian がなく鞍点を昇格させる。xtb の非収束を検出できない | C1-3, C1 見落とし |
| CH-21 | 中 | 確認 | `hfauto/chemistry/geometry_qc.py:42-51` | HF 判定が死んでいて、評価していない項目を「問題なし」と記録 | C1-8 |
| CH-22 | 中 | 一部訂正 | `hfauto/stages/build_complexes.py:52-70,718-778` | アンサンブル選別が GFN2 最安に偏り、重複除去と DFT 再ランクがない | C1-5 |
| CH-23 | 低 | 一部訂正/確認 | `hfauto/chemistry/rdkit_utils.py:79-153`、`stages/detect_sites.py` | RDKit の複数断片が重なる。detect-sites が未配線 | C1-9, C1-10 |
| CH-24 | 中 | 一部訂正 | `hfauto/backends/qm/orca.py`、`configs/stages/sp_orca.yaml` | ORCA バックエンドが古く、未検証 | C2-10, C3-10 |
| CH-25 | 低 | 一部訂正 | `hfauto/backends/qm/nwchem.py:567`、`stages/method_panel.py:83` | xTB/CREST ラッパーの細部、hf_stretch の誤定義、panel の鍵の衝突 | C2-11, C2-12 |
| CH-26 | 致命 | 確認 | `hfauto/backends/reaction_discovery/readuct.py:367-385` | 未射影の固有値で虚振動を数え、有効な TS を全件棄却 | C5-1 |
| CH-27 | 高 | 確認 | `hfauto/backends/reaction_discovery/readuct.py:282-315,445-520` | TS 最適化のモード選択なし、IRC の出発構造一致と端点の極小性を未検証 | C5-2, C5 見落とし |
| CH-28 | 高 | 確認 | `hfauto/stages/explore_reactions.py:400-587`、`chemistry/reaction_trials.py:251-278` | エネルギー窓と妥当性のゲートがなく、TS なし回収・未収束 AFIR を昇格 | C5-3, C5 見落とし |
| CH-29 | 中 | 一部訂正 | `hfauto/chemistry/reaction_trials.py:174-249`、`workflow/discovery_coverage.py` | trial が role を無視、relay の優先度が低い、coverage が切り捨てを示さない | C5-4 |
| CH-30 | 中 | 一部訂正 | `hfauto/backends/reaction_discovery/readuct.py:246-281` | NT2 は線形分子とねじれに不適。B-spline を使っていない | C5-7 |
| CH-31 | 中 | 確認 | `hfauto/backends/reaction_discovery/readuct.py:185-575`、リポジトリ直下のディレクトリ | ReaDuct の出力が CWD に漏れ、一次証拠が失われる。attempt に timeout がない | C5-8 |
| CH-32 | 低 | 一部訂正 | `hfauto/backends/reaction_discovery/readuct.py:74-107` | AFIR の γ=300 と単一ペア | C5-10 |
| CH-33 | 中 | 一部訂正 | `hfauto/backends/ts/base.py:492-575`、`chemistry/saddle_seeds.py:53-103` | saddle の崩壊を分類しない、prominence の下限がない | C3-8 |
| CH-34 | 中 | 確認 | `hfauto/backends/ts/base.py:180-305`、`core/qc.py:154-190` | モードの重なり 0.5 を hard gate にしている | C3-9 |
| CH-35 | 中 | 確認 | `hfauto/chemistry/basin_identity.py:86-181,311-329` | 置換等価な縮退素反応を same_basin に併合 | C4-9 |
| CH-36 | 中 | 一部訂正 | `hfauto/chemistry/reaction_path_qc.py:600-667`、`stages/irc.py:123` | IRC 検証コードの重複・デッドフラグ・HF のハードコード | C4-10 |
| CH-37 | 低 | 確認 | `hfauto/stages/ts_search.py:195,303-311`、`backends/ts/pysisyphus.py:357-374` | attempt の timeout が stage の walltime を無視。search_ts の委譲 | C3-10 |
| CH-38 | 致命 | 確認 | `hfauto/backends/thermo/goodvibes.py:254,520-545`、`backends/qm/nwchem.py:797-866` | GoodVibes が振動ブロックを連結し、HONO の ΔG‡ 24.11(正しくは約 12.2)を production 公開 | C6-1 |
| CH-39 | 高 | 確認 | `hfauto/backends/thermo/goodvibes.py:139-177,461-533`、`stages/thermo.py:926-927` | 内部フォールバック熱化学の誤値。fallback フラグが死んでいる | C6-2 |
| CH-40 | 中 | 一部訂正 | `hfauto/stages/thermo.py:436-575` | エネルギー源と熱補正源で理論レベルが混在しうる(潜在) | C6-3 |
| CH-41 | 中 | 確認 | `hfauto/backends/qm/nwchem.py:753`、`thermo/goodvibes.py:161` | σ=1 固定、標準状態と配座アンサンブルに非対応 | C6-4 |
| CH-42 | 中 | 確認 | `hfauto/backends/thermo/goodvibes.py:365,489-501` | 温度依存を扱えず、温度判定が fail-open | C6-5 |
| CH-43 | 中 | 確認 | `hfauto/backends/kinetics/tst.py:19-28`、`cantera.py:92-151` | kinetics 段が物理的に不正確で未使用 | C6-7 |
| CH-44 | 中 | 確認 | `hfauto/stages/calibrate.py:24-165` | calibrate は手法誤差を評価していない | C6-8 |
| CH-45 | 低 | 確認 | `pyproject.toml:47`、`backends/thermo/goodvibes.py:318-334` | GoodVibes のバージョンが OS で分かれ、スケールを共用。古い CSV が混入しうる | C6-11, C6 見落とし |
| CH-46 | 低 | 確認 | `hfauto/stages/thermo_sensitivity.py:29-54` | 感度解析が剛直な HCN だけ。ソース選択規則が重複 | C6-12 |

### 5.2 参考文献・仕様 URL

**ツール**
- import-linter: https://import-linter.readthedocs.io/en/stable/contract_types.html
- pyrefly: https://pyrefly.org/en/docs/
- radon: https://radon.readthedocs.io/en/latest/intro.html
- Ruff(C901 など): https://docs.astral.sh/ruff/rules/complex-structure/
- vulture: https://github.com/jendrikseipp/vulture
- pydantic: https://docs.pydantic.dev/latest/concepts/models/

**ワークフロー・HPC**
- AiiDA WorkChain: https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/workflows/concepts.html
- jobflow dynamic flows: https://materialsproject.github.io/jobflow/tutorials/5-dynamic-flows.html
- Snakemake rules: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
- Slurm sbatch / job array: https://slurm.schedmd.com/sbatch.html 、https://slurm.schedmd.com/job_array.html
- Open MPI scheduling: https://docs.open-mpi.org/en/main/launching-apps/scheduling.html
- WSL filesystems: https://learn.microsoft.com/en-us/windows/wsl/filesystems

**配座・xTB・CREST・RDKit**
- CREST keywords: https://crest-lab.github.io/crest-docs/page/documentation/keywords.html 、例 3: https://crest-lab.github.io/crest-docs/page/examples/example_3.html
- CREST(JCP 2024): https://pubs.aip.org/aip/jcp/article/160/11/114110/3278084/CREST-A-program-for-the-exploration-of-low-energy
- xtb hessian / optimization: https://xtb-docs.readthedocs.io/en/latest/hessian.html 、https://xtb-docs.readthedocs.io/en/latest/optimization.html
- GFN2-xTB(JCTC 2019): https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176
- CENSO(JPCA 2021): https://pubs.acs.org/doi/10.1021/acs.jpca.1c00971
- RDKit rdDistGeom / rdDetermineBonds: https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html 、https://rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html
- (HF)2 分子線分光: https://pubs.aip.org/aip/jcp/article/81/12/5417/91420/The-molecular-beam-spectrum-and-the-structure-of
- xTB の溶媒・プロトン移動関連: https://pubs.acs.org/doi/10.1021/acs.jctc.1c00471 、https://pubs.acs.org/doi/abs/10.1021/acs.jctc.2c00239

**電子状態**
- NWChem DFT / Geometry Optimization / GEOMETRY / FAQ / Memory / Scratch_Dir / COSMO / Vibration / Hessians:
  https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html 、https://nwchemgit.github.io/Geometry-Optimization.html 、https://nwchemgit.github.io/Keywords-for-the-GEOMETRY-directive.html 、https://nwchemgit.github.io/FAQ.html 、https://nwchemgit.github.io/Memory.html 、https://nwchemgit.github.io/Scratch_Dir.html 、https://nwchemgit.github.io/COSMO-Solvation-Model.html 、https://nwchemgit.github.io/Vibration.html 、https://nwchemgit.github.io/Hessians-and-Vibrational-Frequencies.html
- ORCA 数値積分 / DFT / 3c / NEB: https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/numericalintegration.html 、https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html 、https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html 、https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/neb.html
- ORCA 4 の grid 文法廃止: https://github.com/qcxms/QCxMS/issues/23
- DFT の grid 感度(Bootsma & Wheeler): https://chemrxiv.org/engage/chemrxiv/article-details/60c7436f0f50dba3a3395ef2
- D3(BJ): https://onlinelibrary.wiley.com/doi/10.1002/jcc.21759
- 非局在化誤差と厳密交換(PCCP 2024): https://pubs.rsc.org/en/content/articlehtml/2024/cp/d4cp00907j
- 障壁と diffuse 関数(Zheng/Xu/Truhlar): https://link.springer.com/article/10.1007/s00214-010-0846-z
- GFN2 のアミンのプロトン移動誤差: https://pmc.ncbi.nlm.nih.gov/articles/PMC12288007/
- その他(C2-6 の引用): https://www.sciencedirect.com/science/article/abs/pii/0009261489870605 、https://wires.onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1631 、https://www.sciencedirect.com/science/article/abs/pii/000926149190073I

**TS 探索・IRC**
- NWChem NEB/ZTS マニュアル: https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html
- NWChem 7.2.3 ソース: https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F 、https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/neb/neb_utils.F 、https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/neb/neb_drv.F
- NWChem MEPGS: https://github.com/nwchemgit/nwchem-wiki/blob/master/mepgs.md
- pysisyphus: https://pysisyphus.readthedocs.io/en/latest/tsoptimization.html 、https://pysisyphus.readthedocs.io/en/latest/chainofstates.html 、https://pysisyphus.readthedocs.io/en/latest/irc.html 、https://pysisyphus.readthedocs.io/en/latest/calculators.html 、https://github.com/eljost/pysisyphus/blob/master/pysisyphus/calculators/QCEngine.py 、https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26390
- QCEngine / QCElemental: https://github.com/MolSSI/QCEngine/blob/master/qcengine/programs/nwchem/harvester.py 、https://github.com/MolSSI/QCEngine/blob/master/qcengine/programs/nwchem/runner.py 、https://molssi.github.io/QCElemental/model_molecule.html 、https://github.com/MolSSI/QCElemental/blob/main/qcelemental/models/molecule.py
- optking(座標系の扱い): https://optking.readthedocs.io/en/latest/optimizations.html
- IRC 実務(Maeda らのレビュー / Hratchian & Schlegel / QRC): https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24757 、https://pubs.acs.org/doi/10.1021/ct0499783 、https://www.sciencedirect.com/science/article/abs/pii/S0040403903021798
- 手法論文: CI-NEB(Henkelman 2000, J. Chem. Phys. 113, 9901)、IDPP(Smidstrup 2014, J. Chem. Phys. 140, 214106)、geodesic(Zhu 2019, J. Chem. Phys. 150, 164103)、GSM(Zimmerman 2013, JCTC 9, 3043)、TRIC(https://pubs.aip.org/aip/jcp/article/144/21/214108/313176/Geometry-optimization-made-simple-with-translation)、geomeTRIC engines(https://geometric.readthedocs.io/en/latest/engines.html)、Bofill(J. Comput. Chem. 15, 1 (1994))、dimer(Heyden 2005, J. Chem. Phys. 123, 224101)、ZTS(E, Ren, Vanden-Eijnden 2007, J. Chem. Phys. 126, 164103)
- 低障壁水素結合(Perrin & Nielson 1997): https://www.annualreviews.org/doi/pdf/10.1146/annurev.physchem.48.1.511
- 酸塩基錯体のイオン対: https://pubmed.ncbi.nlm.nih.gov/42490086/ 、https://pubmed.ncbi.nlm.nih.gov/34096305/ 、https://pubs.aip.org/aip/jcp/article-abstract/84/6/2953/793967 、https://arxiv.org/abs/2105.03078

**反応探索**
- ReaDuct 6.1.0: https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/HessianTask.h 、https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/TsOptimizationTask.h 、https://github.com/qcscine/readuct/blob/6.1.0/src/Readuct/App/Tasks/NtOptimization2Task.h 、https://scine.ethz.ch/static/download/manuals/v6.1.0/readuct_manual.pdf
- SCINE Utilities 10.1.0: https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Utils/GeometryOptimization/NtOptimizer2.cpp 、https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Utils/GeometryOptimization/AfirOptimizerBase.h 、https://github.com/qcscine/utilities/blob/10.1.0/src/Utils/Python/NormalModesPython.cpp
- xtb_wrapper 3.0.2: https://github.com/qcscine/xtb_wrapper/blob/3.0.2/src/Xtb/Xtb/Wrapper/XtbSettings.cpp
- Puffin 2.1.0: https://github.com/qcscine/puffin/blob/2.1.0/scine_puffin/jobs/templates/scine_react_job.py 、https://github.com/qcscine/puffin/blob/2.1.0/scine_puffin/jobs/templates/job.py
- Chemoton 2.0: https://arxiv.org/abs/2202.13011 、関連: https://www.nature.com/articles/s43588-021-00101-3
- autodE: https://github.com/duartegroup/autodE/blob/master/autode/config.py 、https://github.com/duartegroup/autodE/blob/master/autode/transition_states/base.py
- GRRM17/AFIR: https://pmc.ncbi.nlm.nih.gov/articles/PMC5765425/ 、https://wires.onlinelibrary.wiley.com/doi/full/10.1002/wcms.1538
- 気相 amine·HX(Legon 1993): https://pubs.rsc.org/en/content/articlelanding/1993/cs/cs9932200153 、NH3·HX と水(JPCA 1999): https://pubs.acs.org/doi/10.1021/jp991918j
- SiO2 の HF エッチング(Knotter 2000): https://pubs.acs.org/doi/abs/10.1021/ja993803z

**熱化学・速度論**
- GoodVibes: https://github.com/patonlab/GoodVibes 、https://doi.org/10.12688/f1000research.22758.1 、https://pypi.org/project/goodvibes/3.2/
- qRRHO(Grimme 2012): https://doi.org/10.1002/chem.201200497 、qh-H: https://doi.org/10.1021/jp509921r 、msRRHO: https://doi.org/10.1039/D1SC00621E 、https://doi.org/10.1021/jp205508z
- 標準状態補正(Bryantsev 2008): https://doi.org/10.1021/jp802665d
- スケール因子: https://comp.chem.umn.edu/freqscale/ 、https://cccbdb.nist.gov/vibscalejust.asp
- Wigner / トンネル: https://reactionmechanismgenerator.github.io/RMG-Py/reference/kinetics/wigner.html 、https://www.osti.gov/biblio/6485793 、Eckart: https://doi.org/10.1021/j150606a003 、レビュー: https://doi.org/10.1021/cr050205w
- Arkane: https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html 、https://reactionmechanismgenerator.github.io/RMG-Py/reference/arkane/essfactory.html
- Cantera species: https://cantera.org/stable/yaml/species.html
- 参照データ: https://doi.org/10.1021/j100007a034(HF2⁻)、https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664393&Mask=8(NIST)

### 5.3 ツール出力ファイルの場所

ディレクトリ: `C:/Users/user/AppData/Local/Temp/claude/c--Users-user-Desktop-HF-reaction-hfauto-final-baseline/425a764f-aa8a-4072-b205-3dc2c3d344ce/scratchpad/review/tools/`

| ファイル | 内容 |
|---|---|
| `versions.txt`、`commands.txt`、`pyproject_after.toml` | ツールのバージョン、実行コマンド、変更後の pyproject |
| `importlinter.txt`、`grimp_matrix.txt` | import-linter の全出力、依存行列・fan-in/out・循環 |
| `pyrefly.txt`、`pyrefly.json`、`pyrefly_min.txt`、`pyrefly_with_warnings_min.txt`、`pyrefly_count_errors.txt`、`pyrefly_stderr.txt` | pyrefly のエラーと警告 |
| `radon_cc.txt/.json`、`radon_cc_worst.txt`、`radon_cc_worst_tests.txt`、`radon_mi.txt/.json`、`radon_raw.txt/.json`、`radon_raw_tests.json`、`radon_summary.txt` | radon の CC・MI・raw と集計 |
| `ruff_default.txt`、`ruff_default_stats.txt`、`ruff_review.txt/.json`、`ruff_review_stats.txt`、`ruff_review_summary.txt`、`ruff_format_check.txt`、`ruff_format_diff.txt`、`ruff_format_files.txt` | Ruff の既定・レビュー用ルール・format の結果 |
| `pytest_collect.txt`、`pytest_collect_per_file.txt`、`test_files_by_size.txt`、`test_reference_map.txt`、`contract_surface.txt` | テストの収集結果、テスト規模、参照マップ、qc/data 契約面の集計 |

補助スクリプト(同じ `review/` 配下): `grimp_matrix.py`、`grimp_fanout.py`、`radon_summary.py`、`ruff_summary.py`、`stage_usage.py`(stage とパイプラインの使用照合)、`stage_cc.py`(stage 別の行数と最大 CC)、および各観点のサブディレクトリ(`a2/`、`a3/`、`a4/`、`a5/`、`c6/` など)。
