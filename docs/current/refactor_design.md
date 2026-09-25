# hfauto 抜本的リファクタリング設計書(2026-09-25、確定版)

- ブランチ: `refactor/fundamental-2026-09`(起点 `5d76201`)
- 根拠: `docs/reviews/2026-09-25_architecture_chemistry_review.md`(以下「レビュー」。指摘 ID の AR-xx / CH-xx / BUG-xx はレビューのもの)、棚卸し M1〜M5、設計案に対する化学批評(chem、24 項目)と工学批評(eng、27 項目)
- 本書は実装エージェント向けの**確定仕様**である。選択肢は示さない。型・関数シグネチャ・既定値は本書のとおりに実装する。変更が必要な場合は、先に本書を改訂する。
- 本書は移行期間中は `docs/current/` に置く。Wave 9 で `docs/reviews/2026-09-25_refactor_design.md` に移して凍結する。完成後の説明は `docs/design.md` に書く。

---

## 0. 批評を受けた改訂の要点

批評の各項目の採否は §13 に記す。本文は採用後の仕様である。

1. **層構造の違反を解消した(eng P0-1)。**
   - system 設定のモデルは `hfauto/core/system.py` に移した。
   - stage は `StageRuntime` Protocol(`stages/spec.py`)だけに依存する。
   - reporting は `render(view, out_dir)` の形にし、run の配置を知らない。
2. **並列化を単純にした(eng P0-2)。**
   - ProcessPool と run をまたぐ Budget 台帳は廃止した。
   - DFT と反応ケースは直列に実行する。スレッドで並列化するのは xTB と CREST だけとする。
   - 同時起動は site 単位の排他ロック 1 つで防ぐ。
3. **判断表を修正した(eng P0-4)。** 17 行にし、multi_max と monotonic の判定を種による分岐より前に置いた。
4. **Hessian を型で受け渡すようにした(eng P0-5、chem 2/20)。**
   - freq の Evidence は、全エンジン共通の正準形式(`.npy`、Eh/bohr²、入力フレーム)で Hessian を持つ。
   - saddle の初期 Hessian には、同じ構造で計算したものであれば、どの Level の Hessian でも使える。1 回目は xTB の Hessian を使う。
   - `frame_residual` は診断値として記録するだけで、ゲートには使わない。
5. **ゲートを opt・saddle の Evidence と結び付けた(chem 3)。**
   - freq の開始構造が opt の最終構造と一致していることを確認する。
   - 両者が同一 PES 上にあることを合否の条件にする。
   - QRC の判定は単調性ではなく、エネルギー差の閾値で行う(chem 7、eng 19)。
   - 変位幅はエネルギー目標から決める。
6. **熱化学は GoodVibes 4.3.0 の Python API に置き換えた(chem 5)。**
   - `goodvibes.api.compute_thermo(qcdata=QCData)` を worker 内で呼ぶ。
   - QCData は、自前で射影した振動数から組み立てる。
   - これでテキストのパース、振動ブロックの連結、古い CSV の混入が構造上起こらなくなる。
7. **振動解析を見直した(chem 6)。** 並進・回転空間は SVD で決め、補空間で対角化する。質量には主同位体の値を使う。
8. **汎用性を高めた(chem 8〜14)。**
   - 全エンジンに電荷と多重度を渡す。
   - trial は極性 H の移動、C–H の移動、重原子の結合組み換え、会合の順に作る。
   - H 結合を作れない錯体は、ランダムな剛体配置で配置する。
   - 開殻系では <S²> を観測する。
   - 会合熱化学を出す。
   - 置換不変 RMSD は scipy で計算し、真性回転だけを許し、初期配向を複数試す。
9. **CLI・設定・契約を減らした(eng 24)。**
   - CLI は 5 コマンドにした。
   - `--set` と、設定の優先順位によるマージを廃止した。
   - キャッシュは常に `<run>/jobs` に置く。
   - AttemptLog は書くだけにした。
10. **Wave を 9 段に組み直した(eng 11〜15)。**
    - 純関数の判断ロジックは W3 に前倒しした。
    - report は W4 に移した。
    - W4 の統合時に、WSL で実エンジンの smoke テストと HCN の実行を行う。
    - 切替と削除は 2 つの Wave に分けた。
    - 文書は最後の Wave で書く。

---

## 1. 目的とスコープ

### 1.1 目的

hfauto は、気相の分子・非共有結合錯体について、次の処理を manifest ベースで自動実行するワークフローである。

1. 候補生成物の発見: ReaDuct NT2/AFIR、無バイアス緩和、CREST のトポロジー停止。
2. 同一 PES 上での極小の検証(opt → 別ジョブの freq → 虚振動の処理)と、組成ごとの極小レジストリ。
3. 反応経路の決定: 障壁の事前判定 → saddle → TS の振動数検証 → QRC による接続確認。
4. 熱化学: GoodVibes の qRRHO を整合ゲート付きで実行する。失敗時は fail-closed とする。会合熱化学も含む。
5. 反応の分類と順位付け: 証拠の階層、不確かさの幅、手法間の符号一致で判断する。

本リファクタリングの目的は次の 4 点である。

1. 使われていない処理・説明・古い内容を削除する。
2. 過剰になった契約・診断・テストを抜本的に減らす。qc 辞書の 287 キー、41 stage、65 種の artifact 型を廃止または集約する。
3. 化学計算部分を作り直す。エンジンにも系にも依存しない汎用の実装にし、目的に適したもの(確実に収束し、科学的に正しく、できるだけ安価)にする。
4. 上記の結果として不要になったものをすべて削除する。

### 1.2 スコープ

| 対象 | 扱い |
|---|---|
| 気相の分子・非共有結合錯体の反応探索(HCN→HNC、HONO 異性化、NH3 反転、TMA/amine·(HF)n、水、ホルムアルデヒド) | 対象。汎用コードで扱い、系固有の分岐は置かない |
| 開殻系(多重度 > 1、NWChem は `odft`、xTB は `--uhf`、<S²> を観測) | 対象 |
| 閉殻一重項のジラジカル TS(broken symmetry) | **対象外**。検出もしない。利用者の責任で扱う |
| 会合熱化学(分離した単量体を基準にした ΔG_assoc と ΔG‡) | 対象(`system.compositions` に構成単量体があるとき) |
| 本番エンジン: NWChem 7.2.3、xTB 6.7.x、CREST 3.0.x、SCINE ReaDuct 6.1.0、pysisyphus 1.0(xTB ネイティブ計算器での低レベル経路のみ)、GoodVibes 4.3.0(Python API) | 対象 |
| 波動関数法の一点計算(NWChem の MP2 / CCSD(T)) | 手法パネルと参照計算に限って対象 |
| 溶媒 | `MethodSpec.solvation`(NWChem COSMO / xTB ALPB)を PES の同一性に含めるだけ。凝縮相キャンペーンの自動化は対象外 |
| wet-etch 反応器・表面モデル、公開 DB による同定、速度論(TST / Cantera / Arkane)、DB 支持スコアの校正、HPC の配列ジョブ、ORCA、NWChem NEB、pysisyphus+QCEngine 経路 | **対象外として削除する** |
| 過去の run(`runs/`) | ファイルとして読めればよい。新コードで読み込む互換層は作らない。golden に必要な出力だけを `tests/golden/data/` に抜粋する |

### 1.3 定量目標

| 指標 | 現状(2026-09-25 実測) | 目標(Wave 8 の受け入れ条件) |
|---|---:|---:|
| 本体の Python 行数(`hfauto` 37,608 + `hfauto_ops` 3,906 + `hfauto_viz` 2,426 + `scripts` 2,844) | 46,784 | **≤ 12,500**(上限 13,000、−73%) |
| テスト行数(`tests/**/*.py`) | 10,533 | ≤ 5,000(golden データは別計上で ≤ 1.2 MB、上限 2 MB) |
| `configs/` 行数 | 1,714 | ≤ 600 |
| 文書(レビューと本書を除く) | 約 1,230 行・18 ファイル | ≤ 600 行・4 ファイル |
| stage の型 | 41 | **8** |
| artifact の型 | 65(runs)/ 44(定数) | **8** |
| `Artifact.qc` のキー | 287 | **0**(フィールドごと廃止) |
| import-linter の契約 | 11(違反 7) | 5(違反 0) |
| 関数の循環的複雑度 | 最大 186 | **最大 15**(`ruff --select C901`、`max-complexity = 15`) |
| 1 ファイルの行数 | 最大 1,296 | ≤ 500 |
| 既定のテスト実行(`-m "not real"`) | 270 passed / 2 skipped、2.1 s | 15 s 未満、失敗 0 |

---

## 2. 設計原則

1. **stage は薄いアダプタにする。**
   - 処理は「`StageRuntime` から型付きの入力とエンジンを得る → driver か純関数を呼ぶ → artifact を返す」の 3 段だけとする。
   - stage は外部出力をパースせず、判定式も持たない。
2. **化学の判定は純関数にする。**
   - `hfauto/chemistry/` は execution・backends・drivers・stages・pipeline を import しない。
   - ファイル IO を行うのは `xyz.py` と `xyz_trajectory.py`(xyz の読み書き)だけとする。その他のモジュールは numpy 配列と型付きレコードを受け取り、値か `Gate` を返す。
3. **証拠は 1 つの型付きモデルで表す。**
   - 外部ジョブ 1 本の観測事実を `Evidence`(frozen な pydantic モデル)で表す。
   - `Evidence` が存在すること自体が、正常終了・収束・観測した理論レベルが要求どおりであること・フレームが保持されていることを保証する。
   - 合否はゲート関数 1 か所で判定し、bool の旗として他へ写さない。
4. **バックエンドは能力 Protocol の実装とする。**
   - 処理は render → run → parse の 3 段で、戻り値は型付きの結果か `Failure` とする。バックエンドは Artifact を作らず、他のバックエンドを生成しない。
   - 型付きの registry(`backends/engines.py`)を import してよいのは pipeline だけである。stage と driver は `StageRuntime.engine()` を通して Protocol 型のエンジンを受け取る。
5. **実行は 1 か所に集める。**
   - subprocess を起動するのは `hfauto/execution/process.py` だけとする。
   - 再利用・再開・再試行は `JobStore` と `JobRunner` に一本化する。
   - 並行実行: DFT のジョブと反応ケースは直列に実行する。スレッドで並列化するのは xTB と CREST だけとする。run の同時起動は `SiteLock` で拒否する。
6. **契約と診断は減らす方向にだけ変える。**
   - 新しい旗・診断キー・監査用 stage を足さない。
   - 人間向けの詳細は stage ごとの `diagnostics.json` に書くだけにする。これを読むコードがないことは grep で検査する。
7. **互換層を作らない。** 旧 manifest・旧設定・旧 stage 名の読み替えはしない。
8. **ラップするより削除する。** 旧コードを包み込まずに消す。新コードは §9 の旧モジュールを import しない。
9. **化学的に価値のある能力は残す**(レビュー 4.5.7)。
   - 対象: minimum-mode-follow、多解像度での障壁なし判定、経路上の中間極小と子反応、端点の対称性の破れ、受容部位に基づく錯体配置、xTB による事前判定。
   - stage としては削除し、統合先で driver の action、純関数、seed の生成器として使う。
10. **fail-closed にする。**
    - ダミーエンジン、内部フォールバックによる熱化学、`fallback_to_dummy`、`global.mode` を本番から消す。
    - テスト用の偽物は `tests/fakes.py` にだけ置く。
11. **エンジンにも系にも依存しないようにする。**
    - 全アダプタに電荷と多重度を必ず渡す。
    - 振動解析と熱化学の入力には、自前で射影した振動数を使う。
    - エンジン名の文字列リテラルは `hfauto/backends/` と `configs/` の外に置かない。
    - HF・amine 固有の名前や分岐を持たない。
12. **層の向きを守る。** 依存の向きは `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core` とする。複数の層で使う型は、使う層のうち最も下の層に置く。

---

## 3. 目標パッケージ構成(完成形)

行数は目安で、合計は約 12,200 行である。

```
hfauto/
├─ core/                               ~1,100   依存は pydantic, numpy, PyYAML のみ
│  ├─ evidence.py        ~200  NEW  FileRef, Level, Geometry, Evidence, PathProfile, FailureKind, Failure
│  ├─ records.py         ~340  NEW  ArtifactType と 8 種の payload、ReactionTrial, CoordinateTerm ほか(§5.3)
│  ├─ manifest.py        ~170  NEW  Artifact, Manifest, load_manifest / save_manifest(原子的に書き込む)
│  ├─ method.py          ~150  NEW  MethodSpec, ExecutionSpec, EngineSite, ThermoSettings, Deadline, level_mismatches()
│  ├─ system.py           ~90  NEW  SpeciesInput, CompositionInput, ReactionInput, SystemConfig, Conditions, load_system()
│  └─ hashing.py / ids.py / units.py / constants.py   ~150  既存(constants に CM1_TO_HARTREE などを追加)
├─ chemistry/                          ~3,300   純関数(IO は xyz だけ)。依存は core, numpy, scipy(rdkit は smiles.py だけで任意)
│  ├─ xyz.py             ~150  既存 + Molecule, geometry_fingerprint(), composition_key()
│  ├─ geometry.py        ~120  既存(Kabsch)
│  ├─ xyz_trajectory.py  ~150  既存(多フレーム xyz、再標本化)
│  ├─ electronic_state.py ~155 既存
│  ├─ topology.py        ~240  NEW  ヒステリシス付き結合判定、断片、labile H、受容原子、ファンデルワールス半径、q、状態ラベル
│  ├─ identity.py        ~230  NEW  置換不変 RMSD(scipy Hungarian と真性回転 Kabsch、複数の初期配向)、比較、一意割付け、縮退、周期最近傍
│  ├─ vibrations.py      ~220  NEW  同位体質量表、並進・回転空間の SVD、補空間での射影振動数、回転定数、曲率、frame_residual(診断用)
│  ├─ modes.py           ~130  NEW  モードに沿った変位、重なり、QRC の変位幅、mode-follow の結果分類
│  ├─ profile.py         ~230  NEW  経路の形状、補間による HEI、gmax による収束・停滞の判定、ノード間隔
│  ├─ interpolation.py   ~120  NEW  写像付きの端点整列と IDPP
│  ├─ gates.py           ~330  NEW  Gate, Policy, zpe_hartree と 10 個のゲート関数(§5.4)
│  ├─ placement.py       ~280  NEW  錯体の配置 seed(H 結合の円錐と、剛体ランダム配置によるフォールバック)
│  ├─ trials.py          ~300  NEW  汎用の反応 trial 生成(4 段の優先順)、緩和による発見
│  ├─ selection.py        ~60  NEW  minima(dft) に回す構造の選択
│  ├─ hypotheses.py      ~200  NEW  反応仮説(極小ペア)の選択と端点構造の選択
│  ├─ classification.py  ~150  NEW  反応ケースの終端分類
│  ├─ thermo.py          ~230  NEW  スケール因子表、composite G、標準状態の換算、アンサンブル、Boltzmann 分布、会合量
│  └─ smiles.py           ~60  NEW  SMILES から 3D 構造(RDKit ETKDG、任意依存)
├─ backends/                           ~2,900
│  ├─ protocols.py       ~220  NEW  Capability, Requirements, 6 つの能力 Protocol, 結果型
│  ├─ engines.py          ~90  NEW  型付きの遅延 registry(8 エントリ)
│  ├─ nwchem/          ~1,150  NEW  input.py(renderer)/ output.py(parser)/ engine.py(QM・string・saddle・WFT の SP)
│  ├─ xtb.py             ~260  NEW  QM(energy / opt / hess)
│  ├─ pysis/             ~280  NEW  engine.py + worker.py(xTB ネイティブ計算器で GS と、同じ入力の中での TSOpt)
│  ├─ crest.py           ~230  NEW
│  ├─ readuct/           ~480  NEW  engine.py + worker.py(NT2 / AFIR を worker サブプロセスで実行)
│  └─ goodvibes/         ~220  NEW  engine.py + worker.py(GoodVibes の API を worker で呼ぶ)
├─ execution/                          ~630
│  ├─ process.py         ~200  run_command(監視フック、プロセスツリーの停止、command_result.json)
│  ├─ jobstore.py        ~170  内容アドレス型キャッシュ(<run>/jobs)、ジョブ単位の排他、終端的な失敗の記憶
│  ├─ jobs.py            ~170  Task, Adapter, JobRunner(attempt ladder とコア数のセマフォ)、thread_map
│  ├─ lock.py             ~40  SiteLock(run の同時起動を拒否する)
│  └─ worker.py           ~50  in-process ライブラリ(SCINE, pysisyphus, GoodVibes)を隔離して実行する
├─ drivers/                            ~1,180
│  ├─ minimum.py         ~280  relax_to_minimum(opt → freq → mode-follow)と Registry
│  └─ reaction_case/     ~900  state.py(CaseState と純粋な decide)/ actions.py / driver.py
├─ stages/                             ~1,600
│  ├─ spec.py ~100(StageConfig, StageSpec, StageRuntime Protocol, Stage)/ catalog.py ~30
│  ├─ structures.py ~150 / conformer_search.py ~230 / minima.py ~230 / explore.py ~250
│  └─ reaction_paths.py ~180 / single_point.py ~90 / thermochemistry.py ~230 / report.py ~110
├─ pipeline/                           ~630   config.py / layout.py / runner.py / preflight.py
├─ reporting/                          ~650   summary.py(順位・被覆率・手法パネル)/ html.py(render(view, out_dir))
└─ cli/main.py                         ~150   doctor / run / status / report / case
tests/   conftest.py, fakes.py, unit/, golden/(data/ + SOURCES.json), integration/, smoke/(-m real)
configs/ sites/, methods/, systems/(+ systems/xyz/), pipelines/(discover, known_endpoints, method_panel)
docs/    README.md(ルート), design.md, environment.md, validation.md, reviews/(凍結)
```

削除するトップレベル:
- `hfauto_ops/`、`hfauto_viz/`、`scripts/`
- `hfauto/hpc/`、`hfauto/workflow/`、`hfauto/data/`、`hfauto/core/schemas/`
- `examples/`(xyz は `configs/systems/xyz/` へ移す)
- リポジトリ直下の ReaDuct 生成物(`afir_biased/` ほか)

---

## 4. stage 一覧(残す・統合・削除)

### 4.1 目標とする stage(8 種)

| # | stage 名(catalog のキー) | モジュール | 役割 | 使う能力 | 入力 → 出力(artifact 型) |
|---|---|---|---|---|---|
| 1 | `structures` | `stages/structures.py` | system ファイルの化学種(xyz / SMILES)を読み、電子状態と、宣言された反応の端点の原子順序を検証する | — | system → `species` |
| 2 | `conformers` | `stages/conformer_search.py` | 単量体の配座探索、錯体の配置 seed、CREST `--nci`、重複の除去、状態ラベル、上位構造の選択 | conformers | `species` → `species` |
| 3 | `minima` | `stages/minima.py` | 構造の選択 → `relax_to_minimum` → レジストリへの登録。`level: screen`(xTB)と `level: dft`(NWChem)の 2 回使う。± 変位が別々の極小に落ちた場合は、`mechanism="mode_follow"` の `discovery` を出す | qm | `species`, `discovery` → `calculation`, `minimum`, `species`, `discovery` |
| 4 | `explore` | `stages/explore.py` | 反応 trial の生成と ReaDuct NT2/AFIR の実行、緩和による発見の抽出、陰性結果の記録 | discovery | `minimum`(screen), `species` → `discovery`, `species` |
| 5 | `reaction-paths` | `stages/reaction_paths.py` | 反応仮説を選び、ReactionCaseDriver で処理する(障壁の事前判定 → saddle → TS 検証 → QRC → 分類) | qm, path, saddle | `minimum`(dft), `discovery`, `species`, system の宣言反応 → `reaction`, `calculation`, `minimum`, `species` |
| 6 | `sp` | `stages/single_point.py` | 固定構造での一点計算(composite 用のエネルギー層、手法パネル) | qm | `minimum`, `reaction` → `calculation` |
| 7 | `thermo` | `stages/thermochemistry.py` | GoodVibes の API、整合ゲート、composite G、感度の幅、会合量、basin 分布 | thermo | `minimum`, `reaction`, `calculation` → `species_thermo`, `reaction_thermo` |
| 8 | `report` | `stages/report.py` | 順位付け、探索の被覆率、手法パネル(符号の一致判定を含む)、表と HTML | — | 全 artifact → `report` |

- stage のディレクトリ名には、パイプライン YAML の安定した `id:` を使う。番号の接頭辞は付けない(AR-07)。
- `reaction-paths` は、デバッグ用に設定 `reaction_ids: [..]` を受け付ける。`scripts/run_*` の代わりはこれで足りる。

### 4.2 旧 41 stage の行き先

| 旧 stage | 判定 | 行き先 |
|---|---|---|
| ingest、enumerate-states(電子状態の検証だけ) | 統合 | `structures`。互変異性体・プロトマーの列挙(`chemistry/microstates.py`)は削除する |
| build-complexes(a: xyz の読み込み) | 統合 | `structures` |
| conformers、build-complexes(b: 配置と CREST)、detect-sites、endpoint-seeds、endpoint-seed-screen | 統合 | `conformers`。配置は `chemistry/placement.py`、受容部位は `topology.acceptor_atoms`、端点の対称性の破れは配置の乱数回転 seed と minima の mode-follow が引き継ぐ |
| preopt、dft-minima、minimum-mode-follow、minimum-mode-assess、minimum-registry | 統合 | `minima`(`drivers/minimum.py`) |
| generate-reactions、explore-reactions、relaxation-discovery、discovery-audit(陰性結果の分類) | 統合 | `explore` |
| connect-minima、reaction-plan、ts-search、path-ensemble、path-ensemble-assess、path-intermediates、irc、reaction-classify、reaction-segments | 統合 | `reaction-paths`(`drivers/reaction_case/`) |
| sp | 書き直し | `sp`(engine を必須にし、`methods:` のリストを受け付ける) |
| thermo、thermo-sensitivity、basin-populations | 統合 | `thermo` |
| reaction-rank、method-panel(集計部分)、discovery-audit(集計部分) | 統合 | `report` |
| recover-path、rank | 削除 | JobStore の継続実行と `report` で不要になる |
| enrich、descriptors、kinetics、calibrate、connector-audit、viz、hpc-plan、ops | 削除 | 目的外(レビュー 4.5.3) |

### 4.3 パイプライン(13 本 → 3 本)

```yaml
# configs/pipelines/discover.yaml
pipeline_id: discover
conditions: {temperatures_K: [298.15], standard_states: [1atm]}
stages:
  - {id: structures, stage: structures}
  - {id: conformers, stage: conformers, engine: crest, method: gfn2}
  - {id: screen,     stage: minima, level: screen, engine: xtb, method: gfn2}
  - {id: explore,    stage: explore, engine: readuct, method: gfn2}
  - {id: dft,        stage: minima, level: dft, engine: nwchem, method: pbe0-d3bj_def2-svpd,
     select: {include: window, per_state: 3, window_kcal: 6.0, rerank_sp: true},
     init_hessian: {engine: xtb, method: gfn2}}
  - {id: paths,      stage: reaction-paths, method: pbe0-d3bj_def2-svpd,
     engines: {qm: nwchem, saddle: nwchem_saddle, path: nwchem_string},
     screen: {method: gfn2, qm: xtb, path: pysis_gs}}
  - {id: thermo,     stage: thermo, engine: goodvibes}
  - {id: report,     stage: report}
```

- `known_endpoints.yaml`:
  - 流れは `structures` → `dft`(minima。`select: {include: all}`、`init_hessian` は discover と同じ)→ `paths` → `thermo` → `report`。
  - `paths` の `screen:` は discover と同じにする。screen では、xTB で DFT の端点を再最適化してから使う。
- `method_panel.yaml`:
  - 流れは `panel_sp`(`stage: sp`、`methods: [...]`、`targets: reaction_stationary_points`)→ `panel_report`(`stage: report`)。
  - 既存の run ディレクトリに `--run-dir` を指定して**追記**する形で実行する(§7.4 の view の定義で成り立つ)。
- 旧 13 本は、system ファイル 7 本(hcn, hono, nh3_inversion, formaldehyde, tma_hf2, amine_hf_panel, water_same_basin)と上の 3 本に置き換える。
- 手法を変える場合は、別の method ファイルを指定する(`--set` は廃止)。xfine を使いたい場合は、利用者が xfine の method ファイルを作る。

---

## 5. 型付き Evidence とゲート関数

### 5.1 `hfauto/core/evidence.py`(Wave 2)

```python
_FROZEN = ConfigDict(frozen=True, extra="forbid")

class FileRef(BaseModel):
    model_config = _FROZEN
    path: str        # run ディレクトリからの相対 POSIX パス(JobStore も <run>/jobs の下にある)
    sha256: str

class Level(BaseModel):
    """出力から観測した理論レベル(要求値ではない)。文字列は小文字に正規化する。"""
    model_config = _FROZEN
    program: str                     # "nwchem" | "xtb" | "fake"
    version: str                     # 出力から観測した版数("7.2.3")
    method: str                      # "pbe0" | "gfn2" | "mp2" | "ccsd(t)"
    basis: str | None = None         # "def2-svpd"。cartesian 基底なら "def2-svpd/cart"
    dispersion: str | None = None    # "d3zero" | "d3bj" | "d4" | None
    solvation: str | None = None     # "cosmo:78.4" | "alpb:water" | None
    charge: int
    multiplicity: int
    grid: str | None = None          # 数値精度の層
    scf_tol: float | None = None     # 数値精度の層
    electronic_temperature_K: float | None = None   # xTB だけ(PES の一部)
    def surface_key(self) -> str: ...  # grid と scf_tol を除いた正規化 JSON の sha256[:16]
    def full_key(self) -> str: ...     # 全項目の sha256[:16]

class Geometry(BaseModel):
    model_config = _FROZEN
    file: FileRef                    # xyz
    fingerprint: str                 # chemistry.xyz.geometry_fingerprint(元素列と、1e-6 Å に丸めた座標)。core では計算しない
    symbols: tuple[str, ...]

Task = Literal["sp", "opt", "freq", "saddle"]

class Evidence(BaseModel):
    """正常終了した外部ジョブ 1 本の観測事実。存在すること自体が次の 6 点を保証する。
    (1) 正常終了し SCF が収束した  (2) opt / saddle なら構造最適化が収束した
    (3) freq なら 3N − n_external 本の振動数がそろっている  (4) 原子の順序が保たれている
    (5) 入力フレームが保たれている(出力に出た初期構造が入力と 1e-4 Å 以内で一致する)
    (6) 観測した Level が、要求した MethodSpec と site の版数 pin に一致する
    どれか 1 つでも満たさなければ、バックエンドは Failure を返す。"""
    model_config = _FROZEN
    kind: Literal["calculation"] = "calculation"
    engine: str
    task: Task
    level: Level
    start: Geometry
    final: Geometry                                   # sp / freq では start と同じ
    energy_hartree: float                             # final での電子エネルギー
    trajectory_energies_hartree: tuple[float, ...] = ()   # opt / saddle の各ステップ(先頭は start のエネルギー)
    frequencies_cm1: tuple[float, ...] | None = None  # freq だけ。final での射影済み振動数(3N − n_external 本)。虚数は負値で表す
    n_external: Literal[5, 6] | None = None           # freq だけ。並進・回転の自由度の数(直線分子なら 5)
    imaginary_modes: tuple[tuple[float, ...], ...] = ()   # ν < 0 のモード。規格化した直交座標の変位(3N 成分、入力フレーム)
    hessian: FileRef | None = None                    # freq だけ。正準形式: .npy、(3N, 3N)、Eh/bohr²、入力フレーム
    s2: float | None = None                           # 観測した <S²>(開殻のときだけ)
    output: FileRef                                   # 主出力ファイル
    job_key: str

class PathProfile(BaseModel):
    model_config = _FROZEN
    kind: Literal["path"] = "path"
    engine: str
    level: Level
    images: FileRef                                   # 多フレームの xyz
    energies_hartree: tuple[float, ...]
    gmax_history: tuple[float, ...] = ()              # 反復ごとの最大勾配(Eh/bohr)
    program_converged: bool                           # プログラム自身の判定。採否は profile.string_converged が決める
    climbing_image: int | None = None
    ts: Geometry | None = None                        # エンジンが同じ入力の中で TS 最適化まで行った場合(pysis_gs)
    ts_energy_hartree: float | None = None
    job_key: str

class FailureKind(StrEnum):
    EXECUTABLE_MISSING = "executable_missing"
    INPUT_INVALID = "input_invalid"           # AUTOZ の失敗、Hessian と構造の不一致など。詳細は reason に書く
    TIMEOUT = "timeout"
    STAGNATED = "stagnated"                   # 監視フックによる打ち切り
    NONZERO_EXIT = "nonzero_exit"
    SCF_NOT_CONVERGED = "scf_not_converged"
    GEOMETRY_MAXITER = "geometry_maxiter"
    INCOMPLETE_OUTPUT = "incomplete_output"   # パースできない、Hessian がない、振動数の本数が合わない
    METHOD_MISMATCH = "method_mismatch"       # 観測した Level が要求と違う(版数の違い、フレームの変化を含む)
    BUDGET_EXHAUSTED = "budget_exhausted"
    GATE_REJECTED = "gate_rejected"           # 化学ゲートで不合格(stage が failed artifact を出すときだけ使う)

class Failure(BaseModel):
    model_config = _FROZEN
    kind: FailureKind
    reason: str
    job_key: str | None = None
```

- `frequencies_cm1` と `hessian`(正準形式の `.npy`)は、どのエンジンでも必ず `chemistry.vibrations` で計算・変換する。エンジン自身が出した振動数は、照合用に `diagnostics.json` に書くだけにする(CH-26 の根本対策)。
- `Level` は出力から観測する。NWChem では次を読む(CH-02)。
  - `NWChem) 7.2.3` の版数
  - `Grid used for XC integration`
  - `Convergence on energy requested`
  - `(spherical)` / `(cartesian)`
  - DFT-D3 ブロック、COSMO ブロック
  - `<S2> =` の行
- `frame_residual`(|H·R_i|)は、診断用に `diagnostics.json` に書くだけにする。非停留点では勾配の寄与によって大きくなるのが正常なので(G08 では 1.8e-2)、ゲートには使わない。フレームが保たれていることは、上の保証 (5) で担保する。

### 5.2 `hfauto/core/method.py` と `hfauto/core/system.py`(Wave 2)

```python
class MethodSpec(BaseModel):       # configs/methods/*.yaml。要求値で、JobStore のキーに入る
    model_config = ConfigDict(frozen=True, extra="forbid")
    id: str
    kind: Literal["dft", "xtb", "wft"]
    functional: str | None = None           # dft
    wft_method: Literal["mp2", "ccsd(t)"] | None = None   # wft(SP だけ)
    basis: str | None = None                # dft / wft
    dispersion: Literal["d3zero", "d3bj", "d4"] | None = None
    gfn: Literal[1, 2] | None = None        # xtb
    electronic_temperature_K: float | None = None   # xtb(None なら 300)
    grid: Literal["coarse", "medium", "fine", "xfine"] | None = None
    scf_energy_tol: float | None = None
    solvation: str | None = None            # "cosmo:<eps>" | "alpb:<solvent>"
    def signature(self) -> dict[str, object]: ...   # id 以外の全項目(キャッシュキー用)

def level_mismatches(requested: MethodSpec, observed: Level, *, version_pin: str) -> list[str]: ...  # 空なら一致

class ExecutionSpec(BaseModel):    # site ファイル由来。JobStore のキーには入れない
    model_config = ConfigDict(frozen=True, extra="forbid")
    ranks: int = 1
    threads: int = 1
    memory_mb_per_rank: int = 1200
    timeout_s: float = 14_400
    maxiter: int | None = None
    env: dict[str, str] = {}
    coordinates: Literal["auto", "cartesian"] = "auto"

class EngineSite(BaseModel):       # site ファイルのエンジン別設定。engines.create に渡す
    model_config = ConfigDict(frozen=True, extra="forbid")
    version: str                           # 版数の pin(必須)。JobStore のキーに入れ、観測値と照合する
    executables: dict[str, str] = {}       # {"nwchem": "...", "mpirun": "..."}。絶対パスは site ファイルにだけ書く
    execution: ExecutionSpec = ExecutionSpec()
    python: str | None = None              # worker 用のインタプリタ(None なら sys.executable)
    scratch_dir: str | None = None         # ext4 上の絶対パス

class ThermoSettings(BaseModel):   # GoodVibes の型付き設定(extra_args は廃止)
    model_config = ConfigDict(frozen=True, extra="forbid")
    temperatures_K: tuple[float, ...] = (298.15,)
    qs: Literal["grimme", "truhlar"] = "grimme"
    cutoff_cm1: float = 100.0
    vib_scale: float | None = None          # None なら chemistry.thermo.scale_factors(level) を使う
    zpe_scale: float | None = None
    symmetry: bool = True
    invert_soft_cm1: float | None = None    # soft 極小だけに使う。TS には付けない
    sensitivity: bool = True                # qs × cutoff{50, 100, 150} の幅を出す

@dataclass
class Deadline:
    end: float                                   # time.monotonic() 基準
    @classmethod
    def after(cls, seconds: float) -> "Deadline": ...
    def remaining(self) -> float: ...
    def expired(self) -> bool: ...
```

```python
# hfauto/core/system.py(system ファイルのモデル。xyz は system ファイルのディレクトリからの相対パス)
class SpeciesInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    xyz: Path | None = None
    smiles: str | None = None
    charge: int = 0
    multiplicity: int = 1
    role: Literal["monomer", "endpoint"] = "monomer"

class CompositionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    components: dict[str, int]          # species id → 個数(例: {"tma": 1, "hf": 2})。電荷と多重度は構成要素から決める

class ReactionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    reactant: str                       # role=endpoint の species id。原子順序は product と同じでなければならない
    product: str
    coordinate: list[CoordinateTerm] = []   # core.records.CoordinateTerm
    torsional: bool | None = None       # None なら、結合変化がないとき True とする

class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system_id: str
    species: list[SpeciesInput]
    compositions: list[CompositionInput] = []
    reactions: list[ReactionInput] = []

class Conditions(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    temperatures_K: tuple[float, ...] = (298.15,)
    standard_states: tuple[Literal["1atm", "1bar", "1M"], ...] = ("1atm",)

def load_system(path: Path) -> SystemConfig: ...   # xyz を絶対パスに解決する
```

版数の pin は site 側(`EngineSite.version`)に置き、`MethodSpec` には持たせない。

### 5.3 `hfauto/core/records.py` と `hfauto/core/manifest.py`(Wave 2)

artifact の型は 8 種に固定する。payload は `kind` で判別する pydantic のタグ付き共用体とし、`data` / `qc` / `method` / `provenance` / `paths` の自由辞書は廃止する。すべて `frozen=True, extra="forbid"` とする。

```python
class ArtifactType(StrEnum):
    SPECIES = "species"; CALCULATION = "calculation"; MINIMUM = "minimum"; DISCOVERY = "discovery"
    REACTION = "reaction"; SPECIES_THERMO = "species_thermo"; REACTION_THERMO = "reaction_thermo"; REPORT = "report"

class SpeciesRecord(BaseModel):
    kind: Literal["species"] = "species"
    species_id: str
    composition_id: str               # chemistry.xyz.composition_key の値
    formula: str
    charge: int
    multiplicity: int
    geometry: Geometry
    source: Literal["input", "conformer", "placement", "crest_topology", "discovery",
                    "mode_follow", "connection", "intermediate"]
    state_label: str                  # topology.state_label(汎用。HF 固有の名前は使わない)
    energy_hartree: float | None = None   # 低レベルでの参考エネルギー
    level_key: str | None = None

class MinimumRecord(BaseModel):
    kind: Literal["minimum"] = "minimum"
    minimum_id: str                   # basin を代表する id
    basin_id: str
    composition_id: str
    species_id: str                   # 代表構造の species
    tier: Literal["screen", "dft"]
    level_key: str                    # Level.full_key()
    opt_calc: str                     # calculation artifact の id
    freq_calc: str
    energy_hartree: float
    state_label: str
    n_fragments: int
    members: tuple[str, ...] = ()     # この basin に落ちた species の id(崩壊した seed を含む)
    notes: tuple[str, ...] = ()       # soft_imaginary_mode, noise_imaginary_mode, spin_contaminated

class ReactionTrial(BaseModel):
    trial_id: str
    source_minimum: str
    kind: Literal["polar_h", "h_shift", "heavy_bond", "association"]
    mechanism: Literal["nt2", "afir"]
    associations: tuple[tuple[int, int], ...] = ()
    dissociations: tuple[tuple[int, int], ...] = ()
    perturbed: bool = False           # 直線分子を曲げた、またはランダムに変位させた

class DiscoveryRecord(BaseModel):
    kind: Literal["discovery"] = "discovery"
    discovery_id: str
    source_minimum: str
    mechanism: Literal["nt2", "afir", "relaxation", "mode_follow"]
    trial: ReactionTrial | None = None
    outcome: Literal["product", "negative", "failed"]
    reason: str | None = None         # negative の理由: monotonic_uphill | collapsed_to:<label> | out_of_window | irc_not_connected_to_source
    product_species: str | None = None
    ts: Geometry | None = None        # 低レベルの TS(射影振動数で虚振動 1 本を確認済み)。REFINE_SADDLE の種になる
    ts_imag_cm1: float | None = None
    barrier_kj_mol: float | None = None      # 300 K での値(SCC 失敗時の再試行でも 300 K で評価し直す)
    reaction_kj_mol: float | None = None
    electronic_temperature_K: float = 300.0

class StoichTerm(BaseModel):
    composition_id: str
    coefficient: int

class CoordinateTerm(BaseModel):
    kind: Literal["distance", "angle", "dihedral"]
    atoms: tuple[int, ...]
    coefficient: float = 1.0

class BarrierVerdict(BaseModel):
    verdict: Literal["proceed", "barrierless", "negative_evidence", "unavailable"]
    max_rel_low_kcal: float | None = None
    max_rel_dft_kcal: float | None = None
    n_dft_points: int = 0
    max_node_spacing_A: float | None = None   # 「解像度」の値として記録する
    below_zpe: bool = False
    seed: Geometry | None = None      # proceed のとき: 低レベルの TS、または DFT プロファイル上の HEI ノード
    reasons: tuple[str, ...] = ()

class SaddleClaim(BaseModel):
    saddle_calc: str
    freq_calc: str
    imag_cm1: float
    energy_hartree: float
    notes: tuple[str, ...] = ()       # soft_secondary_mode, mode_overlap_below_0.3, spin_contaminated

class ConnectionClaim(BaseModel):
    method: Literal["qrc"] = "qrc"
    side_calcs: tuple[str, str]
    minima: tuple[str, str]
    amplitude_A: float
    notes: tuple[str, ...] = ()

class CaseOutcome(StrEnum):
    ELEMENTARY_STEP = "elementary_step"
    DEGENERATE = "degenerate_rearrangement"
    REASSIGNED = "reassigned_step"
    MULTI_STEP = "multi_step"
    BARRIERLESS = "barrierless_at_resolution"
    SAME_BASIN = "same_basin"
    NO_PRODUCT = "no_product_basin"
    OUT_OF_WINDOW = "out_of_window"
    UNRESOLVED = "unresolved_within_budget"
    BLOCKED = "blocked_upstream"

class ReactionRecord(BaseModel):
    kind: Literal["reaction"] = "reaction"
    reaction_id: str
    parent_id: str | None = None
    reactants: tuple[StoichTerm, ...]
    products: tuple[StoichTerm, ...]
    minima: tuple[str, str]           # (反応物側, 生成物側) の minimum_id。縮退反応では同じ id が 2 つ並ぶ
    endpoints: tuple[str, str]        # 経路計算に使う構造の species_id(hypotheses.pick_endpoints で原子写像が物理的に正しいものを選ぶ)
    degenerate: bool = False          # identity.mapped_equivalent(両端) が真(NH3 反転など。CH-35)
    source: Literal["declared", "discovery", "mode_follow", "conformer", "split", "reassigned"]
    coordinate: tuple[CoordinateTerm, ...] = ()
    torsional: bool = False           # 宣言、または結合変化がないときに自動で真
    n_h_transferred: int = 0          # 結合相手が変わる H の数(rankable の許容幅に使う)
    low_level_ts: Geometry | None = None   # discovery / mode_follow が見つけた種
    negative_evidence: tuple[str, ...] = ()
    barrier: BarrierVerdict | None = None
    saddle: SaddleClaim | None = None
    connection: ConnectionClaim | None = None
    outcome: CaseOutcome | None = None
    reasons: tuple[str, ...] = ()
    log: str | None = None            # cases/<reaction_id>/log.jsonl(相対パス、書くだけ)

class SpeciesThermo(BaseModel):
    kind: Literal["species_thermo"] = "species_thermo"
    subject: str                      # minimum_id、または SaddleClaim.freq_calc
    freq_calc: str
    energy_calc: str | None = None    # composite に使った sp
    T_K: float
    G_hartree: float | None           # None は thermo_unavailable を表す
    H_hartree: float | None
    zpe_hartree: float | None
    settings_sha: str
    population: float | None = None   # 同じ組成の中での Boltzmann 重み
    notes: tuple[str, ...] = ()       # scale_factor_from:<level>, scale_factor_unverified, soft_imaginary_mode ...

class ReactionThermo(BaseModel):
    kind: Literal["reaction_thermo"] = "reaction_thermo"
    reaction_id: str
    T_K: float
    standard_state: Literal["1atm", "1bar", "1M"]
    dE_act_kcal: float | None
    dE_rxn_kcal: float | None
    dzpe_act_kcal: float | None                    # ΔE0‡ − ΔE‡(rankable が使う)
    dG_act_kcal: float | None
    dG_rxn_kcal: float | None
    dG_assoc_kcal: float | None = None             # 反応物錯体 − Σ 単量体(標準状態の補正を含む)
    dG_act_vs_separated_kcal: float | None = None  # TS − Σ 単量体
    band_kcal: tuple[float, float] | None = None   # 感度の幅(qs × cutoff)
    blockers: tuple[str, ...] = ()                 # thermo_unavailable, mixed_level_of_theory, spin_contaminated

class RankRow(BaseModel):
    reaction_id: str
    outcome: CaseOutcome
    tier: Literal["screening", "minima", "saddle", "connected"]
    rankable: bool
    rank: int | None
    dG_act_kcal: float | None
    band_kcal: tuple[float, float] | None
    blockers: tuple[str, ...]

class ReportRecord(BaseModel):
    kind: Literal["report"] = "report"
    rows: tuple[RankRow, ...]
    tables: dict[str, FileRef]        # ranking.csv, coverage.csv, method_panel.csv, report.html

Payload = Annotated[Evidence | SpeciesRecord | MinimumRecord | DiscoveryRecord | ReactionRecord
                    | SpeciesThermo | ReactionThermo | ReportRecord, Field(discriminator="kind")]
```

```python
# hfauto/core/manifest.py
class Artifact(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    artifact_id: str
    type: ArtifactType
    parents: tuple[str, ...] = ()
    status: Literal["success", "failed"] = "success"
    payload: Payload | None = None
    failure: Failure | None = None
    # validator: success のとき、payload が None でなく payload.kind == type であること
    #            failed のとき、failure が None でないこと(payload は任意)

class Manifest(BaseModel):
    schema_version: Literal["hfauto.manifest.v2"] = "hfauto.manifest.v2"
    run_id: str
    stage_id: str
    created_at: str
    artifacts: list[Artifact] = []
    def of(self, type: ArtifactType, *, ok_only: bool = True) -> list[Artifact]: ...
    def get(self, artifact_id: str) -> Artifact: ...                  # 同じ id があれば後勝ち
    def records(self, type: ArtifactType, model: type[T]) -> list[T]: ...
    def evidence(self, calc_id: str) -> Evidence: ...
    @classmethod
    def union(cls, manifests: Sequence["Manifest"], *, run_id: str, stage_id: str) -> "Manifest": ...

def load_manifest(path: Path) -> Manifest: ...
def save_manifest(manifest: Manifest, path: Path) -> Path: ...     # tempfile + os.replace で原子的に書き込む
```

- 各 stage の manifest には、**その stage の出力だけ**を入れる。下流の入力は `RunLayout.view(stage_id)` で得る(§7.4)。carry_forward は廃止する。
- `status` は `success | failed` の 2 値とする。`partial` と `skipped` は廃止する。

### 5.4 ゲート関数 `hfauto/chemistry/gates.py`(Wave 2)

import してよいのは標準ライブラリ、numpy、`hfauto.core.*` だけである。

```python
@dataclass(frozen=True)
class Gate:
    ok: bool
    reasons: tuple[str, ...] = ()     # 不合格の理由
    notes: tuple[str, ...] = ()       # 合格したときの注記
    def __bool__(self) -> bool: return self.ok

@dataclass(frozen=True)
class Policy:                          # パイプライン YAML の `gates:` で上書きできる唯一の閾値群
    noise_cm1: float = 10.0            # |ν| < 10 は数値ノイズとみなす
    saddle_cm1: float = 50.0           # ν < −50 は本物の虚振動とみなす
    torsion_saddle_cm1: float = 20.0   # ねじれ反応の TS では ν < −20 を許す
    ts_prominence_hartree: float = 2.0e-5
    scf_noise_factor: float = 20.0     # prominence と降下量の下限は max(既定値, 20 × SCF 閾値)
    qrc_min_drop_hartree: float = 1.0e-5
    barrier_proceed_kcal: float = 1.0
    reaction_window_kcal: float = 40.0
    thermo_zpe_tol_hartree: float = 1.0e-5
    thermo_energy_tol_hartree: float = 1.0e-6
    rank_dzpe_base_kcal: float = 1.0   # |ΔE0‡ − ΔE‡| ≤ base + per_h × max(1, n_h_transferred)
    rank_dzpe_per_h_kcal: float = 4.0
    spin_contamination_tol: float = 0.1

def zpe_hartree(freqs_cm1: Sequence[float], *, scale: float = 1.0, invert_cm1: float | None = None) -> float
    # 0.5 × scale × Σ ν × CM1_TO_HARTREE。対象は ν > 0 と、反転したモード(−invert_cm1 < ν < 0)の |ν|。
    # ZPE の式はここ 1 か所だけに置き、chemistry.thermo もこれを使う
def imaginary_tier(freqs_cm1: Sequence[float], policy: Policy = Policy()) -> Literal["none", "noise", "soft", "saddle"]
def same_pes(*levels: Level, numerics: bool = True) -> Gate
def spin_ok(ev: Evidence, policy: Policy = Policy()) -> Gate          # s2 が None なら合格
def is_minimum(freq: Evidence, *, opt: Evidence, policy: Policy = Policy()) -> Gate
def is_first_order_saddle(freq: Evidence, *, saddle: Evidence, endpoint_energies: Sequence[float] = (),
                          torsional: bool = False, policy: Policy = Policy()) -> Gate
def barrier_verdict(dft_profile: Sequence[float], *, dft_endpoints: tuple[float, float],
                    low_profile: Sequence[float] | None = None, low_endpoints: tuple[float, float] | None = None,
                    negative_evidence: Sequence[str] = (), override: bool = False,
                    tangent_mode_cm1: float | None = None, seed: Geometry | None = None,
                    max_node_spacing_A: float | None = None, policy: Policy = Policy()) -> BarrierVerdict
def connection(ts_freq: Evidence, sides: tuple[Evidence, Evidence], assigned: tuple[str | None, str | None],
               expected: frozenset[str], *, degenerate: bool, sides_distinct: bool = True,
               policy: Policy = Policy()) -> tuple[Gate, Literal["elementary", "degenerate", "reassigned", "failed"]]
def thermo_consistent(freq: Evidence, *, gv_zpe_hartree: float, gv_energy_hartree: float, gv_n_real: int,
                      zpe_scale: float, invert_cm1: float | None, policy: Policy = Policy()) -> Gate
def rankable(reaction: ReactionRecord, thermo: ReactionThermo | None, *, participant_notes: Sequence[str] = (),
             policy: Policy = Policy()) -> Gate
def reaction_tier(reaction: ReactionRecord) -> Literal["screening", "minima", "saddle", "connected"]
```

| ゲート | 判定内容(確定) | 置き換える旧コード | 指摘 |
|---|---|---|---|
| `same_pes` | program・version・method・basis・dispersion・solvation・charge・multiplicity・electronic_temperature_K が一致すること。`numerics=True` なら grid と scf_tol も一致すること。停留点・freq・QRC・熱化学の層は `numerics=True` で判定する。DFT string の結果は種としてだけ使うので `numerics=False` で判定する | `method_lineage` の 5 関数、`path_method_signature`、`minimum_registry._method_signature`、`_validated_minimum_provenance`、`_nwchem_qcengine_evidence` | CH-02、CH-12、CH-16、CH-40、AR-29 |
| `spin_ok` | `s2` が None(閉殻)であるか、\|s2 − S(S+1)\| ≤ 0.1 であること。S = (多重度 − 1)/2。不合格なら注記 `spin_contaminated` とし、rankable の blocker にする | —(新規) | chem 11 |
| `is_minimum` | 次をすべて満たすこと。① `freq.task == "freq"` かつ `opt.task == "opt"`。② `freq.start.fingerprint == opt.final.fingerprint`。③ `same_pes(opt.level, freq.level, numerics=True)`。④ 振動数が 3N − n_external 本ある。⑤ `imaginary_tier` が saddle なら不合格、soft なら合格で注記 `soft_imaginary_mode`、noise なら合格で注記 `noise_imaginary_mode`。saddle 探索が極小に崩壊した構造も登録できるよう、opt の由来(task)は問わない | 約 10 実装(`minimum_promotion_gate`、`is_accepted_optimized_minimum`、`basin_identity._minimum_evidence` の OR 判定ほか) | CH-17、CH-20、AR-21、chem 3 |
| `is_first_order_saddle` | 次をすべて満たすこと。① `freq.start.fingerprint == saddle.final.fingerprint` かつ numerics を含めて同一 PES。② 振動数が 3N − n_external 本ある。③ 最も低いモードが ν < −50(`torsional` なら −20)。④ ほかのモードに ν < −50 がない(あれば不合格 `higher_order`)。⑤ ほかのモードに −50 ≤ ν < −10 があれば合格で、注記 `soft_secondary_mode` を付ける。⑥ `endpoint_energies` が与えられたら E_TS − max(E) ≥ max(2e-5 Eh, 20 × scf_tol)。モードの重なりは合否に使わない | `validate_ts_frequency_calculation`、`ts_qc`、検証済み TS の 11 か所 | CH-26、CH-33、CH-34、chem 3、eng 18 |
| `barrier_verdict` | 相対エネルギーは、両端の高い方(DFT 極小のエネルギー)を基準にする。① 陰性証拠があり、`override=False` なら `negative_evidence`。② `dft_profile` が 3 点未満なら `unavailable`。③ DFT の内部極大が 1.0 kcal/mol 以上、または低レベルの内部極大が 1.0 kcal/mol 以上なら `proceed`。④ どちらも 1.0 kcal/mol 未満なら `barrierless`(`max_node_spacing_A` を解像度として記録する)。`tangent_mode_cm1` が与えられ、DFT の極大 < ½hν なら `below_zpe=True` | `path_diagnostics` の単調 → ブラケット分岐、`recover-path`、`connect-minima` の min_change | CH-06、CH-15、CH-28、chem 16 |
| `connection` | 両側それぞれについて、drop = max(1e-5 Eh, 20 × scf_tol) として次を満たすこと。(a) `side.trajectory_energies_hartree[0]` ≤ E_TS − drop(変位直後の降下)。(b) 軌跡の最大値が E_TS − drop 以下で、最終値が先頭値 − drop 以下。(c) `same_pes(ts, side, numerics=True)`。(d) 非縮退の場合: 割付けが両側とも None でなく互いに異なり、その**集合**が `expected` と一致すれば `elementary`(向きは問わない)。別の登録極小の組につながれば `reassigned`、それ以外は `failed`。縮退の場合(`degenerate=True`、`expected` は 1 要素): 両側ともその basin に割り付き、`sides_distinct` が真なら `degenerate`。単調性の許容値は廃止した | irc の 7 段ゲート、`endpoint_pair_match_qc`、`_has_real_irc`、`_validated_elementary_step` | CH-03、CH-04、CH-35、CH-36、chem 7、eng 19 |
| `thermo_consistent` | 反転したモードの集合を inverted = {ν : −invert < ν < 0} とする。次を満たすこと。GoodVibes の実振動の本数 = count(ν > 0) + \|inverted\|。\|ZPE_GV − zpe_hartree(freqs, scale=zpe_scale, invert_cm1)\| < 1e-5 Eh。\|E_GV − E\| < 1e-6 Eh | `main_values_are_*`、`production_thermo_ready` ほか | CH-38、CH-39、chem 4 |
| `rankable` | 次をすべて満たすこと。outcome が elementary / degenerate / reassigned のいずれか。`thermo.dG_act_kcal` が None でない。blockers と participant_notes に thermo_unavailable、mixed_level_of_theory、spin_contaminated、method_sign_disagreement がない。\|dzpe_act\| ≤ 1 + 4 × max(1, n_h_transferred) kcal/mol | `quality_tier`、`confidence_score`、`scientific_/production_rank_eligible` | AR-27、CH-39、chem 19 |
| `reaction_tier` | minima(両端が dft 極小)→ saddle(`SaddleClaim` がある)→ connected(`ConnectionClaim` がある)の順に上がる。それ以外は screening | Q0〜Q5 と confidence | AR-27 |

**claim を作るのは 1 か所だけにする。**
- `MinimumRecord` を作るのは、`minima` stage と、`reaction-paths` で新しい basin を見つけたときの `Registry` である。
- `SaddleClaim` と `ConnectionClaim` を作るのは、`reaction-paths` の driver だけである。
- 下流(thermo と report)は claim の有無だけを見て、再検証しない。例外は thermo が行う `same_pes`(反応の参加者全体)と `thermo_consistent` の 2 つである。

### 5.5 化学カーネル(Wave 2、`hfauto/chemistry/`)

| モジュール | 公開関数(確定) | 仕様の要点 |
|---|---|---|
| `xyz.py`(追加) | `geometry_fingerprint(symbols, coords) -> str`、`Molecule(xyz: XYZ, charge: int, multiplicity: int)`、`Molecule.fingerprint()`(電荷と多重度も含む)、`Molecule.write(path)`、`composition_key(symbols, charge, multiplicity) -> str` | 既存の `XYZ` はそのまま使う。指紋は、座標を 1e-6 Å で丸めた正規化 JSON の sha256。組成キーは元素数(Hill 順)、電荷、多重度から作る |
| `vibrations.py` | `ISOTOPIC_MASSES`、`external_basis(symbols, coords) -> ndarray`、`projected_frequencies(hessian_eh_bohr2, symbols, coords_A) -> (freqs_cm1, modes, n_external)`、`rotational_constants_ghz(symbols, coords)`、`curvature_along(hessian, mode) -> float`、`frame_residual(hessian, symbols, coords) -> float`、`to_canonical_npy(hessian, path)` | ① 質量を加重する。② 並進 3 本と回転 3 本のベクトルを SVD にかけ、特異値が最大値の 1e-3 倍を超える本数から k(5 または 6)を決める(最適化で得た直線分子の 1e-4 Å 程度の数値ノイズは直線として扱い、179° の屈曲は非直線として扱う)。③ その直交補空間の基底 B で H_int = Bᵀ H_mw B を作り、(3N − k) 次元で対角化する。④ 虚数は負値で返す。⑤ モードは質量加重を外して規格化する。質量は主同位体の値(H 1.007825、C 12.000000、N 14.003074、O 15.994915、F 18.998403 …。H〜Kr と Br、I)を使う。「0 に近い 6 本を捨てる」方式は採らない(chem 6) |
| `topology.py` | `bonds(symbols, coords, previous=None)`、`fragments(...)`、`bond_changes(symbols, a, b)`、`labile_hydrogens(...)`、`acceptor_atoms(...)`、`vdw_radius(symbol)`、`proton_coordinate(coords, d, h, a)`、`state_label(symbols, coords)`、`transferred_hydrogens(symbols, a, b) -> int` | r ≤ 1.15 Σr_cov なら結合、r ≥ 1.45 Σr_cov なら非結合とする。中間は直前の状態を維持する(ヒステリシス)。labile H は N/O/F/S/Cl/Br/I/P に結合した H。q = r(D–H) − r(A–H) とし、プロトンの移動は \|q\| ≥ 0.3 Å で符号が反転した場合だけとする。対称な H 結合(FHF⁻)は 1 断片として扱う。`state_label` は、断片組成のソート列と、元素ラベル付きグラフの WL ハッシュ(3 反復)の先頭 8 桁から作る |
| `identity.py` | `permutation_invariant_rmsd(symbols, a, b) -> (rmsd, perm)`、`compare_minima(symbols, a, b, ea, eb) -> "same" \| "distinct" \| "ambiguous"`、`assign(symbols, coords, energy, candidates) -> str \| None`、`mapped_equivalent(symbols, a, b) -> bool`、`mapped_rmsd(a, b)`、`periodic_nearest(value_deg, targets_deg) -> int` | 同じ元素の中で `scipy.optimize.linear_sum_assignment` による割当てを行い、Kabsch と交互に最大 5 回反復する。Kabsch は **det = +1 の真性回転だけ**を許す(鏡像体を同一視しない)。初期配向は、慣性主軸の符号の 4 通りと、縮退した主軸の面内で 30° 刻みの回転を試し、最小値を採る。閾値: same は RMSD ≤ 0.02 Å かつ \|ΔE\| ≤ 1e-5 Eh。distinct は > 0.05 Å または > 5e-5 Eh。その間は ambiguous。`assign` は ≤ 0.05 Å かつ ≤ 5e-5 Eh で、次点が最良の 3 倍以上、または 0.1 Å 以上離れていること。`mapped_equivalent` は、恒等写像での Kabsch RMSD が 0.05 Å を超え、かつ置換不変 RMSD が 0.05 Å 以下のとき真 |
| `modes.py` | `displace(coords, mode, amplitude_A) -> (plus, minus)`、`overlap(a, b)`、`qrc_amplitude(hessian, mode, *, target_hartree, bounds_A=(0.05, 0.4)) -> float`、`classify_mode_follow(source, plus, minus) -> "replace" \| "ts_candidate" \| "one_side" \| "same_as_source"` | 変位は最大原子変位で規格化する。`qrc_amplitude`: モード u を最大原子変位 1 Å に規格化し、曲率 κ = uᵀHu(Eh/Å²)として s = √(2·target/\|κ\|) を bounds に収める(chem 7) |
| `profile.py` | `interior_maxima(energies, resolution)`、`shape(energies, resolution) -> "monotonic" \| "single_max" \| "multi_max"`、`hei(coords_list, energies) -> (index, energy, coords)`、`string_converged(gmax_history, gmax_tol=1e-3, non_worsening=3)`、`stagnated(gmax_history, xmax_history=(), window=10, min_improvement=0.05, xmax_zero_run=5)`、`max_node_spacing(coords_list)`、`tangent(coords_list, index)` | NWChem が出す "converged" の表示は使わない(CH-08) |
| `interpolation.py` | `align_mapped(a, b) -> b_aligned`、`idpp(symbols, a, b, n_images) -> list[ndarray]` | IDPP の前に、恒等写像での Kabsch 整列を必ず行う(CH-09)。Smidstrup 2014 の目的関数を勾配降下で最大 500 反復する。最小原子間距離が 0.7 Å 未満になれば ValueError |

---

## 6. バックエンド Protocol と実装

### 6.1 `hfauto/backends/protocols.py`(Wave 3)

```python
class Capability(StrEnum):
    QM = "qm"; PATH = "path"; SADDLE = "saddle"; CONFORMERS = "conformers"; DISCOVERY = "discovery"; THERMO = "thermo"

@dataclass(frozen=True)
class Requirements:
    executables: tuple[str, ...] = ()        # site ファイルの executables キー
    python_modules: tuple[str, ...] = ()     # worker 側で必要なモジュール
    version_command: tuple[str, ...] = ()    # preflight で版数を得るコマンド

@runtime_checkable
class Engine(Protocol):
    name: ClassVar[str]
    def requirements(self) -> Requirements: ...
    def supports(self, method: MethodSpec) -> bool: ...

@runtime_checkable
class QMEngine(Engine, Protocol):
    def energy(self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None) -> Evidence | Failure: ...
    def optimize(self, mol: Molecule, method: MethodSpec, *, tight: bool = False, init_hessian: Evidence | None = None,
                 deadline: Deadline | None = None) -> Evidence | Failure: ...
    def frequencies(self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None) -> Evidence | Failure: ...

@runtime_checkable
class PathEngine(Engine, Protocol):
    def find_path(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
                  initial_path: FileRef | None = None, refine_ts: bool = False,
                  deadline: Deadline | None = None) -> PathProfile | Failure: ...

@runtime_checkable
class SaddleRefiner(Engine, Protocol):
    def refine(self, seed: Molecule, method: MethodSpec, *, hessian: Evidence, mode_index: int | None = None,
               deadline: Deadline | None = None) -> Evidence | Failure: ...

@runtime_checkable
class ConformerEngine(Engine, Protocol):
    def search(self, mol: Molecule, method: MethodSpec, settings: ConformerSettings, *,
               deadline: Deadline | None = None) -> ConformerEnsemble | Failure: ...

@runtime_checkable
class DiscoveryEngine(Engine, Protocol):
    def explore(self, source: Molecule, trial: ReactionTrial, method: MethodSpec, settings: DiscoverySettings, *,
                deadline: Deadline | None = None) -> DiscoveryResult | Failure: ...

@runtime_checkable
class ThermoEngine(Engine, Protocol):
    def thermo(self, freq: Evidence, settings: Sequence[ThermoSettings], *,
               deadline: Deadline | None = None) -> list[ThermoResult] | Failure: ...   # 設定 × 温度ごとに 1 件
```

- `Molecule`(`chemistry.xyz`)は電荷と多重度を持つ。**すべてのアダプタは、これを必ずエンジンに渡す**(chem 8)。
  - NWChem: トップレベルの `charge` と `dft` の `mult`。多重度が 1 より大きければ `odft` も付ける。
  - xTB / CREST: `--chrg`、`--uhf`。
  - pysisyphus: `calc: {type: xtb, charge, mult}`。
  - ReaDuct: `molecular_charge`、`spin_multiplicity`、`spin_mode`。
- **初期 Hessian は型で受け渡す**(eng P0-5、chem 20)。
  - `init_hessian` と `hessian` には、task が `"freq"` の Evidence を渡す。その `final.fingerprint` は、対象構造の `geometry_fingerprint` と一致しなければならない。一致しない場合、アダプタは `INPUT_INVALID("hessian_geometry_mismatch")` を返す。
  - Level は問わない。したがって、xTB の Hessian を NWChem の初期 Hessian に使ってよい。初期 Hessian は推定値にすぎず、最終の判定は必ず別ジョブの DFT freq で行うためである。
  - アダプタは、正準形式の `.npy` からエンジン固有の形式(NWChem なら `.hess` の下三角、a.u.)を書き出す。
- 結果型は同じファイルに置く。frozen な pydantic モデルで、JobStore から復元するための `kind` を持つ。
  - `ConformerSettings(nci: bool, quick: bool = True, threads: int = 4, ewin_kcal: float = 6.0, topology: Literal["on", "off", "noref"] = "on", notopo_atoms: tuple[int, ...] = ())`
  - `ConformerEnsemble(kind="conformers", members: tuple[tuple[Geometry, float | None], ...], topology_removed: int, topology_stops: tuple[Geometry, ...], version: str, job_key: str)`
  - `DiscoverySettings(max_scf_iterations: int = 300, electronic_temperature_K: float = 300.0, scc_retry_temperature_K: float = 1000.0, afir_gamma_kj_mol: float = 125.0, afir_gamma_retry_kj_mol: float = 300.0, nt_total_force_norm: float = 0.1, imag_cutoff_cm1: float = 50.0, timeout_s: float = 600.0)`
  - `DiscoveryResult(kind="discovery_result", outcome: Literal["product", "negative"], reason: str | None, product: Geometry | None, ts: Geometry | None, ts_imag_cm1: float | None, barrier_kj_mol: float | None, reaction_kj_mol: float | None, irc_connected_to_source: bool, electronic_temperature_K: float, job_key: str)`
  - `ThermoResult(kind="thermo_result", settings_sha: str, T_K: float, G_hartree: float, H_hartree: float, E_hartree: float, zpe_hartree: float, n_real: int, S_rot: float, notes: tuple[str, ...], job_key: str)`
- エンジンは `Engine(jobs: JobRunner, site: EngineSite)` で生成する。

### 6.2 型付き遅延 registry `hfauto/backends/engines.py`(最終形の表を Wave 3 で先に記入する)

```python
_TABLE: dict[Capability, dict[str, str]] = {
    Capability.QM:         {"nwchem": "hfauto.backends.nwchem.engine:NWChemEngine",
                            "xtb":    "hfauto.backends.xtb:XTBEngine"},
    Capability.PATH:       {"nwchem_string": "hfauto.backends.nwchem.engine:NWChemString",
                            "pysis_gs":      "hfauto.backends.pysis.engine:PysisGrowingString"},
    Capability.SADDLE:     {"nwchem_saddle": "hfauto.backends.nwchem.engine:NWChemSaddle"},
    Capability.CONFORMERS: {"crest": "hfauto.backends.crest:CRESTEngine"},
    Capability.DISCOVERY:  {"readuct": "hfauto.backends.readuct.engine:ReaDuctEngine"},
    Capability.THERMO:     {"goodvibes": "hfauto.backends.goodvibes.engine:GoodVibesEngine"},
}
def create(capability: Capability, name: str, *, jobs: JobRunner, site: EngineSite) -> Engine: ...
def requirements_for(selection: Mapping[Capability, Sequence[str]]) -> dict[str, Requirements]: ...
@contextmanager
def override(capability: Capability, name: str, factory: Callable[..., Engine]) -> Iterator[None]: ...  # テスト専用
```

- エントリは 8 つだけである。`dummy`、`internal`、`quasi_rrho` などの別名は置かない。表にないキーを指定すると `KeyError` になる。
- pysisyphus の RS-I-RFO は SADDLE として登録しない。低レベルの TS は、`pysis_gs` を `refine_ts=True` で呼んで得る。このとき 1 本の pysisyphus 入力(`cos` + `tsopt`)の中で、`run_tsopt_from_cos` が HEI の接線との重なりから root を選ぶ。そのため、root の選択を driver 側で作り直す必要はない(chem 1)。
- この registry を import してよいのは `hfauto.pipeline` だけである(§11.1 の契約 2)。

### 6.3 具象実装(Wave 4)

| エンジン | 実装する能力 | 要点(組み込む修正) |
|---|---|---|
| `nwchem`(`backends/nwchem/`) | QM(energy / optimize / frequencies)、PATH(`nwchem_string`)、SADDLE(`nwchem_saddle`)、WFT の energy | 構造指定・電子状態<br>・geometry は常に `units angstrom nocenter noautosym` とする。<br>・`coordinates: auto` のときは autoz で始め、`AUTOZ failed` が出たら ladder が `cartesian`(`noautoz`)で再投入する(CH-14)。<br>・トップレベルに `charge`、`dft` ブロックに `mult` を書く。多重度が 1 より大きければ `odft` を付ける。基底は spherical で指定する。<br>・分散は、d3bj なら `disp vdw 4`、d3zero なら `disp vdw 3` とする。`solvation: cosmo:<eps>` なら COSMO ブロックを書く(CH-16)。<br>・driver は既定の閾値で maxiter 100 とする(CH-12)。<br><br>ジョブの分割と Hessian<br>・**opt と saddle のジョブには hessian / frequencies タスクを入れない。freq は別ジョブ(`task dft frequencies`)で計算する**(CH-11、CH-38)。<br>・optimize の `init_hessian` と saddle の `hessian` は、正準形式の `.npy` から `<name>.hess` を permanent_dir に書き出し、`inhess 2` で読ませる。<br><br>saddle<br>・`trust 0.1`、`sadstp 0.1`、maxiter 50 とする。<br>・`mode_index` が与えられ、それが最低のモードでない場合は、`noautoz` と `moddir <mode_index+1>` を指定する(chem 17)。<br><br>string(ZTS)<br>・`nbeads`、`maxiter 20`、`stepsize 0.05`、`interpol 3`、`impose`、`freeze1`、`freezeN`、`xyz_path` を指定し、**`tol 1e-5` を明示する**(chem 18)。<br>・`string: gmax,grms,xrms,xmax=` の行をパースして `gmax_history` に入れる(CH-08)。<br><br>WFT<br>・SP だけを扱う。`task mp2 energy` または `task ccsd(t) energy` とし、`freeze atomic` を付ける。<br><br>出力の検査<br>・出力から `Level` と `<S2>` を読み取る。<br>・出力に出た初期構造を入力と照合し、1e-4 Å を超えてずれていれば `METHOD_MISMATCH("frame_changed")` とする。<br>・`level_mismatches` が空でなければ `METHOD_MISMATCH` とする。<br>・振動数は、`.hess` を正準形式の `.npy` に変換し、`projected_frequencies` で求める。<br><br>継続実行と実行環境<br>・継続実行では、最新の `final-NNN.xyz`、`.movecs`、`.drv.hess` を新しい attempt にコピーする(AR-08)。<br>・scratch は site で指定した ext4 上のパスに置き、`OMP_NUM_THREADS=1` とする(CH-13)。<br><br>削除するもの: HF 記述子、スキーマ移行、`_fallback`、backend 内の geometry_qc |
| `xtb`(`backends/xtb.py`) | QM | ・最適化は `--opt vtight`(`tight=True` でも同じ)、Hessian は `--hess` で求める。<br>・`--gfn`、`--chrg`、`--uhf`、`--alpb` を渡す。`electronic_temperature_K` が与えられた場合は `--etemp` も渡す。<br>・収束の条件は、rc=0 であること、出力に `FAILED TO CONVERGE` がないこと、`NOT_CONVERGED` ファイルがないことの 3 つ(BUG-07、CH-20)。<br>・`hessian` ファイルを正準形式の `.npy` に変換し、`projected_frequencies` にかける。単位は G17 と smoke テストで確認する。<br>・`init_hessian` は受け取っても使わない(xTB では不要)。<br>・環境変数は `OMP_NUM_THREADS=n,1`、`OMP_STACKSIZE=4G` とする。<br>・Dummy クラスからの継承は廃止する(AR-16) |
| `pysis`(`backends/pysis/`) | PATH(`pysis_gs`) | ・**xTB ネイティブ計算器でだけ使い、DFT では使わない**。QCEngine は依存ごと削除する(CH-01、CH-02)。<br>・worker サブプロセスの中で `pysisyphus.run.run_from_dict` を呼ぶ。<br>・`calc: {type: xtb, gfn, charge, mult, pal: threads}` とする。worker の環境変数 `PATH` に、site の xtb ディレクトリを入れる。<br>・GS は `max_nodes = images − 2`、`climb true`、`climb_rms 5e-3` とする。<br>・**GS の座標系は、直線分子か 2 断片以上なら `cart`、それ以外は `dlc` とする。`tric` は GS では使えない**(chem 1)。<br>・GS の `interpol` は既定値のままにする。pysisyphus 1.0 が指定を要求する場合だけ `idpp` を指定する(WP4.2 がソースで確認する)。<br>・`refine_ts=True` のときは、同じ入力に `tsopt: {type: rsirfo, thresh: gau}` を加える(2 断片以上なら coord_type は `tric`)。得られた TS を `PathProfile.ts` に入れる |
| `crest`(`backends/crest.py`) | CONFORMERS | ・設定から CLI への対応表: `--gfn2 [--nci] [--quick] -T <threads> --ewin <ewin> --chrg <q> --uhf <m-1> [--notopo <atoms>] [--noreftopo] [--alpb <solvent>] --scratch <ext4>`。<br>・`crest --version` の出力を記録する。<br>・エネルギーが欠けている構造は None のまま扱う(BUG-08)。<br>・CREGEN のトポロジー除去数と、トポロジー停止した構造を結果に含める。<br>・RDKit へのフォールバックは削除する |
| `readuct`(`backends/readuct/`) | DISCOVERY | 実行<br>・1 回の attempt を worker サブプロセスで実行する。cwd は attempt ディレクトリ、timeout は 600 s とする(CH-31)。<br><br>計算器の設定<br>・`molecular_charge` と `spin_multiplicity` を渡す。`spin_mode` は、多重度が 1 なら restricted、それ以外は unrestricted とする。<br>・`max_scf_iterations 300`、`electronic_temperature 300` とする。<br>・SCC が失敗したら、1000 K で 1 回だけ再試行する。その場合は TS と生成物の構造で **300 K の xTB SP を取り直し**、エネルギー窓をその値で評価する(chem 11)。<br><br>TS の最適化と検証<br>・TS 最適化には `automatic_mode_selection = sorted(associations ∪ dissociations)` を渡す(CH-27)。<br>・虚振動は `projected_frequencies` で求め、ν < −50 のものを数える(CH-26)。<br>・IRC の片端が出発構造と一致すること(グラフ同型かつ置換不変 RMSD < 0.1 Å)と、両端の低レベルでの n_imag が 0 であることを確認する。<br><br>AFIR と例外処理<br>・AFIR の γ は 125 kJ/mol から始め、生成物が出なければ 300 で 1 回だけ試す。<br>・`except` で捕まえるのは SCINE の実行時例外だけに限る |
| `goodvibes`(`backends/goodvibes/`) | THERMO | **GoodVibes 4.3.0 の Python API を worker から呼ぶ**(chem 5)。<br><br>呼び出し<br>・worker は Evidence から `goodvibes.io.QCData` を組み立て、`goodvibes.api.compute_thermo(qcdata=..., QS, s_freq_cutoff, temperature, freq_scale_factor, zpe_scale_factor, invert, symm)` を、設定と温度の組ごとに呼ぶ。<br>・QCData に入れるもの: `scf_energy`、`charge`、`multiplicity`、`atom_types`、`cartesians`、`frequency_wn` / `im_frequency_wn`(**自前で射影した振動数**)、`linear_mol`(`n_external == 5` のとき)、`molecular_mass`(同位体質量)、回転定数(`rotational_constants_ghz`)。<br>・`calc_bbe` が実際に使うフィールドを 4.3.0 のソースで確認し、必要なものだけを埋める。<br>・スケール因子は必ず明示して渡す(自動検索させない)。<br><br>検査とキャッシュ<br>・S_rot > 0 であることを確認する。<br>・キャッシュキーは、freq Evidence の job_key と hessian の sha、設定一式の組とする(CH-42、CH-45)。<br><br>持たないもの<br>・テキストのパース、CSV の読み取り、`--help` の解析、内部フォールバック(CH-38、CH-39)。<br><br>W4 の smoke テスト<br>・同じ freq 出力について、CLI で計算した G と 1e-5 Eh 以内で一致することを 1 回だけ確認する |

削除する実装:
- ORCA(`qm/orca.py`、`ts/orca_nebts.py`)
- NWChem NEB
- pysisyphus+QCEngine(IRC、dimer、RS-I-RFO)
- 線形補間ブラケット
- dummy 3 種(`tests/fakes.py` へ移す)
- RDKit の配座エンジン
- DB プロバイダ、kinetics、arkane、simple

### 6.4 テスト用の偽物 `tests/fakes.py`(Wave 3 で API を凍結する)

- `FakeQM(root: Path, pes)` は、解析的な PES から energy、gradient、hessian を計算する。
  - optimize は BFGS で gmax < 1e-5 まで進め、`trajectory_energies_hartree` を記録する。
  - frequencies は `projected_frequencies` で求め、正準形式の `.npy` を書き出す。
  - 観測値の Level は `program="fake"` とする。
  - PES ライブラリには次の 5 種を用意する。
    - `harmonic`
    - `double_well`: 3 原子で、障壁の高さを指定できる
    - `symmetric_double_well`: 縮退あり
    - `triple_well`: 中間体あり
    - `flat_uphill`: 障壁なし
- ほかに次の偽物を置く。
  - `FakePath`: `refine_ts` に対応する
  - `FakeSaddle`: PES 上で固有ベクトルを追跡する
  - `FakeConformers`、`FakeDiscovery`、`FakeThermo`: 台本どおりの結果を返す
- 登録は `engines.override(...)` の fixture だけで行う。本番の registry に `fake` キーは存在しない。
- W4 の各 WP が必要とするドメイン固有の偽の挙動(saddle の崩壊など)は、その WP のテストファイル内でサブクラスとして定義する。`tests/fakes.py` は W3 の後は変更しない。

---

## 7. オーケストレーション

### 7.1 実行層 `hfauto/execution/`(Wave 3)

```python
# process.py(subprocess を import する唯一のモジュール)
@dataclass(frozen=True)
class Command:
    argv: tuple[str, ...]
    cwd: Path
    env: Mapping[str, str]
    stdin: Path | None = None

@dataclass(frozen=True)
class CommandResult:
    returncode: int | None
    timed_out: bool
    stopped: str | None          # monitor が返した停止理由
    duration_s: float
    stdout: Path
    stderr: Path

def run_command(cmd: Command, *, timeout_s: float, monitor: Callable[[Path], str | None] | None = None,
                poll_s: float = 10.0) -> CommandResult     # command_result.json を書く。プロセスツリーを停止できる
def resolve_executable(name: str, explicit: str | None) -> str | None

# jobs.py
@dataclass(frozen=True)
class Task:
    engine: str
    version_pin: str               # EngineSite.version(実行前に分かる。eng 7)
    kind: str                      # "optimize", "frequencies", "string", ...
    key_payload: Mapping[str, Any] # method.signature()、分子の fingerprint、パラメータ、入力 FileRef の sha
    execution: ExecutionSpec
    inputs: Mapping[str, Any] = field(default_factory=dict)   # render 用(キーには入れない)

class Adapter(Protocol[T]):
    result_type: type[T]           # pydantic モデル。JobStore からの復元に使う(文字列からの動的 import はしない)
    def prepare(self, task: Task, workdir: Path) -> Command: ...
    def parse(self, task: Task, workdir: Path, result: CommandResult) -> T | Failure: ...
    def monitor(self, task: Task) -> Callable[[Path], str | None] | None: ...
    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None: ...

class JobRunner:
    def __init__(self, store: JobStore, *, cores: int) -> None: ...
    def run(self, task: Task, adapter: Adapter[T], *, deadline: Deadline | None = None) -> T | Failure: ...
    def stats(self) -> JobStats: ...          # hits, misses, failures_by_kind

def thread_map(fn: Callable[[I], R], items: Sequence[I], *, workers: int) -> list[R]: ...   # 入力順で返す

# lock.py
class SiteLock:            # <scratch_root>/.hfauto_site.lock を O_EXCL で作る(pid、run_dir、時刻を書く)。pid が死んでいれば奪う
    def __enter__(self): ...   # 別の run が保持していれば、保持者を示して直ちに RuntimeError
```

- **JobStore**
  - 置き場所は常に `<run>/jobs/` とし、run をまたいで共有するキャッシュは作らない(eng 6)。
  - キーは、正規化 JSON `{engine, version_pin, kind, key_payload}` の sha256 とする。`ExecutionSpec`(ranks、memory、timeout、実行ファイルのパス、scratch)はキーに**含めない**。
  - レイアウトは `jobs/<k[:2]>/<key>/{job.json, result.json, attempt_NN/}` とする。
  - `result.json` の形式は `{"kind": ..., "data": model_dump(mode="json"), "files": {relpath: sha256}}` とする。再利用するときは、ファイルの sha を照合してから `adapter.result_type.model_validate` で復元する。
  - 終端的な失敗(`INPUT_INVALID`、`METHOD_MISMATCH`、`EXECUTABLE_MISSING`、`INCOMPLETE_OUTPUT`)も保存し、同じ要求は再実行しない。再実行するのは `--retry-failed <kinds>` で指定された種類だけとする。
  - ジョブ単位の排他は、`.lock` を O_EXCL で作って行う。
- **版数の照合**: アダプタは、出力から観測した `Level.version` を `Task.version_pin` と比べ、違えば `METHOD_MISMATCH` を返す。
- **attempt ladder**: `jobs.py` の定数表で定義する。これ以外の自動再試行はしない。

| FailureKind | 対処 | 上限 |
|---|---|---:|
| TIMEOUT、GEOMETRY_MAXITER | `adapter.continuation` で継続する(NWChem は最新構造・movecs・drv.hess、xTB は最終構造。pysis と ReaDuct は None) | 2 |
| INPUT_INVALID(`autoz`) | `coordinates="cartesian"` で再投入する | 1 |
| SCF_NOT_CONVERGED | 前回の vectors を使い、damping と level shift をかける | 1 |
| その他 | 直ちに失敗を返す | 0 |

- attempt の timeout は min(`ExecutionSpec.timeout_s`, `deadline.remaining()`) とする。残り時間が 60 s 未満なら実行せず、`BUDGET_EXHAUSTED` を返す(CH-37)。
- **並行実行**(eng P0-2)
  - `JobRunner` はコア数のセマフォを持ち、各ジョブは ranks × threads 分を確保する。同時に使うコアの合計は site.cores を超えない。
  - DFT と反応ケースは直列に実行する。
  - xTB と CREST だけは、stage が `thread_map(workers = site.cores // threads_per_job)` で並列に実行する。
  - `SiteLock` は `hfauto run` の実行中ずっと保持し、別の run の同時起動を拒否する。CH-13 で報告された 17 倍の遅延を防ぐためである。
- **worker**: `python -m hfauto.execution.worker <module:function> <job.json>` の形で起動する。SCINE、pysisyphus、GoodVibes を import するのは、この子プロセスの中だけとする。cwd は attempt ディレクトリにする。

### 7.2 MinimumDriver `hfauto/drivers/minimum.py`(Wave 3)

```python
@dataclass(frozen=True)
class MinimumPolicy:
    max_mode_follow: int = 2
    amplitude_A: float = 0.1
    gates: Policy = Policy()

@dataclass(frozen=True)
class MinimumOutcome:
    status: Literal["minimum", "soft_minimum", "saddle", "known", "failed"]
    opt: Evidence | None
    freq: Evidence | None
    history: tuple[str, ...]
    known_basin: str | None = None          # status == "known" のとき既存 basin の id
    failure: Failure | None = None
    ts_candidate: tuple[Evidence, Evidence] | None = None   # ± 変位が別々の極小に落ちた場合(追加計算なしで得られる TS 候補)

def relax_to_minimum(mol: Molecule, method: MethodSpec, qm: QMEngine, *, known: Registry | None = None,
                     init_hessian: Evidence | None = None, policy: MinimumPolicy = MinimumPolicy(),
                     deadline: Deadline | None = None) -> MinimumOutcome

class Registry:
    """組成 × level_key ごとの極小の集合。identity.compare_minima と assign で判定する。"""
    def __init__(self, minima: Iterable[MinimumRecord], load_xyz: Callable[[Geometry], XYZ]) -> None: ...
    def find(self, symbols, coords, energy) -> str | None: ...          # identity.assign
    def add(self, outcome: MinimumOutcome, species: SpeciesRecord, *, tier: Literal["screen", "dft"]
            ) -> tuple[MinimumRecord, Literal["new", "same", "ambiguous"]]: ...
    def members(self, basin_id: str) -> tuple[str, ...]: ...
```

アルゴリズム:
1. `opt` を実行する。`init_hessian` が与えられていれば、それを初期 Hessian に使う。
2. `known` が与えられ、opt の最終構造が既存の極小と `same` なら、freq を省いて `status="known"` とする。
3. `freq` を計算し、`gates.is_minimum(freq, opt=opt)` と `imaginary_tier` で判定する。
4. tier が `saddle` の場合: 最低の虚モードに沿って ±0.1 Å 変位し、両方を opt → freq にかける。結果は `modes.classify_mode_follow` で次のように分類する。これを最大 2 サイクル行う。
   - 両側が同じ極小 X に落ちた: 元の構造を X で置き換える。
   - 別々の極小 A と B に落ちた: `ts_candidate` とする。
   - 片側だけが収束した: その側を採る。
5. tier が `soft` の場合: 1 回だけ変位し、それでも残れば `soft_minimum` とする。
6. tier が `noise` の場合: そのまま `minimum` とする。
7. `Registry.add` の結果が ambiguous のとき、呼び出し側(minima stage)は `tight=True` で再最適化し、freq を取り直して、1 回だけ再判定する(CH-17、4.5.7)。
8. Registry への追加は、stage が結果をすべて集めてから、species_id の順に決定的に行う。

### 7.3 ReactionCaseDriver `hfauto/drivers/reaction_case/`(state.py は Wave 3、actions.py と driver.py は Wave 4)

```python
class Action(StrEnum):
    SCREEN = "screen"                    # 低レベル経路と DFT SP による障壁の事前判定
    REFINE_SADDLE = "refine_saddle"      # 初期 Hessian(1 回目は xTB、2 回目は DFT)→ NWChem saddle
    VALIDATE_TS = "validate_ts"          # 別ジョブの DFT freq → is_first_order_saddle
    FIND_PATH = "find_path"              # DFT string(種が得られないとき、または saddle が失敗したときだけ)
    CONNECT = "connect"                  # QRC: エネルギー目標で変位 → opt → 一意割付け → connection
    VALIDATE_INTERMEDIATE = "validate_intermediate"   # relax_to_minimum → Registry
    COMPLETE = "complete"
    BLOCKED = "blocked"

@dataclass(frozen=True)
class CasePolicy:                        # 既定値はここ 1 か所だけに置く
    screen: bool = True
    max_saddle_attempts: int = 2
    string_beads: int = 9
    string_chunks: int = 3
    confirm_barrierless_beads: int = 15  # 多解像度での確認(旧 path-ensemble)
    max_path_runs: int = 2               # DFT string の本計算 1 回と確認 1 回
    screen_images: int = 11              # GS のノード数(両端を含む)
    qrc_target_hartree: float = 3.0e-4   # 変位で下げたいエネルギー(max(3 × drop, この値))
    qrc_bounds_A: tuple[float, float] = (0.05, 0.4)
    qrc_retry_factor: float = 2.0
    walltime_s: float = 6 * 3600
    max_split_depth: int = 2
    override_negative_evidence: bool = False
    gates: Policy = Policy()

@dataclass(frozen=True)
class Seed:
    geometry: Geometry
    source: Literal["discovery_ts", "screen_ts", "screen_hei", "path_hei", "higher_order_retry"]
    tangent: tuple[float, ...] | None    # 経路の接線、または端点の差ベクトル(写像付きで整列したもの)

@dataclass(frozen=True)
class CaseState:                         # driver がメモリ上で積み上げるスナップショット(ログから再生はしない)
    expired: bool = False
    screen: BarrierVerdict | None = None
    seeds: tuple[Seed, ...] = ()         # 未使用の種(先頭から使う)
    saddle_attempts: int = 0
    last_saddle: Literal[None, "converged", "failed"] = None
    ts_check: Literal[None, "ok", "collapsed", "higher_order"] = None
    claim: SaddleClaim | None = None
    connection: Literal[None, "elementary", "degenerate", "reassigned", "failed"] = None
    connection_attempts: int = 0
    path_runs: tuple[Literal["single_max", "multi_max", "monotonic", "failed"], ...] = ()
    intermediate: Literal[None, "distinct", "same_as_endpoint"] = None

@dataclass(frozen=True)
class Decision:
    action: Action
    reason: str
    outcome: CaseOutcome | None = None

def decide(case: ReactionRecord, state: CaseState, policy: CasePolicy) -> Decision     # 純関数(IO なし)
def drive_case(case: ReactionRecord, rt: CaseRuntime, policy: CasePolicy) -> CaseResult
```

**判断表**(全 17 行。上の行から順に評価する)

| # | 条件 | Decision |
|---:|---|---|
| 1 | 両端が、同じ `level_key` をもつ dft 極小になっていない | BLOCKED(`endpoints_not_on_one_pes`) |
| 2 | 両端が同じ basin にあり、`degenerate=False` | COMPLETE(SAME_BASIN) |
| 3 | ΔE_rxn > `reaction_window_kcal`(40) | COMPLETE(OUT_OF_WINDOW) |
| 4 | `state.expired` | COMPLETE(UNRESOLVED, `walltime`) |
| 5 | `connection` が判定済み | elementary → ELEMENTARY_STEP、degenerate → DEGENERATE、reassigned → REASSIGNED で COMPLETE。failed で `connection_attempts < 2` なら CONNECT(振幅 × 2)。それ以外は COMPLETE(UNRESOLVED, `connection_failed`) |
| 6 | `claim` がある(TS 検証済み) | CONNECT |
| 7 | `last_saddle == "converged"` かつ `ts_check is None` | VALIDATE_TS |
| 8 | `ts_check == "collapsed"` かつ `intermediate is None` | VALIDATE_INTERMEDIATE(崩壊した構造) |
| 9 | `intermediate == "distinct"` | COMPLETE(MULTI_STEP)。子反応 R→I と I→P を生成する(`source="split"`、化学量論は親から継承、深さは `max_split_depth` まで) |
| 10 | `ts_check == "higher_order"` かつ `saddle_attempts < max` | REFINE_SADDLE(種: 2 本目の虚モードに沿って 0.1 Å 下った構造、`source="higher_order_retry"`) |
| 11 | `policy.screen` が有効で `screen is None` | SCREEN |
| 12 | `screen.verdict` が `"barrierless"` なら COMPLETE(BARRIERLESS)、`"negative_evidence"` なら COMPLETE(NO_PRODUCT) | — |
| 13 | 直近の DFT 経路が `multi_max` で、`intermediate is None` | VALIDATE_INTERMEDIATE(極大に挟まれた最低ノード) |
| 14 | 直近の DFT 経路が `monotonic` | 直近 2 回とも monotonic なら確認済みとして COMPLETE(BARRIERLESS, `dft_path_monotonic`)。そうでなく `len(path_runs) < max_path_runs` なら FIND_PATH(beads 15、多解像度での確認) |
| 15 | `seeds` が空でなく、`saddle_attempts < max_saddle_attempts` | REFINE_SADDLE(先頭の種) |
| 16 | `len(path_runs) == 0`(DFT string が未実行。screen が無効または unavailable の場合、種を使い切った場合、saddle が失敗した場合を含む) | FIND_PATH(初期経路の優先順: 低レベル GS の経路を再標本化したもの → SCREEN が作った IDPP → 新たに作る IDPP) |
| 17 | 上のどれにも当たらない | COMPLETE(UNRESOLVED, `attempts_exhausted`) |

- DFT 経路の HEI を種に加えるのは、その経路の形が `single_max` のときだけである(FIND_PATH の action が seeds に追加する)。このため、multi_max や monotonic の経路から REFINE_SADDLE が選ばれることはない(eng P0-4)。
- 行 7〜10 で扱った種は消費済みとして扱う。`ts_check == "ok"` のときは、VALIDATE_TS の action が `claim` を設定し、次の評価で行 6 に進む。
- REFINE_SADDLE の action は、開始時に `last_saddle`、`ts_check`、`intermediate` を None に戻す。種ごとに行 7〜10 を評価し直すためである。
- 経路が multi_max で、VALIDATE_INTERMEDIATE の結果が `same_as_endpoint` だった場合は、その action が最も高い極大のノードを `path_hei` の種として追加する。
- `CaseRuntime`(`driver.py` の frozen dataclass)は、`qm`、`saddle`、`path`、`screen_qm`、`screen_path`、`method`、`screen_method`、`registry`、`load_xyz`、`file_ref`、`case_dir`、`deadline` を持つ。stage が `StageRuntime` から組み立てる。driver が import してよいバックエンドは `hfauto.backends.protocols` だけである。

**action の実装**(`actions.py`)
- 1 つの action は 60 行以内とする。すべて JobRunner を経由するので、どの action も冪等である。

- **SCREEN**(chem 16)
  1. **低レベル TS がある場合の近道**: `case.low_level_ts` があれば、次を行う。
     - xTB の freq で、最低モードが ν < −50 であることを確認する。
     - その構造と DFT の両端で DFT SP を計算し、3 点で `barrier_verdict` を判定する。
     - `proceed` なら、この構造を種(`discovery_ts`)にして終える。`proceed` でなければ 2 へ進む。
  2. **xTB による端点の再最適化**: DFT の両端を xTB で再最適化する。恒等写像での RMSD が 0.05 Å 以下(xTB で同じ構造に崩壊した)なら、次を行う。
     - GS は行わない。
     - DFT の端点を `align_mapped` で整列してから `idpp` で 11 点の経路を作り、全ノードで DFT SP を計算する。
     - `barrier_verdict(low_profile=None)` で判定する。
     - `proceed` のときは種を作らない。この IDPP 経路を、行 16 の FIND_PATH の初期経路として使う。
  3. **GS による経路**: 崩壊しなかった場合は `pysis_gs(refine_ts=True, images=11)` を実行し、次を行う。
     - TS が得られたら、xTB の freq で確認する(最低モードが ν < −50 で、ほかに ν < −50 のモードがないこと)。
     - 全ノード(11 点)と TS で DFT SP を計算し、`barrier_verdict` で判定する。
     - `tangent_mode_cm1` には、高い方の端点の DFT freq のうち、経路の接線との重なりが最大のモードの振動数を使う。
     - 種は、xTB の TS があればそれ(`screen_ts`)、なければ DFT の HEI ノード(`screen_hei`)とする。
  4. **失敗時**: GS が失敗したら `unavailable` とし、理由を記録する。

- **REFINE_SADDLE**(chem 20)
  - **1 回目の試行**: 種の構造で xTB の `frequencies` を計算する。負の固有値がちょうど 1 本で、接線との重なりが 0.3 以上なら、これを `hessian` として使う。
  - **それ以外と 2 回目の試行**: 種の構造で DFT の `frequencies` を計算して使う。同じ構造・同じ LOT の Hessian は JobStore で再利用する。
  - `mode_index` には、接線(または宣言された座標)との重なりが最大の負モードを選ぶ。重なりが 0.3 未満なら、合否には使わず注記だけ付ける。
  - 最後に `saddle.refine(seed, hessian=h, mode_index)` を呼ぶ。

- **VALIDATE_TS**
  - `qm.frequencies(saddle.final)` を計算し、`is_first_order_saddle(freq, saddle=..., endpoint_energies, torsional=case.torsional)` と `spin_ok` で判定する。
  - 合格なら `SaddleClaim` を作る。
  - n_imag = 0 なら `collapsed`、2 本目の虚振動が残れば `higher_order` を記録する。

- **CONNECT**
  - 変位の振幅は `modes.qrc_amplitude(TS の Hessian, 虚モード, target=max(3·drop, 3e-4))` で決める。2 回目は 2 倍にする。
  - ± 両方向に変位し、それぞれ `qm.optimize` で最適化してから `Registry.find` で既知の極小に割り付ける。
  - 未知の構造なら `relax_to_minimum` → `Registry.add` で新しい basin として登録する。
  - `gates.connection(degenerate=case.degenerate, sides_distinct=identity.mapped_equivalent(両側の最終構造))` で判定し、`ConnectionClaim` を作る。
  - ねじれ座標で割り付けるときは `periodic_nearest` を使う(CH-04)。

- **FIND_PATH**
  - `path.find_path(endpoints, images=beads)` を、最大 `string_chunks` 回に分けて実行する。各チャンクの `*.string_final.xyz` を次のチャンクの `initial_path` にする。
  - 採否は `profile.string_converged` と `stagnated` で判定し、経路の `shape` を記録する。
  - `single_max` なら、HEI を種(`path_hei`)に追加する。

- **VALIDATE_INTERMEDIATE**
  - `relax_to_minimum` → `Registry.add` を行う。両端のどちらとも別の basin であれば `distinct` とする。

- **再開**
  - 再実行すると、`decide` は決定論的に同じ判断をたどり、完了済みのジョブは JobStore から即座に返る。そのため、中断した地点まで数秒で再生される。
  - `cases/<reaction_id>/log.jsonl` は監査と `hfauto case` 用に書き出すだけで、読み戻さない(eng 24)。

- **出力**
  - `reaction`(子反応を含む)
  - `calculation`: saddle、TS の freq、QRC 両側の opt、SCREEN の SP
  - 新しい basin の `species` と `minimum`

### 7.4 stage の型(`hfauto/stages/spec.py`)と pipeline 層(`hfauto/pipeline/`)(Wave 3)

```python
# stages/spec.py(stage が依存してよい実行時インタフェースはこれだけ。pipeline がこれを実装する)
class StageConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

@dataclass(frozen=True)
class StageSpec:
    name: str
    config: type[StageConfig]
    consumes: tuple[ArtifactType, ...]
    produces: tuple[ArtifactType, ...]

class StageRuntime(Protocol):
    run_id: str
    stage_id: str
    stage_dir: Path
    system: SystemConfig
    conditions: Conditions
    policy: Policy
    def engine(self, capability: Capability, name: str) -> Engine: ...
    def method(self, method_id: str) -> MethodSpec: ...
    def load_xyz(self, geometry: Geometry) -> XYZ: ...
    def file_ref(self, path: Path) -> FileRef: ...                 # run ディレクトリからの相対パスと sha
    def thread_map(self, fn, items, *, threads_per_item: int) -> list: ...
    def deadline(self, seconds: float) -> Deadline: ...

class Stage(Protocol):
    spec: ClassVar[StageSpec]
    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]: ...

# stages/catalog.py
STAGES: dict[str, str] = {
    "structures": "hfauto.stages.structures:StructuresStage",
    "conformers": "hfauto.stages.conformer_search:ConformersStage",
    "minima": "hfauto.stages.minima:MinimaStage",
    "explore": "hfauto.stages.explore:ExploreStage",
    "reaction-paths": "hfauto.stages.reaction_paths:ReactionPathsStage",
    "sp": "hfauto.stages.single_point:SinglePointStage",
    "thermo": "hfauto.stages.thermochemistry:ThermoStage",
    "report": "hfauto.stages.report:ReportStage",
}
def get(name: str) -> type[Stage]: ...
@contextmanager
def override(name: str, cls: type[Stage]) -> Iterator[None]: ...   # テスト専用
```

- **`pipeline/config.py`**
  - 定義するモデル:
    - `SiteConfig(site, scratch_root, cores, memory_mb, engines: dict[str, EngineSite])`
    - `StageEntry(id, stage, extra fields)`
    - `PipelineConfig(pipeline_id, conditions: Conditions, gates: dict[str, float], stages)`
    - `ResolvedConfig(pipeline, system, site, methods: dict[str, MethodSpec], code_version)`
  - 読み込みは `load(pipeline_path, system_path, site_path) -> ResolvedConfig` で行う。method は `pipeline_path.parent.parent / "methods"` から id で読む。
  - 4 つの層は中身が重ならないので、マージや優先順位は設けない(eng 24)。stage 固有の追加項目は、その stage の config モデル(`extra="forbid"`)で検証する。
- **`pipeline/layout.py`**: `RunLayout(run_dir)` を定義する。
  - ディレクトリの配置: `<run>/<stage_id>/manifest.json`、`<run>/<stage_id>/diagnostics.json`、`<run>/<stage_id>/cases/`、`<run>/jobs/`、`<run>/run_state.json`、`<run>/resolved_config.yaml`。
  - `run_state.json` は、stage を**実行した順**に並べたリストである。各要素は {stage_id, pipeline_id, status: pending/running/done/failed/stale, input_sha, config_sha, started, finished, n_ok, n_failed, jobs: {hits, misses, failures_by_kind}}。
  - **`view(stage_id)` は、run_state でその stage より前にあり status が done の stage の manifest をすべて集め、実行順に `Manifest.union` したものとする**(eng 8)。
  - stage_id は 1 つの run の中で一意とする。別の pipeline_id で同じ id を使おうとした場合はエラーにする。
  - 同じ pipeline を再実行したとき、done で input_sha と config_sha が変わっていない stage は飛ばす(再開)。
  - パスはすべて run_dir からの相対パスで扱う。
- **`pipeline/runner.py`**
  - `Runtime` は `StageRuntime` の実装で、エンジンを (capability, name) ごとにキャッシュする。
  - `build_runtime(resolved, layout) -> Runtime`。
  - `execute_stage(entry, resolved, layout, runtime) -> Manifest` を唯一の実行 API とする。処理の順序は次のとおり。
    1. 設定を検証する
    2. consumes を検査する
    3. `stage.run` を呼ぶ
    4. produces を検査する
    5. manifest を原子的に書き込む
    6. run_state を更新する
  - `run_pipeline(resolved, run_dir, *, start=None, stop=None, dry_run=False, retry_failed=()) -> RunLayout`。`--from X` を指定すると、run_state 上の X とそれ以降を stale にする(AR-07)。実行中は `SiteLock` を保持する。
- **`pipeline/preflight.py`** は次を検査する。
  - 実行ファイルがあること
  - 版数の pin(`engines.requirements_for` と version_command)
  - worker 側の Python モジュール
  - `engine.supports(method)`
  - scratch が `/mnt/c` 上にないこと
- **CLI**(`hfauto/cli/main.py`)のコマンドは次の 5 つである。PIPELINE には `configs/pipelines/<name>.yaml` の name、またはパスを指定する。
  - `doctor [--site SITE]`
  - `run PIPELINE --system S --site SITE [--run-dir DIR] [--from ID] [--to ID] [--dry-run] [--retry-failed kinds]`
  - `status RUN_DIR`: run_state、FailureKind 別の件数、JobStore の再利用率を表示する
  - `report RUN_DIR`: view 全体から `reporting.html.render` を呼ぶ
  - `case RUN_DIR REACTION_ID`

### 7.5 設定の 4 分割(AR-10。層どうしをマージしない)

| 層 | 置き場所 | 中身 | 例 |
|---|---|---|---|
| site | `configs/sites/*.yaml` | 実行ファイルのパス、scratch、cores / memory、エンジン別の `EngineSite`(版数の pin、`ExecutionSpec`、env、worker 用の Python) | `wsl_local.yaml`: 4 vCPU / 11,000 MB。nwchem 7.2.3 は ranks 4、1,200 MB/rank、timeout 14,400 s。xtb は threads 1、crest は threads 4。pysis・readuct・goodvibes の python は `/home/user/.venvs/hfauto-prod/bin/python`。`.wslconfig` を 16 vCPU / 48 GB に拡張したときの値をコメントで併記する |
| method | `configs/methods/*.yaml` | `MethodSpec`(ハッシュの対象) | `gfn2`、`pbe0-d3bj_def2-svpd`、`pbe0-d3bj_def2-svp`、`pbe0-d3bj_def2-tzvp`(エネルギー層)、`wb97x-d3_def2-tzvpd`(パネル)、`ccsd-t_def2-tzvp`(参照 SP、小さい分子だけ) |
| system | `configs/systems/*.yaml` と `configs/systems/xyz/<system>/*.xyz` | 化学種(xyz / SMILES、電荷、多重度、役割)、組成、宣言した反応(反応物、生成物、任意の座標、torsional) | hcn、hono、nh3_inversion、formaldehyde、tma_hf2、amine_hf_panel、water_same_basin |
| pipeline | `configs/pipelines/*.yaml` | 化学的な処理の流れ、stage の設定、`conditions`、`gates:`(Policy の上書き) | discover、known_endpoints、method_panel |

- 絶対パスは site ファイルにだけ書く。
- 次の設定は廃止する: `global.mode`、`require_real_qm`、`fallback_to_dummy`、`allow_subprocess`、`extra_args`、`engine_order`、`ts_mode_overlap_threshold`、`--set`。

---

## 8. 化学計算プロトコル

### 8.1 共通方針

- **虚振動の 3 段方針**(`gates.imaginary_tier`)
  - ν < −50 は本物の虚振動、−50 ≤ ν < −10 は soft、−10 ≤ ν < 0 はノイズとして扱う。
  - 振動数は、全エンジンに共通の補空間射影で求めたものだけを使う。
- **PES の同一性は 2 層で判定する。**
  - 停留点・freq・QRC・熱化学: `full_key` の一致を要求する。
  - DFT の経路生成: `surface_key` の一致でよい。結果は種としてしか使わないためである。
- **電荷と多重度**はすべてのエンジンに渡す。開殻系では <S²> を観測し、スピン汚染があれば rankable の blocker にする。
- **資源**
  - 同時に使うコア数の合計は site.cores 以下とする(JobRunner のセマフォで制御する)。
  - メモリは 1 rank あたり 1,000〜1,500 MB とし、Σ(memory × ranks) が 0.8 × MemTotal 以下になるよう site ファイルで設定する。
  - スレッド数: NWChem は `OMP_NUM_THREADS=1`、xTB と CREST は `OMP_NUM_THREADS=n,1` と `OMP_STACKSIZE=4G` とする。
  - scratch は WSL の ext4 上に置く。
- **予算**
  - 1 反応あたり 6 h までとする。saddle は 2 回、DFT string は 1 回(多解像度での確認を 1 回まで追加できる)、QRC は 2 振幅までとする。
  - 各 attempt の timeout は min(設定値, 残り時間) とする。

### 8.2 stage ごとの既定値・停止条件・早期打ち切り

| stage / 段 | アルゴリズムと既定値 | 収束・終了 | 早期打ち切り・失敗時 |
|---|---|---|---|
| structures | ・xyz はそのまま使う。<br>・SMILES は RDKit ETKDG で 1 構造にする(`smiles.py`)。<br>・`resolve_electronic_state` を実行する。<br>・宣言反応の両端は、元素列(原子順序)が一致していなければならない | 全化学種で電子数と多重度が整合していること | ・`'.'` を含む SMILES は `INPUT_INVALID`(錯体は組成で表す。CH-23)。<br>・宣言反応の原子順序が一致しなければ `INPUT_INVALID` |
| conformers(単量体) | `crest --gfn2 --quick -T 4 --ewin 6 --chrg --uhf` | ・rc=0 で、`crest_conformers.xyz` があること。<br>・CREGEN の除去数を記録する | ・重原子が 3 個以下で、回転可能結合が 0 本なら、CREST を省いて入力をそのまま渡す。<br>・トポロジー停止したら `--noreftopo` で 1 回だけ再実行する |
| conformers(錯体) | 1. `placement.seeds` で配置 seed を作る(chem 10)。<br>・**H 結合の配置**: donor は極性 H だけとし、受容原子の lone-pair 円錐(X–H 軸から 110〜120°)上に、方位 0/120/240° で置く。固定 seed の RNG で ±30° 回転させたものを 2 つずつ作り、計 6 seed とする。<br>・**フォールバック**: donor か受容原子がない場合は、剛体をランダムな向きに置く。接触距離はファンデルワールス半径の和 + 0.5 Å とする。<br>・3 断片以上は逐次に積み上げる。<br>・衝突判定: H···H ≥ 1.2 Å、重原子間 ≥ 2.2 Å。棄却サンプリングは最大 50 回。<br>2. 最良の seed から `crest --nci --quick -T 4` を実行する。<br>3. 酸塩基の組成(labile H と受容原子の両方がある)では、`--notopo <labile H, 受容原子>` を自動で付ける(CH-18、CH-19) | 状態ラベルごとに、4 kcal/mol 以内の上位 6 個を species として出力する | ・有効な seed が 0 個なら、その組成を `GATE_REJECTED`(`no_collision_free_seed`)にする。<br>・CREST が失敗したら、seed をそのまま出力する(`source="placement"`) |
| minima(screen) | ・xTB GFN2 で `--opt vtight` → `--hess` を実行する。<br>・`relax_to_minimum` を使う(mode-follow は最大 2 サイクル)。<br>・`thread_map` で並列に実行する | `is_minimum` に合格した構造を Registry に登録する | ・非収束なら、最終構造から 1 回だけ継続する。<br>・状態ラベルが seed と違えば、崩壊として members に記録する(棄却はしない。崩壊した seed は DFT に流さない。CH-06、CH-20) |
| explore | **出発点**: 組成と状態ラベルの組ごとに、screen エネルギーが最も低い 2 極小。<br>**trial の生成**: `trials.generate` で次の優先順に作り、出発点あたり最大 10 件で打ち切る(chem 9)。<br>1. 極性 H の移動(n ≥ 2 では、H を介した relay を単純な移動より優先する)<br>2. C–H を含む 1,2- と 1,3- の H 移動(3 Å 以内)<br>3. 重原子の結合の形成・切断(変化は 2 本まで、3.5 Å 以内、芳香環は除外)<br>4. 断片の会合<br>**直線分子**: 10° 曲げ、RNG で 0.05 Å 変位させてから NT2 を試す(`perturbed=True`、CH-30)。<br>**実行順**: NT2(`nt_total_force_norm 0.1`)を先に試す。生成物がなければ AFIR(γ 125 → 300 kJ/mol) | attempt ごとに product / negative / failed のいずれかになる | ・ΔE‡_xTB > 150 kJ/mol または ΔE_rxn > 100 kJ/mol(300 K で評価)なら `out_of_window`(CH-28)。<br>・IRC が出発構造につながらなければ `irc_not_connected_to_source`。<br>・TS 検証のない NT2 回収構造と、未収束の AFIR 構造は生成物にしない。<br>・緩和で状態ラベルが変わった極小は、`mechanism="relaxation"` の発見として記録する。<br>・ねじれ異性化は trial にせず、`hypotheses` の conformer ペアで扱う |
| minima(dft) | **選択**(`include: window`): 次を対象にする。<br>・組成と状態ラベルの組ごとに、screen エネルギーで 6 kcal/mol 以内の最大 3 個<br>・explore の生成物<br>・入力単量体の最低構造(会合熱化学用)<br>`rerank_sp: true` なら、組成ごとの上位 8 個を DFT SP で並べ直してから選ぶ。`include: all` なら、入力 species をすべて使う。<br>**計算条件**: PBE0-D3BJ/def2-SVPD、grid fine、SCF 1e-7、driver default、maxiter 100、`coordinates auto`、`nocenter noautosym`。opt と freq は別ジョブにする。<br>**初期 Hessian**: 2 断片以上の species では、開始構造での xTB Hessian を `init_hessian` に使う(chem 20) | `is_minimum` と `spin_ok` に合格した構造を Registry に登録する | ・既存の極小と same なら、freq を省く(`known`)。<br>・timeout 時は最大 2 回まで継続実行する(AR-08)。<br>・SCF 非収束なら、damping をかけて 1 回だけ再試行する |
| reaction-paths: 仮説選択(`hypotheses.select`) | **候補**(優先順): ① system の宣言反応、② explore の product(低レベル TS があるものを優先)、③ mode-follow の `ts_candidate`、④ 同じ組成・同じ状態ラベルで、結合変化がなく二面角が 30° 以上違う conformer ペア。<br>**条件**: 別 basin であること、ΔE_rxn ≤ 40 kcal/mol、結合変化・宣言座標・座標変化(0.2 Å または 30° 以上)のいずれかがあること(CH-07)。組成あたり最大 6 件。<br>**端点構造**: `pick_endpoints` で、各 basin の members のうち、相手との恒等写像 RMSD が最小のものを選ぶ。代表構造は使わない(chem 13)。<br>**付随する値**: `torsional` は、宣言があるか、結合変化がない場合に真とする(chem 14)。`n_h_transferred` は `topology.transferred_hydrogens` で求める。<br>**縮退反応**: 宣言反応の両端が同じ basin に割り付いても、`mapped_equivalent` が真なら `degenerate=True`、`minima=(m, m)` として続行する(CH-35)。<br>**陰性証拠**: 同じ組成・同じ結合変化の negative を `negative_evidence` に集める | 仮説ごとに ReactionRecord を初期化する | same_basin なら直ちに COMPLETE |
| reaction-paths: SCREEN | §7.3 のとおり。<br>・低レベル TS がある場合は、3 点で済ませる近道を使う。<br>・ない場合は、xTB で端点を得たあと、GS(11 ノード、CI、同じ入力内で TSOpt)を実行する。xTB で両端が崩壊した場合は、代わりに DFT の IDPP を使う。<br>・全ノードで DFT SP を取り、`barrier_verdict` で判定する。内部極大が 1.0 kcal/mol 以上なら proceed。barrierless と判定するのは、DFT プロファイルが単調で、かつ低レベルも単調な場合だけ | BarrierVerdict を出したら終了する(反復しない) | ・単調なら `barrierless` とし、ΔE_max、`max_node_spacing_A`、`below_zpe` を添える。<br>・陰性証拠がある場合は、override がない限り proceed にしない。<br>・低レベル計算が失敗したら `unavailable` として FIND_PATH へ進む。理由は記録し、fail-open にはしない |
| reaction-paths: REFINE_SADDLE | ・初期 Hessian は、1 回目は xTB(条件付き)、2 回目は DFT のものを使う。<br>・NWChem `saddle` を driver default、`trust 0.1`、`sadstp 0.1`、maxiter 50、`inhess 2` で実行する。<br>・最低モード以外を追う場合は `moddir` を使う | driver が収束すること | ・timeout なら 1 回だけ継続する。<br>・最適化に失敗したら、次の種か FIND_PATH へ進む |
| reaction-paths: VALIDATE_TS | 別ジョブの freq を取り、`is_first_order_saddle` と `spin_ok` で判定する | SaddleClaim を作る | ・n_imag = 0 なら崩壊とみなし、VALIDATE_INTERMEDIATE へ進む(CH-33)。<br>・高次の場合は、2 本目のモードに沿って下った種から saddle をやり直す(試行回数の範囲内)。<br>・soft な 2 本目は注記だけ付ける |
| reaction-paths: FIND_PATH(DFT string) | ・NWChem ZTS: nbeads 9、stepsize 0.05、interpol 3、impose、freeze1 / freezeN、**tol 1e-5**、maxiter 20 × 最大 3 チャンク。<br>・初期経路は、低レベル GS の経路を再標本化したもの(写像付きで整列)。なければ IDPP。<br>・LOT は同じ汎関数・分散補正・電荷とする(`surface_key` の一致) | `string_converged`: 最終 gmax ≤ 1.0e-3 Eh/bohr で、直近 3 反復で悪化していないこと(CH-08) | ・`stagnated`(10 反復で gmax の改善が 5% 未満、または Xmax = 0 が 5 反復続く)で停止する。<br>・未収束でも内部極大が 1 つだけなら、HEI を種としてだけ使う |
| reaction-paths: CONNECT(QRC) | ・振幅はエネルギー目標 max(3 × drop, 3e-4 Eh) から曲率で求め、0.05〜0.4 Å の範囲に収める。2 回目は 2 倍にする。<br>・同じ LOT で opt → `Registry.find`。未知なら `relax_to_minimum` で新規登録してから `gates.connection` で判定する | ConnectionClaim を作る | ・ねじれ座標は、周期を考慮した最近傍で割り付ける(CH-04)。<br>・宣言と別の組につながった場合は、reassigned として採用し、元の仮説は未解決に戻す |
| sp | 指定した `methods:` の各 LOT で、反応に参加する極小と TS の固定構造の SP を計算する(WFT の SP を含む) | Evidence | ・エネルギー層(def2-TZVP)は composite 用、パネルは report 用 |
| thermo | **GoodVibes 4.3.0 の API**: `QS=grimme`、`s_freq_cutoff=100`、`symm=True`。<br>**スケール因子**: 振動用と ZPE 用を分ける(`chemistry.thermo.scale_factors`)。<br>・GoodVibes 同梱の Truhlar DB で、汎関数と基底の完全一致を探す。<br>・なければ同じ汎関数の値を使い、`scale_factor_from:<level>` と注記する。PBE0 では PBE0/MG3S の 0.989 / 0.975(Alecu 2010)を使う。<br>・それもなければ 1.0 とし、`scale_factor_unverified` と注記する。<br>**標準状態**: 気相 1 atm で計算し、1 bar・1 M への換算は純関数で行う(Δn あたり 1.89 kcal/mol、298.15 K・1 M)。<br>**composite**: G = E_SP(エネルギー層)+ (G_GV − E_GV)。<br>**会合量**(組成に構成単量体がある場合。chem 12): dG_assoc = G(反応物錯体) − Σ G(単量体)、dG_act_vs_separated = G(TS) − Σ G(単量体)。いずれも `standard_states` ごとに Δn = 1 − Σn_i の補正をかける。<br>**感度幅**: qs {grimme, truhlar} × cutoff {50, 100, 150}。<br>**basin 分布**: 同一組成内の Boltzmann 重み | `thermo_consistent` に合格すること | ・不合格や GoodVibes の失敗なら G = None(`thermo_unavailable`)とし、再試行しない。<br>・LOT が一致しなければ `mixed_level_of_theory`(CH-40) |
| report | ・**順位付け**: `rankable` を満たす反応だけを、ΔG‡(T) で並べる。感度幅が重なる反応は同順位とする。<br>・**被覆率**: 機構ごとの試行数、生成物数、陰性結果の理由別件数、FailureKind 別件数。<br>・**手法パネル**: `Level.full_key` をキーに ΔE の幅を出す(CH-25)。パネルの手法間で ΔE_rxn の符号が一致しない反応には、blocker `method_sign_disagreement` を付ける | ReportRecord と、`ranking.csv` / `coverage.csv` / `method_panel.csv` / `report.html` を出力する | 疑似的な confidence は出さない(AR-27) |

### 8.3 既定値の一覧(コード上の置き場所)

| 項目 | 既定値 | 置き場所 |
|---|---|---|
| 虚振動の 3 段 | −50 / −10 cm⁻¹(ねじれ TS は −20) | `gates.Policy` |
| レジストリ | same は ≤ 0.02 Å かつ ≤ 1e-5 Eh、distinct は > 0.05 Å または > 5e-5 Eh | `identity.compare_minima` の既定引数 |
| 端点の割付け | ≤ 0.05 Å かつ ≤ 5e-5 Eh で、次点との差が 3 倍以上または 0.1 Å 以上 | `identity.assign` の既定引数 |
| 結合判定 | 1.15 / 1.45 × Σr_cov | `topology` の定数 |
| 仮説 | ΔE_rxn ≤ 40 kcal/mol、変化 ≥ 0.2 Å または ≥ 30°、組成あたり 6 件 | `hypotheses.select` の既定引数(窓は `Policy`) |
| 障壁判定 | 内部極大 ≥ 1.0 kcal/mol、GS 11 ノード、全ノードで DFT SP | `Policy`、`CasePolicy` |
| TS の prominence | max(2e-5 Eh, 20 × SCF 閾値) | `Policy` |
| QRC | 目標 3e-4 Eh、振幅 0.05〜0.4 Å、2 回目は 2 倍、降下量 ≥ max(1e-5 Eh, 20 × SCF 閾値) | `CasePolicy`、`Policy` |
| DFT(停留点) | PBE0-D3BJ/def2-SVPD、fine、1e-7、driver default、maxiter 100 | `configs/methods/pbe0-d3bj_def2-svpd.yaml`、NWChem renderer |
| 初期 Hessian | 2 断片以上の DFT opt と saddle の 1 回目は xTB、saddle の 2 回目は DFT | minima の `init_hessian`、`actions.REFINE_SADDLE` |
| DFT string | tol 1e-5、gmax 1e-3、直近 3 反復で悪化なし | NWChem renderer、`profile` |
| 熱化学 | qRRHO(grimme、100 cm⁻¹)、スケール因子は表から、`symm`、1 atm で計算して換算、ZPE の許容差 1e-5 Eh、E の許容差 1e-6 Eh | `ThermoSettings`、`chemistry.thermo`、`Policy` |
| 順位付けの ZPE 許容差 | 1 + 4 × max(1, n_H) kcal/mol | `Policy` |
| スピン汚染 | \|<S²> − S(S+1)\| ≤ 0.1 | `Policy` |
| 予算 | 1 反応 6 h、saddle 2 回、string 1 回(+ 確認 1 回)、ジョブの timeout 4 h | `CasePolicy`、site の `ExecutionSpec` |

---

## 9. 削除一覧

### 9.1 Wave 1 で削除するもの(新設計に依存せず、残るコードからも参照されないもの)

| 対象 | 行 | 理由 |
|---|---:|---|
| `hfauto_ops/**`、`hfauto_viz/**`。可視化の 3 機能は Wave 5 で `hfauto/reporting/html.py` に作り直す(参照元は `git show 5d76201:<path>`) | 3,906 + 2,426 | 実行実績がほぼ 0。境界違反があり、hpc と三重に実装されている(AR-02〜04) |
| `hfauto/hpc/**`、`hfauto/stages/{hpc_plan,ops,viz}.py`、CLI のコマンド `hpc-plan`、`cache-index`、`retry-plan`、`compare-runs`、`science-status`、`next-actions`、`outputs`、`compare-backends` | 628 + 221 + 約 140 | 同上 |
| `scripts/**` | 2,844 | ゲートを迂回するもの、または起動できないもの(AR-05、CH-03) |
| `hfauto/stages/{enrich,descriptors,kinetics,calibrate,connector_audit,rank,recover_path}.py` | 1,561 | 目的外、legacy、fail-open |
| `hfauto/backends/db/**`、`backends/kinetics/**`、`backends/thermo/{arkane,simple}.py`、`backends/qm/orca.py`、`backends/ts/orca_nebts.py`。あわせて `hfauto/workflow/reaction_state.py:23` の文字列 `"orca_nebts"` を消す | 971 + 365 + 203 + 1,329 | 未使用、または実績 0(CH-24、CH-43) |
| `hfauto/core/{thermochemistry,thermo,reference_data,public_data,artifact_types}.py`、`core/schemas/{calculation,species,site,descriptor,reaction}.py` | 約 1,000 | 死蔵、または重複(AR-25) |
| `hfauto/chemistry/{descriptors,path_validation}.py`(QM バックエンド 3 つの descriptors 呼び出しも消す)、`hfauto/reporting/molecule_dossier.py`、`hfauto/data/**` | 約 320 | 死蔵 |
| `configs/{stages,hpc,ops,viz}/**`(15 本)、`configs/pipelines/{amine_hf_water_extension,amine_hf_pilot_dft_followup,hcn_alternative_pes_validation,hcn_numerical_sensitivity,same_basin_water_nwchem,tma_hf2_recovered_path,tma_hf2_validated_minima_path}.yaml` | 約 900 | 参照 0、または代替がある |
| `examples/{candidates.sdf,db_fixtures/,generic_reaction/,reference_sets/,production/amine_hf_panel.json,production/amine_hf_full.sdf,m3_trimethylamine_hf2/README.md}` | — | 参照 0。panel の SMILES は `git show 5d76201:examples/production/amine_hf_panel.json` から Wave 5 で移す |
| `TECHNICAL_REPORT.md`、`docs/{architecture,developer_guide,stage_contracts,hf1_scan_pilot}.md`、`docs/current/{index,production_science_validation_plan,extension_backlog}.md` | 約 260 | 古い、または実態と矛盾している |
| リポジトリ直下の `afir_biased/`、`nt2_guess/`、`product_unbiased/`、`ts_bofill/`、`ts_hessian/`、`internal_coords.log`、孤立した `__pycache__` | — | CWD への漏れ(CH-31)。`.gitignore` にも追加する |
| テスト: `test_ops_current_stages`、`test_viz_cli`、`test_viz_current_tables`、`test_energy_visualization`、`test_pysisyphus_irc_manifest`、`test_amine_hf_panel`、`test_amine_hf_m3_candidate`、`test_recover_path`、`test_thermo_scientific_tasks`、`test_reaction_inputs` | 約 700 | 削除対象のテスト、設定値を固定しているだけのもの、旗の形しか見ていないもの |
| pyproject: import-linter の契約 10 本(`core-is-foundation` だけを残す)、entry point `hfauto-ops` と `hfauto-viz`、extras の web / viz / viz-prod / hpc / chemiscope | — | Wave 3 で新しい契約を別ファイルに入れ、Wave 6 で最終形の 5 契約に置き換える |

### 9.2 Wave 5 で削除するもの(切替と同時)

- `configs/pipelines/` の残り 6 本と、`examples/` の全体。xyz は `configs/systems/xyz/<system>/` へ `git mv` し、`examples/production/amine_hf_pilot.sdf` は削除する。
- `examples/` と `configs/` を参照している旧テスト 4 本: `test_reaction_case_workflow.py`、`test_reaction_trials.py`、`test_saddle_seeds.py`、`test_workflow_architecture.py`。価値のあるケースは Wave 4 で移植済みである。
- `hfauto/reporting/html_report.py` はここでは消さない。旧 `hfauto/workflow/runner.py` がまだ import しているので、Wave 6 で削除する。

### 9.3 Wave 6 で削除するもの(旧実装の一括削除)

- `hfauto/stages/` の旧モジュールすべて: `base.py`、`registry.py`、`basin_populations`、`build_complexes`、`conformers`、`connect_minima`、`detect_sites`、`dft_minima`、`discovery_audit`、`endpoint_seed_screen`、`endpoint_seeds`、`enumerate_states`、`explore_reactions`、`generate_reactions`、`ingest`、`irc`、`method_panel`、`minimum_mode_assess`、`minimum_mode_follow`、`minimum_registry`、`path_ensemble`、`path_ensemble_assess`、`path_intermediates`、`preopt`、`reaction_classify`、`reaction_plan`、`reaction_rank`、`reaction_segments`、`relaxation_discovery`、`sp`、`thermo`、`thermo_sensitivity`、`ts_search`
- `hfauto/backends/` の旧実装: `base.py`、`registry.py`、`conformer/**`、`qm/**`、`reaction_discovery/**`、`thermo/**`、`ts/**`
- `hfauto/chemistry/` の旧モジュール: `basin_identity`、`complex_seed_ensemble`、`complexes`、`connectivity`、`geometry_qc`、`method_lineage`、`methods`、`microstates`、`minima`、`minimum_connections`、`minimum_recovery`、`mode_following`、`nwchem_evidence`、`path_convergence`、`path_diagnostics`、`path_initialization`、`populations`、`proton_transfer`、`rdkit_utils`、`reaction_classification`、`reaction_path_qc`、`reaction_profile`、`reaction_trials`、`reactions`、`saddle_seeds`、`site_detection`、`stoichiometry`
- `hfauto/core/`: `artifacts.py`、`config.py`、`environment.py`、`executables.py`、`frequency_qc.py`、`qc.py`、`thermo_models.py`、`io.py`、`schemas/**`、および `constants.DEFAULT_HF_BOND_A` と `DEFAULT_HF_STRETCH_CM1`
- `hfauto/workflow/**`、`hfauto/reporting/html_report.py`
- `tests/test_*.py`(直下にある旧テストの残り全部)
- `importlinter-next.toml`(最終契約は pyproject に移す)
- 依存: `qcengine`、`cantera`、`requests`、`beautifulsoup4`、`graphviz`、`networkx`、`jinja2`、`plotly`、`matplotlib`、`Pillow`、`snakemake`、`chemiscope`、`ase`
  - 本体の依存は `pydantic`、`typer`、`PyYAML`、`rich`、`numpy`、`scipy` とする。
  - extras は `chem`(rdkit)、`production`(rdkit、pysisyphus>=1.0,<2、numba、llvmlite、scine-readuct==6.1.0、scine-xtb-wrapper==3.0.2、`goodvibes==4.3.0`、`pymsym`。Linux のみ)、`dev` の 3 つだけにする。

### 9.4 Wave 8〜9 で削除するもの

- 静的解析(vulture、または grimp による到達可能性、radon、grep)で見つかった孤立コード。書き込まれるだけで読まれない payload のフィールドも含む。
- `docs/current/**`。本書は `docs/reviews/2026-09-25_refactor_design.md` に移して凍結する。

---

## 10. テスト戦略

### 10.1 構成と予算

```
tests/
  conftest.py         override_engine、tmp_run
  fakes.py            FakeQM(解析 PES)、FakePath、FakeSaddle、FakeConformers、FakeDiscovery、FakeThermo(W3 で凍結)
  unit/               IO なし、外部プログラムなし。core / chemistry / execution / pipeline / drivers / backends(renderer・parser)/ reporting / cli
  golden/             conftest.py(golden fixture)、data/<engine>/、SOURCES.json、excerpt.py(決定的な抜粋)
  integration/        conftest.py(fake_runtime)。fake エンジンで stage をつなぐテストと、パイプライン全体のテスト
  smoke/              conftest.py(real_engine fixture)。marker `real`。既定では除外し、WSL で HFAUTO_REAL=1 のときだけ実行する
```

- pyproject(Wave 1 で設定する):
  - `markers = ["golden", "integration", "real"]`
  - `pythonpath = ["tests"]`
  - `addopts = "-m 'not real' --strict-markers --import-mode=importlib"`
  - importlib モードなので、テストの basename が一意である必要はない。
  - rdkit、pysisyphus、scine、goodvibes を使うテストは `importorskip` で保護する。
- 目標は、全体で 5,000 行以下、既定の実行で 15 s 未満とする。各 WP のテスト行数の上限は次のとおりで、合計は 4,950 行になる。

| WP | 上限(行) | WP | 上限(行) |
|---|---:|---|---:|
| 1.3 golden 基盤 | 250 | 4.2 xtb・pysis・minima(smoke を含む) | 350 |
| 2.1 core とゲート | 450 | 4.3 reaction-paths | 300 |
| 2.2 化学カーネル(golden を含む) | 500 | 4.4 structures・conformers(smoke を含む) | 250 |
| 3.1 execution | 300 | 4.5 explore(smoke を含む) | 250 |
| 3.2 fakes と engines・MinimumDriver | 500 | 4.6 thermo・sp(smoke を含む) | 250 |
| 3.3 pipeline と integration/smoke の conftest | 300 | 4.7 report summary | 150 |
| 3.4 仮説・分類・decide | 350 | 5.1 と 5.2 の CLI・configs・html | 150 |
| 4.1 NWChem(smoke を含む) | 400 | 6.2 e2e と再開 | 200 |

- monkeypatch は 10 回以下に抑え、`engines.override` と `catalog.override` を使う。
- YAML の設定値を固定するテストは作らない。configs は `integration/test_configs.py` で、モデルの検証と不変条件だけを確認する。
- レビューで誤りとされた挙動(M5 §2.4 の 15 件)を期待値にしている旧テストは、移植しない。期待値を反転した新しいテストを書く。

### 10.2 golden フィクスチャ(Wave 1 で抽出。目標 1.2 MB、上限 2 MB)

- **抜粋の規則**(`tests/golden/excerpt.py`。決定的に抜粋し、使ったコマンドを SOURCES.json に記録する)
  - NWChem の出力で残すもの:
    - 版数、`Grid used`、`Convergence on energy requested`、基底の種別、DFT-D3 / COSMO / `<S2>` の行
    - 入力のエコーと、最初と最後の Geometry ブロック
    - `Total DFT energy` の行
    - すべての振動解析ブロック(`P.Frequency` を含む。G01 ではブロック数も検査に使う)
    - 最適化の収束・失敗メッセージ
    - `string:` と `@zts` の行
    - 末尾の 50 行
  - SCF の反復と勾配の出力は捨てる。
  - 50 KB を超えるファイルは gzip する。
- **ディレクトリ名の制約**: `.gitignore` に当たる名前(`runs/`、`cache/`、`*.log`)を data の配下で使わない。`*.log` のファイルは `.txt` に改名し、元の名前を SOURCES に記録する。
- **xfail(strict) の規則は廃止する**(eng 16)。golden は、それを使う新コードの WP が初めて検査する。

| ID | 元ファイル(`runs/` からの相対パス) | 期待値 | 対応 |
|---|---|---|---|
| G01 | `hono_isomerization_v3/06_ts-search/rxn_candidate_dft_998cc65513/attempt_00_nwchem_neb/saddle_freq/nwchem.out` | 振動ブロックが 2 つある(ZPE 0.019317 と 0.018747)。パーサは最終ブロック(−680.11)だけを採る。この出力は freq ジョブとしては不正(ブロック数 ≠ 1)で、`INCOMPLETE_OUTPUT` になる | CH-38 |
| G02 | `hono_isomerization_v3/09_thermo/backend_runs/*/Goodvibes.csv` | TS の zpe 0.0375115 は 0.985 × 0.018747 と一致しない → `thermo_consistent` は不合格。24.11 は再現しない | CH-38 |
| G03 | `hono_isomerization_v1/02_dft-minima/spc_hono_{trans,cis}/nwchem.out` | G01 と合わせて ΔE‡ 13.62。ΔG‡ ≈ 12.2 kcal/mol(GoodVibes を実行するテストは `real`) | CH-38 |
| G04 | `hcn_readuct_nwchem_validation/04_dft-minima/spc_discovered_product_620510d259663999/nwchem.out`(HNC) | NWChem の RRHO は回転エントロピーが負になる → 使わない。Evidence はパースでき、直線分子として n_external = 5 | CH-39 |
| G05 | 同じ run の `spc_hcn/nwchem.out` | 直線分子で、振動数は 3N−5 本 | CH-39 |
| G06 | `hcn_goodvibes_v4/11_thermo/backend_runs/*/Goodvibes.csv` | ΔG_rxn 12.742、ΔG‡ 42.281 kcal/mol | CH-39 |
| G07 | `hcn_saddle_v3/work/saddle_frequency/{nwchem.out,hfauto_job.hess,final.xyz}` | `.hess` から**同位体質量**で射影した振動数が、NWChem の値(−1131.57)と 1 cm⁻¹ 以内で一致する。`frame_residual` はおよそ 6e-6(診断値の確認だけ)。`is_first_order_saddle` に合格する | CH-26、CH-34 |
| G08 | `m3_tma_p0_seed_hessian_001/frequency/{hfauto_job.hess,nwchem.nw}` | 射影後の虚振動は**本数 4 だけを検査する**(非停留点なので値は NWChem と約 30 cm⁻¹ ずれる)。`frame_residual` は大きい(約 1e-2)が、ゲートには使わない | CH-26、chem 2 |
| G09 | `hcn_irc_v5/work/pysis_irc.out` | 両分岐が TS より上にある → 分岐の軌跡として `connection` は不合格 | CH-03 |
| G10 | `hcn_irc_v5/work/qm_calcs/irc_000.001.qce_nwchem_stdout` | NWChem パーサが grid medium を観測する → xfine を要求した MethodSpec に対して `level_mismatches` が空でない。入力デッキと出力の初期構造が剛体回転でずれている → frame 検査で不合格 | CH-01、CH-02 |
| G11 | `m3_real_qm_ammonia_inversion_004/05_irc/rxn_ammonia_inversion/pysis_irc.out` | 後退側が step 0 で上昇している → 不合格 | CH-03 |
| G12 | `hono_isomerization_v3/07_irc/rxn_candidate_dft_basins_a4130bf1b004a237/{pysis_irc.out,*_endpoint_optimization/final.xyz}` と G03 の最終構造 | 降下の条件は合格するが、割付けは 4 通りとも閾値内 → 一意性で不合格 | CH-03、CH-04 |
| G13 | `m3_tma_path_ensemble_001/rxn_tma_hf2_proton_reorganization/images_11/string/nwchem_string.out`(抜粋) | gmax が 9.14e-4 から 9.15e-3 に悪化している → `string_converged` は False | CH-08 |
| G14 | `hono_isomerization_v3/06_ts-search/.../neb/hfauto_neb.neb_epath`(実在するパスに解決する) | 最後の 58 反復で Xmax = 0 → `stagnated` は True | CH-09 |
| G15 | `hcn_ts_irc_validation/07_ts-search/rxn_candidate_620_a0961220a0/attempt_00_nwchem_neb/neb/{hfauto_neb.neb_final.xyz,neb_final_epath}` | 最小の原子間距離が 0.529 Å → 初期経路の検査(0.7 Å)で不合格。同じ端点から作った IDPP は合格 | CH-09 |
| G16 | `sio2_wet_etch_aniline_001/cases/path00_control/02_xtb_scan/segment_b_departure/{stdout.log,stderr.log,result.json}` | rc=0 でも `FAILED TO CONVERGE` がある → `GEOMETRY_MAXITER` | BUG-07、CH-20 |
| G17 | `amine_hf_pilot_hf1_hf3_v3_discovery/04_preopt/spc_NH3_HF1_seed00_nci_0/attempt_00/{xtb.out,xtbopt.xyz}` | 正常に収束している | — |
| G18 | `amine_hf_pilot_hf1_hf3_v3/03_build-complexes/crest_nci/TMA_HF2/TMA_HF2_seed00/attempt_00/spc_TMA_HF2_seed00/crest/{crest.stdout,crest.energies,crest_conformers.xyz}` | 配座数とエネルギーが読める。欠損は None になる | BUG-08 |
| G19 | `amine_hf_pilot_hf1_hf3_v1/03_build-complexes/crest_nci/ANL_HF1/mol00002_state000/crest/crest.stdout` | トポロジー停止として分類される | CH-18 |
| G20 | `amine_hf_pilot_hf1_hf3_v3_discovery/pipeline.stdout.log`(`Vib. Frequencies` の 8 ブロックを抜粋)と、`07_explore-reactions/reaction_discovery_attempts.jsonl` の該当行 | 8 件とも虚振動は 1 本(−1173.2 〜 −239.2) | CH-26 |
| G21 | `same_basin_water_v1/01_dft-minima/spc_water_{reference,distorted}/nwchem.out` と最終 xyz | `compare_minima` は same(RMSD 8e-6 Å、ΔE 2.6e-10 Eh) | — |
| G22 | `amine_hf_pilot_dft_followup_v1/00_dft-minima/spc_NH3_HF3_seed00_nci_0/{nwchem.out の先頭 400 行と末尾 200 行, command_result.json}` | rc 124 → `TIMEOUT`。継続に使うファイルの一覧が得られる | AR-08 |
| G23 | `m3_tma_mode3_saddle_001/work/saddle_frequency/nwchem.out` の抜粋 | `GEOMETRY_MAXITER`。失敗した saddle の Hessian を TS として公開しない | AR-26 |
| G24 | `hono_isomerization_v1/06_ts-search/.../neb/command_result.json` | rc −15 → 中断(`NONZERO_EXIT`)として分類される | AR-26 |

### 10.3 実エンジンでの確認

- **Wave 4 の統合時**(WSL。eng 14)
  - `uv pip install --python /home/user/.venvs/hfauto-prod/bin/python -e <repo> --no-deps` を実行したうえで、`pytest -m real tests/smoke` を実行する。
  - Python API(`run_pipeline`)で、HCN の known_endpoints を一時設定で 1 回実行する。
  - smoke の内容:
    - xtb: H2O の opt と hess。`--cycles 2` で `GEOMETRY_MAXITER` になること。
    - nwchem: H2O の opt と、別ジョブの freq(Level の一致、振動ブロックが 1 つであること)。平面 NH3 の saddle(`inhess 2`)。`moddir` の経路を 1 件。HCN の MP2 SP。ωB97X-D3 の Level の観測。
    - pysis: HCN↔HNC の GS に `refine_ts` を付けたときに、TS の xTB freq で虚振動が 1 本であること。
    - crest: 小さい錯体で `--chrg` / `--uhf` が効くこと。
    - readuct: HCN の NT2 を曲げてから実行し、射影した振動数で数えること。
    - goodvibes: API と CLI の G が 1e-5 Eh 以内で一致すること。S_rot > 0。
  - 合計で 15 分以内、ranks は 2 以下とする。
- **Wave 7 の受け入れ**(WSL の実計算。結果は `docs/validation.md` に記録する)

| 系 | 期待値 |
|---|---|
| HCN→HNC(known_endpoints) | TS ≈ −1131i、ΔE‡ ≈ 46.7、ΔG‡ ≈ 42.3、ΔG_rxn ≈ 12.7 kcal/mol、`elementary_step`、所要 1 h 未満 |
| HONO trans→cis | TS ≈ −680i、ΔE‡ ≈ 13.6、ΔG‡ ≈ 12.2 kcal/mol(24.11 は再現しない) |
| NH3 反転 | ≈ −757i、≈ 4.25 kcal/mol、`degenerate_rearrangement`。`symm` の分として +0.41 |
| 水(同一 basin) | `same_basin`。TS 探索を起動しない |
| TMA·(HF)₂(discover) | reaction-paths stage が SCREEN で終わる(`barrierless_at_resolution` / `no_product_basin` / `same_basin` のいずれか)。**reaction-paths stage は 1 h 未満**(旧実装では約 30 h) |
| amine·HF pilot の探索(NH3·HF と TMA·HF の 2 組成に限り、`--to explore`) | 低レベルの TS 8 件が、すべて虚振動 1 本と数えられる |
| xTB の初期 Hessian の A/B | 2 断片以上の DFT opt で、ステップ数を初期 Hessian なしの場合と比べて記録する |

---

## 11. 移行手順(Wave)

各 Wave の終わりには、`$PY -m pytest -q -p no:cacheprovider`(既定の marker)が緑であること。

- **Wave 内の WP**: 所有するパスは互いに素とする。共有ファイル(`pyproject.toml`、registry、`__init__.py`、`tests/conftest.py`、`importlinter-next.toml`)は、その Wave の中で 1 つの WP だけが所有するか、統合手順に予約する。
- **新コードの扱い**: Wave 2〜4 で書く新コードは、§9.3 の旧モジュールを import しない。旧コードと並走させ、Wave 5〜6 で一括して切り替える。

**全 WP に共通する規則**

- **実行環境**: `PY=.venv-win-qa/Scripts/python.exe`。静的検査は `.venv-win-qa/Scripts/{ruff,lint-imports,radon,pyrefly}` を使う。
- **WSL の実行環境**: `/home/user/.venvs/hfauto-prod/bin/python`(scipy 1.18、goodvibes 4.3.0、pysisyphus 1.0.0、pymsym、scine)。Git Bash から `wsl.exe` を呼ぶときは `MSYS_NO_PATHCONV=1` を付ける。
- **新コードの禁止 import**(Wave 2〜5 の受け入れ条件で 0 件であること):

```
grep -rnE "hfauto\.(workflow|core\.(schemas|qc|artifacts|executables|environment|frequency_qc|thermo_models|config|io)|backends\.(qm|ts|thermo|conformer|reaction_discovery|registry|base)|chemistry\.(basin_identity|complex_seed_ensemble|complexes|connectivity|geometry_qc|method_lineage|methods|microstates|minima|minimum_connections|minimum_recovery|mode_following|nwchem_evidence|path_convergence|path_diagnostics|path_initialization|populations|proton_transfer|rdkit_utils|reaction_classification|reaction_path_qc|reaction_profile|reaction_trials|reactions|saddle_seeds|site_detection|stoichiometry)|stages\.(base|registry))\b" <新規パス>
```

- **旧コードの扱い**: 旧コードは移植の参考として読むだけにする。旧テストのうち価値のあるケース(M5 で K 判定)は、担当ドメインの新テストに期待値を直したうえで移植する。
- **新テストのデータ**: 新テストは `examples/` を参照しない。小さな xyz はテストの中に書くか、`tests/golden/data/` を使う。
- **サイズと複雑度**: 新しい関数は CC ≤ 15(`ruff check --select C901`)、新しいファイルは 500 行以下とする。
- **コミット**: WP はコミットしない。Wave の統合担当が、統合検査に合格した後に 1 回だけコミットする。メッセージは `refactor(Wn): ...` とし、末尾に規定の Co-Authored-By 行を付ける。
- **editable の再インストール**: pyproject を変えた WP は、`uv pip install --python .venv-win-qa/Scripts/python.exe -e . --no-deps` を実行する。
- **新パッケージの空の `__init__.py`**: `hfauto/execution`、`hfauto/drivers`、`hfauto/drivers/reaction_case`、`hfauto/pipeline` の空の `__init__.py` は、Wave 2 の WP2.1 がまとめて作る。これにより、Wave 3 の WP 間で `__init__.py` の所有が重ならない。

| Wave | 目的 | WP |
|---|---|---|
| 1 | 死蔵・未使用・legacy の削除、golden の抽出、pytest と依存の準備 | 1.1 add-on / hpc / scripts / configs / docs / 衛生 / pyproject、1.2 未使用 stage とバックエンドの連鎖・死蔵モジュール、1.3 golden |
| 2 | 型と純関数の土台 | 2.1 core の型(evidence / records / manifest / method / system)とゲート、2.2 化学カーネル |
| 3 | 実行・能力・配線・判断の土台 | 3.1 execution、3.2 protocols・engines・fakes・MinimumDriver、3.3 pipeline 層・stage spec / catalog・integration と smoke の conftest・importlinter-next、3.4 仮説・分類・decide(純関数) |
| 4 | 化学ドメインの再構築と、実エンジンでの smoke | 4.1 NWChem、4.2 xTB・pysis・minima、4.3 reaction-paths の action・driver・stage、4.4 structures・conformers、4.5 explore、4.6 thermo・sp、4.7 report の集計と stage |
| 5 | 切替 | 5.1 CLI・configs・examples の移動、5.2 HTML レポート |
| 6 | 旧実装の一括削除と e2e | 6.1 旧コードと旧テストの削除・pyproject の最終化、6.2 e2e と再開のテスト |
| 7 | WSL の実計算による受け入れ | 7.1 実計算・最小限の修正・validation.md |
| 8 | 最終整理 | 8.1 孤立コードの削除と定量目標の検証 |
| 9 | 文書 | 9.1 README / design / environment と本書の凍結 |

### 11.1 import-linter の契約

- Wave 1 では、pyproject の契約を `core-is-foundation` だけに減らす(exhaustive の layers 契約と、`exhaustive_ignores` を削除する)。
- Wave 3 では、WP3.3 が `importlinter-next.toml` を作り、新しいモジュールだけを対象にした forbidden 契約を入れる。統合時に `lint-imports --config importlinter-next.toml` がすべて KEPT であることを確認する。
  - `next-core`: 新しい core → chemistry / backends / execution / drivers / stages / pipeline / workflow / core.schemas / core.io / core.artifacts を禁止する。
  - `next-chemistry`: Wave 2 で作った chemistry モジュール → execution / backends / drivers / stages / pipeline / workflow / core.schemas を禁止する。
  - `next-execution`: execution → backends / chemistry / drivers / stages / pipeline / workflow を禁止する。
  - `next-protocols-only`: drivers、stages.spec、stages.catalog → backends.engines / execution / pipeline / workflow を禁止する。
  - `next-pipeline`: pipeline → workflow / core.schemas / core.config / stages.registry / stages.base を禁止する。
- Wave 4 の統合時には、対象を広げる。
  - 新しい stage と reporting.summary を `next-protocols-only` の source に加える。
  - `next-adapters-independent`(独立性: nwchem、xtb、pysis、crest、readuct、goodvibes)を追加する。
- Wave 6 では、pyproject を最終形の 5 契約に置き換え、`importlinter-next.toml` を削除する。
  1. `hfauto-layers`(layers、containers `hfauto`、exhaustive): `cli > pipeline > stages > drivers | reporting > backends > execution | chemistry > core`
  2. `adapters-private`(forbidden): `hfauto.{cli,pipeline,stages,drivers,reporting}` → `hfauto.backends.{nwchem,xtb,pysis,crest,readuct,goodvibes}`
  3. `protocols-only`(forbidden): `hfauto.{stages,drivers,reporting}` → `hfauto.backends.engines`、`hfauto.execution`
  4. `adapters-independent`(independence): 上の 6 アダプタ
  5. `stages-independent`(independence): 8 つの stage モジュール
- ruff については、`C901`(max-complexity 15)と `TID251`(banned-api。`subprocess` は `hfauto/execution/process.py` だけに許可する)を Wave 6 で有効にする。
- エンジン名リテラルの検査は、アーキテクチャテストにせず、統合検査の grep で行う(AR-14、eng 24)。

---

## 12. 既知の致命バグの修正方針

| ID | 重大度 | 根本原因 | 修正(新設計での場所) | WP | 検証 |
|---|---|---|---|---|---|
| CH-01 | 致命 | pysisyphus の QCEngine 計算器が、NWChem の再配向フレームのまま勾配と Hessian を返す | DFT 段での pysisyphus+QCEngine を**削除**し、`qcengine` への依存も外す。DFT の saddle は NWChem のネイティブ実行(`nocenter noautosym`)で行う。接続確認は QRC で行う。pysisyphus は xTB のネイティブ計算器でだけ使う。フレームの保持は Evidence の保証 (5)(出力の初期構造と入力の照合)で担保する | 2.2、4.1、4.2、6.1(`qcengine` の削除) | G07 は合格、G10 は不合格 |
| CH-02 | 高 | QCEngine 経路の grid と SCF が要求と異なり、provenance には要求値が記録される | 経路ごと削除する。一般規則として、`Level` を出力から観測し、`level_mismatches` があれば `METHOD_MISMATCH` とする。停留点の層は `same_pes(numerics=True)` で確認する | 2.1、4.1 | G10 は不合格 |
| CH-03 | 致命 | IRC の合否が分岐の下降を見ておらず、最適化後の端点の照合だけで決まっている | IRC を廃止し、QRC に一本化する。`gates.connection` が、変位直後の降下、軌跡の上限、最終値の降下、一意の割付け、集合としての一致を判定する | 2.1、4.3 | G09、G11、G12 は不合格。fake PES では elementary で合格 |
| CH-05 | 致命 | reaction-plan と classify の間で dict を受け渡しており、必須キーが欠けている | dict の受け渡しを廃止する。`ReactionRecord` と `MinimumRecord` を型のまま渡し、分類は claim の有無だけで決める | 2.1、3.4、4.3 | integration: 宣言反応が fake PES で `elementary_step` になる |
| CH-26 | 致命 | ReaDuct の TS を未射影の固有値で数えており、真の 1 次鞍点をすべて棄却している | 全エンジン共通の `projected_frequencies`(補空間での対角化と同位体質量)、`automatic_mode_selection`、IRC の出発構造との一致確認 | 2.2、4.5 | 合成 Hessian(未射影では 3 本、射影後は 1 本)、G20、G08 |
| CH-38 | 致命 | saddle ジョブ内の 2 つの振動ブロックを GoodVibes が連結し、ZPE を二重に計上する | opt / saddle と freq を別ジョブにする。GoodVibes には**自前の振動数を API で渡し**、テキストは読まない。`thermo_consistent` で検査する。HONO の ΔG‡ 24.11 の撤回を `docs/validation.md` に記録する | 4.1、4.6、7.1 | G01、G02 は不合格。実計算で ≈ 12.2 |
| CH-39 | 高 | GoodVibes が失敗すると内部フォールバックの熱化学(HNC の負の S_rot)で続行し、旗は常に偽のまま | フォールバックを削除し、G = None(`thermo_unavailable`、`rankable` は不合格)とする。回転定数は自前で計算して渡し、S_rot > 0 を確認する。tier・confidence・旗は削除する | 4.6 | G04、G06(12.74)。HCN のフォールバック値 18.05 は出ない |

同時に構造的に解消する指摘:
- CH-04: 周期を考慮した最近傍と、一意の割付け
- CH-06: 線形補間ブラケットの削除と、障壁の事前判定
- CH-07: ヒステリシス付きの結合判定と、min_change 0.2 Å
- CH-08: gmax による判定と、tol 1e-5
- CH-09: NEB の削除と、IDPP の前の写像付き整列
- CH-11: Hessian の型付きの受け渡し
- CH-12: fine / 1e-7
- CH-13: SiteLock、コア数のセマフォ、環境変数
- CH-14: autoz
- CH-17: 3 段の方針と mode-follow
- CH-18 / CH-19: CREST のオプションと配置
- CH-20: xTB の収束判定
- CH-28: エネルギー窓
- CH-31: worker と cwd
- CH-33: prominence の下限と、崩壊時の扱い
- CH-34: 重なりを合否に使わない
- CH-35: 縮退反応
- CH-37: timeout = min(設定値, 残り時間)
- CH-40: 混合 LOT の検出
- CH-41: 標準状態
- CH-42 / CH-45: GoodVibes を API で呼び、キャッシュキーを明示する
- CH-46: 感度の幅
- AR-05: 判断と実行の接続
- AR-07: RunLayout と view
- AR-08: JobStore と継続実行
- AR-16: dummy の除去
- AR-20: thermo の分割
- AR-21 / AR-23 / AR-24: Evidence とゲート
- AR-27: tier と blockers

---

## 13. 批評への対応

**採用**

- 次の項目は、本文の仕様として採用した。
  - chem: 1〜14、16〜19、21〜24
  - eng: P0-1〜10、P1-11〜17、P1-18・19・21〜23、P2-25〜27
- repo の事実を確認した結果:
  - テストはベースラインで 270 passed / 2 skipped。
  - `hfauto/workflow/reaction_state.py:23` に `orca_nebts` がある。
  - QA venv に scipy はないが、`uv pip install` で入る(ネットワークは使える)。
  - WSL の hfauto-prod には goodvibes 4.3.0 の `api.compute_thermo(qcdata=QCData)`、scipy 1.18、pymsym がある。
  - pysisyphus 1.0 の COS は `valid_coord_types = ("cart","cartesian","dlc")`。
  - GoodVibes に同梱された Truhlar DB には PBE0/def2 の項目がなく、`PBE0/MG3S = 0.989 / 0.975` だけがある。

**修正して採用**

| 項目 | 本書での扱い | 理由 |
|---|---|---|
| eng 20(ZPE の許容差は実測してから決める) | 1e-5 Eh に締める | API 経由で自前の振動数を渡すので、GoodVibes の ZPE と参照値の差は本来ほぼ 0 になる。1e-5 でも二重計上(約 2 倍)は確実に検出できる。G06 と smoke で確認する |
| chem 20(xTB の初期 Hessian は A/B 計測の後で既定にする) | 既定で有効にし、Wave 7 で A/B を記録する | 最終判定は常に DFT の freq で行うので、初期 Hessian は結論を変えない。悪化した場合は本書を改訂して既定から外す |
| chem 5(GoodVibes の CLI を使い続ける場合の注意) | API を採用したので不要 | API の採用で、テキストのパースや `--csv` などの問題が構造的に消える |
| eng 16(golden の容量) | 目標 1.2 MB、上限 2 MB。決定的に抜粋する | 抜粋の規則は §10.2 のとおり |

**不採用**

| 項目 | 理由 |
|---|---|
| eng 24(自前の IDPP を pysisyphus の IDPP に置き換える) | Windows の QA 環境でも純関数としてテストでき、Protocol を増やさずに約 120 行で済むため、自前実装を残す |
| chem 15(soft な 2 本目の虚振動で saddle をやり直す) | 反応座標と直交する soft モードが ΔE‡ に与える影響は小さい。TMA 系では saddle と freq を 1 回やり直すだけで数時間かかる。そのため注記 `soft_secondary_mode` を付けるだけにする(eng 18 の案を採用) |
