# hfauto

気相の分子と非共有結合錯体(クラスター)について、反応の発見から熱化学までを自動で実行するワークフロー。外部プログラム(NWChem、xTB、CREST、SCINE ReaDuct、pysisyphus、GoodVibes)の結果は型付きの証拠(`Evidence`)として受け取り、合否は `hfauto/chemistry/gates.py` のゲート関数だけが決める。結果は stage ごとの manifest に残り、同じ run ディレクトリで再実行すると、済んだジョブは再利用される。

## 目的と範囲

- 対象: 気相の分子と非共有結合錯体。電荷を持つ系と開殻系(多重度 > 1。<S²> を観測する)も扱う。溶媒は PES の指定(NWChem の COSMO、xTB の ALPB)として扱うだけである。
- 求めるもの: 候補生成物、同一 PES 上で検証した極小、1 次の鞍点(別ジョブの振動数で検証)と QRC による接続の確認、qRRHO の熱化学(会合量を含む)、証拠の階層と不確かさの幅に基づく順位。
- 対応元素: Z=1〜57、72〜86(計算できる範囲で、遷移金属は未検証)。Z>36 は def2 系の基底だけで扱い、def2-ECP を自動で書く。単原子の化学種も通常の経路を通る(explore の出発点にはならない)。
- 状態: 組成と結合グラフ(r < Σr_cov + 0.4 Å)の状態ラベルで区別する。イオン–双極子錯体は 1 断片になりうる。
- スピン状態: 化学種と組成の多重度は利用者が宣言する。未宣言なら SMILES は不対電子数 + 1(原子は高スピン仮定)、xyz は 1、組成はスピン結合で値が 1 つに決まるときだけその値。d ブロック元素を含む化学種は宣言が必須。
- 範囲外: wet-etch の反応器・表面モデル、開殻一重項(ビラジカル、ラジカル対。BS-UKS は使わず、ラジカル再結合の生成物は単量体として宣言する)と MECP・スピン交差、速度論、公開 DB による同定。
- fail-closed: ダミーエンジンや内部フォールバックの熱化学はない。外部ジョブの失敗は `FailureKind` 付きの failed artifact として残る。

## 処理の流れ

8 種の stage を pipeline YAML でつなぐ。`minima` は xTB(screen)と DFT(dft)の 2 回使う。

```text
structures → conformers → minima(screen) → explore → minima(dft) → reaction-paths → sp(任意) → thermo → report
```

| stage | 役割 |
|---|---|
| `structures` | system ファイルの化学種(xyz / SMILES)を読み、電子状態と宣言反応の原子順序を検査する |
| `conformers` | CREST による配座探索と、錯体の配置 seed(H 結合の円錐、剛体のランダム配置。CREST が失敗した組成(開殻など)は seed を出す) |
| `minima` | opt → 別ジョブの freq → 虚振動に沿った mode-follow → 極小のレジストリ |
| `explore` | 元素に依らない結合変化のテンプレート(移動 / リレー / 形成 / 切断)で反応 trial を作り、ReaDuct の NT2(極大がなければ AFIR)で生成物を探す(陰性結果も記録する) |
| `reaction-paths` | 反応仮説ごとに、障壁の事前判定 → saddle → TS の振動数検証 → QRC → 分類 |
| `sp` | 順位に使う点だけの一点計算(エネルギー層と手法パネル。CCSD(T) は小さい分子で opt-in、開殻は ROHF-CCSD(T)) |
| `thermo` | GoodVibes 4.3.0 の API による qRRHO、キラリティ(m = 2)、整合ゲート、会合量、順位の量 δG_eff と感度の幅 |
| `report` | δG_eff による順位付け(ranking.csv)、探索の被覆率、手法パネルの表と HTML |

パイプラインは 3 本である。

| pipeline | stage の並び | 用途 |
|---|---|---|
| `configs/pipelines/discover.yaml` | structures → conformers → screen → explore → dft → paths → thermo → report | 単量体と組成から反応を探す |
| `configs/pipelines/known_endpoints.yaml` | structures → dft → paths → thermo → report | system に宣言した反応の端点から、経路と熱化学を求める |
| `configs/pipelines/method_panel.yaml` | panel_sp → panel_thermo → panel_report | 既存の run に `--run-dir` で追記し、停留点を PBE0/def2-TZVPD と ωB97X-D3/def2-TZVPD で比べ、参照手法(`energy_method`)のエネルギー層で順位を付け直す(後から追記した層が元の thermo の値を上書きし、ranking.csv の `energy_level` 列に層が出る)。CCSD(T)/def2-TZVPD は小さい分子(重原子約 4 個まで)のときだけ `methods` に足し、`energy_method` にする |

## 設定(4 分割)

4 つの層は中身が重ならないので、マージも優先順位もない。手法を変えるときは別の method ファイルを使う。停留点の手法は minima(dft)と reaction-paths の 2 か所に書き、食い違えば読み込み時に拒否される。宣言反応の端点は xyz で与える(SMILES では原子の対応と配座を決められない)。system と pipeline の YAML は未知のキーを実行前に拒否する。温度と標準状態は thermo stage の `temperatures_K`・`standard_states`(既定 298.15 K・1 atm)だけで指定し、report は ReportConfig の指定がなければ thermo の最初の (T, 標準状態) で並べる。thermo の settings のスケール因子は `vib_scale` の 1 つだけである。

| 層 | 置き場所 | 中身 |
|---|---|---|
| site | `configs/sites/` | 実行ファイルの絶対パス、scratch、コア数とメモリ、エンジンごとの版数の pin と実行設定 |
| method | `configs/methods/` | 汎関数・基底・分散補正・grid・SCF 閾値(xTB は GFN と電子温度) |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | 化学種、組成、宣言反応。hcn、hono、nh3_inversion、formaldehyde、tma_hf2、amine_hf_panel、water_same_basin と、レビュー §6 の検証セット(sn2_cl ほか 17 系。対応は docs/validation.md の r6) |
| pipeline | `configs/pipelines/` | stage の並びと設定、ゲート閾値の上書き(`gates:`) |

## インストール

本番は WSL2 の Ubuntu(Python 3.12)で、NWChem 7.2.3、xTB 6.7、CREST 3.0 を入れた環境で動かす。production extra(pysisyphus、SCINE ReaDuct、GoodVibes など)は Linux 専用である。版数と WSL の設定は [docs/environment.md](docs/environment.md) にある。リポジトリの直下で次を実行する。

```bash
uv venv --python 3.12 /home/user/.venvs/hfauto-prod
uv pip install --python /home/user/.venvs/hfauto-prod/bin/python -e ".[production,dev]"
hfauto doctor --site configs/sites/wsl_local.yaml
```

`hfauto doctor` は、実行ファイル、版数の pin、worker 側の Python モジュール、scratch の場所(ext4 上であること)を検査する。

## CLI

| コマンド | 役割 |
|---|---|
| `hfauto doctor [--site SITE]` | site の全エンジンを検査する(問題があれば終了コード 1) |
| `hfauto run PIPELINE --system SYSTEM --site SITE [--run-dir DIR] [--from ID] [--to ID] [--dry-run] [--retry-failed KINDS]` | パイプラインを実行する。同じ run ディレクトリでは続きから再開する。stage の失敗か failed の artifact があれば終了コード 1 |
| `hfauto status RUN_DIR` | stage の状態、FailureKind 別の失敗数、ジョブの再利用率を表示する |
| `hfauto report RUN_DIR` | run 全体の HTML レポートを書く |
| `hfauto case RUN_DIR REACTION_ID` | 反応ケースの判断ログを表示する |

PIPELINE・SYSTEM・SITE には、YAML のパスか `configs/<kind>/` の下の名前を指定する(名前はカレントディレクトリの `configs/` で解決する)。pipeline が使う method は、pipeline ファイルの隣の `methods/<id>.yaml` があればそれを、なければカレントディレクトリの `configs/methods/<id>.yaml` を読む。

## クイックスタート(HCN → HNC)

```bash
hfauto run known_endpoints --system configs/systems/hcn.yaml --site configs/sites/wsl_local.yaml
```

- run ディレクトリは `--run-dir` で指定する。省略すると `runs/<system_id>_<pipeline_id>`(ここでは runs/hcn_known_endpoints)になる。本番の run は ext4 上に置く。
- 4 vCPU の WSL で約 2 分かかる。結果は `hfauto status <run>` と `hfauto report <run>` で見る。
- 同じコマンドをもう一度実行すると、完了した stage は飛ばされる。`--from paths` で途中の stage から実行し直せる。

## 検証済みのベンチマーク

WSL(4 vCPU / 11 GB)、PBE0-D3BJ/def2-SVPD での実測値(第2ラウンド改良後の再検証 v3)である。9 項目の全体、改良前(W7)・v2 との比較と、撤回した旧値は [docs/validation.md](docs/validation.md) にある。順位の量は δG_eff = max(G_TS, G_R, G_P) − G_R で、G_R・G_P は反応物・生成物と同じ状態の最小 G である。障壁なし(barrierless_at_resolution)と、ZPE で障壁が沈む反応(順方向か逆方向の ΔE‡ + ΔZPE‡ ≤ 0。注記 `submerged_barrier`)は max(ΔG_rxn, 0) で同じ表に並ぶ。ranking.csv の `torsional` 列は結合変化のないねじれの段を示し、手法パネルを追記すると ΔE‡ の最小・最大の列が付く。順位は PBE0-D3BJ/def2-SVPD の δG_eff の序数として読む。PBE0 の障壁の誤差は反応クラスで違い(文献の平均符号付き誤差は H 移動 −4.2、重原子移動 −6.6、求核置換 −1.9 kcal/mol。CCSD(T) との差の実測は HCN −1.2、HONO +2.0、SN2 −3.0)、クラスが混ざるときは差が数 kcal/mol 未満の順序は入れ替わり得る。精度が要るときは method_panel を追記してエネルギー層で並べ直す(根拠は [docs/design.md](docs/design.md) §7)。

| 系 | pipeline | 結果 | 所要時間 |
|---|---|---|---|
| HCN → HNC | known_endpoints | TS −1128.5i cm⁻¹、ΔG‡ 42.18、ΔG_rxn 12.68 kcal/mol、elementary_step | 1 分 24 秒 |
| HONO trans → cis | known_endpoints | TS −681.0i cm⁻¹、ΔG‡ 11.82 kcal/mol(TS がキラルで m = 2。旧値 24.11 は撤回) | 5 分 23 秒 |
| TMA·(HF)₂ | discover | same_basin(プロトン移動の TS はない)。reaction-paths は 0 秒 | 30 分 58 秒(旧実装は約 30 h) |

## 開発

```bash
.venv-win-qa/Scripts/python.exe -m pytest -q -p no:cacheprovider
.venv-win-qa/Scripts/ruff check hfauto tests
.venv-win-qa/Scripts/lint-imports
```

Windows の QA 用 venv では、fake エンジンと golden のテストだけが動く。実エンジンの smoke(`pytest -m real`)は WSL で実行する([docs/environment.md](docs/environment.md))。

## 文書

- [docs/design.md](docs/design.md): 現状の設計(層と import 契約、Evidence とゲート、payload、能力 Protocol、判断表、化学プロトコルの既定値、JobStore)
- [docs/environment.md](docs/environment.md): 版数の固定、WSL の資源、テストの実行方法
- [docs/validation.md](docs/validation.md): WSL の実計算による検証の記録
- [docs/reviews/2026-09-25_architecture_chemistry_review.md](docs/reviews/2026-09-25_architecture_chemistry_review.md): リファクタリングの根拠になったレビュー
- [docs/reviews/2026-09-25_refactor_design.md](docs/reviews/2026-09-25_refactor_design.md): 計画時の設計書(凍結)
