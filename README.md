# hfauto

気相の分子と非共有結合錯体(クラスター)について、反応の発見から熱化学までを自動で実行するワークフロー。外部プログラム(NWChem、xTB、CREST、SCINE ReaDuct、pysisyphus、GoodVibes)の結果は型付きの証拠(`Evidence`)として受け取り、合否は `hfauto/chemistry/gates.py` のゲート関数だけが決める。結果は stage ごとの manifest に残り、同じ run ディレクトリで再実行すると、済んだジョブは再利用される。

## 目的と範囲

- 対象: 気相の分子と非共有結合錯体。電荷を持つ系と開殻系(多重度 > 1。<S²> を観測する)も扱う。溶媒は PES の指定(NWChem の COSMO、xTB の ALPB)として扱うだけである。
- 求めるもの: 候補生成物、同一 PES 上で検証した極小、1 次の鞍点(別ジョブの振動数で検証)と QRC による接続の確認、qRRHO の熱化学(会合量を含む)、証拠の階層と不確かさの幅に基づく順位。
- 範囲外: wet-etch の反応器・表面モデル、broken-symmetry(閉殻一重項ジラジカルの TS。検出もしない)、速度論、公開 DB による同定。
- fail-closed: ダミーエンジンや内部フォールバックの熱化学はない。外部ジョブの失敗は `FailureKind` 付きの failed artifact として残る。

## 処理の流れ

8 種の stage を pipeline YAML でつなぐ。`minima` は xTB(screen)と DFT(dft)の 2 回使う。

```text
structures → conformers → minima(screen) → explore → minima(dft) → reaction-paths → sp(任意) → thermo → report
```

| stage | 役割 |
|---|---|
| `structures` | system ファイルの化学種(xyz / SMILES)を読み、電子状態と宣言反応の原子順序を検査する |
| `conformers` | CREST による配座探索と、錯体の配置 seed(H 結合の円錐、剛体のランダム配置) |
| `minima` | opt → 別ジョブの freq → 虚振動に沿った mode-follow → 極小のレジストリ |
| `explore` | 反応 trial を作り、ReaDuct の NT2 / AFIR で生成物を探す(陰性結果も記録する) |
| `reaction-paths` | 反応仮説ごとに、障壁の事前判定 → saddle → TS の振動数検証 → QRC → 分類 |
| `sp` | 停留点での一点計算(手法パネル。CCSD(T) は小さい閉殻分子で opt-in) |
| `thermo` | GoodVibes 4.3.0 の API による qRRHO、整合ゲート、会合量、感度の幅 |
| `report` | 順位付け、探索の被覆率、手法パネルの表と HTML |

パイプラインは 3 本である。

| pipeline | stage の並び | 用途 |
|---|---|---|
| `configs/pipelines/discover.yaml` | structures → conformers → screen → explore → dft → paths → thermo → report | 単量体と組成から反応を探す |
| `configs/pipelines/known_endpoints.yaml` | structures → dft → paths → thermo → report | system に宣言した反応の端点から、経路と熱化学を求める |
| `configs/pipelines/method_panel.yaml` | panel_sp → panel_report | 既存の run に `--run-dir` で追記し、停留点を PBE0/def2-TZVPD と ωB97X-D3/def2-TZVPD で比べる。CCSD(T)/def2-TZVPD は小さい閉殻分子のときだけ `methods` に足す。追記した後は panel_report の順位が最終 |

## 設定(4 分割)

4 つの層は中身が重ならないので、マージも優先順位もない。手法を変えるときは別の method ファイルを使う。停留点の手法は minima(dft)と reaction-paths の 2 か所に書き、食い違えば読み込み時に拒否される。宣言反応の端点は xyz で与える(SMILES では原子の対応と配座を決められない)。system と pipeline の YAML は未知のキーを実行前に拒否する。温度は `conditions` だけで指定し、thermo の settings のスケール因子は `vib_scale` の 1 つだけである。

| 層 | 置き場所 | 中身 |
|---|---|---|
| site | `configs/sites/` | 実行ファイルの絶対パス、scratch、コア数とメモリ、エンジンごとの版数の pin と実行設定 |
| method | `configs/methods/` | 汎関数・基底・分散補正・grid・SCF 閾値(xTB は GFN と電子温度) |
| system | `configs/systems/`(xyz は `configs/systems/xyz/`) | 化学種、組成、宣言反応。hcn、hono、nh3_inversion、formaldehyde、tma_hf2、amine_hf_panel、water_same_basin |
| pipeline | `configs/pipelines/` | stage の並びと設定、温度と標準状態、ゲート閾値の上書き(`gates:`) |

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
| `hfauto run PIPELINE --system SYSTEM --site SITE [--run-dir DIR] [--from ID] [--to ID] [--dry-run] [--retry-failed KINDS]` | パイプラインを実行する。同じ run ディレクトリでは続きから再開する |
| `hfauto status RUN_DIR` | stage の状態、FailureKind 別の失敗数、ジョブの再利用率を表示する |
| `hfauto report RUN_DIR` | run 全体の HTML レポートを書く |
| `hfauto case RUN_DIR REACTION_ID` | 反応ケースの判断ログを表示する |

PIPELINE・SYSTEM・SITE には、YAML のパスか `configs/<kind>/` の下の名前を指定する(名前はカレントディレクトリの `configs/` で解決する)。

## クイックスタート(HCN → HNC)

```bash
hfauto run known_endpoints --system configs/systems/hcn.yaml --site configs/sites/wsl_local.yaml
```

- run ディレクトリは `--run-dir` で指定する。省略すると `runs/<system_id>_<pipeline_id>`(ここでは runs/hcn_known_endpoints)になる。本番の run は ext4 上に置く。
- 4 vCPU の WSL で約 2 分かかる。結果は `hfauto status <run>` と `hfauto report <run>` で見る。
- 同じコマンドをもう一度実行すると、完了した stage は飛ばされる。`--from paths` で途中の stage から実行し直せる。

## 検証済みのベンチマーク

WSL(4 vCPU / 11 GB)、PBE0-D3BJ/def2-SVPD での実測値である。7 項目の全体と、撤回した旧値は [docs/validation.md](docs/validation.md) にある。順位は PBE0-D3BJ/def2-SVPD の ΔG‡ の序数として読み、約 1.5 kcal/mol 未満の差は手法誤差の範囲と見る(根拠は [docs/design.md](docs/design.md) §7)。

| 系 | pipeline | 結果 | 所要時間 |
|---|---|---|---|
| HCN → HNC | known_endpoints | TS −1128.5i cm⁻¹、ΔG‡ 42.22、ΔG_rxn 12.68 kcal/mol、elementary_step | 1 分 40 秒 |
| HONO trans → cis | known_endpoints | TS −681.5i cm⁻¹、ΔG‡ 12.24 kcal/mol(旧値 24.11 は撤回) | 7 分 22 秒 |
| TMA·(HF)₂ | discover | same_basin(プロトン移動の TS はない)。reaction-paths は 0 秒 | 52 分 33 秒(旧実装は約 30 h) |

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
