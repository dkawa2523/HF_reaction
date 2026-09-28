# hfauto

気相の分子と非共有結合錯体について、反応の発見から熱化学と順位付けまでを自動で行うワークフロー。外部プログラム(NWChem、xTB、CREST、SCINE ReaDuct、pysisyphus)の結果は型付きの証拠(`Evidence`)として受け取り、合否は `hfauto/chemistry/gates.py` のゲート関数だけが決める。熱化学は GoodVibes を同じプロセスで呼ぶ。結果は stage ごとの manifest に残り、同じ run ディレクトリで再実行すると済んだジョブは再利用される。

## 対象範囲

| 軸 | 対象 | 対象外 |
|---|---|---|
| 元素 | Z=1〜57、72〜86。Z>36 は def2 系の基底だけで、def2-ECP を元素ごとに自動で書く。計算できる範囲で、遷移金属は未検証 | Ce〜Lu、Z>86、Z>36 に def2 以外の基底(入口で拒否) |
| 電荷 | 任意の整数(陰イオンは SVPD の拡散関数で扱う) | 気相の多価陰イオンは基底に依存しうる |
| スピン | 閉殻一重項(RKS)と高スピン(UKS、⟨S²⟩ を検査)。多重度は利用者が宣言し、未宣言なら SMILES は不対電子数 + 1、xyz は 1、組成はスピン結合で 1 つに決まるときだけその値。d ブロック元素を含む化学種は宣言が必須 | 開殻一重項(ビラジカル、ラジカル対)、MECP・スピン交差、スピン軌道補正(原子で Cl 0.84、Br 3.5、I 7.3 kcal/mol がそのまま誤差になる) |
| 反応型 | 異性化、H・プロトン移動(リレーを含む)、SN2 型置換、H 引き抜き、付加・会合、開殻・電荷系の結合切断、縮退転位、宣言したねじれ。単原子も通常の経路を通る | 溶液相、速度論(マスター方程式、トンネル補正) |
| 系のサイズ | 4 コア・def2-SVPD で約 35 原子まで(解析 freq が N·nbf^2.8 に比例し、15 原子で約 11 分、40 原子で約 6 h) | freq の途中継続(ないので timeout 4 h で全損する。大きな系は ranks を増やす) |

状態は組成と結合グラフ(r < r_cov,i + r_cov,j + 0.4 Å)の状態ラベルで区別する。化学の規則と精度の目安は [docs/design.md](docs/design.md) §7。

## 処理の流れ

8 種の stage を pipeline YAML でつなぐ。`minima` は xTB(screen)と DFT(dft)の 2 回使い、`sp` は `method_panel` の追記で使う。

```text
structures → conformers → minima(screen) → explore → minima(dft) → reaction-paths → thermo → report
```

| stage | 役割 |
|---|---|
| `structures` | system の化学種(xyz / SMILES)を読み、元素・電荷・多重度と宣言反応の原子順序を検査する |
| `conformers` | 単量体の CREST 配座探索と、組成(錯体)の配置 seed と `--nci` 探索 |
| `minima` | opt → 別ジョブの freq → 虚振動に沿った mode-follow → 極小の登録 |
| `explore` | 元素に依らない結合変化の列挙で反応 trial を作り、ReaDuct の NT2(極大がなければ AFIR)で生成物を探す |
| `reaction-paths` | 反応仮説ごとに、経路の分類 → 鞍点 → TS の振動数検証 → QRC による接続 |
| `sp` | 順位に使う点だけの一点計算(エネルギー層と手法パネル) |
| `thermo` | qRRHO(GoodVibes 4.3.0)、キラリティ、会合量、順位の量 δG_eff |
| `report` | ranking.csv、coverage.csv、method_panel.csv、report.html |

| pipeline | stage の並び | 用途 |
|---|---|---|
| `discover` | structures → conformers → screen → explore → dft → paths → thermo → report | 単量体と組成から反応を探す |
| `known_endpoints` | structures → dft → paths → thermo → report | 宣言した反応の端点から経路と熱化学を求める |
| `method_panel` | panel_sp → panel_thermo → panel_report | 既存の run に `--run-dir` で追記し、PBE0/def2-TZVPD と ωB97X-D3/def2-TZVPD で比べ、ωB97X-D3 のエネルギー層で並べ直す |

## 結果の読み方

- 順位の量は δG_eff = max(G_TS, G_R, G_P) − G_R(kcal/mol)。G_R・G_P は反応物・生成物と同じ状態の DFT 極小の最小 G。障壁なし(`barrierless_at_resolution`)と ZPE で沈む障壁(注記 `submerged_barrier`)は max(ΔG_rxn, 0) で同じ表に並ぶ。qRRHO の扱いによる感度の幅が重なる反応は同順位。
- 既定の順位は PBE0-D3BJ/def2-SVPD の値の序数として読む。PBE0 の障壁の誤差は反応クラスで違うので、クラスが混ざるとき数 kcal/mol 未満の差の順序は入れ替わり得る。精度が要るときは `method_panel` を追記する(ranking.csv の `energy_level` 列に並べた層が出て、手法ごとの ΔE‡ の最小・最大が列に付く)。
- 順位を付けない反応も outcome と blockers 付きで ranking.csv に載る。判断の経過は `hfauto case` で見る。

## インストール

本番は WSL2 の Ubuntu(Python 3.12)に NWChem 7.2.3、xTB 6.7.1、CREST 3.0.2 を入れて動かす。手順と版数は [docs/environment.md](docs/environment.md) にある。

## CLI

| コマンド | 役割 |
|---|---|
| `hfauto doctor [--site SITE]` | site の全エンジン(実行ファイル、版数の pin、worker のモジュール、scratch)を検査する。問題があれば終了コード 1 |
| `hfauto run PIPELINE --system SYSTEM --site SITE [--run-dir DIR] [--from ID] [--to ID] [--dry-run] [--retry-failed KINDS]` | パイプラインを実行する。同じ run ディレクトリでは続きから再開する |
| `hfauto status RUN_DIR` | stage の状態、FailureKind 別の失敗数、ジョブの再利用率 |
| `hfauto report RUN_DIR` | run 全体の report.html を書く |
| `hfauto case RUN_DIR REACTION_ID` | 反応ケースの判断ログ(log.jsonl)を表示する |

- PIPELINE・SYSTEM・SITE には YAML のパスか、カレントディレクトリの `configs/<kind>/` の下の名前を指定する。method は pipeline ファイルの隣の `methods/<id>.yaml`、なければ `configs/methods/<id>.yaml` を読む。
- `run` の終了コード: 0 は成功、1 は failed の artifact か失敗した stage がある、2 は設定か preflight の問題、SIGTERM・SIGHUP では外部プログラムを止めて 128 + signum。failed の artifact は後続の stage を止めず report も出るので、結果は終了コードではなく `hfauto status` と report で判断する。
- `--from ID` はその stage 以降を取り直し、`--retry-failed KINDS`(例 `incomplete_output,timeout`)は JobStore に残った該当の失敗だけを再実行する。

## クイックスタート(HCN → HNC)

```bash
hfauto run known_endpoints --system hcn --site wsl_local --run-dir /home/user/hfauto_runs/hcn
hfauto status /home/user/hfauto_runs/hcn
hfauto report /home/user/hfauto_runs/hcn
```

リポジトリの直下で実行する。`--run-dir` を省くと `runs/<system_id>_<pipeline_id>` になるので、本番の run は ext4 上を明示する。4 vCPU の WSL で約 1.5 分かかる。

## 基準値

WSL(4 vCPU / 11 GB)、PBE0-D3BJ/def2-SVPD、298.15 K・1 atm。全体と所要時間は [docs/validation.md](docs/validation.md) にある。

| 系(system) | pipeline | outcome | δG_eff(kcal/mol) |
|---|---|---|---|
| HCN → HNC(`hcn`) | known_endpoints | elementary_step | 42.18 |
| HONO trans → cis(`hono`) | known_endpoints | elementary_step | 11.82(TS がキラルで m = 2) |
| NH3 の反転(`nh3_inversion`) | known_endpoints | degenerate_rearrangement | 3.84 |
| 水(`water_same_basin`) | known_endpoints | same_basin | — |
| TMA·(HF)₂(`tma_hf2`) | discover | same_basin(プロトン移動の TS はない。全体で約 30 分) | — |

## 文書

- [docs/design.md](docs/design.md): 設計(層と import 契約、処理区分の責務、Evidence とゲート、判断表、化学プロトコル、順位の量、実行基盤、設定)
- [docs/environment.md](docs/environment.md): 版数、インストール、WSL の資源、テスト
- [docs/validation.md](docs/validation.md): WSL の実計算による検証の記録
- [docs/reviews/](docs/reviews/): 設計の根拠になったレビュー(現行の改良は `2026-09-27_platform_review.md`)
