# hfauto

気相の分子と非共有結合錯体について、反応の発見から熱化学と順位付けまでを自動で行うワークフロー。外部プログラム(NWChem、xTB、CREST、SCINE ReaDuct、pysisyphus)の結果は型付きの証拠(`Evidence`)として受け取り、合否は `hfauto/chemistry/gates.py` のゲート関数だけが決める。熱化学は GoodVibes を同じプロセスで呼ぶ。結果は stage ごとの manifest に残り、同じ run ディレクトリで再実行すると済んだジョブは再利用される。

## 対象範囲

| 軸 | 対象 | 対象外 |
|---|---|---|
| 元素 | Z=1〜57、72〜86。Z>36 は def2 系の基底だけで、def2-ECP を元素ごとに自動で書く。計算できる範囲で、遷移金属は未検証 | Ce〜Lu、Z>86、Z>36 に def2 以外の基底(入口で拒否) |
| 電荷 | 任意の整数(陰イオンは SVPD の拡散関数で扱う) | 気相の多価陰イオンは基底に依存しうる |
| スピン | 閉殻一重項(RKS)と高スピン(UKS、⟨S²⟩ を検査)。低スピン結合の組成(CH3·O2 の二重項など)も UKS で受け付け、その会合は分離した単量体からの緩和スキャンで判定する。化学種の多重度は必ず宣言する(省略すると system の読み込みで拒否。SMILES のラジカル電子数も xyz も多重度を決めない)。組成は宣言値か、スピン結合で 1 つに決まるときだけその値 | 開殻一重項(ラジカル対の一重項は `low_spin_singlet_unsupported` で止める。1 分子のビラジカルは RKS で計算され、不安定性を検出しない)、MECP・スピン交差(`spin_crossing_reaction_unsupported`)、2 本以上の結合が同時にできる会合、会合の速度定数(VRC-TST)、スピン軌道補正(原子で Cl 0.84、Br 3.5、I 7.3 kcal/mol がそのまま誤差になる) |
| 反応型 | 異性化、H・プロトン移動(リレーを含む)、SN2 型置換、H 引き抜き、付加・会合、開殻・電荷系の結合切断、縮退転位、宣言したねじれ。単原子も通常の経路を通る | 溶液相、速度論(マスター方程式、トンネル補正) |
| 系のサイズ | 4 コア・def2-SVPD で約 35 原子まで(解析 freq が N·nbf^2.8 に比例し、15 原子で約 11 分、40 原子で約 6 h) | freq の途中継続(ないので timeout 4 h で全損する。大きな系は ranks を増やす) |

状態は組成と結合グラフ(r < r_cov,i + r_cov,j + 0.4 Å)の状態ラベルで区別する。化学の規則と精度の目安は [docs/design.md](docs/design.md) §7、範囲外(多参照、トンネル、溶媒、ビット単位の再現性など)と既知の制限は同 §10。

## 処理の流れ

8 種の stage を pipeline YAML でつなぐ。`minima` は xTB(screen)と DFT(dft)の 2 回使い、`sp` は既定の pipeline のエネルギー層(M06-2X)と `method_panel` の追記で使う。

```text
structures → conformers → minima(screen) → explore → minima(dft) → reaction-paths → sp → thermo → report
```

| stage | 役割 |
|---|---|
| `structures` | system の化学種(xyz / SMILES)を読み、元素・電荷・多重度と宣言反応の原子順序を検査する |
| `conformers` | 単量体の CREST 配座探索と、組成(錯体)の配置 seed と `--nci` 探索 |
| `minima` | opt → 別ジョブの freq → 虚振動に沿った mode-follow → 極小の登録。反応する組成だけを DFT にかけ、厳密な像と既知の basin はジョブなしで登録する。DFT の入口で低レベルの辺を DFT//xTB の高さで 1 回だけ選ぶ |
| `explore` | 結合グラフの編集(形成 1〜2、切断 0〜2、生成物に Lewis 構造)を状態ごとに列挙し、ReaDuct の NT2 で状態の間の辺を探す。届いた状態から世代を重ね、試行は件数の予算 1 つ(`max_trials` 3000)。screen で失われた状態は TS のない辺になる |
| `reaction-paths` | 反応仮説ごとに、経路の分類 → 鞍点 → TS の振動数検証 → QRC による接続。検証済みの鞍点(mode-follow、分割した親)は探し直さず検証から始める |
| `sp` | 順位に使う点だけの一点計算(エネルギー層と手法パネル) |
| `thermo` | qRRHO(GoodVibes 4.3.0)、キラリティ、会合量、順位の量 δG_eff |
| `report` | ranking.csv、coverage.csv、method_panel.csv、report.html |

| pipeline | stage の並び | 用途 |
|---|---|---|
| `discover` | structures → conformers → screen → explore → dft → paths → sp → thermo → report | 単量体と組成から反応を探す |
| `known_endpoints` | structures → dft → paths → sp → thermo → report | 宣言した反応の端点から経路と熱化学を求める |
| `method_panel` | panel_sp → panel_report | 既存の run に `--run-dir` で追記し、PBE0/def2-TZVPD と ωB97X-D3/def2-TZVPD の ΔE を並べて示す(順位は run のエネルギー層のまま) |

## 結果の読み方

- 順位の量は δG_eff(kcal/mol)。反応は共通のゼロを持つ点の鎖 [分離した単量体?, R, TS, P, 分離した単量体?] で、δG_eff = max_k(G_k − それより前の最も低い井戸の G)。単分子の素過程では max(G_TS, G_R, G_P) − G_R で、G_R・G_P は反応物・生成物と同じ状態の DFT 極小の最小 G。分離した点は、端の断片が system で宣言した組成の単量体と一致するときだけ持ち、ゼロ(錯体か分離)は ranking.csv の `reference` に出る。ZPE で沈む障壁(注記 `submerged_barrier`)は TS を外した鎖の値。障壁なし(`barrierless_at_resolution`)は序数の順位に入れず、report.html の「Barrierless steps」(捕獲律速)に ΔG_rxn を出す。qRRHO の扱いによる感度の幅が重なる反応は同順位。
- 既定の順位は M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD の G の序数として読む(ranking.csv の `energy_level` は `m06-2x-d3zero/def2-tzvpd`)。CCSD(T)/def2-TZVPD との差は障壁で平均 0.8、最大 1.8 kcal/mol(5 反応)、ΔE_rxn では平均 1.2、最大 2.3 で停留点レベルより良くない。1 kcal/mol 未満の差は序数を保証しない。鞍点の有無と PES の形は PBE0 で決まる。層の SP がない反応は `energy_layer_missing` で順位から外れ、PBE0 で代用しない。遷移金属・溶媒・多参照性の強い系は未検証。詳しくは [docs/design.md](docs/design.md) §7.3・§7.4。
- 順位を付けない反応も outcome と blockers 付きで ranking.csv に載る。多段反応(`multi_step`)の親は順位を持たず、分割した子反応(`<id>_split<n>`)がそれぞれ並ぶ。判断の経過は `hfauto case` で見る。駆動しなかった子(同じ問いの結論を写した `same_as:<id>` と、分割の深さの上限の `split_depth`)は log を持たないので、`hfauto case` は `no case log` で終わる。

## インストール

本番は WSL2 の Ubuntu(Python 3.12)に NWChem 7.2.3、xTB 6.7.1、CREST 3.0.2 を入れて動かす。手順と版数は [docs/environment.md](docs/environment.md) にある。

## CLI

| コマンド | 役割 |
|---|---|
| `hfauto doctor [--site SITE]` | site の全エンジン(実行ファイル、版数の pin、worker のモジュール、scratch)を検査する。問題があれば終了コード 1 |
| `hfauto run PIPELINE --system SYSTEM --site SITE [--run-dir DIR] [--from ID] [--to ID] [--dry-run] [--retry-failed KINDS]` | パイプラインを実行する。同じ run ディレクトリでは続きから再開する |
| `hfauto status RUN_DIR` | stage の状態、FailureKind 別の失敗数、ジョブの再利用率。終了コードは `run` と同じ規則(1 = failed・incomplete の stage か failed の artifact) |
| `hfauto report RUN_DIR` | run 全体の report.html を書く |
| `hfauto case RUN_DIR REACTION_ID` | 反応ケースの判断ログ(log.jsonl)を表示する |

- PIPELINE・SYSTEM・SITE には YAML のパスか、カレントディレクトリの `configs/<kind>/` の下の名前を指定する。method は pipeline ファイルの隣の `methods/<id>.yaml`、なければ `configs/methods/<id>.yaml` を読む。
- `run` の終了コード: 0 は成功、1 は failed の artifact か失敗した stage がある、または SIGINT・SIGTERM・SIGHUP で止めた(外部プログラムを止め、実行中の stage を incomplete にする。同じコマンドを流し直せば JobStore から続ける)、2 は設定か preflight の問題。failed の artifact は後続の stage を止めず report も出るので、結果は終了コードではなく `hfauto status` と report で判断する。
- `--from ID` はその stage 以降を取り直し、`--retry-failed KINDS`(例 `incomplete_output,nonzero_exit`)は JobStore に記録された該当の種類の失敗を再実行する。timeout は記録しないので、取り直せば常に再実行される。

## クイックスタート(HCN → HNC)

```bash
hfauto run known_endpoints --system hcn --site wsl_local --run-dir $HOME/hfauto_runs/hcn
hfauto status $HOME/hfauto_runs/hcn
hfauto report $HOME/hfauto_runs/hcn
```

リポジトリの直下で実行する。`--run-dir` を省くと `runs/<system_id>_<pipeline_id>` になるので、本番の run は ext4 上を明示する。4 vCPU の WSL で約 1 分 10 秒かかる(QM 242 core 秒。うち M06-2X の層の SP 3 点が 35 core 秒)。記録の型が変わる前の run dir は途中から再開できない(design.md §8)。

## 基準値

WSL(4 vCPU / 11 GB)、M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD、298.15 K・1 atm。全体と所要時間は [docs/validation.md](docs/validation.md) にある。

| 系(system) | pipeline | outcome | δG_eff(kcal/mol) |
|---|---|---|---|
| HCN → HNC(`hcn`) | known_endpoints | elementary_step | 41.57 |
| HONO trans → cis(`hono`) | known_endpoints | elementary_step | 9.88(TS がキラルで m = 2) |
| NH3 の反転(`nh3_inversion`) | known_endpoints | degenerate_rearrangement | 4.38 |
| 水(`water_same_basin`) | known_endpoints | same_basin | — |
| TMA·(HF)₂(`tma_hf2`) | discover | same_basin(プロトン移動の TS はない) | — |

## 文書

- [docs/design.md](docs/design.md): 設計(層と import 契約、処理区分の責務、Evidence とゲート、判断表、化学プロトコル、順位の量、実行基盤、設定、範囲外)
- [docs/environment.md](docs/environment.md): 版数、インストール、WSL の資源、テストと品質ゲート
- [docs/validation.md](docs/validation.md): WSL の実計算による全体の再検証(VAL9)、期待値とその変更の理由、到達表
- [docs/roadmap.md](docs/roadmap.md): round 10 の改良の計画(提案と wave、計画の変更の記録、見送り)。検証の case と道具は `validation/`(cases.yaml、check.py、replay.py、BH76 の参照)
- [docs/reviews/](docs/reviews/): 設計の根拠になったレビューと改良の結果報告(最新は [2026-10-02_round9_result.md](docs/reviews/2026-10-02_round9_result.md)、その根拠の分析は [2026-09-30_remaining_issues_analysis.md](docs/reviews/2026-09-30_remaining_issues_analysis.md))
