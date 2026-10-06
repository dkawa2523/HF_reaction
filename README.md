# hfauto

気相の分子・分子の組合せについて、反応候補の探索、量子化学による経路検証、熱化学の比較を行う研究用ワークフロー。計算結果・未解決の理由・来歴を保存し、停止後に同じ run ディレクトリから再開できる。

目的は、未知の反応候補を広く発見すること、検証済み候補を条件付きで比較すること、初期系から主要生成物・優勢経路を評価すること、および第三者が理解・運用・拡張できるシンプルな基盤を保つことである。現在は個別反応の探索・検証と熱化学比較までを実装している。初期系全体の主要生成物・収率・時間依存の予測は今後の開発対象で、現在の順位から直接判断しない。

## 最初に読む文書

| 文書 | 内容 |
|---|---|
| [docs/roadmap.md](docs/roadmap.md) | 現在の実装・統合・検証状況と、次に行う作業 |
| [docs/improvement-plan.md](docs/improvement-plan.md) | 目的・責務・削除方針・段階別の受入条件を定めた改良計画 |
| [docs/design.md](docs/design.md) | 現行実装の構造と化学的な判断 |
| [docs/environment.md](docs/environment.md) | インストール、実行資源、開発時の検証手順 |
| [docs/validation.md](docs/validation.md) | コードの版と評価段階を明記した実計算の記録 |
| [docs/reviews/README.md](docs/reviews/README.md) | 過去のレビュー・結果報告の位置付け |

## 現在の対象範囲

気相の小〜中規模の主族系が中心。入力には構造、電荷、多重度を指定する。多重度をSMILESやXYZから自動決定しない。高スピン系と一部の低スピン結合錯体はUKSで扱うが、同じ多重度の別の空間電子状態を完全に識別する機能はない。

| 項目 | 現在の扱い |
|---|---|
| 元素 | Z=1〜57、72〜86。Z>36はdef2系の基底とECP。遷移金属の精度は未検証 |
| 反応 | 異性化、H・プロトン移動、置換、H引き抜き、付加・会合、結合切断、縮退転位。自動列挙の被覆は探索ルールと予算に依存 |
| 規模 | 現行4コア環境では約35原子までが実用上の目安。構造・電子状態により費用は大きく変わる |
| 対象外 | 溶液相、スピン交差、開殻一重項の一般的な取扱い、多参照法、捕獲速度・圧力依存・時間依存の反応網予測 |

状態ラベルは結合グラフと立体を用いる。配座、同じ両端を結ぶ別経路、電子状態に関する残る制限は [design.md](docs/design.md) §10にある。単離原子の電子準位・スピン軌道補正の一部は実装済みで、分子・錯体・TS全般の補正ではない。

## 処理の流れ

```text
structures → conformers → minima(screen) → explore → minima(dft) → reaction-paths → sp → thermo → report
```

| stage | 責務 |
|---|---|
| `structures` | 構造と宣言した化学状態を読む |
| `conformers` | 単量体の配座と錯体の開始配置を生成する |
| `minima` | 最適化・振動計算で極小を検証し、登録する。DFT入口の候補選択も行う |
| `explore` | 結合変更を列挙し、低レベルの反応候補を多世代で探索する |
| `reaction-paths` | 仮説からTSを探索し、振動・両側の接続・電子状態を検証する |
| `sp` | 比較に使う計算点の電子エネルギーを求める |
| `thermo` | 熱補正と反応量を計算する |
| `report` | 順位、被覆、未解決の理由をCSV・HTMLで表示する |

`known_endpoints`は宣言した端点から、`discover`は生成物を指定せず開始する。`method_panel`は既存runに別手法のエネルギー比較を追記する。既定の停留点・振動はPBE0-D3BJ/def2-SVPD、エネルギー層はM06-2X-D3(0)/def2-TZVPD、低レベルはGFN2-xTB。

## 実行と再開

本番環境はWSL2 Ubuntu / Python 3.12。[インストールと版数](docs/environment.md)を確認し、リポジトリ直下で実行する。本番runはWSLのext4上に置く。

```bash
hfauto run known_endpoints --system hcn --site wsl_local --run-dir /home/user/hfauto_runs/hcn
hfauto status /home/user/hfauto_runs/hcn
hfauto report /home/user/hfauto_runs/hcn
```

同じ`run`コマンドが停止後の再開になる。設定・記録形式・計算鍵を変えたときの再利用範囲は [design.md](docs/design.md) §8を確認する。

| コマンド | 用途 |
|---|---|
| `hfauto run PIPELINE --system SYSTEM --site SITE --run-dir DIR` | 実行・再開。`--from` / `--to`で処理区間を指定する |
| `hfauto status DIR` | 完了状態、項目の失敗、再利用状況を確認する |
| `hfauto report DIR` | CSV・HTMLを生成する |
| `hfauto case DIR REACTION_ID` | 個別の反応判断ログを見る |
| `hfauto doctor --site SITE` | site全体を点検する。通常のrunは使用するエンジンだけを点検する |

設定はYAMLのパスか`configs/<kind>/`の名前で指定する。`--dry-run`は実行計画の確認、`--retry-failed KINDS`は保存済みの特定の失敗を再試行する指定である。

終了コードは、**0：全stage終了、1：再実行する仕事あり、2：stageまたは設定の失敗**。0でも化学的な未解決や棄却はあり得る。opt・saddleは最後の構造から継続できるが、freq・spの中断はやり直しになる。

## 結果の読み方

- `ranking.csv`の順位は、そのrunの理論・温度・標準状態における個別反応の`δG_eff`による。現在は出発状態の異なる反応も同じ表に並ぶため、初期系からの主要生成物や収率の順位として使わない。
- `δG_eff`は反応の点の鎖の最大の上りで、素反応の速度定数に一律変換する量ではない。分離・錯体の基準は`reference`列にある。定義と不確かさは [design.md](docs/design.md) §7にある。
- 障壁なしの反応は捕獲律速として別表に出る。高精度の一点計算が欠けた反応や認証に落ちた点には順位を付けず、理由を表示する。
- `coverage.csv`は探索・失敗・未試行を示す。計算を終了したこと、期待値との回帰一致、外部参照に対する精度は別々に確認する。

## 開発を再開する

[現在状況](docs/roadmap.md)から未完了項目を確認し、[改良計画](docs/improvement-plan.md)の依存順で進める。実装・運用・計画・過去の証拠を別々の文書で管理し、進捗の正本はroadmap、今後の設計の正本はimprovement-planに置く。[AGENTS.md](AGENTS.md)に責務と変更時の方針、[environment.md](docs/environment.md) §4に検証コマンドがある。
