# 計算済み例題から反応探索を読む

[index.html](index.html)をブラウザーで開く。ネット接続や元のWSL runがなくても、5例の3D構造、科学図、判断ログ、文献と出所を閲覧できる。フォルダー全体を渡せば第三者も同じ資料を開ける。表示にはWebGL対応のブラウザーを使う。

## 例題と読む順番

| 例題 | 示せる処理 | 保存結果 |
|---|---|---|
| NH₃反転 | 平面入力が鞍点であること、虚振動、両側の接続、同じ状態を結ぶTS | W5 fresh。δG_eff 4.38 kcal/mol |
| NH₃・HFの水素交換 | 生成物を指定しない探索からTS認証、原子対応と熱化学比較 | W4のQM、W5の熱化学再評価。δG_eff 36.65 kcal/mol、錯体基準 |
| TMA・(HF)₂ | 別名で宣言した端点と独立した極小の違い | W4でsame_basin。認証TSと反応の順位なし。W5 freshは未完了 |
| BH₃＋NH₃ | 会合の拘束走査、中点確認、障壁なしの扱い | W3 fresh。ΔG_rxn −17.83 kcal/mol。速度定数なし |
| OH＋CH₄の引き抜き | 未知候補のDFT入口、接続と電子状態、分離基準 | W4のQM、W5の熱化学再評価。δG_eff 10.55 kcal/mol、分離基準 |

最初はNH₃反転で「TSを認証する」処理を読み、NH₃・HFで候補探索まで広げる。TMA・(HF)₂とBH₃＋NH₃を対比すると、同じ極小への合流、認証したTS、障壁なしを別々に判断する理由が分かる。探索クエリの内訳と判断チャートで、候補の多さと認証結果の違いを確認する。

これらは今回プッシュしたチェックポイント以前のW3〜W5の保存結果である。現行の個別経路処理の説明に用い、現在版での全系fresh評価や、初期系全体の主要生成物比の実証とは扱わない。資料作成に伴う新規QM、元runの変更、計算の停止・再開は行っていない。

## 図の読み方

- **自由エネルギーGの図**：保存された熱化学の停留点だけを水平線で表示する。破線は順序を示し、連続したポテンシャル曲線、時間、距離、IRCではない。ゼロは各図に明記する。異なる例題の障壁を共通の初期混合物の順位として比較しない。
- **電子エネルギーEの図**：保存されたIRCまたは拘束走査の節点を描く。点の間は案内線で、スプライン補間しない。GFN2-xTBの低レベルIRCとPBE0の走査を、M06-2X層のGやDFTのTS認証と混ぜない。IRCの符号付き反復数は物理的な距離ではない。GFN2の電子温度は保存条件を参照する。
- **反応物・生成物の基準**：NH₃・HFとOH＋CH₄では反応物錯体と分離単量体のGを別点で表示する。OH＋CH₄の錯体基準ΔG‡は6.26、分離基準は10.55 kcal/mol。ΔG_rxn −15.07 kcal/molは反応前後の錯体間の値である。BH₃＋NH₃のG図のRは単量体の和で、3Dの接近入力配置とは別である。
- **3D構造**：XYZの原子順を保つ。NH₃反転の左右は保存されたmode-follow構造。NH₃・HFのPは現在の`identity.carry`が採用するQRC＋側の対称像を再現したものなので、独立したQM端点と区別する。QRC端点と、同じbasinに登録された熱化学の振動計算構造は別の場合がある。表示する結合線は距離による描画で、結合の化学的認証ではない。孤立単量体の原子番号は各ファイル内の番号である。
- **判断**：DFTでの停留性、一次の鞍点、QRCによる両側の最適化と端点照合、電子状態、比較層の不足を保存ログから追う。DFTのQRCを連続したIRCと表記しない。低レベルのIRCがあっても、TMAの宣言反応に認証したDFT TSがあるとは限らない。
- **探索件数**：図は最終report時点の保存分類を示し、探索とDFT入口選抜を集約する。`coverage.csv`の`attempts`は成功artifactに保存されたcoverage分類の総数で、未試行も含む。S10のNH₃・HF／TMA・HF pilot全体の3175分類は試行済み2997＋未試行178、TMA・(HF)₂の3118分類は試行済み2998＋未試行120である。計算失敗は別集計され、実際のNT2試行はそれぞれ3000件（2997分類＋失敗3件、2998分類＋失敗2件）。S10では低レベルで850辺を得て、DFT入口のエネルギー選抜後に1辺を保持した。最終`products`は保持した探索辺、`negatives`は探索で候補化しなかった分類とDFT入口で非採用となった辺を表す。OH＋CH₄の棒も同じ保存分類の定義で読む。排他的なoutcomeを同じ棒に積み、理由と計算失敗は各表で確認する。
- **不確かさ**：`band_kcal`は熱化学設定の感度幅で、電子エネルギー精度の信頼区間ではない。障壁なしは調べた計算面と節点の分解能内の判断。捕獲速度と生成物分布は未計算である。

条件は気相、298.15 K、標準状態1 atm。停留点・振動はPBE0-D3BJ/def2-SVPD、Gの電子エネルギーはM06-2X-D3(0)/def2-TZVPD。各例題の実行版、計算条件、元ファイルのSHA-256は`data/data.json`の`provenance`に保存した。`species_records`の`m`は光学異性体数で、多重度ではない。

## 配布する図とデータ

- [5例の位置付け](figures/examples_overview.svg)：各例で示せる処理の一覧。
- [反応判断チャート](figures/reaction_decision_chart.svg)：極小、TS、接続、比較の判定順。
- [探索クエリの内訳](figures/exploration_coverage.svg)：coverage分類、試行済み分類、未試行、保持した候補と認証の関係。
- `figures/<example_id>.svg`：編集・拡大できる各例の科学図。図中の3D構造はMatplotlibで保存座標を投影したもの。
- `figures/<example_id>.png`：300 dpiの配布用画像。
- `data/data.json`、`data/xyz/`：図の数値、XYZ、保存ログ、出所。元のjobsやrun全体をコピーせず、必要なデータだけを同梱する。
- [references.json](references.json)：手法論文、公式文書、ライセンス、図の解釈。

参考文献は手法・表示の根拠である。論文の障壁や別反応の値を本計算の数値として使っていない。表示には[3Dmol.js公式文書](https://3dmol.org/doc/)とRego・Koesの[原論文](https://doi.org/10.1093/bioinformatics/btu829)、図には[Matplotlib](https://matplotlib.org/stable/)を用いた。NT2、PBE0、M06-2X、GoodVibes、QRCの一次文献は資料内に掲載している。

3Dmol.js 2.5.5は公式npmパッケージから取得し、`vendor/source.json`に取得元、パッケージのintegrity、表示スクリプトのSHA-256を保存した。BSD-3-Clause等の配布条件を`vendor/LICENSE.3Dmol`に同梱する。通常の閲覧と再生成ではネットからスクリプトを読み込まない。

## 再生成

リポジトリ直下で、WSL prod venvを使う。配布済みデータからの検証・再作図は元runを必要としない。Matplotlibは作図用の依存で、本番の反応処理へ追加する依存ではない。

```bash
PY=/home/user/.venvs/hfauto-prod/bin/python
$PY examples/visualization/extract_examples.py --check
$PY examples/visualization/render_figures.py
$PY examples/visualization/build_viewer.py
```

元runからデータを再抽出する場合だけ次を使う。既存runは読み取り専用で扱い、runのreport再生成やQMを起動しない。

```bash
$PY examples/visualization/extract_examples.py --runs-root /home/user/hfauto_r10
```

`extract_examples.py`は保存事実の抽出と整合確認、`render_figures.py`は静的作図、`build_viewer.py`はHTMLへのデータ埋め込みを担う。`viewer.js`と`viewer.css`は閲覧操作と表示だけを担い、反応の合否や順位を再計算しない。例題を追加するときは抽出する記録、比較基準、構造の出所、対応する判定をセットで定める。
