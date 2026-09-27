# hfauto 計算基盤レビュー — 他の分子種・原子種へ汎用的で、化学的に妥当で、有用かつ効率的な反応計算基盤にするための課題と抜本的改良案（2026-09-27）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `07521ed`（v4 の round-4 addendum まで反映済み）。本文は読み取りだけのレビューで、リポジトリのコード・設定・既存文書は変更していません（追加したのはこのファイルだけ）。
- ベースライン: `docs/reviews/2026-09-26_calculation_spec_review_v4.md`（round-4 addendum を含む。以下「v4」）。v4 までの ID（C8、C12 など）は、ここでも同じ意味で参照します。
- 方法（3 段）:
  1. **実計算プローブ**: これまで実計算で一度も通っていなかった経路（手法パネル、電荷のある系、開殻系、単原子、障壁の小さい経路、負モードが 2 本以上ある鞍点、Z>36 の元素、計算時間）を、WSL（4 vCPU、NWChem 7.2.3、xTB 6.7.1、CREST 3.0.2、ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0）で 16 本走らせました。入力・run・ログは `/home/user/hfauto_probe_r5/` にあります（§2）。
  2. **処理区分ごとの精査**: U0〜U9 の 10 区分について、コード、過去の実 run（`/home/user/hfauto_w7`、`hfauto_v2`〜`hfauto_v5`）、各ライブラリのマニュアルとソース（NWChem の `opt_drv.F`・`bas_input.F`・`hess_check.F`、CREST の `confparse.f90`・`cregen.f90`・`setuptest.f90`、GoodVibes の `thermo.py`、pymsym、ReaDuct のマニュアル）、反応理論の文献を照合し、課題と改良案を作りました。
  3. **反証的な検証**: 別の担当が 2 の主張をコードと生データで確かめ、課題ごとに confirmed / partially / refuted、改良案ごとに adopt / adopt_modified / reject を判定しました。
- 記載の規則:
  - 課題の重大度と評価記号は、検証段で訂正した値を使いました。全面的に反証された主張は載せていません。一部だけ正しい主張は、訂正後の内容で書きました。
  - 改良案は「採用」「修正して採用」のものだけを、修正後の形で載せました。却下した案と、評価段で見送った案は、各節の「見送り」と §4.4 に理由付きで載せました。
  - 複数の区分で同じ内容になった課題・改良案は、1 つの**担当区分**にまとめ、他の区分からは「→ 担当 ID」で参照します（統合表は §4.1）。
  - 行番号は HEAD `07521ed` のものです。

凡例:

| 記号 | 意味 |
|---|---|
| critical / high / medium / low | 課題の重大度（critical は「黙って誤った数値を出す」もの） |
| must / should / could | 優先度。must は「結論を誤らせる」「主要な化学で run が止まる」「主要機能が実質働かない」もの |
| 【複雑さ減】/【中立】/【複雑さ増】 | コード・設定・分岐・knob への影響 |
| ◎ / ○ / △ / × | 評価。◎ 問題なし、○ 軽微な課題のみ、△ 改良が必要、× 目的どおりに機能していない |

---

## 0. 要旨

### 0.1 計算基盤としての総合評価

**骨格は健全だが、汎用性は「中性・閉殻・第 2 周期・単分子異性化」の外で急に崩れる。** 今の hfauto は「決められた小分子の単分子異性化を、再現性よく、DFT 停留点レベルで確かめる基盤」としてはよく働いています（HCN・HONO・NH3・水の回帰は v4→v5 で不変）。しかし今回のプローブで未通過経路を通すと、16 本中 9 本で、黙った誤り・run 全体の停止・正しい TS を捨てる結論のいずれかが起きました。

| 観点 | 評価 | 一言 |
|---|---|---|
| 化学的妥当性 | △ | 停留点・LOT 照合・qRRHO は正しい。Z>36 の ECP 欠落（黙った誤り）、縮退判定が入力座標依存、TS/QRC ゲートが正しい TS を捨てる、順位の量が未定義 |
| 汎用性（元素・電荷・スピン・反応型） | × | 単原子で run が落ちる、元素ごとに失敗の場所が違う、組成の多重度は最高スピン固定、探索の drive が極性 H と C–H 移動に偏る |
| 有用性 | △ | 探索は約 150 attempt で生成物 1 件、手法パネルが順位に届かない、障壁なしが順位から消える |
| 効率 | ○ | 1 ジョブの費用は妥当。読まれない診断の DFT freq、誤った方向による DFT Hessian と鞍点探索の空回り、停滞した鞍点の継続が無駄 |
| 簡潔さ | △ | 上書きされない knob 23、同じ判定・定義の重複、閉じた式の外部エンジン化、読まれない記録 |

実装上の土台（内容アドレス型 JobStore、純関数の決定表 `decide`、Evidence と gates の分離、観測した Level の照合、execute_stage を唯一の入口とする pipeline）は、このまま汎用基盤の土台として使えます。**直すべきは土台ではなく、各区分の「化学の規則」の素朴さと、入口の宣言の欠如と、失敗の閉じ込め単位です。**

### 0.2 最重要の課題（8 件）

| # | 課題 | 重大度 | 実証 | 担当 |
|---|---|---|---|---|
| 1 | def2 基底で Z>36 に ECP を書かない。HI は 54 電子の全電子計算で黙って E=−2005.8 Eh、I⁻+CH3I は SCF 未収束で run が落ちる | critical | P6 | U0-I1 |
| 2 | 単原子の化学種（Cl⁻、H•、F⁻、金属イオン）で run 全体が落ち、分離反応物基準も出ない。1 件の parse 例外・1 ケースの例外・上流の全件失敗でも run が止まる | high | P2、P6 | U3-I1、U9-I1 |
| 3 | 縮退反応の判定を入力座標で行い、少しずれた入力の宣言反応を same_basin で黙って捨てる。鏡像を別 basin にし、鏡像化が「反応」として順位 1 位に載る | high | P5、P3 | U5-I8、U3-I2 |
| 4 | 鞍点探索の方向ベクトル（全原子の端点差）が傍観者の運動に支配され、良い xTB Hessian を捨てて DFT Hessian を払い、誤ったモードを追う（acac で 60 分打ち切り） | high | T | U6-I2、U6-I1 |
| 5 | QRC ゲートが、直線変位で TS より上から始まる正しいねじれ型 TS を棄却する。TS ゲートは宣言端点基準の low_prominence と torsional 依存の閾値という TS 固有でない条件を持つ | high | P4 | U6-I5、U6-I4 |
| 6 | 経路の形を 2 つの別定義で判定し、どちらも井戸（中間体）を数えない。障壁なし・中間体の出口は実計算で一度も通っていない | high | — | U5-I1 |
| 7 | 探索の drive が極性 H と C–H 移動に偏り、SN2・H 引き抜き・付加・分子間交換を生成しない。縮退転位を same_as_source として捨てる。実 run の約 150 attempt で生成物は 1 件 | high | 実 run | U4-I1、U4-I5 |
| 8 | 順位の量が定義されていない（1 配座基準、障壁が ZPE で沈んでも数値のまま、障壁なしは順位外）。手法パネルは動くが順位に届かず、代わりに順位の量でない ΔE_rxn の符号で順位を止める | high | P1、P4 | U7-I1、U8-I2、U8-I3 |

### 0.3 抜本的改良の方向性

1. **入口で 1 回だけ宣言し、下流は宣言を信じる。** 元素表を 1 つにして対応範囲を structures で判定する（U1-P1）。ECP は def2 の定義として元素ごとに自動で書く（U0-P1）。スピン状態は利用者が宣言し、一意でなければ宣言を求める（U1-P3）。単原子は特別扱いせず同じ経路に通す（U3-P1）。
2. **失敗は最小の単位で型付きに閉じ込め、run は続ける**（U9-P1）。診断を増やすのではなく、分単位の経路網羅の回帰セット（§6）で検証の空白を塞ぐ。
3. **化学の規則を元素に依らない最小の規則に置き換える。** 探索は「移動・置換／リレー／形成／ラジカル・イオンの切断」の結合変化列挙器 1 つ（U4-P1）。経路の形は山と井戸を同じ関数で数える分類器 1 つ（U5-P1）。反応方向は「低レベル TS の虚モード → 反応中心の接線」の 1 規則（U6-P2）。初期 Hessian は「反応モードだけ負」の 1 関数（U6-P1）。変位は「ν と質量からのエネルギー目標」の 1 規則（U3-P3、U6-P4）。
4. **判定は TS・接続の定義そのものだけにする。** TS は「負の固有値 1 本」、接続は「両側が別々の極小へ下る」（U6-P5）。同一性は basin の最適化構造と置換で決める（U5-P6、U3-P2）。
5. **順位の量を 1 つ定義する**（U7-P1）。`δG_eff = max(G_TS, G_R, G_P) − G_R`、障壁なし・ZPE で沈んだ障壁も同じ量の上で並べる。精度は構造・Hessian（PBE0、解析 Hessian）とエネルギー層（ωB97X-D3/TZVPD、CCSD(T)）の分業で上げ、composite を既定の経路にする（U8-P3）。
6. **削る。** 上書きされない knob 23→8、読まれない記録と、それを作る DFT freq、符号不一致と dzpe の経験則 blocker、GoodVibes の子プロセス・JobStore・実行時照合、溶媒和・d4・MP2 の半実装、自動の配座対仮説（§4.3）。

### 0.4 今のままでよい点

- 内容アドレス型の JobStore（`Task.key_payload`）、型付き `Failure` と上限付き LADDER、`input_sha`・`config_sha` による再開、`execute_stage` を唯一の入口にした pipeline。
- 純関数の決定表 `decide`（`hfauto/drivers/reaction_case/state.py`）、Evidence・gates・Action の分離、ケースごとの `log.jsonl`。
- 停留点を PBE0-D3BJ/def2-SVPD・grid fine・`convergence energy 1e-7` の 1 レベルにそろえ、解析 Hessian の別ジョブ freq で確かめること。観測した Level と要求の照合（`core/method.py:84-94`、`backends/nwchem/output.py:87-111`）と `same_pes`。
- 電荷と多重度の全エンジンへの受け渡し（NWChem の charge/mult/odft、xTB・CREST の `--chrg/--uhf`、pysis、ReaDuct の restricted_open_shell、GoodVibes の S_el）。SN2（電荷 −1）と CH3O•（二重項）で初めて実証しました。
- 探索を GFN2-xTB で行い、DFT の検証を下流に分ける分業。ReaDuct の NT2→Bofill→IRC→両端 opt→結合グラフ照合という標準の流れ。
- 峠の定理に基づく上界の論理（連続経路の最大 ≥ 峠）、エネルギー目標の QRC、composite `G = E_SP + (G_GV − E_GV)`、1 atm からの Δn による標準状態換算、qRRHO の感度幅。

---

## 1. 目標とする計算基盤の姿

### 1.1 処理区分ごとの責務と入出力（改良後）

| 区分 | 責務（1 行） | 入力 | 出力 | 主エンジン |
|---|---|---|---|---|
| U0 プロトコル | どの段をどの理論レベルで計算し、順位をどの面の量で付けるかを決める。LOT をエンジンに正しく渡し観測で照合する（ECP は def2 の定義として自動） | method YAML、pipeline YAML | Level、デッキ | — |
| U1 入力・電子状態 | 化学種と組成を、対応元素・電荷・宣言された多重度・結合グラフ（状態ラベル）を持つ構造に変える。対応範囲外は入口で 1 回だけ拒否する | system YAML（SMILES/xyz、組成、multiplicity） | SpeciesRecord、composition_id、state_label | RDKit |
| U2 配座・錯体配置 | 下流の screen に渡す候補構造を作る（同一性の確定と窓による選抜はしない） | SpeciesRecord、組成 | 候補構造（状態ラベルごと keep_per_state） | CREST |
| U3 極小の確定 | 同じ LOT で確かめた極小を basin に登録する（単原子を含む）。鞍点に落ちた入力は mode-follow で解消し、親の宣言をつなぐ | 候補構造、生成物 | MinimumRecord（basin、members、chiral、state_label） | xTB、NWChem |
| U4 反応探索 | 生成物を仮定しない片端探索で、低レベルの固有な生成物・TS を出す | screen 極小 | DiscoveryRecord（固有な生成物のみ、縮退を含む） | ReaDuct |
| U5 仮説・経路 | 極小対から仮説を作り、DFT 極小を固定端とする経路で「障壁なし／単一の山／中間体」を分類し、種を渡す | MinimumRecord、DiscoveryRecord | ReactionCase、Seed、分類 | pysisyphus、NWChem ZTS |
| U6 鞍点・接続 | 種を一次鞍点に精密化し、別ジョブの freq で TS を確かめ、QRC で両側の極小につなぐ | Seed | SaddleClaim、ConnectionClaim、outcome | NWChem |
| U7 熱化学 | 各停留点の G（qRRHO、σ・m、composite）と、1 つに定義した順位の量 δG_eff を出す | freq Evidence、SP Evidence | SpeciesThermo、ReactionThermo | GoodVibes（同一プロセス） |
| U8 一点計算・報告 | エネルギー層（RSH、CCSD(T)）を順位に使う点にだけ計算し、δG_eff で並べ、パネルの幅を列で示す | 停留点、ReactionThermo | ranking.csv、method_panel.csv、report.html | NWChem |
| U9 実行基盤 | stage を動かし、ジョブを再生・継続し、失敗を最小単位で型付きに閉じ込め、コアと予算を管理する | pipeline、site | manifest、JobStore、終了コード | — |

### 1.2 汎用性の対象範囲（改良後に宣言する範囲）

| 軸 | 対象 | 対象外（文書に明記） |
|---|---|---|
| 元素 | Z=1–57、72–86（GFN2-xTB の Z≤86 と NWChem の def2-SVPD/TZVPD ライブラリの共通部分）。Z>36 は def2-ECP を元素ごとに自動で付ける。「計算できる」範囲で、遷移金属は「検証済み」ではない | Ce〜Lu、Z>86、全電子の相対論計算。Z>36 に def2 以外の基底（入口で拒否） |
| 電荷 | 任意の整数（陰イオンは SVPD の拡散関数で扱う） | 気相の多価陰イオンは基底依存になりうる（注意書き） |
| スピン | 閉殻一重項（RKS）、高スピン・二重項以上（UKS、⟨S²⟩ を検査）。組成の多重度は利用者が宣言し、一意でなければ宣言を求める。d ブロック元素を含む化学種は宣言必須 | 開殻一重項（ビラジカル、ラジカル対の一重項面）、MECP・スピン交差、スピン軌道補正（原子の目安 Cl 0.84、Br 3.5、I 7.3 kcal/mol） |
| 反応型 | 単分子異性化、H・プロトン移動（リレーを含む）、SN2 型の置換、H 引き抜き、付加・会合、ラジカルとイオンの切断（β 開裂）、縮退転位、宣言したねじれ反応。単原子を含む分離反応物基準 | 溶液相（気相専用）、多段の総括速度・マスター方程式、多次元トンネル |
| 系のサイズ | 4 コア・def2-SVPD で約 35 原子まで（freq ∝ N·nbf^2.8、15 原子で約 11 分、40 原子で約 6 h > timeout 4 h）。大きな中性系は def2-SVP の method を選ぶ（U0-P13）。HPC では ranks を増やす | freq の途中継続（ない。timeout で全損） |


---

## 2. 未検証だった経路の実計算プローブ結果

すべて code `07521ed`、PBE0-D3BJ/def2-SVPD（停留点）、WSL 4 コア、逐次実行です。run とログは `/home/user/hfauto_probe_r5/runs/<name>` と `/home/user/hfauto_probe_r5/logs/<name>.log`。

| P | 系（run 名） | 目的 | 結果 | 時間 | 失敗・異常の根本原因 |
|---|---|---|---|---|---|
| P1a | HCN の v5 複製に method_panel 追記（`hcn_panel`） | 手法パネル（C43）の正式な実計算 | 成功。SP 6 件、失敗 0。ΔE‡ は PBE0/SVPD 46.62、PBE0/TZVPD 46.43、ωB97X-D3/TZVPD 46.37 kcal/mol。ranking.csv の ΔG‡ は 42.178（SVPD）のまま | 35 s | `configs/pipelines/method_panel.yaml` に thermo（`energy_method`）がなく、パネルが順位に届かない |
| P1b | Cl⁻···CH3Cl のパネル（CCSD(T) 入り、`sn2_panel`） | 電荷系で手法の幅を見る | 成功。中心障壁 ΔE‡ は PBE0/SVPD 10.45、PBE0/TZVPD 11.15、ωB97X-D3 15.10、CCSD(T)/def2-TZVPD 13.44。ranking は 11.08 のまま | 約 10 分 | 同上。PBE0/SVPD は CCSD(T) より −3.0、ωB97X-D3 は +1.7 |
| P1c | CH3O• のパネル（`ch3o_panel`） | 開殻のパネル | DFT 6 件は成功（33.01/32.71/33.58）。CCSD(T) 3 件が input_invalid で coverage に失敗として数えられた。初回は configs の外に置いたため rc=2 | 1:45 | `backends/nwchem/engine.py:326-328` の `wft_closed_shell_only`。`pipeline/config.py:135` の methods の位置が暗黙 |
| P2a | Cl⁻ + CH3Cl、Cl⁻ を単量体に含む（`sn2_cl`） | 電荷 −1、単原子 | **rc=1**。minima stage で `ValidationError: n_external Input should be 5 or 6, input_value=3` | 35 s | `core/evidence.py:85` の `Literal[5, 6]`。parse の例外が `execution/jobs.py:151` から素通しで run を落とす |
| P2b | 同、Cl⁻ を外した版（`sn2_cl_b`） | 電荷 −1 の同一 SN2 | 成功。degenerate_rearrangement、TS 393.6i、ΔE‡ 10.45、ΔG‡ 11.08（幅 10.53〜11.58） | 10.6 分（paths 491 s、うち QRC 252 s = 51%） | 分離反応物基準は単原子未対応で出ない。手計算では PBE0/SVPD の TS は Cl⁻+CH3Cl より −0.74 kcal/mol 下（参考値） |
| P3a | CH3O• → CH2OH•、生成物を Cs で入力（`ch3o_doublet`） | 二重項 | 宣言反応は **blocked_upstream**（endpoints_not_on_one_pes）。代わりに鏡像対 mf1↔mf2 の「反応」（ΔG‡ 3.46）が ranking 1 位 | 10:09 | 宣言端点が鞍点（533i）に収束し、親に basin が付かない（`stages/minima.py:129`、`chemistry/hypotheses.py:120-126`）。同一性が真の回転のみ（`chemistry/identity.py:4`） |
| P3b | 同、生成物の対称性を崩した入力（`ch3o_doublet_b`） | 二重項 | reassigned_step。ΔE‡ 33.0、ΔG‡ 30.6、TS 2015i、⟨S²⟩ 0.754〜0.759。odft mult 2、xTB `--uhf 1`、pysis mult 2 が初めて通った | 6:13 | QRC の着地点が宣言生成物の鏡像（ΔE 5e-8 Eh）で、別 basin 扱い |
| P3c | H + H2、共線 vdW 錯体を端点（`h3_doublet_b`） | 開殻の二分子引き抜き | **unresolved**（attempts_exhausted） | 42 s | GFN2-xTB が両端を同じ偽の直線 H3 に潰す → Cartesian IDPP（DFT 最大 +81 kcal/mol、⟨S²⟩ 最大 1.19）→ ZTS の bead 4 で SCF が 3 回とも未収束。IDPP は種を返さない（`drivers/reaction_case/actions.py:203`） |
| P4a | マロンアルデヒドの分子内 PT（`p4_malon`） | 低障壁、ZPE で沈む障壁 | degenerate。TS 1100i、ΔE‡ 2.01、ΔZPE −2.40、ΔE0‡ −0.39、ΔG‡ 0.25 が注記なしで ranking 1 位 | 16:05 | 順位の量が未定義（`chemistry/gates.py:297-317` は符号を見ない）。序数（最速）自体は化学的に妥当で、欠けているのは値の意味 |
| P4b | H2O2 gauche(+)→gauche(−)（`p4_h2o2`） | 低障壁、ねじれ | TS（290.6i、1.15 kcal/mol、実験約 1.1）は正しく得た。QRC 4 本が棄却され **unresolved** | 3:51 | QRC の直線変位で出発点が TS より +0.055/+2.65 kcal/mol 上になり、`gates.py:222-225` の no_initial_descent/trajectory_above_ts で棄却（4 本とも 1 歩目から単調に下り gauche に収束）。xTB に gauche 極小がなく IDPP → ZTS（100 s）を余分に回した |
| P5a | DME、両メチルを 120° 回した粗い入力（`p5_dme`） | 縮退（メチル回転） | 宣言反応が **same_basin で黙って捨てられた** | 2:18 | 縮退判定が入力座標（`hypotheses.py:135, 141`）。入力の置換不変 RMSD 0.056 > 0.05（DFT 構造なら 0.0008） |
| P5b | 同、H を正確に置換した入力（`p5_dme_b`） | 同 | 成功。degenerate、226.9i、ΔE‡ 2.33、ΔG‡ 1.96（実験 V3 2.72） | 9:44（QRC 59%） | xTB TSOpt が対称性を破り、xTB Hessian の EF が 4 歩で収束 |
| P5c | C2v の重なり形 DME を種に実エンジンを直接駆動（`ho_harness.py`） | 負モード 2 本以上 | Cartesian moddir 2 は設計どおり動作し、2 次鞍点（−232i/−129i、+4.03 kcal/mol）→ higher_order。再試行の鞍点探索は maxiter で失敗（50 反復×3、495 s） | — | NWChem の Hessian 更新で生じた tiny eigenvalue（−6.0e-5）に対し `opt_drv.F` が |g| 長の最急降下を取り続ける這い（gmax 7.3e-4 一定、1 歩 ΔE −1.6e-6）。maxiter 継続が同じ壊れた Hessian で繰り返す。押す方向（modes[1] か modes[0] か）は原因ではない |
| P6a | I⁻ + CH3I（`p6_sn2_i`） | Z>36、電荷 −1 | **rc=1**。α 電子 58、救済を含め SCF 未収束（E≈−4052 Eh）、paths が「has no input of type minimum」で停止 | 2:37 | `backends/nwchem/input.py:82` が `* library <basis>` だけで ECP を書かない。手で `I library def2-ecp` を足すと 6.7 s で収束、α 30、E=−635.3667 Eh。`* library def2-ecp` は `bas_tag_lib` で abort |
| P6b | H2Te（`p6_h2te`） | 新しい元素 | structures で input_invalid（Te の共有結合半径なし）→ 次の dft stage が ValueError で **rc=1** | — | 元素表が 4 つに分かれ、半径表が事実上の許可リスト。`pipeline/runner.py:150-157` は入力 0 件の stage を例外にする |
| P6c | HI の SP（`ecp_test/hi`） | ECP の有無 | 正常終了、Level も一致、α 27（計 54 電子）、E=−2005.8147 Eh（**黙った誤り**。ECP ありなら 26 電子、≈−298.31 Eh） | 数秒 | ECP ブロックなし |
| T | アセチルアセトンの PT、15 原子（`timing_acac`） | 計算時間、中規模 | 60 分で打ち切り（rc=124）。鞍点探索が +3.25 から +0.04 kcal/mol まで滑落、36 歩で未収束。外部 SIGTERM の後 NWChem の MPI ランクが孤児として残った | >60 分 | 全原子 chord の |c|² の 82% をメチル H が占め、xTB Hessian（負 1 本、−1172 = PT、重なり 0.223）を棄却 → DFT Hessian 703 s（負 3 本）→ メチルねじれのモードを追跡（重なり 0.963）。追わない負モード（PT）は NWChem が trust いっぱいで下らせた |

**まだ実計算で通っていない経路**: SCREEN と FIND_PATH の barrierless（`state.py` 行 12・14）、collapsed → 中間体 → 分割、higher_order 再試行の成功、`_register` による新しい basin、開殻の WFT、単原子の熱化学、複数の温度と 1 M。§6 の検証セットで塞ぎます。

**時間の実測（4 コア、PBE0/def2-SVPD）**:

| 原子数 | SP | opt 1 歩 | 鞍点 1 歩 | 解析 freq | 備考 |
|---|---|---|---|---|---|
| 9（マロンアルデヒド） | 10 s | 11〜13 s | — | 150 s | QRC が paths 651 s の 310 s |
| 13 → 17（TMA·(HF)₂） | — | — | — | 267 → 898 s | 17 原子で opt 377 s、SP 35 s（v3〜v5） |
| 15（acac） | 36 s | 42 s | 33 s | 679〜703 s | freq 2 件でジョブ時間の 57% |
| xTB 全体 | — | — | — | — | run 全体で 1〜16 s（無視できる） |


---

## 3. 処理区分ごとの課題と改良案

各節は「現状 / 評価 / 課題 / 改良案 / 見送り / 維持」の順です。改良案の各項目は、変更・理論的根拠・効果・コスト・複雑さへの影響・検証方法を書きます。他の区分と重複する課題・改良案は担当区分にまとめ、ここでは「→ 担当 ID」で示します。

### 3.0 U0 全体プロトコルと理論レベル

**現状**

- 3 層の分業: 低レベル GFN2-xTB（配座・screen・explore・SCREEN）、停留点レベル PBE0-D3BJ/def2-SVPD・grid fine・`convergence energy 1e-7`（opt・saddle・string・解析 freq。`configs/methods/pbe0-d3bj_def2-svpd.yaml`）、熱補正は GoodVibes の qRRHO（vib_scale 1.0）。
- 既定の pipeline（`configs/pipelines/discover.yaml`、`known_endpoints.yaml`）には sp 層も `energy_method` もなく、順位は PBE0/SVPD の ΔG‡ で付く。composite（`stages/thermochemistry.py:40-45, 115-116`）は実装済みだが、使うには run の複製と追記の手順が要る（`docs/design.md:150`）。
- NWChem のデッキ（`backends/nwchem/input.py:73-104`）は `basis spherical / * library <basis>` だけで ECP を書かない。DFT の iterations は既定の 30、mult>1 は odft（UKS）、多重度 1 は常に RKS、SCF 救済は `damp 40 ncydp 30 lshift 0.5` の 1 回（lshift 0.5 は既定値そのもの）。
- MethodSpec（`core/method.py:18-35`）は、実行できない d4、smoke でしか使わない mp2、実 run のない solvation（cosmo/alpb）を受け付ける。化学の結果を変える maxiter と coordinates が JobStore の鍵の外の ExecutionSpec（`core/method.py:102-113`）にある。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | 中性・閉殻・主族（Z≤36）なら妥当。Z>36 で黙って誤る（唯一の critical）。既定の順位が PBE0/SVPD（重原子移動の障壁 MSE −6.6、H 移動 −4.2 kcal/mol）。順位は序数なので同じ反応クラス内では誤差が一部打ち消すが、クラスが混ざる汎用用途で順位が入れ替わりうる |
| 汎用性 | × | 黙って誤るのは Z>36（I）。単原子・Te・奇数ラジカルの SMILES は大声で落ちる。組成の多重度を宣言できない |
| 有用性 | △ | 順位の量が停留点レベルのまま。composite に手作業が要る |
| 効率 | ○ | 停留点は解析 Hessian、高精度は SP 層という分業は正しい（17 原子で freq 898 s、SP 35 s）。数値 Hessian への暗黙の切り替えは、既定の configs（停留点は PBE0 固定）では起きない |
| 簡潔さ | ○ | MethodSpec は小さい。使われない選択肢、鍵の外の化学 knob、手順書頼みの composite が複雑さを生む |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U0-I1 | def2 で Z>36 に ECP を書かない。HI は全電子（54 電子）で正常終了し E=−2005.81 Eh（黙った誤り）、I⁻+CH3I は SCF 未収束で run が落ちる | correctness | critical | `input.py:82`。NWChem は def2 ライブラリの `ASSOCIATED_ECP` を ecp 指令なしでは読まない（`src/basis/bas_input.F`）。`* library def2-ecp` は ECP のない H・C で `bas_tag_lib` abort（2 回実測）。def2-ecp の対象は Rb–La、Hf–Rn（P6） |
| U0-I3 | 既定の順位が停留点レベル（PBE0/SVPD）の ΔG‡。composite は既定で使われず手作業が要る | validity | high | P1b で SN2 の CCSD(T) 13.44 に対し PBE0/SVPD 10.45。`design.md:149` の「約 3 kcal/mol」は単分子異性化 2 系からの目安 |
| U0-I4 | 停留点に RSH・meta-GGA・COSMO・CD を指定すると NWChem は解析 Hessian を無効にし黙って数値 Hessian になる（約 3〜4 倍）。使われない solvation・d4・mp2 が残る | generality | low | `hess_check.F`、`task_hessian.F`（WSL の NWChem は conda バイナリでソース照合は未了）。既定 configs では起きない |
| U0-I5 | SCF の反復上限が既定 30、救済は実質 damping だけ。WFT の scf maxiter も既定のまま | validity | medium | `input.py:87-104, 101-102, 156-167`。H+H2 の ZTS bead 4 で 3 回未収束（P3c）。I⁻+CH3I の主因は ECP 欠落 |
| U0-I10 | ジョブ間で軌道を引き継がない（opt→freq、経路ノード） | validity | low | `input.py:129-132`、`actions.py:156-158`。開殻で別解に落ちた実例は未観測（仮説段階） |
| U0-I13 | 停留点の def2-SVPD がコストを決め、4 コアで約 35 原子が実用上限（freq に継続なし、timeout 4 h） | efficiency | low | freq ∝ N·nbf^2.8、`configs/sites/wsl_local.yaml:15` |
| 重複 | 元素範囲の不一致 → U1-I1、開殻 WFT → U8-I6、電子状態の宣言 → U1-I4・I5、単原子 → U3-I1、鍵の外の knob → U9-I2、スピン汚染の許容幅 → U1-I9、検証の空白 → U9-I11 | | | |

**改良案**

#### U0-P1 def2 の Z>36 に元素ごとの ECP 行を自動で書き、def2 以外の基底は拒否する（must、【中立】）

- 変更: `input._system`（`input.py:73-84`）で、basis 名が（小文字で）`def2` で始まり構造に Z>36 の元素があれば、`ecp` ブロックに `<El> library def2-ecp` を元素ごとに 1 行書く。`*` は使わない。Z>36 の元素があって basis が def2 系でなければ（WFT パネルの cc-pVnZ など）、全電子入力を書かずに ValueError（INPUT_INVALID）にする。Z は U1-P1 の元素表から取る。knob・MethodSpec・Level に項目は足さない（ECP は def2 の定義の一部で、method と元素から純関数で決まり、ジョブ鍵で区別される）。
- 理論的根拠: def2 の Rb–Rn は Stuttgart–Köln の def2-ECP と組で定義されている（Weigend & Ahlrichs, PCCP 7, 3297 (2005)）。NWChem の ECP 指定（https://nwchemgit.github.io/ECP.html）。UKS と ECP の組でも解析 Hessian は使える。
- 効果: 唯一の critical（黙った誤り）が消える。I・Sn・Te・Xe・Pd などが DFT と CCSD(T) で正しく動く。
- コスト: 約 10 行と単体テスト 2 件。
- 検証: (1) 単体: CH3I のデッキに `I library def2-ecp` が入り `*` が入らない。cc-pVDZ と I の組が拒否される。(2) smoke: HI の SP で α 13（計 26 電子、def2-ECP の芯 28 電子と整合）、E≈−298.31 Eh。(3) CCSD(T)/def2-TZVPD の HI で `freeze atomic` が凍結する軌道数を出力で確かめる（I は ECP 後 4s4p4d の 9 軌道）。(4) ecp ブロックがあっても `output.py` の Level（basis）観測が一致することを golden で確かめる。(5) I⁻+CH3I（P6a）が minima から paths まで通る（中心障壁の文献値約 8 kcal/mol は参考値）。

#### U0-P4 停留点レベルの原則を文書化し、使われない溶媒和・d4・MP2 を削除する（should、【複雑さ減】）

- 変更: `docs/design.md` に「opt・saddle・string・freq は大域混成または GGA + D3、RKS か UKS。RSH・meta-GGA・二重混成・CD・COSMO は SP 層だけ」と 1 行書き、停留点の method は configs で PBE0 に固定する。`MethodSpec.solvation` と各アダプタの COSMO/ALPB 分岐（約 37 行）、dispersion の `d4`、wft_method の `mp2` とその描画・読み取りを削除し、気相専用と明記する（U9-P12 を統合）。`_dft_supported`（`engine.py:58-61`）で PT2 を含む xc（b2plyp など）を拒否する（1 行。NWChem は `dftmp2` がないと PT2 抜きの値を黙って返す）。
- 理論的根拠: 解析 Hessian の可否（https://nwchemgit.github.io/Hessians-and-Vibrational-Frequencies.html）。溶媒和には、溶媒モデル・1 M 標準状態・数値 Hessian の費用までの一貫した設計が要り、半実装は正しさを保証しない。
- 効果: 実行できない・検証されない選択肢とエンジン間の食い違いが消える。約 −50 行。
- コスト: 削除が中心。
- 検証: b2plyp の MethodSpec が INPUT_INVALID。grep で solvation・d4・mp2 の参照 0 件。既存の smoke 12 件（MP2 の smoke は CCSD(T) に差し替え）。

#### U0-P5 SCF 収束の既定値を実用的にする（should、【中立】）

- 変更: すべての DFT デッキに `iterations 100` を書く（Setup 側なので鍵は変わらず、収束する計算の費用も変わらない。`ncydp 30` が「最初の 30 サイクルだけ減衰」として意味を持つ）。render_wft は scf ブロックを常に書き `maxiter 100`。救済策は、前回の movecs から再開して `cgmin`（二次収束、UKS 可）とする案と `convergence rabuck 20` 案を、P3c の bead 4 構造と ECP 修正後の I⁻+CH3I で比べて 1 本に決める。smear・fon は使わない（U9-P10 の SCF 部分を統合）。
- 理論的根拠: NWChem DFT の ITERATIONS 既定 30、SCF の MAXITER 既定 20、cgmin・rabuck（https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html）。smear/FON はエントロピー項で PES を変え、`same_pes` を黙って崩す。
- 効果: 開殻・陰イオン・TS 近傍・重元素で 30 反復打ち切りの偽の失敗が減る。
- コスト: 数行と、救済策を決める probe 1 本。
- 検証: 上記 2 構造で現行と候補を比較。閉殻 golden の反復数とエネルギーが不変。

#### U0-P13 大きな中性系用に def2-SVP の method を用意し、スケール上限を文書化する（could、【中立】）

- 変更: 既定は def2-SVPD のまま（陰イオンと H 結合錯体、ΔG_assoc の LOT を系の中でそろえるため）。`configs/methods/pbe0-d3bj_def2-svp.yaml` を置き、負電荷の化学種を含まない大きな系の選択肢とする（エネルギーは U8-P3 の TZVPD 層で補う）。`design.md`・`environment.md` に「freq ∝ N·nbf^2.8、15 原子で約 11 分、約 35 原子で timeout、freq に継続なし」を書く。
- 理論的根拠: 拡散関数は陰イオンと弱い相互作用のエネルギーに要るが、中性分子の構造・振動数への効果は小さい（Rappoport & Furche, doi:10.1063/1.3484283）。
- 効果: 40 原子級の中性系で freq が約 3〜4 倍速い選択肢（要実測）。
- コスト: yaml 1 本と文書。
- 検証: HCN・HONO・NH3・acac で SVP と SVPD を比較（構造 RMSD、最低振動数、composite の ΔG‡ 差 < 0.5 kcal/mol）、時間比を記録。

他区分へ統合: 元素表 → U1-P1、composite を既定に → U8-P3、電子状態の宣言 → U1-P3、単原子 → U3-P1、開殻 WFT → U8-P7、ExecutionSpec の整理 → U9-P2、符号 blocker 削除とパネル列 → U8-P2、M06-2X → U8-P13、検証セット → §6。

**見送り**

- ECP を `* library def2-ecp` で一括指定: NWChem 7.2.3 が abort する（実測）。
- Level・MethodSpec に ECP 項目を足して `same_pes` で照合: ECP は基底名と元素の純関数。描画の単体テストと HI の golden で足りる。
- freq 出力の `Finite-difference Hessian` を検出して Failure にする（U0-P4 原案）: 検出時点で 3〜4 倍の計算は済んでおり、結果は（ノイズを除けば）有効。原則の文書化と configs の固定で足りる。
- 停留点レベルを RSH・meta-GGA にする／CD fitting で高速化／ROKS を既定に: 解析 Hessian が無効になり freq が 4〜6 倍遅くなる。高精度化は SP 層で行う。
- SCF 救済に smear・FON: PES が変わる。
- BS-UKS 最適化と Yamaguchi 射影を本体に: 複雑さに見合わない（検出だけ U6-P11）。
- CCSD(T) を既定のエネルギー層に／DLPNO: N^7 で重原子約 5 個まで、DLPNO は NWChem にない。
- freq だけ grid・基底を変える、autosym: `same_pes` と入力フレームの契約が崩れる。
- 組成・電荷で停留点の基底を自動で切り替える: 系内で LOT がそろわず ΔG_assoc と `same_pes` が壊れる。
- 同じ基底の後続ジョブを親の movecs から始める（U0-P9）: 結果（SCF 解）を左右しうる入力を鍵の外で渡すことになり、実例もない。代わりに freq と親のエネルギー一致の 1 条件（U3-P7）で検出する。
- 陰イオンの HOMO>0 ゲート: SIE で束縛陰イオンにも誤警報を出す。

**維持**: 停留点を 1 レベル（PBE0-D3BJ/def2-SVPD、grid fine）にそろえること、全 DFT タスクで `convergence energy 1e-7`、freq を独立ジョブにし初期 Hessian は inhess 2、Level の観測照合と `same_pes`、電荷・多重度の全エンジンへの受け渡し、mult>1 の UKS、composite の式、入力フレームのまま扱う `nocenter noautosym`、振動数は未スケール 1 因子。

---

### 3.1 U1 入力構造・電子状態（電荷・多重度・元素）

**現状**

- `core/system.py:24-38`: SpeciesInput は charge=0、multiplicity=1 が既定。CompositionInput は components だけで多重度を宣言できない（extra=forbid）。
- `chemistry/smiles.py:10-38`: 形式電荷の照合、同位体・複数断片の拒否、ETKDG。ラジカル電子数と多重度は照合しない。
- `chemistry/electronic_state.py`: 118 元素の表（:7-14）、パリティ検査（:17-33）、`coupled_multiplicities`（:36-58、2S の整数演算で厳密）。
- `stages/conformer_search.py:206-211`: 組成の多重度は `coupled_multiplicities(...)[-1]`（最高スピン）、`check_electronic_state` を再実行（各成分が妥当なら構造上必ず通る冗長検査）。
- 元素表が 4 つ: `topology.py:29-45`（共有結合半径 H–Kr と I、vdW 半径は Sc–Co 欠落）、`vibrations.py:18-28`（質量 H–Kr と I）、`goodvibes/worker.py:38-39`（dict の挿入順の先頭 36 個に I を足して原子番号を作る）、`electronic_state.py`（118 元素）。
- `chemistry/topology.py`: 結合は 1.15/1.45 Σr_cov の帯、ヒステリシスは片方向（`bond_changes` の b 側だけ previous 付き、:117-124）、状態ラベルは 1-WL 3 反復。宣言座標・プロトン座標・酸塩基の表が同居。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | 電荷・多重度の一貫した受け渡しとパリティ・角運動量の計算は正しい。組成の最高スピン固定（²R+³O2 が四重項）、偶数不対電子の SMILES（[O][O]、[CH2]）が黙って RKS 一重項、1.45 Σr_cov の倍率則が I···N 2.80 Å などを共有結合とみなす、bond_changes が向きに依存 |
| 汎用性 | × | 失敗の場所が元素ごとに違う（Te は structures、Sc–Co は極性 H のない組成の placement、Rb〜Xe は freq）。I は全表を通って ECP なしで誤る |
| 有用性 | △ | スピン状態を指定できず、開殻一重項が対象外であることも入口で示されない |
| 効率 | ○ | U1 自体は無視できる。未対応元素の検出が遅く上流の計算が無駄になる |
| 簡潔さ | ○ | 小さいが、元素表 4 つ、冗長な再検査、片方向の previous、topology への責務外の同居 |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U1-I1 | 元素表が 4 つに分かれ対応範囲が暗黙に食い違う。未対応元素は入口でなく途中で、場所ごとに違う形で失敗する | generality | high | 上記の表。P6b は structures で拒否後、`runner.py:157` が上流の失敗理由を付けて停止（意図した停止だが report が出ない → U9-I1）。`placement.py:171` は `conformer_search.py:212` の try の外 |
| U1-I4 | 組成の多重度を宣言できず最高スピンに固定される | validity | medium | `conformer_search.py:207`。GFN2 はスピン非依存で低レベルでは選べない（https://xtb-docs.readthedocs.io/en/latest/spgfn.html）。configs に開殻組成は 0 件で既存 run への影響はない |
| U1-I5 | SMILES のラジカル電子数と多重度を照合しない。偶数不対電子の開殻（[O][O]、[CH2]）が黙って RKS 一重項になる（奇数はパリティ検査で捕まる） | correctness | low | `system.py:24`、`smiles.py:25-30` |
| U1-I6 | 組成での `check_electronic_state` 再実行は構造上必ず通る | complexity | low | `conformer_search.py:208-211` |
| U1-I7 | 1.45 Σr_cov の倍率則が重原子・イオンの二次的接触を共有結合とみなす。`bond_changes` が向きに依存する | validity | medium | `topology.py:24-25, 79-89, 117-124`。比 1.3→1.6 の対は a→b で broken、b→a で変化なし。`hypotheses.py:136, 191` が使う |
| U1-I8 | 状態ラベルの WL が 3 反復固定で、遠い位置異性体を区別しにくい（8 桁ハッシュの衝突確率は無視できる） | correctness | low | `topology.py:27, 211-239` |
| U1-I9 | spin_ok の許容幅が S に依らない絶対値 0.1。エネルギー層の SP は検査しない | validity | low | `gates.py:106-112`。高スピン UKS は実際には汚染が小さく、不当に止めた実例はない。SP 層の未検査のほうが実質的 |
| U1-I10 | 開殻一重項（ビラジカル性の TS、ラジカル対）が対象外であることを入口で明示しない | validity | medium | `input.py:92-94`、`design.md:138` |
| U1-I11 | topology に宣言座標・プロトン座標・酸塩基の表が同居。`_LABILE_PARTNERS` に P があり Se・Te・As がない | responsibility | low | `topology.py:46-48, 127-208` |
| U1-I13 | 開殻原子・直線ラジカルの電子縮重度（L>0）とスピン軌道補正が g=2S+1 近似のまま、文書にもない | validity | low | GoodVibes `thermo.py:219-229`。SO の安定化 Cl 0.84、Br 3.5、I 7.25 kcal/mol |
| 重複 | ECP → U0-I1、単原子 → U3-I1、検証の空白 → U9-I11 | | | |

**改良案**

#### U1-P1 元素表を 1 つにし、対応範囲を入口で 1 回だけ判定する（must、【複雑さ減】）

U0-P2、U2-P3、U9-P4(3) を統合。

- 変更: `hfauto/chemistry/elements.py` を新設。範囲は Z=1–57 と 72–86（Ce〜Lu を除く）。列は記号、Z、主同位体の質量、共有結合半径、vdW 半径の 5 列だけ。RDKit の PeriodicTable から一度だけ静的に生成し、実行時依存にはしない。**既存の H–Kr と I の半径は今の値（Cordero、Bondi/Mantina）をそのまま使い、追加した元素だけ Cordero 2008／Alvarez 2013 を使う**（C の vdW 1.70→1.77 のように変えると placement と状態ラベルが既存 golden で動く）。ECP の要否は元素表に持たせず、`input.py` の 1 か所（U0-P1）に置く。電気陰性度や最大配位数の列は、U4-P1 で必要になった時点で足す。
- 削除: `electronic_state._ELEMENTS`、`topology._COVALENT_RADII_A`・`_VDW_RADII_A`、`vibrations.ISOTOPIC_MASSES`、`goodvibes/worker._ATOMIC_NUMBERS`（挿入順に頼る作り方）。
- 判定: `check_electronic_state` が表にない元素を `unsupported_element` とし、structures で INPUT_INVALID にする（新しい検査は増えず、既存の検査の範囲を表に合わせるだけ）。対応範囲は「計算できる」範囲として宣言し、遷移金属を「検証済み」とは書かない。
- 理論的根拠: GFN2-xTB の範囲 Z≤86（Bannwarth 2019, doi:10.1021/acs.jctc.8b01176）、NWChem の def2-SVPD/TZVPD ライブラリは H–La、Hf–Rn（検証段で確認）。
- 効果: 失敗が入口の 1 種類になり、上流の CREST・xTB・DFT の無駄と途中の例外が消える。対応範囲が 37 元素から約 72 元素へ。原子番号の隠れた結合（Te で KeyError）が消える。
- コスト: 表のデータ約 90 行、旧表 4 つ（約 40 行）の削除と import の付け替え。定義の重複 4→1。
- 検証: 単体: 全元素で Z と質量が一致、Ce・Fr が INPUT_INVALID、**剛体配置になる組成（FeCl3·CH4、Fe(CO)5·CO2）**で placement が例外を出さない（Fe(CO)5·H2O は H 結合配置の経路を通るので欠陥を再現しない）。WSL: H2Te を structures→dft まで通す。既存 golden（HCN、HONO、NH3、TMA·(HF)₂）で質量・σ・G が 1e-10 Eh 以内で一致。

#### U1-P3 スピン状態を宣言できるようにし、一意でなければ宣言を求める（must、【中立】）

U0-P6、U2-P1 を統合し、C8 を置き換える。

- 変更:
  - `SpeciesInput.multiplicity: int | None`（既定 None）。None のとき SMILES は RDKit の `NumRadicalElectrons+1`、xyz は 1。**d ブロック元素を含む化学種で None なら INPUT_INVALID（宣言必須）**（RDKit は [Fe+2] などで 0 を返しやすい）。明示した値はパリティ検査だけ（一重項カルベンなど閉殻一重項の明示は拒否しない）。
  - `CompositionInput.multiplicity: int | None`。None のとき `coupled_multiplicities` が 1 つならその値、複数なら INPUT_INVALID（`declare_multiplicity: candidates (…)`）。宣言した値は集合に含まれるかを検査。開殻の成分から多重度 1 を作る宣言は INPUT_INVALID（`open_shell_singlet_unsupported`）。エラー文と `design.md` に「ラジカル再結合の生成物は単量体として宣言する。BS-UKS は範囲外」と書く。
  - 純関数 `electronic_state.composition_multiplicity` を置き、`conformer_search.py:207` の `[-1]` と冗長な再検査（:208-211）を削除する。複数のスピン状態は別の組成 id として宣言できる（composition_id が `_m` で区別）。
- 理論的根拠: 角運動量の結合で得られる多重度のどれが重要かは反応による（²R+³O2→²RO2）。GFN-xTB はスピン非依存で低レベルでは選べない。RKS は開殻一重項を記述できない。RDKit のラジカル電子数は高スピン配置の不対電子数（[C] は RDKit で五重項、実際は ³P。原子の SMILES は高スピン仮定が外れうると文書に書く）。
- 効果: スピン状態を黙って選ばなくなり、R•+O2 の二重項面などを宣言できる。閉殻どうし・開殻 1 つの組成（今の検証系すべて）は挙動が変わらない。
- コスト: system.py に 2 項目、conformer_search・smiles に約 10 行（冗長検査の削除で相殺）。テスト 6 件。
- 検証: 単体: HF+HF → 1、OH•+OH• 未宣言 → INPUT_INVALID、宣言 3 → 可、宣言 1 → open_shell_singlet_unsupported、²X+¹Y → 2、[O][O] → 3、[CH3] → 2、[Fe+2] 未宣言 → INPUT_INVALID。WSL: CH3•···O2（m=2）の組成で CREST に `--uhf 1` が渡り、UKS の Level（mult 2）と ⟨S²⟩ を確かめる。

#### U1-P4 結合判定を加算型の許容値にし、`bond_changes` を集合の差にする（should、【複雑さ減】）

U2-P7、U4-P3 を統合。

- 変更: 結合は r < Σr_cov + 0.4 Å の単一構造の規則だけにし、`previous` 引数とヒステリシスを削除する。`bond_changes(a, b)` は `formed = bonds(b) − bonds(a)`、`broken = bonds(a) − bonds(b)` の集合の差にする（構成上対称で、fragments・state_label と一致する。3 状態の定義は、ラベルが違うのに「変化なし」となる不整合を生むので採らない）。H の原子価を 1 とする特例は入れない。イオン–双極子錯体が 1 断片になる場合（K⁺···π など +0.4 Å でも結合に残る）は状態ラベルの限界として `design.md` に 1 行書く。
- 理論的根拠: 加算型の許容値（Meng & Lewis, J. Comput. Chem. 12, 891 (1991)）は SCINE の BondDetector、OpenBabel、RDKit の connect-the-dots が採用。静的な 2 構造を比べる関数は向きに対称でなければならない。
- 効果: ハロゲン結合（I···N 2.8 Å は非結合）、Cl⁻···H–F の断片とラベルが正しくなる。bond_changes が向きに依らなくなり、仮説の torsional/changed と trials の判定が安定する。
- コスト: 約 10 行の変更と previous の配管削除。状態ラベルの値が変わるので新しい run dir で取り直す（DFT ジョブは幾何の鍵で再利用）。
- 検証: 単体: I···N 2.8 Å・I···I 3.5 Å が非結合、FHF⁻ 1.14 Å・I3⁻ 2.92 Å が結合、bond_changes の対称性を性質テストで。**TMA·(HF)₂ の N···H 1.43 Å と閾値 1.42 Å、K⁺···C 3.20 Å と閾値 3.19 Å は差 0.01 Å しかない**ので、過去 run の構造 95 本（v2〜v5、w7）でラベルを再計算し差分を実際に見る。NH3···ICl が 2 断片になる。既知端点 4 系の outcome 不変。

#### U1-P5 WL を分割が安定するまで回し、原子クラスを公開する（could、【中立】）

- 変更: `_wl_hash` を `wl_classes(symbols, bonds) -> tuple[int, ...]` と `state_label` に分け、色の分割数が変わらなくなるまで（最大 N 回）反復。`_WL_ITERATIONS` を削除。表示形式（8 桁）は変えない。U4-P2 の drive 重複除去がこれを使う。U1-P4 と同じ回に入れる。
- 理論的根拠: WL の色の細分化は安定するまで回すのが標準（Shervashidze et al., JMLR 12, 2539 (2011)）。
- 効果: 大きな分子の位置異性体の衝突がなくなり、原子の同値性の定義が 1 つになる。
- コスト: 約 10 行。検証: 1,6- と 1,7- 二置換長鎖が別ラベル、ランダム置換で不変。

#### U1-P6 spin_ok を S(S+1) に対する相対許容幅にする（could、【中立】）

- 変更: `Policy.spin_contamination_tol` の意味を相対値（既定 0.1）に変え、式を `tol·max(S(S+1), 0.75)`（S=0 の BS 一重項で分母 0 を避ける床）にする。knob の数は変えない。エネルギー層 SP への適用は U8-P8 で行う。
- 理論的根拠: スピン汚染は S(S+1) に対する比で見るのが慣用（Young, *Computational Chemistry*, 2001）。
- 効果: 多重度で検査の厳しさが変わる問題を直す。実害はまだ観測されていないので could。
- 検証: OH• 0.7526、CH3O 系 TS 0.759 が通り、H3 の IDPP 点 1.19 が落ちる。³O2 と五重項の境界の単体テスト。

#### U1-P7 topology を「結合グラフと状態ラベル」だけにする（should、【複雑さ減】）

- 変更: `declared_coordinate` とその勾配を `geometry.py` へ、labile/acceptor の表と関数を主な利用者の placement へ移す（U4-P1・U2-P2 の後は trials・conformer_search からの依存が消える）。`transferred_hydrogens`・`proton_coordinate` は U7-P1 で dzpe blocker を削除するのと同時に削除する。残すのは `bonds`、`bond_changes`、`fragments`、`wl_classes`、`state_label`。
- 効果: U1 の責務がはっきりし、約 40〜60 行減る。
- 検証: lint-imports、ruff、pyrefly、全テスト、既知端点 4 系の退行確認。

#### U1-P9 電子状態の対象範囲を文書で明示する（should、【中立】）

- 変更: README と `design.md` に、(1) 対応元素（U1-P1 の表、Z>36 は def2-ECP）、(2) スピン状態は利用者が宣言、d ブロックは宣言必須、開殻一重項・MECP は対象外、原子の SMILES は高スピン仮定が外れうる、(3) 開殻原子と ²Π ラジカルは g=2S+1 近似でスピン軌道補正なし（Cl 0.84、Br 3.5、I 7.3 kcal/mol）、(4) 電荷系の SIE と D3 の電荷非依存の限界（エネルギー層で確かめる）、(5) イオン–双極子錯体が 1 断片になりうる、を書く。U7-P9 もここに統合。
- 理論的根拠: Mori-Sánchez, Cohen, Yang (2008) の非局在化誤差、Bauernschmitt & Ahlrichs (1996) の不安定性。
- 効果: 診断やゲートを足さずに誤用を防ぐ。コスト: 文書 10 行程度。

他区分へ統合: 単原子 → U3-P1、検証セット → §6。

**見送り**

- 組成の既定を最高スピンのまま残す: 化学的に誤った面を黙って選ぶ。開殻組成は configs になく互換性の損失もない。
- 全多重度を DFT で計算して最安を自動選択／xTB（spGFN を含む）で多重度を選ぶ: 計算が数倍になり、反応に関わる面は最安とは限らない。GFN はスピン非依存、spGFN は較正が弱い。
- WBO で結合グラフを作る: 計算器・電荷・スピンに依存し、イオン・配位結合で小さく、純関数性も崩れる。
- 状態ラベルに立体（キラリティ、E/Z）を含める: TS 近傍・錯体で立体認識が不安定。鏡像は U3-P2 の縮重で扱う。
- H の原子価 1 の特例（3 中心扱い）: TMA·(HF)₂ の 2 端点は q=0.22/0.14 Å < 0.3 Å でどちらも 3 中心に残り、区別の効果がない。
- `bond_changes` の 3 状態定義: 集合の差のほうが単純で矛盾しない。
- RDKit の DetermineBonds／xyz2mol: 錯体・イオン対・TS 近傍で失敗しやすく、ラベルに要るのは連結性だけ。
- Ce〜Lu、Z>86、相対論の全電子計算: ライブラリと xTB の範囲外。
- 開殻一重項を BS-UKS で入力から扱う: 複雑さが大きい。
- 電子縮重度とスピン軌道の入力欄（U1-P10、U7-P9 原案）: 対象系がまだない。文書の定量値で足りる。
- NWChem の α/β 電子数を ECP の芯電子の表と照合する実行時ゲート: 決定的に書いたうえで golden で 1 回確かめれば足りる。

**維持**: パリティ検査と `coupled_multiplicities` の整数演算、`composition_id = Hill_q{charge}_m{mult}` を PES の同一性にし宣言反応の両端で一致を要求すること、電荷・多重度の全エンジンへの受け渡しと NWChem 出力での照合（`engine.py:232-248`）、開殻は UKS、SMILES の形式電荷照合・同位体と複数断片の拒否、宣言反応の端点に xyz を要求すること、状態ラベルを元素ラベル付き結合グラフの置換不変ハッシュにすること。

---

### 3.2 U2 配座探索・錯体配置

**現状**

- 単量体: `crest --gfn2 --quick -T <threads> --ewin 6 --chrg q --uhf m-1`（`backends/crest.py:44-60`）。重原子 3 個以下で回転可能結合 0 本（単原子を含む）は探索を省く（`placement.py:257-259`、`conformer_search.py:162`）。トポロジー停止は `crest_topology` として残し 1 回だけ再実行（C34、`conformer_search.py:94-111`）。
- 組成: placement.seeds で seed を作り（受容原子の lone-pair 円錐に H を置く約 110 行か、Σr_vdW+0.5 Å の剛体配置）、CREST に渡すのは seed00 だけ（`--nci --quick --noopt`、酸塩基のときだけ labile H と受容原子に `--notopo`、`conformer_search.py:187-191, 217-219`）。残りの seed は CREST 失敗時の出力にだけ使う。
- 後処理: 独自の重複除去（RMSD<0.1 Å かつ ΔE<0.1 kcal/mol、:114-119）、状態ラベルごとに GFN2 の低い keep_per_state（6）個、誰も読まない `diagnostics.json`（:243-244）。
- 実績: v3 の TMA·(HF)₂ で CREST は run の 1% 未満。開殻の組成（NH3·HF⁺）は `--noopt` で trial MTD が失敗（v2）。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | ○ | iMTD-GC と下流 Registry の分担は標準的。主な欠陥は組成の最高スピン固定（→ U1-P3）。状態ラベルがまとめるのは「強い H 結合の中性と共有プロトン」で、真のイオン対は別ラベルになる（F–H の閾値 1.28 Å） |
| 汎用性 | △ | 元素表の欠落（→ U1-P1）と開殻組成の CREST 失敗。`--notopo` が漏れるのは極性 H のない供与体・受容体の組成（BF3·NMe3、SO2·NMe3、ハロゲン結合、極性 H のない配位子）で、陽イオン–水やハロゲン化物–水では水の O・H が除外されるので漏れない |
| 有用性 | ○ | 最安の配座と状態を安く供給。seeds_per_composition は失敗時しか効かない |
| 効率 | ◎ | run の 1% 未満 |
| 簡潔さ | ○ | 重複除去が 3 重（CREGEN、_same、Registry）、読まれない診断、円錐配置 110 行 |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U2-I2 | 開殻の組成で `--noopt` の未緩和 seed から始めた trial MTD が失敗する。組成の CREST argv が電子状態によらない 1 規則になっていない | untested_path | medium | `crest.py:52-56`。W7 の smoke では preopt+`--noreftopo` で rc 0。CREST `setuptest.f90:352-385` は初期最適化後のトポロジー検査で excludeTOPO を見ず、`reftopo=.false.` のときだけ続行 |
| U2-I4 | `--notopo` の対象と H 結合の配置が labile H・受容原子の元素表に依存。極性 H のない供与体・受容体の組成で、CREGEN が配位結合のできた付加体を捨てうる | generality | medium | `conformer_search.py:187-191`、CREST `ztopology.f90:377-385`（bondtotopo_excl）、`cregen.f90:481`。実 run での事例はまだない |
| U2-I5 | 円錐配置（約 110 行）は seed00 しか CREST に渡らないので価値が限られる | complexity | low | `placement.py:92-163`、`design.md:139` |
| U2-I6 | 重複除去が 3 重。conformers の `_same` は CREGEN 済みの出力に対して実質働かない | complexity | low | `conformer_search.py:114-119`、CREGEN の RTHR 0.125 Å（v3 stdout） |
| U2-I8 | 状態ラベルが強い H 結合の中性と共有プロトンを同じラベルにまとめる（U1 から波及） | validity | low | v3 の 3 配座は同じ結合状態で、keep_per_state による欠落は起きていない |
| U2-I9 | 誰も読まない `diagnostics.json` と `topology_removed`、key_payload の冗長な `"noopt"` | over_diagnostics | low | `conformer_search.py:243-244`、`crest.py:40, 68-70, 169` |
| U2-I10 | 大きく柔らかい単量体で GFN2 の MTD 長が指数的に伸びる | efficiency | low | CREST `choose_settings.f90:74-100`。今の系では 6〜12 s（将来の懸念） |
| 重複 | 最高スピン固定 → U1-I4、元素表 → U1-I1、鏡像 → U3-I2、threads_per_item → U9-I6 | | | |

補足（訂正）: 「CREST が鏡像をまとめ identity は別扱い」という不一致は事実ですが（`cregen.f90:1157` の enantio=.true.、`identity.py:1-6` の det=+1）、U2 への影響はキラル配座の片方しか出力に入らないことだけで、U2 の欠陥ではありません。鏡像の扱いは U3-P2 で決めます。

**改良案**

#### U2-P2 組成の CREST 呼び出しを電子状態によらない 1 規則にする（should、【複雑さ減】）

- 変更: 組成では `_notopo`（`conformer_search.py:187-191`）を削除し、`--notopo` に全原子の添字を渡す（状態は hfauto の state_label が決め、CREST はサンプリングだけ）。初期最適化は (A) 現状の `--noopt` と (B) `--noopt` を外して `--noreftopo` をプローブで比べ、両方通るなら 1 規則を採る。どちらも全系で通らなければ (A) を維持し、失敗時の seed 出力を文書に明記する（分岐は増やさない）。単量体は今のまま CREGEN のトポロジー検査を使う。
- 理論的根拠: 原子リスト付き `--notopo` でも discardbroken（共有結合解離の除去、`cregen.f90:629-751`）は残る。`--noreftopo` は初期最適化後のトポロジーを参照にして続行（`confparse.f90:1580-1581`）。
- 効果: conformers が酸塩基の表に依存しなくなり、配位付加体・ハロゲン結合の状態を CREGEN が捨てなくなる。開殻組成の CREST が通る見込み（W7 で rc 0）。
- コスト: 約 −6 行。プローブ 7 系 × 2 通りで約 10 分。
- 検証: TMA·(HF)₂、F⁻·(HF)₂、NH3·HF⁺（m=2）、OH•·H2O（m=2）、Cl⁻·CH3Cl、Na⁺·H2O、**SO2·NMe3 または BF3·NMe3（極性 H のない型）**で rc、配座数、screen 後の状態ラベル集合と各ラベルの最低 xTB エネルギー、時間を比べる。採用条件は既存の状態を失わず開殻系が rc 0 になること。

#### U2-P4 conformers の独自の重複判定と読まれない診断を削除する（should、【複雑さ減】）

- 変更: `_same`、`DUPLICATE_RMSD_A`・`DUPLICATE_DE_KCAL` を削除し、`_select` を「状態ラベルごとにエネルギーの低い keep_per_state 個（エネルギーのない候補の扱い BUG-08 は維持）」だけにする。conformers の `diagnostics.json`、`_Found.diagnostics`、`ConformerEnsemble.topology_removed`、`crest.py` の `_REMOVALS`・`topology_removed()`、key_payload の `"noopt"` を削除する（鍵が変わるので新しい run dir で取り直す）。placement の `DUPLICATE_RMSD_A` は未最適化配置の多様性フィルタで目的が違うので残す。minima の `diagnostics.json` は mode-follow 経過の唯一の記録なので残す（U3・U9 の判定）。
- 理論的根拠: 同一性は最終レベルで 1 回だけ確定する（CREGEN の出力はすでに重複除去済み、残りは Registry が吸収）。監査情報は JobStore の stdout と `crest_topology` の source に残る。
- 効果: 約 −30 行、定数 2 つ減。
- 検証: pytest。v3 の TMA·(HF)₂ discover を新しい run dir で再計算し、screen の basin 集合と DFT 極小が不変。

#### U2-P6 placement を元素によらない剛体の接触配置 1 つにする（could、条件付き、【複雑さ減】）

- 変更: U2-P2 を採用し、そのプローブで剛体配置の seed でも TMA·(HF)₂ と F⁻·(HF)₂ の最安状態を失わないことが確かめられた場合に限り、円錐配置（`placement.py:92-163`、`CONE_DEG`・`AZIMUTHS_DEG`・`TILT_DEG`・`HBOND_SCALE`）を削除する。距離は今の Σr_vdW + 0.5 Å のまま（新しい係数は足さない）。
- 理論的根拠: iMTD-GC は NCI の壁の中で配置そのものをサンプリングするので、seed に要るのは衝突がないことだけ（placement の docstring 自身の前提）。
- 効果: 約 −110 行、定数 4 つと元素名の列挙表 2 つが消える。
- 検証: U2-P2 と同じ系で円錐と剛体を比べ、各ラベルの最低エネルギー（許容 0.1 kcal/mol）と時間。1 系でも最安状態を失えば採らない。

他区分へ統合: 組成の多重度 → U1-P3、元素表 → U1-P1、結合判定 → U1-P4、鏡像 → U3-P2、並列化 → U9-P6。

**見送り**

- 全 seed で CREST を走らせる: 費用が seed 数倍。多様性の不足は状態ラベルの改善で扱う。
- seed を conformers 内で xTB 緩和してから渡す: 2 つ目のエンジン依存が増える（U2-P2 が全系で失敗した場合だけ再検討）。
- conformers にエネルギー窓を戻す: 緩和前の GFN2 で切ると偽陰性。選抜は screen/dft の窓 1 か所。
- 第 2 の配座エンジン（RDKit ETKDG）: 失敗時は seed と入力で縮退すれば足りる。
- 組成に GFN-FF でサンプリング／単量体に `--gfn2//gfnff` の選択肢（U2-P8）: GFN-FF はトポロジー固定でプロトン移動を表せない。単量体の選択肢は対象系が現れるまで入れない（今は 6〜12 s）。
- identity で鏡像をまとめる案を U2 の側で実装（U2-P5 原案）: 鏡像間の障壁（H2O2 の g+→g−）が同一 basin で消える懸念がある。鏡像の扱いは、縮退判定に掌性を含める U3-P2・U5-P6 で解く。
- 新しい診断ファイル（CREGEN 除去数など）: 読むコードがない。

**維持**: CREST iMTD-GC（`--gfn2 --quick --ewin 6 --chrg --uhf`）を唯一の配座エンジンにすること、小さく剛直な単量体（単原子を含む）の探索省略、単量体のトポロジー停止の `crest_topology` 化と 1 回の再実行（C34、実 CREST で 2/2 成功）、状態ラベルごと keep_per_state でエネルギー窓なし、CREST 失敗時に seed か入力を出して下流を止めないこと、placement の分子固定座標系と決定的 RNG と衝突検査、attempt ディレクトリで `--scratch` なし、版数 pin の照合、CREST の thread_map 並列とコアのセマフォ。

---

### 3.3 U3 極小構造の確定（screen/dft/registry/同一性）

**現状**

- `drivers/minimum.py`（277 行）の `relax_to_minimum`: opt（:186）→ Registry で既知 basin なら freq を省く（:189-194）→ 別ジョブ freq と `is_minimum`（:82-94）→ saddle 級は最低虚モードに沿って ±0.1 Å の mode-follow を最大 2 サイクル（`_follow`、:125-144、`amplitude_A=0.1` :33）→ soft 級は + 側へ 1 回（`_soften`、:147-157）。変位後の再最適化は初期 Hessian なしで trust 0.1。`Registry.find`（:227-236）の composition_id・level_key は任意で、既知判定（:191）は渡さない。
- `chemistry/identity.py`: 真の回転（det=+1）だけの置換不変 RMSD。`same_minimum`（:120-132）と `assign`（:135-162、RMSD ≤0.05 Å、|ΔE| ≤5e-5 Eh、次点規則）が同じ基準を二重実装。
- `chemistry/selection.py`: 組成 × state_label ごとに 6 kcal/mol・3 構造の窓。`rerank`（:57-65）は `window_kcal=inf`。
- `stages/minima.py`: screen（xTB vtight）と dft（NWChem opt + 解析 freq）。species 順に直列で relax・登録。mode-follow の両側は新しい species（:137-161）。宣言端点が ts_candidate に落ちると親は failed にも basin にもならない（:129）。
- `core/evidence.py:85` の `n_external: Literal[5, 6]`。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | ○ | opt と別ジョブ freq の照合、4 段の虚振動分類、mode-follow の TS 候補化、既知 basin での freq 省略は堅い。鏡像の偽反応、鞍点端点で basin が空、固定 0.1 Å の変位（見積もりで柔らかいモードの勾配が収束閾値程度。実 run では未観測） |
| 汎用性 | △ | 単原子で run が落ちる。Z>36・元素表・SCF 上限の欠陥が U3 で表面化 |
| 有用性 | ○ | basin・members・state_label は後段の共通キーとして機能 |
| 効率 | ○ | 既知 basin の freq 省略と錯体の xTB 初期 Hessian（−35%）は有効。rerank 後の窓なし、explore 生成物の未重複除去、変位後の trust 0.1 が無駄 |
| 簡潔さ | ○ | 756 行と小さい。same_minimum と assign の二重実装、入力座標による縮退判定（→ U5-I8）、_soften の別分岐、使われない amplitude_A |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U3-I1 | 単原子（Cl⁻、F⁻、H•、ハロゲン原子）で Evidence の検証が失敗し run 全体が落ちる。熱化学側も単原子を拒否する | generality | high（黙った誤りではなく大声で落ちる） | `evidence.py:85`、`xtb.py:187` の cast、`nwchem/engine.py:289-291`、P2a の `sn2_cl.log`。熱化学側は `goodvibes/worker.py:94-97`（S_rot≤0）、`goodvibes/engine.py:102`（hessian 必須）、GoodVibes 4.3.0 `thermo.py:685`（空の rotemp で計算全体を飛ばす）。NWChem の単原子 freq 自体は完走、identity・imaginary_tier・振動本数の検査も自明に通る |
| U3-I2 | 鏡像体を別 basin として扱い、鏡像化が「反応」として順位に載り、宣言反応が reassigned になる | validity | high | `identity.py:4, 28-33`。P3a の mf1↔mf2（ΔG‡ 3.46、1 位）、P3b の reassigned、p4_h2o2 の h2o2_m/h2o2_p。キラル分子・TS の縮重 m=2 がどこにもない |
| U3-I3 | `same_minimum` と `assign` が同じ基準を二重実装 | complexity | medium | `identity.py:120-162`。placement の 0.1 Å と conformers の `_same` は目的が違う（多様性フィルタと枝刈り）ので対象外 |
| U3-I4 | 既知判定で composition_key と level_key を渡さず、ほぼ縮退した一重項と三重項が別多重度の basin に known として join し freq も省かれうる | validity | low | `minimum.py:191, 227-236` |
| U3-I5 | mode-follow の変位が固定 0.1 Å で、50〜100i の柔らかい saddle 級では変位点の勾配・ΔE が収束閾値と同程度になりうる | validity | medium | 見積もり（ν=50i、換算質量 10 amu で勾配約 1.8e-4 < gmax 4.5e-4）。実 run では未観測 |
| U3-I6 | 変位後の再最適化が、手元の鞍点 DFT Hessian を使わず trust 0.1 で回る | efficiency | medium | `minimum.py:96-103`、`input.py:124`。P3a で片側 30 歩（77 s、74 s）。NWChem の `opt_drv.F` は極小化で負の固有値を |e| として扱うので、生の Hessian をそのまま渡せる（QRC と同じ） |
| U3-I7 | `_soften` は + 側 1 回の別分岐 | complexity | low | 片側で足りるのは化学的に正当（鞍点からはどちらでも下る）。問題は固定 0.1 Å と初期 Hessian なし |
| U3-I8 | 宣言端点が鞍点に収束すると親に basin が付かず、宣言反応が別の理由（endpoints_not_on_one_pes）で blocked_upstream になる | validity | high | P3a。`minima.py:129`、`hypotheses.py:120-126, 147-156`。鏡像対の例は U3-P2 で解消するが、鏡像でない 2 極小の間の対称入力では依然起きる |
| U3-I9 | rerank が DFT の SP 後にエネルギー窓を外し、DFT で高い構造にも opt と freq を払う | efficiency | medium | `selection.py:65`、`minima.py:204-209`。17 原子で 1 構造 opt 377 s + freq 898 s |
| U3-I12 | freq と opt が同じ電子状態かを確かめない | validity | medium | `gates.py:125-128` は構造の指紋と LOT だけ照合。UKS で別解に落ちた freq が黙って組み合わされうる（未観測、閉殻の差は 2.4e-9 Eh） |
| U3-I13 | 崩壊した鞍点を opt からやり直す | efficiency | low | `actions.py:347-348, 426-431`。節約できるのは freq を持つ崩壊鞍点の freq 1 本だけ（QRC 終点の再 opt は 1〜2 歩） |
| U3-I16 | 使われない `MinimumPolicy.amplitude_A` | complexity | low | `minimum.py:33`。`diagnostics.json` は mode-follow 経過の唯一の記録で、R5 の検証でも使ったので残す |
| U3-I17 | 開殻・電荷の golden、DFT の soft_minimum、負モード 2 本以上からの mode-follow、単原子、NWChem の 1 原子 opt/freq が未通過 | untested_path | medium | tests/golden/data/nwchem の 125 本はすべて電荷 0・一重項 |
| 他単位由来 | ECP・元素表・SCF 上限（U3-I14 → U0-I1、U1-I1、U0-I5）、状態ラベル（U3-I15 → U1-I7）、explore 生成物の未重複除去（U3-I10 → U4-I8）、screen の直列（U3-I11 → U9-I6） | | | |

**改良案**

#### U3-P1 単原子を特別な分岐なしで同じ経路に通す（熱化学を含む）（must、【複雑さ減】）

U0-P7、U1-P2、U7-P3 を統合。

- 変更: `Evidence.n_external` を `Literal[3, 5, 6]` にし、`frequencies=()` を許す。`xtb.py:187` の cast を削除。`relax_to_minimum` に単原子用の分岐は作らない（1 原子の opt と freq は数秒で終わり、Registry・known・MinimumRecord はそのまま使える）。**受け入れの最初に WSL で NWChem の 1 原子 optimize・frequencies（Cl⁻、F⁻）を確かめ**、driver が失敗する場合に限りエンジン側で「原子の opt は start=final の自明な Evidence」とする 1 か所の処理を足す。N=1 の化学種は explore の出発点（`explore.py:53-57`）から外す（原子から反応は探索できない。CREST は既存の小分子省略で素通り）。熱化学（U7 / U9-P8 の同一プロセス呼び出し）: N=1 では `goodvibes/engine.py:102` の hessian 必須と worker の S_rot 事前検査を外し、GoodVibes の単原子分岐（`zero_point_corr == 0.0`、振動なし）に入れる。GoodVibes 4.3.0 は空の rotemp で計算全体を飛ばすので（`thermo.py:685`）、空でないダミーの rotemp を渡す。L>0 の開殻原子の縮重度とスピン軌道は範囲外として U1-P9 の文書に書く。
- 理論的根拠: 原子には振動・回転の自由度がなく（3N−3=0）、分配関数は並進（Sackur–Tetrode）と電子だけで閉じる。停留点の判定は自明に成り立つ。
- 効果: X⁻+CH3Y の SN2、X•+RH の引き抜き、M⁺ 錯体の単量体が扱え、ΔG_assoc と `dG_act_vs_separated` が出る。run 全体が止まる欠陥がなくなる。
- コスト: 型 1 行、cast と事前検査の削除、explore の出発点の条件 1 行、テスト 2 件。
- 検証: (1) Ar の S°(298.15 K, 1 bar) = 154.846 J/mol/K（JANAF、https://janaf.nist.gov/tables/Ar-001.html。1 atm との差 R ln 1.01325 に注意）を QM なしで再現。(2) F⁻・Cl⁻・H•（二重項）の xTB と NWChem の opt/freq。(3) P2a の `sn2_cl` を再実行し、minima が通り `dG_act_vs_separated` が出る（PBE0/SVPD で約 −0.74 kcal/mol は参考値で合否基準にしない）。(4) F⁻·HF で ΔG_assoc。

#### U3-P2 basin の同一性に鏡映を含め、キラリティ（m=2）を熱化学に渡す（must、【複雑さ減】）

U7-P4 の m の部分を統合。縮退判定の掌性は U5-P6。

- 変更: identity の公開関数を 2 つに整理する。(a) basin の同一性 `same_basin`／`assign` は、1 組の基準（RMSD ≤ 0.05 Å、|ΔE| ≤ 5e-5 Eh、次点規則）で回転に鏡映（det=−1）も許す。`is_chiral(x)`（x と x·diag(−1,1,1) の置換不変 RMSD が 0.05 Å を超えるか）を返す。(b) ラベル付きの写像 `relabel`・`mapped_equivalent` は真の回転だけのまま（NH3 の反転のような縮退転位を恒等と区別するため）。`MinimumRecord` に `chiral: bool` を 1 つ持たせ、U7 はキラルな構造（極小・TS）の G に −RT ln 2 を入れる（アンサンブルでは鏡像対を 1 回だけ数え m=2 で表す）。`same_minimum` を削除し、`assign` はエネルギーで候補を先に絞る。`Registry.find` の既知判定（`minimum.py:191`、U6 の `_assign`）では template の charge・multiplicity から作る `composition_key` と `level_key` を必ず渡す。**宣言された鏡像化（R→S のラセミ化、H2O2 の gauche 鏡像化）は、U5-P6 で掌性を含めた縮退判定により degenerate（ΔG_rxn=0）として評価され、same_basin として捨てられない。**
- 理論的根拠: 気相のアキラルな環境では鏡像体のエネルギーと振動数は厳密に等しく、分配関数には縮重として入る（経路縮重度 L = σ_R m‡/(σ‡ m_R)、Fernández-Ramos et al., Theor. Chem. Acc. 118, 813 (2007), doi:10.1007/s00214-007-0328-0）。CREGEN の既定（enantio=.true.）とも一致する。DFT 極小の再現性は 0.002 Å 以下なので 0.05 Å は 1 つの基準として十分。
- 効果: 鏡像化の偽反応と reassigned がなくなる（P3a/P3b）。HONO の TS（C1、キラル）の ΔG‡ に RT ln2 = 0.41 kcal/mol が入る（12.23 → 11.82）。design.md:146 の「経路縮重度を含まない」の不正確な記述を「σ の比は種の σ を通して、m は chiral を通して入る」に直す。
- コスト: identity で約 +15/−15 行、record 1 フィールド、U7 で数行。basin 数が変わるので新しい run dir で取り直す。
- 検証: 単体: NH3 の反転は degenerate のまま、CHFClBr の R/S は same_basin で chiral=True、CH4 は chiral=False。WSL: `ch3o_doublet`（Cs 入力）で偽反応が消えて宣言反応が elementary、`p4_h2o2` の gauche 鏡像化が degenerate として評価される。HCN・HONO・NH3・水の ΔG‡ は HONO の RT ln2 以外不変。

#### U3-P3 mode-follow の変位をエネルギー目標にし、変位元の freq を初期 Hessian に渡す（should、【複雑さ減】）

- 変更: 変位の振幅は U6-P4 の共通関数（ν と質量からエネルギー目標で決める）を使い、`MinimumPolicy.amplitude_A` を削除する。各側の再最適化は `relax(coords, init=source.freq)` とし、変位元の freq Evidence の Hessian を**そのまま**渡す（NWChem は極小化で負の固有値を |e| として扱い、その方向の歩幅を trust の 3〜30% に制限する。QRC も生の TS Hessian で動いている）。trust は 0.3。saddle 級（虚振動 1 本）は今のとおり ± 両側で TS 候補を得る。soft 級は片側 1 回のまま（両側にすると雑音の ts_candidate と費用 2 倍）。**虚振動 2 本以上の構造は、全虚モードの合成方向へ片側に 1 回だけ変位して緩和し、ts_candidate は作らない**（高次鞍点から ± に下ると片側が 1 次鞍点になりやすく、TS 候補の条件を満たさず費用だけかかるため）。
- 理論的根拠: 最適化器が下り始めるかは期待されるエネルギー低下（½κs²）で決まる（QRC、Goodman & Silva, Tetrahedron Lett. 44, 8233 (2003)）。正確な Hessian からの準ニュートン法はほぼ 2 次収束。
- 効果: 柔らかいモードでも確実に下り始める。再最適化の歩数が減る見込み（P3a の片側 30 歩、初期 Hessian ありの QRC の実績から半分以下）。
- コスト: minimum.py で約 −20 行（振幅関数は U6-P4 と共有）。
- 検証: `ch3o_doublet`（Cs 入力）の歩数・時間・結果、重なり形 C2v DME（虚振動 2 本）が極小に落ちる、平面 NH3（D3h）から反転の ts_candidate、既知端点 4 系の不変。

#### U3-P4 鞍点に収束した宣言端点を、両側の極小のうち入力に近い方に登録する（should、【中立】）

- 変更: minima の register（`minima.py:122-135`）で、ts_candidate の親 species を mode-follow の両側 basin のうち入力構造との置換不変 RMSD が小さい方に join する（同程度ならエネルギーの低い方）。`MinimumRecord.notes` に `endpoint_was_saddle` を 1 語足す。U5 の見落としとして挙がった `hypotheses.endpoint_basin` 側の代替案はこれで不要。
- 理論的根拠: 利用者が意図した化学種は、入力構造がつながる basin。
- 効果: 対称な入力の宣言反応が、違う理由の blocked_upstream で黙って消えない。
- コスト: 約 6 行とテスト 1 件。
- 検証: fake エンジンで親が近い側に join。鏡像でない 2 極小の間の対称入力（対称に描いた H 結合錯体）の実計算。

#### U3-P5 崩壊した鞍点の freq を極小の freq として使う（could、【中立】）

U6-P9 を統合。

- 変更: `relax_to_minimum` が「収束した saddle の Evidence と、同じ構造・同じ Level で虚振動 0 本の VALIDATE_TS の freq」を受け取れば、opt と freq を省いて `is_minimum` → `Registry.add` に進む。`is_minimum` の `not_opt_task` は収束した saddle ジョブも受け付ける。QRC 終点は対象外（再 opt は 1〜2 歩で、freq はどのみち要る）。単量体への xTB 初期 Hessian の拡張は A/B の結果が出るまで現状（2 断片以上）のまま。
- 理論的根拠: 停留点の定義（勾配 0、Hessian 正定値）はドライバに依存しない。
- 効果: 崩壊鞍点 1 件で freq 1 本（17 原子で約 15 分）減。
- 検証: trans-HONO→HNO2（§6）でジョブ数と結果が不変。

#### U3-P6 rerank の後にもエネルギー窓を適用する（should、【中立】）

- 変更: `selection.rerank`（:57-65）に呼び出し側の `window_kcal` を渡す（inf をやめる）。always の候補は残す。explore 生成物の重複除去は U4-P7。
- 理論的根拠: 窓の目的は Boltzmann 分布に効く構造を残すこと。DFT//xTB は GFN2 より良い推定。
- 効果: DFT//xTB で 6 kcal/mol 以上高い構造（298 K の重み 4×10⁻⁵）に opt と freq を払わない（17 原子で 1 構造約 21 分）。
- コスト: 1 行とテスト 1 件。検証: v3 の tma_hf2_discover と amine パイロットで DFT opt 本数の減少と各 state の最安 basin の不変。

#### U3-P7 freq と親の電子状態の一致を既存のゲートで確かめる（should、【中立】）

U9-P13 を修正して統合。

- 変更: `gates._link_reasons`（:124-128）に、同じ幾何での |E_freq − E_parent| ≤ max(1e-6, 20·scf_tol) の 1 条件を足し、外れたら `state_mismatch`。`is_minimum` と `is_first_order_saddle` の両方に効く。
- 理論的根拠: 同じ構造・同じ電子状態なら、エネルギーは SCF の精度で一致しなければならない。UKS の解は初期推定に依存する。
- 効果: 開殻・電荷系で、別の SCF 解に落ちた freq の Hessian と熱補正が黙って組み合わされる事故をその場で検出する。
- コスト: 数行と単体テスト。検証: OH•・CH3O•・³O2・FHF⁻ の golden（§6）、閉殻の回帰 4 系で state_mismatch が出ない（実測差 2.4e-9 Eh）。

#### U3-P8 使われない MinimumPolicy を削除する（could、【複雑さ減】）

- 変更: `MinimumPolicy` を削除し、`relax_to_minimum` には `max_mode_follow` と gates の Policy だけを渡す（amplitude_A は U3-P3 で消える）。`diagnostics.json` と `_Run.history` は残す。
- 検証: 既存テストと回帰 4 系で minima の artifact 不変。

#### U3-P11 初期 Hessian のある opt で line search を切る（could、プローブ次第、【中立】）

- 変更: 初期 Hessian を渡す opt のデッキに `set driver:linopt 0` を 1 行足す。QRC で上り坂の軌跡が増えるなら QRC だけ除く。
- 理論的根拠: 正確な初期 Hessian では準ニュートンの歩がほぼ制限されず、line search は余分な 1 評価になる（NWChem マニュアル https://nwchemgit.github.io/Geometry-Optimization.html）。
- 効果: 実測の SCF 回数と勾配回数の差（17 原子で 12 対 9、HCN の QRC で 27 対 16）から、opt の SCF が約 2 割減る見込み。
- 検証: HCN・HONO・NH3 の QRC と TMA·(HF)₂ の錯体 opt で A/B。

他区分へ統合: 縮退判定を basin 構造で → U5-P6、explore 生成物の重複除去 → U4-P7、ECP・元素表・SCF・並列の失敗閉じ込め → U0-P1、U1-P1、U0-P5、U9-P1、検証 → §6。

**見送り**

- 同一性の基準を CREGEN の既定（0.125 Å、0.05 kcal/mol）に緩める: NWChem は xmax 1.8e-3 bohr を含む 4 条件をすべて満たすまで止まらず、同じ始点からの再現性は 0.002 Å 以下。緩めると錯体の浅い異性体を誤ってまとめる。
- 単原子を SP だけで登録する専用分岐: 型を広げるだけで同じ経路を通る。
- 鏡像を常に別種にし、立体指定時だけまとめる: 状態ラベルに立体を入れない方針と食い違い、縮重 m=2 で足りる。
- `mapped_equivalent` にも鏡映を許す: 縮退転位が恒等と区別できず消える。
- mode-follow と QRC に渡す Hessian の固有値を |λ| に書き換える（condition_hessian）: NWChem がすでに負の固有値を |e| として扱い歩幅を制限しており、QRC も生の Hessian で動いている。
- soft 級も ± 両側に変位: 雑音の ts_candidate と費用 2 倍。
- 数値 Hessian を記録する欄とゲート: U0 の原則と configs の固定で足りる。
- 一重項の全極小に BS-UKS の安定性解析: 大半は閉殻の安定な分子（TS は U6-P11）。
- DFT tier や screen tier の並列化: DFT は 4 rank で CPU と壁時計がほぼ等しい。screen の `relax_to_minimum` は `known=registry` に依存し、並列にすると決定性か freq 本数が変わる。
- 極小の熱補正に xTB Hessian を使う: LOT が混ざる。
- soft 級の追加検証段: 片側緩和と注記で足りる。
- freq を opt の movecs から始める（U3-P7(3)）: U3-P7 の一致条件で不一致が実際に出たときだけ、既存の `Setup.restart_vectors` で追加する。

**維持**: opt と別ジョブ freq の分離と `_link_reasons`、既知 basin での freq 省略と species 順の直列登録、置換不変 RMSD（Hungarian 法と複数初期向き）と次点規則（0.05 Å / 5e-5 Eh）、虚振動の 4 段階と極小側の閾値（noise 10、saddle 50 cm⁻¹）、mode-follow の ts_candidate の discovery 化と 4 分類、選択窓（組成 × state_label、discovery 端点と単量体最安の常時保持）、錯体の xTB 初期 Hessian と trust 0.3、Evidence の保証、UKS の解析 Hessian、Registry の責務（組成 × level_key の basin 集合）、mode-follow の経過を残す minima の `diagnostics.json`。


---

### 3.4 U4 反応探索（ReaDuct・trial 生成）

**現状**

- `stages/explore.py`（151 行）: 組成 × 状態ラベルごとに screen の最低 2 極小を出発点にし（`_sources` :53-57）、各 trial で NT2 → 生成物がなければ同じ drive で AFIR（:73-80）。窓は `Window(barrier_kj 150, reaction_kj 100)`（:38-41）、判定は `trials.product_verdict`（`trials.py:223-237`、gates の外）。source × trial × 機構の三重ループで直列。
- `chemistry/trials.py`（254 行）: 4 段の経験則を固定順に並べ先頭 10 件で切る。polar_h（labile H を受容原子へ、リレー優先）、h_shift（C を含む 1,2/1,3 移動）、heavy_bond（同じ断片内の形成・形成＋切断・単結合切断、芳香環の結合は除く）、association（断片対ごとに最近接 1 対）。重複除去は原子対集合の完全一致だけ（:210-212）。
- `backends/readuct/worker.py`: source 再最適化 → NT2 → Bofill TSOpt（automatic_mode_selection）→ 射影振動数で ν<−50 が 1 本 → IRC → 両端 opt（n_imag=0）→ 状態ラベルが source と一致する端を除き他方を生成物とする（:50-63、167-192）。AFIR は γ 125 → 300 kJ/mol。spin_mode は restricted_open_shell、電荷と多重度は source から。反復上限は ReaDuct の既定 150。
- 実績: 実 run の explore 5 本（v2/v3/w7 の amine_pilot2・tma_hf2、すべて中性・閉殻）で約 150 attempt、kept の生成物は 1 件（w7 tma、AFIR、ΔE_rxn −2.6 kJ/mol）。v3 tma の NT2 out_of_window 9 件は ΔE‡ がすべて 376.37175 kJ/mol（同じ反応を 9 回）。v2 amine の NT2 で得た TS（−1265i）は 3 回出て、すべて NH3·HF の二重 H 交換という縮退転位で same_as_source として捨てられた。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | ○ | NT2→Bofill→IRC→両端 opt→結合グラフ照合は ReaDuct/Chemoton の標準どおり。損失の主因は縮退転位を捨てること、drive の偏り、同値 drive の重複。1.45 帯によるプロトン移動生成物のラベル合流は主因ではない（分離したイオン対は今も別ラベル） |
| 汎用性 | △ | trial が極性 H と C–H 移動に偏り元素名の表に依存。Cl⁻···CH3Cl で SN2 の交換が出ず、CH4···OH で引き抜きが出ず、NH3···ICl は 1 断片と判定（本統合で `trials.generate` を実行して確認） |
| 有用性 | × | 実績として生成物がほぼ出ない（約 150 attempt で 1 件）。自動の配座対仮説をやめると（U5-P7）自動仮説の源はほぼ explore だけになる |
| 効率 | ○ | ReaDuct 自体は安い（31 attempt で 81 s、run の約 4%）。真の費用は下流で、重複生成物ごとの DFT opt と、仮説ごとの case |
| 簡潔さ | △ | 4 段の経験則、元素名の表 2 つ、自前の芳香環検出（約 30 行）、kJ の独自窓、既定と同じ knob、読まれないフィールド |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U4-I1 | drive 族が極性 H と C–H 移動に偏り、SN2・H 引き抜き・ラジカル付加・分子間交換を生成しない | generality | high | `trials.py:59-163`。swap と formation は同じ断片内だけ（:139、149-150）、断片間は最近接の重原子 1 対（:159-160）なので CH4···OH の O···H 引き抜きは原理的に出ない。Cl⁻···CH3Cl の「Cl→H 1,2 移動」は h_shift の near（:85）に C の隣の Cl が入るため |
| U4-I2 | 重複除去が原子対の完全一致だけで、対称等価な drive が上限を埋める | efficiency | medium | `trials.py:210-212`。v3 の同じ 376.37 kJ/mol × 9、v2 の −1265i × 3。TMA の先頭 10 件のうち固有は 2 件 |
| U4-I3 | 重原子の二次的接触（I···N 2.8 Å、Cl⁻···H–F 1.85 Å）が結合とみなされ 1 断片になる | validity | medium | `topology.py:79-89`。tma_hf2 の 2 端点が同じラベルなのは入力がどちらも共有プロトン型だから（N–H/H–F 1.317/1.097 と 1.271/1.135 Å、DFT でも 1 basin）で、プロトン移動の生成物が捨てられる主因ではない → U1-P4 |
| U4-I4 | IRC 端の opt が ReaDuct の既定の反復上限 150 で止まる。余分な負モードが残ると回復せず陰性 | validity | medium | `worker.py:159-160, 181-184`。amine の irc_end_not_minimum 8/9 件（6 件は上限停止、2 件は −170 cm⁻¹ の停留点、TD-08 の記述に依拠） |
| U4-I5 | 縮退転位（恒等 SN2、プロトン・HF 交換、リレー）を same_as_source として捨てる | generality | high | `worker.py:50-63` は置換不変ラベルで比べる。v2 の −1265i（形成 (0,4),(3,5)、切断 (0,3),(4,5)）と v3 の −289i は実在する素反応。実 run で得た in-window の TS はほぼすべてこの型 |
| U4-I6 | 閉殻一重項で単結合の切断（均一開裂）を drive として作る | validity | low | `trials.py:151`。閉殻の均一開裂はほぼ窓の外なので実害は無駄な attempt。電荷のある閉殻系の不均一開裂（R–OH2⁺→R⁺+H2O）は正当な反応 |
| U4-I7 | explore の窓が kJ/mol で DFT の窓より狭い（23.9 対 40 kcal/mol）、判定が gates の外 | validity | low | `explore.py:38-41`、`gates.py:38` |
| U4-I8 | explore の生成物を重複除去せず DFT opt に回す | efficiency | medium | `minima.py:176-196`。今は kept が 1 件なので理論上の費用。P1/P2/P5 後に顕在化 |
| U4-I11 | 読まれない `ReactionTrial.perturbed`、TS のある陰性でも障壁を残さない、誤った理由ラベル `monotonic_uphill`、既定と同じ knob `nt_total_force_norm` | complexity | low | `records.py:68`、`worker.py:89-92, 173` |
| U4-I12 | 電荷・開殻・重元素の explore が実計算で一度も通っていない | untested_path | medium | 実 run は 5 本とも中性・閉殻 |
| U4-I13 | AFIR へのフォールバックがほぼ無効で、NT2 が TS・IRC まで進んだ陰性でも AFIR を走らせる | efficiency | medium | 約 50 件の AFIR の大半が same_as_source か afir_not_converged、kept 1 件。`explore.py:77-80`（検証段の追加） |
| 重複 | 元素範囲・組成の多重度 → U1-I1・I4、並列 → U9-I6 | | | |

**改良案**

#### U4-P1 4 段の経験則を、元素に依らない結合変化の列挙器 1 つに置き換える（must、【複雑さ減】）

- 変更: `_polar_h`、`_h_shifts`、`_rings`、`aromatic_bonds`、`_breakable`、`_formations`、`_heavy_bonds`、`_associations` と定数群を削除し、次の列挙器 1 つにする。
  - 形成の候補対 (a, b): r_ab ≤ Σr_vdW、トポロジー距離 3 以上（断片の内外を問わない。b が切断を伴って移る場合だけ距離 2 の 1,2 移動を許す）、変化後の配位数が元素ごとの上限を超えない（H 1、B/C/N 4、O 3、ハロゲン 1、第 3 周期以降の非ハロゲン p ブロック 6、金属 9。U1-P1 の元素表にこの時点で列を足す。超原子価ヨウ素は対象外と文書化）。
  - T1 移動・置換: a–b を形成し b の既存結合 b–c（c≠a）を 1 本切る（H 移動、1,2/1,3 シフト、SN2、H 引き抜き、配位子交換）。同じ (a, b) に c の候補が複数あれば、角 a–b–c が直線に近い順（背面 SN2・直線型引き抜きを先頭に）。
  - T2 リレー: H の T1 を 2 つ連鎖（Grotthuss、多重プロトン移動）。
  - T3 形成: a–b だけを形成（付加、会合、環化）。
  - T4 切断: **多重度 > 1 または電荷 ≠ 0 の source でだけ作る**（β 開裂、イオンの不均一開裂）。閉殻・中性の均一開裂は作らない（RKS で記述できない）。
  - 並べ方: 上限 10 件の中で T1→T2→T3→T4 を巡回し、各テンプレート内は r_ab/Σr_cov の小さい順。`ReactionTrial.kind` を (transfer, relay, formation, dissociation) にする。`afir_pair` と NT2 の automatic_mode_selection はそのまま使える。
- 理論的根拠: 素反応の大半は自明でない最小の結合変化で表せる（YARP の b2f2、ZStruct の駆動座標と配位数の上下限、Chemoton の BondBased）。ReaDuct の NT2/AFIR は任意の原子対を元素に依らずに駆動する（ReaDuct マニュアル、https://scine.ethz.ch/download/readuct）。原子価の上限は Lewis 構造を作らずにありえない drive を除く最小の制約。
- 効果: SN2（T1: form C–Cl′、break C–Cl）、引き抜き（T1: form O–H、break C–H）、付加（T3）、β 開裂（T4）が生成される。「Cl→H 1,2 移動」は原子価と距離 2 の規則で消える。trials.py は約 254 行から約 150 行になり、元素名の表と芳香環検出が消える。
- コスト: 約 100 行新規、約 200 行削除、fixture 5〜6 件。U1-P1 と U4-P2 が前提（U1-P4 があれば断片判定が正しくなる）。
- 検証: 単体（純関数）: Cl⁻···CH3Cl に背面 SN2 の T1、CH4···OH と H···CH4 に引き抜きの T1、CH3O• に 1,2-H 移動と T4、閉殻 HONO に T4 なし、TMA·(HF)₂ の先頭に N···H–F の T1 と HF リレー。WSL（explore だけ、xTB、各数分）: Cl⁻+CH3Cl 錯体、OH+CH4 錯体、CH3O•、amine_pilot2 で、固有生成物数・陰性理由の分布・時間を今の版と比べる。受け入れ条件は SN2・引き抜きの生成物が得られることと、amine の −1265i の縮退交換が（U4-P5 と組んで）product になること（v3 tma は kept 0 件なので比較基準にしない）。

#### U4-P2 drive の同値類で重複を除く（must、【中立】）

C12 から診断ファイルを外した形。

- 変更: U1-P5 の `wl_classes` を使い、drive の鍵を「形成対と切断対それぞれの frozenset{(class_i, class_j)}」と「形成対の距離を 0.1 Å に丸めて並べたもの」の組にする（錯体の中でグラフ上同値でも幾何学的に異なる部位、N 側と F 側の HF を残すため）。source をまたいだ除去はしない（幾何が違う）。新しい診断ファイルは作らず、trial_id・スキーマ・JobStore の鍵は変えない。
- 理論的根拠: 分子グラフの自己同型で等価な原子に同じ駆動をかけても同じ反応になる。1-WL の分割は自己同型の軌道と同じか粗い。
- 効果: 「31 attempt 中 21 件が同じ移動」「376.37 kJ/mol × 9」「−1265i × 3」が解消し、上限 10 件が異なる化学の drive で埋まる。
- コスト: 15〜20 行とテスト 2 件。
- 検証: TMA の列挙 33 件がクラスで 5〜6 件、TMA·(HF)₂ で N 側と F 側の HF が別 drive。amine_pilot2 で −1265i の 3 件が 1 件になる。

#### U4-P5 縮退転位を発見として扱う（must、【中立】）

実 run の証拠が強いので should から must に上げました。

- 変更: worker で、source の結合集合と、source の原子添字で比べた端点の結合集合（添字付き）が違えば、状態ラベルが同じでも product とする。生成物の構造は原子の添字を保ったまま残す。`DiscoveryResult` に縮退フラグは持たせず、下流が basin の一致から判定する。U4-P7 の重複除去で縮退の生成物を source の basin に合流させない。`hypotheses._auto` は、発見の両端が同じ basin で端点の置換が恒等でない場合に、宣言反応と同じ degenerate の経路（U5-P6 を直したもの）で受け入れる。立体だけの変化は対象外のまま（`design.md:141`）。
- 理論的根拠: 縮退転位は同位体交換や NMR の線形変化で観測される素反応で、TS も障壁も実在する。反応の同一性は原子の対応を含む結合の変化で定義される。
- 効果: HF 交換、恒等 SN2、多重プロトンのリレーという、主対象の化学で中心となる素過程を自動で拾える（v2 の −1265i、v3 の −289i は実在する素反応）。
- コスト: worker 約 5 行、hypotheses 約 5 行（U5-P6 が前提）。
- 検証: amine_pilot2 の NH3·HF 二重 H 交換（−1265i）が product になり degenerate の case が走る。Cl⁻···CH3Cl で恒等 SN2、ギ酸二量体の二重プロトン移動（T2）が縮退として拾われる。

#### U4-P4 worker の数値処理を直す（should、【中立】）

- 変更: (a) 両端の opt（主な停止箇所）と IRC に `convergence_max_iterations=500` を渡す（C13）。(b) 極小側だけ、ν<−cutoff が残った端を最低のモードに沿って ± に変位し、それぞれ 1 回だけ opt し直して n_imag=0 になった側を採る。変位の振幅は U6-P4 の共通関数を使う。TSOpt 後の n_imag≥2 の再 TSOpt は入れない（実 run で発生 0 件、低レベル探索では陰性で十分、鞍点次数の処理は DFT 側の責務）。
- 理論的根拠: 停留点の分類は Hessian の慣性で決まり、余分な負モードに沿って対称性を破り最小化し直すのが定石（autodE）。ReaDuct の既定の反復上限は 150。
- 効果: amine の irc_end_not_minimum 8 件の多くが回復する見込みで、数値的な陰性と化学的な陰性が区別できる。
- コスト: (a) 1〜3 行、(b) 約 15 行とテスト 1 件。
- 検証: amine_pilot2 の irc_end_not_minimum 8 件を同じ入力で再実行し、回復件数と理由の変化を見る。

#### U4-P6 窓を 1 つにし（kcal/mol）、判定を gates へ移す（should、【複雑さ減】）

- 変更: `explore.Window` と `trials.product_verdict` の既定値を削除し、判定を `gates.discovery_verdict(policy, mechanism, ts_validated, dE_act_kcal, dE_rxn_kcal)` にする。ΔE_rxn の窓は `Policy.reaction_window_kcal`（40）を explore と hypotheses で共用し、低レベル側を DFT より狭くしない。低レベルの障壁の上限は gates のモジュール定数（例 50 kcal/mol）。`DiscoveryRecord` の単位を kcal/mol にし、過去 run の読み込み（ranking.csv と同様）に注意する。`nt_total_force_norm` を削除する（ReaDuct の既定と同値）。AFIR の γ は今の 2 値のまま。
- 理論的根拠: 低レベルのふるいが高レベルの基準より狭いと偽陰性を生む。
- 効果: GFN2 の誤差による偽陰性、単位の混在、gates 外のゲートがなくなる。knob 3 つ減。
- 検証: 窓の境界の単体テスト、U4-P1 の比較 run で out_of_window が減り DFT の窓外の生成物が増えないこと。

#### U4-P11 AFIR へのフォールバックを「NT2 が極大を見つけられなかった場合」に限る（should、【複雑さ減】）

検証段の追加。

- 変更: `explore.py:77-80` の条件を、NT2 の結果が `no_nt2_maximum`（U4-P9 で `monotonic_uphill` から改名）のときだけ AFIR を走らせるようにする。NT2 が TS・IRC まで進んだ陰性（same_as_source、irc_not_connected、irc_end_not_minimum）では AFIR を走らせない。
- 理論的根拠: NT2 が鞍点と IRC まで得た drive では、同じ drive の AFIR は同じ谷に落ちるだけ（v3 のリレーでは NT2 と AFIR がともに same_as_source）。
- 効果: attempt 数がほぼ半分になる（実 run の AFIR 約 50 件の大半が same_as_source か未収束）。
- 検証: amine_pilot2 と tma_hf2 で kept 生成物を失わずに attempt 数が減ること。

#### U4-P7 生成物を explore で一度だけ同一性判定し、固有な basin だけを species にする（should、【中立】）

U3-P6(2) を統合。

- 変更: kept の生成物を、同じ組成の screen 極小とそれまでの kept 生成物に対して、`state_label` の一致と置換不変 RMSD（エネルギーは使わない。ReaDuct の scine-xtb と screen の xtb 6.7.1 のエネルギーは 0.05 kcal/mol の精度で一致する保証がない）で判定する。既存 basin なら新しい SpeciesRecord を作らず、`DiscoveryRecord.product_species` に代表を入れる。縮退の生成物（U4-P5）は添字付きで比べ、source に合流させない。U4-P1・P2・P5 の後に入れる。
- 効果: 同じ生成物の DFT opt（17 原子で約 6 分／件）が trial の数だけ走らない。
- 検証: 同じ生成物を返す 2 trial で species が 1 つ。TMA·(HF)₂ の dft stage の opt 件数。

#### U4-P9 記録を整理する（could、【複雑さ減】）

- 変更: `ReactionTrial.perturbed` を削除、`result_dict` の `relative()` から「product のときだけ」の条件を外し TS のある陰性でも障壁を残す（C14′）、`monotonic_uphill` を `no_nt2_maximum` に改名、`relaxation_discoveries` を explore stage の関数に移す（trials は drive の列挙だけにする）。
- 効果: v2 の −1265i のような陰性を障壁の値で読める。
- 検証: 既存の単体テストと `tests/unit/backends/test_readuct_adapter.py` の更新。

他区分へ統合: 結合判定 → U1-P4、並列化 → U9-P6、検証セット → §6（amine_pilot2 の二重 H 交換を P2・P5 の受け入れ基準にする）。

**見送り**

- WBO で結合グラフ: 計算器・電荷・スピン依存（U1 と同じ判断）。
- YARP の Lewis 構造フィルタや b3f3 列挙の既定化: 組合せ爆発と新しい責務。
- Pauling 電気陰性度による極性の順位付け: 順位の経験則が増える。被覆不足が §6 で分かったら再検討。
- 生成物断片の電荷から開殻一重項の生成物を検出するゲート: T4 を作らない規則と窓で実質同じ効果。
- NT2 の下り坂を低レベルの障壁なし生成物として拾う: AFIR が下り坂の生成物を拾う（改名だけ）。
- 列挙数・固有数の diagnostics.json: 読むコードがない。
- xTB で組成のスピン状態を選ぶ: GFN はスピン非依存。
- 状態ラベルに立体: 不安定。
- source 再最適化のキャッシュ: xTB で 1 秒程度、基準エネルギーの一致に必要。
- Chemoton 全体・多段探索: 目的に対して過大。
- max_trials や sources_per_state を増やす: 同値除去と族の一般化なしでは重複と費用だけが増える。
- GFN-FF で探索: 結合の組み換えを表せない。
- TSOpt 後に n_imag≥2 の再 TSOpt（U4-P4(b) の後半）: 実 run で 0 件。
- AFIR の γ を障壁の上限から導く: 結合を増やすだけで根拠が弱い。

**維持**: NT2 → Bofill（automatic_mode_selection）→ 射影振動数 1 本 → IRC → 両端 opt → 結合グラフ照合の流れ、attempt ごとの JobStore ジョブと worker サブプロセス・版数 pin・continuation なし、SCC 失敗時の高電子温度での再試行と基準温度の SP、spin_mode の restricted_open_shell と電荷・多重度の引き継ぎ、直線分子の曲げと乱し（CH-30）、product/negative/failed の理由付き記録と negative を被覆報告にだけ使うこと、relaxation discovery、screen 極小を出発点にし上限 10 件、低レベル探索と DFT 検証の分業。

---

### 3.5 U5 仮説選択・障壁事前判定・経路（SCREEN / FIND_PATH）

**現状**

- 仮説（`chemistry/hypotheses.py`）: 宣言反応 → discovery → 同じ（組成、state_label、LOT）の DFT 極小の全組（conformer 源、:170-176）。最初の極小対が勝ち（`_lend_ts`）、組成あたり最大 6 件、窓 40 kcal/mol。degenerate・torsional・bond_changes・n_h_transferred は species の**入力座標**から計算（:134-143）。
- SCREEN（`drivers/reaction_case/actions.py:184-258`）: 近道は低レベル TS の xTB freq が負 1 本なら DFT SP 1 点で 3 点プロファイル。本流は DFT 両端を xTB で再最適化 → 写像 RMSD ≤ 0.05 Å（`_COLLAPSE_A`、:48）なら DFT IDPP 11 点（種なし）→ それ以外は pysis GS（DLC、分割系・直線系は Cartesian）+ rsirfo TSOpt → 全ノードの DFT SP（11 本の別ジョブ）。判定 `gates.barrier_verdict`（:168-212）は内部点の最大 − DFT 極小の高い方 ≥ 1 kcal/mol（または低レベル側 ≥ 1）なら proceed。below_zpe のために高い方の端点で DFT freq（`_tangent_mode_cm1`、:161-170）。
- FIND_PATH（:445-489）: NWChem ZTS 9 bead × maxiter 20 × 最大 3 チャンク、bead エネルギーが落ち着いたら止める。`profile.shape`（`profile.py:16-44`）で monotonic / single_max / multi_max。single_max は内部の argmax（`profile.py:66`）を種に。
- 通過状況: SCREEN の proceed、GS 失敗から FIND_PATH、IDPP 切り替え（H2O2、H+H2）は通った。行 12・13・14（障壁なしと中間体）は実 run で一度も通っていない。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | 峠の定理・minimax の上界の骨格は健全。形の判定が 2 系統で井戸を数えない。GS 経路（xTB 端点）からの barrierless は上界として厳密でない（IDPP 経路は両端が DFT 極小なので厳密。H2O2 の 1.28 は正しい上界）。縮退を入力座標で判定 |
| 汎用性 | △ | 電荷・多重度は正しく流れる。DLC の GS が決定的に落ちる（HONO）、xTB が DFT 極小を持たない系（H2O2 gauche、H+H2、イオン対）、傍観者が chord を支配（acac） |
| 有用性 | ○ | うまくいく系（HCN、NH3、DME、SN2、マロンアルデヒド）では SCREEN の種で FIND_PATH を省けている。障壁なしと中間体の出口は実質未使用（主因は試験系がないこと） |
| 効率 | △ | 読まれない below_zpe の DFT freq（17 原子で約 900 s）、落ち着き判定までの string 反復、GS の丸ごと再実行、11 本の独立 SP |
| 簡潔さ | △ | 判定器 2 つ、SCREEN 分岐 3 つ、pysis 側の座標系分岐 5 つ、読まれない BarrierVerdict フィールド 4 つ、conformer 分岐 |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U5-I1 | 経路の形を 2 つの別定義（`_max_rel_kcal` と `profile.shape`）で判定し、どちらも井戸を数えない。両端より深い中間体・吸熱の 2 段経路を BARRIERLESS で閉じる。種は検出した山ではなく内部の argmax | validity | high | `gates.py:168-212`、`profile.py:16-44, 66`、`state.py:184-200`。[0,−1.5,−3,−1,1,3,5] は両方 barrierless、[0,2,0.5,3,6,9,10] は single_max なのに SCREEN は barrierless で種は端点隣の添字 5（手計算で確認） |
| U5-I2 | GS 経路（xTB 端点から始まる）の barrierless が DFT 極小基準の上界として厳密でない。行 14 で閉じても record.barrier は SCREEN の値のまま | validity | medium | `gates.py:168-175`、`driver.py:129`。BarrierVerdict は reporting・summary・rankable から読まれず、記録の食い違いの影響は一貫性だけ |
| U5-I3 | 低レベルの端点の扱いが素朴。xTB が DFT 極小を持たない系で、種のない IDPP から最も高い DFT string へ落ちる | generality | high | `actions.py:203, 230-236`。P4b（IDPP の HEI は良い種なのに ZTS 100 s）、P3c（Cartesian IDPP が 2 本同時に切る +81 kcal/mol → ZTS の SCF 失敗） |
| U5-I4 | pysis GS が DLC の次元変化で決定的に失敗し、worker は同じ GS を丸ごと再実行する。engine は座標系分岐 5 つ、`initial_path` を無視 | efficiency | medium | v5 HONO の stderr（broadcast (55,) と (66,)）、`pysis/worker.py:27-34`、`engine.py:54-105, 188`。HONO は FIND_PATH 175 s（全体 332 s の 53%） |
| U5-I5 | 読まれない below_zpe のために SCREEN のたびに DFT freq を回す | over_diagnostics | medium | `actions.py:125-134, 161-170`、`records.py:104-108`。NH3 v5 の d383f050 は dft stage の freq と同じ Hessian の再計算 |
| U5-I6 | FIND_PATH は落ち着いた後も各チャンクを最後まで回す | efficiency | medium | HONO は 7 反復目で落ち着いた後も 20 反復・172 s。17 原子で 1 チャンク約 86〜108 分 |
| U5-I8 | 縮退判定を入力座標で行い、対称性が厳密でない入力の宣言反応を same_basin として黙って捨てる | validity | high | `hypotheses.py:135, 141`、P5a（0.0563 > 0.05、DFT 構造なら 0.0008） |
| U5-I9 | 自動の配座対仮説が全工程（最大 6 h）を走り、順位を汚す（discovery が先に枠を取るので化学反応を締め出しはしない） | efficiency | medium | `hypotheses.py:170-176, 192-195` |
| U5-I10 | walltime の行 4 が証拠で決まる完了（行 5・9・12・14）より前 | correctness | low | `state.py:217-221` → U9-P5 |
| U5-I11 | 経路ノードの DFT SP が 1 ノード 1 ジョブで原子推測から始まる。GS 経路の両端ノードの SP 2 本は捨てられ、IDPP 経路の両端は極小を再計算 | efficiency | low | `actions.py:156-158` |
| U5-I13 | IDPP の衝突判定 0.7 Å が元素に依らない絶対値 | generality | low | `interpolation.py:16, 73-80` |
| U5-I15 | 行 12・13・14、分割、井戸の経路が実計算で未通過（主因は障壁なし・井戸を持つ試験系がないこと） | untested_path | medium | → §6 |
| 重複 | 鞍点に収束した宣言端点 → U3-I8、種の接線 → U6-I2、予算 → U9-I5 | | | |

**改良案**

#### U5-P1 経路の分類器を 1 つにし（山と井戸を同じ関数で数える）、SCREEN と FIND_PATH で共用する（must、【複雑さ減】）

- 変更: `profile.py` に `classify(energies, resolution) -> Literal['barrierless', 'single', 'intermediate']` を置く。peaks = `interior_maxima(E, res)`、wells = `interior_maxima(−E, res)`。wells があれば intermediate（multi_max は必ず間に井戸を持つのでここに含まれる）、井戸がなく peaks があれば single、どちらもなければ barrierless。種は検出した山の添字で放物線補間する（`hei` に index 引数）。
  - 決定表: 行 12 と 14 を「最新のプロファイルが barrierless → BARRIERLESS」の 1 行、行 13 を「intermediate → VALIDATE_INTERMEDIATE（最も低い井戸のノード）」にする。**井戸の検証結果が same_as_endpoint だったときは「内部の最大 − 高い方の端点 < 解像度なら barrierless、そうでなければ最も高い山を種にする」**（今の multi_max の fallback の一般化）。
  - **SCREEN の barrierless で case を閉じてよいのは、経路の両端が DFT 極小の構造そのもの（IDPP 経路、U5-P2 の後は常に）の場合だけ。** GS 経路（xTB 端点）で出た barrierless は FIND_PATH に回す。
  - BarrierVerdict は verdict、source（screen/string）、max_rel_kcal、max_node_spacing_A（主張の解像度）、reasons だけにする。`finalize`（`driver.py:129`）には決め手のプロファイルの verdict を渡す。
  - 削除: `_max_rel_kcal`、`low_profile`・`low_endpoints` の OR（xTB だけの山による proceed）、`shape`、`max_rel_low_kcal`、`n_dft_points`、`seed`。C16 を包含。
- 理論的根拠: 2 つの厳密な極小を結ぶ連続経路の最大は峠のエネルギーの上界（minimax、Ambrosetti & Rabinowitz, J. Funct. Anal. 14, 349 (1973)）。解像度より深い内部の極小は別の極小の存在を示し、上界の議論では消せない。IRC の定義では多段の反応は TS ごとに素過程に分ける。
- 効果: 中間体を含む経路を BARRIERLESS で誤って閉じない。SCREEN の段階で中間体を見つけ分割へ進める。判定の定義が 1 つになる。関数 1 つ、決定表の行 1 つ、record のフィールド 4 つ減。
- コスト: profile に約 10 行、gates で約 −30 行。テストはモデル曲線 3 本と 3 点プロファイル。
- 検証: 単体: 上記モデル曲線が intermediate、3 点プロファイル、単調曲線 → barrierless。実計算: trans-HONO→HNO2（cis を中間体に持つ 2 段）で行 13 → 分割 → 子 case。(HF)2 の供与体–受容体入れ替え（障壁が解像度未満）で barrierless。HCN・HONO・NH3・水で不変。

#### U5-P3 below_zpe と `_tangent_mode_cm1` を削除する（must、【複雑さ減】）

U9-P3 の該当部分を統合。

- 変更: `actions.py:161-170` の `_tangent_mode_cm1`、Ctx.verdict の tangent_mode_cm1 と frames による ZPE 計算、`gates.py:201-211` の半量子の計算、`records.py:108` の below_zpe を削除する（約 −25 行）。ZPE で障壁が消える反応は U7-P1 で順位の量として扱う。
- 理論的根拠: ゼロ点補正後に障壁が残るかは、TS と極小の freq がそろう熱化学で物理的に定義された量として扱うべき。
- 効果: SCREEN ごとの DFT freq 1 本（17 原子で約 900 s、NH3 で 5 s）がなくなる。
- 検証: NH3 v5 の再計算で d383f050 相当の freq ジョブが消え、ΔG‡ 3.8395 のまま。

#### U5-P6 縮退と結合変化を、basin の最適化構造（各 species の原子順に並べ替えたもの、掌性を含む）で判定する（must、【中立】）

U3-P2 の縮退判定部分を統合。

- 変更: `hypotheses._record` で、`driver._endpoint`（`driver.py:66-74`）と同じ方法で basin の代表構造を各端点 species の原子順に並べ替え（U3-P2 で鏡映を許した場合は、鏡映したかどうか＝掌性も反映して並べ直す）、その 2 つに既存の `mapped_equivalent`（SE(3)、0.05 Å）を適用する。「2 つの置換が異なる」という判定は使わない（対称で等価な置換を偽の縮退とみなしうる）。`bond_changes`・`torsional`・`n_h_transferred` も同じ並べ替え済み座標で計算する。`_endpoint` は共有関数にする。
- 理論的根拠: 縮退転位の定義は「同じ極小が非自明な原子置換で自分自身に写る」ことで、最適化構造の上の置換で決まり、入力座標の距離閾値では決まらない（最適化構造では 0.0008 Å と 0.05 Å で明確に分離）。
- 効果: 手で描いた入力（メチル基、NH2 基、プロトン移動）の縮退反応が黙って捨てられない。宣言された鏡像化も degenerate として評価される。
- コスト: 約 +10/−5 行、テスト 3 件。
- 検証: P5a（`p5_dme`、粗い入力）が degenerate_rearrangement（ΔE‡ 約 2.3 kcal/mol）。NH3・SN2 が不変。H2O2 gauche 鏡像化の宣言が（U3-P2 の後）degenerate。

#### U5-P2 SCREEN の低レベル経路を「DFT 極小を固定端とする IDPP → xTB の climbing-image NEB」の 1 本にする（should、プローブで採否、【複雑さ減】）

- 変更: `PysisGrowingString` を PathEngine の本来の形にし、`initial_path`（hfauto の IDPP、sin(πt) の対称性破りを残す）を受け取って pysis で `{geom: {type: cart, fn: initial.trj}, cos: {type: neb, climb: true}, opt: {type: lbfgs, max_cycles: 上限}}` を回す（端は COS の既定 fix_first/fix_last で固定）。tsopt（rsprfo、xTB Hessian）は CI から try/except で 1 回だけ。TSOpt が例外を出したら pysisyphus の出力（current_geometries.trj か最後の cycle の trj）から COS の最終画像を読む（GS の丸ごとの再実行はしない）。未収束の NEB も経路として受け付ける。SCREEN は近道以外、この経路の**内部ノードだけ**に DFT SP を取り、端点は DFT 極小のエネルギーを使って U5-P1 で分類する（single なら山か xTB TS、intermediate なら井戸を種に。IDPP の HEI も同じ扱い）。FIND_PATH の初期経路にも同じ frames を使う。
  - 削除: `_screen_idpp` と `_screen_path` の二分岐、xTB による端点の再最適化、`_COLLAPSE_A` と写像 RMSD の判定、pysis engine の `_axis`・`_linear`・`off_axis`・`_split`・`_NEAR_LINEAR_A`・`_OFF_AXIS_A`・tric の分岐と DLC、worker の GS 再実行。
  - 採用の条件と後退案: プローブで悪化がないこと。悪化した場合は GS を残し、DLC 失敗時に Cartesian で 1 回だけやり直す。
- 理論的根拠: CI-NEB は与えた両端の間の最小エネルギー経路へ緩和し、climbing image は鞍点へ収束する（Henkelman, Uberuaga, Jónsson, J. Chem. Phys. 113, 9901 (2000)、https://pysisyphus.readthedocs.io/en/latest/）。両端が DFT 極小なので、経路上の DFT エネルギーの最大は（ノード解像度の範囲で）DFT の峠の上界になる。低レベル PES に端点の極小が要らない。
- 効果: xTB が DFT 極小を持たない系（H2O2 gauche、イオン対、H+H2）でも判定でき、偽の barrierless と種のない string への転落がなくなる。HONO の DLC 失敗による FIND_PATH（全体の 53%）を避けられる見込み。SCREEN の barrierless が常に厳密な上界になる。分岐と定数が 10 前後減る。
- コスト: pysis engine と worker で正味約 −60 行、actions で約 −30 行。xTB の NEB 1 本は数秒〜数十秒。プローブ約 1〜2 時間。
- 検証: HCN、HONO、NH3、H2O2、マロンアルデヒド、acac、DME、Cl⁻·CH3Cl、H+H2 か OH+H2 で、分類・種・時間・最終 outcome と ΔG‡ を v5 と比べる。採用基準はどの系でも同じか改善し、HONO と H2O2 で FIND_PATH を回避すること。

#### U5-P4 FIND_PATH は 1 回の呼び出しでチャンク 1 つだけ回し、すぐ分類する（should、【複雑さ減】）

- 変更: 1 チャンク（STRING_MAXITER は 20 のまま始め、短くするかはプローブで決める。IDPP から始める大きな系で 8 反復では HEI の質が落ちうる）だけ走らせ U5-P1 で分類する。single なら検出した山を種に REFINE_SADDLE へ。鞍点が失敗して種が尽き、`len(path_runs) < string_chunks` なら、行 16 の条件を「直前のプロファイルが single か intermediate で、種がない」にして前回の経路から次のチャンクを続ける。内側の for ループ、`energies_settled`、`_SETTLED`、end_drift の note、`PathProfile.energy_history`、`nwchem/output.string_path_energies` を削除する。
- 理論的根拠: 上界の議論は未収束の経路でも成り立ち、井戸は VALIDATE_INTERMEDIATE で確かめ、string の収束は種の質のためだけにある（その質は鞍点探索の成否で直接確かめられる。ZTS: E, Ren, Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007)）。
- 効果: 良い種なら string の後半を省ける（HONO で落ち着いた 7 反復目以降の 13 反復）。最悪の反復数は据え置き。
- コスト: actions で約 −15 行、profile・evidence・output で約 −20 行。
- 検証: HONO（GS 失敗経路のまま）、H2O2、acac で string 反復数、鞍点の成否、ΔG‡（差 0.01 kcal/mol 以内）、時間を v5 と比べる。

#### U5-P7 自動の配座対仮説をやめる（宣言したねじれ反応は残す）（should、【複雑さ減】）

- 変更: `hypotheses._candidates` の conformer 分岐、`_max_torsion_change`、`_dihedrals_deg`、`min_angle_deg`、`_auto` の `source=='conformer'` 条件を削除する（約 −30 行）。配座の寄与は U7-P1 の基準（同じ state_label の最小 G）で扱う。
- 理論的根拠: Curtin–Hammett の原理（Seeman, Chem. Rev. 83, 83 (1983), doi:10.1021/cr00054a001）: 速い配座の前平衡は律速過程の順位に寄与しない。
- 効果: 配座対 case の計算時間（1 件最大 6 h）と順位の汚染がなくなる（discovery が先に枠を取るので、化学反応を締め出してはいない）。
- 検証: 同じ state_label の DFT 極小 3 つから仮説 0 件、宣言した torsional は残る、HONO の cis/trans が不変。

#### U5-P9 経路上の DFT SP を 1 ジョブにまとめる（could、【中立】）

- 変更: エンジンに `energies_along(frames, method)` を足し、1 つのデッキに N 組の geometry と `task dft energy` を並べる（NWChem は同じ入力の後続タスクで前の movecs を既定の初期推定に使う）。鍵は frames の sha と手法。どこかで失敗したら unavailable。内部ノードだけに SP を取る変更は U5-P2 の中ですぐ行う。
- 理論的根拠: 前の点の MO から始めるのは断熱状態を追う標準的なやり方。
- 効果: 小分子で NWChem 起動が 11 回→1 回（HCN で約 15 s）、中分子で SCF 反復約 4 割減。開殻でノードごとに別の SCF 解に落ちて U5-P1 が偽の井戸を見るのを防ぐ。
- 検証: HCN と CH3O• の SCREEN で各ノードのエネルギーが個別ジョブと 1e-6 Eh 以内、⟨S²⟩ がそろう。

#### U5-P11 IDPP の衝突判定を r/Σr_cov にする（could、【中立】）

- 変更: `MIN_DISTANCE_A=0.7` を共有結合半径の和に対する比（< 0.5 で衝突）に替える。U1-P1 の後。
- 検証: H–H、C–C、I–I の対の単体テスト。

他区分へ統合: 種の接線 → U6-P2、決定表の行の順序と子 case の予算 → U9-P5、鞍点に収束した宣言端点 → U3-P4、検証セット → §6。

**見送り**

- 端点の曲率から隠れた障壁の上限（0.074κd²）を見積もり、区間に SP を追加: 極小どうしの barrierless は「解像度未満の障壁＝速い平衡」にしかならず、298 K では 1 と 3 kcal/mol で化学的結論は変わらない。解像度は `max_node_spacing_A` で示すだけ。
- SCREEN の基準を端のノードの DFT SP にする: 上界にならない。U5-P2 で両端を DFT 極小に固定すれば問題自体が消える。
- 二分子の SCREEN で結合に沿って拘束付きで押し込む（AFIR 相当）: explore の責務と重なる。固定端 CI-NEB で同じ目的を果たせる。
- SCREEN ノードの ⟨S²⟩ が外れた点をプロファイルから除く: 経路の連続性（上界の前提）が崩れる。TS のスピン汚染は VALIDATE_TS と熱化学が見る。
- 低レベルだけの山で proceed にする OR 判定: DFT の結論に効くのは DFT//低レベルだけ。
- string の収束を勾配（gmax）で判定: 傾いた経路で射影しない gmax は 0 にならず、結論は収束を要しない。
- SCREEN の上界で仮説を枝刈り: 上界は「低ければ有望」の片側しか言えない。
- IRC や mepgs を後戻り手段に: QRC の数倍かかり、分岐（VRI）も検出できない。
- NWChem string の mode parallel やノードの並列ジョブ（4 コア）: 4 rank の NWChem で CPU と壁時計がほぼ等しい。
- GS を残し DLC 失敗時に Cartesian で再試行を主案に: 対症療法。U5-P2 が不利な場合の後退案としてだけ残す。
- TS のエネルギーが上界 + 解像度を超えたら次の種: U6 の領域で効果も小さい。
- 分割した子に親の経路の区間を初期経路として渡す（U5-P10 の前半）: 分割は一度も通っていない経路で、配線が増えるのに SCREEN の再実行は安い。

**維持**: 峠の定理と minimax の上界の骨格、純関数の決定表と JobStore による再生、SCREEN を FIND_PATH の前に置く分業、低レベル TS の近道と `_xtb_saddle_ok`、`interior_maxima` のヒステリシス付き極大検出、hfauto の IDPP と sin(πt) の対称性破り、NWChem ZTS の後戻り手段としての位置付けと凍結端の継ぎ目処理、宣言座標による方向の指定、仮説の順序（宣言 → discovery → 最初の対が勝つ）と `_lend_ts`、電荷・多重度の受け渡し。

---

### 3.6 U6 鞍点精密化・TS 検証・接続（REFINE_SADDLE / VALIDATE_TS / CONNECT / 中間体・分割）

**現状**

- REFINE_SADDLE（`actions.py:287-310`）: 1 回目だけ種の位置で xTB Hessian を試し、「負モードがちょうど 1 本」かつ「方向ベクトルとの重なり ≥ 0.3」（`_first_hessian` :274-284）でなければ DFT Hessian を新しく計算。方向ベクトル `ctx.direction`（:140-145）は宣言座標の勾配か Seed.tangent（経路の接線か全原子の chord、:134-138）。`_mode_index`（:261-271）で方向と最も重なる負モードを選ぶ。`NWChemSaddle.refine`（`backends/nwchem/engine.py:353-366`）は種とまったく同じ構造の Hessian しか受け付けない。`_moddir`（:87-95）は負モードなしなら moddir 0、1 本なら autoz で moddir 1、2 本以上なら Cartesian で P·H·P の順位。trust 0.1、sadstp 0.1、inhess 2。maxiter で driver の Hessian を引き継いで継続（:181-205、LADDER GEOMETRY_MAXITER 2 回）。
- VALIDATE_TS（:313-340）: `is_first_order_saddle`（`gates.py:143-165`）は最低振動数 < −50（torsional は −20）、2 本目 < −50 なら higher_order、E_TS − max(宣言端点) < 2e-5 Eh なら low_prominence。higher_order は `imaginary_modes[1]` に沿って 0.1 Å（`_RETRY_A`）片側に押す。
- CONNECT（:385-415）: TS の Hessian を `np.load` し、`qrc_amplitude`（`chemistry/modes.py:41-62`）で振幅を決め、± に変位して TS の Hessian をそのまま init_hessian に opt（trust 0.3）。`gates.py:215-273` は no_initial_descent / trajectory_above_ts / no_descent と basin の集合で elementary / degenerate / reassigned。2 回目は振幅 2 倍。
- VALIDATE_INTERMEDIATE（:418-442）→ MULTI_STEP → `classification.split`（:83-102）で子反応。子は親の torsional・n_h_transferred を継承し、新しい 6 h を得る。
- NWChem の saddle の実際の動作（`opt_drv.F` の `driver_sad_search_dir` で検証段が確認）: moddir 0 では 1 歩目を大きな勾配成分の方向へ登る。2 歩目以降、負モードが 1 本ならそれを追い、2 本以上なら直前に追ったモードとの重なりで追跡する。追わない負モードには勾配の大きさに関係なく trust いっぱいの下り歩を取る。|e|<1e-4（smalleig）かつ |g|>gmax のモードには min(trust,|g|) の最急降下の歩を取る。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | 「虚振動 1 本 + QRC」の骨格は標準（Goodman & Silva 2003）。TS 固有でない判定（宣言端点基準の low_prominence、torsional で切り替える閾値）、直線的な QRC 変位が正しいねじれ型 TS を no_initial_descent で棄却 |
| 汎用性 | △ | 電荷 −1 と二重項は通った。傍観者の柔らかい自由度が chord を支配するとき、対称な種、柔らかい反応モードで壊れやすい |
| 有用性 | △ | 未通過の分岐は実行すると 3 例中 3 例失敗（DME、H2O2、acac）。うち H2O2 は正しい TS が手元にあるのに unresolved |
| 効率 | △ | 種の DFT Hessian（15 原子で 703 s）を再計算する経路が多い、停滞した saddle の継続（約 330 s）。QRC が paths の 51〜59% を占めるが、大半は平坦な錯体の尾の収束で、負の固有値による歩幅制限は SN2 で各 3 歩だけ |
| 簡潔さ | ○ | コードは小さく決定表は純関数。moddir の 3 分岐、継続時の moddir 抑制、定数（`_RETRY_A`、torsion_saddle_cm1、ts_prominence）、行 10、変位の規則 2 つ |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U6-I1 | 初期 Hessian と追うモードの指定が NWChem の実際の動作と噛み合わない。負モードが複数あると、追わない負モード（反応モードのこともある）が trust いっぱいで下らされジグザグする | correctness | high | T（acac の stdout に「Forcing downhill step in mode 1 eval=-1.0D-01 grad= 5.4D-05」、負 3 本のまま 36 歩で未収束）。moddir 0 は 1 歩目を勾配方向に登る |
| U6-I2 | 反応方向（全原子の chord）が傍観者の自由度に支配され、良い xTB Hessian を捨てて DFT Hessian を払い、誤ったモードを選ぶ。低レベル TS の種でも手元の虚モードを使わない | naive_solution | high | acac の実データ: 全原子 chord で PT モードとの重なりは xTB 0.223、DFT 0.225、メチル −300 は 0.963。bond_changes の原子 {2,6,10} と隣接に限ると 0.965/0.962/0.010。`actions.py:173-181, 194, 222, 282-284, 294` |
| U6-I3 | higher_order の再試行が DFT Hessian を新しく計算し（0.1 Å 隣の TS freq を捨てる）、試行の予算も共有する | naive_solution | medium | `actions.py:295-300, 333-337`、`engine.py:358`。押す方向（modes[1] の重なり 0.982）は失敗の原因ではない（modes[0] を押しても同じ停滞） |
| U6-I4 | TS の受理が TS 固有でない条件に依存（torsional による閾値の切り替え、宣言端点基準の low_prominence） | validity | medium | `gates.py:153, 161-164`、`actions.py:323`。E_TS < E_B の鞍点は A→I の正しい TS でありうる。実際の極小に対する降下は connection がすでに見ている |
| U6-I5 | QRC の no_initial_descent/trajectory_above_ts が、直線変位で TS より上から始まる正しいねじれ型 TS を棄却する。振幅 2 倍の再試行は逆効果 | validity | high | P4b。4 本とも 1 歩目から単調に下り gauche に収束、棄却理由は `gates.py:222-225` だけ |
| U6-I6 | QRC に負の固有値を含む TS Hessian を渡す | efficiency | low | SN2 で歩幅制限は 23・28 歩のうち各 3 歩。損失は数歩程度 |
| U6-I7 | 変位の規則が 2 つ（固定 0.1 Å と QRC のエネルギー目標） | complexity | low | `actions.py:50, 334`、`minimum.py:33`、`modes.py:41-62` |
| U6-I8 | 縮退の case で、宣言した結合変化が起きたかを確かめない（傍観者のメチル回転の TS でも両側が H を並べ替えた同じ極小なら degenerate として受理されうる） | validity | medium | `gates.py:231-247, 241`。reassigned は rankable で実際の minima を記録するので、配座違い・鏡像の害はほぼラベルだけ |
| U6-I9 | refine が intermediate を戻すので行 13 が再発火し、行 10 はその補正のためだけにある | complexity | low | `actions.py:290-291`、`state.py:174-177` |
| U6-I10 | saddle の maxiter で壊れた driver Hessian を引き継いで停滞を繰り返す | efficiency | medium | P5c: attempt_00〜02 が同じ tiny eigenvalue（−6.0e-5）の這いを繰り返す（1 件約 495 s のうち約 330 s）。過去 run の saddle 継続 4 件で収束 0 件 |
| U6-I11 | 柔らかい反応モードは NWChem の saddle で追いにくい（smalleig 1e-4）。症状と対処が文書にない | generality | medium | `opt_drv.F:3142-3143, 3218`、`design.md:144` |
| U6-I12 | 分割した子が親の torsional・n_h_transferred を継承し、新しい 6 h を得る | responsibility | medium | `classification.py:79-80`、`driver.py:106` → 予算は U9-P5 |
| U6-I13 | 一重項は常に RKS で、開殻一重項の TS の不安定性を検出しない | validity | low | `input.py:93-94`、`design.md:138`（文書どおりの既知の限界） |
| U6-I14 | higher_order の成功、collapsed → 分割、`_register` の新 basin、ねじれ QRC の成功が実計算で 0 件 | untested_path | medium | P4b、P5c |

**改良案**

#### U6-P2 反応方向を 1 つの規則にする（低レベル TS の虚モード → 反応中心に限った接線 → ねじれは二面角の勾配）（must、【中立】）

U5-P5 を統合。

- 変更: `ctx.direction`（`actions.py:140-145`、:194、:222）を置き換える。(1) 種が低レベル TS（screen_ts / discovery_ts）なら、`_xtb_saddle_ok` で計算済みの xTB freq の唯一の虚モードを Seed.tangent にする（bool の代わりにモードを返す）。(2) それ以外は経路の接線か chord のうち、`topology.bond_changes(R, P)` に現れる原子とその隣接原子の成分だけを残して正規化したもの。(3) 結合変化がないとき（ねじれ・反転）は、宣言座標があればその勾配、なければ差の最も大きい二面角の勾配。xTB Hessian の採否は U6-P1 の後は「方向との最大の重なり ≥ 0.3」の 1 条件だけにし、負モードの本数は問わない。Hessian のない新しい種では毎回 xTB を先に試す（`saddle_attempts==1` の条件をやめる。xTB freq は数秒）。
- 理論的根拠: 反応座標は結合が組み替わる原子の局所運動で主に決まる。P-RFO の収束域は追うモードが反応モードであることで決まる（Baker, J. Comput. Chem. 7, 385 (1986)）。pysisyphus の hessian_ref（参照虚モードとの最大重なりで root を選ぶ）と同じ考え方（https://pysisyphus.readthedocs.io/en/latest/tsoptimization.html）。
- 効果: acac で xTB Hessian が採用され DFT Hessian 1 本（703 s）が不要になり、追うモードも正しくなる。元素・反応型に依らず結合グラフだけで決まる。
- コスト: 約 10 行。U1-P4（集合差の bond_changes）があれば望ましい。
- 検証: acac の PT、HONO（ねじれ）、NH3（反転）、マロンアルデヒド、SN2 で選んだモードと chord の重なりを log に出す。xTB Hessian の採用率と DFT Hessian の本数を数える。

#### U6-P1 鞍点探索の初期 Hessian を「反応モードだけ負」に整える純関数 1 つにする（must、【複雑さ減】）

- 変更: `chemistry/vibrations.py` に `shape_hessian(H, coords, reaction_mode, floor=1e-3)`（約 15 行）を置く。質量加重しない P·H·P を対角化し、振動の固有ベクトルのうち reaction_mode と |cos| 最大の r を選び、H' = −max(|λ_r|, f)·u_r u_rᵀ + Σ_{i≠r} max(|λ_i|, f)·u_i u_iᵀ とする。NWChem adapter の saddle の prepare（`engine.py:153-157`）で hess_text の前にかける（job 鍵に mode の sha を入れる）。saddle は常に autoz + moddir 1（autoz 失敗時の Cartesian 継続は残す）。`NWChemSaddle.refine` の引数を mode_index から mode ベクトルに変える。
  - 削除: `_moddir`（`engine.py:87-95`）、`cartesian_mode_number`（`vibrations.py:104-113`）、payload・inputs の cartesian、render_saddle の moddir=0 分岐（`input.py:140`）、継続時の moddir 抑制（`engine.py:203-204`）、mode_index=None の経路と `mode_overlap_below_0.3` の注記（`actions.py:301-303`）、`_first_hessian` の「負モードがちょうど 1 本」の条件（:282-284）。
  - 最小化（QRC、mode-follow）に渡す Hessian は整えない（NWChem が負の固有値を |e| として扱い歩幅を制限するので、効果は数歩程度）。
- 理論的根拠: P-RFO／固有ベクトル追跡（Baker 1986、Bofill, J. Comput. Chem. 15, 1 (1994)）の収束域は初期 Hessian の慣性で決まる。合同変換で慣性は保たれるので（Sylvester）inhess 2 の内部座標でも負は 1 本のまま、moddir 1 = 唯一の負モード。pysisyphus の TS 初期 Hessian も接線方向の曲率を負に書き換える。
- 効果: moddir 0 の上り 1 歩目、負モード複数時の誤追跡と強制下り（acac）、Cartesian 分岐がまとめて消える。xTB Hessian の採用率が上がる。固有値の床 f で初期の数歩の smalleig 脱落を防ぐ（更新後にできる tiny eigenvalue は防げない → U6-P6）。
- コスト: 約 15 行追加、約 30 行削除。単体テスト 3 件（慣性、床、反応モードの選び方）。U6-P2 と一緒に入れる。
- 検証: HCN・HONO・NH3・SN2・CH3O• の既知端点 run で ΔG‡ 不変と saddle の歩数、NWChem 出力の 1 歩目が「negative= 1」。acac が PT の TS に収束（現状 60 分打ち切り）。DME の C2v 種（ho_harness）。

#### U6-P5 TS ゲートを TS 固有の条件だけにし、QRC ゲートを「降下 + 割付け」だけにする（must、【複雑さ減】）

- 変更: (a) `is_first_order_saddle(freq, *, saddle, policy)`: 最低振動数 < −noise_cm1 なら TS 候補（torsional に依らない 1 つの数値雑音の閾値）、2 本目 < −saddle_cm1 なら higher_order、その間は `soft_secondary_mode` の注記。`torsional`・`endpoint_energies` 引数、low_prominence、`Policy.torsion_saddle_cm1`・`ts_prominence_hartree` を削除する。**grid 雑音で出る 10〜50i の崩壊した鞍点が QRC に回る退行を避けるため、|ν| < saddle_cm1 で QRC が失敗した TS は collapsed として行 8 に合流させる（1 行）。** (b) `_side_reasons`（`gates.py:215-228`）から no_initial_descent と trajectory_above_ts を削除し、各側の**終点** < E_TS − drop（no_descent、traj[0] 基準の比較はやめる）と、両側が別々の終点に割り付けられること（`_assignment`）だけにする。TS を越えて戻った側は他方と同じ終点になり sides_same_basin で落ちる。振幅 ×2 の再試行は sides_same_basin のときだけ。
- 理論的根拠: 一次鞍点の定義は Hessian の負の固有値がちょうど 1 本で、大きさは問わない。接続の証明は変位からの降下（QRC/IRC の定義）で行う。峠の定理から、TS の高さを宣言端点と比べる根拠はない。
- 効果: A→I の正しい TS（吸熱の宣言で E_TS < E_B）を捨てず reassigned か MULTI_STEP として残す。ねじれ型 TS の QRC（H2O2）が通る。重い断片の再配向や緩い会合の本物の TS（|ν| 10〜50 cm⁻¹）が偽の中間体にならない。knob 2 つ、ゲート理由 3 つ減。
- コスト: 約 15 行削除、テストの差し替え（「E_B より低い TS が reassigned で残る」「変位点が TS より上でも両側が下れば接続」）。
- 検証: P4b（`p4_h2o2`、JobStore 再利用で数分）が通り ΔE‡ 約 1.15 kcal/mol（今の識別では elementary、U3-P2 の後は degenerate）。HCN・HONO・NH3・SN2 の回帰。水二量体の供与体–受容体交換（柔らかい TS）で QRC が偽の TS を落とすこと。

#### U6-P6 決定表の整理と、停滞した鞍点を新しい Hessian で 1 回だけ再開する（should、【複雑さ減】）

U9-P10 の saddle 部分を統合。

- 変更: (a) `validate_intermediate` が自分のきっかけ（ts_check='collapsed'）を消費して ts_check と last_saddle を None に戻し、refine_saddle は intermediate を戻さない（`actions.py:290-291`）。**行 8 の条件を `ts_check=='collapsed'` だけにする**（`intermediate is None` のままだと、multi_max の検証で same_as_endpoint になった後に崩壊した鞍点が検証されない）。行 10（`state.py:174-177`）を削除。(b) 行の順序は U9-P5。(c) `engine.continuation`（`engine.py:188`）で task.kind=='saddle' の GEOMETRY_MAXITER は drv.hess を引き継いで継続せず、最終フレームを種として先頭に戻す（試行の予算は共有）。次の REFINE_SADDLE は U6-P2 の規則で xTB Hessian を取り U6-P1 で整えて 1 回だけ再開する（新しい層は作らず、既存の Seed と refine の経路で済む）。saddle の TIMEOUT 継続は今のまま。
- 理論的根拠: 決定表の状態は一度判定した事実だけを持つべき。停滞の原因は更新 Hessian にできた tiny eigenvalue で、同じ Hessian で続けても直らず、新しい Hessian で直る見込みがある（`opt_drv.F` の saddle 分岐）。
- 効果: CaseState のリセット項目と行が減り、停滞した saddle の継続（1 件約 330 s）がなくなる。DME の C2v の失敗が回復する見込み。
- コスト: 約 10 行の変更、5 行の削除。
- 検証: 決定表の単体テスト（行の順序、multi_max で鞍点が失敗したとき同じ種を積まない、行 8 の発火）。ho_harness で停滞が 1 回（約 165 s）で返り、再開で一次鞍点（約 226i）に収束すること。

#### U6-P3 検証済みの TS freq を再試行の Hessian に使い、higher_order の押し方を単純にする（should、【複雑さ減】）

- 変更: (a) Seed に任意の `hessian: Evidence | None` を足し、初期 Hessian の順序を「近傍（`HESSIAN_NEAR_A` 0.5 Å 以内）の freq Evidence → 種で xTB → 種で DFT」にする。`NWChemSaddle.refine` の `_hessian_file` に `near_A=HESSIAN_NEAR_A` を渡す（`engine.py:358`）。(b) higher_order では、反応モード＝方向（U6-P2）と最も重なる負モードとし、**それ以外の負モードのうち最も負のものを、片側だけ、U6-P4 の振幅で押す**（対称な種では ± は等価）。その種に検証済みの TS freq を Seed.hessian として付け、U6-P1 で反応モードだけを負にして渡す。(c) `_RETRY_A` と `imaginary_modes[1]` への固定をやめる。
- 理論的根拠: P-RFO は追うモード以外のすべてを最小化するが、勾配 0 の対称方向には変位で対称性を破る必要がある（autodE と同じ）。inhess 2 は近い構造の Cartesian Hessian を変換して使える。
- 効果: 再試行 1 回ごとに DFT Hessian 1 本減（17 原子で約 15 分）。QRC だけにあった近傍 Hessian の特例が一般の規則になる。
- 検証: DME の C2v 種（ho_harness）で一次鞍点（約 226i、ΔE‡ 約 2.3 kcal/mol）。直線 H2O は縮退した 2 次鞍点で近くに一次鞍点がないので検証系にしない。

#### U6-P7 縮退の case で、宣言した結合変化が起きたかを確かめる（should、【中立】）

- 変更: 結合が変わる縮退反応では、`_assignment` の degenerate 分岐に「QRC 両側の原子番号付き結合グラフの組が {bonds(R), bonds(P)} と一致する」の 1 条件を足す。basin による割付けの本体は変えない（配座違い・鏡像の reassigned は rankable で実際の minima を記録しているので、書き換えは優先しない）。子の torsional は親から継承せず子の端点から決め直す（`classification.py:79`）。I→B の子反応を作る部分は could とし、§6 で reassigned が実際に起きてから検討する。
- 理論的根拠: 素過程は TS と、それが IRC でつながる 2 つの極小で定義され、化学種の区別は結合グラフで行う。
- 効果: 縮退 PT の case で傍観者のメチル回転の TS を degenerate として受理しない。
- 検証: acac の縮退 PT にメチル回転の TS を人工的に与える fixture で reassigned か failed、HONO と NH3 の回帰。

#### U6-P4 変位の規則を 1 つにする（could、【複雑さ減】）

- 変更: `modes.py` に `amplitude(nu_cm1, mode, symbols, target_hartree, bounds_A)`（κ = ω²·Σm_i|u_i|²、s = clip(sqrt(2·target/κ), bounds)）を置き、`qrc_amplitude` を置き換え、connect、higher_order の押し（U6-P3）、mode-follow（U3-P3）、ReaDuct の極小側の再最適化（U4-P4）で共通に使う。`amplitude_A`、`_RETRY_A`、connect の `np.load` を削除する。
- 理論的根拠: 調和近似で E = ½ω²Q²（正確な基準モードなら uᵀHu と一致）。QRC と pysisyphus の IRC の displ=energy と同じ考え方。
- 効果: 規則が 1 つになる。QRC の振幅は変わらず、H2O2 の QRC の失敗は直らない（直すのは U6-P5）。`_tangent_mode_cm1` の np.load は U5-P3 で消える。
- 検証: TMA·(HF)₂ の実 DFT Hessian で振幅と勾配が収束閾値を十分に上回る。HCN・NH3・SN2 の QRC の振幅が今と一致。

#### U6-P11 一重項の TS の RKS 不安定性を BS-UKS の SP 1 本で検出する（could、【複雑さ増】）

- 変更: VALIDATE_TS の後、多重度 1 の TS に限り NWChem の SP を 1 本（odft、mult 1、`vectors input <freq の movecs> swap beta h h+1`）。E_BS < E_RKS − max(1e-4 Eh, 20×scf_tol) かつ ⟨S²⟩ > 0.1 なら claim に `rks_unstable` を付け、既存の spin_contaminated と同じ blocker に合流させる。BS 最適化とスピン射影は範囲外と明記。BS 解が RKS に戻れば安全側の偽陰性。
- 理論的根拠: 閉殻 KS の一重項→三重項不安定性（Bauernschmitt & Ahlrichs, J. Chem. Phys. 104, 9047 (1996)）、BS 法（Noodleman 1981）。
- 効果: ビラジカル性の TS で RKS の偽の障壁を順位に載せない。費用は TS 1 点あたり SP 1 本。
- 検証: エチレンのねじれ TS（90°）で rks_unstable、閉殻の回帰系で立たない。

他区分へ統合: 子 case の予算と行の順序 → U9-P5、崩壊鞍点の freq 再利用 → U3-P5、検証セット → §6。

**見送り**

- QRC の後戻りに IRC（mepgs、EulerPC）: QRC の数倍かかり、分岐も 1 本しか示さない。
- 柔らかい TS のために pysisyphus の RS-P-RFO を i-PI で NWChem に駆動: 依存と実行形態が増える。§6 で収束率が低いと分かるまで導入しない。
- QRC の側を loose で止め割付けで打ち切る: identity の許容値と両立せず、新 basin では結局 tight な opt が要る。
- QRC の変位点で SP を取り振幅を半分にしていく: SP とループが増える。no_initial_descent の削除のほうが理論的に正しく計算も増えない。
- TS が経路の上界 + 解像度より高ければ別の峠として棄却: SCREEN の上界は離散的で信頼できない。
- low_prominence を注記として残す、|ν| の閾値を元素や質量で調整: TS の定義に含まれない量。
- 縮退で宣言した置換そのものを記録し照合: 結合グラフの一致で足りる。同値な別の縮退過程（1 メチル回転）は正当なより低い過程。
- 種の Hessian を常に DFT で計算: 15〜17 原子で 700〜900 s。
- saddle の既定を Cartesian: 常に autoz + moddir 1 にできる。
- BS-UKS の鞍点最適化とスピン射影を本体に: 複雑さと費用に見合わない。
- 高次鞍点を xTB TSOpt で種から作り直すことを既定の再試行に: 新しい経路が要る。まず U6-P3 と U6-P6 を試す。
- 最小化側（QRC・mode-follow）の初期 Hessian も正定値に整える（U6-P1 原案の後半）: NWChem がすでに |e| として扱い、損失は数歩。
- higher_order で反応モード以外の負モードを ± に押し SP 2 本で低い方を選ぶ（U6-P3 原案）: 押す方向の選択は DME の失敗の原因ではなかった。
- 接続の割付けを全面的に結合グラフ単位にし、配座違い・鏡像の reassigned を elementary に変える（U6-P7 原案）: 実害はラベルだけで、弱い錯体・共有 H で結合判定が境界になり割付けが不安定になる。

**維持**: 別ジョブ freq による TS 検証と `_link_reasons`、エネルギー目標の QRC と TS Hessian の QRC への再利用、xTB Hessian を先に試し DFT を後戻りにする分業、NWChem saddle の P-RFO + PSB 更新（線探索なし）と trust・sadstp 0.1、純関数の決定表と Evidence・ゲート・Action の分離、collapsed か multi_max → validate_intermediate → split と深さの上限、ねじれの case での周期的二面角による割付けのフォールバック（`actions.py:377-379`）、観測した Level での電荷・多重度の照合、TS freq への spin_ok、縮退の case での `mapped_equivalent` による両側の区別。


---

### 3.7 U7 熱化学・速度論

**現状**

- 構成: `stages/thermochemistry.py`（254 行）、`chemistry/thermo.py`（131 行、純関数）、`backends/goodvibes/engine.py`（126 行）と `worker.py`（108 行）、`gates.thermo_consistent`（`gates.py:276-294`）、`gates.rankable` の dzpe 部分（:297-317）、Policy の 4 knob（:39-42）、`reporting/summary.rank_rows`（:87-136）。
- 流れ: DFT 極小と `SaddleClaim.freq_calc` を subject にし、振動数（極小は |ν|、TS は最低モードを除く。`thermo.py:29-33`）・質量（H–Kr と I だけ）・回転定数を GoodVibes 4.3.0 の API に子プロセスで渡し、JobStore にキャッシュ（`engine.py:93-126`）。QH=False、Grimme の qRRHO（cutoff 100）、symm=True（pymsym の σ）、S_el = R ln(2S+1)。qs × cutoff の 6 変種で感度の幅。実行時ゲート `thermo_consistent` が ZPE・E・n_real を照合。`energy_method` があれば composite。
- 反応: `dG_act = G_TS − G(minima[0])`（`thermo.py:121-131`、基準は 1 構造）。1 atm から Δn で 1 bar / 1 M に換算。会合（`thermochemistry.py:136-176`）は単量体を (Hill, q, m) で束ね、`dG_assoc` と `dG_act_vs_separated` を表示用に出す。Boltzmann の population は誰も読まない（:123-133）。
- 補足: どの反応も係数 1 の組成 1 つから組成 1 つへの変化なので（`hypotheses.py:93`、`classification.py:75`）、反応の Δn は常に 0 で、`dG_act`・`dG_rxn` は標準状態に依存しない。標準状態が意味を持つのは会合量だけ。
- v5 で thermo stage は 0〜1 s。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | qRRHO の式、GoodVibes への渡し方、S_el、Δn 換算、σ の比の自動算入（NH3 で RT ln2 = 0.41）、composite は正しい。順位の量が 1 配座基準で、ZPE で沈んだ障壁の値に意味がなく（序数としてのマロンアルデヒド 1 位は化学的に妥当）、障壁なしが順位外、キラリティ m がない |
| 汎用性 | △（単原子は ×） | 電荷と 2S+1 は正しく流れる。単原子が 4 か所で止まる（→ U3-P1）、原子番号が dict の挿入順 |
| 有用性 | △ | 順位の値が「その配座から見た TST の ΔG‡」で、分離基準・障壁なし・トンネルを含まない |
| 効率 | ◎（計算）/ ○（仕組み） | 計算は 1 秒未満。閉じた式に子プロセス・JobStore・版数 pin・thread_map が付いている |
| 簡潔さ | △ | 教科書の式の外部エンジン化、判定が stage・gates・summary の 3 か所、書くだけの値、効かない分岐、dzpe の連鎖 |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U7-I1 | 順位量の基準が単一配座（minima[0]）で、同じ反応物状態の共通基準になっていない。単量体アンサンブルが異性体を混ぜる | validity | high | `thermo.py:121-131`、`thermochemistry.py:145, 166, 170`、`design.md:146`。explore と配座対仮説は別々の配座から反応を作る |
| U7-I2 | ZPE で沈んだ障壁（ΔE0‡ ≤ 0）、ΔG‡ < 0、ΔG‡ < ΔG_rxn の値を順位の量として定義せず注意書きで扱う。最速過程の BARRIERLESS が順位外 | validity | medium | `gates.py:52-54, 297-317`、`design.md:150`、P4a |
| U7-I4 | σ を pymsym の既定閾値で判定し、失敗すると黙って C1（σ=1）になる | correctness | medium | pymsym `high_level.py:105-108`、GoodVibes `thermo.py:869-900`。プローブでは 1e-3 Å の揺らぎで CH4・NH3 が C1。実 run での誤りは未確認（NH3 は正しく、TMA·(HF)₂ は Cs でも C1 でも σ=1） |
| U7-I5 | キラリティ（光学異性体数 m）の縮重がない。`design.md:146` の「対称反応の因子は含まない」は不正確（σ の比は入っている） | validity | medium | HONO の TS は C1 でキラル（鏡像 RMSD 0.496 Å）、ΔG‡ 12.23 → 11.82 になるべき → U3-P2 |
| U7-I6 | トンネル補正がない | validity | medium | 非対称 Eckart（298 K）で HCN κ=5.4（0.99 kcal/mol）、NH3 1.84、HONO 1.65。典型的な ν‡ で 0.3〜1 kcal/mol。2000i のモデルの κ=421 は 1 次元 Eckart の信頼域外（T_c ≈ 458 K） |
| U7-I7 | dzpe_out_of_tolerance の blocker とその連鎖は、構造的にもう起きない旧バグへの防御で、n_H=0 でも 5 kcal/mol の根拠のない上限を課す | over_diagnostics | medium | `gates.py:41-42, 310-316`、`records.py:155`、`hypotheses.py:143`、`classification.py:80`、`topology.py:165`。二重計上は `thermo_consistent` の zpe_mismatch（1e-5 Eh）がすでに検出 |
| U7-I8 | 1 秒未満の閉じた式に外部エンジンの仕組み（子プロセス、JobStore、版数 pin、実行時照合、S_rot 事前検査、knob 2 つ） | complexity | medium | `goodvibes/engine.py`、`worker.py:1-23, 94-97`、`gates.py:276-294` |
| U7-I9 | 順位の判定が stage・gates・summary の 3 か所で重複 | responsibility | low | `thermochemistry.py:179-187`、`gates.py:305-310`、`summary.py:139-152` → U8-P6 |
| U7-I10 | 書くだけの値（population、`ThermoResult.S_rot`、常に空の `im_frequency_wn`）と効かない分岐（縮退で 0 に強制、`thermo.py:123-124`、`thermochemistry.py:206-207`） | complexity | low | grep で読み手 0 件 |
| U7-I11 | 感度の幅が dG_act にしか付かない（剛直な系で約 5e-4、緩い錯体 SN2 で約 1.05 kcal/mol）。支配的な不確かさは手法誤差 | validity | low | `thermochemistry.py:204`。文書に 1 行 |
| U7-I12 | 電子分配関数が 2S+1 だけ | generality | low | → U1-P9 の文書 |
| U7-I13 | 複数温度、1 M、単原子、σ・m、トンネル、障壁なしの熱化学の検証がない | untested_path | medium | → §6 |
| U7-I14 | パネルの順位に thermo が反映されない。composite の SP 層のスピン状態を確かめない。thermo の artifact id に stage id も energy layer も入らず、追記した composite が元の値を黙って上書きし、どの層で順位を付けたかが残らない | validity | medium | `thermochemistry.py:241-249`、`manifest.py:80-81`（後の同じ id が勝つ）、P1 → U8-P3、U8-P8 |
| 重複 | 単原子 → U3-I1（止まる箇所: `evidence.py:85`、`goodvibes/engine.py:102`、`worker.py:94-97`、GoodVibes `thermo.py:685`） | | | |

**改良案**

#### U7-P1 順位の量 δG_eff を thermo で 1 つだけ定義し、障壁なし・沈んだ障壁も順位に入れる。dzpe の連鎖を削除する（must、【複雑さ減】）

U7-P2、U8-P1、U8-P5 を統合。

- 変更:
  - `chemistry/thermo.py` に `effective_barrier(G_ts, G_R, G_P)` を 1 つ置く（約 10 行）。**δG_eff = max(G_TS, G_R, G_P) − G_R**。G_R・G_P は、反応物・生成物と**同じ state_label・同じ freq LOT・同じ energy LOT の DFT 極小の最小 G**（アンサンブルの G は採らない。DFT に流す配座は窓選択で最大 3 件ほどで、補正 ≤ RT ln n は標本不足で手法誤差より小さい）。
  - **順方向か逆方向の ΔE0‡（ΔE‡ + ΔZPE）が 0 以下なら、障壁なしと同じ扱い**（δG_eff = max(G_R, G_P) − G_R）にし、note `submerged_barrier` を付ける（blocker ではない）。BARRIERLESS は G_TS の項がなく max(ΔG_rxn, 0)。`_RANKABLE_OUTCOMES` に BARRIERLESS を加える（「解像度未満の障壁の下限値」であることは outcome 列で読める）。
  - G_sep（分離反応物）は基準に混ぜない。`dG_act_vs_separated` は今の別列のまま（反応はすべて Δn=0 で順位量は標準状態に依らないが、G_sep を混ぜると 1 atm と 1 M で基準が入れ替わり、単分子と二分子の次元の違う量が 1 列に並ぶため）。単量体アンサンブルは (Hill, q, m) ではなく state_label で束ね、異性体を混ぜない。
  - κ（トンネル）と σ の自前判定は順位に入れない。m=2 は U3-P2 の chiral から G に入る。
  - `ReactionThermo` に `dG_eff_kcal` を加えて順位の量にする。`dzpe_act_kcal` は ΔE0‡ の判定に使うので残す。`dG_act_kcal` はその配座から見た TST の値として表示用に残す。ranking.csv に `dG_eff`・`dG_act`・`dG_rxn`・`dG_act_vs_separated`・`torsional` の列を並べる（`torsional` は、mode-follow から来たねじれの段（P3a の ΔG‡ 3.46）を化学反応と見分けるため。宣言したねじれ反応は順位に残す）。新しい列には既定値を付ける（extra='forbid' の下で既存 run を読むため）。
  - 削除: `gates.rankable` の dzpe 部分、`Policy.rank_dzpe_base_kcal`・`rank_dzpe_per_h_kcal`、`ReactionRecord.n_h_transferred` と `hypotheses.py:143` の設定・`classification.py:80` の継承、`topology.transferred_hydrogens`・`proton_coordinate`、`ReportConfig.metric` と `_value`、`SpeciesThermo.population`・`_populated`・`boltzmann_populations`、縮退で 0 に強制する分岐 2 か所、`design.md:150` の注意書き。
- 理論的根拠: 微視的可逆性: TS の G が反応物・生成物の G より低ければ（調和・熱補正で障壁が沈んだ結果）ボトルネックではなく、振動断熱障壁が 0 以下なら鞍点の TST は意味を持たない（Truhlar, Garrett, Klippenstein, J. Phys. Chem. 100, 12771 (1996), doi:10.1021/jp953748q）。Curtin–Hammett: 速い配座の前平衡の下では反応物は状態として 1 つ。ΔZPE‡ の大きさには理論上の上限がない。
- 効果: 順位が 1 つの定義にそろい、反応ごとの基準のばらつき（+3 kcal/mol の配座から出た反応が見かけで上位になる）が消える。最速過程（障壁なし・沈んだ障壁）が同じ量の上に並ぶ。経験則の blocker が消える。新しい計算ジョブは要らない。
- コスト: 約 +20 行、約 −80 行（dzpe の連鎖、population、0 強制、metric）、knob 3 つ減。
- 検証: 単体（QM なし）: 通常、ΔG‡<0、ΔG‡<ΔG_rxn、ΔE0‡ の逆方向 ≤ 0、障壁なし、複数配座で最小 G が基準、ranking の列。実データ: レコードのフィールドを削除すると既存 run の manifest が extra='forbid' で読めないので、**report の追記ではなく、複製した run を `--from reaction-paths` で再実行（JobStore に当たるので安い）**して確かめる。HCN・HONO・NH3（各状態 1 配座で δG_eff = dG_act、回帰）、p4_malon（submerged_barrier が付き値は max 規則、序数は 1 位のまま）、ch3o_doublet（torsional 列）、TMA·(HF)₂（複数配座）。

#### U7-P5 GoodVibes を同一プロセスの純関数として呼び、外部エンジンの足場を外す（should、【複雑さ減】）

U9-P8 を修正して統合（qRRHO の自前実装は見送り）。

- 変更: GoodVibes の計算を同じプロセスで呼ぶ純関数 `species_thermo(...)` にし、子プロセスの worker、JobStore の保存、`ThermoBatch`、`ThermoEngine`/`ThermoResult`/`Capability.THERMO` の登録、実行時の `thermo_consistent` と `Policy.thermo_zpe_tol`・`thermo_energy_tol`、S_rot の事前検査、挿入順の原子番号表、pipeline の `engine: goodvibes` と site の goodvibes 項目を削除する。版数の固定は preflight で確かめる。単原子は U3-P1 の分岐（`zero_point_corr == 0.0`、非空のダミー rotemp）。`ThermoSettings.symmetry` は σ を切る理由がないので削除。
- 理論的根拠: qRRHO は閉じた式（Grimme, Chem. Eur. J. 18, 9955 (2012), doi:10.1002/chem.201200497）。adapter の仕組み（版数 pin、キャッシュ、失敗の分類、子プロセス）は外部実行ファイルを呼ぶ重い計算のためのもの。不変条件は実行時照合ではなく版数を固定した golden で保証する。確立したライブラリを使い続けるほうが、書き直すより検証対象が少ない。
- 効果: 正味約 −200 行、knob 2 つと Capability 1 つ減。thermo stage のジョブ・キャッシュ・失敗種別が消え、責務が「閉じた式の計算」だけになる。
- 検証: `tests/golden/test_goodvibes_thermo.py` を、同一プロセス版と G06 の CSV（HCN、HNC、TS）の一致（1e-8 Eh）にし、二重項（R ln2）、直線分子、単原子 Ar（S° 154.846 J/mol/K、1 bar）、truhlar 変種、cutoff 50/150 の fixture を足す。v5 の 4 系で ranking が 1e-6 kcal/mol 以内で一致。

#### U7-P6 非対称 Eckart のトンネル係数 κ を列として出す（could、【複雑さ増】）

- 変更: `kappa_eckart(nu_imag_cm1, V0f, V0r, T)`（約 25 行、Johnston–Heicklen の数値積分、Arkane と照合）を置き、V0f・V0r（エネルギー層の ΔE + freq の ΔZPE、両方正のときだけ）から κ を ReactionThermo と report に列で出す。T < T_c = hcν‡/(2πk_B) なら note `deep_tunneling`。**δG_eff（順位）には入れない**（T_c 付近以下では 1 次元モデルの誤差が補正と同程度で、序数をモデル誤差で並べ替えてしまう）。モデルは Eckart 1 つだけ。
- 理論的根拠: Eckart, Phys. Rev. 35, 1303 (1930)、Johnston & Heicklen, J. Phys. Chem. 66, 532 (1962)。
- 効果: H・プロトン移動の反応でトンネルの大きさを利用者が読める。
- 検証: ν‡→0 で κ→1、T≫T_c で Wigner と一致、対称 Eckart の解析解、Arkane の参照値。

他区分へ統合: 単原子 → U3-P1、m=2 → U3-P2、判定の一本化 → U8-P6、電子縮重の文書 → U1-P9、パネル・composite・エネルギー層の来歴 → U8-P3、検証 → §6。

**見送り**

- マスター方程式・微速度論、ネットワーク全体の energetic span: 探索・順位付けの目的外。
- 変分 TST・SCT: 経路の Hessian 列が要り DFT の計算量が桁で増える。
- Wigner・Skodje–Truhlar・Eckart を選べる設定: knob が増え、Wigner は過小、ST は T<T_c で発散的に過大。
- 沈んだ障壁を blocker にして順位から外す: 実際には最速の過程で、表から消すのは有用性を損なう。
- barrierless の両端を 1 つのアンサンブルにまとめて熱化学へ: 極小の集合の意味が変わる。max(ΔG_rxn, 0) で足りる。
- G_ref = min(G_R, G_sep)（分離反応物を順位の基準に混ぜる）: 標準状態で基準が入れ替わり、次元の違う量が 1 列に並ぶ。
- 基準を同じ状態の ensemble_G にする: 標本不足で補正が手法誤差より小さく、min の方が単純で頑健。
- pymsym の閾値を緩める: 3e-3 Å までしか頑健でなく identity の許容値とも一致しない。σ の自前判定（Kabsch + Hungarian）は実 run で σ を誤った例がまだないので、GoodVibes を外すときまで保留。
- qRRHO を hfauto で書き直す（U7-P5 原案）: 電子エントロピー・並進・σ の扱いの利点を捨て、検証対象が増える。GoodVibes にない機能が必要になったときに再検討。
- 配座変化の反応を別の順位表: 配座対仮説の廃止（U5-P7）と torsional 列で足りる。
- 1D 束縛回転子、Head-Gordon の QH エンタルピーを主設定: 計算・実装が大きい／評価が定まっていない。
- 最終の ΔG‡ に xTB の熱補正: LOT が混ざる。
- 実行時の GoodVibes 照合ゲートを残す: 守る失敗が構造的に起きない。
- population を report に表示: 順位にも判定にも使わない。
- 元素ごとの SO 補正を自動の表で／電子縮重とスピン軌道の入力欄（U7-P9 原案）: 分子に一般化せず、対象系がまだない。
- log10 k 列とチャネル単位の TS アンサンブル（U7-P10）: δG_eff の単調変換で情報が増えず、グルーピング規則が増える。
- スケール因子の既定を 0.985 に: ΔG‡ は 0.06 kcal/mol 以下しか変わらない。

**維持**: composite `G = E_SP + (G_thermo − E_thermo)`（同じ構造の SP に限り、欠ければ energy_layer_missing で fail-closed）、`thermo_frequencies` の固定規則、1 atm で計算し Δn で換算する `standard_state_shift`、H は RRHO・S は Grimme の qRRHO（cutoff 100）、vib_scale 1 knob（既定 1.0）、S_el = R ln(2S+1)、qs × cutoff の感度幅と `_tie_ranks`（対象を δG_eff に移す）、`ensemble_G` と `association` の純関数、錯体と単量体の LOT 判定の `same_pes(state=False)` と mixed_level_of_theory の fail-closed、SpeciesThermo と ReactionThermo を artifact にする構造、GoodVibes 4.3.0 の G06 出力を golden に残すこと。

---

### 3.8 U8 一点計算・手法パネル・順位付けと報告

**現状**

- sp stage（`stages/single_point.py:30-78`）: freq の最終構造で methods の手法ごとに一点計算を直列に回す。targets は reaction_stationary_points（outcome を問わず全反応の極小と TS）か all_minima。parents に subject を記録。開殻の WFT は input_invalid。
- report（`stages/report.py`、`reporting/summary.py`）: `method_panel`（:223-253）は Level.full_key ごとに ΔE_rxn・ΔE‡ を並べ、ΔE_rxn が +1 超と −1 未満の手法が両方あれば sign_disagreement。`participant_notes`（:139-152）が notes に method_sign_disagreement を足し、`gates.rankable` の 4 つの blocker（thermo_unavailable、mixed_level_of_theory、spin_contaminated、method_sign_disagreement）と dzpe で順位を止める。`rank_rows`（:87-136、CC 15）は rankable な反応を dG_act で並べ、感度幅が重なれば推移閉包で同順位。出力は ranking.csv、coverage.csv、method_panel.csv、report.html（パネルは HTML に出ない）。
- (T, 標準状態) は ReportConfig か、その pipeline 自身の conditions から取る。method_panel.yaml は元の run の conditions を手で写す必要がある（:5-8）。conditions を使うのは thermo と report の 2 か所だけで、宣言を忘れると既定値（298.15 K、1 atm、`system.py:95-98`）で黙って動く。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | △ | 順位の量が TST の障壁として定義されておらず注記もない（→ U7-P1）。ねじれの段が化学反応と同じ表に並ぶ（P3a で 1 位）。既定の順位は PBE0/SVPD。順位を止める判定が順位の量と無関係 |
| 汎用性 | △ | 開殻の WFT 参照がない。Z>36、二重混成で黙って誤りうる。SP 層に spin_ok がかからない |
| 有用性 | ○ | パネルは実計算で動き CSV も妥当。順位にも HTML にも届かない。composite は同じ run への追記でも成り立つが（id の衝突で view が上書き）、どの層で順位を付けたかが記録されない |
| 効率 | ○ | 1 点の費用は妥当。targets が outcome を見ず、非 rankable の点にも CCSD(T) を払う。逆に既定では単量体が抜ける |
| 簡潔さ | ○ | SP と停留点の対応付けが 3 通り、blocker の判定が 3 か所、使われない metric、CC 15 の rank_rows |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U8-I1 | 順位の量が定義されず、沈んだ障壁も数値のまま並び、BARRIERLESS は順位外 | validity | medium | → U7-I2。p4_malon の値は `thermo/manifest.json` で確認。序数の結論自体は化学的に妥当 |
| U8-I2 | 手法パネルが順位に反映されない。既定の順位は PBE0/SVPD。composite の来歴が残らない | usefulness | high | P1a/P1b。`thermochemistry.py:241, 248`、`manifest.py:80-81` |
| U8-I3 | 順位の量でない ΔE_rxn の符号不一致（±1 kcal/mol）で順位を止める | over_diagnostics | medium | `summary.py:31-32, 149-150, 214-215`、`gates.py:56-58`。縮退反応では常に 0 で無意味。SN2 の ΔE‡ 4.7 kcal/mol のずれでは何も起きない |
| U8-I4 | dzpe_out_of_tolerance の経験則 blocker | over_diagnostics | medium | → U7-I7 |
| U8-I5 | 順位の基準が各反応の 1 配座 | validity | medium | → U7-I1。SN2 は錯体が分離反応物より約 12 kcal/mol 低いので分離基準は値を変えない（分離基準が効くのは ΔG_assoc>0 の弱い錯体） |
| U8-I6 | 開殻種の WFT 参照がなく、パネルでは failed として数えられる | generality | medium | `engine.py:326-328`、`input.py:156-167`、P1c |
| U8-I7 | パネルの手法で黙って誤りうる（ECP なし、二重混成の PT2 抜き、SP 層のスピン汚染を検査しない） | validity | medium | `input.py:82`、`engine.py:58-61`、spin_ok は `minimum.py:93`・`actions.py:327` の freq だけ |
| U8-I8 | sp の targets が outcome を見ず、非 rankable の点にも高価な SP を払い、単量体は既定で抜ける | efficiency | medium | `single_point.py:39-43` |
| U8-I9 | SP と停留点の対応付けが 3 通り（parents、start の指紋、final の指紋）、blocker の判定が 3 か所 | complexity | low | `single_point.py:47-52`、`thermochemistry.py:66-81, 179-187`、`summary.py:139-152, 179-190` |
| U8-I10 | 追記用 pipeline が conditions と methods の位置を暗黙に仮定する | responsibility | low | `config.py:133-135`、`report.py:47-48`、P1c（rc=2） |
| U8-I11 | 使われない metric、HTML が各反応の先頭の ReactionThermo を表示（条件がずれうる）、パネルが HTML に出ない | complexity | low | `report.py:26`、`html.py`。RankRow の互換用既定値（`records.py:207-210`）は既存 run への追記に実際に必要 |
| U8-I13 | 順位の欠けた因子（κ、m、σ の脆さ）が系統誤差になる | validity | medium | → U7-I4・I5・I6、U3-P2 |
| U8-I14 | パネル・順位の経路の多くが未検証（符号 blocker、開殻 composite、分離基準、複数条件、barrierless） | untested_path | medium | → §6 |
| U8-I15 | パネルの交換の軸が 2 点（PBE0 25% と ωB97X-D3 の範囲分離）で、HF 交換の多い meta 混成がない | generality | low | 欠落ではなく補完の余地 |
| U8-I16 | mode-follow から来たねじれの段が化学反応と同じ順位に並ぶ | validity | medium | P3a の `rxn_mode_follow_be1e5e1114`（ΔG_rxn −0.0002、ΔG‡ 3.46）が 1 位 → U7-P1 の torsional 列 |
| 重複 | 上流の失敗で report が出ない（U8-I12）→ U9-I1 | | | |

**改良案**

#### U8-P2 符号不一致の blocker を削除し、パネルの幅を ranking に列で出す（must、【複雑さ減】）

U0-P11 を統合。

- 変更: `SIGN_DISAGREEMENT`、`_SIGN_DEADBAND_KCAL`、`PanelRow.sign_disagreement`、`method_panel` の戻り値 disagreements、`_RANK_BLOCKERS` の `method_sign_disagreement` を削除する。**`participant_notes` 関数そのものも削除する**（符号を除くと rankable に届くのは spin_contaminated だけで、それは `thermo._blockers`（`thermochemistry.py:185`）がすでに立てている。rank_rows の引数も減り CC が下がる）。method_panel.csv の min/max 列は残し、ranking.csv にはパネルがあるときだけ `dE_act_panel_min`・`dE_act_panel_max` の 2 列を足す。
- 理論的根拠: 順位は δG_eff の序数で、手法誤差は ΔE‡ の幅で表すのが直接的。ΔE_rxn の符号は順位の量ではない。
- 効果: 良い手法を足すほど順位が消える逆効果がなくなり、手法の不確かさが見える。
- コスト: 約 −25 行、列 2 つ。
- 検証: `tests/unit/reporting/test_summary.py` の符号テストを列の出力テストに差し替え、hcn_panel と sn2_panel に report を再実行し列の値がパネルの min/max と一致。

#### U8-P3 composite のエネルギー層を順位の正式な経路にし、条件と来歴を 1 か所にする（should、【複雑さ減】）

U0-P3、U9-P9（F12・methods の部分）を統合。

- 変更:
  - 切り替える前に、既存 run の複製（HCN、HONO、SN2、p4_malon、ch3o_doublet_b）で ωB97X-D3/def2-TZVPD と CCSD(T) の差が PBE0/SVPD より小さいこと、所要時間の増分（17 原子の TMA·(HF)₂ で計測、discover では極小が数十点）を実測する。そのうえで既定の discover・known_endpoints に `{id: sp, stage: sp, engine: nwchem, methods: [wb97x-d3_def2-tzvpd]}`（対象は U8-P4 の規則）と thermo の `energy_method` を入れる。構造と Hessian は PBE0-D3BJ/def2-SVPD のまま。
  - エネルギー層で ΔE‡ ≤ 0 になる低障壁の経路（TS は PBE0 の PES 上の停留点なので起こりうる）は、U7-P1 の規則（ΔE0‡ ≤ 0 なら障壁なし扱い、note submerged_barrier）で自然に扱い、新しいゲートは作らない。SCREEN の障壁なし判定は PBE0 の PES 上のものだと明記する。
  - conditions を唯一の消費者である `ThermoConfig`（temperatures_K、standard_states）に移し、pipeline 全体の conditions を削除する。report は ReactionThermo にある (T, state) の最初のもの（または ReportConfig の指定）を使う。method_panel.yaml は sp（パネルの手法）→ thermo（energy_method に参照手法、別の stage id）→ report にし、conditions の書き写しをなくす。
  - `ReactionThermo` に energy level のラベルを 1 つ持たせ（`SpeciesThermo.energy_calc` はすでにある）、ranking.csv に列として出す。同じ run への追記では「後から追記した層が勝つ」（`manifest.py:80-81`）を仕様として残し、どの層で並べたかを列で示す。
  - `load_method` は CLI と同じ規則で探す（pipeline の隣の `methods/` があれば優先、なければ `configs/methods`、`config.py:135`）。
  - `design.md:149` の「約 3 kcal/mol」を反応クラス別の文献値（PBE0 の MSE: H 移動 −4.2、重原子移動 −6.6、SN −1.9 kcal/mol）に書き換え、`design.md:150` の cp -a・--run-dir の手順をこの 1 本の手順に置き換える。
- 理論的根拠: 分業（停留点は解析 Hessian の大域混成、精度はエネルギー層）。範囲分離汎関数は長距離交換を正しく扱い、伸びた結合の TS の過安定化が小さい（Zhao & Truhlar, J. Phys. Chem. A 109, 2012 (2005), doi:10.1021/jp045141s、GMTKN55 doi:10.1039/C7CP04913G、Lin et al., J. Chem. Theory Comput. 9, 263 (2013), doi:10.1021/ct300715s）。
- 効果: 実計算で動いたパネルが順位に届き、既定の順位の系統誤差が小さくなる（SN2 で PBE0/SVPD −3.0 に対し ωB97X-D3 +1.7、HCN ではほぼ一致）。conditions の二重宣言と黙った失敗、configs 外の pipeline の失敗が消える。
- コスト: 設定の変更と report・config で数十行。計算費は停留点ごとに TZVPD の RSH SP 1 本（17 原子で数分と推定、freq 約 15 分より小さい。要実測）。
- 検証: HCN、HONO、NH3、SN2、p4_malon を既定 pipeline と CCSD(T) 入りパネルの追記で回し、ranking が composite の値になること、energy level 列、conditions なしで追記できること、configs 外のパネル（ch3o_panel）が動くこと。

#### U8-P4 sp の targets を「順位に使う点」の 1 規則にする（should、【複雑さ減】）

- 変更: `_targets` を、(1) rankable な outcome（elementary、degenerate、reassigned）の participants と BARRIERLESS の両端、(2) その反応物・生成物と同じ state_label の DFT 極小（U7-P1 の最小 G 用）、(3) その組成の単量体・断片の DFT 極小（分離基準用。U3-P1 までは多原子の単量体だけ）にする。same_basin、blocked、out_of_window、unresolved、outcome None は計算しない。`SinglePointConfig.targets` を削除する。
- 理論的根拠: 費用の大きい層は結論に使う点にだけ払い、基準に必要な点は同じ LOT でそろえる。
- 効果: 非 rankable な点への CCSD(T) がなくなり、energy_layer_missing が消える。
- 検証: tma_hf2 と SN2 で SP の件数が減り、ΔG_assoc と vs_separated が composite で出る。outcome ごとの対象集合の単体テスト。

#### U8-P6 SP と停留点の対応付けを parents の 1 規則にし、blocker の生成元を 1 つにする（should、【複雑さ減】）

U7-P7 を修正して統合。

- 変更: パネルと thermo は、どちらも sp artifact の parents（minimum_id か freq_calc）で対象を引き、参照レベルのエネルギーは各対象の freq Evidence のエネルギーを使う。`summary._stationary_points`・`_deltas` の指紋の突き合わせと、`thermo._subjects` の start の指紋の辞書を削除する。blocker は `thermo._blockers` を事実（thermo_unavailable、mixed_level_of_theory、spin_contaminated）の唯一の生成元とし、`rankable` は outcome と、thermo の記録そのものが欠けている場合だけを足す（`rankable` での thermo_unavailable の再判定を削除）。ReactionThermo を事実だけにして判定を gates に全面移管する改造まではしない（最小の変更）。
- 理論的根拠: 1 つの概念は 1 か所で定義する。sp は計算した点を parents に正確に記録している（`single_point.py:47-52, 74-76`）。
- 効果: 判定の食い違いがなくなり、rank_rows の複雑度が下がる。正味で行数が減る。
- 検証: v5 run と hcn_panel に report を再実行し、method_panel.csv と ranking.csv がバイト単位で一致（U8-P2 の列の変更を除く）。

#### U8-P7 開殻の WFT 参照を ROHF-CCSD(T)（TCE）にし、MP2 を削除する（should、プローブ先行、【中立】）

U0-P8 を統合。

- 変更: まず OH• と CH3O• で TCE の ROHF-CCSD(T)/def2-TZVPD の所要時間と出力の文字列（`CCSD(T) total energy / hartree`、`open shells`）をプローブで確かめる。そのうえで render_wft は mult>1 のときだけ `scf; rohf; nopen <m-1>; maxiter 100; end`、`tce; ccsd(t); freeze atomic; end`、`task tce energy` を書き、閉殻は速い ccsd モジュールのまま。`output.py` の `_WFT_ENERGY` に 1 本足し、`wft_closed_shell_only` を削除する。MP2（configs で未使用）は U0-P4 で削除。T1 診断や多参照性のゲートは足さない。
- 理論的根拠: ROHF 参照の CCSD(T) はラジカル熱化学の標準的な参照（W1/HEAT 系）。NWChem の ccsd モジュールは RHF 専用、TCE は ROHF/UHF に対応（https://nwchemgit.github.io/TCE.html）。UMP2 はスピン汚染を受ける。
- 効果: 開殻反応（引き抜き、ラジカル転位）にも参照値が付き、パネルが開殻で失敗として数えられない。
- コスト: 描画の分岐 1 つと正規表現 1 本（MP2 の削除で相殺）。TCE は ccsd より遅いので、重原子約 5 個までという上限は据え置く。
- 検証: ch3o_panel の再実行で CCSD(T) 行が埋まる。同じ出力を golden に加えて多重度の観測を検査する。

#### U8-P8 エネルギー層の SP にも spin_ok をかける（should、【中立】）

- 変更: composite とパネルの UKS の SP Evidence に spin_ok をかけ、既存の spin_contaminated の blocker に合流させる（許容幅は U1-P6 を採るまでは今の絶対値 0.1）。ROHF-CCSD(T) には ⟨S²⟩ がないので対象外。ECP は U0-P1、二重混成の拒否は U0-P4。
- 理論的根拠: 別ジョブの SP は原子推測から SCF を始め、別の解に落ちうる。
- 効果: パネルと composite の値が黙って別の電子状態になる経路を塞ぐ。
- 検証: ⟨S²⟩ が大きい合成 SP Evidence で blocker になる単体テスト。

#### U8-P11 HTML を順位の条件にそろえる（could、【中立】）

- 変更: `html.render` は RankRow の (T_K, standard_state) に合う ReactionThermo を表示し、表に dG_eff とパネルの ΔE‡ min/max を出す。RankRow の互換用既定値は残し（extra='forbid' の下で既存 run への追記に必要）、`design.md` の「互換層を持たない」の記述をこの事実に合わせる。
- 検証: 複数温度（[250, 298.15, 400] K × [1 atm, 1 M]）の report を追記して HTML と ranking.csv が一致。

#### U8-P13 パネルの交換の軸の補完として M06-2X を任意の手法にする（could、【中立】）

- 変更: `configs/methods/m06-2x-d3zero_def2-tzvpd.yaml`（grid xfine、SP 専用）を足すが、既定の method_panel.yaml には入れない。パネルを「基底（SVPD→TZVPD）・交換（PBE0、ωB97X-D3、任意で M06-2X）・参照（CCSD(T)）」の 3 軸として読むと文書に書く。
- 理論的根拠: 障壁の誤差は HF 交換の割合と長距離補正でほぼ決まる（GMTKN55）。M06-2X は grid に敏感（xfine が要る）。
- 検証: HCN と SN2 のパネルに追加して値が文献の範囲にあること。

他区分へ統合: 順位の量と dzpe → U7-P1、上流の失敗と report → U9-P1、ECP → U0-P1、二重混成 → U0-P4、検証 → §6。

**見送り**

- 手法間の ΔE‡ の幅が閾値を超えたら順位を止める blocker: 符号 blocker と同じく止める理由にならない。
- 信頼度スコアや Bayes 的な誤差モデル: 「rankable か blocker を持つか」の 2 値の方針（AR-27）と較正データの欠如。
- パネルの各レベルで順序が変わらないかを基準に同順位を判定（U8-P10）: 同順位の意味を定義し直すことになり、反応ごとに共通のレベル集合が変わって解釈が難しい。U8-P3 と U8-P2 の列で足りる。
- 組成ごとの順位表: 表の構造が増えるだけ。
- metric の knob を残す: ranking.csv に dG_rxn 列があり、δG_eff が生成物側の max を含む。
- UMP2・DFT+MP2 を開殻の参照に: スピン汚染を受ける。
- 二重混成（dftmp2）を今パネルに入れる: PT2 の段、別の読み取り、解析 Hessian なしへの対応が要る。需要が出てから SP 専用で。
- WFT の費用を予測して自動で飛ばす: 予測の診断と knob が増える。
- T1 診断・(T)/E_corr を多参照性のゲートに: 閾値の較正が弱く、ccsd モジュールは T1 を出さない。
- sp の複数の点を 1 つの NWChem ジョブに: 互いに別の構造で movecs の引き継ぎも効かず、JobStore の粒度が崩れる。
- 配座変化の反応を別表: torsional 列で足りる。
- パネルを HTML に完全な表で: CSV で足りる。
- RankRow の互換用既定値の削除（U8-P11 原案、U9-P3 の一部）: 既存 run への追記ができなくなる。
- thermo の artifact id に stage id を入れて層ごとに共存させる: view に同じ反応の複数の thermo が並び、report がどれを使うかの規則が要る。energy level 列で来歴を示す方が単純。

**維持**: summary の純関数性と「rankable か blocker を持つか」の 2 値、同順位の競争順位（1, 1, 3）、coverage.csv、パネルの鍵を Level.full_key にし method_panel.csv で min/max を並べること、sp を freq の最終構造で計算し parents に subject を記録すること、composite の式と層の欠けた化学種を含む ΔE を None にすること（C28）、温度と標準状態を 1 か所の条件から決め 1 atm から Δn で換算すること、ranking.csv の T_K・standard_state・dG_rxn_kcal 列（C25）、WFT を SP 専用にする分業、失敗を型付きの failed artifact として保存すること、静的な report.html。

---

### 3.9 U9 実行基盤・パイプライン・設定・アーキテクチャ

**現状**

- `pipeline/runner.py`: stage の唯一の入口 `execute_stage`（consumes の検査、run、produces の検査、manifest、状態記録）。`input_sha`・`config_sha` による再開、SiteLock。Runtime は engine キャッシュ、method、thread_map、deadline。consumes が欠けると ValueError（:150-157）。
- `execution/jobs.py`: JobStore の内容アドレス鍵（`Task.key_payload`）で再生。LADDER（TIMEOUT 2、GEOMETRY_MAXITER 2、INPUT_INVALID 1、SCF_NOT_CONVERGED 1）。`_CoreSemaphore`（:63-82）が ranks×threads を予約。`adapter.parse` は例外処理なし（:151）。
- `execution/process.py`: 新しいセッションでプロセスグループごとに起動し、timeout か Python の例外でツリーを kill（:103-113）。SIGTERM のハンドラはない。
- 契約: `Policy` 13 knob（`gates.py:28-43`）、`CasePolicy` 10 knob（`state.py:41-55`）を dict で受け手書きで検査（`reaction_paths.py:30-67`）。configs/ に `gates:` も `policy:` も 0 件。ExecutionSpec の maxiter・coordinates が鍵の外。
- 並列: DFT とケースは直列（`reaction_paths.py:116-125`）。thread_map は CREST と GoodVibes だけ。予算はケースごとに新しい 6 h（`driver.py:106`）、分割した子も新しい予算。

**評価**

| 観点 | 評価 | 理由 |
|---|---|---|
| 妥当性 | ○ | 内容アドレス型 JobStore、型付き Failure と上限付き LADDER、再開、唯一の入口、純関数の decide は健全。鍵の外の化学 knob（設定している site は 0 件で潜在的） |
| 汎用性 | △ | 新しい化学で出る想定外の失敗が 1 件あると run 全体が止まる。元素表の分散（worker の原子番号表で Te が KeyError） |
| 有用性 | ○ | 再生・再開・manifest・log.jsonl・preflight は実務で役立つ。1 件の例外で他の結果も出ない、SIGTERM の後に孤児プロセスと「running」のままの run_state、追記 pipeline の黙った失敗 |
| 効率 | △ | 主因は読まれない below_zpe の DFT freq と saddle の LADDER 継続。explore の直列の損失は約 1 分 |
| 簡潔さ | △ | 上書きされない knob 23、読まれないフィールド、定義の重複、試験用 override、GoodVibes の足場 |

**課題**

| ID | 課題 | 種別 | 重大度 | 根拠 |
|---|---|---|---|---|
| U9-I1 | 失敗を閉じ込める単位が不統一で、想定外の例外 1 件で run 全体が落ちる（job の parse、case の例外、上流 stage の全件失敗で report が出ない） | responsibility | high | `jobs.py:151`（P2a）、`reaction_paths.py:120-125`、`runner.py:150-157`（P6a・P6b）。(3) は下流に入力がないので止まること自体は理にかなうが、report が出ず traceback で終わる |
| U9-I2 | 上書きされない knob が 23。数値照合の許容幅と化学の閾値が混在し、化学の結果を変える設定が鍵の外 | complexity | medium | `gates.py:28-43, 67-69`（scf_noise_factor は既定で効かない）、`state.py:41-55`、`method.py:102-113`、`nwchem/engine.py:186, 218`、`input.py:108, 148`、`xtb.py:113` |
| U9-I3 | 読まれない記録フィールドがあり、その一部は高価な計算（below_zpe の DFT freq）を伴う | over_diagnostics | medium | `records.py:101-111`、`evidence.py:100-103`。diagnostics.json は v2 で人のデバッグ用に意図して残したもので負債は軽い |
| U9-I4 | 同じ定義・判定・プリミティブの重複 | complexity | medium | 接続済み outcome の集合 3 か所（`classification.py:17`、`gates.py:53`、`state.py:135`）、ConnectionLabel 2 か所、StandardState 3 か所、Kabsch・二面角の重複、JobStats と JobCounts、file_ref 2 か所、原子的書き込み 4 か所、`engines._OVERRIDES`（:51-105）、Evidence が保証する本数の gates での再検査 |
| U9-I5 | 予算の単位が閉じていない（子 case が新しい 6 h、walltime の行が証拠による完了より前） | efficiency | medium | `driver.py:106`、`state.py:217-221`、`design.md:144`（最悪 1 反応 7 case = 42 h） |
| U9-I6 | 並列の仕組みが二重（セマフォと threads_per_item）で、explore が直列 | efficiency | low | `runner.py:77-83`、`explore.py:73-80`。minima の screen tier は `known=registry` に依存するので単純には並列にできない |
| U9-I7 | 層の責務の漏れ（stage が run のレイアウトを推測、driver が .npy を読む） | responsibility | low | `reaction_paths.py:75, 90`、`actions.py:169, 393`。F12・methods の位置は U8-P3 |
| U9-I8 | SIGTERM やスケジューラの打ち切りで NWChem の MPI ランクが孤児になり、run_state が「running」のまま残る | correctness | medium | `process.py:103-113`、T（`timing_acac.log`）。子は `start_new_session` で別セッション |
| U9-I9 | LADDER の継続がジョブの種類を区別せず、saddle の maxiter でも壊れた Hessian で 2 回継続する | efficiency | medium | `engine.py:46-47`、P5c → U6-P6。SCF の上限 30 → U0-P5 |
| U9-I10 | 閉じた式（qRRHO）を外部エンジンとして動かす | complexity | medium | → U7-I8 |
| U9-I11 | 検証の空白（回帰セットと golden が中性・閉殻・第 2 周期に偏る） | untested_path | high | configs/systems に charge・multiplicity の指定 0 件、golden 125 本はすべて電荷 0・一重項。今回の critical・high の欠陥はここから見えなかった → §6 |
| U9-I12 | 使われない溶媒和と d4 の半実装 | complexity | low | → U0-P4 |

**改良案**

#### U9-P1 失敗の封じ込めを「最小の単位で型付きに記録し、run は続ける」の 1 規則にする（must、【中立】）

U8-P9 を統合。

- 変更: (a) `JobRunner._attempt` で `adapter.parse` を try で囲み、予期しない例外は `Failure(kind=INCOMPLETE_OUTPUT, reason='parse:<例外型>: <要約>')` にして JobStore に保存する（FailureKind は増やさず、LADDER の対象外）。JobStore の鍵には hfauto のコード版数が入らないので、パーサを直した後は既存の `--retry-failed incomplete_output` で取り直すことを docs に 1 行書く。(b) `reaction_paths` の while ループで `drive_case` を try で囲み、例外が出たケースは UNRESOLVED（`reason='error:<例外型>'`）の record にし、要約をそのケースの `log.jsonl` に 1 行書いて他のケースは続ける。(c) `execute_stage` は consumes が欠けたとき、status='done'・0 artifact とし、上流の失敗の要約を記録して次の stage へ進む（report は同じコマンドで必ず出る）。終了コードは「artifact に failed がある」とき非 0 にするだけで、新しい状態は足さない。テストで例外を再送出する厳格モードは pytest 用の環境変数 1 つに限る。値の個別検査は増やさない。
- 理論的根拠: バッチ計算基盤の障害の閉じ込め。仮説・化学種・ジョブは互いに独立な単位なので、失敗もその単位で閉じる。FailureKind による型付きの分類を parse 経路にも一貫して適用する。
- 効果: 単原子・Z>36・未知の出力形式などの想定外が 1 件あっても、他の化学種と反応の結果と report が出る。
- コスト: 約 15〜20 行とテスト 3 件。
- 検証: 単体で 3 経路。WSL で `sn2_cl`、`p6_h2te`、`p6_sn2_i` を修正前（U3-P1、U1-P1 の前）のコードで再実行し、run が完走して failed artifact の reason と report が出ること。v5 の 4 系の不変。

#### U9-P2 knob を「利用者が変える理由のある化学の閾値」だけにする（should、【複雑さ減】）

U0-P10 を統合。

- 変更: `Policy` を 13 から 5 にする（noise_cm1、saddle_cm1、resolution_kcal（barrier_proceed_kcal の改名）、reaction_window_kcal、spin_tol）。torsion_saddle_cm1・ts_prominence_hartree（U6-P5）、rank_dzpe_*（U7-P1）、thermo_*_tol（U7-P5）を削除し、scf_noise_factor・qrc_min_drop_hartree は gates のモジュール定数にする。`CasePolicy` は型付きの `ReactionPathsPolicy{walltime_h, max_saddle_attempts, max_split_depth}` にし、qrc_*・string_*・screen_images はモジュール定数に、`_POLICY_KEYS`・`_known_keys`・`case_policy` の手書き変換を削除する（screen の有無は ScreenConfig の有無で決める）。`ExecutionSpec.maxiter` を削除し、`coordinates` は autoz 失敗時の継続だけで使う `task.inputs['cartesian']` に移す（scf_rescue と同じ運び方）。
- 理論的根拠: 上書きされず検証もされない knob は保証のない仕様。キャッシュの鍵には結果を決めるすべての入力を入れ、結果を決めない実行設定（ranks、memory、timeout）だけを外に置く。
- 効果: knob 23 → 8。site の設定で同じ鍵の下に違う化学の結果が混ざらない。reaction_paths.py が約 25 行減。
- 検証: `tests/unit/pipeline/test_config_load.py` の更新、v5 の 4 系の ΔG‡ と outcome の一致（既定値は不変）、削除した knob の configs での参照 0 件。

#### U9-P3 読まれないフィールドを削除する（should、【複雑さ減】）

- 変更: 規則を 1 つ「判断・report・下流 stage のどれかが読むものだけを record に置く」と決め、`PathProfile` の gmax_history・program_converged・climbing_image とそのパーサ、`ConnectionClaim` の notes・method、`ThermoResult.S_rot` を削除する（BarrierVerdict は U5-P1、ReactionTrial.perturbed は U4-P9、population は U7-P1 で扱う）。minima の `diagnostics.json` と RankRow の既定値は残す。
- 効果: 型が判断の契約だけを表す。
- 検証: grep で読み手 0 件を確かめてから削除、HONO・HCN の report.html が全行を描画。

#### U9-P4 定義とプリミティブを一本化する（should、【複雑さ減】）

- 変更: (1) `core/records.py` に `CONNECTED_OUTCOMES`、`ConnectionLabel`、`StandardState` を 1 つずつ置き、gates・thermo・classification・state はそれを参照。(2) `chemistry/geometry.py` に kabsch、rotation_about、dihedral、neighbours を集め、read_xyz を trajectory 版の 1 フレーム版にし、zpe_hartree を thermo へ。(3) 元素表 → U1-P1。(4) JobCounts を削除して JobStats を使う。file_ref は JobStore の 1 つ。原子的書き込みを core に 1 つ。`engines._OVERRIDES` と `requirements_for` の特例を削除し、試験は `Runtime.create` に fake を渡す。(5) freq の不変条件（本数、虚振動の本数）は Evidence の model_validator で 1 回だけ保証し、gates の `_mode_count_reasons` を削除。site を program 単位にまとめる件（YAML アンカーの解消）は効果が小さいので後回し。
- 理論的根拠: 1 つの概念には 1 つの定義（K4・K8 の判定の食い違いはその実例）。
- 効果: 片方だけを直して判定が食い違うことが構造的に起きない。約 −100 行。
- 検証: ruff、pyrefly、lint-imports、全テスト、v5 の 4 系のビット一致。

#### U9-P5 予算を「1 仮説 = walltime_h（子 case を含む）」に閉じ、証拠で決まる完了を予算判定より先にする（should、【中立】）

U5-P8、U5-P10（予算の部分）、U6-P6(b)、U6-P8 を統合。

- 変更: `reaction_paths` の queue に (case, depth, deadline) を積み、分割した子は親の Deadline を引き継ぐ（`drive_case` は deadline を引数で受け、`CaseRuntime.deadline` の factory は根の case でだけ使う）。決定表で、手元の証拠だけで決まる**完了**の分岐（行 5 の接続完了、行 9、U5-P1 でまとめた barrierless の行）を行 4（walltime）より前に置く。新しい計算を始める分岐（行 5 の connection_retry、VALIDATE_INTERMEDIATE、行 10 の再試行）は行 4 の後に残す。`design.md:144` の「最悪 7 case = 42 h」を直す。
- 理論的根拠: 予算は新しい計算を始めるかどうかの判断で、得た証拠の分類には関係しない。予算は利用者が数えられる単位で閉じる。
- 効果: 最悪時間が「仮説数 × walltime_h」で読める（1 反応最大 42 h → 6 h）。締め切り直前に確定した接続・中間体を捨てない。knob は増えない。
- コスト: 約 10 行と決定表の順序の変更、テスト 2 件。
- 検証: fake engine で分割した子が親の残り時間で打ち切られること、expired=True でも接続完了が ELEMENTARY になること。trans-HONO→HNO2 で子 case の合計時間が walltime 以内。

#### U9-P7 SIGTERM・SIGHUP で外部プロセスを確実に止める（should、【中立】）

- 変更: `process.py` に実行中のプロセスグループの集合を置き（register/unregister で約 5 行）、`cli/main.py` の入口で SIGTERM・SIGHUP のハンドラを登録する（POSIX のみ）。ハンドラは登録済みのグループを `killpg` してから `sys.exit(128+signum)` を呼ぶ。これで thread_map のワーカー（CREST、U9-P6 の後の explore）が起動したプロセスも止まり、`execute_stage` の `except BaseException` で stage が failed として記録される。SIGKILL は捕まえられないので、スケジューラの KillWait の猶予内で終わることだけを文書に書く。
- 理論的根拠: POSIX のシグナルの既定の処理は Python の例外を起こさず、SystemExit は主スレッドにしか届かない。ThreadPoolExecutor の終了処理はワーカーを待つ。
- 効果: バッチの打ち切りで MPI ランクが孤児として残らず、run_state が正しく failed になる。
- 検証: WSL で `timeout 60 hfauto run` を acac で実行し、終了後の `pgrep nwchem/prterun` が 0 件、stage が failed。子プロセスを起動して自分に SIGTERM を送る単体テスト。

#### U9-P6 並列度の上限をコアのセマフォだけにし、explore の試行を並列にする（could、【複雑さ減】）

U4-P8 を統合。

- 変更: `Runtime.thread_map` と StageRuntime から `threads_per_item` を削除し、workers=site.cores（同時実行数は `_CoreSemaphore` が決める）。explore を (source, trial) 単位の thread_map にし、NT2→AFIR のフォールバックは単位の中に置く（結果は入力順で artifact は決定的、ReaDuct の OMP は threads に固定）。minima の screen tier と DFT tier、reaction-paths の case は直列のまま。HPC での ranks・string の mode parallel は `environment.md` に書くだけ。
- 効果: explore が約 1/4（v3 tma_hf2 で 81 s → 約 25 s）。U4-P1 で試行が一般化するほど効く。並列化の仕組みが 1 つ。
- 検証: tma_hf2_discover で explore の時間、JobStore の鍵の集合と artifact の一致、CREST（threads 4）が同時に 1 本だけ走ること。

#### U9-P9 ファイルの解決を Runtime に戻す（could、【複雑さ減】）

- 変更: StageRuntime に `resolve(FileRef) -> Path` を加え、`reaction_paths` の `stage_dir.parent` の推測をやめる。driver の `np.load` は U5-P3（`_tangent_mode_cm1` の削除）と U6-P4（振幅を ν と質量から）で消える。
- 検証: 既存テストと v5 の 4 系。

他区分へ統合: SCF 既定値 → U0-P5、saddle の継続 → U6-P6、GoodVibes の足場 → U7-P5、溶媒和・d4 → U0-P4、軌道の引き継ぎとエネルギー一致 → U3-P7、below_zpe → U5-P3、検証セット → §6。

**見送り**

- 反応ケース単位の並列化・DAG スケジューラ・非同期実行: Registry の順次更新による決定性が崩れ、4 コアでは NWChem 単体で CPU と壁時計がほぼ等しい。
- 異なる構造の SP を 1 入力にまとめる一般的な仕組み: JobStore の粒度と再生が崩れる（同じ経路のノードだけ U5-P9）。
- ExecutionSpec 全体を鍵に入れる: ranks・memory・timeout の変更でキャッシュが無効になる。
- run 全体を try で囲み全例外を握りつぶす: バグが見えなくなる。
- 列挙数・コスト予測・Hessian の種類の新しい診断ファイルや実行時検査: 読み手がいなければ負債。検証は §6 で行う。
- stage 全体の予算と SCREEN 値順の REFINE: 大きな系を回すまで効果が小さい。
- Policy・CasePolicy のすべての値を YAML から変えられるようにする: 検証していない値を開くことになる。
- 数値 Hessian の変位を hfauto で並列に: 停留点は解析 Hessian の手法に限る方針なので対象がない。
- movecs を Evidence の FileRef として一般化し後続ジョブに引き継ぐ（U9-P13 原案）: SCF 節約は 1 本約 20 s で freq の 900 s に比べ小さく、鍵の外で結果を左右しうる入力を渡すことになる。エネルギー一致の条件（U3-P7）で不一致が出たときだけ既存の restart で追加する。
- minima の screen tier の並列化: `known=registry` に依存し、並列にすると決定性か freq 本数が変わる。xTB の緩和は秒単位。
- diagnostics.json 2 本の削除（U9-P3 原案）: minima のものは mode-follow 経過の唯一の記録（conformers のものは U2-P4 で削除）。

**維持**: `execute_stage` を唯一の入口にし consumes・produces の型を検査すること、`input_sha`・`config_sha` による再開と SiteLock、内容アドレス型 JobStore と restart の入力を鍵に入れない方針、型付き Failure と上限付き LADDER、`_CoreSemaphore` によるコア予約、入力順を保つ thread_map、新しいセッションとプロセスツリーの kill と result.json のサイドカー、site を絶対パスを持つ唯一の層にすることと preflight、純関数の decide と Evidence・gates の分離、ケースごとの `log.jsonl`、DFT tier と反応ケースの直列、Level による電荷・多重度の照合、import-linter による層の規則。


---

## 4. 横断的な改良

### 4.1 統合の対応表

| 担当 ID | 統合した ID（原案） |
|---|---|
| U0-P1（ECP） | U1-P1 の ECP 部分、U3-P10 |
| U0-P4（原則と削除） | U9-P12、U8-P8 の二重混成 |
| U0-P5（SCF 既定値） | U9-P10 の SCF 部分 |
| U1-P1（元素表） | U0-P2、U2-P3、U9-P4(3) |
| U1-P3（スピン状態の宣言） | U0-P6、U2-P1、v4 の C8 |
| U1-P4（結合判定） | U2-P7、U4-P3 |
| U1-P9（対象範囲の文書） | U7-P9 |
| U3-P1（単原子） | U0-P7、U1-P2、U7-P3 |
| U3-P2（鏡映を含む同一性と m） | U7-P4 の m、U2-P5（修正） |
| U3-P5（崩壊鞍点の freq 再利用） | U6-P9 |
| U3-P7（freq と親のエネルギー一致） | U9-P13（修正） |
| U4-P7（生成物の同一性判定） | U3-P6(2) |
| U5-P1（経路の分類器） | v4 の C16 |
| U5-P3（below_zpe 削除） | U9-P3 の該当部分 |
| U5-P6（縮退判定） | U3-P2 の縮退判定部分 |
| U6-P2（反応方向） | U5-P5 |
| U6-P6（決定表と停滞鞍点の再開） | U9-P10 の saddle 部分 |
| U7-P1（順位の量と dzpe 削除） | U7-P2、U8-P1、U8-P5 |
| U7-P5（GoodVibes の同一プロセス化） | U9-P8（修正） |
| U8-P2（符号 blocker 削除） | U0-P11 |
| U8-P3（composite を既定・条件と来歴） | U0-P3、U9-P9 の F12・methods 部分 |
| U8-P6（対応付けと blocker） | U7-P7（修正） |
| U8-P7（ROHF-CCSD(T)） | U0-P8 |
| U9-P1（失敗の閉じ込め） | U8-P9 |
| U9-P2（knob の整理） | U0-P10 |
| U9-P5（予算と行の順序） | U5-P8、U5-P10 の予算部分、U6-P6(b)、U6-P8 |
| U9-P6（並列） | U4-P8 |
| §6（検証セット） | U0-P12、U1-P8、U3-P9、U4-P10、U6-P10、U7-P8、U8-P12、U9-P11 |

### 4.2 共通方針（汎用化・効率化・簡素化）

1. **入口で 1 回だけ宣言し、下流は宣言を信じる。** 元素（U1-P1）、スピン状態（U1-P3）、対象範囲（U1-P9）は structures と文書で 1 回だけ決める。途中の stage に同じ検査を散らさない。
2. **失敗は最小の単位で型付きに閉じ込める**（ジョブ → ケース → stage、U9-P1）。run は続け、report は必ず出す。
3. **化学の規則は元素に依らない最小の規則にする。** 結合は Σr_cov + 0.4 Å（U1-P4）、探索は結合変化の 4 テンプレート（U4-P1）、経路は山と井戸の分類器（U5-P1）、反応方向は反応中心（U6-P2）、変位はエネルギー目標（U6-P4）、初期 Hessian は慣性を整える 1 関数（U6-P1）。元素名の列挙表を持たない。
4. **判定は定義そのものだけにし、場所は 1 か所。** TS は負の固有値 1 本、接続は両側の降下と割付け（U6-P5）、同一性は identity の 1 組の基準（U3-P2）、順位の可否は `thermo._blockers` と `rankable`（U8-P6）。順位の量でない量で順位を止めない（U8-P2、U7-P1）。
5. **費用は結論に使う点にだけ払う。** 読まれない診断のための計算をしない（U5-P3）、エネルギー層は rankable な点と基準の点だけ（U8-P4）、窓は最良の推定で切る（U3-P6）、同じ drive・同じ生成物は 1 回（U4-P2、U4-P7）。
6. **キャッシュの鍵は「結果を決めるもの全部、それ以外は外」**（U9-P2）。
7. **検証は実行時の診断ではなく、分単位の回帰セットと golden で行う**（§6）。
8. **スキーマ変更の検証は `--from` による再実行で行う。** records は extra='forbid' なので、フィールドを削除すると既存 run の manifest が読めなくなる。複製した run を `--from reaction-paths` などで再実行すれば JobStore に当たって安い。新しく足す列には既定値を付ける。

### 4.3 削除できるもの一覧

| 対象 | 場所 | 理由 | 担当 | 概算 |
|---|---|---|---|---|
| 元素表 4 つ（`_ELEMENTS`、半径表、`ISOTOPIC_MASSES`、worker の `_ATOMIC_NUMBERS`） | `electronic_state.py`、`topology.py`、`vibrations.py`、`goodvibes/worker.py` | 1 つの表に置き換え | U1-P1 | −40 行（+データ 90 行） |
| 組成の `[-1]` と冗長な電子状態の再検査 | `conformer_search.py:207-211` | 宣言に置き換え、検査は構造上必ず通る | U1-P3 | −5 行 |
| `bonds` の previous とヒステリシス | `topology.py:79-89, 117-124` | 集合の差で対称に | U1-P4 | −10 行 |
| `MethodSpec.solvation`、`d4`、`mp2` とアダプタの分岐 | `core/method.py`、各 backend | 半実装・未使用 | U0-P4 | −50 行 |
| `ExecutionSpec.maxiter`・`coordinates` | `core/method.py:102-113` | 鍵の外の化学 knob | U9-P2 | −10 行 |
| Policy の 8 knob、CasePolicy の dict と手書き検査 | `gates.py:28-43`、`reaction_paths.py:30-67` | 上書き 0 件・検証なし | U9-P2 | −60 行 |
| conformers の `_same`・定数 2 つ・`diagnostics.json`・`topology_removed`・`_REMOVALS`・`"noopt"` 鍵 | `conformer_search.py`、`crest.py` | 実質働かない・読まれない | U2-P4 | −30 行 |
| 円錐配置（条件付き） | `placement.py:92-163` | CREST がサンプリングする | U2-P6 | −110 行 |
| `same_minimum`、`MinimumPolicy`、`amplitude_A`、`_soften` の固定振幅 | `identity.py:120-132`、`minimum.py` | 二重実装・未使用 knob | U3-P2、U3-P3、U3-P8 | −35 行 |
| 4 段の経験則、芳香環検出、酸塩基の表、`Window`、`product_verdict`、`nt_total_force_norm`、`perturbed` | `trials.py`、`explore.py`、`records.py` | 列挙器 1 つと gates の窓へ | U4-P1、U4-P6、U4-P9 | −120 行（正味） |
| `_max_rel_kcal`、`shape`、低レベル OR 判定、BarrierVerdict の 4 フィールド | `gates.py:168-212`、`profile.py`、`records.py:101-111` | 分類器 1 つ | U5-P1 | −30 行 |
| `_tangent_mode_cm1` と below_zpe（DFT freq 1 本） | `actions.py:125-134, 161-170`、`gates.py:201-211` | 読まれない | U5-P3 | −25 行 |
| GS の二分岐、xTB 端点の再最適化、`_COLLAPSE_A`、pysis の座標系分岐 5 つ、GS 再実行（条件付き） | `actions.py`、`backends/pysis/` | 固定端 CI-NEB 1 本 | U5-P2 | −90 行 |
| `energies_settled`、`_SETTLED`、`energy_history`、`string_path_energies` | `actions.py`、`profile.py`、`output.py` | 1 チャンクで分類 | U5-P4 | −35 行 |
| 自動の配座対仮説 | `hypotheses.py:170-176, 192-195` | Curtin–Hammett | U5-P7 | −30 行 |
| `_moddir`、`cartesian_mode_number`、moddir=0 分岐、継続時の moddir 抑制、`mode_overlap_below_0.3` | `engine.py:87-95, 203-204`、`vibrations.py:104-113`、`input.py:140` | 反応モードだけ負の初期 Hessian | U6-P1 | −30 行 |
| low_prominence、torsional 閾値、no_initial_descent、trajectory_above_ts、`torsion_saddle_cm1`、`ts_prominence_hartree` | `gates.py:143-165, 215-228` | TS・接続の定義外 | U6-P5 | −15 行 |
| 行 10、`_RETRY_A`、`imaginary_modes[1]` 固定 | `state.py:174-177`、`actions.py:50, 333-337` | 状態の消費と一般の規則 | U6-P6、U6-P3 | −10 行 |
| dzpe の連鎖（blocker、knob 2、`n_h_transferred`、`transferred_hydrogens`、`proton_coordinate`）、population、0 強制の分岐、`metric` | `gates.py`、`records.py`、`topology.py`、`thermochemistry.py`、`report.py` | 経験則・読まれない・効かない | U7-P1 | −80 行 |
| GoodVibes の worker・JobStore・`ThermoBatch`・`ThermoEngine`・`Capability.THERMO`・`thermo_consistent`・knob 2・S_rot 検査・`ThermoSettings.symmetry` | `backends/goodvibes/`、`gates.py:276-294` | 閉じた式の外部エンジン化 | U7-P5 | −200 行 |
| 符号不一致の blocker 一式と `participant_notes` | `summary.py`、`gates.py:56-58` | 順位の量でない | U8-P2 | −25 行 |
| 指紋による対応付け 2 つ、`rankable` での thermo_unavailable の再判定 | `summary.py:179-190`、`thermochemistry.py:66-81`、`gates.py:305-310` | parents の 1 規則 | U8-P6 | −30 行 |
| `SinglePointConfig.targets` | `single_point.py` | 1 規則 | U8-P4 | −5 行 |
| pipeline 全体の conditions（ThermoConfig へ移す） | `pipeline/config.py` | 消費者は thermo と report だけ | U8-P3 | 移動 |
| PathProfile の gmax_history・program_converged・climbing_image、ConnectionClaim の notes・method、`ThermoResult.S_rot` | `evidence.py`、`records.py`、`protocols.py` | 読まれない | U9-P3 | −30 行 |
| 定義の重複（outcome 集合、ConnectionLabel、StandardState、geometry、JobCounts、file_ref、原子的書き込み、`_OVERRIDES`、`_mode_count_reasons`） | 各所 | 1 か所に | U9-P4 | −100 行 |
| `threads_per_item` | `runner.py:77-83` と呼び出し 3 か所 | セマフォだけで決まる | U9-P6 | −10 行 |

合計で約 −1,000 行（追加約 +300 行）、knob は 23 → 8 の見込みです（採否が条件付きの U2-P6・U5-P2 を含む）。

**簡素化の候補（could）**: 反応はすべて Δn=0 なので、ReactionThermo の `dG_act`・`dG_rxn` は `standard_states` ごとに同じ値の行が出ます（`thermochemistry.py:249-253`）。標準状態ごとの行は会合量だけに持たせれば足ります。

### 4.4 見送った主な案（横断）

各節の「見送り」に全件を理由付きで載せました。基盤の方向を決める主なものは次のとおりです。

| 案 | 見送る理由 |
|---|---|
| 停留点レベルを RSH・meta-GGA・二重混成にする | NWChem で解析 Hessian が無効になり freq が 4〜6 倍遅く、数値ノイズも乗る。高精度はエネルギー層で |
| CCSD(T)（DLPNO を含む）を既定のエネルギー層 | N^7 で重原子約 5 個まで。DLPNO は NWChem にない |
| 溶媒和を残して拡張 | 溶媒モデル・1 M・数値 Hessian まで端から端の設計が要る。気相専用と明記 |
| 開殻一重項（BS-UKS の最適化とスピン射影） | 複雑さに見合わない。入口で対象外、TS の検出だけ could |
| 全多重度の自動計算・xTB での多重度選択 | 計算が数倍、GFN はスピン非依存 |
| WBO による結合グラフ、状態ラベルに立体 | 計算器依存、TS 近傍で不安定 |
| Chemoton 全体・多段探索・マスター方程式・VTST/SCT | 目的に対して過大 |
| 分離反応物を順位の基準に混ぜる、κ を順位に入れる | 標準状態で基準が入れ替わる／T_c 以下でモデル誤差が補正と同程度 |
| qRRHO を自前で書き直す | 検証対象が増える。GoodVibes を同一プロセスで使う |
| 反応ケース・DFT tier の並列化 | Registry の決定性が崩れ、4 コアの NWChem では得がない |
| 実行時の新しい診断（コスト予測、Hessian 種別、列挙数、多参照性、HOMO>0） | 読み手がいなければ負債。検証は §6 で |
| movecs をジョブ間でキーの外で引き継ぐ | 結果を左右しうる入力を鍵の外に置くことになる。エネルギー一致の条件で検出する |

---

## 5. 計算時間の分析と効率化計画

### 5.1 どこに時間がかかっているか（実測）

| 項目 | 実測 | 占める割合・備考 |
|---|---|---|
| DFT の解析 freq | 9 原子 150 s、13 原子 267 s、15 原子 679〜703 s、17 原子 898 s（∝ N·nbf^2.8） | 極小と TS の freq は熱化学に必須。acac では freq 2 件でジョブ時間の 57%（1 件は捨てた xTB Hessian の代わり） |
| DFT opt | 15 原子で 1 歩 42 s、17 原子で 1 構造 377 s | rerank 後の窓なしで余分な構造にも払う |
| 鞍点探索 | 15 原子で 1 歩 33 s。acac は誤ったモードで 36 歩 1185 s、DME の停滞は 1 件 495 s（うち継続 330 s） | 方向と初期 Hessian で決まる |
| QRC | paths の 51%（SN2）、59%（DME）、310/651 s（マロンアルデヒド） | 大半は平坦な錯体の尾の収束。負の固有値による歩幅制限は SN2 で各 3 歩だけ |
| FIND_PATH（ZTS） | HONO 175 s（全体 332 s の 53%、落ち着いた後の 13 反復を含む）、H2O2 100 s、17 原子で 1 チャンク 86〜108 分 | GS の DLC 失敗・種のない IDPP から入る |
| SCREEN の SP | 15 原子で 12 件 430 s（ジョブ時間の 13〜18%） | 両端ノードの 2 本は捨てられている |
| below_zpe の DFT freq | 17 原子で約 900 s、NH3 で 5 s | 誰も読まない |
| CREST | run の 1% 未満（6〜12 s） | 問題なし |
| ReaDuct explore | 31 attempt で 81 s（run の約 4%）、1 コアで直列 | 下流の重複 DFT が真の費用 |
| xTB 全体 | run 全体で 1〜16 s | 無視できる |
| 熱化学 | 0〜1 s | 仕組みが過剰 |
| 最悪の予算 | 1 反応で 7 case × 6 h = 42 h | 子 case が新しい予算を得る |

### 5.2 何をすればどれだけ減るか

| 施策 | 対象 | 減る量（見込み） | 担当 | 確認方法 |
|---|---|---|---|---|
| below_zpe の DFT freq を削除 | SCREEN ごと | 17 原子で約 900 s／SCREEN | U5-P3 | NH3 v5 の再計算で freq ジョブが消え ΔG‡ 不変 |
| 反応中心の方向で xTB Hessian を採用 | 種の初期 Hessian | DFT Hessian 1 本（15 原子 703 s、17 原子約 15 分）と誤追跡（acac の 60 分） | U6-P2、U6-P1 | acac で DFT Hessian 0 本、PT の TS に収束 |
| 停滞した鞍点を継続せず新しい Hessian で 1 回再開 | saddle の maxiter | 1 件約 330 s | U6-P6 | ho_harness で約 165 s で返る |
| higher_order 再試行で TS freq を再利用 | 再試行 | 1 回あたり DFT Hessian 1 本（17 原子で約 15 分） | U6-P3 | DME C2v |
| 固定端 CI-NEB で GS 失敗と種なし IDPP を回避 | SCREEN → FIND_PATH | HONO で FIND_PATH 175 s（全体の約 53%）、H2O2 で ZTS 100 s | U5-P2 | HONO・H2O2 で FIND_PATH 回避 |
| FIND_PATH を 1 チャンクで分類 | string | 良い種なら後半のチャンク（17 原子で 1 時間以上） | U5-P4 | HONO・H2O2・acac の反復数 |
| rerank 後にも窓 | dft tier | 除いた 1 構造あたり 17 原子で約 21 分 | U3-P6 | tma_hf2_discover の DFT opt 本数 |
| drive の同値除去と AFIR の限定 | explore | attempt 数ほぼ半分（31 中 21 件が同じ移動、AFIR 約 50 件の大半が無効） | U4-P2、U4-P11 | amine_pilot2 の −1265i 3→1 件 |
| explore 生成物の重複除去 | dft tier | 重複 1 件あたり DFT opt 約 6 分（17 原子） | U4-P7 | TMA·(HF)₂ の opt 件数 |
| 自動の配座対仮説の廃止 | reaction-paths | 配座対 1 件あたり最大 6 h の case | U5-P7 | 仮説数 |
| sp の targets を rankable な点に | sp | 非 rankable な点の CCSD(T)（大きな系で 1 点数時間） | U8-P4 | SP 件数 |
| 子 case が親の予算を引き継ぐ | 予算 | 最悪 42 h → 6 h／仮説 | U9-P5 | fake engine の単体テスト |
| 崩壊鞍点の freq 再利用 | 中間体 | freq 1 本（17 原子で約 15 分） | U3-P5 | trans-HONO→HNO2 |
| 経路 SP を 1 ジョブに | SCREEN | 小分子で起動 11→1 回（HCN で約 15 s）、中分子で SCF 約 4 割減 | U5-P9 | HCN・CH3O• |
| `driver:linopt 0` | 初期 Hessian つき opt | opt の SCF 約 2 割 | U3-P11 | A/B |
| explore の並列 | explore | 4 コアで約 1/4（81 s → 約 25 s） | U9-P6 | tma_hf2_discover |
| def2-SVP の method（大きな中性系） | freq・opt | freq 約 3〜4 倍速（要実測） | U0-P13 | SVP と SVPD の比較 |
| SIGTERM で孤児を残さない | 実行基盤 | 次の run のコア占有を防ぐ | U9-P7 | pgrep 0 件 |

**増える計算**:
- composite のエネルギー層（U8-P3）で、停留点ごとに ωB97X-D3/def2-TZVPD の SP が 1 本増えます。17 原子で数分／点と推定し、freq（約 15 分）より小さい見込みです（要実測）。
- 探索の一般化（U4-P1）と縮退の発見（U4-P5）で生成物と仮説が増えます。仮説 1 件の費用は、上の施策で上限が決まります（U9-P5）。

### 5.3 規模の目安（4 コア、PBE0/def2-SVPD）

- 既知端点の小分子（≤ 10 原子）: 1 反応 1〜16 分（HCN 1:22、NH3 1:25、HONO 5:32、SN2 10.6 分、マロンアルデヒド 16 分）。
- 中規模（15〜20 原子）: 1 反応 1〜3 時間。freq が支配的で、方向と初期 Hessian の誤りが最大の損失。
- 約 35 原子で freq が timeout（4 h）に達する。これを超える系は def2-SVP（中性）か HPC の ranks 増で扱う。

---

## 6. 汎用性を担保する検証セット

新しい実行時の診断は足さず、分単位の系と golden 出力と QM なしの fixture だけで、元素・電荷・スピン・反応型と、決定表の各行・各エンジンの分岐を 1 回以上通します。到達表は既存の `log.jsonl` の decide の reason と job.json から作り、`docs/validation.md` に 1 回記録します。SN2 などの文献値は PBE0/SVPD の誤差を含むので参考値とし、合否は rc、電子数、Level の照合、⟨S²⟩、outcome、既知の数値（v5 の値）で判定します。

### 6.1 回帰と実計算の系

| # | 系 | 元素・電荷・スピン | 反応型・経路 | 通す分岐 | 期待値（合否） | 時間 |
|---|---|---|---|---|---|---|
| R1〜R4 | HCN→HNC、HONO、NH3 反転、水 | C・N・O・H、中性、一重項 | 異性化、縮退、同一 basin | 既存の全経路 | v5 の値（42.1783、12.2265、3.8395、same_basin）。HONO は U3-P2 の後 RT ln2 だけ下がる（約 11.82） | 各 0:11〜5:32 |
| S1 | Cl⁻ + CH3Cl（Cl⁻ を単量体に含む） | Cl、−1、一重項、単原子 | 恒等 SN2、分離基準 | 単原子、縮退、vs_separated、パネル | rc=0、degenerate、中心 ΔE‡ 10.45（PBE0/SVPD）、`dG_act_vs_separated` が出る（参考 −0.74） | 約 11 分 |
| S2 | I⁻ + CH3I | I（Z>36）、−1 | 恒等 SN2 | ECP、SCF | rc=0、α 電子数が ECP と整合、Level 一致（中心障壁の文献約 8 kcal/mol は参考） | 数分 |
| S3 | H2Te | Te | 極小のみ | 元素表 | structures→dft→thermo が通る | 数分 |
| S4 | CH3O• → CH2OH•（生成物を Cs で入力） | 二重項 | 1,2-H 移動 | 鞍点端点、mode-follow、鏡像 | 宣言反応が評価され（elementary）、ΔG‡ 約 30.6、鏡像化の偽反応が消える | 約 10 分 |
| S5 | CH3• + O2 の組成（m=2 を宣言） | 二重項（²+³） | 会合 | スピン宣言、CREST `--uhf 1` | 宣言なしは INPUT_INVALID、m=2 で UKS mult 2 | 数分 |
| S6 | OH···CH4 錯体（explore まで） | 二重項 | H 引き抜き | T1、restricted_open_shell | 引き抜きの生成物が出る | 数分 |
| S7 | NH3···ICl | I、中性 | ハロゲン結合 | 結合判定、T1/T3 | 2 断片、無意味な polar_h が出ない | 数分 |
| S8 | SO2·NMe3（または BF3·NMe3） | S、中性 | 配位付加体 | CREST `--notopo` 全原子 | 付加体の状態を CREGEN が捨てない | 数分 |
| S9 | FeCl3·CH4（structures・conformers まで） | Fe（d ブロック、宣言必須） | 剛体配置 | 元素表、placement | 例外なし、多重度未宣言は INPUT_INVALID | 数分 |
| S10 | amine_pilot2（NH3·HF） | 中性、一重項 | 二重 H 交換（縮退） | drive 同値除去、縮退の発見 | −1265i の drive が 3→1 件、product になり degenerate case が走る | 約 30 分 |
| S11 | DME、C2v の重なり形を種（`ho_harness`） | 中性 | メチル回転 | 負モード 2 本、higher_order、停滞の再開 | 一次鞍点 約 226i、ΔE‡ 約 2.3、停滞が再現しない | 約 10 分 |
| S12 | DME、粗い入力（`p5_dme`） | 中性 | 縮退 | 縮退判定 | degenerate_rearrangement（ΔE‡ 約 2.3） | 約 10 分 |
| S13 | H2O2 gauche(+)→gauche(−) | 中性 | ねじれ、低障壁 | ねじれ QRC、IDPP/NEB | 接続が確定、ΔE‡ 約 1.15（今の識別では elementary、U3-P2 の後 degenerate） | 約 4 分 |
| S14 | マロンアルデヒド | 中性 | PT、ZPE で沈む障壁 | δG_eff の submerged 規則 | note `submerged_barrier`、δG_eff は max 規則 | 約 16 分 |
| S15 | trans-HONO → HNO2（cis を中間体に持つ 2 段） | 中性 | 多段 | intermediate、分割、子 case、予算 | multi_step と子反応 2 本、合計が walltime 以内 | 数十分 |
| S16 | アセチルアセトンの PT（15 原子） | 中性 | 縮退 PT、傍観者の回転 | 反応方向、初期 Hessian | xTB Hessian 採用、DFT Hessian 0 本、PT の TS | 約 30〜60 分 |
| S17 | (HF)2 の供与体–受容体入れ替え | 中性 | 解像度未満の障壁 | barrierless（行 12・14） | barrierless、δG_eff = max(ΔG_rxn, 0) | 数分 |
| S18 | OH + H2（または H + H2） | 二重項 | 二分子の引き抜き | 固定端 NEB（U5-P2）、SCF の救済 | unresolved でなく TS が得られる（共線 H3 の真の障壁約 9.6 は参考） | 数分 |
| S19 | TMA·(HF)₂ に conditions [250, 298.15, 400] K × [1 atm, 1 M] を追記 | 中性 | 会合 | 複数温度、1 M | 反応の δG_eff は標準状態で不変、ΔG_assoc だけが動く | 数分（JobStore 再利用） |
| S20 | HCN・SN2・CH3O• のパネル | 中性・−1・二重項 | — | パネル、composite、ROHF-CCSD(T) | 列の値が min/max と一致、ranking が composite、CCSD(T) 行が埋まる | 約 15 分 |

### 6.2 golden 出力（数秒の実出力をテストデータに）

- OH•（二重項）の NWChem UKS opt・freq: Level の mult 2、⟨S²⟩ 約 0.753、spin_ok、freq と opt のエネルギー一致（U3-P7）。
- FHF⁻ の opt・freq・SP: 電荷 −1 の正規表現。
- HI の SP（ECP あり）: α 13、E≈−298.31 Eh、Level の basis 観測が一致。CCSD(T) の `freeze atomic` の凍結軌道数。
- Cl⁻ の xTB と NWChem の opt・freq: 振動数 ()、n_external 3。
- OH•（または CH3O•）の ROHF-CCSD(T)（TCE）: `CCSD(T) total energy` と `open shells`（U8-P7 の後）。
- ωB97X-D3 の freq（数値 Hessian）の既存出力は、原則の文書化の例として残す。

### 6.3 QM なしの fixture

- 単原子 Ar の S°(298.15 K, 1 bar) = 154.846 J/mol/K（JANAF）。
- δG_eff の場合分け（通常、ΔG‡<0、ΔE0‡_rev ≤ 0、生成物が TS より高い、障壁なし、最小 G の基準）。
- 経路の分類器のモデル曲線（[0,−1.5,−3,−1,1,3,5]、[0,2,0.5,3,6,9,10]、単調、3 点）。
- `shape_hessian` の慣性・床・反応モードの選び方。
- 結合判定の距離（I···N 2.8、I···I 3.5、FHF⁻ 1.14、I3⁻ 2.92）と `bond_changes` の対称性。
- 多重度の規則（HF+HF、OH•+OH•、²X+¹Y、[O][O]、[CH3]、[Fe+2]）。
- 鏡映を含む同一性（NH3 反転は degenerate、CHFClBr は chiral、CH4 は achiral）。
- Eckart の極限（could を採る場合）。

合計の実計算時間は約 3〜4 時間（S16 が最大）。must の受け入れに要る部分（R1〜R4、S1〜S4、S10、S12〜S14、S17）は約 2 時間です。

---

## 7. 改良ロードマップ

### 7.1 must（結論を誤らせる・主要な化学で止まる・主要機能が働かない）

| 段 | 内容 | 前提 | 段の検証 |
|---|---|---|---|
| M0 検証の土台 | U9-P1（失敗の閉じ込め）、§6 の golden と configs/systems の追加（まず R1〜R4 の基準値） | なし | 修正前のコードで `sn2_cl`・`p6_h2te`・`p6_sn2_i` が完走し failed と report が出る |
| M1 入口の汎用化 | U1-P1（元素表）→ U0-P1（ECP）、U3-P1（単原子）、U1-P3（スピン状態の宣言） | U0-P1 の Z は U1-P1 の表から | golden（OH•、FHF⁻、HI、Cl⁻）、S1〜S3、S5、S9、R1〜R4 の不変 |
| M2 同一性と経路の形 | U5-P6（縮退を basin 構造と掌性で）、U3-P2（鏡映を含む同一性と m）、U5-P1（分類器。SCREEN の barrierless は両端が DFT 極小のときだけ閉じる）、U5-P3（below_zpe 削除） | U3-P2 と U5-P6 は同時（そうしないと宣言された鏡像化が捨てられる）。U7 の m 反映も同時 | S4、S12、S13、S15、S17、R1〜R4（HONO は RT ln2） |
| M3 鞍点と接続 | U6-P2（反応方向）→ U6-P1（初期 Hessian の慣性）、U6-P5（TS・QRC ゲートを定義だけに） | U6-P2 は U1-P4 があれば望ましい（なくても可） | S13、S16、R1〜R4 の歩数と ΔG‡ |
| M4 探索 | U4-P2（drive の同値除去。`wl_classes` の公開を含む）、U4-P1（列挙器。元素表に最大配位数の列を足す）、U4-P5（縮退の発見） | U1-P1、U5-P6 | S6、S7、S10、S1 の explore |
| M5 順位 | U7-P1（δG_eff と dzpe の連鎖の削除）、U8-P2（符号 blocker の削除とパネル列） | M2（障壁なしの行が実際に発火すること） | S14、S17、ch3o_doublet の torsional 列、R1〜R4（`--from` 再実行） |

### 7.2 should（有用性・効率・簡潔さの大きな改善）

| 群 | 内容 | 前提 | 検証 |
|---|---|---|---|
| S-A 精度と報告 | U8-P3（composite を既定、条件を ThermoConfig へ、energy level 列、methods の探索）、U8-P4（sp の targets）、U8-P6（対応付けと blocker）、U8-P8（SP 層の spin_ok）、U8-P7（ROHF-CCSD(T)、プローブ先行）、U0-P5（SCF 既定値、プローブで救済策を決める） | U8-P3 は M5 の後、切り替え前に差と時間を実測 | S20、S19、S18 |
| S-B 経路・鞍点の効率 | U5-P2（固定端 CI-NEB、プローブで採否）、U5-P4（1 チャンク）、U6-P6（決定表と停滞の再開）、U6-P3（TS freq の再利用）、U6-P7（縮退の結合グラフ確認）、U3-P3（mode-follow の振幅と初期 Hessian）、U3-P4（鞍点端点の登録）、U3-P6（rerank の窓）、U3-P7（エネルギー一致）、U5-P7（配座対仮説の廃止） | U5-P2 は U5-P1 の後 | S11、S13、S16、S18、HONO の FIND_PATH 回避 |
| S-C 探索と配座 | U4-P4（worker の数値処理）、U4-P6（窓の統一）、U4-P11（AFIR の限定）、U4-P7（生成物の同一性判定）、U2-P2（CREST の 1 規則、プローブ）、U2-P4（conformers の整理）、U1-P4（加算型の結合判定）、U1-P7（topology の責務） | U4-P7 は M4 の後、U1-P7 の削除部分は U7-P1 の後 | S6〜S8、S10、amine の irc_end_not_minimum の回復 |
| S-D 基盤と簡素化 | U9-P2（knob）、U9-P3（読まれないフィールド）、U9-P4（定義の一本化）、U9-P5（予算と行の順序）、U9-P7（シグナル）、U7-P5（GoodVibes の同一プロセス化）、U0-P4（原則と削除）、U1-P9（対象範囲の文書） | U9-P2 は U6-P5・U7-P1・U7-P5 の後 | 全テスト、lint、R1〜R4 のビット一致、S15（予算）、timeout 60 の孤児チェック |

### 7.3 could

U0-P13（def2-SVP method）、U1-P5（WL を安定まで）、U1-P6（spin_ok の相対許容幅）、U2-P6（剛体配置、U2-P2 のプローブ次第）、U3-P5（崩壊鞍点の freq 再利用）、U3-P8（MinimumPolicy 削除）、U3-P11（linopt 0）、U4-P9（記録の整理）、U5-P9（経路 SP を 1 ジョブ）、U5-P11（IDPP の衝突判定）、U6-P4（変位の規則を 1 つ）、U6-P11（BS-UKS の検出）、U7-P6（Eckart κ の列）、U8-P11（HTML の条件）、U8-P13（M06-2X 任意）、U9-P6（explore の並列）、U9-P9（Runtime.resolve）。

### 7.4 依存関係の要点

- U1-P1（元素表）→ U0-P1（ECP）、U4-P1（最大配位数の列）、U5-P11。
- U3-P2（鏡映）⇔ U5-P6（掌性を含む縮退判定）⇔ U7 の m：同時に入れる。
- U5-P1（分類器）→ U5-P2（NEB）→ U5-P4（1 チャンク）。U5-P1 → U7-P1（障壁なしの行）→ U8-P3（エネルギー層での低障壁の扱い）。
- U6-P2（方向）→ U6-P1（慣性）→ U6-P3・U6-P6（再試行と再開）。
- U5-P6 → U4-P5（縮退の発見を下流で受ける）→ U4-P7（縮退を合流させない）。
- U7-P1・U6-P5・U7-P5 → U9-P2（削除される knob）。
- フィールドを削除する改良（U7-P1、U9-P3 など）の検証は、複製した run の `--from` 再実行で行う。

---

## 8. 参考文献・マニュアル URL

**ライブラリのマニュアル・ソース**

- NWChem: ECP https://nwchemgit.github.io/ECP.html 、Basis https://nwchemgit.github.io/Basis.html 、DFT https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html 、Geometry optimization（DRIVER、moddir、inhess、linopt）https://nwchemgit.github.io/Geometry-Optimization.html 、Hessians and vibrational frequencies https://nwchemgit.github.io/Hessians-and-Vibrational-Frequencies.html 、TCE https://nwchemgit.github.io/TCE.html 、`src/driver/opt_drv.F`（driver_sad_search_dir、smalleig）https://github.com/nwchemgit/nwchem/blob/master/src/driver/opt_drv.F 、`src/basis/bas_input.F`、`src/hessian/hess_check.F`
- xTB: https://xtb-docs.readthedocs.io/en/latest/ 、spGFN（スピン非依存の GFN）https://xtb-docs.readthedocs.io/en/latest/spgfn.html
- CREST: https://crest-lab.github.io/crest-docs/ 、keywords https://crest-lab.github.io/crest-docs/page/documentation/keywords.html 、ソース `confparse.f90`・`cregen.f90`・`setuptest.f90`・`ztopology.f90`・`choose_settings.f90`（v3.0.2）https://github.com/crest-lab/crest
- pysisyphus: https://pysisyphus.readthedocs.io/en/latest/ 、TS optimization https://pysisyphus.readthedocs.io/en/latest/tsoptimization.html
- SCINE ReaDuct: https://scine.ethz.ch/download/readuct
- GoodVibes: https://github.com/patonlab/GoodVibes
- RDKit: https://www.rdkit.org/docs/
- NIST-JANAF（Ar）: https://janaf.nist.gov/tables/Ar-001.html

**理論・文献**

- 基底と ECP: F. Weigend, R. Ahlrichs, Phys. Chem. Chem. Phys. 7, 3297 (2005). 拡散関数: D. Rappoport, F. Furche, J. Chem. Phys. 133, 134105 (2010), doi:10.1063/1.3484283
- GFN2-xTB: C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 15, 1652 (2019), doi:10.1021/acs.jctc.8b01176
- 障壁の汎関数誤差: Y. Zhao, D. G. Truhlar, J. Phys. Chem. A 109, 2012 (2005), doi:10.1021/jp045141s。GMTKN55: L. Goerigk et al., Phys. Chem. Chem. Phys. 19, 32184 (2017), doi:10.1039/C7CP04913G。ωB97X-D3: Y.-S. Lin et al., J. Chem. Theory Comput. 9, 263 (2013), doi:10.1021/ct300715s
- 非局在化誤差: P. Mori-Sánchez, A. J. Cohen, W. Yang, Phys. Rev. Lett. 100, 146401 (2008)
- 共有結合半径・vdW 半径: B. Cordero et al., Dalton Trans. 2832 (2008)。S. Alvarez, Dalton Trans. 42, 8617 (2013)。加算型の結合判定: E. C. Meng, R. A. Lewis, J. Comput. Chem. 12, 891 (1991)
- WL の色の細分化: N. Shervashidze et al., J. Mach. Learn. Res. 12, 2539 (2011)
- 峠の定理: A. Ambrosetti, P. H. Rabinowitz, J. Funct. Anal. 14, 349 (1973)
- 鞍点探索: J. Baker, J. Comput. Chem. 7, 385 (1986)。J. M. Bofill, J. Comput. Chem. 15, 1 (1994)。CI-NEB: G. Henkelman, B. P. Uberuaga, H. Jónsson, J. Chem. Phys. 113, 9901 (2000)。ZTS: W. E, W. Ren, E. Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007)
- QRC: J. M. Goodman, M. A. Silva, Tetrahedron Lett. 44, 8233 (2003)
- 熱化学: S. Grimme, Chem. Eur. J. 18, 9955 (2012), doi:10.1002/chem.201200497（qRRHO）。R. F. Ribeiro et al., J. Phys. Chem. B 115, 14556 (2011)。経路縮重度: A. Fernández-Ramos et al., Theor. Chem. Acc. 118, 813 (2007), doi:10.1007/s00214-007-0328-0
- 変分 TST と振動断熱障壁: D. G. Truhlar, B. C. Garrett, S. J. Klippenstein, J. Phys. Chem. 100, 12771 (1996), doi:10.1021/jp953748q。Curtin–Hammett: J. I. Seeman, Chem. Rev. 83, 83 (1983), doi:10.1021/cr00054a001。energetic span: S. Kozuch, S. Shaik, Acc. Chem. Res. 44, 101 (2011), doi:10.1021/ar1000956
- トンネル: C. Eckart, Phys. Rev. 35, 1303 (1930)。H. S. Johnston, J. Heicklen, J. Phys. Chem. 66, 532 (1962)
- RKS の不安定性: R. Bauernschmitt, R. Ahlrichs, J. Chem. Phys. 104, 9047 (1996)。BS 法: L. Noodleman, J. Chem. Phys. 74, 5737 (1981)。スピン汚染の目安: D. Young, *Computational Chemistry* (Wiley, 2001)
- 自動反応探索の参考: YARP（Q. Zhao, B. M. Savoie）、ZStruct（P. M. Zimmerman）、Chemoton（J. P. Unsleber et al.）、autodE（T. A. Young et al.）、AFIR（S. Maeda et al.）
- 過去のレビュー: `docs/reviews/2026-09-26_calculation_spec_review_v4.md`（round-4 addendum を含む）、`docs/design.md`、`docs/validation.md`
