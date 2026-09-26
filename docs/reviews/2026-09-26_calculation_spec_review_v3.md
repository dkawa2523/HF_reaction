# hfauto 計算処理の仕様 再評価（v3）— 第2ラウンドの改良（M8、S18〜S25、C2）の後の評価・残る課題・次の改良（2026-09-26）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `d7ff974`。`hfauto/` のコードは `dc82145` と同じです（`d7ff974` は README.md と docs/validation.md だけを変更）。読み取り専用のレビューで、リポジトリのコード・設定・既存の文書は変更していません（追加したのはこのファイルだけ）。
- ベースライン: `docs/reviews/2026-09-26_calculation_spec_review_v2.md`（HEAD `573e80b`、以下「v2」）。その §1 の評価、§2 の課題の状態、§3 の単位ごとの課題と改良案（`U<n>-N<k>`、`U<n>-R<k>`）、§5 のロードマップ（must M8、should S18〜S25・C2、could C8〜C39）を基準にしました。さらに前の評価は `docs/reviews/2026-09-26_calculation_spec_review.md`（以下「v1」）です。
- 評価した改良:

| コミット | 内容 |
|---|---|
| `15b0beb`（R2-W1） | M8（FIND_PATH の端点の復元）、S19（画像の逐次整列）、S18（IDPP の反復上限 5,000）、S20（saddle の moddir の明示）、S21（初期 Hessian を渡す opt の trust 0.3）、C2（CCSD の maxiter 50）、S22（rankable）、S24（SMILES の同位体の拒否）、S25（入力切れの例外に上流の理由） |
| `dc82145`（R2-W2） | 文書: S23（composite の手順、手法誤差の数値、パネルの conditions）、M8・S18〜S25・C2 の記載 |
| `d7ff974` | docs/validation.md v3（WSL での実計算による再検証）と README のベンチマーク表 |

- 方法: v2 と同じく、計算単位 U0〜U8 ごとに 3 段で調べました。
  1. **現状仕様の整理（変更点が中心）**: コード、設定、改良前後の実 run の生ファイル（job.nw、job.json、stdout、result.json、manifest.json、diagnostics.json、log.jsonl）を読みました。改良前は v2（コード `4e972c5`、WSL `/home/user/hfauto_v2/<run>`）、改良後は v3（コード `dc82145`、`/home/user/hfauto_v3/<run>`）です。一部は v1 の前の W7（`/home/user/hfauto_w7/<run>`）とも比べました。
  2. **評価**: 化学的妥当性、有用性、複雑さを評価し、ライブラリの文書とソース（NWChem 7.2.3、CREST 3.0.2、ReaDuct 6.1.0、GoodVibes 4.3.0、RDKit）を調べました。必要に応じて小さなプローブ計算を行いました（WSL `/home/user/hfauto_v3_probe/`）。プローブの一部は他のプローブと同時に走ったので、所要時間は参考値です。
  3. **反証的な検証**: 別の担当が 2 の主張をコードと生データで確かめ、状態・重大度・優先度・評価を確定しました。
- 記載の規則:
  - 課題の状態、重大度、優先度、評価は、検証段の判定を使いました。
  - 全面的に棄却された課題はありませんでした。一部だけ正しいと判定された課題（U7-N11）は、訂正後の内容で書きました。
  - 改良案は、検証で「採用」または「修正して採用」とされたものだけを載せ、修正して採用されたものは修正後の形で書きました。評価段で見送った案は、各節の「見送った案」に理由付きで載せました。
  - 複数の単位で重複した改良は §5 で 1 つにまとめました。v2 のロードマップと同じ内容のものは v2 の ID（C12、C28 など）を使い、新しいものには M9、S26、C40〜C44 を振りました。修正して引き継いだものは `C27′`、`C39′` のように ′ を付けました。
  - ID の付け方: v2 までの課題は元の ID のまま。今回新しく見つかった課題は `U<n>-N<k>` の続き番号（U0-N8、U5-N8 など）と、U5 のテストと文書の課題 `U5-T1`、`U5-D1` です。単位ごとの改良案は `U<n>-R<k>` の続き番号です。

凡例:

| 記号 | 意味 |
|---|---|
| critical / high / medium / low | 課題の重大度 |
| must / should / could | 優先度。must は「結論を誤らせる」または「主要機能が働かない」もの |
| 【複雑さ減】/【複雑さ増】 | コード・設定・分岐が減る / 増える改良 |
| ◎ / ○ / △ / × | 評価。◎ 問題なし、○ 軽微な課題のみ、△ 改良が必要、× 目的どおりに機能していない。「○−」は ○ の下限 |
| resolved / partially / unresolved / not_in_scope | 課題の状態（解消 / 一部解消 / 未解消 / 本単位の範囲外・方針により対応しない） |

エネルギーの単位は、断りがなければ kcal/mol です。「v2 → v3」は第2ラウンドの改良前 → 改良後の実 run を指します。

---

## 0. 要旨

### 0.1 全体評価

1. **第2ラウンドの改良は計画どおりに入りました。** v2 のロードマップの must 1 件（M8）と should 9 件（S18〜S25、C2）はすべて実装され、WSL の実計算（v3）か本レビューのプローブで効果を確かめました。`hfauto/` の差分は 9 ファイル・正味 +33 行で、knob、失敗の種類、設定のフィールドは増えていません。pytest 291 件と smoke 12 件（71 s）が通りました。
2. **評価は 4 か所で上がり、下がった単位はありません（§1）。** △ と × はありません。
   - U0 全体プロトコル: 有用性 ○ → ◎（composite の手順が文書どおりに動き、手法誤差の数値が TZVPD 参照になった）
   - U1 structures: 妥当性 ○ → ◎（S24 で同位体を黙って ¹H にする経路が閉じた。U1 単体での評価）
   - U5 FIND_PATH: ○− → ○（M8 で 2 チャンク以上の string の端点が DFT 極小に戻った。ただし新しい U5-N8 が残る間は ○ の下限）
   - U6 鞍点・接続: 妥当性 ○ → ◎（U6-N1〜N3 の解消を、コード・NWChem のソース・実 run で確認）
3. **実計算での効果**（validation.md v3 と本レビューのプローブ）
   - HONO の known_endpoints: 9:32 → 5:23（−44%）。S18 で FIND_PATH の初期経路がねじれ経路になり、HEI の xTB Hessian が採用され、鞍点は 36 ステップ・106 s → 3 ステップ・6.7 s になった。
   - U6 区間（種の Hessian + saddle + TS freq + QRC）の合計: 438 → 181 s（−59%）。QRC は 3 反応で 204 → 107 ステップ、266 → 128 s（S21）。
   - STO-3G HCN の FIND_PATH: 各チャンクの string の端点が DFT 極小と 2e-9 Eh 以内で一致した（v2 の 2 チャンク目は +8.85 kcal/mol の非極小）。誤った `mode_overlap_below_0.3` の note も消えた（S19）。
   - ΔG の v2 からの差は ≤ 0.002 kcal/mol で、順位は 3 反応とも rank 1 のまま。
   - 主対象（TMA·(HF)₂）と水の報告行は `outcome:same_basin` だけになった（S22）。
   - C2: CCSD(T)/def2-TZVPD で HONO の TS が 18/50 反復で収束（余裕 32）。未収束時は NWChem が errquit し、NONZERO_EXIT として記録されないことをソースで確認した。
   - composite（S23）: v3 の run の複製に TZVPD の SP 層を追記すると、文書どおりに動いた（HCN 1:26、GoodVibes は 3/3 再利用）。
4. **新たに見つかった主な課題**
   - **(medium) 1 チャンク目の凍結端点のゆがみ（U5-N8）**: `_initial_path` が内部の画像を ends[0] へ 1 枚ずつ独立に整列するので、直線・平面の端点では初期経路の最後の区間に剛体の跳びが入ります。NWChem の zts_min_motion がこれに合わせて凍結 bead N を回し、L-BFGS の外挿で構造がゆがみます。端点の実エネルギーは真の極小より HONO で +1.09、STO-3G HCN で +13.6 kcal/mol 高くなりました（DFT SP で確認）。M8 は 2 チャンク目以降にしか効かないので、この問題を覆いません。誤分類の実例はありませんが、行 14（monotonic → BARRIERLESS）の根拠が障壁の小さい領域で崩れます。→ **M9**（2 行の置き換え）。
   - (low) 検証記録と文書のずれ: validation.md:59 の S21 の説明（「trust 0.1 でも歩幅が上限に届いていなかった」）は、v2 の stdout（歩幅の制限が 34 回）と矛盾する。validation.md:55・:137 はずれの原因を IDPP に帰している。validation.md:62・:96・:139 は C2 と手法パネルを「未確認」のまま。README・design.md・method_panel.yaml の「conditions が違えば全反応が thermo_unavailable」は S22 の後は順位外の outcome に当てはまらない。
   - (low) そのほか: 手法誤差の目安「約 2」が 2 反応の順位の比較には甘い（U0-N8）、S24 のメッセージの「natural-abundance masses」が実装（主同位体の質量）と合わない（U1-N8）、moddir の分岐と NWChem の smalleig（U6-N8・N9）、composite の report に energy 層の LOT が出ない（U7-N10）、ΔG_assoc が SVPD と TZVPD//SVPD で 2.49 違うのに検証記録に注記がない（U7-N11）、TSOpt の失敗原因の文字列が残らない（U4-N4）、ladder の継続で trust 0.1 に戻る（U3-N8）。
5. **検証の空白**: 次の経路は、まだ実 run で通っていません。
   - M1、S3（行 14）、行 13（中間体）、低レベル TS の近道、SCREEN の IDPP
   - 負モード 2 本以上の鞍点（S20 の Cartesian 分岐。プローブだけ）、higher_order、分割、QRC の再試行、`_register`、reassigned、periodic_nearest、DFT の mode-follow と soft 押し
   - M6 の副次的な負モード、電荷のある系・開殻系（structures → NWChem odft、会合量）、複数の温度、1 bar / 1 M
   - 手法パネル（sp stage → panel_report の経路）。C2 と composite は、本レビューのプローブでだけ確かめました。
6. **次に有用な改良（最小限）**
   - must 1 件: **M9**（FIND_PATH の初期経路を逐次整列にする。2 行の置き換えとテスト 1 件）。
   - should 2 件: **S26**（検証記録と文書の訂正。プローブの結果の記録、S21 の説明、conditions の文言など）、**C38**（WSL で code_version が常に `-dirty` になる問題。引数 2 つ）。
   - could は、どれも数行〜十数行で、knob は増やしません。複雑さを減らすもの（C34、C37、C27′）も含みます。

### 0.2 次の改良（上位 5 件）

| # | 改良 | 由来 | 期待効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| 1 | **M9** FIND_PATH の初期経路を `align_sequential([ends0, *inner, ends1])` にし、同じコミットで validation.md:55・:137 と design.md:130 を直す | U5-R8（修正して採用） | 1 チャンク目の凍結端点のずれが STO-3G 0.074 → 0.004 Å、HONO 0.022 → 0.0014 Å。端点の持ち上がりが +13.6 → +0.04、+1.09 → +0.009 kcal/mol。M8 と合わせて全チャンクで端点が極小になる。追加の計算はない | 2 行の置き換え、テスト約 10 行、文書 2 か所 | 中立 |
| 2 | **S26** 検証記録と文書の訂正（C2・composite・パネルのプローブ結果、S21 の説明、M6 の帰属、conditions の文言、thermo の根拠と制約） | U0-R3′、U3-R4、U7-R5′、U8-R7 ほか | design.md の数値の出典と、C2・S23 が動くことが文書から追える。trust が錯体に効かないという誤読を防ぐ。主対象の ΔG_assoc の LOT 依存（+2.49）が数値で分かる | 文書 15〜20 行 | 中立 |
| 3 | **C38** code_version の git status に `-c core.autocrlf=input` | U0-R6 | WSL の本番経路でも `-dirty` が本当の未コミット変更だけを示す（v2 の 9 run、v3 の 7 run がすべて -dirty。WSL で 83 → 0 行） | 引数 2 つ | 中立 |
| 4 | **C12 + C14′ + C41** explore の被覆と診断（WL クラスで drive の重複を除く、TS がある negative にも ΔE‡、TSOpt の例外の文字列を stderr に 1 行） | U4-R1、U4-R3、U4-R4 | 31 attempt 中 21 件が同じ C→N 1,2-H 移動に使われている状態が解消し、heavy_bond・association の段に届く。陰性の結論を障壁の値で裏付け、失敗原因をプローブの再実行なしに確かめられる | 20〜35 行、テスト 2〜3 件 | 中立 |
| 5 | **C34 + C37** 使われない分岐の削除（単量体の CREST の `--noreftopo` を「停止構造から既定の設定で 1 回」に置き換え、monitor 経路と FailureKind.STAGNATED を削除） | U2-R1、U8-R5 | 実 CREST で 2/2 失敗していた経路が 2/2 成功（v2 のプローブ）。到達しないコード約 30 行と失敗の種類 1 つ、フィールドと CLI フラグ 1 つずつが減る | 約 40 行（大半は削除）、テスト 4 か所 | 【複雑さ減】 |

### 0.3 第2ラウンドで良くなり、今のままでよい点（主なもの）

- M8 + S19: 成功したチャンクごとに、string の端点を DFT 極小（ctx.ends）に戻し、画像を逐次整列してから次のチャンクと `work.path` に渡す。bead のエネルギーは NWChem の値をそのまま使う。`end_drift` の note は監視だけでゲートにしない。
- S18: IDPP を収束させる（上限 5,000）。平面の端点どうしでも面外のねじれ経路になる。
- S20: 追うモードを `_moddir` 1 か所で決め、負モード 1 本なら autoz のまま `moddir 1`、2 本以上なら noautoz と P·H·P の昇順。NWChem 7.2.3 の射影子と並べ方と一致することを、ソースとプローブの stdout で確かめた。
- S21: 初期 Hessian を渡す opt だけ `trust 0.3` + `inhess 2`（NWChem の最小化の既定値）。Hessian のない opt（1 断片、`_register`、mode-follow、soft 押し）は `trust 0.1` のまま。
- C2: `maxiter 50` は ccsd ブロックだけに書き、キーには入れない。未収束は errquit → NONZERO_EXIT で保存されず、次の実行で計算し直される。
- S22: 順位を付けない outcome には `outcome:<o>` と本当の blocker だけを返し、dzpe の blocker は値があるときだけ付ける（G と ZPE は常に一緒に None になるので抜け道はない）。
- S23: composite は「run を複製して別の stage id で追記する」文書の手順で運用し、コードも knob も足さない。後の artifact が勝つマージと content-addressed な JobStore だけで成り立つ。
- S24 + S25: 同位体の SMILES は INPUT_INVALID（xyz の `D` とそろう）。stage の入力切れでは、上流の失敗の理由が例外の最終行に出る。
- JobStore のキーを化学的な要求に限る方針（ExecutionSpec と描画だけの入力は入れない）。デッキの詳細は `attempt_NN/job.nw` に残るので来歴は追える。

---

## 1. 評価の推移表

評価は「妥当性 / 有用性 / 複雑さ」の順です。複雑さの ◎ は「単純で過不足がない」を意味します。

| 単位 | v1 | v2 | v3 | 主な変化（v2 → v3） |
|---|---|---|---|---|
| U0 全体プロトコル | ○ / ○ / △ | ○ / ○ / ○ | ○ / **◎** / ○ | C2（CCSD maxiter 50、実エンジンで 18/50 反復）、S22（主対象の報告行）、S23（composite の手順が動く、手法誤差 HCN −1.2 / HONO +2.0 を再現）、S21 のデッキ規約。ΔE は SVPD で一致、ΔG‡ の差 ≤ 0.002 |
| U1 structures | ○ / ◎ / ◎ | ○ / ◎ / ◎ | **◎** / ◎ / ◎ | S24（同位体の SMILES を拒否）、S25（入力切れの理由）。SpeciesRecord は v2 と完全一致 |
| U2 conformers | △ / ○ / ○ | ○ / ○ / ◎ | ○ / ○ / ◎ | コード変更なし。tma_hf2 の CREST は揺らぎの範囲で一致（stage 19 → 19 s）。S24・S25 との境界は仕様どおり |
| U3 minima | ○ / ○ / △ | ◎ / ◎ / ○ | ◎ / ◎ / ○ | 本体のコード変更なし。S21 で錯体の opt が trust 0.3（A arm の opt 899 → 854 s、27 → 27 点、basin 3 つ）。S25 |
| U4 explore | △ / △ / ○ | ○ / ○ / ○（有用性は下限） | ○ / ○ / ○（有用性は下限） | コード変更なし。31 job の key が v2 とビット単位で一致、結果も同一 |
| U5 paths（仮説・SCREEN・FIND_PATH） | SCREEN ○ / FIND_PATH × / △ | ○（FIND_PATH ○−）/ ○ / ○ | ○（SCREEN ○、FIND_PATH **○**、下限）/ ○ / ○ | M8・S19（継ぎ目の端点が極小と 2e-9 Eh 以内）、S18（HONO のねじれ経路、U5 区間 208 → 171 s）。新しい U5-N8（1 チャンク目の端点のゆがみ）→ M9 |
| U6 paths（鞍点・検証・接続） | ○ / ◎ / ○ | ○ / ◎ / ○ | **◎** / ◎ / ○ | S20（4 件とも moddir 1、HONO の saddle 70 → 3 ステップ）、S21（QRC 204 → 107 ステップ）、M8・S19 の波及。U6 区間 438 → 181 s（−59%） |
| U7 thermo | ○ / ○ / △ | ◎ / ○ / ○ | ◎ / ○ / ○ | コード変更なし。S9・composite・中性系の S14 が初めて実計算（プローブ）で通った。ΔG_assoc は SVPD −11.53、TZVPD//SVPD −9.04 |
| U8 sp・report・実行基盤 | ○ / ○ / ○ | ◎ / ○ / ○ | ◎ / ○ / ○（有用性は ○ の上側） | C2、S22（HTML の差 43 B が消えた 2 語と一致）、S25。report は 0〜1 s。手法パネルは正式な検証で未実行 |

評価の留保（検証で付いたもの）:

- **U0 有用性 ◎**: 条件付きです。validation.md:62・:96・:139 の古い記述（C2 と S23 は未確認・文書だけ）と、WSL で code_version が常に `-dirty` になること（C38）が残っています。妥当性は、U0-N1（層が欠けた subject の ΔE の混成）、U0-N8（順位差の目安）、電荷のある系・開殻系と手法パネルの実 run がないことから ○ に据え置きました。
- **U1 妥当性 ◎**: U1 単体での評価です。電荷のある系・開殻系を structures → NWChem（odft）まで通した実 run はいまだ 0 本で、組成の多重度の高スピン固定（C8）も U2 に残っています。この 2 点は検証範囲の注記です。
- **U2 複雑さ ◎**: v2 と同じ但し書き（組成では到達せず、単量体では実 CREST で失敗する `--noreftopo` の分岐）が残ります。◎ / ◎ に上げる条件は C34、C8、U2-R2（C39′）の実施です。
- **U3 妥当性 ◎**: DFT の mode-follow（−側、ts_candidate、両側の再 settle、discovery）、DFT の soft 押し、M6 の |ν| 化は、v3 でも実エンジンで発火していません。
- **U4 有用性 ○**: 下限のままです。対象系で採用した生成物は 0 件、31 attempt 中 21 件が同じ C→N 1,2-H 移動、TS がある negative は ΔE‡ を持ちません。C12・C14′ で安定します。
- **U5 妥当性 ○**: FIND_PATH は ○ の下限です。U5-N8 が残る間、1 チャンクで結論が出る string の端点は、実エネルギーで最大 +13.6 kcal/mol ゆがみます。M9 の 1 行修正を入れて、はじめて ○ が確かになります。
- **U6 妥当性 ◎**: 実計算で通ったのは 3〜4 原子の共有結合系だけで、錯体の TS はまだ U6 に入っていません。負モード 2 本以上の分岐はプローブだけで確認しました。柔らかい負モードは NWChem の仕様（smalleig）で追えないことがあります（U6-N9）。どれも正しさではなく効率の問題で、VALIDATE_TS と connection のゲートが結果の誤りを防ぎます。
- **U8 有用性 ○**: 主対象の報告の誤読（U8-N3）はなくなりましたが、C25（ranking.csv の列）、U8-N4（手法パネルの正式な実行）、文書の遅れ（U8-N6、validation.md の 3 か所）が残ります。

---

## 2. 第2ラウンドの改良の実装と実計算での確認

### 2.1 改良ごとの実装と確認

| ID | 内容 | 実装箇所 | 実計算・プローブでの確認 | 判定 |
|---|---|---|---|---|
| M8 | FIND_PATH の成功したチャンクごとに、端点を ctx.ends（DFT 極小）に戻して次の initial と `work.path` にする。note `string{n}:c{i}:end_drift:<Å>`（監視だけ） | `drivers/reaction_case/actions.py:470-477` | find_path_sto3g（v3）: 3 チャンクの報告端点が HCN −92.10839487、HNC −92.06836342 で DFT 極小と 2e-9 Eh 以内（v2 の 2 チャンク目は +8.85 kcal/mol の非極小）。end_drift 0.074 / 0.001 / 0.002 Å。テスト `test_string_chunks_restore_the_frozen_ends_and_align_the_images` | 効果あり。1 チャンク目の中のずれは覆わない（U5-N8 → M9） |
| S19 | 画像を `align_sequential` で逐次整列する | `actions.py:476` | STO-3G の SaddleClaim の notes が `['mode_overlap_below_0.3']` → `[]`。接線と種の虚モードの重なり 0.118 → 0.653 | 効果あり |
| S18 | IDPP の反復上限 500 → 5,000 | `chemistry/interpolation.py:17` | HONO（v3）の初期経路がねじれ経路（二面角 180 → 0° を単調、H–O–N 103〜108°。v2 は面内の反転）。string の bead の最大は trans から 34.37 → 13.47。HEI と TS の写像 RMSD 0.325 → 0.009 Å。回帰テスト `test_idpp_leaves_the_plane_between_planar_cis_trans_hono`（500 反復では落ちる） | 効果あり |
| S20 | saddle で追うモードを明示（`_moddir`: 負モード 1 本なら autoz・moddir 1、2 本以上は noautoz と P·H·P の順、該当なしは moddir 0）。継続では moddir を min(k, 1)。key_payload に {moddir, cartesian} | `backends/nwchem/engine.py:87-95, 206-207, 356-367`、`chemistry/vibrations.py:104-114`、`input.py:135-142` | v3 の saddle 4 件はすべて autoz・`moddir 1`、AUTOZ の再試行と継続はなし。HONO は 34 + 36 ステップ（48.4 + 57.4 s）→ 3 ステップ（6.7 s）。負モード 2 本以上はプローブ `/home/user/hfauto_v3_probe/w1b_moddir` だけ（step 1 の追ったモードが P·H·P の第 k 固有ベクトルと |cos| 1.0000） | 効果あり。2 本以上の分岐は実 pipeline で未通過 |
| S21 | 初期 Hessian を渡す opt は `trust 0.3` + `inhess 2`、渡さない opt は `trust 0.1` | `backends/nwchem/input.py:124` | QRC（3 反応）204 → 107 ステップ、265.9 → 128.2 s。trajectory_above_ts 0 件、最終エネルギーの差 ≤ 8e-8 Eh、割付けは同じ。A arm の錯体の opt 899 → 854 s（27 → 27 点）。v3 の opt デッキは Hessian あり 13 本すべて trust 0.3、なし 12 本すべて trust 0.1 | 効果あり。錯体の opt への効果は小さい（§2.3） |
| C2 | CCSD ブロックに `maxiter 50`（mp2 には書かない） | `input.py:159, 166-167` | プローブ `/home/user/hfauto_v3_probe/u0_c2`（NWChemEngine を実エンジンで呼んだ）: デッキ `ccsd / freeze atomic / maxiter 50 / end`、出力 `maxit = 50`。HNC/TZVPD 11 反復・17.0 s、v3 HONO の TS/TZVPD 18 反復・447 s（同時負荷あり） | 効果あり（余裕 2 → 32 反復） |
| S22 | rankable: 順位を付けない outcome は `outcome:<o>` と本当の blocker だけ。dzpe の blocker は値があるときだけ | `chemistry/gates.py:297-317` | 水と TMA·(HF)₂ の ranking.csv の blockers が `outcome:same_basin;thermo_unavailable;dzpe_out_of_tolerance` → `outcome:same_basin`。report.html も同じ（差 43 B） | 効果あり |
| S23 | composite の手順（複製して別の stage id で追記）、手法誤差の数値（HCN −1.2、HONO +2.0、目安約 2）、パネルの conditions をそろえる注意 | `docs/design.md:147, 149, 150`、`README.md:37, 86`、`configs/pipelines/method_panel.yaml:5` | プローブ `/home/user/hfauto_v3_probe/u7_assess/hcn_tz`: 文書どおりに sp_tz → thermo_tz → report_tz を追記して 1:26、SP 3 本（83 s）、GoodVibes 3/3 再利用、ΔG‡ 41.99、rank 1。HONO の +2.0 は C2 のプローブで +2.02 と再現 | 効果あり（検証記録は未更新 → S26） |
| S24 | SMILES の同位体を拒否 | `chemistry/smiles.py:28-30` | プローブ `/home/user/hfauto_v3_probe/s24_s25`: `[2H]O[2H]` は INPUT_INVALID（「specifies isotopes」）。`[13CH4]`、`[1H]O[1H]` も拒否、`O`・`[OH-]` は通る | 効果あり（メッセージの文言は U1-N8） |
| S25 | stage の入力切れの例外に、上流の failed artifact の理由を最大 3 件 | `pipeline/runner.py:152-157` | 同プローブで最終行が `stage 'dft' (minima) has no input of type ['species']; upstream failures: species_d2o: SMILES '[2H]O[2H]' specifies isotopes; ...`、終了コード 1 | 効果あり（Rich のトレースバックが先に出るのは従来どおり） |

### 2.2 run ごとの比較（v2 → v3）

どれも新しい run dir（`/home/user/hfauto_v3/<run>`）での実行です。

| 系 | 所要時間 | 主な変化 | 値の差 |
|---|---|---|---|
| HCN（known_endpoints） | 1:41 → 1:24 | QRC 44 → 27 ステップ、41.7 → 25.5 s（S21）。鞍点 5 → 4 ステップ（S20） | ΔG‡ 42.1789 → 42.1783（−0.0006） |
| HONO | 9:32 → 5:23（−44%） | S18 で初期経路がねじれ経路に。string の最大 34.4 → 13.5。xTB Hessian が採用され、鞍点 36 ステップ・106 s → 3 ステップ・6.7 s（DFT Hessian と AUTOZ の再試行がなくなった）。QRC 70 → 30 ステップ、162 → 69 s | ΔG‡ 12.2237 → 12.2256（+0.0019）。ΔE‡ 13.67 は同じ |
| NH3 | 1:55 → 1:25 | QRC 90 → 50 ステップ（片側 25）、62 → 34 s | ΔG‡ 3.8407 → 3.8395（−0.0012） |
| 水（same_basin） | 0:11 → 0:12 | ranking.csv の blocker が `outcome:same_basin` だけ（S22） | — |
| STO-3G HCN の FIND_PATH | 1:03 → 1:17 | string 2 → 3 チャンク（+18 s）。各チャンクの端点が DFT 極小と 2e-9 Eh 以内 | ΔG‡ 62.1676 → 62.1678 |
| TMA·(HF)₂（discover） | 31:15 → 30:58 | 錯体の opt 8 → 8 ステップ、400 → 378 s。basin 3 つ、ジョブ 60 件（再利用 2）は同じ | ΔG_assoc −11.5266 → −11.5287（TMA 単量体の NWChem 再実行の数値ノイズ 5e-7 Eh による） |
| A/B の A arm（ab_init_hessian） | 40:34 → 39:51 | opt 27 → 27 点、899 → 854 s。最終エネルギーの差 ≤ 3e-8 Eh | — |

- ΔE（SVPD）は v2 と一致しました（HCN 46.6215 / 13.1797、HONO 13.6707 / 0.1315、NH3 4.2601 / 0.0、STO-3G 67.3273 / 25.1201）。ΔG‡ の差は S20・S18 による TS の構造のわずかな違いで説明できます。
- 接続の判定は 3 反応とも v2 と同じで、`trajectory_above_ts` は出ていません。
- STO-3G の 3 チャンク目は、2 チャンク目が反復 10〜13 で一度落ち着いた後、終わりの 3 反復が 0.103 / 0.109 / 0.084 kcal/mol で閾値 0.1 をわずかに超えたためです（判定がチャンクの終わりだけで行われる U5-N5 の副作用で、M8 の欠点ではありません）。

### 2.3 区間ごとの実測

**鞍点（S20）と TS**

| 系 | 種 | v2 | v3 | E_TS の差（v3 − v2） |
|---|---|---|---|---|
| HCN | screen_ts、xTB Hessian（−1427.2） | 5 ステップ、4.7 s | 4 ステップ、4.2 s | +1.9e-9 Eh |
| NH3 | screen_ts、xTB（−971.5） | 3 ステップ、2.9 s | 2 ステップ、2.5 s | −3.6e-10 Eh |
| STO-3G HCN | path_hei、種で DFT freq | 4 ステップ、1.9 s | 3 ステップ、1.8 s | +5e-10 Eh |
| HONO | path_hei | xTB Hessian 不採用（虚振動 2 本）→ 種で DFT freq 15.1 s、AUTOZ 失敗 34 ステップ 48.4 s、Cartesian 36 ステップ 57.4 s | xTB Hessian 採用（虚振動 1 本、−744.4）、3 ステップ 6.7 s | −2.5e-7 Eh |

HONO の種の高さは、TS から +21.9（v2、種の DFT freq）→ +0.05 kcal/mol（v3、saddle の step 0）になりました。

**QRC（S21）**

| 系 | ステップ（片側 / 片側） | 時間 | 負モードで刻みを制限されたステップ |
|---|---|---|---|
| HCN | 20 / 24 → 12 / 15 | 41.7 → 25.5 s | 13 → 9 |
| HONO | 34 / 36 → 15 / 15 | 162.1 → 69.1 s | 61 → 20 |
| NH3 | 45 / 45 → 25 / 25 | 62.1 → 33.6 s | 76 → 45 |
| STO-3G HCN | 30 / 22 → 15 / 14 | 11.9 → 7.9 s | 22 → 9 |

変位直後の降下は 2.72〜3.45e-4 Eh（目標 3e-4、drop 1e-5）で、8 側すべてで max(traj) = traj[0] でした。

**U6 区間の壁時間**（種の Hessian + saddle + TS freq + QRC）: HCN 52.2 → 35.7 s、HONO 298.3 → 91.2 s、NH3 70.1 → 41.2 s、STO-3G 17.5 → 13.4 s。合計 438 → 181 s（−59%）。

**U5 区間**（SCREEN の xTB opt・GS・DFT SP・tangent freq と FIND_PATH の string）: HCN 27.8 → 26.8 s、NH3 28.2 → 28.1 s、HONO 208.4 → 170.9 s（string 1 チャンク 200.4 → 162.9 s）。HONO の string は 6 反復目で落ち着いたのに、残り 14 反復（約 96 s）が判定されずに走りました（U5-N5）。

**minima（S21 の U3 への波及）**: 錯体の opt 1 本あたり約 20 s（−5%）短くなりました。歩数は変わりません。v2 の stdout では NWChem が 4 錯体の opt で歩幅を 34 回制限し（`Restricting large step in mode`）、全体の縮小も 7 回ありました（v3 はそれぞれ 3 回と 2 回）。歩数が同じなのは、v2 では line search が制限された歩みを 1.2〜1.7 倍に伸ばして補っていたことと、TMA·(HF)₂ の 4〜8 歩目が勾配の収束後に XMAX / XRMS を満たすための小さな歩みだからです（U3-N7）。同じ軌跡をたどった run どうしの差は約 2% 以内で、5% の短縮はそれより大きい観測です（原因の「後半の SCF が軽い」は推定）。dft stage では A arm 2,433 → 2,389 s（−1.8%）、tma_hf2_discover 1,771 → 1,754 s（−1.0%）で、時間の大半は DFT の解析 freq です。

### 2.4 v2 ロードマップの実施状況

| ID | 内容 | 実装 | 実計算での確認 | 備考 |
|---|---|---|---|---|
| M8 | 端点の復元 | 済 | find_path_sto3g、HONO | 1 チャンク目の中は U5-N8 → M9 |
| S18 | IDPP の上限 | 済 | HONO | — |
| S19 | 画像の逐次整列 | 済 | STO-3G | 初期経路（`_initial_path`）は対象外だった → M9 |
| S20 | moddir の明示 | 済 | 4 件（1 本の分岐）、プローブ（2 本以上） | 2 本以上の分岐の回帰テストがない → C42 |
| S21 | trust 0.3 | 済 | QRC 3 反応 + STO-3G、A arm、TMA | validation.md:59 の説明が誤り → S26 |
| S22 | rankable | 済 | 水、TMA | conditions の文言がずれた（U8-N6）→ S26 |
| S23 | composite・手法誤差・conditions の文書 | 済 | プローブ（hcn_tz、tma_a、tma_tz、tma_d） | 検証記録が未更新 → S26 |
| S24 | 同位体の拒否 | 済 | プローブ | メッセージの文言 → C40 |
| S25 | 入力切れの理由 | 済 | プローブ | — |
| C2 | CCSD の maxiter 50 | 済 | プローブ（u0_c2） | sp stage 経由は未通過 → C43 |
| C39 | 文書の訂正 | 一部 | — | validation.md §8 の GS の旧対策案の取り下げ、engine の continuation の docstring、design.md:130 の M8 の記載は済。残りは C39′ と S26 |
| C8、C10、C12〜C16、C20、C21、C25、C27〜C37 | could | 未 | — | §5.3 に引き継ぎ（C21 は文書化に、C27 は削除に置き換え、C35 は見送り） |
| C1、C4 のラジカルの部分、C18、C22 | — | 未 | — | v2 の判断のまま（見送り・据え置き） |

横断課題のうち、X4（手法パネル）と X5（composite G）は、単位の判定（U0-I8・U8-3、U0-I3・U8-5 がいずれも resolved）により実質的に解消しました。X9（開殻と多重度）は S24 で同位体の部分が閉じ、C8 が残ります。

### 2.5 本レビューのプローブ（validation の外で確かめたこと）

すべて WSL の `/home/user/hfauto_v3_probe/` の下で、v3 の run の読み取りか、run dir の複製への追記で行いました。v3 の run そのものは変更していません。

| プローブ | 確かめたこと |
|---|---|
| `u0_c2` | C2 の実デッキと NWChem のエコー（`maxit = 50`）、HNC と v3 HONO の TS の CCSD(T)/def2-TZVPD（11 / 18 反復）。CCSD(T) の ΔE‡ 11.65、SVPD − CCSD(T) = +2.02 |
| `s24_s25` | 同位体の SMILES と xyz の `D` の拒否、入力切れの例外の上流の理由 |
| `u1_review`、`u1_assess` | xyz の欠損が dry-run を素通りし、実行で structures が `running` のまま残ること（U1-N2）。RDKit の RemoveHs が同位体付きの H を残すこと |
| `u2_s24` | 同位体で弾かれた部品を含む組成だけが `missing_component` で止まり、ほかは続くこと |
| `u3_review` | v2 / v3 の錯体の opt の stdout の比較（歩幅の制限の回数、line search、walltime） |
| `u4_spec`、`u4_review` | explore の 31 job の比較と、TSOpt の失敗の再現（NoNegativeEigenValueException） |
| `u5_spec`、`u5_assess` | 初期経路の剛体の跳び、逐次整列した初期経路での string の再実行、凍結端点の DFT SP 8 本（U5-N8） |
| `u7_spec`、`u7_assess`（`hcn_tz`、`tma_d`、`tma_a`、`tma_tz`） | composite の手順（S23）、energy_layer_missing（S9）、energy 層を含む中性系の会合量（S14）、TZVPD//SVPD の ΔG_assoc |
| `w1b_moddir`（第2ラウンドの実装時） | 負モード 2 本以上の saddle の Cartesian 分岐（S20） |

### 2.6 実計算でまだ通っていない経路

| 経路 | 状態 |
|---|---|
| M1（陰性結果ゲートの削除）、S3（行 14）、行 13（中間体） | 今回の系に経路に入るものがない |
| 負モード 2 本以上の saddle（S20 の Cartesian 分岐） | プローブだけ（C42 で回帰テスト） |
| higher_order、分割、QRC の再試行、`_register`、reassigned、periodic_nearest | 未到達 |
| DFT の mode-follow・soft 押し、M6 の副次的な負モード | 未発火 |
| 電荷のある系・開殻系（structures → NWChem odft、会合量） | 実 run 0 本（v3 のデッキ 80 本はすべて charge 0 / mult 1） |
| 手法パネル（sp stage → panel_report） | v2 のプローブだけ（C43） |
| 複数の温度、1 bar / 1 M | 未実行 |

---

## 3. v2 の課題の解消状況

単位ごとに再評価した延べ 121 件（v2 で新しく見つかった課題と、v1 から持ち越した課題。同じ課題が 2 つの単位に載るものは両方で数えた）の内訳は、resolved 54、partially 6、unresolved 57、not_in_scope 4 です。v2 で high・medium とした新しい課題のうち、U5-N1 = U6-N3（high）と U5-N2・U6-N1（medium）は resolved です。U5-N3（medium、DLC の GS の脆さ）だけは unresolved のままですが、S18 の後は影響がコストに限られるようになりました。unresolved の 57 件は、ほとんどが v2 で could とした改良が未実施のものです。

第2ラウンドで状態が変わった課題:

| ID | v2 → v3 | 決め手 |
|---|---|---|
| U0-I3 composite が配線されていない | partially → resolved | S23 の手順が実 run の複製で動いた（C1 は見送りのまま） |
| U0-I8 = U8-3 CCSD の maxiter | unresolved → resolved | C2 |
| U0-N2・N3・N4・N5・N6 | 新規 → resolved | S23、C2、S22 |
| U1-03 SMILES の同位体・ラジカル | unresolved → partially | S24（ラジカルは見送り済み） |
| U1-N1 入力切れの理由が出ない | 新規 → partially | S25（トレースバックが先に出る点は残る） |
| U5-N1 = U6-N3 凍結端点のずれを次のチャンクに渡す（high） | 新規 → resolved | M8 |
| U5-N2 IDPP の打ち切り（medium） | 新規 → resolved | S18 |
| U6-I1 QRC の Hessian と trust | partially → resolved | S21 |
| U6-I2 mode_index と moddir の番号 | unresolved → resolved | S20 |
| U6-N1 moddir 0 の迷走（medium）、U6-N2 画像の未整列、U6-N4 trust 0.1 の退行、U6-N5 docstring | 新規 → resolved | S20、S19、S21、C39 の一部 |
| U8-5 composite と thermo の上書き | partially → resolved | S23 |
| U8-8 TS のない outcome の dzpe blocker、U8-N3 thermo_unavailable | unresolved / 新規 → resolved | S22 |
| U8-N1 composite の手順が動かない | 新規 → resolved | S23 |

### 3.1 U0 全体プロトコル

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U0-I3 composite が配線されていない | partially | resolved | design.md:150 の手順どおりに `u7_assess/hcn_tz` で動いた（sp_tz 3 SP・83 s、thermo_tz は GV 3/3 再利用、report_tz rank 1、ΔG‡ 41.9855） |
| U0-I7 同順位の幅に手法誤差が入らない | resolved | resolved | design.md:149・README.md:86 の数値を TZVPD 参照（HCN −1.2、HONO +2.0、目安約 2）に置き換えた。目安としての甘さは U0-N8 |
| U0-I8 CCSD の maxiter 20 | unresolved | resolved | input.py:166。`u0_c2` の job.nw と stdout（`maxit = 50`、11/50 と 18/50 反復） |
| U0-N1 層が欠けると ΔE が freq 層との混成 | 新規 | unresolved | C28 は未実施。thermochemistry.py:190-193 は v2 のまま |
| U0-N2 composite の文書手順が実行できない | 新規 | resolved | S23。U0-I3 と同じ証拠 |
| U0-N3 手法誤差の数値が TZVP 参照 | 新規 | resolved | (−205.4237916 − (−205.442352)) × 627.5095 = 11.65、SVPD の 13.6707 との差 +2.02 |
| U0-N4 CCSD(T) の opt-in に対して余裕が小さい | 新規 | resolved | HONO の TS/TZVPD で 18/50 反復。未収束は NONZERO_EXIT で保存されない（output.py:231-243、jobstore.py:155-157） |
| U0-N5 パネルの conditions の食い違い | 新規 | resolved | S23 の文書対応（method_panel.yaml:5、README.md:37、design.md:147）。文言は S22 の後に不正確（U8-N6） |
| U0-N6 thermo がないとき dzpe の blocker | 新規 | resolved | gates.py:307-316。v3 の水と TMA の ranking.csv |
| U0-N7 WSL で code_version が常に -dirty | 新規 | unresolved | C38 は未実施。WSL で 83 行、`-c core.autocrlf=input` で 0 行。v3 の 7 run すべて `dc82145…-dirty` |

### 3.2 U1 structures

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U1-01 SMILES の端点で H の写像が決まらない | resolved | resolved | system.py:88-91（M7）は変更なし |
| U1-02 宣言座標の添字 | resolved | resolved | system.py:48-55、structures.py:88-89 は変更なし |
| U1-03 SMILES の同位体・ラジカル | unresolved | partially | S24 で同位体の部分を閉じた（smiles.py:28-30、test_smiles_input.py:17-18）。ラジカルの照合は v2 で見送り済み |
| U1-04 id の衝突 | resolved | resolved | 変更なし |
| U1-05 未知元素 | resolved | resolved | `s24_s25/run_xyz.log` で `unknown element(s): D` |
| U1-06 状態ラベルの限界 | unresolved | unresolved | topology.py は変更なし。design.md:141 に限界の記述がない |
| U1-07 structures の規則の文書 | resolved | resolved | design.md:138 に「同位体を指定した SMILES は INPUT_INVALID」 |
| U1-08 組成キー | resolved | resolved | 変更なし |
| U1-N1 入力切れの理由が出ない | 新規 | partially | S25 で最終行に理由が出る。Rich のトレースバック（約 44 行）が先に出て、終了コード 1、report.html は作られない点は残る |
| U1-N2 xyz の欠損が dry-run を素通り、`running` の残留 | 新規 | unresolved | `u1_review` で再現（run1.log と run2.log がバイト単位で同一）。C31 は未実施 |
| U1-N3 片側の端点だけ失敗したときの cascade | 新規 | unresolved | structures.py は変更なし |
| U1-N4 無駄な DFT（両端の q/m の不一致） | 新規 | unresolved | C31(b) は未実施 |
| U1-N5 状態ラベルの限界が文書にない | 新規 | unresolved | design.md:141 に記述なし（C39′） |
| U1-N6 spec の行番号 | 新規 | resolved | v3 の行番号は HEAD と一致 |
| U1-N7 hypotheses.py:158-159 の到達しない ValueError | 新規 | unresolved | 今回は削除せず防御コードとして残す方針に変更（§4.1.8） |

### 3.3 U2 conformers

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U2-ISS-1 組成の CREST のトポロジー停止 | resolved | resolved | crest.py:56-58 の `--noopt` は変更なし。v3 の組成は rc 0、'Initial Geometry Optimization' 0 回 |
| U2-ISS-2 4 kcal/mol の窓 | resolved | resolved | 窓なし |
| U2-ISS-3 組成の多重度を黙って高スピンに | partially | partially | 文書化（design.md:138 の structures 行とモジュールの docstring:7）だけ。conformer_search.py:202 はそのまま（C8） |
| U2-ISS-4 効かない knob | resolved | resolved | 4 knob、extra=forbid、スレッドは site の 1 か所 |
| U2-ISS-5 `--nci --quick` の実効設定 | resolved | resolved | CREST 3.0.2 の confparse.f90 で再確認（後の `-ewin 6` が quick の 5 を上書き）。v3 の stdout に `EWIN: 6.0000` |
| U2-ISS-6 seed00 だけを渡す | resolved | resolved | 変更なし |
| U2-ISS-7 組成の停止構造が流れる | resolved | resolved | 停止そのものが起きない |
| U2-N1 単量体の `--noreftopo` 再実行が PT 型の停止で失敗 | 新規 | unresolved | C34 は未実施。v3 では未到達 |
| U2-N2 開殻の組成が `--noopt` で失敗 | 新規 | unresolved | design.md:139 に注記なし |
| U2-N3 environment.md:64 の「約 25 s」 | 新規 | unresolved | 残っている（実測は約 1.4 s） |
| U2-N4 ConformersConfig の範囲検査 | 新規 | unresolved | C30 は未実施 |
| U2-N5 crest の threads を省くと `-T 1` | 新規 | unresolved | 既定を変えない判断は v2 で済み。文書の注記（C39′）が残る |

### 3.4 U3 minima

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U3-ISS-1 同じ basin の 2 本目の freq | resolved | resolved | 本体のファイルに差分なし。A arm の tma_hf2_afir は v3 でも `known:basin_tma_hf2_705342ba` |
| U3-ISS-2 explore が入力構造から出発 | resolved | resolved | M3 は変更なし |
| U3-ISS-3 thermo が noise / soft の負モードを捨てる | resolved | resolved | M6 は変更なし（実データでの発火はまだない） |
| U3-ISS-4 同一性の基準が 2 種類 | resolved | resolved | identity.py に差分なし |
| U3-ISS-5 window で screen 極小 0 件なら黙って空 | unresolved | unresolved | S25 が止めるのは consume する型が 0 件のときだけで、SPECIES があれば通る（C10） |
| U3-ISS-6 振動数の射影の約束 | resolved | resolved | 変更なし |
| U3-ISS-7 relax 中の known.find が組成を見ない | partially | partially | 変更なし。エネルギー条件で分離されるので「変更不要」を維持 |
| U3-ISS-8 開殻一重項を RKS で扱う | resolved | resolved | 変更なし |
| U3-N1 rerank_sp の SP | 新規 | unresolved | v3 でも SP 3 本・49.8 s（C32） |
| U3-N2 `_register` が gates を渡さない | 新規 | unresolved | actions.py:347-348 は変更なし（C33） |
| U3-N3 DFT mode-follow の再 settle | 新規 | unresolved | 発火 0 回。発火してから入れる判断を維持 |
| U3-N4 design.md:142 の include: all の文言 | 新規 | unresolved | trust の部分だけ直った（C10） |
| U3-N5 rerank の always の非対称 | 新規 | unresolved | 設計上の性質として見送り |
| U3-N6 assign のエネルギー判定が最近傍だけ | 新規 | unresolved | 設計上の性質として記録 |

S21 の U3 への採用条件（新しい run dir で A arm と known_endpoints を再実行し、basin・エネルギー・同定を確かめて記録する）は満たしています（resolved）。ただし記録の仕組みの説明が誤っています（U3-N7）。

### 3.5 U4 explore

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U4-I1 陰性結果の拒否権（M1） | resolved | resolved | explore のコードに差分なし |
| U4-I2 出発構造（M3） | resolved | resolved | v3 の source の opt は 31/31 件が 2 反復・ΔE 0.000 kJ/mol |
| U4-I3 RMSD による一致判定（S6） | resolved | resolved | 縮退 HF 交換（e30f4cd2）を正しく same_as_source と判定 |
| U4-I4 等価な drive で上限が埋まる | unresolved | unresolved | 31 attempt 中 21 件が同じ C→N 1,2-H 移動（C12） |
| U4-I5 ReaDuct の反復上限 150 | unresolved | unresolved | IRC 22/22 本が 150 で打ち切り。無バイアスの opt は最大 109 反復で、TMA 系では実害なし（C13） |
| U4-I6 陰性ラベルの粗さ | unresolved | unresolved | TS がある negative の barrier_kj_mol は None（C14′） |
| U4-I7 AFIR γ=300 の害 | not_in_scope | not_in_scope | 害の実例なし |
| U4-I8 F→N 移動の drive | not_in_scope | not_in_scope | 設計上の判断 |
| U4-N1 validation.md の −289i の判定変化を S6 の効果としている | 新規 | unresolved | validation.md:91 に残る（C39′） |
| U4-N2 状態ラベルの閾値の近くの H 結合 | 新規 | unresolved | 受け入れたリスク。v3 のラベル反転 0 件 |
| U4-N3 AFIR の reason が γ ごとに上書き | 新規 | unresolved | worker.py の `_afir` は変更なし |

v2 の数値の訂正（「2 反復は v2 の 71 attempt だけ」）は、v3 の仕様書が範囲を正しく書いているので resolved です。

### 3.6 U5 paths（仮説・SCREEN・FIND_PATH）

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U5-I7 両端より低い DFT//xTB ノード | unresolved | unresolved | gates.py は変更なし（C16） |
| U5-I8 tangent 用の DFT freq | unresolved | unresolved | NH3 で 5.1 s（C15） |
| U5-N1 凍結端点のずれを次のチャンクに渡す（high） | 新規 | resolved | M8。1 チャンク目の中のずれは U5-N8 として切り出した |
| U5-N2 IDPP の 500 反復打ち切り（medium） | 新規 | resolved | S18 |
| U5-N3 DLC の GS の脆さ（medium） | 新規 | unresolved | v3 の HONO の GS は同じキーで同じく失敗。S18 で退避経路が正しい機構を捉えるので、影響はコスト（GS 8 s + string 163 s）だけ |
| U5-N4 S4 が GS 自体の失敗でも再実行 | 新規 | unresolved | worker.py は変更なし（C36） |
| U5-N5 判定がチャンクの終わりだけ | 新規 | unresolved | HONO で約 96 s の無駄、STO-3G で余分な 1 チャンク |
| U5-N6 design.md:130 の gmax の記述 | 新規 | unresolved | 残っている（C39′） |
| U5-N7 試行を使い切った後の FIND_PATH | 新規 | unresolved | 見送り（v2） |

v2 の U5 の改良案のうち、U5-R1（M8）と U5-R2（S18）は resolved、U5-R5（C39 の U5 分）は partially（validation.md の旧対策案の取り下げと design.md:130 の M8 の記載は済、gmax の記述は残る）、U5-R3（C35）は今回見送り、U5-R4・R6・R7 は unresolved です。

### 3.7 U6 paths（鞍点・検証・接続）

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U6-I1 QRC に Hessian を渡さない・trust | partially | resolved | S7 + S21。QRC デッキ 8 件すべて `trust 0.3 inhess 2` |
| U6-I2 mode_index と moddir の番号 | unresolved | resolved | S20。NWChem の opt_drv.F の射影子と昇順の対角化と一致。プローブで追ったモードと P·H·P の固有ベクトルが |cos| 1.0000 |
| U6-I3 periodic_nearest に許容幅がない | unresolved | unresolved | C18 は未実施。実 run で未到達 |
| U6-I4 QRC の振幅のクリップ順 | resolved | resolved | v3 はすべて 1 回目の振幅で接続 |
| U6-I5 `_register` がゲートの前 | unresolved | unresolved | 設計どおり、発動なし |
| U6-I6 torsional の閾値 | not_in_scope | not_in_scope | hypotheses の担当 |
| U6-I7 試行を使い切った後の FIND_PATH | unresolved | unresolved | = U5-N7 |
| U6-N1 Cartesian で moddir 0 の迷走（medium） | 新規 | resolved | moddir 0 は mode_index が None のときだけ。v3 の saddle 4 件は moddir 1。元の状況（autoz 失敗後の Cartesian）は実 run で再現しておらず、プローブで確認。柔らかい負モードでは同様の迷走がありうる（U6-N9） |
| U6-N2 string 画像の未整列 | 新規 | resolved | S19。重なり 0.653 |
| U6-N3 凍結端点のずれ | 新規 | resolved | M8 |
| U6-N4 trust 0.1 による退行 | 新規 | resolved | S21。validation.md:69 の「144」の誤記は残る（C39′） |
| U6-N5 continuation の docstring | 新規 | resolved | engine.py:185-187 |
| U6-N6 qrc_bounds_A の上限と HESSIAN_NEAR_A | 新規 | unresolved | 文書に記述なし（C39′） |
| U6-N7 子反応の予算 | 新規 | unresolved | design.md:113・:137 は「1 反応 6 h」のまま（C39′） |

### 3.8 U7 thermo

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U7-I1 負モードを捨てる（M6） | resolved | resolved | コード変更 0 行。v3 の TS はどれも負モード 1 本 |
| U7-I2 会合量の LOT とキー | partially | partially | S14 の state=False が `tma_a` で実 run を通った（−11.52885）。None の理由は今も出ない |
| U7-I3 経路縮重度 | resolved | resolved | design.md:146 |
| U7-I4 zpe_scale（S15） | resolved | resolved | GoodVibes 4.3.0 の calc_bbe で ZPE と U_vib の因子がそろうことを確認 |
| U7-I5 ΔG_assoc の感度幅 | unresolved | unresolved | C21 は文書化に置き換え（§5） |
| U7-I6 QH=False の記載 | resolved | resolved | design.md:146 |
| U7-I7 点群・σ の記録 | unresolved | unresolved | C20 は未実施（NH3 の ΔG‡ の差が −0.0013 で σ は安定、実害なし） |
| U7-I8 基準が反応物の極小 | unresolved | unresolved | 方針どおり見送り |
| U7-I9 書くだけの出力 | partially | partially | population・H・S_rot は書くだけ、ThermoResult.notes は常に () |
| U7-I10 トンネル補正なし | resolved | resolved | design.md:146 |
| U7-I11 単原子種 | unresolved | unresolved | 方針どおり対応しない。文書にもない |
| U7-N1（= U0-N1）混成 ΔE | 新規 | unresolved | S23 の手順（targets: all_minima）では SP が失敗しない限り起きない（C28） |
| U7-N2 ΔG_assoc=None の理由が出ない | 新規 | unresolved | `tma_d` で無言の None を再現（C29） |
| U7-N3 (Hill, 電荷) の後勝ち | 新規 | unresolved | 変更なし（C29） |
| U7-N4 Conditions・ThermoSettings の範囲検査 | 新規 | unresolved | C30 |
| U7-N5 vib_scale の根拠と記録 | 新規 | unresolved | design.md:146 に 0.985 の記述なし（S26） |
| U7-N6 ThermoResult.notes・S_rot が使われない | 新規 | unresolved | C20 |
| U7-N7 validation の帰属「M6・S15」 | 新規 | unresolved | validation.md:76 に残る（S26） |
| U7-N8 反応のない錯体では ΔG_assoc が出ない | 新規 | unresolved | 文書に制約として書かれていない（S26） |
| U7-N9 単量体アンサンブルに異性体 | 新規 | unresolved | 変更なし（S26 で文書化） |

v2 §3.7.2 の「残る弱点 4（M6・S14・S9 が実計算で未通過）」は partially になりました（S9 と中性系の S14 が実 run の複製で通った。M6 の副次モードと電荷・開殻の会合量は未通過）。

### 3.9 U8 sp・report・実行基盤

| ID | v2 | v3 | 根拠 |
|---|---|---|---|
| U8-1 符号ゲートで HONO が外れる | resolved | resolved | summary.py の不感帯は変更なし |
| U8-2 既定パネルの CCSD(T) | resolved | resolved | opt-in のまま |
| U8-3 CCSD の maxiter | unresolved | resolved | C2 |
| U8-4 UHF-MP2 が常に失敗 | resolved | resolved | S17 は変更なし |
| U8-5 composite と thermo の上書き | partially | resolved | S23。Manifest.union は後の artifact を優先（manifest.py:80-91） |
| U8-6 同順位の幅の意味 | resolved | resolved | README と design の数値を更新 |
| U8-7 ranking.csv に T・標準状態・ΔG_rxn がない | unresolved | unresolved | C25 |
| U8-8 TS のない outcome の dzpe blocker | unresolved | resolved | S22 |
| U8-9 site.memory_mb が未使用 | unresolved | unresolved | C27′ |
| U8-10 method_panel.yaml が独自の conditions | unresolved | unresolved | 受容済み。注意書きの文言は S22 の後に不正確（U8-N6） |
| 2.8.3(e) hfauto report の二重出力 | not_in_scope | not_in_scope | 見送りのまま |
| U8-N1 composite の手順が動かない | 新規 | resolved | S23 |
| U8-N2 monitor 経路と STAGNATED | 新規 | unresolved | `def monitor` 7 か所が残る（C37） |
| U8-N3 TS のない outcome に thermo_unavailable | 新規 | resolved | S22 |
| U8-N4 手法パネルが正式な検証で未実行 | 新規 | unresolved | v3 にも panel_sp・panel_report はない（C43） |
| U8-N5 metric=dG_rxn で並べた値が CSV に出ない | 新規 | unresolved | C25 |

---

## 4. 単位ごとの再評価

各単位の「追加の改良案」には、§5 のロードマップ ID を【 】で添えました。

### 4.0 U0 全体プロトコルと理論レベルの選択

評価: v2 ○ / ○ / ○ → **v3 ○ / ◎ / ○**

#### 4.0.1 変更点

- **C2**: input.py:166 で CCSD ブロックに `  maxiter 50` を書く（mp2 には書かない）。JobStore のキー（method・molecule）は変わらない。
- **S22**: gates.py:297-317 の rankable を書き直した。順位を付けない outcome は `outcome:<o>` と本当の blocker（thermo_unavailable、mixed_level_of_theory、spin_contaminated、method_sign_disagreement）だけ。dzpe が None なら dzpe_out_of_tolerance を付けない。
- **S23（文書だけ）**: design.md:149（手法誤差を HCN −1.2 / HONO +2.0、目安約 2 に。NH3・±1.7・約 1.5 の記述を削除）、design.md:150（composite を「cp -a で複製し、sp（targets: all_minima）→ thermo（energy_method）→ report を別の stage id で追記」に）、design.md:147・README.md:37・method_panel.yaml:5（conditions を元の run とそろえる注意）、README.md:86（約 1.5 → 約 2）。
- **S21 の波及**: 初期 Hessian ありの opt は `trust 0.3` + `inhess 2`、なしは `trust 0.1`。trust はジョブキーに入らない（engine.py:342-344）。
- 変わっていないもの: 3 層構成、method ファイル 5 本、ThermoSettings、`_check_stages`、符号ゲートの不感帯、S9、S17、パネルの既定の 2 手法。thermochemistry.py・summary.py・core/method.py・pipeline/config.py・configs/methods は 4e972c5 → dc82145 で差分なし。
- 未実施: C28（U0-R4）、C38（U0-R6）、U0-R3（validation.md へのプローブ記録）。C1 は見送りのまま。

#### 4.0.2 妥当性（○、据え置き）

- 3 層構成と停留点の単一 PES（SVPD）は v2 のままで妥当です。SVPD の ΔE は v2 と一致し、ΔG‡ の差（≤ 0.002）は S20 による TS の構造のわずかな違いで説明できます。
- **S21 のデッキ**は NWChem の公式の既定値（最小化の TRUST 0.3、saddle 0.1）と INHESS 2 の意味（restart データを優先し、なければ freq の Cartesian Hessian を変換）に整合します。継続（maxiter・timeout）では Hessian が None になり `trust 0.1`・inhess なしになりますが、driver は既定の INHESS 0 で drv.hess を読むので、動作に矛盾はありません（U3-N8 で文書化を提案）。
- **C2 を実エンジンで確認しました。** NWChem v7.2.3 の ccsd_input.F は `maxiter` を rtdb の `ccsd:maxiter` に入れ、プローブの出力にも `maxit = 50` と出ました。参照 SCF は ccsd.F が自分で 1e-6 まで締めるので、デッキで指定する必要はありません。未収束のときの流れもソースで確かめました。
  1. ccsd_iterdrv2.F が `****maximum iterations exceeded****` を出し、oconverged は .false. のまま返る。
  2. aoccsd2.F の `if (.not. oconverged) goto 999` で (T) とエネルギーの書き込みを飛ばす。
  3. ccsd.F が `ccsd = oconverged` を返し、task.F が `errquit('ccsd(t) energy failed', CALC_ERR)` で非 0 終了する。
  4. hfauto では NONZERO_EXIT（output.py:231-243）。TERMINAL の外なので保存されず（jobstore.py:155-157）、次の実行で計算し直される。誤った値を採用する経路はありません。v2 のロードマップの「`--retry-failed` で」という注意も実際には不要です。
- **手法誤差の数値（S23）は実測と一致します。** v3 の TS での CCSD(T)/def2-TZVPD は −205.4237916 Eh で、trans の −205.442352 Eh と合わせて ΔE‡ 11.65、SVPD − CCSD(T) = +2.02 です。
- **S22 は正しく動いています**（水と TMA·(HF)₂ の ranking.csv が `outcome:same_basin` だけ）。
- ○ に留める理由: U0-N1（層が欠けた subject の ΔE が LOT の混成になる）、U0-N8（順位差の目安が少し楽観的）、電荷のある系・開殻系と method_panel の実 run がないこと。

#### 4.0.3 有用性（◎、v2 は ○）

v2 で ○ とした前提（composite の手順の穴、甘い手法誤差）が、どちらも解消しました。

- **composite**: design.md:150 の手順どおりに、v3 HCN の run の複製へ sp_tz → thermo_tz → report_tz を追記すると動きました（`u7_assess/hcn_tz`）。SP 3 本（stage 83 s）、GoodVibes は 3 本とも JobStore のヒット、全体 1:26。TZVPD 層の値は ΔE‡ 46.43、ΔE_rxn 13.82、ΔG‡ 41.99 で rank 1 です。v2 の手順では `has no input of type ['minimum']` で止まっていました。
- **CCSD(T) の opt-in**: 重原子 4 個の HONO の TS でも余裕 32 反復で収束し、447 s wall（同時負荷あり）でした。method_panel.yaml:3-4 の目安（重原子 4 個程度）と見合っています。
- **主対象の報告行**: ΔG_assoc −11.53 の隣に thermo_unavailable が並ぶ矛盾がなくなりました（S22）。
- ◎ は条件付きです。validation.md:62・:96・:139 には C2 が「未確認」、S23 が「文書だけ」と書かれたままで、design.md:149 の数値の出典（プローブ）もリポジトリの文書にありません（U0-R3′ → S26）。WSL の本番経路では code_version が常に `-dirty` です（U0-N7 → C38）。method_panel そのものは v3 で実行していません。

#### 4.0.4 複雑さ（○、据え置き）

- コードの変化はごく小さいです（C2 は条件付きの 1 要素、S22 は分岐の整理で行数ほぼ同じ、S21 は 1 行、S23 は文書だけ）。knob も method ファイルも増えていません。
- trust 0.3 と CCSD の maxiter 50 がジョブキーに入らない点は、jobstore.py:45-52 の方針（ExecutionSpec と描画の入力はキーに入れない）どおりです。実際のデッキは `attempt_NN/job.nw` に残ります。
- ◎ にしない理由は、利用者が手で守る規則が多いことです。composite では「複製する、stage id を分ける、targets: all_minima にする、SP 層で ΔE‡ ≤ 0 は読まない」、パネルでは「conditions を元の run とそろえる」。conditions は 3 本の pipeline に重複して書かれ、パネルでの食い違いを止めるのは文書だけです。それぞれをコードで止めると規則や knob が増えるので（§4.0.8）、現状の ○ が釣り合いのとれた点です。

#### 4.0.5 残る課題

- **U0-N1（low）**: energy_method の層が欠けた subject では、ΔE が freq 層との混成で report.html に出る（thermochemistry.py:190-193）。順位からは blocker で外れる。→ C28
- **U0-N7（low）**: WSL の checkout で code_version が常に `-dirty`（v2 の 9 run、v3 の 7 run すべて）。→ C38
- **U0-R3 の未実施**: validation.md にパネルと composite のプローブの記録がない。→ S26

#### 4.0.6 新たな課題

- **U0-N8（low）「約 2 kcal/mol 未満の差は手法誤差」は、2 反応の順位を比べる幅としてはやや楽観的**: design.md:149 と README.md:86 が挙げる誤差（SVPD − CCSD(T)/def2-TZVPD）は HCN −1.19、HONO +2.02 で符号が逆です。2 反応の ΔG‡ の大小に効くのは誤差の差で、測った 2 系でも 3.2 kcal/mol あります。一方、ranking.csv の同順位の幅（qs × cutoff の感度）は HCN・HONO とも 0.001 kcal/mol 未満なので、順位の不確かさを読む手がかりは文書の目安だけです。見積もりは 2 系だけからのもので、実際に順位を付けるのは同じ run の中の反応（誤差の符号が相関しうる）なので、low にとどめます。今のベンチマークは 1 run 1 反応なので実害はありません。

#### 4.0.7 追加の改良案

- **U0-R3′（should）validation.md にプローブの小節を足す** 【S26】: (1) C2: `u0_c2` で実エンジンを呼び、デッキ `maxiter 50`、出力 `maxit = 50`、HNC 11/50 反復・17 s、v3 HONO の TS 18/50 反復・447 s（同時負荷あり）。CCSD(T)/def2-TZVPD の ΔE‡ 11.65、SVPD − CCSD(T) = +2.02（design.md:149 の出典）。(2) composite: `u7_assess/hcn_tz` で design.md:150 の手順どおり、ΔG‡ 41.99、全体 1:26。(3) パネル: v2 のプローブで HONO が rank 1 のまま。:62 の「C2 は未確認」を「v3 の run では未実行。プローブで確認」に変え、:96・:139 にも参照を付ける。v2・v3 の run の記録そのものは書き換えない。v2 のプローブのパスはリポジトリの外なので、場所と要点の数値だけを書く。文書 6〜8 行、中立。
- **C38（should）code_version の git status に `-c core.autocrlf=input`** 【C38】: config.py:162 を `_git(Path(tmp), git, "-c", "core.autocrlf=input", "status", "--porcelain", "--untracked-files=no")` にする。WSL で 83 → 0 行、Windows で 0 → 0 行を再測定済みで、副作用はない。引数 2 つ、中立。
- **C28（could）layer_missing の subject を含む ΔE を None にする** 【C28】: `_kcal`（190-193）で `sa.layer_missing or sb.layer_missing` なら None。test_thermo_sp_stages.py の S9 の節に assert を 1 つ。composite で SP が 1 本だけ失敗した場合の誤読（HCN で 164 kcal/mol）を防ぐ。1 行、中立。
- **U0-R7（could、修正して採用）順位差の目安を誤差の符号の違いまで含めて書く** 【S26】: design.md:149 を「誤差の符号は反応で異なる（HCN −1.2、HONO +2.0。2 系での見積もり）。2 反応の ΔG‡ の差が約 3 kcal/mol 未満なら、順序は手法誤差で入れ替わり得る」とし、README.md:86 も同じ表現にそろえる。数値は断定せず、コードの閾値にはしない。文書 2 か所、中立。

#### 4.0.8 見送った案

- **S21 の trust や CCSD の maxiter を JobStore のキーに入れる**: キーは化学的な要求に限る方針に反し、既存のキャッシュがすべて無効になる。行き着く極小は変わらない（v3 の TMA·(HF)₂ で ΔG_assoc の差 0.002）。デッキは `job.nw` に残り、CCSD の未収束は保存されないので自動で計算し直される。
- **パネルの conditions の一致をコードで検査する、T を既存の ReactionThermo から推定する**: 規則が増える一方、食い違ったときは目に見える形で失敗する。文書の注意書きで足りる。
- **手法誤差（約 2〜3 kcal/mol）を同順位の判定の幅に組み込む**: 小分子ではほとんどの反応が同順位になり、ΔG‡ の序数という読み方が失われる。見積もりは 2 系だけで、コードの閾値の根拠として薄い。
- **v3 の検証 run に method_panel を追加で実行する（U0 の目的で）**: パネルの DFT の経路は v2 のプローブで、C2 は v3 のプローブで、同じエンジンのコードを通した。記録を残せば足りる（U8 の観点での正式な実行は C43 として could）。
- **CCSD の maxiter や thresh を method ファイルの knob にする、MP2 にも maxiter を書く**: HONO の TS で 18/50 と余裕がある。参照 SCF は ccsd.F が自分で締める。MP2（semi-direct）は反復しない。
- **初期 Hessian なしの opt も trust 0.3 にそろえる**: 対角推定の Hessian では、mode-follow の変位から TS より上に戻る例がある（HCN→HNC）。v3 でも 0.1 のままで問題は起きていない。
- **composite.yaml を同梱し、blocker `barrier_vanishes_at_energy_layer` を足す（C1）**: S23 の手順が実際に動いた。同梱すると、元の run に追記したときに thermo の view が置き換わる危険が戻る。

#### 4.0.9 変更不要な点

- 3 層構成（GFN2-xTB は探索と事前判定だけ、停留点は PBE0-D3BJ/def2-SVPD の単一 PES、GoodVibes の qRRHO に自前の振動数）と、`PipelineConfig._check_stages` による読み込み時の method の一致の検査。
- 全 DFT デッキの `convergence energy 1e-7`、d3bj の `disp vdw 4`、ωB97X-D3 に disp 行を書かないこと。
- S21 の trust の使い分け、C2 の ccsd ブロックだけの `maxiter 50`、WFT は閉殻だけ（wft_closed_shell_only）。
- S22 の rankable、S23 の文書による composite の運用、パネルの既定を TZVPD の 2 手法にして CCSD(T) を 1 語の opt-in にすること、符号ゲートの不感帯 ±1 kcal/mol。
- ThermoSettings（grimme、100 cm⁻¹、vib_scale 1.0、symm、感度の幅）と、S9 の energy_layer_missing の fail-closed。
- JobStore のキーを化学的な要求に限る方針。

---

### 4.1 U1 入力構造・電子状態（structures）

評価: v2 ○ / ◎ / ◎ → **v3 ◎ / ◎ / ◎**（妥当性は U1 単体での評価）

#### 4.1.1 変更点

- **S24**: smiles.py:28-30 に 3 行（形式電荷の照合の後、AddHs の前で `GetIsotope()` が 0 でない原子があれば ValueError、`'specifies isotopes; hfauto uses natural-abundance masses'`）。docstring の 15 行目、test_smiles_input.py:17-18、design.md:138（「同位体を指定した SMILES は INPUT_INVALID（xyz の D も未知の元素として拒否）」）。
- **S25**: runner.py:152-157（正味 +4 行）。入力切れの例外に `'; upstream failures: <id>: <reason>, …'`（最大 3 件、重複なく整列）を付ける。テストは test_execute_stage.py の `test_missing_input_names_up_to_three_upstream_failures`。
- 変わっていないもの: core/system.py、stages/structures.py、electronic_state.py、topology.py、hypotheses.py、組成の規則、NWChem・xTB・CREST への電荷と多重度の受け渡し。
- 未実施: C31（U1-R3）、U1-R4（C39 の design.md:141）。

#### 4.1.2 妥当性（◎、v2 は ○）

v2 で妥当性を ○ に留めた理由のうち、U1 のコードの中にあったのは「SMILES の同位体を黙って ¹H にする経路」だけで、S24 で閉じました。v2 §3.1 の「S24 を入れれば U1 単体としては ◎ 相当」とも整合します。

- HEAD（prod venv、RDKit 2026.03.3）で、`[2H]O[2H]`、`[2H]O`、`[2H][2H]`、`[13CH4]`、`[1H]O[1H]` はすべて拒否され、`O`、`[OH-]`（q=−1）、`[H]O[H]`、`[CH4:1]`、`[H][H]` は通ります。
- 検査が漏れない理由: MolFromSmiles が内部で呼ぶ RemoveHs は、既定で同位体付きの H を残します（`RemoveHsParameters.removeIsotopes=False` を実機で確認）。そのため `[2H]` は AddHs の前の検査に必ずかかり、AddHs が付ける H（isotope 0）を誤って拒否することもありません。
- 既存の入力に退行はありません。v3 の実 run 7 本で、structures の config_sha と SpeciesRecord 18 件（composition_id、state_label、電荷、多重度、指紋）が v2 と一致し、structures/xyz の md5 も同じでした。HEAD の smiles_to_molecule は v2 の amine 系 run の xyz をバイト単位で再現しました。
- 受け渡し: v3 の NWChem デッキ 80 本はすべて `charge 0` / `mult 1`、odft 0 本。xtb の argv 26 本と crest 2 本もすべて `--chrg 0 --uhf 0` です。
- 検証範囲の注記: 電荷のある系・開殻系を structures → NWChem（odft）まで通した実 run はいまだ 0 本です。組成の多重度の高スピン固定（C8）は U2 の範囲です。状態ラベルの限界（U1-06 / N5）は文書にありません。

#### 4.1.3 有用性（◎、変わらず）

- v3 の 7 run で structures はすべて done、1 秒未満、外部ジョブ 0 本、18 種すべて success でした。
- S25 で、known_endpoints の典型的な失敗（structures の全種が INPUT_INVALID）の理由がトレースバックの最終行に出るようになりました（`s24_s25/run.log`）。
- 残るもの（どれも fail-closed で low）: U1-N1 の残り（Rich のトレースバック約 44 行が先に出る、終了コード 1、report.html が作られない、run_state に理由が残らない）、U1-N2（xyz の欠損が dry-run を素通りし、実行で structures が `running` のまま残る）、U1-N3・N4。

#### 4.1.4 複雑さ（◎、変わらず）

- 正味の追加は 7 行（smiles.py 35 → 38 行、runner.py 270 → 274 行）で、設定キーは 0 のままです。
- リポジトリの複雑さの関門は ruff の C901（mccabe、max-complexity 15。pyproject.toml:56, 61-62）で、execute_stage は v2 も HEAD も 4、smiles_to_molecule は 5 → 6 です。radon の値（execute_stage 10 → 14、smiles 6 → 8）は内包表記の for と if も数えた値で、関門ではありません（U1-N9）。
- 検査の置き場所の分担（記述だけで分かるものは読み込み時の pydantic、原子数・電子数が要るものは structures、SMILES の解釈に依存する同位体は smiles.py）は一貫しています。

#### 4.1.5 残る課題

- **U1-N1 の残り（low）**: トレースバックが先に出る。見た目の改善だけなので見送り（§4.1.8）。
- **U1-N2（low）**: xyz の欠損が dry-run を素通りし、実行すると hashing.py:11 で FileNotFoundError、run_state の structures が `running`・finished null のまま残る。再実行しても同じ（`u1_review` の run1.log と run2.log がバイト単位で同一）。→ C31
- **U1-N3・N4（low）**: 片側の端点だけ失敗したときの cascade と、両端の q/m の typo による無駄な DFT。→ C31(b)
- **U1-06 / N5（low）**: 状態ラベルの限界（1-WL、立体と E/Z を区別しない、中間帯は結合）が文書にない。→ C39′
- **U1-N7（low）**: hypotheses.py:158-159 の到達しない ValueError。削除せず防御コードとして残す。

#### 4.1.6 新たな課題

- **U1-N8（low）S24 のメッセージが「天然存在比の質量」と言うが、実際は主同位体の質量**: smiles.py:29-30 のメッセージは `hfauto uses natural-abundance masses` ですが、NWChem 7.2.3 の freq 出力の Atomic Mass は C 12.000000、N 14.003070、H 1.007825、GoodVibes 4.3.0 の io.py も H 1.00782503207、vibrations.py:18 のコメントは 'Most abundant isotope masses'、design.md:33 は「主同位体の質量」です。計算結果には影響しませんが、`[1H]O[1H]` の拒否理由が自己矛盾に見え、平均原子量で計算していると誤解させます。
- **U1-N9（low、レビュー入力の訂正）**: 仕様整理の「execute_stage の CC 10 → 14（C）、上限 15 の目前」は radon の値を ruff の関門と混同しています。ruff C901 では 4 のままで、切り出しは不要です。

#### 4.1.7 追加の改良案

- **U1-R5（could、修正して採用）S24 のメッセージを「主同位体の質量」に直す** 【C40】: smiles.py:30 を `hfauto uses most-abundant-isotope masses` にする（1 行）。test_smiles_input.py は `match='isotopes'` なので影響なし。validation.md:61 は実際に観測した出力の記録なので書き換えず、「のちに文言を修正」と注記するか、再プローブ後に更新する。中立。
- **U1-R3（could、v2 から継続）読み込み時に xyz の存在と宣言反応の両端の q/m の一致を検査する** 【C31】: (a) load_system（system.py:106-110）で欠けている xyz があれば FileNotFoundError。cli/main.py:121-126 の try が OSError を `_fail` で受けるので、run-dir を作らずに exit 2 で止まる。load_system を呼ぶのは pipeline/config.py:140 の 1 か所だけで、ほかのコマンドへの退行はない。(b) `_check_references` の反応ループで両端の charge・multiplicity の一致を検査。約 5 行とテスト 2 件、中立。
- **U1-R4（could、v2 から継続）状態ラベルの限界を design.md:141 に 1 文** 【C39′】: 「状態ラベルは結合グラフ（1-WL、立体と E/Z を区別しない、1.15〜1.45 Σr_cov の中間帯は結合）なので、配座・立体・E/Z だけ違う IRC の端は same_as_source として捨てる」。中立。

#### 4.1.8 見送った案

- **execute_stage の理由の収集を `_upstream_failures` に切り出す**: 関門の ruff C901 では 4 で、上限 15 まで大きな余裕がある。
- **S25 の理由を欠けている型の failed artifact だけに絞る**: 枠 3 件が無関係な理由で埋まる状況は実際にはほとんど起きない。
- **主同位体と同じ質量数の同位体指定（`[1H]`、`[12C]`）を許可する**: 元素ごとの表と分岐が要り、利用者が得るものはない（書かなければ同じ結果）。
- **CLI で入力切れの例外を捕まえ、トレースバックなしに 1 行で表示する**: 専用の例外型と handler が要り、ValueError を広く捕まえると本物のバグのトレースバックまで隠す。
- **失敗した stage の例外文を run_state に保存して `hfauto status` に表示する**: スキーマ変更になり、理由は manifest.json とログにすでに残っている。
- **hypotheses.py:158-159 の到達しない ValueError を削除する**: structures の S11b を経ない経路への 2 行の保険で、残す費用はほぼゼロ。
- **S25 の例外文に残りの件数を添える**: 全件は manifest.json にあり、典型は 1〜2 件。
- **runner の config_sha の計算を try の中に移す**: `running` が残る実際の引き金は xyz の欠損だけで、C31(a) により読み込み時に止まる。

#### 4.1.9 変更不要な点

- S24 の実装位置と順序（形式電荷の照合の後、AddHs の前）。
- S25 の形（consumes の型の failed artifact を重複なく整列し最大 3 件、失敗がなければ何も付けない）。
- xyz の `D` を未知の元素として拒否することと、SMILES の同位体の拒否が同じ INPUT_INVALID にそろっていること。
- M7（端点は xyz 必須）、S11b、S12、xyz のバイト単位の複製と 1e-6 Å の指紋、ETKDGv3 を 1 回・seed 20260925・力場なし。
- failed artifact として記録して stage を続ける方針、StructuresConfig に設定キーを持たないこと。

---

### 4.2 U2 配座探索・錯体配置（conformers）

評価: v2 ○ / ○ / ◎ → **v3 ○ / ○ / ◎（据え置き）**

#### 4.2.1 変更点

- U2 本体は変更なし（conformer_search.py 241 行、backends/crest.py 184 行、chemistry/placement.py 259 行、protocols.ConformerSettings は 573e80b と差分なし）。
- 間接（S25）: conformers の species が 0 件のとき、例外に上流の理由が付く（プローブ `u2_s24/run_only_d2o`）。
- 間接（S24）: 同位体の SMILES の species は structures で failed になり、それを含む組成だけが `conformers_<comp>` INPUT_INVALID `missing_component:d2o` になって、ほかの組成と単量体は続く（`u2_s24/run_mixed`: n_ok 3 / n_failed 1）。
- v2 の改良案は 4 件とも未実施（U2-R1/C34、U2-R2、U2-R3/C30、U2-R4/C8）。

#### 4.2.2 妥当性（○、据え置き）

- **CREST 3.0.2 の挙動を公式ソース（v3.0.2 タグの confparse.f90）で再確認しました。** `-nci` は runver=4 と NCI を設定し、後の `-quick` が runver=2・ewin=5 に上書きし、さらに後の `-ewin 6` が 6 に戻します。v3 の組成の stdout にも `sorting energy window (EWIN): 6.0000` が出ています。`-noopt` は preopt=.false. を設定するだけで、`-noreftopo` は reftopo=.false. を設定するだけです。v2 の結論（M4 の `--noopt` は正しい、U2-N1 の推測）と矛盾しません。
- **実 run（v3 tma_hf2_discover）**: 組成の argv は v2 と同一（`--gfn2 --nci --quick -T 4 --ewin 6 --chrg 0 --uhf 0 --notopo 2,14,15,16,17 --noopt`）、rc 0、12.10 s（v2 12.05 s）、3 配座。単量体 TMA は 6.25 s（v2 6.39 s）、1 配座。diagnostics.json は v2 と一致しました。エネルギーの差は ≤ 8e-6 Eh（約 0.005 kcal/mol）で、CREST の実行ごとの揺らぎの範囲です。下流の DFT 極小は W7・v2・v3 で同じ 1 つです。
- **S24・S25 の境界**はプローブで仕様どおりに動き、U2 の側に誤りは入っていません。
- ◎ にしない理由（v2 と同じ）: U2-N1（PT 型のトポロジー停止のあとの `--noreftopo` 再実行が実 CREST で 2/2 失敗）、U2-N2（開殻の組成は `--noopt` で trial MTD が失敗）、C8（組成の多重度を黙って高スピンにする）。v3 の実 run はどれもこの経路を通っていません。

#### 4.2.3 有用性（○、据え置き）

- 組成に対する CREST `--nci --noopt`、placement の seed00、`--notopo` の自動付与、重複除去、版の pin、BUG-08 の規則は v3 でも安定して働きました。conformers stage は 19 s で、DFT の freq（約 15 分）と比べれば無視できます。
- 下流への寄与は v2 と同じく限定的です。TMA·(HF)₂ では CREST の 3 配座と入力の構造が xTB で 1 つの basin にまとまりました。conformers が新しい DFT 極小を生んだ実例はまだありません。
- placement の seed の生座標は run ごとに最大 6〜10 Å ずれますが、permutation-invariant RMSD では ≤ 0.0003 Å で、TMA の CREST 出力の向きの違いだけによるものです。

#### 4.2.4 複雑さ（◎、但し書きは v2 と同じ）

- radon は v2 と同じ（_composition C(14)、ConformersStage C(13)、_select C(12)、ConformersStage.run C(12)。3 ファイルの平均 A(3.76)）。利用者が書ける knob は 4 のまま。
- 残る複雑さは v2 と同じ（`--noreftopo` の分岐、key_payload の冗長な noopt、重複の基準が 2 系統、範囲検査がない）。第2ラウンドで U2 に増えた複雑さはありません。

#### 4.2.5 残る課題

- **U2-N1（low）**: 単量体の停止後の `--noreftopo` 再実行が PT 型の停止で失敗する。二重失敗のとき `_lowest` は未緩和の入力を優先する。→ C34
- **U2-N2（low）**: 開殻の組成は `--noopt` で失敗しうる（文書に注記なし）。→ C39′
- **U2-N3（low）**: environment.md:64 の「約 25 s」（実測約 1.4 s）。→ C39′
- **U2-N4（low）**: ConformersConfig の範囲検査がない（keep_per_state=0 で組成が黙って 0 件）。→ C30
- **U2-N5（low）**: site で crest の threads を省くと `-T 1`。文書の注記だけ残る。→ C39′
- **U2-ISS-3（partially）**: 組成の多重度の高スピン固定（文書は design.md:138 の structures 行とモジュールの docstring:7）。→ C8

#### 4.2.6 新たな課題

- **U2-N6（low）モジュールの docstring が knob の数で実装と食い違う**: conformer_search.py:9-10 は「Users set only `quick` and `ewin_kcal`」ですが、ConformersConfig と design.md:190 は 4 knob です（C7 の 116ee0e から変わっていない）。挙動には影響しないので、U2-R2 の docstring 修正の 1 行にまとめます。
- **U2-N7（low、実質 info）組成の missing_component の reject が上流の原因を運ばない**: `u2_s24/run_mixed` で `conformers_h2o_d2o` は reason `missing_component:d2o`、parents は []。根本原因（同位体）は structures の manifest を見ないと分からず、diagnostics.json にも出ません。ただし reason に component の id があるので 1 段で辿れ、挙動の誤りではありません。

#### 4.2.7 追加の改良案

- **U2-R1（could、v2 から継続）単量体の停止後の再実行を、停止構造から既定の設定で 1 回に替え、`--noreftopo` を廃止する** 【C34】: v2 §3.2.6 のとおり。v2 のプローブで PT 型の停止が 2/2 成功（グリシン 5 配座、β-アラニン 7 配座）。フィールドと CLI フラグが 1 つずつ減り、複雑さの ◎ の但し書きが外れる。約 10 行と fake テスト 1 件、【複雑さ減】。
- **U2-R2（could、v2 から継続、U2-N6 を含める）文書と docstring を実測に合わせる** 【C39′】: environment.md:64 の「約 25 s」を「約 1.4 s（`--noopt` なしでは初期トポロジー検査で止まり約 25 s）」にし、crest の threads の既定が 1 であることを添える。design.md:139 に「開殻の組成は `--noopt` で trial MTD が失敗することがあり、そのときは placement の seed が出る」を 1 文。conformer_search.py:9-10 を「Users set quick, ewin_kcal, seeds_per_composition and keep_per_state」に。中立。
- **U2-R3（could、v2 から継続）ConformersConfig に範囲の制約** 【C30】: keep_per_state・seeds_per_composition に `Field(ge=1)`、ewin_kcal に `Field(gt=0)`。3 行とテスト 1 件、中立。
- **U2-R4（could、v2 から継続）組成の多重度が一意でなければ INPUT_INVALID** 【C8】: conformer_search.py:202 で coupled_multiplicities が 2 つ以上の値を返したら `ambiguous_multiplicity:<候補>`。2〜3 行、中立。

#### 4.2.8 見送った案

- **missing_component の reject に上流の理由を連結する、parents に失敗した部品を入れる（U2-N7）**: reason で 1 段で辿れ、stage 全体が空なら S25 が理由を出す。上流理由の収集を stage 側にも複製すると経路が 2 つになる。
- **diagnostics.json に reject された組成も書く**: reject は failed artifact として manifest に残っている。
- **placement の seed の生座標を決定的にする**: 実質同じ構造で、キャッシュにも影響しない。
- **v3 で (HF)₃ の配座（#9）を再実行する**: コードも argv も v2 と同一で、tma_hf2 で揺らぎの範囲で再現した。
- **全原子を除外する小さな酸塩基錯体で、`--notopo` の代わりにトポロジー検査を丸ごと切る**: 実効は同じで、分岐が増えるだけ。
- **v2 で見送った案（key_payload の noopt の削除、重複判定の統一、開殻の組成だけ preopt に戻す、6 seed すべてから CREST、`-T 1` による決定化、crest の threads の既定を 4 に）**: v3 でも前提は変わらない。

#### 4.2.9 変更不要な点

- 組成に必ず `--noopt` を付けることと argv の順。seed00 の labile H と受容原子への `--notopo` の自動指定。
- エネルギー窓を置かず、状態ラベルごとに keep_per_state 個を残すこと（S13）と BUG-08 の None の規則。
- 利用者が書ける knob を 4 つにし、スレッド数を site の 1 か所から決めること（C7）。
- placement の決定的な rng_seed、CREST に渡すのは seed00 だけで、失敗時に全 seed を出す fallback。小さい剛体の単量体を skip すること。
- 版の pin（METHOD_MISMATCH）、key_payload の noopt、diagnostics.json。
- S24・S25 との境界の扱い（弾かれた部品を含む組成だけを止め、stage 全体が空なら runner が理由付きで止める）。

---

### 4.3 U3 極小構造の確定（minima screen / dft、MinimumDriver、Registry）

評価: v2 ◎ / ◎ / ○ → **v3 ◎ / ◎ / ○（据え置き）**

#### 4.3.1 変更点

- U3 本体（minima.py 235 行、minimum.py 277 行、identity.py 179 行、selection.py 65 行、gates の is_minimum・imaginary_tier）は 1 行も変わっていません。
- **S21**: render_optimize（input.py:116-126）は、初期 Hessian を渡すとき `trust 0.3` + `inhess 2`、渡さないとき `trust 0.1`。U3 では 2 断片以上の species（xTB の init_hessian）だけが対象で、1 断片、mode-follow、soft 押しの opt は trust 0.1 のまま。JobStore のキーは変わらない（a80a4e793d96 は v2 と v3 で同じ）。
- **S25**: minima stage に SPECIES が来ないとき、上流の理由を添えて止まる。
- 文書: design.md:142 の minima(dft) の行を「初期 Hessian なしは trust 0.1、ありは trust 0.3」に。U3-N4（include: all の文言）は未修正。
- 未実施: U3-R1（C32）、U3-R2（C10）、U3-R3（C33）。

#### 4.3.2 妥当性（◎、据え置き）

- **S21 の実装は NWChem 7.2.3 の仕様どおりです。** 公式文書と v7.2.3-release の src/driver/opt_drv.F を照合しました。trust の既定は最小化 0.3、鞍点探索 0.1（opt_drv.F:1069-1074）。trust は 2 か所で歩幅を抑えます（各モードの成分を |dv| ≤ trust に切る :1743-1787、歩み全体を alpha = trust/dsmax で縮める :1905-1931）。line search は制限された歩みを最大 4 倍まで伸ばします（:1874）。
- v3 の全 run の opt デッキは、Hessian ありの 13 本がすべて `trust 0.3` と `inhess 2`、なしの 12 本がすべて `trust 0.1` でした。
- **結果は変わっていません。** A arm は 4 錯体で 27 → 27 点、最終エネルギーの差は 3e-8 Eh 以内、basin は 3 つのまま、tma_hf2_afir は v3 でも known（diagnostics.json は v2 と同一）。TMA·(HF)₂ の最低 ν は 51.836 → 51.827 cm⁻¹。1 断片の系（HCN、HONO、NH3、水、STO-3G）はデッキ・キー・ステップ数が同じです。
- **S25 と S24 を実ログで確認しました**（`s24_s25/run.log`、`run_xyz.log`）。
- **留保は v2 と同じです。** DFT の mode-follow（−側、ts_candidate、両側の再 settle、discovery）、DFT の soft 押し、M6 の |ν| 化は実計算で通っていません。v3 の dft の diagnostics は全 run で `freq:none` か `known` でした（soft が発火したのは xTB の screen の `neutral` だけで、v2 と同じ）。

#### 4.3.3 有用性（◎、据え置き）

- S21 の U3 への効果は小さいものの、観測されました。17 原子の錯体の opt 1 本あたり約 20 s（−5%）短くなり、A arm の opt 合計は 899 → 854 s です。同じ軌跡をたどった run どうしの差（約 2%）より大きい差です。
- 時間の大半はまだ DFT の解析 freq です（tma_hf2_discover の dft stage 1,754 s のうち 898 s）。同じ PES と numerics を保つ限り、削る余地はありません。
- 錯体の opt では、時間の約 3 割が XMAX / XRMS の判定だけのための歩みです（v3 の TMA·(HF)₂ の 5〜8 歩目で約 110 s、エネルギー変化 3.2e-6 Eh）。これを削る案は見送ります（§4.3.8）。
- S1、M3、init_hessian の効果は v2 と同じです（screen の xTB freq の JobStore hit 1、dft は hits 1 / misses 9、known の join も同じ）。

#### 4.3.4 複雑さ（○、据え置き）

- S21 は input.py:124 の 1 式の置き換えだけで、テストは test_nwchem_input.py と test_nwchem_engine.py で trust 0.3 と 0.1 の両方を確かめています。
- 同じ版のコードの初回の attempt では、trust はキーに入る `hessian` の sha が None かどうかだけで決まります（engine.py:342、`_hessian_file` が失敗すれば opt に進まない）。デッキとキーの対応が崩れるのは、版をまたいで再開したとき（v2 のキャッシュを v3 のコードで再利用）と、ladder の継続 attempt（U3-N8）です。前者は JobStore 全体の性質です。
- v2 で挙げた複雑さ（C10 の暗黙の分岐、rerank の always の非対称、`_known` の version の自己比較、known.find と add の絞り込みの非対称）はそのまま残っています。

#### 4.3.5 残る課題

- **U3-ISS-5（low）**: window で screen 極小 0 件なら黙って空の dft stage になる（S25 の対象外）。→ C10
- **U3-N1（low）**: rerank_sp が per_state 以下の組にも SP を出し、失敗すると黙って脱落する（TMA で SP 3 本・49.8 s）。→ C32
- **U3-N2（low）**: reaction-paths の `_register` が pipeline の gates を渡さない。→ C33
- **U3-N4（low）**: design.md:142 の include: all の文言。→ C10
- U3-N3・N5・N6 は v2 の判断（発火してから入れる、設計上の性質）を維持します。

#### 4.3.6 新たな課題

- **U3-N7（low）validation.md:59 の S21 の説明が実ジョブの記録と食い違う**: validation.md:59 は S21 の効果が小さい理由を「trust 0.1 でも歩幅が上限に届いていなかったため」と書いています。しかし v2 の 4 錯体の opt の stdout では `Restricting large step in mode` が 2 / 1 / 17 / 14 回（計 34）、`Restricting overall step due to large component` が計 7 回（a80a4e793d96 の alpha は 0.76 / 0.90 / 0.81）で、v3 は計 3 回と 2 回でした。歩数が 8 → 8 で変わらないのは、(1) v2 では line search が制限後の歩みを 1.28 / 1.67 / 1.19 倍に伸ばして補っていた、(2) v3 の 4〜8 歩目は gmax・grms が閾値未満になった後に XMAX / XRMS だけを満たすための歩み、の 2 つによります。
- **U3-N8（low）初期 Hessian ありの opt を ladder で継続すると trust 0.1（inhess なし）に戻る**: continuation（engine.py:201-203）は継続で inputs の `hessian` を None にし、`_render`（:104-105）経由で input.py:124 が `trust 0.1` を選びます。drv.hess は INHESS 0 の restart で引き継がれるので、trust 0.1 の根拠（対角推定の Hessian）は当てはまりません。docstring にも design.md §8 の ladder の行にも書かれていません。実 run では 0 回（U3 の錯体の opt は最大 8 歩、maxiter 100）で、QRC に対しては安全側に働きます。
- **U3-N9（low、レビュー入力の訂正）**: 仕様整理の「−5% は run ごとの揺らぎ（同じデッキで ±10%）の範囲」は正しくありません。同じ軌跡の run どうしの walltime の差は約 2% 以内（v2 395.4 / 398.8 / 403.2 s、v3 376.6 / 376.7 / 376.9 s）で、対照に使った tma（cfb9b780e34c）の +10.6% は 3〜4 歩目で軌跡が分かれて上り坂の歩みが入ったためです。また「デッキとキーが 1 対 1 でない」は、同じ版の初回の attempt には当たりません（§4.3.4）。

#### 4.3.7 追加の改良案

- **U3-R4（should、修正して採用）validation.md の S21 の説明を実ジョブの記録に合わせる** 【S26】: validation.md:59 を「v2（trust 0.1）では、4 錯体の opt の 1〜5 歩目で NWChem が歩幅を制限していた（`Restricting large step` 34 回、全体の縮小 7 回。v3 は 3 回と 2 回）。歩数が同じなのは、line search が制限された歩みを伸ばして補い、TMA·(HF)₂ の 4〜8 歩目が勾配の収束後に XMAX / XRMS を満たすための小さな歩みだからである。5% の短縮は、同じ軌跡の run どうしの差（約 2%）より大きい」にする。n = 2〜3 の比較なので観測の表現にとどめ、原因（後半の SCF が軽い）は推定と明記する。:124 は手を入れなくてよい。文書 2〜3 行、中立。
- **U3-R5（could）ladder の継続で trust 0.1 に戻ることを明記する** 【C39′】: design.md §8 の ladder の行（「NWChem は最新構造・movecs・drv.hess」）の直後に「初期 Hessian ありの opt の継続は、drv.hess を INHESS 0 の restart で引き継ぎ、trust 0.1 で続ける」を 1 文。コードは変えない。中立。
- **U3-R1（could、v2 から継続）rerank_sp の SP を per_state を超える組だけに出す** 【C32】: v2 §3.3.6 のとおり（4〜5 行、テスト 1 件）。
- **U3-R2（could、v2 から継続）window で screen 極小 0 件なら拒否し、include: all の文言を直す** 【C10】: v2 §3.3.6 のとおり。S25 ではこの経路を塞げないことを確認した。
- **U3-R3（could、v2 から継続）reaction-paths の `_register` に pipeline の gates を渡す** 【C33】: 1 行と import。

#### 4.3.8 見送った案

- **trust（やデッキテンプレートの版）を opt の JobStore のキーに入れる**: 同じ版の初回の attempt では trust はキーの hessian の sha から決まるので重複するだけ。版をまたぐ再開では同じ極小に着き害がない。入れると全 opt のキーが変わり v3 のキャッシュを再利用できない。
- **ladder の継続でも trust 0.3 を保つ**: U3 では 1 度も起きておらず、QRC の継続では trust 0.1 の方が安全側。inputs に新しいキーが増えるので、文書化（U3-R5）で足りる。
- **錯体の opt の収束判定を緩める（gmax が収束したら XMAX / XRMS を待たない）**: 約 110 s 減る見込みだが、この歩みは最低 ν 約 52 cm⁻¹ の柔らかい分子間モードに沿った動きで、freq の tier 判定と同一性（0.05 Å）は収束した構造を前提にしている。
- **1 断片の species の opt も trust 0.3 にする**: 1 断片の opt は 1〜5 歩で、制限による損失は見えない。
- **mode-follow と soft 押しの opt に DFT freq の Hessian を渡し trust 0.3 にする**: v3 でも DFT では発火 0 回で、効果の根拠がない。

#### 4.3.9 変更不要な点

- opt と別ジョブの解析 freq を同じ numerics で計算し、is_minimum と虚振動の 3 段方針（noise 10 / saddle 50 cm⁻¹）で判定すること。質量加重の Eckart 射影。
- S21 の範囲（2 断片以上の xTB 初期 Hessian の opt だけ trust 0.3 と inhess 2）。2 断片以上のときだけ xTB の Hessian を init_hessian にすること。
- S1 の直列の settle とその場の登録、同一性の基準を assign の 1 つにしたこと。M3、window 6 kcal/mol・per_state 3・rerank_top 8。
- soft 押し +0.1 Å を 1 回、±0.1 Å の mode-follow を最大 2 サイクル。xTB の `--opt vtight` と `--hess`。
- S25 の停止メッセージと S24 の上流での拒否。OPT_MAXITER 100、HESSIAN_NEAR_A 0.5 Å、opt のキー {molecule, method, hessian sha}。

---

### 4.4 U4 反応探索（explore）

評価: v2 ○ / ○ / ○ → **v3 ○ / ○ / ○（据え置き、有用性は下限）**

#### 4.4.1 変更点

- explore のコードと設定は変わっていません（explore.py、readuct の engine.py・worker.py、trials.py、topology.py、records.py、discover.yaml の差分は空）。
- 間接的に触れた変更は 3 つだけ: runner.py:152-157（S25。v3 では発火なし）、vibrations.py の `cartesian_mode_number` の追加（worker が使う `projected_frequencies` は不変）、protocols.py の SaddleRefiner のコメント 1 行。
- 未実施: C12、C13、C14′、C39 の U4 分（validation.md:91、design.md:141）。

#### 4.4.2 妥当性（○、据え置き）

- **v2 と v3 の一致**: v2・v3 の tma_hf2_discover で job の sha256 の集合が 31/31 で一致し、outcome と reason も 31/31 で同じ、ts_imag の差は 0.0 cm⁻¹ でした。source の opt は 31/31 件が 2 反復・ΔE 0.000 kJ/mol（M3 は保たれている）。1000 K のやり直し、timeout、failed はいずれも 0 件、採用した生成物 0 件、状態ラベルの反転も 0 件でした。
- **TS がある attempt の生データ**（プローブ `u4_review`）: e30f4cd2（−289.46i）は IRC の両端を opt するとどちらも source と同じラベル・ΔE 0.00 kJ/mol で、縮退した HF 交換を正しく same_as_source と判定しています。ff93257c（−588.98i）は両端とも +314.13 kJ/mol の別ラベルの構造で、irc_not_connected_to_source は正しい。イリドの 8c13ba82 は TS +376.37、product +313.89 kJ/mol で窓外は正しい。
- **NT2 の陰性 2 件**（U4-N4）はどちらも Bofill TSOpt が source の basin へ滑り落ちたもので、陰性の判定は化学的に正しい。0fd43cfb は +2.6 kJ/mol の置換コピーに収束し ν1 −24.5 cm⁻¹（cutoff 50 未満なので TS に数えないのは正しい）、7f7d4c91 は +273 から +10.5 kJ/mol まで下って `NoNegativeEigenValueException` で止まりました（プローブでビット単位で再現）。
- **ライブラリの確認**: ReaDuct 6.1.0 の TsOptimizationTask.h は `catch(...)` で「TS Optimization failed with error!」を出し、stop_on_error なら再送出します（原因の diagnostic を出すのは stop_on_error=False のときだけ）。反復上限は `settings.getInt("convergence_max_iterations")` で、C13 の引数名は正しい。
- 残る制約（v2 と同じ）: 等価な drive が上限を埋める（U4-I4）、F→N 移動の drive は作らない（U4-I8）、状態ラベルの閾値近くの H 結合（U4-N2）。

#### 4.4.3 有用性（○、下限のまま）

- explore stage は 81.0 → 81.0 s、readuct の合計は 79.9 → 79.7 s。下流の DFT に偽の生成物を流すことはありません。
- 限界（v2 と同じ）: 対象系の生成物は 0 件。31 attempt のうち 21 件（68%）が同じ C→N 1,2-H 移動に使われた。TS がある negative は ΔE‡ を持たない（縮退 HF 交換 44.05 kJ/mol、−589i の TS 378.13 kJ/mol は生データからしか得られない）。切り捨てた trial の件数も出ない。
- job key が v2 と v3 で完全に一致したので、同じ run dir で再開すれば explore は全件キャッシュが当たります（再現性の面では有用）。
- NT2 と TSOpt の失敗原因の文字列は result.json にも stderr.txt にも残らず、v2・v3 のレビューで 2 回続けてプローブの再実行が要りました。

#### 4.4.4 複雑さ（○、据え置き）

- radon の SLOC は explore.py 124、worker.py 189、engine.py 93、trials.py 180（計 586）で v2 から増減なし。CC が C の関数は ruff の上限 15 未満。knob も増減なし。
- 残る複雑さは v2 と同じ（ReaDuct の反復上限が暗黙の 150、陰性理由は約 10 種あるが件数にしか使わない、AFIR の reason が γ ごとに上書き）。新しく見えた点として、`_Scine.task` が RuntimeError の文字列を捨てます（U4-N4）。

#### 4.4.5 残る課題

- **U4-I4（low）**: 等価な drive で上限が埋まる。→ C12
- **U4-I5（low）**: ReaDuct の反復上限 150（v3 の TMA 系では無バイアスの opt は最大 109 反復で実害なし）。→ C13
- **U4-I6（low）**: TS がある negative に ΔE‡ がない。→ C14′
- **U4-N1（low）**: validation.md:91 は「W7 で irc_not_connected_to_source だった −289i の TS は、グラフの一致で same_as_source と判定された」のまま（実際は M3 の効果で、縮退した HF 交換）。→ C39′
- **U4-N2（low）**: 状態ラベルの閾値近くの H 結合（受け入れたリスク。v3 のラベル付き端点 42 件で反転 0 件）。
- **U4-N3（low）**: AFIR の reason が γ ごとに上書きされる（318c4965 と 8aed609c は γ125 で same_as_source だったが、最終 reason は afir_not_converged）。

#### 4.4.6 新たな課題

- **U4-N4（low）ts_not_converged が最適化器の例外（負の固有値の喪失）もまとめて表し、原因の文字列がどこにも残らない**: 7f7d4c91 の stderr.txt には「TS Optimization failed with error!」の 1 行しかありません。同じ guess から stop_on_error=False で再実行すると `Scine::Utils::NoNegativeEigenValueException` が出ました（`u4_review/tsopt_7f7d4c91.log`）。worker の `_Scine.task`（worker.py:137-140）は RuntimeError の文字列を捨てて False を返します。陰性の判定は正しく、影響は診断だけです。
- **U4-N5（low、数値の訂正）**: 縮退 HF 交換の ΔE‡ は 44.05 kJ/mol（v2 レビューの 41.5 は誤り。TS −24.392938 Eh と source −24.409714 Eh の差）。v3 の仕様整理の「端点の opt は最大 133」の 133 は AFIR γ300 のバイアス付き最適化の反復数で、無バイアスの opt の最大は 109 です（結論「150 に達した例はない」は変わらない）。

#### 4.4.7 追加の改良案

- **C12（could、v2 から継続）drive の重複を WL の原子クラスで除き、列挙数・固有数・実行数を diagnostics.json に出す** 【C12】: v2 §3.4.6 の修正後の形のまま。スキーマと JobStore の key は変えない。TMA は 2 クラスから 5 クラス、TMA·(HF)₂ は 5 クラスから 10 クラスまで被覆が広がり、有用性の ○ が下限でなくなる条件の 1 つ。15〜25 行、テスト 2 件、中立。
- **C14′（could、v2 から継続）TS がある negative にも ΔE‡ を manifest に載せる** 【C14′】: v2 の修正後の形のまま（SCC のやり直し温度では 300 K で評価し直すか載せない。reaction_kj_mol は product のときだけ）。v2 の数値 41.5 は 44.05 に直す。4〜8 行、中立。
- **U4-R4（could、新規）`_Scine.task` で捨てている例外の文字列を stderr に 1 行書く** 【C41】: worker.py:137-140 の `except RuntimeError as exc:` で、SCC の判定の後に `print(f"{run}: {exc}", file=sys.stderr)`。result.json のスキーマ、陰性理由、job key、knob は変えない。NT2 の no-guess と TSOpt の NoNegativeEigenValueException を、再実行なしに生ファイルから区別できる。C14′ と同じ変更で入れる。1〜2 行、中立。
- **C39 の U4 分（could、v2 から継続）** 【C39′】: validation.md:91 を「M3（screen の最適化構造から出発）の効果。IRC の両端は source の置換コピーで、縮退した HF 交換（ΔE‡ 44.05 kJ/mol）」に。design.md:141 の注記は U1-R4 と同じ 1 文。
- **C13（could、v2 から継続）無バイアスの opt に `convergence_max_iterations=500`** 【C13】: `_Scine.minimum` の run_opt_task だけ。v3 の TMA 系では効かないが、v2 の amine run には数値的な陰性の実例がある。1 行、中立。

#### 4.4.8 見送った案

- **ts_not_converged を mode_lost と max_iterations などに分ける**: 陰性ラベルは件数にしか使わず判定を変えない。原因を知る必要は U4-R4 の stderr で満たせる。
- **n_imag ≠ 1 の TS 候補でも ν1 を記録する**: 判定にも報告にも使わない。
- **Bofill が負の固有値を失ったとき別の optimizer やモードで TSOpt をやり直す**: v3 の 2 件はどちらも source の basin へ戻っており、失われた TS の形跡はない。同じ drive の AFIR も走っている。
- **amine パイロットを v3 で再実行する**: コード・設定・worker が同じで、job key が 31/31 で一致した。
- **IRC の反復上限を上げる**: 端点は必ず opt し直し、18〜109 反復で収束している。
- **S21（trust 0.3）を explore の ReaDuct の opt に合わせる**: trust 0.3 は NWChem の DRIVER の設定で、ReaDuct の最適化器とは別物。source の opt は 2 反復で律速ではない。
- **v2 §3.4.7 の見送り案の見直し**: 第2ラウンドは explore に触れておらず、v3 の実データも v2 と同じ。

#### 4.4.9 変更不要な点

- M3（screen の `--opt vtight` の最終構造から出発）と worker の source の opt（PES の一致の確認にもなる）。
- S6（状態ラベルでの一致判定）、M1（陰性結果は DISCOVERY artifact として件数にだけ使う）。
- NT2 → 同じ drive での AFIR へのフォールバック、NT2 の TS 検証ゲート（射影振動数で ν < −50 cm⁻¹ がちょうど 1 本、IRC、両端の n_imag=0）。
- 窓の fail-closed（ΔE_rxn が None なら窓外）、ΔE‡ 150 / ΔE_rxn 100 kJ/mol。
- job key の設計（timeout_s を除いた settings）、scine-readuct 6.1.0 の版の pin、SPIN_MODE の固定、IRC の stop_on_error=False と端点の opt、SCC 失敗時の 1000 K での 1 回のやり直し、relaxation_discoveries。

---

### 4.5 U5 反応経路：仮説選択と障壁の事前判定（hypotheses.select / SCREEN / FIND_PATH）

評価: v2 ○（SCREEN ○、FIND_PATH ○−）/ ○ / ○ → **v3 ○（SCREEN ○、FIND_PATH ○ の下限）/ ○ / ○**

#### 4.5.1 変更点

- **M8 + S19**（actions.py:459-482）: 成功したチャンクごとに、string_final の両端を ctx.ends に戻し、画像を `align_sequential` で逐次整列して `string{n}_c{i}.xyz` に書き、次のチャンクの initial_path と `ctx.work.path` に使う（以前は `initial = run.images` をそのまま渡していた）。note `string{n}:c{i}:end_drift:<Å>` を追加（ゲートにはしない）。bead のエネルギーは NWChem の値のまま。
- `_initial_frames` → `_initial_path`（actions.py:445-456）: FileRef を返すようにした。初期経路の中身（内部の画像を ends[0] へ 1 枚ずつ独立に整列し、両端を DFT 極小にする）は v2 と同じ。
- **S18**（interpolation.py:17）: IDPP の反復上限 500 → 5,000。
- テスト: DriftingPath（test_reaction_case_actions.py:40-54。initial_path を記録し、毎チャンク最後の bead をずらして画像 3 を回転させる scripted engine）と、それを使う `test_string_chunks_restore_the_frozen_ends_and_align_the_images`、`test_idpp_leaves_the_plane_between_planar_cis_trans_hono`。
- 文書: design.md:130 に M8 と end_drift の note。validation.md §8-1 の旧 GS 対策案（Cartesian、再試行）を取り下げた。
- 変わっていないもの: 仮説選択、決定表（state.py）、SCREEN と barrier_verdict、profile.py、pysis の GS と S4、NWChem の string のデッキ。
- 未実施: C15、C16、C35、C36、U5-N5、U5-N7。

#### 4.5.2 妥当性（○、FIND_PATH は ○− → ○ の下限）

- **M8 + S19 は仕様どおりです。** find_path_sto3g（v3）で、3 チャンクの報告端点エネルギーはどれも HCN −92.10839487、HNC −92.06836342 でした。各チャンクの初期経路と string_final の写像 RMSD は bead 1 が 0.0000、bead N が 0.0737 / 0.0008 / 0.0024 Å で note の値と合います。v2 の 2 チャンク目の +8.85 kcal/mol の非極小は消え、2 チャンク以上の monotonic → BARRIERLESS の前提が戻りました。
- **S18 も効いています。** HONO の FIND_PATH の初期経路はねじれ経路になり（二面角 180/160/141/119/96/80/48/24/0°、H–O–N 103〜108°）、string の最大は trans から 13.47 kcal/mol（v2 は 34.37）、HEI の xTB Hessian は虚振動 1 本（−744.4 cm⁻¹）で採用されました。
- **新しい U5-N8 を裏付けました。** 原因は `_initial_path`（actions.py:455）が内部の画像を ends[0] へ 1 枚ずつ独立に整列することです。
  - IDPP の出力はもともと ends[0] の座標系で連続しています（隣との距離は HONO 0.30 Å、HCN 0.475 Å）。
  - ところが ends[0] への共分散の階数が落ちている、または条件が悪い（HCN では特異値が [s, 0, 0]、HONO では第 2 特異値が 0.42 → 0.04）ので、独立に整列すると各画像が 18〜73°（HONO）、112〜176°（HCN）回転し、最後の区間に剛体の跳びが入ります（隣との距離は生の値と整列後で HONO 1.167 / 0.170 Å、STO-3G 0.94 / 0.475 Å）。
  - NWChem 7.2.3 の string.F:569-571, 1967, 2007 は、zts_min_motion（各軸の並進と回転を距離和の線探索で合わせる）を i = 2..nbeads に適用し、凍結 bead N も対象にします。実 run では bead N が剛体として 1.25 Å（HONO）、1.62 Å（HCN）回り、その後の L-BFGS がこの剛体移動を線形に外挿して構造をゆがめました。
  - ゆがんだ端点の実エネルギー（プローブ `u5_assess` で DFT SP 8 本）は、真の極小より HONO（PBE0-D3BJ/def2-SVPD）で +1.09、STO-3G HCN で +13.57 kcal/mol 高く、NWChem が報告するのは nstep 0 の極小の古い値のままでした（string.F:1830-1831）。同じ初期経路を逐次整列しただけのプローブ（`u5_spec/seqinit*`）では +0.009 と +0.041 kcal/mol でした。
  - 報告されるエネルギーは極小の値なのでプロファイルの数値そのものは正しいのですが、N−1 番目の bead はゆがんだ端点に向かって緩和されるので、行 14 の根拠（連続経路の最大は鞍点の上界）が端点の近く、つまり障壁が小さい領域で崩れます。観測した形の判定は変わっておらず（HONO の最大 13.47 と 13.46）、誤分類の実例はありません。
- **impose は効いていません（U5-N9）。** xyz_path を渡すと string_input.F（v7.2.3-release:148-173）が `bead_list:new=.false.` にし、impose の処理は string.F:207-248 の newchain 分岐の中にしかありません。実 run の stdout にも 'Imposing' は一度も出ず、bead は 'loading geom' で 9 本読み込まれていました。
- SCREEN・仮説選択・barrier_verdict・profile は変わっておらず、妥当なままです。HCN・NH3 の SCREEN は v2 と同じ proceed（max_rel_dft 33.98 / 4.465）。HONO の GS は同じキーで決定的に失敗し、unavailable → FIND_PATH と正しく退避しました。

#### 4.5.3 有用性（○、上向き）

- FIND_PATH の退避経路が HONO で正しい機構を安く捉えるようになりました。U5 区間は 208.4 → 170.9 s、HEI と TS の写像 RMSD は 0.325 → 0.009 Å、U6 の鞍点は 36 ステップ・106 s → 3 ステップ・6.7 s。run 全体でも 9:32 → 5:23 で、GS が成功した W7（7:22）より速くなりました。
- STO-3G HCN は 2 → 3 チャンク（+18 s）になりました。U5-N5 と閾値ノイズの副作用で、M8 の欠点ではありません。
- 無駄は残っています。HONO は 6 反復目（65 s）で落ち着いたのに残りの約 96 s が判定されずに走り（U5-N5）、NH3 では tangent 用の DFT freq（5.1 s）が今も走り（C15）、HONO の GS は 8 s × 2 回同じ例外で終わります（U5-N3・N4）。

#### 4.5.4 複雑さ（○）

- 変更は小さく、find_path は 1 反復あたり約 3 行増え、`_initial_path` は 1 行減り、IDPP は定数 1 つです。knob・policy キー・決定表の行は増えていません。actions.py はちょうど 500 行で上限にあります。U5-N8 の修正（M9）は行数を変えません。
- テストは一部良くなりました。DriftingPath が initial_path を記録するので、チャンクの継ぎ目の不具合は fake で捕まります。ただし FakePath 本体は initial_path を無視し、DriftingPath も 1 チャンク目の初期経路の形を検査しないので、U5-N8 を見逃しました（U5-T1）。
- 文書のずれが残っています（design.md:130 の gmax の記述、design.md の reaction-paths 行の impose、validation.md:55・:137 のずれの原因）。

#### 4.5.5 残る課題

- **U5-N3（medium）DLC の GS の脆さ**: 影響はコスト（GS 8 s + string 163 s）だけになりました。C35 は今回見送ります（§4.5.8）。
- **U5-N4（low）**: S4 が GS 自体の失敗でも再実行する。→ C36
- **U5-N5（low）**: 判定がチャンクの終わりだけ（HONO で約 96 s の無駄、STO-3G で余分な 1 チャンク。0.1 kcal/mol は ZTS の反復ノイズ 0.1〜0.3 に近い）。対策は見送り（§4.5.8）。
- **U5-N6（low）**: design.md:130 の「gmax と NWChem の収束判定は記録するだけ」。→ C39′
- **U5-N7（low）**: 試行を使い切った後の FIND_PATH（見送り）。
- **U5-I7・I8（low）**: 両端より低いノード（→ C16）、tangent 用の DFT freq（→ C15）。

#### 4.5.6 新たな課題

- **U5-N8（medium）`_initial_path` が内部の画像を 1 枚ずつ独立に ends[0] へ整列するので、1 チャンク目の凍結端点が NWChem によって回され、ゆがむ**: §4.5.2 のとおり。1 チャンクで結論が出る string では、行 14 の根拠が解像度（1 kcal/mol）と同じ大きさで崩れます（HONO +1.09、STO-3G +13.57 kcal/mol）。M8 は 2 チャンク目以降にしか効かないので覆いません。誤分類の実例がないので high ではなく medium です。
- **U5-N9（low）string デッキの `impose` は xyz_path を渡すと効かない。公式文書の FREEZEN の説明もソースと違う**: 公式文書の FREEZEN は「the last bead remains fixed」ですが、実際には zts_min_motion で剛体として回され、L-BFGS の外挿でゆがみます。design.md の reaction-paths 行は impose を ZTS の設定として列挙しています。挙動への害はなく、U5-N8 の診断を遅らせました。
- **U5-T1（low）テストの穴**: 1 チャンク目の初期経路の形（隣との剛体の跳び）を検査していない。double_well の端点では Kabsch の条件が悪くならないので、今の実装でも跳びは出ません。
- **U5-D1（low）validation.md:55・:137 がずれの原因と意味を誤って記している**: :55 は「NWChem が端点を大きく動かすのは IDPP から始める 1 チャンク目だけ」、:137 は「1 チャンク目の中の端点のずれ（STO-3G で 0.074 Å）は note による監視だけ（M8 の設計どおり）」。同じ IDPP でも逐次整列すればずれは 0.0037 Å で、原因の帰属が誤りです。+13.6 kcal/mol の端点を「設計どおり」とするのも誤解を招きます。

#### 4.5.7 追加の改良案

- **U5-R8（must、修正して採用）`_initial_path` を逐次整列にする** 【M9】: actions.py:455-456 を `inner = [np.asarray(i.coords, dtype=float) for i in images[1:-1]]` と `return ctx.path_file(name, align_sequential([ctx.ends[0], *inner, ctx.ends[1]]))` に置き換える（行数は変わらない）。ends[1] が剛体として回るのは M8 の継ぎ目と同じ扱いで、c2・c3 で実績がある。
  - テストの判定条件は「隣接する画像の直接の RMS ≈ mapped_rmsd（剛体の跳びがない）」にするか、内部の画像どうしの組だけで align_mapped の不変性を見る。当初案の `align_mapped(a, b) ≈ b` は、直線の ends[0] と最初の内部画像の組で Kabsch が縮退して回転が定まらないので、修正後のコードでも誤って落ちうる。直線（HCN/HNC 型）と平面（HONO 型）の IDPP で検査し、今のコードで落ちることを確かめてから入れる。
  - 同じコミットで validation.md:55・:137（U5-D1）と design.md:130 の「1 チャンク内のずれ」の説明（逐次整列の後は 0.001〜0.004 Å 程度）を直す。
  - 効果: 1 チャンク目の凍結端点のずれ STO-3G 0.074 → 0.004 Å、HONO 0.022 → 0.0014 Å、端点の持ち上がり +13.6 → +0.04、+1.09 → +0.009 kcal/mol。追加の計算はなく、bead のプロファイルはほぼ同じ（HONO の最大 13.47 → 13.46）。GS 経路（xTB の座標系）も同じ処理で扱える。
  - 費用: コード 2 行の置き換え、テスト約 10 行、文書 2 か所。string の JobStore キー（initial_path の sha256）が変わるので、既存 run を再開すると string を 1 回計算し直す（validation に明記する）。中立。
- **U5-R9（could）文書の訂正: impose の実態・freezeN・gmax** 【C39′、M9 と同じコミットで】: design.md の reaction-paths 行に「impose は xyz_path を渡すと NWChem が使わない。実際の整列は NWChem の zts_min_motion（凍結 bead N も含む）」を 1 文。design.md:130 の gmax の記述を「NWChem の string は gmax を記録しない（pysis_gs だけ）。NWChem の converged は使わない」に（U5-N6）。デッキは変えない。
- **C36（could、v2 から継続）S4 の再実行を GS が完了した場合に限る** 【C36】: worker.run_growing_string の except で、tsopt がないか final_geometries.trj がなければ raise。1〜2 行、中立。
- **C15（could、v2 から継続）tangent_mode_cm1 を barrierless のときだけ計算** 【C15】: 17 原子の縮退反応では 1 本約 900 s を省ける。約 5 行、中立。
- **C16（could、修正して採用）両端より低い DFT//xTB ノードがあれば unavailable** 【C16】: fake の E2E テスト（FIND_PATH → multi_max → VALIDATE_INTERMEDIATE）を同時に入れることを条件にし、用意できなければ見送る。M9 を先に入れる。3 行とテスト 1 件、【複雑さ増】。

#### 4.5.8 見送った案

- **C35（pysis_gs の平面端点に off_plane）**: S18 で退避経路の FIND_PATH がねじれ機構を正しく捉えるようになり、結論の誤りはなくなった。得られるのは平面の小分子で約 2 分のコストだけで、標的の錯体（2 断片以上）は cart の GS なので DLC の脆さは関係しない。証拠も 1 系だけで、pysis 専用の摂動を約 8 行足すのに見合わない。
- **string デッキから impose を削る**: 効いていないので挙動は変わらない。文書で実態を書けば足りる。
- **end_drift をゲートにする、ずれた端点で追加の SP をとる**: M9 と M8 で原因を除けばずれは 0.001〜0.004 Å（+0.01〜0.04 kcal/mol）で、追加の計算や閾値に見合わない。
- **_SETTLED を 0.1 → 0.3 に緩める**: 根拠は STO-3G の玩具系の 1 例（+18 s）だけで、形の判定の余裕が減る。
- **ZTS を maxiter 10 × 最大 6 チャンクにする（U5-N5 の対策）**: 再開のたびに過渡が入り（STO-3G の c2 の最初の 2 反復で max|ΔE| 6.2 / 6.0 kcal/mol）、bead の SCF もやり直しになる。HONO での利得は約 80 s。
- **行 16 に saddle の試行数のガードを付ける（U5-N7）**: 実例がなく、予算を使い切った後でも monotonic・multi_max の結論には価値がある。
- **IDPP の収束フラグを返して note に残す**: 5,000 反復で実 run の HONO・HCN はどちらも収束した。

#### 4.5.9 変更不要な点

- M8 + S19 の実装（成功したチャンクごとに ctx.ends へ戻し、内部 bead を逐次整列。bead のエネルギーは NWChem の値のまま使い、凍結 bead の古い値は真の極小のエネルギーに等しい）と、end_drift の note を監視だけに使うこと。
- S18 と、その他の IDPP 定数（_FORCE_TOL 1e-3、キック 0.01 Å、MIN_DISTANCE_A 0.7）。
- M2 の骨格（monitor なし、bead エネルギーの安定だけで打ち切り、形を収束と無関係に記録、種は single_max の HEI だけ）と S3（M9 で 1 チャンク目の前提も満たせる）。
- 決定表 17 行の構造と、FIND_PATH を 1 ケース最大 1 回とすること。
- 仮説選択（優先順、ΔE ≤ 40、0.2 Å / 30°、組成あたり 6 件、F1）、SCREEN（低レベル TS の近道、xTB freq による受理、_COLLAPSE_A 0.05 Å、S2、barrier_verdict の判定式）、S4 の骨格。
- 既定値（string 9 bead × maxiter 20 × 3 チャンク、stepsize 0.05、interpol 3、freeze1 / freezeN、GS 11 ノード、off_axis）と DriftingPath のテスト。

---

### 4.6 U6 反応経路：鞍点精密化・TS 検証・接続確認（REFINE_SADDLE / VALIDATE_TS / CONNECT）

評価: v2 ○ / ◎ / ○ → **v3 ◎ / ◎ / ○**

#### 4.6.1 変更点

- **S20**: `engine._moddir`（engine.py:87-95）と `vibrations.cartesian_mode_number`（vibrations.py:104-114）を新設。mode_index が None なら moddir 0、imaginary_modes が 1 本なら autoz・moddir 1、2 本以上なら noautoz と P·H·P（質量加重しない、剛体運動を射影）の負の固有ベクトルのうち選んだモードとの |cos| が最大のものの 1 始まりの順位。key_payload を {method, molecule, hessian sha, moddir, cartesian} に（v2 の saddle のキャッシュは再利用されない）。TIMEOUT / MAXITER の継続では moddir を min(k, 1) に（golden テスト）。
- **S21**: QRC の両側の opt が `trust 0.3` + `inhess 2`。`_register` と mode-follow の opt は trust 0.1 のまま。
- **C39 の一部**: continuation の docstring を「autoz → 同じ種から Cartesian」に訂正（engine.py:185-187）。
- **波及（M8・S19・S18）**: path_hei の種と接線は、端点を戻して逐次整列した経路から作られる。
- 変わっていないもの: 決定表、driver（予算と分割）、refine_saddle・validate_ts・connect・`_assign`・`_register` の本体、is_first_order_saddle・connection・qrc_drop・qrc_amplitude・periodic_nearest。

#### 4.6.2 妥当性（◎、v2 は ○）

v2 で留保した U6-N1（Cartesian で moddir 0 のときの迷走）、U6-N2（string 画像の未整列）、U6-N3（凍結端点のずれ）の 3 件とも、コードと実計算の両方で解消を確かめました。

- **S20**: v3 の saddle デッキ 4 件（bfa49214 / 5b21f557 / 99b68d44 / ae0c3a2a）はすべて `maxiter 50 trust 0.1 sadstp 0.1 inhess 2 moddir 1`、geometry は autoz で、stdout に 'initial eigen-mode to follow (moddir) = 1' が出ています。ステップ数は HCN 5 → 4、NH3 3 → 2、STO-3G 4 → 3、HONO 34 + 36 → 3、E_TS の差は ≤ 2.5e-7 Eh です。
- **NWChem 7.2.3 のソース（src/driver/opt_drv.F）で番号の整合を確かめました。** Cartesian の射影子は重みなしで重心まわりの並進・回転 6 本を Gram-Schmidt 直交化した P（2760-2850）、Hessian は P·H·P + 1000(1−P)（2948-2963）、対角化は既定で util_jacobi を使い昇順に並べ（1138-1139）、1 手目は `e(moddir)` の固有ベクトルを追います（3270-3274）。`cartesian_mode_number` の `external_basis(["H"]*n)`（等質量で重心まわり、同じ部分空間）と昇順の数え方は、NWChem の番号と一致します。プローブ `w1b_moddir` の stdout では、step 1 の 'The mode being followed to the saddle point' のベクトルが、moddir 1 / 2 のどちらでも P·H·P の第 k 固有ベクトルと |cos| = 1.0000 で一致しました。
- 補足: ofirstneg は既定で .true.（1112-1113、3299-3305）なので、負の固有値が 1 本だけになると NWChem は mode 1 に切り替えます。負モード 2 本以上で明示した moddir が効くのは初めの数ステップだけで、結果の誤りは VALIDATE_TS と connection のゲートが防ぎます。
- **S21**: QRC の v3 デッキ 8 件はすべて `maxiter 100 trust 0.3 inhess 2`。8 側すべてで max(traj) = traj[0]、変位直後の降下は 2.72〜3.45e-4 Eh、最終エネルギーは v2 と 8e-8 Eh 以内、割付けも同じ（elementary 3 件、NH3 は degenerate）でした。
- **S19 と M8 の波及**: STO-3G の SaddleClaim の notes が `['mode_overlap_below_0.3']` → `[]`。
- U6-N1 の元の状況（autoz 失敗後の Cartesian）は v3 では起きておらず、その経路はプローブでだけ確かめました。INPUT_INVALID(autoz) の継続では moddir 1 がそのまま残りますが、質量加重と非加重の Hessian は合同なので（Sylvester の慣性法則）負の固有値の本数は同じで、moddir 1 は Cartesian でも正しい。
- 留保: 実計算で通ったのは 3〜4 原子の共有結合系だけで、錯体の TS はまだ U6 に入っていません（TMA·(HF)₂ の case は U6 に入る前に same_basin で終わる）。負モード 2 本以上の分岐はプローブだけです。柔らかい負モードは NWChem の仕様で追えないことがあります（U6-N9）。どれも正しさではなく効率の問題です。

#### 4.6.3 有用性（◎）

- U6 区間の壁時間は合計 438 → 181 s（−59%）: HONO 298.3 → 91.2 s、HCN 52.2 → 35.7 s、NH3 70.1 → 41.2 s、STO-3G 17.5 → 13.4 s。
- HONO の改善は、S18 と S20 で xTB Hessian が採用されるようになったことによります（v2 は DFT freq 15.1 s、AUTOZ の失敗 48.4 s、Cartesian での迷走 57.4 s）。
- QRC は 3 反応の合計で 204 → 107 ステップ、265.9 → 128.2 s で、v2 のレビューの予測（約 110 s）とほぼ同じでした。下流の ΔG‡ の差は ≤ 0.002 kcal/mol で、結論は変わりません。
- 17 原子（1 勾配約 45 s）なら QRC だけで 1 反応あたり約 30 分の短縮になる見込みです（外挿）。NH3 は trust 0.3 でも 25 ステップ中 22〜23 ステップで負モードの刻みが制限されており、残りのコストは負の曲率の領域の広さによるものです。

#### 4.6.4 複雑さ（○、中立）

- 差分は engine.py +11 行、vibrations.py +12 行、input.py −1 行。radon で CC > 15 の関数はなく、決定表は純関数のままです。
- v2 の複雑さの主因だった「mode_index をそのまま moddir に流用する暗黙の変換」と「`if mode_index` のときだけ noautoz にする分岐」が `_moddir` 1 か所にまとまり、saddle の key_payload も明示的になりました。
- 新しく、モードの表現が 2 つ共存します（質量加重の射影モードは選択とゲートに、重みなしの P·H·P は負モード 2 本以上の番号付けにだけ使う）。使う場所は 1 か所で、docstring と design.md に書かれているので許容できます。本数の数え方が 2 通りある点（U6-N8）が小さな不整合として残ります。
- smoke の `test_moddir_follows_the_second_negative_mode` は Evidence が返ることしか検査していません（U6-N11）。

#### 4.6.5 残る課題

- **U6-I3（low）**: periodic_nearest に許容幅がない（C18、実 run で未到達）。v2 の判断（could のまま据え置き）を維持します。
- **U6-I5（low）**: `_register` を connection のゲートより前に呼ぶ（設計どおり）。
- **U6-I7（low）**: = U5-N7。
- **U6-N6（low）**: qrc_bounds_A の上限と HESSIAN_NEAR_A（0.5 Å）の関係が文書にない。→ C39′
- **U6-N7（low）**: 子反応の予算（case ごとに新しく 6 h、最悪 7 case）が文書にない。→ C39′
- v2 の U6-R4 の残り: validation.md:69 の「127(maxiter で継続)/ 17 = 144」（正しくは 126 / 17 = 143）。→ C39′

#### 4.6.6 新たな課題

- **U6-N8（low）autoz / Cartesian の分岐は ν < 0 のすべてのモード（noise 以下を含む）の本数で決まり、モードの選択（ν < −noise）と数え方が違う**: `_moddir` は `len(hessian.imaginary_modes) == 1` で分岐し（engine.py:93。imaginary_modes は f < 0 のすべて）、mode_index と `_first_hessian` は ν < −noise_cm1 で判定します（actions.py:265-266, 282）。本物の虚振動 1 本に −1〜−10 cm⁻¹ のモードが付いた種は Cartesian の分岐に入ります（番号はほぼ 1）。imaginary_modes は昇順なので添字はずれず、追うモードは正しい。影響は効率（v2 の HONO は Cartesian で 36 ステップ）だけで、design.md:144 の「負モードが 1 本なら autoz」とは食い違います。W7・v2・v3 の freq ジョブ 115 件に −10 ≤ ν < 0 のモードは 0 件で、実 run では起きていません。
- **U6-N9（low）NWChem の saddle は |固有値| < 1e-4 かつ勾配の小さいモードをゼロとして末尾へ回すので、柔らかい負モードは moddir を明示しても追えない**: opt_drv.F の smalleig = 1d-4（3135、ソースの parameter で入力から変えられない）。|e| < smalleig かつ |gv| < gmax_tol（既定 4.5e-4）のモードをゼロとして数え、末尾に回して e = 1000 にします（3210-3240）。負のモードが soft だけなら e(1) は最低の正のモードに置き換わり、moddir 1 はそのモードを追います。ねじれの TS や錯体の柔らかい負モードでは U6-N1 と同じ迷走になりえます。S20 以前からある NWChem の制約で、実例はなく、ゲートが結果の誤りを防ぎます。
- **U6-N10（low）trust がジョブキーに入らないので、古い run dir で再開すると QRC は trust 0.1 の結果を再利用する**: validation の残る課題 2 と同じ。同じ極小に着くので害はありません。saddle はキーが変わったので再計算されます。
- **U6-N11（low）smoke の moddir テストが追ったモードを検査していない**: tests/smoke/test_real_nwchem.py:66-72 は負モードが 2 本あることと Evidence が返ることしか確かめていません。S20 の前提（NWChem の Cartesian での番号 = P·H·P の昇順）を確かめたのは本レビューの手作業の照合だけで、NWChem を更新して前提が崩れても気づけません。

#### 4.6.7 追加の改良案

- **U6-R5（could、修正して採用）文書の訂正の残り（C39 の U6 分）と moddir の分岐の数え方** 【C39′】: (1) design.md の reaction-paths 行に「autoz / Cartesian の分岐は ν < 0 のモード（noise 以下を含む）の本数で決める」「NWChem 7.2.3 の saddle は |固有値| < 1e-4 かつ勾配 < gmax_tol のモードをゼロとして扱う（opt_drv.F の smalleig）ので、柔らかい負モードは追えないことがある」。Cartesian でしか成り立たない振動数の目安（数十 cm⁻¹）は書かない。(2) 分割した子反応はそれぞれ 6 h の予算（最悪 7 case）、qrc_bounds_A の上限は HESSIAN_NEAR_A（0.5 Å）以下。(3) validation.md:69 の 144 → 143。コードは変えない。文書 4〜5 行、中立。
- **U6-R6（could）smoke の moddir テストで、NWChem が追ったモードを P·H·P と照合する** 【C42】: saddle の stdout の最初の 'The mode being followed to the saddle point' ブロックを読み、P·H·P の第 2 固有ベクトルとの |cos| > 0.99 を検査する。stdout を読む処理はテストの中の helper に置き、本体の output パーサには足さない。1 ブロック目のベクトルが nvar = 3N の Cartesian で出力されることを前提にする旨をコメントに書く。smoke に約 8 行（real マーカー、WSL だけ）、中立。

#### 4.6.8 見送った案

- **`_moddir` の本数を ν < −noise で数えるようにコードを変える（U6-N8）**: engine は gates の noise_cm1 を持たないので、Protocol に引数を足すか定数を二重に持つことになる。Cartesian の分岐でも追うモードは正しく、実 run での発生は 0/115。文書で足りる。
- **trust を optimize のジョブキーに入れる（U6-N10）**: v2・v3 のキャッシュがすべて無効になり、再利用される結果は同じ極小でゲートも通る。
- **柔らかい負モード（U6-N9）への対処として vardir を使う、Hessian を尺度変更する、smalleig を変える**: smalleig は入力から変えられず、Hessian の加工は v2 でも見送った。実例がない。
- **QRC の trust をさらに上げる**: NH3 の QRC は 34 s と小さく、平坦な錯体の PES で浅い障壁を越えて別の basin に落ちる危険が増える。
- **autoz をやめ、saddle を常に noautoz と P·H·P の番号にする**: 負モード 1 本の autoz・moddir 1 は 2〜4 ステップで収束しており、Cartesian にすると遅くなりうる。
- **`_register` と mode-follow の opt も trust 0.3 にする**: 対角推定の Hessian では変位点から TS より上へ戻る例がある。
- **負モード 2 本以上の分岐を合成したケースで pipeline 全体に通す**: 番号付けはプローブと stdout の照合で直接確かめられた。U6-R6 の smoke の検査の方が安く、回帰も防げる。

#### 4.6.9 変更不要な点

- 決定表、予算（saddle 2 回、振幅 2 通り、分割の深さ 2）、純関数の decide。
- `_first_hessian` の採否（1 回目だけ、ν < −noise がちょうど 1 本で重なり ≥ 0.3）。重なりは選択と note にだけ使う。
- S20 の `_moddir` 1 か所への集約と key_payload {moddir, cartesian}、継続の規則（autoz の失敗は同じ種と Hessian から Cartesian、TIMEOUT / MAXITER は最新フレーム + drv.hess + 旧 movecs で moddir を min(k, 1)、SCF は damping と level shift）。
- S21 の範囲、QRC に TS の freq Hessian を渡すこと（HESSIAN_NEAR_A 0.5 Å の検査付き）と、エネルギー目標から決める振幅（C19）。
- VALIDATE_TS（別ジョブの DFT freq と hfauto の射影、is_first_order_saddle）、connection のゲートと割付けのラベル、higher_order のときの 0.1 Å の再試行。

---

### 4.7 U7 熱化学（thermo）

評価: v2 ◎ / ○ / ○ → **v3 ◎ / ○ / ○（据え置き）**

#### 4.7.1 変更点

- thermo 単位のコード変更は 0 行です（thermochemistry.py、chemistry/thermo.py、backends/goodvibes/、core/method.py の ThermoSettings、core/system.py の Conditions、core/records.py へのコミットなし）。thermo stage の config_sha も全 run で v2 と同じでした。
- 下流の S22（U8 担当）で、ranking.csv と report.html の ΔG_assoc −11.53 の隣の誤った blocker が消えました。thermo stage の ReactionThermo.blockers は v2・v3 とも [] です。
- S23（文書）で composite G の手順が書き直されました。
- 間接の変化: S20・S18 で TS の構造がわずかに変わり、TS の振動数と E が動きました（HONO の虚振動 −678.1 → −681.0、NH3 −758.0 → −757.7）。M8・S19 で STO-3G の TS の SpeciesThermo.notes も `['mode_overlap_below_0.3']` → `[]` になりました（SaddleClaim.notes をそのまま合流させるため）。
- 未実施: C20、C21、C22、C28、C29、C30、U7-R5 の文書。

#### 4.7.2 妥当性（◎、据え置き。裏付けが強まった）

- thermo_frequencies の固定規則、vib_scale 1 本を freq_scale_factor と zpe_scale_factor の両方に渡すこと、thermo_consistent と fail-closed は v2 と同じです。
- **GoodVibes 4.3.0 の実装を WSL で確認しました**（prod venv の goodvibes/thermo.py L697-747 と api.compute_thermo のシグネチャ）。既定値は QH=False・symm=False で、hfauto は symm=True を明示しています。calc_bbe は ZPE を zpe_scale_fac で、U_vib と S_vib を scale_fac で計算するので、vib_scale を 1 本にすると両者がそろいます。QS='truhlar' が変えるのは S_vib だけで、band は「エントロピーの扱いの感度」だけを表します（設計どおり）。
- **M6 の約束**: S20 は鞍点の最適化中に追うモードを変えるだけで、freq の後の約束（SaddleClaim.imag_cm1 = min、QRC は imaginary_modes[0]）は変えていません。participants() が TS を含めるのは QRC の接続で最低モードの意味が確かめられた outcome だけなので、thermo が捨てるモードは QRC で確かめたモードと同じです。
- **v3 の 15 subject**: 失敗 0、attempt はすべて 1 回。整合ゲートの余裕は |ΔZPE| ≤ 1.63e-8 Eh、|ΔE| = 0、n_real はすべて一致。ΔG の v2 からの差（≤ 0.0021）は、S20・S18 による TS の変化と、TMA 単量体の NWChem 再実行での E のずれ（4.96e-7 Eh。同じ構造・同じ LOT で SP と freq を比べても 5.8e-7 Eh ずれる、SCF 1e-7 での再現性の範囲）で説明できます。NH3 の TS は構造が変わっても ΔG‡ の差が −0.0013 で、σ は安定しています（σ が変われば RT ln 2 = 0.41 ずれる）。
- **新しい実計算の証拠**（`u7_assess`、v3 の run を cp -a で複製して追記）:
  - `hcn_tz`: S23 の手順どおりに sp_tz（TZVPD、all_minima）→ thermo_tz（energy_method）→ report_tz を追記。1:26、SP 3 件 miss（83 s）、GoodVibes 3/3 hit（FileRef の複製も指紋の検査を通った）。ΔE‡ 46.62 → 46.43、ΔG‡ 42.18 → 41.99、ΔG_rxn 12.68 → 13.33、band 幅 0.00057 は同じ、rank 1、blockers []。
  - `tma_d`（既定の targets）: min_hf と min_tma は G=None、notes [thermo_unavailable, energy_layer_missing]（S9 が初めて実 run で通った）。ΔG_assoc は None で理由は出ない。
  - `tma_a`（all_minima、SVPD）: ΔG_assoc −11.52885（元の −11.52869 との差 −0.00016）。energy 層での state=False の LOT 判定が通った。
  - `tma_tz`（all_minima、TZVPD）: ΔG_assoc −9.0394（3:58、rc 0）。
- 残る弱点（low、潜在的）: U7-N1（混成 ΔE）、U7-N3（(Hill, 電荷) の対応）、U7-N9（単量体アンサンブルの異性体）。M6 の副次的な負モード、電荷・開殻の会合量、複数の温度、1 bar / 1 M は実計算で通っていません。

#### 4.7.3 有用性（○、据え置き）

- 良くなった点: S22 で主対象の報告行が `outcome:same_basin` だけになり、STO-3G の TS の SpeciesThermo.notes から誤った note が消えました。composite G の手順が HCN・TMA の複製で動き、追加のコストは SP だけです（HCN 3 点の TZVPD で 83 s、TMA·(HF)₂ 3 点で約 4 分）。1 subject あたり 0.22〜0.28 s で無視できます。
- ○ にとどめる理由: ΔG_assoc の感度の幅が出ない（298 K で 1.155）、ΔG_assoc=None の理由が出ない（`tma_d` で再現）、ReactionRecord のない錯体では ΔG_assoc がどこにも出ない（U7-N8）、点群が記録されない。
- 今回の新しい知見: 主対象の ΔG_assoc は SVPD の −11.53 に対し TZVPD//SVPD の composite で −9.04（+2.49）で、band の幅 1.155 よりも手法誤差の目安（約 2）よりも大きい差です。composite の report.html と ranking.csv には energy 層の LOT が表示されません（U7-N10）。

#### 4.7.4 複雑さ（○、据え置き）

- 622 行（stage 253、chemistry/thermo.py 131、engine 129、worker 108、__init__ 1）と knob は v2 と同じ。radon の CC は最大 C(14) で、すべて ruff の上限 15 未満。隣の rankable も S22 の後で C(11) のまま。
- composite は新しいコードを足さずに、後の artifact が勝つマージ（Manifest.union）と content-addressed な JobStore だけで成り立ち、今回の実 run の複製で動いたので、構造としては簡素です。
- ◎ にしない理由は、使われない経路が残っていることです（ThermoResult.notes は常に ()、S_rot は返すだけで誰も読まない、population は書くだけ。C20 で片づく）。

#### 4.7.5 残る課題

- **U7-N1（= U0-N1、low）**: 層が欠けた subject の混成 ΔE。→ C28
- **U7-N2・N3（low）**: ΔG_assoc=None の理由が出ない、(Hill, 電荷) の後勝ち。→ C29
- **U7-N4（low）**: Conditions・ThermoSettings の範囲検査。→ C30
- **U7-N5・N7・N8・N9（low）**: vib_scale の根拠、validation の帰属「M6・S15」（validation.md:76）、反応のない錯体では ΔG_assoc が出ない、単量体アンサンブルの異性体が文書にない。→ S26
- **U7-N6・U7-I7（low）**: ThermoResult.notes・S_rot が使われない、点群が記録されない。→ C20
- **U7-I5（low）**: ΔG_assoc の感度幅。→ 数値の文書化（S26）で代替
- U7-I8（基準が反応物の極小）と U7-I11（単原子種）は方針どおり見送り。

#### 4.7.6 新たな課題

- **U7-N10（low）composite の結果に energy 層の LOT が表示されず、SVPD の report と見分けにくい**: `hcn_tz/report_tz/report.html` は ΔE‡ 46.43・ΔG‡ 41.99 を「view of stage report_tz」の見出しだけで表示し、'def2' の文字列は 0 件、ranking.csv にも LOT の列がありません（method_panel.csv には level が出る）。LOT は SpeciesThermo.energy_calc → Evidence.level を辿らないと分かりません。stage id が report_tz なので完全に見分けられないわけではありませんが、S23 で composite が公式の手順になったので取り違えのおそれがあります。
- **U7-N11（low、訂正後の内容）主対象の ΔG_assoc の LOT 依存が検証記録に注記されていない**: `tma_tz` で ΔG_assoc は SVPD −11.529 → TZVPD//SVPD −9.039（+2.49）で、差は band の幅 1.155 よりも手法誤差の目安（約 2）よりも大きい。ただし design.md:149 は SVPD の過大評価と BSSE 未補正をすでに定性的に書いているので、欠陥は validation.md:18・:49 が −11.53 を注記なしで載せていることに限られます。TZVPD//SVPD も BSSE を補正していない composite であり、−9.04 を正しい値として扱ってはいけません。
- 検証で追加された注意（missed）: thermo の artifact id に stage id が入らないので、S23 で複製した run では thermo_tz の値が元の thermo の値を上書きします。複製した run で元の report stage を再実行すると、SVPD ではなく composite の値が出ます。design.md:150 はこの点に触れていません。

#### 4.7.7 追加の改良案

- **U7-R5′（should）thermo の文書の訂正と、S23・composite の実測値の記録** 【S26】: (a) validation.md:76 の「M6・S15」を「S15（M6 の寄与は 0。負の副次モードがないため）」に直し、:84 の表の注記（ΔG_assoc +0.07）も S15 だけの帰属に合わせる。(b) design.md:146 の thermo 行に「vib_scale=1.0 は未スケールの調和振動数。PBE0/def2-SVPD の文献の ZPE 因子は 0.985（Kesharwani 2015）で、切り替えても ΔG‡ は ≤ 0.06、ΔG_assoc は ≤ 0.09 しか変わらない」「ΔG_assoc の qRRHO の扱いの幅は 298 K で約 1.2 kcal/mol」「単原子の単量体を含む組成と、反応を持たない錯体では ΔG_assoc は出ない。単量体のアンサンブルには同じ (Hill, q, m) の異性体も入る」を足す。(c) validation.md:62 の「S23 は文書だけ」をプローブの結果に置き換え（HCN の TZVPD composite で ΔG‡ 41.99、GoodVibes 3/3 再利用。TMA·(HF)₂ の ΔG_assoc は TZVPD//SVPD・BSSE 未補正のプローブ値で −9.04。既定の targets では単量体が energy_layer_missing）、表 #5 の −11.53 に「SVPD、BSSE 未補正」と注記する。(d) design.md:150 に「複製した run では thermo_tz が元の thermo を上書きする（元の report を再実行すると composite の値が出る）」を 1 文。文書 6〜10 行、中立。
- **U7-R7（could、U8 側）composite のときは report の見出しに energy 層の LOT を出す** 【C44】: report stage で、使った SpeciesThermo の energy_calc の Evidence.level が freq と違うときだけ集め、report.html の meta 行に「energy: <LOT> // freq: <LOT>」を 1 行。ReactionThermo にフィールドは足さない。約 5 行と html のテストの assert 1 行、中立。
- **U7-R1（could）層が欠けた subject では ΔE も出さない** 【C28】: v2 と同じ。1 行、assert 1 つ。
- **U7-R2（could）ΔG_assoc=None の理由を 1 つ記録し、曖昧な対応を fail-closed にする** 【C29】: v2 の修正版と同じ（ReactionThermo.notes に `assoc_unavailable:<no_composition|monomer_missing|monomer_thermo_unavailable|mixed_level>`、blockers には入れない）。`tma_d` で無言の None が再現したので診断の価値がある。約 10 行、テスト 1 件、【複雑さ増】。
- **U7-R3（could）Conditions と ThermoSettings に宣言的な値の制約** 【C30】: v2 と同じ。約 8 行、テスト 2 件（U2 と合わせて）、中立。
- **U7-R4（could）点群を ThermoResult.notes で記録し、S_rot を protocol から外す** 【C20】: v2 と同じ。約 4 行、中立。
- **U7-R6（could、修正して採用）ΔG_assoc の感度幅**: コードで assoc_band を足すのはやめ、既知の幅（298 K で 1.155）を S26 の文書に数値で書くだけにする。LOT の差（2.49）の方が大きいと分かったので、フィールドを増やす利益は薄い。

#### 4.7.8 見送った案

- **discover.yaml に TZVPD の composite を既定で組み込み、ΔG_assoc を TZVPD//SVPD で出す**: 停留点の LOT を方針ごと変えることになり、validation の値との連続性が崩れる。1 run に置ける熱化学の層は実質 1 つで SVPD の値が見えなくなる。S23 の追記手順で必要な人が約 4 分で得られる。
- **ReactionThermo または ThermoConfig に energy_level のフィールドを足す**: SpeciesThermo.energy_calc から辿れる。スキーマを変えると古い結果の読み込みにも影響する。表示の側で導けば足りる（C44）。
- **thermo の artifact id に stage id を含め、1 run に複数の thermo 層を並べる**: report と panel_report が反応 id と T で ReactionThermo を引く前提が崩れ、選ぶ規則が新たに要る。
- **S21 の trust を thermo の JobStore キーに入れる**: thermo のキーには freq の job_key と Hessian の sha が入っており、thermo の入力はすべて freq の出力で決まる。
- **composite の E_SP と E_freq の差（同じ LOT・同じ構造）をゲートにする**: 実測の差は ≤ 0.0004 kcal/mol で判定に効かない。
- **S20 に合わせて thermo でも反応座標のモードを overlap で選び直す**: thermo が捨てるモードは QRC で確かめたモードと常に同じ。
- **既定の vib_scale を 0.985 にする**: 効果は ΔG‡ ≤ 0.054、ΔG_assoc 0.089 で、LOT の差 2.49 より 1 桁以上小さい。

#### 4.7.9 変更不要な点

- thermo_frequencies の純関数と固定規則、GoodVibes に invert=None で渡すこと、vib_scale 1 本を両方の因子に渡すこと。
- qs × cutoff の 6 variant の band、JobStore のキー {freq の job_key, hessian_sha, 渡す振動数, 温度, settings}（composite の追記で GoodVibes の結果が 3/3 再利用された）。
- 温度は conditions だけから取ること、ThermoSettings の frozen / extra='forbid'、thermo_consistent と fail-closed、失敗しても再試行しないこと。
- S9 の energy_layer_missing、S14 の単量体キー (Hill, q, m) と state=False、composite G = E_SP + (G_GV − E_GV) と S23 の手順。
- 1 atm で計算し 1 bar / 1 M に換算すること、population を組成と LOT ごとに計算すること、SaddleClaim / MinimumRecord の notes の合流、1 subject に 1 worker。

---

### 4.8 U8 一点計算・手法パネル・順位付けと報告・実行基盤（sp・method_panel・report・execution・configs）

評価: v2 ◎ / ○ / ○ → **v3 ◎ / ○ / ○（有用性は ○ の上側）**

#### 4.8.1 変更点

- **C2**: `render_wft` の ccsd ブロックに `  maxiter 50`（input.py:159, 166-167）。テストは test_nwchem_input.py。
- **S22**: gates.rankable を書き換え（gates.py:304-317、333 → 335 行）。test_gates.py に 7〜8 行（same_basin、same_basin + spin_contaminated、thermo 欠落で thermo_unavailable が 1 回だけ）。
- **S23（文書）**: design.md:47・145・147・149・150、README:37・86、method_panel.yaml:5。
- **S25**（U1 担当だがコードは実行基盤の runner.py:142-143, 152-157、+4 行）。
- 変わっていないもの: single_point.py、report.py、summary.py、html.py、execution/*、preflight.py、layout.py、cli。
- 未実施: C25、C27、C37、U8-N4（パネルは正式な検証で未実行）。

#### 4.8.2 妥当性（◎、維持）

- **C2**: 実ジョブ（`u0_c2/run/jobs/15/15c10003…`、`jobs/34/3441d771…`）のデッキに `ccsd / freeze atomic / maxiter 50 / end`、NWChem のエコーは `maxit = 50`・`convi = 0.100E-05`。どちらも `Failed reading restart vector from ./job.t2` → `Using MP2 initial guess vector` から始まり、HNC/TZVPD（83 関数）は 11 反復、HONO TS/TZVPD（126 関数）は 18 反復で収束しました。HONO では rms が 16 反復目 1.9e-6 → 17 反復目 4.5e-6 → 18 反復目 4.3e-7 と非単調で、旧既定の 20 では余裕が 2 しかなかったことが実データで確かめられました。公式文書（CCSD.html）でも既定は MAXITER 20・THRESH 1e-6、閉殻（RHF）専用です。
- **未収束時の挙動**は §4.0.2 のとおり NONZERO_EXIT で、failures_by_kind に数えられ（jobs.py:163-170）、`--retry-failed nonzero_exit` で stage が再実行されます（runner.py:212）。未収束の値が黙って入る経路はありません。maxiter は描画だけの入力で key_payload に入らないので、収束済みの結果を再利用するのは正しい扱いです。
- **S22**: dzpe が None のとき dzpe_out_of_tolerance を付けなくても、順位に抜け道はできません。SpeciesThermo は G と ZPE を常に一緒に None にし（thermochemistry.py:110-114）、reaction_delta は TS があるときだけ dG_act と dzpe を返す（thermo.py:109-131）ので、dG_act がある反応では dzpe もあります。
- **変わっていない部分の数値**（v3 の run）: ranking.csv は HCN 42.17825、HONO 12.22557、NH3 3.83946（degenerate）、STO-3G 62.16778。method_panel.csv は各 run 1 行で sign_disagreement は False。JobStore の hits/misses は HCN 2/27、HONO 0/16、NH3 1/26、水 0/4、STO-3G 0/15、TMA 2/58、A/B 0/11。failures_by_kind は HONO の nonzero_exit 1 件（GS、v2 と同じ）だけで、report はすべて done（0〜1 s）。
- Manifest.union は後の artifact を優先する（manifest.py:80-91）ので、S23 の composite 手順（複製に thermo_tz を追記）では report_tz と html が同じ ReactionThermo を読みます。NWChem の `memory total` はプロセスごとの量（公式文書 Memory.html）で、4 ranks × 1200 MB = 4.8 GB は 11 GB の範囲に収まります。

#### 4.8.3 有用性（○、上側）

- v2 で有用性を下げていた最大の要因（U8-N3）が解消しました。主対象と水の報告行は `outcome:same_basin` だけになり、blockers 列は「順位を付けられない本当の理由」だけを示します（HTML の差 43 B は消えた 2 語の長さと一致）。
- S25 で入力の不正の原因が stage の例外に出るようになり、C2 で opt-in の CCSD(T) 参照点を反復切れで落とす可能性が小さくなりました（余裕 2 → 32 反復）。
- ○ にとどまる理由: (a) ranking.csv に T・standard_state・dG_rxn がなく、metric=dG_rxn のとき並べた値が表に出ない（C25）。(b) 手法パネルは v3 の正式検証でも主対象でも実行されていない（U8-N4）。(c) validation.md:62・:96・:139 が C2 と手法パネルを「実計算で未確認」としたまま。(d) S22 の後の conditions の記述が不正確（U8-N6）。

#### 4.8.4 複雑さ（○、維持）

- 差分は gates.py +2 行、render_wft は条件付きの 1 要素、runner.py +4 行、method_panel.yaml はコメント 1 行。設定項目、失敗の種類、状態は増えていません。
- S22 で分岐が 1 つ増えましたが、規則は「順位を付けない outcome は outcome と本当の blocker だけ」「dzpe は値があるときだけ検査」の 2 文に整理され、blockers 列の意味は単純になりました。
- 削減の余地は v2 と同じく残っています（monitor 経路と FailureKind.STAGNATED、使われない site.memory_mb）。

#### 4.8.5 残る課題

- **U8-7・U8-N5（low）**: ranking.csv に T・standard_state・dG_rxn がない。→ C25
- **U8-9（low）**: site.memory_mb が未使用で、wsl_local.yaml:3-4 のコメントは効かない値の書き換えを勧めたまま。→ C27′
- **U8-10（low）**: method_panel.yaml が独自の conditions（受容済み。文言は U8-N6）。
- **U8-N2（low）**: monitor 経路と STAGNATED が到達しないコード（6 アダプタの monitor() はすべて None、v3 の command_result.json の stopped は常に null）。→ C37
- **U8-N4（low）**: 手法パネルが正式な検証でも主対象でも未実行。→ C43

#### 4.8.6 新たな課題

- **U8-N6（low）S22 の後、「conditions が違えば全反応が thermo_unavailable」という記述が順位外の outcome に当てはまらない**: README.md:37、design.md:147、method_panel.yaml:5 はそう書いていますが、rank_rows は (T, state) が一致しないと t=None にし（summary.py:101-105）、rankable(same_basin, None) は S22 で `("outcome:same_basin",)` だけを返します（gates.py:306-308）。主対象（same_basin のみ）では、conditions の違う panel_report の行がそろえた場合と区別できません。report.html の ΔG_assoc は view 順で最初の ReactionThermo（元の conditions）から取る（html.py:241-255）ので、値そのものは誤りません。S22 と S23 が同じ round に入ったことによる文書のずれです。
- **U8-N7（low、レビュー入力の訂正）CCSD の未収束が NONZERO_EXIT になる根拠の記述**: 仕様整理はバイナリ内の文字列 `ccsd_energy_loc: maxiter exceeded` を根拠に挙げていましたが、これは TCE モジュールのもので、hfauto が使う `task ccsd(t) energy`（非 TCE の src/ccsd）の経路ではありません。実際の経路は §4.0.2 のとおり（ccsd_iterdrv2.F → aoccsd2.F → ccsd.F → task_energy.F → task.F の errquit）。結論（NONZERO_EXIT）は同じで、コードの変更は要りません。
- **U8-N8（low、既知の残課題 7 の再掲）**: 入力の不正で止まったときも Rich のトレースバック（約 44 行）が先に出る（cli/main.py:139 は ValueError を捕まえない、終了コード 1）。回帰ではなく、メッセージの内容は S25 で十分になりました。

#### 4.8.7 追加の改良案

- **U8-R7（could、修正して採用。validation.md の部分は U0-R3′（should）と重なるので同じコミットで）文書の遅れをまとめて直す** 【S26】: (1) validation.md:62・:96・:139 を直す（C2 は U0 のプローブで確認済み、手法パネルは review のプローブ（v2: HCN 6 SP・36 s、HONO で panel_report と report の rank 1 が一致）でだけ確認済み）。(2) README:37、design.md:147、method_panel.yaml:5 の「全反応が thermo_unavailable」を「順位を付けられる反応は thermo_unavailable になる（same_basin などは outcome だけ）」に。(3) CCSD の未収束の根拠の訂正（U8-N7）は本レビューの注記で足りる。文書 5〜7 行、中立。
- **U8-R8（could）method_panel を v3 の run の複製で正式に 1 回実行して記録する** 【C43】: `cp -a /home/user/hfauto_v3/hono_known_endpoints` と hcn の複製に既定の method_panel を `--run-dir` で追記し（HCN では methods に ccsd-t_def2-tzvpd も足す）、panel_report/ranking.csv と report/ranking.csv の一致、method_panel.csv の 3 行、SP の数と時間を validation.md に記録する。元の v3 の run は変えない。コードは変えない。計算は数分。
- **C37（could、継続）使われなくなった monitor 経路と FailureKind.STAGNATED を削除** 【C37】: v2 の U8-R5 のとおり。timeout のときにプロセスグループごと kill する処理は残す（`_supervise` の timeout の部分だけを wait(timeout) にして残す）。約 30 行の削除とテスト 3 か所、【複雑さ減】。
- **C25（could、継続）ranking.csv に T_K・standard_state・dG_rxn_kcal** 【C25】: RankRow に 3 フィールドを足し、`_RANK_HEADER` と `_rank_cells` に列を加える。数行とテスト 1 件、中立。
- **C27′（could、修正して採用）使われない site.memory_mb を削除する** 【C27′】: 検査を足す代わりに、未使用の site.memory_mb（config.py:36）と wsl_local.yaml:3-4 の誤解を招くコメントを削除する（knob の削減）。検査を残すなら、NWChem の memory はプロセスごとの量なので `ranks × memory_mb_per_rank ≤ site.memory_mb` で比べる。【複雑さ減】（検査を残す場合は中立）。

#### 4.8.8 見送った案

- **S21 の trust や C2 の maxiter のような描画入力を JobStore のキーに入れる**: 収束した結果の妥当性は描画入力によらず、CCSD の未収束はもともと記録されない。キーを変えると既存のキャッシュがすべて無効になる。
- **CCSD の未収束を独自の FailureKind にし、job.t2 から続きを計算する**: パーサのパターン、失敗の種類、継続の規則が増える。maxiter 50 で余裕 32 反復。
- **conditions が一致しないとき、順位を付けない outcome にも thermo_unavailable を付ける**: S22 で消した same_basin の行の雑音が戻る。数値は誤らないので、文書の 1 句で足りる。
- **CLI で run_pipeline の ValueError を捕まえ、`_fail` で終える**: ValueError は stage 内のバグでも使われるので、バグのトレースバックまで消える。
- **S25 で添える失敗を欠けた型の artifact だけに絞る**: 他の型の失敗もたいてい同じ上流の原因で、絞っても情報は増えない。
- **主対象 TMA·(HF)₂ で既定パネルを実行して費用を測る**: 反応は same_basin の 1 行だけで順位の情報がなく、17 原子の TZVPD の SP 2 点に数十分かかる。小分子の正式な実行（C43）の方が情報が多い。
- **CCSD の thresh を明示する、maxiter をさらに上げる**: 既定の 1e-6（rms）でエネルギーの精度は 1e-8 Eh 台あり、50 反復で足りないことを示すデータはない。

#### 4.8.9 変更不要な点

- S22 の rankable の規則、C2（ccsd ブロックだけ、キーに入れない）、S25、S17（WFT は閉殻専用、開殻はジョブを作らずに INPUT_INVALID）。
- JobStore のキー = {engine, version_pin, kind, key_payload}。TERMINAL の 4 種だけを記録し、NONZERO_EXIT と TIMEOUT は再実行で計算し直す。
- sp stage（freq.final の構造と freq.level の電荷・多重度、Level と構造の echo 1e-4 Å の照合、calc_id で同じ構造を 1 artifact に）。
- M5 の ±1 kcal/mol の不感帯、S8 の既定パネル、CCSD(T) の 1 語での opt-in、band の推移的な重なりによる同順位、coverage.csv の機構別の内訳。
- Manifest.union で後の artifact を優先すること。実行基盤（キー単位の O_EXCL ロック、コア数セマフォ、SiteLock、/mnt 上の scratch の拒否、ladder）。

---

## 5. 今後の改良候補

有用性が確認できたもの（検証で「採用」「修正して採用」とされたもの）だけを載せました。工数は実装とテストの目安、「複雑さ」はコードと設定の量の増減です。どの改良も新しい knob は増やしません。

### 5.1 must

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|---|
| M9 | FIND_PATH の初期経路を逐次整列にする（`_initial_path` で `align_sequential([ends0, *inner, ends1])`）。同じコミットで validation.md:55・:137 と design.md:130 の「1 チャンク内のずれ」の説明を直す | U5-R8（修正して採用）、U5-N8・T1・D1 | 1 チャンク目の凍結端点のずれ STO-3G 0.074 → 0.004 Å、HONO 0.022 → 0.0014 Å、端点の持ち上がり +13.6 → +0.04、+1.09 → +0.009 kcal/mol。M8 と合わせて全チャンクで端点が極小になり、行 14 と形の判定の前提が 1 チャンク目でも成り立つ。追加の計算はない | 2 行の置き換え、テスト約 10 行、文書 2 か所 | 中立 | テストの判定は「隣接画像の直接の RMS ≈ mapped_rmsd」（`align_mapped(a, b) ≈ b` は直線の端点で縮退する）。直線と平面の IDPP で今のコードで落ちることを確かめてから入れる。string の JobStore キーが変わるので、再開時に 1 回再計算になることを validation に書く。新しい run dir で find_path_sto3g と HONO を再実行して記録する |

### 5.2 should

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|---|
| S26 | 検証記録と文書の訂正（第2ラウンド分）。内訳は下表 | U0-R3′、U3-R4（修正）、U7-R5′、U8-R7（修正）、U0-R7（修正、could）、U7-R6（修正、could）、検証段の指摘 | design.md の数値の出典と C2・S23 が動くことが文書から追える。S21 の仕組みの誤読、M6 の誤った帰属、S22 の後の conditions の文言のずれがなくなる。主対象の ΔG_assoc の LOT 依存が数値で分かる | 文書 15〜20 行 | 中立 | v2・v3 の run の記録そのものは書き換えない（観測値の引用には注記で対応）。v2 のプローブのパスはリポジトリの外なので、場所と要点の数値だけを書く |
| C38 | code_version の git status に `-c core.autocrlf=input` | U0-R6、U0-N7 | WSL の本番経路でも `-dirty` が本当の未コミット変更だけを示す（WSL 83 → 0 行、Windows 0 → 0 行） | 引数 2 つ | 中立 | テストは既存の code_version のテストで足りる |

S26 の内訳:

| 場所 | 内容 | 由来 | 優先度 |
|---|---|---|---|
| validation.md（§3 の後に新しい小節「レビューのプローブ」） | C2（`u0_c2`: デッキ `maxiter 50`、出力 `maxit = 50`、HNC 11/50・17 s、v3 HONO の TS 18/50・447 s、CCSD(T) の ΔE‡ 11.65、SVPD − CCSD(T) = +2.02）、composite（`hcn_tz`: ΔG‡ 41.99、GoodVibes 3/3 再利用、1:26。`tma_tz`: ΔG_assoc −9.04（TZVPD//SVPD、BSSE 未補正）。`tma_d`: 既定の targets では単量体が energy_layer_missing）、手法パネル（v2 のプローブ: HCN 6 SP・36 s、HONO で rank 1 が一致）。負荷時の時間は参考値と断る | U0-R3′、U7-R5′(c)、U8-R7(1) | should |
| validation.md:62・:96・:139 | 「C2 は未確認、S23 は文書だけ、M5・S8 は method_panel を実行していない」を「v3 の run では未実行。プローブで確認（新しい小節）」に | U0-R3′、U8-R7（:96 を追加） | should |
| validation.md:59（:124 は任意） | S21 の説明を実ジョブの記録に合わせる（v2 では歩幅の制限が 34 回・全体の縮小 7 回。歩数が同じなのは line search の補いと、勾配の収束後の XMAX / XRMS のための歩み。5% の短縮は同じ軌跡の run どうしの差（約 2%）より大きい観測で、原因は推定） | U3-R4（修正） | should |
| validation.md:76・:84 | 「M6・S15」→「S15（M6 の寄与は 0）」。:84 の表の注記も同じ帰属に | U7-R5′(a)、U7 の検証の指摘 | should |
| validation.md の表 #5（:18・:49） | ΔG_assoc −11.53 に「SVPD、BSSE 未補正」と注記 | U7-R5′(c)、U7-N11 | should |
| design.md:146（thermo 行） | vib_scale = 1.0 の根拠と文献値（0.985、差は ≤ 0.06 / 0.09）。ΔG_assoc の qRRHO の幅（298 K で約 1.2）。単原子の単量体・反応のない錯体では ΔG_assoc は出ない。単量体のアンサンブルには同じ (Hill, q, m) の異性体も入る | U7-R5′(b)、U7-R6（修正） | should |
| design.md:150 | 複製した run では thermo_tz が元の thermo を上書きする（元の report を再実行すると composite の値が出る） | U7 の検証の指摘 | should |
| README.md:37、design.md:147、method_panel.yaml:5 | 「全反応が thermo_unavailable」→「順位を付けられる反応は thermo_unavailable になる（same_basin などは outcome だけ）」 | U8-R7(2)、U8-N6 | could（同じコミットで） |
| design.md:149、README.md:86 | 「誤差の符号は反応で異なる（HCN −1.2、HONO +2.0。2 系での見積もり）。2 反応の ΔG‡ の差が約 3 kcal/mol 未満なら、順序は手法誤差で入れ替わり得る」 | U0-R7（修正） | could（同じコミットで） |

### 5.3 could

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| C8 | 組成の多重度が一意でなければ INPUT_INVALID | U2-R4 | ラジカル対の組成で高スピンの PES を黙って選ばない | 2〜3 行 | 中立 |
| C10 | window で screen 極小がなければ拒否、design.md:142 の include: all の文言 | U3-R2、U3-ISS-5、U3-N4 | 空の dft stage を「反応なし」と読み違えない（S25 では塞げない） | 数行、テスト 1 件 | 中立 |
| C12 | drive の重複を WL 原子クラスで除き、列挙数・固有数・実行数を diagnostics.json に | U4-R1、U4-I4 | 31 attempt 中 21 件が同じ移動の状態が解消し、heavy_bond・association の段に届く | 15〜25 行、テスト 2 件 | 中立 |
| C13 | 無バイアスの opt に `convergence_max_iterations=500` | U4-R2、U4-I5 | amine 系の数値的な陰性がなくなる（TMA 系では効かない） | 1 行 | 中立 |
| C14′ | TS がある negative にも ΔE‡（v2 の 41.5 は 44.05 に訂正） | U4-R3、U4-I6 | 陰性の結論を障壁の値（44.05 / 378.13 kJ/mol）で裏付けられる | 4〜8 行 | 中立 |
| C15 | tangent_mode_cm1 を barrierless のときだけ計算 | U5-R6、U5-I8 | 17 原子の縮退反応で 1 本約 900 s を省く | 約 5 行 | 中立 |
| C16 | 両端より低い DFT//xTB ノードがあれば unavailable（fake の E2E テストを同時に入れることが条件。M9 の後） | U5-R7（修正）、U5-I7 | SCREEN が中間体のある経路を BARRIERLESS で閉じることを防ぐ | 3 行、テスト 1 件 | 【増】 |
| C20 | 点群を ThermoResult.notes で記録し、S_rot を protocol から外す | U7-R4、U7-I7、U7-N6 | σ を監査でき、使われない経路が 2 つ片づく | 約 4 行 | 中立 |
| C25 | ranking.csv に T_K・standard_state・dG_rxn_kcal | U8-R3、U8-7、U8-N5 | CSV が単独で読め、metric=dG_rxn で並べた値も出る | 数行、テスト 1 件 | 中立 |
| C27′ | 未使用の site.memory_mb と wsl_local.yaml:3-4 の誤解を招くコメントを削除（検査を残すなら ranks × memory_mb_per_rank で比べる） | U8-R4（修正）、U8-9 | 効かない knob がなくなる | 数行 | 【減】 |
| C28 | layer_missing の subject を含む ΔE を None にする | U0-R4 = U7-R1、U0-N1 | 2 つの LOT を混ぜた ΔE が表示されない | 1 行、assert 1 つ | 中立 |
| C29 | ΔG_assoc=None の理由を ReactionThermo.notes に、曖昧な対応は None | U7-R2、U7-I2、U7-N2・N3 | 無言の None（`tma_d` で再現）と原因を区別できる | 約 10 行、テスト 1 件 | 【増】 |
| C30 | 設定値の範囲検査（ConformersConfig、Conditions、ThermoSettings） | U2-R3 + U7-R3、U2-N4、U7-N4 | 設定の誤りを読み込み時に止める | 約 8 行、テスト 2 件 | 中立 |
| C31 | 読み込み時に xyz の存在と宣言反応の両端の q/m の一致を検査 | U1-R3、U1-N2・N4 | `running` の残留と dry-run の素通り、片側の typo による無駄な DFT を防ぐ | 約 5 行、テスト 2 件 | 中立 |
| C32 | rerank_sp の SP を per_state を超える組だけに出し、失敗を diagnostics に | U3-R1、U3-N1 | 単量体の最低構造が SP の失敗で消えない。TMA で約 50 s 減 | 4〜5 行、テスト 1 件 | 中立 |
| C33 | reaction-paths の `_register` に pipeline の gates を渡す | U3-R3、U3-N2 | 1 つの PES の tier の基準が 1 つにそろう | 1 行 | 中立 |
| C34 | 単量体の停止後の再実行を「停止構造から既定の設定で 1 回」にし、`--noreftopo` を廃止 | U2-R1、U2-N1 | 実 CREST で 2/2 失敗していた経路が 2/2 成功（v2 のプローブ）。フィールドと CLI フラグが 1 つずつ減る | 約 10 行、fake テスト 1 件 | 【減】 |
| C36 | S4 の再実行を GS が完了した場合に限る | U5-R4、U5-N4 | GS 自体が落ちたとき同じ GS を繰り返さない | 1〜2 行 | 中立 |
| C37 | 使われなくなった monitor 経路と FailureKind.STAGNATED を削除（timeout の kill は残す） | U8-R5、U8-N2 | 到達しないコード約 30 行と失敗の種類 1 つが減る | 約 30 行の削除、テスト 3 か所 | 【減】 |
| C39′ | 文書の訂正の残り（下表） | U1-R4、U2-R2、U3-R5、U4 の C39 分、U5-R9、U6-R5 | 実装・実測と文書の食い違いをなくす | 文書 10〜15 行と docstring 1 行 | 中立 |
| C40 | S24 のメッセージを `hfauto uses most-abundant-isotope masses` に（validation.md:61 の観測記録は書き換えず注記） | U1-R5（修正）、U1-N8 | 設計書・NWChem・GoodVibes の実際の質量と文言がそろい、`[1H]` の拒否が自己矛盾に見えない | 1 行 | 中立 |
| C41 | `_Scine.task` で捨てている例外の文字列を stderr に 1 行（C14′ と同じ変更で） | U4-R4、U4-N4 | TSOpt・NT2 の失敗原因をプローブの再実行なしに生ファイルから確かめられる | 1〜2 行 | 中立 |
| C42 | smoke の moddir テストで、NWChem が追ったモードを P·H·P の第 2 固有ベクトルと照合（|cos| > 0.99、helper はテスト内） | U6-R6、U6-N11 | 実 pipeline で未通過の S20 の分岐の前提を回帰テストで守る | smoke に約 8 行 | 中立 |
| C43 | method_panel を v3 の run の複製（HONO、HCN + CCSD(T)）で正式に実行し、validation.md に記録 | U8-R8、U8-N4 | sp stage → panel_report の経路が正式な検証の記録として残る | 計算数分、文書数行 | 中立 |
| C44 | composite のとき report の見出しに energy 層の LOT を出す（freq と違うときだけ） | U7-R7、U7-N10 | S23 の手順で作った report を SVPD の report と取り違えない | 約 5 行、assert 1 行 | 中立 |

C39′ で直す文書の内訳:

| 場所 | 内容 | 由来 |
|---|---|---|
| design.md:141（explore 行） | 状態ラベルは結合グラフ（1-WL、立体と E/Z を区別しない、1.15〜1.45 Σr_cov の中間帯は結合）なので、配座・立体・E/Z だけ違う IRC の端は same_as_source として捨てる | U1-R4、U4 の C39 分 |
| validation.md:91 | −289i の TS の判定の変化は M3 の効果で、縮退した HF 交換（ΔE‡ 44.05 kJ/mol） | U4 の C39 分、U4-N1 |
| design.md:139（conformers 行）、environment.md:64、conformer_search.py:9-10 | 開殻の組成は `--noopt` で trial MTD が失敗することがあり、そのときは seed が出る。「約 25 s」→ 約 1.4 s、crest の threads の既定は 1。docstring を 4 knob に | U2-R2、U2-N2・N3・N5・N6 |
| design.md §8（ladder の行） | 初期 Hessian ありの opt の継続は drv.hess を INHESS 0 の restart で引き継ぎ、trust 0.1 で続ける | U3-R5、U3-N8 |
| design.md:130、design.md の reaction-paths 行（M9 と同じコミットで） | NWChem の string は gmax を記録しない（pysis_gs だけ）。impose は xyz_path を渡すと使われず、実際の整列は zts_min_motion（凍結 bead N も含む） | U5-R9、U5-N6・N9 |
| design.md の reaction-paths 行、validation.md:69 | autoz / Cartesian の分岐は ν < 0 の本数で決める。NWChem の smalleig で柔らかい負モードは追えないことがある。子反応はそれぞれ 6 h の予算（最悪 7 case）、qrc_bounds_A の上限は 0.5 Å 以下。「144」→「143」 | U6-R5（修正）、U6-N6・N7・N8・N9 |

v2 の could のうち今回再提案しなかったもの: C1（見送り、§4.0.8）、C4 のラジカルの部分（見送り）、C18（実 run で到達例がなく、v2 の判断のまま据え置き）、C21（数値の文書化 S26 に置き換え）、C22（方針どおり見送り）、C35（見送り、§4.5.8）。

### 5.4 進め方と確認

1. **M9（+ C39′ の U5・U6 分の文書）**: actions.py の 2 行と、直線・平面の IDPP のテスト。確認は、新しい run dir での find_path_sto3g と HONO の再実行（1 チャンク目の end_drift が 0.001〜0.004 Å 程度になり、形の判定と TS が変わらないこと）。string のキーが変わることを validation に明記する。C16 を入れるならこの後。
2. **S26 と C39′ の残り（文書）**: まとめて 1 コミットで行えます。v2・v3 の run の記録は書き換えず、訂正は注記で行います。
3. **小さなコード**: C38、C28、C40、C41（+ C14′）、C25。どれも数行で、報告・診断・来歴の正確さを上げます。
4. **複雑さを減らす**: C34、C37、C27′。C34 は実エンジンのプローブ（約 2 分）で確認します。
5. **検証の空白を埋める最小の実計算**（任意）: C43（手法パネルの正式な実行）と C42（moddir の smoke の検査）。M9 の後、monotonic になる軽い系で行 14 を 1 回通す。

### 5.5 あえて行わないこと

- **既定の流れを重くする変更**: 停留点レベルを TZVPD や ωB97X-D3 に上げる、TZVPD の composite を discover.yaml の既定に組み込む、CCSD(T) を既定のパネルに入れる、主対象 TMA·(HF)₂ で既定パネルを回す、composite.yaml の同梱（C1）。効果に対して費用（17 原子の TZVPD の SP 数分〜数十分、CCSD(T) 1 点数時間）が見合わず、validation の値との連続性も崩れます。
- **キーとスキーマを広げる変更**: trust・maxiter などデッキの詳細を JobStore のキーに入れる、thermo の artifact id に stage id を入れる、ReactionThermo に energy_level のフィールドを足す、失敗した stage の例外文を run_state に保存する、TS 候補の ν1 や陰性理由の細分化を result に足す。
- **knob や閾値を増やす変更**: CCSD の maxiter・thresh の knob 化、手法誤差を同順位の幅に組み込む、end_drift のゲート、_SETTLED の緩和、パネルの conditions の一致検査、錯体の opt の収束判定の緩和、ZTS の小さなチャンク。
- **根拠のない性能改善・再試行**: 全 optimize の trust 0.3、ladder の継続での trust 0.3 の維持、QRC の trust のさらなる引き上げ、mode-follow・soft 押しへの Hessian、柔らかい負モードへの vardir や Hessian の加工、Bofill が失敗したときの別の optimizer での TSOpt、IRC の反復上限の変更、C35（pysis の off_plane）。
- **表示だけのための配管**: CLI で ValueError を捕まえてトレースバックを隠す、missing_component の reject に上流の理由を連結する、S25 の件数の表示、placement の seed の座標の正規化、`_moddir` の本数の数え方をコードで直す（文書で足りる）。
- **対象外の化学への対応**: 開殻 WFT、単原子種の全面対応、SMILES のラジカル電子数の照合、主同位体と同じ質量数の同位体指定の許可、状態ラベルへの立体の追加。

---

## 6. 参考文献・仕様 URL

v1 と v2 の一覧（Bursch/Grimme 2022、Rappoport/Furche 2010、Pracht/Bohle/Grimme 2020、Smidstrup ら 2014 の IDPP、E/Ren/Vanden-Eijnden 2007 の ZTS、Goodman & Silva 2003 の QRC、Grimme 2012 の qRRHO、Luchini ら 2020 の GoodVibes、Kesharwani ら 2015 のスケール因子、各ライブラリの文書とソース）は引き続き有効です。今回あらたに参照したもの、または再確認したものを挙げます。

**文献**

- Kesharwani, Brauer, Martin, J. Phys. Chem. A 2015, 119, 1701 — https://webhome.weizmann.ac.il/home/comartin/OAreprints/260.pdf （PBE0/def2-SVPD の ZPVE 因子 0.9848。U7 の vib_scale の根拠）
- Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. 140, 214106 (2014)（IDPP。目的関数を収束させた経路が前提。S18）

**ライブラリの仕様・ソース**

- NWChem Geometry Optimization（TRUST の既定、INHESS、MODDIR、FIRSTNEG）— https://nwchemgit.github.io/Geometry-Optimization.html
- NWChem CCSD（MAXITER 20、THRESH 1e-6、RHF 専用）— https://nwchemgit.github.io/CCSD.html
- NWChem NEB / ZTS（FREEZE1、FREEZEN、IMPOSE、XYZ_PATH）— https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html
- NWChem Memory（memory はプロセスごと）— https://nwchemgit.github.io/Memory.html
- NWChem 7.2.3 のソース:
  - driver（trust の既定 1069-1074、歩幅の制限 1743-1787・1874・1905-1931、Cartesian の射影子 2760-2850・2948-2963、対角化 1138-1139、moddir 3270-3274、ofirstneg 1112-1113・3299-3305、smalleig 3135・3210-3240）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/driver/opt_drv.F
  - ZTS（impose 207-248、zts_min_motion 569-571・1967・2007、凍結 bead のエネルギー 1830-1831）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F
  - ZTS の入力（xyz_path と bead_list:new 148-173）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string_input.F
  - CCSD（ccsd_input.F の maxiter、ccsd_iterdrv2.F の最大反復、aoccsd2.F の未収束時の分岐、ccsd.F の戻り値と参照 SCF の閾値 34-42）— https://github.com/nwchemgit/nwchem/tree/v7.2.3-release/src/ccsd
  - task（task_energy.F、task.F の errquit）— https://github.com/nwchemgit/nwchem/tree/v7.2.3-release/src/task
- CREST 3.0.2 の引数解析（`-nci`、`-quick`、`-ewin`、`-noopt`、`-notopo`、`-noreftopo`）— https://raw.githubusercontent.com/crest-lab/crest/v3.0.2/src/confparse.f90
- SCINE ReaDuct 6.1.0 TsOptimizationTask（例外の扱い、`convergence_max_iterations`）— https://raw.githubusercontent.com/qcscine/readuct/6.1.0/src/Readuct/App/Tasks/TsOptimizationTask.h
- RDKit rdmolops（RemoveHsParameters.removeIsotopes）— https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html
- GoodVibes 4.3.0（WSL prod venv の `goodvibes/thermo.py` L697-747 の calc_bbe、`goodvibes/api.py` の compute_thermo のシグネチャ、`goodvibes/io.py` L223-225 の原子質量）
- ruff の mccabe（C901、max-complexity）— https://docs.astral.sh/ruff/rules/complex-structure/

**リポジトリ内の関連文書**

- `docs/reviews/2026-09-26_calculation_spec_review_v2.md`（ベースライン）、`docs/reviews/2026-09-26_calculation_spec_review.md`（v1）
- `docs/design.md`、`docs/validation.md`（v3）、`docs/environment.md`、`README.md`

**実 run とプローブ**（WSL）

- 改良前: `/home/user/hfauto_v2/<run>`（コード `4e972c5`）。改良後: `/home/user/hfauto_v3/<run>`（コード `dc82145`）。参考: `/home/user/hfauto_w7/<run>`。
- 本レビューのプローブ: `/home/user/hfauto_v3_probe/{u0_c2, s24_s25, u1_review, u1_assess, u2_s24, u3_review, u4_spec, u4_review, u5_spec, u5_assess, u7_spec, u7_assess}`。第2ラウンドの実装時のプローブ: `/home/user/hfauto_v3_probe/w1b_moddir`。v2 のプローブ: `/home/user/hfauto_review2_probe/`。
