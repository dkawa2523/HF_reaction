# hfauto 計算処理の仕様 再評価（v4）— 第3ラウンドの改良（M9、S26、C38）の後の評価・残る課題・今後の候補（2026-09-26）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `de9daf8`。`hfauto/` のコードは `00680e0` と同じです（`de9daf8` は docs/validation.md と docs/design.md だけを変更）。読み取り専用のレビューで、リポジトリのコード・設定・既存の文書は変更していません（追加したのはこのファイルだけ）。
- ベースライン: `docs/reviews/2026-09-26_calculation_spec_review_v3.md`（以下「v3」）。その §1 の評価、§3 の課題の状態、§4 の単位ごとの課題と改良案、§5 のロードマップ（must M9、should S26・C38、could C8〜C44）を基準にしました。それより前の評価は v2（`..._review_v2.md`）と v1（`..._review.md`）です。
- 評価した改良:

| コミット | 内容 |
|---|---|
| `00680e0`（R3） | M9（FIND_PATH の初期経路の逐次整列）、C38（code_version の git status に `-c core.autocrlf=input`）、S26（検証記録と文書の訂正） |
| `de9daf8`（R3-validate） | docs/validation.md §11（WSL での実計算による v4 の再検証）と、§3・§9・§10 から §11 への参照。design.md:130 の端点のずれの目安を「0.004 Å 以下」に |

- 方法: v3 と同じ 3 段です。今回は、第3ラウンドの改良が通る 4 単位（U0、U5、U6、U8）を再評価しました。残りの 5 単位（U1〜U4、U7）は、差分と実 run の照合で挙動が変わっていないことを確かめ、v3 の評価を引き継ぎました（§4.9）。
  1. **現状仕様の整理**: コード、設定と、改良前後の実 run の生ファイル（job.nw、stdout、result.json、manifest.json、resolved_config.yaml、ranking.csv、report.html）を読み取りだけで比べました。改良前は v3（コード `dc82145`、WSL `/home/user/hfauto_v3/<run>`）、改良後は v4（コード `00680e0`、`/home/user/hfauto_v4/<run>`）です。
  2. **評価**: 化学的妥当性、有用性、複雑さを評価しました。NWChem 7.2.3 の string のソース（string.F、string_input.F）と git の改行変換（convert.c）を調べ、小さなプローブを行いました（WSL `/home/user/hfauto_v4_probe/m9`、`/home/user/hfauto_v4_probe/u8_crlf`）。
  3. **反証的な検証**: 別の担当が 2 の主張をコードと生データで確かめ、状態・重大度・優先度・評価を確定しました。
- 記載の規則:
  - 課題の状態、重大度、優先度、評価は、検証段の判定を使いました。
  - 全面的に棄却された課題はありません。一部だけ正しいと判定された課題（U5-N10）は、訂正後の内容で書きました。
  - 改良案は、検証で「採用」または「修正して採用」とされたものだけを載せました。修正して採用されたもの（U5-R10）は修正後の形で書きました。評価段で見送った案は、各節の「見送った案」と §5.5 に理由付きで載せました。
  - ID の付け方: v3 までの課題と改良は元の ID のままです。今回見つかった課題は `U<n>-N<k>` の続き番号（U0-N9、U5-N10 など）、単位ごとの改良案は `U<n>-R<k>` の続き番号です。§5 では、今回見つかった検証記録の訂正を新しい **C45** にまとめました。v3 の文書訂正 **C39′** は、U5 分の範囲を縮めて引き継ぎました。
  - 対象単位の節番号は v3 と同じにしました（§4.0、§4.5、§4.6、§4.8）。

凡例:

| 記号 | 意味 |
|---|---|
| critical / high / medium / low | 課題の重大度 |
| must / should / could | 優先度。must は「結論を誤らせる」または「主要機能が働かない」もの |
| 【複雑さ減】/【複雑さ増】 | コード・設定・分岐が減る / 増える改良 |
| ◎ / ○ / △ / × | 評価。◎ 問題なし、○ 軽微な課題のみ、△ 改良が必要、× 目的どおりに機能していない。「○−」は ○ の下限 |
| resolved / partially / unresolved / not_in_scope | 課題の状態（解消 / 一部解消 / 未解消 / 本単位の範囲外・方針により対応しない） |

エネルギーの単位は、断りがなければ kcal/mol です。「v3 → v4」は第3ラウンドの改良前 → 改良後の実 run を指します。

---

## 0. 要旨

### 0.1 全体評価

1. **第3ラウンドの改良は、v3 の指定どおりに入りました。** must 1 件（M9）と should 2 件（S26、C38）がすべて実装され、WSL の実計算（v4）で効果を確かめました。
   - `hfauto/` の差分は 2 ファイル・正味 +1 行です（actions.py の 2 行の置き換えと、config.py の引数 2 つ）。
   - knob、設定項目、失敗の種類、出力列は増えていません。
   - pytest 293 件と smoke 12 件（71 s）が通り、ruff・pyrefly・lint-imports も問題ありません。
2. **評価の記号は、9 単位とも v3 と同じです（§1）。** 下がった単位はなく、△ と × もありません。変わったのは留保の中身です。
   - U5 FIND_PATH: 「○ の下限」が外れて ○ になりました（M9 で U5-N8 が解消）。
   - U0 有用性 ◎: v3 で付けた 4 つの条件のうち 3 つが解消しました（S26 と C38）。残る条件は、method_panel を検証の run で実行していないことだけです。
   - U6 ◎ / ◎ / ○ と U8 ◎ / ○（上側）/ ○ は据え置きです。
   - U1〜U4 と U7 は、コードと設定に差分がなく、実 run でも挙動が同じなので、v3 の評価を引き継ぎました。
3. **実計算での効果**（validation.md §11 と本レビューの再計算）
   - **M9**: 1 チャンク目の凍結端点のずれは、STO-3G HCN で 0.0737 → 0.0008 Å、HONO で 0.0220 → 0.0016 Å になりました。その端点の実エネルギーの持ち上がりは +13.57 → +0.002、+1.09 → +0.012 kcal/mol です。STO-3G は v3 の予測（0.004 Å、+0.04）より小さく、HONO は予測（0.0014 Å、+0.009）と同程度でした。初期経路の画像の間の剛体の跳び（0.27 / 0.50 Å）は 1e-13 Å 未満になりました。
   - **結論は変わりません。** チャンク数、経路の形、HEI の位置、鞍点、接続の判定は v3 と同じです。ΔG の差は 0.001 kcal/mol 以下で、順位の付く 4 反応はすべて rank 1 のままです。
   - **C38**: v4 の 5 run の code_version はすべて `00680e0fa525…` で、`-dirty` は付きませんでした（v3 は 7 run すべてに付いていた）。WSL の status の出力は 83 → 0 行です。
   - **S26**: v3 の S26 の内訳 9 行がすべて入り、記述がコードと実測に合うことを確かめました。S24 のメッセージの文言だけは、v3 の割り当てどおり C40 に残っています。
4. **新たに見つかった課題は、すべて low です。** 大半は検証記録の書き方の不正確さです。
   - validation.md:147 の 2 か所: 「HEI は bead の最大」は実装（放物線補間）と違います（U6-N12）。「M9 より前の run dir で再開すると string を 1 回計算し直す」は、paths を再実行するときにしか成り立たず、そのときは U6 のジョブも計算し直します（U6-N13 = U8-N10）。
   - validation.md:170: HONO の string が 12.5 s 遅くなった原因を「特定していない」としていますが、stdout から分かります。NWChem が damped Verlet に切り替えた回数の違いで、M9 の系統的なコストではありません（U5-N10）。
   - validation.md:72: design.md:149 の「HCN −1.2」の出典である CCSD(T) 47.81 が、methods に opt-in で足した別の実行だったことが書かれていません（U0-N9）。
   - code_version は再開のたびに上書きされるので、版をまたいで再開した run の来歴は正確ではありません（U8-N9。以前からの性質で、C38 で値が信頼できるようになったので意味を持つようになった）。
   - README のベンチマーク表は v3 の値のままです。「v3」と明記しているので誤りではなく、対応は不要です（U0-N10 = U8-N11）。
5. **残る課題**: medium は U5-N3（DLC の GS の脆さ。影響はコストだけ）の 1 件で、ほかはすべて low です。実計算で通っていない経路（§2.8）は v3 と同じです。
6. **今後**: must と should に当たる改良は残っていません。could のうち、文書だけで済むもの（C45、C39′）と、評価を上げる条件になっている小さな改良（C25・C43、C37・C27′、C28）が主な候補です。同じ形のラウンド（実装 → 実計算での再検証 → 再評価）をもう一度回す必要はありません（§5.4）。

### 0.2 今後の候補（上位 5 件）

どれも could です。

| # | 改良 | 由来 | 期待効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| 1 | **C45 + C39′** 検証記録と文書の訂正（validation.md:72・:79・:147・:170、design.md:130・:144・:158 など） | U0-R8、U5-R10（修正）、U6-R7、U8-R9、U5-R9′、C39′ の U6 分 | 検証記録と design.md が実装（profile.hei、runner._fresh、NWChem のソース）と一致する。再開の挙動と code_version の意味を誤読しない。「原因不明の +7.7%」がなくなり、M9 に性能劣化の疑いが残らない | 文書 10〜15 行 | 中立 |
| 2 | **C12 + C14′ + C41** explore の被覆と診断（v3 のまま） | U4-R1、U4-R3、U4-R4 | U4 の有用性を下限に留めている「31 attempt 中 21 件が同じ移動」が解消し、陰性の結論を障壁の値で裏付けられる | 20〜35 行、テスト 2〜3 件 | 中立 |
| 3 | **C25 + C43** ranking.csv の列と、手法パネルの正式な実行 | U8-R3、U8-R8 | ranking.csv が単独で読め、「panel_report が最終」という前提が正式な検証の記録に入る。U8 の有用性を ◎ に上げる条件 | 数行とテスト 1 件、計算数分 | 中立 |
| 4 | **C37 + C27′ + C34** 使われない分岐と knob の削除 | U8-R5、U8-R4、U2-R1 | 到達しないコード約 30 行、失敗の種類 1 つ、効かない knob 1 つが減る。実 CREST で失敗していた単量体の再実行が成功する。U8 と U2 の複雑さを ◎ に上げる条件 | 約 45 行（大半は削除） | 【複雑さ減】 |
| 5 | **C28** layer_missing の subject を含む ΔE を None にする | U0-R4 = U7-R1 | 2 つの LOT を混ぜた ΔE（HCN で 164 kcal/mol）が report に出なくなる。U0 の妥当性を ○ に留めているコード側の唯一の理由を消す | 1 行、assert 1 つ | 中立 |

### 0.3 第3ラウンドで良くなり、今のままでよい点

- **M9 の形**: 初期経路、チャンクの継ぎ目、最終経路の 3 か所が、同じ `align_sequential` の規則（ends[0] を固定し、ends[1] を含む各画像を直前の画像へ剛体整列）になりました。ends[1] は剛体として回るだけです。デッキの endgeom は、xyz_path を渡すと NWChem が読まないので、回った ends[1] を渡しても問題ありません。
- **end_drift の扱い**: note で監視するだけで、ゲートにしない（M9 の後は 0.0008〜0.0016 Å）。
- **C38 の形**: git の引数 2 つだけで、.gitattributes も専用テストも足さない。index に CRLF の blob はなく、改行以外の内容の変更は引き続き `-dirty` になる。
- **code_version の位置づけ**: config_sha にも JobStore のキーにも入れない。再開と再利用の価値を守る。
- **順位差の目安**: 「約 3 kcal/mol」は文書の読み方にとどめ、コードの閾値にしない。

---

## 1. 評価の推移表

評価は「妥当性 / 有用性 / 複雑さ」の順です。複雑さの ◎ は「単純で過不足がない」を意味します。「継承」は、本ラウンドで再評価せず v3 の評価を引き継いだ単位です（根拠は §4.9）。

| 単位 | v1 | v2 | v3 | v4 | 主な変化（v3 → v4） |
|---|---|---|---|---|---|
| U0 全体プロトコル | ○ / ○ / △ | ○ / ○ / ○ | ○ / ◎ / ○ | ○ / ◎ / ○ | C38（v4 の 5 run すべて -dirty なし）、S26（プローブの記録、手法誤差の符号と約 3 kcal/mol の目安、conditions の文言）。有用性 ◎ の 4 条件のうち 3 つが解消 |
| U1 structures | ○ / ◎ / ◎ | ○ / ◎ / ◎ | ◎ / ◎ / ◎ | ◎ / ◎ / ◎（継承） | 差分なし。5 run の structures stage の payload が v3 と完全一致 |
| U2 conformers | △ / ○ / ○ | ○ / ○ / ◎ | ○ / ○ / ◎ | ○ / ○ / ◎（継承） | 差分なし。配座探索が通る run は再実行していない |
| U3 minima | ○ / ○ / △ | ◎ / ◎ / ○ | ◎ / ◎ / ○ | ◎ / ◎ / ○（継承） | 差分なし。5 run の dft stage は artifact id と job_key が v3 と同じ |
| U4 explore | △ / △ / ○ | ○ / ○ / ○（有用性は下限） | ○ / ○ / ○（有用性は下限） | ○ / ○ / ○（有用性は下限、継承） | 差分なし。探索が通る run は再実行していない |
| U5 paths（仮説・SCREEN・FIND_PATH） | SCREEN ○ / FIND_PATH × / △ | ○（FIND_PATH ○−）/ ○ / ○ | ○（SCREEN ○、FIND_PATH ○ の下限）/ ○ / ○ | ○（SCREEN ○、FIND_PATH **○**）/ ○ / ○ | M9 で U5-N8 が解消（1 チャンク目の端点のずれ 0.074 → 0.0008 Å、持ち上がり +13.6 → +0.002）。FIND_PATH の「下限」が外れた |
| U6 paths（鞍点・検証・接続） | ○ / ◎ / ○ | ○ / ◎ / ○ | ◎ / ◎ / ○ | ◎ / ◎ / ○ | コード差分なし。M9 で FIND_PATH の種が約 0.004 Å 動いたが、TS は同じ（≤ 1.3e-4 Å）。U6 区間 181 → 178 s |
| U7 thermo | ○ / ○ / △ | ◎ / ○ / ○ | ◎ / ○ / ○ | ◎ / ○ / ○（継承） | コード差分なし。S26 の文書（vib_scale の根拠、ΔG_assoc の LOT の注記など）。極小の ΔG の差 ≤ 5.1e-6 |
| U8 sp・report・実行基盤 | ○ / ○ / ○ | ◎ / ○ / ○ | ◎ / ○ / ○（有用性は ○ の上側） | ◎ / ○ / ○（有用性は ○ の上側） | C38、S26（conditions の文言、プローブの記録）。ranking.csv と report.html に回帰なし |

評価の留保（検証で付いたもの）:

- **U0 妥当性 ○**: U0-N1（層が欠けた subject の ΔE の混成。C28 は未実施）と、電荷のある系・開殻系・method_panel を実 run で通していないことが理由です。U0-N8 は解消しました。
- **U0 有用性 ◎**: 残る条件は、method_panel を検証の run で実行していないことだけです。これは v3 §4.0.8 が「プローブの記録で足りる」とした範囲で、今回その記録が validation.md §3 に入りました。
- **U1 妥当性 ◎**: v3 と同じく U1 単体での評価です（電荷のある系・開殻系の実 run は 0 本、C8 は U2 に残る）。
- **U2 複雑さ ◎**: v3 と同じ但し書きです（C34、C8、C39′ の U2 分）。
- **U3 妥当性 ◎**: v3 と同じです（DFT の mode-follow・soft 押し・M6 の |ν| 化は実エンジンで未発火）。
- **U4 有用性 ○**: 下限のままです（C12・C14′）。
- **U5 妥当性 ○**: ◎ にしない理由は 3 つです。U5-I7（両端より低い DFT//xTB ノード。C16 は未実施）が残ること、行 14（monotonic → BARRIERLESS）と multi_max が実 run で一度も通っていないこと、GS や SCREEN の IDPP を初期経路にする FIND_PATH も実 run で通っていないこと（align_sequential は座標系に依存しないので、危険は小さい）。
- **U6 妥当性 ◎**: v3 と同じです。実計算で通ったのは 3〜4 原子の共有結合系だけで、負モード 2 本以上の分岐はプローブだけ、柔らかい負モードは NWChem の smalleig で追えないことがあります（U6-N9）。FIND_PATH 由来で 10 原子を超える TS も、まだ実 run で通っていません。
- **U8 有用性 ○（上側）**: C25（ranking.csv の列）と C43（手法パネルの正式な実行）が残ります。この 2 つで ◎ にできます。複雑さ ○ は、C37 と C27′ で ◎ にする余地があります。

---

## 2. 第3ラウンドの改良の実装と実計算での確認

### 2.1 改良ごとの実装と確認

| ID | 内容 | 実装箇所 | 実計算・プローブでの確認 | 判定 |
|---|---|---|---|---|
| M9 | `_initial_path` で中間画像を ends[0] へ 1 枚ずつ整列するのをやめ、`align_sequential([ends[0], *inner, ends[1]])` にする | `hfauto/drivers/reaction_case/actions.py:455-456`（2 行の置き換え。ファイルは 500 行のまま、import も不変） | 1 チャンク目の端点のずれ: STO-3G 0.0737 → 0.0008 Å、HONO 0.0220 → 0.0016 Å。持ち上がり: +13.57 → +0.002、+1.09 → +0.012。テスト `test_first_string_chunk_starts_from_an_idpp_without_rigid_jumps`（直線 LINEAR・平面 PLANAR）は現行コードで PASS。M9 前の `_initial_path` にメモリ上で戻すと両方 FAIL する（隣の組 8 組中 7 組が不一致）ことを、本レビューで再現した | 効果あり |
| C38 | code_version の git status に `-c core.autocrlf=input` | `hfauto/pipeline/config.py:162-163`（100 字制限で 2 行に折り返し） | v4 の 5 run は `00680e0fa52527fa1e18ec72a33db492cdb14a0c`（-dirty なし）、v3 の 7 run は `dc82145…-dirty`。WSL の status は 83 → 0 行。code_version() は WSL・Windows とも `de9daf8…`（-dirty なし） | 効果あり。副作用なし（§2.4） |
| S26 | 検証記録と文書の訂正（v3 §5.2 の内訳 9 行） | validation.md:20・:51・:57・:61・:64・:66-72・:86・:88・:106・:147・:149、design.md:130・:146・:147・:149・:150、README.md:37・:86、method_panel.yaml:5-6 | 文書だけ。各記述をコードと v2・v3 レビューの数値に照合し、すべて一致した（§2.5）。プローブの置き場所 5 か所が WSL 上に実在する | 実施済み。新しい不正確さ（low）は C45 |

指定との差（00680e0 のコミットメッセージに記載）:

- 新しいテストのために、DriftingPath が Failure をそのまま返すようにした（2 行）。
- config.py の呼び出しを、100 字制限のため 2 行に折り返した。
- S24 のメッセージの「natural-abundance」は、v3 の割り当てどおり C40 に残した。

README.md:86 の見出し「第2ラウンド改良後の再検証 v3」とベンチマーク表の値は v3 のままです。コミットメッセージには書かれていませんが、S26 の指定に含まれていないので逸脱ではありません（U0-N10）。

### 2.2 run ごとの比較（v3 → v4）

どれも新しい run dir（`/home/user/hfauto_v4/<run>`）での実行です。

| 系 | 所要時間 | ジョブ数（再利用） | 結果（v3 との差） | ΔG‡ |
|---|---|---|---|---|
| HCN→HNC（known_endpoints） | 1:24 → 1:24 | 29(2) → 29(2) | `elementary_step`、TS −1128.5i、鞍点 4 ステップ、QRC 12 / 15。すべて v3 と同じ | 42.1783 → 42.1783 |
| HONO trans→cis | 5:23 → 5:32 | 16(0) → 16(0) | GS は v3 と同じく失敗し、FIND_PATH の 1 チャンク（20 反復）で `single_max`（bead の最大 13.47 で同じ）。鞍点 3 ステップ、TS −681.0i → −681.5i（W7 は −681.5i）、QRC 15 / 15 → 15 / 13 | 12.2256 → 12.2267（+0.001） |
| NH3 反転 | 1:25 → 1:25 | 27(1) → 27(1) | `degenerate_rearrangement`、TS −757.7i、鞍点 2 ステップ、QRC 25 / 25。すべて v3 と同じ | 3.8395 → 3.8395 |
| 水（same_basin） | 0:12 → 0:11 | 4 → 4 | `same_basin`、blocker は `outcome:same_basin` だけ | — |
| STO-3G HCN の FIND_PATH | 1:17 → 1:19 | 15(0) → 15(0) | string 3 チャンクで `single_max`。HEI から鞍点 −1223.9i・3 ステップ、`elementary_step`。QRC 14 / 15 → 14 / 14 | 62.1678 → 62.1677 |
| 実エンジンの smoke | 71 s → 71 s | 12 → 12 件 | 12 passed | — |

- ΔE‡ と ΔE_rxn は v3 と 1e-5 kcal/mol 以内で一致しました。ΔG の差は 0.001 kcal/mol 以下で、順位の付く 4 反応はすべて rank 1 のままです。同順位の幅は、どの反応も 0.001 kcal/mol 未満のままでした。
- resolved_config.yaml は、code_version を除くと v3 と同じです。STO-3G だけは一時 xyz のパスが 2 か所違います（`hfauto_v3` → `hfauto_v4`）。
- JobStore の再利用と失敗の内訳（HONO の GS の nonzero_exit 1 件だけ）は v3 と同じです。
- report.html のバイト数は、HCN・HONO・NH3・水で v3 と同じでした。STO-3G だけ 5181 → 5184 B で、TS のモードの座標の負号が 3 つ増えた分です。STO-3G では QRC の表示順も入れ替わりましたが、固有ベクトルの符号が反転しただけで、接続の判定は同じです。
- TMA·(HF)₂、A/B、amine パイロット、(HF)₃ は、第3ラウンドの改良が通らないので再実行していません。

### 2.3 M9 の実測

1 チャンク目の凍結端点のずれ（string_final と初期経路の端点の写像 RMSD）と、その端点の実エネルギーの持ち上がり（string と同じデッキでの DFT SP の差）を、v3 と v4 の run に同じスクリプト（`/home/user/hfauto_v4_probe/m9/m9_ends.py`）を当てて測りました。v3 の run での値（+13.569、+1.092）は、v3 レビューのプローブ（`u5_assess`）と一致します。

| 量 | STO-3G HCN（HNC 側） | HONO（cis 側） |
|---|---|---|
| 1 チャンク目の端点のずれ | 0.0737 → 0.0008 Å | 0.0220 → 0.0016 Å |
| 端点の持ち上がり（DFT SP） | +13.569 → +0.002 | +1.092 → +0.012 |
| v3 の予測 | 0.004 Å、+0.04 | 0.0014 Å、+0.009 |
| 初期経路の隣の画像の \|直接 RMS − 写像 RMSD\| の最大（剛体の跳び） | 0.269 → 3.4e-14 Å | 0.499 → 2.6e-16 Å |
| note の end_drift（チャンクごと） | 0.074 / 0.001 / 0.002 → 0.001 / 0.000 / 0.000 | 0.022 → 0.002 |
| チャンク数・形・HEI の位置 | 3・`single_max`・bead 4（同じ） | 1・`single_max`・bead 4（同じ） |
| bead の最大 | 1 チャンク目 69.74 → 71.75、最終チャンク 66.73 → 66.81 | 13.47 → 13.47 |

- 始点側（bead 1）のずれと持ち上がりは、v3・v4 とも 0.0000 Å・0.000 kcal/mol です。STO-3G の 2・3 チャンク目のずれも 0.0008 / 0.0024 → 0.0004 / 0.0004 Å に下がりました。
- NWChem が報告する凍結端点のエネルギーは nstep 0 の値のままなので、端点のゆがみは bead のプロファイルには現れません（v3 と v4 の bead 9 の値は 1.2e-9 Eh で一致）。持ち上がりは SP でだけ分かります。
- M9 の効果は、プロファイルの端点の近くに出ました。HONO の cis 側の隣の 2 本（1 始まりで bead 7・8）が、trans から 8.47 → 8.28、2.85 → 2.83 に下がりました。N−1 番目の bead が、ゆがんだ端点に向かって緩和しなくなったためです。
- ずれは原理的に 0 にはなりません。NWChem の zts_distance は原子ごとの距離の和（L1）を、座標ごとの線探索（刻み 0.001）で最小化します。Kabsch は二乗和を最小化するので基準が違い、しかも NWChem は毎反復で整列し直すので、bead N はわずかに回り続けます。その大きさは 0.001〜0.002 Å・+0.01 kcal/mol で、解像度（1 kcal/mol）の 1/100 です。design.md:130 の「0.004 Å 以下」は小分子 2 系の実測で、保証ではありません。note で監視しているので問題はありません。

**U5 区間の時間（HONO の string）**: 1 チャンクの時間は 162.9 → 175.4 s でした（反復数は同じ 20）。validation.md:170 は原因を「特定していない」としていますが、stdout から分かります（U5-N10）。

- NWChem の string（algorithm 3）は、QN fixed point のステップを棄却すると damped Verlet に切り替えます。
- v3 は 4 反復目に 1 回切り替えました（Verlet のステップは 16 回）。v4 は 5 反復目と 11 反復目の 2 回切り替え（10 回）、11 反復目で bead の評価が 1 回り（7 本）増えて 15.7 s（96.1 → 111.8 s。ふだんの反復は約 7 s）かかりました。評価の総数は 149 → 156 です。
- STO-3G でもチャンクごとの評価数は揺れています。初期経路が変わったことによる最適化の軌跡の違いで、M9 の系統的なコストではありません。
- 無駄は v3 から変わらず残っています。HONO の E max は 7 反復目以降 1.44e-5 Eh（0.01 kcal/mol）の幅に収まっているのに、20 反復目まで走りました（7 → 20 反復目で約 98 s。U5-N5）。

**U6 への波及**（FIND_PATH 由来の種）:

| 量 | HONO（1 チャンク） | STO-3G HCN（3 チャンク） |
|---|---|---|
| 種（HEI）の v3 と v4 の差 | 0.0036 Å | 0.0041 Å |
| 種の Hessian の負モード | xTB −744.4 → −736.3 cm⁻¹（採用） | DFT −1130.9 → −1131.6 cm⁻¹ |
| 方向との重なり | 0.5705 → 0.5702（宣言した二面角の勾配） | 0.6534 → 0.6510（最終経路の接線） |
| saddle | 3 → 3 ステップ | 3 → 3 ステップ |
| TS どうしの差 | 1.3e-4 Å、E_TS の差 7e-9 Eh | 6e-6 Å |
| 虚振動 | −681.0i → −681.5i | −1223.88i → −1223.89i |

- U6 区間の壁時間（種の Hessian + saddle + TS freq + QRC）は、HCN 35.7 → 35.1 s、HONO 91.2 → 88.0 s、NH3 41.2 → 41.3 s、STO-3G 13.4 → 13.2 s で、合計 181 → 178 s でした。差は揺らぎの範囲です。QRC は 3 反応の合計で 107 → 105 ステップです。
- saddle の JobStore キーは、HONO（5b21f557 → e423c384）と STO-3G（ae0c3a2a → ec7087ab）で変わりました。種が M9 の後の経路から出ていることの裏付けです。HCN のキーの変化（bfa49214 → 6ee0230e）は上流の xTB opt のキーの変化（run ごとの数値の揺らぎ）によるもので、M9 ではありません。NH3 は 99b68d44 のまま同じでした。

### 2.4 C38 の実測

- v4 の 5 run の resolved_config.yaml は、すべて `00680e0fa52527fa1e18ec72a33db492cdb14a0c`（-dirty なし）でした。v3 は 7 run すべて `dc82145f…-dirty` です。v4 の起動ログも 5 件とも `sha=00680e0 dirty=0` でした。
- Windows 側で git を使った直後の WSL で数え直すと、通常の `git status --porcelain --untracked-files=no` は 83 行、`-c core.autocrlf=input` 付きは 0 行でした（`GIT_OPTIONAL_LOCKS=0` で index を書き換えずに測定）。Windows 側は両方とも 0 行です。
- index の EOL は i/lf 233・i/-text 4・i/none 10 で、CRLF の blob はなく、.gitattributes もありません。input の変換で隠れるのは改行だけの差で、内容の変更は引き続き -dirty になります。
- プローブ（`/home/user/hfauto_v4_probe/u8_crlf`、WSL の git 2.43 と Git for Windows 2.53）でも確かめました。
  - index が LF で作業ツリーが CRLF（Windows の checkout と同じ状態）のファイルは、stat を変えても clean でした。
  - 同じサイズの本当の編集は ` M` と検出されました。
  - `core.safecrlf=true` にしても、status は落ちずに差分を報告しました。
- **レビュー入力の訂正（U8-N12）**: 「CRLF のままコミットしたファイルは、autocrlf=input で常に -dirty に見える」は誤りです。index の blob がすでに CR を含むパスでは、git は CRLF → LF の変換をしません（convert.c の has_crlf_in_index）。プローブでも、通常の status と input 付きのどちらでも clean でした。潜在的な誤判定はなく、C38 に追加の対処は要りません。
- code_version は resolved_config.yaml と doctor の見出しに出るだけで、config_sha（runner.py:110-124）にも JobStore のキーにも入りません。キャッシュと再開の挙動は変わりません。

### 2.5 S26 の実施状況（v3 §5.2 の内訳との対応）

| v3 の内訳 | v4 の記述 | 確認 |
|---|---|---|
| validation.md に小節「レビューのプローブ」 | :66-72。C2（u0_c2）、S23（hcn_tz・tma_tz・tma_d）、手法パネル（v2 のプローブ）。:68 に負荷時の時間は参考値との断り | v2 レビュー（:330・:1344-1345）と v3 レビュー（§2.1、§4.0.2〜4.0.3）の数値と一致。u0_c2 の job.nw と stdout に `maxiter 50` と `maxit = 50` |
| validation.md の「未確認」3 か所 | :64・:106・:149 を「v3 の run では未実行。プローブで確認（§3）」に | 一致 |
| S21 の説明 | :61（歩幅の制限 34 回・全体の縮小 7 回、line search の補い、XMAX / XRMS の歩み） | 実ジョブの記録と一致 |
| M6 の帰属 | :86・:88 を「S15（M6 の寄与は 0）」に | 一致 |
| ΔG_assoc の注記 | :20・:51 に「SVPD、BSSE 未補正（TZVPD//SVPD −9.04）」 | 一致 |
| design.md:146（thermo 行） | vib_scale 1.0 と文献値 0.985（差は ΔG‡ ≤ 0.06、ΔG_assoc ≤ 0.09）、qRRHO の幅約 1.2、単原子の単量体・反応のない錯体では ΔG_assoc が出ない、同じ (Hill, q, m) の異性体 | ΔG_assoc は reaction thermo の中で計算しており（thermochemistry.py:207）、記述と一致 |
| design.md:150 | 複製した run では thermo_tz が元の thermo を上書きする | artifact id に stage id が入らない（thermochemistry.py:242・:248）。Manifest.union は後の artifact を優先する（manifest.py:80-91） |
| conditions の文言 | README.md:37、design.md:147、method_panel.yaml:5-6 を「順位を付けられる反応は thermo_unavailable になる（same_basin などは outcome だけ）」に | gates.py:297-317 と summary.py:100-103 の挙動に一致 |
| 手法誤差の符号と順位の目安 | design.md:149、README.md:86 を「HCN −1.2、HONO +2.0（2 系での見積もり）。ΔG‡ の差が約 3 kcal/mol 未満なら順序は入れ替わり得る」に | −1.19 と +2.02 の差 3.21 と合う。旧版の「約 2」は残っていない。コードの閾値にはしていない |
| （M9 と同じコミットの U5 の文書） | validation.md:57（ずれの原因を中間画像の 1 枚ずつの整列に訂正）、:147（「M8 の設計どおり」を撤回）、design.md:130（剛体の跳びと 3 か所の逐次整列） | M9 と実測に一致。ただし :147 の書き方には U6-N12・U6-N13 が残る |

### 2.6 v3 ロードマップの実施状況

| ID | 内容 | 実装 | 実計算での確認 | 備考 |
|---|---|---|---|---|
| M9 | 初期経路の逐次整列 | 済 | find_path_sto3g、HONO | — |
| S26 | 検証記録と文書の訂正 | 済 | 文書の照合 | 新しく見つかった不正確さ（low）は C45 |
| C38 | code_version の autocrlf | 済 | v4 の 5 run | — |
| C8、C10、C12〜C16、C20、C25、C27′、C28〜C34、C36、C37、C39′、C40〜C44 | could | 未 | — | §5.3 に引き継ぎ。C39′ の U5 分は freezeN の部分が M9 で済み、impose と gmax の 2 文に縮めた（U5-R9′）。C16 は前提の M9 が満たされ、着手できる状態になった |

### 2.7 本レビューで確かめたこと（validation の外）

- `/home/user/hfauto_v4_probe/m9`: `m9_ends.py`（`--sp` なし、読み取りのみ）と `m9_profile.py` を再実行し、保存済みの SP の stdout を読み直した。
- `/home/user/hfauto_v4_probe/u8_crlf`（と Windows 側の一時ディレクトリ）: C38 の改行変換の挙動（§2.4）。
- M9 のテストを、M9 前の `_initial_path` にメモリ上で戻して実行した（リポジトリは変更していない）。
- HONO と STO-3G の string の stdout（damped Verlet への切り替え、評価の回数、E max の推移）。
- NWChem 7.2.3 のソース（v7.2.3-release）:
  - string_input.F:148-173 は xyz_path の全 bead（両端を含む）を rtdb に置き、`bead_list:new` を .false. にする。そのため string.F の newchain 分岐（endgeom の読み込みと impose）は通らず、凍結端点は初期経路ファイルの先頭と末尾そのものになる。
  - zts_min_motion は bead i を bead i−1 へ剛体として合わせ、i = 2..nbeads に初回と毎回の更新で適用される（string.F:569-571、1967、2007 付近）。ソースのコメントにも「Can rotate/translater frozen N bead」とある。
  - freezeN は勾配を 0 にするだけ（1894-1900）で、凍結 bead のエネルギーは nstep 0 のときしか計算しない（1830-1831）。
  - zts_distance（1204-1223）は原子ごとの距離の和（L1）。
- git の convert.c（has_crlf_in_index）。
- v3 と v4 の manifest・job_key・resolved_config の照合（非対象単位の継承根拠、§4.9）。

### 2.8 実計算でまだ通っていない経路

v3 §2.6 から変わっていません。今回、FIND_PATH の初期経路の分岐と、FIND_PATH 由来の大きな TS を明示的に加えました。

| 経路 | 状態 |
|---|---|
| M1（陰性結果ゲートの削除）、S3（行 14 の monotonic → BARRIERLESS）、行 13（中間体）、multi_max | 今回の系に経路に入るものがない |
| GS や SCREEN の IDPP を初期経路にする FIND_PATH | v3・v4 の実 run はどれも新しい IDPP。align_sequential は座標系に依存しないので、危険は小さい |
| FIND_PATH 由来で 10 原子を超える TS | 検証セットに該当する case がない（TMA·(HF)₂ は same_basin で終わる） |
| 負モード 2 本以上の saddle（S20 の Cartesian 分岐） | プローブだけ（C42 で回帰テスト） |
| higher_order、分割、QRC の再試行、`_register`、reassigned、periodic_nearest | 未到達 |
| DFT の mode-follow・soft 押し、M6 の副次的な負モード | 未発火 |
| 電荷のある系・開殻系（structures → NWChem odft、会合量） | 実 run 0 本 |
| 手法パネル（sp stage → panel_report） | v2 のプローブだけ（C43） |
| 複数の温度、1 bar / 1 M | 未実行 |

---

## 3. v3 の課題の解消状況（対象単位）

再評価した 4 単位で、延べ 64 行を判定しました（同じ課題が 2 つの単位に載るものと、複数の ID をまとめた行を含む）。内訳は resolved 30、partially 4、unresolved 25、not_in_scope 5 です。

- v3 で medium だった U5-N8 は resolved です。対象単位に残る medium は U5-N3（DLC の GS の脆さ）だけで、影響はコストに限られます。
- unresolved の 25 行は、ほとんどが v3 で could とした改良が未実施のものです。

第3ラウンドで状態が変わった課題:

| ID | v3 → v4 | 決め手 |
|---|---|---|
| U5-N8（medium）1 チャンク目の凍結端点のゆがみ | unresolved → resolved | M9 |
| U5-T1 1 チャンク目の初期経路の形を検査していない | unresolved → resolved | M9 のテスト（M9 前のコードで FAIL することを再現） |
| U5-D1 validation.md のずれの原因と意味の誤り | unresolved → resolved | validation.md:57・:147、design.md:130 |
| U5-N9 impose が効かないこと・FREEZEN の文書 | unresolved → partially | freezeN の実態は design.md:130 に入った。impose（design.md:144）が残る |
| U0-N7 = C38 WSL で code_version が常に -dirty | unresolved → resolved | C38 |
| U0-N8 順位差の目安「約 2」 | unresolved → resolved | S26（U0-R7） |
| U0-R3′ = U8-R7(1) プローブの記録 | 未実施 → resolved | validation.md:66-72 |
| U8-N6（= U0-N5 の残り）S22 の後の conditions の文言 | unresolved → resolved | S26 |
| U5-N8 の U6 への波及 | — → resolved | M9。種は約 0.004 Å 動いたが TS と判定は同じ |
| U6-N4 trust 0.1 の退行（validation の「144」を含む） | resolved → partially | 退行そのものは解消のまま。「144」の誤記を含めて判定し直した（C39′） |
| U8-10 method_panel.yaml が独自の conditions | unresolved → not_in_scope | 設計上受け入れたもの。注意書きの文言は S26 で正確になった |

表の U5 と U0・U8 の「unresolved」は、v3 で新しく見つかった課題（v3 の時点で未対応）を含みます。

### 3.1 U0 全体プロトコル

| ID | v3 | v4 | 根拠 |
|---|---|---|---|
| U0-N7（WSL で code_version が常に -dirty） | unresolved | resolved | C38（config.py:162-163）。v4 の 5 run は -dirty なし、v3 の 7 run はすべて -dirty。WSL 83 → 0 行 |
| U0-R6（C38 の実施） | 未実施 | resolved | v3 の指定どおり引数 2 つ。2 行への折り返しだけが違う |
| U0-R3 / U0-R3′（validation.md へのプローブの記録） | 未実施 | resolved | validation.md:66-72。:64・:106・:149 からも参照。プローブの置き場所 5 か所が WSL 上に実在 |
| U0-N8（「約 2 kcal/mol」が順位差の目安として楽観的） | 新規 | resolved | design.md:149、README.md:86 が「符号は反応で異なる、約 3 kcal/mol」。コードの閾値にはしていない |
| U0-R7（順位差の目安の書き直し） | 未実施 | resolved | 同上 |
| U0-N5 の残り（S22 の後の conditions の文言） | 文言が不正確 | resolved | README.md:37、design.md:147、method_panel.yaml:5-6。gates.py:297-317 と一致 |
| U0-N1（層が欠けた subject の ΔE が LOT の混成） | unresolved | unresolved | thermochemistry.py:190-193 の _kcal は v3 と同じ（C28） |
| U0-R4（C28） | 未実施 | unresolved | 第3ラウンドの範囲外 |
| v3 §4.0.3 の「◎ は条件付き」の条件 | — | partially | validation.md の古い記述、design.md:149 の出典、WSL の -dirty の 3 つは解消。method_panel を検証の run で実行していないことだけが残る |
| U0-I3・I7・I8・N2・N3・N4・N6（v3 で resolved） | resolved | resolved | 後退なし。input.py・gates.py・thermochemistry.py に変更なし |
| S24 のメッセージの「natural-abundance」（参考、U1-N8） | — | not_in_scope | v3 の割り当てどおり C40。00680e0 のコミットメッセージに明記 |

### 3.2 U5 paths（仮説・SCREEN・FIND_PATH）

| ID | v3 | v4 | 根拠 |
|---|---|---|---|
| U5-N8（medium）1 チャンク目の凍結端点のゆがみ | 新規 | resolved | M9。§2.3 の実測（ずれ 0.0008 / 0.0016 Å、持ち上がり +0.002 / +0.012） |
| U5-R8（M9） | 未実施 | resolved | 指定どおりの 2 行。チャンク数・形・HEI・鞍点・接続の判定は v3 と同じ |
| U5-T1（テストの穴） | 新規 | resolved | 新テストは現行コードで PASS、M9 前で FAIL |
| U5-D1（validation.md の誤り） | 新規 | resolved | validation.md:57・:147、design.md:130 |
| U5-N1（high、v2）凍結端点のずれを次のチャンクに渡す | resolved | resolved | M8 は有効。v4 の c2・c3 のずれは 0.0004 / 0.0004 Å |
| U5-N2（medium、v2）IDPP の打ち切り | resolved | resolved | S18 は有効。HONO はねじれ経路で、最大 13.47・bead 4 |
| U5-N3（medium）DLC の GS の脆さ | unresolved | unresolved | v4 の HONO も `screen:unavailable:screen_path:nonzero_exit`。影響はコストだけ（C35 は見送り） |
| U5-N4 / U5-R4（C36）S4 が GS 自体の失敗でも再実行 | unresolved | unresolved | worker.py は変更なし |
| U5-N5 判定がチャンクの終わりだけ | unresolved | unresolved | v4 の HONO で E max は 7 反復目以降 1.44e-5 Eh の幅なのに 20 反復目まで走った（約 98 s）。対策は v3 で見送り |
| U5-N6 design.md:130 の gmax の記述 | unresolved | unresolved | 「gmax と NWChem の収束判定は記録するだけ」が残る（C39′） |
| U5-N9 impose と FREEZEN の文書 | 新規 | partially | freezeN の実態は design.md:130 に入った。design.md:144 の impose が残る |
| U5-R9（C39′ の U5 分） | 未実施 | partially | freezeN の部分は済み。impose と gmax の 2 文が残る（U5-R9′） |
| U5-N7 試行を使い切った後の FIND_PATH | unresolved | unresolved | 見送り（v2）のまま |
| U5-I7 / U5-R7（C16）両端より低い DFT//xTB ノード | unresolved | unresolved | gates.py は変更なし。前提の M9 は満たされ、着手できる状態になった |
| U5-I8 / U5-R6（C15）tangent 用の DFT freq | unresolved | unresolved | 変更なし（NH3 の max_rel_dft 4.465 は v3 と同じ） |
| U5-R3（C35）pysis_gs の平面端点に off_plane | not_in_scope | not_in_scope | v3 で見送り。v4 でも退避経路で正しい鞍点に到達しており、判断を変える材料はない |

### 3.3 U6 paths（鞍点・検証・接続）

| ID | v3 | v4 | 根拠 |
|---|---|---|---|
| U6-I1 QRC に Hessian を渡さない・trust | resolved | resolved | v4 の QRC デッキ 8 件はすべて `maxiter 100 trust 0.3 inhess 2`。8 側すべてで max(traj) = traj[0] |
| U6-I2 mode_index と moddir の番号 | resolved | resolved | engine._moddir は変更なし。v4 の saddle 4 件は autoz・moddir 1 |
| U6-I3 periodic_nearest に許容幅がない | unresolved | unresolved | C18 は未実施。v4 でも未到達 |
| U6-I4 QRC の振幅のクリップ順 | resolved | resolved | 4 反応とも 1 回目の振幅で接続 |
| U6-I5 `_register` がゲートの前 | unresolved | unresolved | 設計どおり、発動なし |
| U6-I6 torsional の閾値 | not_in_scope | not_in_scope | hypotheses の担当 |
| U6-I7 試行を使い切った後の FIND_PATH | unresolved | unresolved | = U5-N7 |
| U6-N1（medium）Cartesian で moddir 0 の迷走 | resolved | resolved | v4 の saddle は 4/4 が moddir 1・autoz、AUTOZ の再試行なし |
| U6-N2 string 画像の未整列 | resolved | resolved | S19 と M9 で、初期経路を含むすべての経路が逐次整列 |
| U6-N3 凍結端点のずれ（= U5-N1） | resolved | resolved | M8 に加え、M9 で 1 チャンク目も直った |
| U5-N8 の U6 への波及 | — | resolved | M9。種の差 0.0036 / 0.0041 Å、TS の差 1.3e-4 / 6e-6 Å、Hessian の採否・重なり・判定は同じ |
| U6-N4 trust 0.1 の退行（validation の「144」を含む） | resolved | partially | S21 の効果は維持。validation.md:79 の「127 / 17 = 144」（正しくは 126 / 17 = 143）が残る（C39′） |
| U6-N5 continuation の docstring | resolved | resolved | engine.py:185-187 は v3 と同じ |
| U6-N6 qrc_bounds_A の上限と HESSIAN_NEAR_A | unresolved | unresolved | 文書に記述なし（C39′） |
| U6-N7 子反応の予算 | unresolved | unresolved | design.md:113 は「1 反応 6 h」のまま（C39′） |
| U6-N8 autoz / Cartesian の分岐の本数の数え方 | unresolved | unresolved | コードは変更なし。v4 の freq 17 件に −10 ≤ ν < 0 のモードはなく、実 run では起きていない |
| U6-N9 NWChem の smalleig | unresolved | unresolved | 文書に記述なし（C39′） |
| U6-N10 trust がジョブキーに入らない | unresolved | unresolved | 見送り。validation.md §10 の 5 に記載 |
| U6-N11 smoke の moddir テストが追ったモードを検査しない | unresolved | unresolved | tests/smoke は変更なし（C42） |
| U6-R5（C39′ の U6 分） | 未実施 | unresolved | 第3ラウンドの範囲外 |
| U6-R6（C42） | 未実施 | unresolved | 第3ラウンドの範囲外 |

U6-N4 は、v3 では S21 の退行そのものを resolved とし、「144」の誤記を C39′ に残していました。今回の判定は、誤記を含めて partially としたものです。

### 3.4 U8 sp・report・実行基盤

| ID | v3 | v4 | 根拠 |
|---|---|---|---|
| U0-N7 / C38（WSL で code_version が常に -dirty） | unresolved | resolved | §3.1 と同じ。test_config_load.py は 3 passed |
| U8-N6（S22 の後の「全反応が thermo_unavailable」） | 新規 | resolved | S26。summary.py:100-103（T・state が合わなければ t=None）と gates.py:307-310 に一致 |
| U8-R7(1) / U0-R3′（「未確認」の記述とプローブの記録） | 未実施 | resolved | validation.md:66-72、:64・:106・:149 |
| U8-R7(3) / U8-N7（CCSD の未収束の根拠の記述） | 新規 | resolved | v3 の割り当てどおり、レビューの注記だけで対応 |
| U8-3（CCSD の maxiter） | resolved | resolved | 実エンジンでの確認結果（maxit = 50、11/50・18/50 反復）が validation.md:70 に記録された |
| U8-5（composite と thermo の上書き） | resolved | resolved | design.md:150 に上書きの注意を追記 |
| U8-6 / U0-R7（同順位の幅と手法誤差の読み方） | resolved | resolved | design.md:149、README.md:86。ゲートにはしていない |
| S26 の U7・U3 由来の項目（ΔG_assoc の注記、M6 の帰属、thermo 行、S21 の説明） | 未実施 | resolved | §2.5 のとおり |
| U8-1、U8-2、U8-4、U8-8、U8-N1、U8-N3 | resolved | resolved | report.py・summary.py・gates.py・single_point.py は変更なし。v4 の ranking.csv は §2.2 のとおり |
| U8-7・U8-N5（ranking.csv に T・standard_state・dG_rxn がない） | unresolved | unresolved | v4 の見出しは `rank,reaction_id,outcome,tier,rankable,dG_act_kcal,band_low_kcal,band_high_kcal,blockers` のまま（C25） |
| U8-9（site.memory_mb が未使用） | unresolved | unresolved | config.py:36 と wsl_local.yaml:3-4 のコメントが残る（C27′） |
| U8-10（method_panel.yaml が独自の conditions） | unresolved | not_in_scope | 設計上受け入れたもの。注意書きの文言は S26 で正確になった |
| U8-N2（monitor 経路と STAGNATED が到達しない） | unresolved | unresolved | `def monitor` 7 か所、FailureKind.STAGNATED 4 か所が残る（C37） |
| U8-N4（手法パネルを正式な検証で実行していない） | unresolved | unresolved | 文書の空白は埋まった。正式な実行（C43）は未実施 |
| U8-N8（入力の不正で Rich のトレースバック） | unresolved | unresolved | cli は変更なし。見送りの判断は維持 |
| C40（S24 のメッセージの文言） | — | not_in_scope | U1 の担当。validation.md:61 の引用は観測の記録として正しい |

---

## 4. 単位ごとの再評価

各単位の「追加の改良案」には、§5 のロードマップ ID を【 】で添えました。

### 4.0 U0 全体プロトコルと理論レベルの選択

評価: v3 ○ / ◎ / ○ → **v4 ○ / ◎ / ○**（有用性 ◎ の条件がほぼ外れた）

#### 4.0.1 変更点

- **C38**: config.py:162-163 で、code_version() の git status 呼び出しを `_git(Path(tmp), git, "-c", "core.autocrlf=input", "status", "--porcelain", "--untracked-files=no")` にしました。rev-parse（:161）と戻り値の規則（:164-166）は変わっていません。専用テストは足さず、既存の test_config_load.py:27 の正規表現のままです（v3 の指定どおり）。
- **S26**: validation.md の小節「レビューのプローブ」（:66-72）と、:64・:106・:149 からの参照。design.md:149 と README.md:86 の手法誤差と順位の目安。README.md:37・design.md:147・method_panel.yaml:5-6 の conditions の文言。design.md:146・:150 の thermo と composite の追記。validation.md:20・:51・:86・:88 の注記と帰属の訂正。
- **変わっていないもの**: 3 層構成、method ファイル、input.py の CCSD ブロック（`freeze atomic`、`maxiter 50`）、ThermoSettings、パネルの既定の 2 手法。thermochemistry.py:190 の _kcal（U0-N1）も変わっていません。
- v4 の resolved_config.yaml は、code_version と一時 xyz のパスを除けば v3 と同じでした（5 run）。

#### 4.0.2 妥当性（○、据え置き）

C38 と S26 は正しく入っています。一方、○ に留めた v3 の理由のうち、U0-N1 と、電荷のある系・開殻系・method_panel を実 run で通していないことが残っています。U0-N8 は解消しました。

- **理論レベルは変わっていません。** hfauto/ の dcc624d..HEAD の差分は config.py（C38）と actions.py（M9）だけです。v3 と v4 の resolved_config.yaml を code_version を除いて比べると、HCN・HONO・NH3・水は差分 0 行で、STO-3G は一時 xyz のパスが 2 か所違うだけでした。
- **C38 は正しく、副作用もありません。**
  - `_git` の argv は `git -C <repo> -c core.autocrlf=input status --porcelain --untracked-files=no` です。`-c` はサブコマンドより前にあるので、git の大域オプションとして効きます。Windows 側のシステム設定（`C:/Program Files/Git/etc/gitconfig` の autocrlf=true）も、WSL 側の未設定も上書きします。
  - CRLF の blob はなく、.gitattributes もないので、input の変換で隠れるのは改行だけの差です（§2.4）。
- **U0-R7 の数値は整合します。** HCN は SVPD 46.62 − CCSD(T) 47.81 = −1.19、HONO は 13.67 − 11.65 = +2.02 で、両者の差 3.21 が「約 3 kcal/mol」と合います。design.md:149 と README.md:86 は同じ表現で、コードの閾値にはしていません。
- **conditions の新しい文言は実装と一致します。** gates.py:297-317 の rankable では、順位を付けない outcome は `outcome:<o>` と実際に立っている blocker だけを返し、rankable な outcome で dG_act が None なら thermo_unavailable を先頭に入れます。
- **design.md:150 の「複製では thermo_tz が元の thermo を上書きする」** は、thermochemistry.py:242・:248 の artifact id（`species_thermo_{subject}_{T}K`、`reaction_thermo_{id}_{T}K_{state}`）が stage id を含まないことと整合します。
- **結論と順位**: ΔG の v3 からの差は 0.0011 以下で、4 反応とも rank 1 です。同順位の幅（qs × cutoff）は 0.001 未満のままなので、順位の不確かさを読む手がかりが design.md:149 の目安だけという構図は変わりません。検証の run は 1 run 1 反応なので、実害はありません。

#### 4.0.3 有用性（◎、条件がほぼ外れた）

v3 で ◎ に付けた 4 つの条件のうち、3 つが解消しました。

1. **validation.md の古い記述**: :64・:106・:149 が「v3 の run では未実行。プローブで確認（§3）」に直りました。
2. **design.md:149 の数値の出典**: validation.md:66-72 から、C2 の +2.02 と、パネルの HCN −1.2 の出典（CCSD(T) 47.81）を追えるようになりました。
3. **WSL で -dirty になる問題（C38）**: 本番経路の code_version が、本当の未コミット変更だけを示すようになりました。v4 の 5 run は、コミット 00680e0 に一意に結び付きます。

残る条件は、method_panel を検証の run で実行していないことだけです。これは v3 §4.0.8 が「プローブの記録で足りる」とした範囲で、今回その記録が入りました。

利用者向けの読み方も良くなりました。

- 順位の読み方が、誤差の符号が反応で異なることと、約 3 kcal/mol という目安を含む表現になりました。
- 主対象の ΔG_assoc には「SVPD、BSSE 未補正、TZVPD//SVPD では −9.04」という LOT 依存の注記が付き、結論の確からしさを文書だけから判断できます。
- composite の手順に、thermo_tz が上書きするという注意が加わりました。

細かい点が 2 つあります（どちらも有用性を下げるほどではありません）。validation.md:72 は、CCSD(T) 47.81 がパネルに opt-in で足した計算だったことを書いていません（U0-N9）。README のベンチマークは v3 の値と表記のままです（U0-N10）。

#### 4.0.4 複雑さ（○、据え置き）

- **コード**: config.py に git の引数を 2 つ足しただけです（hfauto の正味 +1 行）。knob、設定項目、method ファイル、テストは増えていません。
- **文書**: S26 で増えたのは、既存の挙動の説明と数値の出典です。利用者が手で守る規則は増えていません。
- **◎ にしない理由**: v3 §4.0.4 と同じで、手で守る規則が多いことです（composite の 4 つの規則と、パネルの conditions の一致）。それをコードで止めると規則や knob が増えるので、○ が釣り合いのとれた点だという判断も変わりません。

#### 4.0.5 新たな課題

- **U0-N9（low）validation.md:72 が、HCN の CCSD(T) 47.81（design.md:149 の −1.2 の出典）をパネルへの opt-in の追加で得たことを書いていない**: validation.md:72 は「既定の method_panel を追記した…HCN は 6 SP・36 秒」に続けて、区切りなしに「HCN の CCSD(T)/def2-TZVPD の ΔE‡ は 47.81」と書いています。既定のパネルは 2 手法（6 SP = 3 構造 × 2 手法）で、CCSD(T) は v2 のプローブで methods に 1 語足した別の実行でした（v2 レビュー:343・:1182・:1195。C2 より前の maxiter 20 で、11〜13 反復で収束しているので値は有効）。このままだと既定のパネルに CCSD(T) が入ると読め、README.md:37 と method_panel.yaml:3-4 の「CCSD(T) は小さい閉殻分子のときだけ opt-in」と一見食い違います。値そのもの（47.81、−1.19）は正しいです。
- **U0-N10（low、対応不要）README のベンチマーク（:86-92）は v3 の表記と値のまま**: README.md:86 は「第2ラウンド改良後の再検証 v3」、:91 は HONO −681.0i・5 分 23 秒で、v4 は −681.5i・5 分 32 秒です。v3 の値に v3 と書いているので誤りではなく、差は ±9 s、0.5 cm⁻¹、ΔG 0.001 以下です。validation.md:3 が v4（§11）へ案内しているので、実害はありません（U8-N11 と同じ）。

#### 4.0.6 追加の改良案

- **U0-R8（could）validation.md:72 の CCSD(T) の値が opt-in の別の実行で得たものだと明記する** 【C45】: 「HCN の CCSD(T)/def2-TZVPD の ΔE‡ は 47.81」の前に、「続けて methods に ccsd-t_def2-tzvpd を opt-in で足した別の実行（C2 より前の maxiter 20、11〜13 反復で収束）で」と補う。v2 の run の記録そのものは書き換えない。design.md:149 の「HCN −1.2」の出典が正確に追え、既定のパネルに CCSD(T) が入るという誤読がなくなる。文書 1 か所、十数文字、中立。
- **C28（could、v3 から継続）layer_missing の subject を含む ΔE を None にする** 【C28】: thermochemistry.py:190-193 の _kcal で `sa.layer_missing or sb.layer_missing` なら None にし、test_thermo_sp_stages.py の S9 の節に assert を 1 つ足す。U0 の妥当性を ○ に留めている、コード側で唯一残る理由を解消する。composite で SP が 1 本だけ失敗したときに、2 つの LOT を混ぜた ΔE（HCN で 164 kcal/mol）が report.html に出なくなる。1 行、中立。

#### 4.0.7 見送った案

- **README のベンチマーク（:86-92）を毎ラウンド新しい値と表記に更新する**: 今の表記は v3 の値に v3 と書いていて正確で、validation.md:3 が v4 へ案内している。差は揺らぎの範囲で読み方は変わらず、ラウンドごとに README を書き換えるだけの変更になる。
- **code_version の git 呼び出しに `--no-optional-locks` を足し、WSL から Windows の checkout の index を書き換えないようにする**: v3 以前からの挙動で、index を書き換えても害はない（C38 の後はどちらの順で操作しても -dirty にならない）。lock の競合は観測されておらず、効果を測れない引数が増えるだけ。
- **C38 の argv を固定する専用の単体テストを足す**: subprocess をモックして実装の細部を固定するだけになる。実際の効果は v4 の 5 run と、WSL の 83 → 0 行で確かめられている。
- **U0 の目的で v4 の run に method_panel や composite を追記して実行する**: v3 §4.0.8 のとおり。パネルの DFT の経路は v2 のプローブ、C2 は v3 のプローブで同じエンジンのコードを通しており、その記録が validation.md §3 に入った。正式な実行は U8 の C43 の範囲。
- **手法誤差（約 3 kcal/mol）を同順位の幅に組み込む、パネルの conditions の一致をコードで検査する**: v3 §4.0.8 と同じ。組み込むと小分子ではほとんどの反応が同順位になり、ΔG‡ の序数という読み方が失われる。見積もりは 2 系だけで、閾値の根拠として薄い。conditions が食い違ったときは thermo_unavailable として目に見える形で失敗する。

### 4.5 U5 反応経路：仮説選択と障壁の事前判定（hypotheses.select / SCREEN / FIND_PATH）

評価: v3 ○（SCREEN ○、FIND_PATH ○ の下限）/ ○ / ○ → **v4 ○（SCREEN ○、FIND_PATH ○。下限ではなくなった）/ ○ / ○**

#### 4.5.1 変更点

- **M9（U5-R8）**: actions.py:455-456 の 2 行を、v3 の指定どおりに次のように置き換えました。旧コードは `inner = [align_mapped(ctx.ends[0], ...) for i in images[1:-1]]` で、中間画像を 1 枚ずつ ends[0] に整列していました。新コードは ends[0] を固定し、ends[1] を含む各画像を直前の画像へ剛体整列します（interpolation.py:32-39。Kabsch の固有回転・恒等写像、geometry.py:8-29）。

```python
inner = [np.asarray(i.coords, dtype=float) for i in images[1:-1]]
return ctx.path_file(name, align_sequential([ctx.ends[0], *inner, ctx.ends[1]]))
```

- **テスト（U5-T1）**: test_reaction_case_actions.py に `test_first_string_chunk_starts_from_an_idpp_without_rigid_jumps` を足しました。直線の HCN→HNC（LINEAR）と平面の cis→trans HONO（PLANAR）でパラメータ化し、判定条件は v3 の検証で修正した形です。1 チャンク目の initial_path の隣り合う画像すべてで「直接の RMS = mapped_rmsd」（abs 1e-6）、両端は frames[0] = a、mapped_rmsd(frames[-1], ends[1]) < 1e-6 です。付随して、DriftingPath が Failure をそのまま返すようにしました（2 行）。
- **文書（U5-D1）**: validation.md:57 でずれの原因を「中間画像を 1 枚ずつ整列していたこと（U5-N8）」に訂正し、:147 で「M8 の設計どおり」を撤回して string のキーが変わることを書きました。design.md:130 に、凍結端点が剛体の跳びで回されてゆがむこと、3 か所での逐次整列、「0.004 Å 以下」の目安を足しました。validation.md §11 に M9 の実測を記録しました。
- **変わっていないもの**: 仮説選択、SCREEN と barrier_verdict、決定表（state.py）、profile.py、find_path のチャンクのループ（M8・S19）、NWChem の string のデッキ、pysis の GS と S4、IDPP の定数。

#### 4.5.2 妥当性（○、FIND_PATH の「下限」を外す）

- **M9 は NWChem の実際の挙動に対して正しい修正です。** NWChem 7.2.3 のソースで v3 の診断を確かめ直しました（§2.7）。xyz_path を渡すと endgeom と impose は使われず、zts_min_motion は凍結 bead N も含めて毎反復で剛体として合わせ、freezeN は勾配を 0 にするだけです。したがって凍結端点をゆがめる原因は初期経路の剛体の跳びそのもので、M9 はそこを直しています。
- **ends[1] を回して渡しても問題ありません。** 初期経路の最後の画像（回った ends[1]）はデッキの endgeom と一致しなくなりますが、endgeom は読まれません。これは M8 のチャンク 2・3 と同じ扱いで、テストの mapped_rmsd(frames[-1], ends[1]) < 1e-6 でも、ends[1] が剛体として回るだけであることを確かめました。
- **実 run の値を読み直して再計算すると、報告値とすべて一致しました**（§2.3）。1 チャンク目の端点のずれは 0.0008 / 0.0016 Å、剛体の跳びは 1e-13 Å 未満、始点側は 0.0000 Å、NWChem が報告する bead 9 のエネルギーは v3 と 1.2e-9 Eh で一致しました。
- **行 14 の根拠が端点の近くでも成り立つようになりました。** HONO の cis 側の隣の 2 本が 0.19 / 0.02 kcal/mol 下がり、N−1 番目の bead がゆがんだ端点に向かって緩和しなくなりました。「連続な DFT 経路の最大は鞍点の上界」が、障壁の小さい領域でも成り立ちます。U5-N8 は解消しました。
- **テストは正しい形です。** 現行コードで 3 passed（M8 のテストを含む）でした。M9 前の `_initial_path` をメモリ上で戻すと、LINEAR・PLANAR ともに FAIL しました（隣の組 8 組中 7 組が不一致）。判定条件（直接 RMS = mapped_rmsd）は、直線の ends[0] で Kabsch が縮退しても成り立ちます。
- **ずれは原理的に 0 にはなりません**（L1 と二乗和の基準の違いと、毎反復の再整列。§2.3）。残りは 0.001〜0.002 Å・+0.01 kcal/mol で、design.md:130 の「0.004 Å 以下」は小分子 2 系の実測です。note で監視しているので問題はありません。
- **SCREEN と仮説選択は変わっていません。** HCN・NH3 は v3 と同じく proceed（max_rel_dft は HCN 33.984 → 34.003、NH3 4.465 → 4.465。HCN の 0.02 は実行ごとの揺らぎ）、HONO は unavailable → FIND_PATH です。
- **◎ にしない理由**は 3 つです。U5-I7（C16 は未実施）、行 14 と multi_max が実 run で一度も通っていないこと、GS や SCREEN の IDPP を初期経路にする FIND_PATH も実 run で通っていないこと（align_sequential は座標系に依存しないので、危険は小さい）。

#### 4.5.3 有用性（○、維持）

- 追加の計算なしで、1 チャンク目の凍結端点が真の極小になりました（持ち上がり +13.57 → +0.002、+1.09 → +0.012 kcal/mol）。1 チャンクで結論が出る string でも、障壁の小さい反応の BARRIERLESS / single_max の判定を信頼できます。
- 結論と ΔG は変わっていません（ΔG‡ 62.1678 → 62.1677、12.2256 → 12.2267）。
- HONO の string が 12.5 s 遅くなったのは、NWChem の damped Verlet への切り替えの回数の違いによる軌跡の揺れで、M9 の系統的なコストではありません（§2.3、U5-N10）。
- 無駄は v3 から変わらず残っています。HONO の E max は 7 反復目以降 0.01 kcal/mol の幅に収まっているのに 20 反復目まで走り（約 98 s、U5-N5）、HONO の GS は 8 s で失敗します（U5-N3・N4）。
- M9 で string の JobStore キー（初期経路の sha256、engine.py:388）が変わりました。新しい run dir には影響しません。古い run dir での挙動は U6-N13 のとおりです。

#### 4.5.4 複雑さ（○、維持）

- hfauto の U5 分の増減は 0 行です（actions.py は 500 行で上限ちょうど、import も不変）。knob、policy キー、決定表の行は増えていません。テストの追加は約 20 行です。
- 初期経路、チャンクの継ぎ目、最終経路の 3 か所が同じ align_sequential の規則になり、実装も design.md:130 の説明も単純になりました。
- 文書で残る食い違いは、design.md:144 の impose（U5-N9）と、design.md:130 の「gmax と NWChem の収束判定は記録するだけ」（U5-N6）の 2 か所です。actions.py に docstring を足す余地がないのは v3 と同じです。

#### 4.5.5 新たな課題

- **U5-N10（low）validation.md:170 が、HONO の string の遅れ（162.9 → 175.4 s）を「原因は特定していない」としているが、job の stdout で原因が分かる**: v3 の stdout（`jobs/fb/fb36e6c6…`）でも 4 反復目に「string: switching to damped Verlet」が 1 回出ています（Verlet のステップは 16 回）。v4 の stdout（`/home/user/hfauto_v4/hono_known_endpoints/jobs/4b/4bdb9131…/attempt_00/stdout.txt`）では、5 反復目と 11 反復目の 2 回切り替え（10 回）、11 反復目で bead の評価が 1 回り（7 本）増えて 15.7 s（96.1 → 111.8 s）かかりました。評価の総数は 149 → 156 です。STO-3G でもチャンクごとの評価数は揺れています。これは NWChem の QN fixed point と damped Verlet の軌跡の違いで、M9 の系統的なコストではありません。記録がないと、M9 の性能劣化を疑われます。
  - 評価段の当初の根拠（「v3 では切り替えが起きていない」）は、検証で誤りと分かったので訂正しました。結論は変わりません。

#### 4.5.6 追加の改良案

- **U5-R10（could、修正して採用）validation.md:170 に、HONO の string が遅くなった原因を記録する** 【C45】: 「原因は特定していない」を、たとえば次の 1 文に置き換える。「v3 も 4 反復目で damped Verlet に 1 回切り替えたが、v4 は 5 反復目と 11 反復目の 2 回切り替え、2 回目で bead の評価が 1 回り（7 本、11 反復目は 15.7 s）増えた（評価 149 → 156）。STO-3G でもチャンクごとの評価数は揺れており、M9 の系統的なコストではない」。「v3 には切り替えがない」とは書かない。文書 1 文、中立。
- **U5-R9′（could、v3 の U5-R9 から範囲を縮小）C39′ の U5 分の残りを impose と gmax の 2 か所に絞る** 【C39′】: design.md:144 の reaction-paths 行に「impose は xyz_path を渡すと NWChem が使わない（string.F の newchain 分岐の中だけ）」を 1 文足す。design.md:130 の「gmax と NWChem の収束判定は記録するだけ」を「NWChem の string は gmax を記録しない（pysis_gs だけ）。NWChem の converged は使わない」に直す。freezeN の説明は M9 で済んでいるので含めない。デッキは変えない。U5-N6・U5-N9 を閉じる。文書 2 文、中立。

#### 4.5.7 見送った案

- **align_sequential の整列基準を NWChem の zts_distance（L1）に合わせて、1 チャンク目のずれを 0 に近づける**: 残りは 0.0008〜0.0016 Å・+0.01 kcal/mol で解像度の 1/100。NWChem は毎反復で bead N を整列し直すので、初期経路をどう整列しても 0 にはならない。自前の L1 最適化を足すのに見合わない。
- **design.md:130 の目安を「0.002 Å 以下」に締める、end_drift をゲートにする**: 実測は小分子 2 系（4 チャンク）だけで、17 原子の錯体では bead N−1 の動きが大きく、ずれも大きくなりうる。note で監視すれば足りる（v3 §4.5.8 と同じ判断）。
- **`_initial_path` に、逐次整列の理由を docstring で書く**: actions.py は 500 行の上限ちょうどで、docstring を足すと上限を超える。理由は design.md:130 とテストのコメント（M9: NWChem rotates frozen bead N）にあり、誤って戻すとテストが落ちる。
- **GS 経路を初期経路にする分岐の単体テストを足す**: frames を選んだ後の処理は IDPP のときと同じで、align_sequential は座標系に依存しない。今のパラメータ化したテストで整列の規則は覆えている。
- **NWChem の string の algorithm を固定して、damped Verlet への切り替えによる時間の揺れをなくす**: 切り替えは QN のステップを棄却するための保護で、揺れは両方向に ±10% 程度。系統的な利得の証拠がなく、収束の頑健さを損なうおそれがある。
- **string デッキから impose を削る**: v3 と同じ判断。効いていないので削っても挙動は変わらず、文書に実態を書けば足りる（U5-R9′）。

### 4.6 U6 反応経路：鞍点精密化・TS 検証・接続確認（REFINE_SADDLE / VALIDATE_TS / CONNECT）

評価: v3 ◎ / ◎ / ○ → **v4 ◎ / ◎ / ○（据え置き）**

#### 4.6.1 変更点

- **U6 のコードは変わっていません。** refine_saddle、validate_ts、`_register`、connect、`_mode_index`、`_first_hessian`、Ctx.direction、engine._moddir（engine.py:87-94）、continuation の docstring（engine.py:185-187）は v3 と同じです。
- **M9 の間接的な波及**: 1 チャンク目の NWChem への入力が変わるので、FIND_PATH の種（`seed:path_hei`）と接線（最終経路上の profile.tangent）が変わります。
- **JobStore への影響**: string のキーが変わると HEI の構造が変わり、下流の種の Hessian、saddle（key_payload の molecule と hessian の sha）、TS freq、QRC のキーもすべて変わります（§2.3 で saddle のキーの変化を確認）。
- テスト: U6 の double-well のテストは変わらず、9 passed でした。

#### 4.6.2 妥当性（◎、維持）

- **M9 は U6 の前提を改善しました。** 前提とは、種を出す DFT string の 1 チャンク目が、ゆがんでいない凍結端点で進むことです。NWChem 7.2.3 のソースでは xyz_path の両端がそのまま凍結 bead になり（string_input.F:148-173）、M9 の後は zts_min_motion の回転がほぼ恒等になります（§2.7）。
- **実 run（v3 と v4 を読み取りだけで照合）**:
  - v4 の saddle 4 件は、どれも autoz で `maxiter 50 trust 0.1 sadstp 0.1 inhess 2 moddir 1`。stdout は `moddir = 1`、key_payload は cartesian False・moddir 1、試行は 1 回ずつ。SaddleClaim の notes は 4 件とも [] でした。
  - FIND_PATH 由来の種は約 0.004 Å 動きましたが、Hessian の採否、方向との重なり、saddle のステップ数、到達した TS（≤ 1.3e-4 Å）、接続の判定はすべて v3 と同じでした（§2.3 の表）。
  - FIND_PATH を通らない HCN と NH3 は、TS どうしの差が 7e-7 Å と 0 でした。
  - QRC の 8 側は、どれも `maxiter 100 trust 0.3 inhess 2` で、traj の最大は traj[0]、どの反応も 1 回目の振幅で接続しました。
  - v4 の freq 17 件に −10 ≤ ν < 0 のモードはなく、U6-N8 の条件は起きていません。
- **留保**は v3 と同じです。実 run は 3〜4 原子の共有結合系だけで、負モード 2 本以上の分岐はプローブだけ、柔らかい負モードは smalleig で追えないことがあります（U6-N9）。どれも効率の問題で、結果の誤りは VALIDATE_TS と connection のゲートが防ぎます。

#### 4.6.3 有用性（◎、維持）

- U6 区間の壁時間は合計 181 → 178 s で、差は揺らぎの範囲です。QRC は 3 反応の合計で 107 → 105 ステップ、ΔG‡ の差は 0.001 kcal/mol 以下で、outcome と順位は変わりません。
- M9 の U6 への効果は、数値の差ではなく前提を保証したことにあります。1 チャンクで結論が出る string でも、凍結端点のゆがみ（v3 で +13.6 / +1.09 kcal/mol）に引かれた経路から種を出すことがなくなりました。3〜4 原子の系では、種が約 0.004 Å 動いても同じ TS に同じステップ数で収束するので、差としては表れません。

#### 4.6.4 複雑さ（○、維持）

- U6 のコード、決定表、JobStore のキー構成は変わっていません。
- actions.py はちょうど 500 行で、上限に余裕がありません。U6 を今後変える場合（C18 の許容幅など）は、行数を増やさない変更にする必要があります。
- 負モードの本数の数え方が 2 通りある小さな不整合（U6-N8）と、smoke の moddir テストが弱い点（U6-N11）は v3 と同じです。

#### 4.6.5 新たな課題

- **U6-N12（low）validation.md:147 の「FIND_PATH の HEI は bead の最大」「HONO で HEI 13.5」が実装と違う**: profile.hei（profile.py:54-75）は、bead の最大とその両隣を通る放物線（offset は ±0.5 でクリップ）からエネルギーと分数位置を出し、座標は隣の bead との線形補間にします。v4 の HONO の string のエネルギーから再計算すると、k = 4.301、bead の最大 13.47、HEI 13.63、鞍点（ΔE‡）13.67 kcal/mol でした（13.5 は bead の最大の値）。STO-3G の c3 も k = 3.937 で、HEI は bead の最大より +0.04 kcal/mol 高くなります。結論（鞍点をわずかに下回ることがあり、種としては問題ない）は変わりません。この記述は d7ff974 で入り、v3 のレビューは見落としていました。
- **U6-N13（low）validation.md:147 の、M9 より前の run dir を再開したときの再計算の説明が不正確**（= U8-N10）: validation.md:147 は「M9 より前の run dir で再開すると string を 1 回計算し直す」と書いていますが、実際は次のとおりです。
  - runner._fresh（runner.py:204-215）は、status が done で input_sha と config_sha が同じ stage を飛ばします。config_sha（runner.py:110-124）は code_version を含みません。したがって、完了した run dir を普通に再開しても paths は飛ばされ、string を含めて何も計算し直しません。
  - paths を再実行するとき（`--from paths`、未完了の stage、設定の変更）は、string のキーが変わって HEI の構造が変わるので、FIND_PATH 由来の種の Hessian、saddle、TS freq、QRC の opt、その TS の thermo も 1 回ずつ計算し直します（HONO の saddle 5b21f557 → e423c384、STO-3G ae0c3a2a → ec7087ab。SCREEN 由来の NH3 は 99b68d44 のまま）。着く TS は同じなので害はありません。

#### 4.6.6 追加の改良案

- **U6-R7（could）validation.md:147 の HEI と再開時の説明を直す** 【C45】: コードは変えない。
  1. 「HEI は bead の最大」を「HEI は bead エネルギーの放物線補間の最大（offset は ±0.5 でクリップ、座標は隣の bead との線形補間）」に直し、数値を「HONO で bead の最大 13.47、HEI 13.63、鞍点 13.67 kcal/mol」にする。
  2. 再開の文を「完了した paths stage は、code_version を config_sha に含まないので再開では飛ばされる。`--from paths` などで paths を再実行すると、string と、FIND_PATH の種から始まる U6 のジョブ（種の Hessian、saddle、TS freq、QRC）と、その TS の thermo が 1 回ずつ計算し直される（同じ TS に着く）」にする。U8-R9(b) と同じ箇所なので 1 つにまとめる。
  - 検証記録が profile.hei、runner._fresh、キーの構成と一致し、利用者が再開の挙動を誤解せずに済む。文書 2 行、中立。
- **C39′ の U6 分（could、v3 から継続）** 【C39′】: v3 の U6-R5 のとおり。design.md:144 に、autoz / Cartesian の分岐は ν < 0 の本数（noise 以下を含む）で決まること、NWChem の smalleig で柔らかい負モードを追えないことがあること、子反応はそれぞれ 6 h の予算（最悪 7 case）で qrc_bounds_A の上限は HESSIAN_NEAR_A（0.5 Å）以下であることを書く。validation.md:79 の「127 / 17 = 144」を「126 / 17 = 143」に直す。U6-N6〜N9 と「144」の食い違いがなくなる。文書 4〜5 行、中立。U6-R7 と同じ回に行う。
- **C42（could、v3 から継続）smoke の moddir テストで、NWChem が追ったモードを P·H·P と照合する** 【C42】: v3 の U6-R6 のとおり。stdout の最初の 'The mode being followed to the saddle point' のブロックを読み、P·H·P の第 2 固有ベクトルとの |cos| > 0.99 を検査する。読み取りの helper はテストの中に置く。実 pipeline で通っていない S20 の Cartesian 分岐の前提（U6-N11）を回帰テストで守れる。smoke に約 8 行（real マーカー付きで WSL だけで実行）、中立。優先度は低い。

#### 4.6.7 見送った案

- **end_drift を種の採否のゲートや警告の閾値にする**: M9 の後の 1 チャンク目のずれは 0.0008〜0.0016 Å で、種と TS に差は出ていない。knob が増えるうえ、結果の誤りは VALIDATE_TS と connection がすでに防いでいる。
- **HEI の補間構造で DFT の SP を取り直してから種にする、または bead の最大そのものを種にする**: HEI からの saddle は 3 ステップで収束し、TS までの写像 RMSD は 0.008〜0.025 Å しかない。SP を 1 本足すコストと分岐に見合う効果がない。
- **code_version を config_sha に入れ、コードを変えたら stage を再実行させる**: コミットのたびにすべての再開が無効になる。M9 では U6 の結果は変わらない（同じ TS）。
- **string デッキから使われない endgeom と impose を消す**: xyz_path を渡すと NWChem はどちらも使わないので、消しても結果は変わらない。impose の実態は文書（U5-R9′）で扱う。
- **FIND_PATH の種の接線が、ゆがみのない経路から来ることを検査する U6 のユニットテストを足す**: M9 のテストと既存の double-well の U6 テストで十分。fake の上で接線の重なりを検査しても、実 NWChem の freezeN の挙動は再現できない。
- **M9 の U6 への効果を確かめるため、FIND_PATH 由来の TS を持つ大きな錯体の実 run を足す**: 今の検証セットに該当する case がなく（TMA·(HF)₂ は same_basin で終わる）、1 case でも数時間かかる。M9 の効果は、凍結端点のずれの測定ですでに直接確かめてある。

### 4.8 U8 一点計算・手法パネル・順位付けと報告・実行基盤（sp・method_panel・report・execution・configs）

評価: v3 ◎ / ○ / ○（有用性は ○ の上側）→ **v4 ◎ / ○ / ○（有用性は ○ の上側。◎ に上げる条件は C43 と C25）**

#### 4.8.1 変更点

- **C38**: §4.0.1 のとおり。v3 以降に変わった U8 のコードは config.py だけで、正味 +1 行です。
- **S26**: conditions の文言（README.md:37、design.md:147、method_panel.yaml:5-6）、手法誤差と順位の目安（design.md:149、README.md:86）、プローブの記録（validation.md:66-72）、U7 由来で報告の読み方に効くもの（validation.md:20・:51 の ΔG_assoc の注記、design.md:146・:150、validation.md:86・:88 の帰属、:61 の S21 の説明）。
- **de9daf8**: validation.md:1-3 の見出しと §11 への案内、§11（:152-196。C38 の節は :188-191）、:141 を「v2・v3・v4 で修正なし」に。README は更新していません。
- **変わっていないもの**: single_point.py、input.py、report.py、summary.py、html.py、gates.py、execution/*、preflight.py、layout.py、cli。設定項目、失敗の種類、出力列は増えていません。

#### 4.8.2 妥当性（◎、維持）

- **C38 は正しい形で入っています**（§2.4、§4.0.2）。v4 の 5 run すべてで -dirty が付かないことを WSL で直接確かめ、code_version() は WSL の prod venv でも Windows の .venv-win-qa でも `de9daf8…`（-dirty なし）でした。test_config_load は 3 passed です。
- **S26 の文書はコードと合っています**（§2.5）。conditions の文言は summary.py:100-103 と gates.py:307-310 の挙動に、手法誤差の数値は v2・v3 のレビューに、プローブの小節の数値はすべて出典に一致しました。
- **報告の出力に回帰はありません**（§2.2）。ranking.csv の ΔG‡ は、HCN 42.17825 → 42.17826、HONO 12.22557 → 12.22670、NH3 3.83946 → 3.83946、STO-3G 62.16778 → 62.16771 で、4 反応とも rank 1、blockers は空です。水は `outcome:same_basin` だけで、method_panel.csv の sign_disagreement は各 run とも False です。report stage はどの run でも done です。
- **留保（low）が 2 点あります。** code_version は実行のたびに上書きされる一方、done の stage は版によらず飛ばされるので、版をまたいで再開した run の来歴は正確ではありません（U8-N9）。validation.md:147 の再開の記述は、paths を再実行する場合にしか成り立ちません（U8-N10 = U6-N13）。

#### 4.8.3 有用性（○、上側、維持）

v3 で ○ に留めた理由 4 つのうち、2 つが解消しました。

- (c) validation.md の「C2・パネルは未確認」は、:64・:106・:149 とプローブの小節 :66-72 に置き換わりました。design.md:149 の数値の出典も文書から追えます。
- (d) conditions の文言（U8-N6）は正確になりました。

C38 で、WSL の本番経路でも code_version が本当の未コミット変更だけを示すようになりました（v2 の 9 run と v3 の 7 run はすべて -dirty でした）。来歴の記録として初めて使える値になり、doctor の見出しにも同じ値が出ます。

残る理由は 2 つです。

- (a) ranking.csv に T・standard_state・dG_rxn がなく、CSV だけでは読めません（C25）。
- (b) 手法パネルは、正式な検証の run ではまだ一度も実行されていません。sp stage → panel_report の経路の記録は、v2 のプローブだけです（C43）。

この 2 つが片づけば ◎ にできます。

#### 4.8.4 複雑さ（○、維持）

- U8 のコード差分は正味 +1 行で、設定項目、失敗の種類、出力列、状態は増えていません。文書は増えましたが、コードの分岐には影響しません。
- 削減の余地は v3 と同じく残っています。到達しない monitor 経路（`def monitor` 7 か所）と FailureKind.STAGNATED（crest.py:134、nwchem/output.py:230、pysis/engine.py:148、xtb.py:57）、使われない site.memory_mb と wsl_local.yaml:3-4 の効かない値を勧めるコメントです（C37、C27′）。この 2 つの削減が ◎ に上げる条件です。

#### 4.8.5 新たな課題

- **U8-N9（low）code_version は pipeline_id ごとに毎回上書きされ、done で飛ばした stage を計算した版が残らない**: runner.py:266 は plan の前に毎回 write_resolved_config を呼び、layout.py:195-201 は records[pipeline_id] を丸ごと置き換えます。何もしない再開でも code_version は新しい版になります。一方、config_sha（runner.py:110-124）は code_version を含まないので、done の stage は版が変わっても飛ばされます（_fresh、:204-215）。run_state.json の stage 項目にも版の記録はありません。以前からの性質で回帰ではありませんが、C38 で値が信頼できるようになったので、「全 stage がこの版で計算された」と誤読する余地が意味を持つようになりました。design.md:158 は resolved_config.yaml の名前を挙げるだけで、この性質を書いていません。
- **U8-N10（low）validation.md:147 の再開の記述は、paths を再実行する場合にしか成り立たない**（= U6-N13）: 完了済みの M9 より前の run を再開しても paths は飛ばされ、M9 より前の string の結果が残ります。計算し直すのは、paths が未完了のときか `--from paths` のときだけです。U8-N9 と重なると、「code_version は新しい版なのに string は古い版」という記録になります。
- **U8-N11（low、対応不要）README のベンチマーク表は v3 の値のまま**（= U0-N10）。
- **U8-N12（low、レビュー入力の訂正）「CRLF のままコミットしたファイルは autocrlf=input で常に -dirty」は誤りで、潜在的な誤判定はない**: §2.4 のとおり（convert.c の has_crlf_in_index）。現状の index に CRLF のファイルは 0 件なので、いずれにせよ C38 に追加の対処は要りません。

#### 4.8.6 追加の改良案

- **U8-R9（could）code_version の意味と再開時の上書き、M9 の取り込みに要る操作を文書化する** 【C45】: コードとスキーマは変えない。
  - (a) design.md:158 の run ディレクトリの行に 1 文を足す。「resolved_config.yaml の code_version は git の sha（追跡中のファイルに差分があれば -dirty、改行コードだけの差は数えない、git がなければ unknown）で、pipeline_id ごとに最後に実行したときの値に上書きされる。done で飛ばした stage を計算した版は残らないので、版をまたいだ結果が要るときは新しい run dir か `--from` で取り直す」。
  - (b) validation.md:147 の再開の文を、paths を再実行する場合（未完了か `--from paths`）に限定する。U6-R7 の (2) と同じ箇所なので 1 つにまとめる。
  - U8-N9・U8-N10 を塞ぎ、code_version の仕様の説明がコードの外にも残る。文書 2 文、中立。
- **C43（could、継続）method_panel を v4 の run の複製で正式に 1 回実行し、validation.md に記録する** 【C43】: `cp -a /home/user/hfauto_v4/hono_known_endpoints` と hcn の複製に、既定の method_panel を `--run-dir` で追記する（HCN では methods に ccsd-t_def2-tzvpd も足す）。panel_report/ranking.csv と report/ranking.csv の一致、method_panel.csv の行、SP の数と時間を validation.md に記す。元の v4 の run は変えない。コードは変えない。計算数分、文書数行、中立。有用性を ◎ に上げる条件の 1 つ。
- **C25（could、継続）ranking.csv に T_K・standard_state・dG_rxn_kcal を足す** 【C25】: RankRow に 3 フィールドを足し、`_RANK_HEADER`（summary.py:253）と `_rank_cells` に列を加える。テスト 1 件。数行、中立。有用性を ◎ に上げるもう 1 つの条件。
- **C37（could、継続）使われなくなった monitor 経路と FailureKind.STAGNATED を削除する** 【C37】: v3 §4.8.7 のとおり。6 アダプタの monitor() と STAGNATED の 4 か所の分岐を削除し、`_supervise` は timeout の wait とプロセスグループの kill だけを残す。約 30 行の削除とテスト 3 か所、【複雑さ減】。
- **C27′（could、継続）使われない site.memory_mb と wsl_local.yaml:3-4 の誤解を招くコメントを削除する** 【C27′】: NWChem が実際に使うのは memory_mb_per_rank（method.py:108）だけで、SiteConfig.memory_mb はどこからも読まれていない。検査を残すなら ranks × memory_mb_per_rank で比べる。数行、【複雑さ減】。

#### 4.8.7 見送った案

- **stage ごとに code_version を run_state.json に記録する（StageState にフィールドを足す）**: スキーマを広げる変更で、v3 §5.5 の方針に反する。検証は毎回新しい run dir で行っており、版をまたいだ再開の誤読は U8-R9 の 1 文で防げる。
- **code_version を config_sha に入れ、コードが変わったら done の stage を計算し直す**: コミットのたびに全 stage が stale になり、再開と JobStore の価値がなくなる。化学に効く変更は、すでにジョブのキーか設定の sha に表れる。
- **git status に `--no-optional-locks` を付ける**: 観測された害がない（§4.0.7 と同じ）。
- **C38 の回帰テスト（一時 repo に CRLF の作業ツリーを作り、_REPO_ROOT を差し替えて -dirty にならないことを確かめる）**: テストに git を必須にし、約 10 行と monkeypatch が増える。対象は定数の引数で、v4 の実 run（5/5 が clean）と WSL の 83 → 0 行で確認済み。
- **.gitattributes（`* text=auto eol=lf`）を置き、`-c` の引数の代わりにする**: リポジトリの正規化と全利用者の checkout に影響する広い変更。C38 は引数 2 つで同じ目的を果たしている。
- **report.html と ranking.csv に code_version を表示する**: 表示だけのための配管（v3 §5.5）。resolved_config.yaml は同じ run dir にあり、追記した pipeline ごとに版が違い得るので、1 つの値では表せない。
- **README のベンチマーク表を v4 の値に更新する**: §4.0.7 と同じ。
- **S24 の「natural-abundance masses」の文言を U8 の範囲で直す**: v3 の割り当てどおり、U1 の C40 で扱う。validation.md:61 の引用は観測の記録で、書き換えない。

### 4.9 非対象単位（U1〜U4、U7）の継承根拠

**共通の根拠**

- `git diff dcc624d HEAD --stat` で変わったのは 7 ファイル（+104 / −24）です。コードの変更は `hfauto/drivers/reaction_case/actions.py` の `_initial_path`（M9、U5 の FIND_PATH）と、`hfauto/pipeline/config.py` の `code_version`（C38、U0 の記録）の 2 か所だけです。ほかは README・design.md・validation.md、method_panel.yaml のコメント、reaction_case のテストです。
- `hfauto/stages/`、`hfauto/backends/`、`hfauto/chemistry/`、`hfauto/drivers/minimum.py`、`hfauto/execution/`、`hfauto/core/`、`configs/methods`・`configs/systems`・`configs/sites`、`discover.yaml`、`known_endpoints.yaml` は差分 0 件です。
- `_initial_path` を呼ぶのは `actions.py:464` の `find_path` だけで、reaction_paths stage からしか通りません。
- `code_version` は `resolved_config.yaml` に記録されるだけで、stage 再利用の config_sha（runner.py:110-124）にも JobStore のキーにも入りません。キャッシュと再開の挙動は変わりません。
- M9 の再検証で再実行した 5 run（HCN・HONO・NH3・水・STO-3G）について、v3 と v4 の manifest を照合しました（照合スクリプトはリポジトリの外）。5 単位のどれにも、挙動の変化は見つかりませんでした。

| 単位 | v4 の評価 | 根拠 |
|---|---|---|
| U1 入力構造 | ◎ / ◎ / ◎（継承） | 構造の生成と検証に関わるコード・設定に差分なし。5 run とも structures stage の payload が v3 と完全一致。S24 のメッセージの文言は C40 のまま |
| U2 配座探索 | ○ / ○ / ◎（継承） | `stages/conformer_search.py`、`backends/crest.py`、discover.yaml に差分なし。配座探索が通る run（#5 TMA·(HF)₂、#9 (HF)₃）は再実行していないので、根拠は差分がないことだけ |
| U3 極小の確定 | ◎ / ◎ / ○（継承） | `stages/minima.py`、`drivers/minimum.py`、NWChem backend に差分なし。5 run の dft stage は artifact id と job_key が v3 とすべて同じ。極小の ΔE の差 ≤ 1.0e-6 kcal/mol、構造の RMSD ≤ 1.5e-7 Å、振動数の差 ≤ 0.006 cm⁻¹（4 rank MPI の数値ノイズの大きさ）。validation.md:61 の S21 は説明の補足だけ |
| U4 反応探索 | ○ / ○ / ○（有用性は下限、継承） | `stages/explore.py`、`backends/readuct`、`chemistry/trials.py`・`hypotheses.py`、discover.yaml に差分なし。探索が通る run（#5、#6 amine パイロット）は再実行していないので、根拠は差分がないことだけ |
| U7 熱化学 | ◎ / ○ / ○（継承） | `stages/thermochemistry.py`、`chemistry/thermo.py`、`backends/goodvibes` に差分なし。5 run の極小の species_thermo の ΔG の差 ≤ 5.1e-6 kcal/mol で、差のあったキーは G・H・ZPE・population の数値だけ。HONO の ΔG‡ の +0.001 は TS の G の変化（M9 で鞍点の出発構造が変わった U5/U6 側の差）で、熱化学のコードの差ではない |

補足:

- v3 で S26 に割り当てた U7 と U3 の文書の課題（U7-N5 vib_scale の根拠、U7-N7 の帰属、U7-N8 反応のない錯体、U7-N9 単量体アンサンブルの異性体、U7-N11 ΔG_assoc の LOT の注記、U7-I5 の感度幅の文書化、U3-N7 の S21 の記録の説明）は、記述が入り、コードと合うことを U0・U8 の検証で確かめました（§2.5）。「反応を持たない錯体では ΔG_assoc は出ない」は、ΔG_assoc を reaction thermo の中で計算している既存の挙動（thermochemistry.py:207）と合っています。これらは文書上の対応で、U7 と U3 の評価は動かしません。
- U2 と U4 は、実 run の照合がなく、根拠はコードと設定に差分がないことだけです。挙動を変える経路がないので、継承で足ります。

---

## 5. 残る課題と今後の候補

### 5.1 残る課題

**medium**

| ID | 内容 | 影響 | 対応 |
|---|---|---|---|
| U5-N3 | DLC の GS の脆さ（HONO で決定的に失敗） | コストだけ（GS 8 s。退避先の FIND_PATH が正しい鞍点に到達する） | C35 は見送りのまま |

**low（第3ラウンドで見つかったもの）**

| ID | 内容 | 対応 |
|---|---|---|
| U6-N12 | validation.md:147 の「HEI は bead の最大」「HEI 13.5」が実装（放物線補間、HEI 13.63）と違う | C45（U6-R7） |
| U6-N13 = U8-N10 | validation.md:147 の再開の記述が、条件（paths の再実行時だけ）と範囲（U6 のジョブも）の両方で不正確 | C45（U6-R7、U8-R9） |
| U5-N10 | validation.md:170 の「原因は特定していない」（damped Verlet への切り替えの回数の違いで説明できる） | C45（U5-R10） |
| U0-N9 | validation.md:72 が CCSD(T) 47.81 を opt-in の別の実行で得たことを書いていない | C45（U0-R8） |
| U8-N9 | code_version は再開のたびに上書きされ、版をまたいだ再開の来歴が不正確（以前からの性質） | C45（U8-R9） |
| U0-N10 = U8-N11 | README のベンチマーク表は v3 の値のまま（「v3」と明記しており誤りではない） | 対応不要 |
| U8-N12 | レビュー入力の訂正（CRLF のままコミットしたファイルでも -dirty にはならない） | 対応不要 |

**low（v3 から残るもの、単位別）**

- U0: U0-N1（C28）。
- U1〜U4、U7: v3 §4.1〜§4.4、§4.7 の残る課題のまま（本ラウンドでは再評価していない）。U7 の文書の課題は S26 で対応済み（§4.9）。
- U5: U5-N4（C36）、U5-N5（見送り）、U5-N6・N9（C39′）、U5-N7（見送り）、U5-I7（C16）、U5-I8（C15）。
- U6: U6-I3（C18、据え置き）、U6-I5（設計どおり）、U6-N4 の「144」・N6・N7・N8・N9（C39′）、U6-N10（見送り）、U6-N11（C42）。
- U8: U8-7・N5（C25）、U8-9（C27′）、U8-N2（C37）、U8-N4（C43）、U8-N8（見送り）。

### 5.2 今後の候補

- **must**: ありません。結論を誤らせる課題も、主要機能が働かない課題も残っていません。
- **should**: ありません。
- **could**: 下表のとおりです。どれも新しい knob を増やしません。「v3 のまま」は、本ラウンドで再評価していない単位の改良で、v3 §5.3 の判定（採用・修正して採用）をそのまま引き継いだものです。

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| C45 | 検証記録の訂正（第3ラウンド分）。内訳は下表 | U0-R8、U5-R10（修正）、U6-R7、U8-R9、U0-N9・U5-N10・U6-N12・U6-N13・U8-N9・U8-N10 | 検証記録と design.md が実装と一致する。再開の挙動と code_version の意味を誤読しない | 文書 5〜6 か所 | 中立 |
| C39′ | 文書の訂正の残り（下表。U5 分は impose と gmax に縮小） | v3 の C39′、U5-R9′ | 実装・実測・NWChem のソースと文書の食い違いをなくす | 文書 10〜15 行と docstring 1 行 | 中立 |
| C28 | layer_missing の subject を含む ΔE を None にする | U0-R4 = U7-R1、U0-N1 | 混成 LOT の ΔE が report に出ない。U0 の妥当性を ○ に留めるコード側の理由を消す | 1 行、assert 1 つ | 中立 |
| C25 | ranking.csv に T_K・standard_state・dG_rxn_kcal | U8-R3、U8-7、U8-N5 | CSV が単独で読める。U8 の有用性 ◎ の条件 | 数行、テスト 1 件 | 中立 |
| C43 | method_panel を v4 の run の複製（HONO、HCN + CCSD(T)）で正式に実行し、validation.md に記録 | U8-R8、U8-N4 | sp stage → panel_report の経路が正式な検証の記録に入る。U8 の有用性 ◎ の条件 | 計算数分、文書数行 | 中立 |
| C37 | monitor 経路と FailureKind.STAGNATED を削除（timeout の kill は残す） | U8-R5、U8-N2 | 到達しないコード約 30 行と失敗の種類 1 つが減る | 約 30 行の削除、テスト 3 か所 | 【減】 |
| C27′ | 未使用の site.memory_mb と誤解を招くコメントを削除 | U8-R4、U8-9 | 効かない knob がなくなる | 数行 | 【減】 |
| C42 | smoke の moddir テストで追ったモードを P·H·P と照合 | U6-R6、U6-N11 | S20 の Cartesian 分岐の前提を回帰テストで守る | smoke に約 8 行 | 中立 |
| C16 | 両端より低い DFT//xTB ノードがあれば unavailable（fake の E2E テストを同時に入れることが条件） | U5-R7、U5-I7 | SCREEN が中間体のある経路を BARRIERLESS で閉じることを防ぐ。前提の M9 は満たされた | 3 行、テスト 1 件 | 【増】 |
| C15 | tangent_mode_cm1 を barrierless のときだけ計算 | U5-R6、U5-I8 | 17 原子の縮退反応で 1 本約 900 s を省く | 約 5 行 | 中立 |
| C36 | S4 の再実行を GS が完了した場合に限る | U5-R4、U5-N4 | GS 自体が落ちたとき同じ GS を繰り返さない | 1〜2 行 | 中立 |
| C12 | drive の重複を WL 原子クラスで除き、列挙数・固有数・実行数を diagnostics.json に（v3 のまま） | U4-R1、U4-I4 | 31 attempt 中 21 件が同じ移動の状態が解消する | 15〜25 行、テスト 2 件 | 中立 |
| C14′ | TS がある negative にも ΔE‡（v3 のまま） | U4-R3、U4-I6 | 陰性の結論を障壁の値で裏付けられる | 4〜8 行 | 中立 |
| C41 | `_Scine.task` で捨てている例外の文字列を stderr に 1 行（C14′ と同じ変更で。v3 のまま） | U4-R4、U4-N4 | TSOpt の失敗原因を生ファイルから確かめられる | 1〜2 行 | 中立 |
| C13 | 無バイアスの opt に `convergence_max_iterations=500`（v3 のまま） | U4-R2、U4-I5 | amine 系の数値的な陰性がなくなる | 1 行 | 中立 |
| C34 | 単量体の停止後の再実行を「停止構造から既定の設定で 1 回」にし、`--noreftopo` を廃止（v3 のまま） | U2-R1、U2-N1 | 実 CREST で 2/2 失敗していた経路が 2/2 成功。フィールドと CLI フラグが 1 つずつ減る | 約 10 行、fake テスト 1 件 | 【減】 |
| C8 | 組成の多重度が一意でなければ INPUT_INVALID（v3 のまま） | U2-R4 | ラジカル対の組成で高スピンの PES を黙って選ばない | 2〜3 行 | 中立 |
| C10 | window で screen 極小がなければ拒否、include: all の文言（v3 のまま） | U3-R2、U3-ISS-5、U3-N4 | 空の dft stage を「反応なし」と読み違えない | 数行、テスト 1 件 | 中立 |
| C32 | rerank_sp の SP を per_state を超える組だけに出す（v3 のまま） | U3-R1、U3-N1 | 単量体の最低構造が SP の失敗で消えない | 4〜5 行、テスト 1 件 | 中立 |
| C33 | reaction-paths の `_register` に pipeline の gates を渡す（v3 のまま） | U3-R3、U3-N2 | 1 つの PES の tier の基準が 1 つにそろう | 1 行 | 中立 |
| C20 | 点群を ThermoResult.notes で記録し、S_rot を protocol から外す（v3 のまま） | U7-R4、U7-I7、U7-N6 | σ を監査でき、使われない経路が 2 つ片づく | 約 4 行 | 中立 |
| C29 | ΔG_assoc=None の理由を notes に（v3 のまま） | U7-R2、U7-N2・N3 | 無言の None と原因を区別できる | 約 10 行、テスト 1 件 | 【増】 |
| C44 | composite のとき report の見出しに energy 層の LOT を出す（v3 のまま） | U7-R7、U7-N10 | composite の report を SVPD の report と取り違えない | 約 5 行 | 中立 |
| C30 | 設定値の範囲検査（v3 のまま） | U2-R3 + U7-R3 | 設定の誤りを読み込み時に止める | 約 8 行、テスト 2 件 | 中立 |
| C31 | 読み込み時に xyz の存在と両端の q/m の一致を検査（v3 のまま） | U1-R3、U1-N2・N4 | `running` の残留と無駄な DFT を防ぐ | 約 5 行、テスト 2 件 | 中立 |
| C40 | S24 のメッセージを `hfauto uses most-abundant-isotope masses` に（v3 のまま） | U1-R5、U1-N8 | 実際の質量と文言がそろう | 1 行 | 中立 |

C45 の内訳（すべて文書。コードは変えない）:

| 場所 | 内容 | 由来 |
|---|---|---|
| validation.md:72 | 「HCN の CCSD(T)/def2-TZVPD の ΔE‡ は 47.81」の前に「続けて methods に ccsd-t_def2-tzvpd を opt-in で足した別の実行（C2 より前の maxiter 20、11〜13 反復で収束）で」を補う | U0-R8 |
| validation.md:147（HEI） | 「HEI は bead の最大」→「HEI は bead エネルギーの放物線補間の最大（offset は ±0.5 でクリップ、座標は隣の bead との線形補間）」。数値を「HONO で bead の最大 13.47、HEI 13.63、鞍点 13.67 kcal/mol」に | U6-R7(1) |
| validation.md:147（再開） | 「完了した paths stage は code_version を config_sha に含まないので再開では飛ばされる。`--from paths` などで paths を再実行すると、string と、FIND_PATH の種から始まる U6 のジョブ（種の Hessian、saddle、TS freq、QRC）と、その TS の thermo が 1 回ずつ計算し直される（同じ TS に着く）」 | U6-R7(2)、U8-R9(b) |
| validation.md:170 | 「原因は特定していない」→「v3 も 4 反復目で damped Verlet に 1 回切り替えたが、v4 は 5 反復目と 11 反復目の 2 回切り替え、2 回目で bead の評価が 1 回り（7 本、11 反復目は 15.7 s）増えた（評価 149 → 156）。STO-3G でもチャンクごとの評価数は揺れており、M9 の系統的なコストではない」 | U5-R10（修正） |
| design.md:158（run ディレクトリの行） | 「resolved_config.yaml の code_version は git の sha（追跡中のファイルに差分があれば -dirty、改行コードだけの差は数えない、git がなければ unknown）で、pipeline_id ごとに最後に実行したときの値に上書きされる。done で飛ばした stage を計算した版は残らないので、版をまたいだ結果が要るときは新しい run dir か `--from` で取り直す」 | U8-R9(a) |

C39′ の内訳（v3 からの更新）:

| 場所 | 内容 | 由来 | v4 での変化 |
|---|---|---|---|
| design.md:144（reaction-paths 行） | 「impose は xyz_path を渡すと NWChem が使わない」を 1 文 | U5-R9′、U5-N9 | freezeN の部分は M9 で済んだので外した |
| design.md:130 | 「gmax と NWChem の収束判定は記録するだけ」→「NWChem の string は gmax を記録しない（pysis_gs だけ）。NWChem の converged は使わない」 | U5-R9′、U5-N6 | 変わらず |
| design.md:144、validation.md:79 | autoz / Cartesian の分岐は ν < 0 の本数（noise 以下を含む）で決める。NWChem の smalleig で柔らかい負モードを追えないことがある。子反応はそれぞれ 6 h の予算（最悪 7 case）、qrc_bounds_A の上限は 0.5 Å 以下。「127 / 17 = 144」→「126 / 17 = 143」 | U6-R5、U6-N4・N6〜N9 | validation.md の行番号が :69 → :79 に |
| design.md:141（explore 行）、validation.md:91 | 状態ラベルの限界（配座・立体・E/Z だけ違う IRC の端は same_as_source）、−289i の判定の変化は M3 の効果 | U1-R4、U4 の C39 分 | 変わらず（再評価していない） |
| design.md:139、environment.md:64、conformer_search.py:9-10 | 開殻の組成と `--noopt`、「約 25 s」→ 約 1.4 s、crest の threads、docstring を 4 knob に | U2-R2 | 変わらず（再評価していない） |
| design.md §8（ladder の行） | 初期 Hessian ありの opt の継続は drv.hess を INHESS 0 の restart で引き継ぎ、trust 0.1 で続ける | U3-R5、U3-N8 | 変わらず（再評価していない） |

v3 で見送り・据え置きとした C1、C4 のラジカルの部分、C18、C21（文書化で代替済み）、C22、C35 は、判断を変える材料がないので再提案しません。

### 5.3 評価を動かしうる改良

評価の記号を上げる条件になっているものを、単位ごとにまとめました。

| 単位 | 今の評価 | 上げる条件 | 性格 |
|---|---|---|---|
| U0 | ○ / ◎ / ○ | 妥当性: C28 に加えて、電荷のある系・開殻系の実 run。複雑さ: 手で守る規則をコードで止めると規則や knob が増えるので、○ が釣り合いの点 | 妥当性の残りは主に検証範囲 |
| U2 | ○ / ○ / ◎ | 複雑さの但し書き: C34、C8、C39′ の U2 分 | 小さなコード・文書 |
| U4 | ○ / ○（下限）/ ○ | 有用性の下限: C12、C14′ | 小さなコード |
| U5 | ○ / ○ / ○ | 妥当性: C16（【複雑さ増】）と、行 14・multi_max を実 run で通すこと | 主に検証範囲 |
| U6 | ◎ / ◎ / ○ | 複雑さ: actions.py が 500 行の上限にあり、モードの表現が 2 つ共存する構造。上げる具体的な改良はない | — |
| U7 | ◎ / ○ / ○ | 有用性・複雑さ: C20、C29、C44（v3 のまま） | 小さなコード |
| U8 | ◎ / ○（上側）/ ○ | 有用性: C25 + C43。複雑さ: C37 + C27′ | 小さなコードと計算数分 |

### 5.4 進め方と、改良を続ける価値

- **必須の改良はもうありません。** 3 ラウンドで must・should はすべて片づき、△ と × もなくなりました。残る課題は medium 1 件（コストだけ）と low だけです。今回の M9・S26・C38 で見つかった新しい課題も、すべて検証記録の書き方の不正確さか、以前からの軽い性質です。
- **同じ形のラウンドをもう一度回す必要はありません。** 実装 → 実計算での再検証 → 9 単位の再評価という重い手順は、結論を誤らせる課題があるときに見合うものです。今の候補はどれも could で、効果は文書の正確さか、評価の留保を 1 段外す程度です。
- **続けるなら、次の順が釣り合います。**
  1. **C45 + C39′（文書だけ、1 コミット）**: 実計算は要りません。検証記録と design.md が実装に一致し、利用者が再開の挙動や code_version を誤読しなくなります。
  2. **小さなコードをまとめて 1 回**: C25・C28・C40・C41（+ C14′）と、【複雑さ減】の C37・C27′・C34。どれも数行〜数十行で、確認は pytest と、C34 だけ実 CREST のプローブ（約 2 分）で足ります。U8 の複雑さ、U2 の但し書き、U0 のコード側の理由が片づきます。
  3. **検証の空白を埋める実計算（任意）**: C43（手法パネルの正式な実行、計算数分）で U8 の有用性が ◎ になります。評価の上限を決めているのは、コードよりも実計算で通っていない経路（§2.8。行 14、負モード 2 本以上、電荷のある系・開殻系）です。
- C16 は前提の M9 が満たされましたが、【複雑さ増】で、fake の E2E テスト（FIND_PATH → multi_max → VALIDATE_INTERMEDIATE）を同時に入れることが条件です。用意できなければ見送ってかまいません。

### 5.5 あえて行わないこと

v3 §5.5 の方針はすべて維持します（既定の流れを重くする変更、キーとスキーマを広げる変更、knob や閾値を増やす変更、根拠のない性能改善・再試行、表示だけのための配管、対象外の化学への対応）。第3ラウンドの評価で新たに見送ったものを、同じ分類で加えます。

- **キーとスキーマを広げる変更**: code_version を config_sha に入れて done の stage を計算し直す。stage ごとの code_version を run_state.json に記録する。コミットのたびに再開と JobStore の価値がなくなるか、スキーマが広がります。版をまたいだ再開の誤読は、C45 の 1 文で防げます。
- **knob や閾値を増やす変更**: end_drift を種の採否のゲートにする、「0.004 Å 以下」の目安を締める。M9 の後のずれは解像度の 1/100 で、実測は小分子 2 系だけです。
- **根拠のない性能改善**: align_sequential の基準を NWChem の L1 に合わせる、NWChem の string の algorithm を固定する、HEI の補間構造で SP を取り直す。どれも効果が揺らぎの範囲か、頑健さを損なうおそれがあります。
- **効果を測れない変更・テスト**: git に `--no-optional-locks` を付ける、.gitattributes を置く、C38 の argv をモックで固定する専用テスト、fake の上での接線のテスト。実 run と実エンジンのプローブで確かめた挙動を、モックで固め直すだけになります。
- **表示だけのための配管と定期的な書き換え**: report.html と ranking.csv に code_version を出す、README のベンチマーク表を毎ラウンド更新する。README は「v3」と明記していて正確で、v4 は validation.md §11 にあります。
- **重い実計算の追加**: M9 の U6 への効果を確かめるための大きな錯体の run、U0 の目的での method_panel・composite の実行。前者は該当する case がなく数時間かかり、後者は記録が validation.md §3 に入りました（正式な実行は C43 の範囲）。
- **不要な削除**: string デッキの endgeom と impose。使われないので削っても挙動は変わらず、文書で実態を書けば足ります。
