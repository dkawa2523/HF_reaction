# hfauto 計算処理の仕様 再評価（v2）— 改良 W1〜W4 の後の評価・残る課題・次の改良（2026-09-26）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `56d9cb1`。コードは `4e972c5` と同じです（`56d9cb1` は README.md と docs/validation.md だけを変更）。読み取り専用のレビューで、リポジトリのコード・設定・既存の文書は変更していません（追加したのはこのファイルだけ）。
- ベースライン: `docs/reviews/2026-09-26_calculation_spec_review.md`（HEAD `92ce8d6`、以下「前回」）。その §2.x の課題、§3 の X1〜X12、§4 のロードマップ（must M1〜M7、should S1〜S17、could C1〜C27）を基準にしました。
- 評価した改良:

| コミット | 内容 |
|---|---|
| `116ee0e`（W1） | M1〜M5、M7、S2〜S6、S8、S10〜S13、C5、C7 |
| `9b61e8d`（W2） | M6、S9、S14、S15、S17、C23 |
| `7b9a192`（W3） | S1、S7、C19 |
| `4e972c5`（W4） | 文書: S16、C3、C6、C9、C11、C24、X11 |
| `56d9cb1` | docs/validation.md v2（WSL での実計算による再検証） |

- 方法: 前回と同じく、計算単位 U0〜U8 ごとに 3 段で調べました。
  1. **現状仕様の整理（変更点が中心）**: コード、設定、改良前後の実 run の生ファイルを読みました。改良前は W7（コード `6b87814`、WSL `/home/user/hfauto_w7/<run>`）、改良後は v2（コード `4e972c5`、`/home/user/hfauto_v2/<run>`）です。
  2. **評価**: 化学的妥当性、有用性、複雑さを評価し、ライブラリの文書とソースを調べました。必要に応じて小さなプローブ計算を行いました（WSL `/home/user/hfauto_review2_probe/`）。プローブの多くは他のプローブと同時に走り、load 8〜15 の状態だったので、所要時間は参考値です。
  3. **反証的な検証**: 別の担当が 2 の主張をコードと生データで確かめ、状態・重大度・優先度・評価を確定しました。
- 記載の規則:
  - 課題の状態、重大度、優先度、評価は、検証段の判定を使いました。
  - 全面的に棄却された課題はありませんでした。一部だけ正しいと判定された課題（U2-N3、U3-N4、U4-N3、U6-N2）は、訂正後の内容で書きました。
  - 改良案は、検証で「採用」または「修正して採用」とされたものだけを載せ、修正して採用されたものは修正後の形で書きました。評価段で見送った案は、各節の「見送った案」に理由付きで載せました。
  - 複数の単位で重複した改良は、§5 で 1 つにまとめました。前回のロードマップと同じ内容のものは前回の ID（C2、C8 など）を使い、新しいものには M8、S18〜S25、C28〜C39 を振りました。
  - ID の付け方: 前回の課題は前回の ID（U0-I1、U2-ISS-1、U8-1 など）のまま。今回新しく見つかった課題は `U<n>-N<k>`、単位ごとの改良案は `U<n>-R<k>` です。

凡例:

| 記号 | 意味 |
|---|---|
| critical / high / medium / low | 課題の重大度 |
| must / should / could | 優先度。must は「結論を誤らせる」または「主要機能が働かない」もの |
| 【複雑さ減】/【複雑さ増】 | コード・設定・分岐が減る / 増える改良 |
| ◎ / ○ / △ / × | 評価。◎ 問題なし、○ 軽微な課題のみ、△ 改良が必要、× 目的どおりに機能していない。「○−」は ○ の下限 |
| resolved / partially / unresolved / not_in_scope | 前回の課題の状態（解消 / 一部解消 / 未解消 / 本単位の範囲外・方針により対応しない） |

エネルギーの単位は、断りがなければ kcal/mol です。「W7 → v2」は改良前 → 改良後の実 run を指します。

---

## 0. 要旨

### 0.1 全体評価

1. **改良は計画どおりに入りました。** 前回の must 7 件と should 17 件はすべて実装され、could は 27 件中 9 件（C3、C5、C6、C7、C9、C11、C19、C23、C24）が実装されました。前回「結論を黙って誤らせる」とした 5 件はコードから消えました。
   - negative_evidence の拒否権（M1）
   - FIND_PATH の成り立たない収束判定と停滞 kill（M2）
   - 手法パネルの許容幅 0 の符号ゲート（M5）
   - SMILES で与えた宣言反応の端点（M7）
   - 熱化学が小さな負の振動数を捨てること（M6）
2. **評価は全単位で同じか、良くなりました（§1）。** △ と × は、U5 の「2 チャンク以上の DFT string で monotonic → BARRIERLESS」という条件付きの 1 か所（M8 まで △）を除いて、なくなりました。
   - U5 の FIND_PATH: × → ○−（初めて実計算で機能した）
   - U4 探索: 妥当性・有用性 △/△ → ○/○
   - U3 極小: ○/○/△ → ◎/◎/○
   - U7 熱化学と U8 報告の妥当性: ○ → ◎
   - U2 配座の妥当性: △ → ○、複雑さ ○ → ◎
3. **実計算での効果**（validation.md v2 と本レビューのプローブ）
   - TMA·(HF)₂ discover: 52:33 → 31:15（−41%）。dft stage 3,050 → 1,771 s（−42%）。主因は M3（偽の AFIR 生成物が消えた）。
   - A/B arm A の dft stage: 3,317 → 2,433 s（−27%、S1 で AFIR 生成物が既知の basin に入り 17 原子の freq を省いた）。NH3・水の DFT freq は 2 → 1 本。
   - HONO の QRC: 143 → 70 ステップ、311 → 162 s（S7）。
   - (HF)₃ の CREST: `--noopt` で 2/2 成功（4.8 s / 19.2 s。前回のプローブは 87〜236 s で 1 系失敗）。TMA·(HF)₃ の最安はイオン対。
   - FIND_PATH: SCREEN の GS が失敗した HONO を救い、W7 と同じ ΔE‡ 13.67 の TS に到達した。
   - ΔG の変化は S15（スケール 0.989/0.975 → 1.0）だけによるもので、ΔG‡ は −0.04 以内、ΔG_assoc は +0.07。
   - 手法パネル（プローブ）: 既定の PBE0/TZVPD と ωB97X-D3/TZVPD を HONO に追記しても rank 1（ΔG‡ 12.22）が残り、report と panel_report の順位が一致した。HCN の既定パネルは 6 SP・36 s（前回の 4 手法は 12 SP・63.7 s）。
4. **新たに見つかった主な課題**
   - **(high) NWChem ZTS の凍結端点のずれ（U5-N1 = U6-N3）**: `freezeN` の bead の座標は NWChem が動かします。FIND_PATH はチャンクを継ぐとき、ずれた端点をそのまま次のチャンクに渡します（STO-3G の HCN で、+8.85 kcal/mol の非極小を凍結端点として計算した）。S3 で monotonic な string 1 本を BARRIERLESS の十分な証拠にしたので、誤分類の余地があります。→ M8。
   - **(medium) IDPP の打ち切り（U5-N2）**: 反復 500 で黙って打ち切られ、平面の端点（HONO）では面内の反転経路になります。FIND_PATH の HEI が 34.4（真の TS は 13.7）になり、鞍点の精密化の手間が約 10 倍になりました。
   - **(medium) moddir 0 の迷走（U6-N1）**: 非停留の種に対し、NWChem の saddle が「Cartesian で moddir 0（既定）」の組合せのとき、初手 uphill の後に正の曲率のモードを追い、+82 kcal/mol まで登ってから戻りました。
   - **(medium) DLC の GS の脆さ（U5-N3）**: pysisyphus の DLC の GS は、平面 4 原子で 1e-6 Å の差によって成否が変わります。S4 の再実行は決定的に同じ失敗になります。
   - **(効率) S7 を trust 0.1 のまま入れたこと（U6-I1 の残り）**: TS Hessian の負の固有値のため、刻みが 0.03 Å に制限されます。HCN は効果なし、HONO の速い側は 17 → 36 ステップ。trust 0.3 のプローブでは全側でステップ数がほぼ半分、最終エネルギーの差は 6e-8 Eh 以下でした。
   - (low) 主対象（TMA·(HF)₂、水）の報告行に誤読を招く blocker が付く、composite の文書の手順がそのままでは空の結果になる、手法誤差の目安が甘い、CCSD の maxiter 20 の余裕が小さい、SMILES の同位体を黙って ¹H にする、設定値の範囲を検査しない、など。
5. **検証の空白**: 次の経路は、まだ実 run で通っていません。
   - M1、S3、行 13（中間体）、低レベル TS の近道、SCREEN の IDPP
   - higher_order、分割、QRC の再試行、_register、reassigned、periodic_nearest、moddir ≥ 1、DFT の mode-follow の全体
   - M6、S9、S14 と、電荷のある系・開殻系を structures から NWChem まで通すこと
   - 手法パネルと composite は、本レビューのプローブ（v2 の run の複製への追記）でだけ確かめました。
6. **次に有用な改良（最小限）**
   - must 1 件: **M8**（チャンクの継ぎ目と最終経路で、端点を DFT 極小に戻す。3〜8 行）。
   - should 9 件: IDPP の反復上限（S18）、FIND_PATH の画像の整列（S19）、moddir の明示（S20）、Hessian を渡す opt の trust 0.3（S21）、rankable の blocker の正確化（S22）、composite と手法誤差の文書（S23）、SMILES の同位体の拒否（S24）、stage の入力切れに上流の理由を添える（S25）、CCSD の maxiter 50（C2）。
   - どれも数行〜十数行で、knob は増やしません。

### 0.2 次の改良（上位 5 件）

| # | 改良 | 由来 | 期待効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| 1 | **M8** FIND_PATH のチャンクの継ぎ目と最終経路で、端点を DFT 極小に戻し、ずれを note に残す | U5-R1（= U6-R3 の後半） | 2 チャンク以上の string でも monotonic → BARRIERLESS の根拠（minimax の上界）が成り立つ。STO-3G HCN の +8.85 kcal/mol のような非極小の凍結がなくなる。追加の計算はない | 3〜8 行とテスト約 15 行 | 中立 |
| 2 | **S18** IDPP の反復上限を 500 → 5,000 | U5-R2 | 平面の端点で FIND_PATH がねじれの機構を捉える（HONO の DFT//IDPP の最大 37.0 → 14.8）。U6 の DFT Hessian と AUTOZ の再試行（HONO で約 110 s）を避けられる見込み | 定数 1 つとテスト 1 件 | 中立 |
| 3 | **S20 + S19** 鞍点で追うモードを明示する（常に moddir、負モード 2 本以上は Cartesian と P·H·P の順）。FIND_PATH の画像を逐次整列する | U6-R1、U6-R3 の前半 | HONO の path_hei の種で saddle 70 → 18 ステップ（プローブ）、登り返しがなくなる。STO-3G の重なり 0.118 → 0.651 で、誤った note と xTB Hessian の誤った棄却がなくなる | 約 15 行 + 1 行、テスト 2 件 | 中立 |
| 4 | **S21** 初期 Hessian を渡す optimize は trust 0.3 | U6-R2 | QRC のステップが HCN 44 → 27、HONO 70 → 29、NH3（片側）45 → 24（プローブ）。v2 の QRC 合計 約 266 → 約 110 s | 1 行と文書、再検証 約 1 時間 | 中立 |
| 5 | **S22** rankable で、順位を付けられない outcome には outcome の理由と本当の blocker だけを返す | U8-R1 + U0-R5 | 主対象の報告行が `outcome:same_basin` だけになり、ΔG_assoc −11.53 の隣に `thermo_unavailable` が出る矛盾がなくなる（C26 を置き換える） | 3〜4 行とテスト 1 行 | 中立 |

### 0.3 今回の改良で良くなり、今のままでよい点（主なもの）

- M1: 低レベルの陰性結果は判定に使わず、coverage.csv の件数としてだけ報告する。
- M2: 打ち切りは bead エネルギーの安定だけで決め、経路の形は収束を問わず記録する（連続な DFT 経路の最大値は鞍点の上界）。
- M3 + S6: explore は screen の最適化構造から出発し、生成物の同一性は状態ラベル（結合グラフ）で判定する。v2 の 71 attempt すべてで ReaDuct の source の緩和が 2 反復・ΔE 0.00 kJ/mol になった。
- M4: 組成の CREST に `--noopt`。CREST 3.0.2 のソースで、初期トポロジー検査（trialOPT）そのものが呼ばれないことを確認した。
- M5 + S8: 符号ゲートの ±1 kcal/mol の不感帯（モジュール定数）と、diffuse 付きの TZVPD パネル。CCSD(T) は 1 語の opt-in。
- M6 + S15: 負モードの固定規則（純関数 `thermo_frequencies`）と、スケール 1 本（`vib_scale`）。QRC・SaddleClaim の「最低モード＝反応座標」の約束と一貫し、GoodVibes の U_vib の ZPE と報告される ZPE がそろった。
- M7、S10〜S12、C5: 誤った入力と設定を、読み込み時か structures で計算前に止める。検査の置き場所（記述だけで分かるものは読み込み時、原子数・電子数が要るものは structures）も明確。
- S1: その場の登録と同一性の 1 基準（0.05 Å / 5e-5 Eh。CREGEN の既定より厳しい）。1 basin に freq 1 本。
- S9、S14、S17: energy_method、会合量、開殻 WFT を fail-closed にした。
- 設定の量: StageConfig の knob 59 → 51、method ファイル 6 → 5、ThermoSettings 8 → 5 フィールド、ConformersConfig の利用者が書ける項目 9 → 4。hfauto 全体 10,707 → 10,651 行。

---

## 1. 評価の比較表

評価は「妥当性 / 有用性 / 複雑さ」の順です。複雑さの ◎ は「単純で過不足がない」を意味します。

| 単位 | 前回の評価 | 今回の評価 | 主な変化 | 実測の変化 W7 → v2 |
|---|---|---|---|---|
| U0 全体プロトコル | ○ / ○ / △ | ○ / ○ / ○ | M5 不感帯、S8 TZVPD パネルと CCSD(T) の opt-in、S9 fail-closed、S10 読み込み時の検査、S15 スケール 1 本、C3・C24 の文書 | ΔE は完全一致（HCN 46.62、HONO 13.67、NH3 4.26）。ΔG‡ −0.01〜−0.04、ΔG_assoc +0.07（S15）。パネル: HONO の rank 1 が残る（プローブ） |
| U1 structures | ○ / ◎ / ◎ | ○ / ◎ / ◎ | M7 端点は xyz 必須、S11 宣言座標の検査、S12 静的な検査、C5 未知元素の fail-closed、C6 文書 | 出力（指紋・組成 id・状態ラベル・電荷・多重度）は W7 と完全一致。範囲外の添字は「DFT の後の IndexError」から「structures で 1.56 s」に（プローブ） |
| U2 conformers | △ / ○ / ○ | ○ / ○ / ◎ | M4 組成に `--noopt`、S13 窓の削除、C7 knob 9 → 4 とスレッドの一本化、C9 文書 | tma_hf2 の stage 17 → 19 s（組成の CREST 10.9 → 12.1 s、2 → 3 配座）。(HF)₃ は 2/2 成功（4.8 / 19.2 s） |
| U3 minima | ○ / ○ / △ | ◎ / ◎ / ○ | S1 その場の登録と同一性の 1 基準、tight・rejudge・Verdict の削除。M3・M6 の波及 | dft stage 3,050 → 1,771 s（TMA·(HF)₂）、A/B 3,317 → 2,433 s（−27%）。NH3・水の DFT freq 2 → 1 本。screen の xTB freq 8 → 4 ジョブ |
| U4 explore | △ / △ / ○ | ○ / ○ / ○（有用性は下限） | M1 拒否権の削除、M3 出発構造、S6 状態ラベルでの一致判定 | 偽の生成物 1 → 0。source の緩和 −2.9〜−33.2 kJ/mol → 0.00。attempt 37 → 31（tma）、44 → 40（amine） |
| U5 paths（仮説・SCREEN・FIND_PATH） | SCREEN ○ / FIND_PATH × / △ | ○（FIND_PATH ○−）/ ○ / ○ | M1、M2 bead エネルギーでの打ち切り、S2 GS 画像の整列、S3 monotonic 1 本、S4 GS の回収、S5 FakePath | FIND_PATH が初めて機能（HONO 1 チャンク 200 s、STO-3G 2 チャンク）。HONO の U5 区間は GS の失敗で 48 → 208 s |
| U6 paths（鞍点・検証・接続） | ○ / ◎ / ○ | ○ / ◎ / ○ | S7 QRC に TS の Hessian、C19 振幅の上限と drop の一本化。M2 で path_hei の種が実際に入る | QRC: HONO 143 → 70 ステップ（311 → 162 s）、NH3 100 → 90、HCN 43 → 44。HONO の saddle は path_hei の種で 10 → 121 s |
| U7 thermo | ○ / ○ / △ | ◎ / ○ / ○ | M6 負モードの固定規則、S15 スケール 1 本、C23 温度は conditions だけ、S14 会合量のキー、S9 | ΔG‡ −0.01〜−0.04、ΔG_assoc +0.067（すべて S15）。worker 0.23〜0.28 s/subject、失敗 0 |
| U8 sp・report・実行基盤 | ○ / ○ / ○ | ◎ / ○ / ○ | M5、S8、S17 WFT は閉殻専用、knob 59 → 51、method ファイル 6 → 5 | report 0〜1 s（変化なし）。パネル（プローブ）: HCN 6 SP・36 s、HONO で report と panel_report の rank 1 が一致 |

評価の留保（検証で付いたもの）:

- **U0 有用性 ○**: composite の文書の手順（design.md:150）がそのままでは動かず、手法誤差の目安（約 1.5）が甘いことが前提です。S23 で解消します。
- **U1 妥当性 ○**: SMILES の同位体を黙って ¹H にする経路（S24 で閉じる）に加えて、電荷のある系・開殻系を structures から NWChem まで通した実 run がないこと、組成の多重度が高スピン固定であること（C8）も理由です。
- **U2 複雑さ ◎**: 組成では到達せず、単量体では実 CREST で失敗する `--noreftopo` の分岐が残っています。C34（U2-R1）で名実ともに ◎ になります。
- **U3 妥当性 ◎**: DFT の mode-follow は、実エンジンでは「TS に留まる opt → freq の saddle 判定 → +側の opt」までしか通っていません（−側、ts_candidate、両側の再 settle、discovery、DFT の soft 押し、M6 の |ν| 化は未通過）。
- **U4 有用性 ○**: 下限です。対象系の実 run で採用した生成物は 0 件、NT2 の slot の約 6 割が同じイリド TS、陰性結果は件数しか残りません。C12・C14′ で安定します。
- **U5 妥当性**: SCREEN は ○、FIND_PATH は ○−。M8 を入れるまで、2 チャンク以上の string による monotonic → BARRIERLESS は △ です。M8 の後も、1 チャンク内の端点のずれは note による監視になります。
- **U6 妥当性 ○**: U6-N1 が起きるのは「Cartesian で moddir 0」の組合せに限られ、重なりの過小評価（U6-N2）の主因は string 画像の未整列です。どちらも VALIDATE_TS と QRC のゲートが結果の誤りは防いでいます。

---

## 2. 前回の課題の解消状況

### 2.1 課題ごとの状態

前回の課題 78 件の内訳は、resolved 44、partially 7、unresolved 22、not_in_scope 5 です。前回 critical・high とした課題（X1、U1-01、U2-ISS-1、U3-ISS-1、U5-I2、U5-I3、U0-I1、U8-1）は、すべて resolved です。unresolved の 22 件は、ほとんどが前回 could か「方針により対応しない」としたものです。

| ID | 単位 | 状態 | 根拠 |
|---|---|---|---|
| U0-I1 符号ゲートに許容幅がない | U0 | resolved | `summary.py:32` の `_SIGN_DEADBAND_KCAL = 1.0`、211-212。単体テスト test_summary.py:92-99。プローブで HONO の既定パネルは +0.13 / +0.24 / +0.33 で sign_disagreement=False、rank 1。実データで解消させたのは S8 で、不感帯の分岐は単体テストだけで確認 |
| U0-I2 パネル参照に diffuse がない | U0 | resolved | configs/methods から SVP・TZVP・CCSD(T)/TZVP を削除し、TZVPD の 2 本を追加。method_panel.yaml:9 は TZVPD の 2 手法。tests/smoke にインラインの `MethodSpec(id="pbe0-d3bj_def2-svp")` が 2 か所残るが、config を参照しないので実害なし |
| U0-I3 composite が配線されていない | U0 | partially | S9 と文書（design.md:150）は入った。composite.yaml と blocker は未実施（C1 は見送り）。文書の手順は新しい run-dir では `stage 'comp_sp' (sp) has no input of type ['minimum']` で止まる（本レビューで再現、U0-N2） |
| U0-I4 energy_method が fail-open | U0 | resolved | thermochemistry.py:59, 81, 96-100。SP がなければ GoodVibes を呼ばず G=None、`energy_layer_missing`。ΔE だけは freq 層に戻る（U0-N1） |
| U0-I5 停留点 method の一致をテストでしか見ない | U0 | resolved | `PipelineConfig._check_stages`（config.py:66-78）。test_config_load.py:41-47 |
| U0-I6 amine·(HF)n は PBE0 の PES 上の結論 | U0 | resolved | design.md:149 に明記 |
| U0-I7 同順位の幅に手法誤差が入らない | U0 | resolved | design.md:149 と README:86 に「約 1.5 未満は手法誤差」。数値は TZVP 参照に基づく（U0-N3） |
| U0-I8 CCSD の maxiter 20 | U0 | unresolved | input.py:157-167 は maxiter を書かない。HONO の TS/TZVP は 18/20 反復 |
| U1-01 SMILES の端点で H の写像が決まらない | U1 | resolved | core/system.py:88-91。CLI は exit 2 で run-dir を作らない（0.79 s） |
| U1-02 宣言座標の添字を検査しない | U1 | resolved | 形は system.py:48-55（読み込み時）、上限は structures.py:88-89。範囲外は 1 s 未満・DFT ジョブ 0 で止まる |
| U1-03 SMILES の同位体・ラジカル | U1 | unresolved | smiles.py は無変更。'[2H]O[2H]' → H2O_q0_m1（水と同じラベル）、'[13CH4]' → CH4。xyz の 'D' は拒否されるので入力経路で扱いが違う |
| U1-04 id の衝突で黙って上書き | U1 | resolved | system.py:27-31, 37, 58-59, 71-75。プローブで 5 件とも拒否 |
| U1-05 未知元素で検査を飛ばす | U1 | resolved | electronic_state.py:24-26（spec の 58-61 は誤り）。xyz の 'D' → 'unknown element(s): D' |
| U1-06 状態ラベルの限界 | U1 | unresolved | topology.py は無変更。S6 でラベルが生成物の棄却（same_as_source）にも使われるようになり、「グループ化にしか使わないので対応不要」という前回の根拠は古くなった（U1-N5） |
| U1-07 structures の規則が文書にない | U1 | resolved | design.md:138（C6）と README:41 |
| U1-08 組成キーに多重度がない | U1 | resolved | U7 の S14 で単量体キーが (Hill, 電荷, 多重度) |
| U2-ISS-1 組成の CREST がトポロジー停止する | U2 | resolved | crest.py:56-58 の `--noopt`。v2 の (HF)₃ は 2/2 が rc 0、noreftopo_retry=false、topology_removed=0 |
| U2-ISS-2 4 kcal/mol の窓 | U2 | resolved | conformer_search.py:118-130 に窓がない。window_kcal は ValidationError |
| U2-ISS-3 組成の多重度を黙って高スピンに | U2 | partially | 文書化（design.md:137）だけ。動作は同じ（C8 は未実施） |
| U2-ISS-4 効かない knob | U2 | resolved | ConformersConfig は 4 knob。スレッドは site の 1 か所 |
| U2-ISS-5 `--nci --quick` の実効設定が不明 | U2 | resolved | design.md の conformers 行（C9）。v2 の stdout と一致 |
| U2-ISS-6 seed00 だけを渡すこと・衝突条件 | U2 | resolved | design.md に記載 |
| U2-ISS-7 組成の停止構造が流れる | U2 | resolved | `--noopt` で組成の停止そのものが起きない |
| U3-ISS-1 同じ stage で同じ basin の 2 本目の freq | U3 | resolved | minima.py:92-98, 217-231。A/B で tma_hf2_afir が known（−27%）、NH3・水の freq 2 → 1 |
| U3-ISS-2 explore が入力構造から出発 | U3 | resolved | M3（explore.py:148-150）。dft の always 生成物 0 件 |
| U3-ISS-3 thermo が noise / soft の負モードを捨てる | U3 | resolved | M6。実 run に負の副次モードはなく、実データでの発火はまだない |
| U3-ISS-4 同一性の基準が 2 種類 | U3 | resolved | compare_minima・rejudge・ambiguous・tight は hfauto/ に残っていない |
| U3-ISS-5 window で screen 極小 0 件なら黙って空 | U3 | unresolved | minima.py:199-210 は前回と同じ（C10 は未実施） |
| U3-ISS-6 振動数の射影の約束 | U3 | resolved | design.md:102（C11） |
| U3-ISS-7 Registry.find が組成を見ない | U3 | partially | add は composition_id・level_key で絞る。relax 中の known.find は絞らない。エネルギー条件で分離されるので「変更不要」を維持 |
| U3-ISS-8 開殻一重項を RKS で扱う | U3 | resolved | design.md の structures 行（C11） |
| U4-I1 陰性結果の拒否権 | U4 | resolved | hfauto/・configs/・tests/ に negative_evidence・NO_PRODUCT・override_negative の残りなし。reaction_paths.py:57 は旧キーを拒否。実 run ではこの経路に入っていない |
| U4-I2 出発構造 | U4 | resolved | v2 の 71 attempt すべてで source の opt が 2 反復・ΔE 0.000 kJ/mol（W7 は 14〜38 反復、−2.9〜−33.2） |
| U4-I3 RMSD による一致判定 | U4 | resolved | worker.py:49-51 は状態ラベルの一致。閾値近くの H 結合によるラベル反転のリスクは U4-N2 |
| U4-I4 等価な drive で上限が埋まる | U4 | unresolved | C12 は未実施。v2 で TMA 9/10、TMA·HF 9/10、TMA·(HF)₂ 6/10 の slot が同じ C→N 1,2-H 移動 |
| U4-I5 ReaDuct の反復上限 150 | U4 | unresolved | C13 は未実施。IRC 66/66 本が 150 で打ち切り、amine の irc_end_not_minimum 8 件中 6 件は端点の opt の上限 |
| U4-I6 陰性ラベルが粗い | U4 | unresolved | C14 は未実施。前半（NT2 の失敗の分類）は情報を増やさないと判明、後半（TS がある negative の ΔE‡）が残る |
| U4-I7 AFIR γ=300 の害 | U4 | not_in_scope | 前回に棄却済み。v2 でも害の実例なし |
| U4-I8 F→N 移動の drive が作られない | U4 | not_in_scope | 設計上の判断。trials.py と topology.py は差分なし |
| U5-I1 negative_evidence ゲート | U5 | resolved | state.py の行 12 は barrierless だけ。test_reaction_paths_stage.py:110-128。実 run ではこの経路に入っていない |
| U5-I2 FIND_PATH の収束判定 | U5 | resolved | actions.py:483-485 は形を収束と無関係に記録。v2 の HONO と STO-3G は program_converged=False のまま single_max で TS に到達 |
| U5-I3 収束した string を停滞として kill | U5 | resolved | nwchem/engine.py:170-171 の monitor は None |
| U5-I4 GS 画像の回転の跳び | U5 | resolved | interpolation.py:32-39、actions.py:242。W7 画像のオフライン再計算で 1.149 → 0.190 Å（v2 の HONO は GS 失敗で実 run では未確認） |
| U5-I5 TSOpt の例外で GS を失う | U5 | resolved | worker.py:27-34。プローブ（cart）で ts=None と tsopt_error で回収。GS 自体の失敗でも再実行する副作用は U5-N4 |
| U5-I6 15 bead の確認 | U5 | resolved | state.py:196-200。ただし monotonic → BARRIERLESS の健全性は、端点が真の極小であることが前提で、M8 まで条件付き |
| U5-I7 両端より低いノード | U5 | unresolved | C16 は未実施 |
| U5-I8 tangent 用の DFT freq | U5 | unresolved | C15 は未実施。v2 の NH3 で縮退端点（原子の置換）のため JobStore がミス |
| U5-I9 DFT//xTB の端点のずれ | U5 | not_in_scope | 前回に「変更不要」 |
| U5-I10 ZTS の tol の意味 | U5 | resolved | design.md:144 |
| U6-I1 QRC に Hessian を渡さない | U6 | partially | S7 は入った（HONO 143 → 70 ステップ）。trust 0.1 のままで、負モードの刻みが 0.03 に制限される（HCN 43 → 44、HONO の速い側 17 → 36） |
| U6-I2 mode_index と moddir の番号 | U6 | unresolved | input.py:140-141 は無変更（C17 は未実施）。mode_index=0 では moddir を書かない |
| U6-I3 periodic_nearest に許容幅がない | U6 | unresolved | C18 は未実施。実 run で到達した例なし |
| U6-I4 QRC の振幅のクリップ順と drop の重複 | U6 | resolved | C19。gates.qrc_drop を共有、振幅は min(first×2^(n−1), 0.4) |
| U6-I5 _register がゲートの前 | U6 | unresolved | 設計どおり。実 run で発動なし |
| U6-I6 torsional の閾値 | U6 | not_in_scope | hypotheses の担当 |
| U6-I7 試行を使い切った後の FIND_PATH | U6 | unresolved | 行 16 は無変更。S3 で帰結が変わり、無駄になるのは single_max の場合だけ |
| U7-I1 負モードを捨てる | U7 | resolved | M6。GoodVibes 実物のプローブで G(−x) − G(+x) = 0.0000（旧来は +0.74〜+1.61） |
| U7-I2 会合量の LOT とキー | U7 | partially | S14 は入った。ΔG_assoc=None の理由が出ない。単原子の単量体（F⁻）は常に None |
| U7-I3 経路縮重度 | U7 | resolved | design.md:146、validation.md:24 の ※1 |
| U7-I4 zpe_scale が G に入らない・借用のスケール | U7 | resolved | S15。GoodVibes の U_vib の ZPE と報告値がそろった |
| U7-I5 ΔG_assoc の感度幅 | U7 | unresolved | C21 は未実施。プローブで幅 1.155（298 K） |
| U7-I6 QH=False の未記載 | U7 | resolved | design.md:146 |
| U7-I7 点群・σ の記録 | U7 | unresolved | C20 は未実施 |
| U7-I8 基準が反応物の極小そのもの | U7 | unresolved | C22 は方針どおり未実施（HONO で約 +0.37） |
| U7-I9 効かない設定・書くだけの出力 | U7 | partially | temperatures_K は削除。population・H・S_rot は書くだけのまま。ThermoResult.notes は常に空 |
| U7-I10 トンネル補正なし | U7 | resolved | design.md:146 に明記 |
| U7-I11 単原子種 | U7 | unresolved | 方針どおり。Evidence.n_external が Literal[5, 6] |
| U8-1 符号ゲートで HONO が外れる | U8 | resolved | U0-I1 と同じ。プローブで panel_report と report の ranking.csv が一致（rank 1・12.2237） |
| U8-2 既定パネルの CCSD(T) | U8 | resolved | S8 で opt-in |
| U8-3 CCSD の maxiter | U8 | unresolved | U0-I8 と同じ |
| U8-4 UHF-MP2 が常に失敗 | U8 | resolved | S17。プローブで (1,2)・(0,3) ともジョブ 0 件で拒否 |
| U8-5 composite と thermo の上書き | U8 | partially | S9・C24 は入った。文書の手順が動かない（U8-N1 = U0-N2）。thermo の artifact id に stage id がない点は文書の注意だけ |
| U8-6 同順位の幅の意味 | U8 | resolved | README:86、design.md:149 |
| U8-7 ranking.csv に T・標準状態・ΔG_rxn がない | U8 | unresolved | C25 は未実施。metric=dG_rxn のとき並べた値が CSV に出ない点も加わった（U8-N5） |
| U8-8 TS のない outcome の dzpe blocker | U8 | unresolved | C26 は未実施。範囲は thermo_unavailable にも及ぶ（U8-N3） |
| U8-9 site.memory_mb が未使用 | U8 | unresolved | C27 は未実施 |
| U8-10 method_panel.yaml が独自の conditions | U8 | unresolved | 前回に受容。失敗は thermo_unavailable として明示される（参照行は report.py:47-48） |
| 2.8.3(e) hfauto report の二重出力 | U8 | not_in_scope | 前回に見送り。内容の矛盾は U8-1 の解消でなくなった |

横断課題（前回 §3）:

| ID | 状態 | 根拠 |
|---|---|---|
| X1 陰性結果の拒否権 | partially | ゲートの削除（本体）は resolved。同時に推奨した「TS がある negative にも ΔE‡」（C14 の後半）は未実施 |
| X2 explore の出発構造と一致判定 | resolved | M3・S6。v2 で偽の生成物と配座違いによる非接続判定は 0 件 |
| X3 小さな負の振動数 | resolved | M6 |
| X4 手法パネル | partially | M5・S8・文書は実装済み。C2（CCSD の maxiter）だけが残る |
| X5 composite G | partially | S9・C24 は実装済み。C1 は見送り、文書の手順に穴がある（S23 で直す） |
| X6 同一性の統一と stage 内の freq 省略 | resolved | S1。発火した A/B で −27% |
| X7 計算時間の大部分への対応 | partially | freq の重複（S1）、偽の生成物（M3）、15 bead（S3）、GS の回収（S4）、既定の CCSD(T)（S8）、CREST の停止（M4）は解消。QRC は trust 0.1 のため一部（S21）、tangent の DFT freq（C15）は未実施 |
| X8 入力と設定の誤りは読み込み時に止める | partially | S10・M7・S11・S12・C5・C7・C23 は実装。xyz ファイルの存在、設定値の範囲、C10、C27 が残る |
| X9 開殻と多重度 | partially | S14・S17 は実装。C4（同位体）と C8 が残る。開殻の組成は `--noopt` の CREST で失敗しうる（U2-N2） |
| X10 fake でしか通っていない経路 | partially | FIND_PATH、IDPP、path_hei の種、種での DFT freq、autoz 後の Cartesian 再投入、CREST の組成（`--noopt`）は実計算で通った。§4.5 の経路は未通過 |
| X11 文書化 | resolved | design.md の各行に反映。数値・記述の不正確さは S23・C39 で直す |
| X12 キャッシュの無効化 | resolved | optimize のキーから tight が消え、thermo のキーに振動数と温度が入った。validation.md の値も記録し直した |

### 2.2 must・should・could の実装状況

| ID | 内容 | 実装 | 実計算での確認 | 備考 |
|---|---|---|---|---|
| M1 | negative_evidence ゲートの削除 | 済 | 未（経路に入る系がなかった） | fake テストのみ |
| M2 | FIND_PATH の判定の置き換え | 済 | HONO、STO-3G HCN | チャンク継続の端点のずれ（U5-N1）が新たに判明 → M8 |
| M3 | explore の出発構造 | 済 | tma_hf2、amine | 偽の生成物 1 → 0 |
| M4 | CREST の組成に `--noopt` | 済 | (HF)₃、tma_hf2、smoke | 開殻の組成は失敗しうる（U2-N2） |
| M5 | 符号ゲートの不感帯と SVP の除外 | 済 | プローブのみ（不感帯の分岐は単体テスト） | HONO は S8 だけでも解消 |
| M6 | 熱化学の負モードの固定規則 | 済 | 実 run に負の副次モードなし。GoodVibes のプローブで確認 | — |
| M7 | 宣言反応の端点は xyz | 済 | プローブ（CLI exit 2） | — |
| S1 | その場の登録と同一性の統一 | 済 | A/B、NH3、水、screen | — |
| S2 | GS 画像の逐次整列 | 済 | オフライン（W7 の HONO） | — |
| S3 | 15 bead 確認の撤廃 | 済 | 未 | M8 が前提 |
| S4 | TSOpt 例外時の GS 回収 | 済 | HONO（GS 自体が失敗）、プローブ（cart） | GS 自体の失敗でも再実行（U5-N4） |
| S5 | FakePath と string の smoke | 済 | STO-3G smoke | FakePath は initial_path を無視 |
| S6 | 状態ラベルでの一致判定 | 済 | amine の TMA·HF で判定が変わった | — |
| S7 | QRC に TS Hessian | 済（trust 0.1） | HCN・HONO・NH3・STO-3G | trust 0.3 は見送られた → S21 |
| S8 | パネルの TZVPD 化と CCSD(T) の opt-in | 済 | プローブのみ | — |
| S9 | energy_method の fail-closed | 済 | プローブのみ | ΔE は freq 層に戻る（C28） |
| S10 | 停留点 method の一意性 | 済 | 単体テスト | — |
| S11 | 宣言座標の検査 | 済 | プローブ | — |
| S12 | system の静的な検査 | 済 | プローブ | — |
| S13 | conformers の窓の削除 | 済 | (HF)₃ | — |
| S14 | 会合量の LOT とキー | 済 | 未（fake のみ） | 理由が出ない（C29） |
| S15 | スケールの 1 本化 | 済 | 全 run（ΔG の変化） | 根拠の文書化（C39） |
| S16 | 熱化学の規約の文書 | 済 | — | — |
| S17 | WFT は閉殻専用 | 済 | プローブ | — |
| C1 | composite 用の追記 pipeline | 未 | — | 今回は見送り（§3.0.7） |
| C2 | CCSD の maxiter 50 | 未 | — | 今回 should に格上げ |
| C3 | 停留点レベルの根拠の文書 | 済 | — | 数値が TZVP 参照（S23） |
| C4 | SMILES の同位体とラジカル | 未 | — | 同位体の部分を S24 に。ラジカルの部分は見送り |
| C5 | 未知元素の fail-closed | 済 | プローブ | — |
| C6 | structures の規則の文書 | 済 | — | — |
| C7 | conformers の設定整理 | 済 | — | — |
| C8 | 組成の多重度が一意でなければ拒否 | 未 | — | could のまま |
| C9 | CREST の実効設定の文書 | 済 | — | environment.md:64 の所要時間が不正確（C39） |
| C10 | window で screen 極小なしなら拒否 | 未 | — | could のまま（U3-R2） |
| C11 | 射影の約束と開殻一重項の注記 | 済 | — | — |
| C12 | drive の対称性による重複除去 | 未 | — | could のまま（U4-R1、列挙数の記録を統合） |
| C13 | 無バイアス opt の上限 500 | 未 | — | could のまま |
| C14 | 陰性ラベルの改善 | 未 | — | 後半だけを C14′ として could（U4-R3） |
| C15 | tangent の遅延計算 | 未 | — | could のまま |
| C16 | 両端より低いノードの検出 | 未 | — | could のまま（M8 の後） |
| C17 | moddir の番号を driver 座標で | 未 | — | 拡張して S20 に |
| C18 | periodic_nearest の判定 | 未 | — | 前回の判断（could）のまま据え置き。実 run で到達例なし |
| C19 | QRC の振幅のクリップ順と drop の一本化 | 済 | 単体テスト（再試行は実 run で未発生） | — |
| C20 | 点群の記録 | 未 | — | 最小版を could（U7-R4） |
| C21 | ΔG_assoc の感度の幅 | 未 | — | could のまま |
| C22 | 有効障壁 ΔG‡_eff の表示 | 未 | — | 方針どおり見送り |
| C23 | thermo の temperatures_K の削除 | 済 | — | — |
| C24 | composite の説明を直す | 済 | — | 手順が動かない（S23） |
| C25 | ranking.csv の列 | 未 | — | could のまま |
| C26 | TS のない outcome で dzpe を検査しない | 未 | — | S22 に置き換え |
| C27 | memory の preflight 検査 | 未 | — | could のまま（削除ではなく検査） |

---

## 3. 単位ごとの再評価

各単位の「追加の改良案」には、§5 のロードマップ ID を【 】で添えました。

### 3.0 U0 全体プロトコルと理論レベルの選択

評価: 前回 ○ / ○ / △ → 今回 **○ / ○ / ○**（有用性の ○ は、composite の文書の手順に穴があり、手法誤差の目安が甘いことが前提。S23 で解消）

#### 3.0.1 現状仕様の要点（変更点を中心に）

- **3 層構成は変わっていません。** GFN2-xTB（探索と経路の事前判定）→ PBE0-D3BJ/def2-SVPD（grid fine、全 DFT デッキに `convergence energy 1e-7`、`disp vdw 4`）の単一 PES 上で停留点 → GoodVibes 4.3.0 の qRRHO（grimme、100 cm⁻¹、symm、1 atm から換算）。任意で手法パネルと composite G。
- 【新】**S10**: 読み込み時に `PipelineConfig._check_stages`（`pipeline/config.py:66-78`）が、minima(level=dft) と reaction-paths の method id が 1 種類であることを検査します（比べるのは id だけ）。
- 【新】**S1**: NWChem optimize の JobStore キーとデッキから `tight` が消え、停留点の opt はすべて driver の既定閾値 + `trust 0.1`（v2 HCN `jobs/2d/2d53aa6f…` の deck: `driver; maxiter 100; trust 0.1; xyz final`）。
- 【新】**S7**: QRC の両側の opt に TS の freq Hessian を `inhess 2` で渡します（HONO v2 `jobs/31/31434b81…`、key_payload.hessian に TS freq の sha）。
- 【新】**M6・S15・C23**: thermo は `thermo_frequencies`（`chemistry/thermo.py:29-33`）で負モードを固定の規則で扱い、スケールは `vib_scale` 1 本（既定 1.0、振動数と ZPE の両方）、温度は conditions だけから取ります。
- 【新】**S14**: 会合量の単量体キーを (Hill, 電荷, 多重度) にし、錯体と単量体の LOT 判定を `same_pes(state=False)` で行います。
- 【新】**S9**: `energy_method` を指定したとき、同じ構造の SP がない subject は `energy_layer_missing` で G=None（`stages/thermochemistry.py:59, 81, 96-100`）。
- 【新】**M5・S8**: パネルの既定は `[pbe0-d3bj_def2-tzvpd, wb97x-d3_def2-tzvpd]`（`method_panel.yaml:9`）。CCSD(T)/def2-TZVPD はコメントに従う opt-in（`:3-4`、目安は重原子 4 個程度）。符号ゲートは `_SIGN_DEADBAND_KCAL = 1.0`（`reporting/summary.py:32, 211-212`）。method ファイルは 6 → 5 本。
- 【新】**S17**: WFT は閉殻専用（`backends/nwchem/engine.py:317-318`）。
- 文書（W4）: design.md:147（±1 と「panel_report が最終」）、:149（SVPD の根拠、約 1.5 未満は手法誤差、ΔG_assoc は非 CP、amine·(HF)n は PBE0 の PES 上）、:150（composite は自前の pipeline を別 run-dir で）。
- 実エンジンでの呼び出し（プローブを含む）:
  - パネルの SP: ωB97X-D3 は `xc wb97x-d3`（disp 行なし、分散は汎関数に内蔵）、PBE0 は `xc pbe0` + `disp vdw 4`、どちらも `* library def2-tzvpd`、`task dft energy`。
  - CCSD(T): `ccsd; freeze atomic; end` / `task ccsd(t) energy`。maxiter は書かないので NWChem 既定の `maxit = 20`。scf ブロックは restart のときだけ。
- 未実施: C1（composite.yaml と blocker）、C2（CCSD の maxiter）、C26、U8-10。

#### 3.0.2 妥当性（○）

- **3 層構成は妥当なままです。** W7 と v2 で ΔE は完全に一致しました（HCN 46.621、HONO 13.671、NH3 4.26）。S1 で tight が消えましたが、W7 の本流のデッキにも tight は入っていませんでした。
- **S15 のスケール 1.0**: PBE0/def2-SVPD の文献の ZPVE 因子は 0.9848（Kesharwani, Brauer, Martin 2015）なので、1.0 は ZPE を約 1.5% 過大にします。ただし ΔG‡ で ≤0.054、ΔG_assoc で 0.089 しか変わらず、手法誤差より 1 桁以上小さい値です。借り物の PBE0/MG3S 値をやめて 1 本にしたのは、より正直な選択です。GoodVibes は `freq_scale_factor` と `zpe_scale_factor` を明示すれば DB を引きません（installed `api.py:117-123`）。
- **M5** は物理的に意味のある基準です（|ΔE_rxn| ≤ 1 は符号を問わず、±1 をまたいで割れたときだけ blocker）。**S8** でパネルの参照はすべて diffuse 付きになりました。HCN の CCSD(T)/TZVPD は 47.81 / 15.22 で、文献の CCSD(T)/ANO 48.3 / 14.7 に近い値です。
- **SVPD の障壁と CCSD(T)/def2-TZVPD（SVPD 構造での SP）の差**（前回のプローブの出力 `/home/user/hfauto_review_probe/u0/hono_tzvpd` の再集計と、本プローブ）:

| 系 | SVPD | CCSD(T)/TZVPD | 差 |
|---|---:|---:|---:|
| HCN | 46.62 | 47.81 | −1.19 |
| HONO（TS −205.423787、trans −205.442352 Eh） | 13.67 | 11.65 | +2.02 |

  design.md:149 の「±1.7 以内」「約 1.5 未満は手法誤差」は、S8 で退けた TZVP 参照に基づくので少し甘い値です（U0-N3）。NH3 は TZVPD で再評価していません。
- **残る化学的な限界**（文書化済み、設計の範囲内）: ΔG_assoc は非 CP の SVPD の値（NH3·HF で約 3.7 過大に結合）。amine·(HF)n の結論は PBE0 の PES 上のもの。composite で SP 層の障壁が消える（ΔE‡ ≤ 0）場合の扱いは、コードにも文書にもありません（S23 の文書で扱う）。
- **S9 の fail-closed**: G は正しく None になりますが、ReactionThermo の ΔE は freq 層に黙って戻ります（U0-N1）。blocker で順位からは外れるので、実害は表示だけです。
- ライブラリの仕様:
  - NWChem CCSD: 閉殻 RHF 参照だけ、MAXITER の既定 20、THRESH 1e-6。S17 の閉殻の拒否は公式の制約と整合します。
  - NWChem DFT: `convergence energy` の既定は 1e-6 なので、1e-7 を全デッキで明示するのは正しい。`xc wb97x-d3` は分散を内蔵するので disp 行は不要（method ファイルと一致）。

#### 3.0.3 有用性（○）

- **HONO trans→cis の順位が残るようになりました。** 既定パネルを v2 の run の複製に追記すると、ΔE‡ / ΔE_rxn は SVPD 13.67 / +0.13、PBE0/TZVPD 13.24 / +0.24、ωB97X-D3/TZVPD 12.82 / +0.33 で、sign_disagreement=False、panel_report でも rank 1（ΔG‡ 12.224）でした（`/home/user/hfauto_review2_probe/u0/hono_panel`）。前回は外れていました。ただし、実データで解消させたのは S8（負になるレベルがなくなった）で、M5 の不感帯の分岐は単体テスト（`test_summary.py:92-99`）だけで確かめています。
- **パネルは小分子なら数分で終わります。** HCN の既定パネルは 6 SP・36 s（無負荷に近い条件）。CCSD(T) の opt-in は methods に 1 語足すだけで、DFT の 6 点は JobStore から再利用されました（hits 6 / misses 3）。
- **composite**: run を複製して追記すると、SP と GoodVibes がすべて JobStore のヒットになり、2.2 s で終わりました（HCN、ωB97X-D3/TZVPD 層、ΔE‡ 46.37、ΔG‡ 42.18 → 41.93、rank 1）。一方、design.md:150 の手順（別の run-dir で sp → thermo → report）を文字どおり実行すると、`ValueError: stage 'comp_sp' (sp) has no input of type ['minimum']`（`runner.py:151-153`）で即座に止まりました（`/home/user/hfauto_review2_probe/u0_assess/fresh_composite`）。
- composite の主な価値は ΔG_assoc にあります。障壁は HCN で ωB97X 層のほうが CCSD(T) から遠い（46.37 と 47.81）ので、既定に組み込まない判断は正しいです。
- 手法パネルと composite の実計算の確認は本プローブで補えましたが、validation.md の記録は「未実行」のままです（v2 の検証 run についての記述としては正確）。

#### 3.0.4 複雑さ（○、前回 △）

- 減ったもの: method ファイル 6 → 5 本。ThermoSettings の temperatures_K・zpe_scale・invert_soft_cm1 とスケール因子の DB・別名・解決関数（約 35 行）。energy_method の「実装はあるが黙って誤る」状態。
- 増えたもの: `_check_stages` の 6 行、`_SIGN_DEADBAND_KCAL` 1 つ、`_Subject.layer_missing` 1 つ。knob は増えていません。
- 残る複雑さ（小さい）: composite は利用者が組み立てる経路で、しかも文書の手順に穴があります。conditions を 3 本の pipeline に重複して書いています。WSL の本番経路では code_version が常に `-dirty` になります（U0-N7）。

#### 3.0.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U0-I3（前回） | — | composite の配線は一部だけ（上記）。C1 は今回見送り |
| U0-I8（前回） | low | CCSD の maxiter は 20 のまま。HONO の TS/TZVP は 18/20 反復 |
| U0-N1 | low | energy_method の fail-closed のとき、ReactionThermo の dE_act / dE_rxn は freq 層に戻る（`thermochemistry.py:78` と `_kcal` 190-193）。全部欠けると SVPD の 46.62 が composite 扱いで出る。TS の SP だけ欠けると 2 つの LOT を混ぜた値になり、HCN の実エネルギーで 164.4 kcal/mol（真の値 46.37）。blocker（thermo_unavailable、mixed_level_of_theory）で順位からは外れるが、report.html に表示される（U7-N1 と同じ） |
| U0-N2 | low | design.md:150 の composite の手順は、書かれたとおりには実行できない（新しい run-dir では view が空）。実際に使える「run の複製への追記」は文書にない。sp の既定 targets（reaction_stationary_points）のままだと、反応に入らない単量体が energy_layer_missing になり、composite の主な利点である ΔG_assoc が None になる（U8-N1 と同じ） |
| U0-N3 | low | design.md:149 の手法誤差の数値が TZVP 参照に基づく。TZVPD 参照では HCN −1.19、HONO +2.02 で、順位差の目安は約 2 が正確。NH3 は TZVPD で再評価していないので数値を書くべきでない |
| U0-N4 | low | CCSD(T) の opt-in の目安「重原子 4 個程度」に対し、maxit 20 の余裕が小さい（HONO TS/TZVP 18/20）。未収束は rc≠0 の fail-closed で誤った値にはならないが、そのレベルの行が欠け、再実行しても同じ失敗を繰り返す（U0-I8 と同根） |
| U0-N5 | low | 「追記後は panel_report が最終」（README:37、design.md:147）は、method_panel.yaml の conditions が元の run と同じときだけ成り立つ（`report.py:47-48` は pipeline 自身の conditions を使う）。違えば全反応が thermo_unavailable になる。panel_report に T_K を書けば回避できる |
| U0-N6 | low | thermo がないとき、ranking.csv の blockers に dzpe_out_of_tolerance が付き、本当の理由の energy_layer_missing は出ない（プローブ hcn_missing、v2 の水と TMA·(HF)₂ も同じ形）。dzpe が None なら dG_act も None なので、この blocker は情報を増やさない（U8-8 と同根、S22 で扱う） |
| U0-N7 | low | WSL の `/mnt/c` の checkout では、CRLF の差で `git status --porcelain --untracked-files=no` が 83 件を返し、`code_version` が常に `-dirty` になる。`-c core.autocrlf=input` を付けると 0 件（Windows 側はどちらでも 0 件） |

#### 3.0.6 追加の改良案

- **U0-R1 U0 の文書を直す（composite の手順、手法誤差の数値、パネルの conditions）** — should、修正して採用、中立【S23 に統合】
  - (a) design.md:150 を次の内容にする: 「元の run-dir を複製し（`cp -a`）、複製に `--run-dir` で sp（`targets: all_minima`）→ thermo（`energy_method`）→ report を追記する。stage id は元と別にする（例: `sp_tz` / `thermo_tz` / `report_tz`）。JobStore も複製されるので、計算済みの SP と GoodVibes は再利用される。sp の既定 targets では単量体が energy_layer_missing になり、ΔG_assoc は None になる。SP 層で ΔE‡ ≤ 0 になった反応の ΔG‡ は順位として読まない」。FileRef は run からの相対パス（`runner.py:72-75`）なので、複製は安全です。
  - (b) design.md:149 の数値を、測った系だけにする: 「CCSD(T)/def2-TZVPD（SVPD 構造での SP）と比べて HCN −1.2、HONO +2.0。順位差の目安は約 2 kcal/mol」。NH3 の数値は書かない。
  - (c) method_panel.yaml のコメントと README:37 に「conditions は元の run とそろえる」を 1 句足す。
  - 効果: composite が数秒〜数分で使える経路になる（全段の再計算 約 31 分を避けられる）。SP 層で障壁が消えた反応を誤って読むことを、blocker なしで防げる（C1 の blocker の代わり）。コスト: 文書 4〜5 行。
- **U0-R2 CCSD ブロックに maxiter 50 を書く（C2）** — should【§5 の C2】
  - `backends/nwchem/input.py:166` の ccsd ブロックを `[module, "  freeze atomic", *(["  maxiter 50"] if module == "ccsd" else []), "end"]` にする（mp2 には書かない）。NWChem の `ccsd_input.F` は maxiter（rtdb の `ccsd:maxiter`）を受け付けます。JobStore のキー（method・molecule）は変わらないので、既存の成功結果はそのまま使えます。未収束で失敗した既存のジョブは `--retry-failed` で再実行します。
- **U0-R3 validation.md に、パネルと composite のプローブ記録を足す** — could、修正して採用【C39 に統合】
  - 既存の「M5・S8 は method_panel を実行していない」は v2 の検証 run についての記述として正確なので、書き換えない。「レビューのプローブ（v2 の run の複製に追記、`/home/user/hfauto_review2_probe/u0`）」という短い小節を足し、負荷時の所要時間は参考値と断る。
- **U0-R4 layer_missing の subject を含む ΔE を None にする** — could【C28】
  - `_kcal`（`thermochemistry.py:190-193`）で `sa.layer_missing or sb.layer_missing` なら None を返す（1 行）。test_thermo_sp_stages.py の S9 の節に `dE_act_kcal is None` の assert を 1 つ足す。2 つの LOT を混ぜた ΔE（HCN で 164 kcal/mol）が report.html に出なくなります。
- **U0-R5 dzpe の blocker を、値があるときだけにする** — could【S22 に統合】
  - `chemistry/gates.py:313` を `if dzpe is not None and abs(dzpe) > tolerance:` にする。dzpe が None なら dG_act も None で thermo_unavailable がすでに付くので、順位に入る反応は増えません。U8-R1（修正版）と同じ変更でそろえます。
- **U0-R6 code_version の git に `core.autocrlf=input` を渡す** — could【C38】
  - `pipeline/config.py:162` を `_git(Path(tmp), git, "-c", "core.autocrlf=input", "status", "--porcelain", "--untracked-files=no")` にする（引数 2 つ）。WSL と Windows のどちらでも、`-dirty` が本当の未コミットの変更だけを示すようになります。`.gitattributes` で改行を正規化する方法は作業ツリーを書き換えるので、この最小の修正のほうがよい。

#### 3.0.7 見送った案

| 案 | 理由 |
|---|---|
| composite.yaml を同梱し、blocker `barrier_vanishes_at_energy_layer` を足す（C1） | composite の主な利点は ΔG_assoc で、障壁は HCN で ωB97X 層のほうが CCSD(T) から遠い。同梱すると元の run に追記した場合に thermo の view が置き換わる。U0-R1 の文書（複製して追記し、ΔE‡ ≤ 0 は読まない）で目的は足りる |
| thermo の artifact id に stage id を含め、1 つの run に複数のエネルギー層を並べる | report がどの thermo を使うかの規則が新たに要る。複製した run への追記で同じことが数秒でできる |
| 符号ゲートの不感帯を Policy の knob にする、パネルの幅を同順位の判定に組み込む | 定数 1 kcal/mol で HONO 型を正しく扱えた。knob を増やしても判断材料は改善しない（前回の判断を維持） |
| CCSD(T) を既定のパネルに戻す、系の大きさで自動的に入れたり外したりする | 17 原子で 1 点数時間。opt-in とコメントの目安で足りる |
| S10 の検査を低レベル側や run をまたいだ追記に広げる | 低レベルの結果はすべて DFT で判定し直すので単一 PES は壊れない。run をまたいだ食い違いは mixed_level_of_theory が止める |
| energy_layer_missing を ranking.csv の blockers に伝える | species_thermo の notes にすでに残っている。S22 で blockers を正確にすれば足りる |
| report の T・標準状態を既存の ReactionThermo から推定する（U8-10 をコードで直す） | 推定の規則が増える。失敗は明示されるので、U0-R1(c) の 1 句で十分 |
| PBE0/def2-SVPD 専用のスケール因子をフィットする、0.989/0.975 に戻す | 1.0 にした影響は ΔG‡ ≤0.04、ΔG_assoc ≤0.07 と実測され、手法誤差より 2 桁小さい |
| 停留点を TZVPD か ωB97X-D3 に上げる、会合量の CP 補正を自動化する | 前回の見送り理由（17 原子の freq 約 900 s/本、ωB97X の freq は有限差分、NWChem の bsse の不具合）が変わっていない |

#### 3.0.8 変更不要な点

- GFN2-xTB → PBE0-D3BJ/def2-SVPD（grid fine、`convergence energy 1e-7` を全デッキで明示、`disp vdw 4`）の単一 PES → GoodVibes 4.3.0 qRRHO の 3 層構成。
- 停留点系を numerics=True の same_pes で扱うこと、mixed_level_of_theory の blocker。S14 の state=False は会合量の判定だけに使うこと。
- 既定パネル `[pbe0-d3bj_def2-tzvpd, wb97x-d3_def2-tzvpd]` + 停留点レベルの参照行。CCSD(T)/def2-TZVPD は閉殻の小分子で opt-in、開殻は INPUT_INVALID `wft_closed_shell_only`。
- ±1 kcal/mol の不感帯をモジュール定数にしていること。
- energy_method の fail-closed と composite の式 G = E_SP + (G_GV − E_GV)。composite を既定の流れに組み込まないこと。
- 停留点 method の一意性を読み込み時に id だけで検査すること。
- vib_scale 1 本（既定 1.0）を GoodVibes に明示して DB を引かせないこと、温度を conditions だけから取ること、thermo のキーに振動数・Hessian の sha・温度を入れること。
- method_panel を既存 run への追記型の任意実行にし、Level.full_key をキーにすること。

---

### 3.1 U1 入力構造・電子状態（structures）

評価: 前回 ○ / ◎ / ◎ → 今回 **○ / ◎ / ◎**（妥当性を ○ に留める理由は 3 つ。SMILES の同位体を黙って ¹H にする経路、電荷のある系・開殻系を structures から NWChem まで通した実 run がないこと、組成の多重度が高スピン固定であること（C8、U2 の範囲）。S24 を入れれば U1 単体としては ◎ 相当）

#### 3.1.1 現状仕様の要点（変更点を中心に）

- **読み込み時の検査**（`core/system.py`、pydantic、extra=forbid）。違反はすべて ValidationError で、CLI は exit 2、run-dir は作られません。
  - 【新】`_one_source`（27-31）: xyz と SMILES のちょうど一方（S12）。
  - 【新】`components: dict[str, PositiveInt]`（37）: 個数 ≥ 1（S12）。
  - 【新】`ReactionInput._check_coordinate`（48-55）: 座標項の原子数が distance 2 / angle 3 / dihedral 4、添字は 0 以上で相異なる（S11a）。
  - 【新】`_check_references`（69-92）: species id と composition id を合わせて一意、reaction id も一意（S12）。未知の参照、端点の role=endpoint。【新】**宣言反応の端点は xyz**（88-91、M7。メッセージ 'must be given as xyz; SMILES cannot fix the atom mapping or the conformer'）。
  - xyz ファイルが存在するかは検査しません。
- **structures stage**: `_write`（`stages/structures.py:48-54`）から実行時の排他検査を削除（S12）。SMILES の処理（`chemistry/smiles.py`、35 行、変更なし）は、'.' の拒否、形式電荷の照合、`AddHs`、ETKDGv3（seed 20260925）で 1 回埋め込み、力場なし。
- `check_electronic_state`（`chemistry/electronic_state.py:17-33`）: m < 1 の拒否（21-22）、【新】未知元素は `ValueError('unknown element(s): …')`（24-26、C5）、電子数の偶奇（27-33）。
- 【新】**宣言反応の検査**（`structures.py:81-93`、S11b）: 両端の composition_id（組成式・電荷・多重度）と元素列の一致、座標の添字 < 原子数。違反なら両端を INPUT_INVALID。片方の端点だけが失敗していれば continue。
- **状態ラベル**（`topology.py:231-239`、変更なし）: 結合 ≤ 1.15、非結合 ≥ 1.45 Σr_cov、WL 3 反復、先頭 8 桁。【役割の拡大】S6 以降は ReaDuct の `matches_source`（`readuct/worker.py:48-51`）で生成物の棄却にも使われます。
- 下流への受け渡し（変更なし）: NWChem は `charge q` と `mult m`（m > 1 なら `odft`）、xTB と CREST は `--chrg q --uhf m−1`。組成の電荷は和、多重度は高スピン（`conformer_search.py:201-206`）。
- 文書（C6）: design.md:138 に読み込み時の検査の一覧、SMILES は 1 配座で立体は任意、同位体は扱わない、端点は xyz、組成の規則、m=1 は閉殻 RKS。
- 未実施: C4（SMILES の同位体とラジカル）、C8。

#### 3.1.2 妥当性（○）

- **M7 の前提はライブラリの挙動で裏付けられます。** RDKit 2026.03.3 の `MolFromSmiles` は既定で明示の [H] を除き、`AddHs` は H を重原子の後ろに付け直します（'[H]OC=O' → O,C,O,H,H）。SMILES では原子写像を利用者が決められません。
- S11・S12・C5 で、電荷や多重度の取り違え、範囲外の添字、id の衝突、未知元素が DFT の前に止まるようになりました。v2 の実 run の NWChem デッキ 81 本はすべて `charge 0` / `mult 1`（odft 0 本）、xTB と CREST の argv は `--chrg 0 --uhf 0` でした。
- **まだ黙って誤る経路**（すべて low）:
  - SMILES の同位体: '[2H]O[2H]' は H2O_q0_m1 / H2O_b0287123 で水と同じラベルになり、'[13CH4]' は CH4 になります。RDKit は `GetIsotope()=2`、`GetMass()=2.014` を持っていますが、`GetSymbol()` が 'H' なので、書き出す xyz では ¹H になります。xyz で 'D' と書けば拒否されるので、入力経路によって扱いが違います。design.md:138 の「同位体は扱わない」は、拒否ではなく黙って変換するという意味になっています。
  - '[O][O]' m=1 と '[CH2]' m=1 は通ります（対象化学の外で、[CH2] m=1 は正当な意図でもある）。
  - 組成の多重度は高スピン固定（C8、U2 の範囲）。
- **S6 で状態ラベルの役割が変わりました。** 前回の `matches_source` は「置換不変 RMSD < 0.1 Å かつ同じ結合グラフ」でした。今は状態ラベルの一致だけを見るので、1-WL の衝突、立体・E/Z を区別しないこと、中間帯を結合とみなすことが、same_as_source の判定に直接効きます。HONO の trans と cis はどちらも HNO2_93fd17a7 です。W7 → v2 で、tma_hf2 の explore の唯一の「生成物」（同じラベルの配座違い、ΔE −2.6 kJ/mol）が same_as_source になり、生成物は 1 → 0 件になりました。これは配座を生成物にしないという意図どおりで、真の反応生成物を失った証拠はありません。
- ライブラリの仕様:
  - RDKit の ETKDGv3 の既定値は spec の一覧と一致し、公式の Getting Started も ETKDG の後の最小化は不要としています（1 回埋め込み・力場なし・AddHs・seed 固定は推奨どおり）。
  - Windows の RDKit 2026.03.5 と WSL の 2026.03.3 で、'N' と 'CN(C)C' の座標は同一でした（違いは CRLF だけで、指紋は同じ）。
  - xTB の `--chrg` と `--uhf` は `.CHRG` / `.UHF` ファイルと xcontrol より優先されます。

#### 3.1.3 有用性（◎）

- 1 秒未満・ジョブ 0 本の必須の入口です（W7 は 8 run・23 種、v2 は 9 run・24 種で、すべて done、failed 0）。SpeciesRecord の指紋、composition_id、状態ラベル、電荷、多重度は W7 と v2 で完全に一致しました。SMILES 由来の単量体の xyz も同じバイト列です（tma の md5 bb3c3bd5…）。
- 新しい検査は v2 の実 run では一度も発火していません（同梱の system は正しい入力）。効果はプローブで確かめました。
  - M7: SMILES の端点で CLI exit 2、0.79 s、run-dir は作られない。
  - S11a: angle の原子 2 個、負の添字、重複した添字の 3 件を拒否。S12: 5 件を拒否。
  - S11b: 範囲外の添字（angle [0,1,7]）を `--to structures` で 1.56 s に止めた（前回は両端の DFT opt と freq の後に IndexError）。
  - C5 など: 'Xx' → 'unknown element(s): Xx'。'Xe'・'Sn'・'cl' → 'no covalent radius tabulated …'。
- **有用性を削いでいる点**（どれも fail-closed だが、利用者に理由が届かないか、計算が無駄になる）:
  1. 端点 2 つだけの system（known_endpoints の典型）で S11b が両端を落とすと、次の stage が `has no input of type ['species']` の traceback で止まります。理由は structures/manifest.json にしか残らず、report.html（1.6 KB）にもありません（U1-N1）。
  2. xyz ファイルがないとき、dry-run は通ります。実行すると config_sha で FileNotFoundError になり、structures は 'running' のまま残ります（U1-N2）。
  3. 片側の端点だけが失敗すると（q/m の取り違えで偶奇の検査に落ちるなど）、残った端点が DFT opt と freq まで進み、原因とずれた理由（endpoints_not_on_one_pes）で BLOCKED になります（U1-N4）。

#### 3.1.4 複雑さ（◎）

| ファイル | 前回 | HEAD |
|---|---|---|
| core/system.py | 83 行 | 111 行 |
| stages/structures.py | 90 行 | 94 行 |
| chemistry/electronic_state.py | 59 行 | 58 行 |
| chemistry/smiles.py | 35 行 | 35 行（変更なし） |

- 追加は宣言的な validator 3 本、PositiveInt、structures の条件 1 本だけです。未知元素の fail-open の分岐と、`_write` の実行時の排他検査は削除されました。設定キーは 0 のまま。
- 検査の置き場所が明確になりました（記述だけで分かるもの＝形・一意性・端点の形式は読み込み時、原子数や電子数が要るもの＝範囲・偶奇・組成の一致は structures）。import-linter の layers 契約（core は xyz を読まない）とも合っています。U1 関連のテスト 36 件は pass。

#### 3.1.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U1-03（前回） | low | SMILES の同位体・ラジカルを検査しない（上記）。S24 で同位体の部分を閉じる |
| U1-06（前回） | low | 状態ラベルの限界。S6 で生成物の棄却にも使われるようになり、「グループ化にしか使わない」という前回の根拠は古くなった |
| U1-N1 | low | structures で全種が落ちると、次の stage が理由のない traceback で止まる（`runner.py:150-152`、`cli/main.py:139` は run_pipeline の例外を捕まえない、`reporting/html.py` は failed artifact を描かない）。S11b の改善（早く止まる）の副作用として、known_endpoints では必ずこの経路に入る |
| U1-N2 | low | xyz ファイルの欠損を dry-run が通し、実行すると `runner.py:114` の config_sha が FileNotFoundError を投げる。`runner.py:168-176` の `layout.update(config_sha=…)` が try の外にあるので run_state が 'running' のまま残り、再実行しても `_fresh` が短絡して同じことを繰り返す。preflight も xyz を見ない |
| U1-N3 | low | 反応単位の入力誤りが、共有端点を通じて他の正しい宣言反応も止める（プローブ: r1 A→B が範囲外の添字、r2 B→C が正しいとき、B が落ちて r2 は continue になり、r2 自身の理由はどこにも残らない）。同梱の system は反応 1 本 |
| U1-N4 | low | 片側の端点だけ失敗したとき、残った端点が DFT まで進み、原因とずれた理由で BLOCKED になる（`structures.py:83-84`、`hypotheses.py:147-160`、`state.py:108-111`）。17 原子なら freq だけで約 15 分の無駄 |
| U1-N5 | low | S6 で状態ラベルが生成物の棄却にも使われるようになったが、ラベルの限界（1-WL、立体と E/Z を区別しない、中間帯は結合）が design.md:141 に書かれていない |
| U1-N6 | low | 仕様書の行番号の誤り（文書の正確さの問題）: electronic_state.py の `check_electronic_state` は 17-33、未知元素は 24-26、topology の `_lookup` は 51-55 |
| U1-N7 | low | （検証で追加）`hypotheses._declared`（158-159）の 'endpoints differ in atom order' の ValueError は、S11b が structures で元素列の不一致を先に落とすので、実質的に到達しない防御コード。削除するか S11b への参照コメントにすれば、検査の二重化が解消する（低優先の整理） |

#### 3.1.6 追加の改良案

- **U1-R1 SMILES の同位体を拒否する（C4 の同位体の部分）** — should、中立【S24】
  - `smiles.py` の形式電荷の照合の直後に 1 行: `if any(a.GetIsotope() for a in mol.GetAtoms()): raise ValueError(f"SMILES {smiles!r} specifies isotopes; hfauto uses natural-abundance masses")`。design.md:138 の「同位体は扱わない」を「同位体を指定した SMILES は INPUT_INVALID（xyz の 'D' も未知元素として拒否）」に直す。test_smiles_input.py に '[2H]O[2H]' の行を 1 つ足す。
  - 効果: 対象元素（H）で黙って誤る最後の経路を閉じる（H 移動や HF・DF の KIE を意図した入力が ¹H の質量・ZPE で計算されて気付かれない、をなくす）。`Atom.GetIsotope` は Chem の中にあるので import は増えない。約 2 行。
- **U1-R2 stage の入力がないときの例外に、上流の失敗理由を添える** — should、中立【S25】
  - `pipeline/runner.py` の execute_stage の `if missing:` の分岐で、`reasons = sorted({f"{a.artifact_id}: {a.failure.reason}" for t in spec.consumes for a in inputs.of(t, ok_only=False) if a.failure})[:3]` を作り、例外文の末尾に `; upstream failures: …` を付ける。`Manifest.of(type, ok_only=False)` は `manifest.py:52` にあり、layout.view には上流の failed artifact も含まれる。テストを 1 件足す。
  - 効果: known_endpoints で S11b・偶奇・未知元素のどれかが発火したとき、traceback の最終行に理由（例: 'species_hcn: reaction iso: coordinate atom index out of range'）が出る。conformers や minima など、ほかの stage の入力切れにも同じ形で効く。約 4 行。
- **U1-R3 読み込み時に、xyz の存在と宣言反応の両端の電荷・多重度の一致を検査する** — could、中立【C31】
  - (a) `load_system` でパスを解決した後、`missing = [s.id for s in species if s.xyz is not None and not s.xyz.is_file()]` があれば ValueError（CLI は exit 2 で、dry-run でも止まる）。主な価値はこちらで、'running' の残留と dry-run の素通りを一度に塞ぐ。
  - (b) `_check_references` の反応ループで、両端の charge と multiplicity の一致を検査する（宣言的、ファイルは読まない）。両端が成功すれば S11b が捕まえるので、効くのは片側だけ偶奇の検査に落ちる typo だけだが、U1-N4 の主な引き金を塞ぐ。
  - どちらも core の中で完結し、layers 契約に反しない。約 5 行とテスト 2 件。
- **U1-R4 状態ラベルの限界と same_as_source での使われ方を design に書く** — could【C39 に統合】
  - design.md:141 に 1 文: 「状態ラベルは結合グラフ（1-WL、立体と E/Z を区別しない、中間帯は結合）なので、配座・立体・E/Z だけが違う IRC の端は same_as_source として捨てる」。explore で立体・E/Z 異性化が見つからない理由と、今の挙動が意図であることが記録に残る。

#### 3.1.7 見送った案

| 案 | 理由 |
|---|---|
| 片側の端点が失敗したら相方も INPUT_INVALID にする | 相方を使う他の正しい宣言反応まで止める（U1-N3 の cascade を広げる）。主な引き金は U1-R3(b) で読み込み時に止められる |
| S11b の検査（元素列の一致、添字の上限）を全部読み込み時に移す | xyz を読む必要があり、layers 契約で core は chemistry.read_xyz を使えない。得られるのは理由の可視性だけで、これは U1-R2 で足りる |
| structures で反応単位の failed artifact（reaction 型）を作る | 新しい artifact 型と下流の対応が要る。多反応で端点を共有する system は同梱になく、今も理由付きで fail-closed |
| SMILES のラジカル電子数と多重度を照合する（C4 の後半） | '[CH2]' m=1（一重項カルベン、RKS として妥当）を誤って拒否する。'O=O' m=1 は RDKit でもラジカル 0 で検出できない。奇電子を m=1 で与えた場合は偶奇の検査が止める |
| reactant ≠ product を読み込み時に検査する | 下流で mapped_rmsd=0 → 行 2 の SAME_BASIN になり、端点 1 つの DFT だけで閉じる。無駄が小さく、誤った結果も出ない |
| multiplicity を PositiveInt にする、coefficient ≠ 0 を検査する | m < 1 は structures で INPUT_INVALID になり、U1-R2 を入れれば理由も表示される。coefficient 0 は非現実的 |
| xyz の元素記号の大文字小文字を正規化する | 対象元素（H、C、N、O、F）では起きず、今も fail-closed |
| 予約接尾辞（<id>_c00 など）の species id を拒否する | そのような命名は非現実的で、規則を 1 つ増やす価値がない |
| runner で config_sha の計算を try の中に移す | 'running' のまま残る実際の引き金は xyz の欠損だけで、U1-R3(a) で読み込み時に止まる |
| 状態ラベルに立体（CIP、E/Z）を加える、WL の反復を増やす | 錯体や中間帯で結合次数の推定が不安定になる。対象化学（PT、H 結合錯体）では効果がない |
| RDKit の版数を config_sha に入れる | 2 つの版で座標が同一だった。再開時は保存済みの xyz を使うので正しさに影響しない |

#### 3.1.8 変更不要な点

- 検査の分担（記述だけで分かるものは読み込み時の pydantic validator、原子数・電子数が要るものは structures）。
- M7（端点は xyz 必須）と、そのメッセージ（原子写像と配座という理由を示す）。S11b の composition_id と元素列の両方の比較、reason 2 種類。
- S12（排他、PositiveInt、id の一意性）と、`_write` の if/else + cast。C5 の fail-closed。
- xyz をバイト単位で複製して 1e-6 Å の指紋を取ること。ETKDGv3 で 1 回、seed 固定、力場なし。
- '.' を含む SMILES の拒否、形式電荷の照合、偶奇の検査と m ≥ 1。failed artifact として記録して stage を続ける方針。
- composition_key と指紋が q と m を含むこと。組成の電荷は和。Cordero の 37 元素。StructuresConfig に設定キーを持たないこと。
- S6 で matches_source を状態ラベルにしたこと自体（配座を生成物として登録しない、という目的に合う）。

---

### 3.2 U2 配座探索・錯体配置（conformers）

評価: 前回 △ / ○ / ○ → 今回 **○ / ○ / ◎**（複雑さの ◎ は、組成では到達せず単量体では実 CREST で失敗する `--noreftopo` の分岐が残るという但し書き付き。C34 で名実ともに ◎。妥当性を ◎ にしない理由は、単量体の停止後の再実行の失敗、開殻の組成の `--noopt` 失敗、C8 の未実施。有用性は、下流で新しい極小を生んだ実例がまだないので ○）

#### 3.2.1 現状仕様の要点（変更点を中心に）

- 【新】**M4**: 組成（nci=True）の CREST に必ず `--noopt` を付けます（`backends/crest.py:56-58`）。v2 の TMA·(HF)₂ の argv: `crest input.xyz --gfn2 --nci --quick -T 4 --ewin 6 --chrg 0 --uhf 0 --notopo 2,14,15,16,17 --noopt`。key_payload に `"noopt": settings.nci` を追加（`crest.py:180`。M4 以前の組成結果を黙って再利用させないためのキャッシュの無効化）。
- 【新】**C7**: 利用者が書ける項目は quick、ewin_kcal、seeds_per_composition、keep_per_state の 4 つ（`conformer_search.py:48-54`）。旧キー（settings、window_kcal、threads）は extra=forbid で拒否。ConformerSettings から threads と topology='off' を削除（`protocols.py:44-50`）。スレッド数は site の `engines.crest.execution.threads` だけで決まり（`-T` と `OMP_NUM_THREADS=n,1`、`crest.py:123-125`）、同時実行はコアセマフォ（ranks×threads）に任せる。site で省くと `-T 1`。
- 【新】**S13**: `window_kcal` を削除。状態ラベルごとにエネルギーの低い順で keep_per_state 個を残し、BUG-08 の規則（ラベル内に数値エネルギーがあれば None を捨てる）は残す（`:118-130`）。
- 【新】**C5**: 組成の部品に未知元素があれば INPUT_INVALID。
- 変わらない部分: 単量体は `--gfn2 --quick -T 4 --ewin 6 --chrg --uhf`。小さい剛体の単量体（重原子 ≤ 3、回転可能結合 0）は省略。トポロジー停止では crestopt.log の最終フレームを crest_topology 候補として残し、`--noreftopo` で 1 回だけ再実行（組成では `--noopt` のため到達しない）。placement（`chemistry/placement.py` は差分なし）、CREST に渡すのは seed00、重複は 0.1 Å / 0.1 kcal/mol。
- 文書（C9、C6）: design.md の conformers 行に `--nci --quick` の実効設定、seed00、衝突条件、`--noopt` の理由、窓がないこと、スレッドの出所。

#### 3.2.2 妥当性（○）

- **M4 の修正はソースで裏付けられます**（CREST 3.0.2）。
  - `confparse.f90` の `-noopt` は `env%preopt=.false.` だけを設定し、`crest_main.f90` は `if (env%preopt) call trialOPT(env)` なので、初期最適化とトポロジー検査（`setuptest.f90` の trialOPT）はそもそも呼ばれません。
  - `-nopreopt` は QCG 専用の別フラグで、文書は同義のように並べていますが 3.0.2 では別物です。`--noopt` を選んだのは正解です。
  - `-notopo <list>` は excludeTOPO を確保すると checktopo を `.true.` に戻すので、部分除外として働きます。
  - `--noopt` の下では CREGEN の参照トポロジーは seed の quicktopo による結合になります。v2 の全組成で topology_removed は 0 で、正しい構造を誤って除いた形跡はありません。
- **実 run**: NH3·(HF)₃ は rc 0・4.84 s・10 配座、TMA·(HF)₃ は rc 0・19.17 s・10 配座で、最安のイオン対（C3H10N+F+FH+FH、−29.65447255 Eh）は前回のプローブの `--noopt` の最安と一致しました。F⁻·(HF)₂ の smoke は 1.37 s。閉殻の陽イオン（NH4⁺·HF、NH4⁺·(HF)₂、H3O⁺·H2O）もプローブで `--noopt` のまま rc 0 でした（`/home/user/hfauto_review2_probe/u2_ion`）。tma_hf2 の stdout に 'Initial Geometry Optimization' はありません（W7 には 1 回）。
- **S13 も妥当です。** GFN2 の PT エネルギー誤差は MUE 約 3 kcal/mol なので、4 kcal/mol で先に切る根拠はありませんでした。NH3·(HF)₃ では CREST の 10 配座のうち 6 件が `--ewin` の縁（5.64〜5.72）にあり、ラベル FH+FH+FH4N の 8 件中 6 件が残り、2 件が keep_per_state で落ちました。旧 4 kcal/mol 窓なら全体で 3 件でした。下流の件数は explore の sources_per_state=2 と dft の per_state 3 / 6 kcal/mol で抑えられます。
- **新たに確認した弱点**:
  1. **単量体の停止 → `--noreftopo` 再実行は、PT 型の停止では実 CREST 3.0.2 で失敗します**（U2-N1）。hfauto の `_search` をそのまま通したプローブ（`u2_gly`、`u2_bala`）で、グリシンと β-アラニンの双性イオンは初期最適化で N→O の PT が起きて停止し（検出と crestopt.log の読み取りは正常）、続く `--noreftopo` の再実行は trial MTD の 6 回不収束で rc 1 になりました（グリシンは別ディレクトリの手動実行でも再現）。停止構造を新しい入力にして既定の単量体 argv で走らせると、グリシン 5 配座（66 s）、β-アラニン 7 配座（60 s）で 2/2 成功しました。`setuptest.f90` によれば `--noreftopo` でも env%ref の座標は最適化後の構造に置き換わるので、失敗の原因は座標ではなく、元の入力から作られた参照トポロジー（gfnff_topo / WBO 由来の SHAKE など）が残ることと推測されます（推測）。さらに、失敗時の `_lowest`（`:167-169`）は未緩和の入力（双性イオン）を組成の部品に選びます。
  2. **開殻の組成は `--noopt` で失敗します**（U2-N2）。NH3·HF⁺（m=2）は、W7 の preopt + `--noreftopo` では rc 0・1 配座、v2 の `--noopt` では 0.56 s で trial MTD が 6 回失敗して rc 1。プローブで 3 通りの argv を比べて再現しました（noopt は失敗、preopt は停止、noref は成功）。開殻の例は 1 系だけで「開殻に固有」は傍証です。組成は placement の seed に fallback するので全損ではありません。
  3. 組成の多重度を黙って高スピンにする扱い（C8）は未実施のままです（design.md:137 に記述はある）。
- 参考: CREST の GitHub issue #285（3.0 で trial MTD が収束しない例）と #419（`-noreftopo` が大きな系で 1 構造しか返さない例）は、3.0 系の trial MTD と noreftopo が脆いという観察と矛盾しません（直接の根拠ではない）。

#### 3.2.3 有用性（○）

- 組成に対する CREST `--nci`（`--noopt` 付き）、placement の seed00、`--notopo` の自動付与、重複除去、版数の pin、BUG-08 の None の扱いは、安定して働くようになりました。
- コスト: stage は 13 s（amine_pilot2）、19 s（tma_hf2）、32〜34 s（(HF)₃）。CREST 1 回は 1.2〜19 s。前回のプローブ（`--noopt` なし、87〜236 s）との差の大半は当時の負荷によるもので、energy+grad の回数は同程度でした（22,302 / 17,808）。DFT の freq（17 原子で約 15 分）と比べれば無視できます。
- 下流への寄与は限定的です。TMA·(HF)₂ では、CREST の 3 配座と入力の neutral・shared_proton が xTB で 1 つの basin にまとまり、DFT の極小は W7 と同じ 1 つでした。(HF)₃ のイオン対が DFT で独立な極小として残るかは計算していません（`--to conformers`）。
- CREST は同じ argv でも実行ごとに揺れます（TMA 単量体で 2 配座と 1 配座、energy+grad 10,636 と 10,996 回）。screen が basin にまとめて吸収するので実害はありません。validation.md:55 の「W7 と同じ 2 配座に 1 配座が加わった」は `--noopt` の効果とは切り分けられませんが、因果は主張していないので修正は任意です。

#### 3.2.4 複雑さ（◎）

- 減ったもの: 利用者が書ける項目 9 → 4。エネルギー窓は CREST の `--ewin` と下流（screen・dft）の 2 重になった。スレッド数は site の 1 か所。ConformerSettings は 5 フィールド。
- 残るもの（どれも小さい）:
  1. 単量体の停止 → noref の分岐: 組成では到達せず、単量体では PT 型の停止で失敗する。C34 に置き換えれば、ConformerSettings.topology、`--noreftopo` の argv、TOPOLOGY_CONTINUED が消え、コードはむしろ減る。
  2. key_payload の "noopt" は settings.nci と情報が重複するが、X12 に沿ったキャッシュの無効化として必要（settings.nci は M4 の前後で同じ値のまま argv だけが変わったので、このフィールドがないと旧キーと衝突する）。
  3. 重複判定の基準が 2 つ（conformers の 0.1 Å / 0.1 kcal/mol と same_minimum の 0.05 Å / 5e-5 Eh）。conformers は前探索で screen が再最適化するので、統一する利点はない。
  4. 読むコードのない diagnostics.json（人がデバッグするのに役立つ）。
  5. 設定値の範囲検査がない（U2-N4）。

#### 3.2.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U2-ISS-3（前回） | — | 組成の多重度は文書化（design.md:137）だけで、動作は同じ（C8） |
| U2-N1 | low | 単量体のトポロジー停止 → `--noreftopo` 再実行が、PT 型の停止では実 CREST 3.0.2 で失敗する（2 系、再現あり）。対象の単量体（TMA、NH3、HF、H2O、ピリジン、アニリン）では停止しない |
| U2-N2 | low | 開殻の組成（NH3·HF⁺、m=2）は `--noopt` で trial MTD が発散して失敗する（M4 による挙動の変化）。`test_crest_receives_charge_and_spin` は陽イオンについて入力の echo しか見ないので、テストでは見えない。noref の成功例も topology-based removals 144 があり、得られた 1 配座が化学的に元の状態かは未確認 |
| U2-N3 | low（一部訂正） | environment.md:64 の F⁻·(HF)₂ の CREST「約 25 s」は、v2 の実測 1.37 s と合わない（25 s は `--noopt` なしで停止したときの値）。validation.md:55 の書き方の修正は任意。design.md の `--noreftopo` の文は「単量体:」の節の中にあり、「組成では到達しないことを書いていない」という指摘は成り立たない。smoke のコメントもすでに 'may fail' と書いている |
| U2-N4 | low | ConformersConfig に値の範囲の検査がない。keep_per_state=0 で組成が黙って 0 件になり、failed artifact も出ない。seeds_per_composition=0 だと原因の違う GATE_REJECTED(no_collision_free_seed) になる |
| U2-N5 | low | site で crest の threads を省くと、黙って `-T 1` になる（旧 stage の既定は 4）。全エンジン共通の規則で、同梱の wsl_local は 4 を明示している |

#### 3.2.6 追加の改良案

- **U2-R1 単量体の停止後の再実行を「停止構造から既定の設定で 1 回」に替え、`--noreftopo` を廃止する** — could、【複雑さ減】【C34】
  - `conformer_search._search`: 停止したら `result.topology_stops[0]` を読み込み（`rt.load_xyz` を partial で渡す）、Molecule(停止構造, 電荷, 多重度) と item.settings のまま `engine.search` を 1 回だけ呼ぶ。2 回目も停止したら Failure とし、ループしない。停止の判定は 'Change in topology detected' の有無だけに単純化できる。
  - ConformerSettings.topology（Literal on/noref）、`crest.py` の `--noreftopo` の argv（:59）、TOPOLOGY_CONTINUED（:39、:76）を削除する。
  - 二重失敗のときの `_lowest` は、未緩和の入力より停止構造（GFN2 の緩和構造）を優先する（キーの並びを 1 つ入れ替えるだけ）。
  - fake の integration テストは 2 回目の呼び出しの molecule が停止構造であることを確かめる形にし、design.md:139 の 1 文を直す。実エンジンの確認はプローブ 1 回（約 2 分）で足り、常設の smoke は要らない。
  - 効果: 実 CREST で 2/2 失敗していた経路が、プローブでは 2/2 成功した。CREST 自身の停止メッセージが示す選択肢そのもので、追加の依存はない。組成の部品に未緩和の双性イオンが使われることもなくなる。約 10 行。
- **U2-R2 文書を実測に合わせる** — could、修正して採用【C39 に統合】
  - environment.md:64 の「約 25 s」を実測（約 1.4 s。`--noopt` なしでは初期トポロジー検査で停止した）に直す。design.md の conformers 行に「開殻の組成は `--noopt` で trial MTD が失敗しうり、そのときは placement の seed が出る」を 1 文足す。environment.md に crest の threads の既定が 1 であることを短く添える。validation.md:55、design.md の「単量体だけ」の明記、smoke のコメントは変えなくてよい。
- **U2-R3 ConformersConfig に範囲の制約を付ける** — could【C30 に統合】
  - keep_per_state と seeds_per_composition に `Field(ge=1)`、ewin_kcal に `Field(gt=0)`（3 行）。拒否の assert を 1 つ足す。
- **U2-R4 組成の多重度が一意でなければ INPUT_INVALID にする（C8）** — could【C8】
  - `conformer_search.py:202` で `coupled_multiplicities` が 2 つ以上の値を返したら `reject(INPUT_INVALID, 'ambiguous_multiplicity:<候補>')`（2〜3 行）。明示の指定欄は足さない。design.md の「多重度は高スピン」を「一意でなければ拒否」に直す。副作用として、ラジカル対の組成は実行できなくなる（対象範囲が閉殻中心なので失うものはほぼない）。

#### 3.2.7 見送った案

| 案 | 理由 |
|---|---|
| key_payload から冗長な "noopt" を外す | X12 に沿ったキャッシュの無効化。外せばキーがもう一度変わるだけ |
| conformers の重複判定を same_minimum（0.05 Å / 5e-5 Eh）にそろえる | screen が全件を再最適化して basin にまとめるので、結果は変わらず stage 間の結合が増えるだけ |
| 開殻の組成だけ preopt と `--noreftopo` に戻す | 主対象外（閉殻中心、X9）のために分岐が増える。失敗しても seed が下流に流れる。文書（U2-R2）と C8 で足りる |
| `--noopt` が失敗したら preopt と `--noreftopo` で再試行する | 前回のプローブで noref は NH3·(HF)₃ で失敗し、TMA·(HF)₃ で最安を取り逃した。成功例は開殻の 1 件だけ |
| 単量体の停止時に再実行をやめ、停止構造と入力だけを出す | 柔らかい単量体の配座探索を失う。U2-R1 が同じくらいのコード量で 2/2 成功を示している |
| 双性イオンの常設の実エンジン smoke（約 2 分）を加える | 対象の単量体では起きない経路で、fake で分岐は押さえている。U2-R1 を入れるときの 1 回のプローブで十分 |
| `-T 1` や固定の乱数で CREST を決定的にする | 4 倍前後遅くなる。揺れは screen が吸収する |
| 6 つの seed それぞれから CREST を走らせる | コストが 6 倍。`--noopt` 付きの 1 本でプローブの最安に到達している |
| 組成では到達しない停止 → 再実行の分岐を nci で塞ぐ | 害はなく、分岐を増やすだけ。U2-R1 で単量体と組成に共通の単純な形になる |
| diagnostics.json を削除する | 人が追うのに役立ち、書く処理は 2 行 |
| site で crest の threads を省いたときの既定を 4 にする | エンジンごとに既定を変えると規則が崩れる。同梱の site は 4 を明示している |

#### 3.2.8 変更不要な点

- 組成に `--noopt` を必ず付けること（M4）と argv の順。seed00 の labile H と受容原子の `--notopo` の自動指定。
- エネルギー窓を置かず、状態ラベルごとに keep_per_state 個を残すこと（S13）と BUG-08 の規則。
- 4 つの knob と extra=forbid（C7）。CREST のスレッド数を site だけで決め、同時実行をコアセマフォに任せること。
- key_payload の "noopt"。`--scratch` なしで attempt ディレクトリで実行すること、版数の pin、atom_order_changed の検査。
- CREST が失敗したら、単量体は入力、組成は placement の全 seed を出して下流を止めないこと。placement.py。小さい剛体の単量体の省略。`--quick` と `--ewin 6`。組成の電子状態の検査（C5）。

---

### 3.3 U3 極小構造の確定（minima screen / dft、MinimumDriver、Registry）

評価: 前回 ○ / ○ / △ → 今回 **◎ / ◎ / ○**（妥当性 ◎ の留保: DFT の mode-follow は、実エンジンでは本プローブの「TS に留まる opt → freq の saddle 判定 → +側の opt（26 点）」までしか通っていない。−側の opt、classify_mode_follow、ts_candidate、両側の再 settle、discovery、DFT の soft 押し、M6 の |ν| 化は、fake と部分的なプローブでしか裏付けがない）

#### 3.3.1 現状仕様の要点（変更点を中心に）

- 【新】**S1 その場での登録**（`stages/minima.py:92-98, 217-231`）: species_id 順の直列ループで settle（relax → その場で register → diagnostics）。xTB の thread_map 並列と「全部 relax してから一括登録」は削除。mode-follow の両側は親の直後に settle。
- 【新】**同一性の 1 基準**（`drivers/minimum.py:241-257`、`chemistry/identity.py:120-162`）: `Registry.add` は known なら join、それ以外は `find(composition_id, level_key)`（identity.assign: ≤ 0.05 Å かつ ≤ 5e-5 Eh、次点と分離）で join、なければ新しい basin。compare_minima、Verdict、ambiguous、rejudge、`tight` を削除。mode-follow のラベルは同じ閾値の bool `same_minimum`。
- 【新】tight の削除で NWChem opt のキーが変わった（HCN d4fb7dbddc4a → e16a78b90518、TMA d8762f3fe441 → a80a4e793d96。デッキは scratch_dir の行以外同一）。
- 波及: M3（dft の `_pool` と explore が screen の最適化構造から始まる）、M6（thermo は全モード |ν|、noise / soft の note は注記だけ）、S7（optimize は 0.5 Å 以内の init_hessian を受け付ける。minima は開始構造そのもので計算するので実効差なし、mode-follow と soft 押しには Hessian を渡していない）、S10。
- 変わらない部分: 虚振動の 3 段方針（noise 10 / saddle 50 cm⁻¹）、soft 押し +0.1 Å を 1 回、mode-follow ±0.1 Å を最大 2 サイクル、選択（window 6 kcal/mol、per_state 3、rerank_top 8、rerank_sp）、質量加重の Eckart 射影、NWChem opt（`trust 0.1`、maxiter 100、2 断片以上で `inhess 2`）と別ジョブの解析 freq、xTB の `--opt vtight` と `--hess`、LADDER、`_known` の version の自己比較。
- 未実施: C10。

#### 3.3.2 妥当性（◎）

- 手順（opt → 別ジョブの解析 freq → is_minimum → Eckart 射影 → 3 段の虚振動方針）は変わらず、停留点を確認する正しい方法です。
- 前回化学的な誤りとした 2 点が解消しました: explore の出発構造（M3、`explore.py:150`）と、thermo が noise / soft の負モードを捨てる件（M6）。インストール済みの GoodVibes 4.3.0 の `calc_freerot_entropy` の μ' = μ·Bav/(μ+Bav) は ν→0 でも有界で、`calc_damp` は (ν/100)⁴ で減衰するので、|ν| にした小さな負モードのエントロピーは自由回転子の極限（約 5.5 R、298 K で −TS 最大約 3.3 kcal/mol）にとどまり、発散しません。
- 同一性の基準は ≤ 0.05 Å かつ ≤ 5e-5 Eh（0.031 kcal/mol）の 1 つになりました。CREGEN の既定（rthr 0.125 Å、ethr 0.05 kcal/mol）より厳しく、add の統合帯が 0.02 Å / 1e-5 Eh から広がった後も、誤って統合する危険は増えていません。誤りうる方向は過分割（重複 basin）だけです（U3-N6）。
- 実エンジンの mode-follow（プローブ `/home/user/hfauto_review2_probe/u3_mf`）: HONO の cis-trans TS（v2 の saddle の最終構造）から minima(dft) を実行しました。opt は 2 点で TS に留まり（NWChem は収束済みでも 1 ステップ進む）、freq で saddle（P.Frequency −681.67）、+0.1 Å 側の opt は 26 点で 13.3 kcal/mol 下りました。−側と両側の settle は、負荷（load 約 10）のため 10 分の制限内に終わりませんでした。
- 残る妥当性の課題はすべて low です（rerank の SP 失敗で単量体が黙って脱落しうる、reaction-paths の `_register` が pipeline の `gates:` を無視する、design.md:142 の文言の曖昧さ）。
- ライブラリの仕様:
  - NWChem DRIVER: 閾値の既定は GMAX 0.00045 / GRMS 0.0003 / XMAX 0.0018 / XRMS 0.0012、MAXITER の既定 40、TRUST の既定 0.3。INHESS 0 は restart データがあれば使い、2 は restart データを優先しなければ freq の Cartesian Hessian を変換する。LADDER の継続（`inhess 2` を落として job.drv.hess をコピー）は既定の INHESS 0 で restart データを使うので、更新済みの Hessian は引き継がれる。
  - xTB: vtight は Econv 1e-7 Eh・Gconv 2e-4 Eh/α。未収束は `FAILED TO CONVERGE` と空の NOT_CONVERGED ファイル。ANCopt は入力 Hessian を読まないので、XTBEngine.optimize が init_hessian を捨てるのは仕様どおり。

#### 3.3.3 有用性（◎）

- S1 で、同じ stage の中で同じ basin に落ちた構造は known になり、freq を省きます。
  - A/B arm A: tma_hf2_afir が `known:basin_tma_hf2_705342ba` になり、freq 4 → 3 本（17 原子で約 890 s 節約）、dft stage 3,317 → 2,433 s（−27%）。
  - NH3（dft stage 18 → 12 s）と水（14 → 9 s）: DFT freq 2 → 1 本。
  - TMA·(HF)₂ の screen: xTB freq のジョブ 8 → 4 本（manifest の freq calc では 7 → 3）、known 6 件（hf_c00、shared_proton、tma_c00、tma_hf2_c00/c01/c02）。
- TMA·(HF)₂ の dft stage 3,050 → 1,771 s は、主に M3 によるものです（偽の AFIR 生成物の opt 404.7 s と freq 902.9 s が消えた）。この run の dft では S1 は発火していません。
- screen の直列化による増分は 2 s で、無視できます。
- 残る無駄:
  - rerank_sp が per_state 以下の組にも DFT SP を出す（TMA で 3 本・50.2 s、dft stage の 2.8%。選択を 1 度も変えていない。U3-N1）。
  - DFT の mode-follow で TS 候補ができると、両側を settle し直して freq をもう 1 本ずつ取る（U3-N3、実 run で 0 回）。
- 時間の大半は DFT の解析 Hessian です（17 原子で 902 s、CPHF 5 反復で約 62%。dft stage の 66%、discover 全体の 48%）。同じ PES と numerics を保つ限り、削る余地はありません。

#### 3.3.4 複雑さ（○、前回 △）

- minima.py 248 → 235 行、minimum.py 294 → 277 行、identity.py −7 行。W3 全体の hfauto は −12 行（S7 の +16 行が相殺、見込みの −40〜60 行より小さいが構造は明確に単純になった）。
- 経路は直列の 1 本、基準は 1 つ。compare_minima、Verdict、ambiguous、rejudge、tight 引数は hfauto/ に残っていません（classification.py の Verdict は別の無関係な型）。
- 残る複雑さ: C10（暗黙の分岐）、rerank の経路（always の候補は 1 回目の選択では per_state の枠を使い、rerank では使わない）、`_known` の version の自己比較（同梱 pipeline では常に空）、relax 中の find は組成・level で絞らず add は絞る非対称（エネルギー条件で分離されるので維持）。

#### 3.3.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U3-ISS-5（前回） | low | include=window で screen 極小が 0 件なら黙って空、level=screen で select を黙って無視（C10） |
| U3-ISS-7（前回） | low | relax 中の known.find は組成・level で絞らない。前回の「変更不要」を維持 |
| U3-N1 | low | rerank_sp が per_state 以下の組にも DFT SP を出し、SP が失敗すると非 always の候補（単量体の最低構造を含む）が黙って脱落する（`minima.py:115-120, 204-209`、`selection.py:57-65` の rerank）。`single_points` は失敗を diagnostics に残さない。純関数のプローブで tma の SP を欠かすと rerank の結果は ['hf','cx'] になり、ΔG_assoc が理由なく None になりうる（selection.py の行番号は select_for_refinement が :34-54、rerank が :57-65） |
| U3-N2 | low | reaction-paths の `_register`（`actions.py:348-349`）は policy を渡さず、既定の `MinimumPolicy()` を使う。`gates:` を上書きした pipeline では、QRC や中間体で見つかった basin だけが別の閾値で tier 判定される（1 つの PES で基準が 2 つ）。同梱の pipeline は `gates: {}`。f56534d からある不整合 |
| U3-N3 | low | DFT の mode-follow で TS 候補の両側を settle し直すので、側ごとに opt と freq が 1 本ずつ余分にかかる。NWChem は収束済みの構造からでも必ず 1 ステップ進み（プローブ ae76…: step 0 の gmax 2.7e-4 の後 step 1 で xmax 7.3e-4 Å）、指紋は 1e-6 Å で丸めるので freq は必ずミスする。17 原子で 1 候補あたり約 2 × 950 s。実 run で 0 回 |
| U3-N4 | low（一部訂正） | design.md:142 の「開始構造は、screen 極小があればその最適化構造」は all にもかかるように読めるが、実装（`minima.py:201-202`）は include: all なら入力構造から始める。all を使う known_endpoints と A/B には screen の上流がないので、影響は文言の曖昧さだけ。文書側を直す |
| U3-N5 | low | rerank では always の候補が per_state の枠を使わない（プローブ: rerank なし 3 構造、あり 4 構造）。増えるのは always を含む組ごとに最大 1 構造 |
| U3-N6 | low | assign のエネルギー判定は RMSD 最近傍の 1 候補だけにかかり、次点との分離は RMSD だけで見る（プローブで確認）。誤る方向は過分割だけで、実 run で未観測。設計上の性質として記録する |

#### 3.3.6 追加の改良案

- **U3-R1 rerank_sp の DFT SP を、per_state を超える組だけに出す** — could、修正して採用、中立【C32】
  - `minima._jobs` の中で、非 always の候補が per_state を超える (組成, 状態ラベル) の組だけを rerank に渡し、それ以外の組は素通しにする（4〜5 行。selection.py に `only=` のような引数は足さない）。SP の失敗は `history[sid] = ['rerank_sp_failed:<kind>']` として diagnostics.json に残す。テスト 1 件（1 構造の組には SP を出さず、脱落もしない）。
  - 効果: 1 構造だけの組（単量体の最低構造など）が SP の失敗で黙って消えなくなり、docstring と design.md:142 の約束（入力単量体の最低構造は必ず含まれる）が常に成り立つ。TMA·(HF)₂ の SP 3 本・50 s もなくなる。混んでいる組の再順位付けはそのまま残る。
- **U3-R2 暗黙の分岐を明示し、include: all の記述を直す（C10）** — could【C10】
  - `_jobs` で level=dft かつ include=window のとき、上流に screen の MinimumRecord が 0 件なら `ValueError('include: window needs screen minima; use include: all')`（3 行）。MinimaConfig の model_validator で `level == 'screen' and 'select' in model_fields_set` を読み込み時に拒否する（2〜3 行）。design.md:142 を「all は入力 species 全部を入力構造から」に直す（コードは直さない）。テスト 1 件。
- **U3-R3 reaction-paths の `_register` に pipeline の gates を渡す** — could【C33】
  - `actions.py:348` の呼び出しに `policy=MinimumPolicy(gates=ctx.policy.gates)` を加える（1 行と import）。mode_follow は既定の 2 のまま。1 つの PES の tier 基準が 1 つにそろう。

#### 3.3.7 見送った案

| 案 | 理由 |
|---|---|
| mode-follow と soft 押しの opt に手元の freq の Hessian を渡す | DFT での発火は 0 回。HONO のプローブでは Hessian なしの +側 opt が 26 点で収束し、効果の根拠がない |
| mode-follow の両側を settle し直さず、`_follow` の opt と freq をそのまま登録する（U3-N3） | 17 原子で 1 候補約 30 分の節約だが、発火 0 回で MinimumOutcome の API が変わる。DFT の mode_follow discovery が実 run で出てから入れればよい |
| assign で RMSD の閾値内の全候補にエネルギー条件をかける（U3-N6） | 起きうるのは過分割だけで未観測。唯一の基準に分岐を足すことになる |
| relax 中の known.find も組成・level で絞る | 5e-5 Eh のエネルギー条件で分離される。前回の見送りを維持 |
| `_known` の version の自己比較を直す、`_known` を削除する | 1 つの run の中で version は変わらない。直すには engine の site の公開が要る |
| basin の代表を最低エネルギーの member にする | member どうしの差は 5e-5 Eh 以下 |
| screen の xTB relax の並列化を戻す | stage の壁時計時間は 1〜2 s。S1 の単純さの方が価値がある |
| rerank でも always の候補を per_state の枠に数える（U3-N5） | 増えるのは組ごとに最大 1 構造で、S6 以後は生成物が別の状態ラベルになることが多い |
| DFT freq を速くする（grid・CPHF・RI の変更、movecs の再利用） | same_pes が崩れるか PES が変わる。movecs の再利用は約 4% |
| minima の opt の trust を 0.3 に戻す、閾値を緩める | 初期 Hessian のない opt では根拠がない（初期 Hessian を渡す opt は S21 で扱う） |
| level=screen でも window の選択を効かせる | screen は explore の出発点として全件を緩和する必要がある。U3-R2 の拒否で足りる |

#### 3.3.8 変更不要な点

- opt と別ジョブの解析 freq を同じ numerics の層で計算し、is_minimum だけで判定すること。質量加重の Eckart 射影（ランク判定 1e-3）。
- 虚振動の 3 段方針、soft の片側 0.1 Å 押しを 1 回、±0.1 Å の mode-follow を最大 2 サイクル。
- S1 の直列の settle と、同一性の基準を assign の 1 つにしたこと。
- 2 断片以上のときだけ xTB の Hessian を init_hessian にすること（screen freq の JobStore が hit し、opt の点数が約 −33%）。
- M3: dft と explore を screen の最適化構造から始めること。window 6 kcal/mol・per_state 3・rerank_top 8。
- NWChem の opt デッキと freq デッキ、LADDER。xTB の `--opt vtight` と `--hess`、収束の判定。M6。
- コードに固定された閾値（変位 0.1 Å、同一性、OPT_MAXITER、HESSIAN_NEAR_A）を knob にしないこと。

---

### 3.4 U4 反応探索（explore）

評価: 前回 △ / △ / ○ → 今回 **○ / ○ / ○**（有用性の ○ は下限。対象系の実 run で採用した生成物は 0 件、NT2 の slot の約 6 割が同じイリド TS、陰性結果は件数しか残らない。C12 と C14′ を入れれば ○ が安定する。◎ にはしない）

#### 3.4.1 現状仕様の要点（変更点を中心に）

- 【新】**M3**（`stages/explore.py:66-68, 148-150`）: 出発構造は screen 極小の `inputs.evidence(minimum.opt_calc).final`（xTB `--opt vtight` の最終構造）。`trials.generate` と `perturb_linear` もこの構造に適用します。
- 【新】**S6**（`backends/readuct/worker.py:49-51`）: `matches_source` は状態ラベルの一致（元素付きグラフの WL hash と断片の組成式）。SOURCE_RMSD_A と RMSD の import を削除。NT2 の `irc_product` と AFIR の same_as_source の両方に効きます。
- 【新】**M1**: negative_evidence の一式を削除（`_Pool.negatives`、`ReactionRecord.negative_evidence`、verdict `negative_evidence`、`CaseOutcome.NO_PRODUCT`、`CasePolicy.override_negative_evidence`、screen の早期 return）。explore の出力（DiscoveryRecord のスキーマと理由）は変わらず、negative は `reporting/summary.coverage` の件数に使うだけです。
- 変わらない部分: trial 生成（polar_h の relay → 単独 → h_shift → heavy_bond → association、重複除去は原子番号の集合だけ、先頭 10 件、切り捨て件数は記録しない）。NT2 → Bofill TSOpt（automatic_mode_selection）→ 射影振動数で ν < −50 cm⁻¹ がちょうど 1 本 → IRC → 両端の opt と n_imag=0。AFIR（γ 125 → 300、1 原子対）。窓（ΔE‡ 150、ΔE_rxn 100 kJ/mol、None は窓外）。SCC 失敗の 1000 K 再試行。DiscoverySettings。timeout 600 s。版数の pin 6.1.0。ReaDuct の暗黙の反復上限 150（NT2 は 500）。
- 実物の呼び出し（v2 `jobs/e3/e30f4cd2…`）: `run_nt2_task(..., nt_associations=[14,15,13,16], nt_dissociations=[15,16,1,13], nt_total_force_norm=0.1)` → 'Found TS guess after 37 iterations' → `run_tsopt_task(optimizer='bofill', automatic_mode_selection=[1,13,14,15,16])` → ν1 = −289.46 cm⁻¹ → IRC（両方向とも 150 フレームで打ち切り）→ 両端とも source と同じラベル → same_as_source。
- 未実施: C12、C13、C14。

#### 3.4.2 妥当性（○）

- **M1**: 山越えの定理に反する拒否権がなくなりました。hfauto/・configs/・tests/ に negative_evidence・NO_PRODUCT・no_product_basin・override_negative の残りはなく、`reaction_paths.py:52-57` は旧キーを unknown として拒否します。
- **M3**: v2 の 71 attempt すべてで worker の source の opt が 2 反復・ΔE 0.000 kJ/mol でした（W7 は 14〜38 反復、−2.9〜−33.2 kJ/mol）。screen の xtb 6.7.1 と scine-xtb-wrapper 3.0.2 の GFN2 の PES が一致していることも確かめられました。
- **S6** は構成異性の同一性そのもので、hfauto の basin の定義と揃います。
- 実データの検証:
  - TMA·(HF)₂ の relay の TS（−289i）の IRC の後退端は、source の F14↔F16 などを置換したコピーでした（置換不変 RMSD 0.002 Å、恒等写像 0.894 Å）。つまり縮退した HF 交換で、same_as_source は化学的に正しい判定です。前回の基準（置換後の結合も比べるので置換不変）でも M3 だけで同じ判定になります（validation.md:54 は S6 の効果と書いている、U4-N1）。
  - TMA 単量体の 10 trial の結論は、v2 の 2 run（tma_hf2、amine）で完全に一致しました（9 件が 376.4 kJ/mol の out_of_window、1 件が −589i の irc_not_connected）。W7 では run ごとに変わっていました。
  - S6 が判定を変えた実例: TMA·HF の h_shift（+[1,12] −[3,12]）で、IRC の source 側の端点は N···H 2.94 Å の別配座（RMSD 1.055 Å、同じグラフ）でした。これが source と一致と判定され、もう一方が product → out_of_window（397.9 / 162.3 kJ/mol）になりました。W7 は irc_not_connected でした。
- 残る妥当性の制約: 等価な drive が上限 10 件を埋め、heavy_bond の段に一度も届かない（U4-I4）。F→N 移動の drive が作られない（U4-I8、設計上の判断）。状態ラベルの閾値 1.479 Å（N–H）の近くに H 結合がある（TMA·(HF)₂ の N···H 1.413 Å、差 −0.066 Å。v2 の 107 端点でラベルの反転は 0 件、U4-N2）。
- ライブラリの仕様:
  - `NtOptimization2Task`（ReaDuct 6.1.0）: stop_on_error=True では 'No transition state guess was found in Newton Trajectory scan' を含むすべての例外が再送出され、`.nt.failed.xyz` が書かれる。v2 の monotonic_uphill 4 件をプローブで再実行すると、すべてこのメッセージだった。
  - `NtOptimizer2`（Utilities 10.1.0）: 抽出の `last_maximum_before_first_target` は、Savitzky–Golay 平滑化したエネルギーで、最初の target に達する前の極大を探す。したがって no-guess は「target に達するまで平滑化後に極大がない」ことを意味し、monotonic_uphill の名前は実質的に正しい。
  - `GradientBasedCheck`（Utilities 10.1.0）: 既定 maxIter=150。v2 の IRC 66/66 本と端点の opt 6 件が 150 で打ち切られたことと整合する。

#### 3.4.3 有用性（○、下限）

- この単位の「実害」は下流の DFT に偽物を流すことでしたが、それは消えました。W7 の偽の AFIR 生成物 1 件（DFT の opt と freq で 1,308 s）は v2 で 0 件になりました。
- explore 自体は 80〜100 s、1 attempt 0.4〜8 s で、SCC の再試行・timeout・failed はいずれも 0 件。HCN の smoke（NT2 で HNC）は v2 で PASSED。
- 限界:
  - 対象系（amine·(HF)n）の実 run で採用した生成物は、W7・v2 とも 0 件。
  - v2 では NT2 の 40 slot のうち 24（TMA 9、TMA·HF 9、TMA·(HF)₂ 6）と、それぞれの AFIR のフォールバックが、同じ C→N 1,2-H 移動（イリド、TS −1517〜−1519i、ΔE‡ 376〜480 kJ/mol）に使われた。
  - TS がある negative 14 件（TMA·(HF)₂ の縮退 HF 交換 41.5 kJ/mol、NH3·HF の H 交換 138 kJ/mol ×3 を含む）は ΔE‡ を持たない。切り捨てた trial の件数（23〜41）も出ない。「窓内の低障壁反応はない」という陰性の結論の情報量は小さい。
- trial の網羅性（HEAD の trials.generate を v2 の出発構造に適用したプローブ）:

| source | 列挙した drive（固有） | 実行 | 切り捨て | WL 原子クラスでの固有数 | 実行した 10 件がカバーするクラス |
|---|---:|---:|---:|---:|---:|
| TMA | 33 | 10 | 23 | 5 | 2 |
| TMA·HF | 35 | 10 | 25 | 7 | 2 |
| TMA·(HF)₂ | 51 | 10 | 41 | 15 | 5 |
| NH3·HF | 5 | 5 | 0 | 3 | 3 |

- それでも、安価に「窓内の共有結合の反応がない」ことを示し、真の生成物があれば TS 付きで DFT に渡せます。役割としては妥当です。

#### 3.4.4 複雑さ（○）

- M1 で約 45 行・設定 1 個・verdict 1 個・CaseOutcome 1 個、S6 で定数 1 個と import 2 個が消え、増えたのは M3 の約 3 行だけです。knob（sources_per_state、max_trials_per_source、DiscoverySettings 8 項目、window）は変わりません。explore が CALCULATION を consumes に宣言せずに `inputs.evidence` を読む点は、minima._pool と同じ既存の慣習です。
- 残る複雑さ: ReaDuct の反復上限が暗黙のまま（v2 で IRC 66/66 本が打ち切り）、陰性理由は約 10 種あるが件数にしか使わない、relaxation_discoveries は 4 run とも 0 件（約 15 行、安価なので残す）。

#### 3.4.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U4-I4（前回） | low | 等価な drive で上限 10 件が埋まり、heavy_bond の段に届かない。列挙数・固有数・実行数も記録されない（C12） |
| U4-I5（前回） | low | ReaDuct の反復上限 150。amine の irc_end_not_minimum 8 件中 6 件は端点の opt が 150 に達した数値的な陰性（`all()` の短絡で、片側が失敗するともう片側は opt されない例も確認）。同じ −1519i の TS が、等価な drive ごとに別のラベルに分かれる（C13） |
| U4-I6（前回） | low | 陰性ラベルの粗さ。前半（NT2 の失敗メッセージによる分類）は情報を増やさないと判明した。後半（TS がある negative の ΔE‡）と、AFIR の reason の上書き（U4-N3）が残る。M1 の後は報告の精度だけの問題 |
| U4-N1 | low | validation.md:54 は −289i の TS の判定の変化を「グラフの一致で same_as_source」（S6 の効果）と書いているが、実際には M3 の効果で、しかも縮退反応だった（独自の再計算: irc_f_opt は置換不変 0.001 Å、irc_b_opt は置換不変 0.0016 Å・恒等写像 0.894 Å）。文書の因果の書き方だけの問題 |
| U4-N2 | low | S6 の残存リスク: 状態ラベルの閾値の近くの H 結合の伸縮で、配座違いが「生成物」になりうる（TMA·(HF)₂ の N···H は閾値との差 −0.066 Å、TMA·HF は +0.07 Å）。v2 の 107 端点で反転 0 件。起きた場合も S1 で既知の basin に合流して freq は省かれ、被害は DFT opt 1 本（約 7 分）。前回 P3 の注記がそのまま残っている |
| U4-N3 | low（一部訂正） | AFIR の reason は γ ごとに上書きされ、最後の γ の結果だけが残る。v2 の afir_not_converged 2 件は、どちらも γ=125 が 93 反復で収束して出発構造に戻り、γ=300 が 150 反復で打ち切られたもの。「S6 の後は」という限定は不正確で、同じループは S6 の前からある（S6 が変えたのは same_as_source になる頻度だけ）。U4-I6 の一部として扱う |
| 数値の訂正 | — | 仕様書の「実 run の 142 attempt すべてで source の opt が 2 反復」は誤りで、2 反復になるのは v2 の 71 attempt（W7 は 14〜38 反復） |

#### 3.4.6 追加の改良案

- **U4-R1 drive の重複を WL の原子クラスで除き、列挙数・固有数・実行数を記録する（C12、U4-R2-1 と U4-R2-2 の統合）** — could、修正して採用、中立【C12】
  - topology の内部関数 `_wl_hash` を拡張し、原子ごとのラベル（state_label と同じ反復回数）を返すヘルパーにする。新しいグラフ正準化の仕組みは作らない。
  - `trials.generate` の重複除去の鍵を、原子クラスで正準化した (frozenset(assoc), frozenset(dissoc)) にする。距離の近い順を保ち、同じクラスでは最も近い drive を代表 1 件として残す。WL クラスが同じでも幾何的には等価でない drive（TMA·HF の −1519i と −1395i は別の TS に入った。どちらも窓外の高障壁）を落とすことを許容する、と docstring に書く。
  - 同じ変更で、explore の stage_dir に `diagnostics.json`（{minimum_id: {enumerated, unique, run}}）を書く（conformers と同じ慣習、conformer_search.py:238-239）。スキーマと JobStore のキーは変えない。
  - 効果: v2 の対象系では NT2 の 40 slot のうち 24 とその AFIR 24 件が同じイリド TS に使われていた。1 クラス 1 件にすると TMA 系の explore の attempt は約半分になり、heavy_bond と association の段に届く。クラス数が上限を超える TMA·(HF)₂（15 > 10）では、切り捨て件数が見えることにも意味が残る。15〜25 行とテスト 2 件。
- **U4-R2 無バイアスの opt に `convergence_max_iterations=500` を明示する（C13）** — could、中立【C13】
  - `worker.py` の `_Scine.minimum` の `run_opt_task` だけに渡す（1 行、設定の項目にはしない）。IRC 端点の opt と AFIR 後のバイアスなし opt の両方に効き、source の opt（task を直接呼ぶ）には影響しない。IRC と AFIR は変えない。
  - 効果: amine の irc_end_not_minimum 8 件中 6 件の数値的な陰性と、同じ TS のラベルが分かれる雑音がなくなる。1 attempt の最長が約 8 s → 約 17 s。
- **U4-R3 TS がある negative にも ΔE‡ を載せる（C14 の後半 = C14′）** — could、修正して採用、中立【C14′】
  - `_nt2` で TS が得られた後の negative（irc_end_not_minimum、same_as_source、irc_not_connected）に `energies={'ts': …}` を付け、`result_dict` の `relative('ts')` は outcome によらず ts と source があれば計算する。reaction_kj_mol は product のときだけのまま。
  - 修正点 (1): SCC を 1000 K でやり直した attempt で negative になった場合、`_attempt` は SP をやり直さない（product のときだけ）。そのままでは `records.py:83` の「evaluated at 300 K」に反するので、やり直し温度のときは ts がある場合も 300 K の SP で評価し直すか、基準温度のときだけ載せる。
  - 修正点 (2): 値は manifest の DiscoveryRecord に残るだけで、coverage.csv には出ない。報告に列を足す拡張はしない。
  - 効果: v2 の TS がある negative 14 件（縮退 HF 交換 41.5 kJ/mol、NH3·HF の H 交換 138 kJ/mol ×3、イリド TS 約 398 kJ/mol など）が manifest に定量的に残り、「窓内の反応はない」という陰性の結論を障壁の値で裏付けられる（前回 X1 が推奨していた依存関係）。約 4〜8 行。

#### 3.4.7 見送った案

| 案 | 理由 |
|---|---|
| C14 の前半: monotonic_uphill を no_ts_guess と nt2_failed に分ける | v2 の 4 件の再実行はすべて 'No transition state guess was found…' で、ラベルは実質的に正確。分けても情報は増えない |
| 縮退した転位（IRC の両端が source の置換コピー）を仮説や report の行にする | 新しい極小を生まず、「異なる極小の間の反応」という目的から外れる。値は U4-R3 で manifest に残る |
| NT2 が same_as_source の検証済み TS を得たら AFIR を省く | 1 件約 1 s で、規則が増えるだけ |
| M3 の後、worker の source の opt を削除する | 直線分子の摂動後の start には必要。PES の一致の確認にもなる。0.1 s 未満 |
| matches_source に RMSD かエネルギーの許容幅を併用する | S6 が解消した配座違いの問題を再び持ち込む。反転は 0 件で、起きても S1 で既知の basin に合流する |
| 1.45 Σr_cov を変える、H 結合を考慮したラベルにする | topology・状態ラベル・basin の定義全体に波及する。本題は DFT の same_basin で決着している |
| irc_not_connected の両端を新しい species として登録する（Chemoton 風） | 新しい経路と DFT のコストが要る。対象系ではイリド系の高障壁だけ |
| IRC の反復上限を上げる | 端点は必ず opt し直すので打ち切りに実害はない |
| max_trials_per_source か sources_per_state を増やす | 増えた slot も等価な drive で埋まる。U4-R1 の方が効果的で安い |
| h_shift に原子価のフィルタを付ける | 化学的な規則が増える。U4-R1 で slot の浪費はほぼ解消する |
| AFIR の γ を 1 つにする | γ=300 の害の実証はない（前回の棄却理由のまま） |
| worker と engine の docstring の「design §6.3」を直す | 他のアダプタと共通の、元の設計書の節番号を指す慣習 |

#### 3.4.8 変更不要な点

- M3（screen の最適化構造から出発）と worker の source の opt。S6（状態ラベルで一致判定）。M1（陰性結果は報告だけ）。
- NT2 → Bofill TSOpt → 射影振動数で ν < −50 cm⁻¹ が 1 本 → IRC → 両端の opt と n_imag=0 という検証の流れ。
- product_verdict（NT2 は TS と IRC の検証が必須、ΔE_rxn が None なら窓外、窓 150 / 100 kJ/mol、out_of_window なら AFIR を省く）。
- 単一原子対の AFIR（γ 125 → 300）、SCC 失敗の 1000 K 再試行と 300 K の SP、restricted_open_shell。
- trial の段の順、直線分子の摂動、seed の固定、sources_per_state 2、max_trials_per_source 10、DiscoverySettings の既定値。
- 1 attempt を 1 worker サブプロセスで実行すること、timeout 600 s、版数の pin、continuation なし。relaxation_discoveries。

---

### 3.5 U5 反応経路：仮説選択と障壁の事前判定（hypotheses.select / SCREEN / FIND_PATH）

評価: 前回 SCREEN ○ / FIND_PATH × / 複雑さ △ → 今回 **妥当性 ○（SCREEN ○、FIND_PATH ○−）/ 有用性 ○ / 複雑さ ○**（M8 を入れるまで、2 チャンク以上の string による monotonic → BARRIERLESS は △。M8 の後も、1 チャンク内の端点のずれは note による監視になる）

#### 3.5.1 現状仕様の要点（変更点を中心に）

- **仮説選択**（`chemistry/hypotheses.py:220-249`、変更なし）: 宣言反応（F1 で端点の basin を DFT basin に写す）→ explore の product（TS あり → TS なし → mode_follow）→ conformer ペア。両方 dft・同じ level_key・同じ組成・別 basin、ΔE ≤ 40 kcal/mol、結合変化かねじれ ≥ 30° か距離変化 ≥ 0.2 Å、組成あたり 6 件。【新】M1: 陰性の discovery は仮説にも結論にも影響しない（docstring 7 行目）。
- **決定表 17 行**（`drivers/reaction_case/state.py:217-226`）: 【新】行 12 は barrierless だけで BARRIERLESS（negative_evidence の分岐を削除）。【新】行 14 は monotonic な DFT string 1 本で BARRIERLESS（S3、15 bead の確認と max_path_runs を削除）。
- **SCREEN**（`actions.py:184-259`）: 低レベル TS の近道 → DFT 端点の xTB 再最適化 → 崩壊（≤ 0.05 Å）なら IDPP 11 点 → pysis_gs（11 ノード、climb、rsirfo/gau）→ 【新】S2 `align_sequential`（`interpolation.py:32-39`、`actions.py:242`）→ 全ノードで DFT SP → `barrier_verdict`（内部点の最大 − 両端の高い方。DFT か低レベルのどちらかが ≥ 1 kcal/mol なら proceed）。【新】S4: tsopt 付きの run が例外なら、tsopt を外して GS を 1 回再実行し、ts=None と tsopt_error を返す（`pysis/worker.py:18-35`）。tangent_mode_cm1 のため、高い側の端点の DFT freq を proceed でも毎回要求する（C15 は未実施）。
- **FIND_PATH**（`actions.py:446-489`）: 初期経路は work.initial（整列済みの GS 経路か SCREEN の IDPP）、なければ新しい IDPP 9 bead。NWChem ZTS を最大 3 チャンク（`nbeads 9; maxiter 20; stepsize 0.05; interpol 3; tol 1e-5; freeze1; freezeN; impose`）。
  - 【新】M2: monitor はなく（`nwchem/engine.py:170-171`）、チャンクは maxiter まで走る。打ち切りは `profile.energies_settled`（直近 3 反復で全 bead の |ΔE| < 0.1 kcal/mol、`profile.py:78-87`）。エネルギー履歴は 'string: Path Energy #n' から読む（`output.py:206`）。形は収束を問わず記録し（minimax の上界）、single_max のときだけ HEI を path_hei の種にする。gmax と NWChem の converged は判定に使わない。
  - チャンクの継続は `initial = run.images`（`actions.py:476`）で、**端点を DFT 極小に戻していない**。
- 【新】S5: FakePath は「gmax は横ばい、bead エネルギーは数反復で落ち着く」（`tests/fakes.py:268-304`）。PBE0/STO-3G の HCN string の opt-in smoke（`tests/smoke/test_real_nwchem.py:85-101`）。FakePath は initial_path を無視し、毎回 start から end を線形補間する。
- 削除: negative_evidence の一式、`string_converged`、`stagnated`、`_XMAX_ZERO`、`_string_monitor`、`string_history`、policy の confirm_barrierless_beads・max_path_runs・override_negative_evidence（指定すると読み込み時に拒否）。
- 未実施: C15、C16。

#### 3.5.2 妥当性（○、FIND_PATH ○−）

- **仮説選択と SCREEN の判定式は変わらず、妥当です。** S2 で DLC の画像に残っていた回転の跳びが消えました（W7 の HONO の画像で max_node_spacing 1.149 → 0.190 Å、オフラインの確認）。S4 は、pysisyphus 1.0.0 の `run.py` が tsopt の周りで `HEIIsFirstOrLastException` しか捕まえないことへの正しい対処で、プローブ（cart、`RunAfterCalculationFailedException`）でも ts=None と tsopt_error で回収できました。
- **FIND_PATH の骨格（M2・S3）は理論的に正しい。** 射影しない gmax は MEP 上でも 0 にならないので、打ち切りは bead エネルギーの安定だけで決めます。連続な経路の最大値は鞍点の上界なので、形は収束と無関係に記録してよい。実データでも機能しました（HONO は 9 反復目に max|ΔE| 0.037 で安定。STO-3G の HCN は 1 チャンク目が 0.11〜0.29 で落ち着かず、2 チャンク目の 7 反復目に安定）。
- **【U5-N1、high】minimax の上界は、端点が本物の極小である連続経路にしか成り立ちません。** NWChem 7.2.3 の ZTS は `freezeN` の bead N の座標を動かします。
  - 剛体の整列: `zts_min_motion` を i=2..nbeads に適用する（`string.F:1965-1968, 2005-2008`、コメントは 'Can rotate/translater frozen N bead'）。
  - L-BFGS の第 2 ループ（`neb_utils.F:1302`）は (α_k−β_k)(x_{k+1}−x_k) を足すので、凍結 bead の剛体移動が線形に外挿され、結合長が変わりうる（仕組みは推定）。
  - nstep > 0 では、凍結 bead のエネルギーと勾配を計算し直さない（`string.F:1829-1831, 1886-1900`）。報告される bead N のエネルギーは古い値のまま。
  - 実データ: v2 の find_path_sto3g の 1 チャンク目（19713c6076）で、最後の bead の C–N が 1.2031 → 1.2799 Å、N–H が 1.0333 → 1.0993 Å に変わり、報告エネルギーは −92.068363 のまま。これが 2 チャンク目（33918994c1）の初期経路の端になり、NWChem は −92.054267（+8.85 kcal/mol、隣の bead より 5.0 高い）を凍結端点として計算しました。bead 1 は不変です。
  - S3 で monotonic な string 1 本を BARRIERLESS の十分な証拠にしたので、端点の持ち上がりは、その側にある障壁（ずれの大きさ以下）を隠して誤った最終分類を生みえます。標的系（amine·(HF)n の小さな H 移動の障壁）では、ずれと障壁が同程度になりえます。
  - 1 チャンクの中でも座標はずれます。報告エネルギーは古い真の極小の値なので形の判定には効きませんが、ctx.work.path と端に近い HEI の接線はずれた構造を使います。
- **【U5-N2、medium】IDPP は反復上限 500 で黙って打ち切られ、平面の端点では面内の反転経路になります。** 前回・仕様書の「対称性で閉じ込められる」という診断は訂正が必要で、原因は反復上限です。
  - `_FORCE_TOL` に達するまでの反復数は、HONO の DFT 端点で 9 点 1,479、11 点 1,846、HCN で 2,170。これまでの実計算の IDPP はすべて未収束で返っていました。
  - 500 反復では、HONO（平面性の特異値 6.7e-7、1.2e-8 Å）の経路は H–O–N 最大 174° の面内経路（二面角は 180 → 178 → 8 → 0 と跳ぶ）。5,000 反復（1,479 反復で収束）では、H–O–N 103〜108°、二面角 180 → 0° の滑らかなねじれ経路になります（検証段でも再現）。
  - DFT//IDPP（PBE0-D3BJ/def2-SVPD）の最大値: ねじれ経路 14.82、面内の初期経路 37.04、20 反復後の ZTS 34.37、真の TS 13.67。
  - Cartesian の ZTS は、平面の初期経路を対称性で面内に保つので、IDPP を収束させるかどうかが FIND_PATH の捉える機構の正しさを直接決めます。ただし IDPP が使われるのは退避経路（GS の失敗、SCREEN のない pipeline、SCREEN の崩壊）だけなので medium です。
- **【U5-N3、medium】DLC の GS は平面 4 原子で数値的に脆い。** v2 の HONO（6f3b3760ea）では、pysisyphus.log に 'Rebuilt internal coordinates!' が 2 回出て（internal_coords.log は 'Not enough internal coordinates ... Found only 5'）、H–O–N の変角の primitive が落ちて DLC が 6 本から 5 本になり、`StringOptimizer.get_step` が (66,)/(55,) の形状不一致で例外を出しました。W7 の入力（1e-6 Å の差）では完走しています。S4 の再実行は決定的に同じ例外になりました。結論は誤りません（unavailable → FIND_PATH）が、再現性とコストの問題です。証拠は 1 系なので medium は上限側です。
- **残る制約**: monotonic の判定は 9 個の離散 bead によるもので、「解像度 1 kcal/mol の意味で barrierless」の規約で受け入れます（前回 U5-I6 の結論のまま）。
- U5 のコードに回帰はありません。HONO の GS の失敗は、S1 で DFT 極小の座標が 1e-6 Å の桁で変わったことが引き金になり、既存の DLC の脆さが表に出たものです。

#### 3.5.3 有用性（○）

- **FIND_PATH は × → ○ になりました。** v2 の HONO では SCREEN の GS が失敗しましたが、行 16 の FIND_PATH が救い、W7 と同じ ΔE‡ 13.67 の TS（−678.1i）に到達しました。find_path_sto3g も single_max → elementary_step で閉じました。
- **ただし効率は十分ではありません。** HONO では U5 の区間が 48 → 208 s、U6 の精密化が 10.2 → 121 s に増えました。面内の HEI は 2 次の鞍点の領域（DFT −1327.7 / −1037.9、xTB −1378 / −775 cm⁻¹）にあるので、xTB Hessian は不採用になり、DFT Hessian（15.1 s）と AUTOZ の失敗（48.4 s）からの再試行が入りました。収束させた IDPP のねじれ経路では、HEI の候補（bead 4 と 5）の xTB の虚振動がちょうど 1 本（−745、−794 cm⁻¹）だったので、S18 を入れれば xTB Hessian が使え、この追加コストの大半は消える見込みです。
- **SCREEN** は HCN（28 s）と NH3（29 s）で W7 と同じでした。NH3 で余分に走る tangent 用の DFT freq（5.1 s）の原因も分かりました: 縮退反応の端点は代表構造の原子を置換したもの（H2 と H3 の入れ替え）なので fingerprint が変わり、JobStore がミスします。17 原子の縮退反応では 1 本約 900 s になります（C15 の価値）。
- bead エネルギーの判定はチャンクの境目でしか行わないので、HONO では安定した後の 11 反復（約 110 s、チャンクの 55%）が無駄になりました。ただし、再開のたびに過渡が入る（HCN の 2 チャンク目の最初の反復で max|ΔE| 2.02 kcal/mol）ので、チャンクを小さくする利得ははっきりせず、今は変えません。
- M1 で増えたコストは、仮説ごとの SCREEN 1 回（小分子で 1 分未満）で、組成あたり 6 件の上限で抑えられています。

#### 3.5.4 複雑さ（○、前回 △）

- 削除: negative_evidence の一式、gmax 系の判定と monitor、15 bead の確認、policy の 3 キー（約 60 行）。追加: energies_settled（10 行）、string_path_energies（5 行）、align_sequential（7 行）、S4（8 行）。差し引きで減っています。決定表の 17 行と「FIND_PATH は 1 ケース最大 1 回」は変わりません。
- 残る改良は、どれも定数 1 つか 3〜8 行で、knob は増えません。
- テストの限界: FakePath は initial_path を無視するので、前のチャンクの経路を次へ渡す処理の不具合（U5-N1）は fake では捕まりません（実 run で初めて見つかった）。site の ExecutionSpec.maxiter が opt・saddle・string に共通という設定の罠も前回のままです。

#### 3.5.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U5-I7（前回） | low | 両端より低い DFT//xTB ノードの検出がない（C16） |
| U5-I8（前回） | low | proceed でも tangent 用の DFT freq を要求する（C15）。縮退・置換した端点で JobStore がミスする実例（NH3） |
| U5-N1 | **high** | チャンクを継続するとき、NWChem が動かした凍結端点（bead N）をそのまま次のチャンクの初期経路に渡す。1 チャンクの中でも bead N の座標はずれる（上記）。唯一の実マルチチャンク run で発生し、S3 の BARRIERLESS の前提を壊す |
| U5-N2 | medium | IDPP が反復上限 500 で黙って打ち切られ、平面の端点では面内の反転経路になる（FIND_PATH の機構を誤り、U6 のコストを約 10 倍にする） |
| U5-N3 | medium | DLC の GS は平面の端点で数値的に脆い（1e-6 Å の差で成否が変わり、S4 は効かない） |
| U5-N4 | low | S4 は GS 自体が例外で終わったときも GS を 1 回やり直す（HONO で 8.0 s の 2 回とも同じ例外）。仕様書の「最初の例外は失われる」は誤りで、stderr.txt には 'During handling of the above exception, another exception occurred' の連鎖として両方の ValueError が残る。失われるのは Failure.reason に最終行しか入らないことと、1 回目の cycle_*.trj の上書きだけ。GS が完了したかは final_geometries.trj の有無で判定できる（W7 の 150 サイクル未収束の run にはあり、v2 の落ちた run にはない） |
| U5-N5 | low | bead エネルギーの安定をチャンクの境目でしか判定しないので、安定後の反復が無駄になる（HONO で約 110 s）。対策の利得はデータ 2 点では不明で、今は変えない |
| U5-N6 | low | design.md:130 の「gmax と NWChem の収束判定は記録するだけ」は、NWChem string では gmax_history が空（記録するのは pysis_gs だけ）なので実装と合わない |
| U5-N7 | low | saddle の試行を使い切った後でも FIND_PATH を実行し（行 16 は path_runs だけを見る）、single_max の種は使われず UNRESOLVED になる。17 原子では最大 3 h の浪費になりうる。実例はない（U6-I7 と同じ）。monotonic・multi_max の結論は有益なので修正は見送る |

#### 3.5.6 追加の改良案

- **U5-R1 チャンクの継ぎ目と最終経路で、凍結端点を DFT 極小に戻す** — must、修正して採用、中立【M8】
  - `actions.find_path` のループで、成功したチャンクの後に `frames = ctx.frames(run.images)` とし、`frames[-1]` を `align_mapped(frames[-1], ctx.ends[1])`（真の極小を、ずれた bead の向きに重ねたもの）で置き換え、`ctx.path_file(f"{run_name}_c{n}", frames)` で書き直したものを次の initial にする。bead 1 は NWChem が動かさないが、対称性のため frames[0] も同じように戻してよい。
  - **最終チャンクの frames の端点も同じく ctx.ends に置き換え**、ずれ（mapped RMSD）を note（例 `string0:end_drift:0.081`）に残す。ゲートにはしない。
  - テスト: initial_path を記録する scripted PathEngine で、1 回目に最後の frame をずらした PathProfile（energies_settled は偽）を返し、2 回目に渡された initial_path の最後の frame が ctx.ends[1] と写像 RMSD ≈ 0 になることを確かめる。
  - 文書: 「NWChem の freezeN は凍結 bead の座標を保証しないので、チャンクの継続と最終経路で端点を DFT 極小に戻す。1 チャンク内のずれは、離散 bead の規約と同じ扱いで note によって監視する」。minimax の上界が完全に回復するのはチャンクの開始時について言えることで、1 チャンク内のずれは残ることを明記する。
  - 効果: 2 チャンク目以降の string は常に真の極小どうしを結び、行 14 の monotonic → BARRIERLESS と形の判定の前提が回復する。追加の計算はない（2 チャンク目の端点のエネルギーは NWChem が nstep 0 で真の極小について計算する）。3〜8 行とテスト約 15 行。
- **U5-R2 IDPP の反復上限を 500 → 5,000 にする** — should、中立【S18】
  - `interpolation._MAX_ITERATIONS = 5000`（`_FORCE_TOL`、`_STEP`、`_MAX_STEP_A`、キックは変えない）。単体テストを 1 件足す: 厳密に平面な cis/trans の 4 原子（HONO の DFT 端点）の IDPP が面外へ出る（内部の H–O–N < 120°、二面角が 180 → 0 と単調に変わる）。
  - 効果: HONO の FIND_PATH の初期経路が面内の反転（DFT//IDPP の最大 37.04、HEI の xTB 虚振動 2 本）からねじれ（最大 14.82、HEI の候補の xTB 虚振動 1 本）に変わる。効果が出るのは、IDPP が初期経路になる退避時（GS の失敗、SCREEN のない pipeline）と SCREEN の崩壊時だけで、GS が成功すれば GS の経路が初期経路になる。
  - コスト: 定数 1 つとテスト 1 件。17 原子・11 点で強制的に 5,000 反復しても 1.1 s（負荷下）。NH3 の単体テスト（間隔 0.071 Å、最短距離 0.997 Å）と G15（1.002 Å）の値は変わらない。
- **U5-R3 pysis_gs の平面端点に off_plane を入れる（off_axis の一般化）** — could、修正して採用、【複雑さ増】【C35】
  - 厳密に平面な端点（原子 4 個以上、重心からの座標の最小特異値 < 1e-4 Å）を、法線方向に ±0.01 Å（原子ごとに交互）ずらす。適用は **dlc になる 1 断片の平面分子に限る**。JobStore のキーは今と同じく、ずらす前の fingerprint。テスト 1 件。約 8 行。
  - 証拠は HONO の 1 系だけ（v2 の入力で 3/3 完走、HEI 11.6 と TS を再現）。「cart の平面錯体でも面外経路を選べる」は未検証の推測。**S18 を先に入れてから採否を再判断する**（S18 だけでも FIND_PATH の退避経路は正しいねじれ機構を捉えるので、利得は主に SCREEN を完走させるコストの面に限られる）。
- **U5-R4 S4 の再実行を、GS が完了した場合に限る** — could、中立【C36】
  - `worker.run_growing_string` の except で `if "tsopt" not in run_dict or not (workdir / "final_geometries.trj").is_file(): raise`。pysisyphus のファイル名への依存が増えるので、コメントで `Optimizer.py:385`（COS の final_fn が final_geometries.trj）を参照させる。既存のテスト（always_raises）に呼び出し 1 回の assert を足す。1〜2 行。
- **U5-R5 文書の修正（design.md と validation.md）** — could【C39 に統合】
  - design.md:130 の「gmax は記録するだけ」を「NWChem の string は gmax を記録しない（記録するのは pysis_gs だけ）」に直す。M8 を入れたら、その 1 行を足す。
  - validation.md §8-1 の対策案（小さな分子の GS を Cartesian にする、1 回だけ再試行する）はプローブで無効と分かった（cart は面内の経路で HEI 26.56 になり TSOpt も例外、再試行は決定的に同じ失敗）ので、S18・C35 に置き換える。§8-2 に「最初の例外は stderr.txt に連鎖として残る」を添える。
- **U5-R6 tangent_mode_cm1 を barrierless のときだけ計算する（C15）** — could、中立【C15】
  - tangent なしで barrier_verdict を呼び、barrierless のときだけ `_tangent_mode_cm1` を計算して below_zpe を付け直す（約 5 行）。17 原子の縮退反応や置換した members の端点で 1 本約 900 s を省く。
- **U5-R7 両端より低い DFT//xTB ノードがあれば barrierless にしない（C16）** — could、修正して採用、【複雑さ増】【C16】
  - min(内部) < min(端点) − 解像度 なら `unavailable('dft_node_below_endpoints')` とし、FIND_PATH → multi_max → VALIDATE_INTERMEDIATE へ回す（3 行）。**M8 の後に入れ、FIND_PATH → multi_max → VALIDATE_INTERMEDIATE の経路を fake で通すテストを必須とする**。実例はない。

#### 3.5.7 見送った案

| 案 | 理由 |
|---|---|
| 平面の小分子では GS を Cartesian にする（validation.md §8-1 の案） | プローブで v2 の HONO は cart で 12 サイクルで収束したが、経路は面内のまま（HEI 26.56、ねじれは 11.6）で TSOpt も例外。機構を誤る |
| GS を単純にもう 1 回試す | 決定的。v2 では 2 回とも同じサイクルで同じ例外 |
| IDPP の乱数キックを大きくする | 0.05・0.1 Å では 500 反復で面内のまま、0.3 Å でやっと面外へ出るが全画像を乱す。S18 の方が根本的 |
| IDPP の端点にも off_plane を入れる | S18 で IDPP を収束させれば既存のキックで面外に出る |
| ZTS を maxiter 10 × 最大 6 チャンクにする | 再開のたびに過渡・L-BFGS 履歴の喪失・SCF のやり直しが入る。データ 2 点 |
| 実行中に string を kill する、print_shift で途中の経路を回収する | string_final.xyz はループの後にしか書かれない（前回の見送り） |
| 端点のずれに閾値のゲートを付ける | M8 で各チャンクの開始時の端点は常に正しい。最終チャンクの中のずれは note で十分 |
| 形の判定で string の端点エネルギーを ctx.energies に差し替える | M8 の後は NWChem が計算する端点のエネルギーが真の極小の値になり、同じになる |
| tsopt_error を PathProfile のフィールドにする | 情報は result.json と stderr.txt に残る。モデル変更に見合わない |
| 行 16 に saddle の試行数のガードを付ける、path_hei に追加の試行を与える（U5-N7） | 実例がない。予算を使い切った後でも monotonic・multi_max の結論は有益 |
| monotonic の確認に bead を増やす（旧 15 bead の復活） | 17 原子で 5〜6 h。規約で受け入れる |
| IDPP の収束フラグを返して note に残す | API の変更の割に、S18 の後は未収束はまれ |
| PathProfile.gmax_history を削除する | pysis_gs の監査に使え、害がない |

#### 3.5.8 変更不要な点

- M1（陰性 discovery は仮説にも結論にも影響させず、coverage の報告だけ）。
- 仮説の条件（40 kcal/mol、0.2 Å / 30°、組成あたり 6 件、優先順）、pick_endpoints、F1、_lend_ts。閾値を YAML に出さない方針。
- SCREEN の判定式、低レベル TS の近道、xTB freq による TS の受理（−50 cm⁻¹ 未満がちょうど 1 本）、_COLLAPSE_A 0.05 Å。S2。
- M2（monitor なし、energies_settled、形を収束と無関係に記録、種は single_max の HEI だけ）。S3（M8 を前提とする）。S4 の骨格（GS 自体の失敗は伝える）。
- 決定表 17 行の構造、FIND_PATH は 1 ケース最大 1 回。S5。既定値（string 9 bead × 20 反復 × 3 チャンク、GS 11 ノード、off_axis）。

---

### 3.6 U6 反応経路：鞍点精密化・TS 検証・接続確認（REFINE_SADDLE / VALIDATE_TS / CONNECT）

評価: 前回 ○ / ◎ / ○ → 今回 **○ / ◎ / ○**（妥当性の留保: M2 で未収束 string の HEI が種として入るようになり、「Cartesian で moddir 0」の組合せでの迷走（U6-N1）と、string 画像の未整列による重なりの過小評価（U6-N2）が実計算で表に出た。どちらも VALIDATE_TS と QRC のゲートが結果の誤りは防いでいる）

#### 3.6.1 現状仕様の要点（変更点を中心に）

- 決定表の行 4〜10、15〜17、REFINE_SADDLE、VALIDATE_TS は変わっていません。
  - REFINE_SADDLE: 1 回目だけ種の xTB Hessian（ν < −10 の負モードがちょうど 1 本で、接線との重なり ≥ 0.3）、それ以外は種で DFT freq。`_mode_index`（`actions.py:262-272`）は ν < −10 の負モードの中で direction と |cos| 最大のもの。render_saddle（`input.py:135-143`）は `trust 0.1`、`sadstp 0.1`、`inhess 2`、maxiter 50、mode_index ≥ 1 なら `moddir k+1` と `noautoz`（**mode_index=0 なら moddir は書かれず、NWChem 既定の 0**）。
  - VALIDATE_TS: 別ジョブの DFT freq を hfauto で射影し、`is_first_order_saddle`（指紋、same_pes(numerics)、本数、最低 ν < −50（ねじれは −20）、2 本目に < −50 がない、prominence ≥ max(2e-5, 20×scf_tol)）。higher_order なら 2 本目のモードに +0.1 Å の種。
- 【新】**S7**（`actions.py:397-398`、`nwchem/engine.py:43, 129-138, 321-332`）: QRC の両側の opt に TS の freq Hessian を init_hessian として渡します。optimize は、同じ元素の並び・同じフレームで各原子 ≤ 0.5 Å（HESSIAN_NEAR_A）の Hessian を受け付けます（saddle は指紋の完全一致のまま）。デッキは `inhess 2`・`trust 0.1`・maxiter 100。key_payload は {hessian: TS Hessian の sha, molecule}（W7 は {hessian: None, tight: False}）。trust 0.3 は見送られました（validation.md §8-4）。
- 【新】**C19**: `gates.qrc_drop`（72-74）を connection（260）と CONNECT（393）が共有。振幅は `min(first × qrc_retry_factor^(attempt−1), qrc_bounds_A[1])`（396）。
- 割付け（Registry.find → periodic_nearest → `_register` がゲートより前）、connection ゲート（traj[0] ≤ E_TS − drop、max(traj) ≤ E_TS − drop、traj[-1] ≤ traj[0] − drop）、分割（max_split_depth 2）、予算（6 h/case、子反応もそれぞれ新しく 6 h）、継続（TIMEOUT/MAXITER 2、autoz 1 は種から Cartesian、SCF 1）は変わっていません。
- 【波及】M2 で未収束 string の HEI（path_hei）が種として実際に入るようになりました。S2 で GS の画像は整列されますが、NWChem string の画像（path_hei の接線の元、`actions.py:481`）は整列されていません。
- 未実施: C17、C18。

#### 3.6.2 妥当性（○）

- 骨格（種 → 近似 Hessian で P-RFO → 別ジョブの DFT freq を射影 → 曲率から振幅を決めた QRC → エネルギー降下と厳格な割付けのゲート）は、TS 探索の標準的なやり方です。
- S7 は正しく動いています。QRC の stdout に 'Using Cartesian Hessian from previous frequency calculation' が出ます。v2 の QRC 8 側（HCN・HONO・NH3・STO-3G の各 2 側）はすべて max(traj) = traj[0] で、変位直後の降下は 2.73〜3.46e-4 Eh（目標 3e-4、drop 1e-5）。trajectory_above_ts は 0 件でした。C19 は単体テスト（`qrc_bounds_A=(0.03, 0.05)` で 2 回目も 0.05）で確認できます。
- 誤った値が黙って順位に入る経路は見つかりませんでした。
- **M2 で表に出た点**:
  1. **U6-N1（moddir 0 の迷走）**: NWChem の moddir 0 は、1 手目に勾配の大きな成分の方向へ登り（`opt_drv.F:3277-3294`、コメントは「制約付き最適化の直後なら理にかなう」）、その後は最大重なりでモードを選びます。HONO v2 の attempt_01（noautoz、moddir なし）は、step 1 で 'Maximum overlap of modes: mode= 5 overlap= 8.2D-01' となり、正の曲率の mode 5 を 8 ステップ追い（'Forcing uphill step in mode 5'）、E は −205.2896 → −205.1587（+0.131 Eh、+82 kcal/mol）まで登ってから戻り、36 ステップで収束しました。autoz の attempt_00 は 2 手目から mode 1 を重なり 1.0 で追っていましたが、34 ステップ目に AUTOZ（`geom_binvr`）で失敗しています。つまり正の曲率のモードを追ったのは「Cartesian で moddir 0（既定）」の組合せで、現状では autoz 失敗後の継続か、site の `coordinates: cartesian` のときに起きます。
     - プローブ（同じ種・同じ DFT Hessian、`/home/user/hfauto_review2_probe/u6_assess`）: noautoz + `moddir 1` は step 1 から mode 1 を重なり 0.99〜1.0 で追い、18 ステップで単調に −205.32453416 に収束しました（本番 −205.32453408）。autoz + moddir 1 は rc=1 で 4 ステップまでしか走っておらず、「attempt_00 と同じ這い進み」の根拠は 2〜4 ステップのエネルギーの一致（差約 3e-7 Eh）だけです。
     - 負モード 1 本の xTB 種（autoz + moddir 1）: HCN 5 → 4、NH3 3 → 2、HONO（W7）5 → 4 ステップで、TS のエネルギーは ≤ 1e-7 Eh で一致。
  2. **U6-N2（重なりの過小評価、一部訂正）**: STO-3G の path_hei の種で、string 画像の中心差分の接線と虚モードの重なりが 0.118 になり、誤った `mode_overlap_below_0.3` の note が付きました。隣の画像を Kabsch で整列すると 0.651 です。主因は string 画像の未整列の剛体回転です。質量加重モードの M^−1/2 q の正規化はその基準振動の物理的な Cartesian 変位で、幾何学的な接線と cos で比べることは妥当なので、「計量がそろっていない」という位置づけは言い過ぎでした。P·H·P（質量加重なし、並進・回転を射影）の固有ベクトルが必要になるのは、負モードが 2 本以上のときに moddir の番号を決める場合だけです。
  3. **U6-N3（= U5-N1、U5 担当）**: ZTS の端点のずれは、path_hei の種と形の判定にも及びます。
- ライブラリの仕様:
  - NWChem driver（最小化）: 負の固有値のモードの刻みを [0.03, 0.3] × trust に制限する（'Limiting step in negative mode' / 'Forcing step'、`opt_drv.F:1745-1760`）。trust の既定は最小化 0.3、saddle 0.1（`opt_drv.F:1066-1074`）。したがって TS Hessian を渡した QRC を trust 0.1 で実行すると、反応座標の刻みは 0.03 Å になる。
  - NWChem driver（INHESS）: inhess ≠ 1 ならまず drv.hess（restart データ）を調べて優先し、inhess 2 は restart データがなければ freq の Cartesian Hessian を変換する（`opt_drv.F:373-385`）。S7 と継続の挙動はこの仕様と一致する。
  - NWChem string の文書: IMPOSE は初期構造どうしを整列するだけで、最終の画像が整列している保証はない。
  - hfauto の JobStore のキーには deck のパラメータ（trust、moddir、noautoz）が入らない（`jobstore.py:45-52`）。

#### 3.6.3 有用性（◎）

- この単位は「鞍点 → 接続」の証拠を作る唯一の単位で、U7 と U8 の入力になります。v2 で次の経路が初めて実計算で通りました: path_hei の種（HONO、STO-3G）、xTB の不採用 → 種での DFT freq（HONO 15.1 s、STO-3G 1.8 s）、autoz 失敗後の Cartesian 再投入。GS が失敗した HONO でも、W7 と 2.5e-7 Eh で一致する TS に到達し、elementary と判定されました。
- **S7 の効果**: HONO の QRC は 126+17 = 143 → 34+36 = 70 ステップ（311 → 162 s）、maxiter の継続もなくなりました。NH3 は 100 → 90。一方、HCN は 43 → 44 で効果がなく、HONO の速い側は 17 → 36 と遅くなりました。負モードで刻みを制限されたステップは、HONO 30/34・31/36、NH3 38/45・38/45、HCN 7/24・6/20、STO-3G 16/30・6/22 でした。
- **trust 0.3 のプローブ**（v2 と同じ開始構造・同じ TS Hessian）: QRC のステップは HCN 24/20 → 15/12、HONO 34/36 → 14/15、NH3（片側）45 → 24。最終エネルギーの差はすべて ≤ 6e-8 Eh、全側で max(traj) = 開始点。前回のプローブ（12/24/12/15）とも整合します。v2 の QRC 合計は約 266 s → 約 110 s。17 原子（1 勾配約 45 s）なら 1 反応あたり約 30 分の短縮になると見積もられます（外挿による推定）。
- 鞍点も、moddir 1 を明示すると HONO v2 の種から 70 → 18 ステップ（同じ TS、登り返しなし）になりました。

#### 3.6.4 複雑さ（○）

- W3 で drop の式が 1 か所になり（C19）、`tight` が消え、S7 で約 10 行増えました。決定表は純関数のままで、予算も小さい（saddle 2 回、振幅 2 通り、分割の深さ 2）。
- 残る複雑さは、mode_index（質量加重・射影後の順）をそのまま moddir（driver 座標の順）に流用する暗黙の変換と、mode_index ≥ 1 のときだけ noautoz にする分岐です。S20 で 1 つの関数にまとめられ、行数は中立です。

#### 3.6.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U6-I1（前回） | — | S7 は入ったが trust 0.1 のまま。残りは S21 で解消する |
| U6-I2（前回） | low | mode_index をそのまま moddir にする（C17）。継続（TIMEOUT / MAXITER）で moddir が別の構造の更新済み Hessian に再び適用される点も残る |
| U6-I3（前回） | low | periodic_nearest に許容幅がない（C18）。実 run で到達例なし |
| U6-I5（前回） | low | `_register` が connection ゲートより前（設計どおり）。実 run で発動なし |
| U6-I7（前回） | low | 試行を使い切った後の FIND_PATH（= U5-N7） |
| U6-N1 | medium | moddir 0（既定）の初手 uphill が、非停留の種で正の曲率のモードを追う。起きるのは「Cartesian で moddir 0」の組合せ（autoz 失敗後の継続、site の coordinates: cartesian）。M2 で未収束 string の HEI（負モード 2 本、TS より 20.7 kcal/mol 上）が種になるようになり、実計算で初めて表に出た |
| U6-N2 | low（一部訂正） | NWChem string の画像を整列しないので、path_hei の接線と虚モードの重なりを過小評価する（STO-3G 0.118、整列後 0.651）。影響は誤った note、screen のある pipeline での正しい xTB Hessian の誤った棄却（種で DFT freq、17 原子で約 15 分）、負モード 2 本以上でのモード選択の誤りのおそれ |
| U6-N3 | medium | = U5-N1（U5 担当）。端点のずれで single_max が monotonic に化け、行 14 の BARRIERLESS で黙って誤分類しうる経路 |
| U6-N4 | low | S7 を trust 0.1 のまま入れたので一部の側が遅くなった（HONO の速い側 17 → 36、HCN 43 → 44）。全体では改善しているので軽微な退行。validation.md:32 の「127 / 17 = 144」は軌跡の点数から数えたもので、実際のステップは 100 + 26 + 17 = 143 |
| U6-N5 | low | autoz の継続は種から Cartesian でやり直すが、docstring（`engine.py:174-176` の 'A driver job always resumes from its latest frame'）と食い違う。saddle では Hessian が種でしか正しくないので、種からやり直すのが正しい選択。直すのは docstring だけ（HONO v2 で 34 ステップ・48 s を捨てた） |
| U6-N6 | low | policy の qrc_bounds_A の上限と HESSIAN_NEAR_A（0.5 Å）の関係を検査しない。上限を 0.5 より大きくして振幅が 0.5 を超えると、optimize が INPUT_INVALID(hessian_geometry_mismatch) を返し、connection_failed で終わる。既定（0.4）では起きず、失敗はログに残る |
| U6-N7 | low | 分割した子反応は、それぞれ新しく 6 h の予算を持つ（`driver.py:106`、最悪 1 + 2 + 4 = 7 case × 6 h）。仕様の明確化の問題 |

#### 3.6.6 追加の改良案

- **U6-R1 鞍点で追うモードを明示する（C17 の拡張）** — should、修正して採用、中立【S20】
  1. 負モードを選んだら、k = 0 でも常に `moddir k+1` を書く（現状は `if mode_index` なので 0 のとき書かれない）。これで U6-N1 が直る。n_neg == 0 のときだけ従来の moddir 0。
  2. n_neg ≥ 2 のときは `noautoz` にし、k は P·H·P（Cartesian、並進・回転を射影）の負の固有ベクトルの順で決める（C17 の本体。この場合だけ P·H·P を計算する）。n_neg == 1 は autoz のまま moddir 1（負の固有値の数は座標変換で保存される）。Cartesian にする理由は「C17 の番号の整合」として書く（autoz を避ける根拠は HONO の 1 例だけのため）。
  3. TIMEOUT / MAXITER の継続（`engine.py:180-195`、drv.hess と最新フレーム）では、moddir k+1（k ≥ 1）をそのまま再び適用せず、1 にするか外す。
  4. saddle の key_payload に moddir と座標系を加える（deck は JobStore のキーに入らないため。旧キャッシュの結果は TS 検証を通っているので無害だが、比較のためにキーを分ける）。
  - `_first_hessian` の採否の計量は変えない（U6-N2 の実害は接線の未整列で、S19 で解消する。0.3 の閾値の意味を変えない）。
  - 効果: HONO v2 の path_hei の種で、本番 34（AUTOZ 失敗）+ 36（+0.131 Eh の迷走）= 70 ステップ・106 s が、プローブでは 18 ステップ（TS は ≤ 1e-7 Eh で一致）。負モード 1 本の xTB 種でも 1 ステップずつ減る。別の鞍点に迷い込むリスク（reassigned / failed の原因）が減る。
  - コスト: 約 15 行の置き換え、単体テスト 2 件の更新、known_endpoints の 3 反応と find_path_sto3g の再実行（約 15 分、新しい run ディレクトリで）。smoke の `test_moddir_follows_the_second_negative_mode` はそのまま使える。
- **U6-R2 初期 Hessian を渡す optimize は trust 0.3 にする** — should、中立【S21】
  - `input.render_optimize` の `trust 0.1` を `trust 0.3 if init_hessian else 0.1` にする（0.3 は NWChem の最小化の既定、`opt_drv.F:1071`）。Hessian なしの opt（mode-follow、`_register`、init_hessian のない minima）は 0.1 のまま（登り返しの根拠＝対角推定で TS より上へ戻ることが当てはまるのは、こちらだけ）。docstring、design.md（§6 の CONNECT、§7 の minima(dft) と reaction-paths の行）、golden（`test_nwchem_engine.py:95`）を直す。
  - **単位をまたぐ変更**: minima stage の 2 断片以上の xTB 初期 Hessian の opt（`stages/minima.py:104`、TMA·(HF)₂ など本来の対象系）にも効く。採用の条件として、ab_init_hessian（arm A、約 40 分）と known_endpoints の 3 反応を**新しい run ディレクトリで**再実行し（trust は JobStore のキーに入らないので、既存の run では古い deck の結果が黙って再利用される）、所要時間だけでなく basin の数と同定、最終エネルギー、trajectory_above_ts が変わらないことを確かめる。U3 のレビューと validation にも記録する。
  - 効果: QRC のステップがほぼ半分（上記のプローブ）。確認は 6 側・3 反応だけだが、trajectory_above_ts のゲートと振幅の再試行で守られているので許容できる。1 行と文書・テスト、再検証に約 1 時間。
- **U6-R3 FIND_PATH の画像を整列し、チャンクの継ぎ目で端点を戻す** — should、修正して採用（2 つに分ける）
  - (a) `find_path`（`actions.py:481`）を `frames = align_sequential(ctx.frames(last.images))` にする（S2 と同じ形の 1 行）。path_hei の接線と `profile.hei` の補間から剛体回転がなくなる（STO-3G 0.118 → 0.651。S20 と併用すれば 0.989）。【S19】
  - (b) チャンクの継ぎ目で端点を ctx.ends に戻す部分は U5 に回し、U5-R1（M8）として扱う。黙って分類を誤る経路なので、U6 の検証でも could ではなく should 以上とした。
- **U6-R4 文書の修正（autoz の継続、case ごとの予算、QRC の振幅の上限）** — could【C39 に統合】
  - `engine.continuation` の docstring を「autoz の失敗は種から Cartesian でやり直す（saddle の Hessian は種でしか正しくないため）」に直す。design.md に「分割した子反応はそれぞれ 6 h の予算を持つ（最悪 7 case）」「qrc_bounds_A の上限は HESSIAN_NEAR_A（0.5 Å）以下にする」を 1 行ずつ。validation.md:32 の「127 / 17 = 144」を「126 / 17 = 143」に直す。

#### 3.6.7 見送った案

| 案 | 理由 |
|---|---|
| TS Hessian の負の固有値を正に反転して QRC の opt に渡す | 変位点では本当に負の曲率がある（W7 の対角推定でも NH3 は 40/50 ステップが制限された）ので、刻みの制限の根本は trust 0.1 の上限（0.3 × trust）にある。NH3 の反転プローブ（trust 0.1）は rc=1 で 3 ステップで打ち切られ、1 ステップ後の BFGS 更新で負の固有値が戻る兆候を見ただけ。Hessian を加工するコードが増える |
| すべての optimize を trust 0.3 にする | Hessian なし（対角推定）の変位点からでは、HCN の片側が TS より +0.199 Eh 上まで登った（前回のプローブ）。S21 の範囲に限る |
| autoz の失敗後、saddle を最新のフレームから再開する | 種の Hessian は種でしか正しくない。最新フレームから対角推定の Hessian で P-RFO を始めると、追うべきモードを失う。負モード 2 本での autoz は S20 で避けられる |
| 分割した子反応に親の残り時間を引き継がせる | 分割は実計算で一度も起きておらず、Deadline を queue と CaseResult に通す配線が増える。文書（C39）で十分 |
| qrc_bounds_A の上限を読み込み時に 0.5 Å で検査する | 誤設定はログに残る形で失敗し、上書きした例もない。stage の設定が NWChem 固有の定数に依存することになる |
| `_assign` の Registry.find に level_key / composition_id を渡す | reaction-paths に入る DFT 極小は同梱の pipeline では 1 つの level だけ。エネルギー差 5e-5 Eh 以内という条件が別の level を実質的に排除する |
| mode-follow の opt（relax_to_minimum._follow）にも Hessian を渡す | U6 の実計算では `_register` も mode-follow も走っておらず、効果を確かめる材料がない |
| QRC を NWChem の MEPGS（IRC）に置き換える、振幅を大きくする | 前回と同じ理由。S21 で QRC は片側 12〜24 ステップになり、IRC の端点も結局 opt と割付けが必要 |
| 重なりを合否の判定に使う | CH-34 で正しい TS を棄却した実例がある。S20 の後も、重なりは xTB Hessian の採否とモードの選択にだけ使う |
| QRC_AMPLITUDES・_MIN_OVERLAP・_RETRY_A・HESSIAN_NEAR_A を設定できるようにする | 上書きの需要がなく、新しい knob を増やさない方針（X8）に反する |

#### 3.6.8 変更不要な点

- 別ジョブの DFT freq を hfauto で射影して TS を判定すること（is_first_order_saddle の各条件）。
- S7（QRC の両側に TS freq の Hessian を `inhess 2` で渡す。optimize は 0.5 Å の近接を検査し、saddle は指紋の完全一致）。C19。
- 曲率とエネルギー目標 max(3 × drop, 3e-4 Eh) から QRC の振幅を決め、0.05〜0.4 Å に収めること。
- 接続ゲートと割付けのラベル。Registry.find（identity.assign）の厳格な割付けと mapped_equivalent による縮退の判定（NH3 の反転は degenerate と正しく判定された）。
- 1 回目は xTB の Hessian、使えなければ種で DFT freq というフォールバック（v2 で実際に働いた）。autoz 失敗後に種から Cartesian でやり直すこと。
- 決定表を純関数にすること、小さな予算、higher_order のやり直し、collapse / multi_max の中間体での分割、⟨S²⟩ の注記。重なりを合否に使わないこと。

---

### 3.7 U7 熱化学（thermo）

評価: 前回 ○ / ○ / △ → 今回 **◎ / ○ / ○**（有用性の ○ は、ΔG_assoc の幅と None の理由が出ないこと、ΔG_assoc が ReactionRecord のある錯体にしか出ないことによる）

#### 3.7.1 現状仕様の要点（変更点を中心に）

- 対象: tier=dft の全 MinimumRecord（単量体を含む）と、全 ReactionRecord の SaddleClaim.freq_calc。energy_method があれば、同じ start 指紋の SP をエネルギー層に使います。【新】S9: SP がなければ `layer_missing`、engine を呼ばず G=None、notes は [thermo_unavailable, energy_layer_missing]。
- 【新】**C23**: 温度は pipeline の conditions だけから取ります（stage L232）。variant は main（grimme / 100）に qs{grimme, truhlar} × cutoff{50, 100, 150} を加えた 6 個。
- 【新】**M6**: `thermo_frequencies`（`chemistry/thermo.py:29-33`）。極小は全モードを |ν| に、TS は最低モード（反応座標）を除いて残りを |ν| にし、昇順に並べます。engine が送る振動数（`engine.py:114`）と整合ゲートの参照（stage L101、L106-109）の両方でこの関数を使い、GoodVibes には常に `invert=None` を渡します（`worker.py:230`）。`invert_soft_cm1` と soft 限定の分岐は削除。
- 【新】**S15**: スケール因子は `vib_scale` 1 本（`core/method.py:132`、既定 1.0）で、`freq_scale_factor` と `zpe_scale_factor` の両方に渡します（`worker.py:229-230`）。`_SCALE_DB`、`_ALIASES`、`scale_factors`、`resolve_scales`、`zpe_scale`、注記 `scale_factor_from:PBE0/MG3S` は削除。
- 【新】**X12**: JobStore のキーは {freq の job_key, hessian_sha, 規則を適用した後の frequencies_cm1, temperatures_K, settings}（`engine.py:121-122`）。
- 【新】**S14**: 単量体キー (Hill, 電荷, 多重度)、錯体側は (Hill, 電荷) で引き、LOT は `same_pes(state=False)`、G=None の単量体があれば ΔG_assoc=None（L136-176）。反応の blockers は state=True のまま。
- 変わらない部分: thermo_consistent（n_real、ZPE 1e-5 Eh、E 1e-6 Eh）と fail-closed、composite の式、1 atm 基準と 1 bar / 1 M への換算、band は ΔG‡ だけ、population、QH=False、symmno=1 + pymsym の σ、S_rot の事前検査、thread_map、1 subject 1 worker。
- 実物: v2 HCN の TS（`jobs/20/20b351cf…`）の job.json は frequencies_cm1 [2123.04, 2693.56]（−1128.5 は除去済み）、`compute_thermo(..., QS='grimme', s_freq_cutoff=100.0, temperature=298.15, freq_scale_factor=1.0, zpe_scale_factor=1.0, invert=None, symm=True)`、0.277 s。W7 は同じ TS に [−1128.54, 2123.04, 2693.56] と 0.989 / 0.975 を渡していました。
- 文書: design.md:146（QH=False、スケール 1 つ、負モードの規則、σ だけで L は含まない、トンネル補正なし、1 atm、温度は conditions、単量体のキー、state=False、energy_layer_missing、キャッシュのキー）、:46、:150、README:41、validation.md:24（NH3 の ※1）と :39-47（ΔG の変化表）。
- 未実施: C20、C21、C22。

#### 3.7.2 妥当性（◎）

- **M6**: 負モードの約束は、QRC が変位に使うモード（`actions.py:393` の `imaginary_modes[0]`）とも、`SaddleClaim.imag_cm1 = min(freq)`（`actions.py:330`）とも同じ「最低モード＝反応座標」で一貫しています。MinimumRecord.composition_id は必須なので、極小が TS と誤判定されることはありません。本物の GoodVibes 4.3.0 で確かめました: TMA·(HF)₂ の最低モード 51.8 を ±x に置き換えると、HEAD では G(−x) − G(+x) = 0.0000（x = 20 / 10 / 5 / 1）、旧来の扱いでは +0.736 / +0.940 / +1.145 / +1.615（n_real 45 → 44）。
- **S15**: GoodVibes の `calc_bbe`（`thermo.py:703-715`）は、U_vib の中の ZPE を freq_scale_factor で、報告する zpe だけを zpe_scale_factor で計算します。前回は ZPE が 2 つの値に分かれていましたが、1 本でそろいました。
- **W7 → v2 の差は S15 だけで説明できます。** v2 の job.json を vib_scale = 0.989 で再計算すると、W7 の G を 1e-8 Eh 以内で再現しました（構造が変わった HONO の TS −0.003、TMA −0.0025 kcal/mol を除く、`/home/user/hfauto_review2_probe/u7_assess/w7repro.py`）。v2・W7 のどの DFT freq にも反応座標以外の負モードはなく、M6 の実 run での寄与は 0 でした（validation.md:39 の「M6・S15」は「S15」が正しい、U7-N7）。
- **文献との照合**: Kesharwani, Brauer, Martin 2015 の PBE0/def2-SVPD は調和 0.9930、基音 0.9561、ZPVE 0.9848（Table 1 / 3 / 5）、def2-SVP では 0.9915 / 0.9547 / 0.9817。Minnesota DB v5 には PBE0/MG3S しかありません。既定の 1.0 は SVPD の ZPE を約 1.5% 過大にしますが、0.9848 にしたときの差は HCN の ΔG‡ +0.054、HONO +0.022、NH3 +0.013、ΔG_assoc −0.089、HCN の ΔZPE‡ −3.420 → −3.368 で、無視できます。method に依存しない既定として 1.0 は妥当です。
- 整合ゲートの余裕（v2 の 15 subject）: |ΔZPE| ≤ 1.63e-8 Eh、|ΔE| = 0、n_real はすべて一致。
- **残る弱点**（すべて low で潜在的）:
  1. energy_method の層が欠けた subject でも、ΔE‡ / ΔE_rxn を freq の E との混成で出す（U7-N1 = U0-N1）。
  2. 錯体 → 単量体の対応が (Hill, 電荷) だけで決まり、同じ式・電荷の組成が 2 つあると後勝ちになる（U7-N3）。
  3. 単量体のアンサンブル（L166）は同じ (Hill, q, m) の全 dft 極小を集めるので、配座だけでなく化学的な異性体（HCN に対する HNC、単量体の explore の生成物など）も入る。S14 の fail-closed により、エネルギーの高い異性体 1 つの thermo の失敗だけで ΔG_assoc 全体が None になる（U7-N9、検証で追加）。
  4. M6・S14・S9 は実計算で一度も通っておらず、fake と GoodVibes のプローブでしか確かめていない。

#### 3.7.3 有用性（○）

- 順位付けに必要な ΔG‡・ΔG_rxn・ΔZPE‡・band・blockers は、実 run でそろって出ています。v2 は 15 subject で失敗 0、1 subject 0.23〜0.28 s で、コストは無視できます。
- 錯体系の主な出力である ΔG_assoc（v2 の TMA·(HF)₂ では唯一の熱化学の値、−11.53）には課題が残ります。
  - 6 variant の幅（298 K で −12.089〜−10.934 の幅 1.155、273 K で 1.053、373 K で 1.457）が出力されない（C21）。ΔG‡ の band（≤ 0.0006）とは桁が違う。
  - S14 の fail-closed で None が出る場面が増えたのに、理由が出ない。
  - `Evidence.n_external` が `Literal[5, 6]`（`core/evidence.py:85`）なので、単原子の単量体（F⁻）は dft 極小を持てず、F⁻·(HF)ₙ（W1 で CREST の smoke に加わった組成）の ΔG_assoc は構造的に常に None。
  - （検証で追加）ΔG_assoc は `_reaction`（`stages/thermochemistry.py:207`）からしか計算されず、ReactionThermo にしか載らない。ReactionRecord を 1 件も持たない錯体の組成（discover で仮説が 1 つも残らなかった場合など）では、極小の SpeciesThermo は出ても ΔG_assoc はどこにも出ない。v2 の TMA·(HF)₂ で出たのは、宣言反応 neutral_to_shared_proton（same_basin）があったからにすぎない。
- 点群は GoodVibes 内部の `bbe.point_group` では正しく決まっていますが（v2 の全 subject: C∞v、Cs、C1、C3v、D3h、C2v）、`ThermoResult.point_group` は空（`api.py:296`）で、どこにも記録されません。population・H・S_rot は書くだけのままです。
- 標準状態: 1 M の ΔG_assoc は Δn = −2 の shift −3.789 で −15.315（W7 は約 −15.38）。1 bar の shift は +0.0156。

#### 3.7.4 複雑さ（○、前回 △）

- `9b61e8d` で thermo の 6 ファイルが +89 / −132 行（chemistry/thermo.py は +7 / −40 で 131 行）。stage 253 行、engine 129 行、worker 108 行、合計 621 行。設定キーは 3 つ減り（temperatures_K、zpe_scale、invert_soft_cm1）、増えたのは純関数 1 つと same_pes の state 引数だけで、前回の見込みどおりでした。thermo_frequencies を engine と stage で 2 回呼ぶのは、同じ純関数なので問題ありません。
- 小さな残り: `ThermoResult.notes` は常に `()`（`engine.py:68`）で出力の経路として死んでいる、`ThermoResult.S_rot` は worker の事前検査の値を返すだけで誰も読まない、population は書くだけ。前の 2 つは C20 に使い回せば構造を増やさずに片づきます。

#### 3.7.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U7-I2（前回） | low | 会合量の理由が出ない、単原子の単量体は常に None（下記 U7-N2） |
| U7-I5（前回） | low | ΔG_assoc の感度幅がない（C21） |
| U7-I7（前回） | low | 点群・σ が記録されない（C20） |
| U7-I8（前回） | low | 基準が反応物の極小そのもの（C22 は方針どおり見送り。HONO で約 +0.37） |
| U7-I9（前回） | low | population・H・S_rot は書くだけ |
| U7-I11（前回） | low | 単原子種（方針どおり対応しない） |
| U7-N1 | low | = U0-N1。層が欠けた subject でも ΔE を freq の E との混成で出す（fake で再現: TS の SP だけを外すと dE_act=6.275、dG_act=None、blockers に thermo_unavailable と mixed_level_of_theory）。html.py:192, 209 はこの ΔE‡ を表と図に出す |
| U7-N2 | low | `_association`（L155-176）は、組成が宣言されていない・単量体の極小がない・単量体の G が None・LOT が合わない、のどの場合も無言で (None, None) を返す。ReactionThermo に notes がない。単原子の単量体では常に None で、この制約は design.md にも書かれていない |
| U7-N3 | low | 錯体 → 単量体の対応が (Hill, 電荷) だけで、`_monomers`（L144-146）の後勝ちになる。system の検査は同じ式・同じ総電荷の組成を拒否しない。実害が出るのは多原子イオン対を中性の組と並べて宣言した場合（例の {NH4⁺, F⁻} は F⁻ が単原子なのでどのみち None） |
| U7-N4 | low | Conditions と ThermoSettings に値の範囲の検査がない。温度が空なら全 stage の後に `report.py:48` の `temperatures_K[0]` で IndexError、T ≤ 0 なら DFT の後で全 G が None、vib_scale ≤ 0 も通る。重複した温度は manifest の merge で 1 つにまとまるので害はない。cutoff_cm1 = 0 は RRHO として正当 |
| U7-N5 | low | vib_scale = 1.0 の根拠と文献値が書かれておらず、使った値も run の記録に残らない（v2 HCN の resolved_config.yaml の thermo は `engine: goodvibes` だけ。値は jobs/*/job.json と settings_sha にしかない） |
| U7-N6 | low | `ThermoResult.notes` は常に `()`、S_rot は書くだけ（死んでいる経路） |
| U7-N7 | low | validation.md:39 は ΔG の変化を「M6・S15」によるとしているが、実際は S15 だけ |
| U7-N8 | low〜medium | （検証で追加）ΔG_assoc は ReactionRecord のある錯体にしか出ない。非共有結合錯体を目的とするなら、「錯体の組成ごとに 1 件の会合レコード」にするかを判断するか、少なくとも design.md に制約として書くべき |
| U7-N9 | low | （検証で追加）単量体アンサンブルに化学的な異性体も入り、1 つの thermo 失敗で ΔG_assoc 全体が None になる（異性体が高エネルギーなら値への影響は小さい） |

#### 3.7.6 追加の改良案

- **U7-R1 層が欠けた subject では ΔE も出さない** — could【C28、U0-R4 と同じ】
- **U7-R2 ΔG_assoc=None の理由を 1 つ記録し、組成の曖昧な対応を fail-closed にする** — could、修正して採用、【複雑さ増】【C29】
  - `_association` が理由を返し、ReactionThermo に `notes: tuple[str, ...] = ()` を 1 つ足して `assoc_unavailable:<reason>` を入れる。reason は `no_composition` / `monomer_missing`（単原子の単量体を含む）/ `monomer_thermo_unavailable` / `mixed_level` に絞る。blockers には入れない（順位に使うものではない）。html の thermo 表の blockers 列に並べて表示する。
  - 同じ (Hill, 電荷) に別の部品の組が来たら、`_monomers` で重複キーを集合に記録してそのキーを None にするだけにとどめ、system 側では拒否しない（中性の組とイオン対の組を並べるのは正当な使い方）。
  - あわせて、U7-N8 と U7-N9 を design.md に制約として書く（反応のない錯体では ΔG_assoc は出ない、単量体のアンサンブルには同じ (Hill, q, m) の異性体も入る）。「組成ごとに 1 件の会合レコード」にするかは、この改良と C21 を入れるときに判断する。
  - 効果: F⁻·(HF)ₙ のように単原子の単量体があって出ない場合と、単量体の熱化学が失敗した場合を、利用者が区別できる。同じ式の組成を 2 つ宣言したときに誤った基準の ΔG_assoc を黙って出さなくなる。約 10 行とテスト 1 件。
- **U7-R3 Conditions と ThermoSettings に宣言的な値の制約を付ける** — could【C30 に統合】
  - Conditions: temperatures_K を PositiveFloat の tuple・min_length=1・重複禁止、standard_states も min_length=1・重複禁止。ThermoSettings: vib_scale を PositiveFloat、cutoff_cm1 を NonNegativeFloat（0 は許す）。拒否のテスト 1 件。
- **U7-R4 点群を、今は空の ThermoResult.notes で記録し、S_rot を protocol から外す（C20 の最小版）** — could、中立【C20】
  - worker は main の settings の結果に `"point_group": r.bbe.point_group` を入れ、S_rot は返さない（事前検査は残す）。`engine.parse` は `notes=(f"point_group:{pg}",)`。ThermoResult の S_rot フィールドを削除。stage は `r.notes` を SpeciesThermo.notes にすでに合流させているので変更不要。notes の文字列は `_RANK_BLOCKERS` と照合されないので無害。protocol の形が変わるので古い結果は検証で miss になるが、再計算は 0.3 s。約 4 行とアダプタのテストの修正。
  - 効果: σ を監査できる（pymsym の閾値 1e-3 で対称性が崩れて C1 に落ちたとき、NH3 のような σ = 3 の分子で 0.65 kcal/mol ずれるのを事後に検出できる）。死んでいる経路が 2 つ片づく。
- **U7-R5 スケールの根拠と単原子の制約を文書に書き、validation の帰属を直す** — could【C39 に統合】
  - design.md:146 に: 「vib_scale = 1.0 は未スケールの調和振動数。PBE0/def2-SVPD の文献の ZPE 因子は 0.985（Kesharwani 2015）で、切り替えても ΔG‡ は ≤ 0.06、ΔG_assoc は ≤ 0.09 kcal/mol しか変わらない」「単原子の単量体（F⁻ など）を含む組成、反応を持たない錯体では ΔG_assoc は出ない」。validation.md:39 の「M6・S15」を「S15（M6 の寄与は 0。負の副次モードがないため）」に直す。
- **U7-R6 ΔG_assoc にだけ感度の幅を付ける（C21）** — could、【複雑さ増】【C21】
  - `_reaction` で既存の 6 variant の tables に対して `_association` を回し、`assoc_band_kcal: tuple[float, float] | None` を 1 フィールド出す。shift を足し、表示だけで順位には使わない。vs_separated の幅は付けない。約 8 行とテストの assert 1 行。
  - 効果: 弱い錯体の判断に効く不確かさ（298 K で幅 1.155、373 K で 1.457）が付く。

#### 3.7.7 見送った案

| 案 | 理由 |
|---|---|
| 同梱の pipeline の thermo に `vib_scale: 0.985` を明示する | 効果は ΔG‡ ≤ 0.054、ΔG_assoc 0.089 で手法誤差よりずっと小さい。validation の値との連続性が崩れ、method を替えたときに古い因子が残る危険がある。文書に 1 文書けば足りる |
| スケール因子を MethodSpec に移す | method の層に knob が増え、層の結合を生む。利得は 0.1 kcal/mol 未満 |
| vib_scale を全 SpeciesThermo の notes に記録する | 全 subject で同じ値の雑音になる。既定値を残すなら resolved_config.yaml に既定値込みで dump する（U8 の担当）のが筋 |
| 単原子種を最後まで扱えるようにする（Evidence の n_external=3、SP だけの経路、worker の分岐） | 3 つの単位にまたがる変更で、効くのは F⁻ 系の ΔG_assoc だけ。F⁻·(HF)ₙ を正式な対象にするまでは C29 と文書で足りる |
| M6・S14・S9 を実 DFT の run で通す専用の計算 | M6 は純関数で、GoodVibes の実物のプローブで効果を確認済み。17 原子の freq は 1 本 15 分かかる |
| population・H を削除する | 合わせて約 11 行で負担がない。配座の分布は科学的な情報として残す価値があり、C22 の前提にもなる |
| C22（組成アンサンブル基準の有効障壁）を今回入れる | 新しい根拠がない。効くのは HONO だけで約 +0.37 |
| variant ごとに thermo_consistent を判定する、QH=True・トンネル補正・hindered rotor を入れる | 前回の見送り理由が変わっていない |
| torsional TS で最低モード ≠ 反応座標になりうる場合に、overlap で反応座標のモードを選び直す | QRC と SaddleClaim も同じ約束に立っており、thermo だけ変えると不整合になる。G への影響は qRRHO で 0.3 kcal/mol 未満。直すなら C17（S20）と一緒に paths 側で行う |

#### 3.7.8 変更不要な点

- thermo_frequencies の純関数と固定規則、engine と整合ゲートの両方で使うこと、GoodVibes に invert=None で渡すこと。
- vib_scale 1 本を freq_scale_factor と zpe_scale_factor の両方に渡すこと。
- JobStore のキー（規則を適用した後の振動数と温度を含む）。温度は conditions だけ。ThermoSettings の extra='forbid'。
- thermo_consistent と fail-closed。S9 の energy_layer_missing。S14 の単量体キー・state=False・fail-closed。反応の blockers は state=True。
- variant の格子、1 atm 基準と換算、composite の式、S_rot の事前検査、version pin、1 subject 1 worker。自前の振動数・質量・回転定数を渡し、存在しないファイル名で再パースを防ぐこと。pymsym の σ。
- design.md:146 の熱化学の規約と、validation.md:24 の NH3 の ※1。

---

### 3.8 U8 一点計算・手法パネル・順位付けと報告・実行基盤（sp・report・execution・configs）

評価: 前回 ○ / ○ / ○ → 今回 **◎ / ○ / ○**

#### 3.8.1 現状仕様の要点（変更点を中心に）

- **sp**（`stages/single_point.py`、78 行、変更なし）: 対象は反応に参加する tier=dft の極小と各反応の saddle.freq_calc（既定 targets=reaction_stationary_points、all_minima なら全 dft 極小も）。構造は freq Evidence.final、電荷・多重度は freq.level。method × subject を直列、deadline なし。
- **NWChemEngine.energy**: 【新】S17: `kind == 'wft'` かつ multiplicity > 1 はジョブを作らず INPUT_INVALID `wft_closed_shell_only`（`engine.py:317-318`、MP2 も対象）。`render_wft` から nopen / uhf を削除し、scf ブロックは restart 時の `vectors input` だけ（`input.py:157-167`）。CCSD は maxiter と thresh を書かない（NWChem 既定 20 / 1e-6）。
- **method_panel**: 【新】S8 の既定 `[pbe0-d3bj_def2-tzvpd, wb97x-d3_def2-tzvpd]` + 参照行（SVPD の opt / saddle / freq / paths の SP）。CCSD(T)/def2-TZVPD は opt-in。【新】M5 の ±1 の不感帯（`summary.py:32, 211-212`）。
- **report**（`stages/report.py`、68 行、変更なし）: rankable（outcome ∈ {elementary_step, degenerate_rearrangement, reassigned_step}、dG_act あり、blocker 4 種なし、|ΔZPE‡| ≤ 1 + 4·max(1, n_H)。`gates.py:297-315`）、band が推移的に重なれば同順位（1,1,3）、coverage、ranking.csv・coverage.csv・method_panel.csv・report.html（3Dmol.js 2.5.5）。T は conditions の先頭（`report.py:47-48`）。
- **実行基盤**（execution/* 744 行、runner・preflight・layout・cli は変更なし）: JobStore のキー = sha256({engine, version_pin, kind, key_payload})、ExecutionSpec を含まない。終端的失敗 4 種だけ記録・再利用、ladder 4 種、energy ジョブの TIMEOUT は継続しない、コア数セマフォ、SiteLock、`/mnt` 上の scratch の拒否。
- 【新】S10（読み込み時の停留点 method の一意性）、S9（composite の fail-closed）、C24（文書）。
- 周辺の変更: optimize のキーから tight が消えた（W7 の opt は再利用されない）。energy ジョブのキーは不変（W7 以前の SP は再利用できる）。string の monitor を削除。thermo のキーに振動数と温度。
- 設定の量: StageConfig の knob 59 → 51（conformers 11 → 6、thermo 10 → 7）、ThermoSettings 8 → 5、ConformerSettings 6 → 5、method ファイル 6 → 5、configs の YAML 185 → 179 行、hfauto 全体 10,707 → 10,651 行。
- 未実施: C2、C25、C26、C27。

#### 3.8.2 妥当性（◎）

- **sp**: freq.final の固定構造に、freq.level の電荷・多重度で一点計算し、Level と構造の echo（1e-4 Å）を照合する設計は今も正しいです。S1 の後も MinimumRecord の opt_calc と freq_calc は同じ outcome から作られ、代表は差し替えられない（`minimum.py:241-275`）ので、パネルが照合する opt_calc.final の指紋と SP の構造はずれません。
- **M5**: HONO のプローブ（`/home/user/hfauto_review2_probe/u8/hono`）で ΔE_rxn は SVPD 0.132、PBE0/TZVPD 0.236、ωB97X-D3/TZVPD 0.325 とすべて正（実験値 +0.30 と符号が一致）で、sign_disagreement は False、panel_report/ranking.csv と report/ranking.csv はどちらも rank 1・12.2237 でした。不感帯の分岐（±1 以内の負値）は実データでは通っておらず、単体テストだけが確かめています。|ΔE_rxn| < 1 の符号の反転は ΔE‡ による順位付けの結論に影響しないので、判定の定義は妥当です。
- **S8**: 3 行（参照 → 基底の感度 → 汎関数の感度）の構成で、陰イオン・H 結合系でも「誤った方向の参照」はなくなりました。HCN の幅は ΔE_rxn 13.15〜13.82、ΔE‡ 46.37〜46.62。opt-in の CCSD(T)/TZVPD は 15.22 / 47.81（TZVP では 15.43 / 48.08）。
- **S17**: NWChem の CCSD は RHF 専用です（公式文書）。プローブ（`u8/s17`）では (電荷 1, 多重度 2) と (0, 3) の両方がジョブ 0 件で拒否されました。v2 の smoke の MP2/def2-SVP（HCN）は閉殻、Frozen core 2、wall 0.2 s。
- **CCSD の反復**: HCN/TZVPD は 11〜13 反復（TZVP の 10〜13 と同じ）、HONO/TZVP は最大 18 反復で、20 までの余裕は小さいままです。HONO TS/TZVPD は、負荷下（load 約 11、2 rank）で SCF と 4 添字変換に約 7 分、CCSD の 1 反復目に 95.7 s かかり、制限時間内に反復数を確認できませんでした。未収束は NONZERO_EXIT で artifact が failed になるだけで、誤った値が黙って入ることはありません。
- **順位付け**: 変更なし。小分子の band 幅は 0.0004〜0.0006 kcal/mol で、同順位は実質的に起きません。
- **composite**: S9 は正しい。design.md:150 の手順は動きません（U8-N1 = U0-N2）。
- **実行基盤**: energy ジョブのキーは変わっていません。S10 は、minima と reaction-paths の method が必須なので None との偽の食い違いは起きません。
- ライブラリの仕様:
  - NWChem 7.2.3 の `ccsd_input.F`: ccsd ブロックは maxiter（`ccsd:maxiter`）、thresh、freeze、diisbas などを受け付ける。
  - CCSD は permanent_dir の `./job.t2` を再開用に読みにいき、なければ MP2 の初期推定に戻る。attempt ディレクトリは毎回新しいので古い t2 を読む危険はない（HCN/TZVPD で job.t2 と job.db は約 0.7 MB）。
  - 3Dmol.js 2.5.5 は jsDelivr の latest タグと一致し、build/3Dmol-min.js は HTTP 200（537,792 B）。

#### 3.8.3 有用性（○）

- **費用に見合い、有用なもの**: 既定パネル（HCN 6 SP・36 s。前回は 12 SP・63.7 s）、追記型の method_panel と「panel_report が最終」の明記、rankable と ranking.csv、coverage.csv（tma_hf2 の陰性理由の内訳が出る）、JobStore・ladder・SiteLock・セマフォ（v2 で TIMEOUT も終端的失敗も 0 件）、CCSD(T) の 1 語の opt-in（DFT の 6 点は再利用、hits 6 / misses 3）。
- **有用性を下げているもの**:
  1. 主対象の現在の結果（v2 の tma_hf2_discover と water_same_basin）は、どちらも same_basin の 1 行だけで、その行に `outcome:same_basin;thermo_unavailable;dzpe_out_of_tolerance` が付きます。report.html では thermo_unavailable が ΔG_assoc −11.53 の隣に出ます。熱化学は成功しているので、読者はありもしない失敗を探すことになります（U8-N3）。
  2. ranking.csv に T・標準状態・ΔG_rxn の列がない（U8-7）。さらに metric=dG_rxn のとき、rank_rows は dG_rxn で並べるのに CSV の値の列は dG_act と dG_act の band だけで、並べた値が表に出ません。rankable は metric に関係なく dG_act と dzpe を要求するので、dG_rxn で並べても TS のない反応は入りません（U8-N5、検証で追加）。
  3. method_panel は v2 の検証でも主対象でも実行されておらず、17 原子の ωB97X-D3/TZVPD（約 320 関数と見積もり）の費用は実測していません（1 点十数分以上と推定）。
  4. JobRunner が途中で回収したジョブ失敗（v2 の HONO の GS の nonzero_exit）は `hfauto status` にだけ出ます（設計どおりで問題ない）。
- method_panel.yaml:3-4 の「HCN / HONO は 1 点数秒〜数分」は TZVP での実測値です。TZVPD では関数が HCN 68 → 83、HONO 99 → 126 に増え、HCN 約 20 s、HONO 数分と見積もられ、表現の範囲に収まります（負荷時の実測値は約 10 倍に膨らむので参考にならない）。

#### 3.8.4 複雑さ（○）

- 減ったもの: knob、method ファイル、YAML の行数、UHF の描画（input.py 175 → 173 行）と ccsd(t) 専用の分岐（判定は `kind == 'wft'` の 1 つに一般化）。summary.py の変更は定数 1 つと条件式だけ（290 → 293 行）。
- **新たに見つかった削減の余地**: M2 で NWChem string の monitor が消えた結果、6 つのアダプタ（nwchem、xtb、crest、pysis、readuct、goodvibes）の `monitor()` がすべて None を返します。そのため次のコードには実行時に到達しません（約 30 行）。
  - `Adapter.monitor`（`jobs.py:53, 150`）、`run_command` の monitor / poll_s 引数と監視ループ（`process.py:62-137`）、`CommandResult.stopped`
  - stopped → STAGNATED の 4 分岐（`xtb.py:56-57`、`crest.py:133-134`、`pysis/engine.py:147-148`、`nwchem/output.py:229-230`）
  - `FailureKind.STAGNATED`（`evidence.py:113`、design.md:33 の「11 種」）
- 残るもの: `site.memory_mb` は依然として未使用（`config.py:36` で読むだけ）なのに、`wsl_local.yaml:3-4` のコメントはスケールアップ時に memory_mb を書き換えるよう勧めています。報告物は 3 つ（`<run>/report/`、`<run>/panel_report/`、`hfauto report` の `<run>/report.html`）のままですが、内容の矛盾はなくなったので残してよい。conditions の 3 pipeline への重複も、構造変更になるので扱いません（前回の判断どおり）。

#### 3.8.5 残る課題・新たな課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U8-3（前回） | low | CCSD の maxiter 20（= U0-I8、C2） |
| U8-5（前回） | low | composite の手順と thermo の上書き（一部解消。U8-N1） |
| U8-7（前回） | low | ranking.csv に T・標準状態・ΔG_rxn の列がない（C25）。U8-N5 も加わった |
| U8-8（前回） | low | TS のない outcome に dzpe の blocker が付く（C26）。範囲は thermo_unavailable にも及ぶ（U8-N3） |
| U8-9（前回） | low | site.memory_mb が未使用（C27）。wsl_local.yaml:3-4 のコメントは効かない memory_mb の書き換えを勧めている |
| U8-10（前回） | low | method_panel.yaml が独自の conditions を持つ（受容済み。参照行は report.py:47-48） |
| U8-N1 | low | = U0-N2。composite の手順（design.md:150）はそのとおりに実行すると空の結果になる（`layout.py:176-188` の view は同じ run の done stage の union だけ、JobStore も run ごと）。会合量には targets: all_minima が要る |
| U8-N2 | low | M2 の後、monitor の仕組み全体と FailureKind.STAGNATED が使われないコードになった（動作上の害はなく、複雑さだけの問題） |
| U8-N3 | low | TS のない outcome には thermo_unavailable も誤って付き、主対象の報告で ΔG_assoc と矛盾して見える（`gates.py:304-314` は outcome が rankable でなくても dG_act / dzpe が None なら 2 つを付ける。ReactionThermo.blockers には same_basin の場合に thermo_unavailable が出ないので、この表示は rankable ゲートだけが作っている）。回帰ではなく、前回からの問題の範囲が実データで確かめられたもの |
| U8-N4 | low | 手法パネルは正式な検証でも主対象でも実行されておらず、17 原子での費用も未実測 |
| U8-N5 | low | （検証で追加）report.metric=dG_rxn のとき、ranking.csv に並べた値が出ない。rankable は metric に関係なく dG_act と dzpe を要求する |

#### 3.8.6 追加の改良案

- **U8-R1 rankable: 順位を付けられない outcome には、outcome の理由と本当の blocker だけを返す（C26 の修正版）** — should、修正して採用、中立【S22】
  - `chemistry/gates.py` の rankable で、outcome が `_RANKABLE_OUTCOMES` にないときは `[f"outcome:{o}", *[b for b in _RANK_BLOCKERS if b in flags]]`（flags は thermo.blockers と participant_notes）を返す。dG_act が None であることから作る thermo_unavailable と、dzpe の検査だけを省く。thermo 側の `_blockers` が出す本当の thermo_unavailable（参加者の G が欠けた場合）、spin_contaminated、mixed_level_of_theory、method_sign_disagreement は残る（提案のままの早期 return だと、これらも要約行から消えてしまうので修正した）。
  - rankable な outcome では、U0-R5 のとおり dzpe が None なら dzpe_out_of_tolerance を付けない（`if dzpe is not None and abs(dzpe) > tolerance:`。dzpe が None なら dG_act も None で thermo_unavailable がすでに付く）。
  - test_gates.py に same_basin の 1 行を足す。3〜4 行。
  - 効果: 主対象の報告行（v2 の tma_hf2・水）が `outcome:same_basin` だけになり、report.html で thermo_unavailable が ΔG_assoc の隣に出る矛盾がなくなる。hcn_missing 型（thermo 欠落）の行は thermo_unavailable だけになる。blockers 列が「順位を付けられない本当の理由」だけを示す。
- **U8-R2 CCSD ブロックに maxiter 50 を書く（C2）** — could（U0 の検証は should）【C2、U0-R2 と同じ】
- **U8-R3 ranking.csv に T_K・standard_state・dG_rxn_kcal を加える（C25）** — could、中立【C25】
  - RankRow に 3 つのフィールドを足し、rank_rows で選んだ (T, state) と `t.dG_rxn_kcal` を入れる。`_RANK_HEADER` と `_rank_cells` に列を足す。ranking.csv が単独で読めるようになり、metric=dG_rxn のとき並べた値も表に出る。数行とテスト 1 件。
- **U8-R4 site.memory_mb を上限の検査に使う（C27）** — could、中立【C27】
  - preflight で、NWChem 系の EngineSite ごとに `execution.ranks * execution.memory_mb_per_rank <= site.memory_mb` を検査し、超えたら計算の前に拒否する。係数は新しく作らない。wsl_local.yaml のコメントが memory_mb の書き換えを勧めているので、削除ではなく検査を選ぶ。3〜5 行とテスト 1 件。
- **U8-R5 使われなくなった monitor 経路と FailureKind.STAGNATED を削除する** — could、【複雑さ減】【C37】
  - `Adapter.monitor`、`run_command` の monitor / poll_s と監視ループ（`process.wait(timeout=remaining)` の 1 回にする）、`CommandResult.stopped`、6 アダプタの `monitor()`、4 つの stopped → STAGNATED の分岐、`evidence.py:113` の STAGNATED を削除する。テスト（`test_process.py:33-45`、`test_nwchem_output.py:92-97`、`test_cli_commands.py:37-40`）を直し、design.md:33 を「10 種」にする。STAGNATED は TERMINAL ではないので result.json に残っていることはなく、旧 run の run_state の failures_by_kind は文字列の dict で読み込みに影響しない。約 30 行の削除。
- **U8-R6 composite の手順と method_panel の実計算確認を文書に記録する** — could
  - (a) design.md:150 の手順の修正は U0-R1(a) と同じ【S23】。(b) validation.md への記録は U0-R3 と同じ（既存のプローブの結果を書けば足り、再実行は不要）【C39】。

#### 3.8.7 見送った案

| 案 | 理由 |
|---|---|
| energy ジョブの TIMEOUT を記録する（TERMINAL 化）、継続する | 一時的な失敗なので TERMINAL の定義に合わない。事故の原因だった既定の CCSD(T) は S8 で opt-in になった |
| CCSD(T) に系の大きさ（関数の数）の上限を設ける | 閾値と見積もりのコードが要る。opt-in とコメントの目安で足りる |
| JobRunner が途中で回収したジョブ失敗を coverage.csv に加える | `hfauto status` と run_state.json に出ている。coverage は発見の被覆率の表で、重複になる |
| thermo の artifact id に stage id を入れる | report が ReactionThermo を選ぶ箇所に stage の選択が要り、参照箇所も広く変わる。run の複製の文書化（S23）の方が安い |
| 符号ゲートを ΔE‡ にも広げる、手法間の幅に閾値を設けて順位から外す | 新しい閾値が増える。幅は method_panel.csv に出ており、利用者が判断できる |
| method_panel の幅を順位の band に加える | ほぼすべてが同順位になる（前回の判断を維持） |
| 追記用の composite.yaml を同梱する（C1） | 同じ run には置けず、run の複製と stage id の付け替えが必要。同梱しても手順は減らない |
| 開殻 WFT（UHF-MP2 のパーサ修正、TCE による CCSD(T)） | ⟨S²⟩ を検査できない値がパネルに入る。対象の化学は閉殻が中心（X9 の方針） |
| CCSD の job.t2 / job.db を scratch に移す | HCN で約 0.7 MB、HONO でも約 10 MB で、CCSD(T) は小分子の opt-in |
| 報告物を 1 つにまとめる | 前回の見送り。HONO での矛盾はなくなり、複数条件の run を読み直す用途も残る |
| 主対象（TMA·(HF)₂）で既定パネルの時間を今回実測する | 同時実行のプローブで load 8〜12 あり、壁時計時間に意味がない。17 原子の ωB97X-D3/TZVPD は 10 分の上限を超える見込み |

#### 3.8.8 変更不要な点

- sp の構造を freq Evidence.final に固定し、電荷と多重度を freq.level から取り、同じ構造の subject の parents を連結すること。
- Level の観測と照合、構造 echo の 1e-4 Å 検査、WFT の開殻をジョブの前に拒否すること。
- パネルのキーを Level.full_key にすること、±1 の不感帯をモジュール定数にすること、追記型の method_panel と「panel_report が最終」の規則。既定パネルと CCSD(T) の 1 語の opt-in。
- `freeze atomic`、NWChem 既定の thresh 1e-6、restart 時だけの scf ブロック、sp の直列実行。
- rankable を順位の唯一の関門とし、1,1,3 の同順位と coverage.csv を使うこと（AR-27）。energy_method の fail-closed と、composite を既定の流れに入れないこと。
- JobStore のキーに ExecutionSpec を含めず版数の pin を含めること、終端的失敗の記録と `--retry-failed`、4 種の ladder、energy ジョブの SCF 再投入。SiteLock、コア数セマフォ、`/mnt` 上の scratch の拒否、timeout 14,400 s。
- S10 の読み込み時検査と、`gates:` の未知キーの拒否。自己完結の report.html。

---

## 4. 実計算による確認のまとめ

出典は docs/validation.md v2（WSL2 Ubuntu、4 vCPU / 11 GB、NWChem 4 rank × 1,200 MB）と、各 run の `run_state.json`・`command_result.json`・manifest の再集計です。所要時間は `/usr/bin/time` の壁時計時間、ジョブ数は `hfauto status` の misses + hits（括弧内は再利用）。ΔG は 298.15 K・1 atm。

### 4.1 run 全体（W7 → v2）

| 系 | pipeline | 所要時間 W7 → v2 | ジョブ数（再利用）W7 → v2 | 主な差と原因 |
|---|---|---|---|---|
| HCN→HNC | known_endpoints | 1:40 → 1:41 | 29（2）→ 29（2） | 変化なし。QRC は S7 の効果なし（43 → 44 ステップ） |
| HONO trans→cis | known_endpoints | 7:22 → 9:32 | 29（2）→ 17（0） | SCREEN の GS（DLC）が発散して失敗 → FIND_PATH（200 s）と、path_hei の種からの saddle（AUTOZ 失敗を含め 106 s）が加わった。QRC は 311 → 162 s（S7） |
| NH3 反転 | known_endpoints | 2:08 → 1:55 | 28（1）→ 27（1） | QRC 100 → 90 ステップ（S7）、dft の freq 2 → 1 本（S1） |
| 水（同一 basin） | known_endpoints | 0:15 → 0:11 | 5 → 4 | dft の freq 2 → 1 本（S1） |
| TMA·(HF)₂ | discover | 52:33 → 31:15（−41%） | 72（3）→ 60（2） | dft stage 3,050 → 1,771 s（M3 で偽の AFIR 生成物が消えた）、screen の xTB freq 8 → 4 ジョブ（S1）、explore 37 → 31 attempt |
| amine パイロット | discover `--to explore` | 1:56 → 1:54 | 65（4）→ 56（2） | attempt 44 → 40（M3 で TMA 単量体の trial 14 → 11）、低レベル TS 19 → 22 件 |
| A/B arm A | minima(dft)、`include: all` | 55:18 → 40:34（−27%） | — | S1 で AFIR 生成物の構造が既知の basin に入り、17 原子の freq（約 890 s）を省いた |
| FIND_PATH（新規） | screen なし、PBE0/STO-3G | — → 1:03 | — | IDPP → string 2 チャンク → single_max → 鞍点 −1223.9i、elementary_step |
| (HF)₃ の配座（新規） | discover `--to conformers` | — → 0:34 | — | CREST 2/2 成功（M4） |
| 実エンジンの smoke | pytest -m real | 59 s → 72 s | 10 → 12 件 | 新しい 2 件（STO-3G string 14 s、F⁻·(HF)₂ の `--noopt` 1.4 s） |

### 4.2 処理（stage）ごとの時間とジョブ数

| 処理 | 系 | W7 | v2 | 変化の主因 |
|---|---|---|---|---|
| conformers | TMA·(HF)₂ | stage 17 s。CREST 2 本（組成 10.86 s、2 配座、初期最適化あり） | stage 19 s。CREST 2 本（組成 12.05 s、3 配座、初期最適化なし） | M4。配座数の増加は CREST の実行ごとの揺らぎと切り分けられない |
| conformers | amine パイロット | 12 s、species 6 | 13 s、species 5 | TMA 単量体の CREST が 2 配座 → 1 配座（揺らぎ） |
| conformers | (HF)₃ | —（前回のプローブ: 87〜236 s、NH3·(HF)₃ は失敗） | 32 s（壁時計 34 s）。CREST 4.84 s / 19.17 s、2/2 成功、species 7 / 9 | M4 |
| minima(screen) | TMA·(HF)₂ | 0 s。xTB opt 8・freq 8 ジョブ、known 0 | 2 s。xTB opt 9・freq 4 ジョブ、known 6 | S1（直列化とその場の登録）。species の増加は S13 |
| minima(screen) | amine パイロット | xTB opt 7・freq 7、known 0 | opt 6・freq 5、known 3 | S1 |
| explore | TMA·(HF)₂ | 83 s。37 attempt（NT2 20・AFIR 17）、生成物 1、source の緩和 −4.06 kJ/mol・16 反復 | 81 s。31 attempt（20・11）、生成物 0、source の緩和 0.00・2 反復 | M3、S6 |
| explore | amine パイロット | 102 s、44 attempt | 98 s、40 attempt | M3 |
| minima(dft) | TMA·(HF)₂ | 3,050 s。opt 4 本 907 s、freq 4 本 2,070 s、rerank SP 2 本 15 s | 1,771 s。opt 3 本 505 s、freq 3 本 1,177 s、rerank SP 3 本 50 s | M3（偽の AFIR 生成物の opt 405 s と freq 903 s が消えた）。SP +1 本は neutral が always でなくなったため |
| minima(dft) | A/B arm A | 3,317 s。opt 4 本 900 s（27 点）、freq 4 本 2,355 s | 2,433 s。opt 4 本 899 s（27 点）、freq 3 本 1,467 s | S1 |
| minima(dft) | NH3 / 水 | 18 s / 14 s（freq 2 / 2 本） | 12 s / 9 s（freq 1 / 1 本） | S1 |
| minima(dft) | HCN / HONO | 18 s / 54 s | 18 s / 54 s | 変化なし |
| reaction-paths（全体） | HCN / HONO / NH3 | 79 / 386 / 108 s | 82 / 516 / 101 s | HONO は GS の失敗から FIND_PATH へ |
| SCREEN | HCN / NH3 | 約 28 s / 29 s | 同じ | NH3 の tangent 用 DFT freq（5.1 s）が両方で JobStore をミス（C15 未実施） |
| SCREEN → FIND_PATH | HONO | SCREEN 約 48 s（GS は 150 サイクルで未収束・TSOpt は収束、DFT SP 12 点） | GS 失敗 8.0 s（S4 の再実行も同じ例外）→ FIND_PATH 1 チャンク 200.4 s（9 反復目に安定、残り 11 反復は判定されずに走った） | U5-N3、M2 |
| FIND_PATH | STO-3G HCN | — | 2 チャンク 16.5 s + 17.3 s（2 チャンク目の端点は +8.85 kcal/mol の非極小、U5-N1） | M2 |
| REFINE_SADDLE | HCN / NH3 | 4.6 s・5 ステップ / 2.9 s・3 ステップ | 4.7 s・5 / 2.9 s・3 | 変化なし（xTB Hessian） |
| REFINE_SADDLE | HONO | 10.2 s（xTB Hessian、5 ステップ） | 種で DFT freq 15.1 s + saddle 48.4 s（34 ステップ、AUTOZ 失敗）+ 57.4 s（36 ステップ、noautoz、途中 +82 kcal/mol）= 約 121 s | path_hei の種（2 次鞍点の領域）、U6-N1 |
| VALIDATE_TS | HCN / HONO / NH3 | 5.9 / 15.2 / 5.1 s | 5.8 / 15.3 / 5.1 s | 変化なし（HONO の TS は W7 と 2.5e-7 Eh で一致） |
| CONNECT（QRC） | HCN | 19 / 24 = 43 ステップ、40.1 s | 24 / 20 = 44 ステップ、41.7 s | S7 の効果なし |
| CONNECT（QRC） | HONO | 126（maxiter で継続）/ 17 = 143 ステップ、310.9 s | 34 / 36 = 70 ステップ、162.2 s | S7 |
| CONNECT（QRC） | NH3 | 50 / 50 = 100 ステップ、69.6 s | 45 / 45 = 90 ステップ、62.2 s | S7 |
| thermo | 全 run | 0〜1 s。worker 0.22〜0.32 s/subject | 0〜1 s。0.23〜0.28 s/subject、失敗 0 | キーが変わったので再利用なし（X12） |
| report | 全 run | 0〜1 s、ジョブ 0 | 同じ | — |

時間の大部分は今も DFT の解析 Hessian です（17 原子で 1 本 902 s、CPHF 5 反復で約 62%。TMA·(HF)₂ の dft stage の 66%、discover 全体の 48%）。

### 4.3 値の変化（W7 → v2、kcal/mol）

| 量 | W7 | v2 | 差 | 原因 |
|---|---:|---:|---:|---|
| HCN ΔE‡ | 46.6215 | 46.6215 | 0 | — |
| HCN ΔG‡ / ΔZPE‡ | 42.2178 / −3.3349 | 42.1789 / −3.4204 | −0.039 / −0.086 | S15 |
| HCN ΔG_rxn | 12.6848 | 12.6825 | −0.002 | S15 |
| HONO ΔE‡ | 13.671 | 13.671 | 0 | TS は FIND_PATH 経由（−681.5i → −678.1i） |
| HONO ΔG‡ / ΔG_rxn | 12.2430 / 0.0770 | 12.2237 / 0.0759 | −0.019 / −0.001 | S15 |
| NH3 ΔG‡ | 3.8505 | 3.8407 | −0.010 | S15 |
| TMA·(HF)₂ ΔG_assoc | −11.594 | −11.527 | +0.067 | S15 |

- band の幅は 3 反応とも 6e-4 kcal/mol 未満で変わりません。順位は HCN・HONO・NH3 がすべて rank 1、水と TMA·(HF)₂ は same_basin で rankable=False。
- v2 の job.json を vib_scale = 0.989 で再計算すると W7 の値を ±0.003 以内で再現し、差は全部 S15 によるものと確かめました（M6 の寄与は 0）。
- 文献の ZPE 因子 0.9848 にした場合の差: HCN ΔG‡ +0.054、HONO +0.022、NH3 +0.013、ΔG_assoc −0.089。

### 4.4 本レビューのプローブ（validation の外で確かめたこと）

場所は WSL `/home/user/hfauto_review2_probe/`（v2 の run の複製への追記か、独立した小さな計算）。多くは他のプローブと同時で load 8〜15 だったので、所要時間は参考値です。

| 対象 | 内容 | 結果 |
|---|---|---|
| 手法パネル（HONO） | v2 の複製に HEAD の method_panel を追記（6 SP） | ΔE‡ / ΔE_rxn: SVPD 13.67 / +0.13、PBE0/TZVPD 13.24 / +0.24、ωB97X-D3/TZVPD 12.82 / +0.33。sign_disagreement=False、panel_report と report の ranking.csv がどちらも rank 1（12.2237） |
| 手法パネル（HCN） | 既定パネル 6 SP、次に CCSD(T)/TZVPD を 1 語 opt-in | 6 SP・36 s（load 0.2 → 2.7）。CCSD(T) は hits 6 / misses 3、11〜13 反復、1 点 187〜252 s（load 約 12）。CCSD(T)/TZVPD 47.81 / 15.22（ANO 48.3 / 14.7） |
| composite（HCN） | 複製に comp_sp（all_minima、ωB97X-D3/TZVPD）→ comp_thermo → comp_report を追記 | すべて JobStore のヒットで 2.2 s。ΔE‡ 46.37、ΔG‡ 42.18 → 41.93、rank 1 |
| composite の文書の手順 | 新しい run-dir で sp → thermo → report | `ValueError: stage 'comp_sp' (sp) has no input of type ['minimum']` で即座に停止（U0-N2） |
| S9 | 未計算の pbe0-d3bj_def2-tzvpd を energy_method に指定 | 3 subject とも G=None、[thermo_unavailable, energy_layer_missing]。ただし dE_act には freq 層の 46.62 が入る（U0-N1） |
| S17 | (1,2) と (0,3) の CCSD(T) | どちらもジョブ 0 件で INPUT_INVALID `wft_closed_shell_only` |
| M7・S11・S12・C5 | 不正な system の読み込みと structures | M7 は exit 2・0.79 s、S11a 3 件・S12 5 件を拒否、S11b は 1.56 s で停止。'Xx'・'D' は unknown element |
| SMILES の同位体（C4 未実施） | '[2H]O[2H]'、'[13CH4]' | それぞれ H2O_q0_m1（水と同じラベル）、CH4 として通る（U1-03） |
| CREST 単量体の停止 → 再実行 | グリシン・β-アラニンの双性イオンで hfauto の `_search` をそのまま実行 | 停止の検出は正常、`--noreftopo` の再実行は 2/2 失敗。停止構造からの既定の実行は 2/2 成功（5 / 7 配座、60〜66 s） |
| CREST の電荷・開殻の組成 | NH3·HF⁺（m=2）、NH4⁺·HF、NH4⁺·(HF)₂、H3O⁺·H2O | 開殻は `--noopt` で失敗（preopt は停止、noref は成功）。閉殻の陽イオンは `--noopt` で 3/3 成功 |
| DFT の mode-follow | HONO の TS から minima(dft) | TS に留まる opt（2 点）→ freq で saddle（−681.67）→ +0.1 Å 側の opt 26 点で −13.3 kcal/mol。−側以降は時間内に終わらず |
| IDPP | HONO の DFT 端点で反復上限を変える | 500 反復: 面内の反転（DFT//IDPP 最大 37.04）。5,000 反復（1,479 で収束）: ねじれ（最大 14.82、HEI 候補の xTB 虚振動 1 本） |
| GS | v2 の HONO の入力で cart / 面外ずらし | cart は面内の経路（HEI 26.56）で TSOpt も例外。面外に 0.01〜0.05 Å ずらすと 3/3 完走、HEI 11.6 と TS を再現 |
| QRC の trust 0.3 | v2 と同じ開始構造・TS Hessian | HCN 24/20 → 15/12、HONO 34/36 → 14/15、NH3（片側）45 → 24 ステップ。最終エネルギーの差 ≤ 6e-8 Eh、全側で max(traj) = 開始点 |
| saddle の moddir 1 | HONO の path_hei の種、HCN・NH3・HONO（W7）の xTB 種 | noautoz + moddir 1: 70 → 18 ステップ（同じ TS）。xTB 種: 5 → 4、3 → 2、5 → 4 ステップ |
| M6（GoodVibes 実物） | TMA·(HF)₂ の最低モードを ±x に | HEAD は G(−x) − G(+x) = 0.0000、旧来は +0.74〜+1.61 |
| NT2 の失敗メッセージ | v2 の monotonic_uphill 4 件を再実行 | すべて 'No transition state guess was found in Newton Trajectory scan.' |
| git（WSL の /mnt/c） | `git status --porcelain --untracked-files=no` | 既定 83 件、`-c core.autocrlf=input` で 0 件（U0-N7） |

### 4.5 実計算でまだ通っていない経路

- 決定表: 行 12 の barrierless（SCREEN）、行 13（multi_max → VALIDATE_INTERMEDIATE）、行 14（monotonic → BARRIERLESS、S3）、低レベル TS の近道、SCREEN の崩壊時の IDPP。
- M1 の経路（旧ゲートが発火しうる仮説は v2 で 0 件）。
- U6: higher_order のやり直し、collapsed と分割、QRC の再試行（C19）、`_register`、reassigned、periodic_nearest、moddir ≥ 1（smoke だけ）。
- U3: DFT の mode-follow の全体（プローブは +側の opt まで）、DFT の soft 押し。
- U7: M6 の |ν| 化（実 run に負の副次モードがない）、S14（電荷・開殻の会合）、S9（プローブだけ）。
- U1・U2: 電荷のある系・開殻系を structures から NWChem（odft）まで通すこと、CREST の単量体の停止 → 再実行（プローブで失敗を確認）。
- U8: 手法パネルは v2 の検証と主対象（amine·(HF)n）では未実行（プローブの HCN・HONO だけ）。

---

## 5. 今後の改良候補

有用性が確認できたもの（検証で「採用」「修正して採用」とされたもの）だけを載せました。工数は実装とテストの目安、「複雑さ」はコードと設定の量の増減です。どの改良も新しい knob は増やしません。

### 5.1 must

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|---|
| M8 | FIND_PATH のチャンクの継ぎ目と最終経路で、凍結端点を DFT 極小（ずれた bead に整列した ctx.ends）に戻し、ずれを note に残す | U5-R1（= U6-R3 の後半）、U5-N1 | 2 チャンク以上の string でも monotonic → BARRIERLESS と形の判定の前提（端点が真の極小）が成り立つ。追加の計算はない | 3〜8 行、scripted PathEngine のテスト約 15 行 | 中立 | FakePath は initial_path を無視するので、テストは initial_path を記録する scripted engine で書く。1 チャンク内のずれは note で監視することを文書に書く |

### 5.2 should

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|---|
| S18 | IDPP の反復上限 500 → 5,000 | U5-R2、U5-N2 | 平面の端点で FIND_PATH がねじれの機構を捉える（HONO の DFT//IDPP 最大 37.0 → 14.8）。U6 の DFT Hessian と AUTOZ の再試行を避けられる見込み | 定数 1 つ、テスト 1 件 | 中立 | 効くのは IDPP が初期経路になる退避時と SCREEN の崩壊時だけ |
| S19 | FIND_PATH の画像を `align_sequential` で整列する | U6-R3 の前半、U6-N2 | path_hei の接線と HEI の補間から剛体回転が消える（STO-3G の重なり 0.118 → 0.651）。誤った note と xTB Hessian の誤った棄却（17 原子で約 15 分）がなくなる | 1 行 | 中立 | M8 と同じ箇所なので一緒に入れる |
| S20 | 鞍点で追うモードを明示する（常に moddir、負モード 2 本以上は noautoz と P·H·P の順、継続では moddir を 1 か外す、key_payload に moddir と座標系） | U6-R1（C17 の拡張）、U6-I2、U6-N1 | HONO の path_hei の種で 70 → 18 ステップ（プローブ）、登り返しがなくなる。別の鞍点に迷い込むリスクが減る | 約 15 行、テスト 2 件、再実行約 15 分 | 中立 | `_first_hessian` の計量は変えない |
| S21 | 初期 Hessian を渡す optimize は `trust 0.3`（渡さない opt は 0.1 のまま） | U6-R2、U6-I1、U6-N4 | QRC のステップがほぼ半分（プローブ）。v2 の QRC 合計 約 266 → 約 110 s。17 原子で 1 反応約 30 分の短縮（推定） | 1 行と文書・golden テスト、再検証約 1 時間 | 中立 | minima の 2 断片以上の opt にも効く。ab_init_hessian と known_endpoints を**新しい run ディレクトリで**再実行し、basin の数・同定・最終エネルギー・trajectory_above_ts を確かめる |
| S22 | rankable: 順位を付けられない outcome には outcome の理由と本当の blocker だけを返し、dzpe が None なら blocker にしない | U8-R1（修正版）+ U0-R5、U8-8、U8-N3、U0-N6（C26 を置き換え） | 主対象の報告行が `outcome:same_basin` だけになり、ΔG_assoc の隣の誤った thermo_unavailable がなくなる | 3〜4 行、テスト 1 行 | 中立 | spin_contaminated などの本当の blocker は残す |
| S23 | composite の手順（複製して追記、targets: all_minima、stage id を分ける、SP 層の ΔE‡ ≤ 0 は読まない）、手法誤差の数値（HCN −1.2、HONO +2.0、目安約 2）、パネルの conditions をそろえる旨を文書に書く | U0-R1（修正版）+ U8-R6(a)、U0-N2・N3・N5、U8-N1 | composite が数秒〜数分で使える経路になる（全段の再計算 約 31 分を避ける）。順位差の読み方と「panel_report が最終」の前提が正確になる | 文書 4〜5 行 | 中立 | NH3 の数値は書かない |
| S24 | SMILES の同位体を拒否する | U1-R1、U1-03（C4 の同位体の部分） | 対象元素 H で黙って誤る最後の経路を閉じる。xyz の 'D' の拒否とそろう | 約 2 行、テスト 1 行 | 中立 | design.md:138 も直す |
| S25 | stage の入力切れの例外に、上流の failed artifact の理由（最大 3 件、artifact_id 付き）を添える | U1-R2、U1-N1 | known_endpoints で入力誤りがあったとき、traceback の最終行に理由が出る | 約 4 行、テスト 1 件 | 中立 | — |
| C2 | CCSD ブロックに maxiter 50 を書く（mp2 には書かない） | U0-R2 = U8-R2、U0-I8、U0-N4 | opt-in の目安（重原子 4 個程度）の系でも、参照点を落とさなくなる（HONO TS/TZVP 18/20 反復） | 1 行、テスト 1 行 | 中立 | U0 の検証は should、U8 の検証は could。既存の成功結果は再利用される。未収束で失敗した既存のジョブは `--retry-failed` で |

### 5.3 could

| # | 改良 | 由来 | 効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| C8 | 組成の多重度が一意でなければ INPUT_INVALID | U2-R4 | ラジカル対の組成で高スピンの PES を黙って選ばない | 2〜3 行 | 中立 |
| C10 | window で screen 極小がなければ拒否、level=screen で select の明示を拒否、design.md:142 を「all は入力構造から」に | U3-R2、U3-ISS-5、U3-N4 | 空の dft stage を「反応なし」と読み違えることを防ぐ | 数行、テスト 1 件 | 中立 |
| C12 | drive の重複を WL 原子クラスで除き、列挙数・固有数・実行数を explore の diagnostics.json に出す | U4-R1（U4-R2-1 と R2-2 の統合）、U4-I4 | TMA 系の attempt が約半分、heavy_bond・association の段に届く。打ち切りが見える | 15〜25 行、テスト 2 件 | 中立 |
| C13 | 無バイアスの opt に `convergence_max_iterations=500` | U4-R2、U4-I5 | 端点の opt の打ち切りによる数値的な陰性（amine 8 件中 6 件）とラベルの雑音がなくなる | 1 行 | 中立 |
| C14′ | TS がある negative にも ΔE‡ を manifest に載せる（SCC やり直し温度では 300 K で評価し直すか載せない） | U4-R3、U4-I6、X1 | 陰性の結論を障壁の値で裏付けられる | 4〜8 行 | 中立 |
| C15 | tangent_mode_cm1 を barrierless のときだけ計算 | U5-R6、U5-I8 | 縮退・置換した端点での余分な DFT freq（17 原子で約 900 s）を省く | 約 5 行 | 中立 |
| C16 | 両端より低い DFT//xTB ノードがあれば unavailable にして FIND_PATH へ（M8 の後、fake の経路テスト必須） | U5-R7、U5-I7 | SCREEN が中間体のある経路を BARRIERLESS で閉じることを防ぐ | 3 行、テスト 1 件 | 【増】 |
| C20 | 点群を ThermoResult.notes で記録し、S_rot を protocol から外す | U7-R4、U7-I7、U7-N6 | σ を監査でき、死んでいる経路が 2 つ片づく | 約 4 行 | 中立 |
| C21 | ΔG_assoc の感度の幅（表示だけ） | U7-R6、U7-I5 | 弱い錯体の判断に効く不確かさ（298 K で幅 1.155）が出る | 約 8 行 | 【増】 |
| C25 | ranking.csv に T_K・standard_state・dG_rxn_kcal | U8-R3、U8-7、U8-N5 | CSV が単独で読め、metric=dG_rxn で並べた値も出る | 数行、テスト 1 件 | 中立 |
| C27 | ranks × memory_mb_per_rank ≤ site.memory_mb を preflight で検査 | U8-R4、U8-9 | スケールアップ時の OOM を計算前に止め、効かない knob がなくなる | 3〜5 行、テスト 1 件 | 中立 |
| C28 | layer_missing の subject を含む ΔE を None にする | U0-R4 = U7-R1、U0-N1 | 2 つの LOT を混ぜた ΔE（HCN で 164 kcal/mol）が表示されなくなる | 1 行、assert 1 つ | 中立 |
| C29 | ΔG_assoc=None の理由（no_composition / monomer_missing / monomer_thermo_unavailable / mixed_level）を ReactionThermo.notes に。同じ (Hill, 電荷) の曖昧な対応は None | U7-R2（修正版）、U7-I2、U7-N2・N3 | F⁻ 系と熱化学の失敗を区別できる。誤った基準の ΔG_assoc を出さない | 約 10 行、テスト 1 件 | 【増】 |
| C30 | 設定値の範囲検査（ConformersConfig の ge=1 / gt=0、Conditions の正の温度・非空・重複禁止、ThermoSettings の vib_scale > 0、cutoff ≥ 0） | U2-R3 + U7-R3、U2-N4、U7-N4 | 数時間の run が最後に無駄になる設定の誤りを読み込み時に止める（X8） | 約 8 行、テスト 2 件 | 中立 |
| C31 | 読み込み時に xyz ファイルの存在と、宣言反応の両端の電荷・多重度の一致を検査 | U1-R3、U1-N2・N4 | 'running' の残留と dry-run の素通りを塞ぎ、片側だけの typo による無駄な DFT を防ぐ | 約 5 行、テスト 2 件 | 中立 |
| C32 | rerank_sp の SP を per_state を超える組だけに出し、SP の失敗を diagnostics に残す | U3-R1（修正版）、U3-N1 | 単量体の最低構造が SP の失敗で黙って消えない。TMA で 50 s 減 | 4〜5 行、テスト 1 件 | 中立 |
| C33 | reaction-paths の `_register` に pipeline の gates を渡す | U3-R3、U3-N2 | 1 つの PES の tier の基準が 1 つにそろう | 1 行 | 中立 |
| C34 | 単量体の停止後の再実行を「停止構造から既定の設定で 1 回」にし、`--noreftopo` を廃止 | U2-R1、U2-N1 | 実 CREST で 2/2 失敗していた経路が 2/2 成功（プローブ）。フィールドと CLI フラグが 1 つずつ減る | 約 10 行、fake テスト 1 件 | 【減】 |
| C35 | pysis_gs の off_plane（dlc の 1 断片平面分子だけ） | U5-R3（修正版）、U5-N3 | 1e-6 Å の差で GS の成否が変わる非決定性がなくなる（HONO で 3/3 完走） | 約 8 行、テスト 1 件 | 【増】。S18 を入れた後に採否を再判断 |
| C36 | S4 の再実行を GS が完了した場合（final_geometries.trj あり）に限る | U5-R4、U5-N4 | GS 自体が落ちたとき同じ GS を繰り返さない | 1〜2 行 | 中立 |
| C37 | 使われなくなった monitor 経路と FailureKind.STAGNATED を削除 | U8-R5、U8-N2 | 到達しないコード約 30 行と失敗の種類 1 つが減る | 約 30 行の削除、テスト 3 か所 | 【減】 |
| C38 | code_version の git に `-c core.autocrlf=input` | U0-R6、U0-N7 | WSL でも `-dirty` が本当の未コミット変更だけを示す | 引数 2 つ | 中立 |
| C39 | 文書の訂正（まとめて） | U0-R3 = U8-R6(b)、U1-R4、U2-R2、U5-R5、U6-R4、U7-R5 | 実装・実測と文書の食い違いをなくす | 文書 10〜15 行 | 中立 |

C39 で直す文書の内訳:

| 場所 | 内容 | 由来 |
|---|---|---|
| validation.md（新しい小節） | 手法パネルと composite のプローブ結果（HONO で両 ranking.csv が rank 1、HCN 6 SP・36 s、CCSD(T) 47.81 / 15.22、composite 2.2 s、S9・S17 の拒否）。負荷時の時間は参考値と断る | U0-R3、U8-R6(b) |
| validation.md:32 | 「127 / 17 = 144」→「126 / 17 = 143」 | U6-R4 |
| validation.md:39 | ΔG の変化の原因「M6・S15」→「S15（M6 の寄与は 0）」 | U7-R5 |
| validation.md:54 | −289i の TS の判定の変化は M3 の効果（縮退した HF 交換） | U4-N1 |
| validation.md §8-1・§8-2 | GS の対策案（Cartesian、単純な再試行）は無効 → S18・C35。最初の例外は stderr.txt に連鎖で残る | U5-R5 |
| design.md:130 | NWChem の string は gmax を記録しない。M8 の後は端点を戻す旨の 1 行 | U5-R5 |
| design.md:139（conformers 行） | 開殻の組成は `--noopt` で trial MTD が失敗しうり、そのときは seed が出る | U2-R2 |
| design.md:141（explore 行） | 状態ラベルは結合グラフなので、配座・立体・E/Z だけ違う IRC の端は same_as_source として捨てる | U1-R4 |
| design.md:146（thermo 行） | vib_scale = 1.0 の根拠と文献値（0.985、差は ≤ 0.06 / 0.09）。単原子の単量体・反応のない錯体では ΔG_assoc は出ない | U7-R5 |
| design.md（reaction-paths） | 子反応の予算は case ごと（最悪 7 case × 6 h）、qrc_bounds_A の上限は 0.5 Å 以下 | U6-R4 |
| engine.continuation の docstring | autoz の失敗は種から Cartesian でやり直す | U6-R4 |
| environment.md:64 | F⁻·(HF)₂ の CREST「約 25 s」→ 約 1.4 s。crest の threads の既定は 1 | U2-R2 |

前回の could のうち今回再提案しなかったもの: C1（見送り、§3.0.7）、C4 のラジカルの部分（見送り、§3.1.7）、C18（periodic_nearest。実 run で到達例がなく、前回の判断＝could のまま据え置き）、C22（方針どおり見送り）。

### 5.4 進め方と確認

1. **M8 → S19 → S18（FIND_PATH の一連）**: 同じ `actions.find_path` と `interpolation` を触るので、続けて入れます。確認は、scripted PathEngine のテスト（端点のずれを入れた 2 チャンク）と、新しい run ディレクトリでの find_path_sto3g（2 チャンク目の端点が真の極小になり、note に 1 チャンク目のずれが出ること）と HONO の再実行（GS が失敗した場合の FIND_PATH がねじれの HEI から始まること）。C16 を入れるならこの後。
2. **S20 + S21（NWChem の driver）**: deck のパラメータは JobStore のキーに入らないので、**必ず新しい run ディレクトリで**、known_endpoints の 3 反応・find_path_sto3g・ab_init_hessian（arm A）を再実行します（合計約 1 時間）。確かめること: TS のエネルギー（≤ 1e-7 Eh）、QRC の割付けと trajectory_above_ts、minima の basin の数と同定、最終エネルギー。結果を validation.md と U3 の記録に残します。
3. **S22〜S25 と C2**: どれも数行で、報告・文書・入力検査です。S23 と C39 の validation.md の追記はまとめて行えます。
4. **検証の空白を埋める最小の実計算**（任意）: v2 の HONO の複製で `hfauto run method_panel` を正式に記録する。M8 の後、monotonic になる軽い系（STO-3G など）で行 14 を 1 回通す。

### 5.5 あえて行わないこと

- **既定の流れを重くする変更**: 停留点レベルを TZVPD や ωB97X-D3 に上げる、会合量の CP 補正の自動化、CCSD(T) を既定のパネルに戻す・系の大きさで自動的に入れる、composite を既定の流れに組み込む、composite.yaml の同梱（C1）。効果に対して費用（17 原子の freq 約 900 s/本、CCSD(T) 1 点数時間）が見合いません。
- **knob を増やす変更**: 符号ゲートの不感帯・QRC_AMPLITUDES・_MIN_OVERLAP・HESSIAN_NEAR_A などの設定化、スケール因子を method ファイルへ、既定を 0.985 に変える。新しい knob を増やさない方針（X8）を守ります。
- **閾値や規則を増やす変更**: パネルの幅を同順位の band や blocker に使う、端点のずれにゲートを付ける、qrc_bounds_A の読み込み時検査、CCSD(T) の関数数の上限、h_shift の原子価フィルタ、matches_source に RMSD を併用、assign で全候補にエネルギー条件。
- **根拠のない性能改善**: 全 optimize の trust 0.3（Hessian なしでは TS より上へ戻る）、TS Hessian の負の固有値の反転、mode-follow・soft 押しへの Hessian、IRC の反復上限の変更、ZTS の小さなチャンク、CREST の決定化・6 seed すべての探索、DFT freq の grid・CPHF・RI の変更。
- **対象外の化学への対応**: 開殻 WFT、単原子種の全面対応、開殻の組成だけ preopt に戻す、SMILES のラジカル電子数の照合、状態ラベルへの立体の追加、1.45 Σr_cov の変更。
- **表示の配管を増やす変更**: thermo の artifact id に stage id、energy_layer_missing を blocker に、ジョブ失敗を coverage.csv に、報告物の統合、population・H の削除。

---

## 6. 参考文献・仕様 URL

前回の §5 の一覧（Bursch/Grimme 2022、Rappoport/Furche 2010、Pracht/Bohle/Grimme 2020、Smidstrup ら 2014 の IDPP、E/Ren/Vanden-Eijnden 2007 の ZTS、Goodman & Silva 2003 の QRC、Grimme 2012 の qRRHO、Luchini ら 2020 の GoodVibes、Alecu ら 2010 のスケール因子、Fernández-Ramos ら 2007 の対称数、各ライブラリの文書とソース）は引き続き有効です。今回あらたに参照したもの、または再確認したものを挙げます。

**文献**

- Kesharwani, Brauer, Martin, "Frequency and Zero-Point Vibrational Energy Scale Factors for Double-Hybrid Density Functionals (and Other Selected Methods)", J. Phys. Chem. A 2015, 119, 1701 — https://webhome.weizmann.ac.il/home/comartin/OAreprints/260.pdf （PBE0/def2-SVPD: 調和 0.9930、基音 0.9561、ZPVE 0.9848）
- Tikhonov ら, ChemPhysChem 2024, 25, e202400547（PBE0-D3(BJ)/def2-SVP の基音のスケール因子 0.9555）
- Minnesota Database of Scale Factors v5（2021）— PBE0/MG3S の ZPE 0.975・H 0.989・F 0.950
- HCN/HNC の CCSD(T)/ANO 参照値 — https://cdnsciencepub.com/doi/pdf/10.1139/v96-120
- Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. 140, 214106 (2014)（IDPP。目的関数を収束させた経路が前提）
- Pracht, Bohle, Grimme, PCCP 2020, 22, 7169（CREST / iMTD-GC）— https://pubs.rsc.org/en/content/articlelanding/2020/cp/c9cp06869d

**ライブラリの仕様・ソース**

- NWChem CCSD — https://nwchemgit.github.io/CCSD.html
- NWChem DFT — https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html
- NWChem Geometry Optimization（MODDIR、FIRSTNEG、INHESS、TRUST）— https://nwchemgit.github.io/Geometry-Optimization.html
- NWChem NEB / ZTS（IMPOSE）— https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html
- NWChem 7.2.3 のソース:
  - driver（初手 uphill・モード選択 3255-3360、負モードの刻みの制限 1745-1760、trust の既定 1066-1074、INHESS 373-385）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/driver/opt_drv.F
  - ZTS（凍結 bead 1829-1831・1886-1900、zts_min_motion 1965-1968・2005-2008、string_final の書き出し 897-915）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F
  - L-BFGS（1259-1310）— https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/neb/neb_utils.F
  - CCSD の入力（maxiter）— https://raw.githubusercontent.com/nwchemgit/nwchem/v7.2.3-release/src/ccsd/ccsd_input.F
- CREST 3.0.2 のソース:
  - https://raw.githubusercontent.com/crest-lab/crest/v3.0.2/src/confparse.f90 （`-noopt`、`-nopreopt`、`-notopo`、`-quick`、`-nci`）
  - https://raw.githubusercontent.com/crest-lab/crest/v3.0.2/src/crest_main.f90 （preopt が真のときだけ trialOPT）
  - https://raw.githubusercontent.com/crest-lab/crest/v3.0.2/src/algos/setuptest.f90 （trialOPT、trialMD_calculator）
  - CREST のキーワード — https://crest-lab.github.io/crest-docs/page/documentation/keywords.html
  - GitHub issues — https://github.com/crest-lab/crest/issues/285 、https://github.com/crest-lab/crest/issues/419
- SCINE:
  - ReaDuct 6.1.0 NtOptimization2Task — https://raw.githubusercontent.com/qcscine/readuct/6.1.0/src/Readuct/App/Tasks/NtOptimization2Task.h
  - Utilities 10.1.0 NtOptimizer2 — https://raw.githubusercontent.com/qcscine/utilities/10.1.0/src/Utils/Utils/GeometryOptimization/NtOptimizer2.cpp
  - Utilities 10.1.0 GradientBasedCheck — https://raw.githubusercontent.com/qcscine/utilities/10.1.0/src/Utils/Utils/Optimizer/GradientBased/GradientBasedCheck.h
- xtb — https://xtb-docs.readthedocs.io/en/latest/commandline.html 、https://xtb-docs.readthedocs.io/en/latest/optimization.html
- RDKit Getting Started（ETKDG と最小化、AddHs）— https://www.rdkit.org/docs/GettingStartedInPython.html
- pysisyphus 1.0.0（WSL prod venv の `optimizers/Optimizer.py:385`、`optimizers/StringOptimizer.py:120-139`、`run.py:1505-1527`）
- GoodVibes 4.3.0（WSL prod venv の `goodvibes/api.py:85-160, 296`、`goodvibes/thermo.py:337-376, 695-790, 892`）、pymsym 0.3.5
- Python の例外の連鎖 — https://docs.python.org/3/tutorial/errors.html#exception-chaining
- 3Dmol.js（jsDelivr）— https://data.jsdelivr.com/v1/packages/npm/3dmol 、https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js

**リポジトリ内の関連文書**

- `docs/reviews/2026-09-26_calculation_spec_review.md`（前回、ベースライン）
- `docs/design.md`、`docs/validation.md`（v2）、`docs/environment.md`、`README.md`

**実 run とプローブ**（WSL、読み取り専用の複製）

- 改良前: `/home/user/hfauto_w7/<run>`（コード `6b87814`）。改良後: `/home/user/hfauto_v2/<run>`（コード `4e972c5`）。
- 本レビューのプローブ: `/home/user/hfauto_review2_probe/{u0, u0_assess, u1, u1_review, u2_gly, u2_bala, u2_cation, u2_ion, u3_mf, u4, u4_assess, u4_verify, u5_reassess, hono_gs, hono_gs2, u6, u6_assess, u7, u7_assess, u8, u8_reassess}`
- 前回のプローブ: `/home/user/hfauto_review_probe/`
