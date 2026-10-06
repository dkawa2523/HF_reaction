# round 10 の改良履歴

> 2026年10月6日の整理前の計画と変更履歴を保存した資料。表の完了は当時の記載で、現在の統合・検証の完了を示さない。現在状況は [roadmap.md](../roadmap.md)、次の改良は [improvement-plan.md](../improvement-plan.md) を正とする。本文の旧設計書の節番号は当時のもの。


レビュー [2026-10-04_platform_review_2.md](../reviews/2026-10-04_platform_review_2.md) §7 の計画を、段階 0 のプローブ(P0a〜P0e。証拠は `/home/user/hfauto_r10/probe/`)の結果で決めた形で進める。目標の設計は [design.md](../design.md) §11 で、wave を統合するたびに §0〜§10 に取り込む。実計算の記録は [validation.md](../validation.md) §11。

## 1. 提案と wave

| wave(段階) | 提案 | 状態 |
|---|---|---|
| W1(1、must) | X7-1 予算は件数だけ、X7-2 の一部(外からの停止は incomplete・rc 1)、U3-P1 σ を対称操作から、X4-4 の must 部分(状態の G を閉じる)、U1-P1 多重度の宣言を必須に、X7-7 検証を repo に(cases.yaml、BH76 の参照、compare、replay) | 完了 |
| W2(1、must) | U6-P1 χ の質量加重、U6-P2(1 つの結合変化の共有。付け替えは棄却)、U6-P6(押した探索の guess)、U6-P9、X3 SCF の来歴と救済の 3 段(P0a)、開殻 CCSD(T) は UHF 参照(P0d)。最後に再生の基準 B1 を作る | 完了(U5-P3a は W6、QRC の側の guess は W7 へ) |
| W3(2、must) | X5-1 反応が読む点、X4-3 鎖と共通のゼロ(BH76 の分離のゼロ用に h3_doublet と oh_ch4 の単量体を宣言)、U5-P2 捕獲律速、U5-P1、U5-P7、U2-P2・U2-P4 | 完了(U5-P7 は測って不採用、U2-P4 は円錐と表の削除だけ、S5 split2 は新しい run で未達。変更 28〜37) |
| W4(3、must) | U1-P3・U1-P2 正準のラベル、U4-P1 結合グラフの編集の列挙(電荷分離の規則は P0f)、U4-P2・U4-P4・U4-P7、U4-P3 DFT//xTB の窓、U3-P7 | 完了(P0f は不合格で P0c の規則。gates の帯付きの比較は W2 で済み。変更 38〜50) |
| W5(4、should) | X1-2 停留性の認証、U3-P4、U6-P4、U5-P5、U5-P6、X2-2・X2-3 対称性の受理と電子準位、U7-P6、U1-P4、X7-2 の残り(rc 2、項目の閉じ込め)、X7-3、U9-P3、U9-P2(a)、エネルギー層の選び直し(BH76) | 完了(層は M06-2X のまま、revM06 は測った候補。判断表は 9 行、drivers のコード行は +17。未実証は VAL10 へ: 変更 51〜64) |
| W6(4、should) | X4-2 事実の登録簿(辺は自分のラベル付きの両端を持ち、貸された TS の χ をそれで測る: U5-P3a)、U5-P4(近道の TS を DFT SP の順に)、U8-P10・U8-P7・U8-P6、X5-3 層の経路の判定 | 未着手(W5 で記録の欄 `steps`・`layer_verdict`・`dE_act_path_kcal`・`BarrierVerdict.points` と 1 つの止める規則・継続・下りは済み) |
| W7(最後) | 鍵を変える基盤の変更を 1 回で移行(X7-6、U9-P10、X1-3、U0-P3、QRC の側の opt の `scf_guess`)、分割の仕組みの削除(X4-2 段階 2)、整理、VAL10 と BH76 の表の取り直し | 未着手 |

## 2. 計画の変更の記録

レビュー §7.1 からの変更。影響を受ける wave はそれぞれの作業で扱う。

1. U2-P4(配置の縮小)を U2-P2 とともに段階 2 へ移す。純粋な削除で、P0b の結果に依らない。
2. U1-P3、U1-P2、U3-P7 を段階 3 へ移す。検証はラベルの文字列ではなく状態の分割(全単射)で比べ、cases.yaml の署名は端の断片の組成式の式なのでラベルの変更に影響されない。
3. U4-P4(1 世代の辺と多世代)と U4-P7 を段階 3 へ移す。DiscoveryRecord は名前を保ち、辺の意味を持つ。
4. X5-3 は sp の要求と thermo の判定で作る(`layer_verdict`、`dE_act_path_kcal`)。新しい outcome と case の action は作らない。
5. ジョブ鍵を変える基盤の変更(X1-3、U9-P6、U9-P10、U0-P3)は W7 で 1 回にまとめる。開殻 CCSD(T) の参照の変更(9)だけは手法の定義の変更として W2 の B1 で吸収する。
6. X4-2 の b〜d は stage の中の登録簿にする。ReactionRecord は case の結果のまま。
7. 段階 5(could)の項目は、コードか重複を削るものだけ入れる。
8. P0b: GFN-FF(U2-P1)は採らない。W3-3 は GFN2 の CREST と `--notopo`/`--noopt` を保ち、組成の検査は B1 と比べる。
9. P0d: 開殻 CCSD(T) は UHF 参照にし、W5-4 から W2-2 へ移す。
10. P0c: Lewis の篩の電荷分離の規則は、W4 の前の QM なしのプローブ P0f で決める。explore の予算は件数 1 つ(`max_trials` 3000)で、切った類は `not_attempted` と報告する。
11. P0e: BH76 の分離のゼロには宣言した単量体が要るので、W3-1 が h3_doublet に H + H2、oh_ch4 に CH3 + H2O を足す。bh76_* の系は作らない。W6-3 の合格は GMTKN55 の行 33(2.0)に対して \|ΔE − 2.0\| ≤ 1.0 かつ g ≤ 0.3。
12. P0a: 派生する SP はすべて射影した guess から始める。B1 で変わる値を事前に挙げる(W2)。
13. (W1)`hfauto status` も `run` と同じ規則で rc 1 を返す(`cli._unfinished`)。W5-3 はこの関数だけで rc の体系を変える。
14. (W1)外からの停止は KeyboardInterrupt として扱い、`execute_stage` が `incomplete` を記録する(それ以外の例外は failed)。新しい例外の型は作らない。FailureKind は 9 種(budget_exhausted を削除)。
15. (W1)BH76 の参照は P0e のファイルそのもの `validation/bh76/subset.yaml` で、cases.yaml は `bh76: <反応>.<向き>` で行を指す。W5-4 はここから読む。
16. (W1)閉じた状態の G は、既存の `thermo_unavailable` を状態の極小にも広げて止める。新しい blocker は作らない。
17. (W1)多重度のない化学種は system の読み込みで拒否する(INPUT_INVALID の artifact ではない)。repo の外にあった再生の入力は多重度を書いて `validation/inputs` に写した。
18. (W1)`validation/check.py` は、run dir のない case を superseded でなければ失敗にする。
19. (W2)U6-P2 の等価な原子の付け替え(反応物の端を鞍点に Hungarian で写し、1 つの置換を両端にかける)は実装して R2 で測り、棄却した。付け替えが結合の組を変えた 2 件(S6 の Walden の ts_calc、S10 のリレー)はどちらも誤った等価原子を選んだ(χ 0.940 → 0.726、≥ 0.88 → 0.692)。幾何の写像は大きく動く原子(反転、リレー)を追えない。χ・ρ・高次の鞍点の再試行・縮退の判定は、case のラベル付きの両端の変化 `Ctx.change()` 1 つを共有する。
20. (W2)U5-P3a(貸された TS を貸し手の両端で測る)は W6 へ移す。貸し手の両端は今の記録にないので、種に結合変化を持たせても case の変化と常に等しい(何も決めない欄になる)。W6 の登録簿の辺が自分の両端を持つときに入れる。
21. (W2)U6-P6 のうち QRC の側の opt の `scf_guess` は W7 の鍵の移行へ移す。入れると R2 の paths の QRC の opt 34 case の鍵が変わり、B1 の「sp より前は新しいジョブ 0」に反する。押した探索の `scf_guess`(saddle の鍵は guess を渡したときだけ変わり、R2 に押した探索はない)は W2 で入れた。
22. (W2)CCSD(T) の対象は DFT の freq が `spin_ok` を満たす点だけ。「低スピン結合でない」は置かない: BS の点は必ず spin_ok に落ち、spin_ok を通る低スピン結合の点は結合したラジカル(CH3OO•、0.754)で UHF-CCSD(T) が正しい。
23. (W2)開殻 CCSD(T) のジョブ鍵は `reference: uhf` を持つ。持たないと B1(R2 の JobStore を写す)が ROHF の結果を同じ鍵で黙って再利用する。閉殻の鍵は変わらない。
24. (W2)SCF の救済は駆動系のジョブも attempt の始点からやり直す(vectors だけを捨て、継続の driver の Hessian は保つ)。
25. (W2)χ の事前の期待(接続した TS ≥ 0.88、S5 split2 0.69〜0.72)に対し、実測は接続した TS 0.882〜1.000 のうち 2 行が 0.727(r7・W5 の S6 の多段の親 f9bc3a1cd8 を自分の 2 結合の変化で測ったもの)、S5 split2 は 0.635〜0.714。空白は (0.0653, 0.635) で、`REACTION_MODE_MIN` 0.3 はその中にあり変えない。期待値は動かさず、実測を記録した。
26. (W2)層とパネルのスピンの定義は `gates.energy_spin_ok`(freq の spin_ok ∧ `same_spin_state`)1 つ。
27. (W2)S5 の新しい run は split2(CH2OOH• → CH2O + OH)に elementary 0.13 で答えたので、cases.yaml の期待を unresolved から変えた(理由は χ。境界の系として VAL10 で 2 回確かめる)。R2 の経路の再生(B1)では unresolved のままで、`check.py B1` は s5 だけ FAIL になる。W6 の合格条件「split2 の対に答える」は、新しい run で保つことに読み替える。

28. (W3)`thermo.monomer_states` を `declared_monomers`((Hill 式, 電荷) → 組成のリスト)にした。oh_ch4 が CH5O の組成を 2 つ宣言するため。会合と分離した点は、端の断片の状態ラベルが宣言した組成の単量体と一致するときだけ(`thermo.separated_states`)。
29. (W3)barrierless は δG_eff を持たず、`RANKABLE_OUTCOMES` は connected の outcome だけ。会合でない barrierless(h2o_hf_inversion、water_dimer_as)も序数の外になる。cases.yaml のこれらの 5 行の dG_eff を 0.0 から null にした。
30. (W3)`ReactionThermo.reference` と `RankRow.reference`、ranking.csv の `reference` 列を足した。分離した点がないか δG_eff がない行は空。
31. (W3)手法パネルは電荷と多重度を除いた LOT で点を対応させ、自分の点(`thermo.participants`)を読む。TS は connected の outcome だけが持つので、B1 の S5 split2(unresolved)の dE_act のセルが空になる。
32. (W3)補助の点(ΔG_assoc だけが読む R_sep と前駆錯体)の SP の失敗は artifact にしない(値が空になるだけ)。
33. (W3)種の注記 `thermo_unavailable` を削除した(G = None と `energy_layer_missing` で足りる)。
34. (W3)変更 1 は誤りだった。U2-P4 の接触 seed 1 つは U2-P1 を前提にし(レビュー §6)、接触の seed では S5・S6 の GFN2 の screen が失敗して両系の反応が消える。W3-3 は U2-P2 と円錐・元素の表の削除だけを入れ、`seeds_per_composition` と乱数を残した。U2-I1(錯体が screen の basin なしに explore・DFT に届く経路)ができたら再検討する。W4-2a・W4-2b の S5・S6 の出発点は、乱数の seed が screen で崩れた basin である。
35. (W3)近道は種だけを出し、判定と ⟨S²⟩ の検査を持たない。DFT SP の順の並べ替えは W6-1(U5-P4)へ移した(今並べ替えると B1 の S6 の 2 つの TS の順と鍵が変わる)。種が待っている間は SCREEN を始めない(行 9)。
36. (W3)U5-P7(pysisyphus の IDPP)は QM なしのプローブで測って採らなかった(平面の HONO で面内 177.8° の経路、作業ディレクトリへのログ)。
37. (W3)cases.yaml の期待を理由とともに変えた([validation.md](../validation.md) §11.3): h3_doublet 12.71 → 15.55 と S6 の引き抜き [5.5, 6.0] → [9.95, 10.45] (ゼロが宣言した分離の単量体になった)、S6 split1 [20.55, 21.05] → [21.15, 21.65] (宣言した CH3 + H2O の組成が G で 0.60 kcal/mol 低い錯体の配座を出した)、barrierless の 5 行は null(変更 29)。golden G36・G37・G39 を作り直した。S5 の split2 は W3 の新しい run で unresolved(QRC の側が鞍点に落ちた)だったが、期待は変えず未達として記録し、W5-1 と VAL10 で扱う。

38. (W4)P0f(QM なし)は不合格: 宣言した電荷分離の制限は、P0c が見つけたギ酸エチルの 1,5-H 移動の生成物(1,5-双性イオン)を失い((i))、ibuprofen は 5,200 類で freq 1 本を超えた((iii))。事前の宣言どおり P0c の規則(電荷分離を制限しない)を使い、薬物の大きさの網羅は予算で限られると報告する(ibuprofen 10,708 類、3000 件の予算で f2b2 の 7,708 類が `not_attempted`。report に budget-limited と出す)。
39. (W4)類は編集のラベルを付けた出発のグラフの正準 SMILES(自己同型の軌道そのもの)にした。`topology.atom_classes` の組では、ベンゼン環の ortho と para の C–C 形成のような等価でない編集がまとまる。trials は atom_classes を使わない(identity.basin_coords だけが使う)。
40. (W4)世代の数え方: 出発の状態は screen の状態と入力の化学種の状態。S5 の CH3O2 は screen の状態なので、レビューの第 2 世代(CH2OOH)・第 3 世代(β 開裂)は第 1・第 2 世代になる。
41. (W4)`gates.discovery_verdict` を削除した。worker は射影した虚モードが 1 本の TS と Hessian で確かめた両端の極小のときだけ辺を返すので、呼んでも何も棄却しない。U1-P2 の gates の呼び出し(`_bonds_exchanged`)は W2 で帯付きの `Ctx.change()` の比較になっており、4 つの構造を取る形には変えない(gates の import 規則を破り、両側に共通の傍観の結合変化が打ち消されなくなる)。
42. (W4)DFT の入口の数の上限は minima(dft) の `select.max_edges`(組成あたり 6、端の組の数)1 つで、`hypotheses.select` の `max_per_composition` を削除した。`rerank_sp` の knob を削除した(混んだ状態は常に SP で並べ直す)。minima の relaxation の seed の写し(`spc_<discovery_id>`、`collapsed_at_dft_from_seed`)は、TS のない辺の source が自分の species を持つので不要になった。
43. (W4)古い manifest を読むためだけの欄(`ReactionTrial.source_minimum`、`DiscoveryRecord.source_minimum`・`ts_imag_cm1`・`electronic_temperature_K`)は W7-F3 で、鍵の移行と VAL10 の後に削除する。ReaDuct の pin は 3 つの版の正準な文字列で、辞書にするのは W7-F1。
44. (W4)RDKit を本体の依存にした(状態ラベルはどの pipeline も使う)。structures の rdkit の ImportError の分岐を削除した。RDKit の版は環境の一部として固定し、変えたら再生の基準を作り直す。W4 より前の run dir は最初の stage から流し直す(ラベルの書式が変わった)。
45. (W4)受け入れ条件のうち「trials のコード行が増えない」「explore.py が短くなる」は満たさなかった(trials 124 → 165: Lewis の MILP、正確な類のキー、複数配座の実現と剛体の配置。explore 159 → 234: 世代のループ、到達と向き、ボトルネックの順、TS の重複除去、予算の記録)。gates・hypotheses・worker・topology・identity は短くなった。
46. (W4)未緩和の seed の SP を出発の点にすると DFT の入口の高さが甘くなる(S5 の二重項の seed は CH3 + O2 の漸近より約 19 kcal/mol 上)。記録だけで、案(宣言した単量体の seed は単量体の SP の和から測る)は W6-1 の登録簿の出発の点と一緒に扱う。

47. (W4 統合)新しい run で見つけた誤り 2 つを直した(回帰テストあり)。(a) S6 は同じ原子の組成を 2 つ宣言するので、1 つの状態に原子順の違う構造が入り、片方の元素の並びでもう片方を読んでいた(同じ trial_id の重複、S6 の 276 件中 39 組)。explore は原子順ごとに列挙して類でまとめる。(b) DFT の入口の到達を species で判定していたので、別の配座から出る第 2 世代の辺が unreached になった(fluoroethane 51 本中 12 本)。節点を状態にし、状態のエネルギーはその端の構造の SP の最小にした。
48. (W4 統合)DFT の入口の順は (高さ, 道のエネルギーの幅, 道の辺の数) にした。S5 では未緩和の seed が全体より高く、すべての辺が高さ 0 になり、id の順で会合の辺が `over_cap` になった。窓の規則(出発の状態からの高さ)は変えていない。
49. (W4 統合)効率: (a) DFT の入口は、出発の状態から高さの低い順にたどり、窓の内側で届いた状態を出る辺にだけ SP をとる(窓の判定は全部とるのと同じ。SP は S10 で辺 850 本に 71 点、S19 で 975 本に 54 点、fluoroethane で 51 本に 36 点(全部とると 84 点)。S5 は辺がすべて seed より下なので 220 本に 257 点)。(b) explore は 1 世代の状態をボトルネックの低い順に、予算を満たすまでだけ列挙する。(c) `identity.assign` は、重心からの距離の列の差(RMSD の下限)で候補を先に除く(結果は同じ。S10 の TS 250 個の重複判定が 34 s → 0.3 s。以前は explore の時間の大半がこの判定だった)。
50. (W4 統合)cases.yaml: S6 は辺ごとの case になった(引き抜き [9.95, 10.45] はそのまま、CH3 + H2O → CH4O + H elementary 54.62、CH4O + H → CH3O + H2 elementary 13.10、CH2 + H + H2O → CH4O + H barrierless。VAL9 の多段の親と子 2 つは消えた)。golden G36 を W4 の run から作り直した。S7・S8 の生成物の数は 26・1027(世代を重ねた辺)。S5 は変えない: W4 の run は split2 が unresolved(W3 と同じ境界の揺れ)で、入口が入れた辺から行が増えた。2 回流す VAL10 で決める。fluoroethane の 4 中心の HF 脱離の TS は見つかったが、DFT//xTB の高さ 69.5 kcal/mol で窓 40 の外(事前に宣言した衝突。窓は変えない)。

51. (W5-1)maxiter で止まった探索の継続(最後のフレームの DFT freq から)は、同じ失敗の印の中で行う。新しい種として数えると 1 回の停滞ごとに印を 1 つ使い、W4 の S6 の rxn_discovery_d403cb0271(maxiter → χ の棄却 → barrierless の string)が既定の予算 2 で unresolved になる。高次の鞍点と認証に落ちた TS の継続は、次の種として数える。
52. (W5-1)FIND_PATH が印を使うのは、失敗か unavailable の string だけ。barrierless の string は case を閉じ、intermediate は行 6、single の種は精密化で数える。
53. (W5-1)TS の継続は鞍点 x から始める(その勾配が P-RFO を導く)。\|H\| の信頼領域の歩は反応モードに沿って下るので、x + s は極小にだけ使う。停留した高次の鞍点だけ off_saddle で押す。継続(段数 1)がなお認証に落ちる一次の鞍点は、棄却せず `not_stationary` を付けて受理する(thermo の blocker)。継続の高次の鞍点は棄却する。
54. (W5-1)判断表は約 12 行ではなく 9 行になった(静的な行 1〜3 は `hypotheses.select` の閉じた記録、接続の 3 値の行は 1 つの止める規則に吸収)。drivers のコード行は −50 ではなく +17(1090 → 1107: 下りを共有する Relaxer/Settled、`BarrierVerdict.points`、外への到達)。
55. (W5-1)窓の外の未宣言の候補は、黙って落とさず OUT_OF_WINDOW の閉じた記録にする(coverage に出る)。分割の子は select を通らないので窓をかけない(親の窓が効く)。
56. (W5-1)NWChem の継続が使う movecs と drv.hess は、Fortran の記録の印が全体を覆うもの(途中で切れていない)だけ。golden の STUB は記録の形のファイルを書く。
57. (W5-2)libmsym は VAL9 の S6 の TS にどのしきい値でも候補を出さなかったので、候補は置換の後戻り探索にし、libmsym は受理した群の名前だけ。開殻の縮退した群の例外は identity の雑音の基準(BASIN_*)を保つ。`thermo.electronic`(G_el)と `thermo.spin_orbit`(E_SO)は 2 つの関数で、原子の表は `electronic_state` にある。
58. (W5-2)立体は断片ごとに正準の原子順で電荷 0 として知覚し、`state_label` は電荷を取らない(全呼び出しと分離した状態の電荷が要るため)。規則: 電荷のある原子への二重結合に E/Z なし、四面体中心 1 つは立体なし、TB/OH/SP の印は捨てる。ラベルが変わるのは立体を持つ化学種だけ(DFT 極小と入力 452 構造で分割 0、Z-エノールの 18 件)。再生の基準は environment.md の規則で作り直す。
59. (W5-3)rc 1 は、止まった stage と error の項目に加え、`VOLATILE` のジョブの失敗(timeout、executable_missing、out_of_memory)と `running` の stage も含む(再実行に仕事が残る。その stage は fresh でない)。FailureKind は 11 種(+out_of_memory、+error)。パーサの例外は記録した incomplete_output ではなく項目の error。例外を出した反応ケースは登録を巻き戻して failed の記録にする。minima の登録は閉じ込めず、履歴は `opt:reused` から始まる。sp の artifact は (method, freq) ごとに 1 つ(失敗も)。結果を残さなかった attempt の継続は LADDER の回数に数えない。
60. (W5-4)BH76 の表には W4 のラベルの run が要るので、hcn・sn2_cl・sn2_cl_d3h・h3_doublet の known_endpoints を新しく流した。3 つの汎関数は原子 guess からの 1 つの SP の手順で比べ(run の層との差 0.001 kcal/mol)、規則の (1) は類ごとに当てた。M06-2X を残す(revM06: SCF の健全さと費用 1.12 倍 ↔ HTBH38 の MUE 1.20 → 1.91)。
61. (W5 統合)thermo は自分の点の注記 `not_stationary` を同名の blocker にする(W5-1 の要求)。explore の試行と conformers の探索も `StageRuntime.contain` に通す(パーサの例外が stage 全体を失敗にしないように)。site の scratch のパスは読み込みで 1 回だけ展開し、NWChem の timeout の既定を 86,400 s(ハングの検出)にした。使われなくなった `modes.qrc_step` の `factor` を削除した。
62. (W5 統合)cases.yaml の δG_eff を、W5 で測った値に理由とともに変えた(QH と電子準位による熱化学のモデルの変化。同じ構造で ΔE は不変): s1_sn2_cl・sn2_cl_d3h 14.53 → 15.09(錯体の G が QH で −0.61)、dme_c2v_seed 2.13 → 2.15、S6 の引き抜きの幅 [9.95, 10.45] → [10.10, 10.60] (同じ W4 の構造で +0.16: QH −0.29、OH の ²Π 準位 +0.44)、S6 の 2 つの辺 54.62 → 54.57・13.1 → 13.73。測っていない行(s2_sn2_i、s20_sn2_panel、oxalic、S10 など、B1 での見積もり \|Δ\| ≤ 0.02 か +0.45)は VAL10 の値でだけ変える。
63. (W5 統合)極小の並列(U9-P2(a))は s19 の minima(dft) で速くならなかった(opt 4 本を 1 rank ずつ 1,010 s、freq 3 本 1,270 s。4 rank の直列の見積もりより約 14% 長い。round 7 と同じ向き)。core 秒は半分なので、コードは並列のまま VAL10 の時間の表で直列と比べて決める(W7 の検査に足した)。
64. (W5 統合)acac の TS は継続の後も認証に落ち(ΔE_TR 2.05e-4 Eh、平らなメチルのねじれに沿う傾きが NWChem の saddle の既定の収束の内側に残る)、事前の宣言どおり `not_stationary` で順位から外れた。W5 の受け入れ条件「継続の後 ΔE_TR ≤ 5e-5」「QRC の core 秒が減る」は未達。認証の基準は動かさない。saddle の収束を締めるかは、測ってから決める(見送りに足した)。

## 3. 見送り

仕組みを増やすもの、先に測る必要があるもの: GFN-FF の組成探索(U2-P1。DFT の basin の集合でだけ再評価)、Boltzmann 和(X4-4)、TS の配座(U6-P10)、1D-HR(U7-P8)、CVTST(U7-P10)、速度の量(U8-P11)、基底の規則(U0-P7)、case の並列化と slurm(U9-P2 b・c)、有限差分の Hessian(X7-4)、ZPE のスケール(U7-P7)、WFT のメモリ表(U8-P12)、xTB の PES の仕様(U0-P6)、電子状態の選択(U0-P8)、U0-P9、U1-P6、U1-P7、U2-P3、U2-P5、U2-P6、U3-P5、U6-P8、tier の観測 Level への置き換え(U9-P8)、純粋な開裂の逆向きのスキャン(U4-P5)、断片の DFT 極小の自動化(X4-3 の第 2 段)、U9-P11 と core 秒の台帳、saddle の収束を締めること(acac の平らなモードの TS。変更 64)。範囲外はレビュー §0.5(MECP、多参照、溶媒、VRC-TST とマスター方程式、多次元トンネル、遷移金属の検証、ビット単位の再現性、GPU、ML)。


## 整理前の設計書の未統合案

レビュー [2026-10-04_platform_review_2.md](../reviews/2026-10-04_platform_review_2.md) の計画のうち、段階 0 のプローブ(証拠は `/home/user/hfauto_r10/probe/P0a`〜`P0e`)で採否が決まった設計である。実装は wave の順に進め、統合のたびに該当する §0〜§10 を書き換えて本節から外す。

### 11.1 鍵を変える基盤の変更(W7)

- QRC の側の opt も TS の freq の解から始める(`Relaxer.descend` の opt に `scf_guess=freq`)。paths の全 QRC のジョブ鍵が変わる。
- ジョブ鍵を変える基盤の変更(deck の sha、版数の辞書、`canonical_json` の一本化、密度の収束、`levels`)と合わせ、最後に 1 回でまとめて移行する。

### 11.2 化学の判定

- 貸された TS(他の候補の低レベル TS、mode-follow の `ts_calc`)の χ は、貸し手のラベル付きの両端で測る。貸し手の両端は今の記録にないので、W6 の登録簿の辺に自分の両端を持たせ、そのとき種に結合変化を持たせる。それまでは case の両端の変化 `Ctx.change()` で測る(質量加重の χ は W2 で §3 に取り込んだ。等価な原子の Hungarian の付け替えは測って棄却した: [validation.md](../validation.md) §11.2)。

### 11.3 順位の量と計算点(残り)

W3 で反応が読む点、点の鎖と共通のゼロ、宣言した単量体の分離した点、捕獲律速の表、`profile.judge` を §3・§7 に取り込んだ。残り:
- 捕獲律速の表に、ランキングの LOT で経路の節点をたどった判定(`layer_verdict`、`dE_act_path_kcal`)を足す(X5-3。記録の欄は W5 で足した)。
- 事実(状態の組の結論)は stage の中の登録簿 1 つに置き、種は 1 本の順序付きの流れで渡す(W6)。

### 11.4 配座と錯体(P0b)

- 組成の探索は GFN2 の CREST のまま続ける(`--nci --quick --notopo <全原子> --noopt`)。GFN-FF は採らなかった。s10 の 2 つの組成では、入力状態の候補が 1 つずつしか出なかった(GFN2 も同じ)。S5・S6 は rc 0 になったが、出た構造は断片がほぼ離れていた(最短接触の中央値 5.5 Å、GFN-FF のエネルギーの幅 0.7 kcal/mol 未満)。再検討するときは、ラベルではなく DFT の basin の集合で判定する。
- seed を組成ごとに 1 つ(重心を結ぶ線上の接触)にする案(U2-P4)は採らなかった。前提の U2-P1(GFN-FF)が落ち、接触の seed 1 つでは S5・S6 の GFN2 の screen がどちらも失敗して、両系の反応がすべて消える(W3-3)。錯体が screen の basin を経ずに explore・DFT に届く経路(U2-I1)ができたら再検討する。lone-pair 円錐と元素の表は W3 で削除した。
