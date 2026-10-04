# round 10 のロードマップ

レビュー [2026-10-04_platform_review_2.md](reviews/2026-10-04_platform_review_2.md) §7 の計画を、段階 0 のプローブ(P0a〜P0e。証拠は `/home/user/hfauto_r10/probe/`)の結果で決めた形で進める。目標の設計は [design.md](design.md) §11 で、wave を統合するたびに §0〜§10 に取り込む。実計算の記録は [validation.md](validation.md) §11。

## 1. 提案と wave

| wave(段階) | 提案 | 状態 |
|---|---|---|
| W1(1、must) | X7-1 予算は件数だけ、X7-2 の一部(外からの停止は incomplete・rc 1)、U3-P1 σ を対称操作から、X4-4 の must 部分(状態の G を閉じる)、U1-P1 多重度の宣言を必須に、X7-7 検証を repo に(cases.yaml、BH76 の参照、compare、replay) | 完了 |
| W2(1、must) | U6-P1 χ の質量加重、U5-P3a、U6-P2、U6-P6、U6-P9、X3 SCF の来歴と救済の 3 段(P0a)、開殻 CCSD(T) は UHF 参照(P0d)。最後に再生の基準 B1 を作る | 未着手 |
| W3(2、must) | X5-1 反応が読む点、X4-3 鎖と共通のゼロ(BH76 の分離のゼロ用に h3_doublet と oh_ch4 の単量体を宣言)、U5-P2 捕獲律速、U5-P1、U5-P7、U2-P2・U2-P4 | 未着手 |
| W4(3、must) | U1-P3・U1-P2 正準のラベル、U4-P1 結合グラフの編集の列挙(電荷分離の規則は P0f)、U4-P2・U4-P4・U4-P7、U4-P3 DFT//xTB の窓、U3-P7 | 未着手 |
| W5(4、should) | X1-2 停留性の認証、U3-P4、U6-P4、U5-P5、U5-P6、X2-2・X2-3 対称性の受理と電子準位、U7-P6、U1-P4、X7-2 の残り(rc 2、項目の閉じ込め)、X7-3、U9-P3、U9-P2(a)、エネルギー層の選び直し(BH76) | 未着手 |
| W6(4、should) | X4-2 事実の登録簿、U5-P4、U8-P10・U8-P7・U8-P6、X5-3 層の経路の判定 | 未着手 |
| W7(最後) | 鍵を変える基盤の変更を 1 回で移行(X7-6、U9-P10、X1-3、U0-P3)、分割の仕組みの削除(X4-2 段階 2)、整理、VAL10 と BH76 の表 | 未着手 |

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

## 3. 見送り

仕組みを増やすもの、先に測る必要があるもの: GFN-FF の組成探索(U2-P1。DFT の basin の集合でだけ再評価)、Boltzmann 和(X4-4)、TS の配座(U6-P10)、1D-HR(U7-P8)、CVTST(U7-P10)、速度の量(U8-P11)、基底の規則(U0-P7)、case の並列化と slurm(U9-P2 b・c)、有限差分の Hessian(X7-4)、ZPE のスケール(U7-P7)、WFT のメモリ表(U8-P12)、xTB の PES の仕様(U0-P6)、電子状態の選択(U0-P8)、U0-P9、U1-P6、U1-P7、U2-P3、U2-P5、U2-P6、U3-P5、U6-P8、tier の観測 Level への置き換え(U9-P8)、純粋な開裂の逆向きのスキャン(U4-P5)、断片の DFT 極小の自動化(X4-3 の第 2 段)、U9-P11 と core 秒の台帳。範囲外はレビュー §0.5(MECP、多参照、溶媒、VRC-TST とマスター方程式、多次元トンネル、遷移金属の検証、ビット単位の再現性、GPU、ML)。
