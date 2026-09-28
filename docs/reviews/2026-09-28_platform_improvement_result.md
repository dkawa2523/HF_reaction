# hfauto 計算基盤の抜本的改良 — 結果報告（2026-09-28）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `d5fd24c`。基準はレビュー `docs/reviews/2026-09-27_platform_review.md`（以下「レビュー」、コード `fbf7924`）。
- 範囲: must 段 M0〜M5（`ebafbcf`〜`82c46d4`）、should・could 段 S-A〜S-D と FINAL（`6d41726`〜`496fad0`）、§6 検証セットの再実行（`d5fd24c`）。push はしていない。
- 方法: 区分 U0〜U9 ごとに、コード・テストと WSL の実 run（`/home/user/hfauto_r6/M0`〜`M5`、`SA`〜`SD`、`FINAL`、`VAL`）で評価し、別の担当がコードと生データで反証的に確かめた。**印・状態・重大度は、確かめた段で合意したものを使う。**
- このファイルのほかに、リポジトリは変えていない。

凡例

| 記号 | 意味 |
|---|---|
| ◎ / ○ / △ / × | 問題なし / 軽微な課題のみ / 改良が必要 / 目的どおりに機能していない |
| critical / high / medium / low | 課題の重大度 |
| must / should / could | レビューの優先度 |
| done / partial / not_done / rejected | 実施 / 一部 / 未実施 / 不採用（プローブかレビューで） |

エネルギーの単位は、断りがなければ kcal/mol。「実計算」は WSL で NWChem・xTB などを実際に流した run を指し、「QM なし」は fake エンジンや固定データのテストを指す。

---

## 0. 要旨

### 0.1 改良の総括

- **レビューの重大課題 8 件は、コード上すべて手当てした。** 実計算で解消を確かめたのは 1〜5 と 8 である。
  - 1 ECP: I・Te で def2-ECP が自動で入る。
  - 2 単原子と失敗の閉じ込め: Cl⁻ が同じ経路を通り、run は必ず report まで進む。
  - 3 縮退と鏡像: basin 構造と掌性で判定するようになり、偽反応が消えた。
  - 4 反応方向: acac は 36 歩で未収束だったのが 6 歩で収束した。
  - 5 QRC ゲート: H2O2 のねじれ TS が接続するようになった。
  - 8 順位の量: δG_eff 1 つにした。
- 残るのは 6 と 7 である。6 は分類器を 1 つにしたが、障壁なしと中間体の出口を実計算で一度も通っていない。7 は恒等 SN2 と縮退した H 交換は発見になるが、GFN2 で障壁のない引き抜きは発見にならない。
- **削除:** 次の重複や足場をなくした。
  - knob 23 → 8。元素表 4 → 1。経路の形の判定器 2 → 1。blocker の生成元 3 → 1。SP と停留点の対応付け 3 → 1。重複除去 3 段 → 2 段。
  - GoodVibes の外部エンジン（worker・engine 234 行、JobStore、実行時照合）、GS アダプタの座標系分岐 5 つ、溶媒和・D4・MP2 の半実装、符号不一致と dzpe の経験則 blocker、自動の配座対仮説、読まれない記録のフィールド。
- **コード量:** `hfauto/` は 10,635 → 10,336 行（−299、−2.8%）。M5 から数えると −452 行。レビューの見込み（正味 −700）には届かなかった。主な理由は、失敗の閉じ込め、シグナル処理、ROHF-CCSD(T) のデッキ、生成物の同定、mode-follow の端点登録を足したこと。
- **検証:** R1〜R4 はすべて合格した。S 系は 20 件中 17 件が合格、S6 が条件付き合格、S15 と S17 は期待した分岐に入らなかった（どちらも結果は化学的に正しい）。検証中に見つけた対称数 σ の不具合（S2）は直した。
- **プローブで決めたこと:** 固定端 CI-NEB は採用した（U5-P2）。composite（ωB97X-D3 のエネルギー層）は既定にしなかった（U8-P3。17 原子で 1 点 29 分かかるため）。剛体配置は不採用とした（U2-P6）。SCF の救済は cgmin の 1 本にした（U0-P5）。

### 0.2 評価の変化（全体）

| 観点 | 改良前（レビュー §0.1） | 改良後 | 一言 |
|---|---|---|---|
| 化学的妥当性 | △ | ○ | 黙った誤り（ECP 欠落、鏡像の偽反応、正しい TS の棄却）はなくなった。残るのは既定 LOT の系統誤差（CCSD(T) に対し最大約 3）と σ の揺らぎ |
| 汎用性 | × | ○（条件付き） | 元素表 1 つ、ECP 自動、単原子、スピンの宣言。DFT まで通したのは H・C・N・O・F・Cl・I・Te、電荷 0/−1、多重度 1〜3 だけ |
| 有用性 | △ | ○ | 順位の量が 1 つになり、障壁なし・沈んだ障壁・分離基準も同じ量で並ぶ。探索（U4）は △ のまま |
| 効率 | ○ | ○ | 読まれない DFT freq をやめ、方向を正し、explore を並列化した。一方で QRC の長い尾（15 原子）と mode-follow の重複は残る |
| 簡潔さ | △ | ○ | 規則・定義を 1 つずつにした。ただし行数は −2.8% にとどまり、U3・U4 は横ばいか増加 |

### 0.3 残る課題（重大度 medium 以上）

1. **実計算で通っていない分岐がある**（U5・U6・U9）。
   - 決定表の行 1・3・5〜8・11・13・15〜17。
   - multi_step と split、BARRIERLESS、FIND_PATH（smoke の 1 チャンクを除く）、saddle_restart、失敗の閉じ込めの発動。
2. **既定の順位は PBE0-D3BJ/def2-SVPD のまま**（U0・U8）。CCSD(T) との差は SN2 −3.0、HONO +2.0。エネルギー層は method_panel を手で追記して使う。
3. **二重項（ラジカル）の組成と Na⁺·H2O では CREST 3.0.2 が必ず失敗する**（U2）。未緩和の seed だけが screen に渡る。
4. **対称数 σ が小さな構造の揺らぎで変わりうる**（U7）。実 run で影響したのは 126 本中 1 本（S2）で、pymsym の版も固定していない。
5. **mode-follow から、宣言していない配座変化が仮説になり、順位 1 位を取る**（U5）。S4 で起きた。
6. **mode-follow の効率**（U3）。変位元の Hessian を渡すと唯一の実例（S4）で遅くなった。両側の freq も重複して計算している。
7. **高次 saddle の再試行が停滞すると、再開が試行予算の外になる**（U6）。推論では、DFT string を 1 本余分に払うことになる。
8. **GFN2 で障壁のない反応は、緩和の陰性としてだけ残る**（U4）。S6 の引き抜き、S5 の CH3+O2 がこれに当たる。

---

## 1. 実施した改良

### 1.1 must 段 M0〜M5（S 段の前にコミット済み）

| 段 | コミット | 提案 | 主な内容 | 主な削除 | hfauto 行数の増減 |
|---|---|---|---|---|---|
| M0 | `ebafbcf` | U9-P1、§6 | 失敗を 3 か所で閉じ込める（parse → `incomplete_output`、ケース → UNRESOLVED、入力なしの stage → 0 artifact）。§6 の系の設定と golden を追加 | 入力なしの stage の ValueError | +61 |
| M1 | `9400787` | U1-P1、U0-P1、U3-P1、U1-P3 | 元素表を 1 つにした（Z=1–57、72–86）。def2-ECP を元素ごとに自動で書く。単原子を同じ経路に通す。多重度を宣言制にした | 元素表 4 つ、最高スピンの選択、冗長な再検査 | +141 |
| M2 | `a80b01e` | U3-P2、U5-P6、U5-P1、U5-P3、U7 の m | 鏡映を含む basin の同一性と m=2。縮退を basin 構造と掌性で判定。分類器 1 つ | `same_minimum`、profile.shape、below_zpe と SCREEN の DFT freq、BarrierVerdict の 4 フィールド | +23 |
| M3 | `b530f1d` | U6-P2、U6-P1、U6-P5 | 反応方向の規則を 1 つにした。Hessian を「反応モードだけ負」に整える。TS・QRC ゲートを定義だけにした | moddir の 3 分岐、low_prominence、no_initial_descent、trajectory_above_ts、torsion_saddle_cm1 | −1 |
| M4 | `14b700c` | U4-P2、U4-P1、U4-P5 | 元素に依らない結合変化の列挙器 1 つ。drive の同値類。縮退転位を発見として扱う | `_polar_h`、`_h_shifts`、芳香環検出と経験則の定数 | +6 |
| M5 | `82c46d4` | U7-P1、U8-P2 | δG_eff = max(G_TS, G_R, G_P) − G_R。障壁なしと沈んだ障壁も順位に入る。パネルの幅を列にした | dzpe の連鎖、符号 blocker、population、metric、0 に強制する分岐 | −77 |

### 1.2 S-A（`6d41726`）精度と報告

| 群 | 提案 | 内容 |
|---|---|---|
| SA-A | U8-P3（修正版）、U8-P6、U8-P8 | エネルギー層は method_panel の追記で使う（panel_sp → panel_thermo(energy_method) → panel_report）。条件は ThermoConfig へ移し、`energy_level` を ranking.csv の列にした。load_method は pipeline の隣の methods/ を先に探す。SP の対応付けは sp artifact の parents だけで行う。blocker の生成元は `thermo._blockers` 1 つ。SP 層にも spin_ok をかける |
| SA-B | U8-P4 | SP を計算するのは順位に使う点だけ（rankable な反応の TS、反応物・生成物の状態、単量体の状態） |
| SA-C | U8-P7、U0-P5 | 開殻の CCSD(T) は ROHF + TCE（`2eorb 2emet 13`）。SCF は反復 100、救済は cgmin の 1 本（プローブで決めた） |

- 削除: `core.system.Conditions`、`PipelineConfig.conditions`、`SinglePointConfig.targets`、指紋による突き合わせ、rankable の thermo_unavailable の再判定、`wft_closed_shell_only`、damping・level shift の救済、重複していた単量体の対応表。
- 行数: +30（TCE のデッキと SP の対象規則を足したため）。
- プローブ: ωB97X-D3/TZVPD の CCSD(T) に対する平均絶対誤差は 1.29（PBE0/SVPD は 2.07）。ただし 17 原子で 1 点 1,758 s かかるので、既定にはしなかった。

### 1.3 S-B（`36f5cae`）経路と鞍点の効率

| 群 | 提案 | 内容 |
|---|---|---|
| SB-A | U5-P2（プローブで採用）、U5-P4、U6-P6、U6-P3 | SCREEN は、DFT 極小を固定端とする IDPP → xTB の Cartesian CI-NEB 1 本にした。FIND_PATH は 1 回に 1 チャンクだけ回し、すぐ分類する。maxiter で止まった saddle は継続せず、最終フレームから 1 回だけやり直す。検証済み TS の freq を再試行の Hessian に使う |
| SB-B | U3-P3、U3-P4、U3-P6、U3-P8、U6-P4 | 変位はエネルギー目標の `modes.amplitude` 1 つにした（QRC、高次 saddle の押し、mode-follow、ReaDuct で共用）。鞍点に収束した宣言端点は近い側の basin に入れる。rerank の後にも窓をかける。MinimumPolicy を削除 |
| SB-C | U3-P7、U6-P7、U5-P7 | freq と親のエネルギーの一致を見る（state_mismatch）。縮退 case では結合の交換を確かめる（bond_change_missing）。自動の配座対仮説をやめた |

- 削除: GS アダプタ（`_axis`・`_linear`・off_axis・`_split`・tric・DLC の分岐と GS の再実行）、`refine_ts`、`_COLLAPSE_A`、`_SETTLED`、`_RETRY_A`、`energies_settled`、PathProfile の履歴、string の途中エネルギー、saddle の maxiter 継続、`qrc_amplitude`、配座対の分岐、古い決定表の行。
- 行数: −38。pysis の engine と worker は合わせて 295 → 224 行になった。
- プローブ（11 系）: NEB は全系で 5〜20 サイクル、2〜4 s で収束した。(HF)₂ は未解決から C2h の TS に届き、acac は完走し、HONO と H2O2 は FIND_PATH に入らなくなった。

### 1.4 S-C（`0f1d868`）探索と配座

| 群 | 提案 | 内容 |
|---|---|---|
| SC-A | U4-P4、U4-P6、U4-P11、U4-P7、U4-P9、U9-P6 | IRC と端の opt を 500 反復にし、負モードが残る端は ± に変位する。窓は `gates.discovery_verdict` 1 つ（kcal/mol）。AFIR は NT2 が極大を見つけないときだけ走らせる。生成物は 1 回だけ同定する。explore を thread_map で並列にした |
| SC-B | U2-P2（プローブ）、U2-P4、U2-P6（プローブで不採用）、U1-P4、U1-P7 | 組成の CREST は全原子を `--notopo` に渡し、`--noopt` を付ける 1 規則にした。conformers 独自の重複判定と診断を削除した。結合判定は加算型（Σr_cov + 0.4 Å）で、bond_changes は集合の差。topology は結合とラベルだけ（239 → 112 行） |

- 削除: `explore.Window`、`trials.product_verdict`、`ReactionTrial.perturbed`、`nt_total_force_norm`、`threads_per_item`、`topology_removed`、結合判定のヒステリシスと previous、`_notopo`、`_same`、`DUPLICATE_*`、`_REMOVALS`、conformers の diagnostics.json。
- 行数: −19。trials.py は 254 → 176 行になった。
- U2-P6 を不採用にした根拠: 剛体配置にすると、NH3·HF⁺ と OH·H2O で錯体の状態が消えた。

### 1.5 S-D（`ac8f802`）基盤と簡素化

| 群 | 提案 | 内容 |
|---|---|---|
| SD-A | U9-P2、U9-P5、U0-P4 | knob を 23 → 8 にした（Policy 5 + ReactionPathsPolicy 3）。予算は 1 仮説ごとに閉じ、分割した子は親の締め切りを継ぐ。証拠による完了（行 4〜6）を walltime の行 7 より前に置いた。気相専用とし、溶媒和・D4・MP2 を削除、二重混成は INPUT_INVALID |
| SD-B | U7-P5 | GoodVibes 4.3.0 を同一プロセスの純関数で呼ぶ。版数は preflight が照合する |
| SD-C | U9-P7、U1-P9 | SIGTERM・SIGHUP で、登録した外部プロセスのグループを kill する。対象範囲（元素、スピン、SO、SIE、イオン–双極子）を文書にした |

- 削除: `backends/goodvibes`（engine 126 行 + worker 108 行）、ThermoEngine、ThermoResult、`Capability.THERMO`、`thermo_consistent`、CasePolicy、`_POLICY_KEYS`、`case_policy`、`ExecutionSpec.maxiter`・`coordinates`、`ThermoSettings.symmetry`、`MethodSpec.solvation`・D4・MP2、COSMO・ALPB の分岐。
- 行数: −280。この改良で最大の減少。

### 1.6 FINAL（`496fad0`）定義の一本化と文書

| 群 | 提案 | 内容 |
|---|---|---|
| FINAL-A | U9-P4、U9-P3、U9-P9 | 1 概念 1 定義にした（StandardState・ConnectionLabel・CONNECTED_OUTCOMES は records、幾何のプリミティブは geometry、原子的書き込みは core/files、FileRef は JobStore）。freq の不変条件は Evidence で 1 回だけ保証する。読まれないフィールドを削除した。run 内のファイルは Runtime.resolve が解決する |
| FINAL-B | （提案なし） | README・design.md・environment.md を簡潔に書き直した。§6 の検証セットを 1 回流した |

- 削除: `Evidence.trajectory_energies_hartree` とその生成元、PathProfile の `gmax_history`・`program_converged`・`climbing_image`、`ConnectionClaim.method`・`notes`、JobCounts、engines の `_OVERRIDES`、`gates._mode_count_reasons`、`xyz_trajectory.py`（xyz.py に統合）、`HARTREE_TO_KJ_MOL`、古い `__version__`、G02 の CSV golden。
- 重複していたコピーも 1 つにまとめた（written_geometry 7 つ、species_artifact_id 6 つ）。
- 行数: −170。

### 1.7 検証（`d5fd24c`）

- §6.1 の検証セットを、新しい run dir（`/home/user/hfauto_r6/VAL/runs/`）で 1 本ずつ直列に流し直した（§3）。
- 対称数の修正: `chemistry/thermo.symmetry_number` は libmsym の等価しきい値を 2e-3 にして σ を求め、GoodVibes には `symm=False` で渡す。回帰テストも足した。行数は +25。
- `docs/validation.md` を現在の記録として 542 → 129 行に書き直した。段ごとの記録は git の履歴に残る。

### 1.8 実施しなかった提案・一部だけの提案

| 提案 | 状態 | 理由・現状 |
|---|---|---|
| U0-P3 → U8-P3 composite を順位の正式な経路に | done（修正版） | 既定にすることはプローブ（費用）で見送り、記録した。既定の順位が停留点レベルのままなのは設計上の選択（§4 U0） |
| U0-P9 親の movecs から開始 | rejected（レビューで見送り） | SCF 解を左右する入力を鍵の外で渡すことになる。代わりに U3-P7 を入れた |
| U0-P13 def2-SVP の method | partial | スケール上限の記載だけ。実際の要求が出るまで作らない |
| U1-P5 WL を安定するまで回す / U1-P6 spin_tol を相対値に | partial / not_done | could。衝突や誤判定の実例がない |
| U2-P6 剛体配置 1 つ | rejected（プローブ） | CREST が失敗する系で錯体の状態が消える |
| U3-P3 mode-follow の初期 Hessian | done（効果は否定。差し戻し候補） | S4 では Hessian なしで 31/31 歩、あり（trust 0.3）で 44/37 歩 |
| U3-P5 崩壊した鞍点の freq を再利用 | not_done | 対象の経路が実計算で一度も通っていない |
| U3-P11 linopt 0 | not_done（プローブ未実施） | — |
| U4-P10 → §6 の explore 検証 | partial | CH3O• の T1（1,2-H 移動）と T4 の実行経路を実計算で流していない |
| U5-P7 自動の配座対仮説をやめる | partial | conformer の分岐は消えた。しかし `_degenerate` と `_auto` の 0.2 Å 規則を通って、mode-follow の配座変化が今も入る（§4 U5） |
| U5-P4 FIND_PATH 1 チャンク | done（実計算未検証） | 決定表から配線した FIND_PATH は、実計算で一度も走っていない |
| U5-P9 経路 SP を 1 ジョブに / U5-P11 IDPP の衝突判定 | not_done | 節約が小さい（HCN で 16 s）。DFT 極小どうしでは衝突が起きない |
| U6-P11 RKS 不安定性の検出 | not_done | 開殻一重項は対象外と明記した |
| U7-P6 Eckart κ | not_done | トンネル補正は対象外と明記した |
| U8-P11 HTML の条件 / U8-P13 M06-2X | not_done | §5 に最小形を挙げた / 寄与が小さい |

---

## 2. 評価の変化（区分 × 観点、改良前 → 改良後）

改良前はレビュー §3 の評価、改良後は確かめた段で合意した評価。

| 区分 | 妥当性 | 汎用性 | 有用性 | 効率 | 簡潔さ |
|---|---|---|---|---|---|
| U0 全体プロトコルと理論レベル | △→○ | ×→○（条件付き） | △→○ | ○→○ | ○→◎ |
| U1 入力構造・電子状態 | △→○ | ×→○ | △→○ | ○→○ | ○→◎ |
| U2 配座探索・錯体配置 | ○→◎（条件付き） | △→○ | ○→○ | ◎→◎ | ○→◎ |
| U3 極小構造の確定 | ○→◎（但し書き） | △→○ | ○→◎ | ○→○ | ○→○ |
| U4 反応探索 | ○→○ | △→○（条件付き） | ×→△ | ○→◎ | △→○ |
| U5 仮説選択・障壁事前判定・経路 | △→○ | △→○ | ○→○（但し書き） | △→○ | △→○ |
| U6 鞍点精密化・TS 検証・接続 | △→○ | △→○ | △→○ | △→○（QRC は △） | ○→○ |
| U7 熱化学・速度論 | △→○ | △→○ | △→○ | ◎（計算）/○（仕組み）→◎ | △→○ |
| U8 一点計算・手法パネル・報告 | △→○ | △→○ | ○→○ | ○→◎ | ○→○ |
| U9 実行基盤・パイプライン | ○→◎ | △→○ | ○→◎（留保） | △→○ | △→○ |

評価の根拠と、確かめた段が付けた但し書き:

- **U0**: 唯一の critical だった ECP 欠落が直った。S2 は NWChem の全 16 ジョブに `I library def2-ecp` が入り α 30。VAL の NWChem 620 デッキで SCF 未収束は 0 件。汎用性の ○ は、ECP の実計算が I と Te だけで、陽イオンと遷移金属の DFT は未検証という条件付き。MethodSpec は 11 → 9 フィールド。continuation は CC 13 で上限に近い。
- **U1**: スピンを黙って選ばなくなった（S5・S9 は宣言がないと INPUT_INVALID）。ただし xyz で多重度を書かないと 1 になるので、偶数電子の三重項基底状態（O2、CH2）は黙って一重項になる。これはレビューで検証済みの既定で、README にも書いてある。topology は 239 → 112 行。
- **U2**: 旧規則では CREGEN が SO2·NMe3 の付加体を捨てていた（再現した: 付加体 0 件、最低でも +17.3 の vdW 構造だけが残る）。新規則では N···S 2.17 Å の付加体が残る。◎ の条件は、付加体が残るのは ewin のおかげで、状態ラベルは付加体と vdW 錯体を区別していないこと（U1 の閾値の責務）。構造数と時間は MTD の確率性で変わるので、「約 2 倍速い」とだけ書く。
- **U3**: 鏡像の偽反応、鞍点の端点で basin が空になる問題、固定変位がすべて解消した。ただし一部の分岐は実計算で通っていない。効率の ○ は据え置き。U3-P3 は唯一の実例で遅くなり、改善に数えない。4 ファイルは 756 → 815 行に増えた。
- **U4**: 有用性は × → △ にとどまる。VAL で DFT と順位まで届いた発見は S10 の縮退 H 交換 1 件だけ。汎用性の ○ は、生成物まで確かめたのが charge −1 の SN2、NH3·HF の縮退交換、二重項の CH3···H2O の 3 系だけという条件付き。4 ファイルは 787 → 782 行で横ばい。gates の discovery_verdict、elements の max_coordination 列、placement へ移した表を足すと、区分としては実質増えている。
- **U5**: barrierless は「DFT 極小を両端とする離散経路の最大値」で、上界になるのはノード解像度（`max_node_spacing_A`）の範囲まで。有用性の ○ の但し書きとして、S4 と S20 の順位表の先頭に、宣言していないねじれ（δG_eff 3.877、torsional 列あり）が載る。
- **U6**: 以前失敗した DME・H2O2・acac がすべて解けた。ただし S11 DME C2v（−226.8i）は actions を複写した harness の結果で、pipeline で DME が通ったのは S12 だけ。妥当性と汎用性の ○ は、VAL で通った一本道（screen → refine → validate_ts → connect）についての評価。
- **U7**: 順位の量は δG_eff 1 つになり、S14・S16 は `submerged_barrier` で 0.0。◎ にしない理由は 3 つある。σ がノイズに弱いこと、「状態の G」の定義が 2 つあること、BARRIERLESS と別の配座を基準にする場合を実 run で通していないこと。U7 の中核は 620 → 495 行になった。
- **U8**: 簡潔さは ◎ から ○ に下げた。(T, 標準状態) の選び方が report（ReportConfig か thermo[0]）と HTML（常に thermo[0]）の 2 か所にあり、html.render は CC 15 で上限ちょうど。summary.py は 300 → 239 行、rank_rows の CC は 15 → 12。
- **U9**: VAL の 44 run はすべて traceback なし、孤児プロセスも 0。有用性の ◎ は「SIGTERM と SIGHUP に限れば」という留保付き。SIGINT（Ctrl-C）は処理していない（§4 U9）。簡潔さは knob 23 → 8 だが、実行基盤のファイルは閉じ込めとシグナルの分だけ増えた。

---

## 3. 検証結果（§6 セット、`d5fd24c`）

条件:

- コードは `496fad0` に対称数の修正を加えたもの。環境は WSL2 Ubuntu（4 vCPU / 11 GB）、NWChem 7.2.3 を 4 rank × 1,200 MB。
- 手法は PBE0-D3BJ/def2-SVPD（Te・I は def2-ECP）と GFN2-xTB。298.15 K、1 atm で、順位の量は δG_eff。
- run はすべて新しい run dir で 1 本ずつ直列に流した。smoke から S16 までの合計は 4:11。打ち切った run はない。

### 3.1 合否表

| # | 系 | 期待 | 実測 | 時間 | 判定 |
|---|---|---|---|---|---|
| R1 | HCN → HNC | elementary、42.1783 | elementary（−1128.5i）、δG_eff 42.1783 | 1:17 | 合格 |
| R2 | HONO trans → cis | 約 11.82（m=2） | 11.8158 | 3:03 | 合格 |
| R3 | NH3 反転 | degenerate、3.8395 | 3.8395 | 1:15 | 合格 |
| R4 | 水 | same_basin | same_basin（順位なし） | 0:11 | 合格 |
| S1 | Cl⁻ + CH3Cl | degenerate、ΔE‡ 10.45、分離基準の値が出る | ΔE‡ 10.446、ΔG‡ 11.077、ΔG_assoc −5.08、分離基準 +6.00 | 7:30 | 合格 |
| S2 | I⁻ + CH3I | ECP、α 電子数、Level 一致 | 全ジョブで `I library def2-ecp`、α 30、ΔE‡ 6.553、ΔG‡ 6.882（修正後） | 13:19 | 合格（不具合を修正） |
| S3 | H2Te | 極小まで通る | `Te library def2-ecp`、E −269.20860 Eh | 0:17 | 合格 |
| S4 | CH3O• → CH2OH• | elementary、約 30.6、鏡像の偽反応なし | 30.612。CH2OH• のねじれの縮退 case（3.878、torsional）も加わる | 14:30 | 合格 |
| S5 | CH3• + O2 | 多重度未宣言は INPUT_INVALID、m=2 | 未宣言は `declare_multiplicity: candidates (2, 4)`。CREST に `--uhf 1`（rc −11）、CH3 の ⟨S²⟩ 0.7544 | 0:44 | 合格 |
| S6 | OH + CH4 | 引き抜きの生成物 | GFN2 で障壁がなく、緩和（collapsed_to CH3+H2O）としてだけ出る。生成物は CH3OH+H の 2 件 | 0:07 | 条件付き合格 |
| S7 | NH3···ICl | 2 断片、polar_h なし | 6 配座すべて 2 断片、polar_h の drive なし | 0:10 | 合格 |
| S8 | SO2·NMe3 | 付加体が残る | 付加体 2 配座（N···S 2.17 / 2.21 Å）が残る | 1:04 | 合格 |
| S9 | FeCl3·CH4 | 多重度未宣言は INPUT_INVALID | 未宣言は拒否。宣言 6 で CREST rc 0、6 配座 | 0:30 | 合格 |
| S10 | amine_pilot2 | 二重 H 交換が縮退 case になる | DFT −1081.3i、δG_eff 32.289、分離基準 26.63 | 22:24 | 合格 |
| S11 | DME C2v の種（harness） | 一次鞍点 約 226i | 2 次の saddle から押して −226.8i、ΔE‡ 2.33 | 7:03 | 合格（harness） |
| S12 | DME の粗い入力 | degenerate、約 2.3 | ΔE‡ 2.330 | 10:04 | 合格 |
| S13 | H2O2 | ΔE‡ 約 1.15 | 1.147（degenerate） | 1:38 | 合格 |
| S14 | マロンアルデヒド | submerged_barrier | ΔE‡ 2.011、ΔZPE‡ −2.402、δG_eff 0.0 | 16:52 | 合格 |
| S15 | HONO → HNO2 | multi_step（2 段） | 直接の 1,3-H 移動の TS（−2088.1i）で 1 段、52.016 | 2:46 | 期待の分岐に入らず（化学的には正しい） |
| S16 | acac（15 原子） | xTB Hessian で PT の TS | 重なり 1.00、DFT Hessian 0 本、−986.8i、δG_eff 0.0 | 1:30:57 | 合格 |
| S17 | (HF)₂ 入れ替え | barrierless | 障壁 1.256 が解像度 1.0 を超えるので障壁あり（C2h の TS −227.3i） | 1:35 | 期待の分岐に入らず（化学的には正しい） |
| S18 | H + H2 | TS が得られる | −605.7i、δG_eff 4.091 | 0:48 | 合格 |
| S19 | 複数の温度と標準状態 | δG_eff は標準状態で不変 | ΔG_assoc だけが 2RT ln(RT/P°) 動く。δG_eff は 1 atm と 1 M で同じ | 30:38 | 合格 |
| S20 | 手法パネル | CCSD(T) の行が埋まる | 3 系とも埋まり、min/max も method_panel.csv と一致。δG_eff は ωB97X-D3 の層 | 1:24〜11:14 | 合格 |

- 判定の内訳: R1〜R4 はすべて合格。S 系 20 件中 17 件が合格、S6 は条件付き合格、S15 と S17 は期待した分岐に入らなかった。
- 実エンジンの smoke は 13 件合格した（64 s）。golden（OH•、FHF⁻、HI の ECP、Cl⁻、ROHF-CCSD(T)、I⁻···CH3I の σ）は WSL の prod venv で合格した。Windows の QA venv では GoodVibes と pymsym がないので、σ の golden はスキップされる。

### 3.2 R1〜R4 の推移

| # | v5 / M0 | M2 以後 | VAL |
|---|---|---|---|
| R1 | 42.1783 | 42.1783 | 42.1783 |
| R2 | 12.2265（m=1） | 11.8158（m=2 を数えた） | 11.8158 |
| R3 | 3.8395 | 3.8395 | 3.8395 |
| R4 | same_basin | same_basin | same_basin |

どの段でも R1〜R4 は変わらなかった。R2 の差は TS のキラリティ（−RT ln 2）を数えたためで、12.2265 は撤回する（validation.md §6）。

### 3.3 見つけた不具合と修正（S2）

- **症状:** I⁻···CH3I 錯体（C3v、σ 3）が、この run では Cs（σ 1）と判定された。ΔG‡ は 7.533 になり、FINAL の 6.882 と RT ln 3 = 0.65 ずれた。
- **原因:** 最適化の終点で I⁻ が軸から 0.03° ずれただけで、libmsym の既定しきい値（5e-4）を外れた。
- **修正:** しきい値を 2e-3 にして σ を求め、GoodVibes に渡す。126 構造で試し、σ が変わったのはこの錯体だけだった。修正前の run は複製で thermo を取り直し、値が変わったのは S2 だけ（6.8817）。
- **位置づけ:** この方式は、レビューが見送りの欄で「頑健でない」と退けた案（pymsym のしきい値を緩める）そのもの。実 run で σ を誤る例が出たので採った逸脱である（§4 U7 の課題 1）。

### 3.4 到達表（実計算で通った分岐）

| 対象 | 通った | 通っていない（QM なしのテストだけ） |
|---|---|---|
| 決定表の行 | 2、4、9、10、12、14 | 1、3、5〜8、11、13、15〜17 |
| saddle の種 | screen_ts 14、screen_hei 2、discovery_ts 3、higher_order_retry（S11 の harness だけ） | path_hei、saddle_restart |
| outcome | elementary_step 6、degenerate_rearrangement 13、same_basin 3 | multi_step、barrierless、unresolved_within_budget、reassigned |
| エンジン | CREST、xTB、ReaDuct NT2/AFIR、pysis NEB、NWChem opt/freq/saddle/SP（PBE0、ωB97X-D3、CCSD(T)、ROHF-CCSD(T)） | NWChem string（smoke の 1 チャンクだけ） |

### 3.5 計算時間

| 系 | 以前 | 以後 | 主な変化 |
|---|---|---|---|
| S16 acac | 60 分で打ち切り（saddle 36 歩で未収束、レビュー時） | VAL 1:30:57（saddle 6 歩・227 s、QRC 2 側で 2,705 s） | 方向の規則と Hessian の整形で saddle が収束。QRC の尾が paths の約 65% |
| S10 amine_pilot2 | 23:35（M4） | VAL 22:24 | explore の attempt 29 → 25、AFIR 12 → 4、並列化 |
| S4 CH3O• | 8:53（M2） | VAL 14:30 | mode-follow のねじれの case が加わった（実在の反応） |
| tma_hf2 の explore | 124 s（直列、SC） | 32 s（並列、SC）。VAL の s19 は 33 s | artifact は同じ |
| S17 (HF)₂ の paths | 408 s（GS、SB のプローブ） | 73 s（NEB、SB のプローブ） | 固定端 CI-NEB |

---

## 4. 区分ごとの現状の責務と残る課題

### U0 全体プロトコルと理論レベル

**責務:** 段ごとの LOT を決めてエンジンに渡し、出力で観測した Level で照合する。

- 停留点・振動・経路: PBE0-D3BJ/def2-SVPD、grid fine、`convergence energy 1e-7`。
- 低レベル: GFN2-xTB。ωB97X-D3 と CCSD(T) は SP のエネルギー層だけで使う。気相専用。
- ECP と基底: Z>36 は `<El> library def2-ecp` を元素ごとに自動で書く。def2 以外の基底は INPUT_INVALID。
- SCF: 反復 100、救済は 1 回（DFT は cgmin、WFT は前回の vectors から再開）。二重混成は拒否する。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 既定の順位の LOT が PBE0/SVPD | medium | CCSD(T) との差は HCN −1.2、HONO +2.0、SN2 −3.0（SN2 の δG_eff は既定で 11.08、ωB97X-D3 層で 15.73）。既定にしないのは費用によるプローブでの判断。追記した層が thermo を置き換えるのは文書化した使い勝手の制約で、不具合ではない |
| 救済の記述 | low | design.md:170 の「救済は cgmin の 1 回だけ」は DFT にしか当てはまらない。WFT は前回の vectors から再開する |
| 開殻 DFT の cgmin 救済 | 潜在 | cgmin は ⟨S²⟩ を出さないので、開殻の opt・freq・SP・saddle では救済しても必ず捨てられる。r6 の約 2,400 デッキと VAL の 620 デッキで cgmin は 0 件 |
| `freeze atomic` が ECP 原子を凍結しない | 情報 | I を含む CCSD(T) は実需がなく、順位にも関係しない |
| スケール上限（約 35 原子）、def2-SVP の選択肢なし | low | 17 原子の freq は 1 本約 15 分 |
| 実計算で通した範囲 | low | ECP の実計算は I・Te だけ。陽イオン・−2 以下・遷移金属の DFT は未検証 |

### U1 入力構造・電子状態

**責務:** 化学種と組成を、元素・電荷・多重度・結合グラフ・状態ラベルを持つ構造に変える。

- 範囲外は入口で 1 回だけ INPUT_INVALID にする。元素表は `chemistry/elements.py` の 1 つ（72 元素）。
- 多重度: 未宣言なら SMILES は不対電子 + 1、xyz は 1。d ブロックは宣言が必須。組成の多重度は、宣言値かスピン結合で一意に決まる値だけ。
- topology（112 行）は、結合判定（Σr_cov + 0.4 Å）、bond_changes、fragments、reaction_centre、wl_classes、state_label だけを持つ。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 組成の多重度の判定場所が文書と違う | low | 判定は conformers の中で行う。design.md:171 は structures と書いている。無駄になるのは単量体の計算だけ |
| xyz で多重度を書かないと 1 になる | low | 偶数電子の三重項基底（O2、CH2）は黙って一重項になる。検証済みの既定で README に記載。経験則は足さない |
| 閾値の近くでラベルが分かれる | low | amine·HF の配座 83 本中 9 本。S8 の付加体（N···S 2.17 Å、閾値 2.16 Å）は vdW 錯体と同じラベルになる |
| Windows の QA venv で RDKit の DLL がブロックされる | 環境 | SMILES 系の 7 件が落ちる。WSL の prod venv では 451 件合格 |

### U2 配座探索・錯体配置

**責務:** screen に渡す候補構造を作る。同一性の確定とエネルギー窓による選抜はしない。

- 単量体: `crest --gfn2 --quick --ewin 6 --chrg --uhf` で探す。小さく剛直な分子は探索を省く。
- 組成: 円錐配置か剛体配置の seed から、`--nci --notopo <全原子> --noopt` で探す。CREST が失敗したら seed をそのまま出す。
- 選抜は状態ラベルごとに下位 keep_per_state 個。重複除去は CREGEN と Registry の 2 段。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 二重項の組成と Na⁺·H2O で CREST 3.0.2 が失敗する | medium | S5 は rc −11、S6 は rc 1（MTD が収束しない）。seed も xTB で半数以上が落ちる。U2 のコードで直せる簡単な手段はなく、円錐配置の seed が今の最善 |
| 「開殻の組成で必ず失敗」は言い過ぎ | low | 六重項の FeCl3·CH4 は rc 0 で 9 配座を返した。同じ言い過ぎが validation.md:59・110・124、design.md:172、conformer_search.py・placement.py の docstring にある |
| `--ewin` は組成全体の最低構造から効く | low（理屈上） | ラベルをまたいで状態が消えた実例はない |
| `ConformerEnsemble.version`・`job_key` を読むコードがない | low | 版の照合は parse で済んでいる。削除の候補 |

### U3 極小構造の確定

**責務:** 同じ LOT で確かめた極小を、組成 × level_key ごとの basin として Registry に登録する。

- 流れ: opt → 既知の basin なら freq を省く → 別ジョブの freq と is_minimum（|E_freq − E_opt| ≤ qrc_drop）。
- 虚振動は 4 段に分ける。saddle 級はエネルギー目標の変位で mode-follow する。
- 同一性は `identity.assign` の 1 基準（鏡映を含む）で、chiral（m=2）を付ける。鞍点に収束した宣言端点は、近い側の basin に入れる。

| 課題 | 重大度 | 内容 |
|---|---|---|
| mode-follow に変位元の Hessian を渡すと遅くなった | medium | S4 の CH2OH•: Hessian なし（trust 0.1）で 31/31 歩、あり（trust 0.3）で 44/37 歩。硬いモードでの反例はない |
| mode-follow の両側を stage でもう一度緩和し、freq が重複する | medium | s4 で同じ basin の freq が 2 本（30.7 s と 30.4 s）。17 原子なら約 15 分（外挿）。今は再緩和が soft 級の扱いと notes の付け直しも担っている |
| 実計算で通っていない分岐 | medium | DFT の soft_minimum、虚モード 2 本以上の片側 follow、鏡像でない 2 極小への join、本当の state_mismatch、H• の DFT、崩壊した鞍点の登録 |
| is_chiral のしきい値 0.05 Å の余裕 | low | acac は 0.025 と 0.035、DME は 0.014。境界の近くでは m が run ごとに変わりうる（0.41）。ヒステリシスは足さない |
| 行数が増え、同じ判定を回り道で計算する | low | 4 ファイルで 756 → 815 行。`minimum._same` は置換不変 RMSD を 4 回計算するが、判定は「\|ΔE\| ≤ 5e-5 かつ mapped_rmsd ≤ 0.05」と厳密に同じ |

### U4 反応探索

**責務:** 組成 × 状態ごとに、screen 極小の低い順に出発点を選ぶ。

- 列挙: 元素に依らない列挙器 1 つが、T1 移動・置換、T2 リレー、T3 形成、T4 切断の drive を最大 10 件作る。同値類は 1 件にまとめる。
- 1 drive の流れ: ReaDuct の NT2 → Bofill → IRC → 両端 opt。AFIR は NT2 が極大を見つけないときだけ走らせる。
- 採否は `gates.discovery_verdict` 1 つ（ΔE_rxn ≤ 40、低レベル障壁 ≤ 50）。生成物は 1 回だけ同定し、新しい basin だけを species にする。

| 課題 | 重大度 | 内容 |
|---|---|---|
| GFN2 で障壁のない反応は緩和の陰性としてだけ残る | medium | S6 の引き抜き、S5 の CH3+O2 → CH3O2。どちらも DFT の仮説にならない。宣言すれば評価できる |
| AFIR が一度も生成物を残していない | low | SC と VAL で 24 件すべて陰性（same_as_source 20）。狙いの下り坂の会合は、screen が先に緩和してしまう |
| C–H 同士のリレーが attempt の多くを占める | low | NT2 のリレー 40 件中 27 件が C だけで、ΔE‡ 64〜176。共有プロトン極小では max_coordination で N（4）と F（1）が満杯と判定され、T1 が出ない。費用は run の約 3% |
| IRC の端が source を置換した像だと未接続になる | low | s6 の 2 件。生成物はすでに別の drive で残っており、失われた化学はない |
| 緩く結合した開殻の生成物が人工的な極小を作る | low | s6 の CH4O+H が 2 species に分かれ、一方は GFN2 の人工物（H···H–C 1.36 Å） |
| CH3O• の T1（1,2-H 移動）と T4 の実行経路 | low | 実計算で未通過。SC 以後の U4 の差分は挙動を変えないので、S1 の SC の run は HEAD の証拠として使える |

### U5 仮説選択・障壁事前判定・経路

**責務:** 宣言反応と発見の生成物から、同じ LOT の DFT 極小の対を仮説にする。

- 縮退と結合変化は、basin の最適化構造（端点の原子順と掌性に並べ直したもの）で判定する。
- SCREEN: 低レベル TS があれば 3 点で判定し、なければ固定端の IDPP → xTB CI-NEB の内部 9 点で DFT SP を取る。
- 分類器は 1 つ（`profile.classify` と `gates.barrier_verdict`）で、barrierless / single / intermediate に分ける。種が尽きたら FIND_PATH で DFT string を 1 チャンクずつ回す。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 障壁なし・多段・井戸・FIND_PATH の分岐が実計算で未通過 | medium | VAL の 17 case はすべて screen → single → saddle → connect の一本道。SCREEN の barrierless は 11 画像の離散的な上界 |
| mode-follow の発見から、宣言していない配座変化が仮説になり 1 位を取る | medium | S4・S20 の 1 位は rxn_mode_follow（torsional、3.877）。入口は `_degenerate`（結合変化を求めない）と、`_auto` の 0.2 Å 規則の 2 つ。design.md:176 は同じ文の中で矛盾している |
| 近道の種の鞍点探索が失敗すると、NEB を経ずに string へ進む | low | 17 原子なら 1 チャンク約 1.5 h。実計算では未通過 |
| 井戸の緩和に失敗すると same_as_endpoint になる | low | barrierless（well_is_endpoint）で閉じうる。上界としては正しいが、分類は誤りうる |
| SCREEN の SP が 9 本の別ジョブ、IDPP の衝突判定が 0.7 Å | low | 節約・効果が小さいので着手しない |
| actions.py が 494 / 500 行 | low | resample は SCREEN_IMAGES 11 と STRING_BEADS 9 が違うためだけにある |

### U6 鞍点精密化・TS 検証・接続

**責務:** 種を一次鞍点まで精密化し、TS であることを確かめ、両側の極小につなぐ。

- 反応方向: 低レベル TS の虚モード → 反応中心の接線 → 宣言座標 → 最も変わる二面角 → 全原子の接線、の 1 規則。
- 初期 Hessian: 検証済み TS freq → xTB（重なり 0.3 以上）→ DFT の順に選び、`shape_hessian` で反応モードだけ負にする。
- 鞍点探索: NWChem saddle（moddir 1）。maxiter で止まったら、最終フレームから 1 回だけやり直す。
- TS の判定: 負の固有値 1 本。高次 saddle なら押して再試行し、虚振動がなければ collapsed とする。
- 接続: QRC の両側が下ること、別々の極小に割り付くこと、縮退 case では結合の交換も確かめる。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 高次 saddle の再試行が停滞すると、再開が予算外になる | medium | refine_saddle は試行を先に数えるので、行 14 が再開を止める。推論では、行 15 で DFT string を 1 チャンク余分に払い、行 17 で終わる。SB の harness では停滞し、VAL の harness では 37/50 歩で収束した |
| pipeline で通っていない分岐 | medium | higher_order_retry（S11 は harness）、saddle_restart、柔らかい TS から VALIDATE_INTERMEDIATE、split、QRC の側の新しい basin、reassigned、path_hei。S15 と S17 もこれらを通せなかった |
| 15 原子の QRC の長い尾 | low | acac の 2 側で 36 歩・1,387 s と 29 歩・1,318 s。loose で止める案はレビューで見送り済み |
| 柔らかいモードの症状が文書にない | low | 更新 Hessian に 2 本目の負の固有値が十数歩残る。対処（1 回の再開）は文書化済み |
| U6-P7 の照合範囲 | 情報 | 結合の交換を実際に照合したのは、結合が変わる縮退 case だけ（PT、SN2、H3、NH3·HF）。NH3 反転やねじれでは照合は自明に通る |

開殻一重項（RKS の不安定性）は課題ではなく、宣言した範囲の限界（README.md:11）。

### U7 熱化学・速度論

**責務:** DFT 極小と検証済み TS の freq から、GoodVibes 4.3.0 を同一プロセスで呼び、G・H・ZPE を出す。

- 熱化学の中身: qRRHO、pymsym の σ、キラルなら −RT ln 2、単原子は並進と電子だけ。energy_method があれば composite にする。
- 順位の量は δG_eff = max(G_TS, G_R, G_P) − G_R だけ。順方向か逆方向の ΔE0‡ ≤ 0 なら `submerged_barrier`。BARRIERLESS も max(ΔG_rxn, 0) で順位に入る。
- blocker を作るのは `thermo._blockers` だけ。

| 課題 | 重大度 | 内容 |
|---|---|---|
| σ が構造の揺らぎで変わりうる | medium | 実 run では 126 本中 1 本（S2）が影響した。RMS 1e-3 Å の乱数の揺らぎでは 20 回中 5〜19 回変わる（最悪側の見積もり）。誤差は RT ln σ で 0.41〜1.06。pymsym は版を固定しておらず、preflight も照合しない |
| 基準の状態の G がスピン汚染を見ない | low | `_state_G` と `_association` の単量体は spin_contaminated の配座も候補にする |
| 「状態の G」の定義が 2 つある | low | 反応物・生成物は最小の G、単量体は Boltzmann アンサンブル。ΔG_assoc の基準がそろわない（表示だけで、順位には効かない） |
| BARRIERLESS と別の配座を基準にする場合 | low | QM なしのテストでしか通っていない |
| 縮退反応の RT ln 2 の規約 | low（要確認） | 片方向の TST 値で、交換速度とは 2 倍違う可能性がある。縮退どうしでは序数は変わらない |
| トンネル補正がない | low | 範囲外と明記（HCN の κ は 5.4 の見込み） |
| Windows の QA では golden がスキップされる | low | validation.md:82 の「既定のテストで毎回通る」は Windows には当てはまらない |

### U8 一点計算・手法パネル・順位付けと報告

**責務:** SP は順位に使う点だけを計算し、parents で対象に結び付ける。

- thermo: energy_method があれば composite G にし、SP 層にも spin_ok をかける。
- report: outcome と blocker で rankable を決め、δG_eff で並べる（感度幅が重なれば同順位）。手法の幅は method_panel.csv と ranking.csv の min/max 列で示す。
- CCSD(T): 閉殻は RHF、開殻は ROHF + TCE。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 既定の順位が PBE0/SVPD | medium | U0 と同じ。マロンアルデヒドは SVPD で submerged、ωB97X 層では 1.47 と、定性的な判定まで変わる |
| 複数条件の後にパネルを追記すると、古い層の (T, state) で並びうる | low | 置き換わるのは同じ (T, 標準状態) の記録だけ。method_panel.yaml と design.md:190 の「上書き」は不正確 |
| report.html が順位の条件と層にそろっていない | low | `_summary_cells` と `_section` が thermo[0] を使う。rank と同じ行の δG_eff が別の条件の値になりうる。二重定義で、render は CC 15 |
| パネルのほかのレベルに spin_ok がかからない | low | 範囲は design.md:66 に記載 |
| CCSD(T) の参照は重原子約 4 個まで | low | マロンアルデヒドの (T) で GA の確保に失敗した |
| 実計算で未検証の経路 | low | BARRIERLESS の順位、discover の run へのパネルの追記 |

### U9 実行基盤・パイプライン・設定

**責務:** stage の唯一の入口 `execute_stage`（型検査、manifest、再開、SiteLock）と、JobRunner（内容アドレス型の JobStore、LADDER、コアのセマフォ）。

- 外部プロセスは process.py だけが起動し、SIGTERM・SIGHUP ではグループごとに止める。
- 失敗の閉じ込めはジョブ・ケース・stage の 3 単位。report は必ず出る。
- knob は 8 個。結果を変える入力はすべて鍵に入る。

| 課題 | 重大度 | 内容 |
|---|---|---|
| 実計算で通らない分岐 | medium | U5・U6 と同じ。失敗の閉じ込めと子ケースの締め切りの継承も、VAL では一度も発動していない |
| SIGINT（Ctrl-C）を処理していない | low〜medium（推論） | 外部プロセスは別セッションにいて SIGINT が届かない。並列の stage では ThreadPoolExecutor の終了待ちで、残りの単位がすべて終わるまで止まらない |
| `max_saddle_attempts` の説明が実装と違う | low | 実装は 1 ケースごと（子ケースごとに 0 から数える）。design.md:234 と docstring は 1 仮説と書いている |
| 閉じ込めたケースの極小が manifest に残らない | low | その反応の熱化学が欠けるだけで落ちない。まれな経路なので対処しない |
| SIGTERM の後、未着手の単位が起動しては止められる | low | 孤児は残らず、遅れも小さい |

---

## 5. 今後の候補（有用で単純なものだけ）

確かめた段が採用したものだけを挙げる。どれも削除か数行の修正で、新しい仕組みは足さない。

### 5.1 コード

| # | 区分 | 内容 | 規模 |
|---|---|---|---|
| 1 | U5 | 宣言していない対は、結合変化のあるものだけを仮説にする。`_degenerate` に `not record.torsional` を足し、`_max_distance_change` と `min_distance_A` を削除する。design.md:176 は「結合変化のない未宣言の変化は仮説にしない」に直す。FHN のテストを結合変化の有無で書き直す | 約 −10 行 |
| 2 | U3 | mode-follow の再最適化から、変位元の Hessian を外す（`_Ctx.relax` の source 引数を削除し、trust 0.1 に戻す）。S4 の歩数（31/31 対 44/37）を validation.md に 1 行残す | 削除だけ |
| 3 | U3 | `minimum._same` を identity の小さな述語（\|ΔE\| ≤ 5e-5 かつ mapped_rmsd ≤ 0.05）に置き換え、呼び出し元がなくなる `identity.same_basin` を削除する | 数行の削除 |
| 4 | U3 | ts_candidate に両側の MinimumOutcome（soft 級なら `_soften` 済み）を持たせ、minima stage は再緩和をやめて Registry.add で登録する | ts_candidate 1 件ごとに freq 1 本と opt 2 本が減る |
| 5 | U6 | saddle_restart を試行予算に数えない（種を積むときに `saddle_attempts` を戻す）。refine_saddle の docstring と design.md:158 を同時に直し、decide の単体テストを 1 件足す | 1 行 + テスト 1 件 |
| 6 | U7 | `_state_G` と `_association` の単量体の候補を、`_blockers` と同じスピンの条件で絞る | 数行 + テスト 1 件 |
| 7 | U7 | pyproject で pymsym の版を固定する | 1 行 |
| 8 | U8 | HTML で、RankRow と同じ (T, 標準状態) の ReactionThermo を選ぶ小さな関数を `_summary_cells` と `_section` で共用し、energy_level 列を足す（render の CC は上げない）。RankRow の古いコメントを直す | 数行 |
| 9 | U9 | `stop_on_signals` に SIGINT を加える（終了コード 130） | 1 トークン |
| 10 | U2 | `ConformerEnsemble.version`（と `job_key`）を削除する | 純粋な削除 |
| 11 | U4 | AFIR を削除する（worker の `_afir`・`afir_pair`、γ の 2 設定、explore のフォールバック、discovery_verdict の扱い）。根拠は SC と VAL の 24 件がすべて陰性で、生成物が 0 件であること。validation.md に 1 行書く | 削除 |
| 12 | U5（任意） | 井戸の緩和が失敗したときは barrierless で閉じず、最も高い山を種にする | 約 2 行 |
| 13 | U7（任意） | `thermo.ensemble_G` を削除し、単量体も `_state_G`（最小の G）で基準にする。「状態の G」の定義が 1 つになる。変わるのは表示用の ΔG_assoc と分離基準の値だけ | 削除 + テストの期待値を修正 |

### 5.2 文書（新しい文は足さず、置き換える）

- design.md:170 と §7 の ladder: 救済は「DFT は cgmin、WFT は前回の vectors から再開」。
- design.md:171: 組成の多重度は conformers の入口で判定する。
- CREST が失敗する範囲の記述（validation.md:59・110・124、design.md:172、conformer_search.py と placement.py の docstring）を、「二重項（ラジカル）の組成と Na⁺·H2O（rc −11 か、MTD が収束せず rc 1）。六重項の FeCl3·CH4 は rc 0」にそろえる。
- validation.md の S8 の行: 旧規則（`--notopo` なし）では付加体が残らず、最低でも +17.3 の vdW 構造だけになる。
- design.md:158: 柔らかい反応モードや押した種では、2 本目の負の固有値が十数歩残って停滞しやすい、を括弧で加える。
- design.md:234 と ReactionPathsPolicy の docstring: `max_saddle_attempts` は 1 ケースごとで、子ケースと共有するのは walltime_h の締め切りだけ。
- method_panel.yaml と design.md:190: 置き換わるのは同じ (T, 標準状態) の値だけ。
- 品質ゲートの pytest は WSL の prod venv で回す（451 件が基準）。validation.md:82 もこれに合わせる。Windows の RDKit DLL と σ の golden のスキップを、この 1 つで扱う。
- GFN2 で障壁のない反応は、宣言すれば評価できる（1 文）。

### 5.3 検証 run（WSL、1 本ずつ、コードの追加なし）

1. S4 CH3O• の discover を explore まで流し（`--until explore`）、T1（1,2-H 移動）と T4 の結果を記録する。
2. screen を外した一時 pipeline で HCN を流し、行 15 → 1 チャンク → classify → path_hei → saddle を一度通す。
3. 水二量体の acceptor switching を縮退反応として宣言し、行 6（BARRIERLESS）を 1 回試す。single になっても記録して終える。
4. C2v の DME（または平面 NH3）を宣言端点にして、一次鞍点からの mode-follow と、5.1 の 2 を確かめる。
5. 5.1 の 1 の後で S4 を 1 回流し直し、順位が宣言反応の 1 行になることを確かめる。

### 5.4 採らないもの（理由は §4 と各区分の確認記録）

- 精度と費用: composite を既定にすることは再開しない。def2-SVP の method と SVP/SVPD の比較、ECP 原子の freeze の追加設定は作らない。
- SCF: cgmin 救済を閉殻に限る案は、発生が 0 件で、continuation も CC 13 なので保留する。
- 小さな判定: WL を安定するまで回す案、spin_tol を相対値にする案、is_chiral のヒステリシス、鏡映対称な鞍点で 2 側目を省く案は採らない。
- CREST と探索: CREST の代替サンプラー・リトライ・剛直分子の細かい判定、relaxation を自動で仮説にする機構、C–H リレーを除く元素別フィルタ、IRC の端の添字の付け直しは採らない。
- 経路と鞍点: U5-P9、U5-P11、RKS 不安定性の検出、saddle_hessian の注記の拡充は採らない。
- σ: 自前判定は、実 run で 2 件目の誤りが出るまで保留する。揺らぎを加えた golden も入れない。QA venv への GoodVibes の導入はしない。
- 実行基盤: 閉じ込めたケースの極小を manifest に出す仕組み、thread_map の cancel は作らない。
- その他: M06-2X は入れない。

---

## 6. コード量・品質指標の推移

### 6.1 行数（`hfauto/` の .py、tests の .py）

| コミット | 段 | hfauto 行数 | 増減 | .py ファイル | tests 行数 | pytest 合格（既定） | 最大ファイル |
|---|---|---|---|---|---|---|---|
| `fbf7924` | レビュー時 | 10,635 | — | 80 | 5,259 | — | 500 |
| `ebafbcf` | M0 | 10,696 | +61 | 80 | 5,384 | 318 | — |
| `9400787` | M1 | 10,837 | +141 | 81 | 5,604 | 351 | — |
| `a80b01e` | M2 | 10,860 | +23 | 81 | 5,838 | 371 | ≤ 500 |
| `b530f1d` | M3 | 10,859 | −1 | 81 | 5,916 | 376 | 500 |
| `14b700c` | M4 | 10,865 | +6 | 81 | 6,038 | 385 | 500 |
| `82c46d4` | M5 | 10,788 | −77 | 81 | 6,105 | 391 | 500 |
| `6d41726` | S-A | 10,818 | +30 | 81 | 6,304 | — | — |
| `36f5cae` | S-B | 10,780 | −38 | 81 | 6,602 | — | — |
| `0f1d868` | S-C | 10,761 | −19 | 81 | 6,726 | — | — |
| `ac8f802` | S-D | 10,481 | −280 | 78 | 6,806 | 442 | 499 |
| `496fad0` | FINAL | 10,311 | −170 | 78 | 6,851 | 450 | 494 |
| `d5fd24c` | 検証 | 10,336 | +25 | 78 | 6,868 | 451 | 494 |

- 通算は −299 行（−2.8%）、M5 からは −452 行。レビューの見込み（正味 −700）には届いていない。
- テストは +1,609 行。テスト関数は 225 → 299 本で、ほかに golden の実出力を加えた。
- 「—」はコミットの記録に数がないもの。M0〜M5 の最大ファイルは記録上の値。

### 6.2 区分ごとの主なファイルの行数

| 対象 | 改良前 | 改良後 |
|---|---|---|
| topology.py（U1） | 239 | 112 |
| trials.py（U4） | 254 | 176 |
| summary.py（U8） | 300 | 239 |
| pysis の engine + worker（U5） | 295 | 224 |
| U7 の中核（backends/goodvibes + thermo + thermochemistry） | 620 | 495 |
| U0 の 4 ファイル（input・engine・output・method、`07521ed` 比） | 949 | 922 |
| U3 の 4 ファイル（minimum・minima・identity・selection） | 756 | 815 |
| U4 の 4 ファイル（explore・trials・worker・engine） | 787 | 782 |
| U6 の 9 ファイル | 2,078 | 2,031 |

### 6.3 文書の行数

| ファイル | `fbf7924` | `82c46d4` | `496fad0` | `d5fd24c` |
|---|---|---|---|---|
| README.md | 110 | 112 | — | 93 |
| docs/design.md | 228 | 231 | — | 267 |
| docs/validation.md | 217 | 347 | 542 | 129 |
| docs/environment.md | 70 | 71 | — | 75 |

### 6.4 構造の指標と品質ゲート（`d5fd24c`）

| 指標 | 改良前 | 改良後 |
|---|---|---|
| 上書きされない knob | 23 | 8（Policy 5 + ReactionPathsPolicy 3） |
| 元素表 | 4 | 1 |
| 経路の形の判定器 | 2 | 1 |
| blocker の生成元 / SP と停留点の対応付け | 3 / 3 | 1 / 1 |
| 重複除去の段 | 3 | 2 |
| MethodSpec のフィールド | 11 | 9 |
| 熱化学の外部エンジン（worker・JobStore・Capability） | あり | なし（同一プロセスの純関数） |

| ゲート | 結果 |
|---|---|
| pytest（既定） | WSL の prod venv で 451 件合格。Windows の QA venv では 435 件合格・7 件失敗で、失敗はすべて RDKit の DLL がブロックされる環境要因。基準は WSL の 451 件 |
| 実エンジンの smoke | 13 件合格（64 s） |
| ruff | 問題なし |
| pyrefly | 0 errors |
| lint-imports | 5 契約すべて維持 |
| radon CC | 最大 15（`reporting/html.render`、上限ちょうど）。次は C(14) の `ReactionPathsStage` |
| ファイルの最大行数 | 494（`drivers/reaction_case/actions.py`、上限 500） |

上限に近いのは html.render（CC 15）と actions.py（494 行）の 2 つである。5.1 の 8 は、選択を小さな関数に出して render の CC を上げない形にする。actions.py に何か足すときは、ファイルの分割が先に要る。
