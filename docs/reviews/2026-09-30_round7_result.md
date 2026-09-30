# hfauto 計算基盤の改良 round 7 — 結果報告(2026-09-30)

- 対象: ブランチ `refactor/fundamental-2026-09`、W0 `d549945`、W1 `707e656`、W2 `26389e3`、W3 `2eaaf53`、W4 `c5b515e` と W5(削除と文書)。push はしていない。
- 基準: 前回の結果報告 `docs/reviews/2026-09-28_platform_improvement_result.md`(以下「前回」、`f2309b9`)の §0.3 の課題と §5 の候補、レビュー `docs/reviews/2026-09-27_platform_review.md` の方向性(§0.3)と検証セット(§6)。
- 方法: WSL の実 run(`/home/user/hfauto_r7/W0`〜`W5`・`VAL7`、新しい run dir で 1 本ずつ)。分岐を通したと数えるのは、実 run の `log.jsonl` か diagnostics.json に現れたときだけ。VAL(前回の検証)との比較は生の JSON だけを読む。詳細は [validation.md](../validation.md)。
- W5 の全体の再検証は VAL7(`/home/user/hfauto_r7/VAL7`、35 run)と、行 8 の修正の後の `/home/user/hfauto_r7/W5`(S6 を paths から流し直したものと R1〜R4)。

## 0. 要旨

- **正しさ**: 宣言していないねじれ・配座変化は仮説にしない(S4 の 1 位だった CH2OH• のねじれが消えた)。検証済みの DFT 鞍点(mode-follow、分割した親)は探し直さず検証から始める。障壁なしは中点 2 点の DFT SP を足してから受け入れる。片側だけ端点に届く TS は多段に分割する(GEN-05)。締め切りと freq の電子状態の不具合を実 run で見つけて直した。
- **汎用性**: GFN2 で障壁のない引き抜き(S6)は、screen で失われた seed の状態を DFT に 1 回問うことで、W3 では DFT の TS と順位(δG_eff 1.48)に届いた。ただし VAL7・W5 では向きを変える鞍点に収束し、この TS を再現していない(§5 の 1)。IRC の端が source の添字を付け替えた像のときも接続した段として残す。AFIR は M4 より後の 41 件で生成物 0 だったので削除した。
- **無駄な計算の削減**: VAL と対になる 16 run で QM 時間 12,743 → 8,334 s(−34.6%)、ジョブ 369 → 325。一次の鞍点の Hessian を正定値のモデルにし(acac の QRC 34/27 → 15/14 歩)、厳密な像(対称な TS の − 側、置換・鏡像の端点)を計算せず、反応しない組成を DFT に流さず、mode-follow の側を緩和し直さない。ケースの中の独立な SP と QRC の両側は rank を分けて同時に走らせる(SCREEN の SP の区間 −39〜−41%)。
- **分岐の到達**: 決定表の行 3・5・6・7・13・15・17、multi_step、barrierless、reassigned、FIND_PATH、`path_hei`・`saddle_restart` を実 run で初めて通し、`higher_order_retry` を pipeline で初めて通した。行 8 は VAL7 の S6 で初めて通った(そこで予算を残した打ち切りを見つけて直した)。行 1・11・16 は通っていない。
- **行数**: `hfauto/` は 10,336 → 10,609(W5、+273)。W5 は削除と行 8 の修正で −10。

### 0.1 前回の課題(§0.3)の状態

| # | 前回の課題 | 状態 |
|---|---|---|
| 1 | 実計算で通っていない分岐 | 大半を通した(§3.2)。行 1・11・16、2 回目の SCREEN、ケースの例外の閉じ込めは未達 |
| 2 | 既定の順位は PBE0/SVPD | 変えていない(設計の選択)。H 引き抜きの過小評価(分離基準 0.79、文献 約 5)を加えて記録 |
| 3 | 開殻の組成で CREST が失敗 | 原因を訂正: GFN2 に錯体の極小がない(障壁のない会合・引き抜き)か SCC の失敗で、開殻そのものではない(六重項 FeCl3·CH4 は rc 0)。失われた状態は DFT に問う |
| 4 | σ が揺らぐ | pymsym を 0.3.5 に固定。しきい値 2e-3 の方式は同じ |
| 5 | 未宣言の配座変化が仮説になる | 解消(結合変化のある対だけ) |
| 6 | mode-follow の効率 | 解消(正定値のモデルで 36 → 13 歩、側の再緩和と重複 freq をなくした) |
| 7 | 停滞の再開が予算外 | 解消(再開は数えず、上限でも 1 回通す。試行が残らなければ string を走らせない) |
| 8 | GFN2 で障壁のない反応が陰性だけ | S6 は W3 で DFT の TS と順位に届いたが、VAL7・W5 では再現しない。S5 は UKS のスピン汚染で結論なし |

## 1. 実施した改良

| wave | コミット | 主な内容 | 主な削除 | hfauto 行数 |
|---|---|---|---|---|
| W0 | `d549945` | 経路の action を `paths.py` に分けた。VAL のジョブ時間の基準と事前の数え上げ(付け直しの対象 4 件、R6 のジョブ、cgmin の証拠、pymsym の版) | — | +18 |
| W1 | `707e656` | 結合変化のある対だけを仮説に。`ts_calc` の検証から始める case。行 15・16 は試行が残るときだけ、再開は数えない。近道の種が失敗したら string の前に NEB を 1 回。障壁なしの高密度化。井戸の緩和の失敗は結果にしない。GEN-05。状態の G を 1 定義に。HTML を順位の (T, 状態) に。SIGINT。開殻 DFT の cgmin をやめる。pymsym の固定。分岐用の系 6 つ。締め切りの意味の修正 | `BarrierVerdict.max_rel_kcal`・`max_node_spacing_A`、`ConnectionClaim.amplitude_A`、0.2 Å の距離規則、`thermo.ensemble_G`、2 つ目の `_state_G` | +84 |
| W2 | `26389e3` | 一次の鞍点の Hessian の正定値モデル(K1)。対称な TS の − 側と像の端点を計算しない(`identity.carry`、`IMAGE_A`)。mode-follow の側をジョブなしで登録し、収束した QRC の側を再 opt しない。反応する組成だけを DFT に。rerank は混んだ組だけ。`identity.same_as_labelled` | `identity.same_basin`、`minimum._same`、`Registry.members`、側の 2 回目の緩和、選択の `include` 引数 | +108 |
| W3 | `2eaaf53` | IRC の端の添字の付け直し(R5a)。失われた seed の状態を DFT に 1 回(R6、`collapsed_at_dft_from_seed`)。freq の SCF を opt の vectors から(`scf_guess`) | AFIR 一式、explore の seed ごとの `collapsed_to` の記録、`NO_NT2_MAXIMUM` 定数、`ConformerEnsemble.version`・`job_key` | +28 |
| W4 | `c5b515e` | 独立な NWChem ジョブの同時実行(rank の分け前、timeout の伸長、`--bind-to none`)を SCREEN の SP・中点の SP・非対称な QRC の両側に | 極小の stage の段階実行は実 run で遅くなったので入れなかった | +45 |
| W5 | (この報告) | 使われなくなったフィールドの削除、行 8 の修正(両側が同じ basin の TS で試行が残れば探索を続ける)、文書、VAL7 | `ReactionThermo.dzpe_act_kcal`、`DiscoveryResult.job_key`、`PathProfile.job_key`、`reaction_delta` の ΔZPE‡ | −10 |

## 2. 評価の変化(区分 × 観点、前回 → round 7)

| 区分 | 妥当性 | 汎用性 | 有用性 | 効率 | 簡潔さ |
|---|---|---|---|---|---|
| U0 プロトコル | ○→○ | ○(条件付き)→○(条件付き) | ○→○ | ○→○ | ◎→◎ |
| U1 入力・電子状態 | ○→○ | ○→○ | ○→○ | ○→○ | ◎→◎ |
| U2 配座・錯体配置 | ◎(条件付き)→◎(条件付き) | ○→○ | ○→○ | ◎→◎ | ◎→◎ |
| U3 極小の確定 | ◎(但し書き)→◎(但し書き) | ○→○ | ◎→◎ | ○→◎ | ○→○ |
| U4 反応探索 | ○→○ | ○(条件付き)→○ | △→○(条件付き) | ◎→◎ | ○→○ |
| U5 仮説・経路 | ○→○ | ○→○ | ○(但し書き)→◎ | ○→○ | ○→○ |
| U6 鞍点・接続 | ○→○ | ○→○ | ○→◎ | ○(QRC は △)→○ | ○→○ |
| U7 熱化学 | ○→○ | ○→○ | ○→○ | ◎→◎ | ○→◎ |
| U8 一点計算・報告 | ○→○ | ○→○ | ○→○ | ◎→◎ | ○→○ |
| U9 実行基盤 | ◎→◎ | ○→○ | ◎(留保)→◎ | ○→◎ | ○→○ |

根拠:

- **U3**: 効率 ○→◎。mode-follow の側は 36/36 → 13/13 歩でジョブなしに登録され(S4 の opt 4 → 2、freq 5 → 4)、像の端点はジョブ 0、反応しない組成は DFT に流さない(S10 の dft の freq 861 → 30 s)。但し書きだった未通過の分岐のうち、DFT の `ts_candidate` と `one_side` を実 run で通した。DFT の soft_minimum、崩壊した鞍点の登録などは未達なので但し書きは残る。
- **U4**: 有用性 △→○(条件付き)。S6 の H 引き抜きが W3 で DFT の TS(−481.1i)と 1 位に届き(VAL7・W5 では再現せず)、付け直した IRC の端が S6 で 2 件、S19 で 2 件残った。条件は、GFN2 の探索そのものが生む生成物は少ないままで、ラジカル会合(S5)は UKS のスピン汚染で結論が出ないこと。汎用性の条件(生成物まで確かめたのが 3 系だけ)は、CH3O• の 1,2-H 移動(両向き)と CH3···O2 の生成物が加わったので外した。
- **U5**: 有用性 ○(但し書き)→◎。未宣言のねじれが順位の先頭に載らなくなり、多段(oxalic、S6)と障壁なし(H2O·HF、水二量体)が順位表に並ぶ。妥当性は ○ のまま: 高密度化は対称な経路でしか通っておらず、S6 では緩い錯体の分割が深さの上限に達した。
- **U6**: 有用性 ○→◎。検証済みの鞍点を再利用する case は新しいジョブ 0(nh3_planar_seed、oxalic split1)。効率は QRC の側 −27%(VAL と対の 149 → 109 s/側)と同時実行(−18〜−21%)で △ を外した。行 8・11 が未達なので妥当性と汎用性は ○。
- **U7**: 簡潔さ ○→◎。状態の G は 1 定義(spin_contaminated でない最小の G)になり、`ensemble_G` を削除した。妥当性は ○ のまま: 近直線の錯体の直線性の判定が構造の揺らぎで変わる(S18、δG_eff 4.091 と 6.360)。
- **U9**: 有用性の留保(SIGINT)を外した(rc 130、0.8 s で外部プロセスが消える)。効率 ○→◎: 同時実行は同じ鍵で結果が変わらず(最大 8.9e-9 Eh)、効果のない極小の stage には入れなかった。
- 変えなかった区分: U0 の freq の電子状態の修正は黙った棄却を 1 件防いだが、既定の LOT の系統誤差が残る。U8 の HTML の条件の一致は前回の low の課題の解消で、評価は変わらない。

観点ごとの全体: 化学的妥当性 ○、汎用性 ○(条件付き。陽イオンと遷移金属の DFT は未検証)、有用性 ○(順位表の偽の先頭がなくなり、多段・障壁なしも並ぶ。探索は条件付き)、効率 ○ → ◎、簡潔さ ○(規則は 1 つずつのままだが、行数は W4 まで +283)。

## 3. 検証結果(W0〜W5、VAL7)

### 3.1 合否

| 群 | 系 | 結果 |
|---|---|---|
| 回帰 | R1〜R4 | HCN 42.1784、HONO 11.8159、NH3 3.8395、水 same_basin。VAL との差 ≤ 1.2e-4(freq の SCF 修正による)。VAL7 と W5 は W4 と 1e-6 以内 |
| known_endpoints | S1、S2、S12〜S14、S16、S18 | outcome と basin が VAL と同じ、δG_eff の差 ≤ 0.001。S4 は 1 位が宣言反応 30.612 になった(意図した変化) |
| 回帰だけ | S15、S17 | elementary 52.016、degenerate 1.392。多段・障壁なしの被覆には数えない |
| 新しい系 | formaldehyde、h2o_hf_inversion、water_dimer_as、oxalic_two_step、nh3_planar_seed、dme_c2v_seed、hcn・hf_dimer_swap の noscreen、hono の walltime | すべて期待の分岐に入った。h2o_hf_inversion の行 6 はコードの通過だけ(隠れた障壁の検出は未検証)。hf_dimer_swap_noscreen は higher_order_retry に入らず直接 C2h の TS に届いた |
| discover | S5、S6、S7、S8、S10、S19、S4 の explore | S6 は合格(ΔE‡ 2.40、δG_eff 1.48、1 位)、S5 は結論なし(⟨S²⟩ 1.71)、S7・S8 は生成物 0 で VAL と同じ、S10 は 32.2892、S19 は seed が DFT でも崩れた(`collapsed_at_dft_from_seed`) |
| VAL7 | §6 の全体と round 7 の系(35 run) | known_endpoints・新しい系・S3・S7〜S10・S19・パネル(CCSD(T) は閉殻 HCN −1.19、開殻 CH3O• −0.36)は W4 以前・VAL と同じ(δG_eff の差 ≤ 2e-4)。S5 は別の鞍点で reassigned(順位なし)、S6 は引き抜きの TS を再現せず unresolved |
| 実行基盤 | SIGINT、s20_sn2_panel(400 K・1 M) | rc 130・孤児 0。report.html の δG_eff が ranking.csv と一致 |

### 3.2 到達(round 7 の実 run だけ)

| 対象 | 通った | 通っていない |
|---|---|---|
| 決定表の行 | 2〜10、12〜15、17(行 8 は VAL7/s6_oh_ch4 と W5/s6_oh_ch4_rowfix) | 1、11、16 |
| saddle の種 | `screen_ts`、`screen_hei`、`discovery_ts`、`path_hei`、`higher_order_retry`(TS freq の Hessian)、`saddle_restart` | — |
| outcome | elementary、degenerate、same_basin、out_of_window、multi_step、barrierless_at_resolution、reassigned、unresolved_within_budget | blocked_upstream |
| その他 | `ts_calc` の検証、GEN-05、`minus_is_image`、`image_of`、`not_reacting`、`collapsed_at_dft_from_seed`、付け直した IRC の端、FIND_PATH の string | 2 回目の SCREEN、ケースの例外の閉じ込め、DFT の soft_minimum |

## 4. 計算量

| 比較 | 結果 |
|---|---|
| VAL → W2(対の 16 run、QM 秒) | 12,743 s → 8,334 s(−34.6%)、ジョブ 369 → 325。QRC の側 4,477 → 2,064 s、dft の freq 3,357 → 2,432 s、dft の opt 2,074 → 1,203 s、rerank の SP 90 → 0 s |
| W1 → W2(対の 13 run) | 10,008 s → 6,480 s(−35.2%) |
| W3 の追加(R6) | S19 +2,033 s、S6 +536 s、S5 +208 s(失われた状態を DFT に問う費用) |
| W2・W3 → W4(壁時計) | S16 61:20 → 52:42、S14 10:35 → 9:38、S1 5:42 → 4:19、HCN 1:25 → 1:01、HONO 2:59 → 2:18、oxalic 38:16 → 33:32、S19 65:16 → 64:53 |
| VAL → VAL7(対の 24 run) | QM 14,071 → 13,579 s(−3.5%)、ジョブ 467 → 406、壁時計 14,385 → 11,164 s(−22.4%)。S19 の R6(+2,308 s)を除けば QM −22.9%、壁時計 −43.3% |
| W4 → VAL7(同じ 10 系) | QM 13,627 → 13,440 s、壁時計 10,337 → 10,528 s(+1.8%、S10 の QRC の側の歩数と S19 の seed の 1 歩による)。2 回目の試行は S10 の AUTOZ の失敗 1 件、孤児 0 |
| VAL → W4(壁時計) | S16 1:30:57 → 52:42、S14 16:52 → 9:38、S10 22:24 → 3:05、S1 7:30 → 4:19、HCN 1:17 → 1:01、HONO 3:03 → 2:18 |

## 5. 残る課題

| # | 区分 | 重大度 | 内容 |
|---|---|---|---|
| 1 | U5・U6 | high | S6 の引き抜きの TS が再現しない: VAL7・W5 は OH···CH4 の向きを変える鞍点(−73.0i、−58.6i、両側が同じ basin)に収束して試行を使い切った。虚モードが結合変化を担わない鞍点を QRC の前に見分け、結合変化の方向に拘束して探す必要がある |
| 2 | U5・U6 | medium | 行 1・11・16、2 回目の SCREEN、ケースの例外の閉じ込めが実計算で未達。高密度化は隠れた障壁で確かめていない |
| 3 | U0・U8 | medium | 既定の順位は PBE0/SVPD。SN2 −3.0、HONO +2.0 kcal/mol(対 CCSD(T))、H 引き抜きは分離基準 0.79 と文献 約 5 |
| 4 | U4 | medium | ラジカル会合(S5)は UKS の ⟨S²⟩ 1.71 で判定できない。多参照かスピン射影が要る |
| 5 | U7 | medium | 近直線の錯体の直線性が構造の揺らぎで変わる(S18: 4.091 と 6.360)。対称性で決める必要がある |
| 6 | U3・U4 | low | 失われた seed が DFT の 1 歩目で状態を離れても、既知の basin に入るまで opt を続ける(S19 で 44 歩・約 33 分)。ただし BFGS の途中でラベルが変わることは障壁なしの証明ではない |
| 7 | U6 | low | 柔らかい引き抜き・会合で saddle が 2 回 maxiter で止まってから string に進む(S5 で saddle 826 s) |
| 8 | U3 | low | 対称な鞍点の mode-follow の ± の側は厳密な像なのに両方計算する(S6 で 329 s、NH3 は 21 歩/側) |
| 9 | U5 | low | S6 の CH3OH + H 側は、ほぼ縮退した緩い錯体の分割で `max_split_depth` に達し、順位が付かない |
| 10 | U9 | low | 極小の stage の同時実行は、大きさと歩数に応じた割り当てがなければ遅くなる(S19)。今は直列 |
| 11 | U4 | low | mode-follow の発見の片端の化学種は、鞍点の側の化学種ではなく basin の代表から取る(`DiscoveryRecord` に出所の化学種がない) |
| 12 | U0 | 潜在 | 開殻 DFT には SCF の救済がない(cgmin の後に ⟨S²⟩ を出す 2 回目の SCF が要る)。実 run での発生は 0 |

## 6. 採らなかったもの

| 案 | 理由 |
|---|---|
| QRC の早期停止・チャンク化(CA-3a、R19) | K1 の後の acac の尾は basin に入ってから 7/15 歩と 6/14 歩、ΔE 0.004 / 0.009 kcal/mol。チャンクは再開と複雑さを足す |
| 固有値の下限のプローブ(CA-3b) | K1 が同じことを覆う |
| 負の曲率を持つ Hessian すべてをモデルにする(K1 の計画の規則) | 高次の鞍点からの側が這う(DME C2v で 29 歩 → 207 歩でも未収束)。xTB の錯体の鍵と S18 の熱化学が変わる |
| mode-follow の − 側を像にする(CA-1c) | K1 の後は約 40 s で、走っていないジョブの構造ファイルが要る(S6 では 329 s なので課題 8 に残す) |
| IRC の前の障壁の打ち切り(CA-8) | explore は壁時計の 0.8% |
| SCREEN の SP の重複除去と 1 ジョブ化(CA-9) | 最大 30 s。同時実行で区間は −39〜−41% になった |
| CREST の電子温度 1,000 K(K4) | 組成ごとに低レベルの PES を変える新しい規則になる。失われた状態は DFT に問う |
| linopt 0(K6、U3-P11) | 速くなる系と遅くなる系が混ざる |
| RI-J、ジョブをまたぐ movecs の連鎖(K7) | 解析 Hessian がなく、同一 PES のゲートが崩れ、鍵の外で解を左右する。freq の `scf_guess` は電子状態をそろえるためだけに、元のジョブの鍵を含めて使う |
| 近い basin への known の省略(K8) | 同一性の基準を緩めることになる |
| 掌性を libmsym で決める(K9) | 誤りの実例がなく、HONO の TS の m = 2 を危うくし、許容を 0.05 Å から約 0.006 Å に変える |
| 別の 2 状態を結ぶ IRC の段を発見にする(R5b) | source が drive のものでない発見の記録が要る |
| SCREEN の画像数と string の bead 数の統一 | DFT の string のチャンクが 11 bead になる |
| 行 8・11・16 と 2 回目の SCREEN を強制する系 | 自然に通る系がなく、期待や閾値を動かして発火させることはしない |
| 極小の stage の段階実行と rerank の SP の同時実行(K2 の一部) | S19 で等しい rank の割り当てが遅くなり(TMA 105 → 206 s、2 コアが遊ぶ)、S10 も変わらない。rerank の SP は発生しなかった |
| AFIR を残す | M4 より後の 41 件で生成物 0。M4 の生成物は NT2 か screen が既に持つ状態 |

## 7. コード量と品質指標

| コミット | 段 | hfauto 行数 | 増減 | tests 行数 |
|---|---|---|---|---|
| `f2309b9` | 前回 | 10,336 | — | 6,868 |
| `d549945` | W0 | 10,354 | +18 | 6,869 |
| `707e656` | W1 | 10,438 | +84 | 7,168 |
| `26389e3` | W2 | 10,546 | +108 | 7,396 |
| `2eaaf53` | W3 | 10,574 | +28 | 7,530 |
| `c5b515e` | W4 | 10,619 | +45 | 7,575 |
| W5 | W5 | 10,609 | −10 | 7,591 |

- W5 の時点で最大のファイルは `hfauto/drivers/reaction_case/actions.py` の 461 行、radon の最大 CC は 14。WSL の pytest 498 件(実計算の smoke 13 件も VAL7 で成功)、ruff・pyrefly 0 件・import-linter 5 契約。
- 行数の目標(前回の 10,336 を下回る)は達していない。W1〜W4 で足した機能(R5a、R6、K1、K2、GEN-05 など)の +283 に対し、W5 の走査で呼び手のない定義は見つからず、削除できたのは読まれないフィールドだけだった。
- 文書は design.md 267 → 268 行、validation.md 129 → 159 行(round 6 の記録は git の履歴へ)、environment.md 75 → 76 行、README 93 行。
