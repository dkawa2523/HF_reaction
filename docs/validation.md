# 実計算による検証の記録

計算基盤の改良(`docs/reviews/2026-09-27_platform_review.md`、以下レビュー)の M0〜M5・S-A〜S-D・FINAL を終えたコードで、レビュー §6.1 の検証セットを 1 回ずつ流した記録である。改良の途中の段ごとの実測、v3〜v5 と W7 の記録は git の履歴(`git show 496fad0:docs/validation.md`)にある。

## 1. 条件

- コード: `496fad0` と、この検証で直した対称数の修正(§4)。修正は S2 の後に入れたので、S20 CH3O• 以降の run は修正後のコードで流し、それより前の順位の付く run は修正後のコードで thermo だけを取り直した(§4)。
- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0 と pymsym)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(Te・I は def2-ECP)、低レベルは GFN2-xTB。ΔG は 298.15 K、1 atm(S19 を除く)。順位の量は δG_eff。
- run: WSL の `/home/user/hfauto_r6/VAL/runs/<run>`(すべて新しい run dir を 1 本ずつ直列)。スクリプトは同じ場所の `chain.sh`・`runcase.sh`・`ho_harness.py`(S11)。configs の外の入力(S5・S9 の多重度を外した system、S19 の条件、S20 の CCSD(T) 入りのパネル)は `/home/user/hfauto_r6/VAL/outside`。
- 時間は `/usr/bin/time` の壁時計。上限は `timeout` の 3,600 s(S16・S19 は 5,400 s)で、打ち切った run はない。smoke から S16 までの直列の合計は 4:11(S16 1:31、S19 0:31、S10 0:22)。
- 合否: レビュー §6 のとおり、終了コード、電子数と Level の照合、⟨S²⟩、outcome、既知の数値(以前の段の値)で判定する。文献値は PBE0/SVPD の誤差を含むので参考にとどめる。終了コード 1 は「failed の artifact がある」の意味で、一部の seed の失敗を含む S5・S6 の discover ではそれ自体を不合格にしない。

## 2. 実計算で通した範囲

| 項目 | DFT(PBE0/SVPD)まで通したもの | 低レベル(xTB・CREST・ReaDuct)だけ |
|---|---|---|
| 元素 | H、C、N、O、F、Cl、I(ECP)、Te(ECP) | S(SO2·NMe3)、Fe(FeCl3·CH4) |
| 電荷 | 0、−1(Cl⁻・I⁻ の SN2) | — |
| スピン | 一重項、二重項(CH3O•、H + H2、CH3•)、三重項(O2 の極小) | 二重項の錯体(CH3·O2、OH·CH4)、六重項(Fe、宣言が必須) |
| 反応型 | 1,2-H 移動(HCN、CH3O•)、ねじれ異性化(HONO)、1,3-H 移動(HONO → HNO2)、縮退転位(NH3 反転、H2O2 と DME のねじれ、マロンアルデヒドと acac の PT、Cl・I の恒等 SN2、H + H2 と NH3·HF の H 交換、(HF)₂ の入れ替え)、会合(TMA·(HF)₂、同一 basin) | H 引き抜き(OH + CH4)、ハロゲン結合(NH3···ICl)、配位付加体(SO2·NMe3)、ラジカル会合(CH3 + O2) |
| エネルギー層 | PBE0/TZVPD、ωB97X-D3/TZVPD、閉殻 CCSD(T)、ROHF-CCSD(T)(いずれも def2-TZVPD の SP) | — |

通していないもの: 陽イオン、−2 以下の陰イオン、遷移金属の DFT、開殻一重項(broken-symmetry は対象外)、溶媒(気相だけが対象)。

## 3. 結果

判定の内訳: R1〜R4 はすべて合格。S 系 20 のうち 17 が合格、S6 は条件付き合格(引き抜きは緩和としてだけ出る)、S15 と S17 は期待した分岐(多段、BARRIERLESS)に入らなかった(どちらも結果は化学的に正しい)。S2 で見つけた対称数の不具合は直した(§4)。

### 3.1 R1〜R4(基準値)

| # | 系 | 期待(v5 と M2 以後の値) | 実測 | 時間 | 判定 |
|---|---|---|---|---|---|
| R1 | HCN → HNC(`hcn`) | `elementary_step`、ΔG‡ 42.1783 | `elementary_step`(TS −1128.5i)、ΔE‡ 46.621、δG_eff 42.1783 | 1:17 | 合格 |
| R2 | HONO trans → cis(`hono`) | `elementary_step`、12.2265 − RT ln 2 ≈ 11.82 | `elementary_step`(−681.6i)、ΔE‡ 13.671、δG_eff 11.8158(TS はキラルで m = 2) | 3:03 | 合格 |
| R3 | NH3 反転(`nh3_inversion`) | `degenerate_rearrangement`、3.8395 | `degenerate_rearrangement`(−757.7i)、ΔE‡ 4.260、δG_eff 3.8395 | 1:15 | 合格 |
| R4 | 水(`water_same_basin`) | `same_basin` | `same_basin`、順位なし(blocker `outcome:same_basin`) | 0:11 | 合格 |

### 3.2 known_endpoints の系

| # | 系 | 型 | 期待 | 実測(kcal/mol) | 時間 | 判定 |
|---|---|---|---|---|---|---|
| S1 | Cl⁻ + CH3Cl(`sn2_cl`) | 恒等 SN2、単原子の単量体 | rc 0、degenerate、ΔE‡ 10.45、`dG_act_vs_separated` が出る(参考 −0.74) | rc 0、`degenerate_rearrangement`(−393.6i)、ΔE‡ 10.446、ΔG‡ 11.077、ΔG_assoc −5.08、`dG_act_vs_separated` +6.00 | 7:30 | 合格。解離極限基準の値は 1 atm の G で会合のエントロピー損失を含み、参考値とは基準が違う |
| S2 | I⁻ + CH3I(`sn2_i`) | 恒等 SN2、ECP | rc 0、α 電子数が ECP と整合、Level 一致 | 全 NWChem ジョブで `I library def2-ecp`・α 30、`degenerate_rearrangement`(−307.8i)、ΔE‡ 6.553。ΔG‡ は 7.533 と出て、対称数の修正後に 6.882(§4) | 13:19 | 合格(対称数の不具合を修正) |
| S3 | H2Te(`h2te`) | 極小のみ、Te | structures → dft → thermo が通る | 通過(`Te library def2-ecp`、E −269.20860 Eh) | 0:17 | 合格 |
| S4 | CH3O• → CH2OH•(`ch3o_doublet`、生成物は Cs) | 二重項の 1,2-H 移動 | 宣言反応が elementary、ΔG‡ ≈ 30.6、鏡像化の偽反応がない | 宣言反応 `elementary_step`(−2015.5i)、ΔE‡ 33.012、δG_eff 30.612、⟨S²⟩ ≤ 0.76。Cs の CH2OH• は DFT で saddle になり、mode-follow が CH2OH• のねじれの縮退 case(−533.3i、δG_eff 3.878、torsional)を加えた | 14:30 | 合格。加わった case は実在のねじれで偽反応ではない |
| S12 | DME の粗い入力(`dme_rough`) | 縮退 | `degenerate_rearrangement`、ΔE‡ ≈ 2.3 | `degenerate_rearrangement`(−226.9i)、ΔE‡ 2.330、δG_eff 1.964 | 10:04 | 合格 |
| S13 | H2O2 g+ → g−(`h2o2_gauche`) | ねじれ、低障壁 | 接続が確定、ΔE‡ ≈ 1.15、degenerate | `degenerate_rearrangement`(−290.5i、NEB の HEI から、DFT Hessian の種)、ΔE‡ 1.147、δG_eff 1.180 | 1:38 | 合格 |
| S14 | マロンアルデヒドの PT(`malonaldehyde`) | ZPE で沈む障壁 | note `submerged_barrier`、δG_eff は max 規則 | `degenerate_rearrangement`(−1100.2i)、ΔE‡ 2.011、ΔZPE‡ −2.402、`submerged_barrier`、δG_eff 0.0(表示用の ΔG‡ 0.250) | 16:52 | 合格 |
| S15 | trans-HONO → HNO2(`hono_hno2`) | 多段を想定 | `multi_step` と子反応 2 本 | `elementary_step`(直接の 1,3-H 移動の TS −2088.1i)、ΔE‡ 55.316、δG_eff 52.016。cis の井戸は経路に現れず分割は起きない | 2:46 | 期待の分岐は未達。直接の TS は化学的に正しい。中間体・分割・子 case は QM なしのテストだけ |
| S16 | acac の PT(`acac`、15 原子) | 縮退 PT、傍観者の回転 | xTB Hessian 採用、DFT Hessian 0 本、PT の TS | `degenerate_rearrangement`(−986.8i、`soft_secondary_mode`)。種の Hessian は xTB(方向との重なり 1.00)。ΔE‡ 1.980、ΔZPE‡ −2.607、ΔG‡ −0.93 で `submerged_barrier`、δG_eff 0.0。dft 1,334 s、paths 4,123 s | 1:30:57 | 合格 |
| S17 | (HF)₂ の入れ替え(`hf_dimer_swap`) | 解像度未満の障壁 | barrierless、δG_eff = max(ΔG_rxn, 0) | `degenerate_rearrangement`(C2h の TS −227.3i)、ΔE‡ 1.256、δG_eff 1.392 | 1:35 | 期待の分岐は未達。障壁 1.26 が解像度 1.0 を超えるので BARRIERLESS にならないのが正しい |
| S18 | H + H2(`h3_doublet`) | 二重項の引き抜き | unresolved でなく TS が得られる | `degenerate_rearrangement`(−605.7i)、⟨S²⟩ 0.750〜0.768、ΔE‡ 4.563、δG_eff 4.091(共線 H3 の真の障壁 約 9.6 は参考。PBE0 の過小評価) | 0:48 | 合格 |

### 3.3 discover の系(explore まで、または全体)

| # | 系 | 型 | 期待 | 実測 | 時間 | 判定 |
|---|---|---|---|---|---|---|
| S5 | CH3• + O2(`ch3_o2`) | 二重項の会合 | 多重度の宣言なしは INPUT_INVALID、m = 2 で UKS mult 2 | 宣言なし: conformers が `input_invalid`(`declare_multiplicity: candidates (2, 4)`)。宣言 2: CREST の argv に `--uhf 1`(rc −11 で失敗、開殻の組成の既知の挙動)、placement の 6 seed のうち 1 つが xTB で CH3OO•(1 断片)に結合し、5 つは xTB の収束失敗。`known_endpoints --to dft` で CH3• は m 2・⟨S²⟩ 0.7544、O2 は m 3・⟨S²⟩ 2.0096 | 0:02 / 0:04 / 0:44 | 合格 |
| S6 | OH···CH4(`oh_ch4`) | 二重項の H 引き抜き | 引き抜きの生成物が出る | 引き抜きは GFN2 で障壁がなく、placement の seed が screen で CH3···H2O に落ちる(緩和の発見 `collapsed_to:CH3+H2O` 2 件)。CH3···H2O を出発点に 11 attempt(AFIR 3)、生成物は CH3OH + H の 2 件(ΔE‡ 31.6 / 37.8)。seed 4 本は xTB の失敗 | 0:07 | 条件付き合格。引き抜きの生成物は緩和として出るが、TS 付きの発見にはならない |
| S7 | NH3···ICl(`nh3_icl`) | ハロゲン結合 | 2 断片、無意味な polar_h がない | 錯体 6 配座はすべて 2 断片の状態(ClI+H3N)、polar_h の drive なし。4 attempt で生成物 0 | 0:10 | 合格 |
| S8 | SO2·NMe3(`so2_nme3`) | 配位付加体 | CREGEN が付加体を捨てない | 付加体 2 配座が CREGEN を通り、screen 極小(ラベルは 2 断片)として残る。24 attempt で生成物 0 | 1:04 | 合格 |
| S9 | FeCl3·CH4(`fecl3_ch4`、conformers まで) | Fe、剛体配置 | 例外なし、多重度未宣言は INPUT_INVALID | 宣言なし: `species_fecl3` が `input_invalid`(`declare_multiplicity: d-block element(s) Fe`)。宣言 6: 例外なく CREST 6 配座(`CH4Cl3Fe_q0_m6`) | 0:01 / 0:30 | 合格 |
| S10 | NH3·HF、TMA·HF(`amine_pilot2`、全体) | 二重 H 交換(縮退) | −1265i の drive が 1 件、product になり縮退 case が走る | explore 25 attempt(AFIR 4)、NH3·HF の二重 H 交換(xTB −1265.3i)が生成物 1 件。発見の case が `degenerate_rearrangement`(DFT −1081.3i)、ΔE‡ 32.888、δG_eff 32.289、ΔG_assoc −5.66 | 22:24 | 合格 |
| S11 | DME の C2v 重なり形の種(`ho_harness.py`) | 負モード 2 本 | 一次鞍点 約 226i、ΔE‡ ≈ 2.3、停滞が再現しない | 1 回目(xTB Hessian、重なり 0.30)は 2 次の saddle(−232.1i / −128.6i)→ `higher_order`。−128.6i に 0.23 Å 押し、TS freq を Hessian に 2 回目で一次鞍点 −226.8i、ΔE‡ 2.33。`max_saddle_attempts` 2 の中で終わる | 7:03 | 合格 |
| S19 | TMA·(HF)₂(`tma_hf2`)と条件の追記 | 会合、複数条件 | 反応の δG_eff は標準状態で不変、ΔG_assoc だけが動く | `discover`: explore 15 attempt(すべて NT2)と緩和の発見 3、生成物 0。宣言反応 `neutral_to_shared_proton` は DFT で shared_proton が neutral に落ちて `same_basin`(順位なし)、ΔG_assoc −11.53。追記した thermo → report で ΔG_assoc は 1 atm で −14.52 / −11.53 / −5.11、1 M で −17.52 / −15.32 / −10.66(250 / 298.15 / 400 K、差は 2RT ln(RT/P°))。反応の δG_eff の不変性は S1 の複製への同じ追記で確かめた(250 / 298.15 / 400 K の δG_eff 10.898 / 11.077 / 11.435 が 1 atm と 1 M で同じ。ΔG_assoc は 1 M で 1.50 / 1.90 / 2.78 低い) | 30:38 / 0:01 | 合格 |

### 3.4 手法パネル(S20、known_endpoints の run の複製に CCSD(T) 入りの method_panel を追記)

ΔE‡(kcal/mol)。δG_eff は ωB97X-D3/TZVPD の層。panel の min / max は method_panel.csv の最小・最大と一致した。

| 系 | PBE0/SVPD | PBE0/TZVPD | ωB97X-D3/TZVPD | CCSD(T)/TZVPD | δG_eff | 時間 | 判定 |
|---|---|---|---|---|---|---|---|
| HCN → HNC | 46.621 | 46.429 | 46.369 | 47.807 | 41.926 | 1:24 | 合格 |
| Cl⁻ + CH3Cl | 10.446 | 11.146 | 15.101 | 13.436 | 15.732(`dG_act_vs_separated` 11.43) | 11:14 | 合格 |
| CH3O• → CH2OH• | 33.012 | 32.713 | 33.578 | 33.372(ROHF-CCSD(T)) | 31.178 | 6:19 | 合格 |
| CH2OH• のねじれ | 4.482 | 4.665 | 4.818 | 4.836(ROHF-CCSD(T)) | 4.213 | (同上) | 合格 |

### 3.5 smoke と golden

- 実エンジンの smoke(`HFAUTO_REAL=1 pytest -m real tests/smoke`): 13 件合格(64 s)。水の opt と別ジョブの freq、NH3 の saddle と 2 本目の負モードの追跡、I の def2-ECP と Cl⁻、閉殻 CCSD(T) と ωB97X-D3、ROHF-CCSD(T)、HCN の string 1 チャンク(STO-3G)、ReaDuct の NT2、xTB、pysis NEB、CREST(電荷・スピン、陰イオン)、GoodVibes の CLI との一致。
- golden(レビュー §6.2)は `tests/golden` にあり、既定のテストで毎回通る: OH•(UKS、⟨S²⟩ 0.7526)の opt・freq、FHF⁻ の opt・freq・SP、HI の SP(α 13、E −298.3091 Eh)と CCSD(T) の凍結軌道、Cl⁻ の xTB と NWChem、OH• の ROHF-CCSD(T)、I⁻···CH3I の対称数(§4)。
- QM なしの fixture(レビュー §6.3): Ar の S°(154.850 J/mol/K、JANAF 154.846)、δG_eff の場合分け、経路の分類器のモデル曲線、`shape_hessian`、結合判定の距離、多重度の規則、鏡映を含む同一性。

### 3.6 到達表

`log.jsonl` の decide の reason と job.json から数えた(括弧は通った run の数)。

| 対象 | 実計算で通ったもの | 通らないもの(QM なしのテストで確認) |
|---|---|---|
| 決定表の行 | 2 same_basin(3)、4 connection(17)、9 ts_validated・10 saddle_converged・12 screen・14 seed(各 17) | 1・3・5〜8・11・13・15〜17(17 行すべての reason を `tests/unit/drivers/test_reaction_case_decide.py` が確かめる) |
| saddle の種と Hessian | 種 `screen_ts` 14・`screen_hei` 2・`discovery_ts` 3、S11 の `higher_order_retry`。Hessian は xTB 15・DFT 2・検証済み TS freq(S11) | 種 `path_hei`・`saddle_restart` |
| outcome | `elementary_step` 6、`degenerate_rearrangement` 13、`same_basin` 3 | `multi_step`、`barrierless`、`unresolved_within_budget` |
| エンジン | CREST、xTB opt/freq、ReaDuct の NT2/AFIR、pysis NEB、NWChem の opt/freq/saddle/SP(PBE0 SVPD/TZVPD、ωB97X-D3、閉殻 CCSD(T)、ROHF-CCSD(T))。GoodVibes は同一プロセス | NWChem string(FIND_PATH)は smoke の 1 チャンクだけ |

## 4. この検証で直した不具合

- **対称数の揺らぎ(S2)**: 同じ入力の I⁻···CH3I 錯体(C3v)が、FINAL の run では C3v(σ 3)、この run では Cs(σ 1)と判定され、ΔG‡ が RT ln 3 = 0.65 kcal/mol 違った(6.882 → 7.533)。NWChem の最適化の終点で I⁻ が C3 軸から 0.03° ずれただけで、振動数は 0.1 cm⁻¹ 以内で同じだった。GoodVibes が呼ぶ libmsym(pymsym)の等価判定のしきい値(相対 5e-4)が、最適化の収束の幅より狭いのが原因である。
- 修正: `chemistry/thermo.symmetry_number` が libmsym の等価しきい値を 2e-3 にして点群を求め、σ を GoodVibes に渡す(`symm=False`)。FINAL と VAL の NWChem freq の構造 126 本(21 組成)で、しきい値 2e-3 が既定と違う点群を出すのは S2 のこの構造と H + H2 の TS(Cs から「群なし」、どちらも σ 1)だけだった。S2 の構造はしきい値 5.5e-4 で既に C3v になり、FINAL の構造はしきい値 1e-2 まで σ が変わらない。回帰テストは `tests/golden/test_goodvibes_thermo.py::test_symmetry_number_tolerates_an_optimised_c3v_complex`。
- 修正前に流した run(R1〜R4、S1〜S4、S13、S15、S17、S18)と S12・S14 の複製に、`--from thermo` で thermo と report を取り直した。14 run の複製(`rt_<run>`)の結果: S2 の錯体の G だけが変わり(−635.36312 → −635.36208 Eh)、ΔG‡ は 7.533 → 6.882(FINAL の run の 6.882 と一致)。ほかの 13 run は反応と化学種の thermo が 1e-9 以内で同じ。S20 の HCN・SN2 のパネルの構造も σ は変わらない(上の 126 本に含む)。

## 5. 改良中の実計算プローブで決めたこと

| 案 | プローブ | 決定 |
|---|---|---|
| U8-P3 エネルギー層 | 5 系で PBE0/SVPD・ωB97X-D3/TZVPD・CCSD(T) の ΔE‡ を比較。ωB97X-D3 は CCSD(T) との平均絶対誤差 1.29(PBE0/SVPD 2.07)だが、17 原子で 1 点 29 分 | 既定の pipeline は PBE0/SVPD のまま。層は method_panel の追記で使う |
| U5-P2 SCREEN | GS と固定端 CI-NEB を 11 系で比較。NEB は全系で収束し、(HF)₂ は未解決から C2h の TS に、acac は完走し、HONO と H2O2 は FIND_PATH に入らない | 固定端 CI-NEB を採用 |
| U0-P5 SCF の救済 | 収束しない H3 の bead で、反復増・damping・level shift・rabuck は未収束、`cgmin` だけが 5〜8 反復で収束 | 救済は `cgmin` の 1 本。cgmin は ⟨S²⟩ を出さないので、⟨S²⟩ のない開殻 DFT は `incomplete_output` |
| U8-P7 開殻 CCSD(T) | CH3O• の TCE ROHF-CCSD(T) は `2eorb` なしで GA の確保に失敗、`2eorb 2emet 13` で 63 s | TCE に `2eorb 2emet 13` |
| U2-P2・U2-P6 CREST | `--noopt` 案と初期最適化案、円錐と剛体の配置を 8 組成で比較 | `--noopt` と全原子の `--notopo`、円錐の配置を残す(開殻の組成で CREST は失敗し seed を使う) |
| U7-P5 GoodVibes | 同一プロセスの呼び出しで既存 run の G・H・ZPE を 1.1e-13 Eh 以内で再現 | 同一プロセスで呼ぶ |
| U9-P7 外部プロセス | `timeout` で止めた `mpirun nwchem` が signal handler ありで残らない | handler を入れた |

## 6. 撤回する値

- HONO trans → cis の ΔG‡ **24.11 kcal/mol**(旧 `hono_isomerization_v3`)。ZPE の二重計上と basin 判定の回帰による。正しくは ΔE‡ 13.67、ΔG‡ 12.23(m = 1)、m = 2 を数えた今の δG_eff は 11.82。
- HCN → HNC の ΔG_rxn **18.05 kcal/mol**(内部フォールバックの熱化学)。D3BJ の実測は 12.68。
- I⁻ + CH3I の ΔG‡ **7.53 kcal/mol**(この検証の最初の S2、対称数の誤判定)。正しくは 6.88(§4)。
- M2 以前の R2 の 12.2265(v5)は、TS のキラリティ(m = 2)を数えていない。比べるときは 11.82 を使う。

## 7. 残る課題

1. **実計算で通らない決定表の行**: 1・3・5〜8・11・13・15〜17(QM なしの `tests/unit/drivers/test_reaction_case_decide.py` だけ)。中でも多段(S15 は直接の TS が見つかる)、BARRIERLESS(S17 は障壁が解像度を超える)、FIND_PATH(smoke の 1 チャンクだけ)は実計算の系がない。これらを通す系(cis の井戸を必ず経る多段反応、障壁 1 kcal/mol 未満の入れ替え)は見つかっていない。
2. **開殻の組成の配座探索**: CREST は開殻の組成で必ず失敗し(rc −11)、placement の seed も xTB の SCF・最適化で半数以上が落ちる(S5 5/6、S6 4/5)。S6 の引き抜きは GFN2 で障壁がなく、TS 付きの発見にならない。
3. **範囲外**: 陽イオン、遷移金属の DFT、開殻一重項、溶媒。宣言した範囲(README)どおり。
4. **精度**: PBE0/SVPD の ΔE‡ は CCSD(T) から最大 3 kcal/mol ずれる(SN2 −3.0、HONO +2.0)。高精度が要る反応は method_panel を追記する(17 原子の ωB97X-D3/TZVPD は 1 点約 30 分)。H + H2 の障壁(4.6)は PBE0 の既知の過小評価。
5. **時間**: 15 原子の acac は 1:31(paths 4,123 s の大半が QRC と TS freq)。17 原子の DFT freq は 1 本約 15 分で、discover の dft stage の大半を占める。
6. **対称数のしきい値**: 2e-3 は検証セットの 126 構造で確かめただけである。より大きく柔らかい錯体では、最適化の収束の幅がしきい値を超えうる。
7. **S7**: M4 で得た NH2I + HCl(NT2)は、S-C で状態のまとめ方を変えてから出発点に入らず、再現しない。
