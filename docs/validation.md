# 実計算による検証(r6 の基準と検証セット、v3・v4・v5 の記録)

## r6: 計算基盤の抜本改良

`docs/reviews/2026-09-27_platform_review.md`(以下レビュー)の M0〜M5・S-A〜S-D を段ごとに入れ、各段で下の基準値とレビュー §6 の検証セットを確かめる。run は WSL の `/home/user/hfauto_r6/<段>/<run>`(すべて新しい run dir)に置く。環境と停留点の手法(`pbe0-d3bj_def2-svpd`)は v5 と同じ。

**基準値(v5、コード `13a4340`、pipeline `known_endpoints`)**

| # | 系(config) | v5 の結果 | 所要時間 | r6 で期待する変化 |
|---|---|---|---|---|
| R1 | HCN→HNC(`hcn`) | `elementary_step`、ΔG‡ 42.1783 kcal/mol | 1:22 | なし |
| R2 | HONO trans→cis(`hono`) | `elementary_step`、ΔG‡ 12.2265 kcal/mol | 5:32 | M2(U3-P2)の後、RT ln 2 だけ下がる(約 11.82) |
| R3 | NH3 反転(`nh3_inversion`) | `degenerate_rearrangement`、ΔG‡ 3.8395 kcal/mol | 1:25 | なし |
| R4 | 水(`water_same_basin`) | `same_basin` | 0:11 | なし |

**検証セット(レビュー §6.1)と最初に確かめる段**

| # | 系 | config(`configs/systems/`) | 最初の段 |
|---|---|---|---|
| R1〜R4 | 上の基準値 | `hcn`、`hono`、`nh3_inversion`、`water_same_basin` | M0(以後の各段) |
| S1 | Cl⁻ + CH3Cl(Cl⁻ を単量体に含む) | `sn2_cl` | M0(完走と report)、M1(合否) |
| S2 | I⁻ + CH3I | `sn2_i` | M0(完走と report)、M1(合否) |
| S3 | H2Te | `h2te` | M0(完走と report)、M1(合否) |
| S4 | CH3O• → CH2OH•(生成物は Cs) | `ch3o_doublet` | M2 |
| S5 | CH3• + O2(m=2 を宣言) | `ch3_o2` | M1 |
| S6 | OH• + CH4 | `oh_ch4` | M4 |
| S7 | NH3···ICl | `nh3_icl` | M4 |
| S8 | SO2·NMe3 | `so2_nme3` | S-C |
| S9 | FeCl3·CH4(Fe を m=6 と宣言) | `fecl3_ch4` | M1 |
| S10 | NH3·HF、TMA·HF | `amine_pilot2` | M4 |
| S11 | DME の C2v 重なり形 | config なし(`/home/user/hfauto_probe_r5/ho_harness.py`) | S-B |
| S12 | DME(粗い入力) | `dme_rough` | M2 |
| S13 | H2O2 gauche(+)→gauche(−) | `h2o2_gauche` | M2 |
| S14 | マロンアルデヒドの PT | `malonaldehyde` | M5 |
| S15 | trans-HONO → HNO2 | `hono_hno2` | M2 |
| S16 | アセチルアセトンの PT | `acac` | M3 |
| S17 | (HF)₂ の供与体・受容体の入れ替え | `hf_dimer_swap` | M2 |
| S18 | H + H2(共線の vdW 錯体) | `h3_doublet` | S-A |
| S19・S20 | 複数条件、手法パネル | S-A で pipeline に追記 | S-A |

- 入力の出所: S1〜S4・S12〜S14・S16・S18 は `/home/user/hfauto_probe_r5/systems` の xyz、S10 は W7 の `amine_pilot2`。S15 の HNO2 は `[H][N+](=O)[O-]` の RDKit ETKDG 構造を GFN2-xTB で最適化し、trans-HONO と同じ原子順(H O N O)に並べて重ねた。S17 は文献に近い trans 屈曲の (HF)₂(F···F 2.72 Å)と、2 つの HF の座標を入れ替えた構造(原子順は同じ)。GFN2-xTB では (HF)₂ がほぼ直線になるので、最適化していない。
- golden(レビュー §6.2): G25 は OH•(二重項、UKS)の opt と freq、G26 は FHF⁻ の opt・freq・SP。今のコードの `NWChemEngine` で生成した(`/home/user/hfauto_r6/m0/golden`)。G25 の ⟨S²⟩ は 0.7526、振動は 1 本(3N−5)。opt と freq のエネルギーの差は 9.6e-7 Eh で、U3-P7 の一致判定に使う。

**M0 の実測(U9-P1、`/home/user/hfauto_r6/M0/runs/<系>_known_endpoints`)**

| 系 | 結果 | 終了コード | 所要時間 | ジョブ(ヒット) |
|---|---|---|---|---|
| S3 `h2te` | structures で `species_h2te` が失敗(Te の共有結合半径がない)。dft・paths・thermo は 0 artifact で done、report を出力 | 1 | 0:02 | 0 |
| S2 `sn2_i` | dft の 2 極小が `scf`(ECP なし、SCF 未収束)。paths・thermo は 0 artifact で done、report を出力 | 1 | 2:43 | 4(0) |
| S1 `sn2_cl` | Cl⁻ の極小が `incomplete_output`(`parse:ValidationError: ...`)。`sn2_identity` は `degenerate_rearrangement`、ΔE‡ 10.446、ΔG‡ 11.077 kcal/mol | 1 | 9:11 | 34(1) |
| R1 `hcn` | `elementary_step`、ΔG‡ 42.1783 | 0 | 1:24 | 29(2) |
| R2 `hono` | `elementary_step`、ΔG‡ 12.2275(v5 比 +0.0010。GS の失敗から FIND_PATH の経路は同じで、ΔE‡ 13.671。v4→v5 の −0.0002 と同じ string の揺らぎ) | 0 | 5:40 | 16(0) |
| R3 `nh3_inversion` | `degenerate_rearrangement`、ΔG‡ 3.8395 | 0 | 1:27 | 27(1) |
| R4 `water_same_basin` | `same_basin` | 0 | 0:11 | 4(0) |

- 7 run ともトレースバックなしで終わり、report.html・ranking.csv・coverage.csv がある。入力のない stage の warning は stderr に出る。実エンジンの smoke(`-m real tests/smoke`)は 12 件合格(71 秒)。

**M1 の実測(U1-P1・U1-P3・U0-P1・U3-P1、`/home/user/hfauto_r6/M1/runs/<run>`)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| golden・smoke | WSL の既定 suite 351 件、`-m real tests/smoke` 13 件が合格。G27 は HI の SP(PBE0-D3BJ/def2-SVPD、I に def2-ECP、alpha 13、E −298.3091 Eh)と CCSD(T)/def2-TZVPD(`freeze atomic` は ECP 原子の軌道を凍結せず `number of core 0`)。G28 は Cl⁻ の xTB opt/hess と NWChem opt/freq(n_external 3、振動なし)。Ar の S(298.15 K、1 bar)は 154.850 J/mol/K(JANAF 154.846) | 1:17 |
| S1 `sn2_cl` | Cl⁻ の極小を登録。`sn2_identity` は `degenerate_rearrangement`、ΔE‡ 10.446、ΔG‡ 11.077、ΔG_assoc −5.08、ΔG‡(解離極限基準)+6.00 kcal/mol(1 atm の G で、会合のエントロピー損失を含む。レビューの参考値 −0.74 とは基準が違うとみられ、合否の基準ではない) | 9:25 |
| S2 `sn2_i` | 全ジョブ alpha 30(`I library def2-ecp`)、観測 Level は一致。`degenerate_rearrangement`、ΔE‡ 6.55、ΔG‡ 6.88 kcal/mol(参考の中央障壁 約 8) | 14:23 |
| S3 `h2te` | structures → dft → thermo を通過(`Te library def2-ecp`、G −269.2188 Eh) | 0:17 |
| S5 `ch3_o2` | 組成の多重度を外した複製は conformers で `INPUT_INVALID`(`declare_multiplicity: candidates (2, 4)`)。宣言 2 では `discover --to screen` の CREST の argv に `--uhf 1`、組成 id `CH3O2_q0_m2`(CREST は TRIALMD で signal 11 になり、placement の 6 seed が出た。開殻の組成の既知の挙動)。`known_endpoints --to dft` で CH3 の極小は UKS(`odft`)、Level の多重度 2、⟨S²⟩ 0.7544 | 0:04 / 0:43 |
| S9 `fecl3_ch4` | 多重度を外した複製は `species_fecl3` が `INPUT_INVALID`(`declare_multiplicity: d-block element(s) Fe`)。宣言 6 では `--to conformers` が例外なく終わり、CREST で 6 構造(`CH4Cl3Fe_q0_m6`) | 0:30 |
| R1 `hcn` | `elementary_step`、ΔG‡ 42.1783 | 1:26 |
| R2 `hono` | `elementary_step`、ΔG‡ 12.2263(v5 比 −0.0002、string の揺らぎの範囲)、ΔE‡ 13.671 | 5:51 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、ΔG‡ 3.8395 | 1:27 |
| R4 `water_same_basin` | `same_basin` | 0:12 |

**M2 の実測(U3-P2・U5-P6・U5-P1・U5-P3、`/home/user/hfauto_r6/M2/runs/<run>`、pipeline `known_endpoints`)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| S4 `ch3o_doublet` | 宣言した `ch3o_to_ch2oh` が `elementary_step`(TS −2015.5i)、ΔE‡ 33.01、ΔG‡ 30.61 kcal/mol。順位は 1 行だけで、鏡像化の擬似反応はない | 8:53 |
| S12 `dme_rough` | 粗い入力でも basin 構造で判定して `degenerate_rearrangement`、ΔE‡ 2.33、ΔG‡ 1.96 kcal/mol | 9:50 |
| S13 `h2o2_gauche` | g+ → g− は 1 つの basin の縮退反応(degenerate=True、`same_basin` にならない)。SCREEN は single(1.28 kcal/mol)、鞍点 −290.5i、QRC の接続は失敗して `unresolved_within_budget`(M3 で扱う) | 1:58 |
| S15 `hono_hno2` | SCREEN の IDPP が single(85.1 kcal/mol)で、その HEI から trans-HONO と HNO2 を直接つなぐ TS(−2085.6i)に到達し `elementary_step`、ΔG‡ 52.06 kcal/mol。cis の井戸は経路に現れず、分割は起きなかった | 4:52 |
| S17 `hf_dimer_swap` | 2 つの入力は 1 つの basin(degenerate=True)。ただし写像付きで整列した IDPP・string は H が同じ側に回る経路(最大 7.8 / 7.2 kcal/mol)で single になり、鞍点 2 回が maxiter で `unresolved_within_budget`。C2h の入れ替え経路を初期経路が通らないためで、M1 のコードでも経路は同じ(分類器の問題ではない) | 7:11 |
| R1 `hcn` | `elementary_step`、ΔG‡ 42.1783 | 1:20 |
| R2 `hono` | `elementary_step`、ΔG‡ 11.8158 = 12.2263 − RT ln 2(TS は C1 でキラル、m = 2)。ΔE‡ 13.671 は同じ | 5:50 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、ΔG‡ 3.8395。paths のジョブは 22 → 18(SCREEN の DFT freq、両端ノードの SP 2 本、SCREEN の TS の SP がなくなった) | 1:17 |
| R4 `water_same_basin` | `same_basin` | 0:11 |

**M3 の実測(U6-P2・U6-P1・U6-P5、`/home/user/hfauto_r6/M3/runs/<run>`、pipeline `known_endpoints`)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| S13 `h2o2_gauche` | 鞍点 −290.5i から QRC が両側とも同じ basin の鏡像へ下り `degenerate_rearrangement`、ΔE‡ 1.147、ΔG‡ 1.180 kcal/mol(M2 は `connection_failed`)。xTB の種に負モードがなく DFT Hessian を使った | 1:37 |
| S16 `acac` | SCREEN の xTB TS の虚モードを方向にして xTB Hessian を採用(`saddle_hessian:xtb:overlap:1.00`、種の DFT Hessian 0 本)。NWChem の 1 歩目は `negative= 1`、6 歩(226 s)で PT の TS(−986.75i、ほかに −38 / −25 のメチル、`soft_secondary_mode`)に収束し、開始から約 40 分で TS 検証まで通過(以前は誤ったモードを 36 歩追って 60 分打ち切り)。QRC は 1 側 31 歩(1,390 s)かかり、60 分の上限では接続判定まで届かない(1 側目の終点 −345.191570 Eh は反応物極小 −345.191565 Eh と同じ高さ) | 60:00 で打ち切り(再開 45:00 で打ち切り) |
| S4 `ch3o_doublet` | M2 の run を複製して `--from paths`: 新しいジョブキーの saddle 1 本だけを再計算(hits 17 / misses 1)、6 歩で同じ TS、`elementary_step`、ΔG‡ 30.6117 | 0:26 |
| R1 `hcn` | `elementary_step`、ΔG‡ 42.1783、saddle 5 歩(xTB Hessian、重なり 1.00) | 1:20 |
| R2 `hono` | `elementary_step`、ΔG‡ 11.8160(string の揺らぎ ±0.002 以内)、saddle 4 歩(xTB Hessian、重なり 0.57) | 5:34 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、ΔG‡ 3.8395、saddle 3 歩(xTB Hessian、重なり 1.00) | 1:16 |
| R4 `water_same_basin` | `same_basin` | 0:12 |

**M4 の実測(U4-P2・U4-P1・U4-P5、`/home/user/hfauto_r6/M4/runs/<run>`、S 系は `discover --to explore`、S10 は `discover` 全体)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| S1 `sn2_cl` | 2 つの錯体とも trial 0 が背面 SN2 の T1(C–Cl′ 形成、C–Cl 切断)。NT2 で TS −431.8i(ΔE‡ 44.36 kJ/mol)、IRC の端は状態ラベルが source と同じで C–Cl′ に結合した恒等 SN2 の生成物として残る(以前は same_as_source)。16 attempt で生成物 6(すべて恒等 SN2。重複は U4-P7) | 0:27 |
| S6 `oh_ch4` | GFN2-xTB では引き抜きが障壁なしで、placement の seed は screen で CH3···H2O に崩壊する(`collapsed_to:CH3+H2O` の緩和の発見 2 件)。H2O + CH3• は screen 極小として残り、explore の出発点もこの CH3···H2O なので、13 attempt は逆向き(生成物 0)。引き抜きの T1 は単体テストの CH4···OH で先頭に出るが、ReaDuct の実計算(O···H 2.30 Å)では NT2 が `ts_not_converged`、AFIR が `same_as_source` | 0:14 |
| S7 `nh3_icl` | 2 断片の錯体 2 つで T1 のハロゲン移動(N–I 形成・I–Cl 切断)と N–Cl 形成・I–Cl 切断。polar_h の drive はない。NT2 で NH2I + HCl(−193.6i、124.4 kJ/mol)、AFIR で 2 生成物 | 0:17 |
| S10 `amine_pilot2` | NH3·HF の二重 H 交換(交換リレー、形成 (0,4),(1,5)、切断 (4,5),(0,1))は 1 attempt で −1265.3i の TS に到達し、生成物として残った(W7 は単独の移動 3 件がすべて same_as_source)。`paths` で発見の縮退 case が走り `degenerate_rearrangement`、DFT の TS −1081.2i、ΔE‡ 32.89、ΔG‡ 32.29 kcal/mol。TMA の 1,2-H 移動 ×9 の同じ drive は 1 件になり、376.37 kJ/mol の TS は W7 の 6 回から 3 回(1,2-H 移動 1 件と、同じ TS に落ちるメチル間のリレー 2 件)に減った。explore 全体は 37 attempt・生成物 1(W7 は 44・0) | 23:35 |
| R1 `hcn` | `elementary_step`、ΔG‡ 42.1783(停止した S10 の NWChem が約 40 秒重なった) | 1:54 |
| R2 `hono` | `elementary_step`、ΔG‡ 11.8158 | 3:22 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、ΔG‡ 3.8395 | 1:20 |
| R4 `water_same_basin` | `same_basin` | 0:12 |

- 実 run で見つかった列挙の穴: CREST の NH3·HF では F···H(N) が 3.03 Å で Σr_vdW(2.67 Å)を超え、レビューどおりの接触の規則では二重 H 交換のリレーが作られなかった。2 つ目の H が 1 つ目の受容原子から供与原子へ移る交換は、1 つ目の移動で両原子が近づくので 2 つ目の接触を問わないことにした(単体テストはこの CREST 構造)。
- M3・S6 の −289i の TMA·(HF)₂ の HF 交換(下の §5)は、M4 からは same_as_source ではなく生成物(縮退転位)になる。

**M5 の実測(U7-P1・U8-P2、`/home/user/hfauto_r6/M5/runs/<run>`、以前の run の複製を `--from paths`、プローブの複製は `--from structures`)**

M5 から順位の量は δG_eff、ranking.csv の見出しは `…,dG_eff_kcal,band_low_kcal,band_high_kcal,dG_act_kcal,dG_rxn_kcal,dG_act_vs_separated_kcal,torsional,blockers,notes`(手法パネルがあると `dE_act_panel_min_kcal,dE_act_panel_max_kcal` が付く)。blocker の `dzpe_out_of_tolerance` と `method_sign_disagreement` はなくなった。以下の §1〜§12 に出てくるこの 2 つは当時の記録である。

| 系 | 結果 | 所要時間 |
|---|---|---|
| S14 `malonaldehyde`(p4_malon の複製) | `degenerate_rearrangement`(TS −1100.4i)、ΔE‡ 2.011、ΔZPE‡ −2.401 で ΔE0‡ ≤ 0 → 注記 `submerged_barrier`、δG_eff = max(ΔG_rxn, 0) = 0.0(band 0.0〜0.0)で rank 1。表示用の ΔG‡ 0.250 は残る | 8:46 |
| S17 `hf_dimer_swap` | M2 と同じく初期経路が H の同じ側の回転を通り、鞍点が maxiter で `unresolved_within_budget`(順位なし)。障壁なしの行の順位付けは QM を使わないテストでだけ確認(S-A の固定端 CI-NEB 待ち) | 2:18 |
| S4 `ch3o_doublet` | `elementary_step`、δG_eff = ΔG‡ 30.6117、torsional 列は False(結合変化のある 1,2-H 移動) | 0:01 |
| 手法パネル HCN(R1 の複製に `method_panel` を追記) | ranking.csv の dE_act_panel_min / max 46.369 / 46.621 は method_panel.csv の最小・最大と一致(sign_disagreement 列はない) | 0:30 |
| 手法パネル SN2(sn2_panel の複製、既定のパネル) | dE_act_panel_min / max 10.446 / 15.101 が method_panel.csv と一致。δG_eff 11.077 で rank 1 | 0:27 |
| R1 `hcn` | `elementary_step`、δG_eff = ΔG‡ 42.1783 | 0:02 |
| R2 `hono` | `elementary_step`、δG_eff = ΔG‡ 11.8158(状態の最小 G は trans)、torsional True | 0:01 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff = ΔG‡ 3.8395 | 0:01 |
| R4 `water_same_basin` | `same_basin`、順位なし(blocker は `outcome:same_basin` だけ) | 0:01 |

**S-A の実測(SA-A: U8-P3・U8-P6・U8-P8、`/home/user/hfauto_r6/SA/runs/<run>`)**

U8-P3 のプローブ: M5 の run の複製に、sp(ωB97X-D3/def2-TZVPD と CCSD(T)/def2-TZVPD)→ thermo(`energy_method` ωB97X-D3)→ report を追記した(コードは `82c46d4` の複製)。ΔE‡(kcal/mol、括弧は CCSD(T) との差):

| 系 | PBE0-D3BJ/def2-SVPD | ωB97X-D3/def2-TZVPD | CCSD(T)/def2-TZVPD | δG_eff(SVPD → ωB97X 層) |
|---|---|---|---|---|
| HCN → HNC | 46.62(−1.19) | 46.37(−1.44) | 47.81 | 42.18 → 41.93 |
| HONO trans → cis | 13.67(+2.02) | 12.82(+1.17) | 11.65 | 11.82 → 10.97 |
| SN2 Cl⁻···CH3Cl | 10.45(−2.99) | 15.10(+1.67) | 13.44 | 11.08 → 15.73 |
| マロンアルデヒドの PT | 2.01(−2.1) | 3.23(−0.9) | (T) が失敗。CCSD 4.48、文献の CCSD(T) 約 4.1(Wang et al., J. Chem. Phys. 128, 224314 (2008))を基準にした | 0.00(`submerged_barrier`)→ 1.47(ΔE0‡ = +0.83 で沈まない) |
| CH3O• → CH2OH• | 33.01 | 33.58(⟨S²⟩ 0.754〜0.758) | 開殻は SA-C | 30.61 → 31.18 |

- SP 1 点の時間(4 rank): TMA·(HF)₂(17 原子、`tma_hf2/neutral.xyz`)で PBE0/SVPD 35.5 s、ωB97X-D3/TZVPD 1,758 s(freq 約 900 s の約 2 倍)。マロンアルデヒド(9 原子)の ωB97X-D3 は 577〜750 s、CCSD(T) は CCSD が約 1,070 s/点の後、(T) が GA の配列を確保できず abort した(`ga_create_JKblocked: cannot allocate`、1,200 MB/rank。重原子 5 個)。HCN・HONO・SN2・CH3O の ωB97X-D3 は 6〜115 s(HONO は同時負荷あり)。
- 判定: ωB97X-D3/TZVPD は 4 系中 3 系で CCSD(T) に近い(平均絶対誤差 2.07 → 1.29 kcal/mol。HCN だけ 0.25 遠い)が、17 原子で 1 点 29 分かかり(レビューの見込みは数分)、discover では極小が数十点あるので、**既定の discover・known_endpoints には入れない**。エネルギー層は method_panel の追記(panel_sp → panel_thermo → panel_report)で使う。
- U8-P6: M5 の hcn_panel・sn2_panel の複製で panel_report だけを新しいコードで再実行すると、ranking.csv は `energy_level` 列が加わるだけで、method_panel.csv は参照レベルの行だけが 1e-4 kcal/mol 以下変わった(参照エネルギーを opt / saddle から freq の SCF に変えたため。値は ReactionThermo の ΔE‡ と一致: HCN 46.62147、SN2 10.44594)。
- 新しいコード(作業ツリーの複製)で M5 の hcn_panel に既定の method_panel を追記: conditions を写さずに panel_sp → panel_thermo → panel_report が通り、ranking.csv の `energy_level` は `wb97x-d3/def2-tzvpd`、δG_eff 41.926(上の composite と同じ)。同じ pipeline ファイルを configs/ の外に置き、s4_ch3o_doublet の複製に追記した(method は作業ディレクトリの configs/methods から): 1:18 で完走し(レビュー P1c の rc=2 が解消)、⟨S²⟩ 0.754〜0.758 で `spin_contaminated` なし、δG_eff 31.18。HCN の TS の CCSD(T) は最初の実行で別の run と重なって rc 1 で落ち、`--retry-failed nonzero_exit` で 18 s で通った。

**S-A の SA-B・SA-C の確認(U8-P4・U8-P7・U0-P5)**

- U8-P4(SP の対象、QM なしの数え直し): M5 の run で、順位の付く run の SP 対象数は変わらない(hcn 3、hono 3、nh3 2、malonaldehyde 2、ch3o 3、sn2 2)。`unresolved` の S17 は 1 → 0、`same_basin` の水は 1 → 0。
- U8-P7(`/home/user/hfauto_r6/SA/sac_probe`、別の NWChem 4 rank と同時に走った時間を含む): TCE の ROHF-CCSD(T)/def2-TZVPD の出力は SCF の `open shells = 1` と `CCSD(T) total energy / hartree`。OH• は 1.8 s。CH3O• は `2eorb` なしだと v2 積分の GA 確保(378 MB/rank)に失敗して abort、`2eorb` だけで 407 s、`2eorb 2emet 13` で 63 s(E −114.874239093 Eh、3 通りとも 1e-10 Eh で一致)。golden G29 は NWChemEngine で作った OH•(E −75.640597872 Eh、Level の多重度 2、⟨S²⟩ なし)。
- U0-P5: NWChem 7.2.3 の DFT の反復上限の既定は 50(レビューの 30 ではない)。P3c の ZTS bead 4(H3、二重項)は原子の初期推定からなら既定で 13 反復で収束する(E −1.548056711 Eh、⟨S²⟩ 0.987)が、失敗した attempt の movecs 4 通りから始めると、既定・`iterations 100`・`damp 40 ncydp 30 lshift 0.5`・`convergence rabuck 20` はすべて未収束で、`cgmin` だけが 5〜8 反復で収束した(勾配 5e-4 で止まるので E は 1.5〜3.6e-7 Eh 高い)。I⁻+CH3I(ECP あり)はどの方法でも収束した(E の差 2e-8 Eh 以内)。救済は cgmin の 1 本。cgmin は ⟨S²⟩ を出さないので、⟨S²⟩ のない開殻 DFT の Evidence は `incomplete_output`(`s2_not_reported`)にした。閉殻の G26(FHF⁻)と G27(HI)の SP を新しいデッキで再計算すると E の差は 2e-9 Eh 以内で、G27 の SCF は 7 反復のまま。smoke(`-m real tests/smoke/test_real_nwchem.py`)は OH• の ROHF-CCSD(T) を加えた 7 件が合格。

**S-A の統合後の実測(作業ツリー、`/home/user/hfauto_r6/SA/int/runs/<run>`、R1〜R4・S18・S19 は新しい run dir)**

U8-P3 のプローブで既定の pipeline は PBE0-D3BJ/def2-SVPD のままなので、R1〜R4 の δG_eff は M5 と同じで、ranking.csv に `energy_level` 列(`pbe0-d3bj/def2-svpd`)が加わった。以後の段はこの表の値と比べる。エネルギー層は method_panel の追記で使う。

| 系 | 結果 | 所要時間 |
|---|---|---|
| smoke | `-m real tests/smoke` 14 件合格 | 1:15 |
| R1 `hcn` | `elementary_step`、δG_eff 42.1783 | 1:21 |
| R2 `hono` | `elementary_step`、δG_eff 11.8158、torsional True | 3:22 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff 3.8395 | 1:16 |
| R4 `water_same_basin` | `same_basin`、順位なし | 0:11 |
| S20 HCN(R1 の複製に CCSD(T) 入りのパネルを追記) | configs/ の外の pipeline(method は configs/methods から)で完走。`energy_level` `wb97x-d3/def2-tzvpd`、δG_eff 41.926。ΔE‡: PBE0/SVPD 46.621、PBE0/TZVPD 46.429、ωB97X-D3 46.369、CCSD(T) 47.807。ranking.csv の panel min / max 46.369 / 47.807 は method_panel.csv と一致 | 1:23 |
| S20 SN2(M5 の sn2_panel の複製) | ΔE‡: PBE0/SVPD 10.446、PBE0/TZVPD 11.146、ωB97X-D3 15.101、CCSD(T) 13.436。δG_eff 15.732(ωB97X 層)、panel min / max 10.446 / 15.101 が一致 | 0:01(ジョブはすべて再利用) |
| S20 CH3O•(M5 の s4 の複製) | ROHF-CCSD(T) の行が埋まった(ΔE‡ 33.372、ΔE_rxn −7.802。以前は `wft_closed_shell_only`)。ωB97X-D3 33.578、PBE0/TZVPD 32.713。δG_eff 31.178、`spin_contaminated` なし、panel min / max 32.713 / 33.578 が一致 | 4:38 |
| S18 `h3_doublet` | `degenerate_rearrangement`(TS −605.8i)、ΔE‡ 4.563、δG_eff 4.091。SCREEN の single から saddle が直接収束し、string も SCF の救済も使わなかった(cgmin は SA-C のプローブで P3c の bead を収束させた) | 0:47 |
| S19 `tma_hf2`(`discover`、新しい run) | explore 25 attempt で生成物 0。宣言反応 `neutral_to_shared_proton` は DFT で shared_proton が neutral に落ちて `same_basin`(順位なし)。SP の対象は 0 点(dft 極小 3、旧 `reaction_stationary_points` は 1、`all_minima` は 3)。dft stage が 1,787 s で大半 | 32:02 |
| S19 複数条件(tma_hf2 に thermo_conditions → report_conditions を追記。T 250 / 298.15 / 400 K × 1 atm / 1 M) | ΔG_assoc(TMA + 2 HF → 錯体)は 1 atm で −14.52 / −11.53 / −5.11、1 M で −17.52 / −15.32 / −10.66 kcal/mol(差 2RT ln(RT/P°) = 3.00 / 3.79 / 5.55)。report は最初の (298.15 K、1 atm) で並べた。SN2(sn2_panel の複製に同じ追記)の δG_eff は 1 atm と 1 M で同じ(250 K 10.898、298.15 K 11.077、400 K 11.435) | 0:02 |

**S-B の SB-A のプローブ(U5-P2 の採否、`/home/user/hfauto_r6/SB/sba/runs/<run>`)**

作業ツリーの複製(2026-09-28 06:30、SB-B・SB-C の途中の変更を含む。`driver.py` の `split` にだけ R・I・P の構造を渡す配線を足した)で、基準の run の複製に `known_endpoints --from paths` を 1 本ずつ流した。基準は S-A の統合 run(R1〜R3・S18)か、各系の最新の M 段の run(S-A は経路を変えていない)。時間は paths stage(基準がジョブをほぼ再利用した run は、最初に全部を計算した run の値)。

| 系 | 基準 | SCREEN(基準 → NEB) | outcome・ΔG‡(kcal/mol) | paths(s) |
|---|---|---|---|---|
| HCN | SA int | GS の TS → NEB の TS(13 サイクル)。上界 33.98 → 34.21 | `elementary_step`、42.1783 で同じ | 59 → 58 |
| NH3 | SA int | TS → TS(5)。4.46 → 4.48 | `degenerate_rearrangement`、3.8395 で同じ | 62 → 62 |
| H + H2(S18) | SA int | HEI → HEI(20、CI の TSOpt なし)。81.09 → 7.56。種の Hessian は xTB(重なり 0.91)→ DFT(xTB に方向の負モードなし) | `degenerate_rearrangement`、4.0909 で同じ | 36 → 38 |
| H2O2(S13) | M3 | IDPP の HEI → NEB の HEI(5、TSOpt なし)。1.28 → 1.30 | `degenerate_rearrangement`、1.1802 で同じ。FIND_PATH なし | 74 → 76 |
| HONO | SA int | TS → TS(12)。14.21 → 14.37 | `elementary_step`、11.8158 で同じ。FIND_PATH なし | 146 → 126 |
| (HF)₂(S17) | M2 / M5 | HEI(同じ側の H 回転、7.8)→ NEB の TS(9)。1.74 | `unresolved_within_budget`(saddle 2 回 maxiter、FIND_PATH)→ **`degenerate_rearrangement`(C2h の入れ替え、TS −227.3i)、ΔE‡ 1.256、ΔG‡ 1.392** | 408 → 73 |
| Cl⁻·CH3Cl(S1) | M1 | TS → TS(11)。12.25 → 12.60 | `degenerate_rearrangement`、ΔE‡ 10.4459 で同じ、ΔG‡ 11.0771 → 11.0773 | 441 → 314 |
| DME(S12) | M2 | TS → TS(19) | `degenerate_rearrangement`、1.9638 で同じ | 446 → 446 |
| マロンアルデヒド(S14) | M5 | TS → TS(7)。3.43 → 3.10 | `degenerate_rearrangement`、δG_eff 0.0 で同じ(ΔG‡ 0.2499 → 0.2498) | 651 → 640 |
| trans-HONO → HNO2(S15) | M2 | HEI → TS(14)。85.09 → 55.55 | `elementary_step`、同じ直接の TS(−2085.6i → −2088.1i、E_TS は 7e-5 Eh 低い)、ΔG‡ 52.060 → 52.016 | 240 → 114 |
| acac(S16) | M3(60 分で打ち切り) | TS → TS(13)。3.03 | M3 は QRC の途中で打ち切り → **`degenerate_rearrangement`(TS −986.8i、`soft_secondary_mode`)、ΔE‡ 1.980、ΔG‡ −0.93 で δG_eff 0.0**。縮退の判定を実計算で初めて確認 | 全体 69:31(paths 4,171: SP 9 本 325、saddle 229、TS freq 712、QRC 1,402 + 1,339) |

- NEB は全系で収束し(5〜20 サイクル、1 本 2〜4 s、未収束の NEB はなかった)、CI からの TSOpt は H2O2 と H + H2 以外で収束した。どの系も FIND_PATH に入らなかった(基準では S17 だけが入った)。TS は S15 以外で E_TS が 1e-9 Eh 以内で同じ。
- SCREEN の上界(内部の DFT 最大 − 高い方の極小)は、xTB の端点から始まる GS より DFT 極小を固定した NEB の方が TS に近い(S15 85.1 → 55.6、H + H2 81.1 → 7.6。TS の ΔE‡ は 55.3、4.6)。
- 判定: どの系も同じか改善(S17 は未解決から順位の付く縮退転位になり、acac は完走し、S15・S1・HONO は速くなった)、HONO と H2O2 は FIND_PATH を通らないので、U5-P2 の固定端 CI-NEB を採用した(GS と DLC の後退案は入れない)。
- S17 は障壁 1.26 kcal/mol で解像度(1.0)より高いので BARRIERLESS ではない。実計算で順位の付く BARRIERLESS の行はまだない。S15 の cis-HONO を経る 2 段の経路は実計算では通らない(直接の TS が見つかる)。
- 同じ段の SB-A の QM なしの確認: 決定表(行 8 のきっかけは VALIDATE_INTERMEDIATE が消費する、旧行 10 の削除、行 15 の次チャンク)、multi_max の井戸が端点だったときに同じ種を積み直さないこと、saddle の maxiter で最終フレームが種になり 1 回だけやり直すこと、0.5 Å 以内の TS freq の再利用(NWChem の stub と fake)、pysis の入力。

**S-B の統合後の実測(作業ツリー、`/home/user/hfauto_r6/SB/runs/<run>`、R1〜R4・S4・S13・S15・S18 は新しい run dir)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| smoke | `-m real tests/smoke` 14 件合格(pysis_neb の HCN と NWChem の 1 チャンク string を含む) | 1:12 |
| R1 `hcn` | `elementary_step`、δG_eff 42.1783(S-A と同じ) | 1:17 |
| R2 `hono` | `elementary_step`、δG_eff 11.8158(同じ)。SCREEN は NEB の TS、FIND_PATH なし | 3:00 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff 3.8395(同じ) | 1:14 |
| R4 `water_same_basin` | `same_basin`(同じ) | 0:12 |
| S13 `h2o2_gauche` | `degenerate_rearrangement`、ΔG‡ 1.1802(NEB の HEI、DFT Hessian の種)。FIND_PATH なし | 1:35 |
| S18 `h3_doublet` | `degenerate_rearrangement`(TS が得られる)、δG_eff 4.0909 | 0:46 |
| S15 `hono_hno2` | `elementary_step`、ΔG‡ 52.016(直接の TS。cis の井戸は経路に現れず `multi_step` にはならない) | 2:41 |
| S4 `ch3o_doublet`(Cs 入力) | 宣言反応は `elementary_step`、ΔG‡ 30.6117(同じ)。Cs の CH2OH• は DFT で saddle になり、mode-follow の両側が鏡像になって `ts_candidate` → `rxn_mode_follow_…` が `degenerate_rearrangement`(torsional、δG_eff 3.876)として順位 1 に加わり、ch2oh は注記 `endpoint_was_saddle` 付きで basin に入った(U3-P3・U3-P4 で予期した増加) | 13:41(M2 8:53。増えたのは mode-follow の case) |
| S19 `tma_hf2`(S-A の run の複製を `--from dft`) | dft の DFT ジョブ 10 本(すべて JobStore の再利用、S-A も 10 本)で U3-P6 の窓で増えない。`same_basin` は同じ | 0:02 |
| S16・S17・S1・S12・S14(SB-A プローブの run の複製を `--from paths`) | paths のジョブ 16 本すべて再利用で、outcome と ΔG‡ はプローブと同じ(acac `degenerate_rearrangement` δG_eff 0.0。U6-P7 の結合グラフ確認を通る。種の DFT Hessian 0 本。M3 は 60 分で QRC の途中だったが、プローブでは全体 69:31 で完走し、QRC が paths の 2,741 / 4,171 s) | 各 0:02〜0:10 |
| S11 DME の C2v 種(`/home/user/hfauto_r6/SB/ho_harness.py`、REFINE_SADDLE と VALIDATE_TS の手順を再現) | 1 回目: xTB Hessian(二面角方向との重なり 0.30)で 125 s、C2v の 2 次 saddle(−232.1i / −128.6i、ΔE 4.03)→ `higher_order`。反応モードでない −128.6i に沿って `modes.amplitude` 0.23 Å 押し、TS freq を Hessian に 2 回目 → maxiter で 190 s で停止(以前の drv.hess 継続は 495 s)。最終フレームから xTB Hessian(重なり 0.82)で再開 → 17 歩 160 s で一次鞍点 −226.9i、ΔE‡ 2.33(P5b と同じ) | 8:02 |

- S11 の再開は 3 回目の saddle で、`max_saddle_attempts` 2(再開も同じ予算)の実際の case では FIND_PATH に進む。停滞と再開の動作はレビューの予測どおり。

**S-C の SC-B のプローブ(U2-P2・U2-P6・U1-P4、`/home/user/hfauto_r6/SC/scb/{runs,logs}`、集計は `probe_sum.py`)**

structures → conformers → screen(GFN2)。組成はすべて全原子を `--notopo` に渡す。A は `--noopt`、B は初期最適化あり + `--noreftopo`。CREST の列は rc / 配座数 / 秒。

| 組成 | A | B | screen の状態ラベル(A と B で同じ) |
|---|---|---|---|
| TMA·(HF)₂ | 0 / 3 / 11.3 | 0 / 3 / 10.8 | C3H10FN+FH(最低 E の差 2e-5 kcal/mol) |
| F⁻·(HF)₂ | 0 / 1 / 1.4 | 0 / 1 / 2.8(継続した実行に 'Change in topology detected' が出て、パーサは停止として再実行した) | F+FH+FH |
| Cl⁻·CH3Cl | 0 / 3 / 3.9 | 0 / 3 / 3.8 | CH3Cl+Cl |
| SO2·NMe3 | 0 / 1 / 11.9 | 0 / 2 / 10.0 | C3H9N+O2S(付加体は CREGEN を通る) |
| NH3·HF⁺(m=2)、OH·H2O(m=2)、Na⁺·H2O、CH3·O2(S5、m=2) | rc −11(trial MTD の segfault。`-T 1`・quick なしでも同じ) | rc 1(初期最適化の失敗) | placement の seed を screen が緩和: F+H4N 0.00 / FH+H3N 5.86、H2O+HO、H2NaO、CH3O2 |

- U2-P2: A を残す(B も開殻と Na⁺·H2O で失敗し、F⁻·(HF)₂ で停止の誤検出を起こす)。採ったのは全原子の `--notopo` だけで、CREST が失敗した組成は seed を出す。
- U2-P6(円錐 vs 剛体、A): TMA·(HF)₂ と F⁻·(HF)₂ は同じラベルで ΔE 4e-5 / 2e-5 kcal/mol(剛体は TMA·(HF)₂ で 1 配座 8.9 s、円錐は 3 配座 11.3 s)。CREST が失敗する NH3·HF⁺ と OH·H2O では剛体 seed の xTB 最適化 6 本がすべて失敗し、錯体の状態が消える。円錐を残す。
- U1-P4 のラベル差分(hfauto_v2〜v5・w7 の 204 参照、83 構造): 変わったのは 9 構造で、すべて GFN2 の CREST 配座の H 結合 amine·HF 錯体(v2 nh3_hf3 c02〜c05 の N···H 1.445〜1.448 Å、v2/v3/w7 tma_hf2 c01/c02 の 1.426〜1.464 Å。閾値 1.42 Å)が 1 断片から分かれた(C3H10FN+FH → C3H9N+FH+FH)。tma_hf2 c00(1.413 Å)、DFT・screen・入力の構造(TMA·(HF)₂ の端点 N–H 1.27 / 1.32 Å を含む)は変わらない。閾値の ±0.05 Å に 20 対。K⁺ の構造はこれらの run にない。S7 の N···I 2.767 Å と N···Cl 2.487 Å、SO2·NMe3 の N···S 2.167 Å(閾値 2.160 Å)は 2 断片になる。Na⁺·H2O は 1 断片のまま(旧規則の Na···H の結合が消えハッシュだけ変わる)。

**S-C の SC-A の確認(U4-P4・P6・P7・P9・P11、U9-P6、`/home/user/hfauto_r6/SC/sca`)**

- S6: CH4···OH の接触構造(O···H 2.30 Å)から引き抜きの T1 を実 ReaDuct で実行した。worker の source 再最適化(GFN2)で構造がすでに CH3···H2O(C···H 1.57 Å、O–H 1.02 Å)に落ちるので、NT2 は `ts_not_converged`(Bofill: No more negative eigenvalues)、AFIR(γ 125・300、restricted_open_shell、SCC 収束)は `same_as_source`。原因は GFN2 で引き抜きが障壁なしであることで、反復上限・変位・開殻の AFIR ではない。引き抜きは screen の緩和の発見(`collapsed_to:CH3+H2O`)として残る。
- U4-P4: M4 の amine_pilot2 の `irc_end_not_minimum` 7 件を同じ job.json で再実行すると 7 件すべて終わった(生成物 5、`irc_not_connected_to_source` 2、各 7〜17 s、変位を使ったのは 1 件)。障壁は 102〜175 kcal/mol で窓の外。
- U4-P11(M4 の manifest で数え直し): attempt 数は amine_pilot2 29 → 20(AFIR 12 → 3)、S1 16 → 13、S6 13 → 10、S7 5 → 5。失う生成物の basin はない(S1 の AFIR 生成物 2 件は NT2 と同じ恒等 SN2)。
- U4-P7: S1 の恒等 SN2 の生成物 6 件は species 1 つになり、source の basin(c00・c01)には合流しない。

**S-C の統合後の実測(作業ツリー、`/home/user/hfauto_r6/SC/runs/<run>`、すべて新しい run dir。S 系は `discover --to explore`)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| smoke | `-m real tests/smoke` 14 件合格 | 1:14 |
| R1 `hcn` | `elementary_step`、δG_eff 42.1783(S-B と同じ) | 1:20 |
| R2 `hono` | `elementary_step`、δG_eff 11.8158(同じ) | 3:06 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff 3.8395(同じ) | 1:15 |
| R4 `water_same_basin` | `same_basin`(同じ) | 0:12 |
| S1 `sn2_cl` | 13 attempt(AFIR 3)。恒等 SN2 の生成物 4 件(2 つの錯体から NT2 の T1 と T4 切断)が species 1 つ(M4 は 16 attempt、生成物 6 件がそれぞれ species) | 0:14 |
| S6 `oh_ch4` | 引き抜きの生成物はない(上の SC-A の確認のとおり GFN2 で障壁なし、緩和の発見 2 件)。出発点の CH3···H2O から 11 attempt(M4 は 13)で、CH3OH + H(ΔE‡ 31.6 / 37.8、ΔE_rxn 10.9 / 6.2 kcal/mol)が新しい窓(ΔE‡ ≤ 50、ΔE_rxn ≤ 40)に入って 2 生成物として残る(遊離 H の位置が違い RMSD で別の basin)。rc 1 は CREST の開殻の失敗と placement seed の xTB 失敗(既知) | 0:06 |
| S7 `nh3_icl` | 錯体 5 構造はすべて 2 断片の 1 状態(ClI+H3N)、polar_h の drive はない。出発点(最低 2 つ、N···I 2.77 Å の c00 と c01)の drive は N–I 形成・I–Cl 切断の T1 だけで、NT2 は `no_nt2_maximum`、AFIR は `same_as_source`。M4 の生成物 3 件のうち 2 件(AFIR の ClH3IN)はハロゲン結合の錯体そのもので、新しい結合規則では生成物でない。NT2 の NH2I + HCl(124.4 kJ/mol)は緩い c02 から得たもので、M4 では c00 が別の状態(1 断片)だったので c02 が出発点に入ったが、c00 が同じ状態に入った今は最低 2 つに入らない | 0:13 |
| S8 `so2_nme3` | 付加体(N···S 2.17 Å)は CREGEN を通り screen 極小として残るが、ラベルは 2 断片(C3H9N+O2S)。N–S 形成の T3・T1 は `no_nt2_maximum` → AFIR `same_as_source`。23 attempt で生成物 0 | 1:04 |
| S10 `amine_pilot2`(`discover`) | explore 23 attempt(AFIR 2)、37 s(M4: 36 attempt・AFIR 15・120 s)。`irc_end_not_minimum` 7 → 0。NH3·HF の二重 H 交換は 2 attempt(リレーと T1)が同じ species 1 つに入り、`degenerate_rearrangement`、ΔG‡ 32.29(M4 と同じ)。DFT ジョブ 18(M4 も 18)、dft 1,134 s(M4 1,168)。CREST は同時に 1 本(ReaDuct は最大 4 本) | 21:36(M4 23:35) |
| S19 `tma_hf2`(`discover`、cores 4) | explore 15 attempt(すべて NT2。S-A は 25 で AFIR 10)、32 s。生成物 0、`same_basin` は同じ。dft の DFT ジョブ 9 + 再利用 1(S-A も同じ)、1,778 s(S-A 1,787)。CREST 2 本は同時に 1 本、ReaDuct は最大 4 本 | 0:53(`--to explore`)+ 29:39(残り) |
| 同上の直列(run の複製から ReaDuct のジョブを消し、cores 1 の site で `--from explore`) | explore 124 s(並列はその 1/3.9)。explore・screen・conformers の artifact は run dir を除いて並列と一致し、ReaDuct の JobStore キー 15 個も一致(直列のキー 32 個はすべて並列の run にある) | 2:06 |

- 並列の確認に使った同時実行数は、各ジョブの attempt ディレクトリの最初のファイル(`input.xyz` を除く)から `result.json` までの区間で数えた。job.json はコアの予約の前に書かれるので、その時刻で数えると CREST が 2 本重なって見える。
- 独立に走らせた 2 本の run(並列・直列)は CREST の MTD の揺らぎで配座数が違い(7 / 5 artifact)、artifact の比較には使えない。上の直列は同じ run の複製で比べた。
- 判定: R1〜R4 は S-B と同じ値(結合規則と explore の変更は known_endpoints の順位に影響しない)。AFIR は `no_nt2_maximum` の後だけになり、attempt は S10 で 36 → 23(半分には届かない。減ったのは AFIR 15 → 2)、tma_hf2 で 25 → 15。`irc_end_not_minimum` は S10 で 7 → 0。生成物の同定で S1 は species 6 → 1、S10 は 2 attempt → species 1 で、同定で失った basin はない。S7 の NH2I + HCl は、状態が 1 つにまとまって出発点(状態ごとに最低 2 つ)から外れたために失った(表の S7)。S6 の引き抜きの生成物は得られない(GFN2 で障壁なし。新しい機構は入れない)。

**S-D の確認(U7-P5・U9-P7、実行はすべてリポジトリの外)**

- U7-P5: GoodVibes 4.3.0 を同一プロセスで呼ぶ species_thermo で、既存 run 7 本(SC の HCN・HONO・NH3・水、SB の CH3O・TMA·(HF)₂、M5 の SN2 パネル)の SpeciesThermo 22 件を再計算し、G・H・ZPE が worker 経路の値と 1.1e-13 Eh 以内で一致(順位は不変)。G06 CSV とは H・ZPE が 1e-10 Eh 以内、G は 1.3〜2.0e-8 Eh(CLI が NWChem の 7 桁の回転定数を読んだ差)。1 回 1〜2 ms(17 原子)。Windows の .venv-win-qa には pymsym がビルドできず GoodVibes を入れられないので、golden は WSL で走る。
- U9-P7(`/home/user/hfauto_r6/SD/sdc_signal`): `mpirun -np 2 nwchem`(ベンゼン PBE0/def2-SVPD freq)を run_command で走らせる harness に `timeout 20`。signal handler ありで rc 124、残るプロセス 0。なしでは prterun と 2 rank が残る(U9-I8 の再現)。

**S-D の統合後の実測(作業ツリー、`/home/user/hfauto_r6/SD/runs/<run>`、すべて新しい run dir。MethodSpec と Level から solvation が消えたのでジョブの鍵はすべて変わった)**

| 系 | 結果 | 所要時間 |
|---|---|---|
| smoke | `-m real tests/smoke` 13 件合格(MP2 の smoke は閉殻 CCSD(T) に置換) | 1:07 |
| R1 `hcn` | `elementary_step`、δG_eff 42.1783(同じ) | 1:17 |
| R2 `hono` | `elementary_step`、δG_eff 11.8158(同じ) | 3:02 |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff 3.8395(同じ) | 1:15 |
| R4 `water_same_basin` | `same_basin`(同じ) | 0:11 |
| S15 `hono_hno2` | `elementary_step`、δG_eff 52.0163(直接の TS。分割は起きず、子の予算は fake の stage テストで確認) | 2:43 |
| acac(`timeout 60`) | rc 124。終了後 `pgrep -f nwchem`・`pgrep -f prterun` は空、run_state の dft は `failed` | 1:01 |
| b2plyp | NWChemEngine.energy が `input_invalid`(`unsupported_method:b2plyp_def2-svpd`)、ジョブ 0 本 | — |

- 判定: R1〜R4 の outcome と表示桁の値は S-C と同じ。全桁では δG_eff が 2e-7〜1.4e-6 kcal/mol 違う。同じデッキ(scratch_dir だけ違う)の HCN の最初の opt でも E が 7.7e-10 Eh 違うので、NWChem の実行ごとの揺らぎ(ジョブの鍵が変わって再計算した)で、thermo の変更によるものではない(保存済みの Evidence からの再計算は 1.1e-13 Eh で一致)。種ごとの差は G で 3e-9 Eh 以下。
- 外部の timeout の後に NWChem が残らないので、runcase.sh の pgrep / pkill の回避策は要らない。

**FINAL の実測(レビュー §6.1 の全体を 1 回、作業ツリー、`/home/user/hfauto_r6/FINAL/runs/<run>`)**

FINAL で Evidence の `trajectory_energies_hartree` と ConnectionClaim の `method`・`notes` を消した(`extra="forbid"`)ので、以前の run の Evidence と manifest は読めない。ジョブの鍵は変わらないが、すべて新しい run dir で流した。S19 と S20 は同じ表の run の複製に pipeline を追記し、S11 は S-B の harness を今の API に合わせた複製(`/home/user/hfauto_r6/FINAL/ho_harness.py`)で流した。ジョブは misses + hits(括弧内は hits)。

| 系 | 結果 | 所要時間 | ジョブ |
|---|---|---|---|
| smoke | `-m real tests/smoke` 13 件合格 | 1:04 | — |
| R1 `hcn` | `elementary_step`、δG_eff 42.1783(S-D と同じ) | 1:16 | 20(1) |
| R2 `hono` | `elementary_step`、δG_eff 11.8158(同じ) | 3:03 | 20(1) |
| R3 `nh3_inversion` | `degenerate_rearrangement`、δG_eff 3.8395(同じ) | 1:15 | 19(1) |
| R4 `water_same_basin` | `same_basin`(同じ) | 0:11 | 3(0) |
| S1 `sn2_cl` | `degenerate_rearrangement`、ΔE‡ 10.446、ΔG‡ 11.0773、`dG_act_vs_separated` 6.00 | 7:20 | 25(1) |
| S2 `sn2_i` | 全ジョブ alpha 30(`I library def2-ecp`)、`degenerate_rearrangement`、ΔG‡ 6.88 | 13:11 | 21(1) |
| S3 `h2te` | structures → dft → thermo を通過 | 0:17 | 2(0) |
| S4 `ch3o_doublet` | 宣言反応 `elementary_step` ΔG‡ 30.6117、mode-follow の縮退 case(torsional)3.877 が順位 1(S-B と同じ) | 13:27 | 34(2) |
| S5 `ch3_o2` | 多重度を外した複製は conformers で `INPUT_INVALID`(`declare_multiplicity: candidates (2, 4)`)。宣言 2 では CREST(`--uhf 1`)が rc −11(既知)で placement の 6 seed を screen。`--to dft` で CH3 は m 2・⟨S²⟩ 0.7544、O2 は m 3・⟨S²⟩ 2.0096 | 0:01 / 0:04 / 0:43 | — |
| S6 `oh_ch4` | explore 11 attempt(AFIR 3)、生成物は CH3OH + H の 2 件(ΔE‡ 31.6 / 37.8、S-C と同じ。引き抜きは GFN2 で障壁なし) | 0:06 | 25(2) |
| S7 `nh3_icl` | 4 attempt、生成物 0、polar_h の drive なし(S-C と同じ) | 0:13 | 22(2) |
| S8 `so2_nme3` | 付加体は screen 極小として残る。24 attempt(AFIR 5)、生成物 0 | 1:02 | 35(1) |
| S9 `fecl3_ch4` | 多重度を外した複製は `species_fecl3` が `INPUT_INVALID`(`declare_multiplicity: d-block element(s) Fe`)。宣言 6 では `--to conformers` が例外なく終わり、組成の CREST 配座 6 | 0:01 / 0:31 | 2(0) |
| S10 `amine_pilot2`(`discover`) | explore 24 attempt(AFIR 3)、生成物 1、発見の case が `degenerate_rearrangement` ΔG‡ 32.2887(同じ)。dft 1,195 s | 22:41 | 65(5) |
| S11 DME の C2v 種(harness) | 1 回目: xTB Hessian(重なり 0.30)で 135 s、2 次の saddle(−232.1i / −128.6i、ΔE 4.03)→ `higher_order`。−128.6i に 0.23 Å 押し、TS freq を Hessian に 2 回目が 281 s で一次鞍点 −226.9i、ΔE‡ 2.33。停滞は再現せず、`max_saddle_attempts` 2 の中で終わる(S-B は 2 回目が maxiter で 3 回目の再開が要った) | 7:02 | — |
| S12 `dme_rough` | `degenerate_rearrangement`、ΔG‡ 1.9638(同じ) | 9:49 | 19(1) |
| S13 `h2o2_gauche` | `degenerate_rearrangement`、ΔG‡ 1.1802(同じ) | 1:37 | 19(0) |
| S14 `malonaldehyde` | `degenerate_rearrangement`、注記 `submerged_barrier`、δG_eff 0.0(ΔG‡ 0.2498) | 16:12 | 19(1) |
| S15 `hono_hno2` | `elementary_step`、δG_eff 52.0163(直接の TS、分割なし) | 2:45 | 20(1) |
| S16 `acac` | `degenerate_rearrangement`(xTB Hessian、重なり 1.00、種の DFT Hessian 0 本)、ΔE‡ 1.980、ΔG‡ −0.93、注記 `submerged_barrier`、δG_eff 0.0(S-B と同じ)。dft 1,326 s、paths 4,138 s | 1:31:05 | 19(1) |
| S17 `hf_dimer_swap` | `degenerate_rearrangement`(C2h の入れ替え)、ΔG‡ 1.3917(同じ) | 1:34 | 21(1) |
| S18 `h3_doublet` | `degenerate_rearrangement`、δG_eff 4.0909(同じ) | 0:47 | 21(0) |
| S19 `tma_hf2`(`discover`) | explore 15 attempt(すべて NT2)、生成物 0、宣言反応は `same_basin`。dft 1,733 s | 29:46 | 41(2) |
| S19 複数条件(上の複製に thermo → report を追記) | ΔG_assoc は 1 atm で −14.52 / −11.53 / −5.11、1 M で −17.52 / −15.32 / −10.66 kcal/mol(250 / 298.15 / 400 K、S-A と同じ) | 0:03 | 0 |
| S20 HCN(R1 の複製に CCSD(T) 入りのパネル) | `energy_level` `wb97x-d3/def2-tzvpd`、δG_eff 41.9261。ΔE‡: PBE0/SVPD 46.621、PBE0/TZVPD 46.429、ωB97X-D3 46.369、CCSD(T) 47.807。panel min / max は method_panel.csv と一致 | 1:23 | 9 |
| S20 SN2(S1 の複製) | δG_eff 15.7324、`dG_act_vs_separated` 11.43(ωB97X 層、Cl⁻ と CH3Cl の SP を含む)。ΔE‡ 10.446 / 11.146 / 15.101 / 13.436 | 10:33 | 12 |
| S20 CH3O•(S4 の複製) | ROHF-CCSD(T) の ΔE‡ 33.372、δG_eff 31.1776。mode-follow の case は 4.212 | 6:08 | 12 |

- 判定: R1〜R4 と S 系の outcome・数値は、比べられるものはすべて以前の段と同じ(FINAL は定義の統合と読まれない欄の削除で、化学は変えていない)。orphan の NWChem は全 run で 0。
- 全体の壁時計時間は 4:06(smoke から S16 まで直列。S16 が 1:31、S19 が 0:30、S10 が 0:23)。

到達表(`log.jsonl` の decide の reason と job.json から数えた。数は通った run の数):

| 対象 | 実計算で通ったもの | 実計算で通らないもの(QM なしのテストで確認) |
|---|---|---|
| 決定表の行 | 2 same_basin(3)、4 connection(17)、9 ts_validated・10 saddle_converged・12 screen・14 seed(各 17) | 1・3・5〜8・11・13・15〜17(17 行すべての reason を `tests/unit/drivers/test_reaction_case_decide.py` が確かめる) |
| saddle の種と Hessian | 種 `screen_ts` 14・`screen_hei` 2・`discovery_ts` 3、S11 の `higher_order_retry`。Hessian は xTB 15・DFT 2・TS freq(S11) | `path_hei`、`saddle_restart` |
| エンジン | CREST・xTB opt/freq・ReaDuct の NT2/AFIR・pysis NEB・NWChem の opt/freq/SP(PBE0 SVPD/TZVPD、ωB97X-D3、閉殻 CCSD(T)、ROHF-CCSD(T))・nwchem_saddle。GoodVibes は同一プロセス(ジョブなし) | nwchem_string(FIND_PATH)は smoke の 1 チャンクだけ |

---

以下は v3・v4・v5(第2〜第4ラウンド改良後の再検証、2026-09-26)の記録。第3ラウンドの改良(M9・C38・S26、コード `00680e0`)後の v4 の再検証は §11、第4ラウンドの小さな改良(コード `13a4340`)後の退行確認 v5 は §12 にまとめた。§1〜§10 は v3 の記録で、v4 で数値が変わった箇所には §11 への参照を付けた。

- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0)。
- 手法: DFT は `pbe0-d3bj_def2-svpd`(grid fine、SCF 1e-7)、低レベルは GFN2-xTB。設定は `configs/` のまま使い、amine パイロット・(HF)₃ の配座・FIND_PATH・A/B だけ一時 system / pipeline を使った(v2 の `/home/user/hfauto_v2/configs/` を `/home/user/hfauto_v3/configs/` に複製。内容は同じ)。
- コード: `dc82145`(第2ラウンドの改良 R2-W1 `15b0beb`、R2-W2 `dc82145`)。v3 の再検証でもコードの修正は要らなかった。
- run の置き場所: WSL の ext4 上の `/home/user/hfauto_v3/<run>`(すべて新しい run dir)。v2(コード `4e972c5`)の run は `/home/user/hfauto_v2/<run>`、改良前の W7(コード `6b87814`)の run は `/home/user/hfauto_w7/<run>` に残した。リポジトリの `runs/`(過去の run)には書き込まない。
- 実エンジンの smoke(`HFAUTO_REAL=1 HFAUTO_SITE=configs/sites/wsl_local.yaml pytest -m real tests/smoke`): 12 passed(71 秒。v2 は 72 秒)。
- 所要時間は `/usr/bin/time` の壁時計時間。ジョブ数は `hfauto status` の misses + hits(括弧内は再利用)。opt・saddle のステップ数は Evidence の `trajectory_energies_hartree` の長さ − 1。ΔG は 298.15 K、1 atm。

## 1. ベンチマーク(v3)

| # | 系 | pipeline・手法 | v3 の実測値 | 判定 | 所要時間 v3(v2) | run |
|---:|---|---|---|---|---|---|
| 1 | HCN→HNC | known_endpoints | TS −1128.5i、ΔE‡ 46.62、ΔG‡ 42.18、ΔG_rxn 12.68、ΔE_rxn 13.18 kcal/mol、`elementary_step` | 合格 | 1 分 24 秒(1 分 41 秒) | `hcn_known_endpoints` |
| 2 | HONO trans→cis | known_endpoints | TS −681.0i、ΔE‡ 13.67、ΔG‡ 12.23、ΔG_rxn 0.08 kcal/mol、`elementary_step`。SCREEN の GS は v2 と同じく失敗し、FIND_PATH(DFT string)の経路で到達した(§3) | 合格 | 5 分 23 秒(9 分 32 秒) | `hono_known_endpoints` |
| 3 | NH3 反転 | known_endpoints | TS −757.7i、ΔE‡ 4.26、ΔG‡ 3.84 kcal/mol、`degenerate_rearrangement` | 合格(※1) | 1 分 25 秒(1 分 55 秒) | `nh3_inversion_known_endpoints` |
| 4 | 水(同一 basin) | known_endpoints | `same_basin`。reaction-paths のジョブ 0 件。ranking.csv の blocker は `outcome:same_basin` だけ(S22) | 合格 | 12 秒(11 秒) | `water_same_basin_known_endpoints` |
| 5 | TMA·(HF)₂ | discover | `same_basin`(paths 0 秒・0 件)。explore の生成物 0 件。DFT の極小 3 つ(−374.7574786 Eh ほか)は v2 と同じ。ΔG_assoc −11.53 kcal/mol(SVPD、BSSE 未補正。TZVPD//SVPD では −9.04、§3 のプローブ) | 合格 | 30 分 58 秒(31 分 15 秒)。dft stage 1,754 秒(1,771 秒) | `tma_hf2_discover` |
| 6 | amine パイロット(NH3·HF、TMA·HF) | discover `--to explore` | v3 では再実行していない(第2ラウンドの改良が通らない範囲)。v2: 40 attempt、TS 22 件はすべて虚振動 1 本、生成物 0 件 | — | v2: 1 分 54 秒 | `/home/user/hfauto_v2/amine_pilot2_discover` |
| 7 | xTB 初期 Hessian の A/B | minima(dft)、初期 Hessian は GFN2 | A の opt 27 点・854 s(v2 は 27 点・899 s)。最終エネルギーは v2 と 3e-8 Eh 以内で一致(§8) | 記録済み | A: 39 分 51 秒(40 分 34 秒) | `ab_init_hessian` |
| 8 | FIND_PATH | known_endpoints から screen を外した一時 pipeline、PBE0/STO-3G | HCN→HNC。IDPP から string 3 チャンクで bead エネルギーが落ち着き `single_max`。各チャンクの端点は DFT 極小のまま(M8)。HEI から鞍点 −1223.9i、`elementary_step` | 合格 | 1 分 17 秒(1 分 03 秒) | `find_path_sto3g` |
| 9 | (HF)₃ の配座 | discover `--to conformers` | v3 では再実行していない(CREST に関わる改良はない)。v2: NH3·(HF)₃ 7 配座、TMA·(HF)₃ 9 配座(最安はイオン対) | — | v2: 34 秒 | `/home/user/hfauto_v2/amine_hf3_conformers` |

- ※1: `symm` の寄与 +0.411 kcal/mol は外部回転の対称数の比(NH3 の σ 3、平面の TS の σ 6。RT ln 2)による。ΔG‡ は経路縮重度 L を含まない(design.md §7)ので、L = 1 と数えた有効障壁は `symm` を除いた 3.43 kcal/mol(W7 は 3.44)に当たる。
- 検証の基準値(W7 の D3BJ の値)と、旧 run(D3 zero、G06・G07)との差は W7 で確認済みである。

## 2. 第2ラウンドの比較(v2 → v3)

| 系 | 所要時間 | ジョブ数(再利用) | 主な差と原因 |
|---|---|---|---|
| HCN | 1:41 → 1:24 | 29(2)→ 29(2) | QRC 24 / 20 → 12 / 15 ステップ、41.7 → 25.5 秒(S21)。鞍点 5 → 4 ステップ(`moddir 1`、S20) |
| HONO | 9:32 → 5:23(−44%) | 17(0)→ 16(0) | paths stage 516 → 268 秒。IDPP がねじれ経路になり(S18)、string 1 チャンク 200 → 163 秒。HEI の xTB Hessian が採用され、DFT Hessian(15 秒)と、AUTOZ 失敗からの再試行を含む鞍点(106 秒・36 ステップ)が、鞍点 3 ステップ・6.7 秒になった。QRC 34 / 36 → 15 / 15 ステップ、162 → 69 秒(S21) |
| NH3 | 1:55 → 1:25 | 27(1)→ 27(1) | QRC 45 / 45 → 25 / 25 ステップ、62 → 34 秒(S21)。鞍点 3 → 2 ステップ |
| 水 | 0:11 → 0:12 | 4 → 4 | 計算は同じ。ranking.csv の blocker が `outcome:same_basin;thermo_unavailable;dzpe_out_of_tolerance` → `outcome:same_basin`(S22) |
| TMA·(HF)₂ | 31:15 → 30:58 | 60(2)→ 60(2) | 錯体の opt(初期 Hessian あり、trust 0.3)は 8 → 8 ステップ、400 → 378 秒。blocker は S22 で `outcome:same_basin` だけ。explore の 31 attempt の内訳は v2 と同じ |
| FIND_PATH(STO-3G) | 1:03 → 1:17 | 14(0)→ 15(0) | string 2 → 3 チャンク(+18 秒)。v2 の 2 チャンク目は、ずれた端点(HNC 極小より +8.85 kcal/mol)の上で落ち着いていた(M8、§3)。QRC 22 / 30 → 14 / 15 ステップ |
| A/B の A | 40:34 → 39:51 | 11 → 11 | opt 27 → 27 点、899 → 854 秒(§8) |
| smoke | 72 秒 → 71 秒 | 12 → 12 件 | — |

QRC 3 反応(HCN・HONO・NH3)の合計は 204 → 107 ステップ、266 → 128 秒になった。レビューのプローブの予測(HCN 27、HONO 29、NH3 片側 24)とほぼ一致した。接続の判定は 3 反応とも v2 と同じで、`trajectory_above_ts` は出なかった。

ΔG の変化はノイズの範囲(0.002 kcal/mol 以下)に収まった。

| 量 | v2 | v3 | 差 |
|---|---:|---:|---:|
| HCN ΔG‡ / ΔG_rxn | 42.18 / 12.68 | 42.18 / 12.68 | −0.001 / 0.000 |
| HONO ΔG‡ / ΔG_rxn | 12.22 / 0.08 | 12.23 / 0.08 | +0.002 / 0.000 |
| NH3 ΔG‡ | 3.84 | 3.84 | −0.001 |
| TMA·(HF)₂ ΔG_assoc(SVPD、BSSE 未補正) | −11.53 | −11.53 | 0.00 |

HONO の TS は別の Hessian(v2 は DFT、v3 は xTB)から精密化したので、虚振動は −678.1i → −681.0i(W7 は −681.5i)に変わったが、ΔE‡ は 13.67 で同じだった。

## 3. 第2ラウンドの改良項目ごとの確認

- **M8(FIND_PATH の端点を DFT 極小に戻す)**: STO-3G の HCN で、各チャンクの string の両端のエネルギーは DFT 極小(HCN −92.1083949、HNC −92.0683634 Eh)と 1e-8 Eh 以内で一致した。v2 では 2 チャンク目の HNC 側の端点が −92.0542667 Eh(+8.85 kcal/mol)だった。2 チャンク目の初期経路 `string0_c1.xyz` の両端は DFT 極小である。ずれの note は `string0:c1:end_drift:0.074`、`c2:0.001`、`c3:0.002`(Å)で、端点が大きく動いたのは 1 チャンク目だけだった。HONO は 1 チャンクで `string0:c1:end_drift:0.022`。1 チャンク目のずれの原因は IDPP ではなく、初期経路の中間画像を 1 枚ずつ端点へ整列していたことである(画像の間に剛体の跳びが入り、NWChem が凍結端点を回してゆがめる。ゆがんだ端点は極小より STO-3G で +13.6、HONO で +1.09 kcal/mol 高い)。レビュー v3 で分かり(U5-N8)、M9 で逐次整列にした。v4 では 1 チャンク目のずれが STO-3G 0.0008 Å、HONO 0.0016 Å、端点の持ち上がりが +0.002、+0.012 kcal/mol になった(§11)。
- **S18(IDPP の反復上限 5,000)**: HONO の初期経路は、H–O–N=O の二面角が 180° → 160° → 141° → 119° → 96° → 80° → 48° → 24° → 0° と単調にねじれた。v2 は 180° → 178° → 7° と面内に崩れていた。string は 1 チャンクで落ち着き、bead エネルギーの最大は trans から 13.5 kcal/mol(v2 は面内の経路で 34.4)で、HEI の二面角は 89° だった。HEI の xTB Hessian は虚振動 1 本で接線に沿うので採用され、DFT Hessian と AUTOZ の失敗(v2 で約 110 秒)がなくなった。
- **S19(画像の逐次整列)**: STO-3G の HCN で、v2 の鞍点に付いていた `mode_overlap_below_0.3` の note が消えた(レビューのプローブの重なりは 0.118 → 0.651)。
- **S20(moddir の明示)**: 4 つの鞍点(HCN、HONO、NH3、STO-3G)はどれも虚振動 1 本で、`moddir 1` を書いた。AUTOZ の再試行はどの run にもない。負モード 2 本以上の Cartesian・P·H·P の分岐は、実 run では通っていない(プローブ `/home/user/hfauto_v3_probe/w1b_moddir` でだけ確認)。
- **S21(初期 Hessian ありの opt は trust 0.3)**: job.nw は、QRC の両側と錯体の opt が `trust 0.3` と `inhess 2`、それ以外の opt が `trust 0.1` だった。QRC は §2 のとおり約半分になった。xTB の初期 Hessian を渡す錯体の opt(A/B の A と TMA·(HF)₂)は、点数が 4 構造とも v2 と同じで、時間は 5% 減った。v2(trust 0.1)でも、4 錯体の opt の 1〜5 歩目で NWChem は歩幅を制限していた(`Restricting large step` 計 34 回、全体の縮小 計 7 回。v3 は 3 回と 2 回)。歩数が同じなのは、line search が制限された歩みを伸ばして補っていたことと、TMA·(HF)₂ の 4〜8 歩目が勾配の収束後に XMAX / XRMS を満たすための小さな歩みであることによる。5% の短縮は同じ軌跡の run どうしの差(約 2%)より大きいが、n = 2〜3 の比較での観測で、その原因(後半の歩みの SCF が軽くなった)は推定である。trust はジョブキーに入らないので、v2 のキャッシュが残る run dir で再開すると trust 0.1 の結果が再利用される(同じ極小なので害はない)。
- **S22(rankable の理由)**: 順位を付けない `same_basin`(水、TMA·(HF)₂)の blocker は `outcome:same_basin` だけになった。順位の付く 3 反応は v2 と同じ rank 1。
- **S24・S25(同位体と上流の理由)**: 一時 system(`/home/user/hfauto_v3_probe/s24_s25`、scratch を分けた一時 site)で確かめた。SMILES `[2H]O[2H]` の species は `INPUT_INVALID`(`SMILES '[2H]O[2H]' specifies isotopes; hfauto uses natural-abundance masses`)になり、次の stage が `stage 'dft' (minima) has no input of type ['species']; upstream failures: species_d2o: SMILES '[2H]O[2H]' specifies isotopes; ...` で止まった(終了コード 1)。xyz の `D` は `unknown element(s): D` になった。CLI はこのメッセージの前に Rich のトレースバックを表示する(従来どおり。r6 の M0 以後は、この stage は 0 artifact で done になり、report まで書かれる)。メッセージの文言は、のちに実際の質量(主同位体)に合わせて `hfauto uses most-abundant-isotope masses` に変えた(C40)。
- v3 の run で通っていない改良: C2(CCSD の maxiter 50)と S23(composite G の手順)は v3 の run では実行していない。どちらも下の「レビューのプローブ」で確認した。

### レビューのプローブ

検証の run の外で、レビュー(v2・v3)が v2・v3 の run の複製への追記か、エンジンの直接の呼び出しで確かめた。置き場所はリポジトリの外(WSL の `/home/user/hfauto_v3_probe/`、v2 は `/home/user/hfauto_review2_probe/`)なので、要点の数値だけを記す。多くは同時負荷のもとで測ったので、所要時間は参考値である。

- **C2(CCSD の maxiter 50)**: `u0_c2` で NWChem を実エンジンで呼んだ。デッキは `ccsd / freeze atomic / maxiter 50 / end`、出力のエコーは `maxit = 50`。CCSD(T)/def2-TZVPD は HNC で 11/50 反復・17 秒、v3 の HONO の TS で 18/50 反復・447 秒(同時負荷あり)で収束した。HONO の CCSD(T) の ΔE‡ は 11.65 kcal/mol で、SVPD − CCSD(T) = +2.02(design.md §7 の「HONO +2.0」の出典)。
- **S23(composite G)**: `u7_assess/hcn_tz` で、design.md §7 の手順どおりに v3 の HCN の run の複製へ sp_tz(def2-TZVPD、`all_minima`)→ thermo_tz → report_tz を追記した。全体 1 分 26 秒(SP 3 本・83 秒、GoodVibes は 3/3 再利用)、ΔE‡ 46.43、ΔG‡ 41.99 kcal/mol、rank 1。TMA·(HF)₂ では、`tma_tz`(TZVPD//SVPD)の ΔG_assoc が −9.04 kcal/mol(SVPD の −11.53 より +2.49。どちらも BSSE 未補正)、`tma_d`(sp の既定の targets)では単量体が `energy_layer_missing` になり ΔG_assoc は None だった。
- **手法パネル(M5・S8)**: v2 の run の複製に既定の method_panel を追記した(v2 のプローブ)。HCN は 6 SP・36 秒。HONO は panel_report と report の ranking.csv がどちらも rank 1(ΔG‡ 12.22、`method_sign_disagreement` なし)で一致した。続けて methods に ccsd-t_def2-tzvpd を opt-in で足した別の実行(C2 より前の maxiter 20、11〜13 反復で収束)で、HCN の CCSD(T)/def2-TZVPD の ΔE‡ は 47.81 kcal/mol(SVPD − CCSD(T) = −1.2、design.md §7 の出典)。

## 4. 第1ラウンドの比較(W7 → v2)

| 系 | 所要時間 | ジョブ数(再利用) | 主な差と原因 |
|---|---|---|---|
| HCN | 1:40 → 1:41 | 29(2)→ 29(2) | QRC は 19 / 24 → 24 / 20 ステップで、S7 の効果はほぼない |
| HONO | 7:22 → 9:32 | 29(2)→ 17(0) | SCREEN の GS が発散して失敗し、DFT string(200 秒)と鞍点(NWChem の AUTOZ 失敗からの再試行を含め 106 秒)が加わった。QRC は 126(maxiter で継続)/ 17 = 143 ステップ・311 秒 → 34 / 36 = 70 ステップ・162 秒(S7) |
| NH3 | 2:08 → 1:55 | 28(1)→ 27(1) | QRC 50 / 50 → 45 / 45 ステップ(S7)。dft の freq 2 → 1 本(S1: nh3_up が既知の basin に入る) |
| 水 | 0:15 → 0:11 | 5 → 4 | dft の freq 2 → 1 本(S1) |
| TMA·(HF)₂ | 52:33 → 31:15(−41%) | 72(3)→ 60(2) | dft stage 3,050 → 1,771 秒(−42%): opt 4 → 3、freq 4 → 3 本(M3 で偽の AFIR 生成物が消えた)。screen の xTB freq 7 → 3 本(S1)。explore 37 → 31 attempt |
| amine パイロット | 1:56 → 1:54 | 65(4)→ 56(2) | attempt 44 → 40(TMA 単量体の trial が 14 → 11。出発構造が xTB の最適化構造になったため、M3)。TS 19 → 22 件、すべて虚振動 1 本 |
| smoke | 59 秒 → 72 秒 | 10 → 12 件 | 新しい 2 件の分 |

ΔG の変化(S15: 振動数と ZPE のスケールを 0.989 / 0.975 から 1.0 / 1.0 に 1 本化。M6 の寄与は 0 で、反応座標以外の負モードはどの freq にもなかった)は、予測(ΔG‡ で 0.04 以下、ΔG_assoc で 0.07 以下)の範囲に収まった。

| 量 | W7 | v2 | 差(S15) |
|---|---:|---:|---:|
| HCN ΔG‡ / ΔZPE‡ | 42.22 / −3.34 | 42.18 / −3.42 | −0.04 / −0.09 |
| HCN ΔG_rxn | 12.68 | 12.68 | 0.00 |
| HONO ΔG‡ / ΔG_rxn | 12.24 / 0.08 | 12.22 / 0.08 | −0.02 / 0.00 |
| NH3 ΔG‡ | 3.85 | 3.84 | −0.01 |
| TMA·(HF)₂ ΔG_assoc | −11.59 | −11.53 | +0.07 |

HONO の TS は FIND_PATH の HEI から精密化したので虚振動が −681.5i から −678.1i に変わったが、エネルギーは ΔE‡ 13.67 で W7 と一致した。旧い熱化学と NWChem optimize のキャッシュはキーが変わったので再利用されない。

## 5. 第1ラウンドの改良項目ごとの確認(v2)

- **M2・S5(FIND_PATH)**: 初めて実計算で通った。HONO では GS の失敗から行 16 の FIND_PATH に入り、IDPP から 1 チャンク(20 反復)で bead エネルギーが落ち着いた(最後の 3 反復の HEI の変化 0.004 kcal/mol)。NWChem は未収束を報告したが、形は `single_max`(HEI は面内の経路で trans から 34.4 kcal/mol)で、HEI の種から正しい TS に到達した。HCN(STO-3G、#8)でも 2 チャンクで `single_max` になった(このときの 2 チャンク目は端点がずれていた。M8 で解消、§3)。
- **M3・S6(explore)**: TMA·(HF)₂ の neutral からの AFIR は、W7 の「未収束 7・出発に戻る 2・生成物 1」から「未収束 2・出発に戻る 8・生成物 0」になった。W7 で `irc_not_connected_to_source` だった −289i の TS が `same_as_source` と判定されたのは、M3(screen の最適化構造から出発)の効果である。IRC の両端は source の置換コピーで、縮退した HF 交換(ΔE‡ 44.05 kJ/mol)だった。
- **M4・S13・C7(conformers)**: tma_hf2 は W7 と同じ 2 配座(エネルギー差 2e-5 Eh 以内)に、1.48 kcal/mol 上の 1 配座が加わった(3 配座とも screen で同じ basin に入る)。(HF)₃ は、レビューのプローブ(`--noopt` なし)で失敗していた NH3·(HF)₃ を含め 2 組成とも成功し、TMA·(HF)₃ の最安はイオン対だった(#9)。
- **S1(その場の登録)**: 1 basin に freq 1 本になった。screen の xTB freq(TMA·(HF)₂ 7 → 3、amine 7 → 5)、dft の DFT freq(NH3・水 2 → 1)。A/B(§8)では AFIR 生成物の構造が既知の basin に入り、17 原子の freq(約 15 分)を省いた。
- **S7・C19(QRC)**: 3 反応とも `trajectory_above_ts` は出ず、接続の判定は W7 と同じだった。HONO の maxiter での継続はなくなった。
- **S4(GS の回収)**: HONO で GS 自体が例外で終わったため、GS だけの再実行も同じ例外で終わり、`nonzero_exit` として SCREEN が unavailable になった(GS の失敗は握りつぶさない、の設計どおり)。
- 実計算で通っていない改良: M1(陰性結果ゲートの削除)、S3(行 14 の 15 bead 確認の撤廃)は、今回の系に単調な経路や陰性結果で止まる宣言反応がなく、経路に入らなかった(v3 も同じ)。M5・S8(手法パネル)は `method_panel` を実行していない(v3 も同じ。レビューのプローブでだけ確認した、§3)。

## 6. 撤回する値

- HONO trans→cis の ΔG‡ **24.11 kcal/mol**(`hono_isomerization_v3`、CH-38)は撤回する。原因は ZPE の二重計上と basin 判定の回帰である。実測値は **12.22 kcal/mol**(v2、ΔE‡ 13.67。W7 は 12.24、v3 は 12.23)。
- HCN→HNC の ΔG_rxn **18.05 kcal/mol**(内部フォールバックの熱化学、CH-39)は撤回する。正しい値は **12.74 kcal/mol**(D3 zero、G06)で、D3BJ の実測値は 12.68 kcal/mol(W7・v2・v3)。

## 7. 陰性結果の要約(気相の amine·(HF)n)

- **n = 1**(NH3·HF、TMA·HF): HF から N へのプロトン移動は、NT2 が単調上昇(TS なし)、AFIR は出発構造に戻った。NH3·HF の NT2 の TS(−1265i)3 件は、IRC の両端とも出発構造に戻る H 交換だった。v2 でも同じ。
- **n = 2**(TMA·(HF)₂): 入力の neutral・shared_proton と CREST の 3 配座は、xTB で 1 つの basin にまとまった。DFT の極小は neutral の 1 つ(−374.7574786 Eh、W7・v2・v3 で同じ)。宣言反応 neutral→shared_proton は `same_basin` で、TS 探索は起動しない。旧実装では、この系に約 30 h をかけて TS 0 件だった。
- **n = 3**: CREST の配座は得られるようになった(#9。最安は TMA·(HF)₃ でイオン対)。反応探索(障壁なしという設計時点の結論)は再計算していない。
- 同じ組成の旧 run で見つかった TS(G20)は、虚振動の数え方の誤りで 2〜4 本と数えられていた。正しく数えると 1 本だが、いずれもプロトン移動の生成物にはつながらない。

## 8. xTB 初期 Hessian の A/B

2 断片の錯体 4 構造を、同じ開始構造(GFN2 で最適化した構造)から DFT で最適化した。A は discover と同じ設定(`init_hessian: {engine: xtb, method: gfn2}`、NWChem `inhess 2`。v3 は S21 で trust 0.3)、B は初期 Hessian なし(trust 0.1)。一時 pipeline(structures → minima(dft、`include: all`))で実行した。点数は opt の Evidence の `trajectory_energies_hartree` の長さ(開始点を含む)。B は W7 だけで実行した(改良は B の条件に関係しない)。

| 構造(原子数) | A(W7): 点数 / 時間 | A(v2): 点数 / 時間 | A(v3): 点数 / 時間 | B(W7): 点数 / 時間 |
|---|---|---|---|---|
| NH3·HF(6) | 5 / 9 s | 5 / 9 s | 5 / 8 s | 6 / 10 s |
| TMA·HF(15) | 4 / 90 s | 4 / 89 s | 4 / 89 s | 6 / 148 s |
| TMA·(HF)₂(17) | 9 / 405 s | 9 / 397 s | 9 / 378 s | 14 / 625 s |
| TMA·(HF)₂ の AFIR 生成物(17) | 9 / 413 s | 9 / 405 s | 9 / 378 s | 14 / 627 s |
| 計 | 27 点 / 917 s | 27 点 / 899 s | 27 点 / 854 s | 40 点 / 1,410 s |

- W7 では、初期 Hessian で点数は 33%、opt の時間は 35% 減った。xTB の Hessian は 1 本 0.3 秒未満である。
- v2 の A は、opt の点数が W7 と同じで、最終エネルギーも 1e-8 Eh 以内で一致した。AFIR 生成物の構造はその場で TMA·(HF)₂ の basin に入り(members 2)、freq は 4 → 3 本になった。stage 全体の所要時間は 55 分 18 秒 → 40 分 34 秒(−27%)。
- v3 の A(trust 0.3)は、点数が 4 構造とも v2 と同じで、最終エネルギーは 3e-8 Eh 以内で一致した。途中の歩幅は変わった(TMA·(HF)₂ の 1 歩目 −374.75687 → −374.75717 Eh)が、opt の時間は 5% 減っただけで、stage 全体は 40 分 34 秒 → 39 分 51 秒だった。basin 3 つと freq 3 本は v2 と同じ。
- 残りの時間の大半は DFT freq(17 原子で 893 秒、15 原子で 555 秒)である。

## 9. 修正

- **F1**(W7、`hfauto/chemistry/hypotheses.py`): 宣言反応の端点が screen の basin に崩壊した seed(members の 1 つ)の場合、その端点は screen の極小に割り付き、`blocked_upstream`(`endpoints_not_on_one_pes`)になっていた(TMA·(HF)₂ の shared_proton)。端点の basin がその screen basin の DFT 極小であればそちらを使うように変えた。回帰テストは `tests/unit/chemistry/test_hypotheses.py::test_declared_endpoint_collapsed_at_screen_takes_the_dft_basin`。
- v2・v3・v4 の再検証と v5 の退行確認では、コードの修正が必要な不具合は出なかった。

## 10. 未解決事項

1. **HONO の SCREEN(GS)の発散**: pysisyphus の GS(DLC 座標)が約 10 サイクル目から発散し、H–O–N の変角が内部座標から落ちた後、座標数の不一致(66 → 55)で例外終了する。v3・v4 でも同じだった(8 秒、`screen:unavailable:screen_path:nonzero_exit`)。同じ入力では決定的に再現し、W7 の入力(端点の差 1e-6 Å 程度)では収束したので、平面 4 原子での GS の数値的な不安定であり、コードの回帰ではない。v2 で挙げた対策案(小さな分子の GS を Cartesian 座標にする、1 回だけ再試行する)は、レビューのプローブで無効と分かった(Cartesian では面内の経路で HEI 26.56 kcal/mol になり TSOpt も例外、再試行は同じ失敗)ので取り下げる。S18 で FIND_PATH の退避経路がねじれ機構を捉えるようになり、HONO の所要時間は W7(GS 成功、7 分 22 秒)より短い 5 分 23 秒になった。最初の例外は、ジョブの `stderr.txt` に連鎖として残る。
2. S4 の回収は、TS 最適化ではなく GS 自体が例外で終わった場合にも GS を 1 回再実行する(HONO で 4 秒)。害はないので残した。
3. 17 原子・PBE0/def2-SVPD の解析 Hessian は 4 rank で 1 本約 15 分かかり、discover の所要時間の大半(TMA·(HF)₂ の dft stage 1,754 秒のうち 898 秒)を占める。
4. FIND_PATH の HEI は bead エネルギーの放物線補間の最大(offset は ±0.5 でクリップ、座標は隣の bead との線形補間)なので、鞍点をわずかに下回ることがある(HONO で bead の最大 13.47、HEI 13.63、鞍点 13.67 kcal/mol)。種としては問題ない。v3 の 1 チャンク目の端点のずれ(STO-3G で 0.074 Å)は初期経路の整列の仕方によるもので(§3 の M8)、M9 で初期経路を逐次整列にした(v4 の実 run で STO-3G 0.0008 Å・端点の持ち上がり +0.002 kcal/mol、HONO 0.0016 Å・+0.012 kcal/mol。§11)。ずれは引き続き note で監視するだけで、ゲートにはしない。M9 で string の JobStore キー(初期経路の sha256)が変わったが、完了した paths stage は code_version を config_sha に含まないので再開では飛ばされる。`--from paths` などで paths を再実行すると(未完了の場合も)、string と、FIND_PATH の種から始まる U6 のジョブ(種の Hessian、saddle、TS freq、QRC)と、その TS の thermo が 1 回ずつ計算し直される(同じ TS に着く)。
5. S21 の trust はジョブキーに入らない(§3)。
6. 検証の run で通していない経路: 手法パネル(M5・S8・C2。レビューのプローブでだけ確認、§3)、陰性結果ゲートの削除(M1)、15 bead 確認の撤廃(S3)、負モード 2 本以上の鞍点(S20 の Cartesian 分岐)、電荷のある系・開殻系。
7. 設計(WP7.1)では run を `runs/validation_<system>_<date>` に置くことになっているが、過去の run を守るため、WSL の ext4 上に置いた。

## 11. v4(第3ラウンド改良後の再検証、2026-09-26)

- コード: `00680e0`(第3ラウンドの改良 R3。M9・C38・S26)。環境・手法・一時 system / pipeline は v3 と同じ(`/home/user/hfauto_v3/configs/` を `/home/user/hfauto_v4/configs/` に複製)。run は `/home/user/hfauto_v4/<run>`(すべて新しい run dir)、M9 のプローブは `/home/user/hfauto_v4_probe/m9`。
- `hfauto doctor --site wsl_local`: 8 エンジンすべて ok。実エンジンの smoke: 12 passed(71 秒。v3 も 71 秒)。
- 再実行したのは、M9・C38 が通る 5 run(#1〜#4、#8)である。TMA·(HF)₂(#5)、A/B(#7)、amine パイロット(#6)、(HF)₃(#9)は第3ラウンドの改良が通らないので再実行していない。
- v4 の再検証でもコードの修正は要らなかった。

### v3 → v4

| 系 | 所要時間 v3 → v4 | ジョブ数(再利用) | v4 の結果(v3 との差) | ΔG‡ v3 → v4(kcal/mol) |
|---|---|---|---|---|
| HCN→HNC | 1:24 → 1:24 | 29(2)→ 29(2) | `elementary_step`、TS −1128.5i、鞍点 4 ステップ、QRC 12 / 15。すべて v3 と同じ | 42.1783 → 42.1783(0.000) |
| HONO trans→cis | 5:23 → 5:32 | 16(0)→ 16(0) | GS は v3 と同じく失敗し、FIND_PATH の 1 チャンク(20 反復)で `single_max`(bead の最大 13.47 で同じ)。鞍点 3 ステップ。TS −681.0i → −681.5i(W7 は −681.5i)、ΔE‡ 13.67 は同じ。QRC 15 / 15 → 15 / 13 | 12.2256 → 12.2267(+0.001) |
| NH3 反転 | 1:25 → 1:25 | 27(1)→ 27(1) | `degenerate_rearrangement`、TS −757.7i、鞍点 2 ステップ、QRC 25 / 25。すべて v3 と同じ | 3.8395 → 3.8395(0.000) |
| 水(同一 basin) | 0:12 → 0:11 | 4 → 4 | `same_basin`、blocker は `outcome:same_basin` だけ | — |
| FIND_PATH(STO-3G、HCN→HNC) | 1:17 → 1:19 | 15(0)→ 15(0) | string 3 チャンクで `single_max`(3 チャンク目の最大 66.73 → 66.81)。HEI から鞍点 −1223.9i・3 ステップ、`elementary_step`。QRC 14 / 15 → 14 / 14 | 62.1678 → 62.1677(−0.0001) |

- ΔG の差は 0.001 kcal/mol 以下で、順位の付く 4 反応(STO-3G を含む)はすべて rank 1 のままだった。ΔE‡ と ΔE_rxn は v3 と 1e-5 kcal/mol 以内で一致した。
- HONO の +9 秒は string の 1 チャンクの時間(162.9 → 175.4 秒、反復数は同じ 20)による。v3 も 4 反復目で damped Verlet に 1 回切り替えたが、v4 は 5 反復目と 11 反復目の 2 回切り替え、2 回目で bead の評価が 1 回り(7 本、11 反復目は 15.7 s)増えた(評価 149 → 156)。STO-3G でもチャンクごとの評価数は揺れており、M9 の系統的なコストではない。
- QRC 3 反応(HCN・HONO・NH3)の合計は 107 → 105 ステップ。

### M9(FIND_PATH の初期経路の逐次整列)の確認

1 チャンク目の凍結端点のずれ(`job.string_final.xyz` と初期経路の端点の写像 RMSD)と、その端点の実エネルギーの持ち上がり(初期経路の端点 = DFT 極小と、ずれた端点を、string と同じデッキで DFT SP した差)を、v3 と v4 の run に同じスクリプト(`m9_ends.py`)を当てて測った。

| run(ずれる端点) | 1 チャンク目の端点のずれ v3 → v4(Å) | 端点の持ち上がり v3 → v4(kcal/mol) | レビューのプローブの予測 |
|---|---|---|---|
| STO-3G HCN→HNC(HNC 側) | 0.0737 → 0.0008 | +13.57 → +0.002 | 0.004 Å、+0.04 |
| HONO(cis 側) | 0.0220 → 0.0016 | +1.09 → +0.012 | 0.0014 Å、+0.009 |

- v3 の run での値(+13.569、+1.092)はレビュー v3 のプローブ(`u5_assess`)と一致した。始点側(bead 1)のずれと持ち上がりは v3・v4 とも 0.0000 Å・0.000 kcal/mol。
- 初期経路 `string0_initial.xyz` の隣り合う画像では、直接の RMS と写像付き RMSD の差の最大が v3 の 0.27 Å(STO-3G)・0.50 Å(HONO)から、v4 では 1e-13 Å 未満になった(剛体の跳びがなくなった。M9 のテストの判定と同じ量)。
- note `string0:c<k>:end_drift` は、STO-3G で 0.074 / 0.001 / 0.002 → 0.001 / 0.000 / 0.000、HONO で 0.022 → 0.002(Å)。
- チャンク数(STO-3G 3、HONO 1)、形(`single_max`)、HEI の位置(bead 4)、鞍点と接続の判定は v3 と同じだった。STO-3G の 1 チャンク目の bead の最大は 69.74 → 71.75 kcal/mol で、2・3 チャンク目の最大は v3 と 0.2 kcal/mol 以内だった。NWChem が報告する凍結端点のエネルギーは nstep 0 の値のままなので、端点のゆがみは bead のプロファイルには現れない(持ち上がりは SP でだけ分かる)。
- string の JobStore キー(初期経路の sha256)は変わった。v4 はすべて新しい run dir なので、再計算の影響はない。

### C38(code_version)の確認

- v4 の 5 run の `resolved_config.yaml` の code_version は `00680e0fa52527fa1e18ec72a33db492cdb14a0c`(`-dirty` なし)。v3 は 7 run すべて `dc82145…-dirty` だった。`hfauto doctor` の表示も `-dirty` なし。
- Windows 側で `git status` を実行した直後、WSL の `git status --porcelain --untracked-files=no` は 83 行、`-c core.autocrlf=input` 付きは 0 行だった。WSL で後者が index を更新すると、付けない status も 0 行になるが、Windows 側で git を使うたびに 83 行に戻る。C38 により、どちらの順序でも clean な tree は `-dirty` にならない。

### v4 で残る事項

- HONO の GS の発散(§10 の 1)は v4 でも同じだった(`screen:unavailable:screen_path:nonzero_exit`)。
- 検証の run で通していない経路(§10 の 6)は v4 でも変わらない。

## 12. v5(第4ラウンド後の退行確認、2026-09-26)

- コード: `13a4340`(C25・C28・C40・C41・C37・C27′・C34。文書だけの C45・C39′ は `e73adc9`)。環境・手法は v4 と同じ。run は `/home/user/hfauto_v5/<run>`(すべて新しい run dir)、プローブは `/home/user/hfauto_v5_probe/`。
- 対象は known_endpoints の 4 run(HCN、HONO、NH3、水)。第4ラウンドの項目は FIND_PATH(STO-3G)や discover の数値に効かないので、それらは再実行していない。
- `hfauto doctor --site wsl_local`: 8 エンジンすべて ok、版数は `13a4340…`(`-dirty` なし)。実エンジンの smoke: 12 passed(72 秒)。
- コードの修正は要らなかった。

| 系 | 所要時間 v4 → v5 | ジョブ数(再利用) | v5 の結果(v4 との差) | ΔG‡ v4 → v5(kcal/mol) |
|---|---|---|---|---|
| HCN→HNC | 1:24 → 1:22 | 29(2)→ 29(2) | `elementary_step`、TS −1128.5i。QRC 12 / 15 → 15 / 12(両側の向きが入れ替わっただけ) | 42.1783 → 42.1783(−0.00001) |
| HONO trans→cis | 5:32 → 5:32 | 16(0)→ 16(0) | GS は同じく失敗(`screen_path:nonzero_exit`)し、FIND_PATH から鞍点。TS −681.5i → −681.6i。QRC 15 / 13 → 15 / 15 | 12.2267 → 12.2265(−0.0002) |
| NH3 反転 | 1:25 → 1:25 | 27(1)→ 27(1) | `degenerate_rearrangement`、TS −757.7i、QRC 25 / 25。すべて v4 と同じ | 3.8395 → 3.8395(0.000) |
| 水(同一 basin) | 0:11 → 0:11 | 4 → 4 | `same_basin`、blocker は `outcome:same_basin` だけ | — |

- ΔE‡ と ΔE_rxn は v4 と 1e-5 kcal/mol 以内で一致し、順位の付く 3 反応はすべて rank 1 のままだった。4 run の code_version は `13a4340a…`(`-dirty` なし)。
- **C25**: 4 run の ranking.csv の見出しが `rank,reaction_id,outcome,tier,rankable,T_K,standard_state,dG_act_kcal,band_low_kcal,band_high_kcal,dG_rxn_kcal,blockers` になり、値は T_K 298.15、standard_state `1atm`、dG_rxn_kcal は thermo と同じ(HCN 12.683、HONO 0.076、NH3 0.0、水 0.0)。順位の付かない水の行にも T_K と dG_rxn_kcal が入る。
- **C37・C27′**: 73 個の `command_result.json` に `stopped` のキーはない。`resolved_config.yaml` の site から `memory_mb` がなくなった(v4 は 11000)。NWChem の `memory_mb_per_rank` は従来どおり。
- **C41**: v3 の TMA·(HF)₂ で失敗した ReaDuct の attempt 3 件を、同じ `job.json` から v5 の worker で再実行した(`/home/user/hfauto_v5_probe/c41`)。結果(`ts_not_converged`、`monotonic_uphill`、`afir_not_converged`)は v3 と同じで、`stderr.txt` に原因の 1 行が加わった(`run_tsopt_task: No more negative eigenvalues, optimization failed.`、`run_nt2_task: No transition state guess was found in Newton Trajectory scan.`、`run_afir_task: Problem: AFIR optimization did not converge.`)。v3 では後の 2 件の `stderr.txt` は空だった。
- **C40**: `smiles_to_molecule("[2H]F", 0, 1)` は `SMILES '[2H]F' specifies isotopes; hfauto uses most-abundant-isotope masses` を出した。
- **C28・C34**: 4 run には層の欠けた subject も conformers stage もないので、実 run では通っていない(C28 は S9 のテスト、C34 は第4ラウンドの CREST のプローブで確認済み)。
