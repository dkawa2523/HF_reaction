# 実計算による検証の記録

round 9(r9-PRE〜S6 と FINAL。結果報告 [2026-10-02_round9_result.md](reviews/2026-10-02_round9_result.md))の後の全体の再検証の記録である。今の記録は、コミット済みの HEAD(`4137beb`)で全系を新しい run dir から流し直した **VAL9/R2** で、FINAL の作業ツリーで流した VAL9 と VAL7(round 7)を比べる。期待値は受け入れた最新の再生(CUR)と round 9 の検証 run の値で、VAL7 から変えた期待値は §4 に理由とともに挙げる。wave ごとの詳しい記録は `git show c6eff84:docs/validation.md`(§10〜§25)、round 7 は `git show 1bd2ff5:docs/validation.md` にある。決定表の行番号は今の 14 行の表([design.md](design.md) §6.2)のもの。

## 1. 条件

- コード: VAL9/R2 は `4137beb`(2026-10-02 07:49〜14:22)。最初の run(`hcn`)の `git describe` だけが `-dirty` なのは WSL の git から見た CRLF の差で、`git diff --ignore-cr-at-eol` は空(内容は HEAD と同一)。VAL9 は `c6eff84` に FINAL の作業ツリーを足したもの(01:25〜07:37)。
- 環境: WSL2 Ubuntu 24.04(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 2,000 MB、常に `--bind-to none`。xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0、pymsym 0.3.5)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(I と Te は def2-ECP)、順位のエネルギー層は M06-2X-D3(0)/def2-TZVPD、低レベルは GFN2-xTB。298.15 K・1 atm、順位の量は δG_eff(kcal/mol)。
- run: `/home/user/hfauto_r9/VAL9/R2/<run>`(VAL9 は `/home/user/hfauto_r9/VAL9/<run>`)に新しい run dir で 1 本ずつ。`/home/user/hfauto_r9/runcase.sh`(共有のロック、`/usr/bin/time -v`、CAP は SIGTERM で送り 60 s 後に SIGKILL、孤児の検査、status と report)。手順は `/home/user/hfauto_r9/val9r2_chain.sh`(VAL9 の `val9_chain.sh` と同じ系・同じ CAP。CAP は既定 3,600 s、R1〜R4 は 1,500 s、S10・S16・S19・S6 は 5,400 s、シュウ酸は 7,200 s)。
- 対象(41 run と smoke): VAL7 の 35 run(noscreen とパネルの pipeline は `/home/user/hfauto_r9/pipelines/` の写しで、既定と同じ M06-2X の層を持つ)と、round 9 で足した 6 run(`sn2_cl_d3h`、`malonaldehyde_c2v_undeclared`、`h_c2h4`、`bh3_nh3`、ニトロメタン(`g3_soft.yaml`)、`hono_panel`)。
- 比べ方: 結論の水準で比べる(outcome、生成物の状態の集合、各段の区分、δG_eff)。並列の NWChem と CREST は run ごとに軌跡が分かれる(design.md §10)ので、id ではなく化学的な署名で組む(`/home/user/hfauto_r9/FINAL/concl.py`)。期待の分岐に入らなかった run は未達として記録し、期待値・しきい値・予算は動かさない(分析 §4.4)。
- 時間: 一意のジョブ鍵ごとの core 秒(duration × `-np`、`tools/jobtime.py`)と歩数。壁時計は参考値。比較の出力は `/home/user/hfauto_r9/VAL9/R2/cmp/`。

## 2. 実計算で通した範囲

| 項目 | DFT(PBE0/SVPD)まで通したもの | 低レベルだけ |
|---|---|---|
| 元素 | H、B、C、N、O、F、Cl、I(ECP)、Te(ECP、H2Te) | S(SO2·NMe3)、Fe(FeCl3·CH4) |
| 電荷 | 0、−1(Cl⁻・I⁻ の SN2) | — |
| スピン | 一重項、二重項(CH3O•、H + H2、OH···CH4、H + C2H4)、三重項(O2)、低スピン結合の二重項(CH3• + O2) | 六重項(FeCl3·CH4) |
| 反応型 | 1,2-H 移動(HCN、CH3O•、HCOH → H2CO)、1,3-H 移動、ねじれ(HONO、シュウ酸の 2 段)、縮退転位(NH3、DME、H2O2、PT、恒等 SN2、H 交換、(HF)₂)、未宣言の縮退転位(鞍点だけから)、H 引き抜き(OH + CH4)、障壁なし(H2O·HF の反転、水二量体の受容体交換)、分離した単量体からの会合(CH3 + O2、H + C2H4、BH3 + NH3。いずれも PBE0 で障壁なし) | ハロゲン結合、配位付加体(SO2·NMe3) |

通していないもの: 陽イオン、−2 以下の陰イオン、遷移金属の DFT、開殻一重項、溶媒(design.md §10)。

## 3. 結果(VAL9/R2)

判定は、期待と同じなら「合格」、違えば理由を書いて「未達」とする。δG_eff は 0.01 kcal/mol 以内を同じとみなす(S6 は §3.2 の幅)。「VAL9」欄の「同」は `concl.py` で結論・δG_eff が VAL9 と同一(差 ≤ 0.0001)を示す。時間は壁時計 / QM core 秒。

**まとめ: 41 run と smoke のすべてが期待どおり(合格 41/41)。** VAL9 と比べて 40 run が結論の水準で同一、差は S6 の TS の点群(Cs と C1 の入れ替わり、δG_eff が RT ln 2 だけ動く。期待の幅の中)だけ。rc 1 は期待どおりの 6 本(拒否・止める設計の 3 本、低レベルのジョブの失敗を持つ S5・S6、層の SP が収束しない `h_c2h4`)。

### 3.1 known_endpoints と回帰

| 系 | 型 | 期待(outcome、δG_eff) | VAL9/R2 の実測 | VAL9 | VAL7 | 時間 | 判定 |
|---|---|---|---|---|---|---|---|
| R1 `hcn` | 1,2-H 移動 | elementary 41.57 | 41.57(ν1 −1128.5i、χ 0.644) | 同 | 42.18 | 1:09 / 239 | 合格 |
| R2 `hono` | ねじれ | elementary 9.88 | 9.88(−681.4i) | 同 | 11.82 | 2:47 / 602 | 合格 |
| R3 `nh3_inversion` | 縮退転位 | degenerate 4.38 | 4.38 | 同 | 3.84 | 0:44 / 159 | 合格 |
| R4 `water_same_basin` | 同じ basin | same_basin | same_basin | 同 | 同 | 0:11 / 37 | 合格 |
| `ch3n_spin` | 拒否 | `spin_crossing_reaction_unsupported`、ジョブ 0 | 同、rc 1 | 同 | 旧い文言 | 0:01 / 0 | 合格 |
| `formaldehyde` | 1,2-H 移動 | HCOH → H2CO 31.29、逆は out_of_window | 31.29、out_of_window | 同 | 29.02 | 1:53 / 411 | 合格 |
| `h2o_hf_inversion` | 障壁なし | barrierless 0.0(`screen_midpoints:barrierless`) | 同 | 同 | 同 | 0:49 / 175 | 合格 |
| `water_dimer_as` | 障壁なし | barrierless 0.0 | 同 | 同 | 同 | 1:06 / 239 | 合格 |
| `nh3_planar_seed` | 宣言した鞍点 | `follow1:ts_candidate`、4.38、paths の misses 0 | 同(paths hits 2・misses 0) | 同 | 3.84 | 1:06 / 254 | 合格 |
| `dme_c2v_seed` | 高次の鞍点の種 | `follow1:one_side`、degenerate 2.13 | 同 | 同 | 1.96 | 13:39 / 3,163 | 合格 |
| `hcn_noscreen` | FIND_PATH | `path_hei` から R1 と同じ TS、41.57 | 41.575(−1128.2i) | 同 | 42.18 | 2:34 / 584 | 合格 |
| `hf_dimer_swap_noscreen` | 縮退転位 | degenerate 0.80 | 0.80 | 同 | 1.39 | 1:50 / 425 | 合格 |
| `hono_walltime` | 行 7 | `walltime`、report あり、孤児 0 | 同(unresolved_within_budget) | 同 | 同 | 0:53 / 201 | 合格 |
| S1 `s1_sn2_cl` | 恒等 SN2 | degenerate 14.53 | 14.53(−393.6i、χ 0.930) | 同 | 11.08 | 4:48 / 1,106 | 合格 |
| S2 `s2_sn2_i` | 恒等 SN2(ECP) | degenerate 9.61 | 9.61(χ 0.836) | 同 | 6.88 | 7:01 / 1,625 | 合格 |
| S3 `s3_h2te` | ECP の極小 | Te の DFT 極小 | 同 | 同 | 同 | 0:15 / 56 | 合格 |
| S4 `s4_ch3o_doublet` | 開殻 1,2-H 移動 | elementary 30.26 | 30.26(χ 0.688) | 同 | 30.61 | 6:52 / 1,492 | 合格 |
| S12 `s12_dme_rough` | ねじれの縮退 | degenerate 2.13 | 2.13 | 同 | 1.96 | 7:00 / 1,610 | 合格 |
| S13 `s13_h2o2_gauche` | ねじれの縮退 | degenerate 1.29 | 1.29 | 同 | 1.18 | 1:03 / 234 | 合格 |
| S14 `s14_malonaldehyde` | PT | degenerate 0.93 | 0.93(χ 0.832) | 同 | 0.0 | 12:28 / 2,888 | 合格 |
| S15 `s15_hono_hno2` | 1,3-H 移動 | elementary 56.83 | 56.83(χ 0.640) | 同 | 52.02 | 2:45 / 591 | 合格 |
| S16 `s16_acac` | PT | degenerate 0.0(`submerged_barrier`) | 同(ν2 −36.9i は noise 内) | 同 | 同 | 1:12:12 / 16,823 | 合格 |
| S17 `s17_hf_dimer_swap` | 縮退転位 | degenerate 0.81 | 0.81 | 同 | 1.39 | 1:24 / 309 | 合格 |
| S18 `s18_h3_doublet` | H 交換 | degenerate 12.71 | 12.71(χ 1.000) | 同 | 4.09 | 0:32 / 109 | 合格 |
| `oxalic_two_step` | 2 段のねじれ | multi_step、split1 11.81(`ts_calc`)、split2 12.71 | 同 | 同 | 12.78、13.76 | 40:27 / 9,417 | 合格 |
| `malonaldehyde_c2v_undeclared` | 未宣言の PT | `follow1:ts_candidate` → degenerate 0.93、`qrc1:minus_is_image` | 同 | 同 | — | 16:32 / 3,849 | 合格 |
| `sn2_cl_d3h` | 未宣言の SN2 | degenerate 14.53、`qrc1:minus_is_image`、paths の misses 0 | 同 | 同 | — | 4:58 / 1,150 | 合格 |
| `h_c2h4` | 会合 | `scan:barrierless` 0.0(PBE0 の山 +0.60 < 1.0) | 同。錯体の M06-2X の SP が `scf_not_converged` で rc 1(§10) | 同 | — | 40:31 / 9,485 | 合格 |
| `bh3_nh3` | 会合 | barrierless 0.0(`known:basin_adduct`) | 同 | 同 | — | 4:38 / 1,072 | 合格 |
| ニトロメタン | soft な極小 | `soft:resolved`、minimum | 同(ν1 +36.2) | 同 | — | 3:58 / 921 | 合格 |

### 3.2 discover

| 系 | 期待 | VAL9/R2 の実測 | VAL9 | VAL7 | 時間 | 判定 |
|---|---|---|---|---|---|---|
| S5 `s5_ch3_o2` | 会合 CH3 + O2 → CH3OO•: `scan:barrierless`(AP、`scan:ap:4`)0.0。CH3OO• → CH2O + OH: multi_step、split1 elementary 46.85、split2 は O–O 開裂の TS を期待しない。split2 が unresolved なら親と split2 は `thermo_unavailable` | 同(split1 −1866.7i・χ 0.775、split2 は χ 0.273 の鞍点 2 つを棄却)、rc 1 | 同 | 結論なし | 43:48 / 10,060 | 合格 |
| S6 `s6_oh_ch4` | OH + CH4 → CH3 + H2O: elementary、一次の TS(ν1 約 −490i、ν2 > −50i、χ ≥ 0.5)、δG_eff 5.5〜6.0(TS が C1 なら Cs より RT ln 2 低い)。CH3 + H2O ↔ CH3OH + H: multi_step、split1 `same_as`、split2 約 44.8 | elementary 5.93(−488.5i、ν2 +45.5、χ 0.675、TS は Cs)。multi_step、split1 `same_as` 20.98、split2 44.83(−1347.5i、χ 0.798)、rc 1 | 5.54(C1)、20.60、44.82 | unresolved | 31:00 / 6,713 | 合格 |
| S7 `s7_nh3_icl`(`--to explore`) | 生成物 0 | 生成物 0(trial 2) | 同 | 同 | 0:09 / 24 | 合格 |
| S8 `s8_so2_nme3`(`--to explore`) | trial 19、生成物 0 | 同 | 同 | 同 | 1:09 / 257 | 合格 |
| S9 `s9_undeclared` / `s9_conformers` | `declare_multiplicity` / 六重項で CREST rc 0 | 同(rc 1 / rc 0) | 同 | 同 | 0:01 / 0:30 | 合格 |
| S10 `s10_amine_pilot2` | 二重 H 交換 degenerate 36.64、C3H9N・C3H10FN は `not_reacting` | 同(χ 0.856) | 同 | 32.29 | 3:24 / 778 | 合格 |
| S19 `s19_tma_hf2` | same_basin、R6 の seed なし、層の SP と ΔG_assoc なし | 同(CREST の配座 6) | 同 | seed あり、−11.53 | 27:48 / 6,490 | 合格 |
| `s5_undeclared`(`--to screen`) | CH3·O2 だけが `declare_multiplicity`、単量体は screen まで(rc 1) | 同 | 同 | 同 | 0:01 / 0 | 合格 |

S6 の 2 回の新しい run(VAL9 と R2)は同じ引き抜きの TS(同じ添字で RMSD 0.033 Å、ν1 −485.3i と −488.5i、χ 0.674 と 0.675)に届き、点群の判定だけが C1(m = 2)と Cs に分かれた。これは結果報告の残る課題 8(ほぼ Cs の平らな TS の点群が雑音で入れ替わる)の再現で、δG_eff の差 0.39 は RT ln 2(0.41)と一致する。

### 3.3 手法パネルと精度

- CCSD(T)/def2-TZVPD との差(5 反応、[design.md](design.md) §7.4 の表。round 9 の実測): 障壁は M06-2X で平均 0.83、最大 1.79 kcal/mol(PBE0/SVPD は 2.55、6.2)。ΔE_rxn は M06-2X で平均 1.17、最大 2.30 で、PBE0/SVPD(1.05、2.25)より良くない。
- composite の δG_eff(同じ run に `energy_method: ccsd-t_def2-tzvpd` の thermo を足したもの)との差: SN2 +0.46(M06-2X 14.53、CCSD(T) 14.07)、HCN −1.79、CH3O• −0.72、HONO +0.09。
- マロンアルデヒドの PT の ΔE‡: CCSD(T) 3.47、M06-2X 2.69、PBE0/SVPD 2.01 kcal/mol。文献の CCSD(T)/CBS 約 4 より 0.5 低い。S6 の分離基準の ΔE‡: M06-2X 5.73〜5.81、CCSD(T) 6.96、PBE0/SVPD 0.75〜0.80。
- VAL9/R2 のパネル 4 本(`s20_sn2_panel` 10:51 / 2,533、`s20_hcn_panel` 1:23 / 319、`s20_ch3o_panel` 4:28 / 1,037、`hono_panel` 7:39 / 1,778): `method_panel.csv` は VAL9 と 1e-5 kcal/mol 以内で同一。ΔE‡ の M06-2X − CCSD(T) は HCN −1.79、SN2 +0.46、HONO −0.01、ΔE_rxn は HCN −2.30、HONO −0.50。SN2 の 400 K・1 M は 14.89 で report と ranking.csv が一致。いずれも design.md §7.4 の主張(障壁の差 最大 1.8、ΔE_rxn 最大 2.3)の中。

### 3.4 実行基盤、TS の性格、再現性

- smoke(`pytest -m real`、13 件): 13 件合格(66 s)。孤児(nwchem、prterun)は 41 run すべてで 0。
- TS の性格(`tools/chi.py --src VAL9/R2`): 順位のある一次の TS 23 本の χ は 0.571〜1.000(下限 0.3 を下回る接続した TS は 0)、棄却した鞍点は S5 の split2 の 2 本(χ 0.273)。VAL9 との差は S6 の TS(ν1 −485.3i → −488.5i、ν2 +48.9 → +45.5)、`hcn_noscreen`(−1128.7i → −1128.2i)、S5 の棄却した鞍点(χ 0.273・0.275 → 0.273・0.273)だけ。
- 到達(`F2/reach.py` の決定表の行・理由・注記): VAL9 と同じ。差は S6 の `screen_neb:nonzero_exit` が R2 にないことと、シュウ酸の `ts_calc` の計算 id だけ。
- 再現性: 同じコードの新しい run 2 回(VAL9 と R2)で、ジョブ数(812)、opt・saddle の歩数(差は S6 の DFT opt の 442 → 428 と saddle の 231 → 228 だけ)、SCF 数が一致した。VAL9 の自己一致の再生(`tools/replay.py FINAL --src VAL9 --strict`、40 run)は 40/40 PASS で、CUR はその再生(`FINAL/replay`)を指す。
- 品質ゲート(`4137beb`): WSL の pytest 668 passed、ruff・pyrefly 0 errors・import-linter 5 kept・`radon cc -n D` は空、最大のファイル 481 行。

## 4. 期待値の変更とその理由(VAL7 → VAL9)

| 系 | VAL7 | 今の期待 | 理由 |
|---|---|---|---|
| 順位のあるすべての反応 | PBE0/SVPD の δG_eff | M06-2X 層の δG_eff | 既定の順位を M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD にした(design.md §7.3)。値は PBE0 の非局在化誤差の向き(H 移動・求核置換で上がる)に動く: HCN −0.61、HONO −1.94、NH3 +0.54、SN2 Cl +3.45・I +2.73、S10 +4.35、HONO → HNO2 +4.81、H + H2 +6.35、ねじれと配座は ±1 前後 |
| S18 H + H2 | 4.09 | 12.71 | 上の層(+6.35)に加え、ほぼ直線の H···H2 を点群 C∞v で直線として扱う(+2.27) |
| S17、`hf_dimer_swap_noscreen` | 1.39 | 0.81、0.80 | TS の点群を 1 つの判定で C2h(σ 2)に決め、層で −0.59 |
| S14 マロンアルデヒド | 0.0(`submerged_barrier`) | 0.93 | M06-2X の障壁(ΔE‡ 2.69)は PBE0(2.01)より高く、ZPE で沈まない |
| R2 HONO | 11.82 | 9.88 | PBE0 のねじれの障壁の過大(+2.0)が消える。層で cis が trans より 0.09 低くなり、同じ状態の G_R がその最小になるので δG_eff は ΔG‡(9.79)より大きい。cis/trans の順は期待にしない(差が 1 kcal/mol 未満) |
| S5 CH3• + O2 | 結論なし(BS 錯体を端にした問い) | 会合は barrierless 0.0、順位あり | 障壁のない会合では反応物は漸近で極小ではない。分離した単量体から付加体への緩和スキャンで問い、低スピン結合のスキャンは AP のエネルギーで判定する。ΔH₀ −31.63 は実験の D₀ 約 32 と 0.4 以内 |
| S5 CH3OO• → CH2OOH• | reassigned 32.67 | multi_step(split1 46.85、split2 は結論を期待しない) | 結合が変わる case を化学状態で判定するので、TS の側(CH2OOH•)は中間体になる。O–O 開裂の領域の鞍点(χ 0.27〜0.29)は反応モードでないので棄却する。REACTION_MODE_MIN は変えない |
| S6 OH + CH4 | unresolved(−73.0i の回転子の鞍点) | elementary、5.5〜6.0 | χ で回転子の鞍点(χ ≤ 0.023)を QRC の前に棄却し、ρ だけを負の曲率にした初期 Hessian で引き抜きの TS に届く。TS の点群が Cs か C1 かで RT ln 2 動くので幅で期待する |
| S6 CH3 + H2O ↔ CH3OH + H | reassigned 30.66 | multi_step(CH4 + HO を経る) | 同じ状態の組の仮説を 1 つにまとめ、TS の側を状態で判定する。子 CH3 + H2O → CH4 + HO は引き抜きと同じ問いなので結論を写す |
| S19 TMA·(HF)₂ | seed の DFT opt で `collapsed_at_dft_from_seed`、ΔG_assoc −11.53 | seed なし、ΔG_assoc なし | seed と basin のラベルの違い(N···H)は結合変化の帯(±0.1 Å)の中。順位の点を持たない反応には層の SP を計算しない |
| `h_c2h4`(新) | — | barrierless 0.0 | 計画の期待(約 2 kcal/mol の山)は満たさない。PBE0 の山は +0.60 kcal/mol で解像度未満(参照の障壁 1.72)。resolution とスキャンの点は変えていない |
| `ch3n_spin` | `endpoints differ in charge, multiplicity or atom order` | `spin_crossing_reaction_unsupported` | 多重度だけが違う端点に固有の理由を付けた |

## 5. 到達表(今のコード)

VAL9・VAL9/R2 と round 9 の検証 run の case log と diagnostics から数えた(`/home/user/hfauto_r9/F2/reach.py`)。

| 対象 | 通ったもの | 通っていないもの |
|---|---|---|
| 決定表の行 | 1(宣言した平面 NH3、mode_follow 0 の検証 run)、2(水、S19)、3(formaldehyde)、4(多数)、5(oxalic、S5、S6)、6(SCREEN と中点、会合のスキャン)、7(`hono_walltime`)、9・10・12(多数)、11(S5 の split2)、13(noscreen 2 本、S5)、14(S5 の split2) | 8 |
| 駆動しない子 | `same_as`(S6 の split1)、`split_depth`(`max_split_depth: 0` の oxalic の検証 run) | — |
| saddle の種 | `screen_ts`、`screen_hei`、`discovery_ts`、`path_hei` | `higher_order_retry`(M4 の S6 の実 run だけ)、近道の種がすべて失敗した後の NEB |
| 初期 Hessian | xTB(`xtb:mode`、`xtb:rho`)、DFT(noscreen 2 本) | TS freq(M4 の S6 だけ) |
| 検証 | `ts_calc`、`ts_rejected:not_reaction_mode`(S5 の split2) | `ts_rejected:higher_order`、`higher_order:not_stationary`、`saddle:above_path_bound` |
| QRC | `minus_is_image`、`end<i>_to_new_basin` | 2 倍の振幅の 2 回目の QRC |
| 極小 | `ts_candidate`、`one_side`、`image_of`、`known`、`not_reacting`、`soft:resolved` | 非停留の極小の Newton の押し出し |
| 開殻 | 端の SCF の解からのプロファイルの SP、AP(S5)、SCF の救済(harness) | 2 つの guess の分岐、実計算の `scf_branch_jump` |
| outcome | blocked_upstream 以外のすべて(blocked_upstream は行 1 の検証 run) | — |

## 6. 計算量

| 比較 | 結果 |
|---|---|
| VAL9/R2(41 run) | QM 89,327 core 秒(812 ジョブ)、壁時計の和 23,299 s、連鎖全体 6 h 33 min |
| VAL9(41 run) | QM 84,471 core 秒(812 ジョブ)、壁時計の和 22,067 s。R2 との差 +4,856 のうち +5,530 は `h_c2h4` の収束しない SP(救済の試行が 470 s → 1,851 s。§10)。それを除く 40 run は −0.8% |
| VAL7(35 run) → VAL9/R2 の対の 35 run | 78,752 → 71,071 core 秒(−9.8%、M06-2X の層を含む)。減った主なものは S19(−9,708、seed の DFT がない)、S5(−2,434)、S6(−2,252)、増えたのは層の重い acac(+4,316)とシュウ酸(+1,751) |
| エネルギー層 | VAL9/R2 で M06-2X の SP 82 本・18,257 core 秒(`h_c2h4` の錯体の 7,423 を含む)。1 点 10〜60 s が多く、シュウ酸約 140 s、acac(15 原子)約 600 s |
| HCN のクイックスタート | 壁時計 1:11、QM 242 core 秒(22 ジョブ。`/home/user/hfauto_r9/F2/hcn_quickstart`) |
| 残る重複 | DFT の mode-follow は対称な鞍点でも ± 両側を opt と freq にかける(マロンアルデヒドで freq 4 本が QM 時間の 69%、`sn2_cl_d3h` で像の側が 39%) |

rank 2 の site での `--bind-to none` の修正(M4)より前の rank 2 の時間は、QRC の opt で約 5 倍に膨らんでいるので比べない。

## 7. 再生ゲート

振る舞いを保つ整理の合格基準(分析 §3 X8)。記録済みの run を今の作業ツリーで JobStore から再生し、新しい QM ジョブが 0 件で、結論の記録が変わらないことを確かめる。

- 道具(repo の外、`/home/user/hfauto_r9/`): `tools/replay.py`(本体)、`tools/runs.json`(再生する run と pipeline の連鎖)、`runcase.sh`(1 本ずつの実行)、`tools/jobtime.py`。
- 手順: source を `<wave>/replay/<name>` に複製し、pipeline の最初の stage から再ステージする。パネルの run は続けて `RESUME=1` で panel pipeline を再開する。
- 比べるもの: 再生した stage ごとの `run_state.json`(status、hits、misses、failures_by_kind)、manifest の結論の記録(species、minimum、discovery、reaction、species_thermo、reaction_thermo、report)を欄の単位で、複製の中で走ったジョブ。reason の文字列、FileRef、created_at と壁時計は比べない。
- 合格(`--strict`): 再生したすべての stage で misses 0、source でジョブを走らせた stage はすべて再生して hits > 0、status と記録が同一。振る舞いを変える wave は `--allow-new` で挙げた run の差を化学的な理由で説明する。
- CUR: 受け入れた最新の再生への symlink(`/home/user/hfauto_r9/CUR/<name>`、今は `FINAL/replay`)。次の wave は `--src CUR` で再生する。
- 除外: `hono_walltime`(Deadline は実時間なので、キャッシュの hit で行 7 の発火が変わる)、`s9_*`(遷移金属の DFT は未検証)。
- 限界: 保証するのは同じ鍵に同じ結果が返る範囲だけ。新しいジョブが 1 本でも出ると、その下流は run 間の比較(結論の水準)になる。古い JobStore の Evidence には勾配がなく、Failure にはエネルギーがない。

## 8. プローブと実 run で決めたこと

| 決めたこと | 根拠 |
|---|---|
| 反応モード性の下限 `REACTION_MODE_MIN` = 0.3 | 接続した TS 28 件は 0.571〜1.000(重原子の結合変化を含むものは ≥ 0.776)、接続しない一次の鞍点は ≤ 0.023(S6 の回転子、S5 の −72i)。最も近い棄却は S5 の CH2OOH• の領域の 0.22〜0.29(`M3/chi/`、`tools/chi.py`) |
| ρ だけを負の曲率にした初期 Hessian | VAL7 の S6 の種のモデルは負の固有値 −0.0028、ρ との cos 0.345。新しいモデルは −0.2187、cos 1.000。S6 を 3 種 × rank 2・4 で流した 6 run がすべて引き抜きの TS に届いた(`M4/real/`) |
| 原子の対応(basin の構造を WL クラスで端点の添字に並べる) | S6 の xTB の IRC の端は PBE0 の basin から 0.957 Å 離れ、元素だけで並べると結合変化が 8 本(IRC 自身は 4 本)。端の監査の不一致 6/17 → 0(`M1/audit/`) |
| 結合変化の帯 `RESOLVED_A` = 0.1 Å | S19 の N···H は GFN2 と PBE0 で 0.15 Å 違う。帯でラベルだけが割れる組は 11 組、本物の結合変化は帯の外(最小は acac の O–H の +0.175 Å)(`S3/S3b/`) |
| 停留性の基準 5e-5 Eh(basin のエネルギー幅と同じ) | VAL7 の DFT 極小 64 点の ΔE_N は ≤ 1.2e-5、W3 の肩(−6.5i)だけが 2.66e-4(`tools/stationarity.py`) |
| 点群の受理を basin の基準にする | 1e-3 Å の揺らぎ 20 回で 110 subject の点群・σ・m・直線性が不変、既知の点群 6 つが一致、S18 の 2 run の差 2.27 → 0.006 kcal/mol(`tools/symcheck.py`、`M7/`) |
| 会合のスキャン(付加体 + 1.5 Å から 8 点) | 1.5 Å は Coulson–Fischer 点(CH3·O2 の C–O で約 2.2 Å)の外。S5 の BS の山 +0.10 は解像度未満(`M8/real/s5_ch3_o2`) |
| 開殻のプロファイルの SP を端の SCF の解から直接解く | 端の vectors から解いた点は逐次の継続と 5e-8 Eh 以内で一致し、山は +0.80。atomic guess では +28.9 の偽の山(`probes/G2-1v/`) |
| 枝の跳びの判定に spin_tol を使う | 本物の TS の最大点と隣の ⟨S²⟩ の差は ≤ 0.022、SCF の枝の跳びは 0.77〜0.95(`S4/audit/`) |
| 再開の上限(プロファイルの最大 + resolution) | VAL7 の S5 の再開は最大より 27〜32 kcal/mol 上から、収束した再開は 21〜69 kcal/mol 下から始まっていた(`tools/bound_check.py`) |
| trial の番号付けの不変性 | 原子を無作為に 20 回並べ替えても drive がすべて同じ(`S5/int/trialperm_*.txt`)。S8 の新しい run で trial_id の列と生成物 0 が一致 |
| エネルギー層 M06-2X-D3(0)/def2-TZVPD | 4 反応の ΔE‡ の CCSD(T) との平均絶対誤差 1.03(PBE0/SVPD 2.68、ωB97X-D3 1.56)。grid fine と xfine の差 ≤ 0.01(r8 のプローブ G6-1) |
| WFT のメモリ(2,000 MB/rank、GA に 70%)と `freeze <n>` | `memory total 1200 mb` ではマロンアルデヒドの (T) が配列を確保できなかった。HI は `freeze 4`、I⁻···CH3I は `freeze 9`(`S1/real/`) |
| autoz の失敗は最後の frame から Cartesian で続ける | S10 の QRC の側: 始点からの 19 歩・42 s → 1 歩・3.3 s、同じ basin(`S1/harness/autoz`) |
| NWChem は常に `--bind-to none` | rank 2 の site で QRC の両側が同じコアを取り合い、1 歩 30 s(rank 4 では 6 s)(`M4/real/`) |
| 一次の鞍点の Hessian は正定値のモデル、高次の鞍点はそのまま | 負のまま渡すと下りが遅い(acac の QRC 34/27 → 15/14 歩)。高次の鞍点からの側はモデルにすると這う(DME C2v で 207 歩でも未収束)(round 7) |
| 極小の stage は直列 | 4 つの opt を 1 rank ずつにすると TMA 105 → 206 s(round 7) |
| AFIR を削除 | 41 件すべて陰性で生成物 0(round 7) |
| 採らなかったもの | χ の下限を下げる、2 回目の higher_order_retry、tight の収束、noise_cm1・saddle_cm1 を動かす、ωB97X-D3 の層、CREST の電子温度を上げる、linopt 0 |

## 9. 撤回・訂正した値

- round 7 までに撤回したもの: HONO の ΔG‡ 24.11(ZPE の二重計上)、HCN の ΔG_rxn 18.05(内部フォールバック)、I⁻ + CH3I の ΔG‡ 7.53(対称数)、R2 12.2265(m = 2 を数えていない)、S4 の未宣言のねじれ 3.878。
- W3 の S6 の引き抜き δG_eff 1.48(−481.1i)は、平らな PES の肩(−6.5i、ΔE_N 2.7e-4 Eh)を端にした PBE0 の値で、今のコードの記録にはない。今の値は 5.5〜6.0(M06-2X)。
- S6 の VAL9 の 5.54 を「別の鞍点」とした結果報告の読みは補う: VAL9/R2 の TS は VAL9 の TS と RMSD 0.033 Å の同じ鞍点で、点群が Cs と判定されて 5.93 になった。値の差は点群の判定だけ(§3.2)。
- VAL7 の S5 の「結論なし(多参照性のある二重項の結合)」: 主因は多参照性ではなく問いの形だった(§4)。
- M4 で記録した S6 の「平らな回転子の領域(ν2 −77i〜+24i)」は、停留していない点の Hessian で測った値だった。
- 分析の「既定の mode_follow でも宣言した鞍点の端点は行 1 に落ちる」は成り立たない: 平面 NH3 は mode-follow で NH3 の basin に加わり、same_basin になる。
- 計画の凍結軌道の例(I は 9、I⁻ + CH3I は 19、CH3O は 1)は NWChem の `freeze atomic` の規則と合わない(I は 4、I⁻···CH3I は 9、CH3O は 2)。

## 10. 実証していないこと・残る課題

- 層の SP の SCF(結果報告の残る課題 6): `h_c2h4` の H···C2H4 錯体(二重項)の M06-2X/def2-TZVPD の SP は、新しい run 2 回とも atomic guess と cgmin の救済の両方で収束しない(rc 1)。錯体の構造は 2 回で 1e-5 Å しか違わないが、救済の試行は 470 s と 1,851 s(`ga_iter_lsolve` の停滞が 86 回)で、run の QM 時間が 2.4 倍になった。層の SCF は freq の vectors から始めず、救済の試行には反復数のほかに上限がない。barrierless の順位(0.0)は変わらないが、順位の点で起きれば `thermo_unavailable` になる。
- TS の点群の揺らぎ(残る課題 8): ほぼ Cs の平らな TS の点群が run ごとに Cs と C1 で入れ替わり、δG_eff が RT ln 2 動く(S6: 5.54 と 5.93)。受理の許容は basin の基準のまま。
- 決定表の行 8 は今のコードの実計算で通っていない。近道の種がすべて失敗した後の NEB、2 倍の振幅の 2 回目の QRC、ケースの例外の閉じ込めは単体テストだけ。
- 極小の Newton の押し出しと TS の `higher_order:not_stationary` は単体テストとプローブだけ。再開の上限は QM なしの評価だけ。
- 中点の密化が隠れた山を見つけた例はない。会合のスキャンの山から REFINE_SADDLE に進む分岐は単体テストだけ(DFT の山が 1 kcal/mol を超える会合が検証セットにない)。
- 開殻: 両端の ⟨S²⟩ が違うプロファイル(2 つの guess)、実計算での `scf_branch_jump`、`spin_ok` を通る SCF の救済、四重項の CH3·O2 で AP がかからないこと。
- TS の停留性: 一次の鞍点の受理には停留性をかけていない(40 cm⁻¹ 未満の大振幅のモードでは二次モデルが成り立たない。acac の TS で 8.6e-3 Eh)。
- 3 次以上の鞍点の `off_saddle`、`max_trials_per_source` の境界で同値の類を落とす規則は単体テストだけ。3 原子以上の直線分子の NT2 の開始構造は番号付けに依る。
- 縮退転位の第 2 候補(ギ酸二量体 D2h)、陽イオン(HOC⁺)、遷移金属・Te を含む系の M06-2X の SP、ECP 原子を含む CCSD(T) の参照との比較。
