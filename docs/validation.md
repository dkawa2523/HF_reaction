# 実計算による検証の記録

round 9(r9-PRE〜S6 と FINAL。結果報告 [2026-10-02_round9_result.md](reviews/2026-10-02_round9_result.md))の後の全体の再検証の記録である。今の記録は、コミット済みの HEAD(`4137beb`)で全系を新しい run dir から流し直した **VAL9/R2** で、FINAL の作業ツリーで流した VAL9 と VAL7(round 7)を比べる。期待値は受け入れた最新の再生(CUR)と round 9 の検証 run の値で、VAL7 から変えた期待値は §4 に理由とともに挙げる。wave ごとの詳しい記録は `git show c6eff84:docs/validation.md`(§10〜§25)、round 7 は `git show 1bd2ff5:docs/validation.md` にある。決定表の行番号は今の 13 行の表([design.md](design.md) §6.2)のもの(round 10 の W1 で壁時計の行 7 を削除し、旧 8〜14 を 7〜13 に詰めた)。round 10 の wave ごとの記録は §11 にある。

## 1. 条件

- コード: VAL9/R2 は `4137beb`(2026-10-02 07:49〜14:22)。最初の run(`hcn`)の `git describe` だけが `-dirty` なのは WSL の git から見た CRLF の差で、`git diff --ignore-cr-at-eol` は空(内容は HEAD と同一)。VAL9 は `c6eff84` に FINAL の作業ツリーを足したもの(01:25〜07:37)。
- 環境: WSL2 Ubuntu 24.04(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 2,000 MB、常に `--bind-to none`。xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0、pymsym 0.3.5)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(I と Te は def2-ECP)、順位のエネルギー層は M06-2X-D3(0)/def2-TZVPD、低レベルは GFN2-xTB。298.15 K・1 atm、順位の量は δG_eff(kcal/mol)。
- run: `/home/user/hfauto_r9/VAL9/R2/<run>`(VAL9 は `/home/user/hfauto_r9/VAL9/<run>`)に新しい run dir で 1 本ずつ。`/home/user/hfauto_r9/runcase.sh`(共有のロック、`/usr/bin/time -v`、CAP は SIGTERM で送り 60 s 後に SIGKILL、孤児の検査、status と report)。手順は `/home/user/hfauto_r9/val9r2_chain.sh`(VAL9 の `val9_chain.sh` と同じ系・同じ CAP。CAP は既定 3,600 s、R1〜R4 は 1,500 s、S10・S16・S19・S6 は 5,400 s、シュウ酸は 7,200 s)。
- 対象(41 run と smoke): VAL7 の 35 run(noscreen とパネルの pipeline は `/home/user/hfauto_r9/pipelines/` の写しで、既定と同じ M06-2X の層を持つ)と、round 9 で足した 6 run(`sn2_cl_d3h`、`malonaldehyde_c2v_undeclared`、`h_c2h4`、`bh3_nh3`、ニトロメタン(`g3_soft.yaml`)、`hono_panel`)。
- 比べ方: 結論の水準で比べる(outcome、生成物の状態の集合、各段の区分、δG_eff)。並列の NWChem と CREST は run ごとに軌跡が分かれる(design.md §10)ので、id ではなく化学的な署名で組む(当時は `/home/user/hfauto_r9/FINAL/concl.py`。今は repo の `validation/check.py` と `hfauto/reporting/validation.py` の `compare`。§11)。期待の分岐に入らなかった run は未達として記録し、期待値・しきい値・予算は動かさない(分析 §4.4)。
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
| `hono_walltime`(W1 で置き換え。§11) | 旧行 7(壁時計の予算。削除) | `walltime`、report あり、孤児 0 | 同(unresolved_within_budget) | 同 | 同 | 0:53 / 201 | 合格 |
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
| S9 `s9_undeclared` / `s9_conformers`(今の s9_undeclared は system の読み込みで拒否される。§11) | `declare_multiplicity` / 六重項で CREST rc 0 | 同(rc 1 / rc 0) | 同 | 同 | 0:01 / 0:30 | 合格 |
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

VAL9・VAL9/R2 と round 9 の検証 run の case log と diagnostics から数えた(`/home/user/hfauto_r9/F2/reach.py`。行番号は今の 13 行の表)。

| 対象 | 通ったもの | 通っていないもの |
|---|---|---|
| 決定表の行 | 1(宣言した平面 NH3、mode_follow 0 の検証 run)、2(水、S19)、3(formaldehyde)、4(多数)、5(oxalic、S5、S6)、6(SCREEN と中点、会合のスキャン)、8・9・11(多数)、10(S5 の split2)、12(noscreen 2 本、S5)、13(S5 の split2。子の試行の予算は偽エンジンのテストでも通す) | 7 |
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

- 道具(repo の中。round 9 までは `/home/user/hfauto_r9/tools` の replay.py・runs.json・runcase.sh・jobtime.py と CUR の symlink): `validation/replay.py OUT --src DIR [--runs A,B] [--strict]`。OUT は `/home/user/hfauto_r10` の下で、各 pipeline を QM のロックの下で 1 本ずつ流す。連鎖は `validation/cases.yaml` の `chain`。
- 手順: source の run dir を `OUT/<name>` に複製し、連鎖の最初の pipeline をその最初の stage から再ステージし、後の pipeline(パネル)は続きから再開する。
- 比べるもの: 再生した stage ごとの `run_state.json`(status、hits、misses)と、manifest の calculation 以外の記録を欄の単位で。reason の文字列と report.html は比べない。
- 合格(`--strict`): 再生したすべての stage で misses 0、source でジョブを走らせた stage はすべて再生して hits > 0、status と記録が同一。
- 除外: `superseded` の case(`hono_walltime`。壁時計の予算は W1 で削除)は流さない。`s9_*` は case にしない(遷移金属の DFT は未検証、未宣言の多重度は読み込みで拒否)。
- 限界: 保証するのは同じ鍵に同じ結果が返る範囲だけ。新しいジョブが 1 本でも出ると、その下流は run 間の比較(結論の水準)になる。古い JobStore の Evidence には勾配がなく、Failure にはエネルギーがない。

## 8. プローブと実 run で決めたこと

| 決めたこと | 根拠 |
|---|---|
| 反応モード性の下限 `REACTION_MODE_MIN` = 0.3 | W2 から質量加重の χ(§11.2): 接続した TS 97 行は 0.727〜1.000、回転子・会合の鞍点 9 行は ≤ 0.0653、S5 の split2(O–O 開裂)は 0.635〜0.714。空白 (0.0653, 0.635) の中なので値は変えない(`validation/chi_table.py`)。それより前の Cartesian の χ では接続した TS 0.571〜1.000、S5 split2 は 0.22〜0.29 で棄却されていた(§3 の表の χ は Cartesian の値) |
| ρ だけを負の曲率にした初期 Hessian | VAL7 の S6 の種のモデルは負の固有値 −0.0028、ρ との cos 0.345。新しいモデルは −0.2187、cos 1.000。S6 を 3 種 × rank 2・4 で流した 6 run がすべて引き抜きの TS に届いた(`M4/real/`) |
| 原子の対応(basin の構造を原子クラスで端点の添字に並べる) | S6 の xTB の IRC の端は PBE0 の basin から 0.957 Å 離れ、元素だけで並べると結合変化が 8 本(IRC 自身は 4 本)。端の監査の不一致 6/17 → 0(`M1/audit/`) |
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

- 層の SP の SCF(結果報告の残る課題 6): W2 で解消。`h_c2h4` の H···C2H4 錯体の M06-2X の SP は、射影の guess では発散し、救済の (2)(3) で ⟨S²⟩ 0.7506 に収束する(B1 の sp stage は 2,019 → 311 s。§11.2)。経路の節点の SP の救済は未確認。
- ⟨S²⟩ の照合は、同じ多重度の別の空間状態を区別しない(S5 の qrc1_1。§11.2)。PBE0 の極小が励起の SCF 解にあることも検出しない。
- TS の点群の揺らぎ(残る課題 8): ほぼ Cs の平らな TS の点群が run ごとに Cs と C1 で入れ替わり、δG_eff が RT ln 2 動く(S6 の錯体基準の ΔG‡ で 5.54 と 5.93)。受理の許容は basin の基準のまま。
- S5 の split2 は境界の系: 同じ鞍点の QRC の側が run によって鞍点に落ち(`unassigned_side`、unresolved)、別の run では same_state から string で結ばれる(elementary 0.13)。並列の NWChem の非決定性で、W5-1 の下り(側が鞍点に落ちたとき)と VAL10 の 2 回の run で扱う(§11.3)。
- 決定表の行 7 は今のコードの実計算で通っていない。近道の種がすべて失敗した後の NEB、2 倍の振幅の 2 回目の QRC、ケースの例外の閉じ込めは単体テストだけ。
- 極小の Newton の押し出しと TS の `higher_order:not_stationary` は単体テストとプローブだけ。再開の上限は QM なしの評価だけ。
- 中点の密化が隠れた山を見つけた例はない(W3 から会合のスキャンも 1 回密化する。h_c2h4、bh3_nh3、S5 で山は出なかった)。会合のスキャンの山から REFINE_SADDLE に進む分岐は単体テストだけ(DFT の山が 1 kcal/mol を超える会合が検証セットにない)。
- 開殻: 両端の ⟨S²⟩ が違うプロファイル(2 つの guess)、実計算での `scf_branch_jump`、四重項の CH3·O2 で AP がかからないこと。
- TS の停留性: 一次の鞍点の受理には停留性をかけていない(40 cm⁻¹ 未満の大振幅のモードでは二次モデルが成り立たない。acac の TS で 8.6e-3 Eh)。
- 3 次以上の鞍点の `off_saddle` は単体テストだけ。explore の予算の切り方は偽の NT2 で示しただけ(ibuprofen、§11.4)で、実計算で予算に届いた系はない。
- 縮退転位の第 2 候補(ギ酸二量体 D2h)、陽イオン(HOC⁺)、遷移金属・Te を含む系の M06-2X の SP、ECP 原子を含む CCSD(T) の参照との比較。

## 11. round 10 の記録

計画は [roadmap.md](roadmap.md)、設計の差分は [design.md](design.md)。証拠は `/home/user/hfauto_r10/<wave>/`。

### 11.1 W1(段階 1 の前半: 予算、σ、多重度、検証の repo への取り込み)

- **予算は件数だけ**: 壁時計の予算(`walltime_h`、Deadline、行 7、`budget_exhausted`)を削除した。偽エンジンのテストで、(1) 兄弟の子 rx_split1 が saddle の試行 2 回を使い切って `attempts_exhausted` で閉じた後も、rx_split2 は 0 回から駆動されて elementary に届く(`test_a_split_child_queued_behind_a_sibling_that_used_its_budget_has_its_own`)、(2) reaction-paths の最初の QM ジョブの中で SIGTERM → paths は `incomplete`・rc 1・`status` も rc 1 → 流し直すと rc 0 で、paths の記録は止めなかった run と同一(`test_sigterm_leaves_the_stage_incomplete_and_a_rerun_resumes_it`)。R2 は件数の予算にほとんど届かないので、この 2 つは再生ではなく偽エンジンで示す。
- **σ を対称操作から数える**: 文献値(Fernández-Ramos 2007、Gilson–Irikura 2010)をテストに手で書いた: CH4 12、NH3 3、C6H6 12、SF6 24、アレン 4、S4 の C(OH)4 2、キュバン 24、ドデカヘドラン 60、1e-3 Å の揺らぎの D10h C10H10 20、CO2 2、HCN 1。旧コードは S4 で 4 を返し、D10h で KeyError だった。QM なしの thermo の再生(R2 の thermo 38 stage、103 の (run, subject))で点群・σ・m・直線性・対称化した構造は不変、thermo の記録 159 件はビット単位で同一(`/home/user/hfauto_r10/w1_2/`)。
- **状態の G を閉じる**: 層の SP か G を欠く極小を含む状態の G は None。R2 の `energy_layer_missing` は h_c2h4 1、hono_walltime 2、s19_tma_hf2 3、s3_h2te 1、s5_ch3_o2 1、water_same_basin 1 で、どれも反応に読まれないか、すでに `thermo_unavailable` の反応の参加者なので、R2 の結論は変わらない(計画が例に挙げた h_c2h4 の錯体は、どの反応も読まないので変化しない)。
- **多重度は必ず宣言**: 全 77 の同梱の化学種で宣言値は旧い暗黙の値と同じ。R2 の structures の SpeciesRecord 83 件(組成、多重度、状態ラベル、xyz の sha、指紋)は同一。基底状態の根拠: 三重項は O2(X³Σg⁻)だけ。二重項は CH3(X²A2″)、H(²S)、OH(X²Π)、CH3O/CH2OH、H3、H···C2H4 と C2H5。六重項は FeCl3(高スピン d⁵、⁶A1′)。ほかはすべて閉殻一重項(NH3、NMe3、HF、H2O、BH3、CH4、C2H4、CH3Cl、Cl⁻、I⁻ 錯体、ICl、SO2、H2Te、HCN/HNC、H2CO と trans-HCOH(三重項は約 25〜30 kcal/mol 上)、HONO、H2O2、エノール、シュウ酸、DME、錯体)。旧い暗黙の値が基底状態でなかった化学種はない(カルベン・ナイトレンは同梱していない)。
- **検証を repo に**: `validation/check.py /home/user/hfauto_r9/VAL9/R2` は 41/41 PASS(`hono_walltime` は superseded)と DEVIATION 1 件(h_c2h4 の barrierless は BH76 行 33 の 2.0 kcal/mol(分離のゼロ)と矛盾。合格に数えない)。BH76 の部分集合は P0e のファイルそのもの(sha 一致)。実 HCN のテスト(`tests/smoke/test_real_known_endpoints.py`)は 68 s で合格し、smoke は 14 件。
- **再生ゲート**: `validation/replay.py W1/replay_R2 --src VAL9/R2 --strict` と `W1/replay_VAL9 --src VAL9 --strict` はどちらも 40/40 PASS(新しいジョブ 0、記録は同一)。再生した run にも check.py は 40/40 PASS。統合で見つけた誤り: check.py は run dir のない case を MISSING と出すだけで rc 0 だった。superseded でなければ失敗にした(回帰テストあり)。

### 11.2 W2(段階 1 の後半: χ の質量加重、SCF の来歴、UHF-CCSD(T)、基準 B1)

- **χ(質量加重)**: `validation/chi_table.py /home/user/hfauto_r9/VAL9/R2 /home/user/hfauto_r9/VAL9 /home/user/hfauto_r7/VAL7 /home/user/hfauto_r9/M3/replay /home/user/hfauto_r9/M3/real /home/user/hfauto_r7/W5`(QM なし。証拠 `/home/user/hfauto_r10/W2/chi/all.txt`)で 154 行、一次の鞍点 113。接続した TS 97 行は 0.882〜1.000、ただし 2 行が 0.727(r7・W5 の S6 の多段の親 f9bc3a1cd8 を自分の 2 結合の変化で測ったもの)。回転子・会合の鞍点 9 行は 0.0054〜0.0653、S5 split2(O–O 開裂)は 6 行で 0.635〜0.714。事前の期待(接続 ≥ 0.88、split2 0.69〜0.72)とは 2 点で違うが、期待もしきい値も動かさず実測を記録した。`REACTION_MODE_MIN` 0.3 は空白 (0.0653, 0.635) の中にある。判定が変わるのは S5 split2 だけ(棄却 2 → 0。R2・VAL9・M3)。旧い r7・W5 の S6 の 790468f505(0.065、接続しない鞍点)は QRC でなく χ で棄却されるようになるが、結論は変わらない。縮退の判定(帯付きの変化の比較)は R2 の縮退の接続 16 件すべてで旧い 4 グラフの判定と同じ。
- **等価な原子の付け替え(U6-P2)は棄却**: 実装して R2 で測ると、結合の組を変えた 2 件(S6 の Walden の ts_calc: χ 0.940 → 0.726、S10 のリレー: ≥ 0.88 → 0.692)はどちらも誤った等価原子を選んだ([roadmap.md](roadmap.md) 変更 19)。
- **X3 を hfauto で**(`/home/user/hfauto_r10/W2/x3_p0a`、`W2/tools/x3_accept.py`): P0a の 5 点を R2 の freq を親にして流した。射影の deck(`vectors input project parent`)と出力の 'Orbital projection guess' はすべての (1) にあり、SP の鍵はすべて R2 と違う。⟨S²⟩ は親の類のまま(CH3O 0.7543 → 0.7544、S6 の TS 0.7578 → 0.7587、S5 の BS 錯体 1.7115 → 1.735 / 1.715 / 1.728、H···C2H4 錯体 0.7501 → 0.7506)。エネルギーは P0a の同じ段と ≤ 2.5e-9 Eh、HONO は atomic guess の解と 7.7e-9 Eh。H···C2H4 錯体は (1) が発散し(110 s)、(2) cgmin 31 回 + (3) で E −79.068928358 と P0a の (3) に 3.9e-10 Eh(合計 296 s。以前は救済が約 7.4k core-s 停滞)。ΔE_rxn(b3febc8052) は M06-2X −35.997、PBE0/TZVPD −34.386、ωB97X-D3 −34.714 kcal/mol(P0a −36.00 / −34.39 / −34.71)。CH3O の UHF-CCSD(T) は −114.873944099508 で P0d と 4e-13 Eh(鍵は `reference: uhf`)。
- **基準 B1**(`validation/replay.py B1 --src /home/user/hfauto_r9/VAL9/R2`、`/home/user/hfauto_r10/B1`、比較 `/home/user/hfauto_r10/W2/b1_analysis.txt`): 40 case(`hono_walltime` は superseded)。sp より前の新しいジョブは S5 の paths の 12 件だけ(χ で split2 の鞍点を受理した。事前に挙げた変化)。DFT の SP 121 件はすべて射影の guess で新しい鍵。結論(outcome・順位・blocker・notes)の変化 0、δG_eff の差は最大 2.4e-3 kcal/mol。手法パネルの変化は事前に挙げたものだけ: S5 の BS 錯体の M06-2X ΔE_rxn −83.15 → −36.00(`spin_contaminated:dE_rxn` が付いた。以前は層の SP が ⟨S²⟩ 0.759 の別の解で spin_ok を通っていた)、CH3O の CCSD(T) ΔE‡ 33.37 → 33.47・ΔE_rxn −7.80 → −7.81(UHF)、h_c2h4 の M06-2X の行が加わった(錯体の層の SP が (2)(3) で収束し、sp stage は 2,019 → 311 s)。`check.py B1` は s5 以外 PASS(DEVIATION は h_c2h4 の 1 件で W1 と同じ)。s5 は R2 の経路の再生なので split2 は unresolved のまま: 1 本目の鞍点の QRC は両側が同じ basin(`same_state`)、2 本目は側が鞍点に落ちた(`unassigned_side`)。
  - 事前の期待との違い: (a) h_c2h4 で「dG_assoc が出る」としたが、会合の ΔG_assoc は設計上 ΔG_rxn なので出ない。変わったのは錯体の G が `energy_layer_missing` でなくなったこと。(b) S5 の QRC の側 qrc1_1(CH3OO と同じ状態ラベルで、PBE0 では CH3OO の極小より 22.5 kcal/mol 上)の層の SP は、射影で freq の解に留まり atomic guess の解より 19.0 kcal/mol 上になった(⟨S²⟩ 0.754 と 0.756)。同じ多重度の別の空間状態(CH3OO• の Ã 状態は約 0.9 eV 上で、それと推測する)は ⟨S²⟩ では区別できない。この極小は状態の G の最小ではないので順位には効かない。
- **S5 の新しい run**(`/home/user/hfauto_r10/W2/s5`、1 回): 会合は `scan:barrierless` 0.0、split1 elementary 46.85、**split2 は elementary_step(δG_eff 0.13、M06-2X ΔE‡ 1.93、PBE0 0.81)**で、O–O 開裂の対に答えた(1 本目の QRC は same_state、string の鞍点の QRC が 2 倍の振幅で両端を結んだ)。BS 錯体のパネルは M06-2X −36.00 で `spin_contaminated`。validation/cases.yaml の split2 を unresolved から elementary 0.13 に変えた(理由はこの χ の変更。境界の系なので VAL10 で 2 回流して確かめる)。golden の G37(cases の比較のテストが読む S5 の manifest)はこの run から作り直した。
- **統合で見つけた誤り**: guess を渡したジョブの deck は `vectors input guess.movecs` だけで出力を名指ししておらず、NWChem は収束した vectors を guess.movecs に上書きして job.movecs を残さなかった。そのジョブを親にする次のジョブは guess を見つけられず、黙って atomic guess から始まっていた(鍵には scf_guess が入らない)。最初の S5 の新しい run(`W2/s5_before_vectors_fix`)では会合のスキャンの 1 点おきに guess が落ち、2.50 Å と 2.29 Å の点が二重項の枝(⟨S²⟩ 0.761、BS の枝より 22 kcal/mol 上)に跳んで会合が `unresolved_within_budget` になった。`output job.movecs` を書くように直し、偽の NWChem も NWChem と同じ場所に vectors を書くようにして回帰テストを足した。B1 は影響を受けない(B1 の SP の親は R2 の freq で、SP は末端。鍵も同じ)。直した後の S5 ではスキャンの 7 点のうち最初以外の 6 点が guess を持ち、⟨S²⟩ は 1.72 → 0.76 と単調で R2 と同じ。

### 11.3 W3(段階 2: 反応が読む点、点の鎖と共通のゼロ、捕獲律速、`profile.judge`、組成の seed)

- **B1 の QM なしの再計算**(`/home/user/hfauto_r10/W3-1/replay_points.py`、`replay/replay.json`。対照として HEAD のコードで同じ計算をすると B1 の ReactionThermo を差 0 で再現する: `control/`): 単分子の素過程はすべて 1e-9 以内で不変。S6 の引き抜きは 5.93 → 10.37(幅 10.04〜10.68、分離した OH + CH4 がゼロ)。SN2 14.53 と S10 36.63 は不変(ゼロは錯体)。分割の子の ΔG_assoc と分離基準の ΔG‡(S6 split1 の −10.61 / 10.37、S5 の子)は消えた。S5 の会合、h_c2h4、bh3_nh3、h2o_hf_inversion、water_dimer_as は δG_eff 0.0 → なし(捕獲律速の表へ)。sp の対象は B1 の計算済みの対象と一致した(108、欠け 0、余り 0)。手法パネルは行を失わない(S5 split2 の dE_act のセルだけが空になる。unresolved の鞍点は反応の TS ではない)。`energy_layer_missing` の件数は不変で、値のある反応はどれも読まない。
- **B1 の再生**(`validation/replay.py W3/replay --src /home/user/hfauto_r10/B1`、設定・スキャン・組成の変わらない 31 case): すべての stage で新しいジョブ 0(SCREEN・string の中点の SP も同じ鍵)。記録の差は意図したものだけ: h2o_hf_inversion と water_dimer_as の δG_eff 0.0 → null と順位の削除、SN2(s1、s20)の `reference` complex、種の注記 `thermo_unavailable` の削除(water_same_basin、s3)。ほかは和の順序による 1e-9 未満の差。`check.py` は 31/31 PASS(`W3/replay.check.txt`)。
- **新しい run**(`/home/user/hfauto_r10/W3/fresh`、各 1 回。`check.py` は 8 PASS、s5 FAIL: `W3/fresh.check.txt`):
  - h3_doublet: 鎖に分離した H + H2 の点(min_h、min_h2)があり、ゼロは separated。δG_eff 15.55(錯体ゼロでは 12.71、ΔG_assoc +2.83)。cases.yaml を 15.55 にした(理由: ゼロの定義の変更。BH76 の行 47/48 と同じゼロ)。
  - S6: 引き抜き 10.39(幅 10.04〜10.68、separated、C1 の TS で ΔG‡ 5.95)。生成物側の分離した CH3 + H2O の点(min_ch3、min_h2o)がある。宣言した ch3_h2o の組成は CREST rc 0 で入力状態の候補 6 個を出し、そのうち 1 つが VAL9 の CH3···H2O より G で 0.60 kcal/mol 低い。そのため逆向きの split1 は 20.98 → 21.60 になった(状態の G は最低の配座)。split2 44.83 は不変。cases.yaml の引き抜きを [9.95, 10.45](Cs なら 9.98)、split1 を [21.15, 21.65] にした。2 回目の run は VAL10。
  - h_c2h4、bh3_nh3: 会合のスキャンに中点が 1 回足され(`scan_midpoints:barrierless`)、barrierless のまま。ΔG_rxn −30.65 と −17.83 は B1 と同じ。h_c2h4 の錯体の層の SP はこの run では救済の (2)(3) でも収束せず(`scf_unavailable:scf`)、補助の点なので ΔG_assoc が空になるだけで blocker も rc も出ない(B1 の再計算では R2 の freq から +4.21)。
  - S5: 会合は `scan_mid:ap:4` と `scan_midpoints:barrierless` で barrierless、ΔG_rxn −23.48。split1 46.85。**split2 は unresolved_within_budget**(未達): W2 と同じ鞍点(E −189.8523359 Eh)の QRC の側 1 が鞍点に落ちた(`unassigned_side`)。W2 の run では同じ鞍点の QRC が same_state になり、string の種から結ばれた。QRC は W3 で変えていないので並列の NWChem の非決定性と判断し、期待(elementary 0.13)は変えない。W5-1 の下りと VAL10 の 2 回の run で扱う。
  - S10 36.63(complex)、S19 same_basin、S7・S8(explore まで)は B1 と同じ。近道の case(S5・S6 の親、S10)の記録の判定は None になった。
- **組成**(W3-3。B1・R2 の run と比べた): CREST の rc は B1・R2 と同じ(S5 −11、S6 1、ほかは 0、fecl3_ch4 も 0)。すべての組成に入力状態の構造がある(S5・S6 は seed 6 個ずつ、S7 6、S8 1、S10 1 + 1、S19 3(B1 は 4)、FeCl3·CH4 6、新しい CH3·H2O 6)。(組成, 状態) ごとの DFT の最低エネルギーはすべて B1 + 5e-5 Eh 以下(最大の差は S19 の TMA で +5.5e-7)。DFT の basin の集合は S6 の新しい組成の分を除いて同じ。
- **プローブ**(QM は flock の下): 接触の seed 1 つは S5(scf_not_converged)と S6(geometry_maxiter)の GFN2 の screen で失敗し、両系の反応がすべて消える(`/home/user/hfauto_r10/W3/w3_3_crest_probe/screen_check.log`)。円錐を外した seed は S5・S6・S8・S9 で R2 と同一(差 0)、S7・S10・S19 は新しく、CREST の GFN2 の状態ごとの最低は R2 と 1e-7 Eh 以内(S19 の C3H9N+FH+FH だけ +3.1e-5 Eh、その状態の候補 1 個)。pysisyphus の IDPP は平面の cis/trans HONO で経路を面内に通し(H–O–N 177.8°、二面角が 0 → 180° に跳ぶ。自前は最大 107.8° で単調)、NH3 の最小距離は 0.994 Å(自前 0.997 Å)、作業ディレクトリにログを書く(`/home/user/hfauto_r10/W3/idpp_probe`)。採らなかった。
- **golden**: G36(S6)と G39(h_c2h4)を W3 の新しい run から、G37(S5)を W2 の run を W3 のコードで sp から取り直した `/home/user/hfauto_r10/W3/s5_w2_from_sp`(sp の新しいジョブ 0、split2 elementary 0.13)から作り直した。G38 は VAL9/R2 のまま。

### 11.4 W4(段階 3: 正準のラベル、結合グラフの編集の列挙、辺と世代、DFT//xTB の入口)

証拠は `/home/user/hfauto_r10/W4/`(`w4_1/`、`budget/`、`w42b/`、`w43_rejudge/`、`fresh/`、`fresh.check.txt`、`state_cover.json`)と `/home/user/hfauto_r10/probe/P0f/`。

- **ラベル**(W4-1、QM なし): 正規に近いグラフの構造異性体 3 組(decalin/bicyclopentyl ほか)は旧い WL では同じラベル、新しいラベルでは別。decalin の 100 回の並べ替えでラベルとクラスは不変。B1(40 run、極小 99、species 151)、VAL9/R2、W3/fresh で旧と新の状態の分割は同一(全単射)。次点の規則を消しても B1・R2・W3 の合流は 1 つも変わらない。RDKit は 2026.03.3。
- **P0f**(QM なし、宣言どおり): 不合格。電荷分離の制限はギ酸エチルの 1,5-H 移動の生成物を失い((i))、ibuprofen は 5,200 類で freq 1 本を超えた((iii))。並べ替えの不変性((ii))は 23 系で合格。宣言どおり P0c の規則を使う。repo の列挙器は 23 系すべてで P0c と同じ類の集合(ibuprofen 10,708 類)。
- **予算**(QM なし、偽の NT2): ibuprofen の 1 状態で 3000 件を予算の順に試し、f2b2 の 7,708 類を `not_attempted` にした(試した最大の結合変化 4 = 切った最小 4)。report は `not_attempted` があれば budget-limited と出す。
- **DFT//xTB の再判定**(W4-3、R2 の 74 試行、`w43_rejudge/`): 32 本の辺に SP 60 点。S5 e4ce72fe・S10 7334d674(32.9)は入り、S6 eddde71 は 41.4 で窓の外、316385 が 35.7 で入って同じ CH3OH + H の状態に届く。
- **新しい run**(`fresh/`、各 1 回。統合で直した誤りの後に該当する段から流し直した): R2 の生成物の状態(端の状態ラベル)はどの系でも W4 の集合に含まれる(S5 2 → 35、S6 2 → 17、S7 0 → 10、S8 0 → 436、S10 1 → 204、S19 0 → 161。`state_cover.json`)。`check.py` は S6・S7・S8・S10・S19 が PASS、S5 が FAIL。
  - 予算: S8・S10・S19 は 3000 件に届いた(R2 は 19〜74 件)。`not_attempted` は 527・178・120(うち列挙しないまま残った状態 430・162・112)。生成物の辺は 1027・850・975 で、GFN2 の ΔE‡ が 40 を超えるものが 8 割。S5・S6・S7 と fluoroethane は予算の内(961・267・107・148 件)。受け入れ条件「R2 の系は予算に届かない」は S8・S10・S19 で満たさない。
  - DFT の入口: S10 は 850 本のうち 1 本(32.95、R2 の 7334d674 と同じリレー)だけが入り、縮退転位 36.63 は R2 と同じ。S19 は 975 本のうち 1 本(9.10)で same_basin は R2 と同じ。S6 は 67 本のうち 14 本(引き抜きの辺 0.00、CH3 + H2O → CH4O + H 35.75 など。eddde71 に当たる辺は 41.41 で窓の外)。結果は引き抜き 10.40(帯の中)、CH3 + H2O → CH4O + H elementary 54.62、CH4O + H → CH3O + H2 elementary 13.10、CH2 + H + H2O → CH4O + H barrierless で、cases.yaml を変えた(変更 50)。
  - S5(宣言なし): CH3 + O2 → CH3OO(TS のない辺、第 1 世代)、CH3OO → CH2OOH(第 1 世代、GFN2 ΔE‡ 36.5)は辺として出た。CH2OOH → CH2O + OH の β 開裂は直接の辺としては出ず、CH2O + OH の状態には第 3 世代で CH2O2 + H からの TS のない辺で届いた。入口では未緩和の seed が全体より高く、220 本すべてが高さ 0 になり、順は道のエネルギーの幅で決まった(CH3OO の井戸から出る幅は 51.3)。会合 barrierless、split1 46.85 は VAL9 と同じ。split2 は unresolved(W3 と同じ境界の揺れ。期待は変えない)。入口が入れた辺から CHO + H2O → CH3O2 elementary 38.20・42.48 などの行が増えた。S5 の期待は 2 回流す VAL10 で決める。
  - fluoroethane: 4 中心の 1,2-HF 脱離の TS は第 1 世代で出た(f1b2、GFN2 ΔE‡ 72.2)が、DFT//xTB の高さは 69.5 で窓 40 の外。受け入れ条件「DFT の case が elementary」は満たさない(事前に宣言した衝突。窓は変えない)。
- **統合で直した誤り**(回帰テストあり。[roadmap.md](roadmap.md) 変更 47〜49): (1) S6 の 1 つの状態に原子順の違う構造が入り、別の元素の並びで読んでいた(276 件中 39 組の trial_id が重複)。(2) 入口の到達を species で見ていて、別の配座から出る辺が unreached になった(fluoroethane 12/51)。(3) 入口の同じ高さの順が id の順で、S5 の会合の辺が `over_cap` になった。(4) 入口がすべての辺で SP をとった(S10 で約 2,500 点になる)。(5) 予算の最後の世代で全状態を列挙していた。(6) TS の重複判定が全 TS との RMSD で 2 乗に増えた。(4)〜(6) を直した後、S10・S19 の explore はキャッシュから 42・61 s(以前は 93・103 分)で、辺は同じだった。
- **efficiency の所見**: 世代を重ねる探索は、中くらいの系(16〜17 原子)で予算をすべて使い、辺の大部分は窓のはるか外にある(S10: 849/850 が窓の外)。NT2 の 1 試行は約 10 s・4 並列で、3000 件に約 2 時間かかる。
