# 実計算による検証の記録

round 9(r9-PRE〜S6 と FINAL。結果報告 [2026-10-02_round9_result.md](reviews/2026-10-02_round9_result.md))の全体の再検証 VAL9 の記録である。期待値は、受け入れた最新の再生(CUR = `/home/user/hfauto_r9/S6/replay`、今のコードの結論)と round 9 の検証 run の値で、VAL7(round 7 の全体の再検証)から変えた期待値は §4 に理由とともに挙げる。wave ごとの詳しい記録は `git show c6eff84:docs/validation.md`(§10〜§25)、round 7 は `git show 1bd2ff5:docs/validation.md` にある。決定表の行番号は今の 14 行の表([design.md](design.md) §6.2)のもの。

## 1. 条件

- コード: FINAL のコミット(`c6eff84` に FINAL の作業ツリーを足したもの。run の `git describe` は `c6eff84-dirty`)。2026-10-02 01:25〜07:37。
- 環境: WSL2 Ubuntu 24.04(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 2,000 MB、常に `--bind-to none`。xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0、pymsym 0.3.5)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(I と Te は def2-ECP)、順位のエネルギー層は M06-2X-D3(0)/def2-TZVPD、低レベルは GFN2-xTB。298.15 K・1 atm、順位の量は δG_eff(kcal/mol)。
- run: `/home/user/hfauto_r9/VAL9/<run>` に新しい run dir で 1 本ずつ(`/home/user/hfauto_r9/runcase.sh`: 共有のロック、`/usr/bin/time -v`、CAP は SIGTERM で送り 60 s 後に SIGKILL、孤児の検査、status と report)。
- 対象(41 run と smoke): VAL7 の 35 run(noscreen とパネルの pipeline は `/home/user/hfauto_r9/pipelines/` の写しで、既定と同じ M06-2X の層を持つ。S20 の SN2 の 400 K・1 M の report は、その層の上の検証用の thermo で出す)と、round 9 で足した 6 run(`sn2_cl_d3h`、`malonaldehyde_c2v_undeclared`(`/home/user/hfauto_r9/inputs/malonaldehyde_c2v_undeclared_ab.yaml`)、`h_c2h4`、`bh3_nh3`、ニトロメタン(`/home/user/hfauto_r9/inputs/nitromethane.yaml`、pipeline `g3_soft.yaml`)、`hono_panel`(HONO に CCSD(T) のパネル))。手順は `/home/user/hfauto_r9/val9_chain.sh`。
- 比べ方: 結論の水準で比べる(outcome、生成物の状態の集合、各段の区分、δG_eff が感度の幅の中か)。並列の NWChem と CREST は run ごとに軌跡が分かれる(design.md §10)ので、id ではなく構造の対応で組む。期待の分岐に入らなかった run は未達として記録し、期待値・しきい値・予算は動かさない(分析 §4.4)。
- 時間: 一意のジョブ鍵ごとの core 秒(duration × `-np`、`/home/user/hfauto_r9/tools/jobtime.py`)と歩数。壁時計は参考値。
- 自己一致: VAL9 を再生ゲート(§7)の strict で再生する。

## 2. 実計算で通した範囲

| 項目 | DFT(PBE0/SVPD)まで通したもの | 低レベルだけ |
|---|---|---|
| 元素 | H、B、C、N、O、F、Cl、I(ECP)、Te(ECP、H2Te) | S(SO2·NMe3)、Fe(FeCl3·CH4) |
| 電荷 | 0、−1(Cl⁻・I⁻ の SN2) | — |
| スピン | 一重項、二重項(CH3O•、H + H2、OH···CH4、H + C2H4)、三重項(O2)、低スピン結合の二重項(CH3• + O2) | 六重項(FeCl3·CH4) |
| 反応型 | 1,2-H 移動(HCN、CH3O•、HCOH → H2CO)、1,3-H 移動、ねじれ(HONO、シュウ酸の 2 段)、縮退転位(NH3、DME、H2O2、PT、恒等 SN2、H 交換、(HF)₂)、未宣言の縮退転位(D3h の恒等 SN2 の鞍点だけから)、H 引き抜き(OH + CH4)、障壁なし(H2O·HF の反転、水二量体の受容体交換)、分離した単量体からの会合(CH3 + O2 → CH3OO•、H + C2H4 → C2H5•、BH3 + NH3 → H3B–NH3。いずれも PBE0 で障壁なし) | ハロゲン結合、配位付加体(SO2·NMe3) |

通していないもの: 陽イオン、−2 以下の陰イオン、遷移金属の DFT、開殻一重項、溶媒(design.md §10)。

## 3. VAL9 の結果

判定は、期待と同じなら「合格」、違えば理由を書いて「未達」とする。CUR との結論の比較は `/home/user/hfauto_r9/FINAL/concl.py`(artifact の種類ごとの数、反応の記録(組成、source、outcome、理由、虚振動数)、順位の行)で、結果は `FINAL/concl_val9_vs_cur.txt`。CUR を持つ 36 run のうち 32 本は同一(δG_eff の差 ≤ 0.0004)、差は S6(下)、S19 の CREST の配座の数、パネル 2 本の CUR に残った古い `panel_thermo` の stage(round 7 の写しの残りで、今のパネルには無い)だけである。

### 3.1 known_endpoints と回帰

期待の δG_eff は CUR の値で、差が 0.01 kcal/mol 以内なら同じとみなす(軌跡の分かれる平らな PES の S5・S6 は §3.2 の幅)。

| 系 | 期待(outcome、δG_eff) | VAL7 | VAL9 | 判定 |
|---|---|---|---|---|
| R1 `hcn` | elementary 41.57 | 42.18 | elementary 41.57 | 合格 |
| R2 `hono` | elementary 9.88 | 11.82 | elementary 9.88 | 合格 |
| R3 `nh3_inversion` | degenerate 4.38 | 3.84 | degenerate 4.38 | 合格 |
| R4 `water_same_basin` | same_basin | same_basin | same_basin | 合格 |
| `ch3n_spin` | 両端を `spin_crossing_reaction_unsupported` で拒否、ジョブ 0 | 拒否(旧い文言) | `spin_crossing_reaction_unsupported` で拒否、ジョブ 0(rc 1) | 合格 |
| `formaldehyde` | HCOH → H2CO elementary 31.29、逆向きは out_of_window | 29.02、out_of_window | 31.29、out_of_window | 合格 |
| `h2o_hf_inversion` | barrierless 0.0(`screen_midpoints:barrierless`) | 同じ | barrierless 0.0(`screen_midpoints:barrierless`) | 合格 |
| `water_dimer_as` | barrierless 0.0 | 同じ | barrierless 0.0 | 合格 |
| `nh3_planar_seed` | `follow1:ts_candidate`、degenerate 4.38、paths の新しいジョブ 0 | 3.84 | `follow1:ts_candidate`、4.38、paths の misses 0 | 合格 |
| `dme_c2v_seed` | `follow1:one_side`、degenerate 2.13 | 1.96 | `follow1:one_side`、2.13 | 合格 |
| `hcn_noscreen` | FIND_PATH の `path_hei` から R1 と同じ TS、41.57 | 42.18 | `dft_path` → 41.575(CUR と 0.0004 差) | 合格 |
| `hf_dimer_swap_noscreen` | degenerate 0.80 | 1.39 | 0.80 | 合格 |
| `hono_walltime` | 行 7(`walltime`)、report が出て孤児 0 | 同じ | 行 7(`walltime`)、unresolved_within_budget、report あり、孤児 0 | 合格 |
| S1 `s1_sn2_cl` | degenerate 14.53 | 11.08 | 14.53 | 合格 |
| S2 `s2_sn2_i` | degenerate 9.61 | 6.88 | 9.61 | 合格 |
| S3 `s3_h2te` | Te の ECP で DFT 極小(反応なし) | 同じ | 同じ | 合格 |
| S4 `s4_ch3o_doublet` | elementary 30.26 | 30.61 | 30.26 | 合格 |
| S12 `s12_dme_rough` | degenerate 2.13 | 1.96 | 2.13 | 合格 |
| S13 `s13_h2o2_gauche` | degenerate 1.29 | 1.18 | 1.29 | 合格 |
| S14 `s14_malonaldehyde` | degenerate 0.93 | 0.0(`submerged_barrier`) | 0.93 | 合格 |
| S15 `s15_hono_hno2` | elementary 56.83 | 52.02 | 56.83 | 合格 |
| S16 `s16_acac` | degenerate 0.0(`submerged_barrier`) | 同じ | 0.0(`submerged_barrier`) | 合格 |
| S17 `s17_hf_dimer_swap` | degenerate 0.81 | 1.39 | 0.81 | 合格 |
| S18 `s18_h3_doublet` | degenerate 12.71 | 4.09 | 12.71(CUR との差 0.000) | 合格 |
| `oxalic_two_step` | multi_step(GEN-05)。split1 11.81(`ts_calc` から)、split2 12.71 | 12.78、13.76 | multi_step、11.81、12.71 | 合格 |
| `malonaldehyde_c2v_undeclared` | 宣言なし。C2v の PT の鞍点だけから `follow1:ts_candidate` → degenerate(M1 の実 run)、`qrc1:minus_is_image` | —(round 9 で追加) | degenerate 0.93(S14 と同じ)、`qrc1:minus_is_image` | 合格 |
| `sn2_cl_d3h` | 宣言なし。`follow1:ts_candidate` → `rxn_mode_follow_*` degenerate 14.53、`qrc1:minus_is_image`、paths の新しいジョブ 0 | —(round 9 で追加) | degenerate 14.53、`qrc1:minus_is_image`、paths の misses 0 | 合格 |
| `h_c2h4` | 会合のスキャンで barrierless 0.0(PBE0 の山 +0.60 < 1.0) | — | barrierless 0.0(`scan:barrierless`)。錯体の M06-2X SP が `scf_not_converged`(atomic guess と cgmin の救済の 2 試行、574 s)で rc 1。錯体は会合の順位の点でなく、順位は変わらない(§10、結果報告の残る課題 6) | 合格 |
| `bh3_nh3` | 錯体は付加体に落ち、その入力構造で会合。barrierless 0.0 | — | barrierless 0.0(`scan:barrierless`、`known:basin_adduct`) | 合格 |
| ニトロメタン | 重なり形が soft(停留)、1 回押して `soft:resolved`、status は minimum | — | `soft:resolved`、minimum(ν1 +36.2) | 合格 |

### 3.2 discover

| 系 | 期待 | VAL7 | VAL9 | 判定 |
|---|---|---|---|---|
| S5 `s5_ch3_o2` | 会合 CH3 + O2 → CH3OO•: `scan:barrierless`(AP、注記 `scan:ap:4`)、barrierless 0.0、順位あり。CH3OO• → CH2OOH• → CH2O + OH: multi_step、split1 elementary 46.85、split2 は O–O 開裂の TS を期待しない(CUR は χ 0.271・0.287 の鞍点を棄却して unresolved)。split2 が unresolved なら親と split2 は `thermo_unavailable` | 結論なし | 会合 barrierless 0.0(`scan:ap:4`)。multi_step、split1 46.85、split2 は χ 0.273・0.275 の鞍点を棄却して unresolved、親と split2 は `thermo_unavailable`。CUR と結論・値が同一 | 合格 |
| S6 `s6_oh_ch4` | OH + CH4 → CH3 + H2O: elementary、一次の TS(ν1 約 −490i、ν2 は −50i より上、χ ≥ 0.5)、δG_eff 5.5〜6.0(TS の点群が C1 なら Cs より RT ln 2 低い)、分離基準の ΔE‡(M06-2X)5〜7。CH3 + H2O ↔ CH3OH + H: multi_step(CH4 + HO を経る)、split1 は引き抜きの結論を写す(`same_as`)、split2 elementary 約 44.8 | 引き抜きは unresolved(−73.0i の回転子の鞍点) | elementary 5.54(ν1 −485.3i、ν2 +48.9、χ 0.674、分離基準の ΔE‡ 5.69、TS は C1)。multi_step、split1 `same_as`(20.60)、split2 44.82。CUR(5.95)との差 0.41 は TS の点群(CUR は Cs)の RT ln 2 | 合格 |
| S7 `s7_nh3_icl`(`--to explore`) | 生成物 0 | 同じ | 生成物 0(trial 2) | 合格 |
| S8 `s8_so2_nme3`(`--to explore`) | 付加体が screen の極小に残る。trial 19、生成物 0 | 同じ | trial 19、生成物 0 | 合格 |
| S9 `s9_undeclared` / `s9_conformers` | `declare_multiplicity` / 六重項で CREST rc 0 | 同じ | `declare_multiplicity`(rc 1)/ CREST rc 0 | 合格 |
| S10 `s10_amine_pilot2` | 二重 H 交換 degenerate 36.64。C3H9N と C3H10FN は `not_reacting` | 32.29 | 36.64、`not_reacting` | 合格 |
| S19 `s19_tma_hf2` | same_basin。R6 の seed は作らない(ラベルの違いが帯の中)。順位の点がないので層の SP も ΔG_assoc もない | same_basin、seed は `collapsed_at_dft_from_seed`、ΔG_assoc −11.53 | same_basin、seed なし、ΔG_assoc なし(CREST の配座が 1 つ多い 6 個、結論は同じ) | 合格 |
| `s5_undeclared`(`--to screen`) | 組成 CH3·O2 の多重度が 1 つに決まらず(二重項か四重項)、組成だけが `declare_multiplicity` で止まり、単量体は screen まで進む(rc 1) | 同じ | 同じ(rc 1) | 合格 |

### 3.3 手法パネルと精度

- CCSD(T)/def2-TZVPD との差(5 反応、[design.md](design.md) §7.4 の表。round 9 の実測): 障壁は M06-2X で平均 0.83、最大 1.79 kcal/mol(PBE0/SVPD は 2.55、6.2)。ΔE_rxn は M06-2X で平均 1.17、最大 2.30 で、PBE0/SVPD(1.05、2.25)より良くない。
- composite の δG_eff(同じ run に `energy_method: ccsd-t_def2-tzvpd` の thermo を足したもの)との差: SN2 +0.46(M06-2X 14.53、CCSD(T) 14.07)、HCN −1.79、CH3O• −0.72、HONO +0.09。
- マロンアルデヒドの PT の ΔE‡: CCSD(T) 3.47、M06-2X 2.69、PBE0/SVPD 2.01 kcal/mol(4 ranks × 2,000 MB、1 点約 23 分)。文献の CCSD(T)/CBS 約 4 より 0.5 低い。
- S6 の分離基準の ΔE‡: M06-2X 5.73〜5.81、CCSD(T) 6.96、PBE0/SVPD 0.75〜0.80。
- VAL9 のパネル: `s20_sn2_panel`(400 K・1 M の report と ranking.csv の一致。CUR 14.89)、`s20_hcn_panel` と `s20_ch3o_panel`(CCSD(T)、閉殻と開殻)、`hono_panel`。VAL9: SN2 の 400 K・1 M は 14.89 で report と ranking.csv が一致。ΔE‡ の M06-2X − CCSD(T) は HCN −1.79、SN2 +0.46、HONO −0.01、ΔE_rxn は HCN −2.30、HONO −0.50(M06-2X −0.04、CCSD(T) +0.47 で、cis/trans の順は期待にしない)。いずれも design.md §7.4 の主張(障壁の差 最大 1.8、ΔE_rxn 最大 2.3)の中で、CUR と 1e-4 kcal/mol 以内で同じ。

### 3.4 実行基盤と再生

- smoke(`pytest -m real`、13 件): 13 件合格(66 s)。
- 孤児(nwchem、prterun)が全 run で 0: 合格(41 run)。rc 1 は、拒否・止める設計の 3 本(`ch3n_spin`、`s9_undeclared`、`s5_undeclared`)、低レベルのジョブの失敗を持つ S5・S6(CUR と同じ)、層の SP が失敗した `h_c2h4`。
- 自己一致(`tools/replay.py FINAL --src /home/user/hfauto_r9/VAL9 --strict`、`hono_walltime` と `s9_*` を除く 40 run): 40/40 PASS、misses 0、新しいジョブ 0、記録の差 0(`/home/user/hfauto_r9/FINAL/self_replay.log`)。CUR はこの再生(`FINAL/replay`)を指す。

## 4. 期待値の変更とその理由(VAL7 → VAL9)

| 系 | VAL7 | VAL9 の期待 | 理由 |
|---|---|---|---|
| 順位のあるすべての反応 | PBE0/SVPD の δG_eff | M06-2X 層の δG_eff | 既定の順位を M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD にした(design.md §7.3)。値は PBE0 の非局在化誤差の向き(H 移動・求核置換で上がる)に動く: HCN −0.61、HONO −1.94、NH3 +0.54、SN2 Cl +3.45・I +2.73、S10 +4.35、HONO → HNO2 +4.81、H + H2 +6.35、ねじれと配座は ±1 前後 |
| S18 H + H2 | 4.09 | 12.71 | 上の層(+6.35)に加え、ほぼ直線の H···H2 を点群 C∞v で直線として扱う(+2.27)。2 つの run の差は 2.27 から 0.006 になった |
| S17、`hf_dimer_swap_noscreen` | 1.39 | 0.81、0.80 | TS の点群を 1 つの判定で C2h(σ 2)に決め、層で −0.59 |
| S14 マロンアルデヒド | 0.0(`submerged_barrier`) | 0.93 | M06-2X の障壁(ΔE‡ 2.69)は PBE0(2.01)より高く、ZPE で沈まない |
| R2 HONO | 11.82 | 9.88 | PBE0 のねじれの障壁の過大(+2.0)が消える。δG_eff が ΔG‡(9.79)より大きいのは、層で cis が trans より 0.09 低くなり、同じ状態の G_R がその最小になるからである。cis/trans の順は期待にしない(M06-2X の ΔE_rxn −0.04、CCSD(T) +0.47 で、差が 1 kcal/mol 未満) |
| S5 CH3• + O2 | 結論なし(⟨S²⟩ 1.71 の BS 錯体を端にした問い。−72.0i の鞍点で reassigned、順位なし) | 会合は barrierless 0.0、順位あり | 障壁のない会合では反応物は漸近で極小ではない。分離した単量体から付加体への緩和スキャンで問い(⟨S²⟩ は 1.72 → 0.755 と単調)、低スピン結合のスキャンは AP のエネルギーで判定する。ΔE_e −37.56、ΔH₀ −31.63 は実験の D₀ 約 32 と 0.4 以内 |
| S5 CH3OO• → CH2OOH• | reassigned 32.67 | multi_step(split1 46.85、split2 は結論を期待しない) | 結合が変わる case を化学状態で判定するので、−1866.6i の TS の側(CH2OOH•)は中間体になる。split2 の SCREEN は端の SCF の解から解くと収束し、O–O 開裂の領域の鞍点(χ 0.271〜0.287)は反応モードでないので棄却する。REACTION_MODE_MIN は変えない |
| S6 OH + CH4 | unresolved(−73.0i の回転子の鞍点で試行を使い切る) | elementary、5.5〜6.0 | χ で回転子の鞍点(χ ≤ 0.023)を QRC の前に棄却し、ρ だけを負の曲率にした初期 Hessian で引き抜きの TS に届く(3 種 × rank 2・4 の 6 run すべて)。R6 の seed は xTB Hessian のモデルから直接極小に入るので、仮説の id が変わった。W3 の 1.48 は PBE0 の値 |
| S6 CH3 + H2O ↔ CH3OH + H | reassigned 30.66(W3 は multi_step) | multi_step(CH4 + HO を経る) | 同じ状態の組の仮説を 1 つにまとめ、TS の側を状態で判定する。子 CH3 + H2O → CH4 + HO は引き抜きと同じ問いなので、その結論を写す |
| S19 TMA·(HF)₂ | seed の DFT opt(9,406 core 秒)で `collapsed_at_dft_from_seed`、ΔG_assoc −11.53 | seed を作らない。ΔG_assoc なし | seed と落ちた basin のラベルの違い(N···H)は結合変化の帯(±0.1 Å)の中で、失われた状態ではない。順位の点を持たない反応には層の SP を計算せず、PBE0 で代用しない |
| `h_c2h4`(新) | — | barrierless 0.0 | 計画の期待(約 2 kcal/mol の山から saddle へ)は満たさない。PBE0-D3BJ/def2-SVPD の山は +0.60 kcal/mol(谷からの深さ 0.66)で解像度未満。参照の障壁 1.72 を混成汎関数が過小評価するのは既知で、resolution とスキャンの点は変えていない。層は障壁のない判定をやり直さない |
| `ch3n_spin` | `endpoints differ in charge, multiplicity or atom order` | `spin_crossing_reaction_unsupported` | 多重度だけが違う端点に固有の理由を付けた |

## 5. 到達表(今のコード)

今のコードの記録(VAL9 の 41 run と round 9 の検証 run)の case log と diagnostics から数えた(`/home/user/hfauto_r9/F2/reach.py`)。

| 対象 | 通ったもの | 通っていないもの |
|---|---|---|
| 決定表の行 | 1(宣言した平面 NH3、mode_follow 0 の検証 run)、2(水、S19)、3(formaldehyde)、4(多数)、5(oxalic、S5、S6)、6(SCREEN と中点: H2O·HF、水二量体。会合のスキャン: S5、H + C2H4、BH3 + NH3)、9・10・12(多数)、11(S5 の split2 の `path_intermediate`)、7(VAL9 の `hono_walltime`)、13(noscreen 2 本、S5)、14(S5 の split2) | 8(VAL7 と W5 の S6 は前のコード) |
| 駆動しない子 | `same_as`(S6 の split1)、`split_depth`(`max_split_depth: 0` の oxalic の検証 run) | — |
| saddle の種 | `screen_ts`、`screen_hei`、`discovery_ts`、`path_hei` | `higher_order_retry`(M4 の S6 の実 run 6 本では通ったが、今のコードの S6 は直接一次の TS に届く)、近道の種がすべて失敗した後の NEB |
| 初期 Hessian | xTB(`xtb:mode`、`xtb:rho`)、DFT(noscreen 2 本) | TS freq(M4 の S6 だけ) |
| 検証 | `ts_calc`(nh3_planar_seed、oxalic、S5、S6、`sn2_cl_d3h`)、`ts_rejected:not_reaction_mode`(S5 の split2) | `ts_rejected:higher_order` と `higher_order:not_stationary`(今のコードの記録になし)、`saddle:above_path_bound` |
| QRC | `minus_is_image`、`end<i>_to_new_basin` | 2 倍の振幅の 2 回目の QRC |
| 極小 | `ts_candidate`、`one_side`、`image_of`、`known`、`not_reacting`、`soft:resolved`(S19 の screen、ニトロメタンの DFT) | 非停留の極小の Newton の押し出し |
| 開殻 | 端の SCF の解からのプロファイルの SP、AP(S5)、SCF の救済(harness) | 2 つの guess の分岐、実計算の `scf_branch_jump` |
| outcome | blocked_upstream 以外のすべて(blocked_upstream は行 1 の検証 run) | — |
| その他 | ケースの例外の閉じ込めは単体テストだけ | — |

## 6. 計算量

| 比較 | 結果 |
|---|---|
| VAL7(35 run) | QM 78,752 core 秒(665 ジョブ)、壁時計 20,237 s |
| VAL9(41 run) | QM 84,471 core 秒(812 ジョブ)、壁時計 22,312 s(参考)。VAL7 と対の 35 run は 71,731 core 秒(742 ジョブ、−8.9%)で、うち M06-2X の層が 9,755。層を除くと −21%。減った主なものは S19(−9,658、seed の DFT がない)、S5(−2,520)、S6(−1,819)、増えたのは層の重い acac(+4,558、層 4,848)とシュウ酸(+1,821)。全体の層は SP 82 本・12,940 core 秒(`h_c2h4` の収束しない錯体の 2,299 を含む) |
| HCN のクイックスタート(新しい run、`/home/user/hfauto_r9/F2/hcn_quickstart`) | 壁時計 1:11、QM 242 core 秒(22 ジョブ)。M06-2X の層の SP 3 点で 35 core 秒 |
| エネルギー層(M9 の再生、36 run) | M06-2X の SP 105 本、11,860 core 秒。1 点 10〜60 s が多く、シュウ酸(8 原子)約 140 s、acac(15 原子、337 基底関数)約 600 s |
| R6 の seed(M6) | S19 の seed の DFT opt(9,406 core 秒)がなくなった |
| 会合(M8) | S5 の会合は saddle と string の代わりにスキャンの opt 7 本(687 core 秒) |
| autoz の継続(S1) | S10 の QRC の側の 2 回目の試行が、始点からの 19 歩・42 s から最後の frame からの 1 歩・3.3 s に |
| 残る重複 | DFT の mode-follow は対称な鞍点でも ± 両側を opt と freq にかける(マロンアルデヒドで freq 4 本が 9.1k core 秒の 69%、`sn2_cl_d3h` で像の側が QM 時間の 39%) |

rank 2 の site での `--bind-to none` の修正(M4)より前の rank 2 の時間は、QRC の opt で約 5 倍に膨らんでいるので比べない。

## 7. 再生ゲート

振る舞いを保つ整理の合格基準(分析 §3 X8)。記録済みの run を今の作業ツリーで JobStore から再生し、新しい QM ジョブが 0 件で、結論の記録が変わらないことを確かめる。

- 道具(repo の外、`/home/user/hfauto_r9/`): `tools/replay.py`(本体)、`tools/runs.json`(再生する run と pipeline の連鎖)、`runcase.sh`(1 本ずつの実行)、`tools/jobtime.py`。
- 手順: source を `<wave>/replay/<name>` に複製し、pipeline の最初の stage から再ステージする。パネルの run は続けて `RESUME=1` で panel pipeline を再開する。
- 比べるもの: 再生した stage ごとの `run_state.json`(status、hits、misses、failures_by_kind)、manifest の結論の記録(species、minimum、discovery、reaction、species_thermo、reaction_thermo、report)を欄の単位で、複製の中で走ったジョブ。reason の文字列、FileRef、created_at と壁時計は比べない。
- 合格(`--strict`): 再生したすべての stage で misses 0、source でジョブを走らせた stage はすべて再生して hits > 0、status と記録が同一。`--from <run>=paths` のような浅い再ステージは不合格になる。振る舞いを変える wave は `--allow-new` で挙げた run の差を化学的な理由で説明する。
- CUR: 受け入れた最新の再生への symlink(`/home/user/hfauto_r9/CUR/<name>`)。次の wave は `--src CUR` で再生する。
- 除外: `hono_walltime`(Deadline は実時間なので、キャッシュの hit で行 7 の発火が変わる)、`s9_*`(遷移金属の DFT は未検証)。
- 限界: 保証するのは同じ鍵に同じ結果が返る範囲だけで、新しいジョブが 1 本でも出ると、その下流は run 間の比較(結論の水準)になる。古い JobStore の Evidence には勾配がなく、Failure にはエネルギーがない。再生ではそれぞれ「判定しない」「上限をかけない」として扱う。
- round 9 の結果: 道具の確認(陰性対照 3 本がすべて FAIL)、HEAD の自己一致(36 本、新しいジョブ 0)、振る舞いを保つ wave(M5〜M7、S6 など)の strict PASS。各 wave の数は結果報告 §1。

## 8. プローブと実 run で決めたこと

| 決めたこと | 根拠 |
|---|---|
| 反応モード性の下限 `REACTION_MODE_MIN` = 0.3 | M3 の時点の CUR の χ: 接続した TS 28 件は 0.571〜1.000(重原子の結合変化を含むものは ≥ 0.776)、接続しない一次の鞍点は ≤ 0.023(S6 の回転子 −57.3〜−77.1i、S5 の −72i)。最も近い棄却は S5 の CH2OOH• の領域の 0.221〜0.287 で、0.3 との差は小さい(`/home/user/hfauto_r9/M3/chi/`、`tools/chi.py`) |
| ρ だけを負の曲率にした初期 Hessian | VAL7 の S6 の種のモデルは負の固有値 −0.0028、ρ との cos 0.345。新しいモデルは −0.2187、cos 1.000。S6 を 3 種 × rank 2・4 で流した 6 run がすべて引き抜きの TS に届き、回転子の鞍点には一度も着かなかった(`M4/real/`) |
| 原子の対応(basin の構造を WL クラスで端点の添字に並べる) | S6 の xTB の IRC の端 bfac は PBE0 の basin から 0.957 Å 離れていて、元素だけで並べると結合変化が 8 本(IRC 自身は 4 本)になった。端の監査の不一致 6/17 → 0(`M1/audit/`) |
| 結合変化の帯 `RESOLVED_A` = 0.1 Å | S19 の N···H は GFN2 と PBE0 で 0.15 Å 違う。帯でラベルだけが割れる組は 11 組、本物の結合変化は帯の外(最小は acac の O–H の +0.175 Å)(`S3/S3b/`) |
| 停留性の基準 5e-5 Eh(basin のエネルギー幅と同じ) | VAL7 の DFT 極小 64 点の ΔE_N は ≤ 1.2e-5、W3 の肩(−6.5i)だけが 2.66e-4(`tools/stationarity.py`) |
| 点群の受理を basin の基準にする | 1e-3 Å の揺らぎ 20 回で 110 subject の点群・σ・m・直線性が不変、既知の点群 6 つが一致、S18 の 2 run の差 2.27 → 0.006 kcal/mol(`tools/symcheck.py`、`M7/`) |
| 会合のスキャン(付加体 + 1.5 Å から 8 点) | 1.5 Å は Coulson–Fischer 点(CH3·O2 の C–O で約 2.2 Å)の外。S5 の BS の山 +0.10 は解像度未満で、ΔH₀ は実験の D₀ から 0.4 以内(`M8/real/s5_ch3_o2`) |
| 開殻のプロファイルの SP を端の SCF の解から直接解く | プローブ G2-1v: 端の vectors から直接解いた点は逐次の継続と 5e-8 Eh 以内で一致し、山は +0.80(解像度未満)。atomic guess では +28.9 の偽の山。近い端だけでは中央の点が偽の枝に落ちうるので、両端の ⟨S²⟩ が違えば 2 つの guess の低い方をとる(`probes/G2-1v/`) |
| 枝の跳びの判定に spin_tol を使う | 本物の TS の最大点と隣の ⟨S²⟩ の差は ≤ 0.022、SCF の枝の跳びは 0.77〜0.95(`S4/audit/`) |
| 再開の上限(プロファイルの最大 + resolution) | VAL7 の S5 の再開は最大より 27.2 と 32.0 kcal/mol 上から、収束した再開(W3、W5)は 21〜69 kcal/mol 下から始まっていた(`tools/bound_check.py`) |
| trial の番号付けの不変性 | 開始点 32 + 19 の原子を無作為に 20 回並べ替えても drive がすべて同じ。旧コードは 56 の開始点 × cap で失敗(`S5/int/trialperm_*.txt`)。S8 の新しい run 2 本で trial_id の列と生成物 0 が一致 |
| エネルギー層 M06-2X-D3(0)/def2-TZVPD | 4 反応の ΔE‡ の CCSD(T) との平均絶対誤差 1.03(PBE0/SVPD 2.68、ωB97X-D3 1.56)。TMA·(HF)₂ の 1 点 255 s(ωB97X-D3 1,689 s)。grid fine と xfine の差 ≤ 0.01(r8 のプローブ G6-1) |
| WFT のメモリ(2,000 MB/rank、GA に 70%)と `freeze <n>` | `memory total 1200 mb` ではマロンアルデヒドの (T) が配列を確保できず、新しい配分で 1 点約 23 分で完走。HI は `freeze 4`、I⁻···CH3I は `freeze 9`(`S1/real/`、`S1a/`) |
| autoz の失敗は最後の frame から Cartesian で続ける | S10 の QRC の側: 始点からの 19 歩・42 s → 1 歩・3.3 s、エネルギー差 1.4e-8 Eh、同じ basin(`S1/harness/autoz`) |
| NWChem は常に `--bind-to none` | rank 2 の site で QRC の両側が同じコアを取り合い、1 歩 30 s(rank 4 では 6 s)(`M4/real/`) |
| 一次の鞍点の Hessian は正定値のモデル、高次の鞍点はそのまま | 負のまま渡すと下りが遅い(acac の QRC 34/27 → 15/14 歩)。高次の鞍点からの側はモデルにすると這う(DME C2v で 29 歩 → 207 歩でも未収束)(round 7) |
| 極小の stage は直列 | 4 つの opt を 1 rank ずつにすると TMA 105 → 206 s で 2 コアが遊び、17 原子の seed は約 2 h の見込み(round 7、S19) |
| AFIR を削除 | 41 件すべて陰性で生成物 0(round 7) |
| 採らなかったもの | χ の下限を下げる、2 回目の higher_order_retry、tight の収束、noise_cm1・saddle_cm1 を動かす、ωB97X-D3 の層、CREST の電子温度を上げる、linopt 0(分析と round 7 の理由のとおり) |

## 9. 撤回・訂正した値

- round 7 までに撤回したもの: HONO の ΔG‡ 24.11(ZPE の二重計上)、HCN の ΔG_rxn 18.05(内部フォールバック)、I⁻ + CH3I の ΔG‡ 7.53(対称数)、R2 12.2265(m = 2 を数えていない)、S4 の未宣言のねじれ 3.878。
- W3 の S6 の引き抜き δG_eff 1.48(−481.1i)は、平らな PES の肩(−6.5i、ΔE_N 2.7e-4 Eh)を端にした PBE0 の値で、今のコードの記録にはない。今の値は 5.5〜6.0(M06-2X)。
- VAL7 の S5 の「結論なし(多参照性のある二重項の結合)」: 主因は多参照性ではなく問いの形だった(§4)。
- M4 で記録した S6 の「平らな回転子の領域(ν2 −77i〜+24i)」は、停留していない点の Hessian で測った値だった。停留点では一次の鞍点になる(Newton の歩のプローブ)。
- 分析の「既定の mode_follow でも宣言した鞍点の端点は行 1 に落ちる」は成り立たない: 平面 NH3 は mode-follow で NH3 の basin に加わり、same_basin になる。
- 計画の凍結軌道の例(I は 9、I⁻ + CH3I は 19、CH3O は 1)は NWChem の `freeze atomic` の規則と合わない(I は 4、I⁻···CH3I は 9、CH3O は 2)。
- S17 の M4 での変化(1.39 → 0.98。TS を取り直して点群が Cs と判定された)を当時記録していなかった。M7 で C2h に戻った。

## 10. 実証していないこと

- 決定表の行 8 は今のコードの実計算で通っていない。近道の種がすべて失敗した後の NEB、2 倍の振幅の 2 回目の QRC、ケースの例外の閉じ込めは単体テストだけ。
- 極小の Newton の押し出しと、TS の `higher_order:not_stationary` は単体テストとプローブだけ(実計算に非停留の極小も、勾配を持つ高次の鞍点もなかった)。再開の上限は QM なしの評価だけ(取り直した maxiter の saddle は収束した)。
- 中点の密化が隠れた山を見つけた例はない(対称な経路でしか通っていない)。会合のスキャンの山から REFINE_SADDLE に進む分岐は単体テストだけ(DFT の山が 1 kcal/mol を超える会合が検証セットにない)。
- エネルギー層の SP の SCF: VAL9 の `h_c2h4` で H···C2H4 錯体(二重項)の M06-2X の SP が atomic guess と cgmin の救済の両方で収束しなかった。層の SCF は freq の vectors から始めない(結果報告の残る課題 6)。順位の点ではないので順位は変わらないが、順位の点で起きれば `thermo_unavailable` になる。
- 開殻: 両端の ⟨S²⟩ が違うプロファイル(2 つの guess)、実計算での `scf_branch_jump`、`spin_ok` を通る SCF の救済、四重項の CH3·O2 で AP がかからないこと(dft の段で CAP に達した)。
- TS の停留性: 一次の鞍点の受理には停留性をかけていない。S6 の引き抜きの領域では、ほぼ直線の C–H···O をまたぐ autoz のねじれで内部座標の勾配が縮み、ΔE_N 6.9e-5〜2.2e-4 Eh の点が収束と判定されたことがある。40 cm⁻¹ 未満の大振幅のモードでは二次モデルが成り立たない(acac の TS で 8.6e-3 Eh)。
- 3 次以上の鞍点の `off_saddle`、`max_trials_per_source` の境界で同値の類を落とす規則は単体テストだけ。3 原子以上の直線分子の NT2 の開始構造は番号付けに依る。
- 縮退転位の第 2 候補(ギ酸二量体 D2h)、陽イオン(HOC⁺)、遷移金属・Te を含む系の M06-2X の SP、ECP 原子を含む CCSD(T) の参照との比較。
