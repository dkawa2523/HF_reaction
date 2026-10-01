# 実計算による検証の記録

round 7 の改良(W0〜W5、結果報告 [2026-09-30_round7_result.md](reviews/2026-09-30_round7_result.md))の実計算の記録である。検証セットはレビュー `docs/reviews/2026-09-27_platform_review.md` §6 と、round 7 で足した分岐用の系。round 6 までの記録(VAL を含む)は git の履歴(`git show f2309b9:docs/validation.md`)にある。決定表の行番号は、各節の時点の表のもの(round 7 は 17 行、M5〜S1 は 16 行。今の表は [design.md](design.md) §6.2)。

## 1. 条件

- コード: 各 wave の作業ツリー(W1 = `707e656`、W2 = `26389e3`、W3 = `2eaaf53`、W4 = `c5b515e` でコミットしたもの)。VAL7 は W5 の作業ツリー(W4 に挙動を変えない削除だけを足したもの)で、§6 の全体と round 7 の系を `/home/user/hfauto_r7/VAL7` に流し直した(`/home/user/hfauto_r7/VAL7/chain.sh`、比較は `/home/user/hfauto_r7/VAL7/val7check.py` と `/home/user/hfauto_r7/VAL7/val7walls.py`)。数値の比較先 VAL は round 6 の検証(`/home/user/hfauto_r6/VAL/runs`、コード `d5fd24c`)。
- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB(r9-S1 から 2,000 MB。メモリはジョブ鍵に入らない)、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0、pymsym 0.3.5)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(I は def2-ECP)、低レベルは GFN2-xTB。298.15 K・1 atm、順位の量は δG_eff(kcal/mol)。
- run: `/home/user/hfauto_r7/<wave>/<run>` に新しい run dir で 1 本ずつ直列(`/home/user/hfauto_r7/runcase.sh`。`/usr/bin/time -v`、終了後の孤児プロセス検査、status と report)。記録の型が変わったので round 6 の run は再開せず、VAL との比較は生の JSON と xyz だけを読む。QM 時間は一意のジョブ鍵ごとの秒数(`/home/user/hfauto_r7/tools/jobtime.py`)。GNU timeout は負荷の下で約 5.5% 遅れて発火する。
- 一時 pipeline(`/home/user/hfauto_r7/pipelines/`):
  - `known_endpoints_noscreen`: known_endpoints から reaction-paths の `screen:` を外したもの(行 12 がなく、FIND_PATH が string を走らせる)。
  - `known_endpoints_walltime`: known_endpoints の reaction-paths に `policy: {walltime_h: 0.01}` を足したもの(行 7 の予算の仕組みだけを見る)。
- 合否: 終了コード、outcome、⟨S²⟩、既知の数値(VAL)からの差 0.01 kcal/mol 以内(意図した変化を除く)。新しい系の期待値は、DFT プローブ(`/home/user/hfauto_r7/W1/probe`、pipeline と同じ Level)で状態を確かめてから system の見出しに書いた。期待の分岐に入らなかった run は未達として記録し、期待や閾値は動かしていない。

## 2. 実計算で通した範囲

| 項目 | DFT(PBE0/SVPD)まで通したもの | 低レベルだけ |
|---|---|---|
| 元素 | H、B(M8 の BH3 + NH3)、C、N、O、F、Cl、I(ECP)。Te(ECP、VAL7/s3_h2te) | S(SO2·NMe3)。Fe(FeCl3·CH4、VAL7/s9_conformers) |
| 電荷 | 0、−1(Cl⁻・I⁻ の SN2) | — |
| スピン | 一重項、二重項(CH3O•、H + H2、OH···CH4)、三重項(O2)、低スピン結合の二重項(CH3• + O2、M8) | — |
| 反応型 | 1,2-H 移動(HCN、CH3O•、HCOH → H2CO)、ねじれ(HONO、シュウ酸の 2 段)、1,3-H 移動、縮退転位(NH3、DME、H2O2、PT、恒等 SN2、H 交換、(HF)₂)、H 引き抜き(OH···CH4 → CH3···H2O)、障壁なし(H2O·HF の反転、水二量体の受容体交換)、会合(TMA·(HF)₂)。分離した単量体からの会合(M8: CH3 + O2 → CH3OO•、H + C2H4 → C2H5•、BH3 + NH3 → H3B–NH3。いずれも障壁なし) | ハロゲン結合、配位付加体(SO2·NMe3) |

通していないもの: 陽イオン、−2 以下の陰イオン、遷移金属の DFT、開殻一重項、溶媒。

## 3. 結果

### 3.1 回帰(R1〜R4、known_endpoints)

| # | 系 | VAL | W0〜W2 | W3・W4・VAL7 | 判定 |
|---|---|---|---|---|---|
| R1 | HCN → HNC | 42.1783 | 42.1783(`hcn_noscreen` の string 経由も 42.1786、TS は R1 と 1e-8 Eh 以内) | 42.1784 | 合格 |
| R2 | HONO trans → cis | 11.8158 | 11.8158 | 11.8159 | 合格 |
| R3 | NH3 反転 | 3.8395 | 3.8395(`nh3_planar_seed` 3.8400) | 3.8395 | 合格 |
| R4 | 水 | same_basin | same_basin | same_basin | 合格 |

VAL7 と、行 8 の修正の後の `/home/user/hfauto_r7/W5` の再実行は W4 と 1e-6 以内(W5 の壁時計はホストの別の計算で約 2 倍に伸びたので比べない)。W3 以後の差(≤ 1.2e-4)は freq の SCF を opt の vectors から始める修正(§7)による。

### 3.2 known_endpoints の系(W2、S1・S14・S16 は W4、すべて VAL7 で再現)

VAL7 の δG_eff は round 7 の値と 2e-4 以内で、outcome も同じ。時間は VAL7 の壁時計。

| # | 系 | outcome | VAL → round 7 | 時間 VAL → VAL7 | 判定 |
|---|---|---|---|---|---|
| S1 | Cl⁻ + CH3Cl | degenerate(−393.6i) | 11.077 → 11.0774(ΔG_assoc −5.08、分離基準 +6.00) | 7:30 → 4:18 | 合格 |
| S2 | I⁻ + CH3I | degenerate(−307.8i) | 6.8817 → 6.8816 | 13:19 → 7:05 | 合格 |
| S4 | CH3O• → CH2OH• | elementary(−2015.5i) | 1 位が未宣言のねじれ(3.878)→ 宣言反応 30.612 だけ | 14:30 → 6:27 | 合格(意図した変化) |
| S12 | DME の粗い入力 | degenerate(−226.9i) | 1.964 → 1.9638 | 10:04 → 6:53 | 合格 |
| S13 | H2O2 のねじれ | degenerate(−290.5i) | 1.180 → 1.1802 | 1:38 → 1:05 | 合格 |
| S14 | マロンアルデヒドの PT | degenerate、`submerged_barrier` | 0.0 → 0.0 | 16:52 → 10:10 | 合格 |
| S15 | trans-HONO → HNO2 | elementary(−2088.1i、ΔE‡ 55.3) | 52.016 → 52.016 | 2:46 → 2:15 | 回帰だけ(多段の被覆には数えない) |
| S16 | acac の PT(15 原子) | degenerate、`submerged_barrier` | 0.0 → 0.0 | 1:30:57 → 52:14 | 合格 |
| S17 | (HF)₂ の入れ替え | degenerate(C2h −227.3i、ΔE‡ 1.256) | 1.392 → 1.392 | 1:35 → 1:01 | 回帰だけ(障壁なしの被覆には数えない) |
| S18 | H + H2 | degenerate(−605.7i) | 4.091 → 4.0909 | 0:48 → 0:31 | 合格 |

### 3.3 round 7 で足した系(W1、`oxalic_two_step` は W2)

VAL7 でもすべて同じ outcome で、δG_eff の差は 2e-4 以内(oxalic は 33:21、`hono_walltime` は W1 の `hono_walltime_b` と同じくジョブ 4 本で行 7)。

| 系 | 通す分岐 | 実測 | 判定 |
|---|---|---|---|
| `formaldehyde` | 行 3 と 1,2-H 移動 | 順方向 out_of_window(ΔE_rxn 53.75)。逆方向 elementary(−2084.1i、ΔE‡ 32.96、δG_eff 29.02) | 合格 |
| `h2o_hf_inversion` | 行 6 | `screen_midpoints:barrierless` → barrierless_at_resolution、SP 11 本(9 + 中点 2)、saddle ジョブ 0、δG_eff 0。両端はピラミッド形(DFT プローブで 45.5°) | コードの通過だけ(対称な経路なので、中点の SP が隠れた障壁を見つけることは確かめていない) |
| `water_dimer_as` | 行 6 か 8・11 | 行 6(中点 2 点を加えても barrierless)、δG_eff 0。文献の障壁 0.5〜0.6 は解像度 1.0 未満 | そのまま記録 |
| `oxalic_two_step` | GEN-05 の multi_step | `qrc1:end0_to_new_basin` → multi_step。split1 cTc → tTc は `ts_calc` から検証して新しいジョブ 0、elementary(−666.8i、ΔE‡ 15.08、δG_eff 12.781)。split2 tTc → tTt は elementary(−630.2i、ΔE‡ 13.23、δG_eff 13.763)。W1 は QRC の片側が OH のねじれを 110 歩以上滑って CAP 5,400 s に達し、W2 は 38:16 | 合格 |
| `nh3_planar_seed` | DFT の mode-follow と `ts_calc` | D3h の NH3 が `follow1:ts_candidate` で NH3 の basin に入り、R3 の case は `ts_calc` の validate_ts から始まる(paths の新しいジョブ 0)。3.8400 | 合格 |
| `dme_c2v_seed` | 虚モード 2 本の片側 follow | C2v の DME が `follow1:one_side` で DME の極小へ。宣言反応は 1.9638 | 合格 |
| `hcn_noscreen` | 行 15 → string → `path_hei` | TS は R1 と同じ(1e-8 Eh 以内)、42.1786 | 合格 |
| `hf_dimer_swap_noscreen` | `path_hei` → higher_order_retry | `path_hei` から直接 C2h の TS(−226.9i、ΔE‡ 1.256)、degenerate 1.392 | higher_order_retry はここでは未達(§4 で W3/s6 が通した) |
| `hono_walltime_b` | 行 7 | 最初の decide で `walltime`、report が出て孤児 0 | 予算の仕組みだけ |

宣言反応で行 1 を通す系はない。両端の多重度が違う宣言は structures が拒否する(`/home/user/hfauto_r7/W1/probe/ch3n_spin`: `endpoints differ in charge, multiplicity or atom order`。M8 からは `spin_crossing_reaction_unsupported`、§18)。

### 3.4 discover の系(W3、VAL7)

VAL7 では S7・S8(生成物 0、試行数 2・19)、S10(32.2892)、S19(`collapsed_at_dft_from_seed`、ΔG_assoc −11.5288)が W3・W4 と同じ。S5 と S6 は平らな PES の上で NWChem の軌跡が run ごとに分かれ、次のように変わった。

- S5(`/home/user/hfauto_r7/VAL7/s5_ch3_o2`): CH3···O2 → CH3OO• の case は 2 回目の再開が DFT の Hessian(重なり 0.21)で収束し、2 つの CH3···O2 vdW 極小を結ぶ −72.0i の鞍点で reassigned(`spin_contaminated` で順位なし)。会合の結論がないことは W3 と同じ。
- S6(`/home/user/hfauto_r7/VAL7/s6_oh_ch4`): 同じジョブ鍵の CH3OH + H の opt が W3 では肩(−6.5i の雑音モード、0.10 kcal/mol 上)に止まり、VAL7 では真の極小に入って 2 つの発見が 1 つの basin にまとまった(W3 の multi_step は reassigned 30.659 に、発見の仮説は 1 件減った)。R6 の仮説 OH···CH4 → CH3···H2O は screen の頂点から −73.0i の鞍点に収束したが、QRC の両側が 2 回とも OH···CH4 の basin に戻り(行 8 の初めての到達)、saddle の試行を 1 回残したまま `connection_failed` で終わった。この打ち切りを直し(§7)、同じ run を paths から流し直した `/home/user/hfauto_r7/W5/s6_oh_ch4_rowfix` では FIND_PATH → `path_hei`(maxiter)→ 再開 → −58.6i の鞍点も両側が同じ basin で、試行を使い切って unresolved。W3 の −481.1i の引き抜きの TS(δG_eff 1.48)は VAL7・W5 では再現していない。

| # | 系 | 実測 | 時間 | 判定 |
|---|---|---|---|---|
| S5 | CH3• + O2 | CREST rc −11 で配置の seed を使う。失われた seed の状態を DFT に問うと CH3···O2 の vdW 極小(C···O 2.84 Å)が残るが ⟨S²⟩ 1.71(spin_contaminated)。CH3···O2 → CH3OO• の case は saddle の maxiter 4 回で attempts_exhausted。別の発見の case は reassigned(32.67) | 46:49 | 結論なし(多参照性のある二重項の結合)。M8 で会合は障壁なし・順位あり(§18.4) |
| S6 | OH···CH4 | CREST rc 1。IRC の端の添字の付け直しで未接続 2 → 0。seed の DFT opt は対称な鞍点に止まり mode-follow で 1 つの basin。仮説 OH···CH4 → CH3···H2O は SCREEN single → saddle の maxiter 2 回 → string → `path_hei` で −481.1i、elementary、ΔE‡ 2.40(分離基準 0.79、文献 約 5)、δG_eff 1.48 で 1 位。CH3···H2O → CH3OH + H の発見は multi_step(split1 elementary 48.48、split2 は分割の深さの上限で未解決) | 1:33:13 | W3 は合格(基準どおり)。VAL7・W5 は TS を再現せず(上記) |
| S7 | NH3···ICl(`--to explore`) | 2 attempt、生成物 0(VAL と同じ) | 0:13 | 合格 |
| S8 | SO2·NMe3(`--to explore`) | 付加体(N···S 2.17 Å)が screen の極小に残る(`--notopo` のない旧規則では残らず、最低でも +17.3 の vdW 構造だけ)。19 attempt、生成物 0(VAL と同じ) | 1:10 | 合格 |
| S10 | amine_pilot2 | 二重 H 交換 degenerate 32.2892(W2 32.2887)。C3H9N と C3H10FN は `not_reacting` で DFT に流さない | 3:30 | 合格 |
| S19 | TMA·(HF)₂ | 失われた seed(17 原子)は DFT の 1 歩目で状態を離れ、45 歩・2,033 s で既知の neutral の basin に落ちた(`collapsed_at_dft_from_seed`、freq なし)。宣言反応は same_basin、ΔG_assoc −11.528。付け直した IRC の端 2 件は障壁の上限で棄却(ΔE‡ 115.9 / 139.9) | 1:05:16 | 合格 |
| S4 | CH3O•(`--to explore`) | T1 の 1,2-H 移動が両向きで生成物(GFN2 の ΔE‡ 31.2 / 39.7)、T4 の切断 5 件は `no_nt2_maximum` | 0:05 | 合格 |

### 3.5 手法パネルと条件

- `/home/user/hfauto_r7/W1/s20_sn2_panel`(新しい run に method_panel を追記、report は 400 K・1 M): report.html の δG_eff は全行で ranking.csv と一致(16.09)。298.15 K・1 atm の ωB97X-D3 層は 15.732(VAL と同じ)。
- ΔG_assoc は状態の G を 1 つの定義(spin_contaminated でない最小の G)にしても S19 で −11.5275(VAL −11.5274)。
- VAL7 で流し直した低レベル・パネルの系: S3 H2Te(`Te library def2-ecp`、opt・freq のエネルギーは VAL と 5e-9 Eh 以内)、S9 FeCl3·CH4(宣言がなければ `declare_multiplicity: d-block element(s) Fe`、六重項で CREST rc 0)、CCSD(T) パネル(PBE0/SVPD の ΔE‡ − CCSD(T): HCN −1.19(`s20_hcn_panel`、閉殻)、CH3O• −0.36(`s20_ch3o_panel`、開殻)、どちらも VAL と 1e-3 以内)、`s20_sn2_panel` 16.090(W1 と同じ)。S11 の DME C2v の harness は流していない(VAL だけ)。

### 3.6 テストとシグナル

- 既定のテストの合否は WSL の prod venv で判定する。Windows の QA venv は RDKit の DLL がアプリケーション制御で止まり、GoodVibes と pymsym も入らない。pymsym は 0.3.5 に固定し、既定のテストが版を照合する。
- SIGINT(`/home/user/hfauto_r7/W1/sigint_hcn`): 実行中の HCN の run が rc 130 で終わり、外部プロセスは 0.8 s で消えた。

## 4. 到達表(round 7 の実 run だけ)

`log.jsonl` の decide の reason と注記、minima の diagnostics.json から数えた。run dir は `/home/user/hfauto_r7/` 以下。行番号は round 7 の決定表(17 行)のもの。

| 対象 | 通ったもの | 通っていないもの |
|---|---|---|
| 決定表の行 | 8(VAL7/s6_oh_ch4 の `connection_retry` → `connection_failed`、W5/s6_oh_ch4_rowfix の `connection_retry` → FIND_PATH)、2(W0〜W4/water_same_basin、s19)、3(W1/formaldehyde)、4(R1〜R3 ほか多数)、5(W1・W2・W4/oxalic_two_step、W3/s6_oh_ch4)、6(W1/h2o_hf_inversion、W1/water_dimer_as)、7(W1/hono_walltime_b)、9・10・12・14(R1〜R3 ほか多数)、13(W3/s6_oh_ch4)、15(W1/hcn_noscreen、W1/hf_dimer_swap_noscreen、W3/s5_ch3_o2、W3/s6_oh_ch4)、17(W3/s5_ch3_o2) | 1(端点の DFT 極小が欠けたときだけ発火)、11(崩壊した鞍点)、16(string の次のチャンク)、近道の種が失敗した後の 2 回目の SCREEN、ケースの例外の閉じ込め |
| 検証済みの鞍点(`ts_calc`) | W1・W2/nh3_planar_seed、W1・W2・W4/oxalic_two_step(split1)、W3/s6_oh_ch4(split1) | ゲートを通らず SCREEN へ戻る分岐(ts_calc は order 1 の DFT 鞍点か親が検証した TS なので、χ で拒否されるか freq が失敗したときだけ。専用のコードはなく、ほかの鞍点の失敗と同じ経路) |
| saddle の種 | `screen_ts`(R1〜R3 ほか多数)、`screen_hei`(W2/s13_h2o2_gauche、W2/s18_h3_doublet、W3/s5、W3/s6)、`discovery_ts`(W2〜W4/s10、W3/s5、W3/s6)、`path_hei`(W1 の noscreen 2 本、W3/s5、W3/s6)、`higher_order_retry`(W3/s6 の split2_split2、Hessian は検証済み TS freq)、`saddle_restart`(W3/s5 に 2 回、W3/s6 に 2 回) | — |
| 初期 Hessian | xTB(65 case)、DFT 5 case(W1 の noscreen 2 本、W2/s13、W2/s18、W3/s6)、TS freq 1(W3/s6) | — |
| outcome | elementary、degenerate、same_basin、out_of_window(W1/formaldehyde)、multi_step(oxalic、S6)、barrierless_at_resolution(W1/h2o_hf_inversion、water_dimer_as)、reassigned(W3/s5、W3/s6)、unresolved_within_budget(W1/hono_walltime_b、W3/s5) | blocked_upstream |
| QRC | `minus_is_image`(W2 の NH3、nh3_planar_seed、dme_c2v_seed、S1、S2、S10、S12〜S14、S17、S18。W3・W4 の再実行も)、`end<i>_to_new_basin`(oxalic、S6) | — |
| 極小 | DFT の mode-follow `ts_candidate`(W1・W2/nh3_planar_seed、W1・W2/s4、W3/s6)と `one_side`(W1・W2/dme_c2v_seed)、`image_of`(W2 以後の NH3、nh3_planar_seed、S1、S2、S10、S13、S14、S16〜S18)、`not_reacting`(W2〜W4/s10)、`collapsed_at_dft_from_seed`(W3・W4/s19)、screen の `soft:resolved`(s19)、DFT の `soft:resolved`(S3 のニトロメタン、§22.8) | DFT の Newton の歩による押し出し(非停留の極小が実計算で出ていない) |
| explore | 付け直した IRC の端(W3/s6 2 件、W3/s19 2 件)、relaxation(W3/s5、s6、s19)、CREST 失敗時の seed(W3/s5、s6) | S8・S10 の付け直し(該当する端がない) |
| エンジン | CREST、xTB、ReaDuct NT2、pysis NEB、NWChem の opt・freq・saddle・SP・string(W1 の noscreen 2 本、W3/s5、W3/s6)、ωB97X-D3 の SP(W1/s20_sn2_panel)、CCSD(T) の SP(VAL7/s20_hcn_panel が閉殻、VAL7/s20_ch3o_panel が開殻) | — |

## 5. 計算量

QM 秒は `/home/user/hfauto_r7/tools/jobtime.py` が数える一意のジョブ鍵の合計(複写したジョブを除く)。W4 は同時実行で 1 ジョブの秒数が伸びるので壁時計で比べる。

| 比較 | 対象 | 結果 |
|---|---|---|
| VAL → W2 | 対になる 16 run | 12,743 s / 369 ジョブ → 8,334 s / 325 ジョブ(−34.6%)。QRC の側の opt 4,477 → 2,064 s(529 → 259 歩)、dft の freq 3,357 → 2,432 s、dft の opt 1,153 → 644 s と Hessian 付き 921 → 559 s、rerank の SP 90 → 0 s。経路の SP と saddle は横ばい |
| W1 → W2 | 対になる 13 run | 10,008 s → 6,480 s(−35.2%) |
| W2 → W3 | 失われた seed の DFT(R6) | S19 +2,033 s(17 原子の seed opt、Hessian なし 45 歩)、S6 +536 s(mode-follow の 2 側を含む)、S5 +208 s。予測(3 opt + 3 freq、約 1,371 s)に対して 5 opt + 4 freq |
| VAL → VAL7 | 対になる 24 run | 14,071 s / 467 ジョブ → 13,579 s / 406 ジョブ(−3.5%)、壁時計 14,385 → 11,164 s(−22.4%)。S19 は失われた seed の DFT(R6、+2,308 s)を含む。S19 を除く 23 run は 12,231 → 9,430 s(−22.9%)、壁時計 12,548 → 7,111 s(−43.3%)。種類別: QRC の側 4,477 → 3,051 s(529 → 279 歩)、dft の freq 3,367 → 2,454 s、Hessian 付きの dft の opt 921 → 584 s、rerank の SP 90 → 0 s。経路の SP は 750 → 1,321 s(rank を分けた同時実行でジョブあたりの秒数が伸びる。壁時計は減る) |
| W4 → VAL7 | 同じ 10 系 | QM 13,627 → 13,440 s(−1.4%)、ジョブ 234 → 229、壁時計 10,337 → 10,528 s(+1.8%)。増えたのは S10(185 → 225 s: 上流の入力の違いで QRC の側が 21 → 40 歩、NWChem の AUTOZ の失敗による 2 回目の試行 1 件)と S19(3,893 → 4,053 s: seed の opt 52 → 53 歩)。ほかは ±6% 以内。孤児は全 35 run で 0 |
| W2・W3 → W4 | 壁時計 | S16 61:20 → 52:42、S14 10:35 → 9:38、S1 5:42 → 4:19、HCN 1:25 → 1:01、HONO 2:59 → 2:18、oxalic 38:16 → 33:32、S19 65:16 → 64:53。SCREEN の SP の区間 S14 93 → 55 s(−41%)、S16 315 → 193 s(−39%)、非対称な QRC の対 HCN −21%、acac 1,314 → 1,074 s(−18%)、HONO −19%、oxalic −6%(片側だけが長い)。同じ鍵のエネルギーの差は最大 8.9e-9 Eh、max RSS ≤ 329 MB、孤児 0 |

内訳(VAL → W2): QRC と mode-follow の側の Hessian を正定値のモデルにした(acac の QRC 34/27 → 15/14 歩、S4 の mode-follow の側 36/36 → 13/13 歩)、厳密な像の端点と対称な TS の − 側を計算しない、反応しない組成を DFT に流さない(S10 の dft の freq 861 → 30 s)、mode-follow の側を stage で緩和し直さない(S4 の opt 4 → 2、freq 5 → 4)、rerank の SP を混んだ組だけにした(S10 4 → 0、S19 3 → 0)。

極小の stage を同時に走らせる案は落とした。`/home/user/hfauto_r7/W4/s19_tma_hf2`(途中で止めた run)では 4 つの opt を 1 rank ずつにすると、TMA が 206 s(4 rank で 105 s)、錯体は約 620 s で 4 歩(直列は 383 s で 8 歩)で 2 コアが遊び、43 歩の seed opt は 1 rank で約 2 h の見込みになった。S10 も変わらなかった(dft 45 s 対 48 s)。

## 6. プローブと実 run で決めたこと

| 案 | 根拠 | 決定 |
|---|---|---|
| K1 最小化の初期 Hessian | 計画の規則(射影した固有値が −1e-3 未満なら正定値のモデル)では、DME C2v の 2 本目の鞍点方向(−129.6i)が這い、207 歩で未収束(そのままなら 29 歩。`/home/user/hfauto_r7/W2/superseded_k1neg/dme_c2v_seed`)。xTB の錯体の Hessian も鍵が変わり S18 の熱化学が動いた | 一次の鞍点の Hessian だけをモデルにする。VAL と対の QRC の側は 17.6 → 13.6 歩/側(−23%)。小さな対称の側は速くならない(S10 14.5 → 22 歩、`/home/user/hfauto_r7/knowhow_probe/qrc`) |
| QRC の尾の早期停止 | K1 の後の acac の側は basin に入ってから 7/15 歩と 6/14 歩で、その間の ΔE は 0.004 / 0.009 kcal/mol | 採らない |
| K6 linopt 0(U3-P11) | `/home/user/hfauto_r7/knowhow_probe/opt2.log` と `qrc_all.log`: DME の QRC 143 → 91 s、acac 723 → 662 s は速く、HCN 11 → 15 s、マロンアルデヒド 105 → 148 s、acac の極小 293 → 326 s は遅い | 採らない(閉じる) |
| AFIR | M4 より後の run(SA〜VAL)の AFIR 41 件(一意の discovery)はすべて陰性で生成物 0。M4 の生成物は、今の結合規則ではハロゲン結合錯体 ClI+H3N(W3/s7 の screen basin にある)と、同じ run で NT2 も見つけた恒等 SN2 | W3 で削除 |
| R17 開殻 DFT の SCF 救済 | cgmin の出力は ⟨S²⟩ を出さず、救済したジョブは必ず捨てられる(`/home/user/hfauto_r6/SA/sac_probe/scf/h3p3c/a00_bead_000003_cgmin/stdout.txt`)。実 run での発生は 0 | 開殻 DFT は救済しない(round 7 の run でも発生なし)。r9-S1 で、cgmin の後に通常の SCF を 1 回回して ⟨S²⟩ を読む形にした(§20) |
| K4 CREST の電子温度 | `/home/user/hfauto_r7/knowhow_probe/crest`: OH·CH4 は 1,000 K・3,000 K なら完走するが、組成ごとに低レベルの PES を変える | 採らない。CREST が失敗した組成は seed を使い、失われた状態は DFT に 1 回問う |
| 失われた seed の DFT(R6) | S6 は PBE0 で OH···CH4 が残り、引き抜きの TS と順位に届いた。S19 は DFT でも崩れた(`collapsed_at_dft_from_seed`) | 採用。崩れても障壁なしとは呼ばない |
| 付け直した IRC の端(R5a) | 事前の数え上げ(0.05 Å で 4 件)どおり、未接続が S6 2 → 0、S19 5 → 3 | 採用。許容は 0.05 Å のまま |
| pymsym | prod venv の版は 0.3.5 | 固定。入れ直さない |

## 7. round 7 で直した不具合

- **締め切りの意味**(`/home/user/hfauto_r7/W1/hono_walltime`): 予算を使い切った case が、拒否されるジョブを重ねて行 15・17 まで進んだ。`Deadline.expired()` を「どのジョブも始められない(残り 60 s 未満)」にして、行 7 で閉じるようにした(`hono_walltime_b`)。
- **freq の電子状態**(`/home/user/hfauto_r7/W3/superseded_nofreqguess/s6_oh_ch4`): UKS の OH···CH4 錯体の freq が SCF をゼロから始めて OH の別の π 成分に落ち(opt より 7.2e-5 Eh 高い)、`is_minimum` が `state_mismatch` で seed を捨てた(同じ run の QRC の側も)。freq は opt・saddle の収束 vectors から始める(`scf_guess`、鍵に入る)。golden の回帰テストを足した。

- **予算を残した打ち切り**(`/home/user/hfauto_r7/VAL7/s6_oh_ch4`): 両側が 2 回とも同じ basin に戻る TS(別の過程の鞍点)で、行 8 が saddle の試行を残したまま `connection_failed` で閉じていた。試行が残れば行 12〜17 に進み、次の saddle 探索がその TS の主張を捨てる(`/home/user/hfauto_r7/W5/s6_oh_ch4_rowfix` で FIND_PATH に進むのを確認)。回帰テストを足した。

## 8. 撤回・訂正した値

- HONO trans → cis の ΔG‡ 24.11(旧 `hono_isomerization_v3`、ZPE の二重計上)。HCN → HNC の ΔG_rxn 18.05(内部フォールバック。正しくは 12.68)。I⁻ + CH3I の ΔG‡ 7.53(対称数の誤判定。正しくは 6.8817)。M2 以前の R2 12.2265(キラリティ m = 2 を数えていない。今は 11.8158)。
- S4 の 1 位だった CH2OH• のねじれ(3.878): 宣言していないねじれは仮説にしないので、順位から消える。
- 期待値の訂正(system の見出し): S15 は直接の 1,3-H 移動の TS(ΔE‡ 55.3、δG_eff 52.0)で elementary、S17 は障壁 1.256 が解像度 1.0 を超えるので degenerate、ホルムアルデヒドの順方向は ΔE_rxn 53.75 で窓の外。

## 9. 残る課題

1. **実計算で通っていない分岐**: 2 回目の SCREEN、ケースの例外の閉じ込め(単体テストだけ)、Newton の歩で非停留の極小から抜けること(S3、§22.7。DFT の soft の分岐はニトロメタンで通った、§22.8)と、非停留の高次の鞍点からの `higher_order:not_stationary`(単体テストとプローブだけ)、maxiter からの再開の上限(G1-P3。単体テストと QM なしの評価だけ、§21.6)。行 6 の中点による高密度化は対称な経路でしか通っておらず、隠れた障壁を見つけた例はない。round 7 で未到達だった行 11 は M5 で削除した(§15)。round 7 の行 16(string の次のチャンク)は M4 の S5 split2 で通り(SCREEN の SP が未収束で種がなく、string の種を χ で拒否した後)、S2 で FIND_PATH の 1 行にまとめた。行 1 は S2 で実計算に通った(§21.6)。
2. **既定の順位の精度**: PBE0/SVPD は CCSD(T) から最大 3 kcal/mol ずれ(SN2 −3.0、HONO +2.0)、H 引き抜きは分離基準の ΔE‡ 0.79 と文献 約 5 を大きく下回る(round 8 のプローブの CCSD(T)/def2-TZVPD は 6.96 で、差は −6.2)。M9 で、既定の順位を M06-2X-D3(0)/def2-TZVPD のエネルギー層にした(§19)。実 run の 5 反応で障壁の誤差は平均 0.83、最大 1.79 kcal/mol。ΔE_rxn は改善しない(平均 1.17、最大 2.30。HCN)。
3. **開殻の会合**: S5 の UKS CH3···O2 は ⟨S²⟩ 1.71 で、障壁なしかどうかを判定できない。多参照かスピン射影の扱いが要る。M8 で、会合を分離した単量体からの緩和スキャンで問うようにし、S5 の会合は障壁なし(UKS、AP なし)で順位が付いた(§18.4)。結論が出なかった主因は多参照性ではなく問いの形だった。
4. **失われた seed の費用**: S19 の seed は DFT の 1 歩目で状態を離れたのに、既知の basin に入るまで 44 歩(約 33 分)払った。M6 で、この seed のラベルの違いはしきい値の帯の中(0.006 Å)なので、失われた状態にしなくなった(§16)。本物の R6 の seed の opt(Hessian なし)の費用は残っていた。S3 で、seed も xTB Hessian の正定値のモデルから始めるようにした(G4-P2、§22)。
5. **鞍点探索の空回り**: 柔らかい H 引き抜き・ラジカル会合では string の前に saddle が maxiter で 2 回止まる(S5 で saddle 826 s)。
6. **対称な鞍点の mode-follow**: ± の側は厳密な像なのに 2 本とも opt と freq を払う(S6 で 329 s、NH3 は 21 歩/側)。
7. **近直線の錯体の熱化学**: S18 の H···H2 は、軸から 0.006 Å 外れた極小が非直線(振動 3 本、δG_eff 4.091)、正確に直線の極小が直線(4 本、6.360)と扱われる(`vibrations._EXTERNAL_RANK_TOL` = 1e-3、`/home/user/hfauto_r7/W2/superseded_k1neg/s18_h3_doublet`)。直線性を対称性で決める必要がある。M7 で、点群 1 つから直線性を決めるようにした。2 つの run の差は 2.27 から 0.006 kcal/mol になった(§17)。
8. **分割の深さ**: W3 の S6 の CH3OH + H 側は、ほぼ縮退した緩い錯体(ΔE_rxn ≤ 0.11)の分割で `max_split_depth` に達し、順位が付かない(VAL7 ではその 2 つの極小が 1 つの basin にまとまった)。M1 で偽の結合変化が、M2 で同じ状態の分割がなくなった(W3_s6 の再生、§11.3、§12.3)。
9. **平らな PES の引き抜きの TS の再現性**: S6 の OH···CH4 → CH3···H2O は W3 だけが −481.1i の TS に届き、VAL7・W5 は OH···CH4 の向きを変える鞍点(−73.0i、−58.6i、両側が同じ basin)に収束して試行を使い切った。虚モードが仮説の結合変化を担わない鞍点を QRC の前に見分け、結合変化の方向に拘束して探す必要がある(その鞍点の QRC は 2 振幅 × 2 本)。M3 の χ と M4 の ρ で、7 つの文脈のうち 6 つで同じ引き抜きの TS に届いた(§14.5)。その TS の領域の鞍点は、Cartesian の勾配で見るとどれも停留していない(ΔE_N 1.0e-4〜1.5e-4 Eh、§22.3)。
10. **NWChem の AUTOZ の失敗**: VAL7/s10 の QRC の側の opt が 21 歩目で内部座標の再構築に失敗し(53 s)、2 回目の試行は始点から Cartesian でやり直して 19 歩で収束した(43 s)。VAL7 で 2 回目の試行はこの 1 件。r9-S1 で最後の frame から続けるようにした(§20)。

## 10. 再生ゲート

振る舞いを保つ整理の合格基準で、各 wave の受け入れの土台でもある(分析 [2026-09-30_remaining_issues_analysis.md](reviews/2026-09-30_remaining_issues_analysis.md) §3 X8、§4.0)。記録済みの run を今の作業ツリーで JobStore から再生し、新しい QM ジョブが 0 件で、結論の記録が変わらないことを確かめる。

- **道具**(repo の外、`/home/user/hfauto_r9/`)
  - `tools/replay.py`: ゲート本体。
  - `tools/runs.json`: 再生する run の一覧。
  - `runcase.sh`: 1 本ずつの実行。
  - `tools/jobtime.py`: core 秒、歩数、SCF 回数の集計。
  - `tools/backfill.py`: 旧 JobStore が保存しなかった ladder の最後の失敗を、source の複製に一度だけ保存する(下記)。
- **再生の手順**: source を `/home/user/hfauto_r9/<wave>/replay/<name>` に複製し、`runcase.sh` で pipeline の最初の stage(`structures`)から再ステージする。パネルの run では続けて panel pipeline を再開する(`RESUME=1`)。最初の pipeline を再ステージするとパネルの stage は stale になる。`--from panel_sp` は、実行順でそれより後にある最初の pipeline の stage まで stale にして入力を失うので使わない。ロック、CAP、log は実 run と同じである。
- **比べるもの**
  - stage ごとの `run_state.json`: status、hits、misses、failures_by_kind。再生の開始より後に始まった stage だけを「再生した」と数える。スキップされた stage は hits 0・misses 0 に見えるので、それを合格にしないためである。
  - 再生したすべての stage の manifest にある結論の record(species、minimum、discovery、reaction、species_thermo、reaction_thermo、report)。artifact id ごとに、欄の単位で比べる。
  - 複製の中で走ったジョブ(engine、kind、鍵、core 秒)と、source がその鍵を返せなかった理由(鍵がない、結果を残していない、保存された結果を使えなかった)。
- **比べないもの**
  - calculation の artifact: ジョブの結果そのものなので、misses で見る。
  - reason と reasons の文字列: 結論ではない(§4.4)。
  - report.html の FileRef: run id を含む。
  - manifest の created_at と run_id、壁時計。
- **合格(`--strict`)**: 次をすべて満たすこと。
  - 再生したすべての stage で misses 0。
  - source でジョブを走らせた stage では hits > 0 で、その stage をすべて再生している。
  - status と記録が同一。
  - `--from hcn=paths` のような浅い再ステージは、minima stage を再生していないので不合格になる。ranking.csv だけを比べることはしない。
- **振る舞いを変える wave**: `--allow-new a,b` で挙げた run は報告だけにし、ほかの run は strict で判定する。報告した run の新しいジョブと差分は、その wave が化学的な理由で説明する。
- **BASE と CUR**
  - BASE: `runs.json` の src。r7 の証拠の run で、読むだけである。
  - CUR: `/home/user/hfauto_r9/CUR/<name>`。受け入れた最新の再生の複製への symlink で、wave を受け入れたら `ln -sfn /home/user/hfauto_r9/<wave>/replay/<name> /home/user/hfauto_r9/CUR/<name>` で張り替える。次の wave は `--src CUR` で再生する。
- **対象**: VAL7 の全 run と、S6 の基準 2 本。
  - `W3_s6`: −481.1i の TS に届いた W3 の run。
  - `W5_s6`: 行 8 の修正のあとの S6。
- **除外**(理由は runs.json にも書いてある)
  - `hono_walltime`: Deadline は実時間で判定され、キャッシュの hit は瞬時に返るので、行 7 が発火するかどうかが再生の速さで変わる。
  - `s9_*`: FeCl3·CH4。ジョブの前か、CREST だけで止まる。遷移金属の DFT は未検証と宣言している。
  - `smoke`: run ではない。
- **source の既知の差**
  - 2c3c894 までの JobStore は、ladder の最後の失敗のうち旧 TERMINAL の 4 種以外(geometry_maxiter、scf_not_converged、nonzero_exit)を保存しなかった。該当は VAL7/s5 が 9 件(saddle の maxiter 3、xTB 5、CREST 1)、VAL7/s6 が 5、W3_s6 が 8(saddle 3)、W5_s6 が 6(saddle 1)で、ほかの run は 0 件。
  - そのまま再生すると、これらは再実行される。NWChem の saddle は別の軌跡をたどる(NXTVAL)。実測では VAL7/s5 で saddle 2 本の再実行と新しい鍵の saddle 2 本(計 3,033 core 秒)が走り、`rxn_discovery_b3febc8052` が reassigned_step から unresolved_within_budget に変わった。これは HEAD の自己一致ではない。
  - そこで `tools/backfill.py` で、この 4 本の複製(`/home/user/hfauto_r9/BASE_SRC/<name>`)に、source が実際に受け取った Failure を一度だけ保存した。最後の attempt を adapter 自身の parse で読み直し、鍵を刻印する。ladder が続けるはずの失敗と壁時計の失敗は拒否する。種類ごとの件数は 4 本とも source の run_state の failures_by_kind と一致した。
  - VAL7/s6 は行 8 の修正(1bd2ff5)より前の run である。HEAD では FIND_PATH に進み、新しいジョブが出る。今のコードでの S6 の基準は W5_s6 である。
- **限界**: 再生が保証するのは、同じ鍵に同じ結果が返る範囲の一致だけである。並列の NWChem と CREST はビット単位では再現しないので、新しいジョブが 1 本でも出ると、その下流は run 間の比較(結論の水準、§4.4)になる。
- **時間**
  - core 秒(duration × argv の `-np`。CREST は `-T`)、歩数、SCF 回数で比べる(`jobtime.py`)。
  - 壁時計は参考値にとどめる。`runcase.sh` が開始と終了の load average を log に残す。
  - CAP は SIGTERM で送り、60 s 後に SIGKILL にする(`timeout -k 60`)。
- **runcase.sh の追加**
  - `FRESH_ENGINES=<engine,...>`: SRC の複製から、指定した engine のジョブだけを消す。上流を共有したまま、その engine の軌跡を新しく走らせるため(A/B)。
  - `SITE=`: site を差し替える。
  - `RESUME=1`: 既存の run を `--from` なしで再開する(パネルの連鎖)。
- **r9-PRE の結果**(作業ツリー = 2c3c894 + PRE-1。36 本、1 本ずつ)
  - 道具の確認: hcn を 2 回続けて strict で再生し、2 回とも PASS。陰性対照 3 本はすべて FAIL になった。
    - (a) paths の 1 ランク SP の job dir を 1 つ消す: `paths: 1 misses`、鍵 07c7bd7f が `key not in the source` として出る。
    - (b) source の outcome を書き換える: `paths/hcn_to_hnc.payload.outcome` の差分が出る。
    - (c) `--from hcn=paths`: `dft (minima) ran 4 jobs in the source but was not re-staged`。
  - BASE(`/home/user/hfauto_r9/BASE/replay`、r7 の source と BASE_SRC から): 34 本が strict で PASS(misses 0、記録が同一)。報告にとどめた 2 本は次のとおり。
    - VAL7/s6: 行 8 の修正で FIND_PATH に進み、新しいジョブ 10 本(4,649 core 秒。string 1、saddle 2、opt 2、freq 3、xTB freq 2)。結論は同じで、鞍点が −73.0i から −73.2i に、障壁の出所が screen から string に変わっただけ。
    - W3_s6: W5 の整理(1bd2ff5)で PathProfile と DiscoveryResult から job_key を除いたので、W3 が保存した結果が検証を通らない(`extra_forbidden`)。ReaDuct 8、NEB 4、string 1 の計 13 本(1,839 core 秒)を再実行した。結論の差は、TS のファイルの attempt の番号と、W5 で除いた dzpe_act_kcal 欄だけで、−481.1i の TS を含む分岐はそのまま。
    - パネルの 3 本は `RESUME=1` で再開した(上記)。
  - 自己一致(`tools/replay.py PRE --src /home/user/hfauto_r9/BASE/replay --strict`): 36 本すべて PASS。再生した stage はすべて misses 0 で、hits は BASE の同じ stage のジョブ数と等しい(計 971)。failures_by_kind と記録も同一で、複製の中で走ったジョブは 0 件(jobtime)。壁時計は計 205 秒だった。
  - CUR は `PRE/replay/<name>` を指す。以降の wave は `--src CUR` で再生する。
- **実行例**(`cd /home/user/hfauto_r9`)
  - `tools/replay.py PRE --src /home/user/hfauto_r9/BASE/replay --strict`: 準備した BASE から全体を再生する。
  - `tools/replay.py M1 --src CUR --allow-new W3_s6`: 振る舞いを変える wave。
  - `tools/replay.py PRE-c --runs hcn --from hcn=paths --strict`: 浅い再ステージが不合格になることの確認。
  - 結果は `/home/user/hfauto_r9/<wave>/replay/summary.json` と `summary.txt` に出る。FAIL があれば終了コードは 1。

## 11. M1 原子の対応(X2-1)

分析の §3 X2-1 と §4.1 M1 にあたる。規則は [design.md](design.md) §6.2 の「原子の対応」。コードは r9-M1 のコミット。

### 11.1 変えたこと

- 中間体: 既知の basin に落ちた井戸は、到達した構造をジョブなしの species `spc_<reaction_id>_intermediate` にして basin に加える。分割の子はその添字で続く(`drivers/reaction_case/connection.py`。接続と中間体の action を `actions.py` から移した。移しただけの段階では振る舞いは変わらない)。
- 発見: `DiscoveryRecord.source_species` を足し、仮説の両端を discovery 自身の (`source_species`, `product_species`) にした。`hypotheses.pick_endpoints` と、メンバーの対から端を選ぶ分岐は削除した。explore の生成物は、添字付きの結合が等しい既知の構造にだけ同定する(縮退の生成物だけの特例をやめた。CUR の 14 件の生成物では同定は変わらない)。
- `identity.basin_coords`: 端点の構造が basin と同じ状態ラベルを持つとき、原子を同じ WL クラスの原子にだけ対応させる。統合で選んだ方法である。到達した構造を記録する案もあったが、新しい記録が要り、全 run の MinimumRecord が変わるので採らなかった。
  - 理由: VAL7/s6 と W5_s6 の bfac(xTB の IRC の端)は、それを含む PBE0 の basin から 0.957 Å 離れている。元素だけで並べると H1 が O に、H6 が C に対応し、`rxn_discovery_f9bc3a1cd8` の結合変化が 8 本になった。IRC 自身の変化は +C0–H2 +C0–O5 −C0–H4 −H2–O5 の 4 本である。この構造を回帰テストにした(`test_basin_coords_keep_the_bonds_of_a_member_far_from_the_basin`)。
  - 状態ラベルが違う端(basin への緩和で結合が変わった seed)は、今までどおり元素だけで対応させる。そこでの結合変化は化学的な結果である(S5 の b3e2f366b6: xTB の CH2OOH が PBE0 で CH2O + OH に緩和した)。

### 11.2 監査(QM なし)

道具と出力は `/home/user/hfauto_r9/M1/` にある(`audit/`、`audit_ends/`、`integration/`)。
- `tools/audit_ends.py` は、発見から作った仮説ごとに、case の端の結合変化と discovery 自身の添字付きの結合変化を比べる。
  - HEAD のコードでは 17 件中 6 件が MISMATCH だった(W5_s6 と s6_oh_ch4 の bfac の端)。
  - 作業ツリーでは 0 件である。M1 の再生(記録された `source_species` を使い、補完は 0 件)でも 0 件だった(`integration/audit_ends_M1_replay.txt`)。
  - DFT_STATE の 1 件(S5)は上記の化学的な変化である。
- CUR の 25 の非代表の端のうち、新しい `basin_coords` で変わるのは bfac の 2 件だけだった(`audit_ends/basin_coords_old_vs_proposed.txt`)。

### 11.3 再生(`tools/replay.py M1 --src CUR`、36 本、壁時計 8 分)

- 発見を持たない 25 本は、strict の条件(misses 0、記録が同一)で PASS した。
- 発見を持つ 8 本(nh3_planar_seed、s4、s20_ch3o、s7、s8、s10、s19、s5)の差は、新しい欄 `source_species` だけだった。新しいジョブは 0 件である。
- s6_oh_ch4 と W5_s6: 上の欄のほかに、`rxn_discovery_f9bc3a1cd8` の端が (p02, 05c1) から、discovery 自身の端である (p02, bfac) に変わった。新しいジョブは 0 件で、outcome と数値は同じ。bfac の添字で並べた basin の構造が、05c1 の代表の座標と一致するためである。
- W3_s6(分析の X2 検証 (b)):
  - 中間体が bfac の代表から、到達した QRC の側 `spc_rxn_discovery_f9bc3a1cd8_intermediate` に変わった。
  - split1 は elementary_step のまま(親の TS を JobStore から再検証)。
  - split2 は `torsional=True`(結合変化なし)で、barrierless_at_resolution になった。以前は偽の H の入れ替えで multi_step になり、split2_split1 と split2_split2 を生んでいた。この 2 つは消えた。
  - 新しいジョブは、端が変わった split2 の SCREEN の 12 本(NEB 1、DFT SP 11。計 200 core 秒)。hits は 130 から 80 に減った。
  - 付随する差として、`rxn_discovery_460166732f` の dG_eff が 32.41 から 31.60 kcal/mol に下がった。消えた偽の連鎖が登録していた CH4O+H の状態の極小 2 つがなくなったためである。その G は、この反応の反応物の極小より低かった。状態の G(最小の G)が変わったので dG_eff が変わった。ΔG‡ と ΔE‡ は同じ。
- CUR は `M1/replay/<name>` に張り替えた。

### 11.4 実計算: マロンアルデヒド(反応の宣言なし、C2v の種から)

- 入力は `/home/user/hfauto_r9/inputs/malonaldehyde_c2v_undeclared_{ab,ba}.yaml`(known_endpoints)。
  - 種は VAL7/s14 の TS の freq 構造(−1100.1i)を C2v に対称化したもの。
  - _ab と _ba では、basin の代表になるエノールの登録順を入れ替えてある。
- 結果(`/home/user/hfauto_r9/M1/real/malon_c2v_{ab,ba}`、1 本ずつ):

| run | 代表 | 仮説の端 | outcome | 記録 | ΔE‡ | dG_eff | 壁時計 | QM |
|---|---|---|---|---|---|---|---|---|
| _ab | `min_r1_enol_a_ff60c312` | `ts_c2v_mf1` → `ts_c2v_mf2`(mode_follow) | degenerate_rearrangement | `connection:degenerate`、`qrc1:minus_is_image` | 2.01 | 0.0 | 37:55 | 8 ジョブ、9,084 core 秒 |
| _ba | `min_r1_enol_b_ff60c312` | 同上 | degenerate_rearrangement | 同上 | 2.01 | 0.0 | 37:58 | 8 ジョブ、9,097 core 秒 |

- 2 本とも、ts_c2v が `follow1:ts_candidate` になり、側が 1 つの basin に入った。case は mode-follow の鞍点(`ts_calc`)を検証し、QRC は + 側だけを走らせた。結論は S14(degenerate、dG_eff 0.0)と同じである。
- 見かけの合格でないことの確認:
  - 2 本とも成功した。
  - 代表が違う。
  - case の source は mode_follow で、宣言した反応ではない。
- 実証していないこと: HEAD では、代表の添字が mf2 と同じ run で仮説が消えるはずだった(分析 G8-3 の推論)。HEAD でこの対は走らせていない。
- 費用: DFT の mode-follow は対称な鞍点でも ± 両側を opt と freq にかける(freq 4 本で 69%)。QRC の像の規則は minima の mode-follow には入っていない。

## 12. M2 判定の粒度(X2-2)

分析の §3 X2-2 と §4.1 M2 にあたる(G8-P1、G5-P2 の簡素版)。規則は [design.md](design.md) §6.2 の「判定の粒度」。

### 12.1 変えたこと

- CONNECT と VALIDATE_INTERMEDIATE は、極小を case のキーで比べる。両端の (composition_id, state_label) が異なる case は状態、同じ状態の case は basin がキーである。`gates.connection` には、このキーをそのまま渡す(gates は変えていない)。
- 両側が同じキーで別の basin の TS は `same_state` とし、この case の反応ではない鞍点として棄却する。行 8 は振幅を広げず、saddle の試行が残れば探索を続け、なければ `connection_failed` にする。振幅を広げる再試行は、同じ basin の場合だけに残した。
- 同じ状態の case では、端点と同じ状態の別の basin を、次の 2 条件をともに満たすときだけ、その端点とみなす(G5-P2)。
  - その端との \|ΔE\| < resolution_kcal。
  - それを運んだ経路(プロファイルの井戸は端までのプロファイル、QRC の側は TS)に resolution_kcal 以上の山がない。
- `ConnectionClaim.minima` は実際に結ばれた極小を持ち、elementary の `ReactionRecord.minima` は書き換えない。
- 削除: 「中間体は他方の端点の配座でもよい」という扱い。結合変化の case では、端点と同じ状態の井戸は中間体にならない。

### 12.2 期待値の訂正(QM なしの監査)

CUR のうち paths を持つ 33 本の 41 case について、記録と case の log から QRC の側と中間体の状態を端点の状態と比べた(`/home/user/hfauto_r9/M2/audit/m2_audit.py`、出力は同じ場所の `audit.txt`)。結論が変わるのは、次の 4 本の 6 case だけである。ほかの case は、QRC の両側が端点の basin か 1 つの basin にあるか、QRC に進んでいないので、判定は変わらない。確定は再生(`tools/replay.py M2 --src CUR --allow-new W3_s6,W5_s6,s6_oh_ch4,s5_ch3_o2`)で行う。

| run | case | 今まで | M2 の期待 | 化学的な理由 |
|---|---|---|---|---|
| W3_s6 | `rxn_discovery_f9bc3a1cd8` | multi_step(split1 elementary、split2 barrierless) | elementary(`connection.minima` = (p02, bfac))、子 case なし | QRC の側 bfac は、端点 05c1 と同じ CH3OH + H の状態(ΔE 0.10 kcal/mol の緩い錯体の別 basin)である。分析 X2 の検証 (a) |
| W3_s6 | `rxn_discovery_460166732f` | reassigned_step | multi_step。CH4 + HO を経る 2 段 | −1347.6i の TS の側は、CH3OH + H(端点の状態)と CH4···HO(新しい状態、+14.13)である。TS を渡された子 CH4 + HO → CH3OH + H は再生、CH3 + H2O → CH4 + HO は新たに探索する。検証 (d) |
| W5_s6、s6_oh_ch4 | `rxn_discovery_f9bc3a1cd8` | reassigned_step | multi_step。CH4 + HO を経る 2 段 | 同じ −1347.6i の TS と同じ側である。検証 (c)。子 split1 が 790468f505 の結果を再利用するのは G8-P7 で、M2 では新しく探索する |
| s5_ch3_o2 | `rxn_discovery_b3febc8052` | reassigned_step | unresolved_within_budget(`connection_failed`) | −72.0i の鞍点は、CH3···O2 の 2 つの basin(0.00 と +0.40)を結ぶ。どちらも反応物と同じ CH3 + O2 の状態なので、会合の TS ではない。saddle の試行は使い切っている。検証 (e) |
| s5_ch3_o2 | `rxn_discovery_b3e2f366b6` | reassigned_step | multi_step。CH2OOH を経る 2 段 | −1866.6i の 1,3-H 移動の TS は、CH3OO• の状態(反応物の状態)の別 basin(+22.48)と CH2OOH•(新しい状態、+15.54 の DFT 極小)を結ぶ。TS を渡された子 CH3OO• → CH2OOH• は再生、CH2OOH• → CH2O + OH は新たに探索する |

- oxalic_two_step は変わらない(検証 (h))。cTc → tTt は同じ状態の case なので basin で判定する。中間体 tTc は TS(ΔE‡ 15.08)の向こうにあり、分解能の規則にはかからない。
- 見かけの合格でないこと: 単体テストは、Registry が別々に登録した basin(手で統合していない)で、キー、`ConnectionClaim.minima`、分割の有無、次の決定(探索を続けるか、打ち切るか)を確かめる(`tests/unit/drivers/test_reaction_case_connection.py`)。`BASIN_A` と `_DE_HARTREE` は変えていない。
- 既知の制限: WL の状態ラベルは立体を区別しないので、ジアステレオマーを同じ端点とみなしうる(G8-P6、未対応)。今の検証系に不斉中心はない。

### 12.3 再生と実計算(r9-M2、`/home/user/hfauto_r9/M2/`)

- strict(`tools/replay.py M2 --src CUR --strict`、上の 4 本を除く 32 本、壁時計 7 分): 32 本すべて PASS。misses 0、新しいジョブ 0、記録は同一である。oxalic_two_step は multi_step のまま(検証 (h))。hcn、hono、nh3_inversion、water_*、S1〜S4、S13〜S18 の δG_eff は変わらない(差 0)。
- W3_s6(`--allow-new W3_s6`、壁時計 27 分。新しいジョブ 15 本、5,319 core 秒。すべて下の split1 の探索)
  - `rxn_discovery_f9bc3a1cd8`: multi_step から elementary_step になった(−1417.9i、ΔE‡ 48.10)。`connection.minima` = (p02, bfac)で、`ReactionRecord.minima` は仮説の端 (p02, 05c1) のままである。`f9bc…_split*` の記録と中間体の species は消えた(検証 (a))。
  - `rxn_discovery_460166732f`: reassigned_step から multi_step になり、CH4···HO(+14.13)を経る 2 段になった(検証 (d))。
    - split2(CH4 + HO → CH3OH + H)は、親の −1347.6i を JobStore から再検証して elementary(δG_eff 41.25)。
    - split1(CH3 + H2O → CH4 + HO)は新しい探索で、−485.8i の TS に届いた。QRC の側の一方は、端と同じ CH4 + HO の状態の別の basin(+13.40。端は +14.13)で、状態の単位で elementary(δG_eff 15.92)になった。
  - 順位は 790468f505(1.48)、460166732f_split1(15.92)、460166732f_split2(41.25)、f9bc(48.48)の順。消えた f9bc の split2(barrierless、torsional)は順位から外れた。
  - CUR/W3_s6 を `M2/replay/W3_s6` に張り替えた。
- s5_ch3_o2(`--allow-new s5_ch3_o2`): CAP(3,600 s)で paths が止まり、paths の記録は残らなかった(新しいジョブ 32 本、11,279 core 秒)。case の log から分かることは次のとおり。
  - `rxn_discovery_b3e2f366b6` は multi_step になった(`qrc1:end0_to_new_basin`、CH2OOH•)。split1(CH3OO• → CH2OOH•)は −1866.6i を再検証して elementary である。
  - split2(CH2OOH• → CH2O + OH)の新しい TS は、CH2O + OH(端の状態)と HCO• + H2O(新しい状態。側の opt は maxiter 100 に達し、4,009 s かかった)を結ぶ H 引き抜きの鞍点だった。そのため split2 はさらに分割された。split2_split1(CH2OOH• → HCO• + H2O)の saddle 探索の途中で CAP に達した。
  - この run の順番では b3febc8052 に進まなかった。そこで、paths を `reaction_ids: [rxn_discovery_b3febc8052]` に限った pipeline(`/home/user/hfauto_r9/pipelines/discover_b3febc.yaml`)で CUR の複製を paths から再ステージした(`M2/real/s5_b3febc`)。hits 25、misses 0、5 秒で終わった。注記は `qrc1:same_state:sides_same_basin` で、結果は unresolved_within_budget(`connection_failed`)である。reassigned_step にはならない(検証 (e))。
  - CUR/s5_ch3_o2 は M1 のまま(再生が途中で止まったため)。
- 実証していないこと
  - W5_s6 と s6_oh_ch4 の f9bc(検証 (c))は再生していない。予算(壁時計 2 h、QM 1.5 h)を W3_s6 と s5 で使い切ったためである。W3_s6 の 460166732f は、同じ −1347.6i の TS、同じ両端 (p02, bfac)、同じ中間体の状態を持つので、判定の経路は同じである。ただし、子の新しい探索の結果は run ごとに違いうる。
  - s5 の b3e2f366b6 の子の連鎖の完了。
- 見えた課題(M2 の外): 行 5 の分割は、TS が端の状態と、それより低い別の状態を結ぶときにも中間体を作る。S5 の split2 では、CH2O + OH から 26.1 kcal/mol 下った HCO• + H2O(QRC の側の opt の最終エネルギーの差)を経る R → I → P になる。これは basin の単位の頃からの GEN-05 の規則で、M2 では変えていない。

## 13. M3 反応モード性 χ(X1 / G1-P2)

分析の §2.1 G1-P2、§3 X1 と §4.1 M3 にあたる。規則は [design.md](design.md) §3 の `reaction_mode_character` と §6.2 の VALIDATE_AND_CONNECT(M3 の時点では VALIDATE_TS)。

### 13.1 変えたこと

- `gates.reaction_mode_chi` / `reaction_mode_character`: χ = ‖Q_Bᵀq̂‖。q̂ は質量の重みを外して正規化した Cartesian の虚モード、Q_B は仮説で変わる結合の Wilson 伸縮ベクトルが張る空間の正規直交基底である。結合が変わらない仮説は宣言座標の勾配との \|cos\| を使い、どちらもない仮説(ねじれの自動検出)にはかけない。
- Q_B は QR ではなく SVD で作る。伸縮の組は線形従属になりうる。共線の 3 原子の 3 結合では、QR は余分な方向を足して並進に χ 0.17 を与える。単体テストは、SVD ではこれが 0 になることを確かめる。
- `validate_ts`: `is_first_order_saddle` を通った後にかける。結合は case のラベル付きの両端(`ctx.ends`)から取り、QRC の側からは取らない。不合格なら注記 `ts_rejected:not_reaction_mode` を残し、`last_saddle='failed'` にする。試行は数えたままで、QRC はしない。
- `REACTION_MODE_MIN` = 0.3 は gates のモジュール定数で、YAML では変えない。値は 13.2 の分布の空白から決めた。

### 13.2 分布(QM なし、`/home/user/hfauto_r9/tools/chi.py`、出力は `/home/user/hfauto_r9/M3/chi/`)

- 道具は、case ごとに種(case のフォルダの単一フレームの xyz、screen_idpp から始めた NEB の TS、`low_level_ts`)から始まった DFT saddle ジョブ、その maxiter の最終フレームからの再開、`ts_calc` を集める。そして、それぞれの最終構造での DFT freq の χ を、case の両端(driver の `_endpoint`)の結合変化で測る。
  - reassigned の記録は `minima` を TS の側に書き換えるので、仮説の極小は端点の species の basin から戻す。
  - CUR の TS freq 38 本のうち 35 本を case に割り当てた。残る 3 本は W3_s6 の JobStore にある、M1 で消えた f9bc の分割の子のものである(今の記録に case がない)。
- 結果(CUR、`chi_CUR.txt`)。接続した TS は、いずれかの case が connected の outcome(multi_step を含む)で主張したもの。

| 区分 | 接続した TS(χ) | 接続しなかった一次の鞍点(χ) |
|---|---|---|
| H 移動(HCN、HONO → HNO2、ホルムアルデヒド、CH3O•、S6 の引き抜き −481.1i と −485.8i、PT 2 件、S10、H + H2) | 13 件、0.640〜1.000 | S6 の 790468f505(W5_s6、s6_oh_ch4)の −73.2i、−73.0i、−58.6i: 4 件、0.002〜0.023 |
| 重原子の結合変化を含む(SN2 の Cl と I 3 件、S6 の C–O の置換 5 件、S5 の O–O 開裂を伴う H 移動 −1866.6i) | 9 件、0.776〜0.930 | — |
| 宣言座標だけ(HONO のねじれ、\|cos\|) | 1 件、0.571 | — |
| 会合(S5 の b3febc8052 −72.0i、CH3 + O2 → CH3OO•) | CUR(M1 の記録)では reassigned だが順位なし(spin_contaminated) | M2 の再ステージ(`chi_M2_s5_b3febc.txt`)では same_state で未接続: 0.004 |
| 結合変化も宣言もない(S12、S13、S17、NH3、oxalic など計 10 件) | ゲートはかからない | — |

- 空白: 接続しなかった鞍点は ≤ 0.023、接続した TS は ≥ 0.571 である(結合変化のある TS だけなら ≥ 0.640、重原子を含むものは ≥ 0.776)。0.3 はこの空白の中にあり、§2.1 の予想と一致する。
- 見かけの合格でないこと
  - S6 だけで決めていない。順位のある接続した TS 23 件(重原子の標本 9 件を含む)と、接続しなかった鞍点 5 件(S6 の 4 件と、M2 の判定での S5 の −72i)の分布から決めた。
  - χ は QRC の側から測っていない。縮退反応(SN2、PT、H 交換)は、ラベル付きの両端の結合変化で 0.83〜1.00 になる。
- 順位のある接続反応(CUR の 23 件)で χ < 0.3 のものは 0 件である。
- 空白に最も近い鞍点: M2 の再生で S5 の split2(CH2OOH• → CH2O + OH、O–O 開裂)が見つけた −299.2i は χ 0.257 である(`s5_split2.txt`。記録が書かれる前に CAP で止まった run なので、親の記録から端を組み立てて測った)。
  - 構造(M3 の再生での −299.3i)は CH2OOH• そのもの(C–O 1.32 Å、O–O 1.50 Å。極小の O–O は 1.42 Å)で、CH2OOH• の極小の 0.8 kcal/mol 上にある。虚モードの大半は CH2 のねじれで、O–O の伸縮は一部だけである。
  - M2 の QRC の両側は CH2O + OH と HCO• + H2O で、どちらの側も CH2OOH• に戻らなかった。
  - M3 ではこれを拒否し、split2 は O–O 開裂の探索を続ける(§12.3 の「見えた課題」)。0.3 との差は 0.04 しかないので、この種の鞍点が増えたら分布を測り直す。
- χ はラベル付きの端点に依る(M1 の添字の連続性が前提)。M1 の前の VAL7 の記録(`chi_VAL7.txt`)では、s6 の f9bc の −1347.6i は、端点の添字の誤りのために 0.289 だった。M1 の端点では 0.798 である。

### 13.3 再生と実計算(統合、`/home/user/hfauto_r9/M3/`)

- strict(`tools/replay.py M3 --src CUR --strict`、下の 3 本を除く 33 本、壁時計 7 分): 33 本すべて PASS。misses 0、新しいジョブ 0、記録は同一である。
- 実計算 S6(`M3/real/s6`。CUR/s6_oh_ch4 の複製を paths から再ステージ。G1-P1 なし。壁時計 13.6 分、新しいジョブ 15 本、2,883 core 秒)
  - `rxn_discovery_790468f505`: screen_hei からの −73.0i と、path_hei の再開からの −73.2i が、どちらも `validate_ts` の直後に `ts_rejected:not_reaction_mode` になった。QRC の注記はなく、この case の新しいジョブは 0 件(QRC の opt の秒数も 0)。試行を使い切って unresolved_within_budget(`attempts_exhausted`)。
  - 新しいジョブ 15 本はすべて `rxn_discovery_f9bc3a1cd8` の分割の子 split1 の探索である(M2 の検証 (c) の再生。SCREEN の SP 9、NEB 1(失敗)、xTB freq 1、saddle 1、freq 1、QRC の opt 2)。
  - f9bc は multi_step になった。split1(CH3 + H2O → CH4 + HO)は −479.4i(χ 0.681)で elementary(δG_eff 16.02)、split2 は親の −1347.6i を再検証して elementary(41.25)である。W3_s6 の 460166732f の子(15.92、41.25)と一致する。
- 報告だけの再生(`--allow-new`。すべて CAP の前に完了)

| run | 壁時計 | 新しいジョブ | 結果 |
|---|---|---|---|
| s6_oh_ch4 | 8.2 分 | 15 本、1,871 core 秒 | 実計算と同じ経過。790468f505 は −73.0i と −73.2i を拒否して attempts_exhausted。f9bc は multi_step(split1 −478.3i、δG_eff 16.04) |
| W5_s6 | 10.7 分 | 15 本、2,165 core 秒 | 790468f505 は −73.0i と −58.6i(χ 0.023)を拒否して attempts_exhausted。W5 で 1,346 s かかった QRC の opt 2 本は走らない。f9bc は multi_step(split1 −480.5i、16.02) |
| s5_ch3_o2 | 18.5 分 | 19 本、4,331 core 秒 | 下記 |

- s5_ch3_o2
  - `rxn_discovery_b3febc8052`: 2 回目の試行の再開で得た −72.0i(χ 0.004)を拒否し、attempts_exhausted。QRC の側の極小と species、この鞍点の熱化学の記録は消えた。paths を b3febc8052 に限った CUR の再ステージ(`M3/real/s5_b3febc`、hits 21、misses 0)でも同じだった。
  - `rxn_discovery_b3e2f366b6`: multi_step。split1(CH3OO• → CH2OOH•)は −1866.6i を再検証して elementary(δG_eff 46.34)。
  - split2(CH2OOH• → CH2O + OH): −299.3i(χ 0.257)と −277.4i(χ 0.221)を拒否し、attempts_exhausted。2 つは同じ領域の鞍点である(E の差 0.002 mEh、C–O と O–O の差 0.004 Å 以内)。O–O 開裂の TS には届かなかった。M2 のように HCO• + H2O を経る分割はしない。
- CUR の 36 本すべてを `M3/replay/<name>` に張り替えた。M2 で M1 のままだった s5・W5_s6・s6_oh_ch4 も、M2 と M3 の完全な再生になった。
- M3 の CUR の分布(`M3/int/chi_CUR_after_M3.txt`): 接続した TS 28 件は χ 0.571〜1.000、順位のある接続反応で χ < 0.3 は 0 件。接続しない一次の鞍点 7 件は χ 0.002〜0.257(S6 の 4 件、S5 の −72i、split2 の 2 件)。
- r6・r7 の S6 の JobStore(読むだけ。W3 の manifest は古い schema なので chi.py では読めない)の、ν1 が −40〜−90 cm⁻¹ の DFT freq はすべて χ ≤ 0.023 だった(−77.1i、−73.0i、−70.4i、−60.6i、−58.6i、−57.3i。`M3/int/s6_soft_saddles.txt`)。

### 13.4 実証していないこと

- S6 の OH + CH4 → CH3 + H2O(790468f505)は、どの run でも unresolved のままである。引き抜きの TS には分割の子(split1、逆向き)から届いた。仮説そのものからの到達は G1-P1(M4)で扱う。
- S5 の CH2OOH• → CH2O + OH(split2)の TS は見つからなかった。
- 最も近い拒否(split2 の 0.221 と 0.257)と 0.3 の差は 0.04〜0.08 である。仮説と同じ結合を一部動かす鞍点が増えたら、分布を測り直す。

## 14. M4 ρ による初期 Hessian(X1 / G1-P1)

分析の §2.1 G1-P1、§3 X1 と §4.1 M4 にあたる。規則は [design.md](design.md) §6.2 の REFINE_SADDLE。

### 14.1 変えたこと

- `vibrations.shape_hessian(H, x, direction=None)`: 方向がなければ H₊(剛体運動を除き、固有値を max(\|λ\|, 1e-3) にしたもの)。方向 d があれば P·H₊·P − κd̂d̂ᵀ(d̂ は d の内部運動の成分を正規化したもの、P = I − d̂d̂ᵀ、κ = max(d̂ᵀH₊d̂, `KAPPA_MIN` = 0.05 Eh/bohr²))。「方向と最も平行な固有ベクトル」を選ぶ処理は削除した。
- `Ctx.direction`: TS の種は自分の虚モード、それ以外は ρ = ∇(Σ_切れる r − Σ_できる r)(ラベル付きの両端の `topology.bond_changes` を距離の `CoordinateTerm` にし、`declared_coordinate_gradient` で評価する)。結合変化がなければ宣言座標、最も変わる二面角、両端の差の順。
- 削除したもの: `_MIN_OVERLAP` と overlap テスト、接線を反応中心の行に制限する分岐、`profile.tangent`、`topology.reaction_centre`、`Seed.tangent`(TS の種だけが `Seed.mode` を持つ)。
- `_seed_hessian`: TS freq、なければ xTB freq(あれば常に)、なければ DFT freq。注記は `saddle_hessian:<ts_freq|xtb|dft>:<mode|rho|coordinate|chord>`。
- `NWChemSaddle` のジョブ鍵に `hessian_model: negative_along_mode` を足した。NWChem の saddle の入力(trust 0.1、sadstp 0.1、inhess 2、moddir 1)は変えていない。
  - 書く Hessian が変わるのに鍵が同じだと、TS の種(方向が以前と同じ)で古いモデルの saddle が JobStore から再生されてしまう。
  - そのため M4 以降の再生では、SADDLE を走らせるすべての case で新しい saddle ジョブが出る。`--strict` ではなく `--allow-new` で比べる。

### 14.2 QM なしの確認

- VAL7 の S6 の screen_hei の種(`rxn_discovery_790468f505`、ラベル付きの結合変化は +H2–O5 −C0–H2)と、その種の xTB Hessian(VAL7 の JobStore)を使った。
  - 新しい `shape_hessian` の出力は、プローブ `/home/user/hfauto_r8_probe/G1/rhoonly` の入力 Hessian と最大 1.4e-10 で一致した(最大要素 0.51)。
  - 負の固有値は 1 本で、−0.2187(κ = ρ̂ᵀH₊ρ̂ > 0.05)、その固有ベクトルと ρ の cos は 1.000000 だった。
  - VAL7 で使ったモデルは、負の固有値が −0.0028 で、その固有ベクトルと ρ の cos は 0.345 だった。
  - このプローブは 14 歩で引き抜きの TS の領域(C–H 1.184、O–H 1.382 Å)に着いた。着いた点は −490.0i / −81.9i の二次の鞍点だった(§2.1 の `G1/rhofreq`)ので、higher_order_retry を 1 回経ることが見込まれる。
- 単体テスト
  - 負の固有値はちょうど 1 本で、d̂ に沿って −κ(κ ≥ 0.05)。ほぼ縮退した 2 本の間の方向でも同じ。残りは P·H₊·P に一致し、剛体運動は零空間に入る。
  - 方向なしのモデルは正定値。
  - アセチルアセトンのプロトン移動では、ρ が 2 本の伸縮ベクトルの差(∇r(O2–H10) − ∇r(O6–H10))に 1e-6 で一致し、メチル回転子の成分は 0 だった。両端の差(以前の接線の代わり)は、メチル回転子の成分との cos が 0.96 である。
- golden G08(4 本の虚モード): 形を整えた Hessian の質量加重の虚モードと、指定したモードとの重なりは 0.979 と 0.991 だった。以前の 0.99 は固有ベクトルを反転する方式での値である。新しい方式が保証するのは、質量加重しない Hessian で負の方向が d̂ に一致することである。Cartesian のモードには質量加重しない剛体運動の成分が含まれる(G08 の 2 本では、モードとその内部運動の成分の cos が 0.91 と 0.96)。そのため試験はこちらを確かめる形に改めた。

### 14.3 実計算 S6(統合、`/home/user/hfauto_r9/M4/real/`、1 本ずつ、CAP 5,400 s)

- 種は 3 通りで、それぞれ rank 4 と rank 2(`/home/user/hfauto_r9/sites/wsl_local_r2.yaml`)で流した。
  - (a) CUR/s6_oh_ch4(VAL7 系)の複製を paths から再ステージ。
  - (b) CUR/W5_s6 の複製を paths から再ステージ。
  - (c) 新しい discover の run。
- 引き抜き `rxn_discovery_790468f505`(OH + CH4 → CH3 + H2O)の経過は 6 本とも同じだった。
  - screen_hei の種で `saddle_hessian:xtb:rho`。14 歩で二次の鞍点(−489.2i / −76.6i、χ 0.688)に着いた。G1 のプローブの予想どおりである。
  - その鞍点から higher_order_retry を 1 回(`saddle_hessian:ts_freq:mode`、22〜31 歩)行い、一次の TS を得て QRC で elementary になった。
  - 回転子の鞍点(−73i 前後、χ ≤ 0.023)には一度も着かず、それに QRC をかけることもなかった。

| run | 壁時計 | 最終 TS ν1 | ν2 | χ | E(Eh) | paths の core 秒 |
|---|---|---|---|---|---|---|
| (a) rank 4 | 30.2 分 | −486.4i | −37.5 | 0.682 | −116.031274 | 6,125 |
| (a) rank 2 | 69.5 分 | −489.5i | +24.1 | 0.681 | −116.031314 | 13,445(下記の束縛の不具合を含む) |
| (b) rank 4 | 27.7 分 | −486.0i | −38.1 | 0.682 | −116.031274 | 5,905 |
| (b) rank 2 | 38.5 分 | −495.4i | +17.9 | 0.683 | −116.031302 | 5,552 |
| (c) rank 4 | 62.0 分 | **−522.4i** | −34.6 | 0.691 | −116.031250 | 8,975(run 全体 13,559) |
| (c) rank 2 | 73.2 分 | **−515.2i** | −15.9 | 0.686 | −116.031284 | run 全体 10,259 |

- 合否
  - 6 本とも elementary、注記は `saddle_hessian:xtb:rho`、higher_order_retry は 1 回、χ ≥ 0.5、ν2 は −50i より上(TS の freq で確認)。
  - 虚振動の窓 −480〜−505i は (a)(b) の 4 本が満たし、(c) の 2 本(−522.4i、−515.2i)は外れた。窓は変えていない。
  - 6 つの TS は同じ鞍点である。C–H 1.183〜1.187 Å、O–H 1.376〜1.385 Å、C–H–O 174〜176°、E の幅は 6.4e-5 Eh(0.04 kcal/mol)。W3 の −481.1i(−116.031295)との差も 0.03 kcal/mol 以内である。ν1 の幅は、この平らな領域での収束点の違いによる。
  - 見かけの合格でないこと: W3 の TS を ts_calc として貸していない(6 本の TS は各 run の saddle ジョブの結果)。歩数だけでなく、各 TS の freq の ν2 と χ を見た。
  - ν2 が小さな虚数(−15.9〜−38.1i)の TS が 4 本ある。ゲート(saddle_cm1 = 50)の定義では一次の鞍点だが、メチル回転子の曲率はノイズの水準である。ν2 が正の 2 本のほうが 0.01〜0.04 kcal/mol 低い。停留性の判定(G5、should)の対象として記録する。
- 費用(rank 4 の (a)): 790468f505 の saddle 2 本で 41 歩 602 core 秒、TS の freq 2 本 506、QRC の opt 2 本 889。VAL7 と W5 では、この case は回転子の鞍点で unresolved だった(W5 は QRC に 1,346 s)。
- 分割の子 split1(CH3 + H2O → CH4 + OH)も `xtb:rho` で elementary(−476.7〜−491.7i、χ 0.673〜0.684)。
- W5 の path_hei の種(`string0_hei.xyz`)は、paths からの再ステージでは使われない。SCREEN がキャッシュから同じ screen_hei の種を返し、最初の試行で解けるからである。そこで同じ種から直接試した(`M4/probe_pathhei/`、NWChemSaddle と同じ入力)。
  - 種の xTB Hessian には虚振動がない(最低 138 cm⁻¹)。以前の overlap 検査では DFT Hessian に回り(W5 の `dft:overlap:0.00`)、maxiter のあと −58.6i の回転子の鞍点に着いていた。
  - ρ のモデル(固有値 −0.249 が 1 本、ρ との cos 1.000000)では 18 歩で収束し、C–H 1.183、O–H 1.385 Å、E −116.031245 Eh になった。NWChem の freq は −480.4i / −47.2i で、一次の鞍点である(ν2 は −50i の近く)。
  - 最初の試行は 3 歩目の勾配で NWChem が非零で終了した。同じ入力のやり直しで収束した(hfauto では nonzero_exit として ladder が扱う)。
- 不具合の修正(rank 2 で見つけた): ranks 2 の site(4 コア)では、QRC の両側が全 rank の 2 ジョブとして同時に走る。以前は「site の ranks より少ないときだけ `--bind-to none`」だったので、2 ジョブとも束縛されてコア 0–1 を取り合った。(a) rank 2 の QRC の opt は 1 歩あたり約 30 s(rank 4 では 6 s)、opt の合計は 11,026 core 秒だった。
  - NWChem は常に `--bind-to none` で走らせるようにした(`engine._argv`。単体テスト `test_nwchem_runs_unbound_at_any_rank_count`)。rank はジョブの鍵に入らないので、結果と再生は変わらない。4 rank の単独ジョブの速さへの影響は測っていない。
  - 修正は (b) rank 4 の開始後に入れた。(b) rank 2 以降は修正後のコードで、(b) rank 2 の opt は 3,110 core 秒(歩数は (a) rank 2 とほぼ同じ 208 歩)だった。

### 14.4 回帰と再生(`/home/user/hfauto_r9/M4/replay/`、`summary.txt`)

- strict(saddle ジョブのない 10 本: ch3n_spin、h2o_hf_inversion、nh3_planar_seed、s19_tma_hf2、s3_h2te、s5_undeclared、s7_nh3_icl、s8_so2_nme3、water_dimer_as、water_same_basin。structures から): 10 本すべて PASS(misses 0、記録が同一)。
- 報告(`--allow-new`、paths から再ステージ、24 本): saddle の鍵がモデル名を含むので、すべての saddle ジョブとその下流(TS の freq、QRC)が新しくなる。TS の比較は `M4/ts_compare.py`(記録ごとの outcome、ν1、E)、歩数は `M4/saddle_steps.py`。
  - 同じ TS(\|ΔE\| ≤ 1e-5 Eh、outcome も同じ): hcn、hono、s15、s4、s1、s2、s14、s16、s18、nh3_inversion、oxalic_two_step(子 2 つ)、s10、dme_c2v_seed、hcn_noscreen、hf_dimer_swap_noscreen、s12、s13、s17、s20 の 3 本。\|ΔE\| は最大 5.7e-7 Eh(s17)。
  - formaldehyde(hydroxymethylene_to_formaldehyde): ΔE = −2.5e-5 Eh(−2084.2i → −2084.8i)。M3 の TS は 0.8° ほど非平面の点で、収束判定の内側(Gmax 1.4e-4)で止まっていた。新しい TS は平面(Gmax 3e-5)で、0.016 kcal/mol 低い。同じ鞍点を、よりよく収束させたものである。
  - 歩数(TS のモードの種): 多くは同じか +1〜3 歩。oxalic の子は 7 → 12 と 14 → 23 歩に増えた。モードの種では、負の方向が xTB の Cartesian のモードそのものになり、H の固有ベクトル(cos 0.98 程度)ではなくなったためである。acac の新しいジョブは 7,806 core 秒(saddle、freq、QRC の再計算)。
  - s5_ch3_o2
    - split1(CH3OO• → CH2OOH•)は同じ TS(−1866.7i)。
    - split2(CH2OOH• → CH2O + OH): ρ の種からの鞍点は −314.1i、χ 0.281 で拒否された(M3 では 0.257 と 0.221)。そのあと string が barrierless と分類し、`barrierless_at_resolution` になった(M3 は attempts_exhausted)。SCREEN の SP が 1 点 SCF 未収束で unavailable になったのは、親の TS が新しい saddle ジョブになり、QRC の中間体の構造が少し変わった下流の差である。
    - b3febc8052(CH3 + O2 の会合): ρ の saddle 4 本がすべて maxiter で、unresolved のまま。会合の扱いは M8。
  - W3_s6
    - 460166732f_split1(引き抜きの逆向き)は同じ鞍点で ΔE = −4.8e-5 Eh(0.03 kcal/mol)、−485.8i → −487.0i、χ 0.680。
    - **790468f505 は elementary(W3 の −481.1i)から unresolved_within_budget に変わった。** screen_hei の ρ で −489.6i / −77.1i、retry で −499.0i / −54.8i。2 本とも higher_order で、試行を使い切った。2 本目の ν2 −54.8i はメチル回転子で、−50i のゲートのすぐ外である。上限とゲートは変えていない。
- CUR の 36 本は M4 の結果に張り替えた(strict と報告の 34 本は `M4/replay/<name>`、s6_oh_ch4 と W5_s6 は `M4/real/s6_val7_r4` と `M4/real/s6_w5_r4`)。W3 の −481.1i の記録は `M3/replay/W3_s6` に残る。

### 14.5 実証していないこと

- (c) の 2 本の最終 TS の虚振動(−522.4i、−515.2i)は、窓 −480〜−505i の外である(同じ鞍点であることは構造とエネルギーで確かめた)。
- W3_s6 の再生では、790468f505 が retry のあとも二次の鞍点(ν2 −54.8i)で、未解決になった。7 つの文脈のうち 6 つで解けた。
- 回転子の曲率がノイズの水準の TS(ν2 −15.9〜−47.2i)を一次と数えているのはゲートの定義による。停留性の判定は G5(should)。
- S5 の O–O 開裂の TS は今回も見つからなかった。拒否された鞍点の χ(0.281)は 0.3 に近づいている。

## 15. M5 行 11 の削除(X5 / G3-P1)

分析の §2.3 G3-P1、§3 X3-4 と X5、§4.1 M5 にあたる。規則は [design.md](design.md) §3 の `is_first_order_saddle` と §6.2 の決定表。

### 15.1 変えたこと

- 行 11(崩壊した saddle と、QRC が失敗した柔らかい TS を VALIDATE_INTERMEDIATE にかける `saddle_collapsed`)、`_soft_ts_failed`、`ts_check` の `collapsed` を削除した。決定表は 16 行になり、旧行 12〜17 は 11〜16 になった。§4 などの前の節の行番号は、当時の番号のままである。
- 虚モードのない点(`ts_rejected:no_imaginary_mode`)は、ほかの拒否と同じく失敗した試行になる。試行には数えたままで、QRC はしない。
- 柔らかい TS(最低モードが −saddle_cm1〜−noise_cm1)は TS として受理する。QRC が失敗したときは、ほかの TS と同じく行 8 で扱う。両側が同じ basin か同じキー(`same_state`)で試行が残っていれば探索を続け、残っていなければ `connection_failed` にする。`same_state` は振幅を広げない(M2 のまま)。
- VALIDATE_INTERMEDIATE は、プロファイルの井戸(行 12、`path_intermediate`)だけを扱う。中間体を出すのは、行 12 と QRC の側(GEN-05)だけである。
- ゲートのしきい値(noise_cm1、saddle_cm1)は変えていない。

### 15.2 到達性と検証

- 消した分岐は、実計算では一度も発火していない。分析では r5〜r7 のすべての log で `ts_rejected:no_imaginary_mode` が 0 件、検証済みの TS で \|ν\| < 50 cm⁻¹ のものも 0 件だった。CUR(M4 の 51 case の log と 34 の SaddleClaim)でも、`saddle_collapsed` と `no_imaginary_mode` はともに 0 件、ν1 > −50 cm⁻¹ の SaddleClaim も 0 件だった(`/home/user/hfauto_r9/M5/audit/row11_reach.sh`)。したがって、VAL7 の記録は変わらないはずである。
- 単体テスト: 虚モードのない saddle は QRC も緩和もせず、SCREEN(行 11)、FIND_PATH(行 14)、次の種(行 13)に進む。柔らかい TS は、硬い TS と同じ行 8 の判断になる(`test_reaction_case_actions.py`、`test_reaction_case_decide.py`)。中間体の原子の対応の試験は、崩壊した saddle の代わりにプロファイルの井戸で行う(`test_reaction_case_connection.py`)。
- 再生(`tools/replay.py M5M7 --src CUR`、36 本を最初の stage から、壁時計 4.5 分、`/home/user/hfauto_r9/M5M7/replay/`): 全 run で misses 0(hits 956)。paths の記録の差は、削除した `chiral` 欄(§17.4)のほかに 0 件で、case log の理由と行動の列も W3_s6、W5_s6、s6_oh_ch4、s5_ch3_o2、oxalic_two_step で CUR と同一だった。GEN-05(`end*_to_new_basin`)は W3_s6(3 件)、oxalic(1)、s5_ch3_o2(1)、s6_oh_ch4 と W5_s6(各 1)で前と同じく発火し、行 12(`path_intermediate`)は s5_ch3_o2 の 1 件で発火した。W3_s6 の f9bc の split2 での行 12 は M4 の再生ですでに消えていて(CUR でも 0 件)、この波の差ではない。

## 16. M6 失われた状態の分解能(X3-1 / G4-P1)

分析の §2.4 G4-P1、§3 X3-1、§4.1 M6 にあたる。規則は [design.md](design.md) §7.1 の explore と「結合と状態」。

### 16.1 変えたこと

- `topology.resolved_bond_changes(symbols, x, y)` を足した(S3 で `bond_changes` に統合)。r − r_thr が一方の構造で ≥ +0.1 Å、他方で ≤ −0.1 Å の対だけを、形成と切断に分けて返す。`RESOLVED_A` = 0.1 Å はモジュール定数で、根拠は S19 の N···H の GFN2 と PBE0 の差(0.15 Å)である。
- explore の R6(`_relaxations`)は、ラベルが screen 極小に残らない seed のうち、seed と落ちた basin の構造(`identity.basin_coords` で seed の添字に並べたもの)の間に resolved な変化が 1 つでもあるものだけを、失われた状態とする。すべての変化が resolved であることは求めない。
- `bond_changes` と状態ラベルは変えていない(帯を広く使う G8-P5 は could)。

### 16.2 QM なしの確認(CUR の explore の入力に、作業ツリーの `_relaxations` をかけた。`/home/user/hfauto_r9/M6/audit/r6_lost.py`)

| run | seed → 落ちた basin | ラベルの違う結合 | r − r_thr(seed / basin、Å) | 判定 |
|---|---|---|---|---|
| s19_tma_hf2 | c01、c02(C3H9N+FH+FH)→ neutral(C3H10FN+FH) | N1–H13 | c01 +0.006 / −0.006、c02 +0.043 / −0.006 | 帯の中。`relax_tma_hf2_c01` は作られない |
| s5_ch3_o2 | p00(CH3+O2)→ CH3O2 | C0–O5 | +2.096 / −0.409 | 残る |
| s6_oh_ch4、W3_s6、W5_s6 | p02、p05(CH4+HO)→ CH3+H2O | H2–O5 と C0–H2(p05 は H4 で同じ形) | +2.004 / −0.350、−0.376 / +0.104 | 残る(代表は p02) |

- s5_undeclared、s7、s8、s10 には失われた状態がない(変わらない)。
- 見かけの合格でないこと: S6 の seed を残しているのは H2–O5 の変化(帯の端から 0.25 Å 以上)である。帯のすぐ外の C0–H2(+0.104)だけで残っているのではない。δ は分析の値のままで、run を見て調整していない。

### 16.3 S19 の再生

- 期待される差は次のとおり。explore の `relax_tma_hf2_c01` と、dft の species `spc_relax_tma_hf2_c01` がなくなり、`min_neutral_ff60c312` の members から seed が外れる。DFT の seed の opt 1 本が走らなくなる(VAL7 では 45 歩、2,351.6 s × 4 rank = 9,406 core 秒で、既知の basin に落ちて freq なし)。dft の hits は 8 から 7 になる。宣言反応(same_basin)と ΔG_assoc(−11.53)は変わらない。S5 と S6 の R6 の仮説は残る。
- 結果(同じ再生): s19 の差は期待どおりの 3 件だけである。explore の `relax_tma_hf2_c01` と dft の `species_spc_relax_tma_hf2_c01` がなくなり、`min_neutral_ff60c312` の members から seed が外れた。dft の hits は 8 から 7 になり、misses は 0。CUR でもこの seed は既知の basin に落ちて仮説を作っていなかったので、paths の反応の集合、ΔE‡ と ΔE_rxn は変わらない(ΔG_assoc −11.53 は 1e-9 kcal/mol 内)。s5_ch3_o2 の `relax_ch3_o2_p00`、s6_oh_ch4・W3_s6・W5_s6 の `relax_oh_ch4_p02` は残り、それらの記録とジョブは CUR と同一(misses 0)。

## 17. M7 対称性の判定を 1 つに(X3-2 / G5-P3)

分析の §2.5 G5-P3、§3 X3-2、§4.1 M7 にあたる。規則は [design.md](design.md) §7.1 の thermo と「対称性」。

### 17.1 変えたこと

- `chemistry/symmetry.py` の `analyze(symbols, coords, hessian)` を足した。thermo は、σ、m、直線性、振動と慣性モーメントを評価する構造を、すべてここから取る。
- 削除したもの: `thermo.symmetry_number`(libmsym の equivalence 2e-3。例外なら σ = 1)、熱化学での `identity.is_chiral` の利用、`MinimumRecord.chiral`、回転定数の床による直線の判定。以前は 4 つの検出器(外部自由度の階数、回転定数の床、libmsym、is_chiral)が別々のしきい値で判定していた。S18 の H···H2 では、libmsym が C∞v、前の 2 つが非直線と食い違っていた。
- `_EXTERNAL_RANK_TOL` は、点群のない場面(エンジンが自分の freq を解析するとき、試行の方向、整えた Hessian)だけに残した。
- 受理の基準は identity の basin の基準(5e-5 Eh、0.05 Å)のままで、系ごとの ε はない。分析が退けた受理のしきい値 1e-6 Eh は使っていない(acac のエノールは Cs、m = 1)。
- libmsym の緩い 2 つの設定は、7 つのしきい値をすべて同じ値にする。zero と orthogonalization を既定のまま残した途中の版では、1e-3 Å の揺らぎで 12 subject の点群が部分群か C1 に落ちた(NH3 と SN2 の D3h の TS、シュウ酸の C2h、CH3 の D3h、H3 の D∞h の TS など)。libmsym が 3 つの設定すべてで部分群を返すか、`no primary axis for reorientation` で例外になったためである。

### 17.2 symcheck(QM なし、`/home/user/hfauto_r9/tools/symcheck.py`、CUR の 110 subject。出力は `/home/user/hfauto_r9/M7/symcheck_cur.txt`)

- (a) 揺らぎ: 原子あたり RMS 1e-3 Å の正規乱数を 20 回加えても、110 subject すべてで(点群、σ、直線、m)が変わらない。
- (b) 既知の点群: NH3 C3v(σ 3)、H2O2 gauche C2(m 2)、HCN と HNC C∞v、H···H2 C∞v、TMA C3v、acac のエノール Cs(m 1)の 6 つとも一致した。VAL7 の run(`symcheck_val7.txt`、92 subject)でも (a)(b) は同じく通る。
- 旧検出器からの点群の移り(`/home/user/hfauto_r9/M7/thermo/groupmoves.txt`。run と subject の組で 109 件)
  - libmsym が例外を出して σ = 1 だった 37 件は、C1 12 件、Cs 21 件、C3v 4 件になった。σ・m・直線性が変わったのは C3v の 4 件だけである。C1 の 12 件は、旧 `is_chiral` でもキラル(m = 2)だった。
  - C3 と判定されていた 2 件は C3v になった。σ は 3 のままで、m も 1 のまま(旧 `MinimumRecord.chiral` は False)である。C3 のままなら m = 2 になるはずで、判定が食い違っていた。
  - σ・m・直線性が変わった subject は、合わせて 6 件である(下の表の原因)。ほかは、点群も 3 つの量も変わらない。

### 17.3 δG_eff の変化(QM なし)

作業ツリーの thermo stage を、CUR の各 thermo stage の入力 view にかけ直した(ジョブなし、run には何も書かない。`/home/user/hfauto_r9/M7/thermo/thermo_recompute.py`、出力は同じ場所の `thermo_cur.txt`)。値は kcal/mol で、ここに挙げたもの以外の δG_eff は 1e-4 kcal/mol 以内で変わらない。

| run | 反応 | δG_eff(CUR → M7) | 原因の subject と点群 |
|---|---|---|---|
| s18_h3_doublet | h_exchange | 4.0911 → 6.3656 | 極小 H···H2(軸から 0.006 Å): C∞v なのに非直線(振動 3 本)として扱っていた → 直線(4 本)。G −2.27 |
| s17_hf_dimer_swap | hf_dimer_swap | 0.9848 → 1.3954 | TS: Cs(σ 1)→ C2h(σ 2)。G +0.41 |
| s6_oh_ch4、W5_s6 | 790468f505(OH + CH4 → CH3 + H2O) | 1.4993 → 0.8483、1.5041 → 0.8531 | 反応物の CH4···HO 錯体(OH の H が CH4 の面を向く): libmsym の例外で σ 1 → C3v(σ 3)。G +0.651 = RT ln 3。ΔG_assoc も +0.651 |
| s6_oh_ch4、W5_s6 | f9bc3a1cd8_split2 | 41.2498 → 40.5987 | 同じ錯体が反応物 |
| W3_s6 | 460166732f_split2 | 41.2498 → 40.5987 | 同じ錯体(W3_s6 では C3v の basin が 2 つ) |
| W3_s6 | f9bc3a1cd8 | 48.48002 → 48.48015 | TS: 例外(σ 1)→ Cs。σ と m は同じで、対称化した構造で振動を評価した差(G +2.1e-7 Eh) |

- δG_eff を持たない量の変化: W3_s6 の 790468f505(未解決)の ΔG_rxn −0.651 と ΔG_assoc +0.651、460166732f_split1 の ΔG_rxn +0.651(生成物が同じ錯体)。
- ほかの species の G の変化は、対称化した構造で振動を評価した分だけで、最大 2.1e-7 Eh である。
- S17 は M4 で生じた差を戻すものである。VAL7 と M3 の TS は旧検出器でも C2h(1.3919)だったが、M4 の再生で saddle を取り直した TS(\|ΔE\| 5.7e-7 Eh)を、libmsym の 2e-3 が Cs と判定していた。§14.4 ではこの δG_eff の変化を記録していなかった。M7 では、どちらの TS も C2h になる。
- S18 の一致(検証 (c)): 正確に直線の極小を持つ `/home/user/hfauto_r7/W2/superseded_k1neg/s18_h3_doublet` は、M7 でも 6.3596(変化 1e-4 未満、`thermo_w2_s18.txt`)である。CUR(VAL7 の run)の 6.3656 との差は、M7 の前の 2.27 から 0.006 になった。
- 見かけの合格でないこと: S18 だけを見ていない(全 subject の揺らぎ、既知の点群、全反応の δG_eff の変化)。受理の基準は系ごとに変えていない。

### 17.4 再生

- `MinimumRecord.chiral` を削除したので、再生ではすべての minimum の記録にこの欄の差が出る(結論の値ではない)。r9-M7 より前の run は、途中の stage からは再開できない([design.md](design.md) §8)。
- 結果(`/home/user/hfauto_r9/M5M7/replay/`、分類は `/home/user/hfauto_r9/M5M7/audit/classify.py`、出力 `classify.txt`): 36 本すべてで misses 0。差 732 件はすべて次のどれかで、説明のつかない差は 0 件である。(1) `chiral` 欄の削除 108 件。(2) §17.3 の表に挙げた変化 33 件。17.3 の再計算の値と 1e-9 の相対誤差で一致し、表のすべての変化が再生でも起きた(S18 6.3656、S17 1.3954、S6 0.8483 / 0.8531、split2 40.5987 など)。同じ記録の H・ZPE・帯の変化が 53 件。(3) 対称化した構造で振動を評価した分の差 352 件(1e-4 kcal/mol 以内か、species で 2.1e-7 Eh 以内。最大は S6 の Cs の TS の H と ZPE の 1.2e-7 Eh)。(4) report の行と表の sha256 が thermo の差を写したもの 183 件。(5) §16.3 の S19 の 3 件。thermo と report のほかの stage(structures〜paths)の差は (1) と (5) だけである。
- S18 の検証 (c): W2 の run(`superseded_k1neg`)の写しを最初の stage から再生すると、今のジョブ鍵では極小の opt をやり直すことになり(新しいジョブ 16 本、97 core 秒)、VAL7 と同じく軸から外れた極小(H の軸からの距離 0.0097 Å)に収束した。そのため 2 つの再生の δG_eff は 6.365608 と 6.365609 で一致するが、これは構造の違いを試していない。正確に直線の W2 の構造そのものに作業ツリーの thermo をかけると 6.3596 のまま(`/home/user/hfauto_r9/M5M7/audit/thermo_w2_s18.json`、変化なし)で、CUR の再生の 6.3656 との差は 0.006 kcal/mol(M7 の前は 2.27)である。
- symcheck を再生した 36 本にかけ直した(`/home/user/hfauto_r9/M5M7/audit/symcheck.txt`): (a) 110/110 が不変、(b) 6 つの既知の点群がすべて一致。

## 18. M8 低スピン結合の組成と会合のスキャン(X4-1/2 / G2-P2、G2-P4、G2-P8)

分析の §2.2 G2、§3 X4、§4.1 M8 にあたる。規則は [design.md](design.md) §7.1 の「会合と低スピン結合の組成」と、§6.2 の SCREEN・FIND_PATH の行(M8 の時点の行 11・14)。

### 18.1 変えたこと

- 低スピン結合の組成(`electronic_state.low_spin_coupled`): 開殻の成分が 2 つ以上あり、宣言した多重度が高スピンの結合 Σ(m_i − 1) + 1 より小さい組成。拒否しない。S5 の CH3• + O2 の二重項がこれにあたる。
  - このクラスの一重項は `low_spin_singlet_unsupported` で止める。旧 `open_shell_singlet_unsupported` の名前を一般化しただけで、止まる組成は変わらない(開殻の成分が 1 つなら、スピン結合で一重項は作れない)。
- スピン交差(G2-P8): 多重度だけが違う宣言反応の理由を `spin_crossing_reaction_unsupported` にした。再生では ch3n_spin の失敗理由の文字列だけが変わる。
- 会合(G2-P4)
  - 仮説: 形成 1 本(2 つの断片の間)で切断なし、単量体の DFT 極小が錯体と同じ LOT にある仮説は、反応物側を分離した単量体にする(`ReactionRecord.monomers`)。出所(宣言、explore、R6)には依らない。
  - 経路: SCREEN が、付加体の結合 + 1.5 Å から付加体までの 8 点の緩和スキャンを走らせる(NWChem の `zcoord` の固定結合、前の点の vectors から SCF)。string は走らせない。
  - 判定: [Σ E(単量体)、スキャン、E(付加体)] を `barrier_verdict` で分ける。
  - 熱化学: 反応物側は分離した単量体の状態の G の和で、錯体は対象にしない。標準状態の換算は化学量論から入り(Δn = −1)、barrierless の δG_eff は max(ΔG_rxn, 0)。
- しきい値と予算(resolution_kcal、spin_tol、saddle の試行)は変えていない。

### 18.2 S5 の期待値の訂正

| case | VAL7 | M2 | M3 | M4 | M8 の期待 |
|---|---|---|---|---|---|
| `rxn_discovery_b3febc8052`(CH3···O2 → CH3OO•、R6 の relaxation) | reassigned_step(−72.0i、`spin_contaminated` で順位なし) | unresolved(`connection_failed`。両側が同じ CH3 + O2 の状態) | unresolved(`attempts_exhausted`。−72.0i を χ 0.004 で拒否) | unresolved(ρ の saddle 4 本がすべて maxiter) | 会合。barrierless_at_resolution で順位あり、ΔE_e 約 −37 kcal/mol |
| `rxn_discovery_b3e2f366b6`(CH3OO• → CH2OOH•) | reassigned_step(32.67、順位あり) | multi_step(CH2OOH• を経る) | split1 elementary(46.34)、split2 attempts_exhausted | split2 barrierless_at_resolution(string) | 変わらない |

b3febc8052 の期待を変える理由:
- 問いの形が誤っていた。障壁のない会合では、反応物は漸近であって極小ではない。これまでの case は、R6 の seed が DFT で残した CH3···O2(C···O 2.84 Å)を反応物の極小にしていた。
- その極小の ⟨S²⟩ 1.71 は、分離した低スピンの二重項の BS 行列式から必ず出る値(1.75 = 二重項 2/3 と四重項 1/3 の混合)である。斥力的な四重項が 1/3 混ざるので、極小そのものが汚染による人工物である。AP の面では、この点から C–O の方向に下る(分析 §2.2)。
- M2〜M4 の unresolved は、人工の極小どうしを結ぶ鞍点(両側が同じ CH3 + O2 の状態、χ 0.004)を正しく棄却した結果で、会合の結論ではない。
- M8 では、反応物の端は分離した CH3(⟨S²⟩ 0.754)と O2(2.010)で、付加体 CH3OO•(0.755)とともに汚染がない。プローブの C–O 緩和スキャン(`/home/user/hfauto_r8_probe/G2-1/scan/scan.log`)では、BS の曲線の山は +0.13 kcal/mol で解像度 1.0 より小さく、1.45 Å で −36.8 kcal/mol(AP −37.0)だった。
- 合格条件: barrierless_at_resolution で rankable。ΔE_e は約 −37 kcal/mol。ZPE を含めた ΔH₀ が、実験の D₀(CH3–OO)約 32 kcal/mol から数 kcal/mol 以内。b3e2f366b6 が残ること(結合した異性体どうしの反応で、結合が切れるので会合の形ではない)。
- 見かけの合格の見分け方: 多重度を上げた宣言(四重項)にすり替えていないか。スキャンの開始を Coulson–Fischer 点(約 2.2 Å)より内側に置いていないか。resolution を広げていないか。点を落としていないか。

S5 全体の結論は、§3.4 の「結論なし(多参照性のある二重項の結合)」から「会合は障壁なし(UKS、AP なし)」に変わる見込みである。結論が出なかった主因は、多参照性ではなく問いの形だった(分析 §0.1)。多参照法による定量的なエネルギーと会合の速度定数は、引き続き範囲外である。

### 18.3 検証系

| 系 | 入力 | 確かめること | 期待 | 見かけの合格の見分け方 |
|---|---|---|---|---|
| H + C2H4(二重項) | `configs/systems/h_c2h4.yaml`(単量体、錯体 H···C2H4 2.9 Å、組成。宣言反応 錯体 → C2H5) | 障壁のある会合(反例)。開殻の成分は 1 つなので低スピン結合ではない | スキャンに山があり、saddle に進む。参照の障壁は 1.72 kcal/mol(W1、NHTBH38) | PBE0 の山が解像度より低ければ barrierless として記録し、resolution を変えない |
| BH3 + NH3(一重項、閉殻) | `configs/systems/bh3_nh3.yaml`(単量体、錯体 B···N 3.2 Å、組成。宣言反応 錯体 → H3B–NH3) | 閉殻の障壁のない会合 | barrierless_at_resolution で rankable。ΔG_rxn は分離した単量体基準(付加体の会合の自由エネルギー) | 錯体の極小を反応物にしていないか |
| CH3 + O2(四重項) | `/home/user/hfauto_r9/inputs/ch3_o2_quartet.yaml`(repo の外、discover) | 高スピンの宣言。should の X4-3(AP)で使う | 低スピン結合ではない。四重項の面は斥力的(プローブで 2.84 Å で +0.83、1.45 Å で +115.9 kcal/mol)で、付加体はない | AP がかかっていないか |
| ch3n_spin(再生) | `/home/user/hfauto_r7/VAL7/outside/ch3n_spin.yaml` | スピン交差の理由 | 両端の species が `spin_crossing_reaction_unsupported` で失敗し、ジョブ 0 | — |

- H + C2H4 と BH3 + NH3 は known_endpoints で流す(discover は要らない)。組成と単量体を持つので、discover でも流せる。
- 宣言した錯体は、端点の原子順と会合の形を決める入力である。会合の case は錯体の極小をエネルギーと熱化学の端点にしない。
- BH3 + NH3 の錯体(孤立電子対を B に向けた配置)は、DFT で付加体に落ちうる(H + C2H4 も、PBE0 に障壁がなければ同じ)。そのときは錯体の入力構造で会合の形を判定し、case の反応物端もその構造にする(統合で追加。行 2 の same_basin は会合には効かない)。錯体の配置や resolution を変えて通すことはしない。

### 18.4 実計算と再生(統合、`/home/user/hfauto_r9/M8/`)

- QM なしの監査(`M8/audit/assoc_audit.py`): 作業ツリーの `hypotheses.select` を CUR の 33 本の paths の入力にかけると、会合の仮説になるのは s5_ch3_o2 の b3febc8052 だけだった(単量体 `min_ch3_5a7804cb`、`min_o2_2fbb137e`)。ほかの仮説の数と id は変わらない。ほかの run のジョブ鍵も変わらない(`fixed_bond` と `scf_guess` は渡したときだけ鍵に入る)。
- 統合で足したこと: 錯体が自分の DFT 極小を持たず付加体に落ちた宣言反応は、錯体の入力構造で会合の形を判定する(`minima[0]` は付加体の極小、case の反応物端はその入力構造、行 2 は会合に効かない)。BH3 + NH3 の錯体は実際に付加体へ落ちた(付加体の members に錯体)。これがなければ、この系は same_basin で終わっていた。会合は ΔG_assoc を出さない(ΔG_rxn が分離した単量体基準の値)。

S5(`M8/real/s5_ch3_o2`、CUR の写しを paths から、CAP 2,700 s。壁時計 2 分 55 秒、新しいジョブはスキャンの opt 7 本で 172 s、687 core 秒。VAL7 の S5 は 54 分):

| r(C–O) Å | 分離 | 2.930 | 2.716 | 2.502 | 2.287 | 2.073 | 1.859 | 1.644 | 1.430(付加体) |
|---|---|---|---|---|---|---|---|---|---|
| ΔE kcal/mol | 0 | −1.04 | −1.04 | −0.94 | −1.15 | −3.70 | −13.98 | −28.64 | −37.56 |
| ⟨S²⟩ | 0.754 / 2.010 | 1.723 | 1.687 | 1.614 | 1.462 | 1.136 | 0.772 | 0.755 | 0.755 |

- b3febc8052: `scan:barrierless`、barrierless_at_resolution、順位あり(blocker なし)。2.716 → 2.502 Å の山は +0.10 kcal/mol で解像度 1.0 より小さい(プローブの +0.13 と同じ)。⟨S²⟩ は Coulson–Fischer 領域(2.3 → 1.9 Å)で 1.46 → 0.77 と単調に変わり、枝の跳びはない。
- ΔE_e = −37.56 kcal/mol(合格条件 −37 ± 2)。ΔZPE +5.93 で ΔH₀ = −31.63 kcal/mol、実験の D₀(CH3–OO)約 32 から 0.4 以内。ΔH298 = −33.20、ΔG298(1 atm)= −23.86。δG_eff は 0.0。
- 順位のあるスピン汚染の reassigned_step はない(CUR の b3febc8052 は unresolved で `spin_contaminated`)。
- b3e2f366b6 は CUR(M4〜M7)と同じ: multi_step、split1 elementary 46.34(帯 46.28〜46.37)、split2 barrierless_at_resolution(string)。VAL7 の reassigned 32.67 は M2 の状態単位の判定で multi_step に変わったもので(§12)、M8 では変わらない。
- CUR との差(paths・thermo・report)は b3febc8052 の記録と、それを写した report の行の並びだけである。

反例と閉殻(known_endpoints、新しい run。CAP 2,700 s):

| 系 | 壁時計 | 錯体 | スキャン(分離 = 0、長い側から、kcal/mol) | 結果 |
|---|---|---|---|---|
| H + C2H4(二重項、C–H) | 6 分 55 秒 | DFT 極小あり(−0.44) | −0.06、+0.30、+0.60(2.18 Å)、0.00、−4.23、−16.65、−34.11、付加体 −44.74。⟨S²⟩ 0.752〜0.803 | `scan:barrierless`、barrierless_at_resolution、順位あり。ΔG_rxn −34.89 |
| BH3 + NH3(一重項、B–N) | 4 分 23 秒 | 付加体に落ちた | −5.66、−7.63、−10.42、−14.29、−19.50、−26.13、−33.44、付加体 −37.77 | `scan:barrierless`、barrierless_at_resolution、順位あり。ΔG_rxn −23.56 |

- H + C2H4 は期待(約 2 kcal/mol の山から saddle へ)を満たさない。PBE0-D3BJ/def2-SVPD の山は +0.60 kcal/mol(谷からの深さ 0.66)で、解像度 1.0 より小さい。参照の障壁 1.72(W1)を混成汎関数が過小評価するのは既知で、resolution とスキャンの点は変えていない。山から REFINE_SADDLE に進む分岐は単体テスト(`test_reaction_case_paths.py`)だけで、実計算では通っていない。
- BH3 + NH3 の ΔE_e −37.8 は CCSD(T)/CBS の約 −31 より深い(BSSE を補正しない def2-SVPD。[design.md](design.md) の ΔG_assoc の注意と同じ)。
- 四重項の CH3·O2(`/home/user/hfauto_r9/inputs/ch3_o2_quartet.yaml`)は should の AP(G2-P3)用で、今回は流していない。

再生(`M8F/replay/summary.txt`):
- s5_ch3_o2 を除く 35 本を CUR から strict: 35 本すべて PASS(misses 0、記録の差 0。ch3n_spin の理由の文字列は比較の対象外)。M8 の WP のコードでの先行の再生(`M8/replay`)も 35/35 PASS。
- s5_ch3_o2 を M8 の実 run から strict: PASS(新しいジョブ 0、記録の差 0)。CUR は 36 本とも `M8F/replay` を指す。
- 再生ツールは、元の記録にない欄が空の既定値(`monomers: []`)で現れる差を、結論の差として数えないようにした(`tools/replay.py` の `diff`、`test_tools.py`)。

## 19. M9 既定の順位のエネルギー層(X7-1 / G6-P1)

分析の §2.6 G6、§3 X7、§4.1 M9 にあたる。規則は [design.md](design.md) §7.3・§7.4。

### 19.1 変えたこと

- `configs/methods/m06-2x-d3_def2-tzvpd.yaml` を足した(m06-2x、def2-tzvpd、d3zero、grid fine、scf_energy_tol 1e-7)。NWChem 7.2.3 は M06-2X のゼロ減衰の D3(s6 1、sr6 1.619、s8 0)をかける(golden G30)。
- discover と known_endpoints: paths の後に `{id: sp, stage: sp, engine: nwchem, methods: [m06-2x-d3_def2-tzvpd]}` を置き、thermo に `energy_method: m06-2x-d3_def2-tzvpd` を付けた。コードは変えていない(sp stage と `energy_method` は既存)。
- method_panel から panel_thermo を削除した(panel_sp → panel_report)。順位の層の定義は既定の pipeline の thermo の 1 つだけで、パネルは表示だけになる。
- しきい値、窓、予算は変えていない。

### 19.2 根拠(プローブ `/home/user/hfauto_r8_probe/G6-1_*`、構造は hfauto の PBE0/SVPD の停留点)

- 4 反応の ΔE‡ と CCSD(T)/def2-TZVPD との差は [design.md](design.md) §7.4 の表。平均絶対誤差は PBE0/SVPD 2.68、M06-2X/TZVPD 1.03 kcal/mol。
- S6 の分離基準の ΔE_rxn: CCSD(T) −12.77、M06-2X −14.54(−1.8)、PBE0/SVPD −15.02(−2.2。W3 の freq のエネルギー)。
- M06-2X の grid: fine と xfine で、S6 の障壁は 5.817 と 5.813 kcal/mol(`G6-1_m062x_grid`)。
- golden G30(`tests/golden/test_nwchem_output.py`): S6 の TS の M06-2X の SP の出力は、method ファイルと Level が一致し(m06-2x、d3zero、def2-tzvpd、fine、1e-7)、PBE0/SVPD の freq の Level とは `same_pes` で method・basis・dispersion が違う。

### 19.3 合格条件(統合の実計算)

- VAL7 の全系を新しい既定で流し、rankable な行の `energy_level` がすべて `m06-2x-d3zero/def2-tzvpd` になること。`energy_layer_missing` で順位から外れた反応は、SP の失敗の種類とともに記録する。
- 再生: paths までの stage は misses 0 で記録が同一。新しいジョブは sp の M06-2X の SP だけで、変わる記録は sp・thermo・report だけ。ほかの engine や kind のジョブが出たら不合格。
- S6(CUR の s6_oh_ch4、W5_s6)の分離基準の ΔE‡ が 5〜7 kcal/mol(プローブ 5.82、CCSD(T) 6.96)。
- SN2(Cl⁻ + CH3Cl)の δG_eff が、CCSD(T) の SP を `energy_method` にした composite から ±1 kcal/mol 以内。
- HONO の M06-2X を測り、ΔE‡ と ΔE_rxn を CCSD(T) と比べる。
- R1〜R4 の回帰で順位が変わった箇所は、化学的な理由(PBE0 の非局在化誤差が大きい反応クラスほど δG_eff が上がる)で説明できるものだけを受け入れる。

### 19.4 見かけの合格の見分け方

- `energy_level` の列名だけを見ていないか。各 `SpeciesThermo` の `energy_calc` が M06-2X の sp を指し、その Level が freq の Level と `same_pes` で method が違うこと(`tests/integration/test_pipelines_e2e.py`)。
- モックのエネルギーのテストだけで通していないか。fake の sp は method に依らず PES のエネルギーを返すので、値は golden G30 と実計算で確かめる。
- PBE0 の値を層として再利用していないか。SP がない対象は G = None(`energy_layer_missing`)で、freq のエネルギーで代用しない。
- 窓や resolution を層に合わせて変えていないか。

### 19.5 結果(統合、`/home/user/hfauto_r9/M9/`)

- **再生**(`tools/replay.py M9 --src CUR --allow-new`、36 本、壁時計 54 分): paths までの全 stage で misses 0、記録の差は sp・thermo・report(パネルは panel_*)だけ。新しいジョブは 105 件すべて sp の nwchem `energy`(M06-2X)で、計 11,860 core 秒(多くは 1 点 10〜60 s、シュウ酸 8 原子 約 140 s、acac 15 原子 337 基底関数 約 600 s)。SP の失敗は 0 件で、`energy_layer_missing` で順位から外れた反応はない。rankable な行(29 本で 38 行)の `energy_level` はすべて `m06-2x-d3zero/def2-tzvpd`、各 `SpeciesThermo` の `energy_calc` は M06-2X の sp で、freq と `same_pes` の method が違う(QM なしの照合)。順位の行がない run(水、H2Te、TMA·(HF)₂ の same_basin)は SP の対象がなく、same_basin の行の blockers に `thermo_unavailable` が足される(順位には元々ない)。パネルの pipeline は R9 の写し(`/home/user/hfauto_r9/pipelines/`、panel_thermo なし。S20 SN2 の 400 K・1 M の report は M06-2X 層のまま検証用の thermo で出し、report.html と ranking.csv は 14.89 で一致)。
- **自己一致**: CUR を M9 の再生に向け直し、`tools/replay.py M9S --src CUR --strict` で 36/36 PASS(新しいジョブ 0、記録の差 0)。
- **精度**(CCSD(T)/def2-TZVPD との差、kcal/mol。S6 は W5_s6 と s6_oh_ch4、パネルは `s20_*` と HONO の写し `M9/real/hono_panel`):

  | | ΔE‡ M06-2X | ΔE‡ PBE0/SVPD | ΔE_rxn M06-2X | ΔE_rxn PBE0/SVPD |
  |---|---|---|---|---|
  | HCN → HNC | 46.02(−1.79) | 46.62(−1.19) | 12.92(−2.30) | 13.18(−2.04) |
  | CH3O• → CH2OH• | 32.66(−0.71) | 33.01(−0.36) | −8.41(−0.61) | −7.75(+0.05) |
  | HONO trans → cis | 11.64(−0.01) | 13.67(+2.02) | −0.04(−0.50) | 0.13(−0.34) |
  | Cl⁻ + CH3Cl(錯体基準、ΔE_rxn は分離 → 錯体) | 13.90(+0.46) | 10.45(−2.99) | −11.27(−0.67) | −11.18(−0.58) |
  | OH + CH4(分離基準) | 5.73〜5.81(−1.2) | 0.75〜0.80(−6.2) | −14.54(−1.77) | −15.02(−2.25) |

  障壁は M06-2X で平均 0.83、最大 1.79(PBE0/SVPD は 2.55、6.2)。ΔE_rxn は M06-2X で平均 1.17、最大 2.30 で、PBE0/SVPD(1.05、2.25)より良くない。HONO の trans/cis の順は M06-2X で逆になる(差 0.04、CCSD(T) は trans が 0.47 低い)。
- **SN2 の δG_eff**: M06-2X 層 14.53、CCSD(T) 層の composite(同じ run の写しに `energy_method: ccsd-t_def2-tzvpd` の thermo を足したもの、`M9/real/s20_sn2_panel_ccsdt`)14.07 で、差 +0.46(合格条件 ±1)。同じ比べ方で HCN −1.79、CH3O• −0.72、HONO +0.09。
- **S6**: 分離基準の ΔE‡ は 5.73〜5.81(合格条件 5〜7)。δG_eff は 790468f505 が 0.85 → 5.51、f9bc の split1 が 16.0 → 20.5、split2 が 40.6 → 44.8(W3_s6 も同じ幅で上がる)。
- **順位と値の変化**: run の中の順位は 36 本とも変わらない。R1〜R4 は 1 反応ずつなので順位は 1 のまま、値は R1 42.18 → 41.57(1,2-H 移動。M06-2X の ΔE‡ は PBE0 より 0.6 低く、CCSD(T) からは −1.79 と PBE0 の −1.19 より離れる)、R2 11.82 → 9.88(PBE0 のねじれの障壁の過大 +2.0 が消える)、R3 3.84 → 4.38(ΔE‡ 4.26 → 4.80。実験の反転障壁 約 5.8 にどちらも届かない)、R4 same_basin。ほかに δG_eff が大きく上がったのは H 移動と求核置換で、PBE0 の非局在化誤差(H 移動 −4.2、求核置換 −1.9 の系統誤差)の向きと合う: H + H2 6.37 → 12.71(ΔE‡ 4.56 → 10.91。参照値 9.6、HTBH38)、Cl⁻ SN2 11.08 → 14.53、I⁻ SN2 6.88 → 9.61、二重 H 交換(S10)32.29 → 36.64、マロンアルデヒドの PT 0.00 → 0.93(ΔE‡ 2.01 → 2.69)、HONO → HNO2 52.02 → 56.83(ΔE_rxn 1.94 → 8.03)。ねじれ・配座(DME +0.17、H2O2 +0.11、シュウ酸 −1.0)と (HF)₂ の交換(−0.59)は 1 kcal/mol 前後である。R2 の δG_eff(9.88)が ΔG‡(9.79)より大きいのは、cis と trans が同じ状態で、層では cis が 0.09 低くなり G_R がその最小になるからである。
- **電子状態の違う SP**: S5 の CH3···O2 錯体(PBE0 の freq は BS、⟨S²⟩ 1.71)の M06-2X の SP は ⟨S²⟩ 0.759 の別の解に収束し、分離した単量体より 46 kcal/mol 高い。この極小は `spin_contaminated` で状態の G から外れ、会合の反応物側は単量体なので順位には入らない(ΔE_rxn は分離基準 −37.18、PBE0 −36.50)。層の SP の電子状態を freq と照合するゲートはない(§19.6)。

### 19.6 実証していないこと

- 遷移金属・Te を含む系の M06-2X の SP(H2Te は順位の行がなく SP の対象がない)。I は SN2 で通したが CCSD(T) の参照はない。
- 層の SP が freq と違う電子状態に収束したことの検出。spin_ok は ⟨S²⟩ の汚染だけを見るので、BS の freq 点で純粋な二重項に落ちた SP を止めない(S5 では順位に影響しなかった)。
- PBE0 に鞍点がない(山が解像度未満の)反応の障壁。層では取り戻せない([design.md](design.md) §7.4)。

## 20. S1 バックエンドの正しさ(X5-1 / G4-P4 = G8-P4、G2-P6、G6-P3)

分析の §3 X5-1、§2.2 G2-P6、§3 X7 の 2 にあたる。規則は [design.md](design.md) §7.1・§7.3・§8。

### 20.1 変えたこと

- 失敗の分類: autoz の失敗は致命的な文言(`insufficient internal variables`、`geom_binvr: #indep variables incorrect`、`regeneration of autoz failed`)だけにした。Cartesian に切り替えたという無害な注記は、正常終了なら失敗でなく、SCF 未収束や maxiter と並べばそちらの種類になる。
- 継続: 駆動系ジョブの継続はすべて最後の frame から。opt の autoz も最後の frame から Cartesian で続ける(vectors だけを引き継ぎ、driver の Hessian は捨てる)。frame がまだない autoz と saddle の autoz は従来どおり同じ始点から Cartesian。
- SCF の救済(G2-P6): `_scf_rescuable` を削除した。DFT は開殻も `cgmin` の後に通常の SCF を 1 回回して ⟨S²⟩ とエネルギーを読み、2 つが qrc_drop より離れたら `rescue_solution_changed` にする。⟨S²⟩ のない開殻の出力を Evidence にしない検査(`s2_not_reported`)は残した。
- 後の wave が読む欄: `Failure.energy_hartree`(maxiter で止まった saddle の最後の歩、G1-P3)、`Evidence.gradient`(opt と saddle の最後の DFT の勾配。最終フレームになければ `gradient_not_at_final`、G5-P1)。VAL7 と r9 の NWChem の opt・saddle の Evidence 487 件で、最後の勾配のブロックはすべて最終フレームにあった(5.2e-7 Å 以内)。
- 参照層の容量(G6-P3): wsl_local の `memory_mb_per_rank` を 1,200 → 2,000 MB にした(site の規則 4 × 2,000 MB = 8 GB ≤ 0.8 × MemTotal。新しいつまみはない)。WFT の deck だけ heap 5%・stack 25%・global 70% に分け、凍結する芯を `freeze <n>` で明示する(NWChem の `freeze atomic` の芯 − def2-ECP の電子。I は 4s4p の 4。計画の例の「I は 4s4p4d の 9、I⁻ + CH3I は 19、CH3O は 1」は `freeze atomic` の規則と合わない: CH3O は VAL7 s20 の出力で凍結 2、全電子の Br も 3d を相関させる)。Kr より重い原子を含む CCSD(T) の鍵は `frozen_core` を持つ。手法パネルは `spin_ok` に落ちるエネルギーを使う値を notes に書き、最小・最大から外す。
- しきい値と予算は変えていない。

### 20.2 検証の道具(repo の外、`/home/user/hfauto_r9/tools/`)

- `harness_autoz.py WAVE`: VAL7 の s10 の QRC の側の opt a0a8c704(NH3·HF、鞍点の freq の初期 Hessian)を、元の始点から新しい run で作業ツリーの NWChem エンジンに流す。合格: ジョブ鍵が元と同じ、attempt_00 が frame を残して autoz で止まる、attempt_01 が Cartesian で attempt_00 の最後の frame(1e-4 Å 以内)から始まり 3 歩以内で収束する、元とのエネルギー差 < 1e-6 Eh、元の run の DFT 極小に対する `Registry.find` が両方に同じ basin を返す。
- `harness_rescue.py WAVE`: r6 の h3p3c の bead(P3c の H3 二重項 ZTS の bead 4。失敗した string の試行が残した vectors から始める)を energy ジョブとして流す。合格: attempt_00 が scf_not_converged、attempt_01 が cgmin と通常の SCF の 2 task、結果が ⟨S²⟩ 付きの Evidence で `spin_ok` を通る。同じ bead を atomic guess から解くと ⟨S²⟩ 0.9868 の解に収束する(`/home/user/hfauto_r6/SA/sac_probe/scf/h3/it100`)ので、救済した解が `spin_ok` に落ちることもありうる。そのときは結果として記録し、spin_tol は動かさない。
- どちらも `R9/.runcase.lock` を取り、`R9/WAVE/harness/<name>/summary.json` に書く。

### 20.3 変更前のコード(HEAD 92491f3 の写し)での harness(`/home/user/hfauto_r9/S1/harness/*_head`)

- autoz: attempt_00 は 21 frame の後に autoz で止まった(51 s)。attempt_01 は始点から Cartesian でやり直し(attempt_00 の最後の frame から 0.918 Å)、19 歩(42 s)。エネルギーは元と 2e-9 Eh、basin も同じで、VAL7 の記録を再現した。不合格は continue と steps だけである。
- rescue: attempt_00 は SCF 未収束なのに、ほぼ直線の H3 で出る無害な AUTOZ の注記のために `input_invalid:autoz` に分類され、同じ始点と vectors から Cartesian で再投入された。attempt_01 も SCF 未収束で、開殻なので救済されず `scf_not_converged` で終わった。救済の前に分類の誤りがあった(20.1 の 1 つ目で直す)。

### 20.4 結果(統合、`/home/user/hfauto_r9/S1/`、`/home/user/hfauto_r9/S1a/`)

- 再生(`S1/replay`、CUR の 36 run、strict): 全 stage で misses 0、新しいジョブ 0。記録の差は report の method_panel.csv の sha256 だけ(36 run、パネルの 3 run は 2 つの report)。列 notes が増えたためで、値は s5_ch3_o2 の 1 行を除いて同じ。s5 の b3febc8052 の PBE0 の ΔE_rxn は BS 錯体の freq のエネルギー(spin_ok 不合格)を使うので notes に `spin_contaminated:dE_rxn` と書かれ、幅から外れた(残る幅は M06-2X の −83.15 だけ。その SP が別の電子状態に落ちた件は §19.6 のとおり)。ranking.csv と結論の記録は同じ。
- harness autoz(`S1/harness/autoz`): 6 項目すべて合格。attempt_00 は 21 frame の後に autoz(54 s)、attempt_01 は frame 21 から Cartesian で 1 歩・3.3 s(HEAD は始点から 19 歩・42 s)。元とのエネルギー差 1.4e-8 Eh、basin は両方 basin_nh3_hf1_c00_ff60c312。
- harness rescue(`S1/harness/rescue`): scf・rescue・evidence は合格、`spin_ok` は不合格。attempt_00 は SCF 未収束(今度は autoz に誤分類されない)、attempt_01 は cgmin と通常の SCF の 2 task で ⟨S²⟩ 0.9868 の Evidence。atomic guess の解と 1.6e-10 Eh で同じ解で、この bead の UKS の基底解そのものが汚染されている。spin_tol は動かしていない(S1a のプローブ `S1a/h3p3c` も ⟨S²⟩ 0.9861、cgmin → 通常の SCF の差 3.7e-8 Eh)。
- DFT の回帰(`S1/real/hcn_fresh`、memory 2,000 MB の新しい run): 極小 2 つは BASE と 1e-8 Eh 以内、TS は CUR と 3e-10 Eh(BASE とは M4 の鞍点の Hessian の変更前で 5.5e-8 Eh)、共通の鍵のジョブは 1.3e-10 Eh 以内、δG_eff 41.5739 は CUR と 1e-7 kcal/mol。
- CCSD(T)(`S1/real/s14_ccsdt`、CUR の s14 の写しに CCSD(T) のパネル): 極小と TS が 1,396 s と 1,405 s(4 ranks、heap/stack/global 100/500/1,400 MB、`freeze 5`)で完走した(1,200 MB の `memory total` では (T) が配列を確保できなかった)。ΔE‡ は CCSD(T)/def2-TZVPD 3.47、PBE0/SVPD 2.01、PBE0/TZVPD 1.90、M06-2X/TZVPD 2.69、ωB97X-D3/TZVPD 3.24 kcal/mol。文献の CCSD(T)/CBS 約 4 kcal/mol より 0.5 低い(基底の不完全さの範囲)。
- 凍結する芯(S1a の実計算): HI は `freeze 4` で NWChem の number of core 4、E −297.735051343150(`freeze atomic` は −297.821560434727、G27)。I⁻···CH3I は `freeze 9`、number of core 9、E −634.213649871623、1,024 s(G35)。
- 実証していないこと: 開殻の救済が `spin_ok` を通る例(h3p3c は汚染された解しかない)。ECP 原子を含む CCSD(T) の参照値との比較。`Failure.energy_hartree` と `Evidence.gradient` を使う側(S3 の G1-P3・G5-P1)。それ以前に保存された JobStore の記録ではどちらも None である。

## 21. S2 決定表と仮説の単位の整理(X5-2〜4、X2-3/4、X6-5)

分析の §3 X5 の 2〜4(G1-P3、G3-P2、G3-P3 + G7-P1、G3-P6、G7-P5)、X2 の 3・4(G8-P7、G8-P3)、X6 の 5(G3-P4)にあたる。規則は [design.md](design.md) §6.2・§7.1・§8。

### 21.1 変えたこと

- 行 1 は「端点に reaction-paths の DFT 極小がない」の 1 節にし、理由を `endpoint_without_dft_minimum`(旧 `endpoints_not_on_one_pes`)にした。level_key の節は削除した。両端が 1 つの PES にあることは、設定の読み込みで `PipelineConfig._check_stages` が保証する(minima(dft) と reaction-paths の method が 1 つ)。
- TS の検証と QRC を 1 つの action(VALIDATE_AND_CONNECT: 別ジョブの freq → `is_first_order_saddle` → χ → QRC。両側が同じ basin なら 2 倍の振幅でもう 1 回)にした。旧行 9・10(`ts_validated` → CONNECT、`saddle_converged` → VALIDATE_TS)は行 9 の 1 行になり、`connection_retry`、`connection_attempts`、`ts_check` を削除した。行 8 の判断(same_basin、same_state、それ以外)は表に残る。決定表は 14 行になった。
- FIND_PATH を 1 行(行 13、理由 `dft_path`。旧 `no_dft_path` と `next_chunk`)にし、`STRING_CHUNKS`(3)を削除した。string は、saddle の試行が残り、走らせた string の数が saddle の試行数以下で、まだ string がないか最新のプロファイルが unavailable でない string のときに走る。種を出さなかった string の後は、すぐに打ち切り(`attempts_exhausted`)に進む。以前は種がなくても 3 チャンクまで続けた。
- 続きの段数(G3-P6): `Seed.depth`。maxiter からの再開と高次の鞍点の押し出しは、どちらも段数 + 1 の続きで、2 段(`actions.MAX_DEPTH`)までにした。再開は種の `saddle_restart` ではなく、`refine_saddle` の中で同じ試行として 1 回だけ行う(締め切りの前だけ)。試行の数え方は前と同じで、種 1 つが 1 試行(押した種も数える)、再開は数えない。`saddle_restart` の特別扱い(数えない、上限でも通す、再開の再開をしない)と、行 12 のその例外は削除した。分析は「段数 0 の種だけを数える」としていたが、押した種は試行として数えたままにした(予算の意味を変えない)。
- 再開の上限(G1-P3): maxiter で止まった saddle の最後の歩のエネルギー(`Failure.energy_hartree`、S1)が、最新の DFT プロファイルの最大 + resolution_kcal を超えるなら再開しない(注記 `saddle:above_path_bound:<screen|string|scan>`)。エネルギーのない Failure(S1 より前に保存されたもの)とプロファイルのない case は、今までどおり再開する。
- 仮説の単位(G8-P7): 仮説は case キー(`hypotheses.pair_key`: 両端の状態が違えば 2 つの状態の組、同じなら 2 つの極小の組)で重複を除く(宣言反応はそれぞれ残す)。どの仮説も、そのキーの全候補の低レベル TS を優先順に `low_level_ts`(タプル)に持つ(持っている TS と構造で 1 つの basin に入るものは除く。以前は 1 つだけ)。SCREEN の近道はそれを順に 3 点で分け、single の種を順に積む(記録する判定は最初の single)。分割の子は、同じキーの case が駆動済みか待ち行列にあれば、その駆動を待ち、結論が出ていれば駆動せずにその結果(outcome と claim)を写す(理由 `same_as:<reaction_id>`、ジョブも log もない)。その case が unresolved なら、子は自分の予算で駆動する(統合で追加。§21.6)。
- 分割の深さ(G8-P3): 深さが `max_split_depth` に達した case の子は、駆動せずに unresolved_within_budget(理由 `split_depth`、ジョブも log もない)として記録する(`classification.undriven`)。以前は記録せずに捨てていた。
- 例外の閉じ込め(G3-P4): `drive_case` が decide → action の繰り返しを囲む。例外の case は unresolved(`error:<型>`)で、それまでに登録した極小・species・計算も artifact に出す。`reaction_paths._drive` を削除した。
- しきい値と予算(resolution_kcal、max_saddle_attempts、max_split_depth、walltime_h)は変えていない。

### 21.2 期待される差(QM なし、CUR の記録と case log から)

- 行 1、FIND_PATH の 1 行、VALIDATE_AND_CONNECT、続きの段数: CUR の記録は変わらない見込み。string を走らせた case は hcn_noscreen、hf_dimer_swap_noscreen、s5_ch3_o2 の b3e2f366b6_split2 の 3 つで、どれも種を出し、2 本目の条件も前と同じに満たす。QRC の 2 段目の振幅(`connection_retry`)に進んだ case はない。押した種(W3_s6・W5_s6・s6_oh_ch4 の 790468f505)と再開(W5_s6 の f9bc3a1cd8_split1)は、前と同じ試行の数え方で同じ経過になる。case log の行動と理由の列は変わる(検証と接続が 1 行、再開は `seed:saddle_restart` の行がなくなり注記だけ)が、再生は理由の文字列を比べない。
- 再開の上限: CUR で再開したのは W5_s6 の f9bc3a1cd8_split1 だけ(screen_hei → maxiter → 再開 → elementary)。CUR の Failure にはエネルギーがないので、再生では今までどおり再開する。
- 仮説の単位: W3_s6 の 460166732f と f9bc3a1cd8 は同じ状態の組(CH3 + H2O ↔ CH4O + H)で、1 つの仮説になる。種は d7d2ff63 の xTB TS(C 上の置換、−1347.6i に着いたもの)と b9c080cc の xTB TS(O 上の置換、−1417.9i)の 2 つ。formaldehyde の 2 つの宣言反応は、宣言なのでそれぞれ残る。
- 子の再利用: CH3 + H2O ↔ CH4 + HO の子(s6_oh_ch4 と W5_s6 の f9bc3a1cd8_split1)は 790468f505 と同じ状態の組なので、駆動せずにその結果を写す(分析 X2 の検証 (c))。W3_s6 の同じ組の子(460166732f_split1)は、790468f505 が attempts_exhausted なので写さず、前と同じに駆動する(§21.6 の統合の修正)。
- 分割の深さ: CUR で `max_split_depth` 2 に達した子はない。

### 21.3 合格条件(統合)

1. 再生(`tools/replay.py S2 --src CUR`、全 run): 次の差の run だけを `--allow-new` で報告し、ほかは strict。(i) 種を出さなかった string の後の打ち切り(CUR では 0 件の見込み)、(ii) 状態の組でまとめた仮説とその追加の種(W3_s6)、(iii) 子の結果の再利用(W3_s6、W5_s6、s6_oh_ch4)。
2. W3_s6: CH3 + H2O / CH3OH + H の組で、O 上の置換(b9c080cc)と C 上の置換(d7d2ff63)の両方の種が試されること。
3. s6_oh_ch4 と W5_s6: 子 split1 が 790468f505 の結果を再利用すること。
4. 分割の深さ: `/home/user/hfauto_r9/inputs/oxalic_split0.yaml`(known_endpoints の reaction-paths に `policy: {max_split_depth: 0}`。検証専用で、製品の既定は 2 のまま)で CUR/oxalic_two_step の写しを paths から再ステージする。親は multi_step のまま、子 2 つが `split_depth` で記録され、新しいジョブがないこと(分析 X2 の検証 (g))。
5. 再開の上限
   - QM なし(`tools/bound_check.py`): VAL7/s5 の eb9bb06f(最後の歩 −189.73505)と b0660ccc(−189.73438)は上限を超えて再開しない。W5_s6 の f9bc3a1cd8_split1 の再開のフレームは上限の下にあり、再開が残る。
   - 実計算: CUR/W5_s6 の写しを paths から `--retry-failed geometry_maxiter` で再ステージし(maxiter の saddle を取り直してエネルギーを記録する)、f9bc3a1cd8_split1 が再開を経て収束すること(予算 1.5 h)。すべての仮説を流すと split1 は same_as になって駆動されないので、検証専用の入力 `/home/user/hfauto_r9/inputs/w5_split1_restart.yaml`(paths の `reaction_ids` を f9bc3a1cd8 だけにする)を使う。分析の例の W3 の split2_split2 は、M1 で記録が消えたので CUR/W3_s6 にない。
6. 行 1(各 15 分)
   - `/home/user/hfauto_r8_probe/G3-1_row1`(平面 NH3 を宣言の端点にし、mode_follow 0): blocked_upstream(`endpoint_without_dft_minimum`)。
   - 既定の mode_follow の `/home/user/hfauto_r9/inputs/nh3_planar_row1_default.yaml`(known_endpoints): 結果をそのまま記録する。コードでは、平面の端点は mode-follow の後、入力に近い側の NH3 の basin に加わる(`endpoint_was_saddle`)ので、DFT 極小を持ち、行 1 は発火しない見込みである(分析 G3 の根本原因 1 は「既定でも行 1 に落ちる」としていた)。
7. 単体テスト: 決定表の試験に加え、振る舞いの試験(虚モードのない鞍点は探索の行に進む、W3 の押し出し → maxiter → 再開の連鎖が段数の範囲で再現され交互の無限連鎖がない、上限、種のない string の後の打ち切り)、状態の組の重複除去、複数の低レベル TS、split_depth の記録、子の再利用、例外の閉じ込め(`_register` の後の例外で minimum の artifact が manifest に残る)。

### 21.4 見かけの合格の見分け方

- 決定表の試験を書き換えただけで、振る舞いの試験がない。
- 予算(max_saddle_attempts、max_split_depth、walltime_h)を変えて通す。oxalic の max_split_depth 0 は検証用の入力だけである。
- manifest から黙って消えた子(S2 より前の捨て方)を split_depth と数える。
- 行 1 を、宣言のすり替えや mode_follow・ゲートの変更で発火させる。
- 上限の確認で、SCREEN と string のどちらのプロファイルの最大を使ったかを取り違える(注記にプロファイルの種類が残る)。

### 21.5 結果(QM なし)

再開の上限(`/home/user/hfauto_r9/tools/bound_check.py`、出力は同じ場所の `bound_check.out`): 記録された case の maxiter で止まった探索ごとに、最後の attempt の最後の歩のエネルギー E_last と、その探索の前の最新の DFT プロファイルの最大を比べた(プロファイルは case のジョブから組み立て直す)。

| run | case | 探索(種) | E_last − 最大(kcal/mol) | プロファイル | 判定 |
|---|---|---|---|---|---|
| VAL7/s5_ch3_o2 | b3febc8052 | eb9bb06f(screen_hei) | +27.21 | screen | `above_path_bound:screen`(再開しない) |
| VAL7/s5_ch3_o2 | b3febc8052 | b0660ccc(path_hei) | +31.99 | string0 | `above_path_bound:string`(再開しない) |
| W3/s6_oh_ch4 | f9bc3a1cd8_split2_split2 | e5ef2bd9(higher_order_retry) | −21.29 | screen | 再開する(W3 の経路が残る) |
| CUR/W5_s6 | f9bc3a1cd8_split1 | 204845c7(screen_hei) | −68.51 | screen(NEB 失敗で IDPP) | 再開する |

上限(+1.0 kcal/mol)の近くの値はない。VAL7/s5 の 2 本は HEI の種から 27〜32 kcal/mol 登ったフレームで、上限で切られる。再開して収束した 2 本は上限より十分下にある。

### 21.6 結果(統合)

再生(`tools/replay.py S2 --src CUR --allow-new W3_s6,W5_s6,s6_oh_ch4`。CUR は S1 の再生を指す。36 本、4.2 分): 33 本は strict に合格(misses 0、記録の差 0)。hits 1,006、misses 4。差は 3 本の S6 だけで、すべて許した差である。

- `low_level_ts` は Geometry | null からリストになった。replay.py は null を []、G を [G] と読み替えて比べる。
- (ii) 状態の組でまとめた仮説: W3_s6 では f9bc3a1cd8 が 460166732f にまとまり(記録・その QRC の species と極小・thermo が消える)、460166732f は d7d2ff63(C 上の置換)と b9c080cc(O 上の置換)の 2 つの低レベル TS を持つ。近道は両方を xTB freq と DFT SP の 3 点で分け、どちらも single で種になった(W3_s6 では f9bc3a1cd8 の前のジョブが再利用され、新しいジョブは 0)。最初の種(C 上の置換、−1347.5i)が QRC で CH4 + HO を経る 2 段(multi_step)を示して case を閉じたので、O 上の置換の種は saddle まで試されていない。前の f9bc3a1cd8 の O 上の置換の素反応(−1418.1i、順位付き)は W3_s6 の記録から消えた。W5_s6 と s6_oh_ch4 の f9bc3a1cd8 も 2 つの TS を持ち、2 つ目の xTB freq と DFT SP が新しいジョブ(各 2 件、約 15 core-s)。結論は前と同じ multi_step。
- (iii) 子の再利用: W5_s6 と s6_oh_ch4 の f9bc3a1cd8_split1 は `same_as:rxn_discovery_790468f505` になり、駆動しない(saddle は 790468f505 の −486.0i / −486.4i。自前の −482.8i / −476.7i の代わり。dG_eff の差は +0.05 / −0.04 kcal/mol)。
- (i) 種を出さなかった string の後の打ち切り: 0 件。
- W3_s6 の 790468f505(screen_hei → 押し出し → 高次 → attempts_exhausted)は strict に同じ。

統合で直したこと: 最初の再生で、W3_s6 の 460166732f_split1 が 790468f505 の unresolved(attempts_exhausted)を写し、自前の素反応(−487.0i、CH4 + HO ↔ CH3 + H2O の引き抜き)を失った。その結果、M06-2X の SP が選ばれなくなった p02 の thermo も、460166732f と 790468f505 の thermo も unavailable になった。結論のなさは共有する結果ではないので、子は同じキーの case の結論だけを写し、その case が unresolved なら自分の予算で駆動するようにした(`reaction_paths._Book`。待ち行列にある case を待った子は、待ち行列が空になってから親の締め切りで駆動する)。予算の意味は G8-P7 の前と同じである。修正後の再生で W3_s6 の 460166732f_split1 は前と同じ記録(新しいジョブ 0)。

χ(`tools/chi.py`、S2 の再生の S6 3 本): 状態の組の種から新しい鞍点は生まれていない(新しい saddle ジョブ 0)。接続して順位の付く TS 9 個の χ は 0.680〜0.798 で、0.3 未満はない。ただし、まとめた仮説の端の添字(C 上の置換)に対して、O 上の置換の TS(−1418.1i)の χ は 0.276 で、REACTION_MODE_MIN 0.3 を下回る(自分の端に対しては 0.80)。今の順序では試されないが、試されれば添字だけで not_reaction_mode と棄却される。残る課題(下)。

audit_ends(S2 の再生): 19 件を監査し、不一致 0、別経路の TS 3(上の 05c19 の O 上の置換)、緩和で状態が変わった極小を端に使う仮説 0。M1 の選択肢 (a) は要らない。

実計算(`/home/user/hfauto_r9/S2/real/`、`s2_real.sh`):

| 確認 | 入力 | 結果 |
|---|---|---|
| 分割の深さ(X2 (g)) | oxalic_split0(CUR/oxalic_two_step の写し、paths から) | ctc_to_ttt は multi_step のまま、子 2 つは unresolved_within_budget(`split_depth`)、log なし。paths は hits 17、misses 0 |
| 行 1(mode_follow 0) | `/home/user/hfauto_r8_probe/G3-1_row1` | up_to_planar は blocked_upstream(`endpoint_without_dft_minimum`)。平面 NH3 の端点に DFT 極小がない |
| 行 1(既定の mode_follow) | nh3_planar_row1_default | 平面 NH3 は mode-follow で NH3 の basin に加わり(メンバー planar_nh3、_mf1、_mf2)、行 1 は発火せず same_basin(行 2)。分析 G3 の根本原因 1 の「既定でも行 1 に落ちる」は成り立たなかった |
| 再開の上限(実計算) | w5_split1_restart(CUR/W5_s6 の写し、paths から `--retry-failed geometry_maxiter`) | 取り直した saddle(同じジョブ鍵 204845c7)は今回は maxiter で止まらず収束した(NWChem の並列の非決定性。前回は 50 歩で maxiter)。再開も上限の判定も起きず、split1 は elementary(−484.8i)。実計算での再開の上限は実証できていない(10 分) |

残る課題:

- O 上の置換と C 上の置換の両方の種を saddle まで試すこと(分析の検証)は満たしていない。種は順に試し、最初に結論が出た種で case が閉じる(G8-P7 の実装)。どの種を先にするか(低レベルの障壁 `dE_act_kcal` の順など)と、2 段の結論より直接の素反応を優先するかは未決。
- まとめた仮説の χ は、仮説の端の添字に対して測る。別の添字の経路の TS は、正しい TS でも χ が下がる(上の 0.276)。端を鞍点の添字に合わせ直して測る(`identity.basin_coords`)案があるが、縮退転位では両端が同じ添字に重なるので、そのまま入れられない。REACTION_MODE_MIN は変えない。
- 再開の上限は、S1 より後に maxiter で止まった saddle にしか効かない(それより前の Failure はエネルギーを持たない)。

## 22. S3 極小の停留性、初期 Hessian の 1 行の規則、結合変化の帯(X3-3、X1、X3-1)

分析の §2.5 G5-P1、§2.3 G3-P5、§2.4 G4-P2、§2.8 G8-P5、§3 X1・X3 にあたる。規則は [design.md](design.md) §5、§6.1、§7.1(explore、minima(dft)、「結合と状態」)。

### 22.1 変えたこと

- 極小の停留性(G5-P1): `vibrations.stationarity_gap` は二次モデルの降下量 ΔE_N = Σ g_i²/(2\|λ_i\|)(freq と同じ質量加重の内部モード、\|λ\| に床なし)を返す。`MinimumDriver` は saddle 級でない点で、opt の最終の DFT 勾配(r9-S1 の `Evidence.gradient`)と freq の Hessian から ΔE_N を求め、5e-5 Eh(`identity.BASIN_DE_HARTREE`)を超えれば振動数によらず soft にする(注記 `soft_imaginary_mode`)。基準は 2 つの極小を 1 つの basin とみなすエネルギーの幅で、新しい数ではない。勾配のない Evidence は判定しない。
- soft の押し出し: 停留していない点は Newton の歩 −H₊⁺g(`vibrations.newton_step`。§5 の正定値のモデルなので、虚モードがなくても下る)、停留した soft の点は従来どおり soft の虚モードの和の方向に、1 回だけ押して緩和する(最大の原子変位 0.4 Å)。
- `soft_minimum` の status を削除した(G3-P5)。登録(`Registry.add`、minima の `register`、connection の `_register`)は以前から minimum と同じ扱いだった。結果は常に `minimum` で、残った soft は history の `soft:persisted` と注記 `soft_imaginary_mode` に残る。
- 初期 Hessian の 1 行の規則(X1、G4-P2): NWChem の optimize は、ν < −saddle_cm1 のモードが 2 本以上の Hessian だけをそのまま書き、ほかはすべて正定値のモデル(ジョブ鍵に `hessian_model: positive`)にする。「ちょうど 1 本のときだけモデル」という engine の推測を置き換えた。minima の `init_hessian` から relaxation の seed の除外を外した。
- 結合変化の帯(G8-P5): `topology.bond_changes` が `RESOLVED_A`(0.1 Å)の帯を持ち、`resolved_bond_changes` を削除した。結合変化の定義は 1 つになり、仮説と分割の子のねじれ、会合の形、ρ と χ、R6 が同じ定義を使う。`bonds` と状態ラベルは変えていない。
- 2 本目の虚モードは停留点でだけ数える(持ち越した M4 の問い、22.3): higher_order の saddle が勾配を持ち ΔE_N > 5e-5 Eh なら、符号付きの Newton 歩(二次モデルの停留点へ、最大 0.4 Å。`minimum.newton_push` を極小と共有)の点で freq を 1 本とり、そこで −saddle_cm1 より下が 1 本以下なら、2 本目のモードに沿う押し出しの代わりにその点から `higher_order_retry`(段数 + 1、注記 `higher_order:not_stationary`)。残れば従来どおり押す。2 回目の higher_order_retry は足していない。統合で、S3a の版(Newton 点の freq なしで種にする)を直した: 本物の高次の鞍点(W3_s6 の 975217d3、−1627.8i / −356.9i、ΔE_N 5.7e-5)では二次モデルの停留点がその鞍点自身なので、押さずに戻ると続きを 1 段失う。
- しきい値(noise_cm1、saddle_cm1、BASIN_DE_HARTREE、RESOLVED_A)と予算は変えていない。tight の収束は足していない(分析 G5 で退けた案)。一次の鞍点の受理には停留性をかけていない(分析 X3-3 の「TS への拡張は今回は入れない」)。

### 22.2 期待される差(QM なし)

- 停留性: CUR の opt の Evidence は勾配を持たない(r9-S1 より前の JobStore)ので、再生では判定せず、記録は変わらない。勾配のない記録を停留とはみなさない。
- 初期 Hessian の規則: 以前そのまま書いていた Hessian(低レベルの錯体と R6 の seed の xTB Hessian、停留した soft の点の freq)を渡す opt は、鍵に `hessian_model: positive` が入って新しいジョブになる。最終構造の指紋が変われば、その下流(freq、経路の SP、saddle、QRC、sp)も新しくなる。CUR では W3_s6(4 本。1 本は DFT の soft の点の freq から押した opt)、W5_s6 と s6_oh_ch4(各 3 本)、s1_sn2_cl、s20_sn2_panel、s2_sn2_i、s10_amine_pilot2、s19_tma_hf2(各 1 本)の opt が該当する。一次の鞍点からの側(QRC、mode-follow)の鍵は変わらない。
- R6 の seed: S5 と S6 の seed の opt は、xTB Hessian の正定値のモデルから始まる新しいジョブになる(S19 の seed は M6 で失われた状態でなくなった)。
- 結合変化の帯: r − r_thr が ±0.1 Å の帯の中でしきい値を横切るだけの対は、結合変化でなくなる。ねじれの判定、未宣言の仮説の採否、ρ の項、χ の結合に効く。

### 22.3 QM なしの確認

道具は `/home/user/hfauto_r9/tools/stationarity.py`(出力 `/home/user/hfauto_r9/S3a/stationarity.txt`)。`stationarity_gap` の式で、記録の Evidence に勾配がなければ stdout から読み直し、最終構造にある勾配だけを使う(なければ unknown で、停留とはみなさない)。unknown は 0 点だった。

- 極小: VAL7 の 64 点は 5e-5 Eh を超えるものがなく、最大は S10 の 8e8cc220 の 1.2e-5。W3 の 29 点で超えるのは 1ee97b40(CH3OH···H の肩、ν1 −6.5i)の 2.66e-4 だけで、そのうち 2.5e-4 が −6.5i のモードの項である。次は 67b7f32a(CH3OH···H の領域、ν1 0.7 cm⁻¹)の 2.0e-5。CUR の 36 本(`--minima` に全 run を渡す)でも 75 点のうち同じ 1 点だけである。分析 G5 のプローブの値と一致する。
- 鞍点(S3 では判定しない): CUR の 82 点のうち 21 点が 5e-5 を超える。S6 の引き抜きの領域が 13 点(run ごとに数える)、W3_s6 の高次の鞍点 975217d3(−1627.8i / −356.9i)が 5.7e-5 で、残る 7 点(回転子の鞍点 −58.6〜−73.2i の 5 点と acac の −987i の TS の 2 点)は 5.5e-4〜1.6e-2 Eh である。この 7 点の項のほとんどは 40 cm⁻¹ 未満のモード(6 点は 8〜16 cm⁻¹ のメチルや OH の回転)にある。大振幅の回転には二次モデルが成り立たないので、この値をそのまま停留点までのエネルギーとは読めない(推定)。置換の TS(−1347.5i、−1418.1i)は 5.2e-7 Eh 以下で停留している。
- S6 の引き抜きの TS の領域(M4 から持ち越した問い): CUR/W3_s6 の鞍点は次のとおり。

  | 鞍点 | ν1 / ν2(cm⁻¹) | ΔE_N(Eh) | 322〜331 cm⁻¹ の 2 本の項 | Cartesian の gmax(Eh/bohr) |
  |---|---|---|---|---|
  | 57fc1cfb(W3 の素反応) | −481.1 / +33.4 | 1.3e-4 | 5.8e-5、4.7e-5 | 1.1e-3 |
  | 9f533dcb(460166732f_split1、M4) | −487.0 / +29.7 | 1.0e-4 | 5.5e-5、2.7e-5 | 1.2e-3 |
  | ae9d9169(同、M3) | −485.8 / +25.1 | 1.5e-4 | 1.0e-4、2.3e-5 | 1.4e-3 |
  | f66b9f98(790468f505、screen_hei) | −489.6 / −77.1 | 1.5e-4 | 7.2e-5、5.5e-5 | 1.2e-3 |
  | 57c69460(790468f505、retry) | −499.0 / −54.8 | 1.3e-4 | 8.5e-5、2.3e-5 | 1.1e-3 |

  - 項の大部分は 322〜331 cm⁻¹ の 2 本にある(ν2 の向きの項だけでは説明できないが、ν2 はこの変角と結合しているので、Newton の歩で ν2 も変わる)。移る H と CH3 の H が動き、C–H···O の角(176〜178°)が 1 Å あたり 30〜44° 変わるモード、つまりほぼ直線の C–H···O の変角である。ν2(OH のねじれ。OH の H の振幅が 76%)そのものの項は 3e-6〜7e-6 Eh と小さい。
  - 原因(saddle の出力の autoz の座標表から Wilson の B 行列を組んで確かめた): この 5 本の autoz の座標は、ほぼ直線の C–H···O をまたぐねじれ(X–C–H–O と C–H–O–H)を含み、その B の行の大きさは通常のねじれの 9〜16 倍である(直線からのずれの sin に反比例する)。そのため Cartesian の勾配(gmax 1.1e-3 Eh/bohr)が内部座標の勾配では 1e-4〜2e-4 に縮み、driver の既定の gmax 4.5e-4 で収束と判定された(57fc1cfb の最後の歩の表示は 9e-5)。
  - M4 の実計算 6 本(`M4/real/s6_*`)の引き抜きの領域の鞍点も、一意の 16 本のうち 14 本が ΔE_N 6.9e-5〜2.2e-4 Eh(Cartesian の gmax 7.5e-4〜1.4e-3)で、ν1 −476〜−522i、ν2 −77〜+41 cm⁻¹ と散らばる。残る 2 本は s6_fresh の r2 と r4 の 460166732f_split1 の TS(7ffefe2d、1f9571ef)で、ΔE_N 1.6e-6 と 2.2e-6 Eh で停留している。C–H–O は 172.5° で、autoz の座標に直線をまたぐねじれはない。この 2 本は ν1 −484.6i と −491.7i、ν2 +53.7 と +49.0 cm⁻¹ の一次の鞍点で、E は −116.031418 と −116.031417 Eh である。同じ run の 790468f505 の TS(非停留、ν2 −15.9i と −34.6i)より 1.3e-4 と 1.7e-4 Eh 低く、その TS の ΔE_N(2.0e-4、1.5e-4)と同じ大きさである。
  - Newton の歩のプローブ(S3a、`/home/user/hfauto_r9/S3a/ts_newton_probe.py`、出力は `S3a/ts_newton/`。runcase のロックの下で NWChem の gradient と freq を 1 回ずつ): 790468f505 の 2 本から、Hessian の符号を保った Newton の歩(二次モデルの停留点へ、最大 0.07 Å)を 1 回とると、f66b9f98 は ν1 −480.2i / ν2 +187.0、57c69460 は −482.5i / +117.9 の一次の鞍点になり、ΔE_N は 2.9e-5 と 1.8e-5 Eh(基準の内側)、E は 1.3e-4 と 1.2e-4 Eh 下がった。
  - したがって、M4 が記録した「平らな回転子の領域」(0.04 kcal/mol の中で ν2 が −77i〜+24)は、停留していない点の Hessian で測った値である。W3_s6 の 790468f505 の 2 本目の虚振動(−54.8i)も非停留の点での判定で、停留点では一次の鞍点になる。
  - S3 のコードは、この 2 本のような非停留の高次の鞍点を Newton の歩の点から続ける(22.1)。CUR の saddle は勾配を持たない(r9-S1 より前)ので、再生ではこの分岐は動かず、新しい saddle にだけ効く。
- 結合変化の帯(S3b、`/home/user/hfauto_r9/S3/S3b/band_audit.py` と `split_labels.py`、出力は同じ場所の `*_cur.txt`): CUR の 36 本で、仮説の集合、記録された case、χ の行は帯を入れても変わらない(差 0)。帯の中でラベルだけが割れる組は 11 組で、以前は結合変化 1 本、今は変化なし(ラベルの違いは別の状態としてだけ残る)。10 組は S19 の N1–H13(c01・c02 の seed と neutral・shared_proton の構造、r − r_thr は +0.006〜+0.043 と −0.006〜−0.149 Å。分析の「amine·HF の割れた組」)、1 組は S5 の NT2 の生成物(O4–O5 が −0.079 Å)とその DFT 極小(CH2O + HO)である。帯に最も近い本物の結合変化は acac の PT の O–H(+0.175 Å)、S10 の N–H(+0.215)、マロンアルデヒドの O–H(+0.228)で、帯の外にある。帯を 0.3 Å にすると仮説 3、case 4、χ の行 8 が変わる(感度の確認だけで、定数は 0.1 のまま)。
- ニトロメタンのプローブ: 22.6。

### 22.4 合格条件(統合)

1. QM なし(`/home/user/hfauto_r9/tools/stationarity.py`): W3 と VAL7 の全 DFT 極小で、ΔE_N > 5e-5 Eh は W3 の 1ee97b40(約 2.7e-4)だけ。VAL7 は 64 点で最大 1.2e-5(分析 G5 のプローブの値)。勾配を読めない点(unknown)を停留とみなさない。
2. 再生(`--src CUR`): 停留性は記録を変えない。初期 Hessian の鍵と帯で変わる run は `--allow-new` で報告し、新しいジョブと差分を 22.2 の理由で説明する。ほかは strict。
3. S6 の dft stage × 3(CUR の s6 の写しを `FRESH_ENGINES=nwchem`、`FROM=dft`、`--to dft` で、rank 4 で 2 本と rank 2 で 1 本): すべての DFT 極小で ΔE_N ≤ 5e-5 Eh。basin の数は報告するが、基準にしない。
4. ニトロメタン(22.6): status は minimum で、`soft:resolved` か注記 `soft_imaginary_mode`。`soft_minimum` は出ない。
5. G4-P2(CUR の s5 と s6 の写しを `FROM=dft`): R6 の seed の終点(極小か鞍点か、どの basin か)と下流のジョブ数が CUR と同じ。崩れるかどうかが変われば不合格。
6. 帯: 再生の差(ねじれ、仮説の採否、分割の子、ρ)をすべて挙げて説明する。S19 の amine·HF の割れた組(分析では未検証)の扱いも報告する。

### 22.5 見かけの合格の見分け方

- run 間で basin の数が一致したことを合格にする。判定するのは各点の停留性である。
- tight の収束を足す。noise_cm1、saddle_cm1、BASIN_DE_HARTREE を動かして soft の分岐に入れる、または避ける。系ごとのしきい値や、2 つ目の帯の定数を置く。
- 勾配のない Evidence を停留とみなす。
- ニトロメタンをずれ形の入力にして、soft の分岐を通さずに合格にする。

### 22.6 検証セットに足した系: ニトロメタン

- 入力: `/home/user/hfauto_r9/inputs/nitromethane.yaml`(repo の外)。構造はプローブ `/home/user/hfauto_r8_probe/G3-4_soft` の `ch3no2.xyz`(C–H の 1 本が NO2 の面にある重なり形)、pipeline は同じプローブの `pipeline.yaml`(structures と minima(dft)、`include: all`)。
- 化学: CH3 のねじれは 6 回対称の数 cm⁻¹ の障壁の準自由回転子で、重なり形は停留した柔らかい鞍点である。これまでの検証セットには、DFT の soft の分岐を通す系がなかった(§4 の DFT の soft_minimum)。
- プローブ(S3 より前のコード): 1 本目の opt は重なり形で止まり、ν1 −35.7 cm⁻¹(soft。NWChem の表示は −36.95 で、分析 §2.3 の値)、ΔE_N 7.6e-8 Eh で停留している。虚モードの方向に 1 回押すと、ずれ形の極小(ν1 +36.4 cm⁻¹、ΔE_N 7.3e-8 Eh)に入り、2.8e-5 Eh(6.1 cm⁻¹)低い。history は `freq:soft`、`soft:resolved`。
- S3 で期待すること: 停留しているので、Newton の歩ではなく虚モードの押し出しを使う。status は minimum で、`soft:resolved`(残れば `soft:persisted` と注記 `soft_imaginary_mode`)。押し出しの opt は鍵に `hessian_model: positive` が入るので、プローブのジョブは再利用しない。

### 22.7 実証していないこと

- TS の停留性。S6 の引き抜きの TS は、受理したものも ΔE_N 6.9e-5〜2.2e-4 Eh で、極小の基準では非停留である(22.3)。原因は、ほぼ直線の角をまたぐ autoz のねじれで内部座標の勾配が縮み、driver が収束と判定することである。停留した 2 本の TS は ν2 が正で 0.08〜0.11 kcal/mol 低く、790468f505 の 2 本も Newton の歩 1 回で一次の鞍点になった(プローブ)。2 本目の虚振動は S3 から停留点でだけ数える(22.1)が、実計算の saddle で `higher_order:not_stationary` の分岐が通るかは 22.8 の再生による。一次の鞍点の受理に停留性をかけるには、40 cm⁻¹ 未満の大振幅のモードで ΔE_N が二次モデルの外に出る(acac の TS で 8e-3 Eh)ことを先に扱う必要がある。高次の判定の ΔE_N もこの影響を受けうるが、種を変えるのは Newton 点の freq で 2 本目が消えたときだけである。
- 極小でも、ほぼ 0 の振動数のモードは ΔE_N を大きくしうる(W3 の 67b7f32a で、0.7 cm⁻¹ のモードの項が 2.0e-5)。W3・VAL7・CUR で、それで基準を超えた例はない。
- Newton の歩で肩から抜ける実例。W3 と VAL7 で非停留の極小は W3 の 1ee97b40 だけで(22.3)、S3 の再生ではその opt が正定値のモデルの初期 Hessian で新しく走り、停留した極小に直接入った(22.8)。実計算 4 本(S6 の dft × 3、ニトロメタン)にも非停留の極小は出ず、極小の Newton の押し出しと TS の `higher_order:not_stationary` は単体テストでしか通っていない。
- 点群の判定が、ほぼ Cs の停留した TS で割れる(22.8、M7 の領分): 同じ TS(E の差 3e-7 Eh)を W5_s6 は C1(鏡像 2、R ln 2)、W3_s6 と s6_oh_ch4 は Cs と判定し、dG_eff が 5.53 と 5.95 kcal/mol に分かれる。1e-3 Å の雑音で Cs と C1 が入れ替わる(`tools/symcheck.py` の UNSTABLE)。

### 22.8 結果(統合)

道具の出力は `/home/user/hfauto_r9/S3/`(`stationarity*.txt`、`replay/summary.txt`、`real/logs/`)。

- 統合で直したこと: S3a の版は、非停留の高次の鞍点を Newton の点から続けるだけで、そこの ν2 を見なかった。本物の高次の鞍点(W3_s6 の 975217d3、−1627.8i / −356.9i、ΔE_N 5.7e-5)では二次モデルの停留点がその鞍点自身なので、押さずに戻り、続きを 1 段失う。仕様どおり、Newton の点で freq を 1 本とり、2 本目が消えたときだけその点を種にする(22.1)。単体テスト(消える・残る・停留の 3 通り)を足した。
- QM なし(`stationarity.txt`): 22.3 と同じ。W3 の 1ee97b40 だけが 2.66e-4、VAL7 の 64 点は最大 1.2e-5。
- ニトロメタン(実計算、4 ジョブ、4 分): 1 本目の opt は重なり形で止まり ν1 −35.7 cm⁻¹(soft、停留)。虚モードの方向に 1 回押した opt(鍵に `hessian_model: positive`)でずれ形の極小(ν1 +36.2、ΔE_N 5.3e-10 Eh)に入り、2.79e-5 Eh 低い。status は `minimum`、history は `opt`、`freq:soft`、`soft:resolved`。`soft_minimum` は出ない。合格。再生の集合(`tools/runs.json`)に足し、CUR からの strict の再生も PASS。
- 再生(`--src CUR` = S2/replay、37 本目のニトロメタンを除く 36 本): 27 本は strict で PASS(初期 Hessian に高次の鞍点の Hessian だけを渡す dme_c2v_seed、h2o_hf_inversion、水と HF の二量体、S18 も鍵は変わらない)。`--allow-new` の 9 本の差は次のとおり。帯(S3b)による差は 0 件で、予測どおり。
  - S1、S20 の SN2 パネル、S2、S10、S19: 錯体の opt の鍵が変わり(xTB の Hessian が正定値のモデルになる)、同じ極小に入る(E の差 ≤ 3.4e-8 Eh)。最終構造の指紋が変わるので、下流の経路の SP、saddle、QRC、sp が新しいジョブになる。結論は同じで、数値の差は dE_act で最大 4.4e-4 kcal/mol(S2)。
  - S5: R6 の seed の opt が xTB Hessian のモデルから始まり、CUR と同じ basin の極小に入る(E の差 2.2e-7 Eh)。dft のジョブは +1(その xTB freq)、paths 35 と sp 8 は CUR と同じ。会合の走査の点の opt は始点が変わるので新しいが、経路の記録は同じ。G は低い振動の多い錯体で 2.2e-4 Eh 動く。
  - S6(W3_s6、W5_s6、s6_oh_ch4): R6 の seed の opt は、CUR では初期 Hessian なしで DFT の鞍点に止まり、mode-follow の両側が同じ極小 mf1 に入っていた。S3 では xTB Hessian のモデルから直接その極小に入る(E の差 1.6e-8 Eh)。崩れるかどうかは変わらないが、終点の「極小か鞍点か」は鞍点 → 極小に変わる(G4-P2 の判定を字義どおりには満たさない)。dft のジョブは 18 → 15。極小の id が `..._mf1` から `spc_relax_oh_ch4_p02` に変わるので、その端をもつ仮説の id が 790468f505 から 1f266f1727 に変わる。
  - W3_s6: CUR で W3 だけにあった肩 1ee97b40(05c19 の極小、−6.5i、ΔE_N 2.7e-4)が、正定値のモデルから走り直した opt で停留した極小(ΔE_N 6.7e-10)に入り、bfac3 の seed と同じ basin になった。仮説は 460166732f から、ほかの S6 と同じ f9bc3a1cd8 に変わり、W3 だけの −481.1i の枝はなくなった。CUR で未解決(attempts_exhausted)だった 790468f505 は、1f266f1727 として素反応で結論した。
  - 3 本の S6 は同じ結論になった: 1f266f1727 は素反応、TS は −488.7〜−490.7i、ν2 +48.8〜+50.1 cm⁻¹ の一次の鞍点で、E は −116.031419 Eh(3 本の差 3e-7)、ΔE_N 6.9e-7〜9.0e-7 Eh(停留)。CUR の非停留の TS(−486.4i、ν2 −37.5i、ΔE_N 1.3e-4)より 0.09 kcal/mol 低い。f9bc3a1cd8 は multi_step、split1 と split2 は素反応のまま。どの新しい case でも higher_order の判定は起きなかった(TS の Newton の分岐は通っていない)。
  - dG_eff: 1f266f1727 は W3_s6 と s6_oh_ch4 で 5.95、W5_s6 で 5.53 kcal/mol(CUR の s6 は 5.50)。差 0.42 は RT ln 2 で、W5 の TS だけを C1(鏡像 2)と判定したため(22.7)。CUR の TS は C1 の非停留点だった。
- S6 の dft stage × 3(CUR の s6_oh_ch4 の写し、NWChem を消して `--from dft --to dft`; rank 4 × 2 と rank 2、12〜19 分): 3 本とも DFT 極小は 5 点(CH4、OH、p02、05c19 + bfac3、R6 の seed)で、ΔE_N はすべて ≤ 2.9e-7 Eh。soft の点は出なかった。basin の数は基準にしない。
- 再生の後の CUR: 37 本の DFT 極小 75 点で ΔE_N > 5e-5 Eh は 0(最大 1.2e-5)。受理した TS で超えるのは acac の 1 本だけ(−987.2i、ν2 −36.9i の大振幅のメチル、8.6e-3)。
- CUR は `/home/user/hfauto_r9/S3/replay/<name>`(37 本)を指す。

