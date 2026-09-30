# 実計算による検証の記録

round 7 の改良(W0〜W5、結果報告 [2026-09-30_round7_result.md](reviews/2026-09-30_round7_result.md))の実計算の記録である。検証セットはレビュー `docs/reviews/2026-09-27_platform_review.md` §6 と、round 7 で足した分岐用の系。round 6 までの記録(VAL を含む)は git の履歴(`git show f2309b9:docs/validation.md`)にある。

## 1. 条件

- コード: 各 wave の作業ツリー(W1 = `707e656`、W2 = `26389e3`、W3 = `2eaaf53`、W4 = `c5b515e` でコミットしたもの)。VAL7 は W5 の作業ツリー(W4 に挙動を変えない削除だけを足したもの)で、§6 の全体と round 7 の系を `/home/user/hfauto_r7/VAL7` に流し直した(`/home/user/hfauto_r7/VAL7/chain.sh`、比較は `/home/user/hfauto_r7/VAL7/val7check.py` と `/home/user/hfauto_r7/VAL7/val7walls.py`)。数値の比較先 VAL は round 6 の検証(`/home/user/hfauto_r6/VAL/runs`、コード `d5fd24c`)。
- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0、pymsym 0.3.5)。
- 手法: 停留点と振動は PBE0-D3BJ/def2-SVPD(I は def2-ECP)、低レベルは GFN2-xTB。298.15 K・1 atm、順位の量は δG_eff(kcal/mol)。
- run: `/home/user/hfauto_r7/<wave>/<run>` に新しい run dir で 1 本ずつ直列(`/home/user/hfauto_r7/runcase.sh`。`/usr/bin/time -v`、終了後の孤児プロセス検査、status と report)。記録の型が変わったので round 6 の run は再開せず、VAL との比較は生の JSON と xyz だけを読む。QM 時間は一意のジョブ鍵ごとの秒数(`/home/user/hfauto_r7/tools/jobtime.py`)。GNU timeout は負荷の下で約 5.5% 遅れて発火する。
- 一時 pipeline(`/home/user/hfauto_r7/pipelines/`):
  - `known_endpoints_noscreen`: known_endpoints から reaction-paths の `screen:` を外したもの(行 12 がなく、FIND_PATH が string を走らせる)。
  - `known_endpoints_walltime`: known_endpoints の reaction-paths に `policy: {walltime_h: 0.01}` を足したもの(行 7 の予算の仕組みだけを見る)。
- 合否: 終了コード、outcome、⟨S²⟩、既知の数値(VAL)からの差 0.01 kcal/mol 以内(意図した変化を除く)。新しい系の期待値は、DFT プローブ(`/home/user/hfauto_r7/W1/probe`、pipeline と同じ Level)で状態を確かめてから system の見出しに書いた。期待の分岐に入らなかった run は未達として記録し、期待や閾値は動かしていない。

## 2. 実計算で通した範囲

| 項目 | DFT(PBE0/SVPD)まで通したもの | 低レベルだけ |
|---|---|---|
| 元素 | H、C、N、O、F、Cl、I(ECP)。Te(ECP、VAL7/s3_h2te) | S(SO2·NMe3)。Fe(FeCl3·CH4、VAL7/s9_conformers) |
| 電荷 | 0、−1(Cl⁻・I⁻ の SN2) | — |
| スピン | 一重項、二重項(CH3O•、H + H2、OH···CH4)、三重項(O2) | — |
| 反応型 | 1,2-H 移動(HCN、CH3O•、HCOH → H2CO)、ねじれ(HONO、シュウ酸の 2 段)、1,3-H 移動、縮退転位(NH3、DME、H2O2、PT、恒等 SN2、H 交換、(HF)₂)、H 引き抜き(OH···CH4 → CH3···H2O)、障壁なし(H2O·HF の反転、水二量体の受容体交換)、会合(TMA·(HF)₂) | ラジカル会合(CH3 + O2 → CH3OO•、UKS 参照がスピン汚染で結論なし)、ハロゲン結合、配位付加体 |

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

宣言反応で行 1 を通す系はない。両端の多重度が違う宣言は structures が拒否する(`/home/user/hfauto_r7/W1/probe/ch3n_spin`: `endpoints differ in charge, multiplicity or atom order`)。

### 3.4 discover の系(W3、VAL7)

VAL7 では S7・S8(生成物 0、試行数 2・19)、S10(32.2892)、S19(`collapsed_at_dft_from_seed`、ΔG_assoc −11.5288)が W3・W4 と同じ。S5 と S6 は平らな PES の上で NWChem の軌跡が run ごとに分かれ、次のように変わった。

- S5(`/home/user/hfauto_r7/VAL7/s5_ch3_o2`): CH3···O2 → CH3OO• の case は 2 回目の再開が DFT の Hessian(重なり 0.21)で収束し、2 つの CH3···O2 vdW 極小を結ぶ −72.0i の鞍点で reassigned(`spin_contaminated` で順位なし)。会合の結論がないことは W3 と同じ。
- S6(`/home/user/hfauto_r7/VAL7/s6_oh_ch4`): 同じジョブ鍵の CH3OH + H の opt が W3 では肩(−6.5i の雑音モード、0.10 kcal/mol 上)に止まり、VAL7 では真の極小に入って 2 つの発見が 1 つの basin にまとまった(W3 の multi_step は reassigned 30.659 に、発見の仮説は 1 件減った)。R6 の仮説 OH···CH4 → CH3···H2O は screen の頂点から −73.0i の鞍点に収束したが、QRC の両側が 2 回とも OH···CH4 の basin に戻り(行 8 の初めての到達)、saddle の試行を 1 回残したまま `connection_failed` で終わった。この打ち切りを直し(§7)、同じ run を paths から流し直した `/home/user/hfauto_r7/W5/s6_oh_ch4_rowfix` では FIND_PATH → `path_hei`(maxiter)→ 再開 → −58.6i の鞍点も両側が同じ basin で、試行を使い切って unresolved。W3 の −481.1i の引き抜きの TS(δG_eff 1.48)は VAL7・W5 では再現していない。

| # | 系 | 実測 | 時間 | 判定 |
|---|---|---|---|---|
| S5 | CH3• + O2 | CREST rc −11 で配置の seed を使う。失われた seed の状態を DFT に問うと CH3···O2 の vdW 極小(C···O 2.84 Å)が残るが ⟨S²⟩ 1.71(spin_contaminated)。CH3···O2 → CH3OO• の case は saddle の maxiter 4 回で attempts_exhausted。別の発見の case は reassigned(32.67) | 46:49 | 結論なし(多参照性のある二重項の結合) |
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

`log.jsonl` の decide の reason と注記、minima の diagnostics.json から数えた。run dir は `/home/user/hfauto_r7/` 以下。

| 対象 | 通ったもの | 通っていないもの |
|---|---|---|
| 決定表の行 | 8(VAL7/s6_oh_ch4 の `connection_retry` → `connection_failed`、W5/s6_oh_ch4_rowfix の `connection_retry` → FIND_PATH)、2(W0〜W4/water_same_basin、s19)、3(W1/formaldehyde)、4(R1〜R3 ほか多数)、5(W1・W2・W4/oxalic_two_step、W3/s6_oh_ch4)、6(W1/h2o_hf_inversion、W1/water_dimer_as)、7(W1/hono_walltime_b)、9・10・12・14(R1〜R3 ほか多数)、13(W3/s6_oh_ch4)、15(W1/hcn_noscreen、W1/hf_dimer_swap_noscreen、W3/s5_ch3_o2、W3/s6_oh_ch4)、17(W3/s5_ch3_o2) | 1(端点の DFT 極小が欠けたときだけ発火)、11(崩壊した鞍点)、16(string の次のチャンク)、近道の種が失敗した後の 2 回目の SCREEN、ケースの例外の閉じ込め |
| 検証済みの鞍点(`ts_calc`) | W1・W2/nh3_planar_seed、W1・W2・W4/oxalic_two_step(split1)、W3/s6_oh_ch4(split1) | ゲートを通らず SCREEN へ戻る分岐 |
| saddle の種 | `screen_ts`(R1〜R3 ほか多数)、`screen_hei`(W2/s13_h2o2_gauche、W2/s18_h3_doublet、W3/s5、W3/s6)、`discovery_ts`(W2〜W4/s10、W3/s5、W3/s6)、`path_hei`(W1 の noscreen 2 本、W3/s5、W3/s6)、`higher_order_retry`(W3/s6 の split2_split2、Hessian は検証済み TS freq)、`saddle_restart`(W3/s5 に 2 回、W3/s6 に 2 回) | — |
| 初期 Hessian | xTB(65 case)、DFT 5 case(W1 の noscreen 2 本、W2/s13、W2/s18、W3/s6)、TS freq 1(W3/s6) | — |
| outcome | elementary、degenerate、same_basin、out_of_window(W1/formaldehyde)、multi_step(oxalic、S6)、barrierless_at_resolution(W1/h2o_hf_inversion、water_dimer_as)、reassigned(W3/s5、W3/s6)、unresolved_within_budget(W1/hono_walltime_b、W3/s5) | blocked_upstream |
| QRC | `minus_is_image`(W2 の NH3、nh3_planar_seed、dme_c2v_seed、S1、S2、S10、S12〜S14、S17、S18。W3・W4 の再実行も)、`end<i>_to_new_basin`(oxalic、S6) | — |
| 極小 | DFT の mode-follow `ts_candidate`(W1・W2/nh3_planar_seed、W1・W2/s4、W3/s6)と `one_side`(W1・W2/dme_c2v_seed)、`image_of`(W2 以後の NH3、nh3_planar_seed、S1、S2、S10、S13、S14、S16〜S18)、`not_reacting`(W2〜W4/s10)、`collapsed_at_dft_from_seed`(W3・W4/s19)、screen の `soft:resolved`(s19) | DFT の soft_minimum |
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
| R17 開殻 DFT の SCF 救済 | cgmin の出力は ⟨S²⟩ を出さず、救済したジョブは必ず捨てられる(`/home/user/hfauto_r6/SA/sac_probe/scf/h3p3c/a00_bead_000003_cgmin/stdout.txt`)。実 run での発生は 0 | 開殻 DFT は救済しない(round 7 の run でも発生なし) |
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

1. **実計算で通っていない分岐**: 行 1・11・16、2 回目の SCREEN、ケースの例外の閉じ込め、DFT の soft_minimum。行 6 の中点による高密度化は対称な経路でしか通っておらず、隠れた障壁を見つけた例はない。
2. **既定の順位の精度**: PBE0/SVPD は CCSD(T) から最大 3 kcal/mol ずれ(SN2 −3.0、HONO +2.0)、H 引き抜きは分離基準の ΔE‡ 0.79 と文献 約 5 を大きく下回る。反応クラスが混ざる順位は method_panel を追記して読む。
3. **開殻の会合**: S5 の UKS CH3···O2 は ⟨S²⟩ 1.71 で、障壁なしかどうかを判定できない。多参照かスピン射影の扱いが要る。
4. **失われた seed の費用**: S19 の seed は DFT の 1 歩目で状態を離れたのに、既知の basin に入るまで 44 歩(約 33 分)払った。
5. **鞍点探索の空回り**: 柔らかい H 引き抜き・ラジカル会合では string の前に saddle が maxiter で 2 回止まる(S5 で saddle 826 s)。
6. **対称な鞍点の mode-follow**: ± の側は厳密な像なのに 2 本とも opt と freq を払う(S6 で 329 s、NH3 は 21 歩/側)。
7. **近直線の錯体の熱化学**: S18 の H···H2 は、軸から 0.006 Å 外れた極小が非直線(振動 3 本、δG_eff 4.091)、正確に直線の極小が直線(4 本、6.360)と扱われる(`vibrations._EXTERNAL_RANK_TOL` = 1e-3、`/home/user/hfauto_r7/W2/superseded_k1neg/s18_h3_doublet`)。直線性を対称性で決める必要がある。
8. **分割の深さ**: W3 の S6 の CH3OH + H 側は、ほぼ縮退した緩い錯体(ΔE_rxn ≤ 0.11)の分割で `max_split_depth` に達し、順位が付かない(VAL7 ではその 2 つの極小が 1 つの basin にまとまった)。
9. **平らな PES の引き抜きの TS の再現性**: S6 の OH···CH4 → CH3···H2O は W3 だけが −481.1i の TS に届き、VAL7・W5 は OH···CH4 の向きを変える鞍点(−73.0i、−58.6i、両側が同じ basin)に収束して試行を使い切った。虚モードが仮説の結合変化を担わない鞍点を QRC の前に見分け、結合変化の方向に拘束して探す必要がある(その鞍点の QRC は 2 振幅 × 2 本)。
10. **NWChem の AUTOZ の失敗**: VAL7/s10 の QRC の側の opt が 21 歩目で内部座標の再構築に失敗し、続きからの 2 回目の試行で収束した(+53 s)。VAL7 で 2 回目の試行はこの 1 件。

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
