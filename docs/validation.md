# 実計算による検証(Wave 7、2026-09-25)

- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、GoodVibes 4.3.0)。
- 手法: DFT は `pbe0-d3bj_def2-svpd`(grid fine、SCF 1e-7)、低レベルは GFN2-xTB。設定は `configs/` のまま使い、amine パイロットと A/B だけ一時 system / pipeline を使った。
- run の置き場所: WSL の ext4 上の `/home/user/hfauto_w7/<run>`(下表の run 列はこの `<run>`)。リポジトリの `runs/`(過去の run)には書き込まない。
- コード: `6b87814`。TMA·(HF)₂ の paths 以降と A/B は、修正 F1(§5。W7 のコミットに含まれる)を入れた後に実行した(A/B は F1 の影響を受けない)。`<run>/resolved_config.yaml` の `code_version` に付く `-dirty` は、WSL の git が Windows 側の改行(CRLF)を差分と見なすためで、改行以外の差分は F1 だけである。
- 実エンジンの smoke(WSL、`HFAUTO_REAL=1 HFAUTO_SITE=configs/sites/wsl_local.yaml pytest -m real tests/smoke`): 10 passed(59 秒)。
- 所要時間は `/usr/bin/time` の壁時計時間。ジョブ数は `hfauto status` の misses + hits。

## 1. ベンチマーク(設計書 §10.3 の 7 項目)

| # | 系 | pipeline・手法 | 実測値 | 期待値 | 判定 | 所要時間 | run | コード |
|---:|---|---|---|---|---|---|---|---|
| 1 | HCN→HNC | known_endpoints、PBE0-D3BJ | TS −1128.5i、ΔE‡ 46.62、ΔG‡ 42.22、ΔG_rxn 12.68、ΔE_rxn 13.18 kcal/mol、`elementary_step`、rankable | −1131i、46.7、42.3、12.7、1 h 未満 | 合格(※1) | 1 分 40 秒 | `hcn_known_endpoints` | 6b87814 |
| 2 | HONO trans→cis | known_endpoints、PBE0-D3BJ | TS −681.5i、ΔE‡ 13.67、ΔG‡ 12.24、ΔG_rxn 0.08 kcal/mol、`elementary_step` | −680i、13.6、12.2 | 合格 | 7 分 22 秒 | `hono_known_endpoints` | 6b87814 |
| 3 | NH3 反転 | known_endpoints、PBE0-D3BJ | TS −758.0i、ΔE‡ 4.26、ΔG‡ 3.85 kcal/mol、`degenerate_rearrangement`。`symm` の寄与 +0.411(symm なしの ΔG‡ は 3.44) | −757i、4.25、`symm` の分 +0.41 | 合格 | 2 分 08 秒 | `nh3_inversion_known_endpoints` | 6b87814 |
| 4 | 水(同一 basin) | known_endpoints、PBE0-D3BJ | `same_basin`。reaction-paths のジョブ 0 件 | `same_basin`、TS 探索なし | 合格 | 15 秒 | `water_same_basin_known_endpoints` | 6b87814 |
| 5 | TMA·(HF)₂ | discover、GFN2 → PBE0-D3BJ | `same_basin`。reaction-paths stage は 0 秒・ジョブ 0 件 | SCREEN で終わる(barrierless / no_product / same_basin)、paths は 1 h 未満 | 合格(F1 の後。F1 の前は `blocked_upstream`) | 全体 52 分 33 秒(うち dft stage 50 分 50 秒)、paths 0 秒 | `tma_hf2_discover` | 6b87814、paths 以降は + F1 |
| 6 | amine パイロット(NH3·HF、TMA·HF) | discover `--to explore`、GFN2(ReaDuct NT2 / AFIR) | 44 attempt。収束した低レベル TS 19 件はすべて射影後の虚振動 1 本。本数の違いで棄却された attempt は 0 件。生成物 0 件 | TS がすべて虚振動 1 本 | 合格(※2) | 1 分 56 秒 | `amine_pilot2_discover` | 6b87814 |
| 7 | xTB 初期 Hessian の A/B | minima(dft)、PBE0-D3BJ、初期 Hessian は GFN2 | opt の点数 27(A)/ 40(B)、opt の時間 917 s(A)/ 1,410 s(B)。詳細は §4 | 2 断片以上の DFT opt でステップ数を比べる | 記録済み | A 55 分 18 秒、B 1 時間 4 分 | `ab_init_hessian`、`ab_no_init_hessian` | 6b87814 + F1 |

- ※1: 期待値は D3(zero) の旧 run(G06、G07)の値である。現行の method ファイルは D3BJ なので、−3 cm⁻¹、−0.05〜−0.06 kcal/mol ずれる。W5・W6 の D3BJ の run(1128.5i、46.62、42.22)とは一致した。
- ※2: 旧 run の 8 件(G20。うち NH3·HF と TMA·HF の 2 件がこのパイロットの組成)は golden テストで射影後 1 本と確認済み。今回の 19 件の内訳は NH3·HF 3 件(−1265i)、TMA·HF 9 件(−1395i / −1519i)、TMA 単量体 7 件。

## 2. 撤回する値

- HONO trans→cis の ΔG‡ **24.11 kcal/mol**(`hono_isomerization_v3`、CH-38)は撤回する。原因は ZPE の二重計上と basin 判定の回帰である。実測値は **12.24 kcal/mol**(ΔE‡ 13.67)。
- HCN→HNC の ΔG_rxn **18.05 kcal/mol**(内部フォールバックの熱化学、CH-39)は撤回する。正しい値は **12.74 kcal/mol**(D3 zero、G06)で、D3BJ の今回の実測値は 12.68 kcal/mol。

## 3. 陰性結果の要約(気相の amine·(HF)n)

- **n = 1**(NH3·HF、TMA·HF): HF から N へのプロトン移動(N–H 形成 + H–F 切断)は、NT2 が単調上昇(TS なし)、AFIR は出発構造に戻った。NH3·HF の NT2 の TS(−1265i)3 件は、IRC の両端とも出発構造に戻る H 交換だった。
- **n = 2**(TMA·(HF)₂): 入力の neutral・shared_proton と CREST の 2 配座は、xTB で 1 つの basin にまとまった。DFT でも、neutral と AFIR の生成物が同じ極小(−374.7574786 Eh)に落ちた。極性 H の NT2 は、TS 2 件(−289i)の IRC が出発構造につながらず、1 件は虚振動 0 本、1 件は未収束だった。宣言反応 neutral→shared_proton は `same_basin` で、TS 探索は起動しない。旧実装では、この系に約 30 h をかけて TS 0 件だった。
- **n = 3**: 設計時点の結論(障壁なし)を引き継ぐ。W7 では再計算していない(§6)。
- 同じ組成の旧 run で見つかった TS(G20)は、虚振動の数え方の誤りで 2〜4 本と数えられていた。正しく数えると 1 本だが、いずれもプロトン移動の生成物にはつながらない。

## 4. xTB 初期 Hessian の A/B

2 断片の錯体 4 構造を、同じ開始構造(GFN2 で最適化した構造)から DFT で最適化した。A は discover と同じ設定(`init_hessian: {engine: xtb, method: gfn2}`、NWChem `inhess 2`)、B は初期 Hessian なし(NWChem driver の既定)。どちらも一時 pipeline(structures → minima(dft、`include: all`))で実行した。点数は opt の Evidence の `trajectory_energies_hartree` の長さ(開始点を含む)。

| 構造(原子数) | A: 点数 / opt の時間 | B: 点数 / opt の時間 | 最終エネルギーの差 B − A |
|---|---|---|---|
| NH3·HF(6) | 5 / 9 s | 6 / 10 s | +3e-8 Eh |
| TMA·HF(15) | 4 / 90 s | 6 / 148 s | +3e-8 Eh |
| TMA·(HF)₂(17) | 9 / 405 s | 14 / 625 s | +4.5e-7 Eh |
| TMA·(HF)₂ の AFIR 生成物(17) | 9 / 413 s | 14 / 627 s | +1.8e-6 Eh |
| 計 | 27 点 / 917 s | 40 点 / 1,410 s | — |

- 初期 Hessian で点数は 33%、opt の時間は 35% 減った。xTB の Hessian は 1 本 0.3 秒未満である。極小の割付け(TMA·(HF)₂ の 2 構造が同じ basin になること)は A と B で同じだった。
- stage 全体の所要時間(A 55 分 18 秒、B 1 時間 4 分)は、DFT freq(17 原子で 1 本約 15 分)が大半を占める。

## 5. 修正

- **F1**(`hfauto/chemistry/hypotheses.py`): 宣言反応の端点が screen の basin に崩壊した seed(members の 1 つ)の場合、その端点は screen の極小に割り付き、決定表の行 1 で `blocked_upstream`(`endpoints_not_on_one_pes`)になっていた(TMA·(HF)₂ の shared_proton)。§8.2 では崩壊した seed を DFT に流さない。このため、端点の basin がその screen basin の DFT 極小であればそちらを使うように変えた。回帰テストは `tests/unit/chemistry/test_hypotheses.py::test_declared_endpoint_collapsed_at_screen_takes_the_dft_basin`。

## 6. 旧実装との比較と未解決事項

旧実装の値はレビュー(`docs/reviews/2026-09-25_architecture_chemistry_review.md` §3.0)による。新実装の run では、どの run でも FailureKind 別の失敗は 0 件、`GEOMETRY_MAXITER` も 0 件だった(W4 で入れた `trust 0.1` は、17 原子の錯体でも 9 点で収束した)。

| 系 | 旧実装(run、所要、結果) | 新実装(所要、ジョブ数(うち再利用)、結果) |
|---|---|---|
| HCN | NEB 1,127 s で未収束。string→saddle 849 + 47 s。IRC 107〜818 s/回で判定誤り。ΔG_rxn 18.05 | 1 分 40 秒、29(2)、elementary |
| HONO | NEB 2.5 h で手動停止。TS と IRC があるのに unresolved、ΔG‡ 24.11 | 7 分 22 秒、29(2)、elementary |
| NH3 | IRC の後退側が step 0 で上昇したのに合格(G11) | 2 分 08 秒、28(1)、degenerate |
| TMA·(HF)₂ | 54 run 以上、約 30 h、検証済み TS 0 件 | 52 分 33 秒、72(3)、same_basin。paths 0 秒 |
| amine パイロット | 15 組成で約 4.5 分(`amine_hf_pilot_hf1_hf3_v3_discovery`)。TS 7 件をすべて「虚振動 2〜4 本」として棄却 | 2 組成で 1 分 56 秒、65(4)。TS 19 件がすべて 1 本 |

未解決事項(型・既定値・判断表は変えていない):

1. minima(dft) の freq の省略(§8.2「既存の極小と same なら freq を省く」)は、stage の開始時点の Registry に対してしか効かない。登録が全 relax の後に species id 順で行われるためである。TMA·(HF)₂ では、neutral と同じ極小に落ちる AFIR 生成物にも 17 原子の freq(約 15 分)が走り、dft stage の 30% を占めた。直列の DFT で 1 件ずつ登録するかどうかは、設計の判断として残す。
2. 17 原子・PBE0/def2-SVPD の解析 Hessian は 4 rank で 1 本約 15 分かかり、discover の所要時間の大半を占める。
3. amine·(HF)₃ の「障壁なし」は、W7 では再計算していない(パイロットは設計どおり 2 組成に限った)。
4. 設計(WP7.1)では run を `runs/validation_<system>_<date>` に置くことになっているが、過去の run を守るため、WSL の ext4 上の `/home/user/hfauto_w7/` に置いた。
