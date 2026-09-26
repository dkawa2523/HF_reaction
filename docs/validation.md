# 実計算による検証(v2: 改良後の再検証、2026-09-26)

- 環境: WSL2 Ubuntu(4 vCPU / 11 GB)、`configs/sites/wsl_local.yaml`(NWChem 7.2.3 を 4 rank × 1,200 MB、xTB 6.7.1、CREST 3.0.2、SCINE ReaDuct 6.1.0、pysisyphus 1.0、GoodVibes 4.3.0)。
- 手法: DFT は `pbe0-d3bj_def2-svpd`(grid fine、SCF 1e-7)、低レベルは GFN2-xTB。設定は `configs/` のまま使い、amine パイロット・(HF)₃ の配座・FIND_PATH・A/B だけ一時 system / pipeline(`/home/user/hfauto_v2/configs/`)を使った。
- コード: `4e972c5`(計算仕様レビューの改良 W1〜W4: `116ee0e`、`9b61e8d`、`7b9a192`、`4e972c5`)。今回の再検証でコードの修正は要らなかった。
- run の置き場所: WSL の ext4 上の `/home/user/hfauto_v2/<run>`。改良前(W7、コード `6b87814`)の run は `/home/user/hfauto_w7/<run>` に残した。リポジトリの `runs/`(過去の run)には書き込まない。
- 実エンジンの smoke(`HFAUTO_REAL=1 HFAUTO_SITE=configs/sites/wsl_local.yaml pytest -m real tests/smoke`): 12 passed(72 秒)。W7 の 10 件に、改良で加えた 2 件(PBE0/STO-3G の HCN string で bead エネルギーが落ち着き `single_max` になること 14 秒、F⁻·(HF)₂ の組成探索が `--noopt` で通ること 1.4 秒)が加わった。
- 所要時間は `/usr/bin/time` の壁時計時間。ジョブ数は `hfauto status` の misses + hits(括弧内は再利用)。ΔG は 298.15 K、1 atm。

## 1. ベンチマーク(v2)

| # | 系 | pipeline・手法 | v2 の実測値 | 判定 | 所要時間 v2(W7) | run |
|---:|---|---|---|---|---|---|
| 1 | HCN→HNC | known_endpoints | TS −1128.5i、ΔE‡ 46.62、ΔG‡ 42.18、ΔG_rxn 12.68、ΔE_rxn 13.18 kcal/mol、`elementary_step` | 合格 | 1 分 41 秒(1 分 40 秒) | `hcn_known_endpoints` |
| 2 | HONO trans→cis | known_endpoints | TS −678.1i、ΔE‡ 13.67、ΔG‡ 12.22、ΔG_rxn 0.08 kcal/mol、`elementary_step`。SCREEN の GS が失敗し、FIND_PATH(DFT string)の経路で到達した(§3) | 合格 | 9 分 32 秒(7 分 22 秒) | `hono_known_endpoints` |
| 3 | NH3 反転 | known_endpoints | TS −758.0i、ΔE‡ 4.26、ΔG‡ 3.84 kcal/mol、`degenerate_rearrangement` | 合格(※1) | 1 分 55 秒(2 分 08 秒) | `nh3_inversion_known_endpoints` |
| 4 | 水(同一 basin) | known_endpoints | `same_basin`。reaction-paths のジョブ 0 件。dft の freq は 1 本(W7 は 2 本) | 合格 | 11 秒(15 秒) | `water_same_basin_known_endpoints` |
| 5 | TMA·(HF)₂ | discover | `same_basin`(paths 0 秒・0 件)。explore の生成物 0 件(W7 は偽の AFIR 生成物 1 件)。ΔG_assoc −11.53 kcal/mol | 合格 | 31 分 15 秒(52 分 33 秒)。dft stage 1,771 秒(3,050 秒) | `tma_hf2_discover` |
| 6 | amine パイロット(NH3·HF、TMA·HF) | discover `--to explore` | 40 attempt。収束した低レベル TS 22 件はすべて射影後の虚振動 1 本。生成物 0 件 | 合格 | 1 分 54 秒(1 分 56 秒) | `amine_pilot2_discover` |
| 7 | xTB 初期 Hessian の A/B | minima(dft)、初期 Hessian は GFN2 | A の opt 27 点・899 s(W7 と同じ 27 点)。freq 4 → 3 本(§6) | 記録済み | A: 40 分 34 秒(55 分 18 秒) | `ab_init_hessian` |
| 8 | FIND_PATH(新規) | known_endpoints から screen を外した一時 pipeline、PBE0/STO-3G | HCN→HNC。IDPP から string 2 チャンクで bead エネルギーが落ち着き `single_max`。HEI から鞍点 −1223.9i、`elementary_step` | 合格 | 1 分 03 秒 | `find_path_sto3g` |
| 9 | (HF)₃ の配座(新規) | discover `--to conformers` | NH3·(HF)₃ 7 配座、TMA·(HF)₃ 9 配座(最安はイオン対 TMAH⁺·F⁻·(HF)₂)。CREST はどちらも rc 0 | 合格 | 34 秒 | `amine_hf3_conformers` |

- ※1: `symm` の寄与 +0.411 kcal/mol は外部回転の対称数の比(NH3 の σ 3、平面の TS の σ 6。RT ln 2)による。ΔG‡ は経路縮重度 L を含まない(design.md §7)ので、L = 1 と数えた有効障壁は `symm` を除いた 3.43 kcal/mol(W7 は 3.44)に当たる。
- 検証の基準値(W7 の D3BJ の値)と、旧 run(D3 zero、G06・G07)との差は W7 で確認済みである。

## 2. 改良前後の比較(W7 → v2)

| 系 | 所要時間 | ジョブ数(再利用) | 主な差と原因 |
|---|---|---|---|
| HCN | 1:40 → 1:41 | 29(2)→ 29(2) | QRC は 19 / 24 → 24 / 20 ステップで、S7 の効果はほぼない |
| HONO | 7:22 → 9:32 | 29(2)→ 17(0) | SCREEN の GS が発散して失敗し、DFT string(200 秒)と鞍点(NWChem の AUTOZ 失敗からの再試行を含め 106 秒)が加わった。QRC は 127(maxiter で継続)/ 17 = 144 ステップ・311 秒 → 34 / 36 = 70 ステップ・162 秒(S7) |
| NH3 | 2:08 → 1:55 | 28(1)→ 27(1) | QRC 50 / 50 → 45 / 45 ステップ(S7)。dft の freq 2 → 1 本(S1: nh3_up が既知の basin に入る) |
| 水 | 0:15 → 0:11 | 5 → 4 | dft の freq 2 → 1 本(S1) |
| TMA·(HF)₂ | 52:33 → 31:15(−41%) | 72(3)→ 60(2) | dft stage 3,050 → 1,771 秒(−42%): opt 4 → 3、freq 4 → 3 本(M3 で偽の AFIR 生成物が消えた)。screen の xTB freq 7 → 3 本(S1)。explore 37 → 31 attempt |
| amine パイロット | 1:56 → 1:54 | 65(4)→ 56(2) | attempt 44 → 40(TMA 単量体の trial が 14 → 11。出発構造が xTB の最適化構造になったため、M3)。TS 19 → 22 件、すべて虚振動 1 本 |
| smoke | 59 秒 → 72 秒 | 10 → 12 件 | 新しい 2 件の分 |

ΔG の変化(M6・S15: 振動数と ZPE のスケールを 0.989 / 0.975 から 1.0 / 1.0 に 1 本化)は、予測(ΔG‡ で 0.04 以下、ΔG_assoc で 0.07 以下)の範囲に収まった。

| 量 | W7 | v2 | 差 |
|---|---:|---:|---:|
| HCN ΔG‡ / ΔZPE‡ | 42.22 / −3.34 | 42.18 / −3.42 | −0.04 / −0.09 |
| HCN ΔG_rxn | 12.68 | 12.68 | 0.00 |
| HONO ΔG‡ / ΔG_rxn | 12.24 / 0.08 | 12.22 / 0.08 | −0.02 / 0.00 |
| NH3 ΔG‡ | 3.85 | 3.84 | −0.01 |
| TMA·(HF)₂ ΔG_assoc | −11.59 | −11.53 | +0.07 |

HONO の TS は FIND_PATH の HEI から精密化したので虚振動が −681.5i から −678.1i に変わったが、エネルギーは ΔE‡ 13.67 で W7 と一致した。旧い熱化学と NWChem optimize のキャッシュはキーが変わったので再利用されない。

## 3. 改良項目ごとの実計算での確認

- **M2・S5(FIND_PATH)**: 初めて実計算で通った。HONO では GS の失敗から行 16 の FIND_PATH に入り、IDPP から 1 チャンク(20 反復)で bead エネルギーが落ち着いた(最後の 3 反復の HEI の変化 0.004 kcal/mol)。NWChem は未収束を報告したが、形は `single_max`(HEI は trans から 34.4 kcal/mol で、鞍点の 13.7 より高い上界)で、HEI の種から正しい TS に到達した。HCN(STO-3G、#8)でも 2 チャンクで `single_max` になった。
- **M3・S6(explore)**: TMA·(HF)₂ の neutral からの AFIR は、W7 の「未収束 7・出発に戻る 2・生成物 1」から「未収束 2・出発に戻る 8・生成物 0」になった。W7 で `irc_not_connected_to_source` だった −289i の TS は、グラフの一致で `same_as_source` と判定された。
- **M4・S13・C7(conformers)**: tma_hf2 は W7 と同じ 2 配座(エネルギー差 2e-5 Eh 以内)に、1.48 kcal/mol 上の 1 配座が加わった(3 配座とも screen で同じ basin に入る)。(HF)₃ は、レビューのプローブ(`--noopt` なし)で失敗していた NH3·(HF)₃ を含め 2 組成とも成功し、TMA·(HF)₃ の最安はイオン対だった(#9)。
- **S1(その場の登録)**: 1 basin に freq 1 本になった。screen の xTB freq(TMA·(HF)₂ 7 → 3、amine 7 → 5)、dft の DFT freq(NH3・水 2 → 1)。A/B(§6)では AFIR 生成物の構造が既知の basin に入り、17 原子の freq(約 15 分)を省いた。
- **S7・C19(QRC)**: 3 反応とも `trajectory_above_ts` は出ず、接続の判定は W7 と同じだった。HONO の maxiter での継続はなくなった。
- **S4(GS の回収)**: HONO で GS 自体が例外で終わったため、GS だけの再実行も同じ例外で終わり、`nonzero_exit` として SCREEN が unavailable になった(GS の失敗は握りつぶさない、の設計どおり)。
- 実計算で通っていない改良: M1(陰性結果ゲートの削除)、S3(行 14 の 15 bead 確認の撤廃)は、今回の系に単調な経路や陰性結果で止まる宣言反応がなく、経路に入らなかった。M5・S8(手法パネル)は `method_panel` を実行していない。

## 4. 撤回する値

- HONO trans→cis の ΔG‡ **24.11 kcal/mol**(`hono_isomerization_v3`、CH-38)は撤回する。原因は ZPE の二重計上と basin 判定の回帰である。実測値は **12.22 kcal/mol**(v2、ΔE‡ 13.67。W7 は 12.24)。
- HCN→HNC の ΔG_rxn **18.05 kcal/mol**(内部フォールバックの熱化学、CH-39)は撤回する。正しい値は **12.74 kcal/mol**(D3 zero、G06)で、D3BJ の実測値は 12.68 kcal/mol(W7・v2)。

## 5. 陰性結果の要約(気相の amine·(HF)n)

- **n = 1**(NH3·HF、TMA·HF): HF から N へのプロトン移動は、NT2 が単調上昇(TS なし)、AFIR は出発構造に戻った。NH3·HF の NT2 の TS(−1265i)3 件は、IRC の両端とも出発構造に戻る H 交換だった。v2 でも同じ。
- **n = 2**(TMA·(HF)₂): 入力の neutral・shared_proton と CREST の 3 配座は、xTB で 1 つの basin にまとまった。DFT の極小は neutral の 1 つ(−374.7574786 Eh、W7 と同じ)。宣言反応 neutral→shared_proton は `same_basin` で、TS 探索は起動しない。旧実装では、この系に約 30 h をかけて TS 0 件だった。
- **n = 3**: CREST の配座は得られるようになった(#9。最安は TMA·(HF)₃ でイオン対)。反応探索(障壁なしという設計時点の結論)は再計算していない。
- 同じ組成の旧 run で見つかった TS(G20)は、虚振動の数え方の誤りで 2〜4 本と数えられていた。正しく数えると 1 本だが、いずれもプロトン移動の生成物にはつながらない。

## 6. xTB 初期 Hessian の A/B

2 断片の錯体 4 構造を、同じ開始構造(GFN2 で最適化した構造)から DFT で最適化した。A は discover と同じ設定(`init_hessian: {engine: xtb, method: gfn2}`、NWChem `inhess 2`)、B は初期 Hessian なし。一時 pipeline(structures → minima(dft、`include: all`))で実行した。点数は opt の Evidence の `trajectory_energies_hartree` の長さ(開始点を含む)。B は W7 だけで実行した(改良は B の条件に関係しない)。

| 構造(原子数) | A(W7): 点数 / 時間 | A(v2): 点数 / 時間 | B(W7): 点数 / 時間 |
|---|---|---|---|
| NH3·HF(6) | 5 / 9 s | 5 / 9 s | 6 / 10 s |
| TMA·HF(15) | 4 / 90 s | 4 / 89 s | 6 / 148 s |
| TMA·(HF)₂(17) | 9 / 405 s | 9 / 397 s | 14 / 625 s |
| TMA·(HF)₂ の AFIR 生成物(17) | 9 / 413 s | 9 / 405 s | 14 / 627 s |
| 計 | 27 点 / 917 s | 27 点 / 899 s | 40 点 / 1,410 s |

- W7 では、初期 Hessian で点数は 33%、opt の時間は 35% 減った。xTB の Hessian は 1 本 0.3 秒未満である。
- v2 の A は、opt の点数が W7 と同じで、最終エネルギーも 1e-8 Eh 以内で一致した。AFIR 生成物の構造はその場で TMA·(HF)₂ の basin に入り(members 2)、freq は 4 → 3 本になった。stage 全体の所要時間は 55 分 18 秒 → 40 分 34 秒(−27%)。
- 残りの時間の大半は DFT freq(17 原子で 890 秒、15 原子で 554 秒)である。

## 7. 修正

- **F1**(W7、`hfauto/chemistry/hypotheses.py`): 宣言反応の端点が screen の basin に崩壊した seed(members の 1 つ)の場合、その端点は screen の極小に割り付き、`blocked_upstream`(`endpoints_not_on_one_pes`)になっていた(TMA·(HF)₂ の shared_proton)。端点の basin がその screen basin の DFT 極小であればそちらを使うように変えた。回帰テストは `tests/unit/chemistry/test_hypotheses.py::test_declared_endpoint_collapsed_at_screen_takes_the_dft_basin`。
- v2 の再検証では、コードの修正が必要な不具合は出なかった。

## 8. 未解決事項

1. **HONO の SCREEN(GS)の発散**: pysisyphus の GS(DLC 座標)が約 10 サイクル目から発散し、H–O–N の変角が内部座標から落ちた後、座標数の不一致(66 → 55)で例外終了した。同じジョブを単独で再実行すると、v2 の入力では 2 回とも同じサイクルで同じ例外になり、W7 の入力では収束した。両者の端点の差は 1e-6 Å 程度なので、平面 4 原子での GS の数値的な不安定であり、コードの回帰ではない。結果は FIND_PATH の経路で正しく得られた(+2 分)。再発が続くなら、小さな分子の GS を Cartesian 座標にする、または 1 回だけ再試行する案を検討する。
2. S4 の回収は、TS 最適化ではなく GS 自体が例外で終わった場合にも GS を 1 回再実行する(HONO で 4 秒)。害はないので残した。
3. 17 原子・PBE0/def2-SVPD の解析 Hessian は 4 rank で 1 本約 15 分かかり、discover の所要時間の大半(TMA·(HF)₂ の dft stage 1,771 秒のうち 902 秒)を占める。
4. QRC の trust 0.3: 前提(3 反応で `trajectory_above_ts` が出ない)は満たしたが、trust 0.1 のままにした。
5. 手法パネル(M5・S8)と、陰性結果ゲートの削除(M1)・15 bead 確認の撤廃(S3)の経路は、実計算では通していない(§3)。
6. 設計(WP7.1)では run を `runs/validation_<system>_<date>` に置くことになっているが、過去の run を守るため、WSL の ext4 上に置いた。
