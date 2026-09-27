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
| S5 | CH3• + O2(m=2 を宣言) | M1 で追加 | M1 |
| S6 | OH• + CH4 | `oh_ch4` | M4 |
| S7 | NH3···ICl | `nh3_icl` | M4 |
| S8 | SO2·NMe3 | `so2_nme3` | S-C |
| S9 | FeCl3·CH4 | M1 で追加 | M1 |
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
