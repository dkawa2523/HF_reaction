# hfauto 計算処理の仕様レビュー — 現状・課題・改良案（2026-09-26）

- 対象: ブランチ `refactor/fundamental-2026-09`、HEAD `92ce8d6`。読み取り専用のレビューで、リポジトリのコードと設定は変更していない。
- 範囲: 気相の分子・非共有結合錯体を対象に、「同じ PES 上の別々の極小 → 障壁の事前判定 → 鞍点 → QRC による接続 → 熱化学 → 分類と順位付け」を行う 8 stage と、それを支える実行基盤。
- 方法: 計算単位 U0〜U8 ごとに、次の 3 段で調べた。
  1. 現状仕様の整理。コード、設定、W7 の実計算の生ファイル（WSL `/home/user/hfauto_w7/<run>/jobs/<key>/attempt_*/`）を読んだ。
  2. 評価。化学的妥当性、有用性、複雑さを評価し、各ライブラリの文書とソースを調べた。必要に応じて小さなプローブ計算を行った（WSL `/home/user/hfauto_review_probe/`）。
  3. 反証的な検証。別の担当が 2 の主張をコードと生データで確かめた。
- 記載の規則:
  - 検証で棄却された課題はない。一部だけ正しいと判定された課題は、訂正後の内容で記載した。
  - 重大度と優先度は検証段の判定を使った。
  - 改良案は、検証で「採用」または「修正して採用」とされたものだけを載せた。修正して採用されたものは、修正後の内容を書いた。
  - 棄却された改良案と、評価段で見送った案は、各節の「見送った案」に 1 行の理由付きで載せた。
  - 新しい提案はこのレビューでは加えていない。§3 では、複数の単位で重複した改良案を 1 つにまとめた。

凡例:

| 記号 | 意味 |
|---|---|
| critical / high / medium / low | 課題の重大度 |
| must / should / could | 優先度。must は「結論を誤らせる」または「主要機能が働かない」もの |
| 【複雑さ減】 | コード・設定・分岐が減る改良 |
| 【複雑さ増】 | 増える改良。効果が上回る場合だけ採用した |
| ◎ / ○ / △ / × | 評価。◎ 問題なし、○ 軽微な課題のみ、△ 改良が必要、× 目的どおりに機能していない |

エネルギーの単位は、断りがなければ kcal/mol。

---

## 0. 要旨

### 0.1 全体評価

1. **プロトコルの骨格は化学的に妥当です。** 流れは、GFN2-xTB による探索と経路の事前判定 → PBE0-D3BJ/def2-SVPD（grid fine、SCF 1e-7）の単一 PES 上での極小・鞍点・QRC・解析 Hessian → GoodVibes 4.3.0 の qRRHO、です。これは Bursch/Grimme 2022 の多層プロトコルに沿っています。停留点レベルと CCSD(T)/def2-TZVP の差は次のとおりで、目的（反応の発見と順位付け）には足ります。
   - 障壁: HCN −1.5、HONO +1.7、NH3 −1.2
   - FHF⁻ の De: 47.0（実験値 45.8±1.6）
2. **実計算で確認できたのは 1 経路だけです。** 通った経路は SCREEN → REFINE_SADDLE → VALIDATE_TS → CONNECT で、HCN・HONO・NH3 のすべてが正しい素反応・縮退反応と判定されました。一方、次のフォールバック群は fake エンジンのテストでしか通っていません。重大な欠陥はこちらに集中しています。
   - FIND_PATH（NWChem ZTS）、IDPP、moddir、higher_order、中間体での分割
   - negative_evidence、CREST のトポロジー停止からの再実行
3. **結論を黙って誤らせる欠陥が 5 件あります。**
   - DFT で別 basin の極小がある反応を、低レベルの陰性結果だけで NO_PRODUCT にする（negative_evidence ゲート、critical）。
   - FIND_PATH の収束判定が、障壁のある経路では原理的に成り立たない。そのうえ、収束した string を停滞として kill して捨てる（high）。
   - 手法パネルの符号ゲートに許容幅がなく、検証済みの HONO trans→cis が順位から外れる（high）。
   - SMILES で与えた宣言反応の端点では H の原子写像が決まらず、別の機構を計算する（high）。
   - 熱化学が、極小と TS の小さな負の振動数を捨てる。1 モードあたり 0.5〜1.6 の偏りになる（medium）。
4. **計算時間の大半は DFT の解析 Hessian と QRC の最適化です。** 17 原子の freq は 1 本約 15 分で、dft stage の 68〜71% を占めます。これに次の 3 つの無駄が重なっています。
   - 同じ stage の中で、同じ basin に落ちた 2 本目以降の freq を省けない。
   - explore が screen 極小ではなく入力構造から出発する。その結果、偽の生成物が dft stage の時間の 43% を使った。
   - QRC 側の opt に初期 Hessian を渡していない。HONO では paths stage の 81% がこれだった。
5. **設定の量は、全体としては目的に見合っています。** 利用者が必ず決めるのは pipeline・system・site・run-dir の 4 つだけです。ただし、効かない knob や使われない分岐が各所に残っています。
   - conformers の nci・topology・threads、thermo の temperatures_K・invert_soft_cm1・スケール因子の DB 検索
   - negative_evidence の override、15 bead の確認、ambiguous → tight の再判定
   - 改良の多くは「削る」方向で済みます。
6. **精度の限界は、実装の誤りではなく読み方の問題です。**
   - 順位は「PBE0-D3BJ/def2-SVPD での ΔG‡ の序数」として読んでください。約 1.5 kcal/mol 未満の差は手法誤差の範囲です。
   - ΔG_assoc は非補正の def2-SVPD の値で、数 kcal/mol 過大に結合している可能性があります（NH3·HF で CCSD(T)/aVTZ-CP より 3.7 深い）。
   - 手法パネルの「参照」である CCSD(T)/def2-TZVP には diffuse 関数がないため、陰イオンでは停留点レベルより悪くなります。
7. **大きな追加は、費用に見合わないので見送りました。** 対象は新しいエンジン（ORCA/DLPNO、aISS、MC-AFIR）、停留点レベルの変更（ωB97X-D3、TZ での再最適化）、IRC への置き換え、トンネル補正、hindered rotor、CP 補正の自動化などです。

### 0.2 最優先の改良（上位 5 件）

| # | 改良 | 対象単位 | 期待効果 | 工数 | 複雑さ |
|---|---|---|---|---|---|
| 1 | negative_evidence ゲートを削除する（陰性結果は coverage の報告に残す） | U4-P1 / U5-P1 | DFT で別 basin がある反応や宣言反応（例: HCN→HNC は GFN2 の障壁 306 kJ/mol で窓外になる）を、ジョブなしで NO_PRODUCT にする誤りがなくなる | 約 40 行の削除とテスト 2 件 | 【複雑さ減】 |
| 2 | FIND_PATH の判定を置き換える。gmax による停滞 monitor と、形の記録に付いている収束ゲートを外し、チャンクの継続は bead エネルギーの安定度だけで決める | U5-P2（修正版） | DFT string が初めて機能し、monotonic・multi_max・single_max を返せる。プローブでは 8 反復目以降にエネルギーが安定した | 約 30 行（削除が中心）とテストの差し替え | 【複雑さ減】 |
| 3 | explore の出発構造を screen 極小の最適化構造にする | U4-P2 = U3-IMP-3 | TMA·(HF)₂ の偽 AFIR 生成物（DFT の opt と freq で 1,308 s、dft stage の 43%）が出なくなる。run ごとに結論が変わる問題も解消する | 約 5 行 | 中立 |
| 4 | CREST の組成探索に `--noopt` を付ける | U2-IMP-1 | (HF)₃ 系の成功が 2/3 から 3/3 になる。TMA·(HF)₃ で 1.21 kcal/mol 低いイオン対が見つかる。所要時間は約 1/3 | 1〜5 行 | 【複雑さ減】 |
| 5 | 手法パネルの符号ゲートに ±1 kcal/mol の不感帯を入れ、パネルから def2-SVP を外す | U0-P1 = U8-IMP-1 | 検証済みの HONO trans→cis が rank 1 に戻る。符号が本当に割れている場合の検出は残る | 2 行、テスト 1 件、method ファイル 1 本の削除 | 【複雑さ減】 |

ほかに must が 2 件あります。熱化学の負モードの固定規則（U7-P1 = U3-IMP-4）と、宣言反応の端点に xyz を必須にすること（U1-I1）です。一覧は §4 にあります。

### 0.3 今のままでよい点（主なもの）

- 停留点・freq・QRC・thermo に `same_pes(numerics=True)` を課し、LOT が一致しなければ `mixed_level_of_theory` で順位から外す設計。すべての DFT デッキに `convergence energy 1e-7` を明示している。
- opt と freq を別ジョブにし、振動数は hfauto 自前の質量加重 Eckart 射影で数える。この射影は NWChem の表示より正しい。
- 2 断片以上に限った xTB 初期 Hessian（`inhess 2`）。DFT opt の点数が −33%、時間が −35%。
- SCREEN の骨格（xTB 端点 → pysisyphus GS と TSOpt → 全ノードで DFT SP）。W7 の 4 run すべてで、screen_ts から DFT saddle が一発で収束した。
- QRC の振幅を曲率から決め、接続ゲートを単調性ではなくエネルギー降下で判定する方式。golden の G09/G11/G12 で、登っていく枝を正しく棄却する。
- GoodVibes の API に、自前の振動数・質量・回転定数とスケール因子を明示して渡す方式。thermo_consistent ゲート、fail-closed、pymsym の σ（QCData 側は symmno=1）。
- JobStore のキーに ExecutionSpec を含めず、版数 pin を含める。W7 から HEAD まで 29/29 件を再利用できた。ladder は 4 種に限り、SiteLock と `/mnt` 上の scratch の拒否もある。
- CREST の `--scratch` なし実行、版数照合、組成への `--notopo` の自動付与、placement の決定論的な seed。
- rankable を順位の唯一の関門とし、信頼度スコアを作らない方針（AR-27）。coverage.csv で陰性結果を報告すること。

---

## 1. 計算処理の一覧

評価は「妥当性 / 有用性 / 複雑さ」の順に書いた。複雑さの ◎ は「単純で過不足がない」を意味する。

| 単位 | 目的 | 使用エンジン・手法 | 主要既定値 | 実測コスト（W7 と本レビューのプローブ、WSL 4 vCPU） | 評価 |
|---|---|---|---|---|---|
| U0 全体プロトコル | 理論レベルの層構成を決め、ΔE‡・ΔG‡・ΔG_rxn・ΔG_assoc と順位を出す | GFN2-xTB / PBE0-D3BJ/def2-SVPD（NWChem 7.2.3）/ GoodVibes 4.3.0。任意で手法パネル（PBE0/SVP・TZVP、ωB97X-D3/TZVPD、CCSD(T)/TZVP） | grid fine、SCF 1e-7、D3BJ（`disp vdw 4`）、qRRHO grimme 100 cm⁻¹、298.15 K、1 atm | HCN 1 分 40 秒、HONO 7 分 22 秒、NH3 2 分 08 秒、TMA·(HF)₂ discover 52 分 33 秒（うち freq 2,070 s）。パネルは HCN の 12 SP で 63.7 s | ○ / ○ / △ |
| U1 structures | system の化学種を 3D 構造と SpeciesRecord（電荷・多重度・状態ラベル）にする | RDKit 2026.03.3 ETKDGv3（SMILES のみ）、外部 QM なし | seed 20260925、1 配座、力場による最適化なし。偶奇の検査と元素列の一致 | 1 s 未満、ジョブ 0。埋め込みは 0.1〜30 ms | ○ / ◎ / ◎ |
| U2 conformers | 単量体の配座と、錯体の配置・PT 状態を前探索する | CREST 3.0.2 iMTD（`--gfn2 --quick`、組成は `--nci`）、H 結合円錐の placement | `--ewin 6`、`-T 4`、窓 4 kcal/mol、keep 6 / 状態、seed 6（CREST に渡すのは seed00） | stage 12〜17 s、CREST 1 回 1.2〜10.9 s。(HF)₃ 系のプローブは 87〜236 s（NH3·(HF)₃ は失敗） | △ / ○ / ○ |
| U3 minima | 候補を同じ PES の極小に確定し、basin にまとめる（screen / dft） | xTB `--opt vtight` と `--hess` / NWChem opt と別ジョブの解析 freq / Registry | 虚振動 −10 / −50 cm⁻¹、mode-follow 2 回・±0.1 Å、同一性 0.05 Å・5e-5 Eh、窓 6 kcal/mol・per_state 3・DFT SP で上位 8 を再順位付け | screen 1 s 未満。dft は HCN 18 s、HONO 54 s、TMA·(HF)₂ 3,050 s（freq 68%、17 原子の freq 1 本約 900 s） | ○ / ○ / △ |
| U4 explore | screen 極小から単端の反応 trial を回し、生成物と低レベル TS を得る | SCINE ReaDuct 6.1.0（NT2 → Bofill TSOpt → IRC、次に AFIR）、GFN2 | 出発点 2 / 状態、trial 10 / 出発点、窓 ΔE‡ 150・ΔE_rxn 100 kJ/mol、AFIR γ 125→300 | stage 83〜102 s、1 attempt 0.5〜8 s。偽の生成物 1 件が DFT 1,308 s を消費 | △ / △ / ○ |
| U5 paths（仮説・SCREEN・FIND_PATH） | 反応仮説を作り、DFT saddle 探索の前に障壁を事前判定する。種がなければ DFT string を使う | hypotheses、xTB + pysisyphus 1.0 GS（11 ノード）+ DFT SP、IDPP、NWChem ZTS | 窓 40 kcal/mol、変化 ≥0.2 Å か ≥30°、判定 1 kcal/mol、string 9 bead × 20 反復 × 3 チャンク | SCREEN は HCN 28 s、HONO 48 s、NH3 29 s。FIND_PATH は未実行（17 原子で 1 チャンク約 1 h と推定） | SCREEN ○ / FIND_PATH × / △ |
| U6 paths（鞍点・検証・接続） | 種を 1 次鞍点に精密化し、TS を検証し、QRC で接続を確定する | NWChem `task dft saddle`（xTB Hessian を `inhess 2` で渡す）、別ジョブの DFT freq、QRC の ±opt | saddle は maxiter 50・trust 0.1。TS の閾値 −50（ねじれは −20）cm⁻¹、prominence 2e-5 Eh、QRC 目標 3e-4 Eh・0.05〜0.4 Å | HCN 50.6 s、HONO 336 s（QRC が 81%）、NH3 78 s | ○ / ◎ / ○ |
| U7 thermo | 極小と TS の G(T) を求め、反応量・会合量・感度幅・分布を出す | GoodVibes 4.3.0 API（worker）、pymsym | qRRHO grimme 100 cm⁻¹、symm、スケール 0.989/0.975（PBE0/MG3S の借用）、感度 6 variant、1 atm | stage 0〜1 s、1 subject 0.22〜0.32 s | ○ / ○ / △ |
| U8 sp・report・実行基盤 | 固定構造での多手法 SP（パネル・composite）、順位付けと報告、ジョブ実行・再開 | NWChem DFT / MP2 / CCSD(T)、reporting、JobStore / ladder / SiteLock | 既定 targets は reaction_stationary_points、符号ゲートの許容幅 0、timeout 14,400 s、ladder 4 種 | report 1 s 未満。HONO の CCSD(T)/TZVP は 1 点 282〜322 s（負荷時）。amine·(HF)n の CCSD(T) は 1 点数時間と推定 | ○ / ○ / ○ |

---

## 2. 単位ごとの詳細

### 2.0 U0 全体プロトコルと理論レベルの選択

#### 2.0.1 現状仕様

**流れ**

1. **設定の解決**（`pipeline/config.py:93-138`）
   - stage 設定の `method`・`energy_method`・`methods` から method id を集め、`configs/methods/<id>.yaml` を `MethodSpec`（frozen、extra=forbid）として読む。
   - 解決結果は `<run>/resolved_config.yaml` に pipeline_id ごとに記録する。
2. **preflight**: エンジンと手法の対応を検査する。
   - NWChem: `dft`（d4 は不可、溶媒は `cosmo:<eps>` だけ）と、energy に限った `wft`。
   - xTB・CREST・pysis_gs: `gfn` と任意の `alpb`。
   - ReaDuct: 溶媒なしのときだけ受け付ける。
3. **低レベル（discover だけ）**: すべて GFN2-xTB で行う（電子温度は未指定なので 300 K、溶媒なし）。
   - CREST による配座探索
   - screen 極小（`--opt vtight` → `--hess`）
   - ReaDuct の NT2/AFIR
   - pysisyphus の GS
4. **停留点レベル**: PBE0-D3BJ/def2-SVPD。次をすべてこの単一 PES 上で計算する。
   - minima(dft) の opt と freq（2 断片以上では xTB Hessian を `inhess 2` で渡す）
   - discover では、選択の前に DFT SP で上位 8 構造を並べ直す
   - REFINE_SADDLE、VALIDATE_TS、QRC、ZTS string、SCREEN の DFT SP
5. **thermo**
   - GoodVibes qRRHO（qs grimme、cutoff 100 cm⁻¹、symm）。
   - 感度幅として qs{grimme, truhlar} × cutoff{50, 100, 150} も計算する。
   - スケール因子は表から引く。PBE0/def2 の項目はないので、PBE0/MG3S の 0.989/0.975 を使う（注記 `scale_factor_from:PBE0/MG3S`）。
   - `energy_method` が None（同梱の全 pipeline）なら、エネルギーは freq ジョブの E を使う。
6. **method_panel（任意、既存 run に追記）**
   - 反応の極小と TS の freq 構造で、4 手法の SP を直列に計算する。
   - report で ΔE_rxn・ΔE‡ の表を作る。ΔE_rxn の符号が割れたら blocker `method_sign_disagreement` を付ける。
   - ΔG の値そのものは変えない。
7. **composite G**（G = E_SP + (G_GV − E_GV)）: thermo の `energy_method` として実装されているが、同梱の pipeline では使っていない。

**エンジン呼び出しの実パラメータ（W7 の生ファイルから）**

- NWChem DFT opt（hcn_known_endpoints の `jobs/24/24a9…/attempt_00/job.nw`）。argv は `mpirun -np 4 nwchem job.nw`、site は 4 rank × 1200 MB、timeout 14,400 s。
  ```
  memory total 1200 mb
  geometry units angstrom nocenter noautosym … end
  charge 0
  basis spherical ; * library def2-svpd ; end
  dft ; xc pbe0 ; mult 1 ; grid fine ; convergence energy 1.0e-07 ; disp vdw 4 ; end
  driver ; maxiter 100 ; trust 0.1 ; xyz final ; end
  task dft optimize
  ```
  2 断片以上（NH3·HF など）では、driver に `inhess 2` が入り、`job.hess`（xTB Hessian を変換したもの）を読む。
- NWChem freq: 上と同じ dft ブロックで `task dft frequencies`。出力で次を確認した。
  - `DFT-D3BJ Model`、s6=1.0、s8=1.2177、a1=0.4145、a2=4.8593
  - `Grid used for XC integration: fine`、`Convergence on energy requested: 1.00D-07`
  - AO の数は HCN 48、TMA·(HF)₂ 214。
- NWChem saddle: driver `maxiter 50 / trust 0.1 / sadstp 0.1 / inhess 2`、`task dft saddle`。
- NWChem ZTS: `string: nbeads 9, maxiter 20, stepsize 0.05, interpol 3, tol 1e-5, freeze1/freezeN, impose, xyz_path` と `task dft string`。W7 では一度も実行されていない。
- ωB97X-D3 の SP: `xc wb97x-d3` で disp 行はない。D3 は汎関数に含まれ、出力に `DFT-D3 Model` が出る。
- CCSD(T) の SP: `ccsd ; freeze atomic ; end` と `task ccsd(t) energy`。閉殻では scf ブロックを書かない（SCF 閾値は NWChem 既定の 1e-6）。
- xTB: `xtb input.xyz --opt vtight --gfn 2 --chrg 0 --uhf 0`、`xtb input.xyz --hess --gfn 2 …`。`--etemp` と `--alpb` は付かない。
- pysisyphus（RUN.yaml）の設定:
  - `calc {type: xtb, gfn: 2, pal: 1}`
  - `cos {type: gs, max_nodes: 9, climb: true, climb_rms: 0.005}`
  - `tsopt {type: rsirfo, thresh: gau, max_cycles: 150}`
- GoodVibes の呼び出し: `compute_thermo(qcdata=QCData(自前の射影振動数ほか), QS='grimme', s_freq_cutoff=100, temperature=298.15, freq_scale_factor=0.989, zpe_scale_factor=0.975, invert=None, symm=True)` に、感度用の 5 設定を加える。1 回 0.27 s。

**主な既定値とゲート**

| 項目 | 値 | 場所 |
|---|---|---|
| 停留点レベル | pbe0 / def2-svpd / d3bj / grid fine / scf 1e-7 | `configs/methods/pbe0-d3bj_def2-svpd.yaml` |
| 低レベル | gfn 2、電子温度 None（300 K 扱い）、溶媒なし | `configs/methods/gfn2.yaml`、`core/method.py:14,81-82` |
| 手法パネル | [pbe0-d3bj_def2-svp, pbe0-d3bj_def2-tzvp, wb97x-d3_def2-tzvpd, ccsd-t_def2-tzvp]、targets: reaction_stationary_points | `configs/pipelines/method_panel.yaml:7-8` |
| composite のエネルギー層 | `energy_method` の既定は None。どの同梱 pipeline も設定していない。一方で `pbe0-d3bj_def2-tzvp.yaml` のコメントは「Energy layer of the composite G」になっている | `stages/thermochemistry.py:44,66-76` |
| 要求と観測の Level 照合 | DFT は汎関数・基底・分散・grid・scf_tol・溶媒・版数。WFT は method・基底・版数だけ | `core/method.py:61-99` |
| same_pes | program・版数・method・基底・分散・溶媒・電荷・多重度・電子温度。numerics=True なら grid と scf_tol も | `chemistry/gates.py:47-51,94-100` |
| LOT 不一致 | 参加種の freq 層か energy 層が一致しなければ `mixed_level_of_theory` | `thermochemistry.py:145-181` |
| rankable | outcome が elementary / degenerate / reassigned、ΔG‡ がある、blocker 4 種がない、\|ΔZPE‡\| ≤ 1 + 4·max(1, n_H) | `gates.py:303-321` |
| パネルの符号ゲート | どこかの level で ΔE_rxn > 0、別の level で < 0 なら不一致（許容幅なし） | `reporting/summary.py:210` |
| エネルギー窓 | explore は ΔE‡ 150 / ΔE_rxn 100 kJ/mol、仮説と行 3 は 40 kcal/mol、障壁の判定は 1.0 kcal/mol | `stages/explore.py:37-40`、`gates.py:37-38` |
| 停留点 method の一意性 | 検査は同梱 pipeline に対するテストだけ（`tests/integration/test_configs.py:28-36`） | — |

**利用者が設定する項目**: site（scratch、cores、engines の版数・実行ファイル・ranks・メモリ・timeout）、system（species、compositions、reactions）、pipeline の選択、method ファイル、conditions、gates、thermo、sp、minima、reaction-paths、report。手法を変えるには、`minima(dft).method` と `paths.method` の 2 か所を書き換える。composite G を使うには、利用者が sp → thermo(energy_method) → report の pipeline を自分で書く必要がある。

**実測（W7、すべて PBE0-D3BJ/def2-SVPD）**

| 系 | 所要時間 | 結果（kcal/mol） |
|---|---|---|
| HCN→HNC | 1 分 40 秒 | TS −1128.5i、ΔE‡ 46.62、ΔG‡ 42.22、ΔE_rxn 13.18、ΔG_rxn 12.68、ΔZPE‡ −3.33 |
| HONO trans→cis | 7 分 22 秒 | ΔE‡ 13.67、ΔG‡ 12.24、ΔE_rxn +0.13、ΔG_rxn +0.08 |
| NH3 反転 | 2 分 08 秒 | ΔE‡ 4.26、ΔG‡ 3.85 |
| TMA·(HF)₂ discover | 52 分 33 秒 | freq は 17 原子で 897 s と 903 s。ΔE_assoc −31.90（CP なし）、ΔG_assoc（1 atm）−11.59 |

**手法パネル（本レビューで初めて実行、W7 の JobStore を複製して追記）**: ΔE‡ / ΔE_rxn

| 手法 | HCN | HONO | NH3（ΔE‡） |
|---|---|---|---|
| CCSD(T)/def2-TZVP | 48.08 / 15.43 | 11.96 / +0.22 | 5.50 |
| PBE0/def2-SVPD（基準） | 46.62 / 13.18 | 13.67 / +0.13 | 4.26 |
| PBE0/def2-TZVP | 46.64 / 13.97 | 13.40 / −0.03 | 4.71 |
| ωB97X-D3/def2-TZVPD | 46.37 / 13.15 | 12.82 / +0.33 | 4.51 |
| PBE0/def2-SVP | 47.33 / 13.40 | 14.62 / −0.64 | 5.74 |

HONO では sign_disagreement=True となり、trans_to_cis が rankable=False（blocker は method_sign_disagreement）になった。W7 では rank 1 だった。

**結合エネルギーのプローブ（手書きの NWChem 入力）**

- F⁻ + HF → FHF⁻ の De。実験値は 45.8±1.6（Wenthold & Squires）。
  - PBE0: def2-SVP 66.68、def2-SVPD 47.01、def2-TZVP 53.82、def2-TZVPD 46.39
  - ωB97X-D3/def2-TZVPD 44.75
  - CCSD(T): def2-TZVP 51.00、def2-TZVPD 43.87、aug-cc-pVTZ 44.94
- NH3·HF の相互作用エネルギー（非補正 / CP 補正）
  - PBE0/def2-SVPD −16.79 / −15.74（BSSE 1.05）
  - ωB97X-D3/TZVPD −14.58 / −14.33
  - CCSD(T)/aug-cc-pVTZ −13.81 / −13.13
- 注意: NWChem 7.2.3 の `bsse` ブロックが出力する「BSSE error」は、第 1 モノマーの項しか含んでいなかった。上の値は、5 つのエネルギーから自分で計算したものである。

#### 2.0.2 化学的妥当性

- **骨格は妥当です。** 安価な手法で構造を探索し、ハイブリッド汎関数で停留点を求め、振動数から熱化学を出す多層構成で、Bursch/Grimme 2022 の推奨に沿っています。
  - 実出力で、PBE0 の D3BJ パラメータが Grimme 2011 の値と一致すること、解析 Hessian に vdW の寄与が加算されていることを確認しました。
  - def2-SVPD の diffuse 関数は陰イオンに不可欠です。FHF⁻ の De は、SVP では +21.7 ずれますが、SVPD では +1.2 に収まります。
- **陰イオン・H 結合の系では弱点があります。**
  - (a) パネルの「参照」である CCSD(T)/def2-TZVP と、composite 層とコメントされた PBE0/def2-TZVP には、どちらも diffuse 関数がありません。FHF⁻ ではどちらも停留点レベルより悪い値でした（51.0 と 53.8。停留点レベルは 47.0）。CCSD(T)/def2-TZVPD にすると 43.87 になり、aVTZ との差は −1.07 に縮みます。
  - (b) 会合量は非補正の SVPD の値です。NH3·HF では CCSD(T)/aVTZ-CP より 3.7 深く、内訳は BSSE 1.05 と汎関数の誤差約 2.6 です。TMA·(HF)₂ の ΔG_assoc も過大な結合である可能性がありますが、これは推定で、直接は確かめていません。
  - (c) ΔE_rxn ≈ 0 の反応では、符号の不一致を判定する物理的な意味がありません（§2.0.6 の U0-I1）。
  - (d) amine·(HF)n の「同一 basin」は PBE0 の PES 上での結論です。同梱のエンジンでは、高いレベルで確かめる手段がありません。NWChem の CCSD は RHF 専用で、17 原子では実行できません。ωB97X-D3 の freq は有限差分になります。これは設計の範囲どおりの限界です。
- 開殻の扱い: DFT は `odft` と spin ゲート、CCSD(T) は INPUT_INVALID で処理されます。気相 1 atm で計算し、1 bar / 1 M へ純関数で換算する処理も正しいです。

#### 2.0.3 有用性

- **欠かせず、費用に見合うもの**
  - 単一 PES、別ジョブの freq、LOT の一致ゲート、GoodVibes の qRRHO。
  - xTB の初期 Hessian（opt の時間 −35%）と、低レベル経路を DFT SP で判定すること。どちらも最も費用対効果が高い部分です。
  - 手法パネル。小分子なら数分で手法への依存性を定量できます（HCN は 63.7 s）。
- **パネルの構成に無駄があります。**
  - def2-SVP の行は、diffuse を持たない DZ 基底の既知の欠陥を再確認するだけです。そのうえ、誤った blocker の原因になっています。
  - TZVP（PBE0 と CCSD(T)）は、陰イオンや H 結合ではパネルを誤った方向に引きます。
- **composite G は、配線されていないので現状の価値は 0 です。**
- **qRRHO の感度幅**は、剛直な小分子では約 1e-4 kcal/mol しかなく、同順位の判定には事実上効きません。それでも安価で、軟モードの多い錯体では意味があるので残します。

#### 2.0.4 複雑さ

- 過剰な点、または減らせる点:
  - (1) method ファイル 6 本のうち、`pbe0-d3bj_def2-svp` は削除でき、TZVP 系 2 本は TZVPD 版に置き換えられる。
  - (2) `energy_method` は「実装はあるが使われない」状態で、これが最も分かりにくい。配線して fail-closed にすれば、この状態は解消する。
  - (3) 停留点の method id を 2 か所に書き、一致をテストでしか検査していない。
- 残してよいもの:
  - Policy の 13 閾値。既定のまま動く。
  - explore の窓の単位（kJ/mol）。フィールド名に単位が入っているので問題ない。
  - xfine を FORBIDDEN にしていること。対象は同梱の設定だけなので、実害のある矛盾ではない。

#### 2.0.5 ライブラリ仕様の要点

- **NWChem 7.2.3**
  - `disp vdw 4` は D3BJ です。文書の対応表には PBE0 がありませんが、実装にはパラメータがあり、実出力と Grimme 2011 の値が一致します（https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html）。
  - grid fine の目標精度は 1e-7 です。DFT の `convergence energy` の既定は 1e-6 なので、1e-7 を明示しているのは正しいです。
  - `xc wb97x-d3` は D3 を自動で含みます。一方、その `task dft frequencies` は有限差分の Hessian になり（delta 0.01）、HCN で 26 s かかります（PBE0 の解析 Hessian は 4.6 s）。
  - CCSD モジュールは閉殻 RHF 専用で、MAXITER の既定は 20、THRESH は 1e-6 です（https://nwchemgit.github.io/CCSD.html）。
  - `bsse` ブロックの「BSSE error」は、第 1 モノマーの項しか含みません（本レビューのプローブで再現）。
- **GoodVibes 4.3.0**
  - `concentration=None` なら 1 atm 基準です。`symm` は pymsym を使います。
  - 同梱のスケール因子 DB に、PBE0 は PBE0/MG3S しかありません。
  - `spc` 機能もありますが、hfauto は同じ式を自前で実装しています。
- **文献**
  - Bursch, Mewes, Hansen, Grimme, Angew. Chem. 2022（https://pmc.ncbi.nlm.nih.gov/articles/PMC9826355/）: 反応エネルギーには TZ 以上の基底が要り、陰イオンには diffuse 関数が要り、多層プロトコルを推奨しています。
  - HCN/HNC の CCSD(T)/ANO: 障壁 48.3、吸熱 14.7 kcal/mol（https://cdnsciencepub.com/doi/pdf/10.1139/v96-120）。

#### 2.0.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U0-I1 | high | `method_sign_disagreement` に許容幅がありません（`summary.py:210` の `any(v>0) and any(v<0)`）。HONO では、正しい符号を示す基準レベル（+0.13）、ωB97X-D3（+0.33）、CCSD(T)（+0.22、TZVPD では +0.47）に対し、diffuse のない PBE0/SVP（−0.64）と PBE0/TZVP（−0.03）だけが負になり、検証済みの素反応が順位から外れました。縮退反応は ΔE_rxn が厳密に 0 なので誤判定は起きず、影響するのは \|ΔE_rxn\| ≈ 0 の非縮退反応に限られます。 |
| U0-I2 | medium（検証で high から下げた） | パネルの参照 CCSD(T)/def2-TZVP と PBE0/def2-TZVP に diffuse がなく、陰イオン・H 結合では停留点レベルより悪い値になります。ただし「エネルギー層」はどの pipeline にも配線されていないので、composite に誤差が入る経路は現状ありません。実害は、パネルの「参照」が誤った方向を示すことです。 |
| U0-I3 | medium | composite G が配線されておらず、TZ の SP は ΔG・会合量・順位に入りません。ただし既定の順位指標は停留点レベルの ΔG‡ なので、会合量の誤差は順位に直接は入りません。TMA·(HF)₂ の値が深すぎるという点は NH3·HF からの推定です。障壁の改善も小さく（MAD 1.47 → 1.19）、欠陥というより「使う手段がない」問題です。 |
| U0-I4 | medium | `energy_method` が fail-open です（`thermochemistry.py:75` の `layer.get(fp, (None, freq))`）。全 subject で SP が欠けると、注記も blocker も付かず、SVPD の ΔG が composite と区別できないまま出ます。同梱の設定からは到達しない潜在的な欠陥ですが、composite を配線するならその前に直す必要があります。 |
| U0-I5 | medium | 停留点の method を 2 か所に書き、一致をテストでしか検査していません。PipelineConfig の検証は stage id の重複しか見ず、決定表の行 1 も両端の level_key どうししか比べません。README は手法を変えるときに 2 か所を書き換えるよう案内しているので、実際に起こりえます。 |
| U0-I6 | low（検証で medium から下げた） | amine·(HF)n のプロトン位置と「同一 basin」の結論は、PBE0 の PES 上のものです。TMA·(HF)₂ の極小は N···H 1.274 Å、H–F 1.132 Å の共有プロトン構造です。設計の範囲どおりの限界なので、文書に明記すれば足ります。 |
| U0-I7 | low | 同順位の判定に使う幅は qRRHO の設定差だけで作られ（HCN 42.2177〜42.2183）、手法の誤差（パネルでの ΔE‡ の幅は 1.5〜2.7）は入りません。解釈上の注意であり、コードの修正は不要です。 |
| U0-I8 | low | CCSD の maxiter は NWChem 既定の 20 で、HONO の TS は 18 反復かかりました。maxiter 5 で試すと rc=255 で、エネルギーは出力されませんでした（fail-closed）。 |

#### 2.0.7 改良案

**U0-P1 符号ゲートに ±1 kcal/mol の不感帯を入れ、パネルから def2-SVP を外す** — must、【複雑さ減】
- 変更:
  1. `reporting/summary.py` に定数 `_SIGN_DEADBAND_KCAL = 1.0` を置き、判定を `any(v > 1.0) and any(v < -1.0)` に変える。Policy の knob にはしない。
  2. `method_panel.yaml` から `pbe0-d3bj_def2-svp` を除き、その method ファイルを削除する。
  3. HONO 型（−0.64〜+0.47）で発火しないテストを 1 件加える。
  4. 検証で追加された修正箇所: `tests/integration/test_pipelines_e2e.py:96` の `assert panel >= {"def2-svp", …}` と、`docs/design.md:162`・`refactor_design.md:1277` の method 一覧。
- 効果: HONO trans→cis が rank 1 に戻る。±1 を超える本物の符号反転は、これまでどおり blocker になる。パネルの SP が停留点 1 点あたり 1 ジョブ減る。
- コスト: コード 2 行、テスト 1 件、ファイル 1 本の削除。U8-IMP-1 と同じ内容なので、§3 の X4 で統合した。

**U0-P2 パネルの TZ 層を def2-TZVPD にそろえる** — should、中立
- 変更: `pbe0-d3bj_def2-tzvp.yaml` を `pbe0-d3bj_def2-tzvpd.yaml` に、`ccsd-t_def2-tzvp.yaml` を `ccsd-t_def2-tzvpd.yaml` に置き換える。あわせて、`test_pipelines_e2e.py:96` の期待値、`design.md:162`、`refactor_design.md` の「エネルギー層（def2-TZVP）」の記述も直す。
- 効果: FHF⁻ の De の CCSD(T)/aVTZ との差が、参照 CCSD(T) で +6.06 → −1.07、PBE0 の TZ 層で +8.88 → +1.45 に縮む。中性の小分子ではほとんど変わらない（HONO の CCSD(T) ΔE_rxn は +0.22 → +0.47）。
- コスト: CCSD(T) の関数の数が HONO で 99 → 126 になる（1 点 103〜135 s）。ファイルの本数は変わらない。
- U8-IMP-2（CCSD(T) を既定から外す）と組み合わせた最終形は、§3 の X4 を参照。

**U0-P3 composite G を追記用 pipeline として別に用意する（修正版）** — could、中立
- 変更: method_panel は書き換えない。追記用の pipeline を 1 本だけ新しく作る（例: `configs/pipelines/composite.yaml`）。段は次のとおり。
  - `{comp_sp: sp, methods: [wb97x-d3_def2-tzvpd], targets: all_minima}`
  - → `{comp_thermo: thermo, energy_method: wb97x-d3_def2-tzvpd}`
  - → `{comp_report: report}`
- あわせて、`gates.rankable` に「`dE_act_kcal ≤ 0` なら blocker `barrier_vanishes_at_energy_layer`」を 2 行で加える。
- 動作: artifact id が同じなので `Manifest.union` の後勝ちで機能し、パネルで計算済みの SP は JobStore で再利用される。前提として U0-P4 が必要。
- 効果（CCSD(T) 比の誤差）: NH3·HF の相互作用エネルギー 3.7 → 1.45、FHF⁻ +2.1 → −0.2、HONO の ΔE‡ +1.7 → +0.86。主な利点は ΔG_assoc で、障壁の改善は MAD で −0.3 と小さい。

**U0-P4 energy_method を fail-closed にする** — should、中立
- 変更: `_species` に requested フラグを渡し、`energy_calc` が None なら `Gate(False, ('energy_layer_missing',))` にする。これで既存の `thermo_unavailable` の経路に乗る。3〜5 行とテスト 1 件。
- 補足: P3 を採らないなら、代わりに energy_method と composite を削除する（約 −15 行）。

**U0-P5 停留点 method の一致を PipelineConfig の検証へ移す** — should、中立
- 変更: `_unique_ids` と同じ `model_validator` に、「`minima(level=dft)` と `reaction-paths` の `settings()['method']` が 1 種類であること」の検査を加える（約 6 行）。テストの assert は読み込み時の拒否に置き換える。method_panel のように停留点の段を持たない pipeline には影響しない。
- 効果: 食い違いを、数時間の計算の前に拒否できる。

**U0-P6 CCSD の maxiter を 50 にする** — could、中立
- 変更: `render_wft` で、ccsd ブロックに `maxiter 50` を加える（1 行）。U8-IMP-3 と同じ内容。

**U0-P7 停留点レベルの根拠と結果の読み方を validation.md に 3〜5 行で書く（修正版）** — could、中立
- 書く内容:
  1. def2-SVPD を選んだ根拠。FHF⁻ の De は、SVP で +21.7、SVPD で実験値の誤差内。障壁は CCSD(T) に対して ±1.7 以内。
  2. 順位の差が約 1.5 kcal/mol 未満なら、手法誤差の範囲であること。
  3. ΔG_assoc は非補正の SVPD の値で、過大に結合している可能性があること。
  4. amine·(HF)n の結論は PBE0 の PES 上のものであること（1 文）。
- ωB97X-D3 の freq が数値微分になるという記述は省く。

#### 2.0.8 見送った案

- 停留点の opt と freq を def2-TZVP(D) でやり直す — 17 原子の SVPD でもすでに約 900 s/本かかり、コストが数倍になる。エネルギーは SP 層で補える。
- 停留点の汎関数を ωB97X-D3 に変える — freq が有限差分になり 5.7 倍遅く、ノイズも大きい。3 系の MAD の改善は 1.47 → 1.19 と小さく、W7 の検証をすべてやり直すことになる。
- 会合量の CP 補正を自動化する — ghost 原子の SP が必要で、TS では断片の定義があいまいになる。NWChem の bsse 出力にも不具合がある。TZVPD 層にすれば BSSE は約 0.2 に下がる。
- gCP / DFT-C、r2SCAN-3c、ωB97X-3c — NWChem は D4 と gCP に対応していない。ORCA は削除済み。
- ORCA / DLPNO-CCSD(T) を再び導入する — W6 で意図して削除したエンジンと、その保守が戻ってくる。当面は文書に限界を明記する（P7）。
- パネルの幅や固定の手法誤差を、同順位の判定幅に組み込む — 固定値は任意で、パネルとの結合も増える。
- 符号ゲートの許容幅を Policy の knob にする — 定数 1 kcal/mol で十分。
- PBE0/def2-SVPD 専用のスケール因子をフィットする — 影響は 0.1 kcal/mol 未満。
- xfine を既定にする、または FORBIDDEN から外す — 3.6 倍のコストで結論は変わらない（CH-12）。
- 溶媒を既定で有効にする — 対象は気相で、ReaDuct も溶媒に対応しない。
- Wigner / Eckart のトンネル補正を加える — 順位は障壁の比較が目的で、速度定数ではない。
- MP2 をパネルに加える — CCSD(T)/TZVPD が小分子で安価に実行できる。
- （U0-P3 の元案）method_panel を all_minima と thermo に書き換える — 全手法を単量体まで計算することになり、診断用のパネルと ranking の意味が混ざる。

#### 2.0.9 変更不要な点

- GFN2-xTB で探索・経路を求め、判定は DFT SP で行う 2 層構成。
- PBE0-D3BJ/def2-SVPD、grid fine、SCF 1e-7。D3BJ のパラメータと、Hessian の vdW 寄与は確認済み。
- 全 DFT デッキに convergence energy を明示すること。停留点系は numerics=True、string は numerics=False という二層の数値方針。
- opt と freq を別ジョブにすること、xTB の初期 Hessian、trust 0.1。
- Level の照合と `mixed_level_of_theory` の blocker。
- GoodVibes に自前の振動数とスケール因子を明示して渡すこと、thermo_consistent ゲート。
- 1 atm で計算し、1 bar / 1 M へ純関数で換算すること。
- CCSD(T) の開殻を INPUT_INVALID で拒否すること。
- method_panel を既存 run への追記型の任意実行にし、Level.full_key をキーにすること（CH-25）。
- qRRHO の感度幅と、気相・溶媒なしを既定にしていること。

### 2.1 U1 入力構造・電子状態（structures）

#### 2.1.1 現状仕様

**流れ**（`stages/structures.py`、`chemistry/smiles.py`、`chemistry/electronic_state.py`、`core/system.py`）

1. `load_system` が system YAML を読む（extra=forbid）。`_check_references` が検査するのは、species id の重複、composition が未知の species を参照していないか、反応の両端が role=endpoint か、の 3 点だけ。
2. 再開の判定（`config_sha`）: system の内容（xyz はファイル内容の sha256）で stage を飛ばすかを決める。RDKit の版数は入らない。
3. 各化学種について `_write` を行う。
   - xyz だけのとき: `<run>/structures/xyz/<id>.xyz` にバイト単位で複製する。
   - SMILES だけのとき: `smiles_to_molecule` で 3D 構造を作る。
   - 両方ある、またはどちらもないとき: ValueError。
4. SMILES の処理:
   - `.` を含めば拒否する（CH-23。錯体は compositions で表す）。
   - `MolFromSmiles` → 形式電荷と宣言電荷が一致するか確認 → `AddHs` → ETKDGv3（seed 20260925）で 1 回だけ埋め込む。
   - 力場による最適化はしない。多重度と RDKit のラジカル電子数は照合しない。
5. 書いたファイルを読み直し、Geometry の指紋（1e-6 Å に丸めた座標の sha256）を取る。
6. `check_electronic_state` で、電子数 = ΣZ − q の偶奇と、電子数 ≥ m−1 を確かめる。表にない元素があれば、この検査を飛ばす。
7. `SpeciesRecord` を作る。
   - `composition_key` = `<Hill 式>_q<電荷>_m<多重度>`
   - `state_label` = 断片の Hill 式と WL ハッシュ（1.45 Σr_cov 未満を結合とみなす）
8. 例外の割り付け: ImportError は EXECUTABLE_MISSING、OSError / ValueError は INPUT_INVALID の failed artifact にする。stage は止めずに次へ進む。
9. 宣言反応の両端の `geometry.symbols` が一致しなければ、両端を INPUT_INVALID にする。
10. 下流への受け渡し:
    - 電荷と多重度は、NWChem（`charge`、`mult`、m>1 なら `odft`）、xTB と CREST（`--chrg`、`--uhf m−1`）にそのまま渡す。
    - 組成の電荷と多重度は conformers stage が決める（電荷は和、多重度は高スピン）。

**エンジン呼び出し**: RDKit 2026.03.3 はプロセス内で呼ぶ。ETKDGv3 の既定値のうち seed だけを変えている（useExpTorsionAnglePrefs、useBasicKnowledge、ETversion 2、enforceChirality など）。例: NH3 の SMILES 'N' から N–H 1.036〜1.042 Å、HNH 106.8〜113.0° の構造ができる。

**主な既定値**

| 項目 | 値 |
|---|---|
| seed、構造数、力場 | 20260925、1、なし（設定では変えられない） |
| 複数断片の SMILES | 拒否する |
| 電子状態の検査 | m ≥ 1。偶奇と、電子数 ≥ m−1。表にない元素があれば検査全体を飛ばす |
| 共有結合半径 | Cordero 2008 の 37 元素。大文字小文字を区別する。表にない元素は INPUT_INVALID |
| 状態ラベル | 結合は ≤1.15 Σr_cov、非結合は ≥1.45。前の状態がなければ中間帯は結合とみなす。WL 3 反復、ハッシュの先頭 8 桁 |
| SpeciesInput の既定 | charge 0、multiplicity 1、role monomer |
| 組成の多重度 | coupled_multiplicities の最大値（高スピン）。利用者は上書きできない |

**利用者が設定する項目**: system_id、species（id、xyz か smiles、charge、multiplicity、role）、compositions（id、components）、reactions（id、reactant、product、coordinate、torsional）。stage 自体には設定キーがない。

**実測**: W7 の 8 run と final_check の 1 run で、structures は常に done、failed は 0、1 s 未満、外部ジョブは 0 件だった。状態ラベルは HCN と HNC を区別し、HONO の配座異性体は同じラベルになった。TMA·(HF)₂ の neutral（N–H 1.317 / H–F 1.097 Å）と shared_proton（1.271 / 1.135 Å）は、どちらも中間帯の N–H が結合と数えられ、同じラベルになった。

SMILES のプローブの結果:
- 埋め込みは 0.1〜30 ms で、同じ SMILES からは毎回同じ座標ができる。
- `[CH2]` m=1、`[O][O]` m=1 は合格した。`[CH3]` m=1 は偶奇の検査で不合格になった。
- `[2H]O[2H]` は H2O_q0_m1 になり、同位体は失われた。
- `[Xe]`、`[Cs+]`、SnCl4 は INPUT_INVALID になった。
- decalin と bicyclopentyl は同じラベルになった（1-WL の限界）。

#### 2.1.2 化学的妥当性

- **電子状態の検査**: 閉殻種で電荷を ±1 取り違えた場合と、奇数電子種を既定の m=1 で宣言した場合は、必ず検出します。閉殻が中心の対象化学には十分です。見逃すのは、±2 の電荷の誤りと、不対電子を偶数個持つ種の多重度の誤りです。
- **SMILES の埋め込み**: ETKDG の構造をそのまま使うことは、RDKit の公式見解（そのまま使える品質）に合っています。後段で必ず最適化されるので、初期構造としては十分です。ただし同位体は黙って ¹H になります。
- **宣言反応（最大の問題）**: 元素列が一致することは必要条件ですが、十分条件ではありません。GS・ZTS・IDPP・`align_mapped` は、入力の原子の順番をそのまま反応の原子写像として使います。
  - SMILES では、RDKit が既定で `[H]` を除き、`AddHs` が H を重原子の後ろに付け直すので、H の対応を利用者が決められません。
  - 例（プローブ）: acetaldehyde 'CC=O' → vinyl alcohol 'C=CO' が元素列の検査に通り、写像上は H が 2 個移る過程として扱われました。正しくは 1,3-H 移動 1 回です。
- **宣言座標（coordinate）**: 原子数も添字も検査されません。負の添字は、numpy の規則で黙って末尾の原子を指します。
- **状態ラベル**: 前の状態がないとき中間帯を結合とみなす規則は、CH-07 の意図どおりです。TMA·(HF)₂ の 2 入力が同じラベルになるのは正しい結果で、DFT でも同一 basin でした。1-WL の衝突と、立体を区別しないことは既知の限界ですが、basin の同一性は RMSD で判定するので、実害は小さいです。
- **組成の多重度を高スピンに固定すること**: NWChem のデッキでは m=1 が RKS になるので、単参照で一貫して扱える唯一の選択です。

#### 2.1.3 有用性

- 必須の入口で、コストは 1 秒未満です。
- 次の部分は費用対効果が高いです: xyz の複製と指紋、偶奇の検査（数行で ±1 の電荷の誤りをすべて止める）、形式電荷の照合、`.` の拒否、元素列の検査、`composition_key`（q と m を含む）。
- SMILES の埋め込みは、実際には amine パネルの単量体にしか使われていませんが、有用です。
- 欠けていて価値が大きいのは、宣言反応まわりの入力検査だけです。

#### 2.1.4 複雑さ

- 本体は 90 行で、設定キーはありません。過剰な部分はありません。固定値（seed、1 構造、1.15 / 1.45、WL 3 反復）を設定にしないのは正しい判断です。
- 削れるもの:
  - 未知元素のときに検査を飛ばす fail-open の分岐。
  - `_write` での実行時の排他検査（validator に移せる）。
- 足りないのは複雑さではなく、入力の検査です（端点の形式、座標の添字、個数 ≥ 1、id の一意性）。

#### 2.1.5 ライブラリ仕様の要点

- **RDKit 2026.03.3**（https://www.rdkit.org/docs/GettingStartedInPython.html 、https://www.rdkit.org/docs/RDKit_Book.html）
  - ETKDG の配座は「後の最小化なしで、そのまま使える」。seed を固定すれば再現する。
  - `MolFromSmiles` は既定で `[H]` を除く（removeHs=True）。`AddHs` は H を重原子の後ろに付け足す。atom map 番号を付けても、原子は並び替わらない。
  - 同位体は `GetIsotope()` に入り、`GetSymbol()` は 'H' を返す。
  - `enforceChirality` が保つのは、指定した中心だけ。
- **ETKDG の文献**: Riniker & Landrum, JCIM 2015（https://doi.org/10.1021/acs.jcim.5b00654）、Wang ら JCIM 2020（https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025）。どちらも最低エネルギー配座は保証しないので、最低配座を CREST に任せる現在の分担は妥当。
- **NWChem**: MULT の既定は 1。ODFT が必要なのは一重項を開殻で解く場合だけ（https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html）。
- **xTB / CREST**: `--chrg` と `--uhf` は `.CHRG` / `.UHF` ファイルより優先する（https://xtb-docs.readthedocs.io/en/latest/commandline.html）。
- **Cordero 2008**（https://doi.org/10.1039/b801115j）: 値は原著と一致する。
- **1-WL**（http://proceedings.mlr.press/v139/balcilar21a/balcilar21a-supp.pdf）: decalin と bicyclopentyl は、反復を増やしても区別できない。

#### 2.1.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U1-01 | high | SMILES で与えた宣言反応の端点は、H の原子写像を利用者が制御できず、元素列の検査に通ったまま別の機構として計算されます。同梱の system で SMILES を端点にしているものはありませんが、README と design は SMILES を化学種一般の入力として案内しています。 |
| U1-02 | medium | 宣言座標の原子数と添字を検査していません。負の添字は黙って別の原子を指します。範囲外の添字は、両端の DFT opt と freq が終わった後の reaction-paths で IndexError になり、stage ごと落ちます（actions.py が捕まえるのは ValueError だけ）。 |
| U1-03 | low（検証で medium から下げた） | SMILES の同位体とラジカル電子を黙って捨てます。ただし、`[CH2]` を m=1 で与えることは一重項カルベンとして正当な意図でもあります。O=O の m=1 は RDKit でも検出できず、対象化学の外の経路です。 |
| U1-04 | medium（検証で low から上げた） | system の静的な検査が薄い点です。xyz と SMILES の排他は実行時にしか見ず、個数 ≤ 0 を許し、組成・反応の id が重複しても通ります。さらに composition id と species id が重なると（例: どちらも 'tma'）、conformers で species_id（`<base>_c00`）、失敗 artifact の id、`best` の dict のキーが実際に衝突し、後の artifact が黙って勝ちます（`conformer_search.py:86-87,143,233`）。 |
| U1-05 | low | `check_electronic_state` の未知元素スキップは fail-open です。検証による訂正: この分岐には実際に到達しますが、直後の `state_label` が INPUT_INVALID にするので、結果は fail-closed です。エラーの理由が分かりにくくなるだけです。 |
| U1-06 | low | state_label の既知の限界（1-WL の衝突、E/Z と鏡像体を区別しない）。ラベルはグループ化にしか使わないので、対応は不要です。 |
| U1-07 | low | 組成の多重度を高スピンに固定すること、m=1 が RKS になること、未指定の立体を ETKDG が任意に選ぶことが、design.md と README に書かれていません。 |
| U1-08 | low | 熱化学の組成キーが (Hill 式, 電荷) で、多重度を含みません（U7 で扱う。§3 の X8）。 |

#### 2.1.7 改良案

**U1-I1 宣言反応の端点は xyz で与えることを必須にする** — must、中立
- 変更: `_check_references` の反応ループに、「`species.xyz is None` なら ValueError（'endpoint must be given as xyz; SMILES cannot fix the atom mapping or the conformer'）」を 1 条件足す。design.md §7 と README に 1 文、テストを 1 件加える。
- 効果: SMILES の端点で別の機構を黙って計算する経路がなくなる。検出は読み込み時なので、無駄な DFT は 0 本。同梱の system は影響を受けない。
- コスト: 約 5 行、約 30 分。

**U1-I2 宣言反応の入力検査を structures で完結させる（簡素化版）** — should、中立
- 変更:
  - (a) `CoordinateTerm` に `model_validator` を 1 つ加える。len(atoms) が {distance: 2, angle: 3, dihedral: 4}[kind] と一致し、添字がすべて ≥0 で重複がないこと。
  - (b) `structures.py:86` の条件を `a.composition_id != b.composition_id or a.geometry.symbols != b.geometry.symbols` に広げる（composition_id が q と m を含むので、別に比較する必要はない）。あわせて `max(atoms) < len(symbols)` を検査する。reason の文字列は既存の形式を使う。
- 効果: 負の添字の取り違えがなくなる。範囲外の添字も、DFT の数十分の後ではなく 1 秒未満で止まる。

**U1-I3 SMILES の同位体とラジカル電子を拒否する（修正版）** — could、中立
- 変更: `AddHs` の前に 2 つの検査を置く。
  1. `any(a.GetIsotope())` なら ValueError。
  2. `sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) > multiplicity − 1` なら ValueError（メッセージで、一重項は xyz で与えるよう案内する）。
- `Descriptors` は import しない（重い import を避ける既存方針に合わせる）。テストを 2 件加える。

**U1-I4 system の静的検査を読み込み時に集める** — should、【複雑さ減】
- 変更:
  1. `SpeciesInput` に「xyz と smiles のちょうど一方」を要求する validator を足し、`_write` を if/else にする。
  2. `components` を `dict[str, PositiveInt]` にする。
  3. species id と composition id の和集合の一意性、reaction id の一意性を検査する。
- 効果: conformers での artifact id の黙った上書きを塞ぐ。

**U1-I5 未知元素のスキップをやめる** — could、【複雑さ減】
- 変更: `return` を `ValueError('unknown element …')` に置き換え、テストの 'Xx' 行を `pytest.raises` にする（2 行）。

**U1-I6 structures と組成の規則を design.md §7 に書く** — could、中立
- 書く内容: SMILES は 1 配座で、未指定の立体は任意に選ばれ、同位体は不可。宣言反応の端点は xyz。組成の多重度は高スピンで、m=1 は閉殻 RKS になり、開殻一重項は扱わない。

#### 2.1.8 見送った案

- ETKDG の後に MMFF / UFF で最適化する — 後段で必ず最適化される。MMFF にはパラメータのない種という新しい失敗の仕方が加わる。
- structures で複数の配座を作る、または立体異性体を列挙する — conformers の責務と重なる。立体中心はほとんどない。
- 埋め込みに失敗したら再試行する — 対象の規模では失敗が起きない。
- SMILES の端点の原子写像を自動で推定する（MCS、RXNMapper） — 費用が大きく、誤って推定する危険が残る。xyz に限る方が確実。
- 同位体に対応する — freq と thermo にまたがる変更で、目的から外れる。
- 組成の多重度を上書きできるようにする — m=1 は RKS になり、開殻一重項には BS の初期推定が要る。
- 多重度や電荷の既定値を SMILES から導く — 暗黙の既定値が増える。
- state_label を InChI などに置き換える、または WL の反復を増やす — 1-WL の限界は越えられず、錯体で結合次数の推定に失敗する。
- 中間帯の閾値を中点にする — CH-07（FHF⁻ を 1 断片に保つ）に反する。
- RDKit の版数を config_sha に入れる — 再開時は保存済みの xyz を使うので、正しさに影響しない。
- preflight で RDKit の有無を検査する、INPUT_INVALID が出たら run 全体を止める — 効果が小さい、または artifact 単位の設計に反する。
- xyz の元素記号の大文字小文字を正規化する、三重項基底状態を警告する — 価値が小さい、または一般には判定できない。
- （U1-I2 の元案）電荷と多重度を別々に比べ、reason を 3 種類に分ける — composition_id の比較で足りる。

#### 2.1.9 変更不要な点

- xyz をバイト単位で複製し、読み直して指紋を取ること。ETKDGv3 で 1 回だけ埋め込み、seed を固定し、力場で最適化しないこと。
- `.` を含む SMILES の拒否と、形式電荷の照合。偶奇の検査。
- failed artifact として記録して先に進む方針。宣言反応の元素列の検査と、行 1 の BLOCKED。
- composition_key と Molecule の指紋が q と m を含むこと。
- 中間帯を結合とみなす規則（CH-07）、Cordero の 37 元素、組成の電荷は和・多重度は高スピン。
- StructuresConfig に設定キーを持たないこと。

### 2.2 U2 配座探索・錯体配置（conformers）

#### 2.2.1 現状仕様

**流れ**（`stages/conformer_search.py`、`backends/crest.py`、`chemistry/placement.py`。discover だけで使う）

1. **単量体**（role=monomer）
   - 重原子 ≤3 かつ回転可能結合 0 本なら、CREST を省いて入力構造を出力する（HF、NH3、H2O）。
   - それ以外は CREST（`nci=False`、トポロジー検査は有効）で探索する。
2. **CREST が初期最適化でトポロジーの変化を検出して止まった場合**: 停止時の構造を `crest_topology` 種として残し、`--noreftopo` で 1 回だけ再実行する。
3. **組成**
   - 電荷は部品の和、多重度は高スピン側。
   - 部品を H 結合の円錐配置（`placement.seeds`）で組み合わせ、衝突のない seed を最大 6 個作る。
   - CREST に渡すのは **seed00 だけ**（`--nci`）。酸と塩基の両方を含む組成では、labile H と受容原子を `--notopo` で自動指定する。
   - 他の seed は、CREST が失敗したときの fallback にしか使わない。
4. **選抜**: 次の順に絞り込む。
   1. `state_label`（1.45 Σr_cov の単一閾値）ごとに分ける。
   2. 重複を除く（置換不変 RMSD < 0.1 Å かつ ΔE < 0.1 kcal/mol）。
   3. ラベルの最低値から 4 kcal/mol 以内のものを、6 個まで残す。
   4. エネルギーが None の候補は、ラベルの最低値も None のときだけ残す（BUG-08）。
5. 出力は SpeciesRecord（source は conformer / placement / crest_topology）。下流の screen が全件に xTB の opt と Hessian をかける。

**エンジン呼び出し（W7 の実物）**

- 単量体（TMA）: `crest input.xyz --gfn2 --quick -T 4 --ewin 6 --chrg 0 --uhf 0`。env は `OMP_NUM_THREADS=4,1`、`OMP_STACKSIZE=4G`。`--scratch` は付けず、attempt ディレクトリで直接実行する。6.34 s。
- 組成（TMA·(HF)₂）: `crest input.xyz --gfn2 --nci --quick -T 4 --ewin 6 --chrg 0 --uhf 0 --notopo 2,14,15,16,17`。10.86 s。CREST 内部の実効設定は次のとおり。
  - 楕円体の壁ポテンシャル
  - t(MTD) 2.5 ps × 6 本、2 反復
  - 多段最適化（crude → regular → tight）
  - GC なし、energy+grad の呼び出し 13,952 回
- CREST 3.0.2 の内部仕様（ソースで確認）: argv が `--nci --quick` の順なので、runver は quick（2）に上書きされる。その結果、MTD の bias は quick のもの（k0 0.002·N、α 1.2/0.6/0.3）、長さの係数は 0.5 になる。NCI の壁と「GC なし」はそのまま残る。明示した `--ewin 6` は、quick の 5 を上書きする。

**主な既定値**

| 項目 | 値 |
|---|---|
| seeds_per_composition / keep_per_state / window_kcal | 6 / 6 / 4.0 |
| settings | nci（必須だが stage が上書きする）、quick true、threads 4、ewin 6.0、topology on、notopo_atoms () |
| 重複の判定 | RMSD < 0.1 Å かつ ΔE < 0.1 kcal/mol |
| placement の定数 | 円錐 110〜120°、方位 0/120/240°、傾き ±30°、H 結合距離 0.75×Σr_vdw、H を含む対 ≥1.2 Å、重原子どうし ≥2.2 Å、50 回、seed 20260925 |
| 版数 pin | 3.0.2（stdout のバナーと照合し、不一致なら METHOD_MISMATCH） |
| timeout | 14,400 s（ExecutionSpec の既定）。stage は deadline を渡さず、ladder もない |

**利用者が設定する項目**: seeds_per_composition、keep_per_state、window_kcal、settings（quick、threads、ewin_kcal、topology、notopo_atoms、nci）、method（gfn、alpb）、site（版数、timeout、threads）、system（role、compositions）。

**実測**

- tma_hf2_discover: stage 17 s。
  - TMA の CREST は 1 配座。
  - TMA·(HF)₂ は、seed 5 個（重複除去後）から CREST 10.86 s で 2 配座（ΔE 0.607）。screen でどちらも min_neutral にまとまった。
- amine_pilot2: stage 12 s。TMA·HF は、sp3 N で方位が縮退するので seed が 1 個しかできなかった。
- W7 では、トポロジー停止、noreftopo の再実行、失敗はいずれも 0 件だった。
- 本レビューのプローブ（hfauto の argv で実行）:
  - TMA·(HF)₃: `--notopo 2,14,15,16,17,18,19` を付けても初期最適化で停止した（影響原子は N2・H14・F15 で、すべて notopo リストに含まれる）。noref の再実行は 159 s で 3 配座だった。
  - NH3·(HF)₃: 停止し、noref の再実行も trial MTD の発散で error stop した。hfauto では NONZERO_EXIT となり、energy None の seed しか残らない。
  - `--noopt` を付けた場合: TMA·(HF)₃ は 80 s で 13 配座（最安構造は noref より 1.21 kcal/mol 低いイオン対）。NH3·(HF)₃ は 51 s で 18 配座。F⁻·(HF)₂ は 23 s（noref は 46 s）。

#### 2.2.2 化学的妥当性

- **構成は妥当です。** CREST の NCI-iMTD で非共有結合錯体の配置を安価に前探索し、H 結合の円錐配置で出発構造を作り、`--notopo` で PT 構造とイオン対を CREGEN に残します。下流では screen と DFT SP で精密化します。これは Pracht/Bohle/Grimme 2020 の標準的な使い方です。GFN2 の H 結合の幾何はよく合いますが（RMSD 0.03 Å）、PT の反応エネルギーの誤差は約 3.2 kcal/mol（MUE、JCTC 2025）あるので、GFN2 の段階でエネルギーによる最終判断をしてはいけません。
- **問題点**
  1. CREST の初期トポロジー検査は `--notopo` のリストを参照しません（`setuptest.f90`）。このため、酸塩基やアニオンの錯体は毎回停止して再実行に回り、NH3·(HF)₃ は失敗、TMA·(HF)₃ は最安構造を取り逃しました。
  2. GFN2 の 4 kcal/mol 窓をラベルごとに適用しているので、下流の窓（6 kcal/mol、同じく GFN2 レベル）より先に候補が切られます。
  3. 組成の多重度を黙って高スピンに決めています（現行の閉殻系には影響しない）。

#### 2.2.3 有用性

- **必須**: 組成に対する CREST `--nci`、placement の seed00、`--notopo` の自動付与、重複除去、版数の pin、エネルギーの欠損を None として扱うこと。
- **中程度**: 単量体の CREST。対象の単量体では入力と同じ極小しか出ないが、柔らかい単量体への保険として約 5 s と安い。
- **限界的**: seed01〜05（失敗時の唯一の出力として役立つ）、組成の crest_topology 種（多くは重複）、diagnostics.json（読むコードがない）、4 kcal/mol の窓。
- コストは DFT と比べて無視できます（1 組成あたり数十秒〜数分）。

#### 2.2.4 複雑さ

- 効かない設定や罠になる設定:
  - (a) `settings.nci` は必須フィールドなのに、stage が必ず上書きする。`{ewin_kcal: 4}` だけを書くと ValidationError になる。
  - (b) `notopo_atoms` は、分子によって意味の違う番号を全単量体に共通で渡してしまう。
  - (c) `topology='off'` を使う理由がない。
  - (d) スレッド数が site と stage の 2 か所にあり、効くのは stage 側だけ。
  - (e) エネルギー窓が 3 重になっている（`--ewin 6`、4 kcal/mol、下流の 6 kcal/mol）。
- 残すべきもの: seeds_per_composition と keep_per_state（件数の上限）、`--quick`。

#### 2.2.5 ライブラリ仕様の要点

- **CREST 3.0.2 のソース**（https://github.com/crest-lab/crest/blob/v3.0.2/src/）
  - `confparse.f90`: 引数は出現順に処理される。`-quick` は runver=2・ewin 5・optlev≤1、`-nci` は runver=4・GC なし・rotamer MD なしを設定する。
  - `choose_settings.f90`: NCI の t(MTD) は 0.10(N+0.1N²)（最小 5 ps）に rfac（quick なら 0.5）を掛ける。
  - `algos/setuptest.f90`: 初期最適化（optlev −1）の後の quicktopo の比較は excludeTOPO を参照しない。excludeTOPO を使うのは `cregen.f90:806,843` だけ。
  - `crest_main.f90:226-229`: preopt が偽なら trialOPT を呼ばない（つまり `--noopt` なら初期トポロジー停止は起きない）。
- **CREST docs**（https://crest-lab.github.io/crest-docs/page/documentation/keywords.html 、https://crest-lab.github.io/crest-docs/page/examples/example_3.html）: `--notopo` は「CREGEN のトポロジー検査を切る」、`--noreftopo` は「初期のトポロジー検査だけを切る」。NCI の例は `crest struc.xyz --nci`。
- **xtb dock（aISS）**（https://xtb-docs.readthedocs.io/en/latest/xtb_docking.html）: 2 断片専用で、PT は扱わない。placement の代わりにはならない。
- **GFN2-xTB の PT の精度**（https://pmc.ncbi.nlm.nih.gov/articles/PMC12288007/）: 気相の 30 反応で MUE 13.5 kJ/mol。
- 手法論文: Pracht, Bohle, Grimme, PCCP 2020, 22, 7169（https://pubs.rsc.org/en/content/articlelanding/2020/cp/c9cp06869d）、CREST JCP 2024, 160, 114110。

#### 2.2.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U2-ISS-1 | high | CREST の初期トポロジー検査は `--notopo` を無視するため、酸塩基・アニオンの錯体は毎回停止して再実行に回ります。NH3·(HF)₃ では再実行も失敗し、TMA·(HF)₃ では最安構造を 1.21 kcal/mol 取り逃しました。失敗しても placement の seed が下流に流れるので全損にはなりませんが、CREST による PT・イオン対の標本化は失われます。amine_hf_panel の (HF)₃ 系 5 組成（nh3/tma/anl/pyr_hf3、tma_hf3_h2o）が該当しえます。 |
| U2-ISS-2 | medium | conformers の GFN2 窓（ラベルごとに 4 kcal/mol）が候補を先に切っています。検証による訂正: 下流の `window 6 / rerank_top 8` も GFN2（screen）のエネルギーに対する窓で、DFT SP を使うのはその後の per_state 3 の再順位付けだけです。正しくは「同じ GFN2 レベルの下流の 6 kcal/mol 窓を無効にしている」問題です。粗いラベルのせいで、NH3·(HF)₃ では中性の鎖構造が落ち、TMA·(HF)₃ では残る、という恣意性もあります。 |
| U2-ISS-3 | low（検証で medium から下げた） | 組成の多重度を黙って高スピンに決めています。configs に多重度 ≥2 の種はありません。 |
| U2-ISS-4 | low（検証で medium から下げた） | 効かない設定・罠になる設定（nci が必須なのに上書きされる、notopo_atoms が共通、topology=off、site の threads が効かない）。既定の configs では誰も書いていません。 |
| U2-ISS-5 | low | `--nci --quick` の順で runver が quick になることが、設計書に書かれていません。検証による訂正: quick で増えた配座は `--ewin` の縁（5.4〜6.0）にある同じラベルの近縮退構造なので、「quick の方が広く標本化する」とは言えません。最安の 2 状態はどちらでも一致します。 |
| U2-ISS-6 | low | design.md の記述がコードと合っていません（「最良の seed」は生成順の先頭にすぎない。衝突条件は H を含むすべての分子間ペアに 1.2 Å）。 |
| U2-ISS-7 | low | 組成のトポロジー停止構造は、多くが重複か縮退した付け替えです（単量体の停止構造には意味がある）。IMP-1 を入れれば、組成ではこの経路自体が消えます。 |

#### 2.2.7 改良案

**U2-IMP-1 組成（nci=True）の CREST に `--noopt` を付ける** — must、【複雑さ減】
- 変更: `backends/crest.py` の `command()` で、`settings.nci` が真なら `--noopt` を加える。`ConformerSettings` に preopt フィールドは足さない。
  - コメントに「CREST の初期トポロジー検査は --notopo を無視する（setuptest.f90）」と書く。
  - 単量体の停止 → noreftopo の経路は残す。
  - 試験:
    - argv 表のテストを更新する。
    - smoke に F⁻·(HF)₂ の組成（約 25 s）を 1 件加え、rc 0 と members ≥1 を確認する。
    - W7 の tma_hf2 の conformers だけを再実行し、同じ 2 配座が得られることを確認する。
    - SMILES 由来の部品を含む組成（amine_pilot2 の NH3·HF。既存の seed を使う）でも rc 0 を 1 件確認する。
- 効果: 成功が 2/3 から 3/3 になる。全系で最安構造が得られ（TMA·(HF)₃ では 1.21 kcal/mol 低いイオン対）、時間は約 1/3（87〜236 s → 23〜80 s）。組成では停止 → 再実行の分岐が不要になる。副作用として、CREGEN のトポロジー参照が未最適化の seed になるが、PT に関わる原子は `--notopo` で除外されるので問題ない。
- コスト: 1〜5 行とテスト更新。

**U2-IMP-2 conformers のエネルギー窓を削除する（修正版）** — should、【複雑さ減】
- 変更: `window_kcal` と数値の窓の条件を削除し、ラベルごとに `keep_per_state` 個まで残すだけにする。ただし「ラベルに数値エネルギーの候補があれば、energy None の候補は残さない」規則（BUG-08 系）は残す。`inside = [c for c in group if lowest is None or c.energy is not None]` で足りる。design.md の「4 kcal/mol 以内」も消す。
- 効果: GFN2 の段階での恣意的な切り捨てがなくなる。knob が 1 つ減る。screen の xTB ジョブは 1 組成あたり数件（各数秒）増える。

**U2-IMP-3 利用者向けの設定を整理する（最小版）** — could、【複雑さ減】
- 変更:
  1. `ConformerSettings.nci` に既定値 False を与える（1 行）。
  2. 利用者が書ける項目を quick と ewin_kcal に限り、nci・topology・notopo_atoms は stage の内部値にする。
  3. threads は site の `engines.crest.execution.threads` だけにする（CRESTEngine と thread_map の両方で site の値を使う）。
- 効果: 効かない knob が 3 つ、効かない site 設定が 1 つなくなる。

**U2-IMP-4 組成の多重度が一意でなければ止める（修正版）** — could、中立
- 変更: フィールドは足さない。`coupled_multiplicities` が 2 つ以上の値を返したら、`INPUT_INVALID('ambiguous_multiplicity:<候補>')` にする（2〜3 行）。明示指定の欄は、ラジカル対を扱う必要が出てから足す。

**U2-IMP-5 design.md に CREST の実効設定と seed の扱いを書く（事実だけ）** — could、中立
- 書く内容:
  - `--nci --quick` で runver が quick になること（MTD は 0.5 倍、2.5〜3 ps × 6 本。壁は NCI のもの。`--ewin 6` が quick の 5 を上書きする）。
  - CREST に渡すのは seed00（生成順の先頭）だけであること。
  - 衝突条件は、H を含むすべての分子間ペアに 1.2 Å であること。
  - `--noopt` を付ける理由。
  - quick を使う理由は「同じ最安状態が 1.5〜2.8 倍速く得られる」と書く。「広く標本化する」とは書かない。

#### 2.2.8 見送った案

- 組成では `--quick` を外し、素の `--nci` にする — 同じ最安状態を得るのに 1.5〜2.8 倍かかる。
- xtb で seed を事前に最適化してから CREST に渡す（停止メッセージの選択肢 A） — 依存が増え、NH3·(HF)₃ と同じ MD の発散が起こりうる。
- 6 つの seed それぞれから CREST を走らせる — コストが 6 倍で、取りこぼしの実例がない。
- placement を xtb dock（aISS）や CREST QCG に置き換える — aISS は 2 断片専用で PT を扱えない。QCG は溶媒用。
- conformers で DFT SP による再ランクを行う — minima(dft) の rerank と重なる。
- トポロジー停止構造のエネルギーを crestopt.log から読む — IMP-1 で組成ではこの経路が消える。
- 状態ラベルを細分化する — topology 全体に関わる変更になる。
- メチル回転子を回転可能結合から除く — TMA は重原子が 4 個なので省略対象にならない。
- CREST の timeout を短くし、deadline を渡す — ハングの実例がない。
- 省略した単量体の再出力をやめる — JobStore のヒットで、コストは 0。
- （U2-IMP-3 の元案）site の threads を削除する — 資源の設定は site に置く、という原則に合わせて site 側に一本化した。
- （U2-IMP-4 の元案）CompositionInput に multiplicity フィールドを足す — 推測にもとづく拡張なので、必要になってから足す。

#### 2.2.9 変更不要な点

- `--scratch` なしで attempt ディレクトリで実行すること、版数の照合、エネルギー欠損を None にすること（BUG-08）、原子の並びが変わったら METHOD_MISMATCH にすること。
- 組成への `--notopo` の自動付与（プローブで、イオン対・共有プロトン・中性の 3 状態が残ることを確認）。
- 単量体ではトポロジー検査を有効のままにし、停止したら 1 回だけ noref で再実行すること。
- placement の設計（極性 H だけを donor にする、lone-pair の円錐、固定 seed、衝突条件、置換不変 RMSD による重複除去）と、失敗時の seed の出力。
- `--quick` の併用、小さい剛体単量体の省略、JobStore のキーから threads を除くこと、組成の電荷の事前検査。

### 2.3 U3 極小構造の確定（minima screen / dft、MinimumDriver、Registry）

#### 2.3.1 現状仕様

**流れ**（`stages/minima.py`、`drivers/minimum.py`、`chemistry/gates.py`、`identity.py`、`selection.py`、`vibrations.py`）

1. **既知の極小の Registry**: 上流にある同じ tier・同じ engine・同じ Level の極小から作る。discover と known_endpoints の初回では空になる。
2. **ジョブの選択**
   - screen、または `include: all` のとき: 全 species を `species.geometry` から始める。
   - dft で `include: window` のとき:
     1. (組成, 状態ラベル) ごとに、screen エネルギーの低い順に `per_state` 3 個、窓 6 kcal/mol 以内を選ぶ。
     2. discovery の生成物とその出発極小は、必ず含める（always）。
     3. discover では `rerank_sp: true` なので、上位 8 構造を同じ DFT 法の SP で並べ直してから per_state 個を選ぶ。
   - 開始構造は screen の opt の最終構造にする。
3. **relax**: xTB は `thread_map` で 4 並列、DFT は直列（各ジョブが `mpirun -np 4`）で実行する。
   - 2 断片以上なら、開始構造での xTB の `--hess` を初期 Hessian にする。
   - opt → Registry.find。一意に割り付けば known になり、freq を省く。
   - 割り付かなければ、別ジョブの freq → `is_minimum` で判定する。判定内容は、task の一致、freq.start と opt.final の指紋の一致、`same_pes(numerics)`、振動数の本数 3N−n_ext。
4. **虚振動の 3 段方針**
   - noise（−10 ≤ ν < 0）: 極小として採用し、注記を付ける。
   - soft（−50 ≤ ν < −10）: + 方向に 0.1 Å 押して opt → freq を 1 回だけ行う。
   - saddle（ν < −50）: ±0.1 Å の mode-follow を最大 2 サイクル行う。両側が別の極小に届けば、TS 候補（DiscoveryRecord mode_follow）にする。
5. **登録**: **すべての relax が終わった後で**、species_id 順に Registry.add を呼ぶ。
   - `compare_minima` の判定: same は ≤0.02 Å かつ ≤1e-5 Eh、distinct は >0.05 Å か >5e-5 Eh、その間は ambiguous。
   - ambiguous なら、`tight` で opt → freq をやり直して 1 回だけ再判定する。
6. **出力**: calc artifact、basin ごとの MinimumRecord（tier、level_key、members、notes）、failed、mode_follow の側の species、discovery、diagnostics.json。

**エンジン呼び出し（W7 の実物）**

- xTB の opt: `xtb input.xyz --opt vtight --gfn 2 --chrg 0 --uhf 0`（0.0〜0.2 s）。表示される収束閾値は Econv 1e-7、Gconv 2e-4、max optcycles 340（17 原子）。
- xTB の Hessian: `xtb input.xyz --hess …`（17 原子で 0.3 s 以下）。hessian ファイルを `hessian.npy` にして、自前で射影する。
- NWChem の opt: §2.0.1 のデッキ（maxiter 100、trust 0.1、閾値は既定の gmax 4.5e-4 / grms 3e-4 / xmax 1.8e-3 / xrms 1.2e-3）。2 断片以上なら `inhess 2` と `job.hess` を使う。
- NWChem の freq: `task dft frequencies`。解析 Hessian で、SCF は atomic guess から始まり、CPHF tol 1e-4。
- rejudge: driver に `tight` を足す（gmax 1.5e-5 など）。W7 では一度も実行されていない。
- 継続（LADDER）: TIMEOUT と GEOMETRY_MAXITER は最新の final-NNN.xyz から 2 回まで（job.movecs と job.drv.hess を引き継ぐ）。autoz の失敗は noautoz で 1 回、SCF の不収束は damping で 1 回。

**実測**

- known_endpoints の dft stage: HCN 18 s、HONO 54 s、NH3 18 s、水 14 s。screen は 1 s 未満。
- TMA·(HF)₂ の dft stage は **3,050 s**。内訳は freq の合計 2,070 s（68%）、opt の合計 907 s（30%）。
  - AFIR の生成物は DFT で neutral と同じ basin（RMSD 0.0000 Å、ΔE 1.4e-9 Eh）になったが、その opt 404.7 s と freq 902.9 s（合わせて 43%）がかかった。
- 17 原子の freq の内訳: SCF 約 36 s、2 電子微分 188 s、CPHF 約 553 s（62%）。
- A/B（初期 Hessian の有無）: opt の合計は 900 s と 1,383 s。freq は約 71% を占める。
- 判定の実績（全 run）: mode-follow 0 回、rejudge 0 回、known による freq の省略 0 回、失敗 0 件。soft は xTB の neutral で 1 回だけ（押すと 0.61 kcal/mol 低い真の極小に移った）。
- 振動数の射影: xTB の表示とは 0.29 cm⁻¹ 以内で一致した。NWChem の表示とは最大 3.87 cm⁻¹ ずれた（TMA·(HF)₂ の最低モードで 51.84 と 55.71）。

#### 2.3.2 化学的妥当性

- **妥当です。**
  - opt と別ジョブの解析 freq を同じ数値の層で計算して判定することは、停留点の確認の正しい方法です。
  - 振動数の質量加重 Eckart 射影は標準的な手順（Miller–Handy–Adams 1980）で、NWChem の表示より正しいです。NWChem は非加重の Cartesian で射影してから質量加重します。G への影響は +0.016 kcal/mol です。
  - 虚振動の 3 段方針、soft を片側だけ押すこと（1 次の鞍点ならどちらに押しても下がる）、mode-follow も妥当で、偽の極小を防ぐ効果が実例で確認できます。
  - 同一性の閾値（0.05 Å / 0.03 kcal/mol）は CREGEN の既定より厳しく、安全側です。
  - 選択（6 kcal/mol の窓、DFT SP による再順位付け）は CENSO 型です。
- **化学的な誤りは minima の外との境界にあります。**
  - (a) explore が screen 極小の最適化構造ではなく、入力構造から出発します（U4 で扱う）。
  - (b) thermo が、採用した極小の noise / soft の虚振動を捨てます（U7 で扱う）。
- **開殻一重項（ジラジカル）** は RKS のまま検出されません。対象の極小は閉殻なので、文書に書けば足ります。

#### 2.3.3 有用性

- **必須で費用に見合うもの**: opt → freq → is_minimum、Registry、screen の全件緩和（1 秒未満）、窓と DFT SP による再順位付け（15 s）、2 断片以上の初期 Hessian（0.3 s の追加で opt −35%）、soft の押しと mode-follow（発火したときだけ費用がかかる安全網）。
- **最も効く改善点は「不要な freq を走らせないこと」です。** 現状は、同じ stage の中で freq を省く仕組みが働きません。
- **費用に見合わないもの**: ambiguous → tight の再判定。実 run で 0 回、xTB では何もせず、NWChem では約 20 分かかり、判定帯（0.006〜0.03 kcal/mol）も熱化学的に意味がありません。

#### 2.3.4 複雑さ

- 過剰な点:
  - (1) 同一性の基準が 2 種類ある（find は assign、add は compare_minima の 3 値）。それに ambiguous 専用の rejudge と、`QMEngine.optimize` の `tight` 引数がぶら下がっている。
  - (2) 並列の xTB と直列の DFT が、同じ「全部 relax してから species_id 順に登録」という処理を共有している。このため freq を省けない。
  - (3) 暗黙の分岐がある: level=screen では select が黙って無視される。include=window で screen 極小がないと、stage が黙って空になる。
- trust 0.1、maxiter 100、変位 0.1 Å、同一性の閾値がコードに固定されているのは適切です。

#### 2.3.5 ライブラリ仕様の要点

- **NWChem driver**（https://nwchemgit.github.io/Geometry-Optimization.html 、`src/driver/opt_drv.F`）
  - 既定の閾値は上記のとおり。TIGHT は 1.5e-5 / 1e-5 / 6e-5 / 4e-5。
  - 最小化の trust の既定は 0.3、MAXITER の既定は 40。
  - INHESS 2 は、restart データ（drv.hess）があればそれを優先する。
  - 負の曲率には下り方向のステップを取る。
- **NWChem vib**（`src/vib/vib_eckart.F`、`vib_vib.F`、https://nwchemgit.github.io/Hessians-and-Vibrational-Frequencies.html）: 非加重の Cartesian で射影してから質量加重する。
- **xtb**（https://xtb-docs.readthedocs.io/en/latest/optimization.html 、https://xtb-docs.readthedocs.io/en/latest/hessian.html）: vtight は Econv 1e-7・Gconv 2e-4。未収束なら NOT_CONVERGED ファイルを書く。imagmin は −5 cm⁻¹。
- **CREGEN**（https://crest-lab.github.io/crest-docs/page/documentation/keywords.html）: rthr 0.125 Å、ethr 0.05 kcal/mol、ewin 6。
- **GoodVibes 4.3.0**: invert=None なら虚振動を捨てる（U7 を参照）。
- **ReaDuct**: `matches_source` は ReaDuct 自身が再最適化した source と比べる。このため、explore の出発構造と minima の代表構造がずれると、偽の生成物ができる（U4 を参照）。

#### 2.3.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| ISS-1 | high | 同じ stage の中で、同じ basin に落ちた 2 本目以降の DFT freq を省けません（`minima.py:236-241`）。TMA·(HF)₂ では、AFIR の生成物の freq 902.9 s（30%）が無駄になりました（opt 404.7 s は find の前に走るので省けない）。下流は member の freq を使わないので、省いても結果は変わりません。 |
| ISS-2 | medium（検証で high から下げた） | explore が入力構造から出発し、偽の生成物を作ります（U4-I2 と同じ事象）。soft 押しが働いた source に限られ、化学的な誤りは DFT の合流で吸収されるので、害はコストと偽の discovery です。 |
| ISS-3 | medium | thermo が noise / soft の虚振動を捨てます。1 モードあたり G が +0.5〜1.2 ずれます。検証による補足: `invert_soft_cm1` の既定が None なので、soft の注記がある極小も既定では反転されません（U7-I1 と統合。§3 の X3）。 |
| ISS-4 | low（検証で medium から下げた） | 同一性の基準が 2 種類あり、ambiguous → tight の再判定は xTB では cache hit になるだけで、一度も発火していません。 |
| ISS-5 | low | include=window で screen 極小がないと、黙って空の結果になります。level=screen では select を黙って無視します。 |
| ISS-6 | low | 振動数の射影の約束が NWChem の表示と違います（最大 3.9 cm⁻¹）。hfauto の方が正しいので、注記すれば足ります。 |
| ISS-7 | low | Registry.find は、組成を見ずに元素列の一致だけで候補を作ります。実害はないので変更不要とし、記録だけ残します。 |
| ISS-8 | low | 開殻一重項は RKS で黙って扱われます。対象外なので文書化で足ります。 |

#### 2.3.7 改良案

**U3-IMP-1 species_id 順に relax したら、その場で登録する** — should、【複雑さ減】（IMP-2 と同時に入れること）
- 変更: `MinimaStage.run` と `relax_all` を、`for job in sorted(jobs): done = relax(job); register(done); sides も relax → register` の 1 経路にする。xTB の `thread_map` の分岐と、全件の一括ソートを削除する。
  - 注意: mode-follow の側（_mf1/_mf2）は親の直後に登録する。順番はわずかに変わるが、決定性は保たれる。
  - design.md §5/§11 と validation.md §6-1 の未解決事項を閉じ、回帰テストを 1 件加える。
- 効果: TMA·(HF)₂ の dft stage が 3,050 s から約 2,150 s に減る（−30%）。配座が DFT で合流するたびに、17 原子なら約 15 分節約できる。screen の直列化で増えるのは数秒〜数十秒。
- 検証の指摘: IMP-1 だけを入れると、stage 内の判定に assign（0.05 Å）と compare_minima（0.02 Å）が混在し、どちらが効くかが登録順で変わる。IMP-2 と同時に入れること。

**U3-IMP-2 同一性判定を assign の 1 基準に統一し、ambiguous / tight の再判定を削除する** — should、【複雑さ減】
- 変更: `Registry.add` を「同じ composition_id・level_key・元素列の basin に `identity.assign` が一意に割り付けば join、そうでなければ new」にする。
  - `_labels` では、bool の `same_minimum`（≤0.05 Å かつ ≤5e-5 Eh）を使う。
  - `_Run.rejudge`、Verdict の 'ambiguous'、`QMEngine.optimize` の `tight` 引数（protocols、nwchem engine/input、xtb）を削除する。
- 効果: 実 run での判定は変わらない。約 40〜60 行減る。統合する範囲は CREGEN より厳しいので、別の極小を誤って統合する危険は増えない。
- 注意: NWChem optimize の job key から tight が消えるので、旧 run のキャッシュは使えなくなる（旧 run とは互換を取らない方針の範囲内）。

**U3-IMP-3 explore の出発構造を screen 極小の最適化構造にする** — U4-P2 と同じ（§2.4.7、§3 の X2）。

**U3-IMP-4 採用した極小の小さな虚振動を thermo で反転する** — U7-P1 に統合（§2.7.7、§3 の X3）。検証の修正案は「MinimumRecord なら invert を policy.saddle_cm1（50）に固定し、TS には適用しない」だった。U7-P1 はこれを含み、さらに TS の副次負モードまで扱う形にしたので、そちらを採る。

**U3-IMP-5 暗黙の分岐を 1 か所だけ明示する（修正版）** — could、中立
- 変更: level=dft かつ include=window で screen 極小が 0 件なら、`ValueError('include: window needs screen minima; use include: all')` にする。
- 任意: level=screen で select を明示したら、設定エラーにする（`model_fields_set` を使えば数行）。

**U3-IMP-6 射影の約束と開殻一重項の制限を 1 文ずつ書く** — could、中立
- vibrations.py の docstring か design §5.5 に「NWChem の Projected Frequencies とは低振動モードで数 cm⁻¹ ずれ、G では約 0.02 kcal/mol の差になる」と書く。
- design §7 に「mult 1 は RKS で、開殻一重項は検出しない」と書く。

#### 2.3.8 見送った案

- minima._pool で、explore の生成物を screen 極小と RMSD で照合する — 原因は explore の出発構造にあるので、そちらを直す（U4-P2）方が単純。
- mode-follow の両側を relax し直さない — 実績 0 回で、API の変更に見合わない。
- minima の opt だけ trust を 0.3 に戻す — 失敗 0 件で根拠がない。
- opt の閾値を緩める — 残差勾配で soft 押しが増え、かえって高くつく。
- freq で opt の movecs を使う — 節約は約 4% だけ。
- DFT freq を速くするために grid・CPHF・RI を変える — same_pes が崩れるか、PES そのものが変わる。
- 1 断片の分子にも xTB の初期 Hessian を使う — 効果の根拠がない。
- init_hessian の xTB freq を calc artifact に残す — JobStore からたどれる。下流が取り違えるおそれがある。
- Registry.find で組成も照合する — エネルギー条件で実質的に分離される。
- soft の押しを ± 両方向にする — 本物の鞍点なら片側で十分で、ノイズなら倍の無駄になる。
- DFT のエネルギー窓で freq を省く — freq がないと極小の確認ができない。
- 開殻一重項を自動で検出する — NWChem に簡便な手段がない。
- MinimumRecord に代表構造を持たせる — スキーマ全体に波及する。
- MinimumRecord のエネルギーを freq の値にする — 差は 1.4e-4 kcal/mol しかない。
- （U3-IMP-5 の元案）mode_follow と rerank_top の knob を削除し、rerank_sp を既定 true にする — 効果が小さく、YAML の書き換えが増えるだけ。

#### 2.3.9 変更不要な点

- opt と別ジョブの freq を同じ数値の層で計算し、is_minimum だけで判定すること。
- 質量加重の Eckart 射影、虚振動の 3 段方針、soft の片側押し、±0.1 Å の mode-follow。
- 置換不変 RMSD（proper rotation だけを許す）と、0.05 Å / 5e-5 Eh の閾値。
- 選択（6 kcal/mol、per_state 3、DFT SP で上位 8 を再順位付け、各組の最低構造を必ず残す）と、dft の開始構造を screen の最適化構造にすること。
- 2 断片以上の xTB 初期 Hessian。NWChem の既定の閾値、trust 0.1、maxiter 100、frame_shift の検査。
- xTB の vtight と 3 条件の収束判定。LADDER の継続。odft と spin_ok。

### 2.4 U4 反応探索（explore: ReaDuct NT2/AFIR、trial 生成、緩和による発見）

#### 2.4.1 現状仕様

**流れ**（`stages/explore.py`、`backends/readuct/{engine,worker}.py`、`chemistry/trials.py`）

1. **緩和による発見**: screen 極小の members のうち、状態ラベルが変わった seed を `DiscoveryRecord(mechanism=relaxation, outcome=negative, reason=collapsed_to:<label>)` として記録する。ReaDuct は呼ばない。
2. **出発点**: (組成, 状態ラベル) ごとに、エネルギーの低い screen 極小を 2 つ選ぶ。構造には代表 species の `species.geometry`（入力の xyz、RDKit、CREST 配座）を使い、**screen の最適化構造は使わない**（`explore.py:66`）。
3. **trial の生成**（`trials.generate`）。次の順に列挙する。
   1. polar_h: labile H を受容原子へ移す。2 段の relay を先に並べる。
   2. h_shift: 1,2- と 1,3-H 移動。
   3. heavy_bond: 重原子間の結合の生成と開裂。
   4. association: 断片間の結合。
   - 重複の除去は原子番号の集合でだけ行う。出発点あたり 10 件で打ち切り、切り捨てた件数は記録しない。
   - 直線分子は 10° 曲げ、各原子を 0.05 Å 乱す。
4. **実行**: trial ごとに NT2 を実行し、生成物が得られなかったときだけ、同じ drive で AFIR を実行する。
   - 1 attempt は 1 worker サブプロセスで、timeout は 600 s。
   - SCC が失敗したら、1000 K で attempt 全体をやり直す。
5. **worker の中の処理**
   - 共通: start を緩和して source（基準の構造とエネルギー）を作る。NT2 と AFIR は source ではなく、**未緩和の start** から走らせる。
   - NT2: `run_nt2_task` → Bofill TSOpt → 射影した振動数で虚振動が 1 本かを確かめる → IRC → 両端を opt して n_imag=0 を確かめる → `irc_product`。
   - 出発構造との一致: `matches_source`（置換不変 RMSD < 0.1 Å かつ結合グラフが同じ）で判定する。
   - AFIR: 1 原子対に γ=125 → 300 kJ/mol の順でかけ、無バイアスの opt と Hessian で確かめ、`matches_source` でなければ生成物とする。
6. **判定**（`product_verdict`）
   - NT2 の生成物で TS と IRC の検証がなければ `ts_not_validated`。
   - ΔE_rxn が不明または 100 kJ/mol 超なら、`out_of_window`。
   - ΔE‡ が 150 kJ/mol 超なら、`out_of_window`。
   - 採用した生成物は SpeciesRecord（source=discovery）になり、minima(dft) に必ず送られる。
7. **下流での使われ方**
   - hypotheses が生成物を仮説にする。
   - 同じ組成・同じ結合変化の negative を `negative_evidence` に集め、SCREEN がジョブなしで NO_PRODUCT にする（U5 を参照）。

**エンジン呼び出し（W7 の実物）**

- argv: `python -m hfauto.execution.worker hfauto.backends.readuct.worker:run_attempt <attempt>/job.json`。env は `OMP_NUM_THREADS=1,1`。
- job.json の settings: `max_scf_iterations 300, electronic_temperature_K 300, scc_retry_temperature_K 1000, afir_gamma 125/300, nt_total_force_norm 0.1, imag_cutoff_cm1 50, timeout_s 600`。
- ReaDuct の呼び出し:
  - `run_nt2_task(nt_associations=[1,13], nt_dissociations=[13,14], nt_total_force_norm=0.1)`
  - `run_tsopt_task(optimizer="bofill", automatic_mode_selection=[…])`
  - `run_irc_task(stop_on_error=False)`
  - `run_afir_task(afir_lhs_list=[1], afir_rhs_list=[5], afir_attractive=True, afir_energy_allowance=125.0)`
- ReaDuct と Utilities の暗黙の既定値（hfauto は設定しない）:
  - opt・IRC・AFIR の反復上限は 150
  - NT2 は last_maximum_before_first_target、最大 500 反復
  - AFIR の phase_in は 30

**利用者が設定する項目**: sources_per_state（2）、max_trials_per_source（10）、settings（上記 8 項目）、window（barrier_kj 150、reaction_kj 100）、site の readuct（版数 6.1.0、python、threads）、下流の `override_negative_evidence`。

**実測**

- amine_pilot2:
  - 102 s、44 attempt（NT2 25、AFIR 19）、失敗 0、採用した生成物 0。
  - TMA·HF の F→N プロトン移動の NT2 は、単調に +141.5 kJ/mol 上がり、TS の推定構造は得られなかった。
  - TMA 単量体の C→N 1,2-H 移動（イリド）は ΔE‡ 376 kJ/mol で窓外。
  - TMA·HF の 10 trial のうち 9 件は、等価なメチル H の移動だった。
- tma_hf2:
  - 83 s、37 attempt、失敗 0。
  - AFIR の「生成物」は 1 件。状態ラベルも結合グラフも出発構造と同じで、RMSD 0.586 Å の配座違いにすぎなかった。それでも DFT に送られ、opt と freq で 1,308 s（dft stage の 43%）を使った。
- 反復上限: IRC はすべて両方向とも 150 反復で打ち切られた。irc_end_not_minimum 10 件のうち 8 件は、端点の opt が 150 で打ち切られたもの。
- HCN の smoke: NT2 は 0.97 s、TS −1426i、ΔE‡ 306 kJ/mol（DFT では 195 kJ/mol）なので、既定の窓では窓外になる。

#### 2.4.2 化学的妥当性

- **エンジンの使い方は妥当です。** NT2 → Bofill TSOpt → 虚振動 1 本の確認 → IRC → 両端の opt と n_imag=0 は、SCINE Chemoton/Puffin の標準手順と同じです。
- **欠陥は判定と、下流での使い方にあります。**
  1. **陰性結果を DFT の前のゲートに使っていること（negative_evidence）。** 仮説の両端は DFT の別々の basin なので、山越えの定理により 1 次の鞍点が必ずあります。GFN2 の単端探索（障壁の MAD は約 58 kJ/mol、HCN では +111 kJ/mol 過大）で TS が見つからなかったことは、その反証になりません。
  2. **出発構造が screen の極小ではないこと。** TMA·(HF)₂ では、ReaDuct の source が screen 極小より 2.6 kJ/mol 高い別配座に落ち、AFIR の「生成物」と relay TS の IRC 両端が、実は screen 極小そのものでした。screen 極小から再実行すると、same_as_source になりました。
  3. **一致の判定が RMSD < 0.1 Å に依存すること。** 柔らかい H 結合錯体では、同じ結合グラフの別配座が「生成物」や「非接続」になります。Puffin はグラフで比べています。
- **標的化学との整合**: TMA·HF の F→N PT が +141.5 kJ/mol で単調に上がることは、気相の PA と酸性度の差から見ても妥当で、DFT の same_basin とも一致します。開殻は restricted_open_shell で扱い、黙って誤ることはありません。

#### 2.4.3 有用性

- コストは小さいです（stage 83〜102 s）。実際のコストは下流で出ます。誤った生成物は 1 件ごとに DFT で約 22 分かかり、誤った陰性証拠は正しい反応を捨てます。
- **必須**: polar_h の NT2（relay を先に）、h_shift、TS の検証、CH-28 の窓。
- **有用**: heavy_bond、relaxation の記録。
- **限界的**: AFIR。W7 では 36 件中、有用な生成物が 0 件でした。それでも 1 件 1.5 s で、NT2 の盲点を補う標準の手法なので残します。
- **逆効果**: negative_evidence のゲート。

#### 2.4.4 複雑さ

- 削れる、または置き換えられるもの:
  - (a) negative_evidence の一式（約 30 行と設定 1 個）。
  - (b) RMSD による一致判定。既存の state_label に置き換えられる。
- ReaDuct の反復上限は、設定にせずコードの 1 か所に明示すれば足ります。
- trial の距離の定数、SCC の 1000 K での再試行、直線分子の曲げは、今のままでよいです。

#### 2.4.5 ライブラリ仕様の要点

- **SCINE ReaDuct 6.1.0**（https://raw.githubusercontent.com/qcscine/readuct/6.1.0/src/Readuct/App/Tasks/IrcTask.h 、AfirOptimizationTask.h、NtOptimization2Task.h）
  - opt・IRC・AFIR の反復上限の既定は 150 で、`convergence_max_iterations` で変えられる。
  - NT2 は、失敗時に `.nt.failed.xyz` を書き、軌跡の各フレームのコメント行にエネルギーを書く。
- **SCINE Utilities 10.1.0 の NtOptimizer2**（https://raw.githubusercontent.com/qcscine/utilities/10.1.0/src/Utils/Utils/GeometryOptimization/NtOptimizer2.cpp）
  - "No transition state guess was found" は、平滑化した勾配に極大がないときに出る。平滑化で消えた小さな極大も含む。
  - このほかに、"NT micro iterations failed" と "Calculation in NT optimization failed" という別の失敗がある。
- **scine-xtb-wrapper 3.0.2**: SCC が収束しないと "Self consistent charge iterator did not converge" を出す（プローブで確認）。hfauto はこれを正しく検出できる。
- **SCINE Puffin**（https://raw.githubusercontent.com/qcscine/puffin/master/scine_puffin/jobs/scine_react_complex_nt2.py 、…/templates/scine_react_job.py）
  - NT2 ジョブの既定は、TSopt 200、IRC 50、IRC 端点の opt 200、opt 500。
  - IRC の端点は、molassembler のグラフで出発構造と比べる。
- **Chemoton**（https://scine.ethz.ch/static/download/documentation/chemoton/v4.0.0/py/api/utilities.html）: `get_atom_pairs` の prune で、対称的に等価な原子を除く。
- **AFIR / GRRM17**（https://pmc.ncbi.nlm.nih.gov/articles/PMC5765425/）: γ は越えられる障壁のおおよその上限で、熱反応なら約 107 kJ/mol。
- **GFN2-xTB の障壁の誤差**（https://benchmarks.rowansci.com/methods/gfn2-xtb）: GMTKN55 の障壁群で約 13.9 kcal/mol。

#### 2.4.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U4-I1 | high（U5-I1 と統合すると critical） | negative_evidence のゲートが、正しい反応を DFT の前に捨てます。(1) AFIR の生成物には、必ず同じ drive の NT2 negative が付きます。(2) HCN の NT2 は out_of_window になるので、宣言反応 hcn_to_hnc も止まります。(3) 数値的な不成立も、化学的な陰性と同じ扱いになります。W7 では 0 回しか発火していません。 |
| U4-I2 | high | explore の出発構造が screen の極小ではありません。TMA·(HF)₂ の偽の生成物と「非接続」判定は、これが直接の原因です。amine の TMA（RDKit の構造）は source より 33.2 kJ/mol 高く、同じ trial の結論が run によって変わりました。 |
| U4-I3 | medium | 生成物と IRC 接続の判定が RMSD < 0.1 Å に依存し、配座の違いを別の化学種として扱います。I2 を直しても、錯体ではこの基準の問題が残ります。 |
| U4-I4 | medium | 対称的に等価な drive を区別しないので、上限の 10 件が同じ種類の trial で埋まります。切り捨てた件数も記録されません（CH-29）。実害が大きいのは、等価な基を持つ一般の有機分子です。 |
| U4-I5 | low | 反復上限 150 が、数値的な陰性を作っています。プローブでは、上限を上げても窓内の生成物は増えませんでした。U4-P1 の後は、陰性が判定に使われる経路もなくなります。 |
| U4-I6 | low | NT2 の失敗をすべて `monotonic_uphill` と呼び、定量値を捨てています。TS がある negative の ΔE‡ も記録されません。 |
| U4-I7 | low | AFIR の再試行 γ=300 が障壁の窓 150 と矛盾する、という指摘。検証による訂正: γ は xTB の PES 上での目安にすぎず、GFN2 は障壁を過大に出すので、「窓と矛盾する」は言い過ぎです。害の実証もありません。 |
| U4-I8 | low | 1.45 Σr_cov の結合判定のため、TMA·(HF)₂ では F→N 移動の drive が作られません。設計上の判断で、本題は DFT の same_basin で決着しているので、参考として記録します。 |

#### 2.4.7 改良案

**U4-P1 negative_evidence ゲートを削除する** — must、【複雑さ減】（U5-P1 と同じ。§3 の X1）
- 変更: 次を削除する。
  - `actions.py` の SCREEN の即時 NO_PRODUCT と、その引数
  - `gates.barrier_verdict` の negative_evidence / override の引数と分岐
  - verdict の `"negative_evidence"`
  - `state.py` の `override_negative_evidence` と行 12 への写像
  - `ReactionRecord.negative_evidence`
  - `hypotheses._Pool.negatives` と `negative_evidence()`
  - 関連テスト 1〜2 件、design.md の 44・119・141 行
- 陰性結果は、summary.coverage で機構別・理由別に報告し続ける。
- 効果: DFT で ΔE‡ 195 kJ/mol の素反応 hcn_to_hnc が、discover 系で NO_PRODUCT になることを防ぐ。約 30 行と設定 1 個が減る。W7 の結果は変わらない。
- 検証の指摘: P6 の簡素化版を一緒に入れれば、陰性結果の情報量は落ちない。

**U4-P2 explore の出発構造を screen の最適化構造にする** — must、中立（= U3-IMP-3）
- 変更: `ExploreStage.run` から `_Explorer.source` に `inputs` を渡し、`inputs.evidence(minimum.opt_calc).final` を読むようにする（`minima._pool` と同じ取り方、約 5 行）。trials と perturb_linear もこの構造に適用する。worker は変えない。
- 効果: W7 の偽の生成物（1,308 s）と afir_not_converged 1 件が、どちらも same_as_source になる（プローブで確認）。run ごとに結論が変わる問題もなくなる。
- 注意: JobStore のキーが変わるので、旧 run のキャッシュは再利用できない。

**U4-P3 出発構造との一致を結合グラフ（state_label）で判定する** — should、【複雑さ減】
- 変更: `matches_source` を `state_label(symbols, end) == state_label(symbols, source)` に置き換え、`SOURCE_RMSD_A` と RMSD の import を削る。
- 効果: W7 の 3 件（配座違いの生成物 1 件、irc_not_connected 2 件）が、正しく same_as_source になる。Puffin と同じ基準で、hfauto の basin の定義とも揃う。
- 注意: 短い H 結合（N···H 1.43〜1.48 Å）の伸び縮みでラベルが変わる配座違いは、まだ「生成物」になりうる（下流の DFT の basin 判定で吸収される）。立体異性体は同じラベルになるが、これは conformer ペアの役割なので問題ない。

**U4-P4 drive の重複を WL の原子クラスで除く（簡素化版）** — could、中立
- 変更: `_wl_hash` の反復ラベルを原子ごとに返す関数を公開し、重複除去の鍵を原子クラスで正準化した `(frozenset(assoc), frozenset(dissoc))` にする。列挙数・固有数・実行数は、stage のログか notes に 1 行出すだけにし、スキーマは変えない。
- 効果（プローブ）: 列挙数から固有数への削減は、TMA 30→5、TMA·HF 35→7、NH3·HF 5→3、TMA·(HF)₂ 51→15。効果が出るのは、等価な基を持つ一般の有機分子。

**U4-P5 無バイアスの opt に反復上限 500 を明示する** — could、中立
- 変更: `_Scine.minimum` の `run_opt_task` に `convergence_max_iterations=500` を 1 か所で渡す。設定ファイルの項目にはしない。AFIR の stop_on_error は変えない。
- 効果: 端点の opt に 150 反復以上かかる真の生成物を取り逃がさないための保険になる。1 attempt は最大約 17 s になる。

**U4-P6 陰性結果のラベルを最小限だけ正確にする（修正版）** — could、中立
- 変更:
  1. `run_nt2_task` の RuntimeError のメッセージを見て、"No transition state guess" なら `no_ts_guess`、それ以外は `nt2_failed` にする。
  2. TS がある negative にも、`_nt2` の energies に ts を入れ、既存の `barrier_kj_mol` に ΔE‡ を載せる。
- ΔE_stop を reason 文字列に埋め込む案は採らない。

#### 2.4.8 見送った案

- **U4-P7 AFIR の γ を 1 つにする（棄却）** — γ は xTB の PES 上での目安にすぎず、GFN2 は障壁を過大に出す。γ=300 の害は実証されていない。減るのは設定 1 個と約 1 s/件だけで、既存の決定（CH-32）を 2 系統のデータで覆すには根拠が弱い。
- negative_evidence を理由や機構で絞って残す — 規則が増えるだけで、そもそも拒否権にする根拠がない。
- AFIR を stop_on_error=False にする — 配座違いの「生成物」が 3 件とも出た。
- AFIR のフォールバックを削除する — 安価で、NT2 の盲点を補う標準の手法。
- 複数の原子対を使う AFIR — 誤った引力が入る。複数の結合を同時に動かすのは NT2 の役割。
- 新しい探索エンジン（B-spline、msreact/MTD、グラフ列挙 + GSM、MC-AFIR） — SCREEN や conformers と重なり、重く、不要。
- explore の並列化、deadline の導入 — stage は 100 s 未満。
- 障壁の窓を 200 kJ/mol に広げる — P1 の後は、宣言反応を拒否しない。
- NT2 の抽出基準を設定として公開する — プローブで結果が同じだった。
- IRC の反復を調整する — 端点は後で opt するので問題にならない。
- 1.45 Σr_cov の閾値を変える、片側だけを切る drive を加える — topology 全体に波及する。
- h_shift に原子価のフィルタを付ける — P4 で足りる。
- SCC の 1000 K 再試行を削除する — 標準的な対処で、実装も小さい。
- スピンモードを unrestricted に切り替える — wrapper が受け付けず、GFN2 には効果もない。
- （U4-P6 の元案）guess.nt.trj.xyz の ΔE_stop を reason 文字列に埋め込む — 数値を理由文字列に入れるのは筋が悪い。

#### 2.4.9 変更不要な点

- NT2 → Bofill TSOpt → ν < −50 cm⁻¹ が 1 本であることの確認 → IRC → 両端の opt と n_imag=0、という検証の流れ（CH-26/27 の修正は正しい）。
- 自前の射影振動数による虚振動の計数、CH-28 の窓と TS 検証の要求、out_of_window を found とみなして AFIR を省くこと。
- trial の段の順序（relay を最初に）、直線分子の曲げ、worker サブプロセスと timeout、版数の pin。
- SCC 失敗の検出と、1000 K で再試行したあと 300 K で SP し直すこと。restricted_open_shell。
- sources_per_state 2、max_trials 10、relaxation の記録と coverage の集計、低レベル TS を xTB freq で確かめてから近道に使うこと、単一原子対の AFIR。

### 2.5 U5 反応経路: 仮説の選択と障壁の事前判定（hypotheses.select / SCREEN / FIND_PATH）

#### 2.5.1 現状仕様

**仮説の選択**（`chemistry/hypotheses.py`）

- **宣言反応**: 常に残す。端点の basin は `endpoint_basin` で DFT の basin に写す（W7 の F1 修正）。
- **explore の生成物**: TS のあるもの、ないもの、mode_follow の順に並べる。
- **conformer の組**: 同じ組成・状態ラベル・level_key を持つ DFT 極小の全ペア。
- **採用の条件**: 両方が tier=dft で、level_key と組成が同じで、basin が別であること。さらに ΔE ≤ 40 kcal/mol で、次のどちらかを満たすこと。
  - 結合の変化、ねじれの変化 ≥30°、距離の変化 ≥0.2 Å のいずれかがある。
  - conformer の場合は、結合の変化がなく、ねじれの変化が ≥30° である。
- **端点の組**: members の中で、恒等写像の Kabsch RMSD が最小になる組を使う。
- **件数**: 組成あたり 6 件まで（宣言反応も数える）。
- **degenerate**: 同じ basin で、`mapped_equivalent` が成り立つもの。
- **torsional の既定**: 結合の変化がなければ True。
- **negative_evidence**: 同じ組成で、結合変化が一致する negative discovery を集める。

**決定表**（`drivers/reaction_case/state.py`、17 行を上から評価する純関数）

- 行 11（SCREEN）: 行 1〜10 がどれも当てはまらないときに選ばれるので、実際にはケースの最初の action になる。
- 行 12: barrierless なら BARRIERLESS、negative_evidence なら NO_PRODUCT。
- 行 13: multi_max なら VALIDATE_INTERMEDIATE。
- 行 14: monotonic が 2 回続けば BARRIERLESS。
- 行 15: 種があれば REFINE_SADDLE。
- 行 16: DFT string をまだ実行していなければ FIND_PATH。
- 行 17: UNRESOLVED。

**SCREEN**（`actions.py`）

- 0. negative_evidence があり override がなければ、ジョブなしで negative_evidence を返す。
- 1. **近道**: 低レベル TS があれば、xTB freq（−50 cm⁻¹ 未満が 1 本）と DFT SP 1 点による 3 点で判定する。proceed のときだけ採用する。
- 2. DFT の両端を xTB で再最適化する（`--opt vtight`）。
- 3. xTB で両端が崩壊した場合（写像 RMSD ≤0.05 Å）: DFT の両端を IDPP で 11 点に補間し、全点で DFT SP を計算する。種は作らない。
- 4. **GS**: pysisyphus の GS（11 ノード、climb）と、同じ入力内の TSOpt（rsirfo / gau）。GS が収束したかどうかは判定に使わない。
- 5. 11 ノード全部と TS で DFT SP を計算する。
- 6. **種**: xTB freq を通った GS TS なら `screen_ts`（接線は chord）、通らなければ DFT HEI ノード（`screen_hei`、接線は中心差分）。
- 7. **判定**（`gates.barrier_verdict`）: rel = 内部点の最大値 − 両端の高い方（DFT の基準は DFT 極小のエネルギー）。rel_dft ≥1.0 または rel_low ≥1.0 なら proceed、どちらも <1.0 なら barrierless。
  - below_zpe = rel_dft < ½|ν_tangent|。ν_tangent のために、高い方の端点の DFT freq を要求する。
- 8. unavailable、または種のない proceed は、FIND_PATH に進む。

**FIND_PATH**

- 初期経路: SCREEN の経路（GS か IDPP）があればそれを、なければ新しく IDPP を作る。Cartesian の弧長で bead 数に再標本化し、各点を ends[0] に整列し、両端を DFT の端点に置き換える。
- チャンク: NWChem ZTS を最大 3 回（1 回ごとに別ジョブ）。次のどれかで打ち切る。
  - Failure（monitor が停滞と判定して止めた STAGNATED を含む）
  - `string_converged`: 最終 gmax ≤1e-3 かつ直近 3 反復で悪化なし
  - `stagnated`
- monitor: 実行中に 10 s ごとに gmax / xmax の履歴を見て、停滞ならプロセスを kill する。
- 形の判定: single_max は収束しなくても HEI を種にする。monotonic と multi_max は string_converged のときだけ記録し、それ以外は 'failed' にする。

**エンジン呼び出し（W7 の実物と本レビューのプローブ）**

- xTB: `xtb input.xyz --opt vtight --gfn 2 --chrg 0 --uhf 0`（0.06 s）。
- pysisyphus の run_dict（HCN）: `geom {type: cart}`、`calc {type: xtb, gfn 2, pal 1}`、`cos {type: gs, max_nodes 9, climb true, climb_rms 0.005}`、`opt {type: string}`、`tsopt {type: rsirfo, thresh: gau}`。
  - HONO と NH3 は `geom.type: dlc`。2 断片以上なら、tsopt の geom を tric にする。
  - 厳密に直線の分子では、±0.01 Å ずらす。
- DFT SP: `task dft energy`（1 点 1.7〜2.9 s）。
- ZTS（W7 では未実行。プローブで hfauto の render_string を使って作ったデッキ）:
  ```
  geometry … end
  geometry endgeom … end
  string ; nbeads 9 ; maxiter 20 ; stepsize 0.05 ; interpol 3 ; tol 1e-5
           freeze1 .true. ; freezeN .true. ; impose ; xyz_path initial_path.xyz ; end
  task dft string
  ```

**主な既定値**: 窓 40 kcal/mol、変化 0.2 Å / 30°、6 件/組成、_COLLAPSE_A 0.05 Å、barrier_proceed_kcal 1.0（SCREEN の判定と、形の判定の解像度を兼ねる）、screen_images 11、string_beads 9、confirm_barrierless_beads 15、string_chunks 3、max_path_runs 2、停滞の定義（10 反復で 5% の改善がない、または |xmax| < 1e-6 が 5 反復）、予算は 1 反応 6 h。

**利用者が設定する項目**: paths の method / engines / screen、policy（上記と override_negative_evidence、walltime_h）、reaction_ids、gates（reaction_window_kcal、barrier_proceed_kcal、saddle_cm1）、system の reactions、site の pysis_gs と nwchem_string。

**実測**

- どの run でも、仮説は宣言反応の 1 件だけだった。tma_hf2 と water は行 2 の same_basin になった（paths のジョブ 0 件）。
- SCREEN の GS 経路（すべて proceed で、種は screen_ts）:

| 系 | GS | SCREEN の区間 | max_rel_dft / low | 精密化後の ΔE‡ |
|---|---|---|---|---|
| HCN | cart、44 サイクルで収束、6.1 s | 約 28 s | 33.98 / 53.20 | 46.62 |
| HONO | dlc、150 サイクルで**未収束**（TSOpt は収束）、15.2 s | 約 48 s | 14.21 / 11.63 | 13.67 |
| NH3 | dlc、7 サイクル、2.6 s | 約 29 s（tangent 用の DFT freq が JobStore をミスして 1 本余分） | 4.46 / 6.14 | 4.26 |

- 未実行の経路: 低レベル TS の近道、IDPP、FIND_PATH、行 12・13・14。
- ZTS のプローブ（HCN、GS 経路から 9 bead、負荷のため STO-3G・1 rank）:
  - 20 反復、77.5 s で、NWChem は "failed to converge" を出した。
  - gmax は 0.044〜0.05 で横ばいだった。hfauto の判定では string_converged=False、stagnated=True で、monitor は 14 反復目に kill する。
  - 一方、bead エネルギーの変化は、8 反復目以降 0.07 kcal/mol 未満で安定しており、プロファイルは single_max だった。
- 改修前の実 DFT の HCN string（92 反復）でも、gmax は 0.0369 で横ばいだった。

#### 2.5.2 化学的妥当性

- **仮説の選択は妥当です。** 次の点は、どれも化学的に筋が通っています。
  - 同じ level_key の別 basin の組だけを選ぶ。
  - 40 kcal/mol の窓（K ≈ 1e-29）。
  - 恒等写像の RMSD が最小になる端点の組（置換経路を避ける）。
  - F1 の修正と、degenerate の判定。
- **negative_evidence は論理的に破綻しています。** 行 1 により両端は DFT の極小なので、生成物の basin は必ず存在します。それなのに、低レベルの陰性結果だけで、ジョブなしに NO_PRODUCT（no_product_basin）で終えます。実データでも、TMA·(HF)₂ の唯一の product には、同じ trial の NT2 の陰性記録が付いていました。テスト（`test_hypotheses.py:47-65`）は、TS まで見つかっている HCN→HNC に negative_evidence が付くことを期待値として固定しています。
- **SCREEN は妥当です。** GSM と HEI から TS を最適化する標準的な流れで、判定の基準に DFT 極小を使うので、経路の最大値は minimax の意味で saddle の上界になります（barrierless の判定は保守的）。
  - 錯体では DFT//xTB の端点のずれ（HONO で 0.93〜1.22）が解像度 1.0 と同程度なので、barrierless は出にくく、結論は FIND_PATH に委ねられます。
  - 判定は形を見ないので、両端の極小より低いノード（未知の basin の厳密な証拠）を見逃すことがあります。
- **FIND_PATH は機能していません。** NWChem の gmax は、射影しない全勾配の最大成分です（`string.F:744-750`）。MEP 上でも接線成分が残るので、障壁や傾きのある経路では gmax ≤1e-3 は原理的に満たされません（E, Ren, Vanden-Eijnden 2007）。そのうえ、gmax の横ばい（ほぼ収束した状態）を停滞とみなして kill し、`string_final.xyz` はループの後にしか書かれないので、結果が捨てられます。
- **DLC の GS 画像**: pysisyphus は Cartesian のときしか画像を整列しません。このため、左右の枝の境目で座標系が回転します（HONO では生の距離 1.149 Å、整列後 0.13 Å）。
- **IDPP の実装**は Smidstrup 2014 と一致しており、正しいです。

#### 2.5.3 有用性

- **必須で費用対効果が高いもの**
  - SCREEN の GS 経路（小分子で 30〜50 s、17 原子で約 5 分と推定）。W7 の全 run で、種から接続までが最短で通った。
  - 低レベル TS の近道、IDPP の経路（xTB が崩壊させる場合の唯一の手段）。
  - 仮説の選択（件数の管理に必須）。
- **価値が負のもの**: negative_evidence（節約できるのは SCREEN 1 回分の数分だけなのに、実在する反応を捨てる）。
- **現状の価値が 0 のもの**: FIND_PATH。実行されれば、ほぼ UNRESOLVED になる。17 原子では 1 チャンク約 1 h、15 bead の確認は 1 チャンク 1.7〜2 h で、予算内に終わらない。
- **小さいもの**: below_zpe と tangent_mode_cm1（誰も読まないのに DFT freq を要求する）、max_node_spacing_A（DLC では誤った値）。

#### 2.5.4 複雑さ

- 削れるもの:
  - (a) negative_evidence の一式（約 40 行、verdict の種類が 1 つ減る）。
  - (b) FIND_PATH の gmax 系の判定（string_converged、stagnated、_string_monitor）。bead エネルギーの安定度 1 つで足りる。
  - (c) 行 14 の多解像度での確認（confirm_barrierless_beads、max_path_runs）。
  - (d) tangent_mode_cm1 の即時計算。
- 今のままでよいもの:
  - barrier_proceed_kcal を 2 つの用途で共有すること（同じ「解像度 1 kcal/mol」という概念）。
  - 仮説の閾値などをコードに固定していること。
- 設定の罠: site の ExecutionSpec.maxiter は、nwchem の anchor で opt・saddle・string に共通です。

#### 2.5.5 ライブラリ仕様の要点

- **NWChem 7.2.3 の ZTS**（https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F 、string_input.F）
  - gmax は射影しない（projection1 の既定は .false.）。
  - `tol` は変位（xrms と xmax）の閾値で、勾配の閾値ではない（L23 のコメント、L879）。
  - 初期 gmax ≤0.5 なら LBFGS で始まる。
  - 反復ごとに `string: Path Energy #n` を出力する。
  - `string_final.xyz` はループの後にだけ書かれる（L897）。
  - `impose` は新しく鎖を作るときにしか効かない。
- **NWChem の公式文書**（https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html）: tol を「最大勾配」と説明しているが、ソースとは違う。
- **pysisyphus 1.0.0**（https://pysisyphus.readthedocs.io/en/latest/chainofstates.html と venv のソース）
  - `run.py:437`: splined HEI に負の固有値がないと、汎用の Exception を出す。`run_from_dict` が捕まえるのは HEIIsFirstOrLastException だけ。
  - `optimizers/Optimizer.py:1007`: COS の整列は Cartesian のときだけ。
- **IDPP**: Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. 140, 214106 (2014)。
- **ZTS の理論**: E, Ren, Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007)。CI-NEB: Henkelman, Uberuaga, Jónsson, J. Chem. Phys. 113, 9901 (2000)。
- **GFN2 と PBE0 の構造のずれ**（W7 の実データ）: xTB 構造での DFT エネルギーは、DFT 極小より TMA·(HF)₂ で 1.19、HONO で 0.93/1.22、HCN で 0.52/0.24、NH3 で 0.02 高い。

#### 2.5.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U5-I1 | critical | negative_evidence のゲートが、DFT の basin がある生成物を、ジョブなしで NO_PRODUCT（no_product_basin）にします（U4-I1 と統合。§3 の X1）。AFIR で見つかった生成物は、別 basin に残れば構造的に必ず誤分類されます。 |
| U5-I2 | high | FIND_PATH の `string_converged` は、射影しない gmax ≤1e-3 を要求するので、障壁のある経路では成り立ちません。monotonic と multi_max は 'failed' になり、行 13・14 に実質的に到達しません。FakePath（収束を gmax 3e-4 と仮定）がこの挙動を隠しています。 |
| U5-I3 | high | string の monitor が、gmax の横ばいを停滞とみなして kill し、結果ごと捨てます。検証の追加: チャンクのループにある `string_converged(...) or stagnated(...)` による break（`actions.py:475`）も同じ規則なので、monitor だけを外しても正常終了したチャンクの後で打ち切られます。 |
| U5-I4 | medium | DLC の GS 画像が整列されていません。検証による範囲の訂正: verdict のエネルギーと screen_ts（接線は chord）には影響しません。影響があるのは、max_node_spacing の記録値、FIND_PATH の初期経路（最短原子間距離 0.78 Å）、screen_hei の接線（TS が不合格のときだけ）です。 |
| U5-I5 | medium | pysisyphus の TSOpt が例外で終わると、GS の経路ごと失われ、SCREEN が unavailable になり、高価な FIND_PATH に落ちます。実例はありません（起こりうる）。 |
| U5-I6 | medium | 15 bead の barrierless 確認（行 14）は、17 原子では 6 h の予算内に終わりません。検証による補足: FIND_PATH に来るケースでは、DFT string の monotonic は低レベルとの矛盾を解く証拠になります。それでも、どの DFT 経路の最大値も saddle の上界なので、1 本で結論してよい、という結論自体は正しいです。 |
| U5-I7 | low | SCREEN の判定は形を見ないので、両端より低い DFT//xTB ノードを見逃して BARRIERLESS で閉じることがあります。検証による補足: 両端より低くない凹み（例: 0→5→2→10）も取りこぼしますが、未緩和の経路では厳密な証拠にならないので、対象にしません。実例はありません。 |
| U5-I8 | low | tangent_mode_cm1 用の DFT freq を proceed のときも毎回要求します。below_zpe は下流の誰も読みません。 |
| U5-I9 | low | DFT//xTB の端点のずれが解像度と同程度なので、SCREEN の barrierless は出にくくなっています（安全側）。「錯体ではほぼ出ない」は推測で、変更は不要です。 |
| U5-I10 | low | NWChem の文書とソースで tol の意味が違います。hfauto は NWChem の converged を使わないので、実害はありません。設計書の記述だけを直します。 |

#### 2.5.7 改良案

**U5-P1 negative_evidence ゲートを削除する** — must、【複雑さ減】（U4-P1 と同じ。§3 の X1）
- 削除の範囲は U4-P1 と同じ。テストの期待値は、`test_hypotheses.py` の該当行と `test_reaction_paths_stage.py:107-124` の該当部分を直す。

**U5-P2 FIND_PATH の判定を置き換える（簡素化版）** — must、【複雑さ減】
- 変更:
  1. gmax の窓による停滞判定と、string の monitor（`_string_monitor` と、engine.monitor の string 分岐）を削除する。チャンクのループの `stagnated` による break も外す。
  2. 形の記録に付いている収束ゲート（`actions.py:482-484`）をやめ、single_max と同じく、どの形もそのまま記録する。
     - 根拠: 連続な DFT 経路のエネルギー最大値は、収束していなくても minimax で saddle の上界になる。multi_max は VALIDATE_INTERMEDIATE が DFT 最適化で確かめるので、見かけの凹みでも安全。
  3. チャンクを続けるかどうかだけを、bead エネルギーの安定度で決める。`string: Path Energy #n` のブロックを読む関数（約 10 行）を足し、直近 3 反復の max|ΔE_bead| < 0.1 × resolution（0.1 kcal/mol）なら打ち切る。
  4. `string_converged`、`stagnated`、gmax_history の string 用途を削除する。
  5. FakePath を「gmax は横ばい、エネルギーは安定」に直し、G13/G14 の golden をエネルギー版に置き換えるか、削除する。
- maxiter を 10 に、チャンクを 4 にする変更は行わない（20 × 3 のまま）。
- 効果: FIND_PATH が初めて機能する。プローブでは 8 反復目以降に安定し、1 チャンクで終わる。収束した string を kill して捨てることもなくなる。
- コスト: 約 30 行（削除の方が多い）とテストの差し替え。

**U5-P3 GS の画像を逐次 Kabsch で整列する** — should、中立
- 変更: 置き場所は次のどちらか。GrowingStringAdapter の profile で images.xyz を書く前に整列する方が、成果物自体が整合するのでよい。
  - (a) GrowingStringAdapter の profile で、images.xyz を書く前に隣接画像へ整列する。
  - (b) `actions._screen_path` で、frames を読んだ直後に `align_mapped(prev, f)` を順に適用する。
- 回転の跳びがある経路のテストを 1 件加える。P2 と同時に入れる。
- 効果（HONO）: max_node_spacing が 1.149 → 0.19 Å、初期 bead の原子間距離のずれが 0.18 → 0.05 Å、最短原子間距離が 0.78 → 0.96 Å。判定には影響しない。

**U5-P4 行 14 の 15 bead 確認をやめる** — should、【複雑さ減】（P2 の後に入れる）
- 変更: `_r14_monotonic` を「path_runs の末尾が monotonic なら BARRIERLESS(dft_path_monotonic)」にする。`confirm_barrierless_beads`、`max_path_runs`、confirm の分岐、`_POLICY_KEYS` の該当キーを削除する。
- 効果: 17 原子で約 5〜6 h を省ける。予算切れで UNRESOLVED になっていたケースが BARRIERLESS で閉じ、判定の基準が SCREEN と揃う。

**U5-P5 pysisyphus の TSOpt が例外で終わっても GS の経路を返す** — should、【複雑さ増】
- 変更: worker で `run_from_dict` が例外を出し、かつ tsopt があるときは、tsopt を除いた run_dict でやり直し、ts=None と注記を返す（約 6 行）。GS 自体の失敗まで握りつぶさないこと。例外を注入する単体テストを 1 件加える。
- 効果: screen_hei の種に自然に切り替わり、FIND_PATH（17 原子で 1 h 以上）に落ちない。GS のやり直しは 3〜15 s。

**U5-P6 below_zpe の接線モードを遅延計算にする（修正版）** — could、中立
- 変更: tangent なしで `barrier_verdict` を呼び、barrierless のときだけ `_tangent_mode_cm1` を計算して below_zpe を付け直す（約 5 行）。ステップの取り方の改良は任意。
- 効果: proceed のたびに DFT freq を要求しなくなる（17 原子で JobStore がミスすると約 900 s）。

**U5-P7 両端より低い DFT//xTB ノードがあれば barrierless にしない（P2 の後に限る）** — could、【複雑さ増】
- 変更: `min(内部) < min(端点) − limit` なら `unavailable('dft_node_below_endpoints')` とし、FIND_PATH → multi_max → VALIDATE_INTERMEDIATE に回す（3 行）。P2 の前は入れない。

**U5-P8 FakePath を実エンジンの挙動に合わせ、string の smoke を 1 本加える（修正版）** — should、中立
- 変更:
  - (1) 主な対策として、FakePath の gmax を横ばいにし、bead エネルギーを数反復で安定させる。これで I2・I3 型の不具合を unit テストで捕まえられる。
  - (2) 実 NWChem の smoke は、opt-in で、軽い基底（STO-3G、1 rank、20 反復で約 78 s）か 3-21G にする。def2-SVP で 2〜5 分という見積もりは楽観的。

#### 2.5.8 見送った案

- NWChem の projection1 で gmax を射影させる — 接線への正しい射影ではなく、最適化の力そのものが変わってしまう。
- NWChem の converged（xrms、xmax）を、大きな tol で使う — 揺らぎで横ばいになり、偽の収束も防げない。
- print_shift で途中の xyz を出させて回収する — 処理が増える。チャンクの境目で判定すれば足りる。
- SCREEN の DFT の基準を DFT//xTB の端点にする — 上界の保証がなくなり、barrierless の誤判定を招く。
- 判定に使っていない両端 2 点の SP を省く — 節約は 17 原子で約 40 s だけ。
- xTB で結合トポロジーが変わったら崩壊とみなす — 実例がない（TMA·(HF)₂ で RMSD 0.061 Å、結合の変化なし）。
- screen_hei を 2 本目の種として積む — saddle の試行 2 回のうち 1 回を、未緩和の HEI が消費してしまう。
- 仮説の閾値や GS / IDPP の定数を YAML で変えられるようにする — 根拠となるデータがない。
- H 結合の組み換えを conformer の仮説にする — 件数の予算を圧迫する。
- GS の SCC で電子温度を上げて再試行する、SP の失敗を 1 点だけ許容する — 失敗 0 件で、効果がはっきりしない。
- screen_images を増やす — 整列後の間隔は 0.08〜0.33 Å で、すでに十分。
- 隣のノードの movecs を初期推定に引き継ぐ — JobStore のキーが複雑になる。
- FIND_PATH を pysisyphus と NWChem の DFT GS/NEB に置き換える — CH-01/02 に反し、変更が大きすぎる。
- （U5-P2 の元案）STRING_MAXITER を 10、チャンクを 4 にする — 根拠がない。
- （U5-P6 の元案）below_zpe の一式を削除する — CH-06 の注記を残すなら、遅延計算で十分。

#### 2.5.9 変更不要な点

- 仮説の条件（別 basin・40 kcal/mol・0.2 Å / 30°・6 件・優先の順番）、pick_endpoints、driver._endpoint、endpoint_basin（F1）、degenerate の判定、_lend_ts。
- SCREEN の骨格、DFT 極小を基準にした minimax の判定と低レベルとの OR、barrier_proceed_kcal=1.0。
- xTB freq による低レベル TS の受理と、同じ Hessian の再利用。低レベル TS の近道。GS の収束を判定に使わないこと。
- 座標系の選び方（cart / dlc / tric）、IDPP の実装と崩壊時の DFT IDPP 経路。
- HEI を種にするのは single_max のときだけにすること。xyz_path によるチャンクの継続。決定表の構造。

### 2.6 U6 反応経路: 鞍点の精密化・TS 検証・接続確認（REFINE_SADDLE / VALIDATE_TS / CONNECT(QRC)）

#### 2.6.1 現状仕様

**流れ**（`drivers/reaction_case/actions.py`、`state.py`、`chemistry/gates.py`、`modes.py`、`classification.py`）

1. **決定表の行**
   - 行 5: 接続が確定したら完了。connection が failed なら、2 回目の振幅でもう一度だけ CONNECT。
   - 行 6: SaddleClaim があれば CONNECT。
   - 行 7: saddle が収束して未検証なら VALIDATE_TS。
   - 行 8: collapsed なら VALIDATE_INTERMEDIATE。
   - 行 9: 別の basin なら MULTI_STEP。
   - 行 10: higher_order で試行が残っていれば REFINE_SADDLE。
   - 行 15: 種があり、試行が 2 回未満なら REFINE_SADDLE。
   - 行 16: string が未実行なら FIND_PATH。
   - 各判断の前に、1 反応 6 h の Deadline を調べる。
2. **REFINE_SADDLE**
   - 1 回目は、種での xTB `--hess` を初期 Hessian にする。条件は、ν < −10 の負モードがちょうど 1 本で、direction との |cos| ≥0.3 であること。
   - 2 回目以降、または xTB が使えないときは、DFT の freq を使う。
   - direction は、宣言座標の勾配か種の接線。重なりが最大の虚モードを `mode_index` にする（重なり 0.3 未満は注記だけで、合否には使わない）。
   - NWChem saddle を実行する。`mode_index ≥1` なら `moddir k+1` を指定し、`noautoz` にする。
3. **VALIDATE_TS**: saddle.final で、別ジョブの DFT freq を計算する。`is_first_order_saddle` の条件は次のとおり。
   - 指紋の一致、same_pes(numerics)、振動数の本数。
   - 最低モードが −50（torsional なら −20）未満であること。
   - 2 本目に −50 未満のモードがないこと。
   - prominence ≥ max(2e-5 Eh, 20×scf_tol)。
   - 分岐: higher_order なら、2 本目のモードに沿って +0.1 Å 動かした種でやり直す。no_imaginary_mode なら collapsed。
4. **CONNECT（QRC）**
   - 振幅: κ = uᵀHu から s = √(2·target/κ) を求め、[0.05, 0.4] Å に収める。目標は max(3·drop, 3e-4 Eh)。2 回目は ×2 にする（クリップの後に掛ける）。
   - ±s に変位させた構造を DFT で opt する（**初期 Hessian なし**、trust 0.1、maxiter 100）。
   - 割付けは `Registry.find` で行う。見つからなければ、torsional の反応は `periodic_nearest`、それ以外は `_register`（relax_to_minimum → Registry.add）。
   - 接続ゲート: 両側で same_pes、traj[0] ≤ E_TS−drop、max(traj) ≤ E_TS−drop、traj[-1] ≤ traj[0]−drop を満たすこと。結果は elementary / degenerate / reassigned / failed。
5. **分類と分割**: REASSIGNED は、実際につながった極小の組に置き換える。collapsed と multi_max は中間体を登録し、R→I と I→P に分ける（深さ ≤2）。

**エンジン呼び出し（W7 の実物）**

- xTB の Hessian: `xtb input.xyz --hess --gfn 2 --chrg 0 --uhf 0`（0.016 s）。
- NWChem saddle（HCN、`jobs/75/751e…/job.nw`）: 前述の dft ブロックに `driver ; maxiter 50 ; trust 0.1 ; sadstp 0.1 ; inhess 2 ; xyz final ; end` と `task dft saddle`。job.hess は xTB の Hessian。出力は `moddir 0`、`firstneg T`、`Initial step taken uphill`。
- moddir の実行例（smoke の平面 NH2OH）: `noautoz`、`moddir 2`。17 ステップ、11.8 s。5 ステップ目以降は、負の固有値が 1 本になってモード 1 を追った。
- TS の freq: `task dft frequencies`。
- QRC の opt: `driver ; maxiter 100 ; trust 0.1 ; xyz final ; end`。inhess はない（Using diagonal initial Hessian）。継続するときは `vectors input job.movecs` と drv.hess を引き継ぐ。

**主な既定値**: max_saddle_attempts 2、walltime 6 h、QRC の目標 3e-4 Eh・範囲 0.05〜0.4 Å・再試行 ×2・2 振幅、drop max(1e-5, 20×scf_tol)、noise / saddle / torsion の閾値 10 / 50 / 20 cm⁻¹、prominence 2e-5 Eh、_MIN_OVERLAP 0.3、_RETRY_A 0.1 Å、max_split_depth 2。

**実測**（4 run とも、流れは screen → refine_saddle(screen_ts) → validate_ts → connect → complete の 1 通りだった）

| 系 | saddle | TS freq | QRC | この単位の時間 |
|---|---|---|---|---|
| HCN | 5 ステップ 4.6 s | −1128.5i、5.9 s | 振幅 0.0543 Å、両側 19 / 24 ステップ | 50.6 s（stage の 64%） |
| HONO | 5 ステップ 10.2 s | −681.5i、15.2 s | 振幅 0.0879 Å。片側は 100 ステップで maxiter に達して継続し、合計 126 ステップ 269 s | 336 s（stage の 87%、QRC は 81%） |
| NH3 | 3 ステップ 2.9 s | −758.0i、5.1 s | 振幅 0.05 Å（下限）、両側 50 ステップ | 78 s（stage の 72%）、degenerate |

- 失敗、higher_order、collapsed、QRC の再試行、`_register`、reassigned、分割、moddir、DFT Hessian へのフォールバックは、実計算では一度も起きていない。
- プローブ（同じ開始構造から、TS freq の job.hess を `inhess 2` で渡した場合）:
  - trust 0.3 と TS の Hessian: HONO 12、NH3 24、HCN 12 / 15 ステップ。最終エネルギーは本番と ≤2e-7 Eh で一致し、どれも最大エネルギーは開始点のままだった。
  - trust 0.1 と TS の Hessian: HONO 38、NH3 45、HCN 20 ステップ。
  - trust 0.3 で対角推定のままにすると、HCN の片側は TS より +0.199 Eh 上まで登った。
  - NWChem の MEPGS（NH3）: 前進・後退とも 9〜10 点で収束し、エネルギー評価は合計 50 回だった。

#### 2.6.2 化学的妥当性

- **妥当です。** 種 → 近似 Hessian を使った P-RFO → 別ジョブの DFT freq を自前で射影 → QRC → 厳格な割付けとエネルギー降下のゲート、は TS 探索の標準的な実務です。
  - QRC は IRC の安価な代用として確立しており、pyQRC のベンチマーク（TS 544 件）では IRC との一致率が 98〜99% です。
  - 振幅は二次近似で正確に決めています。変位方向は、質量加重した固有ベクトルを Cartesian に戻したものです。
  - ゲートを単調性ではなくエネルギー差で判定し、prominence の下限を SCF ノイズに連動させ、重なりを合否に使わない（CH-34）点も妥当です。
  - NH3 の反転は、正しく degenerate と判定されました。
- **懸念は、どれも実計算では発動していない分岐にあります。**
  1. `mode_index`（質量加重・射影後の順）と NWChem の `moddir`（driver 座標の Hessian の順）は、番号の意味が違います。
  2. `periodic_nearest` には許容幅がありません。
  3. torsional の既定のせいで、錯体の組み換えにも −20 cm⁻¹ の閾値が使われます（hypotheses の担当）。
- **QRC 固有の限界**（平坦な PES での近道、post-TS の分岐）は、接続を確かめるという目的の範囲内では許容できます。

#### 2.6.3 有用性

- **不可欠**: REFINE_SADDLE と xTB の初期 Hessian（3〜5 ステップ）、VALIDATE_TS（5〜15 s）、CONNECT/QRC（接続を示す唯一の証拠）。
- **費用の大半は QRC 側の DFT opt です**（stage の 51 / 81 / 65%）。原因は、対角推定の Hessian と trust 0.1 です。HONO の遅い側は、100 ステップ中 99 ステップが tiny eigenvalue による steepest descent でした。NH3 は 50 ステップ中 40 ステップが、負のモードの刻みの制限にかかっていました。
- **低頻度の分岐**（higher_order、moddir、collapse の分割、QRC の再試行）はコード量が小さく、化学的な意味もあるので残します。

#### 2.6.4 複雑さ

- 決定表は純関数で、予算も小さく、全体として簡潔です。
- 減らせるもの:
  - drop の式が 2 か所にある。
  - 再試行の倍率をクリップの後に掛けているので、上限 0.4 を超えて 0.8 Å になりうる。
- saddle の `trust 0.1` と `sadstp 0.1` は、NWChem の saddle の既定値と同じです（`opt_drv.F:1069-1074`）。冗長ですが、明示しておくのは妥当です。

#### 2.6.5 ライブラリ仕様の要点

- **NWChem driver 7.2.3**（https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/driver/opt_drv.F 、https://nwchemgit.github.io/Geometry-Optimization.html）
  - trust の既定は、最小化で 0.3、saddle で 0.1。sadstp は 0.1。
  - moddir≠0 のとき、1 ステップ目に driver 座標の Hessian の、昇順で moddir 番目の固有ベクトルを採る（:3270-3274）。以降は最大の重なり（ovtol 0.7）で追い、negeig==1 になれば firstneg でモード 1 を追う。
  - moddir=0 のときは、1 ステップ目を勾配の成分から上り方向に取る（"Initial step taken uphill"）。
  - 追っている moddir は RTDB に保存されるが、hfauto の継続は `start job` で始めるので、入力の moddir が再び適用される（:3374-3377、:927-960）。
  - 最小化では、負のモードの刻みは [0.03, 0.3]×trust に制限される。|固有値| < 1e-8 のモードは steepest descent になる（:1745-1790）。
  - INHESS 2 は、変位させた構造でも使える（プローブの 5 本で確認）。
- **pyQRC / QRC**（https://github.com/patonlab/pyQRC、Goodman & Silva, Tetrahedron Lett. 44, 8233 (2003)）: 既定の振幅は 0.2（正規モードの倍率）。振幅 0.3 で IRC との一致率は 98.9%。
- **IRC の理論**: Fukui 1970/1981、Gonzalez–Schlegel 1989/1990、Hratchian–Schlegel 2004/2005。post-TS の分岐は、Maeda ら IJQC 115, 258 (2015)、Ess ら ACIE 47, 7592 (2008)。

#### 2.6.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U6-I1 | medium | QRC 側の DFT opt が、対角推定の Hessian と trust 0.1 のせいで遅く、この単位の時間の大半を占めます（maxiter に達するリスクもあります）。検証による補足: trust 0.3 で対角推定のときの登り返しは、側によって起きたり起きなかったりします。「trust 0.1 が必要なのは Hessian がない場合だけ」は、プローブ 5 本から言える範囲の推論です。 |
| U6-I2 | low（検証で medium から下げた） | `mode_index` と `moddir` の番号が対応していません。継続すると、入力の moddir が再び適用されます。ただし、次のステップ以降は firstneg と最大重なりの追跡に戻り、誤ったモードを追っても TS 検証と QRC を通るので、誤った値が黙って順位に入ることはありません。 |
| U6-I3 | medium | ねじれ反応の `periodic_nearest` の割付けに許容幅がなく、第三の極小に落ちた側を端点に割り付けることがあります。補足: エネルギーだけの判定では、鏡像の回転異性体は区別できません。 |
| U6-I4 | low | QRC の 2 回目の振幅が上限 0.4 を超えて 0.8 Å になりうること。drop の式が 2 か所にあること（pyQRC の振幅は倍率なので、そのまま比べることはできない）。 |
| U6-I5 | low | CONNECT が、ゲートの判定より前に `_register`（opt と freq）を実行します。登録されるのは本物の極小で、無駄が出るのは失敗の分岐だけです。 |
| U6-I6 | low | torsional の既定のせいで、錯体の組み換えにも −20 cm⁻¹ の閾値が使われます。hypotheses の担当なので、この単位では変更しません。 |
| U6-I7 | low | 行 16 は、saddle の試行を使い切った後にも DFT string を実行します。無駄になるのは single_max の場合だけです。 |

#### 2.6.7 改良案

**U6-P1 QRC の両側の opt に、TS freq の Hessian を渡す（修正版）** — should、中立
- 変更:
  1. 幾何の一致の条件を緩めるのは optimize の経路だけにする。「元素の並びが同じで、同じフレームでの最大原子変位 ≤0.5 Å」とする。saddle の `_hessian_file` は完全一致の検査を残す。あわせて `protocols.py:135` の契約（「final geometry is mol」）を書き換える。
  2. `connect` で `rt.qm.optimize(ctx.mol(y), rt.method, init_hessian=freq, …)` を呼ぶ。
  3. 最初は trust を 0.1 のままにする。trust を 0.3 に戻すのは、W7 の known_endpoints の 3 反応と ab_init_hessian を再実行し、trajectory_above_ts が出ないことを確かめてからにする（minima の init_hessian の opt にも影響するため）。
  - 補足: `drivers/minimum.py:133` の mode-follow の opt にも、同じ関数で追加コストなしに効く（必須ではない）。
- 効果: trust 0.1 のまま Hessian を渡すだけで、HONO は 126 → 38 ステップになり、maxiter に達する継続がなくなる。NH3 は 50 → 45、HCN は 24 → 20。trust 0.3 まで入れると、HONO 12、NH3 24、HCN 12 / 15。paths stage は、HONO が 386 s → 約 120 s、NH3 が 108 s → 約 75 s になると推定される。
- コスト: 約 15 行とテストの更新、実計算での再確認に約 12 分。
- 優先度: 化学的な誤りではなく効率の改善なので should。

**U6-P2 負のモードが 2 本以上のとき、moddir の番号を driver 座標で決める（修正版）** — could、中立
- 変更: 並進と回転を射影した、質量加重なしの Cartesian Hessian P·H·P の負の固有ベクトルの中から、direction との重なりが最大のものの番号を `moddir` にする（約 10 行、`opt_drv.F:3270-3274` と一致する）。1 本しかなければ従来どおり 0 にする。
- 継続時に moddir を外す変更は見送る（効果は 1 ステップ程度）。

**U6-P3 periodic_nearest の割付けに判定を付ける（修正版）** — could、中立
- 変更: 既存の 5e-5 Eh を使った |ΔE| の判定と、「最も近い端点が一意であること（2 番目に近い端点との差が十分あること）」の両方を満たすときだけ採用する。満たさなければ `_register` に回す。新しい knob は作らない。

**U6-P4 振幅は倍率を掛けてからクリップし、drop の式を 1 か所にまとめる** — could、【複雑さ減】
- 変更: `min(qrc_amplitude(...) * factor**(attempt-1), bounds_A[1])` にする。drop は gates の関数で 1 回だけ計算する。
- 効果: 2 回目の振幅が上限 0.4 Å に収まる（1 回目が上限だった場合は同じジョブになり、キャッシュで無償に済む）。

#### 2.6.8 見送った案

- QRC を NWChem の MEPGS（IRC）に置き換える — P1 を入れれば QRC も 24 ステップで済み、差がなくなる。IRC の端点は結局 opt と割付けが必要になる。
- QRC の振幅を大きくする — ベンチマークでの差は 1% 程度。速度の問題は P1 で解決する。
- saddle と QRC を既定で tight にする — 既定の閾値は同一性の判定と両立している。
- TS で SCF の安定性やジラジカル性を検査する — NWChem に簡便な手段がなく、開殻は <S²> で止まる。
- 行 16 を変える — multi_max や monotonic のときに役立つ。
- 他の負のモードを正に反転した Hessian を渡す — 収束が遅くなり、検証もない。
- 常に moddir 1 を指定する — 余分なのは 1 ステップ程度で、利点がはっきりしない。
- _register をゲートの後にする — 失敗の分岐でしか効かない。
- reassigned の I→P を子反応にする — 新しい分類経路が必要になる。explore を再実行して扱う方が目的に沿う。
- QRC の回数や倍率を設定にする — 上書きしている例がない。
- ±2 本の QRC を並列に実行する — 4 vCPU を 4 rank がすでに使い切っている。
- post-TS の分岐を MD で調べる — 目的外。
- 重なりを hard gate に戻す — CH-34 で、正しい TS を棄却した実例がある。
- （U6-P1 の元案）最初から trust 0.3 に戻す — minima にも影響するので、再確認してから入れる。
- （U6-P2 の元案）継続時に moddir を外す — 効果が 1 ステップ程度で、分岐が増える。

#### 2.6.9 変更不要な点

- 別ジョブの DFT freq を自前で射影して TS を検証すること、imaginary_tier の 3 段階、prominence の下限、same_pes(numerics)。
- 初回に xTB Hessian を使い、DFT Hessian に切り替えるフォールバック。saddle の driver 設定。
- 曲率から QRC の振幅を決め、0.05〜0.4 Å に収めること。エネルギー降下による接続ゲート（G09/G11/G12）。
- Registry.find の厳格な同一性と、mapped_equivalent による縮退の判定。
- Hessian がないときの trust 0.1、重なりを合否に使わないこと。
- 決定表を純関数にすること、小さな予算、中間体での分割、higher_order のやり直し、開殻の <S²> の注記と blocker。

### 2.7 U7 熱化学（thermo: GoodVibes API・スケール因子・標準状態・会合量・感度幅・分布）

#### 2.7.1 現状仕様

**流れ**（`stages/thermochemistry.py`、`chemistry/thermo.py`、`backends/goodvibes/{engine,worker}.py`）

1. **対象**: tier=dft の全極小と、SaddleClaim の freq。`energy_method` があれば、構造指紋が一致する SP をエネルギー層に使う。
2. **温度と variant**: 温度は conditions の値で、ThermoSettings.temperatures_K を上書きする。variant は main（grimme / 100）に、qs{grimme, truhlar} × cutoff{50, 100, 150} を合わせた 6 つ。
3. **並列**: subject ごとに `thread_map` で 4 並列にする。
4. **反転**: `soft_imaginary_mode` の注記がある極小にだけ `invert_soft_cm1` を残す。既定は None。
5. **スケール因子**: 次の順に探す。
   1. 完全一致
   2. 同じ汎関数（PBE0 は PBE0/MG3S の 0.989 / 0.975。注記 `scale_factor_from`）
   3. 見つからなければ 1.0 / 1.0（注記 `scale_factor_unverified`）
   - 別名の表で `wb97xd3` を `wb97xd` に写す。
6. **worker**: subject 1 件で 1 プロセス。QCData に次を詰め、settings × T ごとに `compute_thermo` を呼ぶ。
   - 自前の射影振動数（負は `im_frequency_wn`）、主同位体の質量、回転定数、直線性、symmno=1。
   - 存在しないファイル名を入れて、GoodVibes に再パースさせない。
7. **整合ゲート**: main の結果で `thermo_consistent`（n_real、ZPE ≤1e-5、E ≤1e-6）を判定する。不合格なら G=None とし、thermo_unavailable にする。
8. **composite**: G = E_layer + (G_GV − E_GV)。組成内の Boltzmann 重みは main の表にだけ付ける。
9. **反応量**: dG_act、dG_rxn（degenerate なら 0）、dzpe_act、dE_act、dE_rxn。反応はすべて係数 1 の単一組成なので、標準状態による shift は常に 0 になる。band は 6 variant の dG_act の min/max。
10. **会合量**: 錯体（反応物の極小）と、単量体の ensemble_G（キーは (Hill 式, 電荷)）から求める。LOT が same_pes で一致するときだけ出す。
11. **blockers**: thermo_unavailable、mixed_level_of_theory（freq 層と energy 層の両方）、spin_contaminated。
12. **report**: rankable → ΔG‡ の順 → band が重なれば同順位。

**エンジン呼び出し（W7 HCN の実物）**

- argv: `python -m hfauto.execution.worker hfauto.backends.goodvibes.worker:compute <attempt>/job.json`（0.225 s）。
- job.json: symbols、coords、energy_hartree、`frequencies_cm1 [751.51, 751.51, 2243.82, 3462.34]`、n_external 5、settings 6 件（例: `{qs: grimme, cutoff_cm1: 100, vib_scale: 0.989, zpe_scale: 0.975, symmetry: true, invert_soft_cm1: null}`）。
- 呼び出し: `compute_thermo(qcdata=QCData(...), QS='grimme', s_freq_cutoff=100.0, temperature=298.15, freq_scale_factor=0.989, zpe_scale_factor=0.975, invert=None, symm=True)`。
- 渡さないもの（GoodVibes の既定値のまま）: QH=False、concentration=None（1 atm）、inertia 'global'（Bav 1e-44）、solv / spc は None。
- pymsym 0.3.5（symm から呼ばれる）: 点群と σ を求める。どの出力にも記録されない。

**主な既定値**: qs grimme、cutoff 100、vib_scale / zpe_scale None（表から探す）、symmetry True、invert_soft_cm1 None、sensitivity True、temperatures [298.15]、standard_states [1atm]（1M は 298.15 K で Δn あたり 1.894）、ZPE 許容 1e-5 Eh、E 許容 1e-6 Eh、rankable の ΔZPE‡ 許容 1 + 4·max(1, n_H)。

**利用者が設定する項目**: thermo の engine・settings（上記）・energy_method、conditions、gates、site の goodvibes（版数 4.3.0、python）、report の T_K / standard_state / metric、system の compositions。

**実測**

- stage は 0〜1 s、worker は 0.22〜0.32 s/subject、失敗 0。
- HCN: ΔG‡ 42.218、band の幅 0.0006。全 subject に `scale_factor_from:PBE0/MG3S` が付いた。final_check でも同じ値。
- HONO: 12.243（旧実装の 24.11 は再現しない）。NH3: 3.850。symm の寄与は +0.411。
- TMA·(HF)₂: ΔG_assoc −11.594。
  - 錯体の最低振動数は 51.8、84.0、121.3 cm⁻¹。
  - 6 variant で見ると ΔG_assoc は −12.16〜−11.01 だが、この幅は出力されない。
- 整合ゲートの余裕: 12 subject すべてで ZPE の差 ≤1.6e-8 Eh、E の差 0。
- プローブ:
  - QH=True: ΔG‡ ≤0.001、ΔG_assoc −0.675。
  - symm=False: NH3 の ΔG‡ −0.411、ΔG_assoc +0.651。
  - スケールを 1.0/1.0: ΔG‡ −0.01〜−0.04、ΔG_assoc +0.065。zpe_scale は G に影響しない。
  - 273 / 373 K の ΔG_assoc: −13.15 / −6.88。
  - 1 M の ΔG_assoc: 約 −15.38。
  - 単原子（F⁻）: worker が ZeroDivisionError を出す（上流の Evidence が n_external ∈ {5, 6} に限るので、実際には到達しない）。

#### 2.7.2 化学的妥当性

- **中核は正しいです**（WSL にインストール済みの `goodvibes/thermo.py`・`api.py` を読んで確認）。
  - Grimme qRRHO のエントロピー、Sackur–Tetrode の 1 atm 基準、S_elec = R ln(2S+1)、RRHO の H、pymsym の σ（二重計上なし）、標準状態の換算式、composite の式。
  - 自前の振動数を渡すので、パーサ依存の問題（CH-39 の負の S_rot）も排除できています。
- **妥当性の問題**
  1. **負の振動数を捨てること。** noise・soft の極小と、TS の副次モードの負の振動数は、GoodVibes では捨てられます。プローブ（TMA·(HF)₂ の最低モードを置き換え）では、同じ |ν| でも符号だけで G が変わりました（±20 cm⁻¹ で 0.74、±5 で 1.15、±1 で 1.62）。xtb（imagthr −20）も GoodVibes の `invert='auto'` も、反転する扱いです。
  2. **会合量の LOT 判定が電荷と多重度まで要求すること。** このため、陰イオン + 中性分子や、ラジカルを含む会合では、ΔG_assoc と ΔG‡ vs separated が常に None になり、理由も出ません。README は、電荷を持つ系と開殻系を範囲に含めています。
  3. **経路の縮重度**: σ_R/σ‡ しか入らず、キラリティ（HONO の C1 の TS で L=2）と対称反応の因子 2（NH3 の D3h）は入りません。ただし ΔG‡ を種どうしの G の差と定義するなら、σ を入れた値は熱力学量として正しいです（検証の訂正）。規約として明記すれば足ります。
- **数値への影響が小さい点**
  - zpe_scale（0.975）は G と H に入りません（GoodVibes は U_vib を harm 0.989 で計算する）。
  - スケール因子の借用の効果は ≤0.07。
  - Truhlar QS の実装は原論文と違いますが、感度の variant でしか使いません。
  - トンネル補正はありません（障壁の指標としては許容範囲）。
  - pymsym の既定の閾値は厳しいです（0.001 Å 乱すと TMA は C1 になる）。ただし実際の最適化構造はすべて正しく判定されました。

#### 2.7.3 有用性

- **不可欠**: main の GoodVibes の計算、σ の補正（ΔG_assoc で 0.65）、会合量と ΔG‡ vs separated（錯体を対象にするうえでの中核）、ΔZPE‡ の検査、thermo_consistent、fail-closed。計算コストは無視できます。
- **効果が小さい**:
  1. ΔG‡ の band（幅 <0.001。幅 1.16 の ΔG_assoc には付かない）
  2. スケール因子の 3 段検索（効果 ≤0.07）
  3. 1 bar / 1 M（効くのは会合量だけ）
  4. population・H・S_rot（誰も読まない）
  5. composite（使われていない）
  6. temperatures_K（効かない）
  7. invert_soft_cm1（既定 None）
- 数値を最も大きく左右する「負モードを捨てる扱い」に手当てがない一方で、効果の小さい機能に仕組みが割かれています。

#### 2.7.4 複雑さ

- 削れる、または統合できるもの:
  - (a) スケール因子の DB 検索（`_SCALE_DB`・`_ALIASES`・`resolve_scales`）と zpe_scale を、明示値 1 本にまとめる。
  - (b) `invert_soft_cm1` と、soft だけに効かせる分岐を、固定の規則に置き換える。
  - (c) ThermoSettings.temperatures_K を、利用者が設定するキーから外す。
- これらで約 50 行と設定キー 3 個が減り、負モードの偏りという正しさの問題も同時に直ります。
- variant の格子、標準状態の 3 種、thermo_consistent、純関数への分割は、そのままでよいです。

#### 2.7.5 ライブラリ仕様の要点

- **GoodVibes 4.3.0**（WSL の prod venv にある `goodvibes/thermo.py`、`api.py`、`scaling_factors.json`）
  - `calc_bbe`: enthalpy = E + U_trans + U_rot + U_vib + RT。U_vib は ZPE を含み、freq_scale_factor で計算する。zpe_scale_factor は報告される zpe にしか効かない。
  - `im_frequency_wn` は、反転されない限り ZPE・U・S のどれにも入らない。
  - `_apply_frequency_inversion`: 数値を渡すと、x > invert の負モードを反転する。'auto' なら全負モードを反転し、TSFreq のときだけ最も負のモードを残す。
  - Truhlar QS はエントロピーだけに効く。Grimme の減衰関数はスケール前の ν で計算する。QH=False なら qh_G = H(RRHO) − T·qh_S。
  - `concentration=None` なら 1 atm 基準。`sym_correction` は pymsym を使い、点群を `bbe.point_group` に入れるが、`ThermoResult.point_group` は空になる。
  - スケール因子の DB は Truhlar DB v5 の 237 件。PBE0/def2 の項目はない。def2-SVP でも、汎関数によって約 ±1% ばらつく。
- **pymsym / libmsym 0.3.5**: 既定の閾値は 1e-3。例外が出ると 'C1' を返す。閾値の引数はない。
- **xtb**（https://xtb-docs.readthedocs.io/en/latest/xtb_thermo.html 、https://www.mankier.com/7/xcontrol）: imagthr は −20 cm⁻¹（これより上の虚振動を反転する）、sthr は 50 cm⁻¹。
- **経路の縮重度**: Fernández-Ramos, Ellingson, Meana-Pañeda, Marques, Truhlar, Theor. Chem. Acc. 2007, 118, 813（https://link.springer.com/content/pdf/10.1007/s00214-007-0328-0.pdf）。L = σ_R m‡ / (σ‡ m_R)。
- **qRRHO**: Grimme, Chem. Eur. J. 2012, 18, 9955。Luchini ら F1000Research 2020, 9, 291（GoodVibes）。Alecu, Zheng, Zhao, Truhlar, JCTC 2010, 6, 2872（スケール因子）。Skodje & Truhlar, J. Phys. Chem. 1981, 85, 624（T_c）。Kozuch & Shaik, Acc. Chem. Res. 2011, 44, 101（energetic span）。

#### 2.7.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U7-I1 | medium（検証で high から下げた） | 負の振動数（noise・soft の極小、TS の副次モード）を熱化学から捨てます。1 モードあたり 0.5〜1.6 のずれで、しかも符号によって不連続です。検証による訂正: 「注記なし」は言い過ぎで、極小と TS の −50〜−10 には *_imaginary_mode / soft_secondary_mode の注記が付きます。注記がないのは、TS の −10〜0 の副次モードだけです。どの注記も「熱化学から除外した」ことまでは示しません。W7 では 0 件なので、影響は潜在的です。 |
| U7-I2 | medium | 会合量の LOT 判定が、same_pes の電荷と多重度まで要求します。陰イオンやラジカルを含む会合では、ΔG_assoc が常に None で、理由も出ません。G=None の単量体は黙って除外されます。 |
| U7-I3 | low | 経路の縮重度（キラリティ、対称反応の因子 2）が入っていません。検証による訂正: σ を入れた 3.85 は熱力学量として正しく、validation の「+0.411」は σ が適用されたことの確認としては誤りではありません。速度に対応する有効障壁（L=1）は 3.44 です。これは規約の問題です。 |
| U7-I4 | low | zpe_scale が G と H に入りません（GoodVibes 側の意図した仕様）。スケール因子は借用値で、別名の表が ωB97X-D3 と ωB97X-D を同一視しています。効果は ≤0.07。 |
| U7-I5 | low | 感度の幅が ΔG‡ にだけ付き、ΔG_assoc（幅 1.155）には付きません。 |
| U7-I6 | low | QH=False（RRHO の H と qRRHO の S）が設計文書に書かれていません。 |
| U7-I7 | low | pymsym の点群と σ がどこにも記録されません。 |
| U7-I8 | low | 順位の基準が反応物の極小そのものなので、仮説の向き（conformer は低い側から、discovery は出発極小から）に依存します。検証による訂正: ずれの上限は選択窓の 6 kcal/mol ではありません（状態ラベルが違えば、6 を超えうる）。また、出発極小を基準にした障壁にも意味はあります。 |
| U7-I9 | low | 効かない設定（temperatures_K）と、書くだけの出力（population、H、S_rot）があります。 |
| U7-I10 | low | トンネル補正がありません。T_c は HCN 258 K、amine の PT 290〜348 K。速度論は範囲外なので、規約として明記すれば足ります。 |
| U7-I11 | low | 単原子種を扱えません。検証による訂正: 実質的な制約は上流の `Evidence.n_external ∈ {5, 6}` で、worker の ZeroDivisionError には到達しません。 |

#### 2.7.7 改良案

**U7-P1 負モードを固定の規則で扱う（修正版）** — must、【複雑さ減】（= U3-IMP-4。§3 の X3）
- 変更:
  - chemistry に純関数 `thermo_frequencies(freqs, saddle)` を 1 つ作る。規則は「極小は全負モードを |ν| にする。TS は最低モードだけを除き、残りの負モードは |ν| にする」。
  - engine が job.json の `frequencies_cm1` を作るときと、`thermo_consistent` の参照 ZPE・n_real の両方で、この関数を使う。GoodVibes には常に `invert=None` を渡す。
  - `invert_soft_cm1` と、soft 限定の分岐（`thermochemistry.py:91-92`）を削除する。
  - GoodVibes の `job_type='TSFreq'` と `invert='auto'` は使わない。挙動が文字列の一致と浮動小数の等値比較に依存し、Windows のユニットテストで検証できないため。数値の閾値 `invert=-50` も使わない。ねじれの TS で反応座標まで反転してしまうため。
  - 対象には、注記のない TS の −10〜0 の副次モードも含める（検証の追加）。
- キャッシュの注意: thermo の JobStore キー `{freq job_key, hessian_sha, settings}` は振動数リストを含まない。規則を変えたら、settings のスキーマの変更か規則のタグで、キーを必ず変えること。
- 効果: 1 モードあたり +0.48〜+1.6 の偏りと、符号による不連続がなくなる。W7 の数値は変わらない。設定キーが 1 つ減る。
- コスト: 20〜30 行とテスト 3 件。

**U7-P2 会合量の LOT 判定から電荷と多重度を外す（修正版）** — should、中立
- 変更:
  - 錯体と単量体をまたぐ LOT 判定では、same_pes から charge と multiplicity だけを外して比べる（引数を 1 つ足す）。反応の blockers は現行の same_pes のままにする。
  - 欠けた単量体は、fail-closed で (None, None) にする。
  - 追加: 単量体のキーを (Hill 式, 電荷, 多重度) にする。そうしないと、スピン状態の違う極小が同じアンサンブルに混ざる（U1-08 も解消する）。
- 効果: 陰イオン（多原子）と中性分子、ラジカルと閉殻分子の会合量が出るようになる。中性の amine·HF の値は変わらない。

**U7-P3 スケール因子を明示値 1 本（既定 1.0）にまとめる** — should、【複雑さ減】
- 変更: `_SCALE_DB`・`_ALIASES`・`_canonical`・`scale_factors`・`resolve_scales` と zpe_scale を削除する。`vib_scale: float = 1.0` を、freq_scale_factor と zpe_scale_factor の両方に明示して渡す（GoodVibes の自動検索は使わない）。
- 効果: 約 35 行、設定キー 1 個、テスト 4 件が減る。ZPE と G で使うスケールが揃う。数値の変化は ΔG‡ で ≤0.04、ΔG_assoc で ≤0.07（validation の値は小数第 2 位が変わるので、記録を更新する）。

**U7-P4 熱化学の規約を文書に書く（文書のみ、修正版）** — should、中立
- design.md の thermo の行に、次の 1〜2 文を足す。
  - QH=False（RRHO の H と Grimme qRRHO の S）。G に入るスケール因子は 1 つ。
  - σ は外部回転の対称数だけで、経路縮重度 L（キラリティ m、対称反応の因子）は含まない。
  - トンネル補正はない。基準は 1 atm。
- validation.md の NH3 の行に、「+0.411 は種の G に入る σ の寄与である。速度に対応する有効障壁（L=1）では 3.44 になる」と注記する。「誤った解釈を正す」という書き方はしない。

**U7-P5 点群を SpeciesThermo の notes に記録する（修正版）** — could、中立
- 変更: worker で `r.bbe.point_group` を result に入れ、書くだけの `S_rot` を point_group に置き換える。notes に `point_group:<pg>` を 1 つ足す。σ は点群から決まるので記録しない。

**U7-P6 ΔG_assoc にだけ感度の幅を付ける（修正版）** — could、【複雑さ増】
- 変更: 既存の 6 variant の表で `_association` を回し、`assoc_band_kcal` を 1 フィールド表示する。順位には使わない。vs_separated の幅と、表の列の追加は省く。

**U7-P7 組成アンサンブル基準の有効障壁を表示する（修正版）** — could、【複雑さ増】
- 変更: main の表から `dG_act_eff_kcal = G_TS − ensemble_G`（`_populated` と同じ群）を 1 フィールド出し、ranking.csv と html に表示する。順位の既定にするかは、別の設計判断とする。アンサンブルは DFT の極小選択窓の網羅性に依存する、という注記を付ける。

**U7-P8 ThermoSettings.temperatures_K を削除する** — could、【複雑さ減】
- 変更: 温度は conditions だけから決め、engine.thermo の引数として渡す。キャッシュキーには、温度を別のフィールドとして入れる。

#### 2.7.8 見送った案

- トンネル補正（Wigner / Eckart / SCT） — 順位は障壁の指標で決める。T < T_c では Wigner が成り立たず、Eckart と SCT は複雑になる。
- 経路の縮重度 L を自動で計算する — 効果は 0.41 で、判定を一般化しにくい。
- QH=True に変える — どちらが正しいとも言えず、W7 との連続性が失われる。
- hindered rotor や非調和の補正 — 追加の QM 計算が必要で、目的に対して過剰。
- BSSE の補正を thermo に入れる — SP の責務で、ジョブが増える。
- 感度の variant の格子を削除する — 安価で、柔らかい系では意味がある（CH-46）。
- 1 bar / 1 M を削除する — 約 6 行で、1 M は溶液の文献値と比べるのに必要。
- thermo_consistent を削除する — CH-38 の回帰を防ぐ唯一の番人。
- GoodVibes をやめて自前で書く — 再検証が必要で、利得がない。
- variant ごとに整合ゲートを判定する — 同じ入力と処理なので不要。
- worker を単原子種に対応させる — 上流が制約していて、現状は対象がない。
- pymsym の閾値を緩める — σ を過大に判定する危険がある。記録で監査する（P5）。
- 開殻の軌道縮退やスピン軌道 — 効果は ≤R ln 2 で、該当する種は稀。
- Truhlar QS を原論文どおりに直す — 感度でしか使わず、保守が増える。
- （U7-P1 の元案）GoodVibes の `invert='auto'` と `job_type='TSFreq'` を使う — 文字列と浮動小数の等値比較に依存し、Windows で検証できない。
- （U7-P6 の元案）vs_separated の幅と表の列も加える — 省いてよい。
- （U7-P7 の元案）report の metric に 'dG_act_eff' を足す — 過剰。表示だけにする。
- （U7-P5 の元案）σ も notes に記録する — 点群から決まる。

#### 2.7.9 変更不要な点

- GoodVibes の API を worker サブプロセスで呼び、自前の振動数・質量・回転定数を渡し、再パースを避け、版数を pin し、スケール因子を明示して渡すこと。
- Grimme qRRHO（100 cm⁻¹、Bav 1e-44）、1 atm 基準と standard_state_shift、symm=True で QCData の symmno=1。
- thermo_consistent と fail-closed。composite の式と、両方の層での mixed_level_of_theory。
- 会合量の定義と ΔG‡ vs separated、ΔZPE‡ の許容、縮退で ΔG_rxn=0、TS を含める outcome の集合、純関数への分割、S_rot の事前検査。

### 2.8 U8 一点計算・手法パネル・順位付けと報告・実行基盤（sp・report・site / execution / config）

#### 2.8.1 現状仕様

**流れ**

1. **起動**: `hfauto run method_panel --system <S> --site <SITE> --run-dir <既存 run>`。`--run-dir` を省くと `runs/<system_id>_<pipeline_id>` が新しく作られ、既存の run とは別になる。
2. **preflight**: エンジン名、実行ファイル、版数の pin（NWChem は実行後に出力で照合）、worker のモジュール、scratch が `/mnt/<drive>` 上にないこと、`engine.supports(method)` を検査する。
3. **run_pipeline**
   - `<scratch_root>/.hfauto_site.lock` を run の間ずっと保持する。
   - stage id が他の pipeline と衝突したら、ValueError。
   - input_sha と config_sha が変わっていない done の stage は飛ばす。
   - 入力の view は、それより前に done になった全 stage の manifest の union（同じ id は後の方が勝つ）。
4. **SinglePointStage**
   - 対象は tier=dft の極小だけ。`reaction_stationary_points`（反応の極小と TS の freq）か、`all_minima`。
   - 構造は freq Evidence.final、電荷と多重度は freq.level から取る。
   - method × 対象を直列に、deadline なしで `engine.energy` に渡す。同じ構造は 1 つの artifact にまとめ、parents をつなぐ。
5. **NWChemEngine.energy**
   - CCSD(T) で多重度 >1 なら、INPUT_INVALID `ccsd(t)_closed_shell_only`。
   - Level を観測する。WFT では、`Total CCSD(T) energy` / `Total MP2 energy` と、SCF 出力の `open shells =` 行から読む。
   - 入力構造の echo が 1e-4 Å 以内であることを確かめる。
6. **JobStore**
   - キーは正規化した JSON `{engine, version_pin, kind, key_payload}` の sha256。ExecutionSpec は含めない。
   - `.lock` は O_EXCL で取る。attempt_NN の番号を振る。
   - 終端的な失敗（INPUT_INVALID、METHOD_MISMATCH、EXECUTABLE_MISSING、INCOMPLETE_OUTPUT）は記録して再利用する。それ以外の失敗は、次回に再実行する。
   - ladder: TIMEOUT と MAXITER は driver ジョブだけ 2 回、autoz は 1 回、SCF は 1 回。energy ジョブの timeout は継続しない。
7. **report**
   - `rank_rows`: (T, 標準状態) の ReactionThermo を選び、rankable でふるい、ΔG‡ の順に並べ、band が重なれば同順位（1,1,3 方式）にする。
   - `coverage`: 機構別の試行・生成物・陰性理由と、FailureKind 別の件数。
   - `method_panel`: Level.full_key ごとに、構造指紋 → エネルギーの表から ΔE_rxn と ΔE‡ を出し、min/max と sign_disagreement（許容幅 0）を付ける。参照の opt / saddle / freq も 1 行として入る。
   - 出力: ranking.csv・coverage.csv・method_panel.csv と report.html（3Dmol.js 2.5.5 を CDN から読む）。
8. **補助 CLI**: `hfauto report RUN_DIR` は `<run>/report.html` だけを作り直す。`hfauto status` は状態、失敗の種類、再利用率を表示する。

**エンジン呼び出し（実物）**

- DFT SP（本レビューのプローブ hcn_panel、`jobs/49/49f32be6…/job.nw`）: `mpirun -np 4 nwchem job.nw`、OMP_NUM_THREADS=1。
  ```
  memory total 1200 mb          # heap 300 / stack 300 / global 600 MB/rank
  geometry units angstrom nocenter noautosym … end
  charge 0
  basis spherical ; * library def2-tzvp ; end
  dft ; xc pbe0 ; mult 1 ; grid fine ; convergence energy 1.0e-07 ; disp vdw 4 ; end
  task dft energy
  ```
- MP2（W7 の smoke、HCN/def2-SVP）: `mp2 ; freeze atomic ; end` と `task mp2 energy`。Frozen core 2、SCF 閾値 1e-6、1.5 s。
- CCSD(T)（プローブ、HCN の TS）: `ccsd ; freeze atomic ; end` と `task ccsd(t) energy`。maxiter と thresh は書かないので、NWChem の既定（20、1e-6）になる。
  - HCN: 68 関数、CCSD は 10〜13 反復。
  - HONO: 99 関数、18 反復。
  - T1 / D1 は出力されない。
- 開殻 WFT: `scf ; nopen <m−1> ; uhf ; end` を書く（描画のみ。実際には後述のとおり常に失敗する）。
- SCF の継続: DFT は `vectors input job.movecs` と `convergence damp 40 ncydp 30 lshift 0.5`。WFT は `vectors input` だけ。

**主な既定値**

| 項目 | 値 |
|---|---|
| sp の targets | reaction_stationary_points |
| panel の methods | [pbe0-d3bj_def2-svp, pbe0-d3bj_def2-tzvp, wb97x-d3_def2-tzvpd, ccsd-t_def2-tzvp] |
| memory_mb_per_rank | 1200 |
| timeout | 14,400 s |
| 開始に必要な最低予算 | 60 s |
| SiteLock | scratch_root 単位 |
| コア数セマフォ | min(ranks × threads, cores) |
| report | T は conditions の先頭、標準状態は 1atm、metric は dG_act |
| 符号ゲートの許容幅 | 0 |

**利用者が設定する項目**: CLI（`--run-dir`、`--from`、`--to`、`--dry-run`、`--retry-failed`）、sp（engine、methods、targets）、report（T_K、standard_state、metric）、thermo の energy_method、conditions、gates、method ファイル、site（scratch_root、cores、memory_mb、engines.* の execution）。

**実測**

- W7 と final_check では、sp と method_panel は一度も実行されていない。report は 1 s 未満で、ジョブは 0 件。
- ranking（W7）: HCN rank 1（ΔG‡ 42.2178、band の幅 0.0006）、HONO rank 1（12.2430）、NH3 rank 1（3.8505）。water と tma_hf2 は same_basin で、blockers に `outcome:same_basin;thermo_unavailable;dzpe_out_of_tolerance` が付いた。
- coverage（tma_hf2）: attempts は afir 17 / nt2 20。陰性の理由は monotonic_uphill 12、same_as_source 9、afir_not_converged 7 など。
- パネルのプローブ:
  - HCN は 12 SP で 63.7 s。1 点あたり PBE0 1.6〜2.7 s、ωB97X-D3/TZVPD 6.5〜7.0 s、CCSD(T) 8.9〜10.5 s。
  - HONO は負荷が高い状態（load 10〜13）で 1,400 s。CCSD(T) は 1 点 298〜322 s。sign_disagreement で順位から外れた。
- W7 の jobs を HEAD で再生すると、29/29 件がヒットした（ジョブキーは W7 から HEAD まで変わっていない）。
- MP2/def2-TZVP のプローブ: TMA·HF（215 関数）が wall 287 s、TMA·(HF)₂（252 関数）が wall 295 s。amine·(HF)n の CCSD(T) は、o³v⁴ の外挿で 1 点数時間と推定した（実測ではない）。
- 開殻の MP2（OH）: NWChem は `Total MP2 energy` まで正常に計算したが、hfauto は `level_not_observed`（INCOMPLETE_OUTPUT、終端）で失敗した。

#### 2.8.2 化学的妥当性

- **sp は妥当です。** 構造を freq の最終構造に固定し、Level と構造の echo を照合する設計は、固定構造での一点計算として正しいです。パネルの化学的な妥当性も実測で確かめました。
  - HCN の CCSD(T)/def2-TZVP は ΔE_rxn 15.43 で、文献の ΔE0 14.83（van Mourik）とよく合う。
  - HONO の CCSD(T) は ΔE_rxn +0.22 で、実験値 +0.30（99±25 cm⁻¹）と符号が一致する。
- **WFT のデッキ**: `freeze atomic` と thresh 1e-6 は十分です。NWChem 7.2.3 のソースでは、CCSD が収束しないと (T) を飛ばして errquit で終わるので、未収束の値が黙って使われることはありません。ただし、maxiter 20 に対して HONO の TS は 18 反復で、余裕が小さいです。
- **開殻**: DFT は正しく動きます。CCSD(T) は起動前に拒否されます。UHF-MP2 は、UHF の SCF 出力に `open shells =` 行がないので常に失敗します（誤った値を黙って使うことはない）。閉殻の陰イオン（F⁻）は、MP2・CCSD(T)・PBE0 のどれも正常に動きます。
- **パネルの参照**: CCSD(T)/def2-TZVP には diffuse 関数がないので、陰イオン系の参照としては弱いです（U0 と同じ指摘）。
- **符号ゲートの許容幅 0 は過敏です。** HONO では、正しい符号の 3 レベルに対し、diffuse のない 2 行だけが負になり、検証済みの反応が外れました。
- **report**: rankable と、1,1,3 方式の同順位は論理的に正しいです。band は qRRHO の設定差だけから作るので、順位は「参照レベルでの ΔG‡ の順」と読むのが正しいです。
- **実行基盤**: JobStore、ladder、SiteLock、`/mnt` の拒否は、4 vCPU の WSL で安全に再開するという目的に合っています。

#### 2.8.3 有用性

- **必須で費用に見合うもの**: rankable と ranking.csv、coverage.csv（陰性結果を報告する唯一の出力）、JobStore・ladder・SiteLock・セマフォ、sp の本体とパネル（小分子なら 1 分程度）。
- **限界的なもの**:
  - (a) CCSD(T)/def2-TZVP。小分子では価値が最も高いが、amine·(HF)n では 1 点数時間以上かかり、14,400 s を超える。timeout は記録されないので、再実行のたびにやり直す。
  - (b) PBE0/def2-SVP。安いが情報量が小さく、誤った blocker の原因になった。
  - (c) composite G（使われていない）。
  - (d) report の knob。
  - (e) `hfauto report`（成果物が二重になる）。
  - (f) site の memory_mb（読まれていない）。

#### 2.8.4 複雑さ

- 全体として簡素です（sp 78 行、report 68 行、summary / html は純関数）。
- 削れる、または既定で済ませられるもの:
  - (1) 開殻 WFT の描画（実際には常に失敗する）。
  - (2) memory_mb（検査に使うか、削除する）。
  - (3) 既定のパネルに CCSD(T) が入っていること（opt-in にする）。
  - (4) composite のコメントが実態と合わない。
- conditions が 3 本の pipeline に重複して書かれていることは気になりますが、構造変更になるので扱いません。

#### 2.8.5 ライブラリ仕様の要点

- **NWChem の ccsd モジュール**（https://nwchemgit.github.io/CCSD.html）: 閉殻 RHF 専用。MAXITER 20、THRESH 1e-6、DIISBAS 5。T1 / D1 は出力しない。
- **NWChem の CCSD のソース**（https://raw.githubusercontent.com/nwchemgit/nwchem/v7.2.3-release/src/ccsd/ccsd_iterdrv2.F 、…/aoccsd2.F、…/task/task.F）: 収束しないと "maximum iterations exceeded" を出し、(T) を飛ばして errquit で終わる。hfauto ではこれが NONZERO_EXIT（非終端）になり、再実行のたびに同じ失敗を繰り返す。
- **NWChem の SCF / MP2**（https://nwchemgit.github.io/MP2.html）: UHF の出力は `alpha electrons` / `beta  electrons` で、`open shells` 行がない。既定の出力に ⟨S²⟩ はない。
- **NWChem の DFT**: `disp vdw 4` は D3(BJ)。`xc wb97x-d3` は分散を含む。`memory total` は rank ごとの値。
- **GMTKN55 / BH9**（https://doi.org/10.1039/C7CP04913G 、https://pubs.acs.org/doi/10.1021/acs.jpca.2c03922）: 障壁と反応エネルギーでは、範囲分離ハイブリッドがグローバルハイブリッドより良い。
- **HONO の cis/trans の実験値**（https://www.sciencedirect.com/science/article/abs/pii/S0022285209002574）: cis が 99±25 cm⁻¹（約 +0.30 kcal/mol）高い。

#### 2.8.6 課題

| ID | 重大度 | 内容（検証後） |
|---|---|---|
| U8-1 | medium（検証で high から下げた） | method_sign_disagreement の許容幅が 0 で、HONO trans→cis が外れます（U0-I1 と同じ事象）。検証による補足: 影響するのは panel_report の ranking だけで、元の report/ranking.csv は rank 1 のまま残ります。1 つの run に、矛盾する 2 つの ranking.csv ができます。 |
| U8-2 | medium | 既定のパネルに CCSD(T)/def2-TZVP が入っており、主な対象の amine·(HF)n では実行できません（外挿は推定）。energy ジョブの TIMEOUT は継続も記録もされないので、再実行のたびに同じ 4 h をやり直します。系の大きさで止める仕組みもありません。 |
| U8-3 | low（検証で medium から下げた） | CCSD の maxiter が既定の 20 のままです。未収束は NONZERO_EXIT（非終端）になり、再実行しても同じ失敗を繰り返します。IMP-2 を入れれば頻度は低くなります。 |
| U8-4 | low | 開殻 WFT（UHF-MP2）の経路が常に `level_not_observed` で失敗します。描画のコードは使われていません。 |
| U8-5 | low | composite G がどこにも配線されておらず、コメントが実態と合いません。thermo の artifact id は stage id を含まないので、2 つ目の thermo stage を同じ run に置くと、後勝ちで上書きされます。 |
| U8-6 | low | 順位の区間が qRRHO の幅だけで作られます。検証による訂正: これは AR-27 の意図どおりで欠陥ではなく、解釈上の注意です。 |
| U8-7 | low | ranking.csv に T、標準状態、ΔG_rxn、理論レベルの列がありません。 |
| U8-8 | low | TS のない outcome にも `dzpe_out_of_tolerance` が付きます（表示上の誤読の問題だけ）。 |
| U8-9 | low | site の memory_mb はどこからも読まれず、ranks × memory の上限検査もありません。 |
| U8-10 | low | panel_report は、method_panel.yaml 自身の conditions を使います。元の run と条件が違うと、全反応が thermo_unavailable になります（失敗は明示されます）。 |

#### 2.8.7 改良案

**U8-IMP-1 符号ゲートに ±1 kcal/mol の許容幅を入れる** — U0-P1 に統合（must。§3 の X4）。design.md の該当行も同時に直す。

**U8-IMP-2 既定のパネルから CCSD(T) を外し、小分子での opt-in にする** — should、【複雑さ減】
- 変更: `method_panel.yaml` の methods から ccsd-t を外す。コメントに「閉殻で重原子 4 個程度までなら、CCSD(T) を足す（HCN・HONO で 1 点数秒〜数分）」と書く。method ファイルとコードは残す。
- 効果: amine·(HF)n で、4 h × 点数 × 再実行の事故がなくなる。U0-P2（TZVPD 化）との最終形は §3 の X4 を参照。

**U8-IMP-3 CCSD ブロックに maxiter 50 を書く** — could、中立（= U0-P6）
- 変更: `render_wft` で、module が ccsd のときだけ `maxiter 50` を加える（定数 CCSD_MAXITER）。ジョブキーは変わらないので、既存の成功結果は再利用される。

**U8-IMP-4 WFT を閉殻専用と明示し、UHF の描画を削除する** — should、【複雑さ減】
- 変更: `engine.energy` の判定を `kind=='wft' and multiplicity>1` に一般化し、INPUT_INVALID `wft_closed_shell_only` にする。`input.py:163-164` の nopen / uhf を削除し、restart 用の scf ブロック（vectors input）は残す。パーサを直して UMP2 を通す案は、⟨S²⟩ を検査できない値がパネルに入るので採らない。
- 効果: 開殻では、原因が分かる理由ですぐに失敗する。使われないコードも消える。

**U8-IMP-5 composite 層の説明を実態に合わせる（文書のみ）** — could、【複雑さ減】
- 変更: `pbe0-d3bj_def2-tzvp.yaml` のコメントを「Method panel: TZ basis-set sensitivity (no diffuse functions)」に直す。design.md の thermo の行に、次を書く。
  - composite は、利用者が sp → thermo(energy_method) → report の pipeline を書いて使う。
  - 陰イオンや H 結合の系では、TZVPD の層を使う。
  - 同じ run に thermo を 2 つ置くと後勝ちになるので、別の run-dir で行う。

**U8-IMP-6 ranking.csv に T_K・standard_state・dG_rxn_kcal を加える（修正版）** — could、中立
- 変更: RankRow に 3 フィールドを足し、列を追加する。HTML に理論レベルのラベルを出す部分は見送る（method_panel.csv の参照行で代用できる）。

**U8-IMP-7 TS のない outcome では ΔZPE を検査しない** — could、中立
- 変更: `rankable` で、dzpe の判定を `dG_act_kcal is not None` のときだけ行う。dG_act があって dzpe が None の場合は、従来どおり発火させる。

**U8-IMP-8 preflight で memory を検査する（修正版）** — could、中立
- 変更: 新しい係数 0.8 は持ち込まない。NWChem 系のエンジンについて、`ranks × memory_mb_per_rank ≤ site.memory_mb` を検査する（memory_mb を「使ってよい量」と定義し、environment.md の 0.8 × MemTotal は memory_mb を決めるときの指針として残す）。検査を入れないなら、memory_mb を削除する。

**文書（検証の追加）**: README か design に「パネルを追記した後は、panel_report が最終の順位」と 1 行書く。コードの変更は不要。

#### 2.8.8 見送った案

- composite G を既定の流れか method_panel に組み込み、ωB97X-D3/TZVPD の層で順位を付ける — HONO では良くなるが、HCN では遠ざかる。一貫して良くなる根拠がなく、17 原子では 1 点数十分かかり、artifact が上書きされる。追記用の pipeline として別に用意する案は、U0-P3 で採った。
- method_sign_disagreement を完全に外す — 設計原則 5 がなくなる。許容幅で足りる。
- パネルの幅を同順位の区間に加える — ほぼすべてが同順位になる。
- energy ジョブの TIMEOUT を TERMINAL にする — 一時的な失敗なので定義に合わない。根本の原因は IMP-2 で解消する。
- CCSD(T) に系の大きさの上限を設ける — 閾値と見積もりのコードが必要になる。opt-in で十分。
- TCE による開殻 CCSD(T)、T1 / D1、CBS 外挿、F12 — 計算も実装も大きく増える。
- DLPNO / RI-MP2 で大きな系の参照を作る — ORCA は削除済みで、MP2 は参照として弱い。
- sp で CP 補正を行う — 異性化の障壁では打ち消し合う。会合量は thermo の担当。
- ホスト単位の SiteLock、プロセスをまたぐセマフォ — 過負荷は人工的な並列実行のときだけ起きた。
- conditions を run 単位で継承させる — 構造変更になり、影響は既定外の温度に限られる。
- sp の対象を rankable 候補に絞る — 失敗した saddle のパネルにも価値がある。
- report の knob や `hfauto report` を削除する — 数行で、複数条件の run を読み直すのに必要。
- （U8-IMP-6 の元案）HTML に理論レベルのラベルを出す — 前提に頼るロジックが増える。
- （U8-IMP-8 の元案）0.8 × memory_mb の係数を持ち込む — 新しい係数を増やさない。
- （U8-IMP-4 の代替）UHF の出力も読めるようにパーサを直す — ⟨S²⟩ を検査できない UMP2 の値がパネルに入る。

#### 2.8.9 変更不要な点

- sp の構造を freq の最終構造に固定し、電荷と多重度を freq.level から取り、parents をつなぐこと。
- Level の観測と照合（WFT は method・基底・版数）、構造の echo の検査。CCSD(T) の開殻を起動前に拒否すること。
- パネルのキーを Level.full_key にすること（CH-25）、`freeze atomic` と thresh 1e-6、sp を直列に実行すること。
- rankable を唯一の関門にし、信頼度スコアを作らないこと（AR-27）、1,1,3 方式の同順位、coverage.csv。
- JobStore のキーの設計、終端的な失敗の記録と `--retry-failed`、4 種の ladder、SCF の再投入。
- SiteLock、セマフォ、`/mnt` の拒否、timeout。追記型の method_panel、自己完結の report.html。

---

## 3. 横断的な課題と改良案（重複を統合）

### X1 低レベルの陰性結果を、DFT の前の拒否権に使っている（U4-I1、U5-I1）— critical / must

- **問題**: 仮説の両端は、行 1 により DFT の別々の basin です。山越えの定理から、その間には 1 次の鞍点が必ずあります。それなのに、次のような GFN2 の単端探索の陰性結果を、理由も機構も問わずに集め、ジョブを 1 本も出さずに NO_PRODUCT（no_product_basin）で終えます。
  - 同じ組成で drive が一致する NT2 の単調上昇、AFIR の same_as_source、窓外、数値的な不成立
  - 不具合の出方は次の 2 つです。
    - AFIR で見つかった生成物には、同じ trial の NT2 negative が必ず付く。
    - HCN の宣言反応は、GFN2 の障壁が 306 kJ/mol（DFT では 195 kJ/mol）で窓外になるので止まる。
- **統合した改良**: U4-P1 = U5-P1（ゲート一式の削除、【複雑さ減】）。
  - あわせて U4-P6 の簡素化版（`no_ts_guess` と `nt2_failed` を区別し、TS がある negative にも ΔE‡ を載せる）を入れると、coverage で報告される陰性結果の情報量は落ちません（検証で指摘された依存関係）。
- **影響範囲**: W7 ではゲートは一度も発火していないので、既存の結果は変わりません。コストは、該当する仮説ごとに SCREEN 1 回分です（小分子で 1 分未満、17 原子で約 5 分）。

### X2 explore の出発構造と、生成物の判定基準（U3-ISS-2、U4-I2、U4-I3）— high / must + should

- **問題**
  - explore は `species.geometry`（入力の構造）から出発します。ReaDuct は、それを緩和した source を基準にします。soft 押しが働いた系では、source が screen 極小より高い別の配座に落ち、screen 極小そのものが「生成物」になります（TMA·(HF)₂ では、DFT で 1,308 s、dft stage の 43%）。
  - さらに、一致の判定が RMSD < 0.1 Å なので、柔らかい錯体では配座の違いが別の化学種として扱われます。
- **統合した改良**
  - U4-P2 = U3-IMP-3（`inputs.evidence(minimum.opt_calc).final` から出発する。約 5 行）— must。
  - U4-P3（`matches_source` を state_label で判定する。【複雑さ減】）— should。
- 2 つは互いに補い合います。P2 で原因を取り除き、P3 で柔らかい錯体の IRC や緩和の後に起きる配座の違いも吸収します。

### X3 小さな負の振動数を熱化学から捨てている（U3-ISS-3、U7-I1、U7 の検証で追加された TS の −10〜0 cm⁻¹）— medium / must

- **問題**
  - minima は、noise（−10〜0）を極小、押しても残った soft（−50〜−10）を soft_minimum として採用します。TS も、副次の負モード（−50〜0）を許します。
  - ところが thermo は、既定（invert_soft_cm1=None）でこれらをすべて捨てます。GoodVibes は `im_frequency_wn` を ZPE・U・S のどれにも入れないためです。
  - その結果、同じ |ν| でも、符号だけで G が 0.5〜1.6 kcal/mol 変わります。低振動モードの多い amine·(HF)n の錯体は、まさに対象です。
- **統合した改良**: U7-P1 の修正版（`thermo_frequencies(freqs, saddle)` という純関数で、極小は全負モードを |ν| にし、TS は最低モードだけを除く）— 【複雑さ減】。U3-IMP-4 の修正案（極小は invert を saddle_cm1 に固定する）は、これに含まれます。
- **注意**: thermo の JobStore キーは振動数のリストを含みません。規則を変えたら、キーを必ず変えてください（X12）。

### X4 手法パネルの構成と符号ゲート（U0-I1、U0-I2、U0-I8、U8-1、U8-2、U8-3）— high / must + should + could

- **問題**
  - (a) 許容幅 0 の符号ゲートのせいで、検証済みの HONO が外れます。panel_report と report で、矛盾する 2 つの順位ができます。
  - (b) diffuse のない def2-SVP と def2-TZVP は、陰イオンや H 結合の系で、誤った方向の「参照」になります。
  - (c) 既定のパネルに入っている CCSD(T) は、主な対象（15〜17 原子）では timeout を超え、再実行のたびにやり直しになります。
  - (d) CCSD の maxiter 20 は、余裕が小さいです。
- **統合した改良**（どれも個別に採用されたもので、組み合わせた結果を示す）

| 改良 | 内容 | 優先度 |
|---|---|---|
| U0-P1 = U8-IMP-1 | 符号ゲートに ±1 kcal/mol の不感帯（モジュール定数）を入れ、def2-SVP をパネルから外す【複雑さ減】 | must |
| U0-P2 | TZ 層を def2-TZVPD にそろえる（PBE0 と CCSD(T)） | should |
| U8-IMP-2 | CCSD(T) を既定から外し、小さい閉殻分子での opt-in にする【複雑さ減】 | should |
| U0-P6 = U8-IMP-3 | CCSD の maxiter を 50 にする | could |
| 文書（U8 の検証で追加） | 「パネルを追記した後は、panel_report が最終の順位」と 1 行書く | could |

- **組み合わせた後の既定のパネル**: `[pbe0-d3bj_def2-tzvpd, wb97x-d3_def2-tzvpd]`（これに停留点レベルの参照行が加わる）。小さい閉殻分子では、`ccsd-t_def2-tzvpd` を 1 語足して使う。
- **同時に直すもの**: `tests/integration/test_pipelines_e2e.py:96` の期待値、`docs/design.md:162`、`refactor_design.md:1277` と、同書の「エネルギー層（def2-TZVP）」の記述。

### X5 composite G（SP のエネルギー層）が配線されておらず、fail-open でもある（U0-I3、U0-I4、U8-5）— medium / should + could

- **順番**
  1. 先に U0-P4（`energy_layer_missing` で fail-closed にする）を入れる。
  2. 次に U0-P3 の修正版（追記用の `composite.yaml` を別に用意し、`wb97x-d3_def2-tzvpd` を targets: all_minima で使う。`barrier_vanishes_at_energy_layer` の blocker も足す）を必要に応じて入れる。
  3. 説明は U8-IMP-5 で実態に合わせる。
- 既定の流れには組み込みません（U8 の見送り案）。主な利点は ΔG_assoc で、NH3·HF の誤差が 3.7 → 1.45 kcal/mol になります。障壁の改善は小さいです。

### X6 同一性の基準の統一と、stage 内での freq の省略（U3-ISS-1、U3-ISS-4）— high / should

- U3-IMP-1（relax したら即座に登録する）と U3-IMP-2（assign の 1 基準にし、ambiguous・rejudge・tight を削除する）は、**同時に入れます**。片方だけだと、stage 内で 2 つの基準が混在し、どちらが効くかが登録順で変わります。
- 効果: TMA·(HF)₂ の dft stage が −30%（3,050 s → 約 2,150 s）。約 40〜60 行が減ります。

### X7 計算時間の大半を占める部分への対応（U2、U3、U4、U5、U6、U8 の横断）

| 無駄の発生源 | 実測 | 対応する改良 | 見込み |
|---|---|---|---|
| 同じ basin の 2 本目の DFT freq | 17 原子で 1 本 約 900 s | U3-IMP-1 + IMP-2 | TMA·(HF)₂ の dft stage が −30% |
| explore の偽の生成物 | DFT 1,308 s（43%） | U4-P2（+ U4-P3） | 同じ系で opt と freq が 2 本不要になる |
| QRC 側の opt に Hessian がない | HONO で paths の 81%、126 ステップ | U6-P1 | HONO の paths 386 s → 約 120 s（推定） |
| 15 bead の barrierless 確認 | 17 原子で約 5〜6 h（推定） | U5-P4 | 予算内に閉じるようになる |
| TSOpt の例外で GS が失われる | FIND_PATH（1 h 以上/チャンク）に落ちる | U5-P5 | xTB の GS を 3〜15 s でやり直せば済む |
| proceed でも tangent 用の DFT freq を要求する | JobStore がミスすると 1 本 約 900 s | U5-P6 | 不要になる |
| 大きな系での CCSD(T) | 1 点数時間（推定）、再実行のたびに繰り返す | U8-IMP-2 | 既定では走らない |
| CREST のトポロジー停止と再実行 | (HF)₃ で 87〜236 s、1 系は失敗 | U2-IMP-1 | 23〜80 s で 3/3 成功 |

### X8 入力と設定の誤りは、読み込み時に止める（U0-I5、U1-01/02/04/05、U2-ISS-4、U3-ISS-5、U7-I9、U8-9）

- **共通の方針**: 誤った設定は、計算を始める前に理由付きで拒否する。効かない knob は外す。
- **該当する改良**
  - U0-P5: 停留点 method の一意性
  - U1-I1: 宣言反応の端点は xyz（must）
  - U1-I2: 宣言座標の添字と、組成の一致
  - U1-I4: xyz と SMILES の排他、個数 ≥1、id の一意性【複雑さ減】
  - U1-I5: 未知元素は fail-closed【複雑さ減】
  - U2-IMP-3: conformers の設定整理【複雑さ減】
  - U3-IMP-5: include=window で screen 極小がなければ拒否
  - U7-P8: thermo の temperatures_K の削除【複雑さ減】
  - U8-IMP-8: memory の preflight 検査
- 新しい knob は増やさない、という点を、どの修正版でも守っています。

### X9 開殻と多重度は「黙って決めない・黙って誤らない」に統一する（U1-03、U1-08、U2-ISS-3、U3-ISS-8、U7-I2、U8-4）

- 対象の化学は閉殻が中心なので、開殻への対応は広げません。そのかわり、暗黙に決めていた箇所や、常に失敗する箇所を、理由の分かる明示的な扱いにそろえます。
- **該当する改良**
  - U1-I3: SMILES のラジカル電子と同位体の拒否
  - U2-IMP-4: 組成の多重度が一意でなければ INPUT_INVALID
  - U7-P2: 会合量の LOT 判定から電荷と多重度を外し、単量体のキーに多重度を入れる
  - U8-IMP-4: WFT は閉殻専用【複雑さ減】
  - 文書: U1-I6、U3-IMP-6（m=1 は RKS で、開殻一重項は扱わないこと）

### X10 fake エンジンでしか通っていない経路（検証の空白）

- 実計算で通ったのは SCREEN → REFINE_SADDLE → VALIDATE_TS → CONNECT の 1 経路だけです。次の経路は、fake のテストでしか確かめられていません。
  - FIND_PATH、IDPP、低レベル TS の近道、行 12〜14、moddir、higher_order、collapse / 分割、QRC の再試行、_register、reassigned、periodic_nearest、negative_evidence、CREST の停止と再実行、composite G、method_panel（本レビューで初めて実行した）。
- 今回の重大な課題のうち、X1、U5-I2/I3 は、fake が実エンジンと違う挙動をしていた（FakePath は gmax が下がって収束すると仮定していた）ために見逃されていました。
- **該当する改良**
  - U5-P8 の修正版（FakePath を実エンジンの挙動に合わせる。軽い基底で string の smoke を 1 本加える）
  - U2-IMP-1 に付けた smoke（F⁻·(HF)₂ の組成）と、再実行による確認
  - U6-P1 の実計算での再確認

### X11 文書化（design.md と validation.md への追記を 1 か所にまとめる）

| 追記先 | 内容 | 由来 |
|---|---|---|
| validation.md | def2-SVPD を選んだ根拠。約 1.5 kcal/mol 未満の順位差は手法誤差の範囲。ΔG_assoc は非補正の SVPD の値。amine·(HF)n は PBE0 の PES 上の結論 | U0-P7 |
| validation.md | NH3 の +0.411 は種の G に入る σ で、速度に対応する有効障壁では 3.44 | U7-P4 |
| design.md §7 structures / conformers | SMILES は 1 配座で、立体は任意、同位体は不可。端点は xyz。組成の多重度は高スピンで、m=1 は RKS | U1-I6 |
| design.md conformers | `--nci --quick` の実効設定、seed00 だけを CREST に渡すこと、衝突条件、`--noopt` の理由 | U2-IMP-5 |
| design.md §5.5 と §7 | 振動数の射影の約束（NWChem と数 cm⁻¹ ずれる）と、開殻一重項は扱わないこと | U3-IMP-6 |
| design.md thermo | QH=False、G に入るスケールは 1 つ、σ だけで L は含まない、トンネル補正なし、1 atm。composite の使い方 | U7-P4、U8-IMP-5 |
| design.md report | 符号ゲートは \|ΔE_rxn\| > 1 kcal/mol の値どうしで判定する。パネルを追記した後は panel_report が最終 | U0-P1、U8 の検証 |
| design.md paths | ZTS の `tol` は変位の閾値（NWChem の文書の記述とは違う） | U5-I10 |

### X12 キャッシュの無効化と旧 run との互換

次の改良は、JobStore のキーか settings の sha を変えます。U3-IMP-2（NWChem optimize のキーから tight が消える）、U4-P2（explore の出発構造が変わる）、U7-P1 / P3 / P8（thermo の settings）。旧 run とは互換を取らない方針（design）の範囲内ですが、次の点に注意してください。

- U7-P1: 規則を ThermoSettings の外で変えると、既存の結果が黙って再利用されます。キーを確実に変えてください。
- U7-P3: validation.md の値が小数第 2 位で変わるので、記録を更新してください。

---

## 4. 改良ロードマップ

工数は実装とテストの目安です。「複雑さ」は、コードと設定の量の増減を表します。

### 4.1 must（結論の誤り、または機能不全を直す）

| # | 改良 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|
| M1 | negative_evidence ゲートの削除（U4-P1 = U5-P1） | DFT の basin がある反応や宣言反応を、NO_PRODUCT にする誤りがなくなる | 約 40 行の削除、テスト 2 件 | 【減】 | U4-P6 の簡素化版を同時に入れることを推奨 |
| M2 | FIND_PATH の判定の置き換え（U5-P2 の修正版） | DFT string が初めて機能する（monotonic / multi_max / single_max） | 約 30 行（削除が中心）、FakePath と golden の差し替え | 【減】 | monitor とチャンクのループの `stagnated` の両方を外すこと |
| M3 | explore の出発構造を screen の最適化構造にする（U4-P2 = U3-IMP-3） | 偽の生成物（DFT 1,308 s）がなくなり、run ごとの結論のぶれもなくなる | 約 5 行、テスト 1 件 | 中立 | キャッシュが変わる |
| M4 | CREST の組成に `--noopt`（U2-IMP-1） | (HF)₃ 系が 3/3 成功し、より低い構造が見つかり、時間は約 1/3 | 1〜5 行、smoke 1 件 | 【減】 | tma_hf2 の conformers を再実行して確認する |
| M5 | 符号ゲートの ±1 kcal/mol 不感帯と、SVP の除外（U0-P1 = U8-IMP-1） | HONO が rank 1 に戻る | 2 行、テスト 1 件、ファイル 1 本の削除 | 【減】 | e2e テストと design の method 一覧も直す |
| M6 | 熱化学の負モードの固定規則（U7-P1 の修正版 = U3-IMP-4） | 1 モードあたり 0.5〜1.6 kcal/mol の偏りと、不連続がなくなる | 20〜30 行、テスト 3 件 | 【減】 | キャッシュキーを必ず変える |
| M7 | 宣言反応の端点に xyz を必須にする（U1-I1） | SMILES の端点で別の機構を黙って計算することを防ぐ | 約 5 行、テスト 1 件 | 中立 | 同梱の system は影響を受けない |

### 4.2 should（有用性・効率・堅牢性の改善）

| # | 改良 | 効果 | 工数 | 複雑さ | 依存・注意 |
|---|---|---|---|---|---|
| S1 | その場での登録と、同一性基準の統一（U3-IMP-1 + IMP-2） | dft stage が −30%（TMA·(HF)₂） | 約 1 日 | 【減】 | 必ず 2 つ同時に入れる |
| S2 | GS の画像の逐次整列（U5-P3） | FIND_PATH の初期経路と screen_hei の接線が正しくなる | 3〜5 行 | 中立 | M2 と同時に入れる |
| S3 | 行 14 の 15 bead 確認の撤廃（U5-P4） | 17 原子で約 5〜6 h を省ける | 約 15 行（削除が中心） | 【減】 | M2 の後に入れる |
| S4 | TSOpt が例外で終わったときの GS 経路の回収（U5-P5） | SCREEN が unavailable になるのを防ぐ | 約 6 行、テスト 1 件 | 【増】 | GS 自体の失敗は握りつぶさない |
| S5 | FakePath の修正と、string の smoke（U5-P8 の修正版） | I2・I3 型の不具合を回帰テストで捕まえられる | テスト約 40 行 | 中立 | 軽い基底で opt-in にする |
| S6 | `matches_source` をグラフで判定する（U4-P3） | 配座の違いによる誤った生成物と「非接続」がなくなる | 約 5 行 | 【減】 | M3 と組み合わせる |
| S7 | QRC の opt に TS の Hessian を渡す（U6-P1 の修正版） | HONO の QRC が 126 → 38 ステップ（trust 0.1 のまま） | 約 15 行 | 中立 | trust 0.3 は再確認してから入れる |
| S8 | パネルの TZVPD 化と、CCSD(T) の opt-in（U0-P2 + U8-IMP-2） | 陰イオン系でパネルの参照が正しくなり、大きな系での事故を防ぐ | 設定ファイルの置き換え | 【減】 | X4 の最終形を参照 |
| S9 | energy_method の fail-closed（U0-P4） | composite を要求したときの誤報告を防ぐ | 3〜5 行 | 中立 | U0-P3 の前提 |
| S10 | 停留点 method の一意性を読み込み時に検査する（U0-P5） | 設定の食い違いを計算前に拒否できる | 約 6 行 | 中立 | — |
| S11 | 宣言座標の検査（U1-I2 の簡素化版） | 誤った添字を 1 秒未満で止める | 約 12 行 | 中立 | — |
| S12 | system の静的な検査（U1-I4） | id の衝突による黙った上書きを防ぐ | 約 8 行 | 【減】 | — |
| S13 | conformers の窓の削除（U2-IMP-2 の修正版） | GFN2 段階での恣意的な切り捨てがなくなる | 数行 | 【減】 | None の規則は残す |
| S14 | 会合量の LOT 判定と、単量体のキー（U7-P2 の修正版） | 陰イオン・ラジカルの会合量が出るようになる | 約 10 行 | 中立 | キーに多重度を入れる |
| S15 | スケール因子の 1 本化（U7-P3） | 約 35 行と設定 1 つが減り、ZPE と G のスケールが揃う | 約 1 時間 | 【減】 | validation の値を更新する |
| S16 | 熱化学の規約の文書化（U7-P4 の修正版） | ΔG‡ の意味が正しく伝わる | 30 分 | 中立 | — |
| S17 | WFT は閉殻専用にする（U8-IMP-4） | 常に失敗する経路を、理由の分かる拒否にする | 数行 | 【減】 | — |

### 4.3 could（余裕があれば）

| # | 改良 | 効果 | 工数 | 複雑さ |
|---|---|---|---|---|
| C1 | composite 用の追記 pipeline（U0-P3 の修正版） | 会合量を TZVPD の層で出せる | 設定約 2 行、コード約 2 行 | 中立（S9 の後） |
| C2 | CCSD の maxiter 50（U0-P6 = U8-IMP-3） | 収束の遅い系で参照点を落とさない | 1 行 | 中立 |
| C3 | 停留点レベルの根拠と読み方の文書化（U0-P7 の修正版） | 結果の過大解釈を防ぐ | 3〜5 行 | 中立 |
| C4 | SMILES の同位体とラジカルの拒否（U1-I3 の修正版） | 質量や電子状態を黙って取り違えることを防ぐ | 約 6 行 | 中立 |
| C5 | 未知元素の fail-closed（U1-I5） | エラーの理由が正確になる | 2 行 | 【減】 |
| C6 | structures と組成の規則の文書化（U1-I6） | 制約が利用者に見える | 3〜4 行 | 中立 |
| C7 | conformers の設定整理（U2-IMP-3 の修正版） | 効かない knob が 3 つ減る | 数行 | 【減】 |
| C8 | 組成の多重度が一意でなければ止める（U2-IMP-4 の修正版） | ラジカル対で PES を暗黙に選ばなくなる | 2〜3 行 | 中立 |
| C9 | CREST の実効設定の文書化（U2-IMP-5 の修正版） | 再現性が上がる | 文書のみ | 中立 |
| C10 | window で screen 極小がなければ拒否（U3-IMP-5 の修正版） | 空の結果を「反応なし」と読み違えることを防ぐ | 数行 | 中立 |
| C11 | 射影の約束と開殻一重項の注記（U3-IMP-6） | 誤解を防ぐ | 数分 | 中立 |
| C12 | drive の対称性による重複除去（U4-P4 の簡素化版） | 一般の有機分子で、trial の種類が網羅される | 15〜20 行 | 中立 |
| C13 | 無バイアスの opt の上限を 500 にする（U4-P5） | 端点の opt の打ち切りによる取りこぼしを防ぐ | 1 行 | 中立 |
| C14 | 陰性ラベルの最小限の改善（U4-P6 の修正版） | 陰性結果が定量的になる（M1 と同時が望ましい） | 約 10 行 | 中立 |
| C15 | tangent の遅延計算（U5-P6 の修正版） | 不要な DFT freq を省ける | 約 5 行 | 中立 |
| C16 | 両端より低いノードの検出（U5-P7 の修正版） | 中間体の取りこぼしを防ぐ（M2 の後） | 3 行 | 【増】 |
| C17 | moddir の番号を driver 座標で決める（U6-P2 の修正版） | 負のモードが複数あるとき、正しいモードを追う | 約 10 行 | 中立 |
| C18 | periodic_nearest に判定を付ける（U6-P3 の修正版） | 回転異性体の取り違えを防ぐ | 数行 | 中立 |
| C19 | QRC の振幅のクリップ順と drop の一本化（U6-P4） | 上限が仕様どおりになり、重複が 1 つ減る | 数行 | 【減】 |
| C20 | 点群の記録（U7-P5 の修正版） | σ を監査できる | 約 5 行 | 中立 |
| C21 | ΔG_assoc の感度の幅（U7-P6 の修正版） | 弱い錯体の判断材料になる | 約 10 行 | 【増】 |
| C22 | 有効障壁 ΔG‡_eff の表示（U7-P7 の修正版） | 反応の向きに依存しない指標を並べて見られる | 約 10 行 | 【増】 |
| C23 | thermo の temperatures_K の削除（U7-P8） | 効かない設定がなくなる | 約 15 行 | 【減】 |
| C24 | composite の説明を直す（U8-IMP-5） | 誤読と誤用を防ぐ | 数行 | 【減】 |
| C25 | ranking.csv に T・標準状態・ΔG_rxn を加える（U8-IMP-6 の修正版） | CSV が単独で読める | 数行 | 中立 |
| C26 | TS のない outcome で dzpe を検査しない（U8-IMP-7） | blockers 列が正確になる | 1 行 | 中立 |
| C27 | memory の preflight 検査（U8-IMP-8 の修正版） | 設定ミスによる OOM を事前に止める | 3〜5 行 | 中立 |

### 4.4 進め方と確認

1. **must の 7 件**は、どれも数行〜数十行で、多くは削除です。順番は M1 → M2（+ S2）→ M3 → M4 → M5 → M6 → M7 を推奨します。M1 と M2 は reaction_case の同じ範囲を触るので、続けて入れると差分の確認が楽です。
2. **確認の実計算**（どれも本レビューの入力に書かれた手順）
   - M4 の後: tma_hf2 の conformers を再実行し、同じ 2 配座が得られることを確かめる（約 11 s）。F⁻·(HF)₂ の smoke（約 25 s）も通す。
   - M2 と S5 の後: HCN の ZTS を軽い基底で smoke する（STO-3G、1 rank で約 78 s）。
   - S7 の後: W7 の known_endpoints の 3 反応（約 12 分）と ab_init_hessian を再実行し、trajectory_above_ts が出ないことを確かめてから trust 0.3 を入れる。
   - M6 と S15 の後: validation.md の ΔG の値を記録し直す（小数第 2 位が変わる）。
3. **should のうち S1（同一性の統一）と S7（QRC の Hessian）** は、時間の削減効果が最も大きいです。どちらもコードは減るか、増えません。

---

## 5. 参考文献・仕様 URL 一覧

**文献**

- Bursch, Mewes, Hansen, Grimme, "Best-Practice DFT Protocols for Basic Molecular Computational Chemistry", Angew. Chem. Int. Ed. 2022, 61, e202205735 — https://pmc.ncbi.nlm.nih.gov/articles/PMC9826355/
- Rappoport, Furche, J. Chem. Phys. 2010（def2-SVPD / TZVPD）— https://doi.org/10.1063/1.3484283
- HCN/HNC の CCSD(T)/ANO 参照値 — https://cdnsciencepub.com/doi/pdf/10.1139/v96-120
- FHF⁻ の結合エネルギー（Wenthold & Squires）— https://doi.org/10.1021/j100007a034
- HONO cis/trans のエネルギー差（J. Mol. Spectrosc. 2009）— https://www.sciencedirect.com/science/article/abs/pii/S0022285209002574
- Goerigk ら GMTKN55, PCCP 2017 — https://doi.org/10.1039/C7CP04913G
- Prasad ら BH9, J. Phys. Chem. A 2022 — https://pubs.acs.org/doi/10.1021/acs.jpca.2c03922
- GFN2-xTB のベンチマーク — https://benchmarks.rowansci.com/methods/gfn2-xtb 、https://benchmarks.rowansci.com/methods/wb97x-d
- GFN2-xTB の気相 PT の精度, JCTC 2025, 21, 7149 — https://pmc.ncbi.nlm.nih.gov/articles/PMC12288007/
- Pracht, Bohle, Grimme, PCCP 2020, 22, 7169（CREST / iMTD-GC）— https://pubs.rsc.org/en/content/articlelanding/2020/cp/c9cp06869d
- Riniker & Landrum, JCIM 2015（ETKDG）— https://doi.org/10.1021/acs.jcim.5b00654
- Wang, Witek, Landrum, Riniker, JCIM 2020（ETKDG v3）— https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025
- Cordero ら, Dalton Trans. 2008, 2832（共有結合半径）— https://doi.org/10.1039/b801115j
- Balcilar ら, ICML 2021 supp.（1-WL の限界）— http://proceedings.mlr.press/v139/balcilar21a/balcilar21a-supp.pdf
- Maeda ら, J. Comput. Chem. 2018, 39, 233（GRRM17 / AFIR）— https://pmc.ncbi.nlm.nih.gov/articles/PMC5765425/
- Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. 140, 214106 (2014)（IDPP）
- E, Ren, Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007)（ZTS）
- Henkelman, Uberuaga, Jónsson, J. Chem. Phys. 113, 9901 (2000)（CI-NEB）
- Goodman & Silva, Tetrahedron Lett. 44, 8233 (2003)（QRC）、pyQRC — https://github.com/patonlab/pyQRC
- Grimme, Chem. Eur. J. 2012, 18, 9955（qRRHO）
- Luchini ら, F1000Research 2020, 9, 291（GoodVibes）
- Alecu, Zheng, Zhao, Truhlar, JCTC 2010, 6, 2872（スケール因子）
- Fernández-Ramos, Ellingson, Meana-Pañeda, Marques, Truhlar, Theor. Chem. Acc. 2007, 118, 813（対称数と反応速度）— https://link.springer.com/content/pdf/10.1007/s00214-007-0328-0.pdf
- Skodje & Truhlar, J. Phys. Chem. 1981, 85, 624（トンネルの T_c）
- Kozuch & Shaik, Acc. Chem. Res. 2011, 44, 101（energetic span）
- Miller, Handy, Adams, J. Chem. Phys. 1980, 72, 99（Eckart の射影）

**ライブラリの仕様・ソース**

- NWChem DFT — https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html
- NWChem Geometry Optimization — https://nwchemgit.github.io/Geometry-Optimization.html
- NWChem Hessians and Vibrational Frequencies — https://nwchemgit.github.io/Hessians-and-Vibrational-Frequencies.html
- NWChem NEB / ZTS — https://nwchemgit.github.io/Nudged-Elastic-Band-and-Zero-Temperature-String-Methods.html
- NWChem CCSD — https://nwchemgit.github.io/CCSD.html 、MP2 — https://nwchemgit.github.io/MP2.html
- NWChem 7.2.3 のソース:
  - https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/driver/opt_drv.F
  - https://github.com/nwchemgit/nwchem/blob/v7.2.3-release/src/optim/string/string.F
  - https://github.com/nwchemgit/nwchem/blob/master/src/vib/vib_eckart.F 、…/vib_vib.F
  - https://raw.githubusercontent.com/nwchemgit/nwchem/v7.2.3-release/src/ccsd/ccsd_iterdrv2.F 、…/aoccsd2.F、…/src/task/task.F
- xtb — https://xtb-docs.readthedocs.io/en/latest/commandline.html 、…/optimization.html、…/hessian.html、…/xtb_thermo.html、…/xtb_docking.html、xcontrol: https://www.mankier.com/7/xcontrol
- CREST — https://crest-lab.github.io/crest-docs/page/documentation/keywords.html 、https://crest-lab.github.io/crest-docs/page/examples/example_3.html
- CREST 3.0.2 のソース — https://github.com/crest-lab/crest/blob/v3.0.2/src/ （confparse.f90、choose_settings.f90、algos/setuptest.f90、cregen.f90、crest_main.f90）
- RDKit — https://www.rdkit.org/docs/GettingStartedInPython.html 、https://www.rdkit.org/docs/RDKit_Book.html
- pysisyphus — https://pysisyphus.readthedocs.io/en/latest/chainofstates.html （venv 内の run.py、optimizers/Optimizer.py、cos/GrowingString.py）
- SCINE ReaDuct 6.1.0 — https://raw.githubusercontent.com/qcscine/readuct/6.1.0/src/Readuct/App/Tasks/ （IrcTask.h、AfirOptimizationTask.h、NtOptimization2Task.h）
- SCINE Utilities 10.1.0 — https://raw.githubusercontent.com/qcscine/utilities/10.1.0/src/Utils/Utils/GeometryOptimization/NtOptimizer2.cpp
- SCINE Puffin — https://raw.githubusercontent.com/qcscine/puffin/master/scine_puffin/jobs/scine_react_complex_nt2.py 、…/templates/scine_react_job.py
- SCINE Chemoton 4.0.0 — https://scine.ethz.ch/static/download/documentation/chemoton/v4.0.0/py/api/utilities.html
- GoodVibes 4.3.0（WSL prod venv の `goodvibes/thermo.py`、`api.py`、`scaling_factors.json`）、pymsym 0.3.5（`pymsym/high_level.py`）

**リポジトリ内の関連文書**

- `docs/design.md`、`docs/validation.md`、`docs/environment.md`
- `docs/reviews/2026-09-25_architecture_chemistry_review.md`（CH-xx）
- `docs/reviews/2026-09-25_refactor_design.md`

**本レビューのプローブ**（WSL、読み取り専用の複製）

- `/home/user/hfauto_review_probe/{u0, u1, u1b, u2_conformers, u2_assess, u4_nt2, u4_assess, u5, u5_xtb_reopt, u6, u6_assess, u7, u7_assess, u8, u8b}`

