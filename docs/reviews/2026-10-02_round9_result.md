# hfauto round 9 — 結果報告(2026-10-02)

- 対象: ブランチ `refactor/fundamental-2026-09`。r9-PRE `8670200` から S6 `c6eff84` までの 14 コミットと FINAL(整理、文書、VAL9)。push はしていない。
- 基準: 分析 [2026-09-30_remaining_issues_analysis.md](2026-09-30_remaining_issues_analysis.md)(以下「分析」、HEAD `2c3c894`)の §4 のロードマップ(must の M1〜M9 と、should をまとめた 6 つの wave S1〜S6)。
- 方法: 振る舞いを保つ変更は JobStore の再生ゲート(新しい QM ジョブ 0、記録が同一)で、振る舞いを変える変更は再生の差の化学的な説明と、`/home/user/hfauto_r9/<wave>/` の実 run で受け入れた。しきい値・予算・期待値を、分岐を通すために動かしたことはない(分析 §4.4)。変えた期待値は理由とともに [validation.md](../validation.md) §4 にある。
- 全体の再検証は VAL9(`/home/user/hfauto_r9/VAL9`、41 run と smoke。期待の 38 行がすべて合格)と、その strict の自己再生(40/40 PASS)([validation.md](../validation.md) §3)。

## 0. 要旨

- **正しさ**: 仮説の結合変化を探索と判定の一次情報にした。鞍点探索の初期 Hessian は結合変化 ρ だけに負の曲率を置き、虚モードが ρ を担うこと(χ ≥ 0.3)を TS の必要条件にした。S6 の H 引き抜きは、VAL7 の「回転子の鞍点で試行を使い切る」から、3 種 × rank 2・4 の 6 run と S3 以後の再生の 3 本のすべてで引き抜きの TS(M06-2X の分離基準 ΔE‡ 5.7〜5.8、CCSD(T) 6.96)に届くようになった。派生する端点は到達した添字のまま続き(偽の 4 結合変化が消えた)、結合が変わる case は化学状態の粒度で判定する。
- **問いの形**: 障壁のない会合は、分離した単量体から付加体への緩和スキャン(低スピン結合なら Yamaguchi の AP)で問う。S5 CH3• + O2 は「結論なし」から barrierless(ΔH₀ −31.6、実験の D₀ 約 32)になった。開殻のプロファイルは端の SCF の解から解き、偽の山(+28.9 kcal/mol)が消えた。
- **精度**: 既定の順位を M06-2X-D3(0)/def2-TZVPD // PBE0-D3BJ/def2-SVPD にした。CCSD(T) との障壁の差は平均 2.55 → 0.83、最大 6.2 → 1.8 kcal/mol(5 反応)。ΔE_rxn は良くならない(平均 1.2、最大 2.3)。PES の形は PBE0 のまま。
- **再現性と判定の分解能**: 極小の停留性を勾配と Hessian の二次モデルで判定し(平らな PES の肩を極小にしない)、点群を 1 つの判定にし(S18 の 2 run の差 2.27 → 0.006 kcal/mol)、結合変化に ±0.1 Å の帯を付け、trial を原子の番号付けに依らなくした。
- **簡潔さ**: QRC の変位、像の判定、DFT 極小のストア、basin の構造の写像、対称性の判定は 1 か所ずつになった(§3)。決定表は 17 → 14 行。ただし SLOC は `2c3c894` から +538(+7.1%)で、増えた分は新しい振る舞い(χ、ρ、スキャンと AP、停留性、点群、SCF の guess、記録の欄)である(§2)。
- **費用**: エネルギー層は 36 run で SP 105 本・11,860 core 秒を足した。S19 の失われた seed の DFT(9,406 core 秒)は要らなくなった。HCN のクイックスタートは 1 分 10 秒(QM 242 core 秒)。

## 1. wave ごとの結果

SLOC と docstring は `radon raw` の hfauto の合計(docstring は multi の行)。差は直前のコミットから。

| wave | コミット | 主な内容 | 主な削除 | 検証 | SLOC | docstring |
|---|---|---|---|---|---|---|
| PRE | `8670200` | JobStore が壁時計以外の失敗も保存する。再生ゲート(r9 の tools) | — | HEAD の自己一致 36/36、陰性対照 3/3 が FAIL | −7 | +5 |
| M1 | `769cde2` | 派生する端点は到達した添字を保つ。`DiscoveryRecord.source_species`。basin の構造を WL クラスで端点の添字に並べる | 端を basin の代表やメンバーの対から選ぶ処理 | 端の監査の不一致 6/17 → 0。W3_s6 の偽の 4 結合変化が消えた。マロンアルデヒドを登録順を変えて 2 回 | +10 | +29 |
| M2 | `fb7d7ce` | 結合が変わる case は化学状態、同じ状態の case は basin で判定。分解能で端点と同じ井戸 | 「中間体は他方の端点の配座でもよい」 | strict 32/36(差を挙げた 4 本は報告)。S6 の 2 段、S5 の −72i が会合の TS でないこと | +33 | +18 |
| M3 | `aec5cf9` | 反応モード性 χ を QRC の前の必要条件に(0.3) | — | χ の分布の空白(≤ 0.023 と ≥ 0.571)。strict 33/36。S6 の回転子の鞍点を QRC なしで棄却 | +35 | +11 |
| M4 | `eccc4ac` | 鞍点の初期 Hessian は ρ か種の TS モードに負の曲率 1 本。NWChem は常に `--bind-to none` | 重なりの検査、接線、`profile.tangent`、`topology.reaction_centre` | S6 を 3 種 × rank 2・4 で 6 本すべて引き抜きの TS。saddle の鍵が変わるので strict は 10 本、ほかは同じ TS を報告 | −10 | +7 |
| M5〜M7 | `298f6ca` | 行 11 の削除。失われた状態を帯で判定。熱化学の点群を 1 つに(`symmetry.analyze`) | 行 11、`_soft_ts_failed`、`thermo.symmetry_number`、`MinimumRecord.chiral`、回転定数の床 | 再生 misses 0、差 732 件をすべて説明。揺らぎで 110 subject が不変 | +75 | +31 |
| M8 | `a10516f` | 低スピン結合の組成の分類。会合を分離した単量体からの緩和スキャンで問う。スピン交差の理由 | `open_shell_singlet_unsupported` の特例 | S5 barrierless(ΔE_e −37.56)、H + C2H4、BH3 + NH3。strict 36/36(S5 は M8 の実 run から) | +134 | +60 |
| M9 | `92491f3` | 既定の順位を M06-2X の層に。層の定義は thermo の `energy_method` の 1 つ | method_panel の panel_thermo | SP 105 本だけが新しいジョブ、自己一致 36/36、CCSD(T) との比較 | 0 | 0 |
| S1 | `a5c8982` | 致命的な文言だけで autoz を分類。駆動系の継続はすべて最後の frame から。開殻の SCF の救済が ⟨S²⟩ を出す。WFT のメモリ配分と `freeze <n>` | `_scf_rescuable`、autoz の始点からのやり直し | 再生 36 本で新しいジョブ 0、差はパネルの notes 列だけ。autoz の harness 6/6、マロンアルデヒドの CCSD(T) | +67 | +32 |
| S2 | `9a432f4` | 14 行の決定表と続きの段数。状態の組に 1 つの仮説(低レベル TS をすべて種に)。子は結論だけを写す。`split_depth` の記録。例外を `drive_case` で閉じ込める。再開の上限 | 行 1 の level_key の節、行 9・10 の分割、`STRING_CHUNKS`、`saddle_restart` の特別扱い、stage 側の例外処理 | strict 33/36、S6 の 3 本の差を説明。行 1 と `split_depth` を実計算で | +58 | +59 |
| S3 | `1721634` | 極小の停留性(ΔE_N > 5e-5 Eh なら soft)。初期 Hessian の 1 行の規則。結合変化の帯を `bond_changes` に。2 本目の虚モードは停留点で数える | `soft_minimum`、engine の一次の推測、`resolved_bond_changes` | strict 27/36。S6 が停留した一次の TS(ν2 +49)に。ニトロメタンの `soft:resolved` | +54 | +39 |
| S4 | `74f6e4b` | 開殻のプロファイルの SP を端の SCF の解から。枝の跳びは unavailable。低スピン結合のスキャンは AP | — | プローブ G2-1v。strict 30/37。S5 の AP の曲線が単調 | +79 | +38 |
| S5 | `2f8104e` | trial を原子の番号付けに依らなくした。未宣言の縮退転位を D3h の恒等 SN2 の鞍点だけから実証 | 添字の対による順序と同値の分け方、先着の代表 | 置換 20 回 × 51 開始点で同一。strict 30/38(explore の 8 本は trial_id の差を報告)。`sn2_cl_d3h` が degenerate で paths の新しいジョブ 0 | +9 | +5 |
| S6 | `c6eff84` | 変位(`modes.qrc_step`、`off_saddle`)、像の判定(`identity.is_image`)、DFT 極小のストア(`Registry`)、構造の写像(`identity.member_coords`)を 1 か所に。候補選択を `chemistry/selection` へ | `_Point.step`、`mode_amplitude`、`CaseRuntime.minima`、`Ctx.record`、`_end_records`、写像の重複 | strict 38/38、QM 0 | +5 | +25 |
| FINAL | このコミット | §3 の削除の確認(15 項目は削除・統合済み、`identity.is_chiral` は別の問いなので残す)、使われない定義の AST 走査(`tools/deadcode.py`)、docstring の計画・run の ID 77 個の削除、Linux の pyrefly の 7 件、文書、VAL9 | 誰も渡さない `deadline`(CONFORMERS・DISCOVERY のプロトコル)、`modes.amplitude` の `bounds_A`、`perturb_linear` の `seed`、`path_token` の `max_length`、`smiles_to_molecule` の `random_seed`、thermo の「Hessian のない freq」の例外(Evidence が freq の Hessian を保証する) | CUR からの strict 36/36(QM 0)。VAL9 と自己一致(§4) | −4 | 0 |

## 2. 指標

| コミット | 段 | LOC | SLOC | docstring | 空行 | コメント | tests の行数 | 最大のファイル(行) | 最大の CC |
|---|---|---|---|---|---|---|---|---|---|
| `2c3c894` | 分析 | 10,609 | 7,544 | 871 | 1,940 | 339 | 7,591 | `drivers/reaction_case/actions.py` 461 | 14(`minima._pool` など 3 つ) |
| `92491f3` | must の後 | 11,131 | 7,814 | 1,032 | 2,010 | 377 | 8,765 | `backends/nwchem/engine.py` 439 | 16(`reaction_paths.run`) |
| `c6eff84` | S6 | 11,673 | 8,086 | 1,230 | 2,078 | 388 | 10,019 | `backends/nwchem/engine.py` 481 | 18(`hypotheses.select`) |
| FINAL | 最終 | 11,670 | 8,082 | 1,230 | 2,078 | 392 | 10,030 | `backends/nwchem/engine.py` 481 | 18(`hypotheses.select`) |

- 分析 §0.4 のとおり、物理行数の目標は使わず、SLOC と §3 の概念の実装数で見る。`2c3c894` → FINAL で SLOC +538、docstring +359。振る舞いを保つ統合(S6)は SLOC +5、FINAL の整理は −4 で、`2c3c894` を下回らなかった。増えた分は新しい振る舞いで、下回らせるには振る舞いを変える必要がある。docstring の ID は文の中にあったので、消しても行は減らない(行数のために docstring を縮めることはしない)。
- 品質ゲート(FINAL、WSL の prod venv に dev extra を入れて): pytest 668 件、ruff 0、pyrefly 0(1.3.2)、import-linter 5 契約、`radon cc -n D hfauto` が空、500 行を超えるファイルなし(最大 481)。Windows の QA venv でも ruff・pyrefly・import-linter・radon は同じ。

## 3. 概念の実装数

目標は各 1。「今」は S6 の後のコード。

| 概念 | `2c3c894` | 今 | 今の実装と残る別実装 |
|---|---|---|---|
| QRC の変位(1 本のモードに沿う ±、振幅と上限) | 3(`minimum._Point.step`、`actions._amplitude` + `displace`、`actions._pushed`) | 1 | `modes.qrc_step`(`off_saddle` はその和)。ReaDuct の worker の xTB の ± は同じ振幅の式 `modes.amplitude` を使う(Evidence を持たない低レベルの変位) |
| 像の判定(`carry` ≤ `IMAGE_A`) | 2(minima の像のグループ、CONNECT の − 側) | 1 | `identity.is_image` |
| DFT 極小のストア | 2(`Registry._basins`、`CaseRuntime.minima`) | 1 | `Registry`(`Registry.minima`) |
| basin の構造の写像 | 2(`driver._endpoint`、`hypotheses._Pool.coords`) | 1 | `identity.member_coords`。explore の失われた状態は seed(species でない)を `basin_coords` で直接写す |
| 対称性の判定 | 4(外部自由度の階数、回転定数の床、libmsym、`is_chiral`) | 1 | `symmetry.analyze`(熱化学の σ・m・直線性)。外部自由度の階数は点群のない場面だけ、`is_chiral` は写像の掌性だけに使う |
| 結合変化の定義 | 1(帯なし) | 1 | `topology.bond_changes`(±0.1 Å の帯)。縮退の接続の検査 `gates._bonds_exchanged` は帯のない結合グラフを比べる |
| 続きの規則 | 4(ジョブ: timeout と maxiter は最後の frame、autoz は始点から。case: 再開は数えず上限なし、押した種は数える) | 2(層ごとに 1) | ジョブは `LADDER`(駆動系はすべて最後の frame から)、case は `Seed.depth` ≤ `MAX_DEPTH`(再開と押し出しが同じ段数) |

## 4. 検証(VAL9)

[validation.md](../validation.md) §3〜§6。期待値は今のコードの再生(CUR)と round 9 の検証 run の値で、VAL7 から変えた期待値はすべて理由付き(同 §4)。

- 41 run(2026-10-02 01:25〜07:37、smoke 13 件合格): 期待の行 38 がすべて合格、未達 0。孤児 0。
- CUR との結論の比較(36 run): 32 本が同一(δG_eff の差 ≤ 0.0004)。S6 の引き抜きの TS は別の鞍点(−485.3i、ν2 +48.9、χ 0.674)で点群が C1 に決まり、δG_eff は RT ln 2 低い 5.54(CUR 5.95、期待の幅 5.5〜6.0 の中。残る課題 8)。S19 は CREST の配座が 1 つ多く、結論は同じ。パネル 2 本の差は CUR に残った古い stage。
- パネル: M06-2X − CCSD(T) は ΔE‡ で HCN −1.79、SN2 +0.46、HONO −0.01、ΔE_rxn で HCN −2.30、HONO −0.50(design.md §7.4 の主張の中)。
- 費用: VAL7 と対の 35 run で QM 78,750 → 71,731 core 秒(−8.9%、M06-2X の層 9,755 を含む。層を除くと −21%)。41 run で 84,471 core 秒・812 ジョブ。
- 自己一致: `replay.py FINAL --src VAL9 --strict` で 40/40 PASS(misses 0、記録の差 0)。
- 失敗: `h_c2h4` の H···C2H4 錯体の M06-2X の SP が収束せず rc 1(順位の点ではなく、順位は変わらない。残る課題 6)。

## 5. 残る課題

| # | 区分 | 重大度 | 内容 |
|---|---|---|---|
| 1 | U6 | medium | 状態の組の仮説の種は順に試し、最初に結論が出た種で case が閉じる。種の順(低レベルの障壁の順など)と、2 段の結論より直接の素反応を優先するかは未決。S6 では C 上の置換の種が CH4 + HO を経る 2 段を示して閉じ、O 上の置換の素反応(−1418.1i)は試されない |
| 2 | U6 | medium | まとめた仮説の χ は仮説の端の添字で測るので、別の添字の経路の正しい TS でも χ が下がる(O 上の置換の TS で 0.276、自分の端では 0.80)。端を鞍点の添字に合わせて測る方法は、縮退転位で両端が同じ添字に重なるので、そのままでは使えない。下限は下げない |
| 3 | U6 | medium | S5 の CH2OOH• → CH2O + OH(O–O 開裂)の TS が見つからない。棄却した鞍点の χ は 0.271〜0.287 で下限 0.3 に近い。部分的な結合変化の鞍点が増えたら分布を測り直す |
| 4 | U5 | medium | 結論の出ない分割の子は、その極小を層の SP の対象から外す(順位を付けられる反応だけが対象を足す)。S5 では多段の親と split2 が `thermo_unavailable`・`mixed_level_of_theory` になる。多段の親の端を層の対象に残すかは未決 |
| 5 | U5 | medium | GEN-05 の規則は、TS が端の状態とそれより低い別の状態を結ぶときにも中間体を作り、端より低い中間体を経る R → I → P になりうる(S5 で一度起きた。今は χ でその鞍点を棄却している)。規則は変えていない |
| 6 | U8 | medium | 層の SP と freq の電子状態を照合するゲートがなく、層の SCF は freq の vectors から始めない。S5 の BS 錯体の M06-2X の SP は ⟨S²⟩ 0.759 の別の解に収束した(その極小は `spin_contaminated` で順位には影響なし)。VAL9 の `h_c2h4` では H···C2H4 錯体の M06-2X の SP が atomic guess でも cgmin の救済でも収束しなかった(574 s、順位の点ではない) |
| 7 | U8 | low | 手法パネルは会合の BS 錯体の極小にも CCSD(T)(ROHF 参照)の SP をとる。⟨S²⟩ を持たないので `spin_ok` では外れず、その値は BS の点の校正ではない(文書だけで、ゲートはない) |
| 8 | U7 | medium | ほぼ Cs の平らな TS の点群が雑音で Cs と C1 の間で入れ替わり、δG_eff が RT ln 2 だけ動く(S6 の同じ TS で 5.53 と 5.95)。受理の許容は basin の基準のまま |
| 9 | U6 | medium | 一次の鞍点の受理に停留性をかけていない。40 cm⁻¹ 未満の大振幅のモードでは二次モデルが成り立たない(acac の TS で ΔE_N 8.6e-3 Eh)ので、先にこれを扱う必要がある |
| 10 | U5 | low | 会合の分割の子は錯体の極小(`minima[0]`)から始まる。スキャンの井戸から会合の子(単量体 → 井戸)を作るべきかは未決(実例なし)。錯体が DFT で付加体に崩れたときに会合にするのは、宣言した反応だけ |
| 11 | U5 | low | 会合のスキャンの山から REFINE_SADDLE に進む分岐を実計算で確かめるには、DFT の山が 1 kcal/mol を超える会合が要る(H + C2H4 は PBE0 で +0.60) |
| 12 | U3 | low | 費用: DFT の mode-follow は対称な鞍点でも ± 両側を opt と freq にかける(マロンアルデヒドで freq 4 本が 9.1k core 秒の 69%、`sn2_cl_d3h` で像の側が 39%)。QRC の厳密な像の規則(`identity.is_image`)を `relax_to_minimum` の ± にも使える。mode-follow の側は既知の basin に落ちても freq をとる |
| 13 | U6 | low | 費用: 平らな面の QRC の側の opt は maxiter だけで打ち切られる(S5 の側で 100 歩、4,009 s) |
| 14 | U3 | low | 緩和で状態が変わった basin のメンバーを端に使う仮説では、到達した構造でラベルを付け直していない。今の検証セットでは該当 0 |
| 15 | U7 | 潜在 | pymsym の σ の表は回映のない T・O・I で回転対称数の半分(6、12、30)を返し、n > 8 の Cn・Dn・Sn を持たない(KeyError)。GoodVibes と同じ値のために残している。検証セットには出ない |
| 16 | U4 | low | 3 原子以上の直線分子の NT2 の開始構造は、原子ごとの乱数と SVD の軸の向きを通して番号付けに依る。`max_trials_per_source` の境界で同値の類はまとめて落とすので、上限より少ない trial になりうる(該当 0) |
| 17 | U6 | low | string(ZTS)の bead は ⟨S²⟩ を持たず、端の SCF の解を guess にせず、救済の後の通常の SCF も deck にない。端の vectors から解いた開殻の SP は、atomic guess の解より最大 4.7e-6 Eh 高い同じ ⟨S²⟩ の解に落ちることがある(判定への影響なし) |
| 18 | U4 | low | 四重項の CH3·O2 は DFT の極小の段に 30 分以上かかり、高スピンの宣言で AP がかからないことを実計算で示せていない |
| 19 | U9 | low | 古い JobStore の結果は `Evidence.gradient` と `Failure.energy_hartree` を持たず、停留性と再開の上限はそれらに効かない(取り直すと効く)。maxiter の saddle の再開の上限は実計算で発火していない(取り直した saddle は収束した) |
| 20 | U9 | 解決 | WSL だけの品質ゲート: FINAL で prod venv に dev extra(import-linter、pyrefly、radon)を入れ、Linux の pyrefly の 7 件を直した(GoodVibes の Optional は値がなければ例外、libmsym の群が None なら次の設定、SCINE の型、Windows の `ctypes` は `sys.platform` のブロック)。古い ruff(0.15)だけが出すテストの 3 件も直した |
| 21 | U9 | low | `backends/nwchem/engine.py` が 481 行で、500 行の上限に近い。`/dev/shm` は 5.9 GB で、CCSD(T) の GA は 4 × 1,400 MB = 5.6 GB を置くので、それより大きな WFT は `/dev/shm` で先に止まりうる |
| 22 | U4 | low | `gates._bonds_exchanged` は帯のない結合グラフを比べ、`hypotheses` の単量体の判定は帯のない `topology.fragments` を使う(保守的。該当 0) |

実計算で通っていない分岐は [validation.md](../validation.md) §10 にある。

測定の注意: rank 2 の site での `--bind-to none` の修正より前の rank 2 の時間は、QRC の opt で約 5 倍に膨らんでいるので比べていない。S6 の W3 の記録の 790468f505 の結論(2 本目の虚振動 −54.8i での未解決)は、停留点で数え直すと一次の鞍点になる非停留の点の判定で、今のコードの記録にはない。

## 6. 採らなかったもの

| 案 | 理由 |
|---|---|
| χ の下限を下げて S5 の O–O 開裂の鞍点(0.271〜0.287)を通す | 分岐を通すためにしきい値を動かすことになる。分布の空白で決めた値を保つ |
| 2 回目の higher_order_retry | 停留点で数え直すと一次になった(Newton の歩)。根拠なしに予算を足さない |
| 停留点の汎関数を M06-2X か RSH にする | NWChem が解析 Hessian を使わず freq が 1 桁高くつく。SP//PBE0 の構造の誤差は S6 で +0.3 kcal/mol |
| ωB97X-D3 の層 | M06-2X の 6.6 倍の費用で、S6 で誤差が大きい(−2.9) |
| トンネル補正を δG_eff に入れる | ν‡ が PBE0 の誤差を持ち、モデル誤差で序数を並べ替える。文書の目安だけにした(design.md §10) |
| tight の収束、noise_cm1・saddle_cm1 の変更 | 平らな PES でどこに落ちるかは決まらない。停留性は勾配と Hessian で判定する |
| 電荷 ≤ −2 の一律の拒否 | 束縛される大きな多価陰イオンまで拒否する。文書で範囲外とした |
| 自己像の TS の縮退反応に −RT ln 2 を加える | 大きさ(0.41)が層の残差より小さく、外部コードとの突き合わせの検証が要る。入れていない |
