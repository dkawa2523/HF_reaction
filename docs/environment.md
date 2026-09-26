# 実行環境

本番の計算は WSL2 の Ubuntu 24.04(Python 3.12)で行う。Windows では QA 用の venv で fake エンジンと golden のテストだけを動かす。

## 1. 版数の固定

外部プログラムの版数は site ファイル(`configs/sites/wsl_local.yaml`)の `version` に pin する。pin は JobStore のキーに入り、実行後に出力から観測した版数と照合する(違えば `method_mismatch`)。

| 部品 | 版数 | 入れ方(WSL) | 確認 |
|---|---|---|---|
| NWChem | 7.2.3(conda-forge の MPI ビルド、Open MPI 5.0.10) | micromamba で `/home/user/.local/opt/nwchem-7.2.3` に入れ、`/home/user/.local/bin/nwchem` からリンク | 版数のフラグがないので、ジョブごとに出力の版数行を照合する |
| xTB | 6.7.x(6.7.1) | 公式バイナリを `/home/user/.local/opt/xtb-6.7.1` に置く | `xtb --version`(doctor) |
| CREST | 3.0.x(3.0.2) | 公式バイナリを `/home/user/.local/opt/crest-3.0.2` に置く | `crest --version`(doctor) |
| pysisyphus | 1.0(1.0.0) | production extra(`pysisyphus>=1.0,<2`)。xTB ネイティブ計算器でだけ使う | worker の import(doctor)。site の pin は xTB の版数 |
| SCINE ReaDuct | 6.1.0(scine-xtb-wrapper 3.0.2) | production extra | worker の import(doctor) |
| GoodVibes | 4.3.0 | production extra。API(`goodvibes.api.compute_thermo`)を worker で呼ぶ | `goodvibes.__version__`(doctor) |
| pymsym | 0.3.5 | production extra。GoodVibes の対称数に使う | worker の import(doctor) |

production extra(`pyproject.toml`)は Linux / CPython 3.12 のときだけ入る。ほかに rdkit、numba、llvmlite と、scine_utilities が宣言せずに使う setuptools を含む。本体の依存は pydantic、typer、PyYAML、rich、numpy、scipy だけである。

## 2. インストール(WSL)

リポジトリの直下で実行する(CLI の設定名は `configs/` からの相対で解決される)。

```bash
uv venv --python 3.12 /home/user/.venvs/hfauto-prod
uv pip install --python /home/user/.venvs/hfauto-prod/bin/python -e ".[production,dev]"
/home/user/.venvs/hfauto-prod/bin/hfauto doctor --site configs/sites/wsl_local.yaml
```

- リポジトリは Windows 側にあり、WSL からは `/mnt/c/Users/user/Desktop/HF_reaction/hfauto_final_baseline` に見える。editable で入れるので、コードの変更はそのまま WSL に反映される。
- 外部プログラムを conda でそろえる場合は `environment.production.yml` を使う(版数は `hfauto doctor` で pin と照合する)。
- `hfauto doctor` は site の全エンジンについて、実行ファイル、版数の pin、worker 側の Python モジュール、scratch の場所を検査し、問題があれば終了コード 1 を返す。`hfauto run` も実行前に同じ検査(preflight)を行う。

## 3. WSL の資源

| 項目 | 設定 |
|---|---|
| `.wslconfig` | `C:\Users\user\.wslconfig` に `processors=4`、`memory=12GB`、`swap=2GB`。WSL から 4 vCPU・約 11.9 GB に見える。site は `cores: 4` |
| scratch | ext4 上の `/home/user/hfauto_scratch`(NWChem は `/home/user/hfauto_scratch/nwchem`)。`/mnt/c` などの Windows ドライブ上の scratch は preflight が拒否する。`/tmp` も使わない |
| run ディレクトリ | `--run-dir` で ext4 上に置く(例: `/home/user/hfauto_v2/<run>`)。省略するとリポジトリ内の `runs/<system_id>_<pipeline_id>` になり、過去の run と混ざる |
| スレッド | NWChem は `OMP_NUM_THREADS=1`(MPI rank で並列化)、xTB・CREST・ReaDuct は `OMP_NUM_THREADS=<n>,1` と `OMP_STACKSIZE=4G`。アダプタが設定する。CREST のスレッド数(`-T` と OMP)は site の `engines.crest.execution.threads` だけから決まる(省くと既定の 1 で `-T 1`) |
| メモリ | 1 rank あたり 1,000〜1,500 MB、ranks × memory の合計を MemTotal の 0.8 倍以下にする。`configs/sites/wsl_local.yaml` の NWChem は 4 rank × 1,200 MB、timeout 14,400 s |
| コア数 | `JobRunner` のセマフォで、実行中のジョブの ranks × threads の合計を `cores` 以下に保つ。極小(xTB と DFT)と反応ケースは直列、CREST と GoodVibes だけが並列 |
| SiteLock | `hfauto run` は実行中ずっと `<scratch_root>/.hfauto_site.lock`(ここでは /home/user/hfauto_scratch の下)を保持する。別の run は保持者(pid、run ディレクトリ、時刻)を示して直ちに失敗する。pid が死んでいるロックは奪う |

`.wslconfig` を 16 vCPU / 48 GB に広げた場合の site の値(cores 16、NWChem 16 rank × 2,000 MB など)は `configs/sites/wsl_local.yaml` のコメントにある。

## 4. テスト

### Windows の QA 用 venv

`.venv-win-qa`(Python 3.12)には dev と chem(rdkit)の extra だけを入れ、production extra は入らない(Linux 専用)。既定の `pytest` は `-m 'not real'` で、unit、golden(実出力の抜粋)、integration(fake エンジンで stage をつなぐ)だけが動く。pysisyphus、SCINE、GoodVibes と外部プログラムを要するのは `tests/smoke/` だけである。rdkit がない環境では SMILES のテストが `importorskip` で飛ばされる。

```bash
.venv-win-qa/Scripts/python.exe -m pytest -q -p no:cacheprovider
.venv-win-qa/Scripts/ruff check hfauto tests
.venv-win-qa/Scripts/lint-imports
.venv-win-qa/Scripts/pyrefly check
```

### 実エンジンの smoke(`pytest -m real`、WSL)

`tests/smoke/` のテストは marker `real` が付き、`HFAUTO_REAL=1` と `HFAUTO_SITE=<site ファイル>` がそろったときだけ動く(そろわなければ skip)。ジョブは pytest の一時ディレクトリの下に作られるので、`--basetemp` で ext4 上の `/home/user/hfauto_v2` の下を指定する。W7 では 10 passed、59 秒だった。その後、HCN → HNC の ZTS(PBE0/STO-3G、9 beads、1 rank で約 80 s、4 rank で約 15 s。bead エネルギーが落ち着いて single_max になること)と、F⁻·(HF)₂ の組成の CREST(`--noopt` で初期トポロジー検査を通ること、約 1.4 s。`--noopt` なしでは初期トポロジー検査で止まり約 25 s)を加えた。Git Bash から WSL を呼ぶときは、/home で始まる引数が Windows のパスに書き換えられないように `MSYS_NO_PATHCONV=1` を付ける。

```bash
MSYS_NO_PATHCONV=1 wsl -e bash -lc 'cd /mnt/c/Users/user/Desktop/HF_reaction/hfauto_final_baseline && HFAUTO_REAL=1 HFAUTO_SITE=configs/sites/wsl_local.yaml /home/user/.venvs/hfauto-prod/bin/python -m pytest -m real tests/smoke --basetemp=/home/user/hfauto_v2/smoke'
```

WSL のシェルの中では `MSYS_NO_PATHCONV` は要らない。実計算の run(`hfauto run`)も同じ venv と site ファイルで実行する。
