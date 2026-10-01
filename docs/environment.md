# 実行環境

本番の計算は WSL2 の Ubuntu 24.04(Python 3.12)で行う。Windows の QA 用 venv は編集と静的検査に使う。

## 1. 版数

外部プログラムの版数は site ファイル(`configs/sites/wsl_local.yaml`)の `version` に pin する。pin は JobStore の鍵に入り、ジョブごとに出力から観測した版数と照合する(違えば `method_mismatch`)。

| 部品 | 版数 | 入れ方(WSL) | 確認 |
|---|---|---|---|
| NWChem | 7.2.3(conda-forge の MPI ビルド、Open MPI 5) | micromamba で `/home/user/.local/opt/nwchem-7.2.3` に入れ、`/home/user/.local/bin/nwchem` からリンク | 出力の版数行(ジョブごと) |
| xTB | 6.7.1 | 公式バイナリを `/home/user/.local/opt/xtb-6.7.1` に置く | `xtb --version`(doctor) |
| CREST | 3.0.2 | 公式バイナリを `/home/user/.local/opt/crest-3.0.2` に置く | `crest --version`(doctor) |
| pysisyphus | 1.0 | production extra。xTB ネイティブ計算器の CI-NEB(`pysis_neb`)だけに使う | worker の import(doctor)。pin は xTB の版数 |
| SCINE ReaDuct | 6.1.0(scine-xtb-wrapper 3.0.2) | production extra | worker の import(doctor) |
| GoodVibes | 4.3.0 | production extra。thermo stage が同じプロセスで呼ぶ | `hfauto run` の preflight が thermo を含む pipeline で照合する |
| pymsym | 0.3.5(対称数 σ の libmsym) | production extra で `pymsym==0.3.5` に固定。入れ直さない(別の libmsym のビルドは σ を変えうる) | `importlib.metadata.version('pymsym')`(WSL の既定テストが照合) |

production extra(`pyproject.toml`)は Linux・CPython 3.12 のときだけ入り、rdkit、numba、llvmlite と、scine_utilities が宣言せずに使う setuptools を含む。本体の依存は pydantic、typer、PyYAML、rich、numpy、scipy だけである。外部プログラムを conda でそろえる場合は `environment.production.yml` を使う。

## 2. インストール(WSL)

リポジトリの直下で実行する(CLI の設定名は `configs/` からの相対で解決される)。

```bash
uv venv --python 3.12 /home/user/.venvs/hfauto-prod
uv pip install --python /home/user/.venvs/hfauto-prod/bin/python -e ".[production,dev]"
/home/user/.venvs/hfauto-prod/bin/hfauto doctor --site configs/sites/wsl_local.yaml
```

- リポジトリは Windows 側にあり、WSL からは `/mnt/c/Users/user/Desktop/HF_reaction/hfauto_final_baseline` に見える。editable で入れるので、コードの変更はそのまま WSL に反映される。
- `hfauto doctor` は site の全エンジンの実行ファイル、版数の pin、worker の Python モジュール、scratch の場所を検査する。`hfauto run` も実行前に同じ検査(preflight)を行い、thermo を含む pipeline では GoodVibes の版数も照合する。

## 3. WSL の資源と site の値

| 項目 | 設定 |
|---|---|
| `.wslconfig` | `C:\Users\user\.wslconfig` に `processors=4`、`memory=12GB`、`swap=2GB`(WSL から 4 vCPU・約 11.9 GB)。site は `cores: 4` |
| scratch | ext4 上の `/home/user/hfauto_scratch`(NWChem は `/home/user/hfauto_scratch/nwchem`)。`/mnt/<ドライブ>` 上の scratch は preflight が拒否する。`/tmp` も使わない |
| run ディレクトリ | `--run-dir` で ext4 上に置く(例 `/home/user/hfauto_runs/<run>`)。省略するとリポジトリ内の `runs/<system_id>_<pipeline_id>` になる |
| NWChem | 4 rank × `memory_mb_per_rank` 2,000 MB、`timeout_s` 14,400。ranks × メモリの合計は MemTotal の 0.8 倍以下にする(4 × 2,000 MB = 8 GB)。`OMP_NUM_THREADS=1`(MPI rank で並列化)。反応ケースの中で同時に走る SP・QRC の側は rank を分け合い(cores // 同時の数)、timeout を同じ倍率で延ばす。mpirun には常に `--bind-to none` を付ける。CCSD(T) の deck は GA の global にメモリの 70%(4 × 1,400 MB = 5.6 GB)を置き、これは `/dev/shm`(WSL で 5.9 GB)に載るので、それより大きな WFT のジョブは `/dev/shm` で先に止まりうる |
| xTB・CREST・ReaDuct | `OMP_NUM_THREADS=<threads>,1`、`OMP_STACKSIZE=4G`(アダプタが設定する)。CREST の `-T` は site の `engines.crest.execution.threads`(4)から決まる |
| コア数 | 実行中のジョブの ranks × threads の合計を `cores` 以下に保つ(`JobRunner` のセマフォ)。極小の stage は全 rank で直列 |
| 系のサイズ | 4 コア・def2-SVPD では解析 freq が約 35 原子でジョブの timeout(4 h)に達し、途中継続がないので全損する。大きな系は ranks を増やす |
| シグナル | `hfauto run` は SIGINT・SIGTERM・SIGHUP で外部プログラム(NWChem の MPI ランクを含むプロセスグループ)を止め、128 + signum で終わる。SIGKILL は捕まえられないので、`timeout` は既定の SIGTERM で使い(`-s KILL` にしない)、スケジューラの SIGTERM から SIGKILL までの猶予(Slurm の KillWait)は数秒以上にする |
| SiteLock | `hfauto run` は実行中ずっと `<scratch_root>/.hfauto_site.lock` を保持し、別の run は保持者(pid、run ディレクトリ、時刻)を示して直ちに失敗する。pid が死んだロックは奪う。run を止めるシグナルは hfauto の python プロセスに送る(`time` などの親だけに送ると python が残り、ロックを持ち続ける) |

`.wslconfig` を 16 vCPU / 48 GB に広げたときの site の値(cores 16、NWChem 16 rank × 2,000 MB、CREST threads 4)は `configs/sites/wsl_local.yaml` のコメントにある。

## 4. テストと品質ゲート

品質ゲートの合否は WSL の prod venv(dev extra を入れたもの、§2)の結果だけを正とする(Windows のアプリケーション制御は範囲外)。ゲートは、pytest が通ること、ruff の指摘 0、import-linter の 5 契約、pyrefly の指摘 0、`radon cc -n D hfauto` が何も出さないこと(CC 21 以上がない)、`hfauto/` と `tests/` に 500 行を超える `.py` がないことの 6 つである。

既定の `pytest` は `-m 'not real'` で、unit、golden(実出力の抜粋)、integration(fake エンジンで stage をつなぐ)だけが動く。GoodVibes の golden は GoodVibes がなければ、SMILES のテストは rdkit がなければ飛ばされる。

```bash
# WSL(Git Bash から呼ぶときは、/home で始まる引数を書き換えられないように MSYS_NO_PATHCONV=1 を付ける)
MSYS_NO_PATHCONV=1 wsl -e bash -lc 'cd /mnt/c/Users/user/Desktop/HF_reaction/hfauto_final_baseline && P=/home/user/.venvs/hfauto-prod/bin && $P/python -m pytest -q -p no:cacheprovider && $P/ruff check . && $P/lint-imports && $P/pyrefly check --python-interpreter-path $P/python && $P/radon cc -n D hfauto'
```

- pyrefly の設定(`pyproject.toml` の `python-interpreter-path`)は Windows の venv を指すので、WSL では prod venv を渡す。Linux では production extra(GoodVibes、pymsym、SCINE)の型も読むので、Windows の venv では出ない指摘が出うる。Windows 専用の `ctypes` の呼び出しは `if sys.platform == "win32":` のブロックに置く。
- Windows の `.venv-win-qa`(dev と chem extra)は編集中の手早い確認用で、合否には使わない。GoodVibes(pymsym)が入らず、この PC ではアプリケーション制御が RDKit の DLL を止め、`pyrefly.exe` を止めたこともある。

### 実エンジンの smoke(`pytest -m real`、WSL)

`tests/smoke/` のテストは marker `real` が付き、`HFAUTO_REAL=1` と `HFAUTO_SITE=<site ファイル>` がそろったときだけ動く。各エンジンを小さな系で 1 回ずつ動かし(NWChem の opt・freq・saddle・ECP・CCSD(T)・ZTS、xTB と `pysis_neb`、CREST、ReaDuct の NT2、GoodVibes)、約 1 分かかる。`--basetemp` は実行のたびに消されるので、専用の ext4 のディレクトリを指定する。

```bash
MSYS_NO_PATHCONV=1 wsl -e bash -lc 'cd /mnt/c/Users/user/Desktop/HF_reaction/hfauto_final_baseline && HFAUTO_REAL=1 HFAUTO_SITE=configs/sites/wsl_local.yaml /home/user/.venvs/hfauto-prod/bin/python -m pytest -m real tests/smoke --basetemp=/home/user/hfauto_smoke'
```

実計算の run(`hfauto run`)も同じ venv と site ファイルで実行する。
