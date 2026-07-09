from __future__ import annotations

"""Portable scheduler script generation for SLURM/PBS/LSF/local review.

The scripts are templates for production review.  They do not submit jobs during
bundle generation.
"""

import re
from pathlib import Path
from typing import Any

import pandas as pd


def _safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _mem_mb(memory_gb: Any) -> int:
    try:
        return int(float(memory_gb) * 1024)
    except Exception:
        return 4096


def _script_for_row(row: dict[str, Any], scheduler: str, run_id: str, project_root: str) -> str:
    stage = str(row.get("stage", "stage"))
    job_name = _safe(f"hfauto_{run_id}_{stage}")
    cmd = str(row.get("command", "echo missing command")).replace("${RUN_ID}", run_id)
    ncores = int(row.get("ncores") or 1)
    mem_gb = int(float(row.get("memory_gb") or 4))
    queue = row.get("queue") or row.get("partition") or "short"
    wall = row.get("time_slurm") or row.get("walltime") or "01:00:00"
    prolog = ["#!/usr/bin/env bash"]
    if scheduler == "slurm":
        prolog += [
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --cpus-per-task={ncores}",
            f"#SBATCH --mem={mem_gb}G",
            f"#SBATCH --time={wall}",
            f"#SBATCH --partition={queue}",
            "#SBATCH --output=logs/%x-%j.out",
            "#SBATCH --error=logs/%x-%j.err",
        ]
    elif scheduler == "pbs":
        prolog += [
            f"#PBS -N {job_name}",
            f"#PBS -l select=1:ncpus={ncores}:mem={mem_gb}gb",
            f"#PBS -l walltime={wall if '-' not in str(wall) else '24:00:00'}",
            f"#PBS -q {queue}",
            "#PBS -j oe",
        ]
    elif scheduler == "lsf":
        prolog += [
            f"#BSUB -J {job_name}",
            f"#BSUB -n {ncores}",
            f"#BSUB -M {_mem_mb(mem_gb)}",
            f"#BSUB -q {queue}",
            "#BSUB -oo logs/%J.out",
            "#BSUB -eo logs/%J.err",
        ]
    prolog += [
        "set -euo pipefail",
        "mkdir -p logs",
        f"export PYTHONPATH={project_root}:${{PYTHONPATH:-}}",
        cmd,
        "",
    ]
    return "\n".join(prolog)


def render_scheduler_bundle(
    plan: pd.DataFrame,
    out_dir: str | Path,
    scheduler: str = "slurm",
    run_id: str = "run",
    project_root: str = "$PWD",
    prefix: str = "submit",
) -> dict[str, Path]:
    out = Path(out_dir) / scheduler
    out.mkdir(parents=True, exist_ok=True)
    scripts: list[Path] = []
    for i, row in plan.iterrows():
        stage = str(row.get("stage", "stage"))
        path = out / f"{prefix}_{i:04d}_{_safe(stage)}.{scheduler}.sh"
        path.write_text(_script_for_row(row.to_dict(), scheduler, run_id, project_root), encoding="utf-8")
        path.chmod(0o755)
        scripts.append(path)

    submit_all = out / "submit_all.sh"
    if scheduler == "slurm":
        submit_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", "mkdir -p logs", ""]
        prev = ""
        for i, path in enumerate(scripts):
            var = f"jid_{i:04d}"
            dep = f" --dependency=afterok:${{{prev}}}" if prev else ""
            submit_lines.append(f"{var}=$(sbatch{dep} --parsable {path.name})")
            submit_lines.append(f"echo submitted {path.name} ${{{var}}}")
            prev = var
    elif scheduler == "pbs":
        submit_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", ""]
        for path in scripts:
            submit_lines.append(f"qsub {path.name}")
    elif scheduler == "lsf":
        submit_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", ""]
        for path in scripts:
            submit_lines.append(f"bsub < {path.name}")
    else:
        submit_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", ""] + [f"bash {p.name}" for p in scripts]
    submit_all.write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
    submit_all.chmod(0o755)

    cancel = out / "cancel_jobs.sh"
    if scheduler == "slurm":
        cancel.write_text("#!/usr/bin/env bash\nset -euo pipefail\nsqueue -u $USER -h -o '%i %j' | awk '/hfauto_/ {print $1}' | xargs -r scancel\n", encoding="utf-8")
    elif scheduler == "pbs":
        cancel.write_text("#!/usr/bin/env bash\nset -euo pipefail\nqstat -u $USER | awk '/hfauto_/ {print $1}' | xargs -r qdel\n", encoding="utf-8")
    elif scheduler == "lsf":
        cancel.write_text("#!/usr/bin/env bash\nset -euo pipefail\nbjobs | awk '/hfauto_/ {print $1}' | xargs -r bkill\n", encoding="utf-8")
    else:
        cancel.write_text("#!/usr/bin/env bash\nset -euo pipefail\necho local scheduler has no central cancel operation\n", encoding="utf-8")
    cancel.chmod(0o755)

    status = out / "status_snapshot_command.sh"
    if scheduler == "slurm":
        status.write_text("#!/usr/bin/env bash\nset -euo pipefail\nsqueue -u $USER -o '%i|%j|%T|%M|%D|%R'\n", encoding="utf-8")
    elif scheduler == "pbs":
        status.write_text("#!/usr/bin/env bash\nset -euo pipefail\nqstat -u $USER\n", encoding="utf-8")
    elif scheduler == "lsf":
        status.write_text("#!/usr/bin/env bash\nset -euo pipefail\nbjobs\n", encoding="utf-8")
    else:
        status.write_text("#!/usr/bin/env bash\nset -euo pipefail\necho local planned\n", encoding="utf-8")
    status.chmod(0o755)
    return {"scripts_dir": out, "submit_all": submit_all, "cancel_jobs": cancel, "status_command": status}
