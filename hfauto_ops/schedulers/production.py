from __future__ import annotations

"""Scheduler script renderers for production-oriented operations bundles.

The generated scripts are intentionally reviewable shell scripts. They do not
submit jobs unless a user explicitly executes the generated submit script.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, write_json


def _safe(value: str) -> str:
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in value)[:120]


def _wall(row: dict[str, Any]) -> str:
    return str(row.get("time_slurm") or row.get("walltime") or "01:00:00")


def render_script(row: dict[str, Any], scheduler: str, project_root: str = "$PWD") -> str:
    job = _safe(str(row.get("array_item_id") or row.get("job_id") or row.get("stage", "hfauto")))
    ncores = int(row.get("ncores") or 1)
    mem = int(float(row.get("memory_gb") or 4))
    queue = str(row.get("queue") or "")
    command = str(row.get("command") or row.get("worker_command") or "echo no command")
    scheduler = scheduler.lower()
    if scheduler == "pbs":
        header = [
            "#!/usr/bin/env bash",
            f"#PBS -N {job}",
            f"#PBS -l select=1:ncpus={ncores}:mem={mem}gb",
            f"#PBS -l walltime={_wall(row)}",
        ]
        if queue:
            header.append(f"#PBS -q {queue}")
        header += ["set -euo pipefail", "cd \"$PBS_O_WORKDIR\""]
    elif scheduler == "lsf":
        header = [
            "#!/usr/bin/env bash",
            f"#BSUB -J {job}",
            f"#BSUB -n {ncores}",
            f"#BSUB -M {mem * 1024}",
            f"#BSUB -W {_wall(row)}",
        ]
        if queue:
            header.append(f"#BSUB -q {queue}")
        header += ["set -euo pipefail"]
    elif scheduler == "local":
        header = ["#!/usr/bin/env bash", "set -euo pipefail"]
    else:
        header = [
            "#!/usr/bin/env bash",
            f"#SBATCH --job-name={job}",
            f"#SBATCH --cpus-per-task={ncores}",
            f"#SBATCH --mem={mem}G",
            f"#SBATCH --time={_wall(row)}",
        ]
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        header += ["#SBATCH --output=logs/%x-%j.out", "#SBATCH --error=logs/%x-%j.err", "set -euo pipefail", "mkdir -p logs"]
    body = [f"export PYTHONPATH={project_root}:${{PYTHONPATH:-}}", "echo '[hfauto] start ' $(date)", command, "echo '[hfauto] end ' $(date)", ""]
    return "\n".join(header + body)


def render_array_script(array_items: str | Path, out_dir: str | Path, scheduler: str = "slurm", project_root: str = "$PWD", job_name: str = "hfauto_array") -> Path:
    out = ensure_dir(out_dir)
    items_path = Path(array_items)
    try:
        n_items = sum(1 for line in items_path.read_text(encoding="utf-8").splitlines() if line.strip())
    except Exception:
        n_items = 0
    max_idx = max(0, n_items - 1)
    scheduler = scheduler.lower()
    script = out / f"submit_array_{scheduler}.sh"
    if scheduler == "pbs":
        header = ["#!/usr/bin/env bash", f"#PBS -N {_safe(job_name)}", f"#PBS -J 0-{max_idx}", "#PBS -l select=1:ncpus=1:mem=4gb", "#PBS -l walltime=02:00:00", "set -euo pipefail", "cd \"$PBS_O_WORKDIR\"", "TASK_ID=${PBS_ARRAY_INDEX:-0}"]
    elif scheduler == "lsf":
        header = ["#!/usr/bin/env bash", f"#BSUB -J {_safe(job_name)}[0-{max_idx}]", "#BSUB -n 1", "#BSUB -M 4096", "#BSUB -W 02:00", "set -euo pipefail", "TASK_ID=${LSB_JOBINDEX:-0}"]
    elif scheduler == "local":
        header = ["#!/usr/bin/env bash", "set -euo pipefail", "TASK_ID=${1:-0}"]
    else:
        header = ["#!/usr/bin/env bash", f"#SBATCH --job-name={_safe(job_name)}", f"#SBATCH --array=0-{max_idx}", "#SBATCH --cpus-per-task=1", "#SBATCH --mem=4G", "#SBATCH --time=02:00:00", "#SBATCH --output=logs/%x-%A_%a.out", "#SBATCH --error=logs/%x-%A_%a.err", "set -euo pipefail", "mkdir -p logs", "TASK_ID=${SLURM_ARRAY_TASK_ID:-0}"]
    body = [f"export PYTHONPATH={project_root}:${{PYTHONPATH:-}}", f"python -m hfauto_ops.cli.main run-array-task --array-items {items_path} --task-id $TASK_ID"]
    script.write_text("\n".join(header + body) + "\n", encoding="utf-8")
    script.chmod(0o755)
    return script


def write_submit_status_cancel_scripts(out_dir: str | Path, scheduler: str = "slurm") -> dict[str, Path]:
    out = ensure_dir(out_dir)
    scheduler = scheduler.lower()
    if scheduler == "pbs":
        submit = "qsub $1"
        status = "qstat ${1:-}"
        cancel = "qdel $1"
    elif scheduler == "lsf":
        submit = "bsub < $1"
        status = "bjobs ${1:-}"
        cancel = "bkill $1"
    elif scheduler == "local":
        submit = "bash $1"
        status = "echo local scheduler has no status API"
        cancel = "echo local scheduler has no cancel API"
    else:
        submit = "sbatch $1"
        status = "squeue ${1:+-j $1}"
        cancel = "scancel $1"
    paths = {}
    for name, cmd in {"submit_one.sh": submit, "status.sh": status, "cancel.sh": cancel}.items():
        p = out / name
        p.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + cmd + "\n", encoding="utf-8")
        p.chmod(0o755)
        paths[name.replace(".sh", "")] = p
    return paths


def write_multi_scheduler_bundle(
    resource_plan_csv: str | Path | None,
    array_items: str | Path | None,
    out_dir: str | Path,
    schedulers: list[str] | None = None,
    project_root: str = "$PWD",
) -> dict[str, Path]:
    schedulers = schedulers or ["slurm", "pbs", "lsf", "local"]
    out = ensure_dir(out_dir)
    paths: dict[str, Path] = {}
    for sched in schedulers:
        sched_dir = ensure_dir(out / sched)
        if resource_plan_csv and Path(resource_plan_csv).exists():
            df = pd.read_csv(resource_plan_csv)
            scripts = []
            for _, row in df.iterrows():
                script = sched_dir / f"submit_{int(row.get('stage_index', len(scripts))):02d}_{row.get('stage')}.{sched}.sh"
                script.write_text(render_script(row.to_dict(), scheduler=sched, project_root=project_root), encoding="utf-8")
                script.chmod(0o755)
                scripts.append(script)
            submit_all = sched_dir / "submit_all.sh"
            if sched == "pbs":
                lines = [f"qsub {s.name}" for s in scripts]
            elif sched == "lsf":
                lines = [f"bsub < {s.name}" for s in scripts]
            elif sched == "local":
                lines = [f"bash {s.name}" for s in scripts]
            else:
                lines = [f"sbatch {s.name}" for s in scripts]
            submit_all.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(lines) + "\n", encoding="utf-8")
            submit_all.chmod(0o755)
            paths[f"{sched}_stage_submit_all"] = submit_all
        if array_items and Path(array_items).exists():
            paths[f"{sched}_array_submit"] = render_array_script(array_items, sched_dir, scheduler=sched, project_root=project_root)
        control = write_submit_status_cancel_scripts(sched_dir, scheduler=sched)
        paths.update({f"{sched}_{k}": v for k, v in control.items()})
    write_json(out / "scheduler_bundle_manifest.json", {k: str(v) for k, v in paths.items()})
    paths["scheduler_manifest"] = out / "scheduler_bundle_manifest.json"
    return paths
