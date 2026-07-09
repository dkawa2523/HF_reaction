from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_ops.schedulers.base import SchedulerAdapter


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


class SlurmAdapter(SchedulerAdapter):
    name = "slurm"

    def submit_command(self, script: str | Path, dependency: str | None = None) -> str:
        dep = f" --dependency=afterok:{dependency}" if dependency else ""
        return f"sbatch{dep} --parsable {script}"

    def status_command(self, job_id: str | None = None, user: str | None = None) -> list[str]:
        if job_id:
            return ["squeue", "-j", str(job_id), "-o", "%i %T %M %R %j"]
        if user:
            return ["squeue", "-u", str(user), "-o", "%i %T %M %R %j"]
        return ["squeue", "-u", "$USER", "-o", "%i %T %M %R %j"]

    def cancel_command(self, job_id: str) -> list[str]:
        return ["scancel", str(job_id)]

    def list_user_jobs_command(self, user: str = "$USER") -> str:
        return f"squeue -u {user} -o '%i %T %M %R %j'"

    def script_header(self, row: dict[str, Any], job_name: str) -> str:
        queue = row.get("queue") or row.get("partition") or "short"
        return f"""#!/usr/bin/env bash
#SBATCH --job-name={_safe_name(job_name)}
#SBATCH --cpus-per-task={int(float(row.get('ncores', 1) or 1))}
#SBATCH --mem={int(float(row.get('memory_gb', 4) or 4))}G
#SBATCH --time={row.get('time_slurm') or row.get('walltime') or '01:00:00'}
#SBATCH --partition={queue}
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
"""

    def array_header(self, row: dict[str, Any], job_name: str, n_tasks: int) -> str:
        header = self.script_header(row, job_name)
        return header + f"#SBATCH --array=1-{max(1, int(n_tasks))}\n"

    def task_index_expr(self) -> str:
        return "${SLURM_ARRAY_TASK_ID:-1}"

    def render_stage_script(self, row: dict[str, Any], path: Path, run_id: str, project_root: str = "$PWD") -> Path:
        stage = str(row.get("stage", "stage"))
        command = str(row.get("command", "echo missing command")).replace("${RUN_ID}", run_id)
        content = self.script_header(row, f"hfauto_{run_id}_{stage}") + f"""
set -euo pipefail
mkdir -p logs
export PYTHONPATH={project_root}:${{PYTHONPATH:-}}
{command}
"""
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755)
        return path


def render_slurm_stage_script(row: dict[str, Any], out_dir: Path, run_id: str, project_root: str = "$PWD") -> Path:
    stage = str(row["stage"])
    path = out_dir / f"submit_{int(row['stage_index']):02d}_{stage}.slurm"
    return SlurmAdapter().render_stage_script(row, path, run_id=run_id, project_root=project_root)


def render_slurm_bundle(plan: pd.DataFrame, out_dir: str | Path, run_id: str, project_root: str = "$PWD") -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scripts_dir = out / "slurm"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    adapter = SlurmAdapter()
    paths: list[Path] = []
    for _, row in plan.iterrows():
        paths.append(render_slurm_stage_script(row.to_dict(), scripts_dir, run_id=run_id, project_root=project_root))

    submit_all = scripts_dir / "submit_all.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")\"", "mkdir -p logs", ""]
    previous_var = ""
    for path, (_, row) in zip(paths, plan.iterrows()):
        var = f"jid_{int(row['stage_index']):02d}"
        rel = path.name
        if previous_var:
            lines.append(f"{var}=$({adapter.submit_command(rel, dependency='${' + previous_var + '}')})")
        else:
            lines.append(f"{var}=$({adapter.submit_command(rel)})")
        lines.append(f"echo submitted {row['stage']} ${{{var}}}")
        previous_var = var
    lines.append("")
    submit_all.write_text("\n".join(lines), encoding="utf-8")
    submit_all.chmod(0o755)

    local_dry_run = scripts_dir / "run_local_dry_run.sh"
    dry = ["#!/usr/bin/env bash", "set -euo pipefail", "cd \"$(dirname \"$0\")/../..\"", ""]
    for _, row in plan.iterrows():
        dry.append(str(row["command"]).replace("${RUN_ID}", run_id))
    local_dry_run.write_text("\n".join(dry), encoding="utf-8")
    local_dry_run.chmod(0o755)
    return {"scripts_dir": scripts_dir, "submit_all": submit_all, "local_dry_run": local_dry_run}


def render_slurm_artifact_array(
    artifact_job_plan_csv: str | Path,
    out_dir: str | Path,
    run_id: str,
    project_root: str = "$PWD",
    chunk_size: int = 1,
) -> dict[str, Path]:
    """Render a conservative SLURM array wrapper for artifact-level jobs."""
    import math
    out = Path(out_dir)
    scripts_dir = out / "slurm"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    plan = pd.read_csv(artifact_job_plan_csv) if Path(artifact_job_plan_csv).exists() else pd.DataFrame()
    items = scripts_dir / "artifact_job_items.tsv"
    if plan.empty:
        items.write_text("array_index\tjob_id\tcommand\n", encoding="utf-8")
        n_tasks = 0
        ncores = 1
        mem = 4
        time_slurm = "01:00:00"
        queue = "short"
    else:
        cols = ["array_index", "job_id", "command"]
        plan[cols].to_csv(items, index=False, sep="\t")
        n_tasks = int(math.ceil(len(plan) / max(1, int(chunk_size))))
        ncores = int(plan.get("ncores", pd.Series([1])).max())
        mem = int(plan.get("memory_gb", pd.Series([4])).max())
        time_slurm = str(plan.get("time_slurm", pd.Series(["01:00:00"])).iloc[0])
        queue = str(plan.get("queue", pd.Series(["short"])).iloc[0])
    array_script = scripts_dir / "artifact_array.slurm"
    upper = max(0, n_tasks - 1)
    content = f"""#!/usr/bin/env bash
#SBATCH --job-name=hfauto_{run_id}_artifact_array
#SBATCH --array=0-{upper}
#SBATCH --cpus-per-task={ncores}
#SBATCH --mem={mem}G
#SBATCH --time={time_slurm}
#SBATCH --partition={queue}
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
mkdir -p logs
export PYTHONPATH={project_root}:${{PYTHONPATH:-}}
ITEM_FILE=\"$(dirname \"$0\")/artifact_job_items.tsv\"
python -m hfauto_ops.cli.main run-artifact-job --items "$ITEM_FILE" --array-index "${{SLURM_ARRAY_TASK_ID:-0}}" --chunk-size {int(chunk_size)}
"""
    array_script.write_text(content, encoding="utf-8")
    array_script.chmod(0o755)
    submit = scripts_dir / "submit_artifact_array.sh"
    submit.write_text("#!/usr/bin/env bash\nset -euo pipefail\nsbatch artifact_array.slurm\n", encoding="utf-8")
    submit.chmod(0o755)
    return {"artifact_items_tsv": items, "artifact_array_script": array_script, "artifact_array_submit": submit}
