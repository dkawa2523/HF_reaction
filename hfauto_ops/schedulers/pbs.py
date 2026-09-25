from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir
from hfauto_ops.schedulers.base import SchedulerAdapter


class PBSAdapter(SchedulerAdapter):
    name = "pbs"
    submit_executable = "qsub"
    status_executable = "qstat"
    cancel_executable = "qdel"

    def submit_command(self, script: str | Path) -> list[str]:
        return [self.submit_executable, str(script)]

    def status_command(self, job_id: str | None = None, user: str | None = None) -> list[str]:
        cmd = [self.status_executable]
        if job_id:
            cmd.append(str(job_id))
        elif user:
            cmd += ["-u", str(user)]
        return cmd

    def cancel_command(self, job_id: str) -> list[str]:
        return [self.cancel_executable, str(job_id)]

    def parse_submit_job_id(self, stdout: str) -> str | None:
        return stdout.strip().splitlines()[-1].strip() if stdout.strip() else None


def render_pbs_stage_script(row: dict[str, Any], out_dir: Path, run_id: str, project_root: str = "$PWD") -> Path:
    stage = str(row["stage"])
    path = out_dir / f"submit_{int(row['stage_index']):02d}_{stage}.pbs.sh"
    command = str(row["command"]).replace("${RUN_ID}", run_id)
    queue = row.get("queue")
    queue_directive = f"#PBS -q {queue}\n" if queue else ""
    content = f"""#!/usr/bin/env bash
#PBS -N hfauto_{run_id}_{stage}
#PBS -l select=1:ncpus={int(row['ncores'])}:mem={int(row['memory_gb'])}gb
#PBS -l walltime={row.get('time_slurm', '01:00:00').replace('-', ':')}
{queue_directive}set -euo pipefail
cd "$PBS_O_WORKDIR"
export PYTHONPATH={project_root}:${{PYTHONPATH:-}}
{command}
"""
    path.write_text(content, encoding="utf-8")
    return path


def render_pbs_bundle(plan: pd.DataFrame, out_dir: str | Path, run_id: str, project_root: str = "$PWD") -> dict[str, Path]:
    out = ensure_dir(Path(out_dir) / "pbs")
    scripts = [render_pbs_stage_script(row.to_dict(), out, run_id=run_id, project_root=project_root) for _, row in plan.iterrows()]
    submit_all = out / "submit_all.sh"
    submit_all.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(f"qsub {p.name}" for p in scripts) + "\n", encoding="utf-8")
    submit_all.chmod(0o755)
    return {"scripts_dir": out, "submit_all": submit_all}
