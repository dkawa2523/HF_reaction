from __future__ import annotations

from hfauto_ops.schedulers.base import SchedulerAdapter
from hfauto_ops.schedulers.lsf import LSFAdapter
from hfauto_ops.schedulers.pbs import PBSAdapter
from hfauto_ops.schedulers.slurm import SlurmAdapter


def get_scheduler(name: str | None) -> SchedulerAdapter:
    key = (name or "slurm").lower().replace("_", "-")
    if key == "slurm":
        return SlurmAdapter()
    if key in {"pbs", "torque"}:
        return PBSAdapter()
    if key == "lsf":
        return LSFAdapter()
    if key in {"local", "bash"}:
        return SchedulerAdapter()
    raise KeyError(f"Unknown scheduler: {name}")


from pathlib import Path

import pandas as pd


def _ids(job_table: pd.DataFrame) -> list[str]:
    if job_table is None or job_table.empty:
        return []
    for col in ["scheduler_job_id", "job_id", "slurm_job_id", "pbs_job_id", "lsf_job_id"]:
        if col in job_table.columns:
            vals = [str(v) for v in job_table[col].dropna().tolist() if str(v).strip() and str(v).lower() != "nan"]
            if vals:
                return vals
    return []

def _write(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(0o755)
    return path

class _CommandBundle:
    def __init__(self, scheduler: str, action: str, command_file: str):
        self.scheduler=scheduler; self.action=action; self.command_file=command_file; self.dry_run=True; self.executed=False
    def to_dict(self):
        return {"scheduler": self.scheduler, "action": self.action, "command_file": self.command_file, "dry_run": True, "executed": False}

def render_scheduler_status(scheduler: str, job_table: pd.DataFrame, out_dir: str | Path):
    sched=(scheduler or "slurm").lower(); ids=_ids(job_table); lines=["#!/usr/bin/env bash", "set -euo pipefail"]
    if sched == "pbs":
        lines += [f"qstat -f {x} || true" for x in ids] or ["qstat -u $USER || true"]
    elif sched == "lsf":
        lines += [f"bjobs -l {x} || true" for x in ids] or ["bjobs || true"]
    elif sched == "local":
        lines += ["echo local scheduler has no status API"]
    else:
        lines += [f"sacct -j {x} --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode || squeue -j {x} || true" for x in ids] or ["squeue -u $USER || true"]
    p=_write(Path(out_dir)/f"{sched}_status.sh", lines)
    return _CommandBundle(sched, "status", str(p))

def render_scheduler_cancel(scheduler: str, job_table: pd.DataFrame, out_dir: str | Path):
    sched=(scheduler or "slurm").lower(); ids=_ids(job_table); lines=["#!/usr/bin/env bash", "set -euo pipefail"]
    if sched == "pbs":
        lines += [f"qdel {x} || true" for x in ids] or ["echo 'No job ids available for qdel' >&2"]
    elif sched == "lsf":
        lines += [f"bkill {x} || true" for x in ids] or ["echo 'No job ids available for bkill' >&2"]
    elif sched == "local":
        lines += ["echo local scheduler has no cancel API"]
    else:
        lines += [f"scancel {x} || true" for x in ids] or ["echo 'No job ids available for scancel' >&2"]
    p=_write(Path(out_dir)/f"{sched}_cancel.sh", lines)
    return _CommandBundle(sched, "cancel", str(p))
