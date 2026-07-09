from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, write_json
from hfauto_ops.schedulers.base import SchedulerAdapter, SchedulerResult
from hfauto_ops.schedulers.slurm import SlurmAdapter
from hfauto_ops.schedulers.pbs import PBSAdapter
from hfauto_ops.schedulers.lsf import LSFAdapter


def get_scheduler(name: str) -> SchedulerAdapter:
    n = (name or "slurm").lower()
    if n == "slurm":
        return SlurmAdapter()
    if n == "pbs":
        return PBSAdapter()
    if n == "lsf":
        return LSFAdapter()
    raise KeyError(f"Unknown scheduler: {name}")


def _script_candidates(ops_dir: Path, scheduler: str) -> list[Path]:
    candidates = []
    for sub in [scheduler, "slurm", "pbs", "lsf"]:
        d = ops_dir / sub
        if d.exists():
            candidates.extend(sorted(d.glob("submit_*.slurm")))
            candidates.extend(sorted(d.glob("submit_*.pbs.sh")))
            candidates.extend(sorted(d.glob("submit_*.lsf.sh")))
            candidates.extend(sorted(d.glob("artifact_array*.slurm")))
            candidates.extend(sorted(d.glob("artifact_array*.pbs.sh")))
            candidates.extend(sorted(d.glob("artifact_array*.lsf.sh")))
    return candidates


def submit_scripts(
    ops_dir: str | Path,
    scheduler: str = "slurm",
    dry_run: bool = True,
    allow_execute: bool = False,
    limit: int | None = None,
) -> dict[str, Path]:
    out = ensure_dir(Path(ops_dir))
    adapter = get_scheduler(scheduler)
    records = []
    for script in _script_candidates(out, scheduler)[: limit or None]:
        cmd = adapter.submit_command(script)
        result = adapter.run_scheduler_command("submit", cmd, out, dry_run=dry_run, allow_execute=allow_execute)
        rec = result.to_dict() | {"script": str(script)}
        records.append(rec)
    df = pd.DataFrame(records)
    csv = out / "scheduler_submit_results.csv"
    json = out / "scheduler_submit_results.json"
    df.to_csv(csv, index=False)
    write_json(json, {"scheduler": scheduler, "dry_run": dry_run, "n_results": len(records), "results": records})
    return {"csv": csv, "json": json}


def scheduler_status(
    ops_dir: str | Path,
    scheduler: str = "slurm",
    job_id: str | None = None,
    user: str | None = None,
    dry_run: bool = True,
    allow_execute: bool = False,
) -> dict[str, Path]:
    out = ensure_dir(Path(ops_dir))
    adapter = get_scheduler(scheduler)
    result = adapter.run_scheduler_command("status", adapter.status_command(job_id=job_id, user=user), out, dry_run=dry_run, allow_execute=allow_execute)
    path = out / "scheduler_status.json"
    write_json(path, result.to_dict())
    csv = out / "scheduler_status.csv"
    pd.DataFrame([result.to_dict()]).to_csv(csv, index=False)
    return {"json": path, "csv": csv}


def scheduler_cancel(
    ops_dir: str | Path,
    scheduler: str = "slurm",
    job_id: str | None = None,
    job_file: str | Path | None = None,
    dry_run: bool = True,
    allow_execute: bool = False,
) -> dict[str, Path]:
    out = ensure_dir(Path(ops_dir))
    adapter = get_scheduler(scheduler)
    ids: list[str] = []
    if job_id:
        ids.append(str(job_id))
    if job_file and Path(job_file).exists():
        for line in Path(job_file).read_text(encoding="utf-8").splitlines():
            val = line.strip()
            if val:
                ids.append(val)
    records = []
    for jid in ids:
        result = adapter.run_scheduler_command("cancel", adapter.cancel_command(jid), out, dry_run=dry_run, allow_execute=allow_execute)
        records.append(result.to_dict() | {"job_id": jid})
    csv = out / "scheduler_cancel_results.csv"
    json = out / "scheduler_cancel_results.json"
    pd.DataFrame(records).to_csv(csv, index=False)
    write_json(json, {"scheduler": scheduler, "dry_run": dry_run, "n_results": len(records), "results": records})
    return {"csv": csv, "json": json}
