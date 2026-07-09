from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.resource import DEFAULT_STAGE_RESOURCES, minutes_to_slurm_time
from hfauto_ops.core.run_index import artifact_record

HEAVY_ARTIFACT_TYPES = {
    "species",
    "species_preopt",
    "species_optimized",
    "reaction",
    "reaction_validated",
    "reaction_path_validated",
    "calculation",
}

STAGE_BY_ARTIFACT_TYPE = {
    "species": "preopt",
    "species_preopt": "dft-minima",
    "species_optimized": "sp",
    "reaction": "ts-search",
    "reaction_validated": "irc",
    "reaction_path_validated": "thermo",
    "calculation": "sp",
}


def _safe_stage(stage: str) -> str:
    return str(stage).replace("_", "-")


def _resource_for_stage(stage: str, ops_config: dict[str, Any] | None = None) -> dict[str, Any]:
    ops_config = ops_config or {}
    profile = ops_config.get("resource_profile") or {}
    base = dict(DEFAULT_STAGE_RESOURCES.get(_safe_stage(stage), {"ncores": 1, "memory_gb": 4, "time_min": 60, "queue": "short"}))
    overrides = (profile.get("stages") or {}).get(_safe_stage(stage), {})
    base.update(overrides)
    return base


def build_artifact_task_plan(run_dir: str | Path, ops_config: dict[str, Any] | None = None) -> pd.DataFrame:
    """Create conservative artifact-level tasks for HPC arrays.

    The task plan is not automatically submitted.  It is a reviewable mapping from
    chemistry artifacts to single-item stage commands.  Production deployments can
    convert this table to SLURM/PBS/LSF arrays after validating stage-level input
    manifests and backend licenses.
    """
    run = Path(run_dir)
    manifest = read_manifest(latest_manifest_path(run))
    rows: list[dict[str, Any]] = []
    for art in manifest.artifacts:
        if art.artifact_type not in HEAVY_ARTIFACT_TYPES:
            continue
        rec = artifact_record(art)
        stage = _safe_stage((art.method or {}).get("stage") or rec.get("stage") or STAGE_BY_ARTIFACT_TYPE.get(art.artifact_type, "review"))
        res = _resource_for_stage(stage, ops_config)
        selector = art.data.get("species_id") or art.data.get("reaction_id") or art.artifact_id
        task_id = f"task_{len(rows):05d}_{stage.replace('-', '_')}"
        rows.append(
            {
                "task_id": task_id,
                "run_id": manifest.run_id,
                "stage": stage,
                "artifact_id": art.artifact_id,
                "artifact_type": art.artifact_type,
                "species_id": rec.get("species_id"),
                "reaction_id": rec.get("reaction_id"),
                "mol_id": rec.get("mol_id"),
                "site_id": rec.get("site_id"),
                "hf_n": rec.get("hf_n"),
                "state": rec.get("state"),
                "status": rec.get("status"),
                "fallback_dummy": rec.get("fallback_dummy"),
                "quality_tier": rec.get("quality_tier"),
                "ncores": int(res.get("ncores", 1)),
                "memory_gb": int(res.get("memory_gb", 4)),
                "time_min": int(res.get("time_min", 60)),
                "time_slurm": minutes_to_slurm_time(int(res.get("time_min", 60))),
                "queue": res.get("queue", "short"),
                "selector": selector,
                "input_manifest": str(latest_manifest_path(run)),
                "command": f"python -m hfauto.cli.main run-stage {stage} --in {latest_manifest_path(run)} --out {run}/array_{stage}/${{SLURM_ARRAY_TASK_ID:-0}} --select {selector}",
                "review_note": "single-artifact selection requires stage support; otherwise use as review/task list",
            }
        )
    return pd.DataFrame(rows)


def write_artifact_task_plan(df: pd.DataFrame, out_dir: str | Path) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv = out / "artifact_task_plan.csv"
    json = out / "artifact_task_plan.json"
    df.to_csv(csv, index=False)
    json.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    return {"artifact_task_csv": csv, "artifact_task_json": json}


def write_array_submit_scripts(df: pd.DataFrame, out_dir: str | Path, scheduler: str = "slurm", chunk_size: int = 1) -> dict[str, Path]:
    out = Path(out_dir) / "artifact_arrays"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    if df.empty:
        marker = out / "README.txt"
        marker.write_text("No artifact-level tasks were generated.\n", encoding="utf-8")
        return {"artifact_array_readme": marker}
    scheduler = scheduler.lower()
    for stage, group in df.groupby("stage", dropna=False):
        stage_s = str(stage).replace("/", "_")
        table = out / f"tasks_{stage_s}.csv"
        group.to_csv(table, index=False)
        n = len(group)
        if scheduler == "slurm":
            script = out / f"array_{stage_s}.slurm"
            script.write_text(
                f"""#!/usr/bin/env bash
#SBATCH --job-name=hfauto_array_{stage_s}
#SBATCH --array=0-{max(0, n-1)}%{max(1, int(chunk_size))}
#SBATCH --cpus-per-task={int(group['ncores'].max())}
#SBATCH --mem={int(group['memory_gb'].max())}G
#SBATCH --time={group['time_slurm'].iloc[0]}
#SBATCH --partition={group['queue'].iloc[0]}
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
mkdir -p logs
TASK_TABLE={table}
python -m hfauto_ops.cli.main execute-task --task-table "$TASK_TABLE" --task-index "${{SLURM_ARRAY_TASK_ID}}" --dry-run
""",
                encoding="utf-8",
            )
            script.chmod(0o755)
            paths[f"array_{stage_s}"] = script
        elif scheduler == "pbs":
            script = out / f"array_{stage_s}.pbs"
            script.write_text(
                f"""#!/usr/bin/env bash
#PBS -N hfauto_array_{stage_s}
#PBS -J 0-{max(0, n-1)}
#PBS -l select=1:ncpus={int(group['ncores'].max())}:mem={int(group['memory_gb'].max())}gb
#PBS -l walltime={group['time_slurm'].iloc[0].replace('-', ':')}
#PBS -q {group['queue'].iloc[0]}
set -euo pipefail
TASK_TABLE={table}
python -m hfauto_ops.cli.main execute-task --task-table "$TASK_TABLE" --task-index "${{PBS_ARRAY_INDEX}}" --dry-run
""",
                encoding="utf-8",
            )
            script.chmod(0o755)
            paths[f"array_{stage_s}"] = script
        else:
            script = out / f"array_{stage_s}.sh"
            script.write_text(
                f"""#!/usr/bin/env bash
set -euo pipefail
TASK_TABLE={table}
idx=${{1:-0}}
python -m hfauto_ops.cli.main execute-task --task-table "$TASK_TABLE" --task-index "$idx" --dry-run
""",
                encoding="utf-8",
            )
            script.chmod(0o755)
            paths[f"array_{stage_s}"] = script
    submit_all = out / f"submit_all_{scheduler}.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p logs"]
    for key, path in sorted(paths.items()):
        if scheduler == "slurm" and path.suffix == ".slurm":
            lines.append(f"sbatch {path.name}")
        elif scheduler == "pbs" and path.suffix == ".pbs":
            lines.append(f"qsub {path.name}")
        else:
            lines.append(f"echo review {path.name}")
    submit_all.write_text("\n".join(lines) + "\n", encoding="utf-8")
    submit_all.chmod(0o755)
    paths["artifact_array_submit_all"] = submit_all
    paths["artifact_array_dir"] = out
    return paths


def execute_task_from_table(task_table: str | Path, task_index: int, dry_run: bool = True) -> dict[str, Any]:
    df = pd.read_csv(task_table)
    if df.empty:
        raise IndexError("Task table is empty")
    idx = int(task_index)
    if idx < 0 or idx >= len(df):
        raise IndexError(f"Task index {idx} outside 0..{len(df)-1}")
    row = df.iloc[idx].to_dict()
    cmd = str(row.get("command"))
    return {"task_index": idx, "task_id": row.get("task_id"), "dry_run": dry_run, "command": cmd, "row": row}
