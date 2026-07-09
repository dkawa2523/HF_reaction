from __future__ import annotations

"""Create reviewable job plans from a run manifest or pipeline config."""

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, write_json
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.hpc.resources import estimate_artifact_resources, estimate_stage_resources


@dataclass
class JobRecord:
    job_id: str
    run_id: str
    stage: str
    target_artifact_id: str | None
    artifact_type: str | None
    backend: str
    command: str
    input_manifest: str | None
    output_dir: str
    ncores: int
    memory_gb: float
    walltime: str
    partition: str | None
    queue: str | None
    gpus: int
    priority: int
    status: str = "planned"
    retry_of: str | None = None
    reason: str = ""
    cache_key: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _stage_command(stage: str, input_manifest: str | None, out_dir: str, stage_config: str | None = None) -> str:
    if stage == "ingest":
        return f"hfauto ingest --sdf $HFAUTO_INPUT_SDF --out {out_dir} --run-id $HFAUTO_RUN_ID"
    cmd = f"hfauto run-stage {stage} --in {input_manifest or '$HFAUTO_INPUT_MANIFEST'} --out {out_dir}"
    if stage_config:
        cmd += f" --config {stage_config}"
    return cmd


def planned_stage_jobs(
    run_id: str,
    stages: list[dict[str, Any]],
    run_dir: str | Path,
    input_manifest: str | None = None,
    global_config: dict[str, Any] | None = None,
) -> list[JobRecord]:
    jobs: list[JobRecord] = []
    current_manifest = input_manifest
    run_dir = Path(run_dir)
    enabled = [s for s in stages if s.get("enabled", True)]
    for idx, stage_cfg in enumerate(enabled):
        stage = stage_cfg["name"]
        out_dir = str(run_dir / f"{idx:02d}_{stage.replace('_','-')}")
        req = estimate_stage_resources(stage, stage_cfg | {"resources": (global_config or {}).get("resources", {})})
        command = _stage_command(stage, current_manifest, out_dir, stage_cfg.get("stage_config"))
        payload = {"run_id": run_id, "stage": stage, "idx": idx, "command": command, "resources": req.to_dict()}
        job_id = f"job_{idx:02d}_{stage.replace('-','_')}_{fingerprint_dict(payload)[:8]}"
        jobs.append(JobRecord(
            job_id=job_id,
            run_id=run_id,
            stage=stage,
            target_artifact_id=None,
            artifact_type=None,
            backend=req.backend,
            command=command,
            input_manifest=current_manifest,
            output_dir=out_dir,
            ncores=req.ncores,
            memory_gb=req.memory_gb,
            walltime=req.walltime,
            partition=req.partition,
            queue=req.queue,
            gpus=req.gpus,
            priority=req.priority,
            reason=req.reason,
            cache_key=fingerprint_dict(payload),
        ))
        current_manifest = str(Path(out_dir) / "manifest.json")
    return jobs


def retry_jobs_from_manifest(
    manifest: Manifest,
    run_dir: str | Path,
    retry_root: str | Path,
    config: dict[str, Any] | None = None,
) -> list[JobRecord]:
    jobs: list[JobRecord] = []
    retry_root = Path(retry_root)
    for idx, artifact in enumerate(manifest.artifacts):
        if artifact.status.status != "failed":
            continue
        stage = artifact.method.get("stage") if artifact.method else None
        stage = stage or artifact.data.get("stage") or artifact.artifact_type
        req = estimate_artifact_resources(artifact, stage=stage, config=config)
        out_dir = str(retry_root / artifact.artifact_id)
        command = f"hfauto run-stage {stage} --in $HFAUTO_INPUT_MANIFEST --out {out_dir}"
        payload = {"retry_of": artifact.artifact_id, "stage": stage, "reason": artifact.status.reason, "resources": req.to_dict()}
        jobs.append(JobRecord(
            job_id=f"retry_{idx:04d}_{fingerprint_dict(payload)[:8]}",
            run_id=manifest.run_id,
            stage=str(stage),
            target_artifact_id=artifact.artifact_id,
            artifact_type=artifact.artifact_type,
            backend=req.backend,
            command=command,
            input_manifest=None,
            output_dir=out_dir,
            ncores=req.ncores,
            memory_gb=req.memory_gb,
            walltime=req.walltime,
            partition=req.partition,
            queue=req.queue,
            gpus=req.gpus,
            priority=req.priority,
            status="retry_planned",
            retry_of=artifact.artifact_id,
            reason=artifact.status.reason or "failed artifact",
            cache_key=fingerprint_dict(payload),
        ))
    return jobs


def write_job_plan(jobs: list[JobRecord], out_dir: str | Path, prefix: str = "job_plan") -> dict[str, Path]:
    out = ensure_dir(out_dir)
    rows = [j.to_dict() for j in jobs]
    csv_path = out / f"{prefix}.csv"
    json_path = out / f"{prefix}.json"
    pd.DataFrame(rows, columns=list(JobRecord.__dataclass_fields__.keys())).to_csv(csv_path, index=False)
    write_json(json_path, {"schema_version": "hfauto.job_plan.v1", "n_jobs": len(rows), "jobs": rows})
    return {"csv": csv_path, "json": json_path}


def job_artifact(job_paths: dict[str, Path], n_jobs: int, kind: str = "job_plan") -> Artifact:
    return Artifact(
        artifact_id=kind,
        artifact_type="hpc_plan",
        paths={k: str(v) for k, v in job_paths.items()},
        data={"n_jobs": int(n_jobs), "plan_type": kind},
    )
