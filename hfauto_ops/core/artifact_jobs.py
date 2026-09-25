from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, read_manifest, write_json
from hfauto.hpc.resources import estimate_artifact_resources
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.run_index import artifact_record

ARTIFACT_STAGE_HINTS: dict[str, str] = {
    "species": "preopt",
    "species_preopt": "dft-minima",
    "species_optimized": "sp",
    "calculation": "dft-minima",
    "reaction": "ts-search",
    "reaction_validated": "irc",
    "reaction_path_validated": "irc",
    "ts_path": "ts-search",
    "irc": "irc",
    "irc_attempt": "irc",
    "thermo": "thermo",
    "kinetics": "kinetics",
    "molecule_enriched": "detect-sites",
    "site": "conformers",
    "conformer": "build-complexes",
}

HEAVY_ARTIFACT_TYPES = {
    "species",
    "species_preopt",
    "species_optimized",
    "calculation",
    "reaction",
    "reaction_validated",
    "reaction_path_validated",
    "ts_path",
    "irc",
    "irc_attempt",
    "thermo",
    "kinetics",
}


def _is_truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.lower() in {"true", "1", "yes", "y"}
    return bool(value)


def infer_artifact_target_stage(artifact: Any) -> str:
    method = artifact.method or {}
    data = artifact.data or {}
    stage = method.get("stage") or data.get("stage")
    if (
        stage in {"dft-minima", "ts-search", "irc", "sp", "thermo", "kinetics", "preopt"}
        and (
            artifact.status.status == "failed"
            or artifact.qc.get("fallback_dummy")
            or data.get("fallback_dummy")
        )
    ):
        return str(stage)
    # Calculation artifacts can represent SP or opt/freq depending on task.
    if artifact.artifact_type == "calculation":
        task = str(method.get("task") or data.get("task") or "").lower()
        if "sp" in task or "single" in task:
            return "sp"
        if "ts" in task or "neb" in task or "optts" in task:
            return "ts-search"
        if "irc" in task:
            return "irc"
        return "dft-minima"
    return ARTIFACT_STAGE_HINTS.get(artifact.artifact_type, "unknown")


def _should_plan_artifact(artifact: Any, include_successful: bool = False, include_lightweight: bool = False) -> bool:
    data = artifact.data or {}
    qc = artifact.qc or {}
    if artifact.status.status == "failed":
        return True
    if _is_truthy(qc.get("fallback_dummy")) or _is_truthy(data.get("fallback_dummy")):
        return True
    if not _is_truthy(data.get("production_thermo_ready", True)) and artifact.artifact_type in {"thermo", "kinetics"}:
        return True
    if not _is_truthy(data.get("real_orca_executed", True)) and artifact.artifact_type == "calculation":
        return True
    if include_successful and artifact.artifact_type in HEAVY_ARTIFACT_TYPES:
        return True
    return bool(include_lightweight and artifact.artifact_type in ARTIFACT_STAGE_HINTS)


def artifact_job_fingerprint(artifact: Any, target_stage: str, method_config: dict[str, Any] | None = None) -> str:
    data = artifact.data or {}
    method = artifact.method or {}
    payload = {
        "target_stage": target_stage,
        "artifact_type": artifact.artifact_type,
        "species_id": data.get("species_id") or data.get("source_species_id"),
        "reaction_id": data.get("reaction_id"),
        "mol_id": data.get("mol_id"),
        "site_id": data.get("site_id"),
        "hf_n": data.get("hf_n"),
        "state": data.get("state"),
        "engine": method.get("engine"),
        "method_id": method.get("method_id") or method.get("model") or method.get("method"),
        "task": method.get("task"),
        "method_config": method_config or {},
    }
    return fingerprint_dict(payload)


def build_artifact_job_plan(
    run_dir: str | Path,
    out_dir: str | Path | None = None,
    pipeline_config: str | Path | None = None,
    ops_config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    ops_config = ops_config or {}
    manifest_path = latest_manifest_path(run_dir)
    manifest = read_manifest(manifest_path)
    include_successful = bool(ops_config.get("include_successful_artifact_jobs", False))
    include_lightweight = bool(ops_config.get("include_lightweight_artifact_jobs", False))
    root = Path(out_dir or Path(run_dir) / "artifact_jobs")
    rows: list[dict[str, Any]] = []
    for idx, artifact in enumerate(manifest.artifacts):
        if not _should_plan_artifact(artifact, include_successful=include_successful, include_lightweight=include_lightweight):
            continue
        target_stage = infer_artifact_target_stage(artifact)
        if target_stage == "unknown":
            continue
        req = estimate_artifact_resources(artifact, stage=target_stage, config=ops_config)
        fp = artifact_job_fingerprint(artifact, target_stage, ops_config.get("method_overrides"))
        job_id = f"ajob_{idx:05d}_{target_stage.replace('-', '_')}_{fp[:8]}"
        job_out = root / "jobs" / job_id
        command = (
            f"HFAUTO_TARGET_ARTIFACT_ID={artifact.artifact_id} "
            f"hfauto run-stage {target_stage} --in {manifest_path} --out {job_out}"
        )
        rec = artifact_record(artifact)
        rows.append(
            {
                "artifact_job_id": job_id,
                "target_stage": target_stage,
                "target_artifact_id": artifact.artifact_id,
                "artifact_type": artifact.artifact_type,
                "run_id": manifest.run_id,
                "latest_manifest": str(manifest_path),
                "fingerprint": fp,
                "cache_key": fp,
                "scheduler_group": f"{target_stage}:{(artifact.method or {}).get('engine') or 'unknown'}:{(artifact.method or {}).get('method_id') or 'unknown'}",
                "engine": rec.get("engine") or req.backend,
                "method_id": rec.get("method_id"),
                "task": rec.get("task"),
                "mol_id": rec.get("mol_id"),
                "site_id": rec.get("site_id"),
                "species_id": rec.get("species_id"),
                "reaction_id": rec.get("reaction_id"),
                "state": rec.get("state"),
                "hf_n": rec.get("hf_n"),
                "status": rec.get("status"),
                "failure_category": rec.get("failure_category"),
                "failure_reason": rec.get("failure_reason"),
                "fallback_dummy": rec.get("fallback_dummy"),
                "real_orca_executed": rec.get("real_orca_executed"),
                "quality_tier": rec.get("quality_tier"),
                "confidence_score": rec.get("confidence_score"),
                "ncores": req.ncores,
                "memory_gb": req.memory_gb,
                "walltime": req.walltime,
                "partition": req.partition,
                "queue": req.queue,
                "gpus": req.gpus,
                "priority": req.priority,
                "reason": req.reason,
                "output_dir": str(job_out),
                "command": command,
                "production_reason": _production_reason(artifact),
                "pipeline_config": str(pipeline_config) if pipeline_config else None,
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        # Scheduler/reporting columns are normalized here once.
        df["stage"] = df["target_stage"]
        df["job_id"] = df["artifact_job_id"]
        df["target_artifact_type"] = df["artifact_type"]
        df["array_index"] = range(len(df))
        df["time_slurm"] = df.get("walltime", "01:00:00")
        df = df.sort_values(["priority", "target_stage", "artifact_job_id"], ascending=[False, True, True]).reset_index(drop=True)
        df["array_index"] = range(len(df))
    return df


def _production_reason(artifact: Any) -> str:
    reasons: list[str] = []
    if artifact.status.status == "failed":
        reasons.append("failed")
    if artifact.qc.get("fallback_dummy") or artifact.data.get("fallback_dummy"):
        reasons.append("fallback_dummy")
    if artifact.artifact_type == "calculation" and not artifact.data.get("real_orca_executed", True):
        reasons.append("needs_real_qm")
    if artifact.artifact_type in {"thermo", "kinetics"} and not artifact.data.get("production_thermo_ready", True):
        reasons.append("needs_production_connector")
    return ";".join(reasons) or "production_campaign"


def summarize_artifact_arrays(job_plan: pd.DataFrame, array_chunk: int = 1) -> pd.DataFrame:
    if job_plan.empty:
        return pd.DataFrame(columns=["array_id", "scheduler_group", "target_stage", "n_jobs", "array_chunk", "n_array_tasks", "ncores", "memory_gb", "walltime", "queue", "job_ids"])
    rows: list[dict[str, Any]] = []
    for idx, (group, gdf) in enumerate(job_plan.groupby("scheduler_group", dropna=False)):
        n = len(gdf)
        rows.append(
            {
                "array_id": f"array_{idx:03d}_{str(group).replace(':','_').replace('-','_')[:60]}",
                "scheduler_group": group,
                "target_stage": ";".join(sorted(set(gdf["target_stage"].astype(str)))),
                "n_jobs": n,
                "array_chunk": max(1, int(array_chunk)),
                "n_array_tasks": int((n + max(1, int(array_chunk)) - 1) / max(1, int(array_chunk))),
                "ncores": int(gdf["ncores"].max()),
                "memory_gb": float(gdf["memory_gb"].max()),
                "walltime": str(gdf.sort_values("priority", ascending=False)["walltime"].iloc[0]),
                "queue": str(gdf["queue"].dropna().iloc[0]) if gdf["queue"].notna().any() else "",
                "job_ids": ";".join(gdf["artifact_job_id"].astype(str).tolist()),
            }
        )
    return pd.DataFrame(rows)


def write_artifact_job_plan(
    run_dir: str | Path,
    out_dir: str | Path,
    pipeline_config: str | Path | None = None,
    ops_config: dict[str, Any] | None = None,
) -> dict[str, Path]:
    ops_config = ops_config or {}
    out = ensure_dir(out_dir)
    jobs = build_artifact_job_plan(run_dir, out_dir=out, pipeline_config=pipeline_config, ops_config=ops_config)
    arrays = summarize_artifact_arrays(jobs, array_chunk=int(ops_config.get("artifact_array_chunk", ops_config.get("array_chunk", 1)) or 1))
    jobs_csv = out / "artifact_job_plan.csv"
    jobs_json = out / "artifact_job_plan.json"
    arrays_csv = out / "artifact_job_arrays.csv"
    arrays_json = out / "artifact_job_arrays.json"
    commands = out / "artifact_job_commands.sh"
    jobs.to_csv(jobs_csv, index=False)
    arrays.to_csv(arrays_csv, index=False)
    jobs_json.write_text(jobs.to_json(orient="records", indent=2), encoding="utf-8")
    arrays_json.write_text(arrays.to_json(orient="records", indent=2), encoding="utf-8")
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Artifact-level rerun commands. Review before production use.",
        *jobs.get("command", pd.Series(dtype=str)).dropna().astype(str).tolist(),
    ]
    commands.write_text("\n".join(lines) + "\n", encoding="utf-8")
    commands.chmod(0o755)
    manifest = write_json(
        out / "artifact_job_plan_manifest.json",
        {
            "schema_version": "hfauto.ops.artifact_job_plan.v1",
            "n_artifact_jobs": len(jobs),
            "n_array_groups": len(arrays),
            "pipeline_config": str(pipeline_config) if pipeline_config else None,
            "outputs": {
                "jobs_csv": str(jobs_csv),
                "jobs_json": str(jobs_json),
                "arrays_csv": str(arrays_csv),
                "arrays_json": str(arrays_json),
                "commands": str(commands),
            },
            "note": "Artifact-level execution is planned conservatively. Commands set HFAUTO_TARGET_ARTIFACT_ID for future fine-grained stage filters.",
        },
    )
    return {
        "artifact_jobs_csv": jobs_csv,
        "artifact_job_csv": jobs_csv,
        "artifact_jobs_json": jobs_json,
        "artifact_job_json": jobs_json,
        "artifact_arrays_csv": arrays_csv,
        "artifact_arrays_json": arrays_json,
        "artifact_commands": commands,
        "artifact_plan_manifest": manifest,
        "artifact_job_plan_manifest": manifest,
    }
