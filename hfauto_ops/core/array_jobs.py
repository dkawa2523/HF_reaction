from __future__ import annotations

"""Artifact-level job-array planning and execution helpers.

The helpers in this module are operations-layer utilities. They do not alter a
completed scientific run. They create target-specific manifests and run a single
workflow stage against that manifest when explicitly invoked by an HPC array
worker.
"""

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, read_jsonl, read_manifest, write_json, write_jsonl, write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.stages.base import StageContext
from hfauto_ops.core.resource import DEFAULT_STAGE_RESOURCES, minutes_to_slurm_time

# The target artifact types are selected so that the existing stage code loops
# over one target when we create a sliced manifest. Dependencies needed for the
# stage are retained by ``slice_manifest_for_target`` below.
STAGE_TARGET_TYPES: dict[str, tuple[str, ...]] = {
    "preopt": ("species",),
    "dft-minima": ("species_preopt", "species"),
    "ts-search": ("reaction",),
    "irc": ("reaction_validated",),
    "sp": ("species_optimized", "species_preopt", "species"),
}

STAGE_ORDER = ["preopt", "dft-minima", "ts-search", "irc", "sp"]


def canonical_stage(name: str) -> str:
    return str(name).replace("_", "-")


def canonical_species_id(artifact: Artifact) -> str:
    data = artifact.data or {}
    return str(data.get("species_id") or data.get("source_species_id") or artifact.artifact_id)


def _enabled_stage_configs(pipeline_config: str | Path | None) -> list[dict[str, Any]]:
    if not pipeline_config:
        return []
    cfg = load_yaml(pipeline_config)
    return [s for s in cfg.get("stages", []) if s.get("enabled", True)]


def _stage_config_map(pipeline_config: str | Path | None) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for stage in _enabled_stage_configs(pipeline_config):
        out[canonical_stage(stage.get("name"))] = {k: v for k, v in stage.items() if k not in {"enabled"}}
    return out


def _target_artifacts(manifest: Manifest, stage: str) -> list[Artifact]:
    types = STAGE_TARGET_TYPES.get(stage, ())
    rows: list[Artifact] = []
    seen_keys: set[str] = set()
    for typ in types:
        for artifact in manifest.latest_artifacts(typ):
            if artifact.status.status != "success":
                continue
            if stage in {"dft-minima", "sp"}:
                key = canonical_species_id(artifact)
            elif stage in {"ts-search", "irc"}:
                key = str(artifact.data.get("reaction_id") or artifact.artifact_id)
            else:
                key = artifact.artifact_id
            if key in seen_keys:
                continue
            seen_keys.add(key)
            rows.append(artifact)
    return rows


def _resource_for_stage(stage: str, ops_config: dict[str, Any] | None = None) -> dict[str, Any]:
    ops_config = ops_config or {}
    profile = ops_config.get("resource_profile") or {}
    defaults = DEFAULT_STAGE_RESOURCES.get(stage, {"ncores": 1, "memory_gb": 4, "time_min": 60, "queue": "short"})
    res = dict(defaults)
    res.update((profile.get("stages", {}) or {}).get(stage, {}) or {})
    res.setdefault("time_slurm", minutes_to_slurm_time(int(res.get("time_min", 60))))
    return res


def build_array_items(
    run_dir: str | Path,
    pipeline_config: str | Path | None,
    ops_config: dict[str, Any] | None = None,
    stages: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Create one target item per eligible artifact for heavy stages.

    The resulting table is a plan. Each item can be executed by ``run_array_task``.
    The stage config is stored per item so that different backends/methods can be
    preserved from the original pipeline YAML.
    """
    run = Path(run_dir)
    manifest_path = latest_manifest_path(run)
    manifest = read_manifest(manifest_path)
    stage_cfgs = _stage_config_map(pipeline_config)
    selected_stages = [canonical_stage(s) for s in (stages or STAGE_ORDER)]
    rows: list[dict[str, Any]] = []
    idx = 0
    for stage in selected_stages:
        if stage not in STAGE_TARGET_TYPES:
            continue
        # Skip stages not present in the pipeline config unless no config was supplied.
        if pipeline_config and stage not in stage_cfgs:
            continue
        stage_cfg = stage_cfgs.get(stage, {"name": stage})
        res = _resource_for_stage(stage, ops_config)
        for target in _target_artifacts(manifest, stage):
            target_key = (
                canonical_species_id(target)
                if target.artifact_type in {"species", "species_preopt", "species_optimized"}
                else str(target.data.get("reaction_id") or target.artifact_id)
            )
            payload = {"run_id": manifest.run_id, "stage": stage, "target_artifact_id": target.artifact_id, "target_key": target_key, "stage_config": stage_cfg}
            item_id = f"arr_{idx:05d}_{stage.replace('-', '_')}_{fingerprint_dict(payload)[:8]}"
            out_dir = run / "array_runs" / stage / item_id
            rows.append(
                {
                    "array_index": idx,
                    "array_item_id": item_id,
                    "run_id": manifest.run_id,
                    "stage": stage,
                    "target_artifact_id": target.artifact_id,
                    "target_artifact_type": target.artifact_type,
                    "target_key": target_key,
                    "mol_id": target.data.get("mol_id"),
                    "site_id": target.data.get("site_id"),
                    "species_id": target.data.get("species_id") or target.data.get("source_species_id"),
                    "reaction_id": target.data.get("reaction_id"),
                    "state": target.data.get("state"),
                    "hf_n": target.data.get("hf_n"),
                    "input_manifest": str(manifest_path),
                    "output_dir": str(out_dir),
                    "target_manifest": str(out_dir / "target_manifest.json"),
                    "stage_config_json": json.dumps(stage_cfg, ensure_ascii=False, sort_keys=True, default=str),
                    "ncores": int(res.get("ncores", 1)),
                    "memory_gb": int(res.get("memory_gb", 4)),
                    "time_min": int(res.get("time_min", 60)),
                    "time_slurm": res.get("time_slurm") or minutes_to_slurm_time(int(res.get("time_min", 60))),
                    "queue": res.get("queue", "short"),
                    "target_filter_mode": "manifest_slice",
                    "cache_key": fingerprint_dict(payload),
                }
            )
            idx += 1
    return pd.DataFrame(rows)


def _is_dependency_for_target(artifact: Artifact, target: Artifact, stage: str) -> bool:
    # Always keep non-heavy metadata and all molecules/sites/conformers because they are small.
    if artifact.artifact_type in {"molecule", "molecule_enriched", "site", "conformer", "ranking", "table", "report"}:
        return True
    if artifact.artifact_id == target.artifact_id:
        return True

    target_sid = canonical_species_id(target)
    target_rxn = str(target.data.get("reaction_id") or "")
    data = artifact.data or {}

    if stage in {"preopt"}:
        return False

    if stage in {"dft-minima", "sp"}:
        if artifact.artifact_type in {"species", "species_preopt", "species_optimized"}:
            return canonical_species_id(artifact) == target_sid
        # Keep calculations already attached to the target species for lineage.
        if artifact.artifact_type == "calculation":
            return str(data.get("species_id")) == target_sid
        return False

    if stage == "ts-search":
        if artifact.artifact_type == "reaction":
            return artifact.artifact_id == target.artifact_id
        # Keep all species because TS endpoints point to species ids; this is still
        # small compared with the DFT/TS work and keeps slicing robust.
        return artifact.artifact_type in {"species", "species_preopt", "species_optimized", "calculation"}

    if stage == "irc":
        if artifact.artifact_type == "reaction_validated":
            return artifact.artifact_id == target.artifact_id
        if artifact.artifact_type in {"reaction", "species", "species_preopt", "species_optimized", "calculation", "ts_path"}:
            # Keep TS and endpoint species plus calculations. Reaction-level filtering is
            # conservative because different backends encode endpoint IDs differently.
            if artifact.artifact_type == "reaction" and target_rxn:
                return str(data.get("reaction_id")) == target_rxn
            return True
    return False


def slice_manifest_for_target(manifest: Manifest, target_artifact_id: str, stage: str) -> Manifest:
    target = manifest.find(target_artifact_id)
    if target is None:
        raise KeyError(f"Target artifact not found: {target_artifact_id}")
    out = Manifest.new(run_id=manifest.run_id, stage=f"array-input-{stage}", parents=[manifest.manifest_id], metadata={"target_artifact_id": target_artifact_id, "target_stage": stage})
    kept: list[Artifact] = []
    seen: set[tuple[str, str]] = set()
    for artifact in manifest.artifacts:
        if _is_dependency_for_target(artifact, target, stage):
            key = (artifact.artifact_type, artifact.artifact_id)
            if key not in seen:
                kept.append(artifact)
                seen.add(key)
    out.artifacts = kept
    return out


def write_array_plan(
    run_dir: str | Path,
    out_dir: str | Path,
    pipeline_config: str | Path | None,
    ops_config: dict[str, Any] | None = None,
    stages: Iterable[str] | None = None,
) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    df = build_array_items(run_dir, pipeline_config=pipeline_config, ops_config=ops_config, stages=stages)
    csv = out / "array_job_plan.csv"
    jsonl = out / "array_items.jsonl"
    summary = out / "array_job_summary.json"
    df.to_csv(csv, index=False)
    write_jsonl(df.to_dict(orient="records"), jsonl)
    by_stage = df.groupby("stage").size().to_dict() if not df.empty else {}
    write_json(summary, {"schema_version": "hfauto.array_plan.v1", "n_array_items": int(len(df)), "by_stage": {str(k): int(v) for k, v in by_stage.items()}, "pipeline_config": str(pipeline_config) if pipeline_config else None})
    worker = out / "run_array_task.py"
    worker.write_text(
        "#!/usr/bin/env python\n"
        "from hfauto_ops.cli.main import app\n"
        "if __name__ == '__main__':\n"
        "    app()\n",
        encoding="utf-8",
    )
    worker.chmod(0o755)
    return {"array_csv": csv, "array_jsonl": jsonl, "array_summary": summary, "array_worker": worker}


def run_array_task(array_items: str | Path, task_id: int, dry_run: bool = False) -> dict[str, Any]:
    items = read_jsonl(array_items)
    if task_id < 0 or task_id >= len(items):
        raise IndexError(f"array task id {task_id} outside 0..{len(items)-1}")
    item = items[task_id]
    stage = str(item["stage"])
    input_manifest = read_manifest(item["input_manifest"])
    target_manifest = slice_manifest_for_target(input_manifest, item["target_artifact_id"], stage)
    out_dir = ensure_dir(item["output_dir"])
    target_manifest_path = write_manifest(target_manifest, out_dir)
    stage_cfg = json.loads(item.get("stage_config_json") or "{}")
    if dry_run:
        result = {"status": "dry_run", "array_item": item, "target_manifest": str(target_manifest_path)}
        write_json(out_dir / "array_task_result.json", result)
        return result
    from hfauto.stages.registry import get_stage

    stage_obj = get_stage(stage)
    context = StageContext(out_dir=out_dir / "stage_output", run_id=str(item.get("run_id") or input_manifest.run_id), global_config={})
    out_manifest = stage_obj.run(target_manifest, stage_cfg, context)
    output_manifest_path = write_manifest(out_manifest, context.out_dir)
    result = {"status": "success", "array_item": item, "target_manifest": str(target_manifest_path), "output_manifest": str(output_manifest_path)}
    write_json(out_dir / "array_task_result.json", result)
    return result
