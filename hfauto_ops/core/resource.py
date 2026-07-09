from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path


DEFAULT_STAGE_RESOURCES: dict[str, dict[str, Any]] = {
    "ingest": {"ncores": 1, "memory_gb": 1, "time_min": 10, "queue": "short"},
    "enrich": {"ncores": 1, "memory_gb": 2, "time_min": 30, "queue": "short"},
    "detect-sites": {"ncores": 1, "memory_gb": 2, "time_min": 10, "queue": "short"},
    "conformers": {"ncores": 4, "memory_gb": 8, "time_min": 120, "queue": "medium"},
    "build-hf": {"ncores": 1, "memory_gb": 4, "time_min": 30, "queue": "short"},
    "preopt": {"ncores": 4, "memory_gb": 8, "time_min": 240, "queue": "medium"},
    "dft-minima": {"ncores": 16, "memory_gb": 32, "time_min": 1440, "queue": "long"},
    "ts-search": {"ncores": 16, "memory_gb": 48, "time_min": 2880, "queue": "long"},
    "irc": {"ncores": 16, "memory_gb": 48, "time_min": 1440, "queue": "long"},
    "sp": {"ncores": 16, "memory_gb": 64, "time_min": 720, "queue": "long"},
    "thermo": {"ncores": 1, "memory_gb": 4, "time_min": 30, "queue": "short"},
    "descriptors": {"ncores": 1, "memory_gb": 4, "time_min": 30, "queue": "short"},
    "kinetics": {"ncores": 1, "memory_gb": 4, "time_min": 60, "queue": "short"},
    "calibrate": {"ncores": 1, "memory_gb": 4, "time_min": 60, "queue": "short"},
    "rank": {"ncores": 1, "memory_gb": 4, "time_min": 30, "queue": "short"},
    "viz": {"ncores": 1, "memory_gb": 8, "time_min": 120, "queue": "short"},
    "ops": {"ncores": 1, "memory_gb": 4, "time_min": 30, "queue": "short"},
}

STAGE_OUTPUT_HINTS: dict[str, tuple[str, ...]] = {
    "conformers": ("molecule",),
    "build-hf": ("conformer", "site"),
    "preopt": ("species",),
    "dft-minima": ("species_preopt", "species"),
    "ts-search": ("reaction",),
    "irc": ("reaction_validated",),
    "sp": ("species_optimized", "species_preopt"),
    "thermo": ("calculation",),
    "kinetics": ("thermo",),
}


def _enabled_stages(pipeline_config: dict[str, Any]) -> list[dict[str, Any]]:
    return [s for s in pipeline_config.get("stages", []) if s.get("enabled", True)]


def _stage_name(stage_cfg: dict[str, Any]) -> str:
    return str(stage_cfg.get("name", "unknown")).replace("_", "-")


def _artifact_counts(run_dir: str | Path | None) -> dict[str, int]:
    if not run_dir:
        return {}
    try:
        manifest = read_manifest(latest_manifest_path(run_dir))
    except Exception:
        return {}
    counts: dict[str, int] = {}
    for a in manifest.artifacts:
        counts[a.artifact_type] = counts.get(a.artifact_type, 0) + 1
    return counts


def estimate_stage_item_count(stage_name: str, artifact_counts: dict[str, int]) -> int:
    hints = STAGE_OUTPUT_HINTS.get(stage_name, ())
    vals = [artifact_counts.get(h, 0) for h in hints if artifact_counts.get(h, 0)]
    if vals:
        return max(1, max(vals))
    if stage_name in {"ingest", "enrich", "detect-sites", "rank", "viz", "ops", "calibrate"}:
        return 1
    return 1


def _merge_resource(defaults: dict[str, Any], stage_cfg: dict[str, Any], profile: dict[str, Any]) -> dict[str, Any]:
    resources = dict(defaults)
    stage_overrides = profile.get("stages", {}).get(_stage_name(stage_cfg), {}) if profile else {}
    resources.update(stage_overrides)
    if "resources" in stage_cfg and isinstance(stage_cfg["resources"], dict):
        resources.update(stage_cfg["resources"])
    settings = stage_cfg.get("settings") or {}
    if "ncores" in settings:
        resources["ncores"] = settings["ncores"]
    if "memory_gb" in settings:
        resources["memory_gb"] = settings["memory_gb"]
    elif "memory_mb" in settings:
        resources["memory_gb"] = max(1, math.ceil(float(settings["memory_mb"]) / 1024))
    if "time_min" in settings:
        resources["time_min"] = settings["time_min"]
    return resources


def build_resource_plan(
    pipeline_config_path: str | Path,
    run_dir: str | Path | None = None,
    ops_config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    cfg = load_yaml(pipeline_config_path)
    ops_config = ops_config or {}
    profile = ops_config.get("resource_profile") or {}
    counts = _artifact_counts(run_dir)
    rows: list[dict[str, Any]] = []
    prev_job_id: str | None = None
    for i, stage_cfg in enumerate(_enabled_stages(cfg)):
        name = _stage_name(stage_cfg)
        defaults = DEFAULT_STAGE_RESOURCES.get(name, {"ncores": 1, "memory_gb": 4, "time_min": 60, "queue": "short"})
        res = _merge_resource(defaults, stage_cfg, profile)
        item_count = estimate_stage_item_count(name, counts)
        array_chunk = int(res.get("array_chunk", profile.get("array_chunk", 1)) or 1)
        n_array_tasks = int(math.ceil(item_count / max(1, array_chunk)))
        job_id = f"job_{i:02d}_{name.replace('-', '_')}"
        rows.append(
            {
                "job_id": job_id,
                "stage_index": i,
                "stage": name,
                "depends_on": prev_job_id,
                "run_id_placeholder": "${RUN_ID}",
                "pipeline_config": str(pipeline_config_path),
                "ncores": int(res.get("ncores", 1)),
                "memory_gb": int(res.get("memory_gb", 4)),
                "time_min": int(res.get("time_min", 60)),
                "time_slurm": minutes_to_slurm_time(int(res.get("time_min", 60))),
                "queue": res.get("queue", "short"),
                "item_count_estimate": item_count,
                "array_chunk": array_chunk,
                "n_array_tasks_estimate": n_array_tasks,
                "command": f"python -m hfauto.cli.main pipeline --config {pipeline_config_path} --run-id ${{RUN_ID}} --from {name} --to {name}",
            }
        )
        prev_job_id = job_id
    return pd.DataFrame(rows)


def minutes_to_slurm_time(minutes: int) -> str:
    minutes = max(1, int(minutes))
    days, rem = divmod(minutes, 24 * 60)
    hours, mins = divmod(rem, 60)
    if days:
        return f"{days}-{hours:02d}:{mins:02d}:00"
    return f"{hours:02d}:{mins:02d}:00"


def write_resource_plan(df: pd.DataFrame, out_dir: str | Path) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "resource_plan.csv"
    json_path = out / "resource_plan.json"
    df.to_csv(csv_path, index=False)
    json_path.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    return {"csv": csv_path, "json": json_path}
