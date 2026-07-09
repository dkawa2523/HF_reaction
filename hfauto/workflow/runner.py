from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.core.io import read_manifest, write_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.stages.base import StageContext
from hfauto.stages.registry import get_stage


def stage_dir(run_dir: Path, index: int, name: str) -> Path:
    clean = name.replace("_", "-")
    return run_dir / f"{index:02d}_{clean}"


def _slice_stages(stages: list[dict[str, Any]], from_stage: str | None, to_stage: str | None) -> list[dict[str, Any]]:
    enabled = [s for s in stages if s.get("enabled", True)]
    names = [s["name"] for s in enabled]
    start = names.index(from_stage) if from_stage in names else 0
    end = names.index(to_stage) + 1 if to_stage in names else len(enabled)
    return enabled[start:end]


def run_pipeline(
    config: dict[str, Any],
    run_id: str,
    from_stage: str | None = None,
    to_stage: str | None = None,
    start_manifest: str | Path | None = None,
) -> Path:
    run_root = Path(config.get("run_root", "runs"))
    run_dir = run_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    global_config = dict(config.get("global", {}))
    if "__config_path" in config:
        global_config["pipeline_config_path"] = config["__config_path"]
    current_manifest = read_manifest(start_manifest) if start_manifest else None
    current_manifest_path = Path(start_manifest) if start_manifest else None

    selected = _slice_stages(config.get("stages", []), from_stage, to_stage)
    if selected and selected[0]["name"] != "ingest" and current_manifest is None:
        try:
            current_manifest_path = latest_manifest_path(run_dir)
            current_manifest = read_manifest(current_manifest_path)
        except Exception as exc:
            raise ValueError("Partial pipeline requires --start-manifest or an existing run manifest.path") from exc

    # Preserve original stage index for stable directory names.
    all_enabled = [s for s in config.get("stages", []) if s.get("enabled", True)]
    stage_index = {id(s): i for i, s in enumerate(all_enabled)}

    for local_idx, stage_cfg in enumerate(selected):
        name = stage_cfg["name"]
        stage = get_stage(name)
        idx = stage_index.get(id(stage_cfg), local_idx)
        out_dir = stage_dir(run_dir, idx, name)
        context = StageContext(out_dir=out_dir, run_id=run_id, global_config=global_config)
        cfg = dict(stage_cfg)
        if name == "ingest":
            cfg["sdf"] = config["input"]["sdf"]
            manifest = stage.run(None, cfg, context)
        else:
            if current_manifest is None and current_manifest_path is not None:
                current_manifest = read_manifest(current_manifest_path)
            manifest = stage.run(current_manifest, cfg, context)
        manifest.metadata.setdefault("global_config", global_config)
        manifest.metadata.setdefault("pipeline_id", config.get("pipeline_id"))
        manifest.metadata.setdefault("run_mode", global_config.get("mode", config.get("mode", "development")))
        current_manifest = manifest
        current_manifest_path = write_manifest(manifest, out_dir)

    if current_manifest_path:
        (run_dir / "manifest.path").write_text(str(current_manifest_path), encoding="utf-8")
    return run_dir
