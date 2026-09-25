from __future__ import annotations

from typing import Any

from hfauto.core.config import load_yaml
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto_ops.bundle import build_operations_bundle, operations_artifact


class OpsStage(Stage):
    """Read-only operations/HPC planning stage.

    This stage creates indexes, resource plans, artifact-level job plans,
    scheduler templates, retry/reuse/QCArchive plans, and an operations HTML
    report. It does not modify scientific artifacts.
    """

    name = "ops"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out = manifest.carry_forward(self.name)
        run_dir = config.get("run_dir") or str(context.out_dir.parent)
        pipeline_config = config.get("pipeline_config") or config.get("pipeline_config_path")
        # Accept both flat ops config and stage.settings for pipeline readability.
        ops_config: dict[str, Any] = {}
        if isinstance(config.get("ops_config"), dict):
            ops_config.update(config["ops_config"])
        elif isinstance(config.get("ops_config"), str):
            ops_config.update(load_yaml(config["ops_config"]))
        for key, value in config.items():
            if key not in {"name", "enabled", "pipeline_config", "pipeline_config_path", "ops_config", "settings"}:
                ops_config[key] = value
        if isinstance(config.get("settings"), dict):
            ops_config.update(config["settings"])
        paths = build_operations_bundle(run_dir, context.out_dir, pipeline_config=pipeline_config, ops_config=ops_config)
        out.add_artifact(operations_artifact(run_dir, context.out_dir, paths))
        return out
