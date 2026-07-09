from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.schemas.manifest import Manifest
from hfauto_ops.core.run_index import failure_rows

DEFAULT_RETRY_POLICY: dict[str, dict[str, Any]] = {
    "dft_scf_failed": {"stage": "dft-minima", "priority": 90, "max_retries": 2, "action": "tighten_scf_or_change_initial_guess"},
    "dft_opt_failed": {"stage": "dft-minima", "priority": 85, "max_retries": 2, "action": "restart_from_last_geometry_or_loosen_then_tighten"},
    "freq_imaginary_minimum": {"stage": "dft-minima", "priority": 80, "max_retries": 1, "action": "reoptimize_with_tighter_thresholds"},
    "neb_failed": {"stage": "ts-search", "priority": 95, "max_retries": 2, "action": "fallback_scan_optts_or_zoom_neb"},
    "ts_search_failed": {"stage": "ts-search", "priority": 95, "max_retries": 2, "action": "try_next_ts_backend"},
    "ts_wrong_mode": {"stage": "ts-search", "priority": 90, "max_retries": 1, "action": "constrained_scan_near_expected_coordinate"},
    "irc_failed": {"stage": "irc", "priority": 88, "max_retries": 1, "action": "shorter_step_or_restart_from_ts"},
    "endpoint_mismatch": {"stage": "build-hf", "priority": 100, "max_retries": 1, "action": "rebuild_endpoints_and_atom_order"},
    "missing_input": {"stage": "inspect", "priority": 100, "max_retries": 0, "action": "run_upstream_stage"},
    "xtb_opt_failed": {"stage": "preopt", "priority": 75, "max_retries": 2, "action": "try_rdkit_or_xtb_relaxed_settings"},
    "unknown": {"stage": "inspect", "priority": 10, "max_retries": 0, "action": "manual_inspection"},
}


def build_retry_plan(manifest: Manifest, run_dir: str | Path, config: dict[str, Any] | None = None) -> pd.DataFrame:
    cfg = config or {}
    policy = {**DEFAULT_RETRY_POLICY, **(cfg.get("retry_policy", {}) or {})}
    rows = []
    for idx, row in enumerate(failure_rows(manifest)):
        cat = row.get("category") or "unknown"
        spec = policy.get(cat, policy["unknown"])
        target_stage = row.get("recommended_fallback") or spec.get("stage", "inspect")
        if target_stage.startswith("run_") or "before" in str(target_stage):
            target_stage = spec.get("stage", "inspect")
        command = "manual_inspection_required"
        if target_stage not in {"inspect", "manual"}:
            command = (
                f"hfauto pipeline --config configs/pipelines/phase9_hpc_template.yaml "
                f"--run-id {manifest.run_id} --from {target_stage} --to {target_stage} "
                f"--start-manifest $(cat {Path(run_dir) / 'manifest.path'})"
            )
        rows.append(
            {
                "retry_id": f"retry_{idx:04d}",
                "artifact_id": row["artifact_id"],
                "artifact_type": row["artifact_type"],
                "category": cat,
                "reason": row.get("reason", ""),
                "target_stage": target_stage,
                "priority": int(spec.get("priority", 10)),
                "max_retries": int(spec.get("max_retries", 0)),
                "action": spec.get("action", "manual_inspection"),
                "recoverable": row.get("recoverable", True),
                "command": command,
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["priority", "retry_id"], ascending=[False, True])
    return df


def write_retry_shell(retry_plan: pd.DataFrame, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "# Auto-generated retry commands. Review before running.", ""]
    if retry_plan.empty:
        lines.append("echo 'No retryable failures found.'")
    else:
        for _, row in retry_plan.iterrows():
            lines.append(f"# {row['retry_id']} {row['category']} {row['artifact_id']} action={row['action']}")
            cmd = str(row["command"])
            if cmd == "manual_inspection_required":
                lines.append(f"echo 'Manual inspection required for {row['artifact_id']}'")
            else:
                lines.append(cmd)
            lines.append("")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    p.chmod(0o755)
    return p
