from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path

ARTIFACT_TO_STAGE = {
    "molecule": "ingest",
    "molecule_enriched": "enrich",
    "site": "detect-sites",
    "conformer": "conformers",
    "species": "build-hf",
    "species_preopt": "preopt",
    "calculation": "dft-minima",
    "ts_path": "ts-search",
    "reaction_validated": "ts-search",
    "irc": "irc",
    "irc_attempt": "irc",
    "thermo": "thermo",
    "kinetics": "kinetics",
    "ranking": "rank",
    "visualization_bundle": "viz",
}


def infer_retry_stage(artifact_type: str, category: str | None, fallback: str | None) -> str:
    if fallback:
        f = fallback.lower()
        for key in ["scan", "optts", "neb", "ts"]:
            if key in f:
                return "ts-search"
        if "irc" in f:
            return "irc"
        if "orca" in f or "dft" in f:
            return "dft-minima"
        if "xtb" in f or "preopt" in f:
            return "preopt"
    if category:
        c = category.lower()
        if "scf" in c or "opt" in c or "freq" in c:
            return "dft-minima"
        if "ts" in c or "neb" in c:
            return "ts-search"
        if "irc" in c:
            return "irc"
        if "visual" in c:
            return "viz"
    return ARTIFACT_TO_STAGE.get(artifact_type, "unknown")


def build_retry_dataframe(run_dir: str | Path, pipeline_config: str | Path | None = None) -> pd.DataFrame:
    run = Path(run_dir)
    manifest = read_manifest(latest_manifest_path(run))
    rows: list[dict[str, Any]] = []
    for a in manifest.artifacts:
        include = a.status.status == "failed"
        # Soft retry targets: fallback calculations are scientifically weak and should be rerun in production.
        if not include and (a.qc.get("fallback_dummy") or a.data.get("fallback_dummy")):
            include = True
        if not include:
            continue
        retry_stage = infer_retry_stage(a.artifact_type, a.status.category, a.status.recommended_fallback)
        cmd = None
        if pipeline_config and retry_stage != "unknown":
            cmd = f"python -m hfauto.cli.main pipeline --config {pipeline_config} --run-id {manifest.run_id} --from {retry_stage} --to {retry_stage}"
        rows.append(
            {
                "artifact_id": a.artifact_id,
                "artifact_type": a.artifact_type,
                "status": a.status.status,
                "category": a.status.category,
                "reason": a.status.reason,
                "recommended_fallback": a.status.recommended_fallback,
                "retry_stage": retry_stage,
                "mol_id": a.data.get("mol_id"),
                "site_id": a.data.get("site_id"),
                "reaction_id": a.data.get("reaction_id"),
                "species_id": a.data.get("species_id") or a.data.get("source_species_id"),
                "is_fallback_dummy": bool(a.qc.get("fallback_dummy") or a.data.get("fallback_dummy")),
                "retry_command": cmd,
            }
        )
    return pd.DataFrame(rows)


def write_retry_plan(run_dir: str | Path, out_dir: str | Path, pipeline_config: str | Path | None = None) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = build_retry_dataframe(run_dir, pipeline_config)
    csv = out / "retry_plan.csv"
    json = out / "retry_plan.json"
    shell = out / "retry_commands.sh"
    df.to_csv(csv, index=False)
    json.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    commands = [c for c in df.get("retry_command", pd.Series(dtype=str)).dropna().unique() if str(c).strip()]
    shell.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands) + "\n", encoding="utf-8")
    shell.chmod(0o755)
    return {"csv": csv, "json": json, "shell": shell}
