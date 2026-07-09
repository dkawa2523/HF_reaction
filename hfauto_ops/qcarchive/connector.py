from __future__ import annotations

"""QCArchive/QCFractal integration planning skeleton.

This module deliberately avoids importing qcelemental/qcportal unless production
systems add them.  It produces JSONL payloads and a connection/readiness report
that can be reviewed before pushing calculations to a QCArchive server.
"""

import json
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.run_index import artifact_record, calculation_signature

QCARCHIVE_TYPES = {"calculation", "species", "species_preopt", "species_optimized"}


def qca_payload_for_artifact(artifact, run_id: str) -> dict[str, Any]:
    data = artifact.data or {}
    method = artifact.method or {}
    return {
        "schema_name": "hfauto_qcarchive_payload",
        "schema_version": 1,
        "run_id": run_id,
        "artifact_id": artifact.artifact_id,
        "artifact_type": artifact.artifact_type,
        "signature": calculation_signature(artifact),
        "molecule_hint": {
            "mol_id": data.get("mol_id"),
            "species_id": data.get("species_id") or data.get("source_species_id"),
            "reaction_id": data.get("reaction_id"),
            "state": data.get("state"),
            "hf_n": data.get("hf_n"),
            "xyz_path": data.get("xyz_path") or artifact.paths.get("xyz") or artifact.paths.get("final_xyz"),
        },
        "driver": method.get("task") or method.get("driver") or "unknown",
        "model": {
            "engine": method.get("engine"),
            "method_id": method.get("method_id") or method.get("method") or method.get("model"),
            "basis": method.get("basis"),
            "settings": method.get("settings") or {},
        },
        "qc": artifact.qc,
        "status": artifact.status.model_dump(),
        "paths": artifact.paths,
        "note": "Review and convert to native QCSchema before submission to QCArchive/QCFractal.",
    }


def build_qcarchive_payloads(run_dir: str | Path) -> list[dict[str, Any]]:
    manifest = read_manifest(latest_manifest_path(run_dir))
    payloads = []
    for artifact in manifest.artifacts:
        if artifact.artifact_type in QCARCHIVE_TYPES:
            payloads.append(qca_payload_for_artifact(artifact, manifest.run_id))
    return payloads


def write_qcarchive_plan(run_dir: str | Path, out_dir: str | Path, dataset_name: str | None = None, server_url: str | None = None) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    payloads = build_qcarchive_payloads(run_dir)
    payload_path = out / "qcarchive_payloads.jsonl"
    with payload_path.open("w", encoding="utf-8") as f:
        for payload in payloads:
            f.write(json.dumps(payload, ensure_ascii=False) + "\n")
    rows = []
    for p in payloads:
        rows.append(
            {
                "artifact_id": p["artifact_id"],
                "artifact_type": p["artifact_type"],
                "signature": p["signature"],
                "driver": p["driver"],
                "engine": p["model"].get("engine"),
                "method_id": p["model"].get("method_id"),
                "species_id": p["molecule_hint"].get("species_id"),
                "reaction_id": p["molecule_hint"].get("reaction_id"),
                "status": p["status"].get("status"),
                "ready_for_native_qcschema": bool(p["molecule_hint"].get("xyz_path") and p["model"].get("method_id")),
            }
        )
    df = pd.DataFrame(rows)
    csv_path = out / "qcarchive_import_plan.csv"
    df.to_csv(csv_path, index=False)
    summary = {
        "dataset_name": dataset_name or "hfauto_dataset_review",
        "server_url": server_url,
        "n_payloads": len(payloads),
        "n_ready_for_native_qcschema": int(df.get("ready_for_native_qcschema", pd.Series(dtype=bool)).sum()) if not df.empty else 0,
        "external_qcarchive_executed": False,
        "production_ready": False,
        "note": "This is an integration plan. Install/configure qcportal and approve QCSchema conversion before production import.",
    }
    summary_path = out / "qcarchive_connection_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {"qcarchive_payloads": payload_path, "qcarchive_import_plan": csv_path, "qcarchive_summary": summary_path}
