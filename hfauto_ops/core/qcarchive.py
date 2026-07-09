from __future__ import annotations

"""QCArchive/QCFractal integration boundary.

This module intentionally does not require qcportal at import time.  It exports a
reviewable submit payload and, if qcportal is installed and explicitly enabled,
can be extended to submit records.  The default behavior is offline-safe.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, write_json, write_jsonl


def build_qcarchive_payload(artifact_job_plan_csv: str | Path) -> list[dict[str, Any]]:
    try:
        df = pd.read_csv(artifact_job_plan_csv)
    except Exception:
        df = pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        rows.append(
            {
                "schema_name": "hfauto_qcarchive_submit_v1",
                "job_id": row.get("job_id"),
                "stage": row.get("stage"),
                "artifact_id": row.get("artifact_id"),
                "species_id": row.get("species_id"),
                "reaction_id": row.get("reaction_id"),
                "engine": row.get("engine"),
                "method_id": row.get("method_id"),
                "cache_key": row.get("cache_key"),
                "command": row.get("command"),
                "note": "Review and map to QCArchive dataset/specification before live submit.",
            }
        )
    return rows


def write_qcarchive_bundle(artifact_job_plan_csv: str | Path, out_dir: str | Path, allow_submit: bool = False) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    payload = build_qcarchive_payload(artifact_job_plan_csv)
    jsonl = out / "qcarchive_submit_payload.jsonl"
    summary_path = out / "qcarchive_summary.json"
    write_jsonl(payload, jsonl)
    summary = {
        "schema_version": "hfauto.qcarchive_bundle.v1",
        "n_records": len(payload),
        "allow_submit": bool(allow_submit),
        "submitted": False,
        "status": "offline_payload_only" if not allow_submit else "submit_not_implemented_without_qcportal_config",
        "notes": [
            "No live submission is performed by default.",
            "Use this payload to review dataset/specification mapping before QCArchive/QCFractal integration.",
        ],
    }
    write_json(summary_path, summary)
    return {"payload_jsonl": jsonl, "summary": summary_path}
