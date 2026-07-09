from __future__ import annotations

"""QCArchive/QCFractal export planning.

This module avoids requiring qcelemental/qcportal in normal installations. It
writes portable JSONL records that can be reviewed or converted to QCArchive
submissions in a production environment.
"""

from pathlib import Path
from typing import Any

from hfauto.core.io import ensure_dir, read_manifest, write_json, write_jsonl
from hfauto.reporting.html_report import latest_manifest_path


def _qca_record(artifact) -> dict[str, Any] | None:
    if artifact.artifact_type not in {"species", "species_preopt", "species_optimized", "calculation"}:
        return None
    data = artifact.data or {}
    method = artifact.method or {}
    return {
        "record_type": artifact.artifact_type,
        "artifact_id": artifact.artifact_id,
        "species_id": data.get("species_id") or data.get("source_species_id"),
        "mol_id": data.get("mol_id"),
        "site_id": data.get("site_id"),
        "state": data.get("state"),
        "hf_n": data.get("hf_n"),
        "charge": data.get("charge"),
        "multiplicity": data.get("multiplicity"),
        "xyz_path": data.get("xyz_path") or artifact.paths.get("xyz") or artifact.paths.get("final_xyz"),
        "engine": method.get("engine"),
        "method_id": method.get("method_id") or method.get("model"),
        "task": method.get("task"),
        "status": artifact.status.status,
        "fallback_dummy": bool(artifact.qc.get("fallback_dummy") or data.get("fallback_dummy")),
        "real_orca_executed": bool(artifact.qc.get("real_orca_executed") or data.get("real_orca_executed")),
    }


def write_qcarchive_export(run_dir: str | Path, out_dir: str | Path, dataset_name: str | None = None) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    manifest = read_manifest(latest_manifest_path(run_dir))
    records = [r for a in manifest.artifacts if (r := _qca_record(a)) is not None]
    jsonl = out / "qcarchive_records.jsonl"
    summary = out / "qcarchive_export_summary.json"
    submit = out / "qcarchive_submit_plan.py"
    write_jsonl(records, jsonl)
    dataset_name = dataset_name or f"hfauto_{manifest.run_id}"
    write_json(summary, {"schema_version": "hfauto.qcarchive_export.v1", "run_id": manifest.run_id, "dataset_name": dataset_name, "n_records": len(records), "requires": ["qcportal", "qcelemental"], "status": "export_plan_only"})
    submit.write_text(
        "\"\"\"Review-only QCArchive submission skeleton.\n"
        "Install qcportal/qcelemental and map hfauto records to your QCFractal server before use.\n\"\"\"\n"
        "from pathlib import Path\nimport json\n\n"
        f"DATASET_NAME = {dataset_name!r}\nRECORDS = Path({str(jsonl)!r})\n\n"
        "def main():\n"
        "    rows = [json.loads(line) for line in RECORDS.read_text().splitlines() if line.strip()]\n"
        "    print(f'Would submit {len(rows)} hfauto records to QCArchive dataset {DATASET_NAME}')\n"
        "    print('This is a skeleton; production submission requires site-specific credentials and molecule mapping.')\n\n"
        "if __name__ == '__main__':\n    main()\n",
        encoding="utf-8",
    )
    return {"qca_jsonl": jsonl, "qca_summary": summary, "qca_submit_plan": submit}
