from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.hashing import sha256_text
from hfauto.core.io import read_manifest, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.reporting.html_report import latest_manifest_path


def load_latest_manifest(run_dir: str | Path):
    """Compatibility API returning the latest manifest and its path."""
    path = latest_manifest_path(Path(run_dir))
    return read_manifest(path), path


def latest_table_path(
    manifest,
    artifact_id: str,
    preferred: tuple[str, ...] = ("parquet", "csv", "json", "jsonl"),
) -> Path | None:
    """Resolve a table path from the latest revision of an artifact."""
    artifact = manifest.find(artifact_id)
    if artifact is None:
        return None
    for key in preferred:
        raw = artifact.paths.get(key)
        if raw:
            return Path(raw)
    return Path(next(iter(artifact.paths.values()))) if artifact.paths else None


def read_table(path: str | Path | None) -> pd.DataFrame:
    """Read a supported table, returning an empty frame for missing outputs."""
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    if p.suffix.lower() in {".json", ".jsonl"}:
        return pd.read_json(p, lines=p.suffix.lower() == ".jsonl")
    return pd.read_csv(p)


def failure_rows(manifest) -> list[dict[str, Any]]:
    """Return only unresolved failures from the latest artifact revisions."""
    return [
        {
            "artifact_id": artifact.artifact_id,
            "artifact_type": artifact.artifact_type,
            "category": artifact.status.category,
            "reason": artifact.status.reason,
            "recoverable": artifact.status.recoverable,
            "recommended_fallback": artifact.status.recommended_fallback,
            "stage": (artifact.method or {}).get("stage") or artifact.data.get("stage"),
        }
        for artifact in manifest.latest_artifacts()
        if artifact.status.status == "failed"
    ]


def run_summary(manifest) -> dict[str, Any]:
    latest = manifest.latest_artifacts()
    return {
        "run_id": manifest.run_id,
        "latest_stage": manifest.stage,
        "n_artifacts": len(latest),
        "n_failures": sum(a.status.status == "failed" for a in latest),
        "n_fallback_dummy": sum(bool(a.qc.get("fallback_dummy") or a.data.get("fallback_dummy")) for a in latest),
    }


def _safe_get(d: dict[str, Any] | None, *keys: str, default: Any = None) -> Any:
    cur: Any = d or {}
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def artifact_record(artifact: Artifact) -> dict[str, Any]:
    method = artifact.method or {}
    data = artifact.data or {}
    status = artifact.status
    paths = artifact.paths or {}
    return {
        "artifact_id": artifact.artifact_id,
        "artifact_type": artifact.artifact_type,
        "status": status.status,
        "failure_category": status.category,
        "failure_reason": status.reason,
        "recommended_fallback": status.recommended_fallback,
        "recoverable": status.recoverable,
        "engine": method.get("engine"),
        "method_id": method.get("method_id") or method.get("model") or method.get("method"),
        "task": method.get("task"),
        "stage": method.get("stage") or artifact.provenance.get("stage"),
        "mol_id": data.get("mol_id"),
        "site_id": data.get("site_id"),
        "species_id": data.get("species_id") or data.get("source_species_id"),
        "reaction_id": data.get("reaction_id"),
        "state": data.get("state"),
        "hf_n": data.get("hf_n"),
        "quality_tier": data.get("quality_tier") or artifact.qc.get("quality_tier"),
        "confidence_score": data.get("confidence_score") or artifact.qc.get("confidence_score"),
        "fallback_dummy": bool(artifact.qc.get("fallback_dummy") or data.get("fallback_dummy")),
        "real_orca_executed": bool(artifact.qc.get("real_orca_executed") or data.get("real_orca_executed")),
        "has_paths": bool(paths),
        "n_paths": len(paths),
        "primary_path": next(iter(paths.values()), None) if paths else None,
        "parents": ";".join(artifact.parents),
    }


def build_artifact_dataframe(manifest_path: str | Path) -> pd.DataFrame:
    manifest = read_manifest(manifest_path)
    rows = [artifact_record(a) for a in manifest.artifacts]
    return pd.DataFrame(rows)


def stage_manifests(run_dir: str | Path) -> list[Path]:
    run = Path(run_dir)
    return sorted(run.glob("[0-9][0-9]_*/*manifest.json"))


def build_stage_dataframe(run_dir: str | Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for p in stage_manifests(run_dir):
        try:
            man = read_manifest(p)
        except Exception as exc:
            rows.append({"stage_dir": str(p.parent), "stage": p.parent.name, "manifest": str(p), "read_ok": False, "error": str(exc)})
            continue
        counts: dict[str, int] = {}
        failures = 0
        fallback = 0
        for a in man.artifacts:
            counts[a.artifact_type] = counts.get(a.artifact_type, 0) + 1
            failures += int(a.status.status == "failed")
            fallback += int(bool(a.qc.get("fallback_dummy") or a.data.get("fallback_dummy")))
        rows.append(
            {
                "stage": man.stage,
                "stage_dir": str(p.parent),
                "manifest": str(p),
                "manifest_id": man.manifest_id,
                "created_at": man.created_at,
                "n_artifacts": len(man.artifacts),
                "n_failures": failures,
                "n_fallback_dummy": fallback,
                "artifact_type_counts_json": json.dumps(counts, sort_keys=True),
                "read_ok": True,
            }
        )
    return pd.DataFrame(rows)


def calculation_signature(artifact: Artifact) -> str:
    """Stable fingerprint for duplicate/reuse checks.

    The signature intentionally excludes transient paths and timestamps and focuses on
    chemical identity + computational model + task. It is conservative: two artifacts
    with the same signature should be safe to inspect for possible reuse, not blindly
    substituted without geometry/QC checks.
    """
    data = artifact.data or {}
    method = artifact.method or {}
    payload = {
        "artifact_type": artifact.artifact_type,
        "task": method.get("task"),
        "engine": method.get("engine"),
        "method_id": method.get("method_id") or method.get("model") or method.get("method"),
        "basis": method.get("basis"),
        "settings": method.get("settings"),
        "mol_id": data.get("mol_id"),
        "site_id": data.get("site_id"),
        "species_id": data.get("species_id") or data.get("source_species_id"),
        "reaction_id": data.get("reaction_id"),
        "state": data.get("state"),
        "hf_n": data.get("hf_n"),
        "charge": data.get("charge"),
        "multiplicity": data.get("multiplicity"),
    }
    return sha256_text(json.dumps(payload, sort_keys=True, default=str))


def duplicate_dataframe(manifest_path: str | Path) -> pd.DataFrame:
    manifest = read_manifest(manifest_path)
    rows: list[dict[str, Any]] = []
    for a in manifest.artifacts:
        if a.artifact_type not in {"calculation", "thermo", "kinetics"}:
            continue
        rows.append({**artifact_record(a), "signature": calculation_signature(a)})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    counts = df.groupby("signature")["artifact_id"].transform("count")
    return df[counts > 1].sort_values(["signature", "artifact_id"]).reset_index(drop=True)


def write_index_bundle(run_dir: str | Path, out_dir: str | Path) -> dict[str, Path]:
    run = Path(run_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = latest_manifest_path(run)
    artifacts = build_artifact_dataframe(manifest_path)
    stages = build_stage_dataframe(run)
    duplicates = duplicate_dataframe(manifest_path)

    artifact_csv = out / "artifact_index.csv"
    stage_csv = out / "stage_index.csv"
    duplicate_csv = out / "duplicate_candidates.csv"
    artifact_jsonl = out / "artifact_index.jsonl"
    sqlite_path = out / "ops_index.sqlite"

    artifacts.to_csv(artifact_csv, index=False)
    stages.to_csv(stage_csv, index=False)
    duplicates.to_csv(duplicate_csv, index=False)
    write_jsonl(artifacts.to_dict(orient="records"), artifact_jsonl)

    with sqlite3.connect(sqlite_path) as con:
        artifacts.to_sql("artifact_index", con, if_exists="replace", index=False)
        stages.to_sql("stage_index", con, if_exists="replace", index=False)
        duplicates.to_sql("duplicate_candidates", con, if_exists="replace", index=False)

    summary = {
        "run_dir": str(run),
        "latest_manifest": str(manifest_path),
        "n_artifacts": len(artifacts),
        "n_stage_manifests": len(stages),
        "n_failures": int((artifacts.get("status") == "failed").sum()) if not artifacts.empty else 0,
        "n_fallback_dummy": int(artifacts.get("fallback_dummy", pd.Series(dtype=bool)).fillna(False).sum()) if not artifacts.empty else 0,
        "n_duplicate_candidate_rows": len(duplicates),
    }
    summary_path = write_json(out / "ops_index_summary.json", summary)
    return {
        "artifact_csv": artifact_csv,
        "artifact_jsonl": artifact_jsonl,
        "stage_csv": stage_csv,
        "duplicate_csv": duplicate_csv,
        "sqlite": sqlite_path,
        "summary": summary_path,
    }
