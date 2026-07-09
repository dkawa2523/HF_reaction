from __future__ import annotations

"""Deduplication-backed calculation reuse planning.

The module compares conservative calculation signatures across one or more run
roots.  It never replaces artifacts automatically; it emits reviewable reuse
candidates with QC flags so a production operator can approve reuse.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.run_index import artifact_record, calculation_signature

REUSE_TYPES = {"calculation", "thermo", "kinetics", "species_preopt", "species_optimized"}


def _manifest_for_run(run_dir: str | Path):
    return read_manifest(latest_manifest_path(run_dir))


def build_reuse_catalog(run_dirs: list[str | Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for run_dir in run_dirs:
        run = Path(run_dir)
        try:
            manifest = _manifest_for_run(run)
        except Exception:
            continue
        for artifact in manifest.artifacts:
            if artifact.artifact_type not in REUSE_TYPES:
                continue
            rec = artifact_record(artifact)
            rec.update(
                {
                    "run_id": manifest.run_id,
                    "run_dir": str(run),
                    "signature": calculation_signature(artifact),
                    "eligible_for_reuse": bool(
                        artifact.status.status == "success"
                        and not artifact.qc.get("fallback_dummy")
                        and not artifact.data.get("fallback_dummy")
                    ),
                }
            )
            rows.append(rec)
    return pd.DataFrame(rows)


def build_reuse_plan(target_run: str | Path, cache_roots: list[str | Path] | None = None) -> pd.DataFrame:
    runs = [Path(target_run)] + [Path(p) for p in (cache_roots or [])]
    catalog = build_reuse_catalog(runs)
    if catalog.empty or "signature" not in catalog.columns:
        return pd.DataFrame()
    target = catalog[catalog["run_dir"] == str(Path(target_run))].copy()
    pool = catalog[(catalog["run_dir"] != str(Path(target_run))) & (catalog["eligible_for_reuse"] == True)].copy()
    rows: list[dict[str, Any]] = []
    for _, t in target.iterrows():
        matches = pool[pool["signature"] == t["signature"]]
        for _, m in matches.iterrows():
            rows.append(
                {
                    "target_run_id": t.get("run_id"),
                    "target_artifact_id": t.get("artifact_id"),
                    "target_artifact_type": t.get("artifact_type"),
                    "target_status": t.get("status"),
                    "target_fallback_dummy": t.get("fallback_dummy"),
                    "source_run_id": m.get("run_id"),
                    "source_run_dir": m.get("run_dir"),
                    "source_artifact_id": m.get("artifact_id"),
                    "source_artifact_type": m.get("artifact_type"),
                    "signature": t.get("signature"),
                    "engine": t.get("engine") or m.get("engine"),
                    "method_id": t.get("method_id") or m.get("method_id"),
                    "species_id": t.get("species_id") or m.get("species_id"),
                    "reaction_id": t.get("reaction_id") or m.get("reaction_id"),
                    "reuse_action": "review_reuse_candidate",
                    "safety_note": "Conservative signature match only. Verify geometry, charge/multiplicity, method settings, and QC before substitution.",
                }
            )
    return pd.DataFrame(rows)


def write_reuse_plan(target_run: str | Path, out_dir: str | Path, cache_roots: list[str | Path] | None = None) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    catalog = build_reuse_catalog([Path(target_run)] + [Path(p) for p in (cache_roots or [])])
    plan = build_reuse_plan(target_run, cache_roots)
    catalog_path = out / "reuse_catalog.csv"
    plan_path = out / "reuse_plan.csv"
    json_path = out / "reuse_plan.json"
    catalog.to_csv(catalog_path, index=False)
    plan.to_csv(plan_path, index=False)
    json_path.write_text(plan.to_json(orient="records", indent=2), encoding="utf-8")
    return {"reuse_catalog": catalog_path, "reuse_plan": plan_path, "reuse_json": json_path}
