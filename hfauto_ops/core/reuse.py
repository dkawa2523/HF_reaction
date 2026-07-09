from __future__ import annotations

"""Conservative calculation reuse planning.

The output is intentionally a plan, not automatic substitution. It identifies
successful prior artifacts with the same conservative signature and asks a user
or a later policy engine to decide whether reuse is acceptable.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, read_manifest, write_json
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.run_index import calculation_signature, artifact_record


def _safe_read_csv(path: str | Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def _artifact_rows(run_dir: str | Path) -> pd.DataFrame:
    man = read_manifest(latest_manifest_path(run_dir))
    rows = []
    for a in man.artifacts:
        if a.artifact_type not in {"calculation", "thermo", "kinetics", "reaction_validated", "reaction_path_validated"}:
            continue
        rec = artifact_record(a)
        rec["calculation_signature"] = calculation_signature(a)
        rec["run_dir"] = str(run_dir)
        rows.append(rec)
    return pd.DataFrame(rows)


def build_reuse_plan(target_run: str | Path, cache_roots: list[str | Path] | None = None) -> pd.DataFrame:
    target = _artifact_rows(target_run)
    cache_frames = [_artifact_rows(root) for root in (cache_roots or [])]
    cache = pd.concat(cache_frames, ignore_index=True, sort=False) if cache_frames else pd.DataFrame()
    if target.empty:
        return pd.DataFrame()
    if cache.empty:
        target["reuse_status"] = "no_cache"
        target["reuse_candidate_artifact_id"] = None
        target["reuse_candidate_run_dir"] = None
        target["review_action"] = "run_calculation"
        return target
    cache_success = cache[cache.get("status", "success") == "success"].copy() if "status" in cache.columns else cache
    rows: list[dict[str, Any]] = []
    for _, row in target.iterrows():
        sig = row.get("calculation_signature")
        hit = cache_success[cache_success["calculation_signature"] == sig] if "calculation_signature" in cache_success.columns else pd.DataFrame()
        rec = row.to_dict()
        if not hit.empty:
            h = hit.iloc[0]
            rec.update({
                "reuse_status": "hit",
                "reuse_candidate_artifact_id": h.get("artifact_id"),
                "reuse_candidate_run_dir": h.get("run_dir"),
                "review_action": "review_reuse_candidate",
            })
        else:
            rec.update({
                "reuse_status": "miss",
                "reuse_candidate_artifact_id": None,
                "reuse_candidate_run_dir": None,
                "review_action": "run_calculation",
            })
        rows.append(rec)
    return pd.DataFrame(rows)


def write_reuse_plan(target_run: str | Path, out_dir: str | Path, cache_roots: list[str | Path] | None = None, cache_runs: list[str | Path] | None = None) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    roots = cache_roots if cache_roots is not None else cache_runs
    df = build_reuse_plan(target_run, cache_roots=roots)
    csv = out / "reuse_plan.csv"
    json_path = out / "reuse_plan.json"
    policy = out / "reuse_decisions_template.yaml"
    df.to_csv(csv, index=False)
    write_json(json_path, {"schema_version": "hfauto.reuse_plan.v1", "n_candidates": int(len(df)), "rows": df.to_dict(orient="records")})
    policy.write_text("# Review and edit decisions before applying reuse.\n# decision: accept | reject | defer\nreuse_decisions: []\n", encoding="utf-8")
    return {"reuse_csv": csv, "reuse_json": json_path, "reuse_policy_template": policy}



def write_reuse_bundle(
    target_run: str | Path,
    out_dir: str | Path,
    artifact_job_plan_csv: str | Path | None = None,
    registry_db: str | Path | None = None,
    cache_runs: list[str | Path] | None = None,
) -> dict[str, Path]:
    """Phase 12 compatibility wrapper.

    Writes the conservative reuse plan plus a lightweight registry/review file.
    `artifact_job_plan_csv` and `registry_db` are accepted so production callers
    can keep a stable API; the current implementation does not mutate either.
    """
    out = ensure_dir(out_dir)
    paths = write_reuse_plan(target_run, out, cache_roots=cache_runs)
    review = out / "calculation_reuse_review.csv"
    try:
        df = pd.read_csv(paths["reuse_csv"]) if Path(paths["reuse_csv"]).exists() else pd.DataFrame()
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if artifact_job_plan_csv and Path(artifact_job_plan_csv).exists():
        try:
            jobs = pd.read_csv(artifact_job_plan_csv)
        except pd.errors.EmptyDataError:
            jobs = pd.DataFrame()
        except Exception:
            jobs = pd.DataFrame()
        if not jobs.empty:
            cols = [c for c in ["target_artifact_id", "artifact_job_id", "target_stage", "command"] if c in jobs.columns]
            if cols:
                if not df.empty and "target_artifact_id" in df.columns and "target_artifact_id" in cols:
                    df = df.merge(jobs[cols], on="target_artifact_id", how="outer")
                else:
                    df = jobs[cols]
    df.to_csv(review, index=False)
    registry = out / "calculation_registry_plan.json"
    write_json(registry, {
        "schema_version": "hfauto.reuse_registry_plan.v1",
        "registry_db": str(registry_db) if registry_db else None,
        "artifact_job_plan_csv": str(artifact_job_plan_csv) if artifact_job_plan_csv else None,
        "n_review_rows": int(len(df)),
        "note": "Review-only reuse planning; no scientific artifact substitution performed.",
    })
    paths.update({"reuse_plan": paths["reuse_csv"], "reuse_review_csv": review, "registry_plan": registry})
    return paths
