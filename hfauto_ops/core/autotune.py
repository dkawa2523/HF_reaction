from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from hfauto.core.io import ensure_dir, write_json


def _read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def collect_history(run_dirs: list[str | Path]) -> pd.DataFrame:
    frames = []
    for rd in run_dirs:
        run = Path(rd)
        for rel in ["16_ops/job_history.csv", "16_ops/scheduler_status.csv", "16_ops/resource_plan.csv", "15_hpc-plan/job_plan.csv", "resource_plan.csv"]:
            df = _read_optional(run / rel)
            if not df.empty:
                df["source_run_dir"] = str(run)
                frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def propose_resources(history: pd.DataFrame, safety_factor: float = 1.25) -> dict[str, Any]:
    proposals: dict[str, Any] = {"resource_profile": {"stages": {}}}
    if history.empty or "stage" not in history.columns:
        return proposals
    for stage, grp in history.groupby("stage"):
        rec: dict[str, Any] = {}
        if "elapsed_min" in grp.columns:
            vals = pd.to_numeric(grp["elapsed_min"], errors="coerce").dropna()
            if not vals.empty:
                rec["time_min"] = int(max(10, vals.quantile(0.90) * safety_factor))
        elif "time_min" in grp.columns:
            vals = pd.to_numeric(grp["time_min"], errors="coerce").dropna()
            if not vals.empty:
                rec["time_min"] = int(max(10, vals.median()))
        if "max_rss_gb" in grp.columns:
            vals = pd.to_numeric(grp["max_rss_gb"], errors="coerce").dropna()
            if not vals.empty:
                rec["memory_gb"] = int(max(1, vals.quantile(0.90) * safety_factor + 0.999))
        elif "memory_gb" in grp.columns:
            vals = pd.to_numeric(grp["memory_gb"], errors="coerce").dropna()
            if not vals.empty:
                rec["memory_gb"] = int(max(1, vals.median()))
        if "ncores" in grp.columns:
            vals = pd.to_numeric(grp["ncores"], errors="coerce").dropna()
            if not vals.empty:
                rec["ncores"] = int(max(1, vals.median()))
        if rec:
            proposals["resource_profile"]["stages"][str(stage)] = rec
    return proposals


def build_autotune_table(job_history_csv: str | Path | None, resource_plan_csv: str | Path | None = None) -> pd.DataFrame:
    history = _read_optional(Path(job_history_csv)) if job_history_csv else pd.DataFrame()
    proposal = propose_resources(history)
    rows = []
    for stage, rec in proposal.get("resource_profile", {}).get("stages", {}).items():
        rows.append({"stage": stage, **rec})
    base = pd.DataFrame(rows)
    if resource_plan_csv and Path(resource_plan_csv).exists():
        plan = _read_optional(Path(resource_plan_csv))
        if not plan.empty and "stage" in plan.columns:
            base = plan[[c for c in ["stage", "time_min", "memory_gb", "ncores"] if c in plan.columns]].drop_duplicates().merge(base, on="stage", how="left", suffixes=("_current", "_recommended"))
    return base


def write_autotune_report(run_dirs_or_out_dir, out_dir: str | Path | None = None, safety_factor: float = 1.25, job_history_csv: str | Path | None = None, resource_plan_csv: str | Path | None = None) -> dict[str, Path]:
    # Backward compatible: write_autotune_report(out_dir, job_history_csv=..., resource_plan_csv=...)
    if out_dir is None:
        out = ensure_dir(run_dirs_or_out_dir)
        history = _read_optional(Path(job_history_csv)) if job_history_csv else pd.DataFrame()
    else:
        out = ensure_dir(out_dir)
        dirs = list(run_dirs_or_out_dir) if isinstance(run_dirs_or_out_dir, (list, tuple)) else [run_dirs_or_out_dir]
        history = collect_history(dirs)
    proposal = propose_resources(history, safety_factor=safety_factor)
    hist_csv = out / "resource_history.csv"
    prop_yaml = out / "resource_autotune_proposal.yaml"
    prop_json = out / "resource_autotune_proposal.json"
    csv = out / "resource_autotune_recommendations.csv"
    html = out / "resource_autotune_report.html"
    history.to_csv(hist_csv, index=False)
    rows = [{"stage": stage, **rec} for stage, rec in proposal.get("resource_profile", {}).get("stages", {}).items()]
    pd.DataFrame(rows).to_csv(csv, index=False)
    prop_yaml.write_text(yaml.safe_dump(proposal, sort_keys=False), encoding="utf-8")
    write_json(prop_json, proposal)
    body = ["<html><head><meta charset='utf-8'><title>Resource auto-tune</title></head><body>", "<h1>Resource auto-tune proposal</h1>"]
    body.append("<p>Transparent median/p90 based proposal. Review before production use.</p>")
    body.append(f"<p>History rows: {len(history)}</p>")
    if rows:
        body.append(pd.DataFrame(rows).to_html(index=False, escape=False))
    else:
        body.append("<p>No history rows; defaults retained.</p>")
    body.append("<h2>YAML proposal</h2><pre>" + prop_yaml.read_text(encoding="utf-8") + "</pre></body></html>")
    html.write_text("\n".join(body), encoding="utf-8")
    return {"history_csv": hist_csv, "autotune_csv": csv, "proposal_yaml": prop_yaml, "proposal_json": prop_json, "autotune_html": html, "report_html": html}


def build_resource_autotune_suggestions(run_dir: str | Path, pipeline_config: str | Path | None = None, ops_config: dict | None = None) -> pd.DataFrame:
    from hfauto_ops.core.resource import DEFAULT_STAGE_RESOURCES, build_resource_plan
    rows = []
    plan = build_resource_plan(pipeline_config, run_dir=run_dir, ops_config=ops_config or {}) if pipeline_config else pd.DataFrame()
    if plan.empty:
        for stage, res in DEFAULT_STAGE_RESOURCES.items():
            rows.append({
                "stage": stage,
                "current_ncores": res["ncores"],
                "current_memory_gb": res["memory_gb"],
                "current_time_min": res["time_min"],
                "suggested_ncores": res["ncores"],
                "suggested_memory_gb": res["memory_gb"],
                "suggested_time_min": res["time_min"],
                "item_count_estimate": 1,
                "suggestion": "keep",
                "note": "no pipeline_config provided",
            })
        return pd.DataFrame(rows)
    for _, row in plan.iterrows():
        stage = str(row.get("stage"))
        item_count = int(row.get("item_count_estimate") or 1)
        time_min = int(row.get("time_min") or 60)
        suggestion = "keep"
        note = "default"
        if item_count > 50 and stage in {"preopt", "dft-minima", "ts-search", "sp"}:
            suggestion = "use_artifact_job_array"
            note = "large task count; prefer artifact-level arrays"
        if stage == "ts-search" and time_min < 2880:
            suggestion = "increase_walltime"
            time_min = 2880
            note = "TS search often needs long walltime"
        rows.append({
            "stage": stage,
            "current_ncores": int(row.get("ncores") or 1),
            "current_memory_gb": int(row.get("memory_gb") or 4),
            "current_time_min": int(row.get("time_min") or 60),
            "suggested_ncores": int(row.get("ncores") or 1),
            "suggested_memory_gb": int(row.get("memory_gb") or 4),
            "suggested_time_min": time_min,
            "item_count_estimate": item_count,
            "suggestion": suggestion,
            "note": note,
        })
    return pd.DataFrame(rows)


def write_resource_autotune(df: pd.DataFrame, out_dir: str | Path) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    csv = out / "resource_autotune_suggestions.csv"
    html = out / "resource_autotune_report.html"
    json_path = out / "resource_autotune_suggestions.json"
    df.to_csv(csv, index=False)
    write_json(json_path, {"schema_version": "hfauto.resource_autotune.suggestions.v1", "n_rows": len(df), "rows": df.to_dict(orient="records")})
    html.write_text("<html><body><h1>Resource auto-tune suggestions</h1>" + (df.to_html(index=False) if not df.empty else "<p>No suggestions.</p>") + "</body></html>", encoding="utf-8")
    return {"resource_autotune_csv": csv, "resource_autotune_json": json_path, "resource_autotune_html": html}
