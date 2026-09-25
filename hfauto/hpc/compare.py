from __future__ import annotations

"""Run and backend comparison utilities."""

from pathlib import Path

import pandas as pd

from hfauto.core.io import ensure_dir, read_manifest
from hfauto.reporting.html_report import latest_manifest_path


def _latest(run_dir: str | Path):
    return read_manifest(latest_manifest_path(Path(run_dir)))


def _artifact_table(manifest, artifact_types: set[str] | None = None) -> pd.DataFrame:
    rows = []
    for a in manifest.artifacts:
        if artifact_types and a.artifact_type not in artifact_types:
            continue
        row = {
            "run_id": manifest.run_id,
            "artifact_id": a.artifact_id,
            "artifact_type": a.artifact_type,
            "status": a.status.status,
            "engine": (a.method or {}).get("engine"),
            "method_id": (a.method or {}).get("method_id"),
            "task": (a.method or {}).get("task"),
        }
        for key in [
            "species_id", "reaction_id", "mol_id", "site_id", "hf_n", "state", "quality_tier",
            "electronic_energy_hartree", "gibbs_298K_hartree", "delta_G_assoc_kcal_mol",
            "delta_G_act_kcal_mol", "k_corrected_s-1", "scavenger_score", "activation_score",
        ]:
            if key in a.data:
                row[key] = a.data.get(key)
        rows.append(row)
    return pd.DataFrame(rows)


def compare_runs(run_a: str | Path, run_b: str | Path, out_dir: str | Path) -> dict[str, Path]:
    ma = _latest(run_a)
    mb = _latest(run_b)
    a = _artifact_table(ma)
    b = _artifact_table(mb)
    out = ensure_dir(out_dir)
    summary = {
        "run_a": ma.run_id,
        "run_b": mb.run_id,
        "n_artifacts_a": len(a),
        "n_artifacts_b": len(b),
        "n_failures_a": int((a.get("status") == "failed").sum()) if not a.empty else 0,
        "n_failures_b": int((b.get("status") == "failed").sum()) if not b.empty else 0,
    }
    common_cols = [c for c in ["artifact_id", "artifact_type", "species_id", "reaction_id", "mol_id", "site_id", "hf_n"] if c in a.columns and c in b.columns]
    merged = a.merge(b, on=common_cols, how="outer", suffixes=("_a", "_b"), indicator=True) if common_cols else pd.DataFrame()
    if not merged.empty:
        for metric in ["delta_G_assoc_kcal_mol", "delta_G_act_kcal_mol", "k_corrected_s-1", "scavenger_score", "activation_score"]:
            ca, cb = f"{metric}_a", f"{metric}_b"
            if ca in merged.columns and cb in merged.columns:
                merged[f"diff_{metric}"] = pd.to_numeric(merged[cb], errors="coerce") - pd.to_numeric(merged[ca], errors="coerce")
    summary_path = out / "run_comparison_summary.json"
    diff_path = out / "run_comparison.csv"
    pd.Series(summary).to_json(summary_path, force_ascii=False, indent=2)
    merged.to_csv(diff_path, index=False)
    return {"summary": summary_path, "csv": diff_path}


def compare_backends(run_dir: str | Path, out_dir: str | Path) -> dict[str, Path]:
    m = _latest(run_dir)
    df = _artifact_table(m, {"calculation"})
    out = ensure_dir(out_dir)
    if df.empty or "species_id" not in df.columns:
        result = pd.DataFrame()
    else:
        value_cols = [c for c in ["electronic_energy_hartree", "gibbs_298K_hartree"] if c in df.columns]
        keys = ["species_id", "task"] if "task" in df.columns else ["species_id"]
        rows = []
        for key, group in df.groupby(keys, dropna=False):
            if group["engine"].nunique(dropna=True) + group["method_id"].nunique(dropna=True) < 2:
                continue
            records = group.to_dict(orient="records")
            for i in range(len(records)):
                for j in range(i + 1, len(records)):
                    row = {"comparison_key": str(key), "artifact_a": records[i].get("artifact_id"), "artifact_b": records[j].get("artifact_id"), "engine_a": records[i].get("engine"), "engine_b": records[j].get("engine"), "method_a": records[i].get("method_id"), "method_b": records[j].get("method_id")}
                    for val in value_cols:
                        try:
                            row[f"diff_{val}"] = float(records[j].get(val)) - float(records[i].get(val))
                        except Exception:
                            row[f"diff_{val}"] = None
                    rows.append(row)
        result = pd.DataFrame(rows)
    path = out / "backend_comparison.csv"
    result.to_csv(path, index=False)
    return {"csv": path}
