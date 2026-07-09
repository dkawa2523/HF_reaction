from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto_ops.core.run_index import artifact_record


def _calculation_rows(run_dir: str | Path) -> pd.DataFrame:
    m = read_manifest(latest_manifest_path(run_dir))
    rows = []
    for a in m.artifacts:
        if a.artifact_type != "calculation":
            continue
        rec = artifact_record(a)
        rec.update({k: a.data.get(k) for k in ["electronic_energy_hartree", "gibbs_298K_hartree", "n_imag", "hf_stretch_cm1"]})
        rows.append(rec)
    return pd.DataFrame(rows)


def build_backend_comparison(run_dir: str | Path) -> pd.DataFrame:
    df = _calculation_rows(run_dir)
    if df.empty or "species_id" not in df.columns:
        return pd.DataFrame()
    out_rows: list[dict[str, Any]] = []
    keys = ["species_id", "task"] if "task" in df.columns else ["species_id"]
    for key, group in df.groupby(keys, dropna=False):
        if len(group) < 2:
            continue
        records = group.to_dict(orient="records")
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                a, b = records[i], records[j]
                if (a.get("engine"), a.get("method_id")) == (b.get("engine"), b.get("method_id")):
                    continue
                row = {"comparison_key": str(key), "artifact_a": a.get("artifact_id"), "artifact_b": b.get("artifact_id"), "engine_a": a.get("engine"), "engine_b": b.get("engine"), "method_a": a.get("method_id"), "method_b": b.get("method_id"), "quality_tier_a": a.get("quality_tier"), "quality_tier_b": b.get("quality_tier")}
                for metric in ["electronic_energy_hartree", "gibbs_298K_hartree", "hf_stretch_cm1"]:
                    try:
                        row[f"delta_{metric}"] = float(b.get(metric)) - float(a.get(metric))
                    except Exception:
                        row[f"delta_{metric}"] = None
                    if row.get(f"delta_{metric}") is not None and metric.endswith("hartree"):
                        row[f"delta_{metric}_kcal_mol"] = row[f"delta_{metric}"] * 627.509474
                out_rows.append(row)
    return pd.DataFrame(out_rows)


def write_backend_comparison(run_dir: str | Path, out_dir: str | Path) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    raw = _calculation_rows(run_dir)
    df = build_backend_comparison(run_dir)
    raw_csv = out / "backend_calculations.csv"
    csv = out / "backend_comparison.csv"
    html = out / "backend_comparison_dashboard.html"
    raw.to_csv(raw_csv, index=False)
    df.to_csv(csv, index=False)
    html.write_text("<html><body><h1>Backend comparison</h1>" + (df.head(100).to_html(index=False) if not df.empty else "<p>No comparable backend pairs.</p>") + "</body></html>", encoding="utf-8")
    return {"backend_raw_csv": raw_csv, "backend_comparison_csv": csv, "backend_comparison_html": html, "raw_csv": raw_csv, "comparison_csv": csv, "dashboard_html": html}


def compare_backends(run_dir: str | Path, out_dir: str | Path) -> dict[str, Path]:
    return write_backend_comparison(run_dir, out_dir)
