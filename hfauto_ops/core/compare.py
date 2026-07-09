from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def _read_optional_csv(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def _table(run_dir: str | Path, name: str) -> pd.DataFrame:
    run = Path(run_dir)
    matches = sorted(run.glob(f"*/{name}.csv"))
    if not matches:
        return pd.DataFrame()
    return pd.read_csv(matches[-1])


def _numeric_delta(df: pd.DataFrame, suffix_a: str = "_a", suffix_b: str = "_b") -> pd.DataFrame:
    out = df.copy()
    for col in list(df.columns):
        if not col.endswith(suffix_a):
            continue
        base = col[: -len(suffix_a)]
        other = f"{base}{suffix_b}"
        if (
            other in df.columns
            and pd.api.types.is_numeric_dtype(df[col])
            and pd.api.types.is_numeric_dtype(df[other])
            and not pd.api.types.is_bool_dtype(df[col])
            and not pd.api.types.is_bool_dtype(df[other])
        ):
            out[f"delta_{base}"] = df[other] - df[col]
    return out


def compare_runs(run_a: str | Path, run_b: str | Path, out_dir: str | Path) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cand_a = _table(run_a, "candidate_summary")
    cand_b = _table(run_b, "candidate_summary")
    rxn_a = _table(run_a, "reaction_results")
    rxn_b = _table(run_b, "reaction_results")

    paths: dict[str, Path] = {}
    if not cand_a.empty and not cand_b.empty and "mol_id" in cand_a.columns and "mol_id" in cand_b.columns:
        merged = cand_a.merge(cand_b, on="mol_id", how="outer", suffixes=("_a", "_b"), indicator=True)
        merged = _numeric_delta(merged)
        p = out / "candidate_comparison.csv"
        merged.to_csv(p, index=False)
        paths["candidate_comparison"] = p
    else:
        p = out / "candidate_comparison.csv"
        pd.DataFrame().to_csv(p, index=False)
        paths["candidate_comparison"] = p

    keys = [k for k in ["reaction_id", "mol_id", "site_id", "hf_n"] if k in rxn_a.columns and k in rxn_b.columns]
    if not rxn_a.empty and not rxn_b.empty and keys:
        merged = rxn_a.merge(rxn_b, on=keys, how="outer", suffixes=("_a", "_b"), indicator=True)
        merged = _numeric_delta(merged)
        p = out / "reaction_comparison.csv"
        merged.to_csv(p, index=False)
        paths["reaction_comparison"] = p
    else:
        p = out / "reaction_comparison.csv"
        pd.DataFrame().to_csv(p, index=False)
        paths["reaction_comparison"] = p

    summary = {
        "run_a": str(run_a),
        "run_b": str(run_b),
        "candidate_rows_a": int(len(cand_a)),
        "candidate_rows_b": int(len(cand_b)),
        "reaction_rows_a": int(len(rxn_a)),
        "reaction_rows_b": int(len(rxn_b)),
    }
    summary_path = out / "run_comparison_summary.json"
    summary_path.write_text(__import__("json").dumps(summary, indent=2), encoding="utf-8")
    paths["summary"] = summary_path
    html = out / "run_comparison_report.html"
    html.write_text(
        "<html><body><h1>hfauto run comparison</h1>"
        f"<p>Run A: {run_a}</p><p>Run B: {run_b}</p>"
        f"<p>Candidate rows: {summary['candidate_rows_a']} → {summary['candidate_rows_b']}</p>"
        f"<p>Reaction rows: {summary['reaction_rows_a']} → {summary['reaction_rows_b']}</p>"
        "<ul><li><a href='candidate_comparison.csv'>candidate_comparison.csv</a></li>"
        "<li><a href='reaction_comparison.csv'>reaction_comparison.csv</a></li></ul>"
        "</body></html>",
        encoding="utf-8",
    )
    paths["html"] = html
    return paths
