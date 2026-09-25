from __future__ import annotations

from pathlib import Path

import pandas as pd

from hfauto_ops.core.run_index import latest_table_path, load_latest_manifest, read_table


def _ranking_table(run_dir: str | Path, table_id: str) -> pd.DataFrame:
    man, _ = load_latest_manifest(run_dir)
    path = latest_table_path(man, table_id, preferred=("csv", "parquet"))
    df = read_table(path)
    if not df.empty:
        df["source_run_dir"] = str(run_dir)
    return df


def compare_runs(run_a: str | Path, run_b: str | Path, ranking: str = "candidate_summary") -> pd.DataFrame:
    """Compare candidate-level outputs between two runs.

    Designed for method/backend sensitivity checks, e.g. xTB vs ORCA or
    different DFT methods. Missing columns are tolerated to support partial runs.
    """
    df_a = _ranking_table(run_a, ranking)
    df_b = _ranking_table(run_b, ranking)
    if df_a.empty or df_b.empty:
        return pd.DataFrame(columns=["mol_id", "compare_status"])
    keys = ["mol_id"]
    for extra in ["site_id", "hf_n"]:
        if extra in df_a.columns and extra in df_b.columns:
            keys.append(extra)
    keep = list(dict.fromkeys(keys + [
        "scavenger_score", "activation_score", "quality_tier", "confidence_score",
        "delta_G_assoc_pressure_corrected_kcal_mol", "delta_G_act_kcal_mol", "k_corrected_s-1",
        "recommended_next_action",
    ]))
    a = df_a[[c for c in keep if c in df_a.columns]].copy()
    b = df_b[[c for c in keep if c in df_b.columns]].copy()
    out = a.merge(b, on=keys, how="outer", suffixes=("_a", "_b"), indicator=True)
    for col in ["scavenger_score", "activation_score", "confidence_score", "delta_G_assoc_pressure_corrected_kcal_mol", "delta_G_act_kcal_mol", "k_corrected_s-1"]:
        ca, cb = f"{col}_a", f"{col}_b"
        if ca in out.columns and cb in out.columns:
            out[f"delta_{col}"] = pd.to_numeric(out[cb], errors="coerce") - pd.to_numeric(out[ca], errors="coerce")
    out["compare_status"] = out["_merge"].map({"both": "matched", "left_only": "only_a", "right_only": "only_b"})
    return out.drop(columns=["_merge"])


def compare_backend_calculations(run_a: str | Path, run_b: str | Path) -> pd.DataFrame:
    man_a, _ = load_latest_manifest(run_a)
    man_b, _ = load_latest_manifest(run_b)
    rows_a = []
    rows_b = []
    for label, man, rows in [("a", man_a, rows_a), ("b", man_b, rows_b)]:
        for art in man.artifacts:
            if art.artifact_type != "calculation":
                continue
            rows.append({
                "species_id": art.data.get("species_id") or art.data.get("source_species_id"),
                "reaction_id": art.data.get("reaction_id"),
                "stage": (art.method or {}).get("stage"),
                "engine": (art.method or {}).get("engine"),
                "method_id": (art.method or {}).get("method_id"),
                "electronic_energy_hartree": art.data.get("electronic_energy_hartree"),
                "gibbs_298K_hartree": art.data.get("gibbs_298K_hartree"),
                "n_imag": art.data.get("n_imag"),
            })
    a = pd.DataFrame(rows_a)
    b = pd.DataFrame(rows_b)
    if a.empty or b.empty:
        return pd.DataFrame(columns=["species_id", "stage", "compare_status"])
    keys = [c for c in ["species_id", "reaction_id", "stage"] if c in a.columns and c in b.columns]
    out = a.merge(b, on=keys, how="outer", suffixes=("_a", "_b"), indicator=True)
    if "electronic_energy_hartree_a" in out.columns and "electronic_energy_hartree_b" in out.columns:
        out["delta_E_kcal_mol"] = (pd.to_numeric(out["electronic_energy_hartree_b"], errors="coerce") - pd.to_numeric(out["electronic_energy_hartree_a"], errors="coerce")) * 627.509474
    if "gibbs_298K_hartree_a" in out.columns and "gibbs_298K_hartree_b" in out.columns:
        out["delta_G_298K_kcal_mol"] = (pd.to_numeric(out["gibbs_298K_hartree_b"], errors="coerce") - pd.to_numeric(out["gibbs_298K_hartree_a"], errors="coerce")) * 627.509474
    out["compare_status"] = out["_merge"].map({"both": "matched", "left_only": "only_a", "right_only": "only_b"})
    return out.drop(columns=["_merge"])
