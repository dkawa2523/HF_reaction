from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _safe_nested(d: Any, *keys: str, default=None):
    cur = d
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _identity_fields(mol: dict[str, Any]) -> dict[str, Any]:
    ident = mol.get("identity") or {}
    public = mol.get("public_data") or {}
    gas = public.get("gas_process") or {}
    nist = public.get("nist_webbook") or {}
    comp = public.get("comptox") or {}
    niosh = public.get("niosh") or {}
    cas = public.get("cas_common_chemistry") or {}
    pubchem = public.get("pubchem") or {}
    support = 0.0
    notes: list[str] = []
    if ident.get("identity_confidence") == "high" or pubchem.get("matched"):
        support += 0.25
        notes.append("identity")
    if nist.get("proton_affinity_kj_mol") is not None or nist.get("gas_basicity_kj_mol") is not None:
        support += 0.30
        notes.append("NIST_PA_GB")
    if public.get("atct", {}).get("matched"):
        support += 0.10
        notes.append("ATcT")
    if public.get("cccbdb", {}).get("matched"):
        support += 0.10
        notes.append("CCCBDB")
    if comp.get("matched") or gas.get("gas_process_feasibility"):
        support += 0.15
        notes.append("process_properties")
    if niosh.get("matched"):
        support += 0.07
        notes.append("NIOSH")
    if cas.get("matched"):
        support += 0.03
        notes.append("CAS")
    support = min(1.0, support)
    if ident.get("identity_conflict"):
        support = max(0.0, support - 0.4)
        notes.append("identity_conflict")
    pa = nist.get("proton_affinity_kj_mol")
    gb = nist.get("gas_basicity_kj_mol")
    return {
        "mol_id": mol.get("mol_id"),
        "name": mol.get("name"),
        "canonical_smiles": mol.get("canonical_smiles"),
        "inchikey": mol.get("inchikey"),
        "formula": mol.get("formula"),
        "identity_confidence": ident.get("identity_confidence"),
        "identity_conflict": ident.get("identity_conflict"),
        "pubchem_cid": ident.get("pubchem_cid"),
        "cas_rn": ident.get("cas_rn"),
        "proton_affinity_kj_mol": pa,
        "gas_basicity_kj_mol": gb,
        "gas_process_feasibility": gas.get("gas_process_feasibility") or comp.get("gas_process_feasibility"),
        "ehs_review_flag": gas.get("ehs_review_flag") or comp.get("ehs_review_flag"),
        "vapor_pressure_Pa_25C": gas.get("vapor_pressure_Pa_25C") or comp.get("vapor_pressure_Pa_25C"),
        "boiling_point_C": gas.get("boiling_point_C") or comp.get("boiling_point_C"),
        "manual_ehs_review_required": bool(niosh.get("manual_ehs_review_required") or str(gas.get("ehs_review_flag") or comp.get("ehs_review_flag") or "").lower().find("manual") >= 0),
        "nist_matched": bool(nist.get("matched")),
        "atct_matched": bool(_safe_nested(public, "atct", "matched", default=False)),
        "cccbdb_matched": bool(_safe_nested(public, "cccbdb", "matched", default=False)),
        "comptox_matched": bool(comp.get("matched")),
        "niosh_matched": bool(niosh.get("matched")),
        "cas_matched": bool(cas.get("matched")),
        "public_db_support_score": round(support, 3),
        "reference_support_level": "basicity_reference" if pa is not None or gb is not None else ("process_identity_reference" if support > 0 else "none"),
        "calibration_action": "use_for_basicity_trend_validation" if pa is not None or gb is not None else "no_public_basicity_reference",
        "support_notes": ";".join(notes),
        "db_provider_status": mol.get("db_provider_status", {}),
    }


def _flatten_reference_values(mol: dict[str, Any]) -> list[dict[str, Any]]:
    public = mol.get("public_data") or {}
    rows: list[dict[str, Any]] = []
    nist = public.get("nist_webbook") or {}
    for prop, key, units, db in [
        ("proton_affinity_kj_mol", "proton_affinity_kj_mol", "kJ/mol", "NIST WebBook"),
        ("gas_basicity_kj_mol", "gas_basicity_kj_mol", "kJ/mol", "NIST WebBook"),
    ]:
        if nist.get(key) is not None:
            rows.append({"mol_id": mol.get("mol_id"), "name": mol.get("name"), "inchikey": mol.get("inchikey"), "property": prop, "reference_value": nist.get(key), "units": units, "source": nist.get("source") or db, "reference_status": "available"})
    atct = public.get("atct") or {}
    if atct.get("delta_f_H_298_kj_mol") is not None:
        rows.append({"mol_id": mol.get("mol_id"), "name": mol.get("name"), "inchikey": mol.get("inchikey"), "property": "delta_f_H_298_kj_mol", "reference_value": atct.get("delta_f_H_298_kj_mol"), "units": "kJ/mol", "source": atct.get("source") or "ATcT", "reference_status": "available"})
    cccbdb = public.get("cccbdb") or {}
    if cccbdb.get("dipole_D") is not None:
        rows.append({"mol_id": mol.get("mol_id"), "name": mol.get("name"), "inchikey": mol.get("inchikey"), "property": "dipole_D", "reference_value": cccbdb.get("dipole_D"), "units": "D", "source": cccbdb.get("source") or "CCCBDB", "reference_status": "available"})
    freq = cccbdb.get("hf_stretch_cm1") or (cccbdb.get("vibrational_frequencies_cm1") or [None])[0]
    if freq is not None:
        rows.append({"mol_id": mol.get("mol_id"), "name": mol.get("name"), "inchikey": mol.get("inchikey"), "property": "hf_stretch_cm1", "reference_value": freq, "units": "cm^-1", "source": cccbdb.get("source") or "CCCBDB", "reference_status": "available"})
    return rows


def _to_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


class CalibrateStage(Stage):
    """Audit public DB enrichment and build method/reference validation outputs."""

    name = "calibrate"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)

        mols = [a.data for a in list(manifest.iter_artifacts("molecule_enriched")) or list(manifest.iter_artifacts("molecule"))]
        coverage_rows = [_identity_fields(m) for m in mols]
        coverage_df = pd.DataFrame(coverage_rows)
        coverage_path = out_dir / "public_data_coverage.csv"
        coverage_df.to_csv(coverage_path, index=False)

        audit_path = out_dir / "public_db_audit.csv"
        coverage_df.to_csv(audit_path, index=False)

        ref_rows: list[dict[str, Any]] = []
        for m in mols:
            ref_rows.extend(_flatten_reference_values(m))
        refs_df = pd.DataFrame(ref_rows, columns=["mol_id", "name", "inchikey", "property", "reference_value", "units", "source", "reference_status"])
        refs_path = out_dir / "calibration_reference_records.csv"
        refs_jsonl = out_dir / "calibration_reference_records.jsonl"
        refs_df.to_csv(refs_path, index=False)
        write_jsonl(ref_rows, refs_jsonl)

        thermo_df = pd.DataFrame([a.data for a in manifest.iter_artifacts("thermo")])
        best_thermo = pd.DataFrame()
        if not thermo_df.empty and "mol_id" in thermo_df.columns:
            sort_cols = [c for c in ["T_K", "delta_G_assoc_pressure_corrected_kcal_mol"] if c in thermo_df.columns]
            best_thermo = thermo_df.sort_values(sort_cols).groupby("mol_id", as_index=False).first()
        metrics: list[dict[str, Any]] = []
        if not coverage_df.empty and not best_thermo.empty and "proton_affinity_kj_mol" in coverage_df.columns:
            merged = coverage_df.merge(best_thermo[[c for c in ["mol_id", "delta_G_assoc_pressure_corrected_kcal_mol", "delta_G_act_kcal_mol"] if c in best_thermo.columns]], on="mol_id", how="inner")
            pa = pd.to_numeric(merged.get("proton_affinity_kj_mol"), errors="coerce")
            assoc = -pd.to_numeric(merged.get("delta_G_assoc_pressure_corrected_kcal_mol"), errors="coerce")
            valid = pa.notna() & assoc.notna()
            if valid.sum() >= 2:
                corr = float(pa[valid].rank().corr(assoc[valid].rank(), method="pearson"))
            else:
                corr = None
            metrics.append({"metric": "basicity_rank_spearman_proxy", "value": corr, "n": int(valid.sum()), "note": "Spearman-like rank correlation between NIST PA and stronger HF association proxy; screening diagnostic only"})
        else:
            metrics.append({"metric": "basicity_rank_spearman_proxy", "value": None, "n": 0, "note": "No matched PA/thermo records"})
        metrics.extend([
            {"metric": "public_db_mean_support", "value": float(coverage_df["public_db_support_score"].mean()) if not coverage_df.empty else 0.0, "n": len(coverage_df), "note": "Average candidate public DB support score"},
            {"metric": "identity_conflict_count", "value": int(coverage_df["identity_conflict"].fillna(False).astype(bool).sum()) if not coverage_df.empty and "identity_conflict" in coverage_df else 0, "n": len(coverage_df), "note": "Identity conflicts should be reviewed before production ranking"},
            {"metric": "manual_ehs_review_count", "value": int(coverage_df["manual_ehs_review_required"].fillna(False).astype(bool).sum()) if not coverage_df.empty and "manual_ehs_review_required" in coverage_df else 0, "n": len(coverage_df), "note": "Candidates requiring process/EHS review"},
        ])
        metrics_df = pd.DataFrame(metrics)
        metrics_path = out_dir / "method_validation_metrics.csv"
        metrics_df.to_csv(metrics_path, index=False)

        residual_rows: list[dict[str, Any]] = []
        if not coverage_df.empty and not best_thermo.empty:
            cols = [c for c in ["mol_id", "delta_G_assoc_pressure_corrected_kcal_mol", "delta_G_act_kcal_mol", "quality_tier", "confidence_score"] if c in best_thermo.columns]
            merged_res = coverage_df.merge(best_thermo[cols], on="mol_id", how="left") if cols else coverage_df.copy()
            for _, row in merged_res.iterrows():
                pa = _to_float(row.get("proton_affinity_kj_mol"))
                assoc = _to_float(row.get("delta_G_assoc_pressure_corrected_kcal_mol"))
                residual_rows.append({
                    "mol_id": row.get("mol_id"),
                    "name": row.get("name"),
                    "reference_property": "proton_affinity_kj_mol",
                    "reference_value": pa,
                    "computed_proxy": "-delta_G_assoc_pressure_corrected_kcal_mol",
                    "computed_proxy_value": None if assoc is None else -assoc,
                    "proxy_residual_note": "PA/GB are trend-validation proxies, not direct HF reaction energies",
                    "quality_tier": row.get("quality_tier"),
                    "confidence_score": row.get("confidence_score"),
                })
        residuals_df = pd.DataFrame(residual_rows)
        residuals_path = out_dir / "reference_residuals.csv"
        residuals_df.to_csv(residuals_path, index=False)

        support_by_mol = {str(r["mol_id"]): float(r.get("public_db_support_score") or 0.0) for r in coverage_rows}
        support_notes = {str(r["mol_id"]): str(r.get("support_notes") or "").split(";") if r.get("support_notes") else [] for r in coverage_rows}
        calibration_summary = {
            "status": "reference_values_available" if len(ref_rows) else "no_public_reference_values",
            "n_molecules": len(mols),
            "n_reference_values": len(ref_rows),
            "calibration_support_score": float(coverage_df["public_db_support_score"].mean()) if not coverage_df.empty else 0.0,
            "db_calibration_support_by_mol": support_by_mol,
            "support_notes_by_mol": support_notes,
            "scientific_note": "PA/GB validate gas-phase basicity trends; they are not direct substitutes for HF association, ion-pair, or TS free energies.",
        }
        write_jsonl([calibration_summary], out_dir / "calibration_summary.jsonl")
        report_path = self._write_report(out_dir / "method_validation_report.html", calibration_summary, coverage_df, refs_df, metrics_df)
        dashboard_path = self._write_dashboard(out_dir / "method_calibration_dashboard.html", calibration_summary, coverage_df, refs_df, metrics_df, residuals_df)

        out.add_artifact(Artifact(artifact_id="public_data_coverage", artifact_type="table", paths={"csv": str(coverage_path)}, data={"n_rows": len(coverage_df), "table_type": "public_data_coverage"}))
        out.add_artifact(Artifact(artifact_id="public_db_audit", artifact_type="table", paths={"csv": str(audit_path)}, data={"n_rows": len(coverage_df), "table_type": "public_db_audit"}))
        out.add_artifact(Artifact(artifact_id="calibration_reference_records", artifact_type="table", paths={"csv": str(refs_path), "jsonl": str(refs_jsonl)}, data={"n_rows": len(refs_df), "table_type": "calibration_reference_records"}))
        out.add_artifact(Artifact(artifact_id="method_validation_metrics", artifact_type="table", paths={"csv": str(metrics_path)}, data={"n_rows": len(metrics_df), "table_type": "method_validation_metrics"}))
        out.add_artifact(Artifact(artifact_id="reference_residuals", artifact_type="table", paths={"csv": str(residuals_path)}, data={"n_rows": len(residuals_df), "table_type": "reference_residuals"}))
        out.add_artifact(Artifact(artifact_id="method_validation_report", artifact_type="report", paths={"html": str(report_path), "dashboard_html": str(dashboard_path)}, data=calibration_summary))
        out.add_artifact(Artifact(artifact_id="calibration_summary", artifact_type="calibration", paths={"jsonl": str(out_dir / "calibration_summary.jsonl")}, data=calibration_summary, qc={"calibration_support_score": calibration_summary["calibration_support_score"]}))
        out.add_artifact(Artifact(artifact_id="method_validation", artifact_type="method_validation", paths={"html": str(report_path), "dashboard_html": str(dashboard_path), "metrics_csv": str(metrics_path), "coverage_csv": str(coverage_path), "residuals_csv": str(residuals_path)}, data=calibration_summary, qc={"calibration_support_score": calibration_summary["calibration_support_score"]}))
        return out

    @staticmethod
    def _write_report(path: Path, summary: dict[str, Any], coverage: pd.DataFrame, refs: pd.DataFrame, metrics: pd.DataFrame) -> Path:
        def table(df: pd.DataFrame, n: int = 30) -> str:
            if df.empty:
                return "<p>No rows.</p>"
            return df.head(n).to_html(index=False, escape=True)
        html = [
            "<html><head><meta charset='utf-8'><title>hfauto method validation</title>",
            "<style>body{font-family:Arial,sans-serif;margin:2rem;} table{border-collapse:collapse;font-size:0.9rem;} th,td{border:1px solid #ccc;padding:4px 6px;} th{background:#eee;}</style>",
            "</head><body>",
            "<h1>Method validation and public DB audit</h1>",
            f"<p>Status: <b>{escape(str(summary.get('status')))}</b></p>",
            "<h2>Summary</h2><table>",
        ]
        for k, v in summary.items():
            if k.endswith("_by_mol"):
                continue
            html.append(f"<tr><th>{escape(str(k))}</th><td>{escape(str(v))}</td></tr>")
        html.extend(["</table>", "<h2>Method validation metrics</h2>", table(metrics), "<h2>Public data coverage</h2>", table(coverage), "<h2>Reference records</h2>", table(refs), "</body></html>"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(html), encoding="utf-8")
        return path

    @staticmethod
    def _write_dashboard(path: Path, summary: dict[str, Any], coverage: pd.DataFrame, refs: pd.DataFrame, metrics: pd.DataFrame, residuals: pd.DataFrame) -> Path:
        def table(df: pd.DataFrame, n: int = 30) -> str:
            if df.empty:
                return "<p>No rows.</p>"
            return df.head(n).to_html(index=False, escape=True)
        cards = [
            ("Molecules", summary.get("n_molecules")),
            ("Reference values", summary.get("n_reference_values")),
            ("Mean DB support", round(float(summary.get("calibration_support_score") or 0.0), 3)),
            ("Status", summary.get("status")),
        ]
        html = [
            "<html><head><meta charset='utf-8'><title>hfauto calibration dashboard</title>",
            "<style>body{font-family:Arial,sans-serif;margin:2rem;} .grid{display:flex;gap:1rem;flex-wrap:wrap}.card{border:1px solid #ccc;border-radius:8px;padding:1rem;min-width:160px;background:#f8f8f8} table{border-collapse:collapse;font-size:0.9rem;} th,td{border:1px solid #ccc;padding:4px 6px;} th{background:#eee;} .note{background:#eef6ff;padding:0.6rem}</style>",
            "</head><body>",
            "<h1>Method calibration dashboard</h1>",
            f"<p class='note'>{escape(str(summary.get('scientific_note')))}</p>",
            "<div class='grid'>",
        ]
        for k, v in cards:
            html.append(f"<div class='card'><b>{escape(str(k))}</b><br>{escape(str(v))}</div>")
        html.extend(["</div>", "<h2>Validation metrics</h2>", table(metrics), "<h2>Reference residual/proxy table</h2>", table(residuals), "<h2>Public DB coverage</h2>", table(coverage), "<h2>Reference values</h2>", table(refs), "</body></html>"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(html), encoding="utf-8")
        return path
