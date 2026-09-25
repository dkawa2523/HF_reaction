from __future__ import annotations

"""Production connector audit stage.

This stage is intentionally read-only: it summarizes whether external connectors
(GoodVibes, Arkane, Cantera, public DB providers) were actually used and whether
outputs are production-ready or still screening/fallback data.
"""

from collections import Counter, defaultdict
from html import escape
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class ConnectorAuditStage(Stage):
    name = "connector-audit"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)

        rows: list[dict[str, Any]] = []
        counters: Counter[str] = Counter()

        for st in manifest.iter_artifacts("species_thermo"):
            data = st.data
            rows.append({
                "connector": "goodvibes",
                "artifact_id": st.artifact_id,
                "species_id": data.get("species_id"),
                "status": data.get("external_goodvibes_status") or "not_requested",
                "executed": bool(data.get("external_goodvibes_executed")),
                "production_ready": bool(data.get("production_thermo_ready")),
                "quality": data.get("connector_quality"),
                "path": data.get("goodvibes_csv_path"),
                "note": "species thermochemistry",
            })
            counters[f"goodvibes_status_{data.get('external_goodvibes_status') or 'not_requested'}"] += 1

        for art in manifest.iter_artifacts("connector_validation"):
            data = art.data
            rows.append({
                "connector": data.get("connector") or art.qc.get("connector") or "cantera",
                "artifact_id": art.artifact_id,
                "species_id": None,
                "status": data.get("cantera_validation_status") or data.get("status"),
                "executed": bool(data.get("cantera_available") and data.get("cantera_validation_status") not in {"not_requested", None}),
                "production_ready": bool(data.get("production_cantera_ready")),
                "quality": "production" if data.get("production_cantera_ready") else "draft_or_smoke_test",
                "path": (art.paths or {}).get("json"),
                "note": data.get("reason") or "Cantera validation",
            })

        for art in manifest.iter_artifacts("arkane_input"):
            data = art.data
            rows.append({
                "connector": "arkane",
                "artifact_id": art.artifact_id,
                "species_id": None,
                "status": data.get("arkane_status") or "skeleton",
                "executed": bool(data.get("arkane_external_executed")),
                "production_ready": bool(data.get("production_arkane_ready")),
                "quality": "production" if data.get("production_arkane_ready") else "skeleton_or_not_run",
                "path": (art.paths or {}).get("python"),
                "note": data.get("readiness_note") or "Arkane connector",
            })

        # Public DB provider status from enriched molecule artifacts.
        provider_counts: dict[str, Counter[str]] = defaultdict(Counter)
        for mol in list(manifest.iter_artifacts("molecule_enriched")) or list(manifest.iter_artifacts("molecule")):
            for name, status in (mol.data.get("db_provider_status") or {}).items():
                s = str((status or {}).get("status") or "unknown")
                matched = bool((status or {}).get("matched"))
                rows.append({
                    "connector": f"db:{name}",
                    "artifact_id": mol.artifact_id,
                    "species_id": None,
                    "status": s,
                    "executed": s not in {"offline_fixture", "skipped_missing_api_key", "network_disabled"},
                    "production_ready": matched and s in {"success", "offline_fixture"},
                    "quality": "matched" if matched else "not_matched",
                    "path": None,
                    "note": (status or {}).get("reason") or (status or {}).get("source"),
                })
                provider_counts[str(name)][s] += 1

        df = pd.DataFrame(rows, columns=["connector", "artifact_id", "species_id", "status", "executed", "production_ready", "quality", "path", "note"])
        csv_path = out_dir / "production_connector_audit.csv"
        df.to_csv(csv_path, index=False)
        write_jsonl(rows, out_dir / "production_connector_audit.jsonl")

        summary_rows: list[dict[str, Any]] = []
        if not df.empty:
            for connector, sub in df.groupby("connector"):
                summary_rows.append({
                    "connector": connector,
                    "n_records": len(sub),
                    "n_executed": int(sub["executed"].fillna(False).astype(bool).sum()),
                    "n_production_ready": int(sub["production_ready"].fillna(False).astype(bool).sum()),
                    "statuses": ";".join(f"{k}:{v}" for k, v in Counter(sub["status"].fillna("unknown")).items()),
                })
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = out_dir / "production_connector_summary.csv"
        summary_df.to_csv(summary_csv, index=False)
        summary = {
            "n_rows": len(df),
            "n_connectors": int(summary_df["connector"].nunique()) if not summary_df.empty else 0,
            "n_production_ready": int(df["production_ready"].fillna(False).astype(bool).sum()) if not df.empty else 0,
            "all_production_ready": bool((not df.empty) and df["production_ready"].fillna(False).astype(bool).all()),
            "scientific_note": "Production-ready connector status only means the external tool/provider ran or matched successfully; chemical validity still requires method validation and expert review.",
        }
        write_json(out_dir / "production_connector_summary.json", summary)
        report = self._write_report(out_dir / "production_connector_audit.html", summary, df, summary_df)

        out.add_artifact(Artifact(
            artifact_id="production_connector_audit",
            artifact_type="connector_audit",
            paths={"csv": str(csv_path), "jsonl": str(out_dir / "production_connector_audit.jsonl"), "html": str(report), "summary_csv": str(summary_csv)},
            data=summary,
            qc={"connector_audit_status": "success", "all_production_ready": summary["all_production_ready"]},
        ))
        out.add_artifact(Artifact(
            artifact_id="production_connector_audit_table",
            artifact_type="table",
            paths={"csv": str(csv_path), "summary_csv": str(summary_csv)},
            data={"n_rows": len(df), "table_type": "production_connector_audit"},
        ))
        return out

    @staticmethod
    def _write_report(path: Path, summary: dict[str, Any], df: pd.DataFrame, summary_df: pd.DataFrame) -> Path:
        def table(x: pd.DataFrame, n: int = 60) -> str:
            if x.empty:
                return "<p>No rows.</p>"
            return x.head(n).to_html(index=False, escape=True)
        html = [
            "<html><head><meta charset='utf-8'><title>hfauto production connector audit</title>",
            "<style>body{font-family:Arial,sans-serif;margin:2rem;} table{border-collapse:collapse;font-size:0.9rem;} th,td{border:1px solid #ccc;padding:4px 6px;} th{background:#eee;} .warn{background:#fff2cc;padding:0.6rem;}</style>",
            "</head><body>",
            "<h1>Production connector audit</h1>",
            f"<p class='warn'>{escape(str(summary.get('scientific_note')))}</p>",
            "<h2>Summary</h2><table>",
        ]
        for k, v in summary.items():
            html.append(f"<tr><th>{escape(str(k))}</th><td>{escape(str(v))}</td></tr>")
        html.extend(["</table>", "<h2>By connector</h2>", table(summary_df), "<h2>Records</h2>", table(df), "</body></html>"])
        path.write_text("\n".join(html), encoding="utf-8")
        return path
