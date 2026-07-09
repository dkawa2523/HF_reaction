from __future__ import annotations

from collections import Counter
from html import escape
from pathlib import Path

import pandas as pd

from hfauto.core.io import read_manifest


def latest_manifest_path(run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    pointer = run_dir / "manifest.path"
    if pointer.exists():
        return Path(pointer.read_text(encoding="utf-8").strip())
    manifests = sorted(run_dir.glob("*/manifest.json"))
    if not manifests:
        raise FileNotFoundError(f"No manifest found in {run_dir}")
    return manifests[-1]


def _html_table_from_csv(path: str | Path, max_rows: int = 20) -> str:
    p = Path(path)
    if not p.exists():
        return f"<p>Missing table: {escape(str(p))}</p>"
    try:
        df = pd.read_csv(p).head(max_rows)
        return df.to_html(index=False, escape=True)
    except Exception as exc:
        return f"<p>Could not read {escape(str(p))}: {escape(str(exc))}</p>"


def render_report(run_dir: str | Path, out_path: str | Path) -> Path:
    manifest = read_manifest(latest_manifest_path(run_dir))
    counts = Counter(a.artifact_type for a in manifest.artifacts)
    failures = [a for a in manifest.artifacts if a.status.status == "failed"]
    failure_counts = Counter(f.status.category or "unknown" for f in failures)
    rankings = [a for a in manifest.artifacts if a.artifact_type == "ranking"]
    tables = {a.artifact_id: a for a in manifest.artifacts if a.artifact_type == "table"}
    quality_counts = Counter()
    for t in manifest.iter_artifacts("thermo"):
        quality_counts[t.data.get("quality_tier", "unknown")] += 1
    kinetics_count = sum(1 for _ in manifest.iter_artifacts("kinetics"))
    mechanism_count = sum(1 for _ in manifest.iter_artifacts("mechanism"))
    method_validations = [a for a in manifest.iter_artifacts("method_validation")] + [a for a in manifest.iter_artifacts("calibration")]
    method_summary = tables.get("method_validation_summary") or tables.get("method_validation_metrics")
    db_summary = tables.get("public_db_enrichment_summary") or tables.get("db_enrichment_summary") or tables.get("public_data_coverage") or tables.get("public_db_coverage_calibration")

    html = [
        "<html><head><meta charset='utf-8'><title>hfauto report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:2rem;} table{border-collapse:collapse;font-size:0.9rem;} th,td{border:1px solid #ccc;padding:4px 6px;} th{background:#eee;} code{background:#f5f5f5;padding:2px 4px;}</style>",
        "</head><body>",
    ]
    html.append(f"<h1>hfauto run report: {escape(manifest.run_id)}</h1>")
    html.append(f"<p>Latest stage: <b>{escape(manifest.stage)}</b></p>")
    html.append("<h2>Artifact counts</h2><table><tr><th>type</th><th>count</th></tr>")
    for key, value in sorted(counts.items()):
        html.append(f"<tr><td>{escape(key)}</td><td>{value}</td></tr>")
    html.append("</table>")

    html.append("<h2>Quality tiers</h2><table><tr><th>tier</th><th>count</th></tr>")
    for key, value in sorted(quality_counts.items()):
        html.append(f"<tr><td>{escape(str(key))}</td><td>{value}</td></tr>")
    html.append("</table>")
    html.append("<h2>Thermo / kinetics outputs</h2>")
    html.append(f"<p>Kinetics records: <b>{kinetics_count}</b>; mechanism artifacts: <b>{mechanism_count}</b></p>")

    html.append("<h2>Public DB / method validation</h2>")
    html.append(f"<p>Method validation records: <b>{len(method_validations)}</b></p>")
    if method_summary is not None:
        html.append("<table><tr><th>metric</th><th>value</th></tr>")
        for key, value in method_summary.data.items():
            html.append(f"<tr><td>{escape(str(key))}</td><td>{escape(str(value))}</td></tr>")
        html.append("</table>")
        preview_path = method_summary.paths.get("csv") or method_summary.paths.get("records_csv") or method_summary.paths.get("summary_csv")
        if preview_path:
            html.append("<h3>Method validation records preview</h3>")
            html.append(_html_table_from_csv(preview_path, max_rows=30))
    if db_summary is not None:
        html.append("<h3>DB enrichment summary</h3>")
        html.append(_html_table_from_csv(db_summary.paths.get("csv", ""), max_rows=20))

    html.append(f"<h2>Failures</h2><p>{len(failures)} failure artifacts</p>")
    if failure_counts:
        html.append("<table><tr><th>category</th><th>count</th></tr>")
        for key, value in sorted(failure_counts.items()):
            html.append(f"<tr><td>{escape(str(key))}</td><td>{value}</td></tr>")
        html.append("</table>")
    if failures:
        html.append("<ul>")
        for f in failures[:50]:
            html.append(f"<li><code>{escape(f.artifact_id)}</code>: {escape(str(f.status.category))} / {escape(str(f.status.reason))}</li>")
        html.append("</ul>")

    html.append("<h2>Output tables</h2><ul>")
    for r in rankings:
        path = r.paths.get("csv", "")
        html.append(f"<li>{escape(r.artifact_id)}: <code>{escape(path)}</code></li>")
    for table_id in ["candidate_summary", "reaction_results", "kinetics_table", "arrhenius_fits", "reactor_results", "species_thermo_table", "method_validation_summary", "method_validation_metrics", "public_data_coverage", "public_db_enrichment_summary", "identity_conflicts", "reference_matches", "calibration_comparisons", "method_validation_metrics", "calibration_gaps", "failure_report"]:
        if table_id in tables:
            html.append(f"<li>{escape(table_id)}: <code>{escape(str(tables[table_id].paths))}</code></li>")
    html.append("</ul>")

    if "candidate_summary" in tables:
        html.append("<h2>Candidate summary</h2>")
        html.append(_html_table_from_csv(tables["candidate_summary"].paths.get("csv", ""), max_rows=30))
    if "public_db_enrichment_summary" in tables:
        html.append("<h2>Public DB enrichment summary</h2>")
        html.append(_html_table_from_csv(tables["public_db_enrichment_summary"].paths.get("csv", ""), max_rows=50))
    if "reference_matches" in tables:
        html.append("<h2>Reference matches</h2>")
        html.append(_html_table_from_csv(tables["reference_matches"].paths.get("csv", ""), max_rows=50))
    if "method_validation_metrics" in tables:
        html.append("<h2>Method validation metrics</h2>")
        html.append(_html_table_from_csv(tables["method_validation_metrics"].paths.get("csv", ""), max_rows=50))
    for rep in [a for a in manifest.artifacts if a.artifact_type == "report" and a.artifact_id == "method_validation_report"]:
        html.append("<h2>Method validation report</h2>")
        html.append(f"<p><code>{escape(rep.paths.get('html', ''))}</code></p>")
    for r in rankings:
        html.append(f"<h2>{escape(r.artifact_id)}</h2>")
        html.append(_html_table_from_csv(r.paths.get("csv", ""), max_rows=20))

    html.append("</body></html>")
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(html), encoding="utf-8")
    return out
