from __future__ import annotations

from pathlib import Path

import pandas as pd


def _table_html(path: Path, max_rows: int = 30) -> str:
    if not path.exists():
        return "<p>not generated</p>"
    try:
        df = pd.read_csv(path)
    except Exception as exc:
        return f"<p>could not read {path.name}: {exc}</p>"
    if df.empty:
        return "<p>No rows.</p>"
    return df.head(max_rows).to_html(index=False, escape=False)


def write_operations_report(out_dir: str | Path, title: str = "hfauto operations report") -> Path:
    out = Path(out_dir)
    html = out / "operations_report.html"
    sections = [
        ("Stage index", out / "stage_index.csv"),
        ("Resource plan", out / "resource_plan.csv"),
        ("Artifact-level array plan", out / "array_job_plan.csv"),
        ("Scheduler status dry run", out / "scheduler_status.csv"),
        ("Retry plan", out / "retry_plan.csv"),
        ("Reuse plan", out / "reuse_plan.csv"),
        ("Registry duplicates", out / "registry_duplicate_report.csv"),
        ("Backend comparison", out / "backend_comparison" / "backend_comparison.csv"),
        ("Resource autotune", out / "resource_autotune_recommendations.csv"),
        ("Duplicate/reuse candidates", out / "duplicate_candidates.csv"),
        ("Artifact index", out / "artifact_index.csv"),
    ]
    body = [
        "<html><head><meta charset='utf-8'><title>hfauto ops</title>",
        "<style>body{font-family:Arial,sans-serif;margin:24px} table{border-collapse:collapse;font-size:12px} td,th{border:1px solid #ddd;padding:4px} th{background:#f0f0f0}.warn{background:#fff3cd;padding:8px}.ok{background:#e9f7ef;padding:8px}</style>",
        "</head><body>",
    ]
    body.append(f"<h1>{title}</h1>")
    body.append("<p class='warn'>Operations planning artifacts are read-only. They schedule, retry, deduplicate, and compare workflow runs but do not modify scientific results.</p>")
    body.append("<p class='ok'>The operations layer provides scheduler dry runs, artifact-level job arrays, calculation reuse planning, QCArchive payload export, backend comparison, and resource auto-tuning proposals.</p>")
    for label, path in sections:
        body.append(f"<h2>{label}</h2>")
        body.append(f"<p><a href='{path.relative_to(out) if path.exists() else path.name}'>{path.name}</a></p>")
        body.append(_table_html(path))
    qcarchive_summary = out / "qcarchive" / "qcarchive_summary.json"
    if qcarchive_summary.exists():
        body.append("<h2>QCArchive payload</h2><pre>" + qcarchive_summary.read_text(encoding="utf-8") + "</pre>")
    body.append("</body></html>")
    html.write_text("\n".join(body), encoding="utf-8")
    return html
