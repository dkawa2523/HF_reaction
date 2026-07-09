from __future__ import annotations

from pathlib import Path

import pandas as pd

from hfauto.core.io import ensure_dir


def render_hpc_dashboard(out_dir: str | Path, job_plan_csv: str | Path | None = None, retry_plan_csv: str | Path | None = None, cache_index_csv: str | Path | None = None) -> Path:
    out = ensure_dir(out_dir)
    sections: list[str] = ["<h1>HF Auto HPC Operations Dashboard</h1>"]
    for title, path in [("Job plan", job_plan_csv), ("Retry plan", retry_plan_csv), ("Cache index", cache_index_csv)]:
        sections.append(f"<h2>{title}</h2>")
        if path and Path(path).exists():
            try:
                df = pd.read_csv(path)
            except pd.errors.EmptyDataError:
                df = pd.DataFrame()
            sections.append(f"<p>Rows: {len(df)}</p>")
            if not df.empty:
                sections.append(df.head(50).to_html(index=False, escape=False))
        else:
            sections.append("<p>Not available.</p>")
    html = "<!doctype html><html><head><meta charset='utf-8'><title>HPC Dashboard</title><style>body{font-family:system-ui;margin:2rem}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:4px 8px}th{background:#f2f2f2}</style></head><body>" + "\n".join(sections) + "</body></html>"
    path = out / "hpc_dashboard.html"
    path.write_text(html, encoding="utf-8")
    return path
