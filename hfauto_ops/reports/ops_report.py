from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.io import ensure_dir


def _table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df is None or df.empty:
        return "<p><em>No records.</em></p>"
    return df.head(max_rows).to_html(index=False, escape=False)


def render_ops_report(
    out_path: str | Path,
    summary: dict[str, Any],
    job_plan: pd.DataFrame,
    retry_plan: pd.DataFrame,
    duplicate_report: pd.DataFrame,
    comparison: pd.DataFrame | None = None,
) -> Path:
    p = Path(out_path)
    ensure_dir(p.parent)
    counts = summary.get("artifact_counts", {}) or {}
    quality = summary.get("quality_counts", {}) or {}
    html = f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>HF Auto Ops Report</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 24px; line-height: 1.45; }}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.9rem; }}
th,td {{ border: 1px solid #ddd; padding: 0.35rem; vertical-align: top; }}
th {{ background: #f2f4f7; }}
.card {{ border: 1px solid #ddd; border-radius: 8px; padding: 1rem; margin: 1rem 0; }}
.warn {{ background: #fff7e6; border-left: 4px solid #faad14; padding: 0.8rem; }}
code {{ background: #f6f8fa; padding: 0.1rem 0.25rem; }}
</style></head><body>
<h1>HF Auto Phase 9 Operations Report</h1>
<div class=\"card\">
<h2>Run summary</h2>
<p><b>run_id:</b> {summary.get('run_id')}<br>
<b>latest stage:</b> {summary.get('stage')}<br>
<b>artifacts:</b> {summary.get('n_artifacts')}<br>
<b>failures:</b> {summary.get('n_failures')}</p>
<h3>Artifact counts</h3><pre>{counts}</pre>
<h3>Quality counts</h3><pre>{quality}</pre>
</div>
<div class=\"card\"><h2>HPC job plan</h2>
<p>Rows: {len(job_plan)}. This is a scheduler-neutral plan; review commands/resources before submission.</p>
{_table_html(job_plan)}
</div>
<div class=\"card\"><h2>Retry plan</h2>
<p>Rows: {len(retry_plan)}. Failures are prioritized by category and recoverability.</p>
{_table_html(retry_plan)}
</div>
<div class=\"card\"><h2>Duplicate / cache report</h2>
<p>Rows: {len(duplicate_report)}. Fingerprints identify repeated calculation-like artifacts.</p>
{_table_html(duplicate_report)}
</div>
"""
    if comparison is not None:
        html += f"<div class=\"card\"><h2>Run/backend comparison</h2>{_table_html(comparison)}</div>"
    html += """
<div class=\"warn\"><b>Operational note:</b> generated scheduler scripts are scaffolds. Validate account, partition, modules, environment activation, licenses, scratch policy, and safety constraints before production submission.</div>
</body></html>
"""
    p.write_text(html, encoding="utf-8")
    return p
