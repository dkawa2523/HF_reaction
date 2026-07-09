from __future__ import annotations
from html import escape as _esc
from pathlib import Path
from typing import Any


def esc(v: Any) -> str:
    return _esc("" if v is None else str(v), quote=True)


def write_html(path: str | Path, title: str, body: str, extra_head: str = "") -> Path:
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(title)}</title>
<style>
body{{font-family:system-ui,-apple-system,Segoe UI,sans-serif;margin:0;background:#fafafa;color:#222}}header{{background:#17202a;color:white;padding:1rem 1.4rem}}main{{max-width:1500px;margin:0 auto;padding:1rem}}.card{{background:white;border:1px solid #ddd;border-radius:10px;padding:1rem;margin:1rem 0;box-shadow:0 1px 2px rgba(0,0,0,.04)}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:1rem}}.warn{{border-left:5px solid #d08b00;background:#fff8e6;padding:.7rem}}.data-table{{border-collapse:collapse;width:100%;font-size:.9rem}}.data-table th,.data-table td{{border:1px solid #e0e0e0;padding:.35rem .5rem;vertical-align:top}}.data-table th{{background:#f1f4f6}}.viewer{{width:100%;height:520px;border:1px solid #ddd;border-radius:8px;background:white}}.svg-box{{overflow:auto;background:white;border:1px solid #ddd;border-radius:8px;padding:.5rem}}.muted{{color:#666}}.badge{{display:inline-block;border:1px solid #bbb;border-radius:999px;padding:.15rem .5rem;margin:.1rem;background:#f7f7f7}}
</style>{extra_head}</head><body><header><h1>{esc(title)}</h1></header><main>{body}</main></body></html>"""
    p.write_text(html, encoding="utf-8"); return p


def dataframe_to_html(df, max_rows: int = 30) -> str:
    if df is None or getattr(df, "empty", True): return "<p>No rows.</p>"
    return df.head(max_rows).to_html(index=False, escape=True)


def table_html(rows: list[dict[str, Any]], max_rows: int = 30) -> str:
    import pandas as pd
    return dataframe_to_html(pd.DataFrame(rows), max_rows=max_rows)

# Phase 8 compatibility alias.
def df_to_html(df, max_rows: int = 30) -> str:
    return dataframe_to_html(df, max_rows=max_rows)

# Phase 8 compatibility alias.
def df_to_html(df, max_rows: int = 30) -> str:
    return dataframe_to_html(df, max_rows=max_rows)
