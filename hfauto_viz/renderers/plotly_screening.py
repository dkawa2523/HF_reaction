from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd

from hfauto_viz.core.html import esc, write_html


def write_scatter(df: pd.DataFrame, out_path: str | Path, title: str, x: str, y: str, color: str | None=None, size: str | None=None) -> Path:
    if df is None or df.empty or x not in df.columns or y not in df.columns:
        return write_html(out_path,title,f'<div class="warn">Missing columns for {esc(title)}: {esc(x)}, {esc(y)}</div>')
    try:
        import plotly.express as px
        fig=px.scatter(df,x=x,y=y,color=color if color in df.columns else None,size=size if size in df.columns else None,hover_data=[c for c in df.columns[:20]],title=title,template='plotly_white')
        Path(out_path).parent.mkdir(parents=True,exist_ok=True); fig.write_html(out_path,include_plotlyjs='cdn',full_html=True); return Path(out_path)
    except Exception as exc:
        return write_html(out_path,title,f'<p>Plotly unavailable: {esc(exc)}</p>{df.head(30).to_html(index=False, escape=True)}')

def write_quality_funnel(counts: dict[str,int], out_path: str | Path) -> Path:
    try:
        import plotly.graph_objects as go
        fig=go.Figure(go.Funnel(y=list(counts.keys()), x=list(counts.values()), textinfo='value+percent initial')); fig.update_layout(title='Calculation campaign quality funnel',template='plotly_white')
        Path(out_path).parent.mkdir(parents=True,exist_ok=True); fig.write_html(out_path,include_plotlyjs='cdn',full_html=True); return Path(out_path)
    except Exception:
        return write_html(out_path,'Quality funnel','<table class="data-table">'+''.join(f'<tr><td>{esc(k)}</td><td>{v}</td></tr>' for k,v in counts.items())+'</table>')

def write_failure_heatmap(failures: pd.DataFrame | Counter, out_path: str | Path) -> Path:
    df = pd.DataFrame([{'category':k,'count':v} for k,v in failures.items()]) if isinstance(failures,Counter) else failures
    if df is None or df.empty: return write_html(out_path,'Failure heatmap','<p>No failure artifacts.</p>')
    try:
        import plotly.express as px
        if 'stage' in df.columns and 'category' in df.columns:
            pivot=df.pivot_table(index='stage',columns='category',values='artifact_id',aggfunc='count',fill_value=0); fig=px.imshow(pivot,text_auto=True,title='Failure heatmap',template='plotly_white')
        else: fig=px.bar(df,x=df.columns[0],y=df.columns[-1],title='Failure summary',template='plotly_white')
        Path(out_path).parent.mkdir(parents=True,exist_ok=True); fig.write_html(out_path,include_plotlyjs='cdn',full_html=True); return Path(out_path)
    except Exception:
        return write_html(out_path,'Failure heatmap',df.head(50).to_html(index=False, escape=True))
