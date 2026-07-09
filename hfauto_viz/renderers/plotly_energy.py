from __future__ import annotations
from pathlib import Path
from html import escape
from typing import Any
import pandas as pd
from hfauto_viz.core.html import write_html

def _n(row: dict[str, Any], key: str, default=0.0) -> float:
    try:
        v=row.get(key)
        if v is None or v=='': return default
        return float(v)
    except Exception: return default

def energy_points(row: dict[str, Any], process: bool | None = None, mode: str | None = None):
    # Compatibility: mode="process" returns list[{state,G}], while process=False
    # returns (labels, values) for older callers.
    use_process = True if process is None else bool(process)
    if mode is not None:
        use_process = (mode == "process")
    assoc = _n(row, 'delta_G_assoc_pressure_corrected_kcal_mol' if use_process else 'delta_G_assoc_standard_kcal_mol', _n(row, 'delta_G_assoc_kcal_mol', 0.0))
    act = _n(row, 'delta_G_act_process_kcal_mol' if use_process else 'delta_G_act_kcal_mol', _n(row, 'delta_G_act_kcal_mol', 0.0))
    ion = _n(row, 'delta_G_ionpair_kcal_mol', 0.0)
    labels=['B + (HF)n','B···(HF)n','TS','BH+···F(HF)n-1−']
    vals=[0.0, assoc, assoc+act, assoc+ion]
    if mode is not None:
        return [{'state': lab, 'G': val} for lab, val in zip(labels, vals)]
    return labels, vals

def _fallback_svg(labels, y):
    width,height=760,330; mn=min(y+[0]); mx=max(y+[0])
    if mn==mx: mn-=1; mx+=1
    def xx(i): return 70+i*(width-140)/max(len(y)-1,1)
    def yy(v): return height-50-(v-mn)*(height-100)/(mx-mn)
    parts=[f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}'><rect width='100%' height='100%' fill='white'/>"]
    for i,(lab,v) in enumerate(zip(labels,y)):
        x=xx(i); yv=yy(v); parts.append(f"<line x1='{x-30}' x2='{x+30}' y1='{yv}' y2='{yv}' stroke='#145da0' stroke-width='4'/><text x='{x}' y='{yv-9}' text-anchor='middle' font-size='12'>{v:.2f}</text><text x='{x}' y='{height-20}' text-anchor='middle' font-size='11'>{escape(lab)}</text>")
        if i < len(labels)-1: parts.append(f"<line x1='{x+30}' y1='{yv}' x2='{xx(i+1)-30}' y2='{yy(y[i+1])}' stroke='#999' stroke-dasharray='4 4'/>")
    return ''.join(parts)+'</svg>'

def plotly_energy_div(row: dict[str, Any], title: str='Energy profile') -> str:
    labels, yp = energy_points(row, process=True); _, ys = energy_points(row, process=False)
    try:
        import plotly.graph_objects as go, plotly.io as pio
        fig=go.Figure(); fig.add_trace(go.Scatter(x=labels,y=yp,mode='lines+markers+text',name='process',text=[f'{v:.2f}' for v in yp],textposition='top center')); fig.add_trace(go.Scatter(x=labels,y=ys,mode='lines+markers',name='standard'))
        fig.update_layout(title=title,yaxis_title='Relative G / kcal mol⁻¹',template='plotly_white')
        return pio.to_html(fig, include_plotlyjs='cdn', full_html=False)
    except Exception:
        return _fallback_svg(labels, yp)

def write_energy_profile(row: dict[str, Any], out_path: str | Path, title: str | None=None) -> Path:
    title=title or f"Energy profile: {row.get('reaction_id','reaction')}"
    fields=['reaction_id','mol_id','site_id','hf_n','T_K','quality_tier','confidence_score','delta_G_assoc_pressure_corrected_kcal_mol','delta_G_act_kcal_mol','k_corrected_s-1']
    table='<table class="data-table">'+''.join(f"<tr><th>{escape(k)}</th><td>{escape(str(row.get(k,'')))}</td></tr>" for k in fields)+'</table>'
    return write_html(out_path,title,f"<div class='card'>{plotly_energy_div(row,title)}</div><div class='card'><h2>Key values</h2>{table}</div>")

def write_energy_comparison(df: pd.DataFrame, out_path: str | Path, title: str='Top candidate energy profiles', top_n: int=12) -> Path:
    if df.empty: return write_html(out_path,title,"<div class='warn'>No reaction table available.</div>")
    rows=df.head(top_n).to_dict('records')
    try:
        import plotly.graph_objects as go, plotly.io as pio
        fig=go.Figure()
        for row in rows:
            labels,y=energy_points(row, process=True); fig.add_trace(go.Scatter(x=labels,y=y,mode='lines+markers',name=str(row.get('mol_id',row.get('reaction_id')))))
        fig.update_layout(title=title,yaxis_title='Relative G / kcal mol⁻¹',template='plotly_white')
        div=pio.to_html(fig,include_plotlyjs='cdn',full_html=False)
    except Exception:
        labels,y=energy_points(rows[0],process=True); div=_fallback_svg(labels,y)
    return write_html(out_path,title,f"<div class='card'>{div}</div><div class='card'>{df.head(top_n).to_html(index=False, escape=True)}</div>")

# Compatibility extension for Phase 8 tests and public API.
_energy_points_tuple = energy_points

def energy_points(row: dict[str, Any], process: bool = True, mode: str | None = None):
    labels, vals = _energy_points_tuple(row, process=(mode == 'process') if mode is not None else process)
    if mode is not None:
        return [{'state': s, 'G': v} for s, v in zip(labels, vals)]
    return labels, vals

def write_top_energy_overlay(df: pd.DataFrame, out_path: str | Path, top_n: int = 12) -> Path:
    return write_energy_comparison(df, out_path, top_n=top_n)

def write_energy_profile_static_svg(row: dict[str, Any], out_path: str | Path, process: bool = True) -> dict[str, Path]:
    """Write a lightweight static SVG and CSV next to an interactive Plotly energy profile."""
    import pandas as pd
    p = Path(out_path); p.parent.mkdir(parents=True, exist_ok=True)
    labels, vals = energy_points(row, process=process)
    p.write_text(_fallback_svg(labels, vals), encoding='utf-8')
    csv_path = p.with_suffix('.csv')
    pd.DataFrame({'state': labels, 'relative_G_kcal_mol': vals}).to_csv(csv_path, index=False)
    return {'energy_profile_svg': p, 'energy_profile_csv': csv_path}
