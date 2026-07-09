from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_viz.core.assets import script_tag
from hfauto_viz.core.html import dataframe_to_html, esc, write_html
from hfauto_viz.structure.geometry import coordinate_values, reaction_coordinate_atoms
from hfauto_viz.structure.xyz import Frame


def reaction_coordinate_dataframe(frames: list[tuple[str, Frame]] | list[Frame], reaction_data: dict[str, Any] | None = None) -> pd.DataFrame:
    reaction_data = reaction_data or {}
    first = frames[0][1] if frames and isinstance(frames[0], tuple) else (frames[0] if frames else [])
    atoms = reaction_coordinate_atoms(reaction_data, first)
    rows=[]
    for i,item in enumerate(frames):
        label, frame = item if isinstance(item, tuple) else (f'frame_{i:03d}', item)
        vals=coordinate_values(frame, atoms)
        rows.append({'frame': i, 'label': label, **vals, 'base_atom': atoms.get('base_atom'), 'transfer_h': atoms.get('transfer_h'), 'leaving_f': atoms.get('leaving_f')})
    return pd.DataFrame(rows)


def write_reaction_coordinate_plot(
    frames: list[tuple[str, Frame]] | list[Frame],
    reaction_data: dict[str, Any] | None,
    out_path: str | Path,
    library_mode: str = 'cdn',
    asset_prefix: str = '../../assets',
) -> dict[str, Path]:
    out=Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    df=reaction_coordinate_dataframe(frames, reaction_data)
    csv_path=out.with_suffix('.csv')
    df.to_csv(csv_path,index=False)
    if df.empty:
        html=write_html(out,'Reaction coordinate','<div class="warn">No frames for reaction coordinate.</div>')
        return {'html':html,'csv':csv_path}
    try:
        import plotly.graph_objects as go, plotly.io as pio
        fig=go.Figure()
        if 'r_BH_A' in df:
            fig.add_trace(go.Scatter(x=df['frame'], y=df['r_BH_A'], mode='lines+markers', name='r(B-H) / Å'))
        if 'r_HF_A' in df:
            fig.add_trace(go.Scatter(x=df['frame'], y=df['r_HF_A'], mode='lines+markers', name='r(H-F) / Å'))
        if 'q_BH_minus_HF_A' in df:
            fig.add_trace(go.Scatter(x=df['frame'], y=df['q_BH_minus_HF_A'], mode='lines+markers', name='q=r(B-H)-r(H-F) / Å'))
        fig.update_layout(title='Reaction coordinate along path frames', xaxis_title='Frame', yaxis_title='Distance / Å', template='plotly_white')
        div = script_tag('plotly', library_mode, prefix=asset_prefix) + pio.to_html(fig, include_plotlyjs=False, full_html=False)
    except Exception as exc:
        div = f"<div class='warn'>Plotly unavailable: {esc(exc)}</div>" + dataframe_to_html(df)
    body = f"<div class='card'>{div}</div><div class='card'><h2>Coordinate data</h2>{dataframe_to_html(df)}</div>"
    html=write_html(out,'Reaction coordinate',body)
    return {'html':html,'csv':csv_path}
