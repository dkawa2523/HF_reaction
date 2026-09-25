from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

from hfauto_viz.core.html import esc, write_html
from hfauto_viz.structure.xyz import Frame


def _frame_payload(frame: Frame, name: str):
    return {
        'size': len(frame),
        'names': [a.symbol for a in frame],
        'x': [a.x for a in frame],
        'y': [a.y for a in frame],
        'z': [a.z for a in frame],
        'name': name,
    }


def _safe_json_value(v):
    try:
        if pd.isna(v):
            return None
    except Exception:
        pass
    return v


def write_chemiscope_like(df: pd.DataFrame, structures: dict[str, Frame], out_path: str | Path) -> Path:
    p = Path(out_path); p.parent.mkdir(parents=True, exist_ok=True)
    rows = df.to_dict('records') if df is not None and not df.empty else []
    struct_list = []
    if structures:
        for sid, frame in structures.items():
            struct_list.append(_frame_payload(frame, sid))
    else:
        struct_list = [{'size': 0, 'names': [], 'x': [], 'y': [], 'z': [], 'name': str(r.get('reaction_id', i))} for i, r in enumerate(rows)]
    props = {}
    if rows:
        limit = len(struct_list)
        for k in rows[0]:
            props[k] = {'target': 'structure', 'values': [_safe_json_value(r.get(k)) for r in rows[:limit]]}
    payload = {
        'meta': {
            'name': 'HF reactivity screening',
            'description': 'hfauto Chemiscope-compatible fallback export. Native Chemiscope can consume the same structures and properties when installed.',
        },
        'structures': struct_list,
        'properties': props,
    }
    if p.suffix == '.gz':
        with gzip.open(p, 'wt', encoding='utf-8') as f:
            json.dump(payload, f, ensure_ascii=False)
    else:
        p.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    return p


def write_chemiscope_index(bundle_path: str | Path, out_path: str | Path | None = None) -> Path:
    bundle = Path(bundle_path)
    out = Path(out_path) if out_path else bundle.with_suffix('').with_suffix('.html') if bundle.name.endswith('.json.gz') else bundle.with_suffix('.html')
    body = f"""
<div class='card'>
  <h2>Chemiscope export</h2>
  <p>Structure-property bundle: <a href='{esc(bundle.name)}'>{esc(bundle.name)}</a></p>
  <p class='muted'>This export is intentionally viewer-independent. If the optional <code>chemiscope</code> Python package is available, the JSON bundle can be converted to a native Chemiscope HTML viewer in a production deployment.</p>
  <pre>chemiscope view {esc(bundle.name)}</pre>
</div>
"""
    return write_html(out, 'Chemiscope export', body)


def write_chemiscope_bundle(df: pd.DataFrame, structures: dict[str, Frame], out_path: str | Path) -> dict[str, Path]:
    bundle = write_chemiscope_like(df, structures, out_path)
    index = write_chemiscope_index(bundle)
    return {'chemiscope_json_gz': bundle, 'chemiscope_html': index}
