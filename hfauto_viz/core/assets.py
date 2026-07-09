from __future__ import annotations

import json
from pathlib import Path

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
THREEDMOL_CDN = "https://3Dmol.org/build/3Dmol-min.js"
CYTOSCAPE_CDN = "https://unpkg.com/cytoscape@3.30.2/dist/cytoscape.min.js"

_ASSETS = {
    "plotly": ("plotly.min.js", PLOTLY_CDN, "Plotly"),
    "3dmol": ("3Dmol-min.js", THREEDMOL_CDN, "$3Dmol"),
    "cytoscape": ("cytoscape.min.js", CYTOSCAPE_CDN, "cytoscape"),
}


def normalize_mode(mode: str | None) -> str:
    mode = (mode or "cdn").lower().strip()
    if mode in {"bundle", "offline"}:
        return "local"
    if mode not in {"cdn", "local", "none", "auto"}:
        return "cdn"
    return mode


def asset_filename(kind: str) -> str:
    return _ASSETS[kind][0]


def script_tag(kind: str, mode: str = "cdn", prefix: str = "assets") -> str:
    mode = normalize_mode(mode)
    if kind not in _ASSETS:
        return f"<!-- unknown hfauto-viz asset kind: {kind} -->"
    fname, cdn, _global = _ASSETS[kind]
    if mode == "none":
        return f"<!-- {kind} library intentionally not loaded -->"
    if mode == "local":
        return f'<script src="{prefix.rstrip("/")}/{fname}"></script>'
    if mode == "auto":
        # Local first, then CDN fallback. The onerror fallback keeps reviewed HTML useful
        # when a site has not yet vendored third-party JavaScript assets.
        return (
            f'<script src="{prefix.rstrip("/")}/{fname}" '
            f'onerror="var s=document.createElement(\'script\');s.src=\'{cdn}\';document.head.appendChild(s);"></script>'
        )
    return f'<script src="{cdn}"></script>'


def css_tag(prefix: str = "assets") -> str:
    return f'<link rel="stylesheet" href="{prefix.rstrip("/")}/hfauto_viz.css">'


def _placeholder_js(kind: str) -> str:
    fname, cdn, global_name = _ASSETS[kind]
    return f"""// hfauto-viz placeholder for {fname}\n// Production/offline deployments may replace this file with the real library.\n// Expected global: {global_name}\nconsole.warn('hfauto-viz placeholder loaded for {fname}; interactive {kind} viewer will use HTML fallback unless the real library is installed. CDN reference: {cdn}');\n"""


def write_placeholder_assets(out_dir: str | Path) -> Path:
    assets = Path(out_dir) / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    (assets / "hfauto_viz.css").write_text(
        """
body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;margin:0;background:#fafafa;color:#222}
.card{background:#fff;border:1px solid #ddd;border-radius:10px;padding:1rem;margin:1rem 0;box-shadow:0 1px 2px rgba(0,0,0,.04)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:1rem}.viewer{width:100%;height:520px;border:1px solid #ddd;border-radius:8px;background:white}.warn{border-left:5px solid #d08b00;background:#fff8e6;padding:.7rem}.data-table{border-collapse:collapse;width:100%;font-size:.9rem}.data-table th,.data-table td{border:1px solid #e0e0e0;padding:.35rem .5rem;vertical-align:top}.thumbnail{max-width:140px;max-height:90px}.viz-split{display:grid;grid-template-columns:2fr 1fr;gap:1rem}.badge{display:inline-block;border:1px solid #bbb;border-radius:999px;padding:.15rem .5rem;margin:.1rem;background:#f7f7f7}.svg-box{overflow:auto;background:white;border:1px solid #ddd;border-radius:8px;padding:.5rem}.muted{color:#666}.toolbar{display:flex;gap:.5rem;flex-wrap:wrap;margin:.5rem 0}.toolbar button{padding:.25rem .6rem;border:1px solid #aaa;border-radius:6px;background:#f9f9f9}
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (assets / "hfauto_viz.js").write_text(
        """
function hfautoVizEscapeHtml(s){return String(s ?? '').replace(/[&<>\"']/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}[m]));}
function hfautoVizShowFallback(id, text){const el=document.getElementById(id); if(el){el.innerHTML='<pre>'+hfautoVizEscapeHtml(text)+'</pre>';}}
""".strip()
        + "\n",
        encoding="utf-8",
    )
    for kind in _ASSETS:
        (assets / asset_filename(kind)).write_text(_placeholder_js(kind), encoding="utf-8")
    manifest = {
        "schema_version": "hfauto_viz.assets.v1",
        "note": "Placeholder assets are written for offline packaging. Replace Plotly/3Dmol/Cytoscape files with vendor-approved production copies for full interactivity without CDN.",
        "assets": {k: {"filename": v[0], "cdn": v[1], "expected_global": v[2]} for k, v in _ASSETS.items()},
    }
    (assets / "asset_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (assets / "README.txt").write_text(
        "hfauto-viz asset bundle. Placeholder JS files are intentionally small and safe for repository storage.\n"
        "For a fully offline production report, replace plotly.min.js, 3Dmol-min.js and cytoscape.min.js with approved vendor copies.\n",
        encoding="utf-8",
    )
    (assets / "asset_status.html").write_text(
        "<html><body><h1>hfauto-viz assets</h1><p>Placeholder asset bundle generated. Replace third-party JS files for fully offline interactivity.</p></body></html>",
        encoding="utf-8",
    )
    return assets


write_default_assets = write_placeholder_assets
