from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_viz.core.assets import script_tag
from hfauto_viz.core.html import esc, write_html


def cytoscape_elements(reaction_df: pd.DataFrame, top_n: int = 50) -> list[dict[str, Any]]:
    rows = reaction_df.head(top_n).to_dict('records') if reaction_df is not None and not reaction_df.empty else []
    elements: list[dict[str, Any]] = []
    seen: set[str] = set()
    def node(node_id: str, label: str, state: str, row: dict[str, Any]):
        if node_id in seen:
            return
        seen.add(node_id)
        elements.append({'data': {'id': node_id, 'label': label, 'state': state, 'mol_id': row.get('mol_id'), 'quality_tier': row.get('quality_tier')}})
    for row in rows:
        rid = str(row.get('reaction_id') or f"rxn_{len(elements)}")
        mid = str(row.get('mol_id') or 'candidate')
        rc = str(row.get('reactant_species_id') or f'{rid}_rc')
        ts = str(row.get('ts_species_id') or f'{rid}_ts')
        ip = str(row.get('product_species_id') or f'{rid}_ip')
        node(mid, mid, 'candidate', row)
        node(rc, f"B···(HF){row.get('hf_n','')}", 'reactant_complex', row)
        node(ts, 'TS', 'transition_state', row)
        node(ip, 'ion pair', 'ion_pair', row)
        for source,target,label in [(mid,rc,'association'),(rc,ts,f"ΔG‡ {row.get('delta_G_act_kcal_mol','')}") ,(ts,ip,'PT')]:
            elements.append({'data': {'id': f'{rid}_{source}_{target}', 'source': source, 'target': target, 'label': label, 'reaction_id': rid, 'delta_G_act_kcal_mol': row.get('delta_G_act_kcal_mol')}})
    return elements


def write_cytoscape_network(reaction_df: pd.DataFrame, out_dir: str | Path, top_n: int = 50, library_mode: str = 'cdn', fallback_svg_rel: str = 'reaction_network.svg') -> dict[str, Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    elements = cytoscape_elements(reaction_df, top_n=top_n)
    data_path = out / 'reaction_network_interactive.json'
    data_path.write_text(json.dumps({'elements': elements}, ensure_ascii=False, indent=2), encoding='utf-8')
    asset_prefix = '../assets'
    script = script_tag('cytoscape', library_mode, prefix=asset_prefix)
    # Cytoscape is optional. If absent, a static Graphviz SVG iframe and JSON preview still help reviewers.
    body = f"""
<div class='card'>
  <h2>Interactive reaction network</h2>
  <p class='muted'>Click nodes/edges to inspect metadata. If Cytoscape.js is not available, the static Graphviz network below is shown.</p>
  <div id='cy' style='width:100%;height:680px;border:1px solid #ddd;border-radius:8px;background:white'></div>
  <div id='cyinfo' class='card'></div>
  <details><summary>Network JSON</summary><pre>{esc(json.dumps({'elements': elements[:20]}, ensure_ascii=False, indent=2))}</pre></details>
</div>
{script}
<script src='{asset_prefix}/hfauto_viz.js'></script>
<script>
const hfautoNetwork={json.dumps(elements)};
function initNetwork(){{
  const cyEl=document.getElementById('cy');
  if(typeof cytoscape==='undefined'){{
    cyEl.innerHTML='<iframe src="{fallback_svg_rel}" style="width:100%;height:650px;border:0"></iframe><p class="warn">Cytoscape.js not loaded; showing static Graphviz fallback.</p>';
    return;
  }}
  const cy=cytoscape({{container:cyEl,elements:hfautoNetwork,layout:{{name:'breadthfirst',directed:true,padding:20}},
    style:[
      {{selector:'node',style:{{'label':'data(label)','background-color':'#e7f1ff','border-color':'#4b83c4','border-width':1,'text-valign':'center','text-wrap':'wrap','text-max-width':120,'font-size':10}}}},
      {{selector:'node[state="transition_state"]',style:{{'shape':'hexagon','background-color':'#fff3cd','border-color':'#d08b00'}}}},
      {{selector:'node[state="ion_pair"]',style:{{'background-color':'#e9f7ef','border-color':'#2e8b57'}}}},
      {{selector:'edge',style:{{'curve-style':'bezier','target-arrow-shape':'triangle','line-color':'#999','target-arrow-color':'#999','label':'data(label)','font-size':9,'text-background-color':'white','text-background-opacity':0.8}}}}
    ]}});
  cy.on('tap','node,edge',evt=>{{document.getElementById('cyinfo').innerHTML='<h3>Selected</h3><pre>'+hfautoVizEscapeHtml(JSON.stringify(evt.target.data(),null,2))+'</pre>';}});
}}
if(document.readyState==='loading'){{document.addEventListener('DOMContentLoaded',initNetwork);}}else{{initNetwork();}}
</script>
"""
    html_path = write_html(out / 'reaction_network_interactive.html', 'Interactive reaction network', body)
    return {'cytoscape_html': html_path, 'cytoscape_json': data_path}


def write_cytoscape_network_linked(run, reaction_df: pd.DataFrame, out_dir: str | Path, top_n: int = 50, library_mode: str = 'cdn', fallback_svg_rel: str = 'reaction_network.svg') -> dict[str, Path]:
    """Write an interactive network with a linked 3Dmol structure pane.

    This remains robust in offline reports: if Cytoscape or 3Dmol vendor JS is
    missing, the report falls back to static SVG/metadata/XYZ text rather than
    failing.
    """
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    elements = cytoscape_elements(reaction_df, top_n=top_n)
    sp_idx = run.species_index()
    xyz_by_node: dict[str, str] = {}
    for row in (reaction_df.head(top_n).to_dict('records') if reaction_df is not None and not reaction_df.empty else []):
        for sid in [row.get('reactant_species_id'), row.get('ts_species_id'), row.get('product_species_id')]:
            if not sid:
                continue
            art = sp_idx.get(str(sid))
            p = run.artifact_xyz_path(art)
            if p and p.exists():
                try:
                    xyz_by_node[str(sid)] = p.read_text(encoding='utf-8')
                except Exception:
                    pass
    data = {'elements': elements, 'xyz_by_node': xyz_by_node}
    data_path = out / 'reaction_network_linked.json'
    data_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    asset_prefix = '../assets'
    body = f"""
<div class='card'>
  <h2>Linked reaction network + 3D structure viewer</h2>
  <p class='muted'>Click species nodes to load the corresponding XYZ into the 3D viewer. Click edges to inspect reaction metadata.</p>
  <div class='viz-split'>
    <div><div id='cy_linked' style='width:100%;height:720px;border:1px solid #ddd;border-radius:8px;background:white'></div></div>
    <div><h3>Selected structure</h3><div id='linked3d' class='viewer' style='height:430px'></div><div id='linkedinfo' class='card'></div></div>
  </div>
  <details><summary>Network JSON preview</summary><pre>{esc(json.dumps({'elements': elements[:20]}, ensure_ascii=False, indent=2))}</pre></details>
</div>
{script_tag('cytoscape', library_mode, prefix=asset_prefix)}
{script_tag('3dmol', library_mode, prefix=asset_prefix)}
<script src='{asset_prefix}/hfauto_viz.js'></script>
<script>
const linkedNetwork={json.dumps(elements)};
const linkedXYZ={json.dumps(xyz_by_node)};
let linkedViewer=null;
function loadLinked3D(nodeId){{
  const xyz=linkedXYZ[nodeId];
  const info=document.getElementById('linkedinfo');
  if(info) info.innerHTML='<h3>'+hfautoVizEscapeHtml(nodeId)+'</h3><pre>'+hfautoVizEscapeHtml(JSON.stringify((linkedNetwork.find(e=>e.data && e.data.id===nodeId)||{{}}).data||{{}},null,2))+'</pre>';
  if(!xyz){{hfautoVizShowFallback('linked3d','No XYZ geometry linked for '+nodeId); return;}}
  if(typeof $3Dmol==='undefined'){{hfautoVizShowFallback('linked3d','3Dmol.js not loaded\\n'+xyz); return;}}
  const el=document.getElementById('linked3d'); el.innerHTML='';
  linkedViewer=$3Dmol.createViewer(el,{{backgroundColor:'white'}});
  linkedViewer.addModel(xyz,'xyz'); linkedViewer.setStyle({{}},{{stick:{{radius:0.16}},sphere:{{scale:0.27}}}}); linkedViewer.zoomTo(); linkedViewer.render();
}}
function initLinkedNetwork(){{
  const cyEl=document.getElementById('cy_linked');
  if(typeof cytoscape==='undefined'){{
    cyEl.innerHTML='<iframe src="{fallback_svg_rel}" style="width:100%;height:690px;border:0"></iframe><p class="warn">Cytoscape.js not loaded; showing static Graphviz fallback.</p>';
    return;
  }}
  const cy=cytoscape({{container:cyEl,elements:linkedNetwork,layout:{{name:'breadthfirst',directed:true,padding:20}},
    style:[
      {{selector:'node',style:{{'label':'data(label)','background-color':'#e7f1ff','border-color':'#4b83c4','border-width':1,'text-valign':'center','text-wrap':'wrap','text-max-width':120,'font-size':10}}}},
      {{selector:'node[state="transition_state"]',style:{{'shape':'hexagon','background-color':'#fff3cd','border-color':'#d08b00'}}}},
      {{selector:'node[state="ion_pair"]',style:{{'background-color':'#e9f7ef','border-color':'#2e8b57'}}}},
      {{selector:'edge',style:{{'curve-style':'bezier','target-arrow-shape':'triangle','line-color':'#999','target-arrow-color':'#999','label':'data(label)','font-size':9,'text-background-color':'white','text-background-opacity':0.8}}}}
    ]}});
  cy.on('tap','node',evt=>loadLinked3D(evt.target.data().id));
  cy.on('tap','edge',evt=>{{document.getElementById('linkedinfo').innerHTML='<h3>Selected edge</h3><pre>'+hfautoVizEscapeHtml(JSON.stringify(evt.target.data(),null,2))+'</pre>';}});
  const first=linkedNetwork.find(e=>e.data && !e.data.source && linkedXYZ[e.data.id]); if(first) loadLinked3D(first.data.id);
}}
if(document.readyState==='loading'){{document.addEventListener('DOMContentLoaded',initLinkedNetwork);}}else{{initLinkedNetwork();}}
</script>
"""
    html_path = write_html(out / 'reaction_network_linked.html', 'Linked reaction network and 3D viewer', body)
    return {'cytoscape_linked_html': html_path, 'cytoscape_linked_json': data_path}
