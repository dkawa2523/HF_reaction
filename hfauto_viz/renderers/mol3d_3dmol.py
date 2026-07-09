from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from hfauto_viz.core.assets import script_tag
from hfauto_viz.core.html import esc, write_html
from hfauto_viz.structure.xyz import Frame, frame_to_xyz, frames_to_xyz, read_xyz_frames


def _viewer_script(
    div_id: str,
    xyz_var: str,
    is_multiframe: bool,
    n_frames: int,
    annotations: list[dict[str, Any]] | None = None,
) -> str:
    annotations = annotations or []
    return f"""
<script>
let viewer_{div_id}=null, model_{div_id}=null, timer_{div_id}=null, cur_{div_id}=0;
const annotations_{div_id}={json.dumps(annotations)};
function annList_{div_id}(){{
  return annotations_{div_id}.filter(a=>a.frame===undefined || a.frame===null || a.frame==='all' || a.frame===cur_{div_id});
}}
function annotate_{div_id}(){{
  if(!viewer_{div_id}) return;
  for(const ann of annList_{div_id}()){{
    if(ann.type==='line' && ann.start && ann.end){{viewer_{div_id}.addLine({{start:ann.start,end:ann.end,color:ann.color||'black',dashed:!!ann.dashed,linewidth:ann.linewidth||2}});}}
    if(ann.type==='arrow' && ann.start && ann.end){{viewer_{div_id}.addArrow({{start:ann.start,end:ann.end,color:ann.color||'purple',radius:ann.radius||0.06,radiusRatio:ann.radiusRatio||1.7,mid:ann.mid||0.75}});}}
    if(ann.type==='sphere' && ann.center){{viewer_{div_id}.addSphere({{center:ann.center,radius:ann.radius||0.30,color:ann.color||'orange',alpha:ann.alpha===undefined?0.5:ann.alpha}});}}
    if(ann.type==='label' && ann.position){{viewer_{div_id}.addLabel(ann.text||'',{{position:ann.position,backgroundColor:ann.backgroundColor||'white',fontColor:ann.color||'black',fontSize:ann.fontSize||12}});}}
  }}
}}
function style_{div_id}(mode){{
  if(!viewer_{div_id}) return;
  viewer_{div_id}.setStyle({{}},{{}});
  if(mode==='stick') viewer_{div_id}.setStyle({{}},{{stick:{{radius:0.18}}}});
  else if(mode==='sphere') viewer_{div_id}.setStyle({{}},{{sphere:{{scale:0.42}}}});
  else viewer_{div_id}.setStyle({{}},{{stick:{{radius:0.16}},sphere:{{scale:0.27}}}});
  viewer_{div_id}.removeAllShapes(); viewer_{div_id}.removeAllLabels(); annotate_{div_id}(); viewer_{div_id}.render();
}}
function render_{div_id}(){{
  if(typeof $3Dmol==='undefined'){{hfautoVizShowFallback('{div_id}', '3Dmol.js not loaded\\n'+{xyz_var}); return;}}
  viewer_{div_id}=$3Dmol.createViewer(document.getElementById('{div_id}'),{{backgroundColor:'white'}});
  model_{div_id}=viewer_{div_id}.addModel({xyz_var},'xyz'{', {multimodel:true}' if is_multiframe else ''});
  viewer_{div_id}.setStyle({{}},{{stick:{{radius:0.16}},sphere:{{scale:0.27}}}});
  viewer_{div_id}.zoomTo(); annotate_{div_id}(); viewer_{div_id}.render();
}}
function show_{div_id}(i){{
  cur_{div_id}=Math.max(0,Math.min(i,{max(n_frames-1,0)}));
  if(model_{div_id} && model_{div_id}.setFrame){{model_{div_id}.setFrame(cur_{div_id}); viewer_{div_id}.removeAllLabels(); viewer_{div_id}.removeAllShapes(); annotate_{div_id}(); viewer_{div_id}.render();}}
  const lab=document.getElementById('{div_id}_label'); if(lab) lab.innerText='frame '+cur_{div_id};
  const sl=document.getElementById('{div_id}_slider'); if(sl) sl.value=cur_{div_id};
}}
function play_{div_id}(){{pause_{div_id}(); timer_{div_id}=setInterval(()=>show_{div_id}(cur_{div_id}+1>{max(n_frames-1,0)}?0:cur_{div_id}+1),650);}}
function pause_{div_id}(){{if(timer_{div_id}) clearInterval(timer_{div_id}); timer_{div_id}=null;}}
if(document.readyState==='loading'){{document.addEventListener('DOMContentLoaded',render_{div_id});}}else{{render_{div_id}();}}
</script>"""


def _toolbar(div_id: str, multiframe: bool = False, n_frames: int = 1) -> str:
    buttons = (
        f"<button onclick=\"style_{div_id}('ballstick')\">ball+stick</button>"
        f"<button onclick=\"style_{div_id}('stick')\">stick</button>"
        f"<button onclick=\"style_{div_id}('sphere')\">sphere</button>"
    )
    if multiframe:
        buttons += (
            f"<button onclick='play_{div_id}()'>Play</button>"
            f"<button onclick='pause_{div_id}()'>Pause</button>"
            f"<input id='{div_id}_slider' type='range' min='0' max='{max(n_frames-1,0)}' value='0' oninput='show_{div_id}(parseInt(this.value))'>"
            f"<span id='{div_id}_label' class='badge'></span>"
        )
    return f"<div class='toolbar'>{buttons}</div>"


def html_for_structure(
    frame: Frame | str | None,
    div_id: str = 'viewer',
    title: str | None = None,
    library_mode: str = 'cdn',
    asset_prefix: str = 'assets',
    annotations: list[dict[str, Any]] | None = None,
) -> str:
    xyz = frame_to_xyz(frame, title or '') if isinstance(frame, list) else (frame or '0\nmissing\n')
    var = f"xyz_{div_id}"
    body = [
        script_tag('3dmol', library_mode, prefix=asset_prefix),
        f"<script src='{asset_prefix}/hfauto_viz.js'></script>" if library_mode in {'local','auto'} else "",
        _toolbar(div_id),
        f"<div id='{div_id}' class='viewer'></div>",
        f"<script>const {var}={json.dumps(xyz)};</script>",
        _viewer_script(div_id, var, False, 1, annotations=annotations),
        f"<details><summary>XYZ source</summary><pre>{esc(xyz[:20000])}</pre></details>",
    ]
    return "\n".join(body)


def write_structure_html(
    frame: Frame | str | None,
    out_path: str | Path,
    title: str,
    library_mode: str = 'cdn',
    asset_prefix: str = 'assets',
    annotations: list[dict[str, Any]] | None = None,
) -> Path:
    return write_html(out_path, title, html_for_structure(frame, div_id='viewer', title=title, library_mode=library_mode, asset_prefix=asset_prefix, annotations=annotations))


def write_animation_html(
    frames: list[tuple[str, Frame]] | list[Frame],
    out_path: str | Path,
    title: str = 'Path animation',
    library_mode: str = 'cdn',
    asset_prefix: str = 'assets',
    annotations: list[dict[str, Any]] | None = None,
) -> Path:
    xyz = frames_to_xyz(frames)
    p = Path(out_path); p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix == '.xyz':
        p.write_text(xyz, encoding='utf-8'); html_path = p.with_suffix('.html')
    else:
        html_path = p
    div_id='pathviewer'; var='xyz_pathviewer'; n_frames=len(frames)
    body = (
        f"{script_tag('3dmol', library_mode, prefix=asset_prefix)}"
        + (f"<script src='{asset_prefix}/hfauto_viz.js'></script>" if library_mode in {'local','auto'} else "")
        + _toolbar(div_id, multiframe=True, n_frames=n_frames)
        + f"<div id='{div_id}' class='viewer'></div><script>const {var}={json.dumps(xyz)};</script>"
        + _viewer_script(div_id, var, True, n_frames, annotations=annotations)
        + f"<details><summary>Multi-frame XYZ</summary><pre>{esc(xyz[:40000])}</pre></details>"
    )
    return write_html(html_path, title, body)


def write_xyz_viewer(xyz_path: str | Path, out_path: str | Path, title: str, subtitle: str | None = None, library_mode: str = 'cdn') -> Path:
    frames = read_xyz_frames(xyz_path)
    return write_structure_html(frames[0] if frames else None, out_path, title, library_mode=library_mode)


def write_xyz_animation(xyz_path: str | Path, out_path: str | Path, title: str = 'Path animation', subtitle: str | None = None, library_mode: str = 'cdn') -> Path:
    frames = read_xyz_frames(xyz_path)
    return write_animation_html(frames, out_path, title=title, library_mode=library_mode)
