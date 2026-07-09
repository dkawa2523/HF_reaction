
from __future__ import annotations
from typing import Any

def atom_highlights_for_reaction(reaction_data: dict[str, Any] | None) -> dict[str,int]:
    if not reaction_data: return {}
    coord = reaction_data.get('reaction_coordinate') or {}; atoms = coord.get('atoms', {}) if isinstance(coord, dict) else {}
    out={}
    for key in ['base_atom','transfer_h','leaving_f']:
        if atoms.get(key) is not None:
            try: out[key]=int(atoms[key])
            except Exception: pass
    return out
