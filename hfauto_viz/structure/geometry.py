from __future__ import annotations

from math import sqrt
from typing import Any

from hfauto_viz.structure.xyz import Atom, Frame


def distance(a: Atom, b: Atom) -> float:
    return sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


def safe_atom(frame: Frame, idx: int | None) -> Atom | None:
    if idx is None:
        return None
    try:
        if 0 <= int(idx) < len(frame):
            return frame[int(idx)]
    except Exception:
        return None
    return None


def detect_pt_atoms(frame: Frame) -> dict[str, int | None]:
    """Heuristic fallback for B-H-F proton-transfer visualization.

    It finds the shortest H-F pair and then the nearest heavy non-F atom to H.
    This is only used for visualization if the ReactionRecord lacks atom indices.
    """
    if not frame:
        return {"base_atom": None, "transfer_h": None, "leaving_f": None}
    h_indices=[i for i,a in enumerate(frame) if a.symbol.upper()=='H']
    f_indices=[i for i,a in enumerate(frame) if a.symbol.upper()=='F']
    best_h=None; best_f=None; best=1e9
    for hi in h_indices:
        for fi in f_indices:
            d=distance(frame[hi],frame[fi])
            if d<best:
                best=d; best_h=hi; best_f=fi
    base=None; best_b=1e9
    if best_h is not None:
        for i,a in enumerate(frame):
            if i in {best_h,best_f} or a.symbol.upper() in {'H','F'}:
                continue
            d=distance(frame[i],frame[best_h])
            if d<best_b:
                best_b=d; base=i
    return {"base_atom": base, "transfer_h": best_h, "leaving_f": best_f}


def reaction_coordinate_atoms(reaction_data: dict[str, Any], frame: Frame | None = None) -> dict[str, int | None]:
    rc = reaction_data.get('reaction_coordinate') or {}
    atoms = rc.get('atoms') or {}
    out = {
        'base_atom': atoms.get('base_atom') or atoms.get('B') or atoms.get('base') or reaction_data.get('base_atom'),
        'transfer_h': atoms.get('transfer_h') or atoms.get('H') or atoms.get('proton') or reaction_data.get('transfer_h'),
        'leaving_f': atoms.get('leaving_f') or atoms.get('F') or atoms.get('fluoride') or reaction_data.get('leaving_f'),
    }
    cleaned={k:(int(v) if v is not None and str(v) not in {'','nan','None'} else None) for k,v in out.items()}
    if frame and any(v is None for v in cleaned.values()):
        detected=detect_pt_atoms(frame)
        for k,v in detected.items():
            if cleaned.get(k) is None:
                cleaned[k]=v
    return cleaned


def coordinate_values(frame: Frame, atoms: dict[str, int | None]) -> dict[str, float | None]:
    b=safe_atom(frame, atoms.get('base_atom'))
    h=safe_atom(frame, atoms.get('transfer_h'))
    f=safe_atom(frame, atoms.get('leaving_f'))
    r_bh=distance(b,h) if b and h else None
    r_hf=distance(h,f) if h and f else None
    q=(r_bh-r_hf) if r_bh is not None and r_hf is not None else None
    return {'r_BH_A': r_bh, 'r_HF_A': r_hf, 'q_BH_minus_HF_A': q}
