from __future__ import annotations

from dataclasses import asdict
from math import sqrt
from pathlib import Path
from typing import Any
import csv
import json

from hfauto_viz.structure.xyz import Atom, Frame


def distance(a: Atom, b: Atom) -> float:
    return sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def atom_position(atom: Atom) -> dict[str, float]:
    return {"x": atom.x, "y": atom.y, "z": atom.z}


def midpoint(a: Atom, b: Atom) -> dict[str, float]:
    return {"x": (a.x + b.x) / 2.0, "y": (a.y + b.y) / 2.0, "z": (a.z + b.z) / 2.0}


def _safe_atom(frame: Frame | None, idx: Any) -> Atom | None:
    try:
        i = int(idx)
    except Exception:
        return None
    if frame is None or i < 0 or i >= len(frame):
        return None
    return frame[i]


def reaction_coordinate_atoms(reaction_data: dict[str, Any] | None) -> dict[str, int]:
    """Extract B/H/F atom indexes from a ReactionRecord-like dictionary.

    The canonical HF proton-transfer coordinate used throughout hfauto is
    q = r(B-H) - r(H-F), with keys base_atom, transfer_h and leaving_f.
    """
    if not reaction_data:
        return {}
    coord = reaction_data.get("reaction_coordinate") or {}
    atoms = coord.get("atoms", {}) if isinstance(coord, dict) else {}
    out: dict[str, int] = {}
    aliases = {
        "base_atom": ["base_atom", "B", "donor", "acceptor_atom"],
        "transfer_h": ["transfer_h", "H", "proton", "moving_h"],
        "leaving_f": ["leaving_f", "F", "fluoride", "hf_f"],
    }
    for key, names in aliases.items():
        for name in names:
            if atoms.get(name) is not None:
                try:
                    out[key] = int(atoms[name])
                    break
                except Exception:
                    pass
    return out


def infer_bhf_atoms(frame: Frame | None) -> dict[str, int]:
    """Conservative fallback for B/H/F annotation when ReactionRecord indexes are absent."""
    if not frame:
        return {}
    f_indices = [i for i, a in enumerate(frame) if a.symbol.upper() == "F"]
    h_indices = [i for i, a in enumerate(frame) if a.symbol.upper() == "H"]
    base_indices = [i for i, a in enumerate(frame) if a.symbol.upper() in {"N", "O", "S", "P"}]
    if not f_indices or not h_indices:
        return {}
    # Find the closest H-F pair, then closest base atom to that H.
    best_hf = None
    for hi in h_indices:
        for fi in f_indices:
            d = distance(frame[hi], frame[fi])
            if best_hf is None or d < best_hf[0]:
                best_hf = (d, hi, fi)
    if best_hf is None:
        return {}
    _, hi, fi = best_hf
    base = None
    for bi in base_indices:
        d = distance(frame[bi], frame[hi])
        if base is None or d < base[0]:
            base = (d, bi)
    out = {"transfer_h": hi, "leaving_f": fi}
    if base:
        out["base_atom"] = base[1]
    return out


def bond_change_table(
    reactant: Frame | None,
    ts: Frame | None,
    product: Frame | None,
    reaction_data: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    atoms = reaction_coordinate_atoms(reaction_data)
    if not atoms:
        atoms = infer_bhf_atoms(ts or reactant or product)
    rows: list[dict[str, Any]] = []
    pairs = [
        ("forming_B_H", atoms.get("base_atom"), atoms.get("transfer_h"), "forming"),
        ("breaking_H_F", atoms.get("transfer_h"), atoms.get("leaving_f"), "breaking"),
    ]
    for label, i, j, change_type in pairs:
        if i is None or j is None:
            continue
        row: dict[str, Any] = {"label": label, "change_type": change_type, "atom_i": i, "atom_j": j}
        for name, frame in [("reactant", reactant), ("ts", ts), ("product", product)]:
            ai, aj = _safe_atom(frame, i), _safe_atom(frame, j)
            row[f"r_{name}_A"] = distance(ai, aj) if ai and aj else None
        if row.get("r_reactant_A") is not None and row.get("r_product_A") is not None:
            row["delta_product_minus_reactant_A"] = row["r_product_A"] - row["r_reactant_A"]
        rows.append(row)
    return rows


def annotation_shapes_for_frame(
    frame: Frame | None,
    reaction_data: dict[str, Any] | None = None,
    include_labels: bool = True,
    include_bond_lines: bool = True,
) -> list[dict[str, Any]]:
    if not frame:
        return []
    atoms = reaction_coordinate_atoms(reaction_data) or infer_bhf_atoms(frame)
    annotations: list[dict[str, Any]] = []
    role_style = {
        "base_atom": {"color": "#1f77b4", "label": "B site"},
        "transfer_h": {"color": "#d62728", "label": "H transfer"},
        "leaving_f": {"color": "#2ca02c", "label": "F"},
    }
    for role, idx in atoms.items():
        atom = _safe_atom(frame, idx)
        if not atom:
            continue
        annotations.append({
            "type": "sphere",
            "center": atom_position(atom),
            "radius": 0.33 if role == "transfer_h" else 0.28,
            "color": role_style.get(role, {}).get("color", "#444"),
            "alpha": 0.55,
        })
        if include_labels:
            annotations.append({
                "type": "label",
                "position": {"x": atom.x + 0.18, "y": atom.y + 0.18, "z": atom.z + 0.18},
                "text": f"{role}:{idx}",
                "color": role_style.get(role, {}).get("color", "#111"),
            })
    if include_bond_lines:
        b, h, f = atoms.get("base_atom"), atoms.get("transfer_h"), atoms.get("leaving_f")
        for label, i, j, color, dashed in [
            ("forming B-H", b, h, "#d62728", True),
            ("breaking H-F", h, f, "#2ca02c", False),
        ]:
            ai, aj = _safe_atom(frame, i), _safe_atom(frame, j)
            if ai and aj:
                annotations.append({
                    "type": "line",
                    "start": atom_position(ai),
                    "end": atom_position(aj),
                    "color": color,
                    "dashed": dashed,
                })
                if include_labels:
                    annotations.append({
                        "type": "label",
                        "position": midpoint(ai, aj),
                        "text": f"{label} {distance(ai, aj):.2f} Å",
                        "color": color,
                    })
    return annotations


def write_bond_change_outputs(
    rows: list[dict[str, Any]],
    out_prefix: str | Path,
) -> dict[str, Path]:
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    if rows:
        fields = sorted({k for row in rows for k in row})
        with csv_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    else:
        csv_path.write_text("label,change_type,atom_i,atom_j\n", encoding="utf-8")
    return {"bond_change_json": json_path, "bond_change_csv": csv_path}
