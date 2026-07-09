from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Any
import csv
import json

from hfauto_viz.structure.xyz import Atom, Frame
from hfauto_viz.structure.bond_changes import reaction_coordinate_atoms, infer_bhf_atoms, atom_position


@dataclass
class ModeVector:
    atom_index: int
    symbol: str
    x: float
    y: float
    z: float
    dx: float
    dy: float
    dz: float
    role: str = ""

    @property
    def norm(self) -> float:
        return sqrt(self.dx * self.dx + self.dy * self.dy + self.dz * self.dz)


def _sub(a: Atom, b: Atom) -> tuple[float, float, float]:
    return (a.x - b.x, a.y - b.y, a.z - b.z)


def _normalize(v: tuple[float, float, float]) -> tuple[float, float, float]:
    n = sqrt(sum(x * x for x in v))
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    return tuple(x / n for x in v)  # type: ignore[return-value]


def _safe(frame: Frame | None, idx: Any) -> Atom | None:
    try:
        i = int(idx)
    except Exception:
        return None
    if not frame or i < 0 or i >= len(frame):
        return None
    return frame[i]


def approximate_proton_transfer_mode(
    ts: Frame | None,
    reaction_data: dict[str, Any] | None,
    reactant: Frame | None = None,
    product: Frame | None = None,
    scale: float = 0.65,
) -> list[ModeVector]:
    """Return a chemically interpretable approximate imaginary-mode vector.

    Production runs should eventually read true normal modes from ORCA Hessian or
    frequency output. This fallback is deliberately labeled as approximate and is
    generated from the canonical B-H-F reaction coordinate so reviewers can see
    the expected proton-transfer direction.
    """
    frame = ts or reactant or product
    if not frame:
        return []
    atoms = reaction_coordinate_atoms(reaction_data) or infer_bhf_atoms(frame)
    b = _safe(frame, atoms.get("base_atom"))
    h = _safe(frame, atoms.get("transfer_h"))
    f = _safe(frame, atoms.get("leaving_f"))
    if not h:
        return []
    vectors: list[ModeVector] = []
    if b and f:
        # Proton motion roughly from F toward base atom in PT coordinate.
        vh = _normalize(_sub(b, f))
    elif b:
        vh = _normalize(_sub(b, h))
    elif f:
        vh = _normalize(_sub(h, f))
    else:
        vh = (1.0, 0.0, 0.0)
    # Add H vector strongest; base/F small opposite components for context.
    for idx, role, coeff in [
        (atoms.get("transfer_h"), "transfer_h", 1.0),
        (atoms.get("base_atom"), "base_atom", -0.22),
        (atoms.get("leaving_f"), "leaving_f", -0.18),
    ]:
        atom = _safe(frame, idx)
        if atom:
            vectors.append(ModeVector(int(idx), atom.symbol, atom.x, atom.y, atom.z, vh[0] * scale * coeff, vh[1] * scale * coeff, vh[2] * scale * coeff, role))
    # If reactant/product frames are available, add small displacement vectors for atoms that move most.
    if reactant and product and len(reactant) == len(product):
        existing = {v.atom_index for v in vectors}
        disps = []
        for i, (ra, pa) in enumerate(zip(reactant, product)):
            if i in existing:
                continue
            dx, dy, dz = pa.x - ra.x, pa.y - ra.y, pa.z - ra.z
            norm = sqrt(dx*dx + dy*dy + dz*dz)
            if norm > 0.15:
                disps.append((norm, i, dx, dy, dz))
        for norm, i, dx, dy, dz in sorted(disps, reverse=True)[:6]:
            atom = frame[i]
            unit = _normalize((dx, dy, dz))
            vectors.append(ModeVector(i, atom.symbol, atom.x, atom.y, atom.z, unit[0]*scale*0.25, unit[1]*scale*0.25, unit[2]*scale*0.25, "endpoint_displacement"))
    return vectors


def annotations_for_mode(vectors: list[ModeVector], color: str = "#8e44ad") -> list[dict[str, Any]]:
    ann: list[dict[str, Any]] = []
    for v in vectors:
        if v.norm < 1e-6:
            continue
        start = {"x": v.x, "y": v.y, "z": v.z}
        end = {"x": v.x + v.dx, "y": v.y + v.dy, "z": v.z + v.dz}
        ann.append({"type": "arrow", "start": start, "end": end, "color": color, "radius": 0.055, "radiusRatio": 1.8, "mid": 0.72})
        ann.append({"type": "label", "position": end, "text": f"ν‡ {v.role or v.atom_index}", "color": color})
    return ann


def mode_vectors_table(vectors: list[ModeVector]) -> list[dict[str, Any]]:
    return [{**v.__dict__, "norm": v.norm} for v in vectors]


def write_mode_vectors(vectors: list[ModeVector], out_prefix: str | Path) -> dict[str, Path]:
    prefix = Path(out_prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    rows = mode_vectors_table(vectors)
    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["atom_index", "symbol", "role", "x", "y", "z", "dx", "dy", "dz", "norm"]
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})
    return {"mode_vectors_json": json_path, "mode_vectors_csv": csv_path}
