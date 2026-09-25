from __future__ import annotations

from math import acos, degrees
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.xyz import read_xyz


def _distance(coords: np.ndarray, i: int, j: int) -> float | None:
    if i < 0 or j < 0 or i >= len(coords) or j >= len(coords):
        return None
    return float(np.linalg.norm(coords[i] - coords[j]))


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float | None:
    """Angle ABC in degrees."""
    v1 = a - b
    v2 = c - b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-12 or n2 < 1e-12:
        return None
    cosang = float(np.dot(v1, v2) / (n1 * n2))
    cosang = max(-1.0, min(1.0, cosang))
    return float(degrees(acos(cosang)))


def _min_interatomic_distance(coords: np.ndarray) -> float | None:
    if len(coords) < 2:
        return None
    best = 999.0
    for i in range(len(coords)):
        delta = coords[i + 1 :] - coords[i]
        if len(delta):
            best = min(best, float(np.min(np.linalg.norm(delta, axis=1))))
    return best if best < 999.0 else None


def _component_hf_bonds(species_data: dict[str, Any], coords: np.ndarray) -> list[float]:
    bonds: list[float] = []
    for comp in species_data.get("components", []) or []:
        name = str(comp.get("name", ""))
        indices = comp.get("atom_indices", []) or []
        if name.startswith("HF_") and len(indices) >= 2:
            d = _distance(coords, int(indices[0]), int(indices[1]))
            if d is not None:
                bonds.append(d)
    return bonds


def geometry_qc_from_xyz(species_data: dict[str, Any], xyz_path: str | Path | None) -> dict[str, Any]:
    """Return chemistry-aware geometry QC for HF reaction artifacts.

    The thresholds are deliberately conservative. They are not meant to replace
    manual inspection or IRC validation; they catch common automation failures:
    atom collisions, HF drifting away, and unintended proton transfer during
    low-level preoptimization.
    """
    if not xyz_path:
        return {"geometry_sane": False, "geometry_qc_status": "missing_xyz_path"}
    try:
        xyz = read_xyz(xyz_path)
    except Exception as exc:
        return {"geometry_sane": False, "geometry_qc_status": "xyz_parse_failed", "reason": str(exc)}

    coords = xyz.coords
    state = species_data.get("state")
    expected_atoms = None
    if "components" in species_data:
        try:
            expected_atoms = max(max(c.get("atom_indices", [-1])) for c in species_data.get("components", []) or []) + 1
        except Exception:
            expected_atoms = None
    min_d = _min_interatomic_distance(coords)
    hf_bonds = _component_hf_bonds(species_data, coords)
    qc: dict[str, Any] = {
        "geometry_qc_status": "ok",
        "n_atoms_observed": len(xyz.symbols),
        "n_atoms_expected": expected_atoms,
        "atom_count_ok": expected_atoms is None or expected_atoms == len(xyz.symbols),
        "min_interatomic_distance_A": min_d,
        "collision_detected": bool(min_d is not None and min_d < 0.55),
        "hf_bond_distances_A": hf_bonds,
        "hf_bond_min_A": min(hf_bonds) if hf_bonds else None,
        "hf_bond_max_A": max(hf_bonds) if hf_bonds else None,
    }

    rc = species_data.get("reaction_coordinate", {}) or {}
    atoms = rc.get("atoms", {}) or {}
    base_atom = atoms.get("base_atom")
    transfer_h = atoms.get("transfer_h")
    leaving_f = atoms.get("leaving_f")
    if base_atom is not None and transfer_h is not None and leaving_f is not None:
        base_atom = int(base_atom)
        transfer_h = int(transfer_h)
        leaving_f = int(leaving_f)
        r_bh = _distance(coords, base_atom, transfer_h)
        r_hf = _distance(coords, transfer_h, leaving_f)
        angle_bhf = None
        if r_bh is not None and r_hf is not None:
            angle_bhf = _angle(coords[base_atom], coords[transfer_h], coords[leaving_f])
        qc.update(
            {
                "B_H_distance_A": r_bh,
                "H_F_distance_A": r_hf,
                "r_HF_A": r_hf,
                "B_H_F_angle_deg": angle_bhf,
            }
        )
        if state == "reactant_complex":
            qc["proton_transferred_unintentionally"] = bool(
                r_bh is not None and r_hf is not None and r_bh < 1.25 and r_hf > 1.15
            )
            qc["hf_dissociated"] = bool(
                (r_bh is not None and r_bh > 3.5)
                or (r_hf is not None and r_hf > 1.45 and not qc["proton_transferred_unintentionally"])
            )
        elif state == "ion_pair":
            qc["proton_transfer_confirmed"] = bool(
                r_bh is not None and r_hf is not None and r_bh < 1.30 and r_hf > 1.10
            )
            qc["ion_pair_dissociated"] = bool(r_bh is not None and r_bh > 1.80)
        elif state == "transition_state":
            canonical_q = None if r_bh is None or r_hf is None else float(r_hf - r_bh)
            qc["pt_coordinate_A"] = canonical_q
            qc["q_HF_minus_BH_A"] = canonical_q
            qc["q_BH_minus_HF_A"] = None if canonical_q is None else -canonical_q

    # General HF sanity for clusters and complexes. Ion-pair products may contain
    # FHF- motifs, so not all H-F distances are expected to be covalent there.
    if state in {"hf_cluster", "reactant_complex"} and hf_bonds:
        qc["hf_fragment_bonds_ok"] = all(0.65 <= d <= 1.35 for d in hf_bonds)
    elif hf_bonds:
        qc["hf_fragment_bonds_ok"] = all(0.55 <= d <= 1.80 for d in hf_bonds)
    else:
        qc["hf_fragment_bonds_ok"] = None

    qc["geometry_sane"] = bool(
        qc["atom_count_ok"]
        and not qc["collision_detected"]
        and qc.get("hf_fragment_bonds_ok", True) is not False
        and not qc.get("hf_dissociated", False)
    )
    return qc
