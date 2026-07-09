from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.hf_builder import read_xyz
from hfauto.core.constants import DEFAULT_HF_BOND_A, DEFAULT_HF_STRETCH_CM1


def distance(a, b) -> float:
    return float(np.linalg.norm(np.array(a) - np.array(b)))


def angle_deg(a, b, c) -> float:
    ba = np.array(a) - np.array(b)
    bc = np.array(c) - np.array(b)
    nba = np.linalg.norm(ba)
    nbc = np.linalg.norm(bc)
    if nba < 1e-12 or nbc < 1e-12:
        return float("nan")
    cosang = float(np.dot(ba, bc) / (nba * nbc))
    cosang = max(-1.0, min(1.0, cosang))
    return float(np.degrees(np.arccos(cosang)))


def hf_descriptors_from_species(species_data: dict[str, Any], calc_data: dict[str, Any] | None = None) -> dict[str, Any]:
    xyz_path = species_data.get("xyz_path")
    components = species_data.get("components", [])
    if not xyz_path or not Path(xyz_path).exists():
        return {}
    xyz = read_xyz(xyz_path)
    hf_components = [c for c in components if str(c.get("name", "")).startswith("HF")]
    if not hf_components:
        return {}
    h_idx, f_idx = hf_components[0]["atom_indices"]
    if h_idx >= len(xyz.coords) or f_idx >= len(xyz.coords):
        return {}
    r_hf = distance(xyz.coords[h_idx], xyz.coords[f_idx])
    atoms = (species_data.get("reaction_coordinate") or {}).get("atoms") or {}
    geom_qc = species_data.get("geometry_qc") or {}
    base_atom = species_data.get("base_atom")
    if base_atom is None:
        base_atom = atoms.get("base_atom")
    if base_atom is None:
        base_atom = geom_qc.get("base_atom")
    b_h = None
    b_h_f = None
    if base_atom is not None and int(base_atom) < len(xyz.coords):
        b_h = distance(xyz.coords[int(base_atom)], xyz.coords[h_idx])
        b_h_f = angle_deg(xyz.coords[int(base_atom)], xyz.coords[h_idx], xyz.coords[f_idx])
    hf_stretch = None
    if calc_data:
        hf_stretch = calc_data.get("hf_stretch_cm1") or calc_data.get("frequencies", {}).get("hf_stretch_cm1")
    return {
        "r_HF_A": round(float(r_hf), 6),
        "delta_r_HF_A": round(float(r_hf - DEFAULT_HF_BOND_A), 6),
        "B_H_distance_A": None if b_h is None else round(float(b_h), 6),
        "B_H_F_angle_deg": None if b_h_f is None else round(float(b_h_f), 3),
        "nu_HF_cm1": hf_stretch,
        "delta_nu_HF_cm1": None if hf_stretch is None else round(float(hf_stretch) - DEFAULT_HF_STRETCH_CM1, 3),
        "HF_fragment_charge": None,
        "dipole_D": None,
    }
