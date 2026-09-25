from __future__ import annotations

"""Reaction-path validation helpers for HF proton-transfer workflows.

These utilities deliberately avoid graph-toolkit dependencies so they can be
used by ORCA, pysisyphus, dummy, and offline tests.  They provide conservative,
machine-readable QC for two common questions:

1. Does a TS-like geometry sit between reactant complex and ion-pair product
   along q = r(H-F) - r(B-H)?
2. Do IRC endpoints match the intended reactant/product endpoint geometries?

The endpoint match is not a replacement for full chemical inspection, but it
catches the automation failures that most often invalidate high-throughput TS
jobs: wrong atom ordering, unrelated endpoint, or IRC path returning to the same
minimum on both sides.
"""

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.geometry import xyz_rmsd
from hfauto.chemistry.xyz import XYZ, read_xyz


def _distance(coords: np.ndarray, i: int, j: int) -> float | None:
    if i is None or j is None:
        return None
    if i < 0 or j < 0 or i >= len(coords) or j >= len(coords):
        return None
    return float(np.linalg.norm(coords[i] - coords[j]))


def _atoms(reaction_or_species_data: dict[str, Any]) -> dict[str, int]:
    rc = reaction_or_species_data.get("reaction_coordinate", {}) or {}
    atoms = rc.get("atoms", {}) or {}
    return {k: int(v) for k, v in atoms.items() if v is not None}


def proton_transfer_coordinate(xyz: XYZ, atoms: dict[str, int]) -> dict[str, float | None]:
    base = atoms.get("base_atom")
    h = atoms.get("transfer_h")
    f = atoms.get("leaving_f")
    r_bh = _distance(xyz.coords, base, h) if base is not None and h is not None else None
    r_hf = _distance(xyz.coords, h, f) if h is not None and f is not None else None
    q = None if r_bh is None or r_hf is None else float(r_hf - r_bh)
    return {
        "r_BH_A": r_bh,
        "r_HF_A": r_hf,
        "q_HF_minus_BH_A": q,
        "q_BH_minus_HF_A": None if q is None else -q,
    }


def estimate_ts_mode_overlap_from_geometry(ts_xyz_path: str | Path, reaction_data: dict[str, Any]) -> dict[str, Any]:
    """Return a transparent proxy for TS-mode relevance.

    ORCA normal-mode eigenvectors are not reliably available in all text outputs
    unless extra files are retained.  Until a full Hessian parser is added, this
    geometry proxy is used only as an automation screen and is labelled as such.
    It is high when the TS geometry lies near the proton-transfer midpoint.
    """
    try:
        xyz = read_xyz(ts_xyz_path)
    except Exception as exc:
        return {
            "mode_overlap_score": None,
            "mode_overlap_model": "geometry_proxy_failed",
            "mode_overlap_reason": str(exc),
        }
    atoms = _atoms(reaction_data)
    q = proton_transfer_coordinate(xyz, atoms)
    val = q.get("q_HF_minus_BH_A")
    if val is None:
        return {**q, "mode_overlap_score": None, "mode_overlap_model": "geometry_proxy_missing_reaction_atoms"}
    # |q| near zero resembles a centered proton-transfer TS.  Decay slowly so
    # asymmetric but chemically relevant TSs are not rejected too aggressively.
    score = max(0.0, min(1.0, 1.0 - abs(float(val)) / 1.5))
    return {**q, "mode_overlap_score": float(round(score, 4)), "mode_overlap_model": "geometry_proxy_q=rHF-rBH"}


def endpoint_match(reference_xyz: str | Path, endpoint_xyz: str | Path, reaction_data: dict[str, Any], expected_state: str) -> dict[str, Any]:
    """Compare an IRC endpoint with an intended endpoint.

    ``expected_state`` is either ``reactant_complex`` or ``ion_pair`` and adds a
    reaction-coordinate sanity check on top of atom-order and RMSD checks.
    """
    try:
        ref = read_xyz(reference_xyz)
        end = read_xyz(endpoint_xyz)
    except Exception as exc:
        return {"endpoint_match": False, "endpoint_match_reason": f"xyz_parse_failed: {exc}"}
    symbols_ok = ref.symbols == end.symbols
    rmsd = xyz_rmsd(ref, end) if symbols_ok else None
    atoms = _atoms(reaction_data)
    q = proton_transfer_coordinate(end, atoms)
    r_bh = q.get("r_BH_A")
    r_hf = q.get("r_HF_A")
    coord_ok = None
    if expected_state == "reactant_complex" and r_bh is not None and r_hf is not None:
        coord_ok = bool(r_bh > 1.25 and r_hf < 1.35)
    elif expected_state == "ion_pair" and r_bh is not None and r_hf is not None:
        coord_ok = bool(r_bh < 1.35 and r_hf > 1.05)
    rmsd_ok = bool(rmsd is not None and rmsd <= 1.25)
    # For floppy HF clusters the RMSD threshold can be too strict.  A correct
    # reaction-coordinate endpoint with matching atoms is enough for an automated
    # pass, while exact RMSD is still reported for review.
    ok = bool(symbols_ok and (rmsd_ok or coord_ok is True))
    return {
        "endpoint_match": ok,
        "symbols_match": symbols_ok,
        "endpoint_rmsd_A": rmsd,
        "endpoint_rmsd_ok": rmsd_ok,
        "endpoint_coordinate_ok": coord_ok,
        **q,
    }


def split_irc_endpoint_blocks(blocks: list[XYZ]) -> tuple[XYZ | None, XYZ | None]:
    """Best-effort extraction of backward/forward endpoints from ORCA blocks."""
    if not blocks:
        return None, None
    if len(blocks) == 1:
        return blocks[0], blocks[0]
    return blocks[0], blocks[-1]
