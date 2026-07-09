from __future__ import annotations

"""Geometry helpers for proton-transfer TS and IRC validation.

These utilities intentionally operate on plain XYZ files and atom-order matched
endpoints. They do not try to infer chemistry from connectivity, which keeps
TS/IRC validation auditable and independent from any specific QM package.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.hf_builder import XYZ, read_xyz, write_xyz


@dataclass
class EndpointMatch:
    ok: bool
    rmsd_A: float | None
    orientation: str | None
    forward_matches: str | None
    backward_matches: str | None
    reason: str | None = None


def kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Return RMSD after optimal rigid alignment for same atom ordering."""
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape or len(a) == 0:
        return float("inf")
    ac = a - a.mean(axis=0)
    bc = b - b.mean(axis=0)
    h = ac.T @ bc
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    aligned = ac @ r
    return float(np.sqrt(np.mean(np.sum((aligned - bc) ** 2, axis=1))))


def same_atom_order(a: XYZ, b: XYZ) -> bool:
    return list(a.symbols) == list(b.symbols) and len(a.symbols) == len(b.symbols)


def interpolate_xyz(reactant_xyz: str | Path, product_xyz: str | Path, out_path: str | Path, fraction: float = 0.5, comment: str = "state=transition_state_guess") -> Path:
    r = read_xyz(reactant_xyz)
    p = read_xyz(product_xyz)
    if not same_atom_order(r, p):
        raise ValueError("Cannot interpolate endpoints with different atom order")
    frac = float(fraction)
    coords = (1.0 - frac) * r.coords + frac * p.coords
    return write_xyz(XYZ(symbols=list(r.symbols), coords=coords, comment=comment), out_path)


def proton_transfer_coordinate(xyz_path: str | Path, reaction_coordinate: dict[str, Any] | None) -> dict[str, float | None]:
    """Return q=r(B-H)-r(H-F) and component distances in Angstrom."""
    rc = reaction_coordinate or {}
    atoms = (rc.get("atoms") or {}) if isinstance(rc, dict) else {}
    base = atoms.get("base_atom")
    h = atoms.get("transfer_h")
    f = atoms.get("leaving_f")
    if base is None or h is None or f is None:
        return {"pt_coordinate_A": None, "B_H_distance_A": None, "H_F_distance_A": None}
    xyz = read_xyz(xyz_path)
    base_i, h_i, f_i = int(base), int(h), int(f)
    if max(base_i, h_i, f_i) >= len(xyz.coords):
        return {"pt_coordinate_A": None, "B_H_distance_A": None, "H_F_distance_A": None}
    r_bh = float(np.linalg.norm(xyz.coords[base_i] - xyz.coords[h_i]))
    r_hf = float(np.linalg.norm(xyz.coords[h_i] - xyz.coords[f_i]))
    return {"pt_coordinate_A": r_bh - r_hf, "B_H_distance_A": r_bh, "H_F_distance_A": r_hf}


def proton_transfer_mode_proxy_score(ts_xyz: str | Path, reaction_coordinate: dict[str, Any] | None, midpoint_width_A: float = 0.55) -> dict[str, Any]:
    """Estimate whether a TS geometry lies on the intended PT coordinate.

    A true mode-overlap requires normal-mode eigenvectors. In production ORCA
    runs this function is a fallback only. It grades the TS geometry by how close
    q=r(B-H)-r(H-F) is to zero. The output records the source explicitly so the
    score is not confused with an eigenvector overlap.
    """
    qdata = proton_transfer_coordinate(ts_xyz, reaction_coordinate)
    q = qdata.get("pt_coordinate_A")
    if q is None:
        return {**qdata, "mode_overlap_score": None, "mode_overlap_source": "missing_reaction_coordinate"}
    width = max(float(midpoint_width_A), 1.0e-6)
    score = max(0.0, min(1.0, 1.0 - abs(float(q)) / width))
    return {**qdata, "mode_overlap_score": float(score), "mode_overlap_source": "geometry_proxy_q_near_zero"}


def parse_mode_overlap_from_text(text: str) -> float | None:
    import re

    patterns = [
        r"MODE\s+OVERLAP\s+SCORE\s*[:=]\s*([0-9]*\.?[0-9]+)",
        r"mode_overlap_score\s*[:=]\s*([0-9]*\.?[0-9]+)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.I)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None


def _rmsd_pair(a: str | Path, b: str | Path) -> float | None:
    xa = read_xyz(a)
    xb = read_xyz(b)
    if not same_atom_order(xa, xb):
        return None
    return kabsch_rmsd(xa.coords, xb.coords)


def endpoint_pair_match(
    forward_endpoint: str | Path,
    backward_endpoint: str | Path,
    reactant_ref: str | Path,
    product_ref: str | Path,
    threshold_A: float = 0.75,
) -> EndpointMatch:
    """Match IRC endpoints to reactant/product references.

    ORCA's forward/backward direction is arbitrary. We therefore test both
    assignments and choose the smaller maximum RMSD.
    """
    try:
        fr = _rmsd_pair(forward_endpoint, reactant_ref)
        fp = _rmsd_pair(forward_endpoint, product_ref)
        br = _rmsd_pair(backward_endpoint, reactant_ref)
        bp = _rmsd_pair(backward_endpoint, product_ref)
    except Exception as exc:
        return EndpointMatch(False, None, None, None, None, reason=f"endpoint_parse_failed: {exc}")
    if None in {fr, fp, br, bp}:
        return EndpointMatch(False, None, None, None, None, reason="atom_order_mismatch")
    # assignment A: forward->product, backward->reactant
    a = max(float(fp), float(br))
    # assignment B: forward->reactant, backward->product
    b = max(float(fr), float(bp))
    if a <= b:
        ok = a <= float(threshold_A)
        return EndpointMatch(ok, a, "forward_product_backward_reactant", "product", "reactant", None if ok else "rmsd_above_threshold")
    ok = b <= float(threshold_A)
    return EndpointMatch(ok, b, "forward_reactant_backward_product", "reactant", "product", None if ok else "rmsd_above_threshold")



def endpoint_match(reference_xyz: str | Path, observed_xyz: str | Path, threshold_A: float = 0.75) -> dict[str, Any]:
    """Backward-compatible two-structure endpoint match used by dummy backend."""
    try:
        rmsd = _rmsd_pair(reference_xyz, observed_xyz)
    except Exception as exc:
        return {"endpoint_match": False, "endpoint_rmsd_A": None, "reason": str(exc)}
    ok = bool(rmsd is not None and float(rmsd) <= float(threshold_A))
    return {"endpoint_match": ok, "endpoint_rmsd_A": rmsd, "threshold_A": threshold_A}


def graph_match(reference_xyz: str | Path, observed_xyz: str | Path) -> dict[str, Any]:
    """Lightweight graph proxy: atom count and symbols must match in order."""
    try:
        ref = read_xyz(reference_xyz)
        obs = read_xyz(observed_xyz)
    except Exception as exc:
        return {"graph_match": False, "reason": str(exc)}
    return {"graph_match": same_atom_order(ref, obs), "same_atom_order": same_atom_order(ref, obs), "n_atoms_reference": len(ref.symbols), "n_atoms_observed": len(obs.symbols)}


def reaction_coordinate_value(xyz: XYZ | str | Path, atoms: dict[str, Any] | None) -> float | None:
    """Backward-compatible q=r(B-H)-r(H-F) for XYZ objects or files."""
    atoms = atoms or {}
    base = atoms.get("base_atom")
    h = atoms.get("transfer_h")
    f = atoms.get("leaving_f")
    if base is None or h is None or f is None:
        return None
    obj = xyz if isinstance(xyz, XYZ) else read_xyz(xyz)
    base_i, h_i, f_i = int(base), int(h), int(f)
    if max(base_i, h_i, f_i) >= len(obj.coords):
        return None
    return float(np.linalg.norm(obj.coords[base_i] - obj.coords[h_i]) - np.linalg.norm(obj.coords[h_i] - obj.coords[f_i]))


def write_irc_endpoints_from_expected(reactant_xyz: str | Path, product_xyz: str | Path, forward_out: str | Path, backward_out: str | Path) -> tuple[Path, Path]:
    """Write deterministic fallback IRC endpoints from expected product/reactant."""
    prod = read_xyz(product_xyz)
    reac = read_xyz(reactant_xyz)
    fwd = write_xyz(XYZ(list(prod.symbols), prod.coords.copy(), comment="state=irc_forward fallback=expected_product"), forward_out)
    bwd = write_xyz(XYZ(list(reac.symbols), reac.coords.copy(), comment="state=irc_backward fallback=expected_reactant"), backward_out)
    return fwd, bwd

# Phase 5 public convenience wrappers used by tests and reports.
from types import SimpleNamespace


def proton_transfer_geometry_metrics(xyz_path: str | Path, species_data: dict[str, Any]) -> dict[str, Any]:
    """Return PT coordinate metrics plus the geometry-proxy mode-overlap score."""
    rc = species_data.get("reaction_coordinate") or {}
    metrics = proton_transfer_coordinate(xyz_path, rc)
    proxy = proton_transfer_mode_proxy_score(xyz_path, rc)
    return {**metrics, **proxy}


def match_irc_endpoints(endpoint_a_xyz: str | Path, endpoint_b_xyz: str | Path, reactant_xyz: str | Path, product_xyz: str | Path, threshold_A: float = 0.75):
    """Match two unordered IRC endpoints against reactant/product references."""
    ep = endpoint_pair_match(endpoint_a_xyz, endpoint_b_xyz, reactant_xyz, product_xyz, threshold_A)
    return SimpleNamespace(
        ok=ep.ok,
        forward_ok=ep.ok,
        backward_ok=ep.ok,
        endpoint_graph_match=ep.ok,
        endpoint_assignment=ep.orientation or "unmatched",
        endpoint_rmsd_A=ep.rmsd_A,
        reason=ep.reason,
    )
