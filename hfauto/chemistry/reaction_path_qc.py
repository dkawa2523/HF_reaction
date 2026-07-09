from __future__ import annotations

"""Reaction-path QC helpers for proton-transfer TS/IRC automation.

These helpers are deliberately small and dependency-light.  They do not replace
manual inspection of important candidates; they provide machine-readable gates
that prevent proxy/fallback geometries from being silently treated as validated
activation barriers.
"""

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.hf_builder import XYZ, read_xyz, write_xyz


def _kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float | None:
    """Return atom-order RMSD after optimal rigid alignment."""
    if len(a) != len(b) or len(a) == 0:
        return None
    aa = a - a.mean(axis=0)
    bb = b - b.mean(axis=0)
    h = aa.T @ bb
    try:
        u, _s, vt = np.linalg.svd(h)
    except np.linalg.LinAlgError:
        return float(np.sqrt(np.mean(np.sum((aa - bb) ** 2, axis=1))))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    aa_rot = aa @ r
    return float(np.sqrt(np.mean(np.sum((aa_rot - bb) ** 2, axis=1))))


def xyz_rmsd(path_a: str | Path, path_b: str | Path) -> float | None:
    try:
        a = read_xyz(path_a)
        b = read_xyz(path_b)
    except Exception:
        return None
    if a.symbols != b.symbols:
        return None
    return _kabsch_rmsd(a.coords, b.coords)


def reaction_coordinate_value(species_data: dict[str, Any], xyz_path: str | Path | None) -> float | None:
    """Return q = r(B-H) - r(H-F) for the configured proton-transfer coordinate."""
    if not xyz_path:
        return None
    qc = geometry_qc_from_xyz(species_data, xyz_path)
    r_bh = qc.get("B_H_distance_A")
    r_hf = qc.get("H_F_distance_A")
    if r_bh is None or r_hf is None:
        return None
    return float(r_bh) - float(r_hf)


def make_midpoint_ts_xyz(reactant_xyz: str | Path, product_xyz: str | Path, out_xyz: str | Path) -> Path:
    """Write a midpoint TS guess from atom-matched reactant/product endpoints."""
    rc = read_xyz(reactant_xyz)
    ip = read_xyz(product_xyz)
    if rc.symbols != ip.symbols:
        raise ValueError("reactant/product atom symbols are not identical")
    coords = 0.5 * (rc.coords + ip.coords)
    return write_xyz(XYZ(list(rc.symbols), coords, comment="state=transition_state generated_by=hfauto midpoint_ts"), out_xyz)


def estimate_pt_mode_overlap(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
    n_imag: int | None = None,
    imag_freq_cm1: float | None = None,
    output_text: str | None = None,
) -> tuple[float, str]:
    """Estimate whether the imaginary mode is likely the proton-transfer mode.

    ORCA normal-mode displacement parsing differs by print level and version.
    Phase 5 therefore provides a robust geometry-based gate that can be replaced
    later by a true displacement projection parser.  It is intentionally labelled
    as a heuristic in the output.
    """
    # Allow tests or high-print ORCA parsers to inject a literal score.
    if output_text:
        import re

        m = re.search(r"HFAUTO_MODE_OVERLAP\s*[:=]\s*(0(?:\.\d+)?|1(?:\.0+)?)", output_text, flags=re.I)
        if m:
            return float(m.group(1)), "explicit_hfauto_marker"

    if n_imag != 1 or imag_freq_cm1 is None or imag_freq_cm1 >= -100.0:
        return 0.0, "frequency_gate_failed"
    q = reaction_coordinate_value(species_data, xyz_path)
    if q is None:
        return 0.72, "frequency_only_no_coordinate"
    # A proton-transfer TS should be near the crossing region.  q exactly zero is
    # ideal; q around +/-0.6 A is still plausible for a loose acid-base TS.
    score = 0.90 - min(abs(float(q)), 1.0) * 0.40
    return max(0.0, min(0.95, score)), "geometry_coordinate_heuristic"


def endpoint_match_qc(
    observed_xyz: str | Path | None,
    expected_xyz: str | Path | None,
    species_data: dict[str, Any],
    rmsd_threshold_A: float = 0.75,
) -> dict[str, Any]:
    if not observed_xyz or not expected_xyz:
        return {"endpoint_match": False, "endpoint_match_status": "missing_endpoint_xyz", "endpoint_rmsd_A": None}
    rmsd = xyz_rmsd(observed_xyz, expected_xyz)
    if rmsd is None:
        return {"endpoint_match": False, "endpoint_match_status": "atom_order_or_xyz_mismatch", "endpoint_rmsd_A": None}
    q_obs = reaction_coordinate_value(species_data, observed_xyz)
    q_exp = reaction_coordinate_value(species_data, expected_xyz)
    return {
        "endpoint_match": bool(rmsd <= rmsd_threshold_A),
        "endpoint_match_status": "ok" if rmsd <= rmsd_threshold_A else "rmsd_above_threshold",
        "endpoint_rmsd_A": rmsd,
        "endpoint_q_observed_A": q_obs,
        "endpoint_q_expected_A": q_exp,
        "rmsd_threshold_A": rmsd_threshold_A,
    }


def endpoint_pair_match_qc(
    forward_xyz: str | Path | None,
    backward_xyz: str | Path | None,
    reactant_xyz: str | Path | None,
    product_xyz: str | Path | None,
    species_data: dict[str, Any],
    rmsd_threshold_A: float = 0.75,
) -> dict[str, Any]:
    """Match two IRC endpoint files to reactant/product, accepting either direction."""
    f_to_p = endpoint_match_qc(forward_xyz, product_xyz, species_data, rmsd_threshold_A)
    b_to_r = endpoint_match_qc(backward_xyz, reactant_xyz, species_data, rmsd_threshold_A)
    f_to_r = endpoint_match_qc(forward_xyz, reactant_xyz, species_data, rmsd_threshold_A)
    b_to_p = endpoint_match_qc(backward_xyz, product_xyz, species_data, rmsd_threshold_A)

    direct_ok = bool(f_to_p["endpoint_match"] and b_to_r["endpoint_match"])
    swapped_ok = bool(f_to_r["endpoint_match"] and b_to_p["endpoint_match"])
    if direct_ok or not swapped_ok:
        orientation = "forward_product_backward_reactant"
        chosen_f = f_to_p
        chosen_b = b_to_r
    else:
        orientation = "forward_reactant_backward_product"
        chosen_f = f_to_r
        chosen_b = b_to_p
    ok = bool(direct_ok or swapped_ok)
    return {
        "irc_validated": ok,
        "forward_ok": bool(chosen_f["endpoint_match"]),
        "backward_ok": bool(chosen_b["endpoint_match"]),
        "reactant_endpoint_match": ok,
        "product_endpoint_match": ok,
        "endpoint_orientation": orientation if ok else "unmatched",
        "endpoint_rmsd_A": max(
            x for x in [chosen_f.get("endpoint_rmsd_A"), chosen_b.get("endpoint_rmsd_A")] if x is not None
        ) if chosen_f.get("endpoint_rmsd_A") is not None and chosen_b.get("endpoint_rmsd_A") is not None else None,
        "forward_endpoint_qc": chosen_f,
        "backward_endpoint_qc": chosen_b,
    }
