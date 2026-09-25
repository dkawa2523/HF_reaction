"""Canonical geometry checks for B...H-F proton-transfer paths.

The whole project uses one signed coordinate: ``q = r(H-F) - r(B-H)``.
Consequently a neutral B...H-F complex has negative ``q`` and an ion pair
has positive ``q``.  Older manifests may carry the opposite definition as
metadata; atom indices remain compatible and are normalised here.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.xyz import XYZ, read_xyz

NEUTRAL_MAX_Q_A = -0.15
ION_PAIR_MIN_Q_A = 0.15
NEUTRAL_BASE_H_MAX_A = 2.50
ION_PAIR_PROTON_F_MAX_A = 3.00
REACTION_COORDINATE_DEFINITION = "r(H-F)-r(B-H)"


def proton_transfer_metrics_from_xyz(
    xyz: XYZ,
    base_atom: int,
    transfer_h: int,
    leaving_f: int,
) -> dict[str, float | str]:
    """Return canonical proton-transfer distances for an in-memory geometry."""
    indices = [int(base_atom), int(transfer_h), int(leaving_f)]
    if min(indices) < 0 or max(indices) >= len(xyz.coords):
        raise IndexError(f"Reaction-coordinate atom index outside geometry: {indices}")
    if xyz.symbols[indices[1]].upper() != "H" or xyz.symbols[indices[2]].upper() != "F":
        raise ValueError(
            "Reaction-coordinate atom types are inconsistent: "
            f"transfer_h={xyz.symbols[indices[1]]}, leaving_f={xyz.symbols[indices[2]]}"
        )
    r_bh = float(np.linalg.norm(xyz.coords[indices[0]] - xyz.coords[indices[1]]))
    r_hf = float(np.linalg.norm(xyz.coords[indices[1]] - xyz.coords[indices[2]]))
    q = r_hf - r_bh
    return {
        "r_BH_A": r_bh,
        "r_HF_A": r_hf,
        "q_HF_minus_BH_A": q,
        "endpoint_class": classify_proton_location(q),
    }


def proton_transfer_metrics(
    xyz_path: str | Path,
    base_atom: int,
    transfer_h: int,
    leaving_f: int,
) -> dict[str, float | str]:
    """Return distances and q = r(H-F) - r(B-H), in Angstrom."""
    return proton_transfer_metrics_from_xyz(
        read_xyz(xyz_path), base_atom, transfer_h, leaving_f
    )


def classify_proton_location(
    q_HF_minus_BH_A: float,
    neutral_max_q_A: float = NEUTRAL_MAX_Q_A,
    ion_pair_min_q_A: float = ION_PAIR_MIN_Q_A,
) -> str:
    """Classify q with a conservative, inclusive shared-proton band.

    ``q < -0.15 A`` means H is measurably closer to F, while ``q > +0.15 A``
    means it is measurably closer to the base.  Both exact boundaries remain
    shared/unresolved so a cutoff equality cannot promote an endpoint claim.
    """
    q = float(q_HF_minus_BH_A)
    if q < float(neutral_max_q_A):
        return "neutral_complex"
    if q > float(ion_pair_min_q_A):
        return "ion_pair"
    return "shared_proton"


def reaction_coordinate_atoms(reaction: Any) -> dict[str, int]:
    data = getattr(reaction, "data", reaction) or {}
    coordinate = data.get("reaction_coordinate", {}) or {}
    atoms = coordinate.get("atoms", {}) or {}
    required = ("base_atom", "transfer_h", "leaving_f")
    if not all(key in atoms for key in required):
        raise KeyError("reaction_coordinate.atoms must define base_atom, transfer_h, leaving_f")
    return {key: int(atoms[key]) for key in required}


def canonical_reaction_coordinate(reaction: Any) -> dict[str, Any]:
    """Return coordinate metadata with a stable definition and legacy provenance."""
    data = getattr(reaction, "data", reaction) or {}
    coordinate = data.get("reaction_coordinate", {}) or {}
    original = coordinate.get("definition")
    result = {
        "type": "distance_difference",
        "definition": REACTION_COORDINATE_DEFINITION,
        "atoms": reaction_coordinate_atoms(reaction),
        "index_base": 0,
    }
    if coordinate.get("hf_relay"):
        result["hf_relay"] = coordinate["hf_relay"]
    if original and str(original).replace(" ", "") != REACTION_COORDINATE_DEFINITION:
        result["source_definition"] = str(original)
        result["definition_normalized"] = True
    return result


def infer_spectator_hf_pairs(
    xyz_path: str | Path,
    atoms: dict[str, int],
    max_initial_distance_A: float = 1.35,
) -> list[dict[str, int | float]]:
    """Identify intact non-transferring HF pairs from the starting atom order.

    This is intentionally conservative and dependency-free.  Each fluorine is
    paired with its nearest unused hydrogen only when their initial distance is
    chemically plausible.  The resulting fixed atom pairs can then be followed
    throughout a relaxed scan to detect HF dissociation or proton identity swaps.
    """
    xyz = read_xyz(xyz_path)
    transfer_h = int(atoms["transfer_h"])
    leaving_f = int(atoms["leaving_f"])
    hydrogens = [i for i, symbol in enumerate(xyz.symbols) if symbol.upper() == "H" and i != transfer_h]
    fluorines = [i for i, symbol in enumerate(xyz.symbols) if symbol.upper() == "F" and i != leaving_f]
    unused_h = set(hydrogens)
    pairs: list[dict[str, int | float]] = []
    for f_index in fluorines:
        if not unused_h:
            break
        h_index = min(
            unused_h,
            key=lambda index: float(np.linalg.norm(xyz.coords[index] - xyz.coords[f_index])),
        )
        distance = float(np.linalg.norm(xyz.coords[h_index] - xyz.coords[f_index]))
        if distance <= float(max_initial_distance_A):
            pairs.append(
                {"h_atom": int(h_index), "f_atom": int(f_index), "initial_distance_A": distance}
            )
            unused_h.remove(h_index)
    return pairs


def spectator_hf_distances(
    xyz: XYZ,
    pairs: list[dict[str, int | float]],
) -> list[dict[str, int | float]]:
    """Measure the fixed spectator HF pairs in one scan geometry."""
    measured: list[dict[str, int | float]] = []
    for pair in pairs:
        h_index = int(pair["h_atom"])
        f_index = int(pair["f_atom"])
        if min(h_index, f_index) < 0 or max(h_index, f_index) >= len(xyz.coords):
            raise IndexError(f"Spectator HF atom index outside geometry: {(h_index, f_index)}")
        measured.append(
            {
                **pair,
                "distance_A": float(np.linalg.norm(xyz.coords[h_index] - xyz.coords[f_index])),
            }
        )
    return measured


def hf_chain_metrics(
    xyz: XYZ,
    atoms: dict[str, int],
    spectator_pairs: list[dict[str, int | float]],
) -> list[dict[str, int | float]]:
    """Track the F1...H2-F2 relay coordinate for every spectator HF."""
    base = int(atoms["base_atom"])
    leaving_f = int(atoms["leaving_f"])
    rows: list[dict[str, int | float]] = []
    for pair in spectator_pairs:
        h2 = int(pair["h_atom"])
        f2 = int(pair["f_atom"])
        r_f1_h2 = float(np.linalg.norm(xyz.coords[leaving_f] - xyz.coords[h2]))
        r_h2_f2 = float(np.linalg.norm(xyz.coords[h2] - xyz.coords[f2]))
        r_b_h2 = float(np.linalg.norm(xyz.coords[base] - xyz.coords[h2]))
        rows.append(
            {
                "h_atom": h2,
                "f_atom": f2,
                "r_F1_H2_A": r_f1_h2,
                "r_H2_F2_A": r_h2_f2,
                "r_B_H2_A": r_b_h2,
                "q2_H2F2_minus_F1H2_A": r_h2_f2 - r_f1_h2,
            }
        )
    return rows


def hf_endpoint_metrics_from_xyz(
    xyz: XYZ,
    atoms: dict[str, int],
    spectator_pairs: list[dict[str, int | float]],
    base_h_max_A: float = 1.30,
    separated_h_f_min_A: float = 1.10,
    intact_h_f_max_A: float = 1.15,
    neutral_base_h_max_A: float = NEUTRAL_BASE_H_MAX_A,
    ion_pair_proton_f_max_A: float = ION_PAIR_PROTON_F_MAX_A,
    counterion_f_h_max_A: float = 2.30,
) -> dict[str, Any]:
    """Classify an HFn endpoint without assuming which proton moved.

    Atom labels remain useful for following a scan, but an unconstrained endpoint
    optimization may exchange equivalent HF protons.  This classifier therefore
    considers every mapped HF-chain hydrogen and evaluates q for the H nearest
    the base.  The same inclusive ``[-0.15, +0.15] A`` shared band used by mapped
    scans is authoritative; absolute bond/contact cutoffs only reject implausible
    or disconnected geometries.  A minimum-cost H/F assignment recognizes an
    intact neutral HF cluster.
    """

    base = int(atoms["base_atom"])
    mapped_transfer_h = int(atoms["transfer_h"])
    mapped_leaving_f = int(atoms["leaving_f"])
    chain_h = [mapped_transfer_h, *(int(pair["h_atom"]) for pair in spectator_pairs)]
    chain_f = [mapped_leaving_f, *(int(pair["f_atom"]) for pair in spectator_pairs)]
    chain_h = list(dict.fromkeys(chain_h))
    chain_f = list(dict.fromkeys(chain_f))
    indices = [base, *chain_h, *chain_f]
    if min(indices) < 0 or max(indices) >= len(xyz.coords):
        raise IndexError(f"HF endpoint atom index outside geometry: {indices}")
    if len(chain_h) != len(chain_f):
        raise ValueError("HF endpoint classification requires equal H/F chain counts")
    if any(xyz.symbols[index].upper() != "H" for index in chain_h):
        raise ValueError("Mapped HF-chain hydrogen has a non-H element label")
    if any(xyz.symbols[index].upper() != "F" for index in chain_f):
        raise ValueError("Mapped HF-chain fluorine has a non-F element label")

    b_h_distances = {
        h_index: float(np.linalg.norm(xyz.coords[base] - xyz.coords[h_index]))
        for h_index in chain_h
    }
    proton_candidates: list[dict[str, int | float]] = []
    for h_index in chain_h:
        nearest_f = min(
            chain_f,
            key=lambda index: float(
                np.linalg.norm(xyz.coords[h_index] - xyz.coords[index])
            ),
        )
        r_hf = float(np.linalg.norm(xyz.coords[h_index] - xyz.coords[nearest_f]))
        proton_candidates.append(
            {
                "h_atom": int(h_index),
                "nearest_f_atom": int(nearest_f),
                "r_BH_A": b_h_distances[h_index],
                "r_HF_A": r_hf,
                "q_HF_minus_BH_A": r_hf - b_h_distances[h_index],
            }
        )
    # The H nearest the base is the protonation candidate, independent of its
    # original atom label.  Its nearest-F distance then defines the same q used
    # for mapped scans, without allowing a remote dissociated H to dominate.
    active_proton = min(
        proton_candidates,
        key=lambda row: (float(row["r_BH_A"]), int(row["h_atom"])),
    )
    protonated_h = int(active_proton["h_atom"])
    nearest_f_to_proton = int(active_proton["nearest_f_atom"])
    proton_f_distance = float(active_proton["r_HF_A"])
    identity_invariant_q = float(active_proton["q_HF_minus_BH_A"])
    proton_location_class = classify_proton_location(identity_invariant_q)

    best_assignment = min(
        permutations(chain_f),
        key=lambda perm: sum(
            float(np.linalg.norm(xyz.coords[h_index] - xyz.coords[f_index]))
            for h_index, f_index in zip(chain_h, perm)
        ),
    )
    matched_hf = [
        {
            "h_atom": int(h_index),
            "f_atom": int(f_index),
            "distance_A": float(
                np.linalg.norm(xyz.coords[h_index] - xyz.coords[f_index])
            ),
        }
        for h_index, f_index in zip(chain_h, best_assignment)
    ]
    protonated_contact = bool(
        b_h_distances[protonated_h] <= float(base_h_max_A)
        and proton_f_distance >= float(separated_h_f_min_A)
    )
    protonated_base = bool(
        proton_location_class == "ion_pair" and protonated_contact
    )
    ion_pair_cluster_contact = bool(
        protonated_base and proton_f_distance <= float(ion_pair_proton_f_max_A)
    )
    neutral_base_contact = bool(
        b_h_distances[protonated_h] <= float(neutral_base_h_max_A)
    )
    neutral_hf_cluster = bool(
        neutral_base_contact
        and all(row["distance_A"] <= float(intact_h_f_max_A) for row in matched_hf)
    )
    counterion_h = [index for index in chain_h if index != protonated_h]
    counterion_f_h_distances = [
        {
            "f_atom": int(f_index),
            "nearest_h_atom": (
                int(
                    min(
                        counterion_h,
                        key=lambda h_index: float(
                            np.linalg.norm(xyz.coords[f_index] - xyz.coords[h_index])
                        ),
                    )
                )
                if counterion_h
                else None
            ),
            "distance_A": (
                min(
                    float(np.linalg.norm(xyz.coords[f_index] - xyz.coords[h_index]))
                    for h_index in counterion_h
                )
                if counterion_h
                else None
            ),
        }
        for f_index in chain_f
    ]
    counterion_connected = bool(
        (len(chain_f) == 1 and not counterion_h)
        or (
            counterion_h
            and all(
                float(row["distance_A"]) <= float(counterion_f_h_max_A)
                for row in counterion_f_h_distances
            )
        )
    )
    shared_proton_geometry = bool(
        proton_location_class == "shared_proton"
        and min(b_h_distances[protonated_h], proton_f_distance)
        <= float(base_h_max_A)
        and counterion_connected
    )
    neutral_cluster_connected = bool(neutral_hf_cluster and counterion_connected)
    endpoint_class = (
        "ion_pair"
        if proton_location_class == "ion_pair"
        and protonated_base
        and ion_pair_cluster_contact
        and counterion_connected
        else "ion_pair_dissociated"
        if proton_location_class == "ion_pair" and protonated_base
        else "neutral_complex"
        if proton_location_class == "neutral_complex" and neutral_cluster_connected
        else "shared_proton"
        if shared_proton_geometry
        else "shared_or_rearranged"
    )
    endpoint_geometry_sane = endpoint_class in {
        "neutral_complex",
        "shared_proton",
        "ion_pair",
    }
    return {
        "endpoint_class": endpoint_class,
        "classification_scope": "permutation_invariant_hfn_geometry",
        "proton_location_class": proton_location_class,
        "q_HF_minus_BH_identity_invariant_A": identity_invariant_q,
        "proton_candidates": proton_candidates,
        "protonated_base_geometry": protonated_base,
        "protonated_contact_geometry": protonated_contact,
        "ion_pair_cluster_contact_geometry": ion_pair_cluster_contact,
        "neutral_base_contact_geometry": neutral_base_contact,
        "neutral_hf_cluster_geometry": neutral_hf_cluster,
        "neutral_cluster_connected_geometry": neutral_cluster_connected,
        "shared_proton_geometry": shared_proton_geometry,
        "counterion_connected": counterion_connected,
        "counterion_f_h_distances": counterion_f_h_distances,
        "mapped_transfer_h": mapped_transfer_h,
        "protonated_h_atom": int(protonated_h),
        "proton_identity_swapped": bool(protonated_h != mapped_transfer_h),
        "r_BH_min_A": b_h_distances[protonated_h],
        "r_proton_nearest_F_A": proton_f_distance,
        "nearest_f_to_proton": int(nearest_f_to_proton),
        "chain_h_atoms": chain_h,
        "chain_f_atoms": chain_f,
        "matched_hf_pairs": matched_hf,
        "endpoint_geometry_sane": endpoint_geometry_sane,
        "neutral_base_h_max_A": float(neutral_base_h_max_A),
        "ion_pair_proton_f_max_A": float(ion_pair_proton_f_max_A),
    }


def tilted_hf_chain_seeds(
    xyz: XYZ,
    atoms: dict[str, int],
    spectator_pairs: list[dict[str, int | float]],
    tilt_degrees: float = 12.0,
) -> dict[str, XYZ]:
    """Return deterministic bent orientations of the mapped HF chain.

    The transferring proton is held fixed while the leaving fluorine and the
    downstream HF atoms are rotated about it.  This preserves B--H, H--F, and all
    internal distances on the F side while making B--H--F non-linear.  Two
    orthogonal bend planes and both signs are included because an exactly axial
    seed can have zero transverse gradient on a ridge of the potential surface.
    """

    base = int(atoms["base_atom"])
    transfer_h = int(atoms["transfer_h"])
    leaving_f = int(atoms["leaving_f"])
    required = [base, transfer_h, leaving_f]
    if min(required) < 0 or max(required) >= len(xyz.coords):
        raise IndexError(f"HF-chain atom index outside geometry: {required}")

    moving = {leaving_f}
    for pair in spectator_pairs:
        moving.add(int(pair["h_atom"]))
        moving.add(int(pair["f_atom"]))
    if base in moving:
        raise ValueError("The base atom cannot be part of the rotated HF chain")

    origin = np.asarray(xyz.coords[transfer_h], dtype=float)
    direction = np.asarray(xyz.coords[leaving_f], dtype=float) - origin
    norm = float(np.linalg.norm(direction))
    if norm < 1.0e-10:
        raise ValueError("Cannot bend an HF chain with a zero H--F direction")
    direction /= norm

    axes = np.eye(3, dtype=float)
    reference = min(axes, key=lambda axis: abs(float(np.dot(axis, direction))))
    transverse_1 = reference - float(np.dot(reference, direction)) * direction
    transverse_1 /= float(np.linalg.norm(transverse_1))
    transverse_2 = np.cross(direction, transverse_1)
    transverse_2 /= float(np.linalg.norm(transverse_2))

    def rotation_matrix(axis: np.ndarray, angle_degrees: float) -> np.ndarray:
        axis = np.asarray(axis, dtype=float)
        axis /= float(np.linalg.norm(axis))
        x, y, z = axis
        angle = np.deg2rad(float(angle_degrees))
        c, s = float(np.cos(angle)), float(np.sin(angle))
        cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
        return c * np.eye(3) + (1.0 - c) * np.outer(axis, axis) + s * cross

    variants: list[tuple[str, np.ndarray]] = [
        ("axial", np.eye(3)),
        (
            "tilt_plane1_plus",
            rotation_matrix(np.cross(direction, transverse_1), tilt_degrees),
        ),
        (
            "tilt_plane1_minus",
            rotation_matrix(np.cross(direction, transverse_1), -tilt_degrees),
        ),
        (
            "tilt_plane2_plus",
            rotation_matrix(np.cross(direction, transverse_2), tilt_degrees),
        ),
        (
            "tilt_plane2_minus",
            rotation_matrix(np.cross(direction, transverse_2), -tilt_degrees),
        ),
    ]
    seeds: dict[str, XYZ] = {}
    for label, rotation in variants:
        coords = np.asarray(xyz.coords, dtype=float).copy()
        for index in moving:
            coords[index] = origin + rotation @ (coords[index] - origin)
        seeds[label] = XYZ(
            symbols=list(xyz.symbols),
            coords=coords,
            comment=f"bend_seed={label}; bend_degrees={float(tilt_degrees):.3f}",
        )
    return seeds


def bifluoride_seed_geometry(
    xyz: XYZ,
    atoms: dict[str, int],
    spectator_pairs: list[dict[str, int | float]],
    fraction_from_f1: float = 0.5,
) -> XYZ:
    """Place the second HF proton between F1/F2 to probe an [F--H--F]- branch."""

    if len(spectator_pairs) != 1:
        raise ValueError(
            "A bifluoride seed requires exactly one mapped spectator HF pair"
        )
    fraction = float(fraction_from_f1)
    if not 0.25 <= fraction <= 0.75:
        raise ValueError("fraction_from_f1 must keep H2 between the two fluorines")
    leaving_f = int(atoms["leaving_f"])
    h2 = int(spectator_pairs[0]["h_atom"])
    f2 = int(spectator_pairs[0]["f_atom"])
    coords = np.asarray(xyz.coords, dtype=float).copy()
    coords[h2] = coords[leaving_f] + fraction * (coords[f2] - coords[leaving_f])
    return XYZ(
        symbols=list(xyz.symbols),
        coords=coords,
        comment=f"bifluoride_seed=true; fraction_from_f1={fraction:.3f}",
    )


def artifact_xyz_path(species: Any) -> Path:
    data = getattr(species, "data", {}) or {}
    paths = getattr(species, "paths", {}) or {}
    value = data.get("xyz_path") or paths.get("xyz") or paths.get("final_xyz")
    if not value or not Path(str(value)).exists():
        raise FileNotFoundError(f"Missing endpoint XYZ for {getattr(species, 'artifact_id', 'species')}")
    return Path(str(value))


def reaction_spectator_hf_pairs(
    reaction: Any,
    reference_xyz_path: str | Path,
) -> list[dict[str, int | float]]:
    """Resolve spectator H/F atom identities from reaction metadata.

    Current HF2 reactions carry an explicit ``hf_relay`` mapping.  The
    geometry-based inference is retained only as a compatibility fallback for
    older manifests that lack that mapping.
    """

    data = getattr(reaction, "data", reaction) or {}
    coordinate = data.get("reaction_coordinate", {}) or {}
    relay_value = coordinate.get("hf_relay")
    relays = relay_value if isinstance(relay_value, list) else [relay_value]
    pairs: list[dict[str, int | float]] = []
    for relay in relays:
        relay_atoms = (relay or {}).get("atoms", {}) if relay else {}
        if "spectator_h" in relay_atoms and "terminal_f" in relay_atoms:
            pairs.append(
                {
                    "h_atom": int(relay_atoms["spectator_h"]),
                    "f_atom": int(relay_atoms["terminal_f"]),
                }
            )
    if pairs:
        return pairs
    return infer_spectator_hf_pairs(
        reference_xyz_path,
        reaction_coordinate_atoms(reaction),
    )


def validate_endpoint_pair(reaction: Any, reactant: Any, product: Any) -> dict[str, Any]:
    """Apply the same proton-identity-invariant endpoint gate used upstream."""
    try:
        atoms = reaction_coordinate_atoms(reaction)
        reactant_path = artifact_xyz_path(reactant)
        product_path = artifact_xyz_path(product)
        spectator_pairs = reaction_spectator_hf_pairs(reaction, reactant_path)
        rc = hf_endpoint_metrics_from_xyz(
            read_xyz(reactant_path), atoms, spectator_pairs
        )
        ip = hf_endpoint_metrics_from_xyz(
            read_xyz(product_path), atoms, spectator_pairs
        )
    except Exception as exc:
        return {"accepted": False, "reasons": [f"endpoint_geometry_unreadable:{exc}"]}
    reasons: list[str] = []
    if rc["endpoint_class"] != "neutral_complex":
        reasons.append(f"reactant_is_{rc['endpoint_class']}")
    if ip["endpoint_class"] != "ion_pair":
        reasons.append(f"product_is_{ip['endpoint_class']}")
    rc_q = float(rc["q_HF_minus_BH_identity_invariant_A"])
    ip_q = float(ip["q_HF_minus_BH_identity_invariant_A"])
    if abs(ip_q - rc_q) < 0.30:
        reasons.append("endpoint_coordinate_separation_too_small")
    return {
        "accepted": not reasons,
        "reasons": reasons,
        "classification_scope": "permutation_invariant_hfn_geometry",
        "spectator_pairs": spectator_pairs,
        "reactant": rc,
        "product": ip,
    }
