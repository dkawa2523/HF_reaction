"""Connect validated minima using chemistry-neutral reaction trials.

Single-ended xTB exploration can miss a DFT minimum when the two methods have
different basin topology.  This module compares *validated DFT minima* against
the atom-pair changes already declared by a reaction trial.  It never invents
a barrier. Isomorphic covalent graphs may be linked only as a conservative
continuous-coordinate reorganization, without asserted bond changes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.basin_identity import (
    permutation_invariant_graph_rmsd,
)
from hfauto.chemistry.connectivity import covalent_adjacency
from hfauto.chemistry.reactions import (
    validate_reaction_coordinate_between_geometries,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.schemas.chemistry import ReactionTrialRecord


def _pair_change_satisfied(
    source_graph: Any,
    target_graph: Any,
    atoms: tuple[int, int],
    *,
    form: bool,
) -> bool:
    source_bonded = bool(source_graph[atoms])
    target_bonded = bool(target_graph[atoms])
    return (
        not source_bonded and target_bonded
        if form
        else source_bonded and not target_bonded
    )


def assess_trial_between_minima(
    trial: ReactionTrialRecord,
    source_xyz_path: str | Path,
    target_xyz_path: str | Path,
    *,
    minimum_coordinate_change_A: float = 0.05,
    minimum_bond_distance_change_A: float = 0.15,
) -> dict[str, Any]:
    """Test both directions of one trial between two atom-mapped minima.

    The reverse direction matters when the low-level trial source relaxes into
    the opposite DFT basin.  Reversing a connection changes the reported bond
    changes and reaction-coordinate sign; it does not alter either geometry.
    """

    source = read_xyz(source_xyz_path)
    target = read_xyz(target_xyz_path)
    if source.symbols != target.symbols:
        return {
            "accepted": False,
            "reasons": ["minimum_atom_order_or_composition_mismatch"],
        }

    atom_count = len(source.symbols)
    requested_pairs = [
        *(pair.atoms for pair in trial.associations),
        *(pair.atoms for pair in trial.dissociations),
    ]
    if any(max(pair) >= atom_count for pair in requested_pairs):
        return {
            "accepted": False,
            "reasons": ["trial_atom_index_outside_minimum_geometry"],
        }

    if trial.reaction_coordinate.terms:
        explicit_coordinate = trial.reaction_coordinate.model_dump(mode="json")
        forward_coordinate = validate_reaction_coordinate_between_geometries(
            {"reaction_coordinate": explicit_coordinate},
            source_xyz_path,
            target_xyz_path,
        )
        reversed_coordinate = {
            **explicit_coordinate,
            "terms": [
                {**term, "coefficient": -float(term.get("coefficient", 1.0))}
                for term in explicit_coordinate["terms"]
            ],
        }
        reverse_coordinate = validate_reaction_coordinate_between_geometries(
            {"reaction_coordinate": reversed_coordinate},
            source_xyz_path,
            target_xyz_path,
        )
        if not forward_coordinate["accepted"] and not reverse_coordinate["accepted"]:
            return {
                "accepted": False,
                "reasons": ["trial_coordinates_do_not_connect_minima"],
                "forward_coordinate_evidence": forward_coordinate,
                "reverse_coordinate_evidence": reverse_coordinate,
                "requested_change_count": len(trial.reaction_coordinate.terms),
            }
        direction = "forward" if forward_coordinate["accepted"] else "reverse"
        return {
            "accepted": True,
            "reasons": [],
            "direction": direction,
            "connection_kind": "coordinate_reorganization",
            "requested_change_count": len(trial.reaction_coordinate.terms),
            "forward_coordinate_evidence": forward_coordinate,
            "reverse_coordinate_evidence": reverse_coordinate,
            "bond_changes": [],
            "reaction_coordinate": (
                explicit_coordinate
                if direction == "forward"
                else reversed_coordinate
            ),
        }

    source_graph = covalent_adjacency(source)
    target_graph = covalent_adjacency(target)
    graph_rmsd = permutation_invariant_graph_rmsd(
        source_xyz_path, target_xyz_path
    )

    def distance(atoms: tuple[int, int], *, target_geometry: bool) -> float:
        geometry = target if target_geometry else source
        return float(np.linalg.norm(geometry.coords[atoms[0]] - geometry.coords[atoms[1]]))

    forward_deltas = [
        *(
            distance(pair.atoms, target_geometry=False)
            - distance(pair.atoms, target_geometry=True)
            for pair in trial.associations
        ),
        *(
            distance(pair.atoms, target_geometry=True)
            - distance(pair.atoms, target_geometry=False)
            for pair in trial.dissociations
        ),
    ]
    reverse_deltas = [-value for value in forward_deltas]

    forward = [
        *(
            _pair_change_satisfied(
                source_graph, target_graph, pair.atoms, form=True
            )
            for pair in trial.associations
        ),
        *(
            _pair_change_satisfied(
                source_graph, target_graph, pair.atoms, form=False
            )
            for pair in trial.dissociations
        ),
    ]
    reverse = [
        *(
            _pair_change_satisfied(
                target_graph, source_graph, pair.atoms, form=True
            )
            for pair in trial.associations
        ),
        *(
            _pair_change_satisfied(
                target_graph, source_graph, pair.atoms, form=False
            )
            for pair in trial.dissociations
        ),
    ]
    forward_bond_change = bool(
        forward
        and all(forward)
        and all(
            value >= float(minimum_bond_distance_change_A)
            for value in forward_deltas
        )
    )
    reverse_bond_change = bool(
        reverse
        and all(reverse)
        and all(
            value >= float(minimum_bond_distance_change_A)
            for value in reverse_deltas
        )
    )
    forward_coordinate_change = bool(
        forward_deltas
        and all(value > 0.0 for value in forward_deltas)
        and sum(forward_deltas) >= float(minimum_coordinate_change_A)
    )
    reverse_coordinate_change = bool(
        reverse_deltas
        and all(value > 0.0 for value in reverse_deltas)
        and sum(reverse_deltas) >= float(minimum_coordinate_change_A)
    )
    forward_accepted = forward_bond_change or forward_coordinate_change
    reverse_accepted = reverse_bond_change or reverse_coordinate_change
    if not forward_accepted and not reverse_accepted:
        return {
            "accepted": False,
            "reasons": ["trial_coordinates_do_not_connect_minima"],
            "forward_changes_satisfied": sum(forward),
            "reverse_changes_satisfied": sum(reverse),
            "requested_change_count": len(forward),
            "forward_coordinate_progress_A": sum(forward_deltas),
            "reverse_coordinate_progress_A": sum(reverse_deltas),
            "permutation_invariant_graph_rmsd_A": graph_rmsd,
        }

    direction = "forward" if forward_accepted else "reverse"
    sign = 1.0 if direction == "forward" else -1.0
    connection_kind = (
        "bond_rearrangement"
        if (forward_bond_change if direction == "forward" else reverse_bond_change)
        else "coordinate_reorganization"
    )
    declared_bond_changes = [
        *(
            {
                "kind": "form" if direction == "forward" else "break",
                "atoms": list(pair.atoms),
                "label": pair.label,
            }
            for pair in trial.associations
        ),
        *(
            {
                "kind": "break" if direction == "forward" else "form",
                "atoms": list(pair.atoms),
                "label": pair.label,
            }
            for pair in trial.dissociations
        ),
    ]
    bond_changes = (
        declared_bond_changes if connection_kind == "bond_rearrangement" else []
    )
    coordinate_terms = [
        *(
            {
                "kind": "distance",
                "atoms": list(pair.atoms),
                "coefficient": -1.0 * sign,
                "label": pair.label,
            }
            for pair in trial.associations
        ),
        *(
            {
                "kind": "distance",
                "atoms": list(pair.atoms),
                "coefficient": 1.0 * sign,
                "label": pair.label,
            }
            for pair in trial.dissociations
        ),
    ]
    return {
        "accepted": True,
        "reasons": [],
        "direction": direction,
        "connection_kind": connection_kind,
        "requested_change_count": len(forward),
        "forward_changes_satisfied": sum(forward),
        "reverse_changes_satisfied": sum(reverse),
        "forward_coordinate_progress_A": sum(forward_deltas),
        "reverse_coordinate_progress_A": sum(reverse_deltas),
        "minimum_coordinate_change_A": float(minimum_coordinate_change_A),
        "minimum_bond_distance_change_A": float(
            minimum_bond_distance_change_A
        ),
        "permutation_invariant_graph_rmsd_A": graph_rmsd,
        "bond_changes": bond_changes,
        "reaction_coordinate": {
            "terms": coordinate_terms,
            "min_change": float(minimum_coordinate_change_A),
        },
    }
