"""Bounded, product-free reaction trial generation.

The module proposes internal-coordinate drives only.  It never constructs or
labels a product geometry; that belongs to a single-ended discovery backend.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np

from hfauto.chemistry.connectivity import covalent_adjacency
from hfauto.chemistry.xyz import XYZ
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.chemistry import (
    DrivingAtomPairRecord,
    ReactionCoordinateRecord,
    ReactionTrialRecord,
)

_ACCEPTOR_ELEMENTS = {"N", "O", "F", "P", "S", "Cl", "Br", "I"}
_ELECTROPHILE_ELEMENTS = {"B", "C", "Si", "P", "S"}
_LEAVING_GROUP_ELEMENTS = {"N", "O", "F", "S", "Cl", "Br", "I"}


def _component_ids(
    atom_count: int, components: Iterable[Mapping[str, Any]]
) -> tuple[list[str], list[str]]:
    atom_components = ["system"] * atom_count
    ordered: list[str] = []
    for component_index, component in enumerate(components):
        component_id = str(component.get("component_id") or f"fragment_{component_index}")
        ordered.append(component_id)
        for atom_index in component.get("atom_indices", []):
            index = int(atom_index)
            if not 0 <= index < atom_count:
                raise IndexError("component atom index is outside the source geometry")
            atom_components[index] = component_id
    return atom_components, ordered or ["system"]


def _pair(first: int, second: int, label: str) -> DrivingAtomPairRecord:
    return DrivingAtomPairRecord(atoms=(first, second), label=label)


def _trial(
    *,
    source_species_id: str,
    associations: list[DrivingAtomPairRecord],
    dissociations: list[DrivingAtomPairRecord],
    mechanism_hint: str,
    priority: float,
    rationale: list[str],
    component_ids: list[str],
    charge: int,
    multiplicity: int,
    driver_order: list[str],
    max_attempts: int,
    reaction_coordinate: ReactionCoordinateRecord | None = None,
    settings: Mapping[str, Any] | None = None,
) -> ReactionTrialRecord:
    definition = {
        "source_species_id": source_species_id,
        "associations": [pair.atoms for pair in associations],
        "dissociations": [pair.atoms for pair in dissociations],
        "reaction_coordinate": (
            reaction_coordinate.model_dump(mode="json")
            if reaction_coordinate is not None
            else None
        ),
        "mechanism_hint": mechanism_hint,
        "charge": int(charge),
        "multiplicity": int(multiplicity),
    }
    return ReactionTrialRecord(
        trial_id="trial_" + fingerprint_dict(definition)[:20],
        source_species_id=source_species_id,
        associations=associations,
        dissociations=dissociations,
        reaction_coordinate=reaction_coordinate or ReactionCoordinateRecord(),
        driver_order=driver_order,
        mechanism_hint=mechanism_hint,
        priority=float(priority),
        rationale=rationale,
        component_ids=component_ids,
        charge=int(charge),
        multiplicity=int(multiplicity),
        max_attempts=int(max_attempts),
        settings=dict(settings or {}),
    )


def explicit_reaction_trials(
    source_species_id: str,
    atom_count: int,
    definitions: Iterable[Mapping[str, Any]],
    *,
    component_ids: list[str],
    charge: int,
    multiplicity: int,
    default_driver_order: list[str],
    default_max_attempts: int,
) -> list[ReactionTrialRecord]:
    """Normalize reviewed driving coordinates into the same trial contract."""

    trials: list[ReactionTrialRecord] = []
    for definition in definitions:
        selector = definition.get("source_species_id") or definition.get("system_id")
        if selector and str(selector) not in {source_species_id, source_species_id.removeprefix("spc_")}:
            continue
        associations = [
            _pair(int(atoms[0]), int(atoms[1]), "explicit_association")
            for atoms in definition.get("associations", [])
        ]
        dissociations = [
            _pair(int(atoms[0]), int(atoms[1]), "explicit_dissociation")
            for atoms in definition.get("dissociations", [])
        ]
        if any(max(pair.atoms) >= atom_count for pair in [*associations, *dissociations]):
            raise IndexError("explicit driving atom is outside the source geometry")
        reaction_coordinate = ReactionCoordinateRecord.model_validate(
            definition.get("reaction_coordinate") or {}
        )
        if any(
            max(term.atoms) >= atom_count
            for term in reaction_coordinate.terms
        ):
            raise IndexError(
                "explicit reaction-coordinate atom is outside the source geometry"
            )
        trials.append(
            _trial(
                source_species_id=source_species_id,
                associations=associations,
                dissociations=dissociations,
                mechanism_hint=str(definition.get("mechanism_hint", "reviewed")),
                priority=float(definition.get("priority", 100.0)),
                rationale=["Reviewed explicit driving-coordinate definition."],
                component_ids=component_ids,
                charge=charge,
                multiplicity=multiplicity,
                driver_order=list(definition.get("driver_order") or default_driver_order),
                max_attempts=int(definition.get("max_attempts", default_max_attempts)),
                reaction_coordinate=reaction_coordinate,
                settings=dict(definition.get("settings") or {}),
            )
        )
    return trials


def automatic_reaction_trials(
    source_species_id: str,
    xyz: XYZ,
    components: Iterable[Mapping[str, Any]],
    *,
    charge: int,
    multiplicity: int,
    driver_order: list[str],
    max_attempts: int,
    max_relay_steps: int = 2,
    relay_contact_cutoff_A: float = 3.2,
) -> list[ReactionTrialRecord]:
    """Propose bounded transfer, relay, substitution, and dissociation drives."""

    adjacency = covalent_adjacency(xyz)
    atom_components, component_ids = _component_ids(len(xyz.symbols), components)
    acceptors = [
        index
        for index, symbol in enumerate(xyz.symbols)
        if symbol in _ACCEPTOR_ELEMENTS
    ]
    donor_bonds: list[tuple[int, int]] = []
    for hydrogen, symbol in enumerate(xyz.symbols):
        if symbol != "H":
            continue
        donor_bonds.extend(
            (heavy, hydrogen)
            for heavy in range(len(xyz.symbols))
            if adjacency[heavy, hydrogen]
            and xyz.symbols[heavy] in _ACCEPTOR_ELEMENTS
        )

    trials: list[ReactionTrialRecord] = []
    for donor, hydrogen in donor_bonds:
        for acceptor in acceptors:
            if acceptor == donor or adjacency[acceptor, hydrogen]:
                continue
            intercomponent = atom_components[acceptor] != atom_components[donor]
            distance = float(np.linalg.norm(xyz.coords[acceptor] - xyz.coords[hydrogen]))
            trials.append(
                _trial(
                    source_species_id=source_species_id,
                    associations=[_pair(acceptor, hydrogen, "acceptor-H")],
                    dissociations=[_pair(donor, hydrogen, "donor-H")],
                    mechanism_hint="proton_transfer",
                    priority=(20.0 if intercomponent else 12.0) - min(distance, 8.0),
                    rationale=[
                        "A heteroatom-bound hydrogen can transfer to an unbound acceptor.",
                        "The drive is a search prior, not a product or barrier claim.",
                    ],
                    component_ids=component_ids,
                    charge=charge,
                    multiplicity=multiplicity,
                    driver_order=driver_order,
                    max_attempts=max_attempts,
                )
            )

    if max_relay_steps >= 2:
        for first_donor, first_hydrogen in donor_bonds:
            for second_donor, second_hydrogen in donor_bonds:
                if len({first_donor, first_hydrogen, second_donor, second_hydrogen}) < 4:
                    continue
                relay_distance = float(
                    np.linalg.norm(xyz.coords[first_donor] - xyz.coords[second_hydrogen])
                )
                if relay_distance > float(relay_contact_cutoff_A):
                    continue
                for acceptor in acceptors:
                    if acceptor in {first_donor, second_donor}:
                        continue
                    if atom_components[acceptor] == atom_components[first_donor]:
                        continue
                    trials.append(
                        _trial(
                            source_species_id=source_species_id,
                            associations=[
                                _pair(acceptor, first_hydrogen, "relay_acceptor-H"),
                                _pair(first_donor, second_hydrogen, "relay_bridge-H"),
                            ],
                            dissociations=[
                                _pair(first_donor, first_hydrogen, "relay_first_donor-H"),
                                _pair(second_donor, second_hydrogen, "relay_second_donor-H"),
                            ],
                            mechanism_hint="proton_relay",
                            priority=18.0 - relay_distance,
                            rationale=[
                                "Two hydrogen-donor bonds form a geometrically connected relay.",
                                "All association and dissociation pairs are driven together by NT2.",
                            ],
                            component_ids=component_ids,
                            charge=charge,
                            multiplicity=multiplicity,
                            driver_order=["nt2"],
                            max_attempts=max_attempts,
                        )
                    )

    for electrophile, symbol in enumerate(xyz.symbols):
        if symbol not in _ELECTROPHILE_ELEMENTS:
            continue
        leaving_atoms = [
            other
            for other in range(len(xyz.symbols))
            if adjacency[electrophile, other]
            and xyz.symbols[other] in _LEAVING_GROUP_ELEMENTS
        ]
        for leaving in leaving_atoms:
            for nucleophile in acceptors:
                if nucleophile in {electrophile, leaving} or adjacency[nucleophile, electrophile]:
                    continue
                trials.append(
                    _trial(
                        source_species_id=source_species_id,
                        associations=[_pair(nucleophile, electrophile, "Nu-E")],
                        dissociations=[_pair(electrophile, leaving, "E-LG")],
                        mechanism_hint="substitution",
                        priority=(15.0 if atom_components[nucleophile] != atom_components[electrophile] else 9.0),
                        rationale=["A nucleophilic heteroatom is paired with a bonded leaving group."],
                        component_ids=component_ids,
                        charge=charge,
                        multiplicity=multiplicity,
                        driver_order=driver_order,
                        max_attempts=max_attempts,
                    )
                )

    for first in range(len(xyz.symbols)):
        for second in range(first):
            if not adjacency[first, second]:
                continue
            if xyz.symbols[first] in _LEAVING_GROUP_ELEMENTS or xyz.symbols[second] in _LEAVING_GROUP_ELEMENTS:
                trials.append(
                    _trial(
                        source_species_id=source_species_id,
                        associations=[],
                        dissociations=[_pair(first, second, "polar_bond_dissociation")],
                        mechanism_hint="dissociation",
                        priority=2.0,
                        rationale=["A heteroatom-containing covalent bond is a bounded dissociation trial."],
                        component_ids=component_ids,
                        charge=charge,
                        multiplicity=multiplicity,
                        driver_order=driver_order,
                        max_attempts=max_attempts,
                    )
                )
    return trials


def deduplicate_trials(trials: Iterable[ReactionTrialRecord]) -> list[ReactionTrialRecord]:
    """Keep the highest-priority instance of each normalized drive."""

    unique: dict[tuple[Any, ...], ReactionTrialRecord] = {}
    for trial in trials:
        key = (
            trial.source_species_id,
            tuple(sorted(pair.atoms for pair in trial.associations)),
            tuple(sorted(pair.atoms for pair in trial.dissociations)),
            tuple(
                (term.kind, term.atoms, float(term.coefficient))
                for term in trial.reaction_coordinate.terms
            ),
        )
        current = unique.get(key)
        if current is None or trial.priority > current.priority:
            unique[key] = trial
    return sorted(unique.values(), key=lambda item: (-item.priority, item.trial_id))
