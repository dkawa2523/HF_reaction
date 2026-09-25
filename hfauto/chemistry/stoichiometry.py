"""Chemical-equation validation and linear energy bookkeeping."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable
from typing import Any

from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import StoichiometricTermRecord


def molecular_surface_key(species: Artifact) -> str:
    """Return a composition/electronic-state key for energy comparisons.

    Absolute electronic energies are comparable only for structures with the
    same elemental composition, charge, and multiplicity.  The key is kept
    independent of component labels so that isomers and association topologies
    on the same molecular potential-energy surface remain comparable.
    """

    raw_counts = species.data.get("element_counts") or {}
    counts = {
        str(symbol): int(count)
        for symbol, count in dict(raw_counts).items()
        if int(count) != 0
    }
    if not counts:
        try:
            xyz = read_xyz(species_xyz_path(species))
        except (OSError, TypeError, ValueError):
            identity = str(
                species.data.get("species_id") or species.artifact_id
            )
            return f"unknown:{identity}"
        for symbol in xyz.symbols:
            counts[symbol] = counts.get(symbol, 0) + 1

    chemical_state = species.data.get("chemical_state") or {}
    charge = int(species.data.get("charge", chemical_state.get("charge", 0)))
    multiplicity = int(
        species.data.get(
            "multiplicity",
            chemical_state.get("multiplicity", 1),
        )
    )
    formula = "".join(
        f"{symbol}{counts[symbol]}" for symbol in sorted(counts)
    )
    return f"{formula}|q={charge}|m={multiplicity}"


def normalize_terms(
    terms: Iterable[StoichiometricTermRecord | dict[str, Any]],
) -> list[StoichiometricTermRecord]:
    """Validate terms and merge repeated species while preserving first order."""

    order: list[str] = []
    coefficients: dict[str, float] = defaultdict(float)
    roles: dict[str, str | None] = {}
    for raw in terms:
        term = (
            raw
            if isinstance(raw, StoichiometricTermRecord)
            else StoichiometricTermRecord.model_validate(raw)
        )
        if term.species_id not in coefficients:
            order.append(term.species_id)
            roles[term.species_id] = term.role
        coefficients[term.species_id] += term.coefficient
    return [
        StoichiometricTermRecord(
            species_id=species_id,
            coefficient=coefficients[species_id],
            role=roles[species_id],
        )
        for species_id in order
    ]


def reaction_side_terms(
    reaction: Artifact,
    side: str,
) -> list[StoichiometricTermRecord]:
    """Return normalized terms from a reaction, including the legacy single id."""

    if side not in {"reactants", "products"}:
        raise ValueError("reaction side must be 'reactants' or 'products'")
    terms = list((reaction.data.get("stoichiometry") or {}).get(side) or [])
    if terms:
        return normalize_terms(terms)
    key = "reactant_species_id" if side == "reactants" else "product_species_id"
    species_id = reaction.data.get(key)
    return (
        [StoichiometricTermRecord(species_id=str(species_id), coefficient=1.0)]
        if species_id
        else []
    )


def stoichiometric_sum(
    terms: Iterable[StoichiometricTermRecord | dict[str, Any]],
    value_for_species: Callable[[str], float | None],
) -> float | None:
    """Return a side total, or ``None`` when any required value is missing."""

    total = 0.0
    for term in normalize_terms(terms):
        value = value_for_species(term.species_id)
        if value is None:
            return None
        total += term.coefficient * float(value)
    return total


def reaction_delta(
    reactants: Iterable[StoichiometricTermRecord | dict[str, Any]],
    products: Iterable[StoichiometricTermRecord | dict[str, Any]],
    value_for_species: Callable[[str], float | None],
) -> float | None:
    """Return products minus reactants in the source value's unit."""

    reactant_value = stoichiometric_sum(reactants, value_for_species)
    product_value = stoichiometric_sum(products, value_for_species)
    if reactant_value is None or product_value is None:
        return None
    return product_value - reactant_value


def validate_reaction_balance(
    reactants: Iterable[StoichiometricTermRecord | dict[str, Any]],
    products: Iterable[StoichiometricTermRecord | dict[str, Any]],
    species: dict[str, Artifact],
    *,
    tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    """Validate elemental and charge balance from the actual species geometries."""

    reasons: list[str] = []

    def side_totals(
        terms: Iterable[StoichiometricTermRecord | dict[str, Any]],
    ) -> tuple[dict[str, float], float]:
        elements: dict[str, float] = defaultdict(float)
        charge = 0.0
        for term in normalize_terms(terms):
            artifact = species.get(term.species_id)
            if artifact is None:
                reasons.append(f"stoichiometric_species_missing:{term.species_id}")
                continue
            try:
                xyz = read_xyz(species_xyz_path(artifact))
            except (OSError, TypeError, ValueError) as exc:
                reasons.append(f"stoichiometric_geometry_unreadable:{term.species_id}:{exc}")
                continue
            for symbol in xyz.symbols:
                elements[symbol] += term.coefficient
            charge += term.coefficient * float(artifact.data.get("charge", 0))
        return dict(elements), charge

    reactant_elements, reactant_charge = side_totals(reactants)
    product_elements, product_charge = side_totals(products)
    element_delta = {
        symbol: product_elements.get(symbol, 0.0) - reactant_elements.get(symbol, 0.0)
        for symbol in sorted(set(reactant_elements) | set(product_elements))
    }
    if any(abs(value) > tolerance for value in element_delta.values()):
        reasons.append("elemental_composition_not_balanced")
    charge_delta = product_charge - reactant_charge
    if abs(charge_delta) > tolerance:
        reasons.append("total_charge_not_balanced")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "element_delta": element_delta,
        "charge_delta": charge_delta,
    }
