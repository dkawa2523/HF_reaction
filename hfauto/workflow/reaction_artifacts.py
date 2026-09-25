"""Build and validate backend-independent reaction artifacts."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.reactions import validate_reaction_endpoint_pair
from hfauto.chemistry.stoichiometry import normalize_terms, validate_reaction_balance
from hfauto.core.ids import reaction_id as make_reaction_id
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import (
    ReactionHypothesisRecord,
    StoichiometricTermRecord,
)

SCHEMA_VERSION = "hfauto.define_reactions.v1"


def _reaction_side(
    definition: dict[str, Any],
    side: str,
    default_reference: str,
    states: dict[str, Artifact],
) -> list[StoichiometricTermRecord]:
    raw_side = definition.get(side)
    if raw_side is None:
        raw_side = [{"species": default_reference, "coefficient": 1.0}]
    if isinstance(raw_side, dict):
        raw_side = [
            {"species": reference, "coefficient": coefficient}
            for reference, coefficient in raw_side.items()
        ]
    terms: list[StoichiometricTermRecord] = []
    for raw in raw_side:
        if isinstance(raw, str):
            raw = {"species": raw}
        reference = str(raw.get("species") or raw.get("species_id") or "")
        if reference not in states:
            raise ValueError(f"{side} species not found: {reference!r}")
        artifact = states[reference]
        terms.append(
            StoichiometricTermRecord(
                species_id=str(
                    artifact.data.get("species_id") or artifact.artifact_id
                ),
                coefficient=raw.get("coefficient", 1.0),
                role=raw.get("role"),
            )
        )
    return normalize_terms(terms)


def build_reaction_hypothesis(
    definition: dict[str, Any], states: dict[str, Artifact]
) -> Artifact:
    """Create one balanced, atom-mapped reaction hypothesis."""

    reactant_ref = str(definition["reactant"])
    product_ref = str(definition["product"])
    if reactant_ref not in states or product_ref not in states:
        raise ValueError(
            "reaction endpoints not found: "
            f"reactant={reactant_ref!r}, product={product_ref!r}"
        )
    reactant = states[reactant_ref]
    product = states[product_ref]
    reactants = _reaction_side(definition, "reactants", reactant_ref, states)
    products = _reaction_side(definition, "products", product_ref, states)
    reaction_id = str(
        definition.get("reaction_id")
        or make_reaction_id(reactant_ref, product_ref)
    )
    mechanism_family = str(
        definition.get("mechanism_family", "unspecified")
    )
    hypothesis = ReactionHypothesisRecord(
        reaction_id=reaction_id,
        reactant_state_id=str(reactant.data["species_id"]),
        product_state_id=str(product.data["species_id"]),
        mechanism_family=mechanism_family,
        bond_changes=list(definition.get("bond_changes", []) or []),
        reactants=reactants,
        products=products,
        reaction_coordinate=dict(
            definition.get("reaction_coordinate", {}) or {}
        ),
        rationale=list(definition.get("rationale", []) or []),
        environment=dict(definition.get("environment", {}) or {}),
    )
    data = {
        "reaction_id": reaction_id,
        "reaction_type": str(
            definition.get("reaction_type", mechanism_family)
        ),
        "mechanism_family": mechanism_family,
        "reactant_species_id": reactant.data["species_id"],
        "product_species_id": product.data["species_id"],
        "ts_species_id": None,
        "bond_changes": [
            change.model_dump() for change in hypothesis.bond_changes
        ],
        "stoichiometry": {
            "reactants": [
                term.model_dump() for term in hypothesis.reactants
            ],
            "products": [term.model_dump() for term in hypothesis.products],
        },
        "reaction_coordinate": hypothesis.reaction_coordinate.model_dump(),
        "hypothesis": hypothesis.model_dump(),
        "environment": hypothesis.environment,
        "schema_version": SCHEMA_VERSION,
    }
    for key in ("candidate_id", "comparison_group"):
        if definition.get(key) is not None:
            data[key] = definition[key]
    endpoint_qc = validate_reaction_endpoint_pair(data, reactant, product)
    if not endpoint_qc["accepted"]:
        raise ValueError(
            f"reaction {reaction_id!r} has invalid endpoints: "
            + ", ".join(endpoint_qc["reasons"])
        )
    balance_qc = validate_reaction_balance(
        hypothesis.reactants, hypothesis.products, states
    )
    if not balance_qc["accepted"]:
        raise ValueError(
            f"reaction {reaction_id!r} is unbalanced: "
            + ", ".join(balance_qc["reasons"])
        )
    return Artifact(
        artifact_id=reaction_id,
        artifact_type="reaction",
        parents=[reactant.artifact_id, product.artifact_id],
        data=data,
        qc={
            "endpoint_hypothesis_validated": True,
            "endpoint_validation": endpoint_qc,
            "stoichiometry_validated": True,
            "stoichiometry_validation": balance_qc,
        },
    )
