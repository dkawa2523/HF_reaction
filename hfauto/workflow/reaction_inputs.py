"""Resolve registry-promoted reactions and their validated workflow inputs."""

from __future__ import annotations

from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


def resolve_reaction_endpoints(
    source: Manifest,
    reaction: Artifact,
) -> tuple[Artifact | None, Artifact | None, dict]:
    """Resolve the representative minima recorded on a promoted reaction."""

    species = preferred_species_by_id(source)
    reactant_id = str(reaction.data.get("reactant_species_id") or "")
    product_id = str(reaction.data.get("product_species_id") or "")
    reactant = species.get(reactant_id)
    product = species.get(product_id)
    reasons: list[str] = []
    if reaction.qc.get("distinct_registry_basins") is not True:
        reasons.append("reaction_was_not_promoted_from_distinct_registry_basins")
    if reactant is None:
        reasons.append("declared_reactant_species_missing")
    if product is None:
        reasons.append("declared_product_species_missing")
    return reactant, product, {
        "accepted": not reasons,
        "reasons": reasons,
        "source": "minimum_registry_reaction",
        "selection_artifact_id": None,
    }


def reaction_and_endpoints(
    source: Manifest, reaction_id: str
) -> tuple[Artifact, Artifact, Artifact]:
    reaction = next(
        (
            item
            for item in source.latest_artifacts("reaction")
            if str(item.data.get("reaction_id") or item.artifact_id)
            == reaction_id
        ),
        None,
    )
    if reaction is None:
        raise KeyError(f"Reaction not found: {reaction_id}")
    reactant, product, gate = resolve_reaction_endpoints(source, reaction)
    if not gate["accepted"] or reactant is None or product is None:
        raise KeyError(
            "Validated endpoint species missing or rejected: "
            + ", ".join(gate["reasons"])
        )
    return reaction, reactant, product


def basin_assessment(source: Manifest, reaction_id: str) -> dict:
    reaction = next(
        (
            item
            for item in source.latest_artifacts("reaction")
            if str(item.data.get("reaction_id") or item.artifact_id)
            == reaction_id
        ),
        None,
    )
    if reaction is not None:
        assessment = dict(reaction.data.get("basin_assessment") or {})
        if assessment.get("accepted") is True:
            return assessment
    for item in reversed(source.latest_artifacts("reaction_case")):
        if str(item.data.get("reaction_id") or "") != reaction_id:
            continue
        assessment = dict(item.data.get("basin_assessment") or {})
        if assessment.get("accepted") is True:
            return assessment
    raise KeyError(f"Accepted registry basin assessment not found: {reaction_id}")


def saddle_seed_candidate(
    source: Manifest, reaction_id: str
) -> tuple[dict, Artifact]:
    for attempt in reversed(source.latest_artifacts("saddle_attempt")):
        candidate = dict(attempt.data.get("candidate") or {})
        accepted = bool(
            attempt.status.status == "success"
            and str(attempt.data.get("reaction_id") or "") == reaction_id
            and attempt.data.get("diagnosis") == "resolved_saddle_candidate"
            and attempt.qc.get("saddle_seed_resolved") is True
            and attempt.qc.get("method_evidence_validated") is True
            and candidate.get("xyz_path")
        )
        if accepted:
            return candidate, attempt
    raise KeyError(f"Validated saddle seed not found: {reaction_id}")


def validated_ts_and_endpoints(
    source: Manifest, reaction_id: str
) -> tuple[Artifact, Artifact, Artifact, Artifact, Artifact]:
    reaction, reactant, product = reaction_and_endpoints(source, reaction_id)
    validated = next(
        (
            item
            for item in reversed(source.latest_artifacts("reaction_validated"))
            if str(item.data.get("reaction_id") or "") == reaction_id
            and item.status.status == "success"
            and item.qc.get("ts_validated_by_frequency") is True
            and item.qc.get("n_imag") == 1
        ),
        None,
    )
    if validated is None:
        raise KeyError(f"Frequency-validated TS not found: {reaction_id}")
    ts_species_id = str(validated.data.get("ts_species_id") or "")
    ts_species = source.find(ts_species_id)
    if (
        ts_species is None
        or ts_species.artifact_type != "species"
        or ts_species.status.status != "success"
        or ts_species.qc.get("method_evidence_validated") is not True
    ):
        raise KeyError(f"Validated TS species missing: {ts_species_id}")
    return reaction, ts_species, reactant, product, validated
