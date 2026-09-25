"""Promote only same-PES pairs of distinct validated minima to reactions."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.minima import is_accepted_optimized_minimum
from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.ids import path_token
from hfauto.core.ids import reaction_id as make_reaction_id
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_artifacts import build_reaction_hypothesis
from hfauto.workflow.reaction_evidence import (
    artifacts_by_reaction,
    attempts_by_reaction,
)
from hfauto.workflow.reaction_state import decide_reaction_case


def make_reaction_case_artifact(
    case: dict[str, Any],
    *,
    parents: list[str],
    created_by: str,
) -> Artifact:
    """Publish one reaction-case revision with a common status contract."""

    token = path_token(str(case["reaction_id"]), max_length=48)
    ready = bool(
        case["workflow_complete"]
        or case["ts_search_allowed"]
        or case["irc_allowed"]
    )
    return Artifact(
        artifact_id=f"reaction_case_{token}",
        artifact_type="reaction_case",
        parents=list(dict.fromkeys(parents)),
        data=case,
        qc={
            "ts_search_allowed": case["ts_search_allowed"],
            "irc_allowed": case["irc_allowed"],
            "workflow_complete": case["workflow_complete"],
            "scientific_conclusion_supported": case[
                "scientific_conclusion_supported"
            ],
        },
        provenance={"created_by": created_by},
        status=ArtifactStatus(
            status="success" if ready else "partial",
            category=None if ready else "reaction_case_not_ready",
            reason="; ".join(case["reasons"]),
        ),
    )


def _registry(manifest: Manifest) -> Artifact | None:
    registries = manifest.latest_artifacts("minimum_registry")
    return registries[-1] if registries else None


def _basin_assessment(
    candidate: Artifact,
    registry: Artifact | None,
    basins: dict[str, Artifact],
) -> tuple[dict[str, Any], list[str]]:
    """Resolve a candidate through the registry without repeating basin logic."""

    candidate_id = str(candidate.data.get("candidate_id") or candidate.artifact_id)
    parents = [candidate.artifact_id]
    if registry is None:
        return (
            {
                "candidate_id": candidate_id,
                "status": "invalid_endpoints",
                "accepted": False,
                "reasons": ["minimum_registry_missing"],
            },
            parents,
        )

    parents.append(registry.artifact_id)
    species_to_basin = dict(registry.data.get("species_to_basin") or {})
    reactant_species_id = str(candidate.data["reactant_species_id"])
    product_species_id = str(candidate.data["product_species_id"])
    reactant_basin_id = species_to_basin.get(reactant_species_id)
    product_basin_id = species_to_basin.get(product_species_id)
    common = {
        "candidate_id": candidate_id,
        "reactant_species_id": reactant_species_id,
        "product_species_id": product_species_id,
        "reactant_basin_id": reactant_basin_id,
        "product_basin_id": product_basin_id,
        "registry_id": registry.artifact_id,
    }
    missing_roles = [
        role
        for role, basin_id in (
            ("reactant", reactant_basin_id),
            ("product", product_basin_id),
        )
        if basin_id is None
    ]
    if missing_roles:
        return (
            {
                **common,
                "status": "invalid_endpoints",
                "accepted": False,
                "same_basin": False,
                "distinct_basin": False,
                "reasons": [
                    f"{role}_frequency_validated_minimum_missing"
                    for role in missing_roles
                ],
            },
            parents,
        )

    reactant_basin = basins.get(str(reactant_basin_id))
    product_basin = basins.get(str(product_basin_id))
    if reactant_basin is None or product_basin is None:
        return (
            {
                **common,
                "status": "invalid_endpoints",
                "accepted": False,
                "same_basin": False,
                "distinct_basin": False,
                "reasons": ["minimum_basin_artifact_missing"],
            },
            parents,
        )
    parents.extend([reactant_basin.artifact_id, product_basin.artifact_id])

    if reactant_basin_id == product_basin_id:
        return (
            {
                **common,
                "status": "same_basin",
                "accepted": True,
                "same_basin": True,
                "distinct_basin": False,
                "reasons": ["registry_assigns_both_endpoints_to_one_basin"],
            },
            parents,
        )

    signatures_match = (
        reactant_basin.data.get("signature")
        == product_basin.data.get("signature")
    )
    ambiguous = bool(
        product_basin_id
        in (reactant_basin.data.get("unresolved_against") or [])
        or reactant_basin_id
        in (product_basin.data.get("unresolved_against") or [])
    )
    endpoint_evidence = {
        "reactant_evidence": {
            "basin_id": reactant_basin.artifact_id,
            "electronic_energy_hartree": reactant_basin.data.get(
                "electronic_energy_hartree"
            ),
            "calculation_artifact_ids": list(
                reactant_basin.data.get("calculation_artifact_ids") or []
            ),
        },
        "product_evidence": {
            "basin_id": product_basin.artifact_id,
            "electronic_energy_hartree": product_basin.data.get(
                "electronic_energy_hartree"
            ),
            "calculation_artifact_ids": list(
                product_basin.data.get("calculation_artifact_ids") or []
            ),
        },
    }
    if not signatures_match or ambiguous:
        reasons = []
        if not signatures_match:
            reasons.append("candidate_pair_method_lineage_mismatch")
        if ambiguous:
            reasons.append("registry_basin_identity_ambiguous")
        return (
            {
                **common,
                **endpoint_evidence,
                "status": "unresolved",
                "accepted": False,
                "same_basin": False,
                "distinct_basin": False,
                "reasons": reasons,
            },
            parents,
        )

    return (
        {
            **common,
            **endpoint_evidence,
            "status": "distinct_basin",
            "accepted": True,
            "same_basin": False,
            "distinct_basin": True,
            "same_pes": True,
            "reasons": [],
        },
        parents,
    )


def _assessment_artifact(
    candidate: Artifact,
    assessment: dict[str, Any],
    parents: list[str],
) -> Artifact:
    candidate_id = str(candidate.data.get("candidate_id") or candidate.artifact_id)
    return Artifact(
        artifact_id=f"basin_pair_{path_token(candidate_id, max_length=48)}",
        artifact_type="basin_pair_assessment",
        parents=list(dict.fromkeys(parents)),
        data=assessment,
        qc={
            "registry_resolved": assessment.get("accepted") is True,
            "same_basin": assessment.get("same_basin") is True,
            "distinct_basin": assessment.get("distinct_basin") is True,
            "same_pes": assessment.get("same_pes") is True,
        },
        provenance={"created_by": "reaction-plan"},
        status=ArtifactStatus(
            status="success" if assessment.get("accepted") else "partial",
            category=(
                None
                if assessment.get("accepted")
                else "basin_identity_unresolved"
            ),
            reason="; ".join(assessment.get("reasons") or []),
        ),
    )


def _promote_candidate(
    candidate: Artifact,
    assessment: dict[str, Any],
    basins: dict[str, Artifact],
    minima: dict[str, Artifact],
) -> Artifact:
    reactant_basin = basins[str(assessment["reactant_basin_id"])]
    product_basin = basins[str(assessment["product_basin_id"])]
    reactant_species_id = str(reactant_basin.data["representative_species_id"])
    product_species_id = str(product_basin.data["representative_species_id"])
    reactant = minima[reactant_species_id]
    product = minima[product_species_id]
    states = {
        reactant_species_id: reactant,
        product_species_id: product,
        reactant.artifact_id: reactant,
        product.artifact_id: product,
    }
    candidate_id = str(candidate.data.get("candidate_id") or candidate.artifact_id)
    definition = {
        "reaction_id": make_reaction_id(candidate_id),
        "reactant": reactant_species_id,
        "product": product_species_id,
        "mechanism_family": candidate.data.get("mechanism_family", "discovered"),
        "reaction_type": candidate.data.get("mechanism_family", "discovered"),
        "bond_changes": list(candidate.data.get("bond_changes") or []),
        "reaction_coordinate": dict(
            candidate.data.get("reaction_coordinate") or {}
        ),
        "rationale": list(candidate.data.get("rationale") or []),
        "candidate_id": candidate_id,
    }
    reaction = build_reaction_hypothesis(definition, states)
    reaction.parents = list(
        dict.fromkeys(
            [
                candidate.artifact_id,
                reactant.artifact_id,
                product.artifact_id,
                reactant_basin.artifact_id,
                product_basin.artifact_id,
            ]
        )
    )
    reaction.data.update(
        {
            "reaction_trial_id": candidate.data.get("trial_id"),
            "discovery_evidence_artifact_id": candidate.data.get(
                "evidence_artifact_id"
            ),
            "reactant_basin_id": reactant_basin.artifact_id,
            "product_basin_id": product_basin.artifact_id,
            "basin_assessment": assessment,
        }
    )
    reaction.qc.update(
        {
            "dft_minima_validated": True,
            "distinct_registry_basins": True,
            "same_pes": True,
            "reaction_promoted": True,
        }
    )
    reaction.provenance["created_by"] = "reaction-plan"
    return reaction


class ReactionPlanStage(Stage):
    """Classify candidates and authorize path search only for true reactions."""

    name = "reaction-plan"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("reaction-plan requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        registry = _registry(manifest)
        basins = {
            artifact.artifact_id: artifact
            for artifact in manifest.latest_artifacts("minimum_basin")
        }
        minima = preferred_species_by_id(
            manifest,
            artifact_types=("species_optimized",),
            accept=is_accepted_optimized_minimum,
        )
        attempts = attempts_by_reaction(manifest)
        validated_ts = artifacts_by_reaction(manifest, "reaction_validated")
        validated_paths = artifacts_by_reaction(
            manifest, "reaction_path_validated"
        )
        existing_by_candidate = {
            str(reaction.data.get("candidate_id")): reaction
            for reaction in manifest.latest_artifacts("reaction")
            if reaction.data.get("candidate_id")
        }
        cases: list[dict[str, Any]] = []

        for candidate in manifest.latest_artifacts("reaction_candidate"):
            assessment, parents = _basin_assessment(candidate, registry, basins)
            candidate_id = str(
                candidate.data.get("candidate_id") or candidate.artifact_id
            )
            reaction: Artifact | None = existing_by_candidate.get(candidate_id)
            if assessment.get("distinct_basin") is True and reaction is None:
                try:
                    reaction = _promote_candidate(
                        candidate, assessment, basins, minima
                    )
                except (KeyError, ValueError) as exc:
                    assessment.update(
                        {
                            "status": "invalid_endpoints",
                            "accepted": False,
                            "distinct_basin": False,
                            "reasons": [f"reaction_promotion_failed:{exc}"],
                        }
                    )
                else:
                    out.add_artifact(reaction)
                    existing_by_candidate[candidate_id] = reaction
            assessment_artifact = _assessment_artifact(
                candidate, assessment, parents
            )
            out.add_artifact(assessment_artifact)

            case_reaction_id = (
                str(reaction.data.get("reaction_id") or reaction.artifact_id)
                if reaction is not None
                else candidate_id
            )
            path_attempts = attempts.get(case_reaction_id, [])
            ts_ok = any(
                item.qc.get("ts_validated_by_frequency") is True
                for item in validated_ts.get(case_reaction_id, [])
            )
            irc_ok = any(
                item.qc.get("irc_validated") is True
                for item in validated_paths.get(case_reaction_id, [])
            )
            case = decide_reaction_case(
                reaction_id=case_reaction_id,
                basin_assessment=assessment,
                path_attempts=path_attempts,
                ts_validated=ts_ok,
                irc_validated=irc_ok,
                max_path_attempts=int(config.get("max_path_attempts", 3)),
            )
            case.update(
                {
                    "candidate_id": candidate_id,
                    "reaction_promoted": reaction is not None,
                }
            )
            case_artifact = make_reaction_case_artifact(
                case,
                parents=[
                    assessment_artifact.artifact_id,
                    *(
                        [reaction.artifact_id]
                        if reaction is not None
                        else []
                    ),
                    *[item.artifact_id for item in path_attempts],
                ],
                created_by=self.name,
            )
            out.add_artifact(case_artifact)
            cases.append(case_artifact.model_dump(mode="json"))

        write_jsonl(cases, out_dir / "reaction_cases.jsonl")
        out.metadata.update(
            {
                "reaction_candidate_count": len(
                    manifest.latest_artifacts("reaction_candidate")
                ),
                "reaction_case_count": len(cases),
                "reaction_promoted_count": sum(
                    bool(row["data"].get("reaction_promoted"))
                    for row in cases
                ),
                "minimum_registry_present": registry is not None,
            }
        )
        return out
