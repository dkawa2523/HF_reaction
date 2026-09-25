"""Split an evidenced multistep path into independently validated segments."""

from __future__ import annotations

from itertools import pairwise
from typing import Any

from hfauto.chemistry.reactions import endpoint_derived_coordinate_basis
from hfauto.core.artifacts import canonical_species_id, species_xyz_path
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_inputs import resolve_reaction_endpoints


def _coordinate_from_endpoints(
    reactant: Artifact,
    product: Artifact,
    *,
    minimum_distance_change_A: float,
    maximum_terms: int,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    evidence = endpoint_derived_coordinate_basis(
        {},
        species_xyz_path(reactant),
        species_xyz_path(product),
        minimum_distance_change_A=minimum_distance_change_A,
        maximum_terms=maximum_terms,
    )
    ranked = evidence.get("recommended_coordinate_terms") or []
    if evidence.get("accepted") is not True or not ranked:
        return None, evidence
    coordinate = {
        "terms": [
            {
                "kind": "distance",
                "atoms": term["atoms"],
                "coefficient": term["coefficient_for_positive_progress"],
                "label": term["interpretation"],
            }
            for term in ranked
        ],
        "min_change": min(
            0.05,
            max(1.0e-6, min(abs(float(term["delta_A"])) for term in ranked)),
        ),
    }
    return coordinate, evidence


class ReactionSegmentsStage(Stage):
    """Create child reaction hypotheses only from a validated third minimum."""

    name = "reaction-segments"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("reaction-segments requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        reactions = {
            str(item.data.get("reaction_id") or item.artifact_id): item
            for item in manifest.latest_artifacts("reaction")
        }
        records: list[dict[str, Any]] = []
        child_count = 0

        for classification in manifest.latest_artifacts(
            "reaction_classification"
        ):
            reaction_id = str(classification.data.get("reaction_id") or "")
            reaction = reactions.get(reaction_id)
            if reaction is None or reaction.data.get("parent_reaction_id"):
                continue
            if not (
                classification.status.status == "success"
                and classification.data.get("scientific_conclusion_supported") is True
                and classification.data.get("classification")
                == "multistep_with_intermediate"
            ):
                continue

            reactant, product, endpoint_gate = resolve_reaction_endpoints(
                manifest, reaction
            )
            intermediate_records = (
                (classification.data.get("evidence") or {}).get(
                    "validated_intermediates"
                )
                or []
            )
            intermediate_records = sorted(
                intermediate_records,
                key=lambda item: (
                    int(item.get("source_image_index", -1)),
                    str(item.get("artifact_id") or ""),
                ),
            )
            intermediates = [
                manifest.find(str(item.get("artifact_id") or ""))
                for item in intermediate_records
            ]
            reasons = list(endpoint_gate.get("reasons") or [])
            if reactant is None or product is None:
                reasons.append("validated_parent_endpoints_missing")
            if not intermediate_records:
                reasons.append("validated_intermediate_list_empty")
            if any(item is None for item in intermediates):
                reasons.append("validated_intermediate_artifact_missing")
            nodes = (
                [reactant, *intermediates, product]
                if not reasons
                else []
            )
            segment_specs: list[dict[str, Any]] = []
            for index, (left, right) in enumerate(pairwise(nodes)):
                assert left is not None and right is not None
                coordinate, coordinate_evidence = _coordinate_from_endpoints(
                    left,
                    right,
                    minimum_distance_change_A=float(
                        config.get("minimum_distance_change_A", 0.03)
                    ),
                    maximum_terms=int(config.get("maximum_coordinate_terms", 8)),
                )
                if coordinate is None:
                    reasons.append(
                        f"segment_{index}_independent_coordinate_missing"
                    )
                    continue
                segment_specs.append(
                    {
                        "index": index,
                        "reactant": left,
                        "product": right,
                        "reaction_coordinate": coordinate,
                        "coordinate_evidence": coordinate_evidence,
                    }
                )
            expected_count = max(0, len(nodes) - 1)
            ready = not reasons and len(segment_specs) == expected_count
            child_ids: list[str] = []
            if ready:
                for spec in segment_specs:
                    index = int(spec["index"])
                    left = spec["reactant"]
                    right = spec["product"]
                    child_reaction_id = (
                        f"{reaction_id}::segment_{index:02d}"
                    )
                    child = Artifact(
                        artifact_id="reaction_segment_"
                        + fingerprint_dict(
                            {
                                "parent_reaction_id": reaction_id,
                                "segment_index": index,
                                "reactant_artifact_id": left.artifact_id,
                                "product_artifact_id": right.artifact_id,
                            }
                        ),
                        artifact_type="reaction",
                        parents=[
                            reaction.artifact_id,
                            classification.artifact_id,
                            left.artifact_id,
                            right.artifact_id,
                        ],
                        data={
                            "reaction_id": child_reaction_id,
                            "parent_reaction_id": reaction_id,
                            "segment_index": index,
                            "segment_count": expected_count,
                            "reaction_type": "path_segment",
                            "mechanism_family": "molecular_elementary_segment",
                            "reactant_species_id": canonical_species_id(left),
                            "product_species_id": canonical_species_id(right),
                            "reaction_coordinate": spec["reaction_coordinate"],
                            "bond_changes": [],
                            "minimum_endpoint_rmsd_A": float(
                                config.get("minimum_endpoint_rmsd_A", 0.02)
                            ),
                            "source_classification_artifact_id": (
                                classification.artifact_id
                            ),
                            "source_path_attempt_ids": list(
                                (
                                    classification.data.get("evidence") or {}
                                ).get("multistep_path_attempt_ids")
                                or []
                            ),
                            "coordinate_evidence": spec[
                                "coordinate_evidence"
                            ],
                            "requires_independent_path_validation": True,
                        },
                        qc={
                            "validated_third_minimum_source": True,
                            "elementary_step_claimed": False,
                            "segment_ready_for_minimum_registry": True,
                        },
                        provenance={"created_by": self.name},
                    )
                    out.add_artifact(child)
                    child_ids.append(child.artifact_id)
                    child_count += 1

            plan = Artifact(
                artifact_id="reaction_segmentation_"
                + fingerprint_dict(
                    {
                        "reaction_id": reaction_id,
                        "classification_id": classification.artifact_id,
                        "intermediate_ids": [
                            item.get("artifact_id")
                            for item in intermediate_records
                        ],
                    }
                ),
                artifact_type="reaction_segmentation_plan",
                parents=[
                    reaction.artifact_id,
                    classification.artifact_id,
                    *[
                        str(item.get("artifact_id"))
                        for item in intermediate_records
                        if item.get("artifact_id")
                    ],
                ],
                data={
                    "reaction_id": reaction_id,
                    "classification_artifact_id": classification.artifact_id,
                    "segment_reaction_artifact_ids": child_ids,
                    "segment_count": len(child_ids),
                    "expected_segment_count": expected_count,
                    "ready": ready,
                    "reasons": list(dict.fromkeys(reasons)),
                    "next_stage": "dft-minima" if ready else None,
                },
                qc={
                    "validated_intermediate_count": len(intermediate_records),
                    "all_segments_have_independent_coordinates": ready,
                    "elementary_step_claimed": False,
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if ready else "partial",
                    category=None if ready else "reaction_segmentation_incomplete",
                    reason=None if ready else "; ".join(dict.fromkeys(reasons)),
                ),
            )
            out.add_artifact(plan)
            records.append(plan.model_dump(mode="json"))

        write_jsonl(records, out_dir / "reaction_segmentation_plans.jsonl")
        out.metadata["reaction_segmentation_plan_count"] = len(records)
        out.metadata["reaction_segment_count"] = child_count
        return out
