"""Assess the two frequency calculations produced by minimum-mode-follow."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.basin_identity import (
    assess_minimum_pair_basin_identity,
)
from hfauto.chemistry.method_lineage import evaluate_endpoint_pair_lineage
from hfauto.chemistry.minimum_recovery import (
    assess_minimum_mode_following_source,
)
from hfauto.chemistry.reaction_path_qc import estimate_reaction_mode_overlap
from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _calculation_for_minimum(
    manifest: Manifest,
    minimum: Artifact | None,
) -> Artifact | None:
    if minimum is None:
        return None
    if minimum.artifact_type != "species_optimized":
        return None
    calc_id = minimum.qc.get("dft_calc_id") or (
        minimum.data.get("dft_minima") or {}
    ).get("calc_id")
    calculation = manifest.find(str(calc_id)) if calc_id else None
    return (
        calculation
        if calculation is not None
        and calculation.artifact_type == "calculation"
        else None
    )


def _successor_lineage(
    seed: Artifact | None,
    minimum: Artifact | None,
    calculation: Artifact | None,
    source_calculation_id: str,
) -> dict[str, Any]:
    reasons: list[str] = []
    if seed is None:
        reasons.append("mode_following_seed_missing")
    if minimum is None:
        reasons.append("frequency_validated_successor_minimum_missing")
    if calculation is None:
        reasons.append("successor_minimum_calculation_missing")
    if seed is not None:
        if seed.data.get("mode_following_generation") != 1:
            reasons.append("mode_following_generation_is_not_one")
        if seed.data.get("source_stationary_point_calculation_id") != (
            source_calculation_id
        ):
            reasons.append("seed_source_stationary_point_mismatch")
    if minimum is not None and seed is not None:
        if minimum.artifact_type != "species_optimized":
            reasons.append("successor_is_not_an_optimized_species")
        if seed.artifact_id not in minimum.parents:
            reasons.append("successor_minimum_not_descended_from_seed")
        if minimum.data.get("mode_following_generation") != 1:
            reasons.append("successor_mode_following_generation_is_not_one")
        if minimum.data.get(
            "source_stationary_point_calculation_id"
        ) != source_calculation_id:
            reasons.append("successor_source_stationary_point_mismatch")
        if minimum.qc.get("minimum_accepted") is not True:
            reasons.append("successor_minimum_not_accepted")
        if minimum.qc.get("n_imag") != 0:
            reasons.append("successor_minimum_does_not_have_zero_imaginary_modes")
        if minimum.qc.get("real_qm_executed") is not True:
            reasons.append("successor_real_qm_not_executed")
        if minimum.qc.get("state_method_evidence_validated") is not True:
            reasons.append("successor_method_evidence_not_validated")
    if calculation is not None and seed is not None:
        if calculation.status.status != "success":
            reasons.append("successor_calculation_not_successful")
        if calculation.qc.get("minimum_accepted") is not True:
            reasons.append("successor_calculation_minimum_not_accepted")
        if calculation.qc.get("real_qm_executed") is not True:
            reasons.append("successor_calculation_real_qm_not_executed")
        if calculation.data.get("n_imag") != 0:
            reasons.append("successor_calculation_n_imag_is_not_zero")
        if calculation.data.get("source_geometry_artifact_id") != seed.artifact_id:
            reasons.append("successor_calculation_source_seed_mismatch")
        try:
            observed_seed_hash = sha256_file(seed.paths.get("xyz", ""))
        except (OSError, TypeError, ValueError):
            observed_seed_hash = None
        declared_seed_hash = (seed.data.get("mode_following") or {}).get(
            "seed_xyz_sha256"
        )
        if not declared_seed_hash or observed_seed_hash != declared_seed_hash:
            reasons.append("mode_following_seed_hash_mismatch")
        if calculation.data.get("input_xyz_sha256") != observed_seed_hash:
            reasons.append("successor_calculation_input_seed_hash_mismatch")
        if calculation.data.get("frequency_count_complete") is not True:
            reasons.append("successor_frequency_count_not_complete")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "seed_artifact_id": seed.artifact_id if seed else None,
        "minimum_artifact_id": minimum.artifact_id if minimum else None,
        "calculation_artifact_id": calculation.artifact_id if calculation else None,
    }


def _stationary_point_relative_energies(
    source: Artifact | None,
    calculations: dict[str, Artifact | None],
) -> dict[str, Any]:
    source_energy = (
        source.data.get("electronic_energy_hartree") if source else None
    )
    result: dict[str, Any] = {
        "source_electronic_energy_hartree": source_energy,
        "interpretation": (
            "candidate electronic energy differences only; not activation "
            "barriers until TS and bidirectional IRC validation"
        ),
    }
    for direction, calculation in calculations.items():
        successor_energy = (
            calculation.data.get("electronic_energy_hartree")
            if calculation is not None
            else None
        )
        result[direction] = {
            "successor_electronic_energy_hartree": successor_energy,
            "stationary_point_minus_successor_kcal_mol": (
                (float(source_energy) - float(successor_energy))
                * HARTREE_TO_KCAL_MOL
                if source_energy is not None and successor_energy is not None
                else None
            ),
        }
    return result


def _declared_reaction_mode_projections(
    manifest: Manifest,
    reaction_ids: list[str],
    source: Artifact | None,
    *,
    threshold: float,
) -> list[dict[str, Any]]:
    if source is None:
        return []
    output_path = source.paths.get("output")
    try:
        output_text = (
            Path(output_path).read_text(encoding="utf-8", errors="ignore")
            if output_path
            else None
        )
    except OSError:
        output_text = None
    projections: list[dict[str, Any]] = []
    for reaction_id in reaction_ids:
        reaction = manifest.find(reaction_id)
        if reaction is None or reaction.artifact_type != "reaction":
            projections.append(
                {
                    "reaction_id": reaction_id,
                    "accepted": False,
                    "reasons": ["declared_reaction_artifact_missing"],
                }
            )
            continue
        try:
            score, method = estimate_reaction_mode_overlap(
                reaction.data,
                source.paths.get("final_xyz"),
                source.data.get("n_imag"),
                source.data.get("imag_freq_cm1"),
                output_text,
                mode_displacements=source.data.get(
                    "imaginary_mode_displacements"
                ),
                mode_component_units=(
                    source.data.get("projected_imaginary_mode") or {}
                ).get("component_units"),
                allow_legacy_geometry_heuristic=False,
                imaginary_frequency_cutoff_cm1=float(
                    source.data.get("imaginary_frequency_cutoff_cm1", -1.0)
                ),
            )
        except (OSError, TypeError, ValueError):
            score, method = None, "unavailable"
        accepted = bool(score is not None and float(score) >= threshold)
        projections.append(
            {
                "reaction_id": reaction_id,
                "accepted": accepted,
                "mode_overlap_score": score,
                "mode_overlap_method": method,
                "minimum_overlap_threshold": threshold,
                "reasons": (
                    []
                    if accepted
                    else ["imaginary_mode_does_not_match_declared_reaction_coordinate"]
                ),
            }
        )
    return projections


class MinimumModeAssessStage(Stage):
    """Resolve only the local basin outcome; TS/IRC remain separate gates."""

    name = "minimum-mode-assess"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        species_by_id = preferred_species_by_id(
            manifest,
            artifact_types=("species_optimized",),
        )
        records: list[dict[str, Any]] = []

        for plan in manifest.latest_artifacts("minimum_mode_following_plan"):
            if plan.data.get("accepted") is not True:
                continue
            source_calculation_id = str(
                plan.data.get("source_calculation_artifact_id") or ""
            )
            source_calculation = manifest.find(source_calculation_id)
            source_reaudit = (
                assess_minimum_mode_following_source(source_calculation)
                if source_calculation is not None
                and source_calculation.artifact_type == "calculation"
                else {
                    "accepted": False,
                    "reasons": ["source_stationary_point_calculation_missing"],
                }
            )
            planned_source = plan.data.get("source_eligibility")
            if not isinstance(planned_source, dict):
                planned_source = {}
            if (
                source_reaudit.get("source_file_sha256")
                != planned_source.get("source_file_sha256")
                or source_reaudit.get("imaginary_mode_fingerprint")
                != planned_source.get("imaginary_mode_fingerprint")
            ):
                source_reaudit = {
                    **source_reaudit,
                    "accepted": False,
                    "reasons": [
                        *source_reaudit.get("reasons", []),
                        "source_evidence_changed_after_seed_generation",
                    ],
                }
            seed_ids = [str(value) for value in plan.data.get("seed_artifact_ids", [])]
            seeds = {identifier: manifest.find(identifier) for identifier in seed_ids}
            directions: dict[str, str] = {}
            for identifier, seed in seeds.items():
                if seed is not None:
                    direction = str(seed.data.get("mode_following_direction") or "")
                    if direction in {"plus", "minus"}:
                        directions[direction] = identifier

            minima = {
                direction: species_by_id.get(identifier)
                for direction, identifier in directions.items()
            }
            calculations = {
                direction: _calculation_for_minimum(manifest, minimum)
                for direction, minimum in minima.items()
            }
            lineages = {
                direction: _successor_lineage(
                    seeds.get(identifier),
                    minima.get(direction),
                    calculations.get(direction),
                    source_calculation_id,
                )
                for direction, identifier in directions.items()
            }
            for direction in ("plus", "minus"):
                if direction not in lineages:
                    lineages[direction] = {
                        "accepted": False,
                        "reasons": [f"{direction}_direction_seed_missing"],
                    }

            plus = minima.get("plus")
            minus = minima.get("minus")
            if plus is not None and minus is not None:
                basin = assess_minimum_pair_basin_identity(
                    plus,
                    minus,
                    calculations.get("plus"),
                    calculations.get("minus"),
                    same_basin_distance_rmsd_A=float(
                        config.get("same_basin_distance_rmsd_A", 0.02)
                    ),
                    distinct_basin_distance_rmsd_A=float(
                        config.get("distinct_basin_distance_rmsd_A", 0.05)
                    ),
                    same_basin_energy_tolerance_hartree=float(
                        config.get(
                            "same_basin_energy_tolerance_hartree", 1.0e-5
                        )
                    ),
                )
                method_lineage = evaluate_endpoint_pair_lineage(plus, minus)
            else:
                basin = {
                    "status": "incomplete",
                    "accepted": False,
                    "reasons": ["both_successor_minima_not_available"],
                }
                method_lineage = {
                    "accepted": False,
                    "reasons": ["both_successor_minima_not_available"],
                }

            local_supported = bool(
                source_reaudit.get("accepted") is True
                and all(value.get("accepted") is True for value in lineages.values())
                and basin.get("accepted") is True
                and method_lineage.get("accepted") is True
            )
            reaction_ids = [
                str(value) for value in plan.data.get("reaction_ids", [])
            ]
            reaction_mode_projections = _declared_reaction_mode_projections(
                manifest,
                reaction_ids,
                source_calculation,
                threshold=float(
                    config.get("reaction_mode_overlap_threshold", 0.70)
                ),
            )
            reaction_mode_supported = bool(
                reaction_mode_projections
                and any(
                    projection.get("accepted") is True
                    for projection in reaction_mode_projections
                )
            )
            if local_supported and basin.get("same_basin") is True:
                outcome = "same_basin_descents"
                next_action = "classify_declared_endpoint_change_as_same_basin_if_applicable"
            elif (
                local_supported
                and basin.get("distinct_basin") is True
                and reaction_mode_supported
            ):
                outcome = "distinct_basin_descents"
                next_action = "validate_source_as_ts_and_run_bidirectional_irc"
            elif local_supported and basin.get("distinct_basin") is True:
                outcome = "distinct_basin_descents_off_declared_coordinate"
                next_action = "revise_declared_reaction_coordinate_or_endpoint_hypothesis"
            else:
                outcome = "unresolved"
                next_action = "complete_or_repair_both_successor_minimum_calculations"

            reasons = [
                *source_reaudit.get("reasons", []),
                *[
                    f"{direction}_{reason}"
                    for direction, lineage in lineages.items()
                    for reason in lineage.get("reasons", [])
                ],
                *basin.get("reasons", []),
                *method_lineage.get("reasons", []),
            ]
            token = fingerprint_dict(
                {
                    "plan_id": plan.artifact_id,
                    "successor_calculation_ids": {
                        direction: calculation.artifact_id if calculation else None
                        for direction, calculation in calculations.items()
                    },
                }
            )
            assessment = Artifact(
                artifact_id=(
                    f"minimum_mode_assessment_{path_token(plan.artifact_id, 48)}_{token}"
                ),
                artifact_type="minimum_mode_following_assessment",
                parents=[
                    plan.artifact_id,
                    *[
                        minimum.artifact_id
                        for minimum in minima.values()
                        if minimum is not None
                    ],
                    *[
                        calculation.artifact_id
                        for calculation in calculations.values()
                        if calculation is not None
                    ],
                ],
                data={
                    "plan_artifact_id": plan.artifact_id,
                    "source_calculation_artifact_id": source_calculation_id,
                    "source_species_id": plan.data.get("source_species_id"),
                    "reaction_ids": plan.data.get("reaction_ids", []),
                    "outcome": outcome,
                    "local_basin_conclusion_supported": local_supported,
                    "reaction_classification_supported": False,
                    "source_stationary_point_reaudit": source_reaudit,
                    "successor_lineage": lineages,
                    "successor_basin_assessment": basin,
                    "successor_method_lineage": method_lineage,
                    "declared_reaction_mode_projections": (
                        reaction_mode_projections
                    ),
                    "declared_reaction_mode_supported": (
                        reaction_mode_supported
                    ),
                    "reaction_gate_reasons": list(
                        dict.fromkeys(
                            reason
                            for projection in reaction_mode_projections
                            for reason in projection.get("reasons", [])
                        )
                    ),
                    "stationary_point_relative_energies": (
                        _stationary_point_relative_energies(
                            source_calculation, calculations
                        )
                    ),
                    "next_action": next_action,
                    "reasons": list(dict.fromkeys(reasons)),
                },
                qc={
                    "local_basin_conclusion_supported": local_supported,
                    "source_promoted_to_transition_state": False,
                    "declared_reaction_mode_supported": reaction_mode_supported,
                    "bidirectional_irc_validated": False,
                    "reaction_classification_supported": False,
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if local_supported else "partial",
                    category=(
                        None if local_supported else "mode_following_outcome_unresolved"
                    ),
                    reason=None if local_supported else ";".join(dict.fromkeys(reasons)),
                ),
            )
            out.add_artifact(assessment)
            records.append(assessment.model_dump(mode="json"))

        write_jsonl(records, out_dir / "minimum_mode_following_assessments.jsonl")
        out.metadata["minimum_mode_following_assessment_count"] = len(records)
        out.metadata["minimum_mode_following_supported_count"] = sum(
            record["qc"]["local_basin_conclusion_supported"] for record in records
        )
        return out
