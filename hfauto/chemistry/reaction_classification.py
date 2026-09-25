"""Evidence-based classification of molecular reaction paths.

The module combines evidence; it does not run quantum chemistry.  Each
positive class has an explicit sufficient condition.  Incomplete or mutually
inconsistent evidence stays unresolved instead of being forced into a class.
"""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.basin_identity import match_geometry_to_basin
from hfauto.chemistry.method_lineage import evaluate_endpoint_pair_lineage
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.path import ReactionPathRecord

REACTION_CLASSES = {
    "elementary_first_order_saddle": "第一階鞍点を持つ素反応",
    "same_basin_relaxation": "同一basin内の緩和",
    "effectively_barrierless": "実質的に障壁なし（指定分解能内）",
    "multistep_with_intermediate": "別中間体を含む多段階反応",
    "unresolved": "未確定",
}


def _positive_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number > 0.0 else None


def _reaction_id(artifact: Artifact) -> str | None:
    value = artifact.data.get("reaction_id")
    return str(value) if value else None


def _path_record(attempt: Artifact) -> ReactionPathRecord | None:
    value = attempt.data.get("path")
    if not isinstance(value, dict):
        return None
    try:
        return ReactionPathRecord.model_validate(value)
    except ValueError:
        return None


def validate_path_attempt_evidence(
    attempt: Artifact,
    path_artifacts: dict[str, Artifact],
    *,
    require_real_qm: bool,
) -> tuple[ReactionPathRecord | None, list[str]]:
    record = _path_record(attempt)
    reasons: list[str] = []
    if attempt.status.status != "success":
        reasons.append("path_attempt_not_successful")
    if record is None:
        reasons.append("normalized_path_record_missing")
        return None, reasons
    if not record.converged:
        reasons.append("path_not_converged")
    if record.classification == "missing_profile":
        reasons.append("complete_path_energy_profile_missing")

    parent = next(
        (path_artifacts[parent_id] for parent_id in attempt.parents if parent_id in path_artifacts),
        None,
    )
    if require_real_qm:
        if parent is None:
            reasons.append("reaction_path_parent_missing")
        else:
            real_path = bool(
                parent.qc.get("real_path_executed") is True
                or parent.qc.get("real_neb_executed") is True
            )
            if not real_path:
                reasons.append("real_quantum_path_not_executed")
            if parent.qc.get("path_geometry_validated") is not True:
                reasons.append("path_geometry_not_validated")
            if parent.qc.get("method_evidence_validated") is not True:
                reasons.append("path_method_evidence_not_validated")
            if parent.qc.get("fallback_dummy") is True:
                reasons.append("path_uses_dummy_fallback")
    return record, list(dict.fromkeys(reasons))


def _validated_elementary_step(
    ts_artifacts: list[Artifact],
    irc_artifacts: list[Artifact],
    *,
    require_real_qm: bool,
) -> dict[str, Any]:
    valid_ts = [
        artifact
        for artifact in ts_artifacts
        if artifact.status.status == "success"
        and artifact.qc.get("ts_validated_by_frequency") is True
        and artifact.qc.get("n_imag") == 1
        and (
            not require_real_qm
            or (
                artifact.qc.get("real_ts_search_executed") is True
                and artifact.qc.get("real_qm_executed") is True
                and artifact.qc.get("fallback_dummy") is False
                and artifact.qc.get("method_evidence_validated") is True
            )
        )
    ]
    valid_irc = [
        artifact
        for artifact in irc_artifacts
        if artifact.status.status == "success"
        and artifact.qc.get("irc_validated") is True
        and (
            not require_real_qm
            or (
                artifact.qc.get("real_irc_executed") is True
                and artifact.qc.get("fallback_dummy") is False
            )
        )
    ]
    connected_pairs = [
        {"ts_artifact_id": ts.artifact_id, "irc_artifact_id": irc.artifact_id}
        for ts in valid_ts
        for irc in valid_irc
        if ts.artifact_id in {str(parent) for parent in irc.parents}
    ]
    return {
        "accepted": bool(connected_pairs),
        "ts_artifact_ids": [item.artifact_id for item in valid_ts],
        "irc_artifact_ids": [item.artifact_id for item in valid_irc],
        "connected_ts_irc_pairs": connected_pairs,
        "reasons": [
            *([] if valid_ts else ["frequency_validated_first_order_saddle_missing"]),
            *([] if valid_irc else ["bidirectional_irc_endpoint_connection_missing"]),
            *(
                []
                if connected_pairs or not (valid_ts and valid_irc)
                else ["irc_not_bound_to_frequency_validated_ts"]
            ),
        ],
    }


def _validated_intermediates(
    candidates: list[Artifact],
    attempts: dict[str, Artifact],
    path_artifacts: dict[str, Artifact],
    reactant: Artifact | None,
    product: Artifact | None,
    *,
    require_real_qm: bool,
    basin_match_rmsd_A: float,
) -> list[dict[str, Any]]:
    if reactant is None or product is None:
        return []
    try:
        reactant_xyz = species_xyz_path(reactant)
        product_xyz = species_xyz_path(product)
    except (OSError, TypeError, ValueError):
        return []

    accepted: list[dict[str, Any]] = []
    endpoint_state = (
        int(reactant.data.get("resolved_charge", reactant.data.get("charge", 0)) or 0),
        int(reactant.data.get("resolved_multiplicity", reactant.data.get("multiplicity", 1)) or 1),
    )
    for candidate in candidates:
        source_attempt_id = str(candidate.data.get("source_path_attempt_id") or "")
        attempt = attempts.get(source_attempt_id)
        if attempt is None:
            continue
        record, path_reasons = validate_path_attempt_evidence(
            attempt, path_artifacts, require_real_qm=require_real_qm
        )
        image_index = candidate.data.get("source_image_index")
        well_indices = {
            well.image_index
            for well in (record.intermediate_wells if record is not None else [])
            if well.resolved
        }
        state = (
            int(candidate.data.get("resolved_charge", candidate.data.get("charge", 0)) or 0),
            int(candidate.data.get("resolved_multiplicity", candidate.data.get("multiplicity", 1)) or 1),
        )
        reasons = list(path_reasons)
        if candidate.status.status != "success":
            reasons.append("intermediate_optimization_not_successful")
        if candidate.artifact_type != "species_optimized":
            reasons.append("intermediate_not_optimized")
        if candidate.qc.get("minimum_accepted") is not True:
            reasons.append("intermediate_minimum_not_accepted")
        if candidate.qc.get("is_minimum") is not True or candidate.qc.get("n_imag") != 0:
            reasons.append("intermediate_zero_imaginary_frequency_evidence_missing")
        if require_real_qm and candidate.qc.get("real_qm_executed") is not True:
            reasons.append("intermediate_real_quantum_optimization_missing")
        if (
            require_real_qm
            and candidate.qc.get("state_method_evidence_validated") is not True
        ):
            reasons.append("intermediate_method_evidence_not_validated")
        if candidate.qc.get("fallback_dummy") is True:
            reasons.append("intermediate_uses_dummy_fallback")
        if state != endpoint_state:
            reasons.append("intermediate_electronic_state_mismatch")
        if image_index not in well_indices:
            reasons.append("intermediate_not_linked_to_resolved_path_well")
        try:
            candidate_xyz = species_xyz_path(candidate)
            reactant_match = match_geometry_to_basin(
                candidate_xyz,
                reactant_xyz,
                ordered_rmsd_threshold_A=basin_match_rmsd_A,
                permutation_rmsd_threshold_A=basin_match_rmsd_A,
            )
            product_match = match_geometry_to_basin(
                candidate_xyz,
                product_xyz,
                ordered_rmsd_threshold_A=basin_match_rmsd_A,
                permutation_rmsd_threshold_A=basin_match_rmsd_A,
            )
        except (OSError, TypeError, ValueError):
            reasons.append("intermediate_geometry_unreadable")
            reactant_match = {"accepted": False}
            product_match = {"accepted": False}
        if reactant_match.get("accepted") is True:
            reasons.append("intermediate_relaxed_to_reactant_basin")
        if product_match.get("accepted") is True:
            reasons.append("intermediate_relaxed_to_product_basin")
        if not reasons:
            accepted.append(
                {
                    "artifact_id": candidate.artifact_id,
                    "source_path_attempt_id": source_attempt_id,
                    "source_image_index": image_index,
                    "reactant_basin_match": reactant_match,
                    "product_basin_match": product_match,
                }
            )
    return accepted


def classify_reaction(
    reaction: Artifact,
    basin_assessment: dict[str, Any],
    *,
    path_attempts: list[Artifact] = (),
    reaction_paths: list[Artifact] = (),
    ts_artifacts: list[Artifact] = (),
    irc_artifacts: list[Artifact] = (),
    intermediate_minima: list[Artifact] = (),
    path_ensemble_assessments: list[Artifact] = (),
    minimum_mode_following_assessments: list[Artifact] = (),
    reactant: Artifact | None = None,
    product: Artifact | None = None,
    require_real_qm: bool = True,
    barrier_resolution_kcal_mol: float = 0.05,
    intermediate_basin_match_rmsd_A: float = 0.20,
) -> dict[str, Any]:
    """Return one supported reaction class or an auditable unresolved result."""

    reaction_id = str(reaction.data.get("reaction_id") or reaction.artifact_id)
    path_artifacts = {artifact.artifact_id: artifact for artifact in reaction_paths}
    attempts = {artifact.artifact_id: artifact for artifact in path_attempts}
    paths: list[tuple[Artifact, ReactionPathRecord]] = []
    path_rejections: dict[str, list[str]] = {}
    for attempt in path_attempts:
        record, reasons = validate_path_attempt_evidence(
            attempt, path_artifacts, require_real_qm=require_real_qm
        )
        if record is not None and not reasons:
            paths.append((attempt, record))
        else:
            path_rejections[attempt.artifact_id] = reasons

    endpoint_evidence = [
        basin_assessment.get("reactant_evidence") or {},
        basin_assessment.get("product_evidence") or {},
    ]
    real_endpoint_pair = bool(
        all(item.get("real_qm_executed") is True for item in endpoint_evidence)
        and all(item.get("fallback_dummy") is False for item in endpoint_evidence)
        and all(item.get("method_evidence_validated") is True for item in endpoint_evidence)
    )
    endpoint_pair_lineage = (
        evaluate_endpoint_pair_lineage(reactant, product)
        if reactant is not None and product is not None
        else {
            "accepted": False,
            "reasons": ["endpoint_pair_missing"],
        }
    )
    endpoint_species_ids = {
        str(endpoint.data.get("species_id") or "")
        for endpoint in (reactant, product)
        if endpoint is not None
    }
    active_mode_assessments = [
        assessment
        for assessment in minimum_mode_following_assessments
        if not endpoint_species_ids
        or str(assessment.data.get("source_species_id") or "")
        in endpoint_species_ids
    ]
    latest_mode_assessment = (
        active_mode_assessments[-1] if active_mode_assessments else None
    )
    mode_following_evidence = (
        dict(latest_mode_assessment.data)
        if latest_mode_assessment is not None
        else {}
    )
    mode_following_outcome = str(
        mode_following_evidence.get("outcome") or ""
    )
    all_declared_mode_projections = mode_following_evidence.get(
        "declared_reaction_mode_projections"
    )
    declared_mode_projections = [
        projection
        for projection in (
            all_declared_mode_projections
            if isinstance(all_declared_mode_projections, list)
            else []
        )
        if isinstance(projection, dict)
        and str(projection.get("reaction_id") or "") == reaction_id
    ]
    mode_following_evidence["classification_reaction_id"] = reaction_id
    mode_following_evidence["selected_reaction_mode_projections"] = (
        declared_mode_projections
    )
    declared_mode_mismatch = bool(
        declared_mode_projections
        and not any(
            projection.get("accepted") is True
            for projection in declared_mode_projections
        )
    )
    reaction_mode_reasons = list(
        dict.fromkeys(
            reason
            for projection in declared_mode_projections
            for reason in projection.get("reasons", [])
        )
    )
    off_declared_coordinate = bool(
        declared_mode_mismatch
        or (
            not declared_mode_projections
            and
            mode_following_evidence.get("local_basin_conclusion_supported") is True
            and mode_following_outcome
            == "distinct_basin_descents_off_declared_coordinate"
        )
    )
    basin_supported = bool(
        basin_assessment.get("accepted") is True
        and (
            not require_real_qm
            or (
                real_endpoint_pair
                and endpoint_pair_lineage["accepted"] is True
            )
        )
    )

    elementary = _validated_elementary_step(
        list(ts_artifacts), list(irc_artifacts), require_real_qm=require_real_qm
    )
    intermediates = _validated_intermediates(
        list(intermediate_minima),
        attempts,
        path_artifacts,
        reactant,
        product,
        require_real_qm=require_real_qm,
        basin_match_rmsd_A=intermediate_basin_match_rmsd_A,
    )
    multi_attempt_ids = {
        attempt.artifact_id
        for attempt, record in paths
        if len(record.internal_maximum_indices) >= 2
        and any(well.resolved for well in record.intermediate_wells)
    }
    validated_intermediates = [
        item
        for item in intermediates
        if item["source_path_attempt_id"] in multi_attempt_ids
    ]

    latest_assessment = (
        path_ensemble_assessments[-1] if path_ensemble_assessments else None
    )
    latest_ensemble = (
        latest_assessment
        if latest_assessment is not None
        and latest_assessment.status.status == "success"
        and latest_assessment.data.get("accepted") is True
        else None
    )
    superseded_attempt_ids = {
        str(identifier)
        for identifier in (
            latest_ensemble.data.get("supersedes_path_attempt_ids", [])
            if latest_ensemble is not None
            else []
        )
    }
    active_paths = [
        (attempt, record)
        for attempt, record in paths
        if attempt.artifact_id not in superseded_attempt_ids
    ]
    monotonic_paths = [
        (attempt, record)
        for attempt, record in active_paths
        if record.classification == "monotonic_no_internal_maximum"
        and not record.internal_maximum_indices
    ]
    barrierless_paths = [
        (attempt, record)
        for attempt, record in monotonic_paths
        if record.barrier_threshold_kcal_mol
        <= float(barrier_resolution_kcal_mol)
    ]
    resolved_barrier_paths = [
        (attempt, record)
        for attempt, record in active_paths
        if record.classification == "resolved_internal_maximum"
        and record.internal_maximum_indices
    ]
    ensemble_barrierless = bool(
        latest_ensemble is not None
        and latest_ensemble.data.get("barrierless_at_resolution") is True
    )
    ensemble_path_resolution = _positive_float(
        (latest_ensemble.data or {}).get(
            "path_barrier_resolution_kcal_mol"
        )
        if latest_ensemble is not None
        else None
    )
    ensemble_profile_tolerance = _positive_float(
        (latest_ensemble.data or {}).get(
            "profile_energy_tolerance_kcal_mol"
        )
        if latest_ensemble is not None
        else None
    )
    ensemble_resolution_compatible = bool(
        ensemble_path_resolution is not None
        and ensemble_profile_tolerance is not None
        and ensemble_path_resolution <= float(barrier_resolution_kcal_mol)
        and ensemble_profile_tolerance <= float(barrier_resolution_kcal_mol)
    )
    barrierless_accepted = bool(
        basin_supported
        and basin_assessment.get("distinct_basin") is True
        and ensemble_barrierless
        and ensemble_resolution_compatible
        and not resolved_barrier_paths
        and not off_declared_coordinate
    )
    barrierless = {
        "accepted": barrierless_accepted,
        "scope": "electronic_energy_path_at_declared_resolution",
        "barrier_resolution_kcal_mol": float(barrier_resolution_kcal_mol),
        "ensemble_path_barrier_resolution_kcal_mol": (
            ensemble_path_resolution
        ),
        "ensemble_profile_tolerance_kcal_mol": (
            ensemble_profile_tolerance
        ),
        "ensemble_resolution_compatible": ensemble_resolution_compatible,
        "minimum_independent_resolutions": (
            (latest_ensemble.data or {}).get("minimum_resolutions")
            if latest_ensemble is not None
            else None
        ),
        "path_attempt_ids": [item.artifact_id for item, _record in barrierless_paths],
        "path_ensemble_assessment_id": (
            latest_assessment.artifact_id if latest_assessment else None
        ),
        "path_ensemble": latest_assessment.data if latest_assessment else None,
        "reasons": [
            *(
                []
                if basin_supported
                and basin_assessment.get("distinct_basin") is True
                else ["frequency_validated_distinct_endpoint_basins_missing"]
            ),
            *(
                []
                if ensemble_barrierless
                else ["accepted_path_ensemble_evidence_missing"]
            ),
            *(
                []
                if ensemble_resolution_compatible
                else [
                    "path_ensemble_resolution_coarser_than_classification_policy"
                ]
            ),
            *(
                []
                if len(barrierless_paths) == len(monotonic_paths)
                else ["path_energy_resolution_too_coarse"]
            ),
            *(
                []
                if not resolved_barrier_paths
                else ["converged_path_profiles_disagree"]
            ),
            *(
                []
                if not off_declared_coordinate
                else [
                    "stationary_point_mode_does_not_match_declared_reaction_coordinate"
                ]
            ),
        ],
        "conflicting_resolved_path_attempt_ids": [
            attempt.artifact_id for attempt, _record in resolved_barrier_paths
        ],
    }

    if basin_supported and basin_assessment.get("same_basin") is True:
        classification = "same_basin_relaxation"
        accepted = True
        next_action = "stop_same_basin_duplicate"
        reasons: list[str] = []
    elif (
        basin_supported
        and basin_assessment.get("distinct_basin") is True
        and elementary["accepted"]
    ):
        classification = "elementary_first_order_saddle"
        accepted = True
        next_action = "compute_activation_thermochemistry"
        reasons = []
    elif (
        basin_supported
        and basin_assessment.get("distinct_basin") is True
        and validated_intermediates
    ):
        classification = "multistep_with_intermediate"
        accepted = True
        next_action = "split_path_into_elementary_segments"
        reasons = []
    elif barrierless["accepted"]:
        classification = "effectively_barrierless"
        accepted = True
        next_action = "use_capture_or_variational_kinetics"
        reasons = []
    else:
        classification = "unresolved"
        accepted = False
        if off_declared_coordinate:
            candidate = "unresolved"
            next_action = "revise_declared_reaction_coordinate_or_endpoint_hypothesis"
            reasons = [
                *list(basin_assessment.get("reasons") or []),
                "stationary_point_mode_does_not_match_declared_reaction_coordinate",
                *reaction_mode_reasons,
            ]
        elif not basin_supported:
            if basin_assessment.get("same_basin") is True:
                candidate = "same_basin_relaxation"
            elif multi_attempt_ids:
                candidate = "multistep_with_intermediate"
            elif monotonic_paths or ensemble_barrierless:
                candidate = "effectively_barrierless"
            else:
                candidate = "elementary_first_order_saddle"
            endpoint_reasons = list(basin_assessment.get("reasons") or [])
            refinement_rejected = any(
                str(reason).endswith("latest_minimum_refinement_rejected")
                for reason in endpoint_reasons
            )
            next_action = (
                "follow_endpoint_imaginary_mode_in_both_directions"
                if refinement_rejected
                else "complete_endpoint_frequency_and_basin_identity"
            )
            reasons = [
                *endpoint_reasons,
                "frequency_validated_distinct_endpoint_basins_missing",
                *(
                    ["endpoint_pair_method_lineage_not_validated"]
                    if require_real_qm
                    and endpoint_pair_lineage["accepted"] is not True
                    else []
                ),
            ]
        elif multi_attempt_ids and not validated_intermediates:
            candidate = "multistep_with_intermediate"
            next_action = "optimize_and_frequency_check_path_well"
            reasons = ["resolved_path_well_has_no_distinct_frequency_validated_minimum"]
        elif monotonic_paths or ensemble_barrierless:
            candidate = "effectively_barrierless"
            next_action = "repeat_path_at_independent_tighter_resolution"
            reasons = list(barrierless["reasons"])
        elif basin_assessment.get("same_basin") is True:
            candidate = "same_basin_relaxation"
            next_action = "complete_endpoint_frequency_and_basin_identity"
            reasons = ["same_basin_evidence_not_accepted"]
        else:
            candidate = "elementary_first_order_saddle"
            next_action = "continue_saddle_and_irc_validation"
            reasons = list(elementary["reasons"])
        reasons = list(dict.fromkeys(reasons))
        return {
            "schema_version": "hfauto.reaction_classification.v1",
            "reaction_id": reaction_id,
            "classification": classification,
            "classification_label_ja": REACTION_CLASSES[classification],
            "candidate_classification": candidate,
            "candidate_label_ja": REACTION_CLASSES[candidate],
            "scientific_conclusion_supported": False,
            "next_action": next_action,
            "reasons": reasons,
            "evidence": {
                "basin_assessment": basin_assessment,
                "real_endpoint_pair_validated": real_endpoint_pair,
                "endpoint_pair_lineage": endpoint_pair_lineage,
                "elementary_step": elementary,
                "validated_intermediates": intermediates,
                "multistep_path_attempt_ids": sorted(multi_attempt_ids),
                "superseded_path_attempt_ids": sorted(superseded_attempt_ids),
                "barrierless": barrierless,
                "minimum_mode_following": mode_following_evidence,
                "rejected_paths": path_rejections,
            },
        }

    return {
        "schema_version": "hfauto.reaction_classification.v1",
        "reaction_id": reaction_id,
        "classification": classification,
        "classification_label_ja": REACTION_CLASSES[classification],
        "candidate_classification": classification,
        "candidate_label_ja": REACTION_CLASSES[classification],
        "scientific_conclusion_supported": accepted,
        "next_action": next_action,
        "reasons": reasons,
        "evidence": {
            "basin_assessment": basin_assessment,
            "real_endpoint_pair_validated": real_endpoint_pair,
            "endpoint_pair_lineage": endpoint_pair_lineage,
            "elementary_step": elementary,
            "validated_intermediates": validated_intermediates,
            "multistep_path_attempt_ids": sorted(multi_attempt_ids),
            "superseded_path_attempt_ids": sorted(superseded_attempt_ids),
            "barrierless": barrierless,
            "minimum_mode_following": mode_following_evidence,
            "rejected_paths": path_rejections,
        },
    }
