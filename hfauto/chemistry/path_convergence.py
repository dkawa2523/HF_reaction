"""Discretization-convergence assessment for molecular reaction paths."""

from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy as np

from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.path import ReactionPathRecord


def path_method_signature(path_artifact: Artifact) -> dict[str, Any]:
    """Return method fields that must remain fixed across path resolutions."""

    method = path_artifact.method or {}
    settings = method.get("settings") or {}
    return {
        "engine": method.get("engine"),
        "backend": method.get("backend"),
        "functional": method.get("functional"),
        "basis": method.get("basis"),
        "disp_vdw": method.get("disp_vdw"),
        "program_version": path_artifact.data.get("program_version"),
        "charge": method.get("charge"),
        "multiplicity": method.get("multiplicity"),
        "grid": settings.get("grid", "fine"),
        "scf_energy_tolerance": float(
            settings.get("scf_energy_tolerance", 1.0e-7)
        ),
        "path_tolerance": settings.get(
            "string_tolerance", settings.get("neb_tolerance")
        ),
        "path_stepsize": settings.get(
            "string_stepsize", settings.get("neb_stepsize")
        ),
        "path_maxiter": settings.get(
            "string_maxiter", settings.get("neb_maxiter")
        ),
        "path_history": settings.get("string_nhist"),
        "path_interpolation": settings.get("string_interpol"),
    }


def _profile_kcal_mol(
    member: Artifact,
    path: ReactionPathRecord,
    points: int = 201,
) -> np.ndarray:
    energies = np.asarray(
        [image.energy_hartree for image in path.images], dtype=float
    )
    relative = (energies - energies[0]) * HARTREE_TO_KCAL_MOL
    source = np.asarray(
        member.data.get("normalized_path_coordinate"), dtype=float
    )
    if (
        source.shape != relative.shape
        or not np.isfinite(source).all()
        or abs(float(source[0])) > 1.0e-12
        or abs(float(source[-1]) - 1.0) > 1.0e-12
        or np.any(np.diff(source) <= 0.0)
    ):
        raise ValueError("validated normalized path coordinate is invalid")
    return np.interp(np.linspace(0.0, 1.0, points), source, relative)


def _member_record(member: Artifact) -> ReactionPathRecord | None:
    try:
        return ReactionPathRecord.model_validate(member.data.get("path"))
    except (TypeError, ValueError):
        return None


def _has_valid_path_coordinate(
    member: Artifact, path: ReactionPathRecord | None
) -> bool:
    if path is None:
        return False
    try:
        coordinate = np.asarray(
            member.data.get("normalized_path_coordinate"), dtype=float
        )
    except (TypeError, ValueError):
        return False
    return bool(
        len(path.images) >= 2
        and coordinate.shape == (len(path.images),)
        and np.isfinite(coordinate).all()
        and abs(float(coordinate[0])) <= 1.0e-12
        and abs(float(coordinate[-1]) - 1.0) <= 1.0e-12
        and np.all(np.diff(coordinate) > 0.0)
    )


def assess_path_ensemble(
    members: list[Artifact],
    *,
    minimum_resolutions: int = 2,
    profile_energy_tolerance_kcal_mol: float = 0.05,
) -> dict[str, Any]:
    """Assess image-discretization convergence without inventing path accuracy."""

    if minimum_resolutions < 2:
        raise ValueError("minimum_resolutions must be at least two")
    if profile_energy_tolerance_kcal_mol <= 0.0:
        raise ValueError("profile_energy_tolerance_kcal_mol must be positive")

    eligible: list[tuple[Artifact, ReactionPathRecord]] = []
    rejected: dict[str, list[str]] = {}
    for member in members:
        path = _member_record(member)
        reasons: list[str] = []
        if member.status.status != "success":
            reasons.append("ensemble_member_not_successful")
        if member.qc.get("real_path_executed") is not True:
            reasons.append("real_path_not_executed")
        if member.qc.get("path_geometry_validated") is not True:
            reasons.append("path_geometry_not_validated")
        if member.qc.get("method_evidence_validated") is not True:
            reasons.append("method_evidence_not_validated")
        if member.qc.get("input_method_evidence_validated") is not True:
            reasons.append("rendered_input_method_not_validated")
        if member.qc.get("endpoint_method_lineage_validated") is not True:
            reasons.append("endpoint_path_method_lineage_not_validated")
        if (
            member.qc.get("path_coordinate_validated") is not True
            or not _has_valid_path_coordinate(member, path)
        ):
            reasons.append("path_coordinate_not_validated")
        if member.qc.get("fallback_dummy") is not False:
            reasons.append("dummy_or_unknown_path_backend")
        if path is None:
            reasons.append("normalized_path_missing")
        elif not path.converged:
            reasons.append("path_not_converged")
        elif len(path.images) < 3 or any(
            image.energy_hartree is None for image in path.images
        ):
            reasons.append("complete_image_energy_profile_missing")
        if reasons:
            rejected[member.artifact_id] = reasons
        else:
            eligible.append((member, path))

    signatures = {
        repr(sorted((member.data.get("method_signature") or {}).items()))
        for member, _path in eligible
    }
    image_counts = {len(path.images) for _member, path in eligible}
    input_hashes = {
        str(member.data.get("input_sha256") or "")
        for member, _path in eligible
    }
    output_hashes = {
        str(member.data.get("output_sha256") or "")
        for member, _path in eligible
    }
    path_hashes = {
        str(member.data.get("path_sha256") or "")
        for member, _path in eligible
    }
    pairwise: list[dict[str, Any]] = []
    for (left_member, left), (right_member, right) in combinations(eligible, 2):
        difference = np.abs(
            _profile_kcal_mol(left_member, left)
            - _profile_kcal_mol(right_member, right)
        )
        maximum_index = int(np.argmax(difference))
        deviation = float(difference[maximum_index])
        pairwise.append(
            {
                "left_member_id": left_member.artifact_id,
                "right_member_id": right_member.artifact_id,
                "maximum_profile_deviation_kcal_mol": deviation,
                "maximum_deviation_normalized_string_coordinate": (
                    maximum_index / (len(difference) - 1)
                ),
                "within_tolerance": deviation
                <= float(profile_energy_tolerance_kcal_mol),
            }
        )

    classifications = {path.classification for _member, path in eligible}
    profile_converged = bool(pairwise and all(item["within_tolerance"] for item in pairwise))
    independent_discretizations = bool(
        len(image_counts) >= minimum_resolutions
        and len(input_hashes - {""}) >= minimum_resolutions
        and len(output_hashes - {""}) >= minimum_resolutions
        and len(path_hashes - {""}) >= minimum_resolutions
    )
    reasons: list[str] = []
    rejected_reason_set = {
        reason
        for member_reasons in rejected.values()
        for reason in member_reasons
    }
    if len(eligible) < minimum_resolutions:
        reasons.append("eligible_path_count_insufficient")
    if "endpoint_path_method_lineage_not_validated" in rejected_reason_set:
        reasons.append("endpoint_path_method_lineage_not_validated")
    if "rendered_input_method_not_validated" in rejected_reason_set:
        reasons.append("rendered_input_method_not_validated")
    if len(signatures) > 1:
        reasons.append("path_methods_or_numerical_settings_differ")
    if not independent_discretizations:
        reasons.append("independent_image_discretizations_insufficient")
    if len(classifications) > 1:
        reasons.append("path_profile_classifications_disagree")
    if len(eligible) >= minimum_resolutions and not profile_converged:
        reasons.append("path_energy_profiles_not_converged")

    consensus = next(iter(classifications)) if len(classifications) == 1 else None
    accepted = not reasons
    if accepted and consensus == "monotonic_no_internal_maximum":
        next_action = "classify_effectively_barrierless_at_resolution"
    elif accepted and consensus == "resolved_internal_maximum":
        next_action = "refine_first_order_saddle"
    elif "endpoint_path_method_lineage_not_validated" in reasons:
        next_action = "harmonize_endpoint_and_path_methods"
    elif "rendered_input_method_not_validated" in reasons:
        next_action = "repair_rendered_path_method_evidence"
    elif "path_methods_or_numerical_settings_differ" in reasons:
        next_action = "harmonize_path_methods_and_repeat"
    elif "path_energy_profiles_not_converged" in reasons:
        next_action = "refine_near_maximum_profile_disagreement"
    else:
        next_action = "complete_missing_path_ensemble_members"
    return {
        "schema_version": "hfauto.path_ensemble_assessment.v1",
        "accepted": accepted,
        "consensus_classification": consensus if accepted else None,
        "barrierless_at_resolution": bool(
            accepted and consensus == "monotonic_no_internal_maximum"
        ),
        "resolved_saddle_profile": bool(
            accepted and consensus == "resolved_internal_maximum"
        ),
        "eligible_member_ids": [member.artifact_id for member, _path in eligible],
        "rejected_members": rejected,
        "image_counts": sorted(image_counts),
        "method_signature_count": len(signatures),
        "profile_coordinate": (
            "normalized_cumulative_kabsch_aligned_cartesian_arc_length"
        ),
        "minimum_resolutions": int(minimum_resolutions),
        "profile_energy_tolerance_kcal_mol": float(
            profile_energy_tolerance_kcal_mol
        ),
        "pairwise_profile_comparisons": pairwise,
        "next_action": next_action,
        "reasons": reasons,
    }
