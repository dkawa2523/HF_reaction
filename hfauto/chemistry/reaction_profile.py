"""Method-neutral analysis of one-dimensional reaction-path energies."""

from __future__ import annotations

import re
from itertools import pairwise
from typing import Any

from hfauto.chemistry.xyz import XYZ
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.schemas.path import (
    IntermediateWellRecord,
    PathImageRecord,
    ReactionPathRecord,
)

HARTREE_TO_KJ_MOL = 2625.4996394799

_ENERGY_LABEL = re.compile(
    r"(?:^|[\s,;])(?:energy(?:_hartree)?|e)\s*(?:=|:)?\s*"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)",
    re.IGNORECASE,
)


def _internal_peak_candidates(
    energies: list[float], *, minimum_prominence_hartree: float
) -> tuple[list[int], list[dict[str, Any]]]:
    """Return strict internal maxima and basin-referenced peak evidence."""

    raw_maxima = [
        index
        for index in range(1, len(energies) - 1)
        if energies[index] > energies[index - 1]
        and energies[index] > energies[index + 1]
    ]
    candidates: list[dict[str, Any]] = []
    for index in raw_maxima:
        left = min(range(index + 1), key=energies.__getitem__)
        right = min(range(index, len(energies)), key=energies.__getitem__)
        left_rise = energies[index] - energies[left]
        right_rise = energies[index] - energies[right]
        prominence = min(left_rise, right_rise)
        candidates.append(
            {
                "point_index": index,
                "left_minimum_index": left,
                "right_minimum_index": right,
                "left_rise_hartree": left_rise,
                "right_rise_hartree": right_rise,
                "prominence_hartree": prominence,
                "resolved": prominence >= float(minimum_prominence_hartree),
            }
        )
    return raw_maxima, candidates


def find_intermediate_wells(
    energies: list[float],
    maximum_indices: list[int],
    *,
    minimum_depth_hartree: float,
) -> list[IntermediateWellRecord]:
    """Find resolved wells between adjacent path maxima.

    A geometric image is only a *candidate* intermediate.  The caller must
    still optimize it and confirm zero imaginary frequencies before assigning
    a distinct intermediate basin.
    """

    wells: list[IntermediateWellRecord] = []
    for left, right in pairwise(sorted(set(maximum_indices))):
        if right - left < 2:
            continue
        index = min(range(left + 1, right), key=energies.__getitem__)
        left_drop = energies[left] - energies[index]
        right_drop = energies[right] - energies[index]
        wells.append(
            IntermediateWellRecord(
                image_index=index,
                left_maximum_index=left,
                right_maximum_index=right,
                left_drop_kcal_mol=left_drop * HARTREE_TO_KCAL_MOL,
                right_drop_kcal_mol=right_drop * HARTREE_TO_KCAL_MOL,
                resolved=min(left_drop, right_drop)
                >= float(minimum_depth_hartree),
            )
        )
    return wells


def parse_energy_comment(comment: str) -> float | None:
    """Parse an explicitly labelled Hartree energy from an XYZ comment."""

    match = _ENERGY_LABEL.search(comment)
    return float(match.group(1)) if match else None


def select_path_ts_guess(path: str) -> tuple[XYZ, int, str]:
    """Select only an internal path image; prefer explicit energy evidence."""

    images = read_xyz_trajectory(path)
    if len(images) < 3:
        raise ValueError("reaction path must contain at least three images")
    internal = list(enumerate(images[1:-1], start=1))
    labelled = [
        (index, image, energy)
        for index, image in internal
        if (energy := parse_energy_comment(image.comment)) is not None
    ]
    if labelled:
        index, image, _energy = max(labelled, key=lambda item: item[2])
        return image, index, "highest_energy_internal_image"
    index, image = internal[len(internal) // 2]
    return image, index, "middle_internal_image_no_energy_labels"


def analyze_reaction_path(
    *,
    reaction_id: str,
    engine: str,
    comments: list[str],
    converged: bool,
    barrier_threshold_kcal_mol: float = 0.5,
    energy_source: str = "xyz_comment",
) -> ReactionPathRecord:
    """Normalize a path profile without inventing missing image energies."""

    images = [
        PathImageRecord(index=index, energy_hartree=parse_energy_comment(comment), comment=comment)
        for index, comment in enumerate(comments)
    ]
    energies = [image.energy_hartree for image in images]
    reasons: list[str] = []
    if len(images) < 3:
        reasons.append("fewer_than_three_path_images")
    if any(energy is None for energy in energies):
        reasons.append("image_energy_missing_or_unlabelled")
    if reasons:
        return ReactionPathRecord(
            reaction_id=reaction_id,
            engine=engine,
            converged=bool(converged),
            images=images,
            classification="missing_profile",
            barrier_threshold_kcal_mol=float(barrier_threshold_kcal_mol),
            energy_source=energy_source,
            unresolved_reasons=reasons,
        )

    values = [float(energy) for energy in energies if energy is not None]
    threshold_hartree = float(barrier_threshold_kcal_mol) / HARTREE_TO_KCAL_MOL
    internal_index = max(range(1, len(values) - 1), key=values.__getitem__)
    internal_energy = values[internal_index]
    _raw_maxima, peak_candidates = _internal_peak_candidates(
        values,
        minimum_prominence_hartree=threshold_hartree,
    )
    internal_maxima = [
        int(candidate["point_index"])
        for candidate in peak_candidates
        if candidate["resolved"]
    ]
    intermediate_wells = find_intermediate_wells(
        values,
        internal_maxima,
        minimum_depth_hartree=threshold_hartree,
    )
    resolved = internal_energy > max(values[0], values[-1]) + threshold_hartree
    deltas = [right - left for left, right in pairwise(values)]
    monotonic = all(delta >= -threshold_hartree for delta in deltas) or all(
        delta <= threshold_hartree for delta in deltas
    )
    classification = (
        "resolved_internal_maximum"
        if resolved
        else "monotonic_no_internal_maximum"
        if monotonic
        else "unresolved_internal_profile"
    )
    return ReactionPathRecord(
        reaction_id=reaction_id,
        engine=engine,
        converged=bool(converged),
        images=images,
        classification=classification,
        highest_internal_image_index=internal_index,
        internal_maximum_indices=internal_maxima,
        intermediate_wells=intermediate_wells,
        ts_guess_image_index=internal_index if resolved else None,
        barrier_from_reactant_kcal_mol=(
            (internal_energy - values[0]) * HARTREE_TO_KCAL_MOL
            if resolved
            else None
        ),
        maximum_internal_rise_from_reactant_kcal_mol=(
            internal_energy - values[0]
        )
        * HARTREE_TO_KCAL_MOL,
        reaction_energy_kcal_mol=(values[-1] - values[0]) * HARTREE_TO_KCAL_MOL,
        barrier_threshold_kcal_mol=float(barrier_threshold_kcal_mol),
        energy_source=energy_source,
    )


def reanalyze_reaction_path(
    path: ReactionPathRecord | dict[str, Any],
    *,
    barrier_threshold_kcal_mol: float,
) -> ReactionPathRecord:
    """Reapply one common peak threshold to stored image energies.

    This changes only profile interpretation.  It cannot improve the original
    electronic-structure or path-optimization accuracy.
    """

    record = (
        path
        if isinstance(path, ReactionPathRecord)
        else ReactionPathRecord.model_validate(path)
    )
    comments = [
        (
            f"energy_hartree={image.energy_hartree:.15f}"
            if image.energy_hartree is not None
            else image.comment or ""
        )
        for image in record.images
    ]
    return analyze_reaction_path(
        reaction_id=record.reaction_id,
        engine=record.engine,
        comments=comments,
        converged=record.converged,
        barrier_threshold_kcal_mol=barrier_threshold_kcal_mol,
        energy_source=f"reanalyzed:{record.energy_source}",
    )


def classify_energy_profile(
    energies: list[float], min_peak_prominence_kj_mol: float = 2.0
) -> dict[str, Any]:
    """Classify a Hartree profile without promoting noise or endpoints to TS."""

    if not energies:
        return {
            "profile_classification": "missing_profile",
            "maximum_point_index": None,
            "minimum_point_index": None,
            "internal_maximum_indices": [],
            "raw_internal_maximum_indices": [],
            "internal_peak_candidates": [],
        }
    maximum = max(range(len(energies)), key=energies.__getitem__)
    minimum = min(range(len(energies)), key=energies.__getitem__)
    raw_maxima, candidates_hartree = _internal_peak_candidates(
        energies,
        minimum_prominence_hartree=float(min_peak_prominence_kj_mol)
        / HARTREE_TO_KJ_MOL,
    )
    candidates = [
        {
            "point_index": candidate["point_index"],
            "left_minimum_index": candidate["left_minimum_index"],
            "right_minimum_index": candidate["right_minimum_index"],
            "left_rise_kj_mol": candidate["left_rise_hartree"]
            * HARTREE_TO_KJ_MOL,
            "right_rise_kj_mol": candidate["right_rise_hartree"]
            * HARTREE_TO_KJ_MOL,
            "prominence_kj_mol": candidate["prominence_hartree"]
            * HARTREE_TO_KJ_MOL,
            "resolved": candidate["resolved"],
        }
        for candidate in candidates_hartree
    ]
    resolved = [int(item["point_index"]) for item in candidates if item["resolved"]]
    if resolved:
        classification = "resolved_internal_maximum"
    elif raw_maxima:
        classification = "unresolved_internal_ripple"
    elif maximum == len(energies) - 1 and 0 < minimum < len(energies) - 1:
        # Keep the established public value used by path-planning policy. The
        # name predates this method-neutral module but remains schema data.
        classification = "neutral_side_minimum_then_endpoint_rise"
    else:
        classification = "endpoint_extremum_without_internal_maximum"
    return {
        "profile_classification": classification,
        "maximum_point_index": maximum,
        "minimum_point_index": minimum,
        "internal_maximum_indices": resolved,
        "raw_internal_maximum_indices": raw_maxima,
        "internal_peak_candidates": candidates,
        "min_peak_prominence_kj_mol": float(min_peak_prominence_kj_mol),
    }
