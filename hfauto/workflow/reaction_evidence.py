"""Small manifest joins shared by reaction planning and classification."""

from __future__ import annotations

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


def calculations_by_species(manifest: Manifest) -> dict[str, Artifact]:
    selected: dict[str, Artifact] = {}
    for calculation in manifest.latest_artifacts("calculation"):
        species_id = calculation.data.get("species_id")
        if species_id and calculation.status.status == "success":
            selected[str(species_id)] = calculation
    return selected


def latest_minimum_attempts_by_species(
    manifest: Manifest,
) -> dict[str, Artifact]:
    """Return the latest explicit DFT-minimum attempt, accepted or rejected."""

    selected: dict[str, Artifact] = {}
    for calculation in manifest.latest_artifacts("calculation"):
        method = calculation.method or {}
        if (
            method.get("stage") != "dft-minima"
            or calculation.data.get("task") != "opt_freq"
        ):
            continue
        species_id = calculation.data.get("species_id")
        if species_id:
            selected[str(species_id)] = calculation
    return selected


def artifacts_by_reaction(
    manifest: Manifest,
    artifact_type: str,
) -> dict[str, list[Artifact]]:
    grouped: dict[str, list[Artifact]] = {}
    for artifact in manifest.latest_artifacts(artifact_type):
        reaction_ids = artifact.data.get("reaction_ids")
        identifiers = (
            [str(value) for value in reaction_ids]
            if isinstance(reaction_ids, list)
            else []
        )
        reaction_id = artifact.data.get("reaction_id")
        if reaction_id:
            identifiers.append(str(reaction_id))
        for identifier in dict.fromkeys(identifiers):
            grouped.setdefault(identifier, []).append(artifact)
    return grouped


def attempts_by_reaction(manifest: Manifest) -> dict[str, list[Artifact]]:
    """Return latest path/saddle revisions in their scientific order."""

    grouped: dict[str, list[Artifact]] = {}
    latest = {
        artifact.artifact_id: artifact
        for artifact_type in ("path_attempt", "saddle_attempt")
        for artifact in manifest.latest_artifacts(artifact_type)
    }
    seen: set[str] = set()
    for original in manifest.artifacts:
        if original.artifact_id not in latest or original.artifact_id in seen:
            continue
        artifact = latest[original.artifact_id]
        seen.add(artifact.artifact_id)
        reaction_id = artifact.data.get("reaction_id")
        if reaction_id:
            grouped.setdefault(str(reaction_id), []).append(artifact)
    return grouped
