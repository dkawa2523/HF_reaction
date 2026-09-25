"""Artifact lookup helpers shared by stages and calculation backends."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest

SPECIES_ARTIFACT_PRIORITY = ("species", "species_preopt", "species_optimized")


def artifact_data_matches(artifact: Artifact, filters: dict | None) -> bool:
    """Return whether an artifact matches a small exact-value data selector."""

    return all(
        artifact.data.get(str(key)) == value
        for key, value in (filters or {}).items()
    )


def canonical_species_id(species: Artifact) -> str:
    """Return the stable species id across constructed and optimized revisions."""

    return str(
        species.data.get("species_id")
        or species.data.get("source_species_id")
        or species.artifact_id
    )


def species_xyz_path(species: Artifact, *, must_exist: bool = True) -> Path:
    """Resolve a species geometry using the one canonical path precedence."""

    value = (
        species.data.get("xyz_path")
        or species.paths.get("xyz")
        or species.paths.get("final_xyz")
    )
    if not value:
        raise ValueError(f"Species has no XYZ geometry: {species.artifact_id}")
    path = Path(str(value))
    if must_exist and not path.is_file():
        raise FileNotFoundError(f"XYZ file does not exist for {species.artifact_id}: {path}")
    return path


def preferred_species_by_id(
    manifest: Manifest,
    *,
    states: Iterable[str] | None = None,
    artifact_types: Iterable[str] = SPECIES_ARTIFACT_PRIORITY,
    accept: Callable[[Artifact], bool] | None = None,
) -> dict[str, Artifact]:
    """Select the highest-priority successful geometry for each chemical state.

    Artifact types are ordered from lowest to highest priority. Later revisions
    replace earlier ones without changing the canonical species id.
    """

    allowed_states = set(states) if states is not None else None
    selected: dict[str, Artifact] = {}
    for artifact_type in artifact_types:
        for artifact in manifest.iter_artifacts(artifact_type):
            if artifact.status.status != "success":
                continue
            if allowed_states is not None and artifact.data.get("state") not in allowed_states:
                continue
            if accept is not None and not accept(artifact):
                continue
            selected[canonical_species_id(artifact)] = artifact
    return selected


def preferred_species(
    manifest: Manifest,
    *,
    states: Iterable[str] | None = None,
    artifact_types: Iterable[str] = SPECIES_ARTIFACT_PRIORITY,
    accept: Callable[[Artifact], bool] | None = None,
) -> list[Artifact]:
    """Return preferred species in stable canonical-id order."""

    return list(
        preferred_species_by_id(
            manifest,
            states=states,
            artifact_types=artifact_types,
            accept=accept,
        ).values()
    )
