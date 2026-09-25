from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto_viz.core.paths import latest_manifest_path, resolve_artifact_path

TABLE_ARTIFACT_ALIASES: dict[str, tuple[str, ...]] = {
    "reaction_results": ("reaction_ranking", "reaction_thermo_table"),
    "candidate_summary": ("reaction_ranking",),
    "rank_activation": ("reaction_ranking",),
    "reaction_thermo": ("reaction_thermo_table",),
}


class RunData:
    def __init__(self, run_dir: str | Path):
        self.run_dir = Path(run_dir)
        self.manifest_path = latest_manifest_path(self.run_dir)
        self.manifest: Manifest = read_manifest(self.manifest_path)
        self._by_type: dict[str, list[Artifact]] = defaultdict(list)
        for artifact in self.manifest.latest_artifacts():
            self._by_type[artifact.artifact_type].append(artifact)

    @property
    def run_id(self) -> str:
        return self.manifest.run_id

    def artifacts(self, artifact_type: str) -> list[Artifact]:
        return list(self._by_type.get(artifact_type, []))

    def latest_artifacts(self, artifact_type: str) -> list[Artifact]:
        return self.manifest.latest_artifacts(artifact_type)

    def find(self, artifact_id: str | None) -> Artifact | None:
        return self.manifest.find(artifact_id) if artifact_id else None

    def table_artifact(self, table_id: str) -> Artifact | None:
        accepted_ids = {table_id, *TABLE_ARTIFACT_ALIASES.get(table_id, ())}
        artifacts = self.manifest.latest_artifacts("table") + self.manifest.latest_artifacts(
            "ranking"
        )
        for a in artifacts:
            if (
                a.artifact_id in accepted_ids
                or a.data.get("table_type") in accepted_ids
                or a.data.get("ranking_type") in accepted_ids
            ):
                return a
        return None

    def table_path(self, table_id: str) -> Path | None:
        art = self.table_artifact(table_id)
        if art is None:
            return None
        raw_path = art.paths.get("csv") or art.paths.get("parquet")
        return resolve_artifact_path(raw_path, self.run_dir)

    def read_table(self, table_id: str) -> pd.DataFrame:
        path = self.table_path(table_id)
        if path is None or not path.exists():
            return pd.DataFrame()
        return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)

    def species_index(self) -> dict[str, Artifact]:
        index: dict[str, Artifact] = {}
        for artifact_type in ("species", "species_preopt", "species_optimized"):
            for species in self.manifest.latest_artifacts(artifact_type):
                species_id = (
                    species.data.get("species_id")
                    or species.data.get("source_species_id")
                    or species.artifact_id
                )
                index[str(species_id)] = species
        return index

    def reaction_index(self) -> dict[str, Artifact]:
        index: dict[str, Artifact] = {}
        for artifact_type in ("reaction", "reaction_validated", "reaction_path_validated"):
            for reaction in self.manifest.latest_artifacts(artifact_type):
                reaction_id = reaction.data.get("reaction_id") or reaction.artifact_id
                index[str(reaction_id)] = reaction
        return index

    def molecule_index(self) -> dict[str, Artifact]:
        index: dict[str, Artifact] = {}
        for artifact_type in ("molecule", "molecule_enriched"):
            for molecule in self.manifest.latest_artifacts(artifact_type):
                molecule_id = molecule.data.get("mol_id") or molecule.artifact_id
                index[str(molecule_id)] = molecule
        return index

    def resolve_path(self, raw: str | Path | None) -> Path | None:
        return resolve_artifact_path(raw, self.run_dir)

    def artifact_xyz_path(self, artifact: Artifact | None) -> Path | None:
        if artifact is None:
            return None
        raw_path = (
            artifact.paths.get("xyz")
            or artifact.paths.get("final_xyz")
            or artifact.data.get("xyz_path")
            or artifact.paths.get("source_xyz")
        )
        return self.resolve_path(raw_path)

    def has_dummy_or_fallback(self, parents: list[str] | None = None) -> bool:
        artifacts = (
            self.manifest.latest_artifacts()
            if not parents
            else [self.find(parent_id) for parent_id in parents]
        )
        for artifact in artifacts:
            if artifact is None:
                continue
            if (
                artifact.qc.get("fallback_dummy")
                or artifact.qc.get("engine_is_dummy")
                or artifact.data.get("fallback_dummy")
                or artifact.data.get("real_orca_executed") is False
            ):
                return True
        return False
