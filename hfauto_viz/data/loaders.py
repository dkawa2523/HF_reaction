from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import pandas as pd
from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto_viz.core.paths import latest_manifest_path, resolve_artifact_path
from hfauto_viz.core.html import dataframe_to_html

class RunData:
    def __init__(self, run_dir: str | Path):
        self.run_dir = Path(run_dir)
        self.manifest_path = latest_manifest_path(self.run_dir)
        self.manifest: Manifest = read_manifest(self.manifest_path)
        self._by_type: dict[str, list[Artifact]] = defaultdict(list)
        self._latest: dict[str, Artifact] = {}
        for a in self.manifest.artifacts:
            self._by_type[a.artifact_type].append(a); self._latest[a.artifact_id] = a
    @property
    def run_id(self): return self.manifest.run_id
    def artifacts(self, artifact_type: str) -> list[Artifact]: return list(self._by_type.get(artifact_type, []))
    def latest_artifacts(self, artifact_type: str) -> list[Artifact]: return self.manifest.latest_artifacts(artifact_type)
    def find(self, artifact_id: str | None) -> Artifact | None: return self.manifest.find(artifact_id) if artifact_id else None
    def table_artifact(self, table_id: str) -> Artifact | None:
        for a in self.artifacts("table") + self.artifacts("ranking"):
            if a.artifact_id == table_id or a.data.get("table_type") == table_id or a.data.get("ranking_type") == table_id:
                return a
        return None
    def table_path(self, table_id: str) -> Path | None:
        art = self.table_artifact(table_id)
        if art:
            return resolve_artifact_path(art.paths.get("csv") or art.paths.get("parquet"), self.run_dir)
        fallback = {
            "reaction_results": ["14_rank/reaction_results.csv", "13_rank/reaction_results.csv"],
            "candidate_summary": ["14_rank/candidate_summary.csv", "13_rank/candidate_summary.csv"],
            "rank_scavenger": ["14_rank/rank_scavenger.csv", "13_rank/rank_scavenger.csv"],
            "rank_activation": ["14_rank/rank_activation.csv", "13_rank/rank_activation.csv"],
            "cluster_risk": ["14_rank/cluster_risk.csv", "13_rank/cluster_risk.csv"],
            "reaction_thermo": ["10_thermo/reaction_thermo.csv"],
        }.get(table_id, [])
        for rel in fallback:
            p = resolve_artifact_path(rel, self.run_dir)
            if p and p.exists(): return p
        return None
    def read_table(self, table_id: str) -> pd.DataFrame:
        p = self.table_path(table_id)
        if p is None or not p.exists(): return pd.DataFrame()
        return pd.read_parquet(p) if str(p).endswith(".parquet") else pd.read_csv(p)
    def species_index(self) -> dict[str, Artifact]:
        idx={}
        for t in ["species", "species_preopt", "species_optimized"]:
            for sp in self.manifest.iter_artifacts(t):
                sid=sp.data.get("species_id") or sp.data.get("source_species_id") or sp.artifact_id
                idx[str(sid)] = sp
        return idx
    def reaction_index(self) -> dict[str, Artifact]:
        idx={}
        for t in ["reaction", "reaction_validated", "reaction_path_validated"]:
            for r in self.manifest.iter_artifacts(t): idx[str(r.data.get("reaction_id") or r.artifact_id)] = r
        return idx
    def molecule_index(self) -> dict[str, Artifact]:
        idx={}
        for t in ["molecule", "molecule_enriched"]:
            for m in self.manifest.iter_artifacts(t): idx[str(m.data.get("mol_id") or m.artifact_id)] = m
        return idx
    def resolve_path(self, raw): return resolve_artifact_path(raw, self.run_dir)
    def artifact_xyz_path(self, artifact: Artifact | None) -> Path | None:
        if artifact is None: return None
        return self.resolve_path(artifact.paths.get("xyz") or artifact.paths.get("final_xyz") or artifact.data.get("xyz_path") or artifact.paths.get("source_xyz"))
    def has_dummy_or_fallback(self, parents: list[str] | None = None) -> bool:
        arts = self.manifest.artifacts if not parents else [self.find(p) for p in parents]
        for a in arts:
            if a and (a.qc.get("fallback_dummy") or a.qc.get("engine_is_dummy") or a.data.get("fallback_dummy") or a.data.get("real_orca_executed") is False): return True
        return False
