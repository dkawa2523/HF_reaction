from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _preferred_species_inputs(manifest: Manifest, states: set[str]) -> list[Artifact]:
    optimized_by_id: dict[str, Artifact] = {}
    for art in manifest.iter_artifacts("species_optimized"):
        if art.status.status == "success" and art.data.get("state") in states:
            optimized_by_id[_canonical_species_id(art)] = art
    preopt_by_id: dict[str, Artifact] = {}
    for art in manifest.iter_artifacts("species_preopt"):
        if art.status.status == "success" and art.data.get("state") in states:
            preopt_by_id[_canonical_species_id(art)] = art
    selected: list[Artifact] = []
    seen: set[str] = set()
    for species in manifest.latest_artifacts("species"):
        if species.data.get("state") not in states:
            continue
        cid = _canonical_species_id(species)
        if cid in seen:
            continue
        seen.add(cid)
        selected.append(optimized_by_id.get(cid) or preopt_by_id.get(cid) or species)
    return selected


class SinglePointStage(Stage):
    name = "sp"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine = get_qm_engine(config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {}))
        method = {
            "method_id": config.get("method", "sp"),
            **config.get("settings", {}),
            **config.get("method_settings", {}),
        }
        states = set(config.get("states", ["candidate", "isolated_candidate", "bare_candidate", "hf_cluster", "reactant_complex", "ion_pair", "transition_state"]))
        records: list[dict] = []
        for species in _preferred_species_inputs(manifest, states):
            canonical_id = _canonical_species_id(species)
            workdir = out_dir / canonical_id
            calc = engine.single_point(species, method, str(workdir))
            calc.method = {**(calc.method or {}), "stage": self.name}
            calc.data["species_id"] = canonical_id
            calc.data["source_geometry_artifact_id"] = species.artifact_id
            out.add_artifact(calc)
            records.append(calc.model_dump())
        write_jsonl(records, out_dir / "single_point_records.jsonl")
        return out
