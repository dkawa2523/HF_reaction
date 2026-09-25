from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.core.artifacts import canonical_species_id, preferred_species
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _preferred_species_inputs(manifest: Manifest, states: set[str]) -> list[Artifact]:
    return preferred_species(manifest, states=states)


def _validated_transition_state_inputs(manifest: Manifest) -> list[Artifact]:
    """Return one real, frequency-validated TS geometry per canonical ID."""

    by_species: dict[str, Artifact] = {}
    for calculation in manifest.latest_artifacts("calculation"):
        if (
            calculation.status.status != "success"
            or calculation.data.get("state") != "transition_state"
            or calculation.qc.get("ts_validated_by_frequency") is not True
            or calculation.qc.get("real_qm_executed") is not True
            or not calculation.paths.get("final_xyz")
        ):
            continue
        by_species[canonical_species_id(calculation)] = calculation
    return [by_species[key] for key in sorted(by_species)]


def _single_point_inputs(
    manifest: Manifest,
    states: set[str],
    *,
    include_validated_transition_states: bool,
) -> list[Artifact]:
    inputs = _preferred_species_inputs(manifest, states)
    if include_validated_transition_states and "transition_state" in states:
        inputs = [
            species
            for species in inputs
            if species.data.get("state") != "transition_state"
        ]
        inputs.extend(_validated_transition_state_inputs(manifest))
    return inputs


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
        states = set(
            config.get(
                "states",
                [
                    "candidate",
                    "isolated_candidate",
                    "bare_candidate",
                    "hf_cluster",
                    "encounter_complex",
                    "reactant",
                    "reactant_complex",
                    "product",
                    "product_complex",
                    "adsorbed_reactant",
                    "adsorbed_product",
                    "ion_pair",
                    "transition_state",
                ],
            )
        )
        records: list[dict] = []
        for species in _single_point_inputs(
            manifest,
            states,
            include_validated_transition_states=bool(
                config.get("include_validated_transition_states", False)
            ),
        ):
            canonical_id = canonical_species_id(species)
            workdir = out_dir / canonical_id
            calc = engine.single_point(species, method, str(workdir))
            calc.method = {**(calc.method or {}), "stage": self.name}
            calc.data["species_id"] = canonical_id
            calc.data["source_geometry_artifact_id"] = species.artifact_id
            out.add_artifact(calc)
            records.append(calc.model_dump())
        write_jsonl(records, out_dir / "single_point_records.jsonl")
        return out
