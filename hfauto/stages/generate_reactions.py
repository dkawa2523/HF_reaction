"""Generate bounded single-ended reaction trials from molecular complexes."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.reaction_trials import (
    automatic_reaction_trials,
    deduplicate_trials,
    explicit_reaction_trials,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import (
    artifact_data_matches,
    canonical_species_id,
    preferred_species,
)
from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class GenerateReactionsStage(Stage):
    """Propose search coordinates without constructing a product endpoint."""

    name = "generate-reactions"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("generate-reactions requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        states = {
            str(value)
            for value in config.get(
                "states", ["encounter_complex", "candidate", "isolated_candidate"]
            )
        }
        sources = preferred_species(
            manifest,
            states=states,
            artifact_types=("species", "species_preopt"),
        )
        source_filters = dict(config.get("source_data_filters") or {})
        required_source_keys = {
            str(key) for key in config.get("required_source_data_keys", [])
        }
        sources = [
            source
            for source in sources
            if artifact_data_matches(source, source_filters)
            and all(source.data.get(key) is not None for key in required_source_keys)
        ]
        sources.sort(
            key=lambda source: (
                source.data.get("relative_energy_kcal_mol") is None,
                float(source.data.get("relative_energy_kcal_mol", float("inf"))),
                canonical_species_id(source),
            )
        )
        default_driver_order = [
            str(value) for value in config.get("driver_order", ["nt2", "afir"])
        ]
        max_attempts = int(config.get("max_attempts_per_trial", 2))
        max_per_species = int(config.get("max_trials_per_species", 12))
        max_total = int(config.get("max_total_trials", 100))
        if min(max_attempts, max_per_species, max_total) < 1:
            raise ValueError("reaction-trial budgets must be positive")
        explicit_definitions = list(config.get("explicit_trials", []) or [])
        automatic = bool(config.get("automatic", True))
        records: list[dict[str, Any]] = []
        batches: list[tuple[Artifact, list[Any]]] = []

        for source in sources:
            source_id = canonical_species_id(source)
            xyz = read_xyz(source.data.get("xyz_path") or source.paths["xyz"])
            components = list(source.data.get("components", []) or [])
            if not components:
                components = [
                    {
                        "component_id": source_id,
                        "atom_indices": list(range(len(xyz.symbols))),
                    }
                ]
            component_ids = [
                str(component.get("component_id") or f"fragment_{index}")
                for index, component in enumerate(components)
            ]
            charge = int(source.data.get("charge", 0) or 0)
            multiplicity = int(source.data.get("multiplicity", 1) or 1)
            trials = explicit_reaction_trials(
                source_id,
                len(xyz.symbols),
                explicit_definitions,
                component_ids=component_ids,
                charge=charge,
                multiplicity=multiplicity,
                default_driver_order=default_driver_order,
                default_max_attempts=max_attempts,
            )
            if automatic:
                trials.extend(
                    automatic_reaction_trials(
                        source_id,
                        xyz,
                        components,
                        charge=charge,
                        multiplicity=multiplicity,
                        driver_order=default_driver_order,
                        max_attempts=max_attempts,
                        max_relay_steps=int(config.get("max_relay_steps", 2)),
                        relay_contact_cutoff_A=float(
                            config.get("relay_contact_cutoff_A", 3.2)
                        ),
                    )
                )
            batches.append(
                (source, deduplicate_trials(trials)[:max_per_species])
            )

        # Apply a global budget fairly across conformers.  Filling the budget
        # from the first source silently defeats an upstream NCI ensemble.
        max_batch_size = max((len(trials) for _, trials in batches), default=0)
        for trial_index in range(max_batch_size):
            for source, trials in batches:
                if len(records) >= max_total:
                    break
                if trial_index >= len(trials):
                    continue
                trial = trials[trial_index]
                data = trial.model_dump(mode="json")
                data.update(
                    {
                        "source_geometry_artifact_id": source.artifact_id,
                        "source_xyz_path": str(
                            source.data.get("xyz_path") or source.paths.get("xyz")
                        ),
                        "bond_changes": [
                            change.model_dump(mode="json")
                            for change in trial.bond_changes()
                        ],
                        "reaction_coordinate": trial.reaction_coordinate.model_dump(
                            mode="json"
                        ),
                        "product_geometry_claimed": False,
                        "barrier_claimed": False,
                    }
                )
                artifact = Artifact(
                    artifact_id=trial.trial_id,
                    artifact_type="reaction_trial",
                    parents=[source.artifact_id],
                    paths={"reactant": data["source_xyz_path"]},
                    data=data,
                    qc={
                        "bounded_trial": True,
                        "atom_mapping_validated": True,
                        "product_geometry_claimed": False,
                    },
                    provenance={"created_by": self.name},
                )
                out.add_artifact(artifact)
                records.append(data)
            if len(records) >= max_total:
                break

        summary = {
            "source_species_count": len(sources),
            "reaction_trial_count": len(records),
            "automatic_proposals": automatic,
            "max_trials_per_species": max_per_species,
            "max_total_trials": max_total,
            "driver_order": default_driver_order,
            "source_data_filters": source_filters,
            "required_source_data_keys": sorted(required_source_keys),
            "budget_distribution": "round_robin_by_source",
        }
        if not records:
            out.add_artifact(
                Artifact.failure(
                    "reaction_trial_empty",
                    "reaction_trial",
                    "no bounded reaction trial was generated",
                    category="no_reaction_trials",
                    recoverable=True,
                    recommended_fallback="provide explicit_trials or broader compositions",
                )
            )
        summary_path = write_json(out_dir / "reaction_trial_summary.json", summary)
        out.add_artifact(
            Artifact(
                artifact_id="reaction_trial_summary",
                artifact_type="table",
                paths={"json": str(summary_path)},
                data=summary,
                qc={"bounded_generation": len(records) <= max_total},
            )
        )
        write_jsonl(records, out_dir / "reaction_trials.jsonl")
        return out
