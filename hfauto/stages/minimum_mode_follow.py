"""Create auditable +/- mode-following seeds from a rejected minimum attempt."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.minimum_recovery import (
    assess_minimum_mode_following_source,
)
from hfauto.chemistry.mode_following import (
    DEFAULT_MAXIMUM_ATOM_DISPLACEMENT_A,
    build_cartesian_mode_following_seeds,
)
from hfauto.chemistry.xyz import read_xyz, write_xyz
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _reaction_ids_for_species(manifest: Manifest, species_id: str) -> list[str]:
    matches: list[str] = []
    for reaction in manifest.latest_artifacts("reaction"):
        identifiers: set[str] = set()
        for key in (
            "reactant_species_id",
            "product_species_id",
            "reactant_species_ids",
            "product_species_ids",
        ):
            value = reaction.data.get(key)
            if isinstance(value, list):
                identifiers.update(str(item) for item in value)
            elif value:
                identifiers.add(str(value))
        if species_id in identifiers:
            matches.append(
                str(reaction.data.get("reaction_id") or reaction.artifact_id)
            )
    return matches


def _latest_minimum_attempts(manifest: Manifest) -> dict[str, Artifact]:
    attempts: dict[str, Artifact] = {}
    for calculation in manifest.latest_artifacts("calculation"):
        if calculation.data.get("task") != "opt_freq":
            continue
        if (calculation.method or {}).get("stage") != "dft-minima":
            continue
        species_id = str(calculation.data.get("species_id") or "")
        if species_id:
            attempts[species_id] = calculation
    return attempts


def _selected_sources(
    manifest: Manifest,
    config: dict[str, Any],
) -> list[Artifact]:
    requested_calculations = {
        str(identifier)
        for identifier in config.get("source_calculation_ids", [])
    }
    if requested_calculations:
        selected: list[Artifact] = []
        for identifier in sorted(requested_calculations):
            calculation = manifest.find(identifier)
            if calculation is None or calculation.artifact_type != "calculation":
                raise ValueError(
                    f"source calculation does not exist: {identifier}"
                )
            selected.append(calculation)
        return selected

    attempts = _latest_minimum_attempts(manifest)
    requested_species = {
        str(identifier) for identifier in config.get("species_ids", [])
    }
    if requested_species:
        missing = sorted(requested_species - attempts.keys())
        if missing:
            raise ValueError(
                "latest dft-minima attempt does not exist for: "
                + ", ".join(missing)
            )
        return [attempts[key] for key in sorted(requested_species)]
    return [
        attempts[key]
        for key in sorted(attempts)
        if attempts[key].qc.get("minimum_accepted") is False
    ]


def _source_species(manifest: Manifest, calculation: Artifact) -> Artifact | None:
    identifiers = [
        calculation.data.get("source_geometry_artifact_id"),
        calculation.data.get("source_species_artifact_id"),
        *(calculation.parents or []),
    ]
    for identifier in identifiers:
        if (
            identifier
            and (artifact := manifest.find(str(identifier))) is not None
            and artifact.artifact_type.startswith("species")
        ):
            return artifact
    return None


def _seed_data(
    source: Artifact | None,
    calculation: Artifact,
    *,
    seed_id: str,
    seed_state: str,
    xyz_path: str,
    direction: str,
    reaction_ids: list[str],
    generation: dict[str, Any],
) -> dict[str, Any]:
    source_data = source.data if source is not None else {}
    data = {
        key: source_data[key]
        for key in (
            "mol_id",
            "site_id",
            "hf_n",
            "system_id",
            "atom_order_key",
            "components",
            "element_counts",
            "environment",
            "reaction_coordinate",
        )
        if key in source_data
    }
    data.update(
        {
            "species_id": seed_id,
            "source_species_id": calculation.data.get("species_id"),
            "source_state": source_data.get("state"),
            "source_stationary_point_calculation_id": calculation.artifact_id,
            "state": seed_state,
            "xyz_path": xyz_path,
            "charge": calculation.data.get(
                "resolved_charge", data.get("charge")
            ),
            "multiplicity": calculation.data.get(
                "resolved_multiplicity", data.get("multiplicity")
            ),
            "electron_count": calculation.data.get(
                "electron_count", data.get("electron_count")
            ),
            "reaction_ids": reaction_ids,
            "mode_following_generation": 1,
            "mode_following_direction": direction,
            "mode_following": generation,
        }
    )
    return data


class MinimumModeFollowStage(Stage):
    """Generate two bounded seeds; never promote a saddle or a minimum."""

    name = "minimum-mode-follow"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        displacement = float(
            config.get(
                "maximum_atom_displacement_A",
                DEFAULT_MAXIMUM_ATOM_DISPLACEMENT_A,
            )
        )
        seed_state = str(
            config.get("seed_state", "endpoint_minimum_candidate")
        )
        if not seed_state.strip():
            raise ValueError("seed_state must be a non-empty string")
        plans: list[dict[str, Any]] = []
        seed_records: list[dict[str, Any]] = []

        for calculation in _selected_sources(manifest, config):
            source = _source_species(manifest, calculation)
            eligibility = assess_minimum_mode_following_source(calculation)
            if source is not None and int(
                source.data.get("mode_following_generation", 0) or 0
            ) >= 1:
                eligibility["accepted"] = False
                eligibility["reasons"] = [
                    *eligibility["reasons"],
                    "recursive_mode_following_is_not_allowed",
                ]

            source_species_id = str(
                calculation.data.get("species_id") or calculation.artifact_id
            )
            reaction_ids = _reaction_ids_for_species(
                manifest, source_species_id
            )
            plan_token = fingerprint_dict(
                {
                    "source_calculation_artifact_id": calculation.artifact_id,
                    "source_file_sha256": eligibility.get(
                        "source_file_sha256"
                    ),
                    "maximum_atom_displacement_A": displacement,
                }
            )
            plan_id = f"mode_follow_plan_{path_token(source_species_id, 48)}_{plan_token}"
            seed_ids: list[str] = []
            if eligibility["accepted"]:
                try:
                    geometry = read_xyz(calculation.paths["final_xyz"])
                    generated = build_cartesian_mode_following_seeds(
                        geometry,
                        calculation.data["imaginary_mode_displacements"],
                        maximum_atom_displacement_A=displacement,
                    )
                except (KeyError, OSError, TypeError, ValueError) as exc:
                    generation_error = f"mode_seed_generation_failed:{exc}"
                    eligibility["accepted"] = False
                    eligibility["reasons"] = [
                        *eligibility["reasons"], generation_error
                    ]
                else:
                    for direction in ("plus", "minus"):
                        seed_geometry, generation = generated[direction]
                        seed_id = (
                            f"{path_token(source_species_id, 48)}"
                            f"__mode_follow_{direction}_{plan_token[:10]}"
                        )
                        xyz_path = write_xyz(
                            seed_geometry,
                            out_dir / f"{seed_id}.xyz",
                        )
                        generation = {
                            **generation,
                            "source_calculation_artifact_id": calculation.artifact_id,
                            "source_species_artifact_id": (
                                source.artifact_id if source is not None else None
                            ),
                            "source_file_sha256": eligibility[
                                "source_file_sha256"
                            ],
                            "seed_xyz_sha256": sha256_file(xyz_path),
                        }
                        seed = Artifact(
                            artifact_id=seed_id,
                            artifact_type="species",
                            parents=[
                                identifier
                                for identifier in (
                                    source.artifact_id if source else None,
                                    calculation.artifact_id,
                                    plan_id,
                                )
                                if identifier
                            ],
                            paths={
                                "xyz": str(xyz_path),
                                "source_xyz": calculation.paths["final_xyz"],
                            },
                            data=_seed_data(
                                source,
                                calculation,
                                seed_id=seed_id,
                                seed_state=seed_state,
                                xyz_path=str(xyz_path),
                                direction=direction,
                                reaction_ids=reaction_ids,
                                generation=generation,
                            ),
                            method={"stage": self.name},
                            qc={
                                "minimum_optimization_required": True,
                                "frequency_validation_required": True,
                                "is_minimum": False,
                                "is_transition_state": False,
                                "scientific_conclusion_supported": False,
                            },
                            provenance={
                                "created_by": self.name,
                                "source_eligibility": eligibility,
                            },
                        )
                        out.add_artifact(seed)
                        seed_ids.append(seed_id)
                        seed_records.append(seed.model_dump(mode="json"))

            plan = Artifact(
                artifact_id=plan_id,
                artifact_type="minimum_mode_following_plan",
                parents=[calculation.artifact_id],
                data={
                    "accepted": eligibility["accepted"],
                    "source_species_id": source_species_id,
                    "source_calculation_artifact_id": calculation.artifact_id,
                    "reaction_ids": reaction_ids,
                    "imaginary_frequency_cm1": eligibility.get(
                        "imaginary_frequency_cm1"
                    ),
                    "maximum_atom_displacement_A": displacement,
                    "seed_state": seed_state,
                    "seed_artifact_ids": seed_ids,
                    "source_eligibility": eligibility,
                    "interpretation": "unresolved_stationary_point",
                    "required_downstream_steps": [
                        "optimize_and_frequency_validate_both_seeds",
                        "compare_frequency_validated_endpoint_basins",
                        "project_source_mode_onto_declared_reaction_coordinate",
                        "validate_source_as_first_order_saddle_only_if_basin_and_mode_gates_pass",
                        "run_bidirectional_irc_before_elementary_reaction_classification",
                    ],
                },
                qc={
                    "source_eligible": eligibility["accepted"],
                    "both_direction_seeds_created": len(seed_ids) == 2,
                    "source_promoted_to_minimum": False,
                    "source_promoted_to_transition_state": False,
                    "scientific_conclusion_supported": False,
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if eligibility["accepted"] else "failed",
                    category=None if eligibility["accepted"] else "scientific_gate",
                    reason=(
                        None
                        if eligibility["accepted"]
                        else ";".join(eligibility["reasons"])
                    ),
                    recoverable=True,
                    recommended_fallback=(
                        None
                        if eligibility["accepted"]
                        else "repair or rerun the source opt_freq calculation"
                    ),
                ),
            )
            out.add_artifact(plan)
            plans.append(plan.model_dump(mode="json"))

        write_jsonl(plans, out_dir / "minimum_mode_following_plans.jsonl")
        write_jsonl(seed_records, out_dir / "minimum_mode_following_seeds.jsonl")
        out.metadata["minimum_mode_following_plan_count"] = len(plans)
        out.metadata["minimum_mode_following_seed_count"] = len(seed_records)
        return out
