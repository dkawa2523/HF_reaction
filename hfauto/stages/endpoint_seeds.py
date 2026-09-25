"""Generate bounded rigid-fragment seeds after endpoint-basin collapse."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.complex_seed_ensemble import generate_rigid_fragment_seeds
from hfauto.chemistry.xyz import read_xyz, write_xyz
from hfauto.core.artifacts import preferred_species_by_id, species_xyz_path
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import latest_minimum_attempts_by_species


def _selections(manifest: Manifest) -> dict[str, Artifact]:
    return {
        str(artifact.data.get("reaction_id")): artifact
        for artifact in manifest.latest_artifacts("endpoint_pair_selection")
        if artifact.data.get("reaction_id")
    }


def _source_geometry(
    root_species_id: str,
    species: dict[str, Artifact],
    minimum_attempts: dict[str, Artifact],
) -> tuple[Artifact | None, Artifact | None, Path | None]:
    source = species.get(root_species_id)
    calculation = minimum_attempts.get(root_species_id)
    calculation_path = (
        calculation.paths.get("final_xyz")
        if calculation is not None
        else None
    )
    if calculation_path and Path(calculation_path).is_file():
        return source, calculation, Path(calculation_path)
    if source is None:
        return None, calculation, None
    try:
        return source, calculation, species_xyz_path(source)
    except (OSError, TypeError, ValueError):
        return source, calculation, None


class EndpointSeedsStage(Stage):
    """Create geometry hypotheses only; optimization belongs to dft-minima."""

    name = "endpoint-seeds"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("endpoint-seeds requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        selections = _selections(manifest)
        species = preferred_species_by_id(manifest)
        minimum_attempts = latest_minimum_attempts_by_species(manifest)
        roles = tuple(str(value) for value in config.get("roles", ["reactant"]))
        if not roles or any(role not in {"reactant", "product"} for role in roles):
            raise ValueError("endpoint-seeds roles must contain reactant and/or product")
        requested = {str(value) for value in config.get("reaction_ids", [])}
        angle_degrees = tuple(
            float(value)
            for value in config.get("angle_degrees", [25.0, -25.0, 50.0, -50.0])
        )
        radial_shifts = tuple(
            float(value) for value in config.get("radial_shifts_A", [0.0])
        )
        maximum_seeds = int(config.get("maximum_seeds_per_role", 4))
        generation = int(config.get("generation", 1))
        if not 1 <= maximum_seeds <= 24:
            raise ValueError("maximum_seeds_per_role must be in [1, 24]")
        if generation != 1:
            raise ValueError("endpoint-seeds supports exactly one bounded generation")
        records: list[dict[str, Any]] = []
        plans: list[dict[str, Any]] = []

        for reaction in manifest.latest_artifacts("reaction"):
            reaction_id = str(
                reaction.data.get("reaction_id") or reaction.artifact_id
            )
            if requested and reaction_id not in requested:
                continue
            selection = selections.get(reaction_id)
            if selection is None or selection.qc.get("endpoint_pair_ready") is True:
                continue
            selection_reasons = set(selection.data.get("reasons") or [])
            if not selection_reasons.intersection(
                {
                    "all_endpoint_candidates_collapse_to_same_basin",
                    "frequency_validated_reactant_candidates_missing",
                    "frequency_validated_product_candidates_missing",
                    "chemically_valid_distinct_endpoint_pair_missing",
                }
            ):
                continue
            for role in roles:
                root_key = f"{role}_root_species_id"
                root_species_id = str(
                    selection.data.get(root_key)
                    or reaction.data.get(f"{role}_species_id")
                    or ""
                )
                source, calculation, xyz_path = _source_geometry(
                    root_species_id, species, minimum_attempts
                )
                plan_token = fingerprint_dict(
                    {
                        "reaction_id": reaction_id,
                        "role": role,
                        "root_species_id": root_species_id,
                        "generation": generation,
                        "angles": angle_degrees,
                        "radial_shifts": radial_shifts,
                        "maximum_seeds": maximum_seeds,
                    }
                )
                plan_id = f"endpoint_seed_plan_{path_token(reaction_id, 36)}_{role}_{plan_token}"
                if manifest.find(plan_id) is not None:
                    continue
                reasons: list[str] = []
                if source is None:
                    reasons.append("endpoint_root_species_missing")
                if xyz_path is None:
                    reasons.append("endpoint_source_geometry_missing")
                generated: list[tuple[Any, dict[str, Any]]] = []
                if not reasons:
                    try:
                        generated = generate_rigid_fragment_seeds(
                            read_xyz(xyz_path),
                            reaction.data,
                            angle_degrees=angle_degrees,
                            radial_shifts_A=radial_shifts,
                            minimum_interfragment_distance_A=float(
                                config.get(
                                    "minimum_interfragment_distance_A", 0.70
                                )
                            ),
                            maximum_seeds=maximum_seeds,
                        )
                    except (OSError, TypeError, ValueError) as exc:
                        reasons.append(f"endpoint_seed_generation_failed:{exc}")
                if not generated and not reasons:
                    reasons.append("no_disconnected_rigid_fragments_found")
                seed_ids: list[str] = []
                source_hash = sha256_file(xyz_path) if xyz_path is not None else None
                for seed, seed_evidence in generated:
                    seed_token = fingerprint_dict(
                        {
                            "plan": plan_token,
                            "template": seed_evidence["template_index"],
                            "source_sha256": source_hash,
                        }
                    )
                    identifier = (
                        f"{path_token(root_species_id, 44)}"
                        f"__fragment_seed_{seed_token}"
                    )
                    seed_path = out_dir / f"{identifier}.xyz"
                    write_xyz(seed, seed_path)
                    data = {
                        "species_id": identifier,
                        "state": "endpoint_minimum_candidate",
                        "endpoint_role": role,
                        "source_species_id": root_species_id,
                        "source_geometry_artifact_id": (
                            source.artifact_id if source is not None else None
                        ),
                        "source_stationary_point_calculation_id": (
                            calculation.artifact_id
                            if calculation is not None
                            else None
                        ),
                        "reaction_ids": [reaction_id],
                        "endpoint_seed_generation": generation,
                        "charge": int((source.data if source else {}).get("charge", 0)),
                        "multiplicity": int(
                            (source.data if source else {}).get("multiplicity", 1)
                        ),
                        "resolved_charge": int(
                            (source.data if source else {}).get(
                                "resolved_charge",
                                (source.data if source else {}).get("charge", 0),
                            )
                        ),
                        "resolved_multiplicity": int(
                            (source.data if source else {}).get(
                                "resolved_multiplicity",
                                (source.data if source else {}).get(
                                    "multiplicity", 1
                                ),
                            )
                        ),
                        "atom_order_key": (
                            (source.data if source else {}).get("atom_order_key")
                        ),
                        "xyz_path": str(seed_path),
                        "fragment_seed": seed_evidence,
                    }
                    artifact = Artifact(
                        artifact_id=identifier,
                        artifact_type="species",
                        parents=[
                            reaction.artifact_id,
                            selection.artifact_id,
                            *(
                                [source.artifact_id]
                                if source is not None
                                else []
                            ),
                            *(
                                [calculation.artifact_id]
                                if calculation is not None
                                else []
                            ),
                        ],
                        paths={"xyz": str(seed_path)},
                        data=data,
                        method={"stage": self.name},
                        provenance={
                            "created_by": self.name,
                            "source_xyz_sha256": source_hash,
                            "seed_xyz_sha256": sha256_file(seed_path),
                        },
                        qc={
                            "minimum_optimization_required": True,
                            "frequency_validation_required": True,
                            "minimum_claimed": False,
                            "transition_state_claimed": False,
                            "chemical_state_assigned": False,
                        },
                    )
                    out.add_artifact(artifact)
                    records.append(artifact.model_dump(mode="json"))
                    seed_ids.append(identifier)
                plan_data = {
                    "reaction_id": reaction_id,
                    "endpoint_role": role,
                    "root_species_id": root_species_id,
                    "source_species_artifact_id": (
                        source.artifact_id if source is not None else None
                    ),
                    "source_calculation_artifact_id": (
                        calculation.artifact_id if calculation is not None else None
                    ),
                    "source_xyz_path": str(xyz_path) if xyz_path else None,
                    "source_xyz_sha256": source_hash,
                    "generation": generation,
                    "seed_artifact_ids": seed_ids,
                    "reasons": reasons,
                    "next_action": (
                        "preopt_then_dft_minima"
                        if seed_ids
                        else "provide_fragmented_endpoint_geometry"
                    ),
                }
                plan = Artifact(
                    artifact_id=plan_id,
                    artifact_type="endpoint_seed_ensemble_plan",
                    parents=[reaction.artifact_id, selection.artifact_id],
                    data=plan_data,
                    provenance={"created_by": self.name},
                    qc={
                        "bounded_generation": len(seed_ids) <= maximum_seeds,
                        "seed_count": len(seed_ids),
                        "minimum_claimed": False,
                    },
                    status=ArtifactStatus(
                        status="success" if seed_ids else "partial",
                        category=None if seed_ids else "endpoint_seed_unavailable",
                        reason="; ".join(reasons),
                    ),
                )
                out.add_artifact(plan)
                plans.append(plan.model_dump(mode="json"))
        write_jsonl(records, out_dir / "endpoint_seed_records.jsonl")
        write_jsonl(plans, out_dir / "endpoint_seed_plans.jsonl")
        out.metadata["endpoint_seed_count"] = len(records)
        out.metadata["endpoint_seed_plan_count"] = len(plans)
        return out
