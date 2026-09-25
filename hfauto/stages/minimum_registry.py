"""Build one same-PES registry of frequency-validated minimum basins."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from hfauto.chemistry.basin_identity import assess_minimum_pair_basin_identity
from hfauto.chemistry.minima import (
    has_frequency_validated_minimum_evidence,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import canonical_species_id, preferred_species_by_id, species_xyz_path
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import calculations_by_species


def _method_signature(species: Artifact, calculation: Artifact | None) -> dict[str, Any]:
    provenance = species.provenance or {}
    validated = dict(provenance.get("validated_qm_method") or {})
    numerical = dict(provenance.get("validated_numerical_settings") or {})
    method = calculation.method if calculation is not None else species.method or {}
    for key in ("engine", "functional", "basis", "disp_vdw", "program_version"):
        if validated.get(key) is None:
            validated[key] = (method or {}).get(key)
    return {"method": validated, "numerical_settings": numerical}


def _composition(species: Artifact) -> dict[str, int]:
    xyz = read_xyz(species_xyz_path(species))
    return dict(sorted(Counter(xyz.symbols).items()))


def _energy(calculation: Artifact | None) -> float:
    if calculation is None:
        return float("inf")
    value = calculation.data.get(
        "electronic_energy_hartree",
        calculation.data.get("optimized_electronic_energy_hartree"),
    )
    return float(value) if value is not None else float("inf")


class MinimumRegistryStage(Stage):
    """Cluster accepted DFT minima without assigning reaction roles."""

    name = "minimum-registry"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("minimum-registry requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        calculations = calculations_by_species(manifest)
        minima = list(
            preferred_species_by_id(
                manifest,
                artifact_types=("species_optimized",),
                accept=lambda species: has_frequency_validated_minimum_evidence(
                    species,
                    calculations.get(canonical_species_id(species)),
                    require_real_qm=bool(config.get("require_real_qm", True)),
                ),
            ).values()
        )
        minima.sort(
            key=lambda item: _energy(calculations.get(canonical_species_id(item)))
        )
        thresholds = {
            "same_basin_distance_rmsd_A": float(
                config.get("same_basin_distance_rmsd_A", 0.02)
            ),
            "distinct_basin_distance_rmsd_A": float(
                config.get("distinct_basin_distance_rmsd_A", 0.05)
            ),
            "same_basin_energy_tolerance_hartree": float(
                config.get("same_basin_energy_tolerance_hartree", 1.0e-5)
            ),
        }
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        species_to_basin: dict[str, str] = {}

        for minimum in minima:
            species_id = canonical_species_id(minimum)
            calculation = calculations.get(species_id)
            signature = {
                "composition": _composition(minimum),
                "charge": int(minimum.data.get("charge", 0) or 0),
                "multiplicity": int(minimum.data.get("multiplicity", 1) or 1),
                "method": _method_signature(minimum, calculation),
            }
            group_key = fingerprint_dict(signature)
            matched: dict[str, Any] | None = None
            unresolved_against: list[str] = []
            for basin in groups[group_key]:
                representative = basin["representative"]
                representative_id = canonical_species_id(representative)
                assessment = assess_minimum_pair_basin_identity(
                    representative,
                    minimum,
                    calculations.get(representative_id),
                    calculation,
                    **thresholds,
                )
                if assessment.get("same_basin") is True:
                    matched = basin
                    basin["assessments"].append(assessment)
                    break
                if assessment.get("status") == "unresolved":
                    unresolved_against.append(basin["basin_id"])
            if matched is None:
                basin_id = "basin_" + fingerprint_dict(
                    {
                        **signature,
                        "representative_species_id": species_id,
                        "xyz_sha256": sha256_file(species_xyz_path(minimum)),
                    }
                )[:24]
                matched = {
                    "basin_id": basin_id,
                    "group_key": group_key,
                    "signature": signature,
                    "representative": minimum,
                    "members": [],
                    "assessments": [],
                    "unresolved_against": unresolved_against,
                }
                groups[group_key].append(matched)
            matched["members"].append(minimum)
            species_to_basin[species_id] = matched["basin_id"]

        rows: list[dict[str, Any]] = []
        basin_ids: list[str] = []
        for group_basins in groups.values():
            for basin in group_basins:
                representative = basin["representative"]
                representative_id = canonical_species_id(representative)
                member_ids = [
                    canonical_species_id(member) for member in basin["members"]
                ]
                row = {
                    "basin_id": basin["basin_id"],
                    "representative_species_id": representative_id,
                    "representative_artifact_id": representative.artifact_id,
                    "member_species_ids": member_ids,
                    "member_artifact_ids": [
                        member.artifact_id for member in basin["members"]
                    ],
                    "calculation_artifact_ids": [
                        calculations[member_id].artifact_id
                        for member_id in member_ids
                        if member_id in calculations
                    ],
                    "signature": basin["signature"],
                    "unresolved_against": basin["unresolved_against"],
                    "same_basin_assessments": basin["assessments"],
                    "electronic_energy_hartree": (
                        None
                        if _energy(calculations.get(representative_id)) == float("inf")
                        else _energy(calculations.get(representative_id))
                    ),
                }
                artifact = Artifact(
                    artifact_id=basin["basin_id"],
                    artifact_type="minimum_basin",
                    parents=[
                        *row["member_artifact_ids"],
                        *row["calculation_artifact_ids"],
                    ],
                    paths={"representative_xyz": str(species_xyz_path(representative))},
                    data=row,
                    qc={
                        "frequency_validated_minimum_members": len(member_ids),
                        "same_pes_grouped": True,
                        "ambiguous_against_other_basins": bool(
                            basin["unresolved_against"]
                        ),
                    },
                )
                out.add_artifact(artifact)
                rows.append(row)
                basin_ids.append(artifact.artifact_id)

        registry = {
            "registry_id": "minimum_registry",
            "minimum_count": len(minima),
            "basin_count": len(rows),
            "species_to_basin": species_to_basin,
            "basin_ids": basin_ids,
            "thresholds": thresholds,
            "same_pes_partition_required": True,
        }
        registry_path = write_json(out_dir / "minimum_registry.json", registry)
        out.add_artifact(
            Artifact(
                artifact_id="minimum_registry",
                artifact_type="minimum_registry",
                parents=basin_ids,
                paths={"json": str(registry_path)},
                data=registry,
                qc={
                    "accepted_minima_only": True,
                    "global_across_reaction_trials": True,
                },
            )
        )
        write_jsonl(rows, out_dir / "minimum_basins.jsonl")
        return out
