"""Create reaction candidates from trial-compatible DFT minimum basins."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.minimum_connections import assess_trial_between_minima
from hfauto.core.artifacts import preferred_species_by_id, species_xyz_path
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import ReactionCandidateRecord, ReactionTrialRecord
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class ConnectMinimaStage(Stage):
    """Connect distinct same-PES basins without relying on xTB basin topology."""

    name = "connect-minima"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("connect-minima requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        registries = manifest.latest_artifacts("minimum_registry")
        if not registries:
            raise ValueError("connect-minima requires minimum-registry evidence")
        registry = registries[-1]
        species_to_basin = dict(registry.data.get("species_to_basin") or {})
        basins = {
            basin.artifact_id: basin
            for basin in manifest.latest_artifacts("minimum_basin")
        }
        minima = preferred_species_by_id(
            manifest, artifact_types=("species_optimized",)
        )
        maximum = int(config.get("max_connections", 20))
        if maximum < 1:
            raise ValueError("connect-minima max_connections must be positive")

        trials = [
            (artifact, ReactionTrialRecord.model_validate(artifact.data))
            for artifact in manifest.latest_artifacts("reaction_trial")
            if artifact.status.status == "success"
        ]
        trials.sort(key=lambda item: (-item[1].priority, item[1].trial_id))
        existing_pairs: set[tuple[str, str, str]] = set()
        for candidate in manifest.latest_artifacts("reaction_candidate"):
            reactant_basin = species_to_basin.get(
                str(candidate.data.get("reactant_species_id"))
            )
            product_basin = species_to_basin.get(
                str(candidate.data.get("product_species_id"))
            )
            if reactant_basin and product_basin:
                existing_pairs.add(
                    (
                        *sorted((str(reactant_basin), str(product_basin))),
                        str(candidate.data.get("mechanism_family") or "discovered"),
                    )
                )

        records: list[dict[str, Any]] = []
        assessments: list[dict[str, Any]] = []
        for trial_artifact, trial in trials:
            if len(records) >= maximum:
                break
            source_basin_id = species_to_basin.get(trial.source_species_id)
            source_basin = basins.get(str(source_basin_id))
            if source_basin is None:
                continue
            source_id = str(source_basin.data.get("representative_species_id"))
            source = minima.get(source_id)
            if source is None:
                continue
            for target_basin in basins.values():
                if len(records) >= maximum:
                    break
                if target_basin.artifact_id == source_basin.artifact_id:
                    continue
                if target_basin.data.get("signature") != source_basin.data.get(
                    "signature"
                ):
                    continue
                if (
                    target_basin.artifact_id
                    in (source_basin.data.get("unresolved_against") or [])
                    or source_basin.artifact_id
                    in (target_basin.data.get("unresolved_against") or [])
                ):
                    continue
                target_id = str(
                    target_basin.data.get("representative_species_id")
                )
                target = minima.get(target_id)
                if target is None:
                    continue
                evidence = assess_trial_between_minima(
                    trial,
                    species_xyz_path(source),
                    species_xyz_path(target),
                    minimum_coordinate_change_A=float(
                        config.get("minimum_coordinate_change_A", 0.05)
                    ),
                    minimum_bond_distance_change_A=float(
                        config.get("minimum_bond_distance_change_A", 0.15)
                    ),
                )
                assessments.append(
                    {
                        "trial_id": trial.trial_id,
                        "source_basin_id": source_basin.artifact_id,
                        "target_basin_id": target_basin.artifact_id,
                        **evidence,
                    }
                )
                if evidence.get("accepted") is not True:
                    continue
                mechanism_family = (
                    "coordinate_reorganization"
                    if evidence.get("connection_kind")
                    == "coordinate_reorganization"
                    else trial.mechanism_hint
                )
                pair_key = (
                    *sorted((source_basin.artifact_id, target_basin.artifact_id)),
                    mechanism_family,
                )
                if pair_key in existing_pairs:
                    continue
                existing_pairs.add(pair_key)
                token = fingerprint_dict(
                    {
                        "trial_id": trial.trial_id,
                        "source_basin_id": source_basin.artifact_id,
                        "target_basin_id": target_basin.artifact_id,
                        "direction": evidence["direction"],
                    }
                )[:20]
                candidate_id = f"candidate_dft_basins_{token}"
                record = ReactionCandidateRecord(
                    candidate_id=candidate_id,
                    trial_id=trial.trial_id,
                    reactant_species_id=source_id,
                    product_species_id=target_id,
                    driver="dft_minimum_ensemble",
                    bond_changes=evidence["bond_changes"],
                    composition_preserved=True,
                    charge=trial.charge,
                    multiplicity=trial.multiplicity,
                    endpoint_evidence_level="dft_minimum_ensemble",
                    low_level_endpoint_distinct=None,
                    evidence_artifact_id=registry.artifact_id,
                )
                data = record.model_dump(mode="json")
                data.update(
                    {
                        "mechanism_family": mechanism_family,
                        "rationale": [
                            *trial.rationale,
                            "The pair was recovered from distinct frequency-validated DFT minima.",
                            "Bond changes are reported only when every requested distance change exceeds the bond-change threshold.",
                            "No low-level or biased energy is used as an activation barrier.",
                        ],
                        "reaction_coordinate": evidence["reaction_coordinate"],
                        "minimum_connection_evidence": evidence,
                        "reactant_basin_id": source_basin.artifact_id,
                        "product_basin_id": target_basin.artifact_id,
                    }
                )
                candidate = Artifact(
                    artifact_id=candidate_id,
                    artifact_type="reaction_candidate",
                    parents=[
                        trial_artifact.artifact_id,
                        registry.artifact_id,
                        source_basin.artifact_id,
                        target_basin.artifact_id,
                        source.artifact_id,
                        target.artifact_id,
                    ],
                    data=data,
                    qc={
                        "composition_preserved": True,
                        "structural_change_detected": True,
                        "structural_change_basis": (
                            "dft_minimum_coordinate_reorganization"
                            if evidence.get("connection_kind")
                            == "coordinate_reorganization"
                            else "dft_minimum_connectivity_change"
                        ),
                        "dft_minima_validated": True,
                        "distinct_registry_basins": True,
                        "same_pes": True,
                        "biased_energy_used_as_barrier": False,
                        "reaction_promoted": False,
                    },
                    provenance={"created_by": self.name},
                )
                out.add_artifact(candidate)
                records.append(data)

        summary = {
            "trial_count": len(trials),
            "minimum_basin_count": len(basins),
            "connection_candidate_count": len(records),
            "max_connections": maximum,
            "requires_distinct_frequency_validated_dft_basins": True,
            "low_level_energy_used_as_barrier": False,
        }
        summary_path = write_json(out_dir / "minimum_connection_summary.json", summary)
        out.add_artifact(
            Artifact(
                artifact_id="minimum_connection_summary",
                artifact_type="table",
                parents=[registry.artifact_id],
                paths={"json": str(summary_path)},
                data=summary,
                qc={"bounded_connections": len(records) <= maximum},
            )
        )
        write_jsonl(records, out_dir / "minimum_connection_candidates.jsonl")
        write_jsonl(assessments, out_dir / "minimum_connection_assessments.jsonl")
        return out
