"""Extract only resolved multi-well path images as minimum-search seeds."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.reaction_classification import (
    validate_path_attempt_evidence,
)
from hfauto.chemistry.xyz import write_xyz
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.ids import path_token, species_id
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import attempts_by_reaction


class PathIntermediatesStage(Stage):
    """Turn path-local wells into auditable, backend-neutral species seeds."""

    name = "path-intermediates"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        attempts = attempts_by_reaction(manifest)
        path_artifacts = {
            artifact.artifact_id: artifact
            for artifact in manifest.latest_artifacts("reaction_path")
        }
        require_real_qm = (
            str(context.global_config.get("mode", "")).lower() == "production"
        )
        superseded_attempt_ids = {
            str(identifier)
            for assessment in manifest.latest_artifacts(
                "path_ensemble_assessment"
            )
            if assessment.status.status == "success"
            and assessment.data.get("accepted") is True
            for identifier in assessment.data.get(
                "supersedes_path_attempt_ids", []
            )
        }
        endpoints = preferred_species_by_id(
            manifest,
            artifact_types=("species_optimized", "species_preopt", "species"),
        )
        records: list[dict[str, Any]] = []

        for reaction in manifest.latest_artifacts("reaction"):
            reaction_id = str(reaction.data.get("reaction_id") or reaction.artifact_id)
            reactant_id = str(reaction.data.get("reactant_species_id") or "")
            reactant = endpoints.get(reactant_id)
            for attempt in attempts.get(reaction_id, []):
                if attempt.artifact_id in superseded_attempt_ids:
                    continue
                path, rejection_reasons = validate_path_attempt_evidence(
                    attempt,
                    path_artifacts,
                    require_real_qm=require_real_qm,
                )
                if path is None or rejection_reasons:
                    continue
                wells = [well for well in path.intermediate_wells if well.resolved]
                if not path.converged or not wells:
                    continue
                source = attempt.paths.get("path_xyz") or (
                    attempt.data.get("evidence") or {}
                ).get("path_xyz")
                if not source:
                    continue
                try:
                    images = read_xyz_trajectory(source)
                except (OSError, TypeError, ValueError):
                    continue
                if len(images) != len(path.images):
                    continue

                charge = int(
                    (reactant.data if reactant else {}).get(
                        "resolved_charge",
                        (reactant.data if reactant else {}).get("charge", 0),
                    )
                    or 0
                )
                multiplicity = int(
                    (reactant.data if reactant else {}).get(
                        "resolved_multiplicity",
                        (reactant.data if reactant else {}).get("multiplicity", 1),
                    )
                    or 1
                )
                for well in wells:
                    identifier = species_id(
                        "intermediate",
                        path_token(reaction_id, 36),
                        path_token(attempt.artifact_id, 28),
                        f"image{well.image_index}",
                    )
                    xyz_path = out_dir / f"{identifier}.xyz"
                    image = images[well.image_index]
                    image.comment = (
                        f"state=reaction_intermediate_candidate reaction_id={reaction_id} "
                        f"source_image_index={well.image_index}"
                    )
                    write_xyz(image, xyz_path)
                    artifact = Artifact(
                        artifact_id=identifier,
                        artifact_type="species",
                        parents=[reaction.artifact_id, attempt.artifact_id],
                        paths={"xyz": str(xyz_path), "source_path_xyz": str(source)},
                        data={
                            "species_id": identifier,
                            "state": "reaction_intermediate_candidate",
                            "reaction_id": reaction_id,
                            "source_path_attempt_id": attempt.artifact_id,
                            "source_image_index": well.image_index,
                            "left_maximum_index": well.left_maximum_index,
                            "right_maximum_index": well.right_maximum_index,
                            "charge": charge,
                            "multiplicity": multiplicity,
                            "xyz_path": str(xyz_path),
                        },
                        qc={
                            "path_well_resolved": True,
                            "minimum_optimization_required": True,
                            "scientific_conclusion_supported": False,
                        },
                        provenance={"created_by": self.name},
                    )
                    out.add_artifact(artifact)
                    records.append(artifact.model_dump(mode="json"))

        write_jsonl(records, out_dir / "path_intermediate_seeds.jsonl")
        out.metadata["path_intermediate_seed_count"] = len(records)
        return out
