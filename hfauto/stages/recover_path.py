"""Recover an interrupted double-ended path as a geometry-only saddle seed."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.ts.nwchem_path_support import validate_path_trajectory
from hfauto.chemistry.path_diagnostics import (
    make_path_attempt_record,
    parse_path_optimization_history,
)
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.chemistry.xyz import read_xyz, write_xyz
from hfauto.chemistry.xyz_trajectory import (
    reaction_mode_reference,
    read_xyz_trajectory,
)
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.stages.reaction_plan import make_reaction_case_artifact
from hfauto.workflow.reaction_evidence import attempts_by_reaction
from hfauto.workflow.reaction_inputs import reaction_and_endpoints
from hfauto.workflow.reaction_state import decide_reaction_case


def _reaction_case_by_id(manifest: Manifest) -> dict[str, Artifact]:
    return {
        str(item.data.get("reaction_id")): item
        for item in manifest.latest_artifacts("reaction_case")
        if item.data.get("reaction_id")
    }


class RecoverPathStage(Stage):
    """Translate preserved raw paths into seeds without publishing path energies."""

    name = "recover-path"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("recover-path requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        existing_attempts = attempts_by_reaction(manifest)
        existing_cases = _reaction_case_by_id(manifest)
        records: list[dict[str, Any]] = []

        for request in config.get("paths", []) or []:
            reaction_id = str(request["reaction_id"])
            raw_path = Path(str(request["path_xyz"]))
            if not raw_path.is_file():
                raise FileNotFoundError(f"recovery path does not exist: {raw_path}")
            reaction, reactant, product = reaction_and_endpoints(
                manifest, reaction_id
            )
            images = read_xyz_trajectory(raw_path)
            validate_path_trajectory(
                images,
                read_xyz(species_xyz_path(reactant)),
                read_xyz(species_xyz_path(product)),
                endpoint_tolerance_A=float(
                    request.get("endpoint_rmsd_tolerance_A", 0.10)
                ),
            )
            path_record = analyze_reaction_path(
                reaction_id=reaction_id,
                engine=str(request.get("engine", "recovered_path")),
                comments=[image.comment for image in images],
                converged=False,
                barrier_threshold_kcal_mol=float(
                    request.get("barrier_threshold_kcal_mol", 0.5)
                ),
                energy_source="recovered_unconverged_path_seed_only",
            )
            if (
                path_record.classification != "resolved_internal_maximum"
                or len(path_record.internal_maximum_indices) != 1
                or path_record.highest_internal_image_index is None
            ):
                raise ValueError(
                    f"recovery path for {reaction_id} has no unique internal maximum"
                )
            image_index = path_record.highest_internal_image_index
            token = fingerprint_dict(
                {
                    "reaction_id": reaction_id,
                    "path_sha256": sha256_file(raw_path),
                    "image_index": image_index,
                }
            )
            optimization_path = Path(str(request.get("optimization_history") or ""))
            engine = str(request.get("engine", "recovered_path"))
            optimization = parse_path_optimization_history(
                engine,
                optimization_path if optimization_path.is_file() else None,
                converged=False,
            )
            attempt_record = make_path_attempt_record(
                attempt_id=f"path_attempt_recovered_{token}",
                reaction_id=reaction_id,
                strategy="recovered_unconverged_path",
                engine=engine,
                endpoint_basin_status="distinct_basin",
                path=path_record,
                optimization=optimization,
                evidence={
                    "path_xyz": str(raw_path),
                    "path_sha256": sha256_file(raw_path),
                    "geometry_seed_only": True,
                },
            )
            path_artifact = Artifact(
                artifact_id=f"reaction_path_recovered_{token}",
                artifact_type="reaction_path",
                parents=[
                    reaction.artifact_id,
                    reactant.artifact_id,
                    product.artifact_id,
                ],
                paths={
                    "path_xyz": str(raw_path),
                    "optimization_history": (
                        str(optimization_path) if optimization_path.is_file() else ""
                    ),
                },
                data={
                    "reaction_id": reaction_id,
                    "reaction_path": path_record.model_dump(mode="json"),
                    "scientific_role": "saddle_seed_source_only",
                },
                qc={
                    "real_path_started": True,
                    "path_converged": False,
                    "path_geometry_validated": True,
                    "path_energy_publishable": False,
                    "unconverged_path_used_as_seed_only": True,
                    "fallback_dummy": False,
                },
                provenance={
                    "created_by": self.name,
                    "raw_path_sha256": sha256_file(raw_path),
                },
                status=ArtifactStatus(
                    status="partial",
                    category="interrupted_path_recovered_as_seed",
                    reason="path energy is not publishable; geometry is seed-only",
                ),
            )
            attempt = Artifact(
                artifact_id=attempt_record.attempt_id,
                artifact_type="path_attempt",
                parents=[path_artifact.artifact_id],
                paths={"path_xyz": str(raw_path)},
                data=attempt_record.model_dump(mode="json"),
                qc={
                    "path_converged": False,
                    "path_energy_publishable": False,
                    "saddle_refinement_allowed": True,
                    "unconverged_path_used_as_seed_only": True,
                },
                provenance={"created_by": self.name},
            )
            seed_id = f"saddle_seed_recovered_{token}"
            seed_path = out_dir / f"{path_token(reaction_id, 40)}_seed.xyz"
            write_xyz(images[image_index], seed_path)
            mode_reference = reaction_mode_reference(
                images,
                image_index,
                source="recovered_unconverged_path_seed_only",
            )
            seed = Artifact(
                artifact_id=seed_id,
                artifact_type="species",
                parents=[path_artifact.artifact_id, attempt.artifact_id],
                paths={"xyz": str(seed_path)},
                data={
                    **reactant.data,
                    "species_id": seed_id,
                    "state": "saddle_seed",
                    "xyz_path": str(seed_path),
                    "reaction_id": reaction_id,
                    "reaction_mode_reference": mode_reference,
                },
                qc={
                    "scientific_role": "saddle_seed_only",
                    "path_converged": False,
                    "path_energy_publishable": False,
                    "path_geometry_validated": True,
                },
                provenance={
                    "created_by": self.name,
                    "source_path_sha256": sha256_file(raw_path),
                },
            )
            candidate = {
                "xyz_path": str(seed_path),
                "species_artifact_id": seed_id,
                "path_artifact_id": path_artifact.artifact_id,
                "path_image_index": image_index,
                "source": "recovered_unconverged_path_seed_only",
                "reaction_mode_reference": mode_reference,
            }
            saddle_attempt = Artifact(
                artifact_id=f"saddle_attempt_recovered_{token}",
                artifact_type="saddle_attempt",
                parents=[attempt.artifact_id, seed.artifact_id],
                paths={"seed_xyz": str(seed_path)},
                data={
                    "reaction_id": reaction_id,
                    "strategy": "recovered_unconverged_path",
                    "engine": engine,
                    "diagnosis": "resolved_saddle_candidate",
                    "next_action": "refine_saddle",
                    "candidate": candidate,
                    "search_evidence_validated": True,
                },
                qc={
                    "path_converged": False,
                    "path_energy_publishable": False,
                    "saddle_seed_resolved": True,
                    "path_geometry_validated": True,
                },
                provenance={"created_by": self.name},
            )
            previous_case = existing_cases.get(reaction_id)
            basin_assessment = dict(
                (previous_case.data if previous_case else reaction.data).get(
                    "basin_assessment"
                )
                or {"status": "distinct_basin", "accepted": True}
            )
            case = decide_reaction_case(
                reaction_id=reaction_id,
                basin_assessment=basin_assessment,
                path_attempts=[
                    *existing_attempts.get(reaction_id, []),
                    attempt,
                    saddle_attempt,
                ],
                max_path_attempts=int(request.get("max_path_attempts", 4)),
            )
            case_artifact = make_reaction_case_artifact(
                case,
                parents=[
                    *([previous_case.artifact_id] if previous_case else []),
                    saddle_attempt.artifact_id,
                ],
                created_by=self.name,
            )
            out.extend(
                [path_artifact, attempt, seed, saddle_attempt, case_artifact]
            )
            records.append(
                {
                    "reaction_id": reaction_id,
                    "path_artifact_id": path_artifact.artifact_id,
                    "saddle_seed_artifact_id": seed.artifact_id,
                    "saddle_attempt_artifact_id": saddle_attempt.artifact_id,
                    "reaction_case_artifact_id": case_artifact.artifact_id,
                    "path_energy_publishable": False,
                }
            )

        write_jsonl(records, out_dir / "recovered_paths.jsonl")
        return out
