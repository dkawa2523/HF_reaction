"""Run path calculations at multiple image resolutions."""

from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_ts_engine
from hfauto.chemistry.method_lineage import (
    evaluate_endpoint_method_lineage,
    evaluate_endpoint_path_numerical_lineage,
    evaluate_endpoint_source_integrity,
)
from hfauto.chemistry.path_convergence import path_method_signature
from hfauto.chemistry.reaction_profile import reanalyze_reaction_path
from hfauto.chemistry.xyz_trajectory import (
    normalized_xyz_arc_lengths,
    read_xyz_trajectory,
)
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl, write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.schemas.path import ReactionPathRecord
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import attempts_by_reaction
from hfauto.workflow.reaction_inputs import reaction_and_endpoints


def _path_parent(manifest: Manifest, attempt: Artifact) -> Artifact | None:
    return next(
        (
            parent
            for parent_id in attempt.parents
            if (parent := manifest.find(parent_id)) is not None
            and parent.artifact_type == "reaction_path"
        ),
        None,
    )


def _seed_attempt(manifest: Manifest, reaction_id: str) -> tuple[Artifact, Artifact]:
    attempts = attempts_by_reaction(manifest).get(reaction_id, [])
    for attempt in reversed(attempts):
        if attempt.artifact_type != "path_attempt":
            continue
        if (attempt.data.get("evidence") or {}).get("path_ensemble"):
            continue
        try:
            path = ReactionPathRecord.model_validate(attempt.data.get("path"))
        except (TypeError, ValueError):
            continue
        source = attempt.paths.get("path_xyz") or (
            attempt.data.get("evidence") or {}
        ).get("path_xyz")
        parent = _path_parent(manifest, attempt)
        if path.converged and source and parent is not None:
            return attempt, parent
    raise KeyError(f"Converged source path missing: {reaction_id}")


def _member_artifact(
    *,
    study_id: str,
    role: str,
    attempt: Artifact,
    path_artifact: Artifact,
    peak_threshold_kcal_mol: float,
    include_in_assessment: bool,
    refines_attempt_id: str | None,
    study_definition: dict[str, Any],
) -> Artifact:
    path = reanalyze_reaction_path(
        attempt.data.get("path"),
        barrier_threshold_kcal_mol=peak_threshold_kcal_mol,
    )
    token = fingerprint_dict(
        {"study": study_id, "attempt": attempt.artifact_id, "role": role}
    )
    path_xyz = attempt.paths.get("path_xyz", "")
    normalized_coordinate: list[float] = []
    try:
        normalized_coordinate = normalized_xyz_arc_lengths(
            read_xyz_trajectory(path_xyz)
        )
    except (OSError, TypeError, UnicodeError, ValueError):
        pass
    coordinate_validated = len(normalized_coordinate) == len(path.images)
    return Artifact(
        artifact_id=f"path_ensemble_member_{token}",
        artifact_type="path_ensemble_member",
        parents=[attempt.artifact_id, path_artifact.artifact_id],
        paths={
            "path_xyz": path_xyz,
            "output": path_artifact.paths.get("output", ""),
        },
        data={
            "reaction_id": path.reaction_id,
            "study_id": study_id,
            "study_definition": study_definition,
            "role": role,
            "source_path_attempt_id": attempt.artifact_id,
            "source_path_artifact_id": path_artifact.artifact_id,
            "refines_path_attempt_id": refines_attempt_id,
            "image_count": len(path.images),
            "path": path.model_dump(mode="json"),
            "method_signature": path_method_signature(path_artifact),
            "input_sha256": path_artifact.provenance.get("input_sha256"),
            "output_sha256": path_artifact.provenance.get("output_sha256"),
            "path_sha256": path_artifact.provenance.get(
                "path_sha256",
                path_artifact.provenance.get("neb_path_sha256"),
            ),
            "normalized_path_coordinate": (
                normalized_coordinate if coordinate_validated else []
            ),
            "path_coordinate_definition": (
                "normalized_cumulative_kabsch_aligned_cartesian_arc_length"
                if coordinate_validated
                else None
            ),
        },
        qc={
            "include_in_assessment": include_in_assessment,
            "real_path_executed": bool(
                path_artifact.qc.get("real_path_executed") is True
                or path_artifact.qc.get("real_neb_executed") is True
            ),
            "path_geometry_validated": path_artifact.qc.get(
                "path_geometry_validated"
            ),
            "method_evidence_validated": path_artifact.qc.get(
                "method_evidence_validated"
            ),
            "input_method_evidence_validated": path_artifact.qc.get(
                "input_method_evidence_validated"
            ),
            "endpoint_method_lineage_validated": path_artifact.qc.get(
                "endpoint_method_lineage_validated"
            ),
            "fallback_dummy": path_artifact.qc.get("fallback_dummy", False),
            "path_converged": path.converged,
            "path_coordinate_validated": coordinate_validated,
        },
        provenance={"created_by": "path-ensemble"},
    )


class PathEnsembleStage(Stage):
    """Generate discretization evidence; never promote a TS or mechanism."""

    name = "path-ensemble"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine_name = str(config.get("engine", "nwchem_string"))
        engine = get_ts_engine(engine_name, **(config.get("engine_settings") or {}))
        image_counts = sorted({int(value) for value in config.get("image_counts", [])})
        if len(image_counts) < 2 or image_counts[0] < 3:
            raise ValueError("path-ensemble requires at least two image counts >= 3")
        peak_threshold = float(config.get("peak_threshold_kcal_mol", 0.05))
        if peak_threshold <= 0.0:
            raise ValueError("peak_threshold_kcal_mol must be positive")
        method = {
            "method_id": config.get("method", "path_discretization_ensemble"),
            **(config.get("settings") or {}),
            **(config.get("method_settings") or {}),
        }
        selected_ids = {str(value) for value in config.get("reaction_ids", [])}
        checkpoint_each = bool(config.get("checkpoint_each_path", True))
        member_records: list[dict[str, Any]] = []

        for reaction in manifest.latest_artifacts("reaction"):
            reaction_id = str(reaction.data.get("reaction_id") or reaction.artifact_id)
            if selected_ids and reaction_id not in selected_ids:
                continue
            try:
                reaction, reactant, product = reaction_and_endpoints(
                    manifest, reaction_id
                )
                seed_attempt, seed_path_artifact = _seed_attempt(
                    manifest, reaction_id
                )
            except KeyError as exc:
                out.add_artifact(
                    Artifact.failure(
                        f"path_ensemble_unavailable_{path_token(reaction_id, 48)}",
                        "path_ensemble_result",
                        str(exc),
                        category="path_ensemble_seed_missing",
                        parents=[reaction.artifact_id],
                        reaction_id=reaction_id,
                    )
                )
                continue
            endpoint_lineage = evaluate_endpoint_method_lineage(
                reactant,
                product,
                method,
                engine_name,
            )
            endpoint_numerical_lineage = (
                evaluate_endpoint_path_numerical_lineage(
                    reactant,
                    product,
                    method,
                )
            )
            endpoint_source_lineage = {
                "reactant": evaluate_endpoint_source_integrity(reactant),
                "product": evaluate_endpoint_source_integrity(product),
            }
            endpoint_sources_accepted = all(
                audit["accepted"] is True
                for audit in endpoint_source_lineage.values()
            )
            if (
                endpoint_lineage["accepted"] is not True
                or endpoint_numerical_lineage["accepted"] is not True
                or not endpoint_sources_accepted
            ):
                out.add_artifact(
                    Artifact.failure(
                        f"path_ensemble_unavailable_{path_token(reaction_id, 48)}",
                        "path_ensemble_result",
                        "Path method does not match the validated endpoint PES",
                        category="method_lineage_mismatch",
                        parents=[
                            reaction.artifact_id,
                            reactant.artifact_id,
                            product.artifact_id,
                        ],
                        reaction_id=reaction_id,
                        method_lineage=endpoint_lineage,
                        numerical_lineage=endpoint_numerical_lineage,
                        source_lineage=endpoint_source_lineage,
                    )
                )
                continue
            seed_xyz = seed_attempt.paths.get("path_xyz") or (
                seed_attempt.data.get("evidence") or {}
            ).get("path_xyz")
            study_definition = {
                "reaction": reaction_id,
                "engine": engine_name,
                "seed_attempt": seed_attempt.artifact_id,
                "image_counts": image_counts,
                "method": method,
                "endpoint_method_lineage": endpoint_lineage,
                "endpoint_numerical_lineage": endpoint_numerical_lineage,
                "endpoint_source_lineage": endpoint_source_lineage,
            }
            study_id = "path_ensemble_" + fingerprint_dict(study_definition)
            members = [
                _member_artifact(
                    study_id=study_id,
                    role="seed",
                    attempt=seed_attempt,
                    path_artifact=seed_path_artifact,
                    peak_threshold_kcal_mol=peak_threshold,
                    include_in_assessment=bool(
                        config.get("include_seed_in_assessment", False)
                    ),
                    refines_attempt_id=None,
                    study_definition=study_definition,
                )
            ]
            out.add_artifact(members[0])

            existing_members = {
                int(member.data.get("image_count")): member
                for member in manifest.latest_artifacts("path_ensemble_member")
                if member.data.get("study_id") == study_id
                and member.data.get("role") == "discretization_refinement"
                and member.data.get("image_count") is not None
                and member.status.status == "success"
                and member.qc.get("real_path_executed") is True
                and member.qc.get("path_converged") is True
            }

            for image_count in image_counts:
                if image_count in existing_members:
                    members.append(existing_members[image_count])
                    continue
                ensemble_context = {
                    "study_id": study_id,
                    "role": "discretization_refinement",
                    "seed_path_attempt_id": seed_attempt.artifact_id,
                    "image_count": image_count,
                    "peak_threshold_kcal_mol": peak_threshold,
                }
                run_method = {
                    **method,
                    "path_strategy": "reparameterized_double_ended_path",
                    "endpoint_basin_status": "distinct_basin",
                    "path_image_count": image_count,
                    "path_initial_trajectory": str(seed_xyz),
                    "resample_path_initial_trajectory": True,
                    "path_barrier_threshold_kcal_mol": peak_threshold,
                    "path_ensemble": ensemble_context,
                }
                result = engine.search_ts(
                    reaction,
                    reactant,
                    product,
                    run_method,
                    out_dir / path_token(reaction_id, 42) / f"images_{image_count}",
                )
                for artifact in result.artifacts:
                    artifact.qc["endpoint_method_lineage_validated"] = True
                    artifact.provenance["endpoint_method_lineage"] = endpoint_lineage
                    artifact.provenance["endpoint_numerical_lineage"] = (
                        endpoint_numerical_lineage
                    )
                    artifact.provenance["endpoint_source_lineage"] = (
                        endpoint_source_lineage
                    )
                out.extend(result.artifacts)
                attempt = next(
                    (
                        item
                        for item in result.artifacts
                        if item.artifact_type == "path_attempt"
                    ),
                    None,
                )
                path_artifact = next(
                    (
                        item
                        for item in result.artifacts
                        if item.artifact_type == "reaction_path"
                    ),
                    None,
                )
                if attempt is not None and path_artifact is not None:
                    member = _member_artifact(
                        study_id=study_id,
                        role="discretization_refinement",
                        attempt=attempt,
                        path_artifact=path_artifact,
                        peak_threshold_kcal_mol=peak_threshold,
                        include_in_assessment=True,
                        refines_attempt_id=seed_attempt.artifact_id,
                        study_definition=study_definition,
                    )
                    members.append(member)
                    out.add_artifact(member)
                if checkpoint_each:
                    write_manifest(out, out_dir)

            member_records.extend(
                member.model_dump(mode="json") for member in members
            )

        write_jsonl(member_records, out_dir / "path_ensemble_members.jsonl")
        out.metadata["path_ensemble_member_count"] = len(member_records)
        return out
