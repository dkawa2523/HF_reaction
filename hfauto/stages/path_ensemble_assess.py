"""Audit and assess completed multi-resolution reaction paths."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from hfauto.backends.ts.nwchem_path_support import (
    assess_nwchem_rendered_method,
    existing_file_hash,
)
from hfauto.chemistry.method_lineage import (
    evaluate_endpoint_method_lineage,
    evaluate_endpoint_path_numerical_lineage,
    evaluate_endpoint_source_integrity,
)
from hfauto.chemistry.path_convergence import (
    assess_path_ensemble,
    path_method_signature,
)
from hfauto.chemistry.reaction_profile import reanalyze_reaction_path
from hfauto.chemistry.xyz_trajectory import (
    normalized_xyz_arc_lengths,
    read_xyz_trajectory,
)
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_inputs import reaction_and_endpoints


def _input_method_audit(manifest: Manifest, member: Artifact) -> dict[str, Any]:
    source_id = str(member.data.get("source_path_artifact_id") or "")
    source = manifest.find(source_id)
    if source is None:
        return {"accepted": False, "reasons": ["source_path_artifact_missing"]}
    input_path = source.paths.get("input")
    output_path = source.paths.get("output")
    reasons: list[str] = []
    expected_input_hash = source.provenance.get("input_sha256")
    expected_output_hash = source.provenance.get("output_sha256")
    expected_path_hash = source.provenance.get(
        "path_sha256", source.provenance.get("neb_path_sha256")
    )
    if not input_path or not Path(input_path).is_file():
        reasons.append("rendered_path_input_missing")
        rendered = None
    else:
        rendered = Path(input_path).read_text(encoding="utf-8", errors="strict")
        if existing_file_hash(input_path) != expected_input_hash:
            reasons.append("rendered_path_input_hash_mismatch")
    if not output_path or existing_file_hash(output_path) != expected_output_hash:
        reasons.append("raw_path_output_hash_mismatch")

    path_xyz = source.paths.get("path_xyz") or source.paths.get("neb_path")
    normalized_coordinate: list[float] = []
    if not path_xyz or existing_file_hash(path_xyz) != expected_path_hash:
        reasons.append("raw_path_geometry_hash_mismatch")
    else:
        try:
            normalized_coordinate = normalized_xyz_arc_lengths(
                read_xyz_trajectory(path_xyz)
            )
        except (OSError, TypeError, UnicodeError, ValueError):
            reasons.append("raw_path_geometry_coordinate_unreadable")

    method = {**(source.method or {}), **((source.method or {}).get("settings") or {})}
    engine = str((source.method or {}).get("engine") or "").lower()
    if rendered is None:
        method_audit = {"accepted": False, "checks": {}}
    elif engine == "nwchem":
        method_audit = assess_nwchem_rendered_method(method, rendered)
    else:
        method_audit = {
            "accepted": source.qc.get("input_method_evidence_validated") is True,
            "checks": {},
        }
    if method_audit["accepted"] is not True:
        reasons.append("rendered_input_method_not_validated")
    return {
        "accepted": not reasons,
        "reasons": reasons,
        "source_path_artifact_id": source_id,
        "input_sha256": expected_input_hash,
        "output_sha256": expected_output_hash,
        "path_sha256": expected_path_hash,
        "normalized_path_coordinate": normalized_coordinate,
        "method_checks": method_audit["checks"],
    }


def _endpoint_method_audit(
    manifest: Manifest, member: Artifact
) -> dict[str, Any]:
    source = manifest.find(str(member.data.get("source_path_artifact_id") or ""))
    reaction_id = str(member.data.get("reaction_id") or "")
    if source is None:
        return {"accepted": False, "reasons": ["source_path_artifact_missing"]}
    try:
        _reaction, reactant, product = reaction_and_endpoints(
            manifest, reaction_id
        )
    except KeyError as exc:
        return {"accepted": False, "reasons": [str(exc)]}
    configured = {
        **(source.method or {}),
        **((source.method or {}).get("settings") or {}),
    }
    electronic = evaluate_endpoint_method_lineage(
        reactant,
        product,
        configured,
        str((source.method or {}).get("backend") or ""),
    )
    numerical = evaluate_endpoint_path_numerical_lineage(
        reactant,
        product,
        configured,
    )
    endpoint_input_audits: dict[str, dict[str, Any]] = {}
    endpoint_source_audits: dict[str, dict[str, Any]] = {}
    for role, endpoint in (("reactant", reactant), ("product", product)):
        evidence = (endpoint.provenance or {}).get(
            "validated_numerical_input"
        )
        reasons: list[str] = []
        if not isinstance(evidence, dict):
            evidence = {}
            reasons.append("validated_numerical_input_missing")
        path = evidence.get("path")
        expected_hash = evidence.get("sha256")
        if not path or existing_file_hash(path) != expected_hash:
            reasons.append("validated_numerical_input_hash_mismatch")
        endpoint_input_audits[role] = {
            "accepted": not reasons,
            "reasons": reasons,
            "path": path,
            "sha256": expected_hash,
        }
        endpoint_source_audits[role] = (
            evaluate_endpoint_source_integrity(endpoint)
        )
    endpoint_inputs_accepted = all(
        audit["accepted"] is True
        for audit in endpoint_input_audits.values()
    )
    endpoint_sources_accepted = all(
        audit["accepted"] is True
        for audit in endpoint_source_audits.values()
    )
    return {
        "accepted": electronic["accepted"] is True
        and numerical["accepted"] is True
        and endpoint_inputs_accepted
        and endpoint_sources_accepted,
        "reasons": [
            *electronic["reasons"],
            *numerical["reasons"],
            *[
                f"{role}_{reason}"
                for role, audit in endpoint_input_audits.items()
                for reason in audit["reasons"]
            ],
            *[
                f"{role}_{reason}"
                for role, audit in endpoint_source_audits.items()
                for reason in audit["reasons"]
            ],
        ],
        "electronic_method": electronic,
        "numerical_settings": numerical,
        "numerical_input_audits": endpoint_input_audits,
        "source_integrity_audits": endpoint_source_audits,
    }


class PathEnsembleAssessStage(Stage):
    """Assess existing members without rerunning electronic structure."""

    name = "path-ensemble-assess"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        peak_threshold = float(
            config.get("peak_threshold_kcal_mol", 0.05)
        )
        if peak_threshold <= 0.0:
            raise ValueError("peak_threshold_kcal_mol must be positive")
        studies: defaultdict[tuple[str, str], list[Artifact]] = defaultdict(list)
        for member in manifest.latest_artifacts("path_ensemble_member"):
            if member.qc.get("include_in_assessment") is not True:
                continue
            key = (
                str(member.data.get("reaction_id") or ""),
                str(member.data.get("study_id") or ""),
            )
            if all(key):
                studies[key].append(member)

        records: list[dict[str, Any]] = []
        for (reaction_id, study_id), members in studies.items():
            audited_members: list[Artifact] = []
            audits: dict[str, dict[str, Any]] = {}
            for member in members:
                audit = _input_method_audit(manifest, member)
                endpoint_audit = _endpoint_method_audit(manifest, member)
                audits[member.artifact_id] = {
                    "rendered_input": audit,
                    "endpoint_method_lineage": endpoint_audit,
                }
                revised = member.model_copy(deep=True)
                revised.qc["input_method_evidence_validated"] = audit["accepted"]
                revised.qc["endpoint_method_lineage_validated"] = endpoint_audit[
                    "accepted"
                ]
                source = manifest.find(
                    str(member.data.get("source_path_artifact_id") or "")
                )
                if source is not None:
                    revised.data["input_sha256"] = audit["input_sha256"]
                    revised.data["output_sha256"] = audit["output_sha256"]
                    revised.data["path_sha256"] = audit["path_sha256"]
                    revised.data["method_signature"] = path_method_signature(source)
                path = revised.data.get("path")
                if isinstance(path, dict):
                    revised.data["path"] = reanalyze_reaction_path(
                        path,
                        barrier_threshold_kcal_mol=peak_threshold,
                    ).model_dump(mode="json")
                coordinate = audit["normalized_path_coordinate"]
                coordinate_validated = bool(
                    isinstance(path, dict)
                    and len(coordinate) == len(path.get("images") or [])
                )
                revised.data["normalized_path_coordinate"] = (
                    coordinate if coordinate_validated else []
                )
                revised.data["path_coordinate_definition"] = (
                    "normalized_cumulative_kabsch_aligned_cartesian_arc_length"
                    if coordinate_validated
                    else None
                )
                revised.qc["path_coordinate_validated"] = (
                    coordinate_validated
                )
                revised.provenance["rendered_input_method_audit"] = audit
                revised.provenance["endpoint_method_lineage"] = endpoint_audit
                audited_members.append(revised)
                out.add_artifact(revised)

            assessment = assess_path_ensemble(
                audited_members,
                minimum_resolutions=int(config.get("minimum_resolutions", 2)),
                profile_energy_tolerance_kcal_mol=float(
                    config.get("profile_energy_tolerance_kcal_mol", 0.05)
                ),
            )
            seed_ids = {
                str(member.data.get("refines_path_attempt_id") or "")
                for member in audited_members
                if member.data.get("refines_path_attempt_id")
            }
            assessment.update(
                {
                    "reaction_id": reaction_id,
                    "study_id": study_id,
                    "member_audits": audits,
                    "path_barrier_resolution_kcal_mol": peak_threshold,
                    "supersedes_path_attempt_ids": sorted(
                        seed_ids
                        | {
                            str(value)
                            for value in config.get(
                                "supersedes_path_attempt_ids", []
                            )
                        }
                    ),
                }
            )
            accepted = assessment["accepted"] is True
            artifact = Artifact(
                artifact_id=(
                    f"path_ensemble_assessment_{path_token(study_id, 48)}"
                ),
                artifact_type="path_ensemble_assessment",
                parents=[member.artifact_id for member in audited_members],
                data=assessment,
                qc={
                    "path_discretization_converged": accepted,
                    "barrierless_at_resolution": assessment[
                        "barrierless_at_resolution"
                    ],
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if accepted else "partial",
                    category=None if accepted else "path_ensemble_unresolved",
                    reason="; ".join(assessment["reasons"]),
                ),
            )
            out.add_artifact(artifact)
            records.append(artifact.model_dump(mode="json"))

        write_jsonl(records, out_dir / "path_ensemble_assessments.jsonl")
        out.metadata["path_ensemble_assessment_count"] = len(records)
        return out
