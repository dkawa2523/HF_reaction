"""Run and validate bidirectional IRC for registry-promoted reactions."""

from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_ts_engine
from hfauto.chemistry.method_lineage import (
    evaluate_artifact_method_lineage,
    evaluate_endpoint_method_lineage,
)
from hfauto.chemistry.reaction_path_qc import endpoint_pair_match_qc
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.stages.reaction_plan import make_reaction_case_artifact
from hfauto.workflow.reaction_inputs import validated_ts_and_endpoints
from hfauto.workflow.reaction_state import decide_reaction_case
from hfauto.workflow.ts_execution import normalize_artifact_result


def _reaction_cases_by_id(manifest: Manifest) -> dict[str, Artifact]:
    return {
        str(artifact.data.get("reaction_id")): artifact
        for artifact in manifest.latest_artifacts("reaction_case")
        if artifact.data.get("reaction_id")
    }


def _endpoint_pair_qc(
    irc: Artifact,
    ts_species: Artifact,
    reactant: Artifact,
    product: Artifact,
    method: dict[str, Any],
) -> dict[str, Any]:
    """Recheck backend endpoints against the two registered DFT minima."""

    forward = (
        irc.paths.get("forward_xyz")
        or irc.paths.get("endpoint_a_xyz")
        or irc.data.get("forward_endpoint_xyz")
    )
    backward = (
        irc.paths.get("backward_xyz")
        or irc.paths.get("endpoint_b_xyz")
        or irc.data.get("backward_endpoint_xyz")
    )
    reactant_xyz = (
        reactant.data.get("xyz_path")
        or reactant.paths.get("xyz")
        or reactant.paths.get("final_xyz")
    )
    product_xyz = (
        product.data.get("xyz_path")
        or product.paths.get("xyz")
        or product.paths.get("final_xyz")
    )
    return endpoint_pair_match_qc(
        forward,
        backward,
        reactant_xyz,
        product_xyz,
        ts_species.data,
        float(method.get("endpoint_rmsd_threshold_A", 0.75)),
        q_tolerance_A=float(method.get("endpoint_q_tolerance_A", 0.30)),
        require_identity_invariant_geometry=bool(
            method.get("require_identity_invariant_endpoint_match", True)
        ),
        permutation_rmsd_threshold_A=float(
            method.get("endpoint_permutation_rmsd_threshold_A", 0.20)
        ),
        distance_spectrum_threshold_A=float(
            method.get("endpoint_distance_spectrum_threshold_A", 0.08)
        ),
    )


def _failure(
    reaction: Artifact,
    category: str,
    reason: str,
    *,
    parents: list[str] | None = None,
    **evidence: Any,
) -> Artifact:
    return Artifact.failure(
        f"irc_failed_{reaction.artifact_id}",
        "irc",
        reason,
        category=category,
        parents=parents or [reaction.artifact_id],
        recommended_fallback=(
            "verify the registered minima, TS method lineage, and both IRC branches"
        ),
        **evidence,
    )


class IRCStage(Stage):
    name = "irc"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        production = str(context.global_config.get("mode", "")).lower() == "production"
        require_case = bool(config.get("require_reaction_case", production))
        cases = _reaction_cases_by_id(manifest)
        engine_spec = config.get("engine", "pysisyphus")
        method = {
            "method_id": config.get("method", "irc"),
            **(config.get("settings", {}) or {}),
            **(config.get("method_settings", {}) or {}),
        }
        if production:
            method["allow_legacy_rmsd_only_irc_validation"] = False
            method["require_identity_invariant_endpoint_match"] = True
        records: list[dict[str, Any]] = []

        reactions = list(manifest.latest_artifacts("reaction"))
        known_ids = {
            str(item.data.get("reaction_id") or item.artifact_id)
            for item in reactions
        }
        for validated in manifest.latest_artifacts("reaction_validated"):
            validated_id = str(
                validated.data.get("reaction_id") or validated.artifact_id
            )
            if validated_id not in known_ids:
                reactions.append(validated)
                known_ids.add(validated_id)

        for reaction in reactions:
            reaction_id = str(
                reaction.data.get("reaction_id") or reaction.artifact_id
            )
            case = cases.get(reaction_id)
            if require_case and (
                case is None or case.data.get("irc_allowed") is not True
            ):
                fail = _failure(
                    reaction,
                    "reaction_case_not_irc_ready",
                    "The current reaction-case state does not authorize IRC",
                    parents=[
                        reaction.artifact_id,
                        *([case.artifact_id] if case is not None else []),
                    ],
                    reaction_case=(case.data if case is not None else None),
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue
            try:
                (
                    _,
                    ts_species,
                    reactant,
                    product,
                    validated_ts,
                ) = validated_ts_and_endpoints(manifest, reaction_id)
            except KeyError as exc:
                fail = _failure(
                    reaction,
                    "validated_ts_or_endpoint_missing",
                    str(exc),
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            engine_name = (
                str(engine_spec.get("engine") or engine_spec.get("name"))
                if isinstance(engine_spec, dict)
                else str(engine_spec)
            )
            endpoint_lineage = evaluate_endpoint_method_lineage(
                reactant,
                product,
                method,
                engine_name,
            )
            ts_lineage = evaluate_artifact_method_lineage(
                validated_ts,
                endpoint_lineage,
            )
            if not endpoint_lineage["accepted"] or not ts_lineage["accepted"]:
                fail = _failure(
                    reaction,
                    "method_lineage_mismatch",
                    "IRC calculator, registered minima, and TS are not on one PES",
                    parents=[
                        reaction.artifact_id,
                        reactant.artifact_id,
                        product.artifact_id,
                        validated_ts.artifact_id,
                    ],
                    endpoint_method_lineage=endpoint_lineage,
                    ts_method_lineage=ts_lineage,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            engine = (
                get_ts_engine(engine_spec)
                if isinstance(engine_spec, dict)
                else get_ts_engine(engine_spec, **(config.get("engine_settings") or {}))
            )
            result = engine.run_irc(
                reaction,
                ts_species,
                reactant,
                product,
                method,
                out_dir / reaction_id,
            )
            artifacts, record = normalize_artifact_result(result)
            for artifact in artifacts:
                out.add_artifact(artifact)
            if record is not None:
                records.append(record)

            successful = [
                artifact
                for artifact in artifacts
                if artifact.artifact_type == "irc"
                and artifact.status.status == "success"
            ]
            primary = successful[-1] if successful else None
            if primary is None:
                continue
            independent_qc = _endpoint_pair_qc(
                primary,
                ts_species,
                reactant,
                product,
                method,
            )
            irc_validated = bool(
                independent_qc.get("irc_validated") is True
                and primary.qc.get("irc_validated") is True
                and primary.qc.get("method_evidence_validated") is True
                and primary.qc.get("fallback_dummy") is False
                and (
                    not production
                    or primary.qc.get("real_irc_executed") is True
                )
            )
            primary.qc["independent_endpoint_pair_qc"] = independent_qc
            primary.qc["irc_validated"] = irc_validated
            if not irc_validated:
                fail = _failure(
                    reaction,
                    "irc_path_not_validated",
                    "Both IRC branches did not recover the registered endpoint basins",
                    parents=[reaction.artifact_id, primary.artifact_id],
                    independent_endpoint_pair_qc=independent_qc,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            validated_path = Artifact(
                artifact_id=f"validated_path_{reaction_id}",
                artifact_type="reaction_path_validated",
                parents=[
                    reaction.artifact_id,
                    validated_ts.artifact_id,
                    primary.artifact_id,
                    reactant.artifact_id,
                    product.artifact_id,
                ],
                data={
                    "reaction_id": reaction_id,
                    "ts_species_id": ts_species.artifact_id,
                    "irc_artifact_id": primary.artifact_id,
                    "reactant_species_id": reaction.data["reactant_species_id"],
                    "product_species_id": reaction.data["product_species_id"],
                },
                method=primary.method,
                qc={
                    "ts_validated_by_frequency": True,
                    "irc_validated": True,
                    "real_irc_executed": primary.qc.get("real_irc_executed"),
                    "method_evidence_validated": True,
                    "fallback_dummy": False,
                    "independent_endpoint_pair_qc": independent_qc,
                },
                provenance={"created_by": self.name},
            )
            out.add_artifact(validated_path)
            if case is not None:
                revised_case = decide_reaction_case(
                    reaction_id=reaction_id,
                    basin_assessment=dict(case.data.get("basin_assessment") or {}),
                    path_attempts=[],
                    ts_validated=True,
                    irc_validated=True,
                    max_path_attempts=int(case.data.get("max_path_attempts", 3)),
                )
                out.add_artifact(
                    make_reaction_case_artifact(
                        revised_case,
                        parents=[case.artifact_id, validated_path.artifact_id],
                        created_by=self.name,
                    )
                )

        write_jsonl(records, out_dir / "irc_records.jsonl")
        return out
