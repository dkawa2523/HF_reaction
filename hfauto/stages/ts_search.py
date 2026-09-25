from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from hfauto.backends.registry import get_ts_engine
from hfauto.chemistry.method_lineage import evaluate_endpoint_method_lineage
from hfauto.chemistry.minima import (
    is_accepted_optimized_minimum,
)
from hfauto.chemistry.reactions import validate_reaction_endpoint_pair
from hfauto.core.artifacts import canonical_species_id, preferred_species_by_id
from hfauto.core.ids import path_token
from hfauto.core.io import (
    ensure_dir,
    read_manifest,
    write_jsonl,
    write_manifest,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.stages.reaction_plan import make_reaction_case_artifact
from hfauto.workflow.reaction_evidence import attempts_by_reaction
from hfauto.workflow.reaction_state import decide_reaction_case
from hfauto.workflow.ts_execution import (
    configured_ts_attempts,
    method_for_reaction_case,
    next_authorized_attempt,
    normalize_ts_result,
)


def _reaction_cases_by_id(manifest: Manifest) -> dict[str, Artifact]:
    return {
        str(artifact.data.get("reaction_id")): artifact
        for artifact in manifest.latest_artifacts("reaction_case")
        if artifact.data.get("reaction_id")
    }


def _next_ts_workdir(base: Path) -> Path:
    """Never overwrite raw evidence from an interrupted path attempt."""

    if not base.exists() or not any(base.iterdir()):
        return base
    retry = 2
    while True:
        candidate = base.with_name(f"{base.name}_retry_{retry:02d}")
        if not candidate.exists() or not any(candidate.iterdir()):
            return candidate
        retry += 1


def _merge_checkpoint(source: Manifest, checkpoint: Manifest | None) -> Manifest:
    """Merge a stage checkpoint with fresh upstream evidence by artifact ID."""

    if checkpoint is None:
        return source
    return Manifest.merge(
        [source, checkpoint],
        run_id=source.run_id,
        stage="ts-search-resume",
        metadata={**source.metadata, **checkpoint.metadata},
    )


def _validated_ts_reaction_ids(
    manifest: Manifest, *, require_real: bool
) -> set[str]:
    """Identify reactions whose TS gate is already complete in a checkpoint."""

    return {
        str(artifact.data.get("reaction_id") or artifact.artifact_id)
        for artifact in manifest.latest_artifacts("reaction_validated")
        if artifact.status.status == "success"
        and artifact.qc.get("ts_validated_by_frequency") is True
        and (
            not require_real
            or (
                artifact.qc.get("real_ts_search_executed") is True
                and artifact.qc.get("real_qm_executed") is True
                and artifact.qc.get("method_evidence_validated") is True
                and artifact.qc.get("fallback_dummy") is False
            )
        )
    }


class TSSearchStage(Stage):
    name = "ts-search"

    def run(
        self, manifest: Manifest | None, config: dict[str, Any], context: StageContext
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        production_mode = str(context.global_config.get("mode", "")).lower() == "production"
        checkpoint_each_attempt = bool(
            config.get("checkpoint_each_attempt", production_mode)
        )
        resume_completed_reactions = bool(
            config.get("resume_completed_reactions", production_mode)
        )
        checkpoint_path = out_dir / "manifest.json"
        checkpoint = (
            read_manifest(checkpoint_path)
            if resume_completed_reactions and checkpoint_path.is_file()
            else None
        )
        working = _merge_checkpoint(manifest, checkpoint)
        out = working.carry_forward(self.name)
        require_validated_dft_minima = bool(
            config.get("require_validated_dft_minima", production_mode)
        )
        require_reaction_case = bool(config.get("require_reaction_case", False))
        reaction_cases = _reaction_cases_by_id(working)
        path_attempts = attempts_by_reaction(working)
        artifact_types = (
            ("species_optimized",)
            if require_validated_dft_minima
            else ("species", "species_preopt", "species_optimized")
        )
        species_by_id = preferred_species_by_id(
            working,
            artifact_types=artifact_types,
            accept=lambda artifact: (
                artifact.artifact_type != "species_optimized"
                or is_accepted_optimized_minimum(artifact)
            ),
        )
        method_base = {
            "method_id": config.get("method", "ts-search"),
            **(config.get("settings", {}) or {}),
            **(config.get("method_settings", {}) or {}),
        }
        engine_settings = config.get("engine_settings", {}) or {}
        records: list[dict] = []
        completed_reactions = _validated_ts_reaction_ids(
            working, require_real=production_mode
        )
        resumed_reactions = 0
        attempts_started = 0
        reactions_started = 0
        stage_started = time.monotonic()
        max_total_attempts = int(config.get("max_total_attempts", 100))
        max_reactions = int(config.get("max_reactions", 100))
        max_stage_walltime_s = float(
            config.get("max_stage_walltime_s", float("inf"))
        )
        if min(max_total_attempts, max_reactions, max_stage_walltime_s) <= 0:
            raise ValueError("TS execution budgets must be positive")

        def checkpoint_now() -> None:
            out.metadata.update(
                {
                    "ts_attempts_started": attempts_started,
                    "ts_reactions_started": reactions_started,
                    "ts_reactions_resumed": resumed_reactions,
                    "ts_execution_budget": {
                        "max_total_attempts": max_total_attempts,
                        "max_reactions": max_reactions,
                        "max_stage_walltime_s": (
                            None
                            if max_stage_walltime_s == float("inf")
                            else max_stage_walltime_s
                        ),
                    },
                }
            )
            write_jsonl(records, out_dir / "ts_search_records.jsonl")
            write_manifest(out, out_dir)

        requested_reaction_ids = {
            str(value) for value in config.get("reaction_ids", [])
        }
        reactions = [
            reaction
            for reaction in working.latest_artifacts("reaction")
            if not requested_reaction_ids
            or str(reaction.data.get("reaction_id") or reaction.artifact_id)
            in requested_reaction_ids
        ]
        for rxn in reactions:
            reaction_id = str(rxn.data.get("reaction_id") or rxn.artifact_id)
            if resume_completed_reactions and reaction_id in completed_reactions:
                resumed_reactions += 1
                continue
            budget_reason = None
            if reactions_started >= max_reactions:
                budget_reason = "max_reactions_reached"
            elif attempts_started >= max_total_attempts:
                budget_reason = "max_total_attempts_reached"
            elif time.monotonic() - stage_started >= max_stage_walltime_s:
                budget_reason = "max_stage_walltime_reached"
            if budget_reason:
                fail = Artifact.failure(
                    f"ts_deferred_{rxn.artifact_id}",
                    "ts_result",
                    "TS calculation was deferred by the declared execution budget",
                    category="execution_budget_exhausted",
                    parents=[rxn.artifact_id],
                    recoverable=True,
                    recommended_fallback=(
                        "resume the same stage with a larger explicit budget; this is "
                        "not evidence about reaction existence"
                    ),
                    reaction_id=reaction_id,
                    budget_reason=budget_reason,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue
            reactions_started += 1
            reaction_case = reaction_cases.get(reaction_id)
            rc_id = rxn.data.get("reactant_species_id")
            ip_id = rxn.data.get("product_species_id")
            rc = species_by_id.get(str(rc_id)) if rc_id else None
            ip = species_by_id.get(str(ip_id)) if ip_id else None

            if rc is None or ip is None:
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "missing_endpoint",
                    category="missing_input",
                    parents=[rxn.artifact_id],
                    recommended_fallback="run dft-minima, minimum-registry, and reaction-plan",
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            execution_reaction = rxn.model_copy(deep=True)
            execution_reaction.data["reactant_species_id"] = (
                canonical_species_id(rc)
            )
            execution_reaction.data["product_species_id"] = (
                canonical_species_id(ip)
            )
            endpoint_qc = validate_reaction_endpoint_pair(
                execution_reaction, rc, ip
            )
            if not endpoint_qc["accepted"]:
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "Reaction endpoints failed the declared chemistry gate",
                    category="invalid_or_collapsed_endpoints",
                    parents=[rxn.artifact_id, rc.artifact_id, ip.artifact_id],
                    recommended_fallback=(
                        "rebuild distinct endpoint minima and verify the declared bond changes"
                    ),
                    endpoint_validation=endpoint_qc,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            if reaction_case is None and require_reaction_case:
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "A reaction-case decision is required before TS search",
                    category="reaction_case_required",
                    parents=[rxn.artifact_id, rc.artifact_id, ip.artifact_id],
                    recommended_fallback="run reaction-plan before ts-search",
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue
            if reaction_case is not None and not bool(
                reaction_case.data.get("ts_search_allowed")
            ):
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "The current reaction-case state does not authorize TS search",
                    category="reaction_case_not_ts_ready",
                    parents=[rxn.artifact_id, reaction_case.artifact_id],
                    recommended_fallback=str(
                        reaction_case.data.get("next_action") or "run reaction-plan"
                    ),
                    reaction_case=reaction_case.data,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            success = False
            last_failure: Artifact | None = None
            attempted_engines: list[str] = []
            attempts = configured_ts_attempts(config, method_base)
            completed_attempts: set[tuple[str, str]] = set()
            execution_index = 0
            attempt_budget_reason: str | None = None
            while selected := next_authorized_attempt(
                attempts, reaction_case, completed_attempts
            ):
                if attempts_started >= max_total_attempts:
                    attempt_budget_reason = "max_total_attempts_reached"
                    break
                if time.monotonic() - stage_started >= max_stage_walltime_s:
                    attempt_budget_reason = "max_stage_walltime_reached"
                    break
                _, attempt, attempt_key = selected
                completed_attempts.add(attempt_key)
                attempts_started += 1
                eng_name = attempt.name
                attempted_engines.append(eng_name)
                method = method_for_reaction_case(
                    attempt.method, reaction_case
                )
                require_method_lineage = bool(
                    require_validated_dft_minima or production_mode
                )
                if require_method_lineage:
                    method_lineage = evaluate_endpoint_method_lineage(
                        rc, ip, method, str(eng_name) if eng_name is not None else None
                    )
                    if not method_lineage["accepted"]:
                        parents = [rxn.artifact_id, rc.artifact_id, ip.artifact_id]
                        fail = Artifact.failure(
                            f"ts_failed_{rxn.artifact_id}",
                            "ts_result",
                            "TS method does not match the validated endpoint PES",
                            category="method_lineage_mismatch",
                            parents=parents,
                            recommended_fallback=(
                                "configure TS engine, functional, basis, dispersion, and "
                                "required program version to match both validated endpoints"
                            ),
                            method_lineage=method_lineage,
                        )
                        out.add_artifact(fail)
                        records.append(fail.model_dump())
                        last_failure = fail
                        continue
                if isinstance(attempt.spec, dict):
                    engine = get_ts_engine(attempt.spec)
                else:
                    engine = get_ts_engine(attempt.spec, **engine_settings)
                # Geometry-only imaginary-mode scoring is a legacy diagnostic,
                # never a production TS validation path.  Non-production callers
                # must opt in explicitly.
                method["allow_legacy_mode_overlap_heuristic"] = (
                    bool(method.get("allow_legacy_mode_overlap_heuristic", False))
                    and not production_mode
                )
                reaction_dir = path_token(
                    rxn.data.get("reaction_id", rxn.artifact_id), max_length=28
                )
                engine_dir = path_token(eng_name, max_length=16)
                workdir = _next_ts_workdir(
                    out_dir
                    / reaction_dir
                    / f"attempt_{execution_index:02d}_{engine_dir}"
                )
                execution_index += 1
                result = engine.search_ts(
                    execution_reaction, rc, ip, method, workdir
                )
                artifacts, record, result_success = normalize_ts_result(
                    result, require_real=production_mode
                )
                for art in artifacts:
                    out.add_artifact(art)
                new_path_attempts = [
                    artifact
                    for artifact in artifacts
                    if artifact.artifact_type
                    in {"path_attempt", "saddle_attempt"}
                ]
                if reaction_case is not None and (
                    new_path_attempts or result_success
                ):
                    path_attempts.setdefault(reaction_id, []).extend(
                        new_path_attempts
                    )
                    revised_case = decide_reaction_case(
                        reaction_id=reaction_id,
                        basin_assessment=dict(
                            reaction_case.data.get("basin_assessment") or {}
                        ),
                        path_attempts=path_attempts[reaction_id],
                        ts_validated=result_success,
                        irc_validated=False,
                        max_path_attempts=int(
                            reaction_case.data.get("max_path_attempts", 3)
                        ),
                    )
                    revised_case_artifact = make_reaction_case_artifact(
                        revised_case,
                        parents=[
                            reaction_case.artifact_id,
                            *[
                                artifact.artifact_id
                                for artifact in artifacts
                                if artifact.artifact_type
                                in {
                                    "path_attempt",
                                    "saddle_attempt",
                                    "reaction_validated",
                                }
                            ],
                        ],
                        created_by=self.name,
                    )
                    out.add_artifact(revised_case_artifact)
                    reaction_case = revised_case_artifact
                if record is not None:
                    records.append(record)
                if checkpoint_each_attempt:
                    checkpoint_now()
                failures = [a for a in artifacts if a.status.status == "failed"]
                if result_success:
                    success = True
                    break
                if failures:
                    last_failure = failures[-1]
            if not success and attempt_budget_reason:
                fail = Artifact.failure(
                    f"ts_deferred_{rxn.artifact_id}",
                    "ts_result",
                    "TS calculation was deferred by the declared execution budget",
                    category="execution_budget_exhausted",
                    parents=[rxn.artifact_id],
                    recoverable=True,
                    recommended_fallback=(
                        "resume the same stage with a larger explicit budget; this is "
                        "not evidence about reaction existence"
                    ),
                    reaction_id=reaction_id,
                    budget_reason=attempt_budget_reason,
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
            configured_names = {attempt.name for attempt in attempts}
            authorized_names = (
                {
                    str(name)
                    for name in reaction_case.data.get(
                        "authorized_path_engines", []
                    )
                }
                if reaction_case is not None
                else configured_names
            )
            missing_current_strategy = not bool(
                configured_names & authorized_names
            )
            workflow_complete = bool(
                reaction_case is not None
                and reaction_case.data.get("workflow_complete") is True
            )
            ts_search_pending = bool(
                reaction_case is None
                or reaction_case.data.get("ts_search_allowed") is True
            )
            if not success and (
                last_failure is None or missing_current_strategy
            ) and not workflow_complete and ts_search_pending and not attempt_budget_reason:
                strategy = (
                    reaction_case.data.get("next_strategy")
                    if reaction_case is not None
                    else None
                )
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "No configured TS engine is authorized for the current reaction state",
                    category="path_strategy_not_supported",
                    parents=[rxn.artifact_id, rc.artifact_id, ip.artifact_id],
                    recommended_fallback=(
                        f"configure an engine for {strategy}"
                        if strategy
                        else "configure one TS-search engine"
                    ),
                    path_strategy=strategy,
                    attempted_engines=attempted_engines,
                    configured_engines=sorted(configured_names),
                    authorized_engines=sorted(authorized_names),
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
            if checkpoint_each_attempt:
                checkpoint_now()
        checkpoint_now()
        return out
