"""Evidence-driven workflow state for one molecular reaction hypothesis."""

from __future__ import annotations

from typing import Any

_STRATEGY_FOR_DIAGNOSIS = {
    "endpoint_path_inconsistency": "bracketed_saddle_search",
    "path_stagnation": "reparameterized_double_ended_path",
    "path_unresolved": "adaptive_double_ended_path",
    "saddle_optimizer_failed": "adaptive_double_ended_path",
    "final_frequency_missing": "adaptive_double_ended_path",
    "saddle_order_invalid": "adaptive_double_ended_path",
    "saddle_mode_evidence_missing": "adaptive_double_ended_path",
    "saddle_mode_mismatch": "adaptive_double_ended_path",
    "saddle_energy_order_invalid": "adaptive_double_ended_path",
    "saddle_seed_bracket_failed": "adaptive_double_ended_path",
    "saddle_seed_hessian_rejected": "adaptive_double_ended_path",
    "saddle_frequency_rejected": "adaptive_double_ended_path",
}

_ENGINES_FOR_STRATEGY = {
    "double_ended_path": ["nwchem_neb", "orca_nebts"],
    "adaptive_double_ended_path": ["nwchem_string"],
    "reparameterized_double_ended_path": ["nwchem_string"],
    "bracketed_saddle_search": ["nwchem_saddle"],
    "endpoint_local_saddle_search": ["pysisyphus"],
}


def _attempt_data(attempt: Any) -> dict[str, Any]:
    if hasattr(attempt, "data"):
        return dict(attempt.data or {})
    if hasattr(attempt, "model_dump"):
        return attempt.model_dump()
    return dict(attempt or {})


def _endpoint_next_action(reasons: list[str]) -> str:
    reason_set = set(reasons)
    if "all_endpoint_candidates_collapse_to_same_basin" in reason_set:
        return "expand_discovery_and_rebuild_minimum_registry"
    if "candidate_pair_method_lineage_mismatch" in reason_set:
        return "recompute_endpoints_on_one_pes"
    if "distinct_minima_do_not_match_declared_reaction_chemistry" in reason_set:
        return "revise_reaction_hypothesis"
    if any("frequency_validated" in reason for reason in reasons):
        return "frequency_validate_endpoint_candidates"
    return "audit_basins"


def decide_reaction_case(
    *,
    reaction_id: str,
    basin_assessment: dict[str, Any],
    path_attempts: list[Any] | None = None,
    ts_validated: bool = False,
    irc_validated: bool = False,
    max_path_attempts: int = 3,
) -> dict[str, Any]:
    """Return one deterministic next action; never repeat an identical strategy."""

    attempts = [_attempt_data(attempt) for attempt in (path_attempts or [])]
    strategies = list(
        dict.fromkeys(
            strategy
            for attempt in attempts
            if (strategy := str(attempt.get("strategy") or ""))
        )
    )
    base = {
        "reaction_id": reaction_id,
        "basin_status": basin_assessment.get("status", "unresolved"),
        "basin_assessment": basin_assessment,
        "path_attempt_count": len(strategies),
        "attempted_strategies": strategies,
        "max_path_attempts": int(max_path_attempts),
        "ts_validated": bool(ts_validated),
        "irc_validated": bool(irc_validated),
        "ts_search_allowed": False,
        "irc_allowed": False,
        "workflow_complete": False,
        "scientific_conclusion_supported": False,
        "next_strategy": None,
        "authorized_path_engines": [],
        "reasons": [],
    }
    if irc_validated:
        return {
            **base,
            "state": "validated_elementary_step",
            "next_action": "complete",
            "workflow_complete": True,
            "scientific_conclusion_supported": True,
        }
    if ts_validated:
        return {
            **base,
            "state": "ts_validated",
            "next_action": "validate_connectivity",
            "irc_allowed": True,
        }

    basin_status = str(basin_assessment.get("status") or "unresolved")
    if basin_status == "same_basin":
        return {
            **base,
            "state": "same_basin",
            "next_action": "complete",
            "workflow_complete": True,
            "scientific_conclusion_supported": True,
            "reasons": ["endpoint_pair_does_not_define_an_elementary_step"],
        }
    if basin_status != "distinct_basin":
        basin_reasons = list(basin_assessment.get("reasons") or [])
        return {
            **base,
            "state": "endpoint_evidence_incomplete",
            "next_action": _endpoint_next_action(basin_reasons),
            "reasons": basin_reasons,
        }

    if not attempts:
        strategy = "double_ended_path"
        return {
            **base,
            "state": "path_search_ready",
            "next_action": "search_path",
            "next_strategy": strategy,
            "authorized_path_engines": _ENGINES_FOR_STRATEGY[strategy],
            "ts_search_allowed": True,
        }

    last = attempts[-1]
    diagnosis = str(last.get("diagnosis") or "path_unresolved")
    if diagnosis == "resolved_saddle_candidate":
        candidate = last.get("candidate") or {}
        evidence_validated = bool(
            candidate.get("xyz_path")
            and last.get(
                "search_evidence_validated",
                last.get("all_scan_points_valid", True),
            )
        )
        if evidence_validated:
            return {
                **base,
                "state": "saddle_refinement_ready",
                "next_action": "search_path",
                "next_strategy": "saddle_refinement",
                "authorized_path_engines": [
                    "pysisyphus_saddle",
                    "nwchem_saddle",
                ],
                "ts_search_allowed": True,
                "saddle_seed_candidate": candidate,
                "saddle_seed_evidence_validated": True,
                "reasons": [
                    "resolved_internal_maximum_requires_saddle_refinement"
                ],
            }
        # A profile shape alone is not a geometry seed.  Replan the path once
        # instead of retrying an empty saddle-refinement request.
        diagnosis = "path_unresolved"
    if diagnosis == "multiple_step_candidate":
        return {
            **base,
            "state": "intermediate_basin_validation_required",
            "next_action": "validate_intermediate_basins",
            "reasons": ["path_contains_multiple_internal_maxima"],
        }
    recovery = dict(last.get("saddle_recovery") or {})
    if (
        diagnosis == "saddle_optimizer_failed"
        and recovery.get("next_action") == "refresh_hessian_and_restart"
        and recovery.get("restart_allowed") is True
    ):
        candidate = {
            "xyz_path": recovery.get("restart_seed_xyz"),
            "source": "failed_saddle_final_geometry",
        }
        return {
            **base,
            "state": "saddle_hessian_refresh_ready",
            "next_action": "search_path",
            "next_strategy": "saddle_refinement",
            "authorized_path_engines": ["nwchem_saddle"],
            "ts_search_allowed": True,
            "saddle_seed_candidate": candidate,
            "saddle_seed_evidence_validated": bool(candidate["xyz_path"]),
            "method_updates": dict(recovery.get("method_updates") or {}),
            "reasons": ["nonconverged_saddle_requires_fresh_hessian"],
        }
    if diagnosis == "final_frequency_missing":
        return {
            **base,
            "state": "saddle_frequency_completion_ready",
            "next_action": "run_fixed_geometry_frequency",
            "fixed_geometry_xyz": recovery.get("restart_seed_xyz"),
            "reasons": ["stationary_geometry_requires_final_frequency_analysis"],
        }
    if diagnosis == "same_basin":
        return {
            **base,
            "state": "same_basin",
            "next_action": "complete",
            "workflow_complete": True,
            "scientific_conclusion_supported": True,
            "reasons": ["path_attempt_confirms_same_basin"],
        }

    strategy = _STRATEGY_FOR_DIAGNOSIS.get(diagnosis, "adaptive_double_ended_path")
    if strategy in strategies:
        strategy = None
    budget_exhausted = len(strategies) >= int(max_path_attempts)
    if strategy is None or budget_exhausted:
        return {
            **base,
            "state": "unresolved_within_budget",
            "next_action": "complete",
            "workflow_complete": True,
            "reasons": [
                "path_attempt_budget_exhausted"
                if budget_exhausted
                else "no_untried_strategy_for_diagnosis"
            ],
        }
    authorized_engines = _ENGINES_FOR_STRATEGY[strategy]
    if not authorized_engines:
        return {
            **base,
            "state": "path_strategy_implementation_required",
            "next_action": "implement_path_strategy",
            "next_strategy": strategy,
            "reasons": [
                f"previous_path_diagnosis:{diagnosis}",
                f"no_engine_implements:{strategy}",
            ],
        }
    return {
        **base,
        "state": "path_replan_ready",
        "next_action": "search_path",
        "next_strategy": strategy,
        "authorized_path_engines": authorized_engines,
        "ts_search_allowed": True,
        "adaptive_path_seed_candidate": last.get("candidate"),
        "initial_path_xyz": (last.get("evidence") or {}).get("path_xyz"),
        "reasons": [
            f"previous_path_diagnosis:{diagnosis}",
            *(
                ["resolved_maximum_missing_validated_geometry_seed"]
                if str(last.get("diagnosis")) == "resolved_saddle_candidate"
                else []
            ),
        ],
    }
