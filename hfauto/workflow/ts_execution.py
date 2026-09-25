"""Small orchestration helpers shared by transition-state workflow stages."""

from __future__ import annotations

from dataclasses import dataclass
from json import dumps
from typing import Any

from hfauto.core.schemas.artifact import Artifact


@dataclass(frozen=True)
class TSEngineAttempt:
    """One configured backend and its backend-independent method settings."""

    name: str
    spec: str | dict[str, Any]
    method: dict[str, Any]


def configured_ts_attempts(
    config: dict[str, Any], method_base: dict[str, Any]
) -> list[TSEngineAttempt]:
    """Normalize string/dict engine syntax once, preserving configured order."""

    specs = list(config.get("engine_order") or [config.get("engine", "dummy")])
    attempts: list[TSEngineAttempt] = []
    for spec in specs:
        if isinstance(spec, dict):
            name = spec.get("engine") or spec.get("name")
            method = {**method_base, **(spec.get("settings", {}) or {})}
        else:
            name = spec
            method = dict(method_base)
        if not name:
            raise ValueError("TS engine specification requires a name")
        attempts.append(TSEngineAttempt(str(name), spec, method))
    return attempts


def attempt_is_authorized(
    attempt: TSEngineAttempt, reaction_case: Artifact | None
) -> bool:
    """Select only engines authorized by the current scientific state."""

    if reaction_case is None:
        return True
    strategy = reaction_case.data.get("next_strategy")
    authorized = {
        str(name) for name in reaction_case.data.get("authorized_path_engines", [])
    }
    return bool(strategy and attempt.name in authorized)


def ts_attempt_key(
    attempt: TSEngineAttempt, reaction_case: Artifact | None
) -> tuple[str, str]:
    """Identify an exact engine/scientific-state input to prevent retry loops."""

    if reaction_case is None:
        return attempt.name, dumps(
            {"state": "unplanned", "method": attempt.method},
            sort_keys=True,
            default=str,
        )
    data = reaction_case.data
    state_input = {
        "state": data.get("state"),
        "next_strategy": data.get("next_strategy"),
        "path_attempt_count": data.get("path_attempt_count"),
        "method_updates": data.get("method_updates"),
        "saddle_seed_candidate": data.get("saddle_seed_candidate"),
        "adaptive_path_seed_candidate": data.get(
            "adaptive_path_seed_candidate"
        ),
        "initial_path_xyz": data.get("initial_path_xyz"),
        "method": attempt.method,
    }
    return attempt.name, dumps(state_input, sort_keys=True, default=str)


def next_authorized_attempt(
    attempts: list[TSEngineAttempt],
    reaction_case: Artifact | None,
    completed: set[tuple[str, str]],
) -> tuple[int, TSEngineAttempt, tuple[str, str]] | None:
    """Pick the first configured engine valid for the *current* state."""

    for configured_index, attempt in enumerate(attempts):
        key = ts_attempt_key(attempt, reaction_case)
        if key not in completed and attempt_is_authorized(attempt, reaction_case):
            return configured_index, attempt, key
    return None


def method_for_reaction_case(
    base: dict[str, Any], reaction_case: Artifact | None
) -> dict[str, Any]:
    """Translate one reaction-case decision into backend method inputs."""

    method = dict(base)
    if reaction_case is None:
        return method
    data = reaction_case.data
    method.update(
        {
            "path_strategy": data.get("next_strategy"),
            "endpoint_basin_status": data.get("basin_status", "unresolved"),
            "basin_assessment": dict(data.get("basin_assessment") or {}),
            "saddle_seed_candidate": data.get("saddle_seed_candidate"),
            "saddle_seed_evidence_validated": data.get(
                "saddle_seed_evidence_validated"
            ),
            "adaptive_path_seed_candidate": data.get(
                "adaptive_path_seed_candidate"
            ),
            "initial_path_xyz": data.get("initial_path_xyz"),
            "reaction_case_artifact_id": reaction_case.artifact_id,
        }
    )
    method.update(dict(data.get("method_updates") or {}))
    return method


def has_validated_ts(
    artifacts: list[Artifact], *, require_real: bool = False
) -> bool:
    """Backend completion is not equivalent to a frequency-validated TS."""

    return any(
        artifact.artifact_type == "reaction_validated"
        and artifact.status.status == "success"
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
        for artifact in artifacts
    )


def normalize_ts_result(
    result: Any, *, require_real: bool = False
) -> tuple[list[Artifact], Any, bool]:
    """Normalize legacy backend result shapes at one compatibility boundary."""

    artifacts, record = normalize_artifact_result(result)
    return artifacts, record, has_validated_ts(
        artifacts, require_real=require_real
    )


def normalize_artifact_result(result: Any) -> tuple[list[Artifact], Any]:
    """Normalize legacy result containers until every backend uses dataclasses."""

    if hasattr(result, "artifacts"):
        artifacts = list(result.artifacts)
        record = getattr(result, "record", None)
    elif isinstance(result, dict):
        artifacts = list(result.get("artifacts", []) or [])
        record = result.get("record")
    else:
        artifacts = list(result or [])
        record = None
    return artifacts, record
