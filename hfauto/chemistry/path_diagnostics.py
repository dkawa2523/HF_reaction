"""Method-neutral diagnosis and routing for reaction-path attempts."""

from __future__ import annotations

import re
from pathlib import Path
from statistics import median
from typing import Any

from hfauto.core.schemas.path import (
    PathAttemptRecord,
    PathOptimizationRecord,
    ReactionPathRecord,
)

_VALUES = {
    "iterations": re.compile(r"#\s*NEB Path iteration\s*=\s*(\d+)", re.IGNORECASE),
    "max_gradient": re.compile(r"#\s*Gmax\s*=\s*([-+0-9.Ee]+)", re.IGNORECASE),
    "rms_gradient": re.compile(r"#\s*Grms\s*=\s*([-+0-9.Ee]+)", re.IGNORECASE),
    "max_displacement": re.compile(r"#\s*Xmax\s*=\s*([-+0-9.Ee]+)", re.IGNORECASE),
    "rms_displacement": re.compile(r"#\s*Xrms\s*=\s*([-+0-9.Ee]+)", re.IGNORECASE),
}

_STRING_STEP = re.compile(
    r"^\s*@zts\s+(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+",
    re.IGNORECASE | re.MULTILINE,
)
_STRING_TOLERANCE = re.compile(
    r"^\s*@zts\s+Co(?:n)?vergence\s+Tolerance\s*=\s*([-+0-9.Ee]+)",
    re.IGNORECASE | re.MULTILINE,
)


def assess_optimizer_stagnation(
    rms_step_history: list[float],
    *,
    convergence_target: float,
    window: int = 20,
    far_from_target_factor: float = 50.0,
    minimum_fractional_improvement: float = 0.20,
) -> dict[str, Any]:
    """Diagnose a stalled optimizer from unit-consistent RMS step evidence.

    Two complete windows are required.  A run is considered stalled only when
    every recent step remains far from the convergence target and the recent
    median improved by less than the requested fraction.  The function is
    independent of reaction coordinates and electronic-structure programs.
    """

    if convergence_target <= 0.0:
        raise ValueError("convergence_target must be positive")
    if window < 2:
        raise ValueError("window must be at least 2")
    if far_from_target_factor <= 1.0:
        raise ValueError("far_from_target_factor must be greater than 1")
    if not 0.0 <= minimum_fractional_improvement < 1.0:
        raise ValueError("minimum_fractional_improvement must be in [0, 1)")

    values = [float(value) for value in rms_step_history]
    enough_history = len(values) >= 2 * window
    if not enough_history:
        return {
            "stagnant": False,
            "enough_history": False,
            "sample_count": len(values),
            "window": window,
        }

    previous = values[-2 * window : -window]
    recent = values[-window:]
    previous_median = median(previous)
    recent_median = median(recent)
    fractional_improvement = (
        (previous_median - recent_median) / previous_median
        if previous_median > 0.0
        else 0.0
    )
    far_from_target = min(recent) > convergence_target * far_from_target_factor
    stagnant = bool(
        far_from_target
        and fractional_improvement < minimum_fractional_improvement
    )
    return {
        "stagnant": stagnant,
        "enough_history": True,
        "sample_count": len(values),
        "window": window,
        "convergence_target": convergence_target,
        "far_from_target_factor": far_from_target_factor,
        "minimum_fractional_improvement": minimum_fractional_improvement,
        "previous_median": previous_median,
        "recent_minimum": min(recent),
        "recent_median": recent_median,
        "recent_maximum": max(recent),
        "fractional_improvement": fractional_improvement,
        "far_from_target": far_from_target,
    }


def parse_neb_optimization_history(
    source: str | Path | None,
    *,
    converged: bool,
    stagnation_window: int = 3,
    stagnation_relative_range: float = 0.10,
    stagnation_displacement_A: float = 1.0e-4,
) -> PathOptimizationRecord:
    """Parse NWChem's NEB history without treating a small RMS as convergence."""

    if source is None:
        text = ""
    else:
        candidate = Path(source)
        try:
            text = candidate.read_text(encoding="utf-8", errors="ignore") if candidate.is_file() else str(source)
        except OSError:
            text = str(source)
    iterations = [int(value) for value in _VALUES["iterations"].findall(text)]
    gradients = [float(value) for value in _VALUES["max_gradient"].findall(text)]
    rms_gradients = [float(value) for value in _VALUES["rms_gradient"].findall(text)]
    displacements = [float(value) for value in _VALUES["max_displacement"].findall(text)]
    rms_displacements = [
        float(value) for value in _VALUES["rms_displacement"].findall(text)
    ]
    recent = gradients[-max(2, int(stagnation_window)) :]
    scale = max(abs(sum(recent) / len(recent)), 1.0e-15) if recent else 0.0
    relative_range = (max(recent) - min(recent)) / scale if recent else None
    displacement_small = bool(
        displacements and displacements[-1] <= float(stagnation_displacement_A)
    )
    stagnant = bool(
        not converged
        and len(recent) >= max(2, int(stagnation_window))
        and relative_range is not None
        and relative_range <= float(stagnation_relative_range)
        and displacement_small
    )
    return PathOptimizationRecord(
        iterations=max(iterations) if iterations else None,
        max_gradient=gradients[-1] if gradients else None,
        rms_gradient=rms_gradients[-1] if rms_gradients else None,
        max_displacement=displacements[-1] if displacements else None,
        rms_displacement=rms_displacements[-1] if rms_displacements else None,
        max_gradient_history=gradients,
        converged=bool(converged),
        stagnant=stagnant,
    )


def parse_string_optimization_history(
    source: str | Path | None,
    *,
    converged: bool,
) -> PathOptimizationRecord:
    """Parse NWChem zero-temperature-string step evidence.

    NWChem reports string ``xrms`` and ``xmax`` but not the projected maximum
    gradient in the compact ``@zts`` rows.  Missing quantities stay ``None``;
    they are never inferred from displacements.
    """

    if source is None:
        text = ""
    else:
        candidate = Path(source)
        try:
            text = (
                candidate.read_text(encoding="utf-8", errors="ignore")
                if candidate.is_file()
                else str(source)
            )
        except OSError:
            text = str(source)
    rows = [
        (int(step), float(xrms), float(xmax))
        for step, xrms, xmax in _STRING_STEP.findall(text)
    ]
    tolerance_match = _STRING_TOLERANCE.search(text)
    tolerance = (
        float(tolerance_match.group(1)) if tolerance_match else 4.5e-4
    )
    xrms_history = [xrms for _step, xrms, _xmax in rows]
    window = 3
    stagnant = False
    if len(xrms_history) >= 2 * window:
        previous_best = min(xrms_history[-2 * window : -window])
        recent_best = min(xrms_history[-window:])
        best_improvement = (
            (previous_best - recent_best) / previous_best
            if previous_best > 0.0
            else 0.0
        )
        stagnant = bool(
            recent_best > tolerance * 5.0 and best_improvement < 0.20
        )
    return PathOptimizationRecord(
        iterations=rows[-1][0] if rows else None,
        max_displacement=rows[-1][2] if rows else None,
        rms_displacement=rows[-1][1] if rows else None,
        converged=bool(converged),
        stagnant=bool(not converged and stagnant),
    )


def parse_path_optimization_history(
    engine: str,
    source: str | Path | None,
    *,
    converged: bool,
) -> PathOptimizationRecord:
    """Dispatch preserved raw history without leaking backend logic to stages."""

    if "string" in str(engine).lower():
        return parse_string_optimization_history(source, converged=converged)
    return parse_neb_optimization_history(source, converged=converged)


def diagnose_path_attempt(
    path: ReactionPathRecord | dict[str, Any],
    endpoint_basin_status: str,
    optimization: PathOptimizationRecord | dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Choose the next scientific action from path and endpoint evidence."""

    path_record = (
        path
        if isinstance(path, ReactionPathRecord)
        else ReactionPathRecord.model_validate(path)
    )
    optimization_record = (
        optimization
        if isinstance(optimization, PathOptimizationRecord)
        else PathOptimizationRecord.model_validate(optimization or {})
    )
    reasons: list[str] = []
    if not path_record.converged:
        reasons.append("path_not_converged")
    if optimization_record.stagnant:
        reasons.append("optimizer_stagnant")

    if endpoint_basin_status == "same_basin":
        return {
            "diagnosis": "same_basin",
            "next_action": "close_without_elementary_step",
            "reasons": [*reasons, "endpoints_assign_to_same_basin"],
        }
    if (
        path_record.classification == "monotonic_no_internal_maximum"
        and endpoint_basin_status == "distinct_basin"
    ):
        return {
            "diagnosis": "endpoint_path_inconsistency",
            "next_action": "bracket_narrow_saddle",
            "reasons": [
                *reasons,
                "frequency_validated_distinct_minima_but_path_is_monotonic",
            ],
        }
    if len(path_record.internal_maximum_indices) > 1:
        return {
            "diagnosis": "multiple_step_candidate",
            "next_action": "locate_and_validate_intermediate_basins",
            "reasons": [*reasons, "multiple_resolved_internal_maxima"],
        }
    if path_record.classification == "resolved_internal_maximum":
        return {
            "diagnosis": "resolved_saddle_candidate",
            "next_action": "refine_saddle",
            "reasons": (
                reasons
                if path_record.converged
                else [
                    *reasons,
                    "unconverged_path_energy_is_seed_only_not_a_barrier",
                ]
            ),
        }
    if not path_record.converged:
        return {
            "diagnosis": (
                "path_stagnation" if optimization_record.stagnant else "path_unresolved"
            ),
            "next_action": (
                "reparameterize_path"
                if optimization_record.stagnant
                else "adaptive_path_refinement"
            ),
            "reasons": reasons,
        }
    return {
        "diagnosis": "path_unresolved",
        "next_action": "adaptive_path_refinement",
        "reasons": [*reasons, "path_has_no_resolved_saddle_candidate"],
    }


def make_path_attempt_record(
    *,
    attempt_id: str,
    reaction_id: str,
    strategy: str,
    engine: str,
    endpoint_basin_status: str,
    path: ReactionPathRecord | dict[str, Any],
    optimization: PathOptimizationRecord | dict[str, Any] | None = None,
    evidence: dict[str, Any] | None = None,
) -> PathAttemptRecord:
    path_record = (
        path
        if isinstance(path, ReactionPathRecord)
        else ReactionPathRecord.model_validate(path)
    )
    optimization_record = (
        optimization
        if isinstance(optimization, PathOptimizationRecord)
        else PathOptimizationRecord.model_validate(optimization or {})
    )
    decision = diagnose_path_attempt(
        path_record, endpoint_basin_status, optimization_record
    )
    return PathAttemptRecord(
        attempt_id=attempt_id,
        reaction_id=reaction_id,
        strategy=strategy,
        engine=engine,
        endpoint_basin_status=endpoint_basin_status,
        path=path_record,
        optimization=optimization_record,
        diagnosis=decision["diagnosis"],
        next_action=decision["next_action"],
        reasons=decision["reasons"],
        evidence=dict(evidence or {}),
    )
