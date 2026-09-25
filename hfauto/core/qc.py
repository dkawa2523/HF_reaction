from __future__ import annotations

from collections.abc import Iterable

from hfauto.core.frequency_qc import (
    DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
    is_significant_imaginary_frequency,
    resolve_imaginary_frequency_cutoff,
)


def quality_tier(has_dft: bool, has_ts: bool, has_irc: bool, has_high_level_sp: bool) -> str:
    if has_irc and has_high_level_sp:
        return "Q5"
    if has_irc:
        return "Q4"
    if has_ts:
        return "Q3"
    if has_dft:
        return "Q2"
    return "Q1"


def confidence_score(
    method_quality: float,
    structure_validity: float,
    path_validity: float,
    ensemble_coverage: float,
    db_calibration_support: float,
) -> float:
    score = (
        0.25 * method_quality
        + 0.25 * structure_validity
        + 0.25 * path_validity
        + 0.15 * ensemble_coverage
        + 0.10 * db_calibration_support
    )
    return max(0.0, min(1.0, float(score)))




TIER_CONFIDENCE_CAP = {
    "Q0": 0.05,
    "Q1": 0.30,
    "Q2": 0.50,
    "Q3": 0.65,
    "Q4": 0.85,
    "Q5": 1.00,
}


def cap_confidence(score: float | None, tier: str | None, has_dummy_or_fallback: bool = False) -> float:
    """Cap confidence by scientific evidence tier.

    This prevents development/dummy screening values from looking as reliable
    as validated DFT/TS/IRC results.  The function is intentionally simple so
    ranking tables remain easy to review.
    """
    try:
        value = float(score if score is not None else 0.0)
    except Exception:
        value = 0.0
    cap = TIER_CONFIDENCE_CAP.get(str(tier or "Q0"), 0.05)
    if has_dummy_or_fallback:
        cap = min(cap, 0.20)
    return max(0.0, min(1.0, min(value, cap)))


def scientific_rank_eligible(tier: str | None, has_dummy_or_fallback: bool = False) -> bool:
    """Eligible for scientific screening tables: real DFT minima or better."""
    return tier_numeric(tier) >= tier_numeric("Q2") and not has_dummy_or_fallback


def production_rank_eligible(
    tier: str | None,
    has_dummy_or_fallback: bool = False,
    production_thermo_ready: bool = False,
    real_irc_executed: bool = False,
) -> bool:
    """Eligible for production ranking: validated TS/IRC plus production thermo."""
    return (
        tier_numeric(tier) >= tier_numeric("Q4")
        and not has_dummy_or_fallback
        and bool(production_thermo_ready)
        and bool(real_irc_executed)
    )


def tier_numeric(tier: str | None) -> float:
    table = {"Q0": 0.0, "Q1": 0.2, "Q2": 0.4, "Q3": 0.6, "Q4": 0.8, "Q5": 1.0}
    return table.get(str(tier or "Q0"), 0.0)


def quality_value(tier: str | None) -> float:
    return tier_numeric(tier)


def summarize_failures(artifacts: Iterable) -> dict[str, int]:
    counts: dict[str, int] = {}
    for artifact in artifacts:
        if getattr(getattr(artifact, "status", None), "status", None) == "failed":
            cat = artifact.status.category or "unknown"
            counts[cat] = counts.get(cat, 0) + 1
    return counts


def minimum_qc_from_frequencies(n_imag: int | None) -> dict:
    if n_imag is None:
        return {"is_minimum": None, "minimum_qc_status": "missing_frequency"}
    return {"is_minimum": int(n_imag) == 0, "minimum_qc_status": "ok" if int(n_imag) == 0 else "imaginary_frequency"}


def minimum_qc_from_freq(n_imag: int | None, state: str | None = None) -> dict:
    if state == "transition_state":
        ok = n_imag == 1
        return {"is_valid_ts_frequency": ok, "expected_n_imag": 1, "observed_n_imag": n_imag, "state": state}
    out = minimum_qc_from_frequencies(n_imag)
    out["state"] = state
    return out


def minimum_promotion_gate(calc, require_real_qm: bool = False) -> tuple[bool, list[str]]:
    """Decide whether a completed optimization is safe to expose as a minimum.

    Backend execution status and scientific interpretation intentionally remain
    separate: a normally completed frequency job with an imaginary mode is a
    valid calculation, but it is not a minimum.
    """
    reasons: list[str] = []
    if getattr(getattr(calc, "status", None), "status", None) != "success":
        reasons.append("calculation_not_successful")
    qc = getattr(calc, "qc", {}) or {}
    data = getattr(calc, "data", {}) or {}
    if qc.get("scf_converged") is not True:
        reasons.append("scf_not_converged")
    if qc.get("geometry_converged") is not True:
        reasons.append("geometry_not_converged")
    n_imag = data.get("n_imag", qc.get("n_imag"))
    if n_imag is None:
        reasons.append("frequency_missing")
    elif int(n_imag) != 0:
        reasons.append("imaginary_frequency")
    if qc.get("geometry_sane") is False:
        reasons.append("geometry_not_sane")
    if require_real_qm:
        if qc.get("fallback_dummy") or qc.get("engine_is_dummy"):
            reasons.append("dummy_or_fallback_result")
        if qc.get("real_qm_executed", qc.get("real_orca_executed")) is not True:
            reasons.append("real_qm_not_executed")
    return not reasons, reasons


def ts_qc(
    n_imag: int | None,
    imag_freq_cm1: float | None,
    mode_overlap: float | None,
    *,
    mode_overlap_threshold: float = 0.50,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> dict:
    """Validate TS frequencies and target-mode projection.

    ``n_imag`` must have been counted with the same numerical-noise cutoff.
    IRC connectivity remains a separate, mandatory production path gate.
    """

    cutoff = resolve_imaginary_frequency_cutoff(
        imaginary_frequency_cutoff_cm1
    )
    frequency_valid = bool(
        n_imag == 1
        and is_significant_imaginary_frequency(imag_freq_cm1, cutoff)
    )
    threshold = float(mode_overlap_threshold)
    if not 0.0 < threshold <= 1.0:
        raise ValueError("mode_overlap_threshold must be in (0, 1]")
    projection_valid = bool((mode_overlap or 0.0) >= threshold)
    valid = frequency_valid and projection_valid
    return {
        "ts_validated_by_frequency": bool(valid),
        "ts_imaginary_frequency_gate_passed": frequency_valid,
        "ts_target_mode_projection_gate_passed": projection_valid,
        "n_imag": n_imag,
        "imag_freq_cm1": imag_freq_cm1,
        "imaginary_frequency_cutoff_cm1": cutoff,
        "mode_overlap_score": mode_overlap,
        "mode_overlap_threshold": threshold,
    }


def process_penalty_from_public_data(public_data: dict | None) -> float:
    public_data = public_data or {}
    gas = public_data.get("gas_process", {}) or {}
    penalty = 0.0
    if str(gas.get("gas_process_feasibility", "")).startswith("poor"):
        penalty += 0.5
    if gas.get("ehs_review_flag") in {True, "manual_review_required", "yes"}:
        penalty += 0.3
    return min(1.0, penalty)
