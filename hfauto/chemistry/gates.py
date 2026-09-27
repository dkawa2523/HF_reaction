"""Chemical pass/fail gates (design §5.4): the only place that turns Evidence into verdicts.

Allowed imports: the standard library, numpy, hfauto.core and chemistry.profile.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

from hfauto.chemistry.profile import classify
from hfauto.core.constants import CM1_TO_HARTREE, HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Level
from hfauto.core.records import BarrierVerdict, CaseOutcome, ReactionRecord, ReactionThermo


@dataclass(frozen=True)
class Gate:
    ok: bool
    reasons: tuple[str, ...] = ()  # why it failed
    notes: tuple[str, ...] = ()  # remarks attached to a pass

    def __bool__(self) -> bool:
        return self.ok


@dataclass(frozen=True)
class Policy:
    """The only thresholds a pipeline YAML may override (``gates:``)."""

    noise_cm1: float = 10.0
    saddle_cm1: float = 50.0
    scf_noise_factor: float = 20.0
    qrc_min_drop_hartree: float = 1.0e-5
    resolution_kcal: float = 1.0  # hills and wells of a DFT path count from this depth
    reaction_window_kcal: float = 40.0
    thermo_zpe_tol_hartree: float = 1.0e-5
    thermo_energy_tol_hartree: float = 1.0e-6
    spin_contamination_tol: float = 0.1


_DEFAULT = Policy()
_PES_FIELDS = (
    "program", "version", "method", "basis", "dispersion", "solvation",
    "charge", "multiplicity", "electronic_temperature_K",
)
_STATE_FIELDS = ("charge", "multiplicity")
_NUMERICS_FIELDS = ("grid", "scf_tol")
RANKABLE_OUTCOMES = frozenset({
    CaseOutcome.ELEMENTARY_STEP, CaseOutcome.DEGENERATE, CaseOutcome.REASSIGNED,
    CaseOutcome.BARRIERLESS,
})

ConnectionLabel = Literal["elementary", "degenerate", "reassigned", "failed"]


def _gate(reasons: Sequence[str], notes: Sequence[str] = ()) -> Gate:
    return Gate(ok=not reasons, reasons=tuple(reasons), notes=tuple(notes))


def qrc_drop(level: Level, policy: Policy = _DEFAULT) -> float:
    """Smallest energy drop a QRC side must show below the TS: max(minimum, scf_noise_factor x
    scf_tol); a missing scf_tol counts as 0."""
    return max(policy.qrc_min_drop_hartree, policy.scf_noise_factor * (level.scf_tol or 0.0))


def zpe_hartree(freqs_cm1: Sequence[float], *, scale: float = 1.0) -> float:
    """0.5 x scale x sum(nu) over nu > 0."""
    return 0.5 * scale * sum(nu for nu in freqs_cm1 if nu > 0.0) * CM1_TO_HARTREE


def imaginary_tier(
    freqs_cm1: Sequence[float], policy: Policy = _DEFAULT
) -> Literal["none", "noise", "soft", "saddle"]:
    if not freqs_cm1:
        return "none"
    lowest = min(freqs_cm1)
    if lowest < -policy.saddle_cm1:
        return "saddle"
    if lowest < -policy.noise_cm1:
        return "soft"
    return "noise" if lowest < 0.0 else "none"


def same_pes(*levels: Level, numerics: bool = True, state: bool = True) -> Gate:
    """state=False skips charge and multiplicity (a complex against its monomers)."""
    fields = _PES_FIELDS + (_NUMERICS_FIELDS if numerics else ())
    reasons = [
        f"pes_mismatch:{name}" for name in fields
        if (state or name not in _STATE_FIELDS)
        and len({getattr(level, name) for level in levels}) > 1
    ]
    return _gate(reasons)


def spin_ok(ev: Evidence, policy: Policy = _DEFAULT) -> Gate:
    if ev.s2 is None:
        return Gate(True)
    spin = (ev.level.multiplicity - 1) / 2
    if abs(ev.s2 - spin * (spin + 1)) <= policy.spin_contamination_tol:
        return Gate(True)
    return Gate(False, ("spin_contaminated",))


def _mode_count_reasons(freq: Evidence) -> list[str]:
    if freq.frequencies_cm1 is None or freq.n_external is None:
        return ["missing_frequencies"]
    expected = 3 * len(freq.final.symbols) - freq.n_external
    if len(freq.frequencies_cm1) != expected:
        return [f"frequency_count:{len(freq.frequencies_cm1)}!={expected}"]
    return []


def _link_reasons(freq: Evidence, parent: Evidence) -> list[str]:
    """The freq job must sit on the parent's final geometry, on the same PES incl. numerics."""
    reasons = [] if freq.start.fingerprint == parent.final.fingerprint else ["geometry_mismatch"]
    return reasons + list(same_pes(parent.level, freq.level, numerics=True).reasons)


def is_minimum(freq: Evidence, *, opt: Evidence, policy: Policy = _DEFAULT) -> Gate:
    reasons = [f"not_{want}_task:{ev.task}" for ev, want in ((freq, "freq"), (opt, "opt"))
               if ev.task != want]
    reasons += _link_reasons(freq, opt) + _mode_count_reasons(freq)
    tier = imaginary_tier(freq.frequencies_cm1 or (), policy)
    notes: tuple[str, ...] = ()
    if tier == "saddle":
        reasons.append("imaginary_mode")
    elif tier != "none":
        notes = (f"{tier}_imaginary_mode",)
    return _gate(reasons, notes)


def is_first_order_saddle(freq: Evidence, *, saddle: Evidence, policy: Policy = _DEFAULT) -> Gate:
    """One negative eigenvalue of any size: the lowest mode below -noise_cm1, the second not
    below -saddle_cm1 (higher_order); a second between the two is the note soft_secondary_mode."""
    reasons = _link_reasons(freq, saddle) + _mode_count_reasons(freq)
    lowest, second = (*sorted(freq.frequencies_cm1 or ()), 0.0, 0.0)[:2]  # missing: 0
    if lowest >= -policy.noise_cm1:
        reasons.append("no_imaginary_mode")
    if second < -policy.saddle_cm1:
        reasons.append("higher_order")
    soft = -policy.saddle_cm1 <= second < -policy.noise_cm1
    return _gate(reasons, ("soft_secondary_mode",) if soft else ())


def barrier_verdict(
    energies: Sequence[float],
    *,
    source: Literal["screen", "string"],
    max_node_spacing_A: float | None = None,
    policy: Policy = _DEFAULT,
) -> BarrierVerdict:
    """Class of a DFT profile with the two DFT minima energies at its ends, at resolution_kcal:
    a continuous path's maximum bounds the saddle from above."""
    if len(energies) < 3:
        return BarrierVerdict(verdict="unavailable", source=source, reasons=("too_few_points",))
    rel = (max(energies[1:-1]) - max(energies[0], energies[-1])) * HARTREE_TO_KCAL_MOL
    verdict = classify(energies, policy.resolution_kcal / HARTREE_TO_KCAL_MOL)
    return BarrierVerdict(verdict=verdict, source=source, max_rel_kcal=rel,
                          max_node_spacing_A=max_node_spacing_A)


def _side_reasons(index: int, side: Evidence, ts: Evidence, drop: float) -> list[str]:
    """One QRC side on the TS's PES ends below E_TS - drop; its path there is not judged."""
    reasons = [f"side{index}:{r}" for r in same_pes(ts.level, side.level, numerics=True).reasons]
    if side.energy_hartree >= ts.energy_hartree - drop:
        reasons.append(f"side{index}:no_descent")
    return reasons


def _assignment(
    assigned: tuple[str | None, str | None],
    expected: frozenset[str],
    *,
    degenerate: bool,
    sides_distinct: bool,
) -> tuple[ConnectionLabel, str | None]:
    first, second = assigned
    if first is None or second is None:
        return "failed", "unassigned_side"
    if degenerate and first == second and {first} == expected:
        return ("degenerate", None) if sides_distinct else ("failed", "sides_not_distinct")
    if first == second:
        return "failed", "sides_same_basin"
    if not degenerate and {first, second} == expected:
        return "elementary", None
    return "reassigned", None


def connection(
    ts_freq: Evidence,
    sides: tuple[Evidence, Evidence],
    assigned: tuple[str | None, str | None],
    expected: frozenset[str],
    *,
    degenerate: bool,
    sides_distinct: bool = True,
    policy: Policy = _DEFAULT,
) -> tuple[Gate, ConnectionLabel]:
    drop = qrc_drop(ts_freq.level, policy)
    reasons = [r for i, side in enumerate(sides) for r in _side_reasons(i, side, ts_freq, drop)]
    label, why = _assignment(
        assigned, expected, degenerate=degenerate, sides_distinct=sides_distinct
    )
    if why is not None:
        reasons.append(why)
    if reasons:
        return _gate(reasons), "failed"
    return Gate(True), label


def thermo_consistent(
    freq: Evidence,
    *,
    frequencies_cm1: Sequence[float],
    gv_zpe_hartree: float,
    gv_energy_hartree: float,
    gv_n_real: int,
    scale: float,
    policy: Policy = _DEFAULT,
) -> Gate:
    """GoodVibes' output against the frequencies it was given (chemistry.thermo_frequencies)."""
    n_real = sum(1 for nu in frequencies_cm1 if nu > 0.0)
    reasons = [] if gv_n_real == n_real else [f"n_real:{gv_n_real}!={n_real}"]
    reference = zpe_hartree(frequencies_cm1, scale=scale)
    if not abs(gv_zpe_hartree - reference) < policy.thermo_zpe_tol_hartree:
        reasons.append("zpe_mismatch")
    if not abs(gv_energy_hartree - freq.energy_hartree) < policy.thermo_energy_tol_hartree:
        reasons.append("energy_mismatch")
    return _gate(reasons)


def rankable(reaction: ReactionRecord, thermo: ReactionThermo | None) -> Gate:
    """A rankable outcome and a thermo record without blockers. The thermo stage is the only
    source of the blockers (thermo_unavailable, mixed_level_of_theory, spin_contaminated)."""
    outcome = () if reaction.outcome in RANKABLE_OUTCOMES else (f"outcome:{reaction.outcome}",)
    facts = thermo.blockers if thermo is not None else ("thermo_record_missing",)
    return _gate([*outcome, *facts])


Tier = Literal["screening", "minima", "saddle", "connected"]


def reaction_tier(reaction: ReactionRecord) -> Tier:
    """connected > saddle > minima > screening.

    Both endpoints are DFT minima on one PES once the case driver has produced any outcome
    other than BLOCKED (decision row 1), so the record alone decides the minima tier.
    """
    if reaction.connection is not None:
        return "connected"
    if reaction.saddle is not None:
        return "saddle"
    if reaction.outcome is not None and reaction.outcome != CaseOutcome.BLOCKED:
        return "minima"
    return "screening"
