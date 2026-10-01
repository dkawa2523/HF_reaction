"""Chemical pass/fail gates (design §5.4): the only place that turns Evidence into verdicts.

Allowed imports: the standard library, numpy, hfauto.core and chemistry.profile.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from hfauto.chemistry.profile import classify
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Level
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    BarrierVerdict,
    CaseOutcome,
    ConnectionLabel,
    ProfileSource,
    ReactionRecord,
    ReactionThermo,
)


@dataclass(frozen=True)
class Gate:
    ok: bool
    reasons: tuple[str, ...] = ()  # why it failed
    notes: tuple[str, ...] = ()  # remarks attached to a pass

    def __bool__(self) -> bool:
        return self.ok


@dataclass(frozen=True)
class Policy:
    """The chemical thresholds a pipeline YAML may override (``gates:``); nothing else."""

    noise_cm1: float = 10.0
    saddle_cm1: float = 50.0
    resolution_kcal: float = 1.0  # hills and wells of a DFT path count from this depth
    reaction_window_kcal: float = 40.0
    spin_tol: float = 0.1  # |<S2> - S(S+1)|


_DEFAULT = Policy()
LOW_LEVEL_BARRIER_MAX_KCAL = 50.0  # explore's xTB barrier cap; the DFT stages judge the rest
SCF_NOISE_FACTOR = 20.0  # x scf_tol: the SCF noise of an energy difference at one geometry
QRC_MIN_DROP_HARTREE = 1.0e-5  # the floor of that noise (H2Te freq vs opt: 3.1e-6 Eh)
REACTION_MODE_MIN = 0.3  # in the empty χ gap: saddles no case connected ≤ 0.02, TSs ≥ 0.57
_PES_FIELDS = (
    "program", "version", "method", "basis", "dispersion",
    "charge", "multiplicity", "electronic_temperature_K",
)
_STATE_FIELDS = ("charge", "multiplicity")
_NUMERICS_FIELDS = ("grid", "scf_tol")
RANKABLE_OUTCOMES = frozenset({*CONNECTED_OUTCOMES.values(), CaseOutcome.BARRIERLESS})
Bonds = frozenset[tuple[int, int]]  # a labelled bond graph (topology.bonds)


def _gate(reasons: Sequence[str], notes: Sequence[str] = ()) -> Gate:
    return Gate(ok=not reasons, reasons=tuple(reasons), notes=tuple(notes))


def qrc_drop(level: Level) -> float:
    """The SCF noise floor max(QRC_MIN_DROP_HARTREE, SCF_NOISE_FACTOR x scf_tol); a missing
    scf_tol counts as 0. A QRC side must end this far below the TS, and a freq job may differ
    this much from its parent at one geometry."""
    return max(QRC_MIN_DROP_HARTREE, SCF_NOISE_FACTOR * (level.scf_tol or 0.0))


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
    if abs(ev.s2 - spin * (spin + 1)) <= policy.spin_tol:
        return Gate(True)
    return Gate(False, ("spin_contaminated",))


def _link_reasons(freq: Evidence, parent: Evidence) -> list[str]:
    """The freq job sits on the parent's final geometry, on the same PES incl. numerics, in the
    same SCF solution: an energy within qrc_drop of the parent's (a UKS freq may find another
    solution; NWChem's freq tightens the grid and screening, H2Te 3.1e-6 Eh apart)."""
    if freq.start.fingerprint != parent.final.fingerprint:
        reasons = ["geometry_mismatch"]
    elif abs(freq.energy_hartree - parent.energy_hartree) > qrc_drop(parent.level):
        reasons = ["state_mismatch"]
    else:
        reasons = []
    return reasons + list(same_pes(parent.level, freq.level, numerics=True).reasons)


def is_minimum(freq: Evidence, *, opt: Evidence, policy: Policy = _DEFAULT) -> Gate:
    reasons = [f"not_{want}_task:{ev.task}" for ev, want in ((freq, "freq"), (opt, "opt"))
               if ev.task != want]
    reasons += _link_reasons(freq, opt)
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
    reasons = _link_reasons(freq, saddle)
    lowest, second = (*sorted(freq.frequencies_cm1 or ()), 0.0, 0.0)[:2]  # missing: 0
    if lowest >= -policy.noise_cm1:
        reasons.append("no_imaginary_mode")
    if second < -policy.saddle_cm1:
        reasons.append("higher_order")
    soft = -policy.saddle_cm1 <= second < -policy.noise_cm1
    return _gate(reasons, ("soft_secondary_mode",) if soft else ())


def reaction_mode_chi(mode: Sequence[float], coords: np.ndarray, bonds: Bonds,
                      gradient: np.ndarray | None = None) -> float | None:
    """χ = ‖Q_Bᵀq̂‖ of the unit Cartesian imaginary mode q̂ (mass weighting removed): Q_B is an
    orthonormal basis of the Wilson stretch vectors ∂r_ij/∂x of the hypothesis' changed bonds at
    the saddle ``coords`` (SVD: a set of stretches can be linearly dependent). Without a changed
    bond, |cos(q̂, ∇q)| of a declared coordinate q; None with neither."""
    q = np.ravel(np.asarray(mode, dtype=float))
    q = q / np.linalg.norm(q)
    if bonds:
        x, b = np.reshape(coords, (-1, 3)), np.zeros((q.size, len(bonds)))
        for k, (i, j) in enumerate(sorted(bonds)):
            e = (x[i] - x[j]) / np.linalg.norm(x[i] - x[j])
            b[3 * i:3 * i + 3, k], b[3 * j:3 * j + 3, k] = e, -e
        u, s, _ = np.linalg.svd(b, full_matrices=False)
        return float(np.linalg.norm(u[:, s > 1e-8 * s[0]].T @ q))
    if gradient is None or not np.any(gradient):
        return None
    return float(abs(q @ np.ravel(gradient)) / np.linalg.norm(gradient))


def reaction_mode_character(mode: Sequence[float], coords: np.ndarray, bonds: Bonds,
                            gradient: np.ndarray | None = None) -> Gate:
    """A first-order saddle is this hypothesis' TS only when its imaginary mode carries the
    hypothesis' bond change (or declared coordinate): χ ≥ REACTION_MODE_MIN, else
    ``not_reaction_mode`` (a reorientation or rotor saddle). Not applied without either."""
    chi = reaction_mode_chi(mode, coords, bonds, gradient)
    return _gate(("not_reaction_mode",) if chi is not None and chi < REACTION_MODE_MIN else ())


def barrier_verdict(
    energies: Sequence[float], *, source: ProfileSource, policy: Policy = _DEFAULT
) -> BarrierVerdict:
    """Class of a DFT profile with the two DFT minima energies at its ends (an association's
    scan: the separated monomers' sum first), at resolution_kcal: a continuous path's maximum
    bounds the saddle from above."""
    if len(energies) < 3:
        return BarrierVerdict(verdict="unavailable", source=source, reasons=("too_few_points",))
    verdict = classify(energies, policy.resolution_kcal / HARTREE_TO_KCAL_MOL)
    return BarrierVerdict(verdict=verdict, source=source)


def _side_reasons(index: int, side: Evidence, ts: Evidence, drop: float) -> list[str]:
    """One QRC side on the TS's PES ends below E_TS - drop; its path there is not judged."""
    reasons = [f"side{index}:{r}" for r in same_pes(ts.level, side.level, numerics=True).reasons]
    if side.energy_hartree >= ts.energy_hartree - drop:
        reasons.append(f"side{index}:no_descent")
    return reasons


def _bonds_exchanged(bond_sets: tuple[Bonds, ...] | None) -> bool:
    """A degenerate step that changes bonds must show the two endpoints' labelled bond graphs on
    its QRC sides; one without a bond change (inversion, torsion) is not judged by bonds."""
    if bond_sets is None:
        return True
    reactant, product, *sides = bond_sets
    return reactant == product or set(sides) == {reactant, product}


def _assignment(
    assigned: tuple[str | None, str | None],
    expected: frozenset[str],
    *,
    degenerate: bool,
    sides_distinct: bool,
    bonds_exchanged: bool,
) -> tuple[ConnectionLabel, str | None]:
    first, second = assigned
    if first is None or second is None:
        return "failed", "unassigned_side"
    if degenerate and first == second and {first} == expected:
        if not sides_distinct:
            return "failed", "sides_not_distinct"
        return ("degenerate", None) if bonds_exchanged else ("failed", "bond_change_missing")
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
    bond_sets: tuple[Bonds, ...] | None = None,
) -> tuple[Gate, ConnectionLabel]:
    """Both QRC sides end below the TS in assigned basins. ``bond_sets`` are the labelled bond
    graphs (R, P, side0, side1): a degenerate case whose ends differ in bonds needs exactly those
    two graphs on its sides (a methyl rotation TS leaves the transferred proton in place)."""
    drop = qrc_drop(ts_freq.level)
    reasons = [r for i, side in enumerate(sides) for r in _side_reasons(i, side, ts_freq, drop)]
    label, why = _assignment(
        assigned, expected, degenerate=degenerate, sides_distinct=sides_distinct,
        bonds_exchanged=_bonds_exchanged(bond_sets),
    )
    if why is not None:
        reasons.append(why)
    if reasons:
        return _gate(reasons), "failed"
    return Gate(True), label


def discovery_verdict(
    *,
    ts_validated: bool,
    dE_act_kcal: float | None,
    dE_rxn_kcal: float | None,
    policy: Policy = _DEFAULT,
) -> str | None:
    """None when a low-level NT2 product is kept, else the negative reason. It needs its
    validated TS and IRC; the reaction window is the DFT one (a narrower low-level sieve
    only adds false negatives) and an unevaluated reaction energy lies outside it."""
    if not ts_validated:
        return "ts_not_validated"
    if dE_rxn_kcal is None or dE_rxn_kcal > policy.reaction_window_kcal:
        return "out_of_window"
    if dE_act_kcal is not None and dE_act_kcal > LOW_LEVEL_BARRIER_MAX_KCAL:
        return "out_of_window"
    return None


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
