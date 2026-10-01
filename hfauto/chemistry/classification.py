"""Terminal classification of a reaction case, its split into child steps and the record of a
child that is not driven (design §7.3)."""

from __future__ import annotations

from typing import NamedTuple, Protocol

import numpy as np

from hfauto.chemistry.topology import bond_changes
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    BarrierVerdict,
    CaseOutcome,
    ConnectionClaim,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)


class Verdict(Protocol):
    """What finalize reads from ``drivers.reaction_case.state.Decision`` (chemistry sits below
    drivers, so it cannot import the class)."""

    @property
    def outcome(self) -> CaseOutcome | None: ...

    @property
    def reason(self) -> str: ...


class _Undriven(NamedTuple):
    outcome: CaseOutcome
    reason: str


def _oriented(case: ReactionRecord, connected: tuple[str, str]) -> tuple[str, str]:
    """Connected minima, reactant side first when one of them is a hypothesis endpoint."""
    a, b = connected
    return (b, a) if b == case.minima[0] or a == case.minima[1] else (a, b)


def finalize(
    case: ReactionRecord,
    decision: Verdict,
    *,
    barrier: BarrierVerdict | None,
    claim: SaddleClaim | None,
    connection: ConnectionClaim | None,
) -> ReactionRecord:
    """The case with its outcome, reasons and claims attached.

    A connected outcome needs both claims.  A reassigned step adopts the minima the TS actually
    connects (source ``reassigned``), so the declared pair itself stays unanswered.
    """
    outcome = decision.outcome
    if outcome is None:
        raise ValueError(f"{case.reaction_id}: decision {decision.reason!r} is not terminal")
    if outcome in CONNECTED_OUTCOMES.values() and (claim is None or connection is None):
        raise ValueError(f"{case.reaction_id}: {outcome.value} needs saddle and connection claims")
    update: dict[str, object] = {
        "outcome": outcome,
        "reasons": tuple(dict.fromkeys((*case.reasons, decision.reason))),
        "barrier": barrier,
        "saddle": claim,
        "connection": connection,
    }
    if outcome is CaseOutcome.REASSIGNED and connection is not None:
        update |= {"source": "reassigned", "degenerate": False,
                   "minima": _oriented(case, connection.minima)}
    return case.model_copy(update=update)


def undriven(case: ReactionRecord, reason: str, like: ReactionRecord | None = None
             ) -> ReactionRecord:
    """A split child recorded without a drive, so without a job or a log: with ``like``, a
    driven case of its case key (``same_as:<its id>``), that case's outcome and claims, a
    reassigned one adopting the minima its TS connects (``finalize``); else UNRESOLVED
    (``split_depth``)."""
    if like is None or like.outcome is None:
        return finalize(case, _Undriven(CaseOutcome.UNRESOLVED, reason), barrier=None,
                        claim=None, connection=None)
    return finalize(case, _Undriven(like.outcome, reason), barrier=like.barrier,
                    claim=like.saddle, connection=like.connection)


def _child(parent: ReactionRecord, index: int, minima: tuple[str, str],
           endpoints: tuple[str, str], torsional: bool, ts_calc: str | None) -> ReactionRecord:
    return ReactionRecord(
        reaction_id=f"{parent.reaction_id}_split{index}",
        source="split",
        reactants=parent.products,  # a child of an association starts at its complex
        products=parent.products,
        minima=minima,
        endpoints=endpoints,
        torsional=torsional,
        ts_calc=ts_calc,
    )


def split(
    parent: ReactionRecord,
    intermediate: MinimumRecord,
    intermediate_species: SpeciesRecord,
    structures: tuple[np.ndarray, np.ndarray, np.ndarray],
    *,
    ts: tuple[int, str] | None = None,
) -> tuple[ReactionRecord, ReactionRecord]:
    """Fresh child cases R→I and I→P of a multi-step case.

    ``structures`` are the basin structures of R, I and P in one atom order. Stoichiometry is
    the parent's product side (one composition; an association's monomers are not a child's
    end); a child is torsional when its own ends differ in no bond. The
    declared coordinate and low-level TS describe the whole step and are not inherited; ``ts``
    (child 1 or 2, saddle calc) is a TS the parent validated between that child's ends, which
    the child validates first (``ts_calc``) instead of searching again.
    """
    compositions = {t.composition_id for t in (*parent.reactants, *parent.products)}
    if intermediate.composition_id not in compositions:
        raise ValueError(
            f"{parent.reaction_id}: intermediate {intermediate.minimum_id} has composition "
            f"{intermediate.composition_id}, not one of {sorted(compositions)}"
        )
    symbols, (a, i, b) = intermediate_species.geometry.symbols, structures
    first, second = (not any(bond_changes(symbols, x, y)) for x, y in ((a, i), (i, b)))
    m, s = intermediate.minimum_id, intermediate_species.species_id
    calc = dict([ts]) if ts else {}
    return (
        _child(parent, 1, (parent.minima[0], m), (parent.endpoints[0], s), first, calc.get(1)),
        _child(parent, 2, (m, parent.minima[1]), (s, parent.endpoints[1]), second, calc.get(2)),
    )
