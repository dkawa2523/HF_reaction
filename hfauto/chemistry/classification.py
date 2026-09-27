"""Terminal classification of a reaction case and its split into child steps (design §7.3)."""

from __future__ import annotations

from typing import Protocol

from hfauto.core.records import (
    BarrierVerdict,
    CaseOutcome,
    ConnectionClaim,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)

_CONNECTED = frozenset(
    {CaseOutcome.ELEMENTARY_STEP, CaseOutcome.DEGENERATE, CaseOutcome.REASSIGNED}
)


class Verdict(Protocol):
    """What finalize reads from ``drivers.reaction_case.state.Decision`` (chemistry sits below
    drivers, so it cannot import the class)."""

    @property
    def outcome(self) -> CaseOutcome | None: ...

    @property
    def reason(self) -> str: ...


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
    if outcome in _CONNECTED and (claim is None or connection is None):
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


def _child(parent: ReactionRecord, index: int, minima: tuple[str, str],
           endpoints: tuple[str, str]) -> ReactionRecord:
    return ReactionRecord(
        reaction_id=f"{parent.reaction_id}_split{index}",
        source="split",
        reactants=parent.reactants,
        products=parent.products,
        minima=minima,
        endpoints=endpoints,
        torsional=parent.torsional,
    )


def split(
    parent: ReactionRecord, intermediate: MinimumRecord, intermediate_species: SpeciesRecord
) -> tuple[ReactionRecord, ReactionRecord]:
    """Fresh child cases R→I and I→P of a multi-step case.

    Stoichiometry and ``torsional`` are inherited from the parent; the declared coordinate
    and low-level TS describe the whole step and are not.
    """
    compositions = {t.composition_id for t in (*parent.reactants, *parent.products)}
    if intermediate.composition_id not in compositions:
        raise ValueError(
            f"{parent.reaction_id}: intermediate {intermediate.minimum_id} has composition "
            f"{intermediate.composition_id}, not one of {sorted(compositions)}"
        )
    m, s = intermediate.minimum_id, intermediate_species.species_id
    return (
        _child(parent, 1, (parent.minima[0], m), (parent.endpoints[0], s)),
        _child(parent, 2, (m, parent.minima[1]), (s, parent.endpoints[1])),
    )
