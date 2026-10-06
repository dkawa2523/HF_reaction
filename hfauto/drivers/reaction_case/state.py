"""Reaction-case state and the pure decision table (design §7.3).

``decide`` evaluates the rows of ``ROWS`` from the top and returns the first decision. The
driver accumulates ``CaseState`` in memory; actions only change the state, never the table.
hypotheses.select closes a hypothesis that its static checks decide (no DFT minimum at an end,
one basin, out of the window), so every row reads the case's own evidence: rows 1-3 complete it
from evidence in hand, every row after them starts a computation or gives up. One stop rule:
every search that gives no answer is one failure token (``CaseState.attempts``), and the case
moves on to its next seed or path while ``max_saddle_attempts`` lasts. An intermediate comes only
from a profile's well (row 6) or a QRC side (row 2). An association
(``ReactionRecord.monomers``) is asked from its separated monomers: SCREEN runs its relaxed scan,
with or without a low-level engine (row 5), and no string runs (row 8: a string joins two
minima, and the monomers' end is none).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, ConfigDict

from hfauto.chemistry.gates import Policy
from hfauto.core.evidence import Evidence, Geometry
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    BarrierVerdict,
    CaseOutcome,
    ConnectionLabel,
    ReactionRecord,
    SaddleClaim,
)


class Action(StrEnum):
    SCREEN = "screen"  # low-level TSs or path and DFT single points, or an association's scan
    REFINE_SADDLE = "refine_saddle"  # initial Hessian (xTB first, then DFT) -> saddle
    VALIDATE_AND_CONNECT = "validate_and_connect"  # TS freq -> gates -> QRC
    FIND_PATH = "find_path"  # DFT string, only without a usable seed or after saddle failure
    VALIDATE_INTERMEDIATE = "validate_intermediate"  # relax_to_minimum -> Registry
    COMPLETE = "complete"


class ReactionPathsPolicy(BaseModel):
    """The budget (pipeline ``policy:``), counts only: max_saddle_attempts failure tokens per
    case (``CaseState.attempts``; a split child starts from 0) and the split depth below a
    hypothesis."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    max_saddle_attempts: int = 2
    max_split_depth: int = 2


@dataclass(frozen=True)
class CaseRules:
    """What decide() and the actions read: the gates, the budget and whether SCREEN runs (the
    stage has a ScreenConfig)."""

    gates: Policy = field(default_factory=Policy)
    budget: ReactionPathsPolicy = field(default_factory=ReactionPathsPolicy)
    screen: bool = True


@dataclass(frozen=True)
class Seed:
    """A saddle search's start. A TS seed (screen_ts, discovery_ts) and a continuation carry
    their own imaginary ``mode``; ``hessian`` is a continuation's DFT freq within HESSIAN_NEAR_A,
    also the search's SCF guess. ``depth`` 1: a continuation, which is not continued again."""

    geometry: Geometry
    source: Literal["discovery_ts", "screen_ts", "screen_hei", "path_hei", "continuation"]
    mode: tuple[float, ...] | None = None
    hessian: Evidence | None = None
    depth: int = 0


@dataclass(frozen=True)
class CaseState:
    """Snapshot accumulated by the driver (never replayed from the log).

    ``screen`` is the verdict of the latest DFT profile (SCREEN's path or scan, or a string),
    which decides rows 3, 6 and 8 and goes into the record. ``shortcut_done``: SCREEN tried its
    shortcut from the low-level TSs (once; it gives seeds, never a verdict). ``neb_done``: the
    low-level path of SCREEN ran (after the shortcut, row 5 runs it once its seeds have failed),
    or an association's scan. ``saddle_pending``: a converged saddle awaits its checks.
    ``connection`` holds only a connected label: a failed action leaves neither claim nor
    connection.
    """

    screen: BarrierVerdict | None = None
    shortcut_done: bool = False
    neb_done: bool = False
    seeds: tuple[Seed, ...] = ()  # unused seeds, consumed from the front
    attempts: int = 0  # failure tokens: every saddle search, and every string without an answer
    saddle_pending: bool = False
    claim: SaddleClaim | None = None
    connection: ConnectionLabel | None = None
    path_runs: int = 0  # DFT strings run, failed ones included
    intermediate: Literal["distinct", "same_as_endpoint", "relax_failed"] | None = None


@dataclass(frozen=True)
class Decision:
    action: Action
    reason: str
    outcome: CaseOutcome | None = None


def record_profile(state: CaseState, verdict: BarrierVerdict, seeds: tuple[Seed, ...] = ()
                   ) -> CaseState:
    """A new DFT profile becomes the latest verdict; its seeds (its peak or SCREEN's NEB TS)
    join only for a single-step profile (design §7.3); its well, if any, is judged anew."""
    seeds = (*state.seeds, *seeds) if verdict.verdict == "single" else state.seeds
    return replace(state, screen=verdict, seeds=seeds, intermediate=None,
                   path_runs=state.path_runs + int(verdict.source == "string"))


Row = Callable[[ReactionRecord, CaseState, CaseRules], Decision | None]


def _complete(outcome: CaseOutcome, reason: str) -> Decision:
    return Decision(Action.COMPLETE, reason, outcome)


def _r01_connected(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.connection in CONNECTED_OUTCOMES:
        return _complete(CONNECTED_OUTCOMES[s.connection], f"connection:{s.connection}")
    return None


def _r02_distinct(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.intermediate == "distinct":  # the driver splits the case (classification.split)
        return _complete(CaseOutcome.MULTI_STEP, "intermediate_distinct")
    return None


def _r03_barrierless(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # A continuous DFT path between the two DFT minima bounds the saddle from above.
    if s.screen is not None and s.screen.verdict == "barrierless":
        return _complete(CaseOutcome.BARRIERLESS, f"{s.screen.source}:barrierless")
    return None


def _r04_saddle(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # a converged saddle (a ts_calc case: from the start) is checked and, as a TS, connected
    return Decision(Action.VALIDATE_AND_CONNECT, "saddle_converged") if s.saddle_pending else None


def _r05_screen(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # the cheap low-level path comes before any string, also once the shortcut's seeds have
    # failed; an association's scan needs no low-level engine
    fire = (p.screen or bool(case.monomers)) and not s.neb_done and not s.seeds
    return Decision(Action.SCREEN, "screen") if fire else None


def _r06_intermediate(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.screen is not None and s.screen.verdict == "intermediate" and s.intermediate is None:
        return Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")  # the lowest well
    return None


def _attempts_left(s: CaseState, p: CaseRules) -> bool:
    return s.attempts < p.budget.max_saddle_attempts


def _r07_seed(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.seeds and _attempts_left(s, p):
        return Decision(Action.REFINE_SADDLE, f"seed:{s.seeds[0].source}")
    return None


def _r08_find_path(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    """A DFT string (none for an association): the first one, or the next chunk from where the
    latest string stopped once every earlier string's seed was tried (a string without an
    answer is a token itself), so a string that left nothing to try ends the case (row 9)."""
    fire = _attempts_left(s, p) and s.path_runs <= s.attempts and not case.monomers
    return Decision(Action.FIND_PATH, "dft_path") if fire else None


def _r09_exhausted(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "attempts_exhausted")


ROWS: tuple[Row, ...] = (
    _r01_connected, _r02_distinct, _r03_barrierless, _r04_saddle, _r05_screen,
    _r06_intermediate, _r07_seed, _r08_find_path, _r09_exhausted,
)


def decide(case: ReactionRecord, state: CaseState, rules: CaseRules) -> Decision:
    """First matching row of the table (pure, no IO); row 9 always matches."""
    return next(d for row in ROWS if (d := row(case, state, rules)) is not None)
