"""Reaction-case state and the pure decision table (design §7.3).

``decide`` evaluates the 13 rows of ``ROWS`` from the top and returns the first decision.  The
driver accumulates ``CaseState`` in memory; actions only change the state, never the table.
Rows 4-6 complete a case from evidence already in hand; every row after them starts a
computation or gives up, bounded by counts only (``ReactionPathsPolicy``). A saddle without an
imaginary mode is a failed attempt like any rejected saddle: an intermediate comes only from a
profile's well (row 10) or a QRC side (row 5). An association (``ReactionRecord.monomers``) is
asked from its separated monomers: SCREEN runs its relaxed scan, with or without a low-level
engine (row 9), and no string runs (row 12: a string joins two minima, and the monomers' end is
none).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, ConfigDict

from hfauto.chemistry.gates import Policy
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Geometry
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    BarrierVerdict,
    CaseOutcome,
    ConnectionLabel,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
)


class Action(StrEnum):
    SCREEN = "screen"  # low-level path and DFT single points, or an association's scan
    REFINE_SADDLE = "refine_saddle"  # initial Hessian (xTB first, then DFT) -> saddle
    VALIDATE_AND_CONNECT = "validate_and_connect"  # TS freq -> gates -> QRC (two amplitudes)
    FIND_PATH = "find_path"  # DFT string, only without a usable seed or after saddle failure
    VALIDATE_INTERMEDIATE = "validate_intermediate"  # relax_to_minimum -> Registry
    COMPLETE = "complete"
    BLOCKED = "blocked"


class ReactionPathsPolicy(BaseModel):
    """The budget (pipeline ``policy:``), counts only: max_saddle_attempts per case (a split
    child starts from 0; a stalled search's restart is part of its attempt) and the split depth
    below a hypothesis."""

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
    """A saddle search's start. Only a TS seed (screen_ts, discovery_ts, higher_order_retry)
    carries its own imaginary ``mode``; ``hessian`` is a verified TS freq within HESSIAN_NEAR_A,
    also the search's SCF guess. ``depth``: the continuations since a fresh seed (a push, or a
    restart from the last frame, which keeps the source but neither mode nor Hessian, adds
    one)."""

    geometry: Geometry
    source: Literal["discovery_ts", "screen_ts", "screen_hei", "path_hei", "higher_order_retry"]
    mode: tuple[float, ...] | None = None
    hessian: Evidence | None = None
    depth: int = 0


@dataclass(frozen=True)
class CaseState:
    """Snapshot accumulated by the driver (never replayed from the log).

    ``minima`` (not in the §7.3 listing) holds the registry records of ``case.minima``, None when
    an endpoint has none: rows 1-3 need their basin and energy, which a ReactionRecord does not
    carry. ``screen`` is the verdict of the latest DFT profile (SCREEN or string), which decides
    rows 6, 10 and 12 and goes into the record. ``neb_done``: the low-level path of SCREEN ran
    (after a shortcut, row 9 runs it once its seeds have failed), or an association's scan.
    ``last_saddle``: a converged saddle awaits its checks, or the latest search or its checks
    failed (None once a TS is accepted: its claim stands for it).
    """

    minima: tuple[MinimumRecord | None, MinimumRecord | None] = (None, None)
    screen: BarrierVerdict | None = None
    neb_done: bool = False
    seeds: tuple[Seed, ...] = ()  # unused seeds, consumed from the front
    saddle_attempts: int = 0
    last_saddle: Literal["converged", "failed"] | None = None
    claim: SaddleClaim | None = None
    # a rejected QRC: sides in one basin at both amplitudes, or in one case key but two basins
    connection: ConnectionLabel | Literal["same_basin", "same_state"] | None = None
    path_runs: int = 0  # DFT strings run, failed ones included
    intermediate: Literal["distinct", "same_as_endpoint", "relax_failed"] | None = None


@dataclass(frozen=True)
class Decision:
    action: Action
    reason: str
    outcome: CaseOutcome | None = None


def record_profile(state: CaseState, verdict: BarrierVerdict, seeds: tuple[Seed, ...] = ()
                   ) -> CaseState:
    """A new DFT profile becomes the latest verdict; its seeds (its peak, or the low-level TSs
    of SCREEN's shortcut) join only for a single-step profile (design §7.3), and the checks of an
    earlier saddle no longer apply."""
    seeds = (*state.seeds, *seeds) if verdict.verdict == "single" else state.seeds
    return replace(state, screen=verdict, seeds=seeds, last_saddle=None, intermediate=None,
                   path_runs=state.path_runs + int(verdict.source == "string"))


Row = Callable[[ReactionRecord, CaseState, CaseRules], Decision | None]


def _complete(outcome: CaseOutcome, reason: str) -> Decision:
    return Decision(Action.COMPLETE, reason, outcome)


def _r01_endpoints(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    """An end without a DFT minimum (a declared end that relaxed to a saddle) blocks the case.
    Both ends lie on one PES by construction: PipelineConfig._check_stages gives the DFT minima
    and reaction-paths stages one method, and the stage reads only DFT minima."""
    a, b = s.minima
    if a is None or b is None:
        return Decision(Action.BLOCKED, "endpoint_without_dft_minimum", CaseOutcome.BLOCKED)
    return None


def _r02_same_basin(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    a, b = s.minima
    # an association's reactant end is its separated monomers, never the adduct's basin
    if a and b and a.basin_id == b.basin_id and not case.degenerate and not case.monomers:
        return _complete(CaseOutcome.SAME_BASIN, "same_basin")
    return None


def _r03_window(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    a, b = s.minima
    if a and b and (b.energy_hartree - a.energy_hartree) * HARTREE_TO_KCAL_MOL > (
        p.gates.reaction_window_kcal
    ):
        return _complete(CaseOutcome.OUT_OF_WINDOW, "out_of_window")
    return None


def _r04_connected(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.connection in CONNECTED_OUTCOMES:
        return _complete(CONNECTED_OUTCOMES[s.connection], f"connection:{s.connection}")
    return None


def _r05_distinct(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.intermediate == "distinct":  # the driver splits the case (classification.split)
        return _complete(CaseOutcome.MULTI_STEP, "intermediate_distinct")
    return None


def _r06_barrierless(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # A continuous DFT path between the two DFT minima bounds the saddle from above.
    if s.screen is not None and s.screen.verdict == "barrierless":
        return _complete(CaseOutcome.BARRIERLESS, f"{s.screen.source}:barrierless")
    return None


def _r07_connection(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.connection in ("same_basin", "same_state") and _attempts_left(s, p):
        return None  # a saddle of another process (sides in one basin or key): search on (9-12)
    return None if s.connection is None else _complete(CaseOutcome.UNRESOLVED, "connection_failed")


def _r08_saddle(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # a converged saddle (a ts_calc case: from the start) is checked and, as a TS, connected
    if s.last_saddle == "converged":
        return Decision(Action.VALIDATE_AND_CONNECT, "saddle_converged")
    return None


def _r09_screen(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # the cheap low-level path comes before any string, also once a shortcut's seeds have
    # failed; an association's scan needs no low-level engine
    fire = ((p.screen or bool(case.monomers)) and not s.neb_done
            and (s.screen is None or not s.seeds))
    return Decision(Action.SCREEN, "screen") if fire else None


def _r10_intermediate(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.screen is not None and s.screen.verdict == "intermediate" and s.intermediate is None:
        return Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")  # the lowest well
    return None


def _attempts_left(s: CaseState, p: CaseRules) -> bool:
    return s.saddle_attempts < p.budget.max_saddle_attempts


def _r11_seed(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.seeds and _attempts_left(s, p):
        return Decision(Action.REFINE_SADDLE, f"seed:{s.seeds[0].source}")
    return None


def _r12_find_path(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    """A DFT string (none for an association): the first one, or the next chunk from where the
    latest string stopped once its seeds have failed. A string runs only while its peak seed can
    still be refined, and only after every earlier string's seed was tried: a string that left
    no seed ends the case (row 13)."""
    v, runs = s.screen, s.path_runs
    fire = (_attempts_left(s, p) and runs <= s.saddle_attempts and not case.monomers
            and (not runs or (v is not None and v.source == "string"
                              and v.verdict != "unavailable")))
    return Decision(Action.FIND_PATH, "dft_path") if fire else None


def _r13_exhausted(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "attempts_exhausted")


ROWS: tuple[Row, ...] = (
    _r01_endpoints, _r02_same_basin, _r03_window, _r04_connected, _r05_distinct,
    _r06_barrierless, _r07_connection, _r08_saddle, _r09_screen, _r10_intermediate,
    _r11_seed, _r12_find_path, _r13_exhausted,
)


def decide(case: ReactionRecord, state: CaseState, rules: CaseRules) -> Decision:
    """First matching row of the table (pure, no IO); row 13 always matches."""
    return next(d for row in ROWS if (d := row(case, state, rules)) is not None)
