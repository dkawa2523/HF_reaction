"""Reaction-case state and the pure decision table (design §7.3).

``decide`` evaluates the 16 rows of ``ROWS`` from the top and returns the first decision.  The
driver accumulates ``CaseState`` in memory; actions only change the state, never the table.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from enum import StrEnum
from typing import Literal

from hfauto.chemistry.gates import Policy
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Geometry
from hfauto.core.records import (
    BarrierVerdict,
    CaseOutcome,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
)

ConnectionLabel = Literal["elementary", "degenerate", "reassigned", "failed"]
QRC_AMPLITUDES = 2  # QRC tries at most two amplitudes (§8.1)


class Action(StrEnum):
    SCREEN = "screen"  # low-level path and DFT single points: barrier pre-check
    REFINE_SADDLE = "refine_saddle"  # initial Hessian (xTB first, then DFT) -> saddle
    VALIDATE_TS = "validate_ts"  # separate DFT freq -> is_first_order_saddle
    FIND_PATH = "find_path"  # DFT string, only without a usable seed or after saddle failure
    CONNECT = "connect"  # QRC: displace, optimize, assign, connection gate
    VALIDATE_INTERMEDIATE = "validate_intermediate"  # relax_to_minimum -> Registry
    COMPLETE = "complete"
    BLOCKED = "blocked"


@dataclass(frozen=True)
class CasePolicy:
    """Budget and defaults of one reaction case; the only place they are set."""

    screen: bool = True
    max_saddle_attempts: int = 2
    string_beads: int = 9
    string_chunks: int = 3
    screen_images: int = 11  # GS nodes including both ends
    qrc_target_hartree: float = 3.0e-4  # displacement energy target: max(3 x drop, this)
    qrc_bounds_A: tuple[float, float] = (0.05, 0.4)
    qrc_retry_factor: float = 2.0
    walltime_s: float = 6 * 3600
    max_split_depth: int = 2
    gates: Policy = field(default_factory=Policy)


@dataclass(frozen=True)
class Seed:
    geometry: Geometry
    source: Literal["discovery_ts", "screen_ts", "screen_hei", "path_hei", "higher_order_retry"]
    tangent: tuple[float, ...] | None  # path tangent or mapped endpoint difference


@dataclass(frozen=True)
class CaseState:
    """Snapshot accumulated by the driver (never replayed from the log).

    ``minima`` (not in the §7.3 listing) holds the registry records of ``case.minima``, None when
    an endpoint has none: rows 1-3 need their tier, level, basin and energy, which a
    ReactionRecord does not carry. ``screen`` is the verdict of the latest DFT profile (SCREEN
    or string), which decides rows 12-13 and goes into the record.
    """

    minima: tuple[MinimumRecord | None, MinimumRecord | None] = (None, None)
    expired: bool = False
    screen: BarrierVerdict | None = None
    seeds: tuple[Seed, ...] = ()  # unused seeds, consumed from the front
    saddle_attempts: int = 0
    last_saddle: Literal["converged", "failed"] | None = None
    ts_check: Literal["ok", "collapsed", "higher_order"] | None = None
    claim: SaddleClaim | None = None
    connection: ConnectionLabel | None = None
    connection_attempts: int = 0
    path_runs: int = 0  # DFT strings run, failed ones included
    intermediate: Literal["distinct", "same_as_endpoint"] | None = None


@dataclass(frozen=True)
class Decision:
    action: Action
    reason: str
    outcome: CaseOutcome | None = None


def record_profile(state: CaseState, verdict: BarrierVerdict, seed: Seed | None = None
                   ) -> CaseState:
    """A new DFT profile becomes the latest verdict; its peak seeds only a single-step profile
    (eng P0-4), and the checks of an earlier saddle no longer apply."""
    seeds = (*state.seeds, seed) if verdict.verdict == "single" and seed else state.seeds
    return replace(state, screen=verdict, seeds=seeds, last_saddle=None, ts_check=None,
                   intermediate=None, path_runs=state.path_runs + int(verdict.source == "string"))


Row = Callable[[ReactionRecord, CaseState, CasePolicy], Decision | None]


def _complete(outcome: CaseOutcome, reason: str) -> Decision:
    return Decision(Action.COMPLETE, reason, outcome)


def _r01_one_pes(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    a, b = s.minima
    if a is None or b is None or a.tier != "dft" or b.tier != "dft" or a.level_key != b.level_key:
        return Decision(Action.BLOCKED, "endpoints_not_on_one_pes", CaseOutcome.BLOCKED)
    return None


def _r02_same_basin(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    a, b = s.minima
    if a and b and a.basin_id == b.basin_id and not case.degenerate:
        return _complete(CaseOutcome.SAME_BASIN, "same_basin")
    return None


def _r03_window(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    a, b = s.minima
    if a and b and (b.energy_hartree - a.energy_hartree) * HARTREE_TO_KCAL_MOL > (
        p.gates.reaction_window_kcal
    ):
        return _complete(CaseOutcome.OUT_OF_WINDOW, "out_of_window")
    return None


def _r04_walltime(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "walltime") if s.expired else None


_CONNECTED = {
    "elementary": CaseOutcome.ELEMENTARY_STEP,
    "degenerate": CaseOutcome.DEGENERATE,
    "reassigned": CaseOutcome.REASSIGNED,
}


def _r05_connection(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.connection is None:
        return None
    if s.connection in _CONNECTED:
        return _complete(_CONNECTED[s.connection], f"connection:{s.connection}")
    if s.connection_attempts < QRC_AMPLITUDES:
        return Decision(Action.CONNECT, "connection_retry")  # amplitude x qrc_retry_factor
    return _complete(CaseOutcome.UNRESOLVED, "connection_failed")


def _r06_claim(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    return Decision(Action.CONNECT, "ts_validated") if s.claim is not None else None


def _r07_saddle(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.last_saddle == "converged" and s.ts_check is None:
        return Decision(Action.VALIDATE_TS, "saddle_converged")
    return None


def _r08_collapsed(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.ts_check == "collapsed" and s.intermediate is None:
        return Decision(Action.VALIDATE_INTERMEDIATE, "saddle_collapsed")
    return None


def _r09_distinct(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.intermediate == "distinct":  # the driver splits the case (classification.split)
        return _complete(CaseOutcome.MULTI_STEP, "intermediate_distinct")
    return None


def _r10_higher_order(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.ts_check == "higher_order" and s.saddle_attempts < p.max_saddle_attempts:
        return Decision(Action.REFINE_SADDLE, "higher_order_retry")
    return None


def _r11_screen(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    return Decision(Action.SCREEN, "screen") if p.screen and s.screen is None else None


def _r12_barrierless(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    # A continuous DFT path between the two DFT minima bounds the saddle from above.
    if s.screen is not None and s.screen.verdict == "barrierless":
        return _complete(CaseOutcome.BARRIERLESS, f"{s.screen.source}:barrierless")
    return None


def _r13_intermediate(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.screen is not None and s.screen.verdict == "intermediate" and s.intermediate is None:
        return Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")  # the lowest well
    return None


def _r14_seed(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    if s.seeds and s.saddle_attempts < p.max_saddle_attempts:
        return Decision(Action.REFINE_SADDLE, f"seed:{s.seeds[0].source}")
    return None


def _r15_no_path(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    return Decision(Action.FIND_PATH, "no_dft_path") if not s.path_runs else None


def _r16_exhausted(case: ReactionRecord, s: CaseState, p: CasePolicy) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "attempts_exhausted")


ROWS: tuple[Row, ...] = (
    _r01_one_pes, _r02_same_basin, _r03_window, _r04_walltime, _r05_connection, _r06_claim,
    _r07_saddle, _r08_collapsed, _r09_distinct, _r10_higher_order, _r11_screen,
    _r12_barrierless, _r13_intermediate, _r14_seed, _r15_no_path, _r16_exhausted,
)


def decide(case: ReactionRecord, state: CaseState, policy: CasePolicy) -> Decision:
    """First matching row of the table (pure, no IO); row 16 always matches."""
    return next(d for row in ROWS if (d := row(case, state, policy)) is not None)
