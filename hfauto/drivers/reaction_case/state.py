"""Reaction-case state and the pure decision table (design §7.3).

``decide`` evaluates the 16 rows of ``ROWS`` from the top and returns the first decision.  The
driver accumulates ``CaseState`` in memory; actions only change the state, never the table.
Rows 4-6 complete a case from evidence already in hand and so come before the walltime (row 7);
every row after it starts a computation or gives up. A saddle without an imaginary mode is a
failed attempt like any rejected saddle: an intermediate comes only from a profile's well
(row 12) or a QRC side (GEN-05). An association (``ReactionRecord.monomers``) is asked from its
separated monomers: SCREEN runs its relaxed scan, with or without a low-level engine (row 11),
and no string runs (rows 14-15: a string joins two minima, and the monomers' end is none).
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

QRC_AMPLITUDES = 2  # QRC tries at most two amplitudes (design §7.3, CONNECT)
STRING_CHUNKS = 3  # DFT string chunks (maxiter 20 each) per case


class Action(StrEnum):
    SCREEN = "screen"  # low-level path and DFT single points, or an association's scan
    REFINE_SADDLE = "refine_saddle"  # initial Hessian (xTB first, then DFT) -> saddle
    VALIDATE_TS = "validate_ts"  # separate DFT freq -> is_first_order_saddle
    FIND_PATH = "find_path"  # DFT string, only without a usable seed or after saddle failure
    CONNECT = "connect"  # QRC: displace, optimize, assign, connection gate
    VALIDATE_INTERMEDIATE = "validate_intermediate"  # relax_to_minimum -> Registry
    COMPLETE = "complete"
    BLOCKED = "blocked"


class ReactionPathsPolicy(BaseModel):
    """The budget (pipeline ``policy:``). max_saddle_attempts counts per case (a split child
    starts from 0; a stalled search's restart is not counted); a hypothesis and its split
    children share only the walltime_h deadline."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    walltime_h: float = 6.0
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
    carries its own imaginary ``mode``; ``hessian`` is a verified TS freq within HESSIAN_NEAR_A."""

    geometry: Geometry
    source: Literal["discovery_ts", "screen_ts", "screen_hei", "path_hei", "higher_order_retry",
                    "saddle_restart"]
    mode: tuple[float, ...] | None = None
    hessian: Evidence | None = None


@dataclass(frozen=True)
class CaseState:
    """Snapshot accumulated by the driver (never replayed from the log).

    ``minima`` (not in the §7.3 listing) holds the registry records of ``case.minima``, None when
    an endpoint has none: rows 1-3 need their tier, level, basin and energy, which a
    ReactionRecord does not carry. ``screen`` is the verdict of the latest DFT profile (SCREEN
    or string), which decides rows 6, 12 and 15 and goes into the record. ``neb_done``: the
    low-level path of SCREEN ran (after a shortcut, row 11 runs it once its seed has failed), or
    an association's scan.
    """

    minima: tuple[MinimumRecord | None, MinimumRecord | None] = (None, None)
    expired: bool = False
    screen: BarrierVerdict | None = None
    neb_done: bool = False
    seeds: tuple[Seed, ...] = ()  # unused seeds, consumed from the front
    saddle_attempts: int = 0
    last_saddle: Literal["converged", "failed"] | None = None
    ts_check: Literal["ok"] | None = None
    claim: SaddleClaim | None = None
    # a rejected QRC: sides in one basin (retried wider) or in one case key but two basins
    connection: ConnectionLabel | Literal["same_basin", "same_state"] | None = None
    connection_attempts: int = 0
    path_runs: int = 0  # DFT strings run, failed ones included
    intermediate: Literal["distinct", "same_as_endpoint", "relax_failed"] | None = None


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


Row = Callable[[ReactionRecord, CaseState, CaseRules], Decision | None]


def _complete(outcome: CaseOutcome, reason: str) -> Decision:
    return Decision(Action.COMPLETE, reason, outcome)


def _r01_one_pes(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    a, b = s.minima
    if a is None or b is None or a.tier != "dft" or b.tier != "dft" or a.level_key != b.level_key:
        return Decision(Action.BLOCKED, "endpoints_not_on_one_pes", CaseOutcome.BLOCKED)
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


def _r07_walltime(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "walltime") if s.expired else None


def _r08_connection(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.connection is None:
        return None
    if s.connection == "same_basin" and s.connection_attempts < QRC_AMPLITUDES:
        return Decision(Action.CONNECT, "connection_retry")  # a wider displacement
    if s.connection in ("same_basin", "same_state") and _attempts_left(s, p):
        return None  # a saddle of another process (sides in one basin or key): search on (11-16)
    return _complete(CaseOutcome.UNRESOLVED, "connection_failed")


def _r09_claim(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    fresh = s.claim is not None and s.connection is None
    return Decision(Action.CONNECT, "ts_validated") if fresh else None


def _r10_saddle(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.last_saddle == "converged" and s.ts_check is None:
        return Decision(Action.VALIDATE_TS, "saddle_converged")
    return None


def _r11_screen(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # the cheap low-level path comes before any string, also once a shortcut's seed has failed;
    # an association's scan needs no low-level engine
    fire = ((p.screen or bool(case.monomers)) and not s.neb_done
            and (s.screen is None or not s.seeds))
    return Decision(Action.SCREEN, "screen") if fire else None


def _r12_intermediate(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    if s.screen is not None and s.screen.verdict == "intermediate" and s.intermediate is None:
        return Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")  # the lowest well
    return None


def _attempts_left(s: CaseState, p: CaseRules) -> bool:
    return s.saddle_attempts < p.budget.max_saddle_attempts


def _r13_seed(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # a stalled search's restart (never itself restarted) is not counted
    if s.seeds and (s.seeds[0].source == "saddle_restart" or _attempts_left(s, p)):
        return Decision(Action.REFINE_SADDLE, f"seed:{s.seeds[0].source}")
    return None


# A string runs only while its peak seed can still be refined, never for an association.
def _r14_no_path(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    fire = not s.path_runs and _attempts_left(s, p) and not case.monomers
    return Decision(Action.FIND_PATH, "no_dft_path") if fire else None


def _r15_next_chunk(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    # The seeds of the latest string failed: its next chunk starts where it stopped.
    v = s.screen
    if (v is not None and v.source == "string" and v.verdict in ("single", "intermediate")
            and not s.seeds and s.path_runs < STRING_CHUNKS and _attempts_left(s, p)):
        return Decision(Action.FIND_PATH, "next_chunk")
    return None


def _r16_exhausted(case: ReactionRecord, s: CaseState, p: CaseRules) -> Decision | None:
    return _complete(CaseOutcome.UNRESOLVED, "attempts_exhausted")


ROWS: tuple[Row, ...] = (
    _r01_one_pes, _r02_same_basin, _r03_window, _r04_connected, _r05_distinct,
    _r06_barrierless, _r07_walltime, _r08_connection, _r09_claim, _r10_saddle, _r11_screen,
    _r12_intermediate, _r13_seed, _r14_no_path, _r15_next_chunk, _r16_exhausted,
)


def decide(case: ReactionRecord, state: CaseState, rules: CaseRules) -> Decision:
    """First matching row of the table (pure, no IO); row 16 always matches."""
    return next(d for row in ROWS if (d := row(case, state, rules)) is not None)
