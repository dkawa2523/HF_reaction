"""decide(): the 17 rows of design §7.3 (ports tests/test_reaction_classification.py cases)."""

from __future__ import annotations

from dataclasses import replace
from typing import Literal

import pytest

from hfauto.core import records as r
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import FileRef, Geometry
from hfauto.drivers.reaction_case import state as st
from hfauto.drivers.reaction_case.state import Action, CaseState, Decision, decide, record_profile

A, C, BV = Action, r.CaseOutcome, r.BarrierVerdict
TERM = r.StoichTerm(composition_id="CHN_q0_m1", coefficient=1)
CASE = r.ReactionRecord(reaction_id="rxn", reactants=(TERM,), products=(TERM,), minima=("ma", "mb"),
                        endpoints=("sa", "sb"), source="declared")
DEGENERATE = CASE.model_copy(update={"minima": ("ma", "ma"), "degenerate": True})
GEO = Geometry(file=FileRef(path="seed.xyz", sha256="0" * 64), fingerprint="f", symbols=("H",))
SEED = st.Seed(geometry=GEO, source="screen_ts", tangent=None)
HIGHER = st.Seed(geometry=GEO, source="higher_order_retry", tangent=None)
CLAIM = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-1131.0, energy_hartree=-93.0)
SOFT = CLAIM.model_copy(update={"imag_cm1": -30.0})  # |nu| < saddle_cm1 = 50
RULES = st.CaseRules()


def minimum(mid: str, basin: str, kcal: float = 0.0, *,
            tier: Literal["screen", "dft"] = "dft", level: str = "L") -> r.MinimumRecord:
    return r.MinimumRecord(minimum_id=mid, basin_id=basin, composition_id="CHN_q0_m1",
                           species_id="s" + mid, tier=tier, level_key=level, opt_calc="o",
                           freq_calc="f", energy_hartree=-93.0 + kcal / HARTREE_TO_KCAL_MOL,
                           state_label="x")


MA, MB = minimum("ma", "A"), minimum("mb", "B", 10.0)
BASE = CaseState(minima=(MA, MB))  # nothing done yet
S = replace(BASE, screen=BV(verdict="single", source="screen"))  # screened, one peak found
BLOCKED = Decision(A.BLOCKED, "endpoints_not_on_one_pes", C.BLOCKED)
NO_PATH = Decision(A.FIND_PATH, "no_dft_path")
NEXT_CHUNK = Decision(A.FIND_PATH, "next_chunk")
EXHAUSTED = Decision(A.COMPLETE, "attempts_exhausted", C.UNRESOLVED)
ELEMENTARY = Decision(A.COMPLETE, "connection:elementary", C.ELEMENTARY_STEP)
DISTINCT = Decision(A.COMPLETE, "intermediate_distinct", C.MULTI_STEP)
WALLTIME = Decision(A.COMPLETE, "walltime", C.UNRESOLVED)
COLLAPSED = Decision(A.VALIDATE_INTERMEDIATE, "saddle_collapsed")


ROW_CASES = [  # (row, id, case, state, rules, expected decision)
    (1, "missing", CASE, replace(BASE, minima=(MA, None)), RULES, BLOCKED),
    (1, "screen_tier", CASE, replace(BASE, minima=(MA, minimum("mb", "B", tier="screen"))),
     RULES, BLOCKED),
    (1, "two_levels", CASE, replace(S, minima=(MA, minimum("mb", "B", level="L2")),
                                    path_runs=2), RULES, BLOCKED),
    (2, "same_basin", CASE, replace(BASE, minima=(MA, minimum("mb", "A")), claim=CLAIM), RULES,
     Decision(A.COMPLETE, "same_basin", C.SAME_BASIN)),
    (3, "uphill", CASE, replace(BASE, minima=(MA, minimum("mb", "B", 41.0))), RULES,
     Decision(A.COMPLETE, "out_of_window", C.OUT_OF_WINDOW)),
    (4, "elementary", CASE, replace(S, claim=CLAIM, connection="elementary"), RULES, ELEMENTARY),
    (4, "at_the_deadline", CASE, replace(S, expired=True, claim=CLAIM, connection="elementary"),
     RULES, ELEMENTARY),  # U9-P5: evidence in hand completes the case
    (4, "degenerate", DEGENERATE, replace(S, minima=(MA, MA), claim=CLAIM,
                                          connection="degenerate"), RULES,
     Decision(A.COMPLETE, "connection:degenerate", C.DEGENERATE)),
    (4, "reassigned", CASE, replace(S, claim=CLAIM, connection="reassigned"), RULES,
     Decision(A.COMPLETE, "connection:reassigned", C.REASSIGNED)),
    (5, "distinct", CASE, replace(S, intermediate="distinct"), RULES, DISTINCT),
    (5, "at_the_deadline", CASE, replace(S, expired=True, intermediate="distinct"), RULES,
     DISTINCT),
    (6, "screen", CASE, replace(S, screen=BV(verdict="barrierless", source="screen"),
                                seeds=(SEED,)), RULES,
     Decision(A.COMPLETE, "screen:barrierless", C.BARRIERLESS)),
    (6, "string_at_the_deadline", CASE, replace(S, screen=BV(verdict="barrierless",
                                                             source="string"),
                                                path_runs=1, expired=True), RULES,
     Decision(A.COMPLETE, "string:barrierless", C.BARRIERLESS)),
    (7, "walltime", CASE, replace(S, expired=True, claim=CLAIM), RULES, WALLTIME),
    (7, "no_retry", CASE, replace(S, expired=True, claim=CLAIM, connection="same_basin",
                                  connection_attempts=1), RULES, WALLTIME),
    (8, "same_basin_retry", CASE, replace(S, claim=CLAIM, connection="same_basin",
                                          connection_attempts=1), RULES,
     Decision(A.CONNECT, "connection_retry")),
    (8, "failed", CASE, replace(S, claim=CLAIM, connection="failed", connection_attempts=1),
     RULES, Decision(A.COMPLETE, "connection_failed", C.UNRESOLVED)),
    (9, "claim", CASE, replace(S, claim=CLAIM, last_saddle="converged", ts_check="ok"), RULES,
     Decision(A.CONNECT, "ts_validated")),
    (10, "converged", CASE, replace(S, last_saddle="converged", saddle_attempts=1), RULES,
     Decision(A.VALIDATE_TS, "saddle_converged")),
    (11, "collapsed", CASE, replace(S, last_saddle="converged", ts_check="collapsed",
                                    seeds=(SEED,)), RULES, COLLAPSED),
    (11, "after_an_endpoint_well", CASE, replace(S, last_saddle="converged", ts_check="collapsed",
                                                 intermediate="same_as_endpoint"), RULES,
     COLLAPSED),
    (11, "soft_ts_qrc_failed", CASE, replace(S, claim=SOFT, last_saddle="converged",
                                             ts_check="ok", connection="failed",
                                             connection_attempts=1), RULES, COLLAPSED),
    (12, "screen", CASE, BASE, RULES, Decision(A.SCREEN, "screen")),
    (12, "degenerate_goes_on", DEGENERATE, replace(BASE, minima=(MA, MA)), RULES,
     Decision(A.SCREEN, "screen")),
    (13, "intermediate", CASE, replace(S, screen=BV(verdict="intermediate", source="screen"),
                                       seeds=(SEED,)), RULES,
     Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")),
    (14, "seed", CASE, replace(S, seeds=(SEED,)), RULES,
     Decision(A.REFINE_SADDLE, "seed:screen_ts")),
    (14, "higher_order_retry", CASE, replace(S, seeds=(HIGHER,), last_saddle="failed",
                                             saddle_attempts=1), RULES,
     Decision(A.REFINE_SADDLE, "seed:higher_order_retry")),
    (15, "no_path", CASE, replace(S, last_saddle="failed", saddle_attempts=1), RULES, NO_PATH),
    (15, "screen_off", CASE, BASE, replace(RULES, screen=False), NO_PATH),
    (15, "unavailable", CASE, replace(S, screen=BV(verdict="unavailable", source="screen")),
     RULES, NO_PATH),
    (16, "seeds_failed", CASE, replace(S, screen=BV(verdict="single", source="string"),
                                       path_runs=1, saddle_attempts=1, last_saddle="failed"),
     RULES, NEXT_CHUNK),
    (16, "after_an_endpoint_well", CASE, replace(S, screen=BV(verdict="intermediate",
                                                              source="string"),
                                                 intermediate="same_as_endpoint", path_runs=2,
                                                 saddle_attempts=2), RULES, NEXT_CHUNK),
    (17, "exhausted", CASE, replace(S, path_runs=1, seeds=(SEED,), saddle_attempts=2,
                                    last_saddle="failed"), RULES, EXHAUSTED),
    (17, "chunks_used", CASE, replace(S, screen=BV(verdict="single", source="string"),
                                      path_runs=3, saddle_attempts=1), RULES, EXHAUSTED),
    (17, "string_failed", CASE, replace(S, path_runs=1, saddle_attempts=1), RULES, EXHAUSTED),
]


def test_table_has_17_rows_and_every_row_is_covered():
    assert len(st.ROWS) == 17
    assert {row for row, *_ in ROW_CASES} == set(range(1, 18))


@pytest.mark.parametrize(("row", "case", "state", "rules", "expected"),
                         [c[:1] + c[2:] for c in ROW_CASES],
                         ids=[f"row{c[0]}-{c[1]}" for c in ROW_CASES])
def test_row(row, case, state, rules, expected):
    assert decide(case, state, rules) == expected


def test_a_well_is_validated_before_any_seed_then_the_seeds_go_on():
    """eng P0-4: an intermediate profile adds no peak seed and validates its well first."""
    well = BV(verdict="intermediate", source="string")
    state = record_profile(replace(S, seeds=(SEED,)), well, replace(SEED, source="path_hei"))
    assert state.seeds == (SEED,) and state.path_runs == 1 and state.screen == well
    assert decide(CASE, state, RULES).action is A.VALIDATE_INTERMEDIATE
    after = replace(state, intermediate="same_as_endpoint")
    assert decide(CASE, after, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")


def test_a_peak_seeds_only_a_single_step_profile():
    hei = replace(SEED, source="path_hei")
    assert record_profile(S, BV(verdict="single", source="string"), hei).seeds == (hei,)
    for verdict in ("barrierless", "intermediate", "unavailable"):
        assert record_profile(S, BV(verdict=verdict, source="string"), hei).seeds == ()


def test_a_new_profile_reopens_the_checks_of_an_earlier_saddle():
    """A string after a collapsed saddle: its own well is validated, not the old saddle."""
    old = replace(S, saddle_attempts=2, last_saddle="converged", ts_check="collapsed",
                  intermediate="same_as_endpoint")
    state = record_profile(old, BV(verdict="intermediate", source="string"))
    assert (state.last_saddle, state.ts_check, state.intermediate) == (None, None, None)
    assert decide(CASE, state, RULES) == Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert record_profile(BASE, BV(verdict="single", source="screen")).path_runs == 0


def test_seeds_share_one_attempt_budget():
    """U6-P3 / U6-P6: a higher-order push or a stalled search's restart is an ordinary seed."""
    first = replace(S, seeds=(HIGHER,), last_saddle="failed", saddle_attempts=1, path_runs=1)
    assert decide(CASE, first, RULES) == Decision(A.REFINE_SADDLE, "seed:higher_order_retry")
    assert decide(CASE, replace(first, saddle_attempts=2), RULES) == EXHAUSTED


def test_a_soft_ts_whose_qrc_failed_is_validated_as_a_collapsed_saddle():
    """U6-P5: a TS with |nu| < saddle_cm1 is a TS, but once its QRC fails for good (the one
    same-basin retry included) it joins row 11, whose action withdraws the claim."""
    failed = replace(S, claim=SOFT, last_saddle="converged", ts_check="ok", connection="same_basin",
                     connection_attempts=1)
    assert decide(CASE, failed, RULES) == Decision(A.CONNECT, "connection_retry")
    failed = replace(failed, connection_attempts=2)
    assert decide(CASE, failed, RULES) == COLLAPSED
    validated = replace(failed, claim=None, connection=None, ts_check=None, last_saddle=None)
    assert decide(CASE, replace(validated, intermediate="distinct"), RULES) == DISTINCT
    assert decide(CASE, replace(validated, intermediate="same_as_endpoint"), RULES) == NO_PATH
    unresolved = Decision(A.COMPLETE, "connection_failed", C.UNRESOLVED)
    assert decide(CASE, replace(failed, claim=CLAIM), RULES) == unresolved  # a hard TS
