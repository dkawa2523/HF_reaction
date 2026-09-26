"""decide(): the 17 rows of design §7.3 (ports tests/test_reaction_classification.py cases)."""

from __future__ import annotations

from dataclasses import replace
from typing import Literal

import pytest

from hfauto.core import records as r
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import FileRef, Geometry
from hfauto.drivers.reaction_case import state as st
from hfauto.drivers.reaction_case.state import Action, CaseState, Decision, decide, record_path

A, C, BV = Action, r.CaseOutcome, r.BarrierVerdict
TERM = r.StoichTerm(composition_id="CHN_q0_m1", coefficient=1)
CASE = r.ReactionRecord(reaction_id="rxn", reactants=(TERM,), products=(TERM,), minima=("ma", "mb"),
                        endpoints=("sa", "sb"), source="declared")
DEGENERATE = CASE.model_copy(update={"minima": ("ma", "ma"), "degenerate": True})
GEO = Geometry(file=FileRef(path="seed.xyz", sha256="0" * 64), fingerprint="f", symbols=("H",))
SEED = st.Seed(geometry=GEO, source="screen_ts", tangent=None)
CLAIM = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-1131.0, energy_hartree=-93.0)
POLICY = st.CasePolicy()


def minimum(mid: str, basin: str, kcal: float = 0.0, *,
            tier: Literal["screen", "dft"] = "dft", level: str = "L") -> r.MinimumRecord:
    return r.MinimumRecord(minimum_id=mid, basin_id=basin, composition_id="CHN_q0_m1",
                           species_id="s" + mid, tier=tier, level_key=level, opt_calc="o",
                           freq_calc="f", energy_hartree=-93.0 + kcal / HARTREE_TO_KCAL_MOL,
                           state_label="x")


MA, MB = minimum("ma", "A"), minimum("mb", "B", 10.0)
BASE = CaseState(minima=(MA, MB))  # nothing done yet
S = replace(BASE, screen=BV(verdict="proceed"))  # screened, barrier found
BLOCKED = Decision(A.BLOCKED, "endpoints_not_on_one_pes", C.BLOCKED)
NO_PATH = Decision(A.FIND_PATH, "no_dft_path")


ROW_CASES = [  # (row, id, case, state, policy, expected decision)
    (1, "missing", CASE, replace(BASE, minima=(MA, None)), POLICY, BLOCKED),
    (1, "screen_tier", CASE, replace(BASE, minima=(MA, minimum("mb", "B", tier="screen"))),
     POLICY, BLOCKED),
    (1, "two_levels", CASE, replace(S, minima=(MA, minimum("mb", "B", level="L2")),
                                    path_runs=("monotonic", "monotonic")), POLICY, BLOCKED),
    (2, "same_basin", CASE, replace(BASE, minima=(MA, minimum("mb", "A")), claim=CLAIM), POLICY,
     Decision(A.COMPLETE, "same_basin", C.SAME_BASIN)),
    (3, "uphill", CASE, replace(BASE, minima=(MA, minimum("mb", "B", 41.0))), POLICY,
     Decision(A.COMPLETE, "out_of_window", C.OUT_OF_WINDOW)),
    (4, "walltime", CASE, replace(S, expired=True, claim=CLAIM), POLICY,
     Decision(A.COMPLETE, "walltime", C.UNRESOLVED)),
    (5, "elementary", CASE, replace(S, claim=CLAIM, connection="elementary"), POLICY,
     Decision(A.COMPLETE, "connection:elementary", C.ELEMENTARY_STEP)),
    (5, "degenerate", DEGENERATE, replace(S, minima=(MA, MA), claim=CLAIM,
                                          connection="degenerate"), POLICY,
     Decision(A.COMPLETE, "connection:degenerate", C.DEGENERATE)),
    (5, "reassigned", CASE, replace(S, claim=CLAIM, connection="reassigned"), POLICY,
     Decision(A.COMPLETE, "connection:reassigned", C.REASSIGNED)),
    (5, "retry", CASE, replace(S, claim=CLAIM, connection="failed", connection_attempts=1), POLICY,
     Decision(A.CONNECT, "connection_retry")),
    (5, "failed", CASE, replace(S, claim=CLAIM, connection="failed", connection_attempts=2),
     POLICY, Decision(A.COMPLETE, "connection_failed", C.UNRESOLVED)),
    (6, "claim", CASE, replace(S, claim=CLAIM, last_saddle="converged", ts_check="ok"), POLICY,
     Decision(A.CONNECT, "ts_validated")),
    (7, "converged", CASE, replace(S, last_saddle="converged", saddle_attempts=1), POLICY,
     Decision(A.VALIDATE_TS, "saddle_converged")),
    (8, "collapsed", CASE, replace(S, last_saddle="converged", ts_check="collapsed",
                                   seeds=(SEED,)), POLICY,
     Decision(A.VALIDATE_INTERMEDIATE, "saddle_collapsed")),
    (9, "distinct", CASE, replace(S, ts_check="collapsed", intermediate="distinct"), POLICY,
     Decision(A.COMPLETE, "intermediate_distinct", C.MULTI_STEP)),
    (10, "higher_order", CASE, replace(S, last_saddle="converged", ts_check="higher_order",
                                       saddle_attempts=1), POLICY,
     Decision(A.REFINE_SADDLE, "higher_order_retry")),
    (11, "screen", CASE, BASE, POLICY, Decision(A.SCREEN, "screen")),
    (11, "degenerate_goes_on", DEGENERATE, replace(BASE, minima=(MA, MA)), POLICY,
     Decision(A.SCREEN, "screen")),
    (12, "barrierless", CASE, replace(S, screen=BV(verdict="barrierless"), seeds=(SEED,)), POLICY,
     Decision(A.COMPLETE, "screen:barrierless", C.BARRIERLESS)),
    (13, "multi_max", CASE, replace(S, path_runs=("multi_max",), seeds=(SEED,)), POLICY,
     Decision(A.VALIDATE_INTERMEDIATE, "path_multi_max")),
    (14, "monotonic", CASE, replace(S, path_runs=("monotonic",)), POLICY,
     Decision(A.COMPLETE, "dft_path_monotonic", C.BARRIERLESS)),
    (15, "seed", CASE, replace(S, seeds=(SEED,)), POLICY,
     Decision(A.REFINE_SADDLE, "seed:screen_ts")),
    (16, "no_path", CASE, replace(S, last_saddle="failed", saddle_attempts=1), POLICY, NO_PATH),
    (16, "screen_off", CASE, BASE, replace(POLICY, screen=False), NO_PATH),
    (16, "unavailable", CASE, replace(S, screen=BV(verdict="unavailable")), POLICY, NO_PATH),
    (17, "exhausted", CASE, replace(S, path_runs=("single_max",), seeds=(SEED,),
                                    saddle_attempts=2, last_saddle="failed"), POLICY,
     Decision(A.COMPLETE, "attempts_exhausted", C.UNRESOLVED)),
]


def test_table_has_17_rows_and_every_row_is_covered():
    assert len(st.ROWS) == 17
    assert {row for row, *_ in ROW_CASES} == set(range(1, 18))


@pytest.mark.parametrize(("row", "case", "state", "policy", "expected"),
                         [c[:1] + c[2:] for c in ROW_CASES],
                         ids=[f"row{c[0]}-{c[1]}" for c in ROW_CASES])
def test_row(row, case, state, policy, expected):
    assert decide(case, state, policy) == expected


def test_multi_max_validates_the_well_even_with_seeds_left():
    """eng P0-4: a multi_max path adds no HEI seed and never leads to REFINE_SADDLE."""
    state = record_path(replace(S, seeds=(SEED,)), "multi_max", replace(SEED, source="path_hei"))
    assert state.seeds == (SEED,) and state.path_runs == ("multi_max",)
    assert decide(CASE, state, POLICY).action is A.VALIDATE_INTERMEDIATE
    after = replace(state, intermediate="same_as_endpoint")
    assert decide(CASE, after, POLICY) == Decision(A.REFINE_SADDLE, "seed:screen_ts")


def test_hei_seed_only_from_a_single_max_path():
    hei = replace(SEED, source="path_hei")
    assert record_path(S, "single_max", hei).seeds == (hei,)
    for shape in ("multi_max", "monotonic", "failed"):
        assert record_path(S, shape, hei).seeds == ()  # type: ignore[arg-type]


def test_higher_order_retry_until_attempts_run_out():
    first = replace(S, last_saddle="converged", ts_check="higher_order", saddle_attempts=1,
                    path_runs=("single_max",))
    assert decide(CASE, first, POLICY) == Decision(A.REFINE_SADDLE, "higher_order_retry")
    second = replace(first, saddle_attempts=2)
    assert decide(CASE, second, POLICY) == Decision(A.COMPLETE, "attempts_exhausted",
                                                    C.UNRESOLVED)
