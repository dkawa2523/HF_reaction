"""Split children of the reaction-paths stage (G8-P7, G8-P3): reuse of a conclusion only."""

from types import SimpleNamespace
from typing import cast

from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.method import Deadline
from hfauto.core.records import CaseOutcome as O
from hfauto.core.records import MinimumRecord, ReactionRecord, StoichTerm
from hfauto.drivers.reaction_case.driver import CaseResult, CaseRuntime
from hfauto.stages.reaction_paths import _Book

GEO = Geometry(file=FileRef(path="x.xyz", sha256="0" * 64), fingerprint="x", symbols=("H",))
TERM = StoichTerm(composition_id="H_q0_m2", coefficient=1)
DEADLINE = Deadline.after(3600.0)


def _minimum(mid: str, state: str) -> MinimumRecord:
    return MinimumRecord(minimum_id=mid, basin_id=mid, composition_id="H_q0_m2", species_id=mid,
                         tier="dft", level_key="L", opt_calc="o", freq_calc="f",
                         energy_hartree=-1.0, state_label=state)


def _case(rid: str, minima: tuple[str, str], outcome: O | None = None) -> ReactionRecord:
    return ReactionRecord(reaction_id=rid, reactants=(TERM,), products=(TERM,), minima=minima,
                          endpoints=minima, source="declared", outcome=outcome)


def _book(max_depth: int) -> _Book:
    states = {"a": "A", "a2": "A", "b": "B", "c": "C"}
    minima = {m: (_minimum(m, s), GEO) for m, s in states.items()}
    return _Book(cast(CaseRuntime, SimpleNamespace(minima=minima, calcs={})), max_depth)


def _parent_done(book: _Book) -> list:
    """p (A -> C) splits into A -> B, the key of the queued hypothesis t, and B -> C."""
    book.queue([_case("p", ("a", "c")), _case("t", ("a2", "b"))])
    children = (_case("p_split1", ("a", "b")), _case("p_split2", ("b", "c")))
    result = CaseResult(_case("p", ("a", "c"), O.MULTI_STEP), children, ())
    return book.emit(result, 0, DEADLINE)


def test_a_child_waits_for_a_queued_case_of_its_key_and_takes_its_conclusion() -> None:
    book = _book(max_depth=2)
    assert [c.reaction_id for c, _, _ in _parent_done(book)] == ["p_split2"]
    assert book.waited() == []  # t not driven yet: still waiting
    book.emit(CaseResult(_case("t", ("a2", "b"), O.BARRIERLESS), (), ()), 0, DEADLINE)
    assert book.waited() == []
    child = book.artifacts["p_split1"].payload
    assert isinstance(child, ReactionRecord)
    assert (child.outcome, child.reasons, child.log) == (O.BARRIERLESS, ("same_as:t",), None)


def test_a_failure_to_conclude_is_not_shared_the_child_is_driven() -> None:
    """W3_s6: 790468f505 ends attempts_exhausted, while the split child of the same key finds
    its elementary TS from its own ends; copying UNRESOLVED would lose that TS."""
    book = _book(max_depth=2)
    _parent_done(book)
    book.emit(CaseResult(_case("t", ("a2", "b"), O.UNRESOLVED), (), ()), 0, DEADLINE)
    assert [(c.reaction_id, depth) for c, depth, _ in book.waited()] == [("p_split1", 1)]
    assert "p_split1" not in book.artifacts


def test_a_child_past_the_split_depth_is_recorded_not_driven() -> None:
    book = _book(max_depth=0)
    assert _parent_done(book) == []
    capped = book.artifacts["p_split2"].payload
    assert isinstance(capped, ReactionRecord)
    assert (capped.outcome, capped.reasons) == (O.UNRESOLVED, ("split_depth",))
    book.emit(CaseResult(_case("t", ("a2", "b"), O.UNRESOLVED), (), ()), 0, DEADLINE)
    assert book.waited() == []  # nothing to reuse, and too deep to drive
    split1 = book.artifacts["p_split1"].payload
    assert isinstance(split1, ReactionRecord) and split1.reasons == ("split_depth",)
