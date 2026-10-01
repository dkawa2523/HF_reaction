"""A reaction case's connection and intermediates (§3 X2): the atom labelling a split continues
(X2-1) and the granularity a case is judged at (X2-2). On a stepwise fake surface every named
structure is its own optimum at a chosen energy, so Registry.find alone decides the basins.
The O·H3 doublet has three chemical states here: HO + H2 (reactant), H2O + H (product; "side"
and "far" are two more basins of it, H3 beyond H1 and H3 far off) and the chain O-H1-H3
(bridged). HO + H + H (apart) is a fourth."""

from dataclasses import replace
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry.classification import finalize
from hfauto.chemistry.topology import bond_changes, bonds, state_label
from hfauto.chemistry.xyz import composition_key
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence
from hfauto.core.ids import species_artifact_id, species_id
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import (
    BarrierVerdict,
    CaseOutcome,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
    StoichTerm,
)
from hfauto.drivers.minimum import Registry, relax_to_minimum
from hfauto.drivers.reaction_case import connection, driver
from hfauto.drivers.reaction_case.actions import Profile
from hfauto.drivers.reaction_case.state import (
    Action,
    CaseRules,
    CaseState,
    Decision,
    decide,
    record_profile,
)

K = 1.0 / HARTREE_TO_KCAL_MOL
DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
SYMBOLS = ("O", "H", "H", "H")
COMPOSITION = composition_key(SYMBOLS, 0, 2)
WATER = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]
POINTS = {
    "reactant": np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [-0.24, 0.93, 0.0],
                          [3.74, 0.0, 0.0]]),  # O-H2 and H1-H3
    "product": np.array([*WATER, [0.0, 0.0, 3.0]]),  # H2O with H3 above O
    "side": np.array([*WATER, [3.0, 0.0, 0.0]]),  # the product's state, H3 beyond H1
    "far": np.array([*WATER, [-3.0, -3.0, 0.0]]),  # the product's state, H3 far off
    "bridged": np.array([*WATER, [1.76, 0.0, 0.0]]),  # H1-H3 0.80 A: one fragment O-H1-H3
    "apart": np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [-0.24, 0.93, 0.0],
                       [0.0, 0.0, 3.0]]),  # HO + H + H
}
POINTS["representative"] = POINTS["bridged"][[0, 3, 2, 1]]  # that basin as registered first


def surface(kcal: dict[str, float]) -> fakes.PES:
    """Every point takes the energy of the nearest named structure (kcal/mol, default 0)."""
    names = list(POINTS)

    def energy(x: np.ndarray) -> float:
        x = np.reshape(x, (-1, 3))
        return kcal.get(min(names, key=lambda n: np.linalg.norm(POINTS[n] - x)), 0.0) * K

    return fakes.PES(SYMBOLS, energy, POINTS)


def case_ctx(root: Path, ends=("reactant", "product"), known=(), kcal=None):
    """The case ends[0] -> ends[1], with the DFT basins of ``known`` in the Registry first."""
    pes = surface(kcal or {})
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, minima, species, ids = Registry([], load), {}, {}, {}
    for name in dict.fromkeys((*known, *ends)):
        geo = fakes.write_geometry(root, f"in/{name}.xyz", SYMBOLS, POINTS[name])
        species[name] = SpeciesRecord(
            species_id=name, composition_id=COMPOSITION, charge=0, multiplicity=2, geometry=geo,
            source="input", state_label=state_label(SYMBOLS, POINTS[name]))
        out = relax_to_minimum(pes.molecule(name, multiplicity=2), DFT, qm, load_xyz=load)
        record = registry.add(out, species[name], tier="dft")
        minima[record.minimum_id], ids[name] = (record, out.opt.final), record.minimum_id
    rt = driver.CaseRuntime(
        qm=qm, saddle=fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes),
        screen_qm=None, screen_path=None, method=DFT, screen_method=None, registry=registry,
        load_xyz=load, case_dir=root / "cases", file_ref=lambda p: fakes._ref(root, p),
        resolve=lambda ref: root / ref.path, map=lambda fn, items: [fn(x) for x in items],
        minima=minima, species=species, calcs={})
    terms = (StoichTerm(composition_id=COMPOSITION, coefficient=1),)
    case = ReactionRecord(reaction_id="rx", reactants=terms, products=terms,
                          minima=(ids[ends[0]], ids[ends[1]]), endpoints=ends,
                          source="discovery")
    ctx = driver.open_case(case, rt, CaseRules(screen=False), Deadline.after(600), root,
                           lambda _: None)
    return ctx, CaseState(minima=tuple(minima[m][0] for m in case.minima))


def qrc(ctx, state, monkeypatch, names, ts_kcal=20.0):
    """CONNECT from a validated TS ``ts_kcal`` above the reactant whose QRC sides optimize to
    the named structures (the QRC optimization itself: test_reaction_case_actions)."""
    sides = tuple(ctx.rt.qm.optimize(ctx.mol(POINTS[n]), DFT) for n in names)
    ctx.work.ts_freq = Evidence.model_validate({
        **sides[0].model_dump(), "task": "freq", "energy_hartree": ts_kcal * K, "n_external": 6,
        "frequencies_cm1": (-1500.0, 100.0, 200.0, 300.0, 400.0, 500.0),
        "imaginary_modes": (tuple(np.eye(12)[3]),)})
    monkeypatch.setattr(connection, "_sides", lambda *_: (
        sides, tuple(ctx.coords(s.final) for s in sides)))
    claim = SaddleClaim(saddle_calc="ts", freq_calc="freq", imag_cm1=-1500.0,
                        energy_hartree=ts_kcal * K)
    return connection.connect(ctx, replace(state, claim=claim), Decision(Action.CONNECT, "ts"))


def minimum_of(ctx, name: str) -> str:
    """The minimum id of the basin the named structure lies in."""
    opt = ctx.rt.qm.optimize(ctx.mol(POINTS[name]), DFT)
    return ctx.record(ctx.rt.registry.find(opt)).minimum_id


# X2-2: the granularity -------------------------------------------------------------------

def test_a_bond_changing_case_is_judged_by_chemical_state(tmp_path, monkeypatch):
    """HO + H2 -> H2O + H: a TS from the reactant into another basin of the product's state
    (the exit-channel complex, 5 kcal/mol above the product) is the elementary step. The claim
    keeps the minima actually reached; the record keeps the hypothesis' ends."""
    ctx, state = case_ctx(tmp_path, kcal={"product": -10.0, "side": -5.0})
    state = qrc(ctx, state, monkeypatch, ("reactant", "side"))
    side = minimum_of(ctx, "side")
    assert state.connection == "elementary" and side not in ctx.case.minima
    assert ctx.work.connection.minima == (ctx.case.minima[0], side)
    assert ctx.work.intermediate is None and state.intermediate is None
    decision = decide(ctx.case, state, ctx.rules)
    assert decision == Decision(Action.COMPLETE, "connection:elementary",
                                CaseOutcome.ELEMENTARY_STEP)
    record = finalize(ctx.case, decision, barrier=None, claim=state.claim,
                      connection=ctx.work.connection)
    assert record.minima == ctx.case.minima and record.connection.minima[1] == side


def test_a_ts_whose_sides_join_one_state_is_rejected_and_the_search_goes_on(
        tmp_path, monkeypatch):
    """Both sides in the product's state, in two basins: a saddle of another process (a
    conformer change of H2O···H) at any energy. No wider displacement, no claim, no split; the
    search goes on while saddle attempts are left, else connection_failed."""
    ctx, state = case_ctx(tmp_path, kcal={"product": -10.0, "side": -2.0})
    jobs = len(ctx.rt.qm.calls)
    state = qrc(ctx, state, monkeypatch, ("side", "product"), ts_kcal=15.0)
    assert state.connection == "same_state" and ctx.work.connection is None
    assert ctx.work.intermediate is None and ctx.work.split_ts is None
    assert ctx.rt.qm.calls[jobs:] == ["optimize", "optimize", "frequencies"]  # its new basin
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "no_dft_path")
    spent = replace(state, saddle_attempts=ctx.rules.budget.max_saddle_attempts)
    assert decide(ctx.case, spent, ctx.rules) == Decision(Action.COMPLETE, "connection_failed",
                                                          CaseOutcome.UNRESOLVED)


def test_sides_in_one_basin_are_retried_wider(tmp_path, monkeypatch):
    ctx, state = case_ctx(tmp_path)
    state = qrc(ctx, state, monkeypatch, ("product", "product"))
    assert state.connection == "same_basin"
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.CONNECT, "connection_retry")


def test_a_new_state_on_one_side_splits_and_none_at_an_end_is_reassigned(tmp_path,
                                                                        monkeypatch):
    ctx, state = case_ctx(tmp_path)
    split = qrc(ctx, state, monkeypatch, ("bridged", "reactant"))
    well, _ = ctx.work.intermediate
    assert split.intermediate == "distinct" and ctx.work.split_ts == (1, "ts")
    assert well.minimum_id == minimum_of(ctx, "bridged")
    ctx, state = case_ctx(tmp_path / "b")
    state = qrc(ctx, state, monkeypatch, ("bridged", "apart"))
    assert state.connection == "reassigned" and ctx.work.intermediate is None
    decision = decide(ctx.case, state, ctx.rules)
    record = finalize(ctx.case, decision, barrier=None, claim=state.claim,
                      connection=ctx.work.connection)
    assert record.minima == (minimum_of(ctx, "bridged"), minimum_of(ctx, "apart"))


@pytest.mark.parametrize("side_kcal,ts_kcal,expected", [
    (0.5, 0.8, "same_state"),  # G5-P2: the well and the TS lie within the resolution
    (0.5, 2.0, "distinct"),  # a hill of 1.5 kcal/mol between the well and the end
    (1.5, 2.0, "distinct"),  # 1.5 kcal/mol from the end
])
def test_a_same_state_case_is_judged_by_basin(tmp_path, monkeypatch, side_kcal, ts_kcal,
                                              expected):
    """A same-state case (H3 moved from above O to far off, as a declared torsion is): a new
    basin of that state is a distinct intermediate, unless it is an end at the resolution."""
    ctx, state = case_ctx(tmp_path, ("product", "far"), kcal={"side": side_kcal})
    assert qrc(ctx, state, monkeypatch, ("product", "far"), ts_kcal).connection == "elementary"
    ctx, state = case_ctx(tmp_path / "b", ("product", "far"), kcal={"side": side_kcal})
    state = qrc(ctx, state, monkeypatch, ("product", "side"), ts_kcal)
    assert (state.connection or state.intermediate) == expected
    if expected == "distinct":
        assert ctx.work.intermediate[0].minimum_id == minimum_of(ctx, "side")
    else:
        assert ctx.work.intermediate is None and ctx.work.connection is None


# G5-P2 on a profile: the lowest well of the latest DFT profile ---------------------------

@pytest.mark.parametrize("ends,well,kcal,expected", [
    (("product", "far"), "side", (0, 0.8, -0.5, 3, 0.2), "same_as_endpoint"),  # as end 0
    (("product", "far"), "side", (0, 3, -0.5, 0.8, 0.2), "same_as_endpoint"),  # as end 1
    (("product", "far"), "side", (0, 1.5, -0.5, 3, 0.2), "distinct"),  # a hill to each end
    (("product", "far"), "side", (0, 0.8, -1.2, 3, 0.2), "distinct"),  # too deep
    (("reactant", "product"), "side", (0, 1.5, -3, 3, 0.2), "same_as_endpoint"),  # its state
    (("reactant", "product"), "bridged", (0, 0.8, -0.5, 3, 0.2), "distinct"),  # a new state
])
def test_a_profile_well_is_an_end_by_state_or_at_the_resolution(tmp_path, ends, well, kcal,
                                                                 expected):
    """A bond-changing case: a well in an end's state is that end, a new state an
    intermediate, whatever the energies. A same-state case: a new basin of the state is that
    end only within a resolution of it with no hill of a resolution between them (G5-P2)."""
    names = (ends[0], "apart", well, "apart", ends[1])  # "apart" frames are never relaxed
    ctx, state = case_ctx(tmp_path, ends, kcal={ends[0]: kcal[0], well: kcal[2], ends[1]: kcal[4]})
    frames = [POINTS[n] for n in names]
    ctx.work.path = Profile(frames, tuple(e * K for e in kcal))
    state = record_profile(state, BarrierVerdict(verdict="intermediate", source="string"))
    state = connection.validate_intermediate(
        ctx, state, Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate"))
    assert state.intermediate == expected
    assert (ctx.work.intermediate is not None) is (expected == "distinct")
    if expected == "same_as_endpoint":  # the profile goes on from its highest peak
        assert [s.source for s in state.seeds] == ["path_hei"]


# X2-1: the atom labelling of an intermediate -----------------------------------------------

def by_qrc(ctx, state, monkeypatch):
    """A validated TS whose QRC sides reach the reactant and the H3O chain (two steps)."""
    return qrc(ctx, state, monkeypatch, ("reactant", "bridged"))


def by_profile(ctx, state, monkeypatch):
    """A string whose lowest well is the H3O chain."""
    names = ("reactant", "apart", "bridged", "apart", "product")
    ctx.work.path = Profile([POINTS[n] for n in names], tuple(e * K for e in (0, 3, -1, 3, 0)))
    state = record_profile(state, BarrierVerdict(verdict="intermediate", source="string"))
    decision = Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    return connection.validate_intermediate(ctx, state, decision)


@pytest.mark.parametrize("reach", [by_qrc, by_profile])
def test_an_intermediate_in_a_known_basin_continues_the_case_labelling(tmp_path, monkeypatch,
                                                                       reach):
    """G8-P2: the well joins the known basin as the side's own structure, a member without a
    job; the split child I -> P only breaks H1-H3, where the representative's labels would
    also move H3 onto O."""
    ctx, state = case_ctx(tmp_path, known=("representative",))
    basin, geometry = next(v for v in ctx.rt.minima.values() if v[0].species_id ==
                           "representative")
    jobs = len(ctx.rt.qm.calls)
    state = reach(ctx, state, monkeypatch)
    well, member = ctx.work.intermediate
    assert state.intermediate == "distinct" and well.basin_id == basin.basin_id
    assert "frequencies" not in ctx.rt.qm.calls[jobs:] and list(ctx.work.species) == [
        member.species_id]  # no species for the side found in the reactant's basin
    # (b) the member is the structure reached, as labelled; the representative stays
    assert np.allclose(ctx.coords(member.geometry), POINTS["bridged"])
    assert (member.source, member.state_label) == ("intermediate",
                                                   state_label(SYMBOLS, POINTS["bridged"]))
    assert well.species_id == "representative" != member.species_id
    # (c) the member enters the emitted MinimumRecord and the runtime
    record = ctx.case.model_copy(update={"outcome": CaseOutcome.MULTI_STEP})
    emitted = {a.artifact_id: a.payload for a in driver._artifacts(record, ctx.work)}
    assert emitted[well.minimum_id].members == ("representative", member.species_id)
    assert emitted[species_artifact_id(member.species_id)] == member
    assert ctx.rt.minima[well.minimum_id] == (emitted[well.minimum_id], geometry)
    # (a) the children continue the case labelling
    first, second = driver._children(record, ctx.rt, ctx)
    assert first.endpoints == ("reactant", member.species_id) and not first.torsional
    assert second.endpoints == (member.species_id, "product") and not second.torsional
    middle = driver._endpoint(ctx.rt, well.minimum_id, member.species_id)
    assert bonds(SYMBOLS, middle) == bonds(SYMBOLS, POINTS["bridged"])
    assert bond_changes(SYMBOLS, middle, POINTS["product"]) == (frozenset(), {(1, 3)})
    formed, broken = bond_changes(SYMBOLS, POINTS["representative"], POINTS["product"])
    assert len(formed | broken) == 3


@pytest.mark.parametrize("reach,name", [(by_qrc, "qrc1_1"), (by_profile, "int0_1")])
def test_an_intermediate_in_a_new_basin_is_the_cases_own_species(tmp_path, monkeypatch, reach,
                                                                 name):
    ctx, state = case_ctx(tmp_path)
    state = reach(ctx, state, monkeypatch)
    well, own = ctx.work.intermediate
    assert state.intermediate == "distinct"
    assert well.species_id == own.species_id == species_id("rx", name)
    assert well.members == (own.species_id,) and list(ctx.work.species) == [own.species_id]
    assert np.allclose(driver._endpoint(ctx.rt, well.minimum_id, own.species_id),
                       POINTS["bridged"])
