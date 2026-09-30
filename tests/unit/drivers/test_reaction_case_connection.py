"""Atom-index continuity of a case's intermediate (§3 X2-1). On a flat fake surface every
structure is its own optimum, so Registry.find alone decides the basins. A TS of HO + H2 ->
H2O + H reaches the exit-channel complex O-H1···H3, whose basin another case registered first
with H1 and H3 exchanged."""

from dataclasses import replace
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry.topology import bond_changes, bonds, state_label
from hfauto.chemistry.xyz import composition_key
from hfauto.core.ids import species_artifact_id, species_id
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import (
    CaseOutcome,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
    StoichTerm,
)
from hfauto.drivers.minimum import Registry, relax_to_minimum
from hfauto.drivers.reaction_case import connection, driver
from hfauto.drivers.reaction_case.state import Action, CaseRules, CaseState, Decision

DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
SYMBOLS = ("O", "H", "H", "H")
COMPOSITION = composition_key(SYMBOLS, 0, 2)
WATER = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]
POINTS = {
    "reactant": np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [-0.24, 0.93, 0.0],
                          [3.74, 0.0, 0.0]]),  # O-H2 and H1-H3
    "product": np.array([*WATER, [0.0, 0.0, 3.0]]),  # H2O with H3 above O
    "side": np.array([*WATER, [3.0, 0.0, 0.0]]),  # the same bonds, H3 beyond H1: another basin
}
POINTS["representative"] = POINTS["side"][[0, 3, 2, 1]]  # that basin as registered first


def case_ctx(root: Path, known: tuple[str, ...]):
    """The case reactant -> product, with the DFT basins of ``known`` in the Registry."""
    pes = fakes.PES(SYMBOLS, lambda x: 0.0, POINTS)
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, minima, species = Registry([], load), {}, {}
    for name in known:
        geo = fakes.write_geometry(root, f"in/{name}.xyz", SYMBOLS, POINTS[name])
        species[name] = SpeciesRecord(
            species_id=name, composition_id=COMPOSITION, charge=0, multiplicity=2, geometry=geo,
            source="input", state_label=state_label(SYMBOLS, POINTS[name]))
        out = relax_to_minimum(pes.molecule(name, multiplicity=2), DFT, qm, load_xyz=load)
        record = registry.add(out, species[name], tier="dft")
        minima[record.minimum_id] = (record, out.opt.final)
    ids = {record.species_id: m for m, (record, _) in minima.items()}
    rt = driver.CaseRuntime(
        qm=qm, saddle=fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes),
        screen_qm=None, screen_path=None, method=DFT, screen_method=None, registry=registry,
        load_xyz=load, case_dir=root / "cases", file_ref=lambda p: fakes._ref(root, p),
        resolve=lambda ref: root / ref.path, map=lambda fn, items: [fn(x) for x in items],
        minima=minima, species=species, calcs={})
    terms = (StoichTerm(composition_id=COMPOSITION, coefficient=1),)
    case = ReactionRecord(reaction_id="rx", reactants=terms, products=terms,
                          minima=(ids["reactant"], ids["product"]),
                          endpoints=("reactant", "product"), source="discovery")
    ctx = driver.open_case(case, rt, CaseRules(screen=False), Deadline.after(600), root,
                           lambda _: None)
    return ctx, CaseState(minima=tuple(minima[m][0] for m in case.minima))


def by_qrc(ctx, state):
    """A validated TS whose QRC sides reach the reactant and the complex (reassigned)."""
    claim = SaddleClaim(saddle_calc="ts", freq_calc="freq", imag_cm1=-1500.0, energy_hartree=0.0)
    sides = [ctx.rt.qm.optimize(ctx.mol(POINTS[n]), DFT) for n in ("reactant", "side")]
    reached = tuple(connection._assign(ctx, s, ctx.coords(s.final), f"qrc1_{i}")
                    for i, s in enumerate(sides))
    state = replace(state, claim=claim, connection_attempts=1)
    return connection._two_steps(ctx, state, "reassigned", reached)


def by_relaxation(ctx, state):
    """A saddle search that collapsed into the complex."""
    ctx.work.saddle = ctx.rt.qm.energy(ctx.mol(POINTS["side"]), DFT)
    decision = Decision(Action.VALIDATE_INTERMEDIATE, "saddle_collapsed")
    return connection.validate_intermediate(ctx, state, decision)


@pytest.mark.parametrize("reach", [by_qrc, by_relaxation])
def test_an_intermediate_in_a_known_basin_continues_the_case_labelling(tmp_path, reach):
    """G8-P2: the well joins the known basin as the side's own structure, a member without a
    job; the split child I -> P only moves H3 around the water (torsional), where the
    representative's labels would make a fake H scramble."""
    ctx, state = case_ctx(tmp_path, ("reactant", "product", "representative"))
    basin, geometry = next(v for v in ctx.rt.minima.values() if v[0].species_id ==
                           "representative")
    jobs = len(ctx.rt.qm.calls)
    state = reach(ctx, state)
    well, member = ctx.work.intermediate
    assert state.intermediate == "distinct" and well.basin_id == basin.basin_id
    assert "frequencies" not in ctx.rt.qm.calls[jobs:] and list(ctx.work.species) == [
        member.species_id]  # no species for the side found in the reactant's basin
    # (b) the member is the structure reached, as labelled; the representative stays
    assert np.allclose(ctx.coords(member.geometry), POINTS["side"])
    assert (member.source, member.state_label) == ("intermediate",
                                                   state_label(SYMBOLS, POINTS["side"]))
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
    assert second.endpoints == (member.species_id, "product") and second.torsional
    middle = driver._endpoint(ctx.rt, well.minimum_id, member.species_id)
    assert bonds(SYMBOLS, middle) == bonds(SYMBOLS, POINTS["product"])
    assert any(bond_changes(SYMBOLS, POINTS["representative"], POINTS["product"]))


@pytest.mark.parametrize("reach,name", [(by_qrc, "qrc1_1"), (by_relaxation, "int0_0")])
def test_an_intermediate_in_a_new_basin_is_the_cases_own_species(tmp_path, reach, name):
    ctx, state = case_ctx(tmp_path, ("reactant", "product"))
    state = reach(ctx, state)
    well, own = ctx.work.intermediate
    assert state.intermediate == "distinct"
    assert well.species_id == own.species_id == species_id("rx", name)
    assert well.members == (own.species_id,) and list(ctx.work.species) == [own.species_id]
    assert np.allclose(driver._endpoint(ctx.rt, well.minimum_id, own.species_id), POINTS["side"])
