"""reaction-paths stage on fake engines (§7.3, §8.2), including the CH-05 regression."""

import json
from collections import Counter
from functools import partial
from pathlib import Path

import fakes
import pytest

from hfauto.backends.protocols import Capability as Cap
from hfauto.chemistry.xyz import composition_key
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType as T
from hfauto.core.records import CaseOutcome as O
from hfauto.core.records import (
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    ReactionTrial,
    SpeciesRecord,
)
from hfauto.core.system import ReactionInput, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import Registry, calc_id, relax_to_minimum
from hfauto.drivers.reaction_case import driver
from hfauto.drivers.reaction_case.state import Action
from hfauto.stages.reaction_paths import ReactionPathsConfig, ReactionPathsStage

pytestmark = pytest.mark.integration
DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
XTB = MethodSpec(id="gfn2", kind="xtb", gfn=2)
ENDS = {"reactant": "reactant", "product": "product"}  # species id -> PES point
SYSTEM = SystemConfig(  # declared endpoints name an xyz, which this stage never reads
    system_id="t", reactions=[ReactionInput(id="rx", reactant="reactant", product="product")],
    species=[SpeciesInput(id=n, role="endpoint", xyz=Path(f"{n}.xyz"), multiplicity=1)
             for n in ENDS])
FAKES = {(Cap.QM, "nwchem"): fakes.FakeQM, (Cap.PATH, "nwchem_string"): fakes.FakePath,
         (Cap.QM, "xtb"): fakes.FakeQM, (Cap.PATH, "pysis_neb"): partial(fakes.FakePath, tsopt=True)}


class CollapsingSaddle(fakes.FakeSaddle):  # every saddle search falls into the intermediate
    def refine(self, seed, method, *, hessian, mode):
        key = self._key("collapse", seed.fingerprint())
        return self._evidence("saddle", seed, method, key, self._start(seed, key),
                              self.pes.points["intermediate"])


class KeyedQM(fakes.FakeQM):  # the calculation id of every result, as a JobStore sees it
    def __init__(self, root, pes):
        super().__init__(root, pes)
        self.jobs = Counter()

    def frequencies(self, mol, method, **kw):
        ev = super().frequencies(mol, method, **kw)
        self.jobs[calc_id(ev)] += 1
        return ev

    def optimize(self, mol, method, *, init_hessian=None):
        ev = super().optimize(mol, method, init_hessian=init_hessian)
        self.jobs[calc_id(ev)] += 1
        return ev


class FirstStepSaddle(fakes.FakeSaddle):
    """Every search that starts with H on the reactant's side of the triple well's intermediate
    (the R -> I step) stops without a last frame (no restart); the I -> P searches converge."""

    def refine(self, seed, method, *, hessian, mode):
        a, h, b = seed.xyz.coords
        if (h - 0.5 * (a + b)) @ (b - a) < 0:
            self.calls.append("refine:failed")
            return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="scripted")
        return super().refine(seed, method, hessian=hessian, mode=mode)


def dft_view(root, pes, points=ENDS):
    """The view after structures and minima(dft): species, calculations and minima."""
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, arts, basins = Registry([], load), [], {}
    for sid, point in points.items():
        geo = fakes.write_geometry(root, f"in/{sid}.xyz", pes.symbols, pes.points[point])
        s = SpeciesRecord(species_id=sid, state_label="x", charge=0, multiplicity=1,
                          geometry=geo, source="input",
                          composition_id=composition_key(pes.symbols, 0, 1))
        out = relax_to_minimum(pes.molecule(point), DFT, qm, load_xyz=load)
        record = registry.add(out, s, tier="dft")
        basins[record.basin_id] = record
        arts += [(sid, T.SPECIES, s), (calc_id(out.opt), T.CALCULATION, out.opt),
                 (calc_id(out.freq), T.CALCULATION, out.freq)]
    arts += [(m.minimum_id, T.MINIMUM, m) for m in basins.values()]
    source = next(m.minimum_id for m in basins.values() if m.species_id == "reactant")
    artifacts = [Artifact(artifact_id=i, type=t, payload=p) for i, t, p in arts]
    return Manifest(run_id="run", stage_id="dft", created_at="now", artifacts=artifacts), source


def run_stage(fake_runtime, root, pes, view, *, screen=True, saddle=None, qm=None,
              system=SYSTEM, **policy):
    engines = {key: cls(root, pes) for key, cls in FAKES.items()}
    engines[Cap.QM, "nwchem"] = qm or engines[Cap.QM, "nwchem"]
    engines[Cap.SADDLE, "nwchem_saddle"] = saddle or fakes.FakeSaddle(root, pes)
    rt = fake_runtime(system, engines, methods={"pbe0": DFT, "gfn2": XTB})
    config = ReactionPathsConfig(method="pbe0", policy=policy, engines={
        "qm": "nwchem", "saddle": "nwchem_saddle", "path": "nwchem_string"},
        screen={"method": "gfn2", "qm": "xtb", "path": "pysis_neb"} if screen else None)
    arts = ReactionPathsStage().run(view, config, rt)
    return {a.artifact_id: a.payload for a in arts if a.type == T.REACTION}, arts, engines


def test_declared_reaction_becomes_a_typed_elementary_step(tmp_run, fake_runtime) -> None:
    view, _ = dft_view(tmp_run, fakes.double_well())  # CH-05: typed records end to end
    reactions, arts, _ = run_stage(fake_runtime, tmp_run, fakes.double_well(), view)
    rx = reactions["rx"]
    assert isinstance(rx, ReactionRecord) and rx.outcome is O.ELEMENTARY_STEP
    assert (rx.barrier.verdict, rx.barrier.source) == ("single", "screen")
    assert rx.saddle.imag_cm1 < -50
    calcs = {a.artifact_id for a in arts if a.type == T.CALCULATION}
    assert {rx.saddle.saddle_calc, rx.saddle.freq_calc, *rx.connection.side_calcs} <= calcs
    assert (tmp_run / "stage" / rx.log).read_text().count("\n") >= 5


def test_a_case_that_raises_keeps_what_it_registered_and_the_next_case_runs(
        tmp_run, fake_runtime, monkeypatch) -> None:
    """G3-P4: an exception right after _register (rx's well on the triple-well screen path)
    leaves rx UNRESOLVED (error:<type>) and still emits the minimum, species and calculations it
    registered (rx2 joins that basin without a freq, so the species and freq are rx's own); the
    next case runs. HFAUTO_STRICT=1 re-raises."""
    validate = driver.HANDLERS[Action.VALIDATE_INTERMEDIATE]

    def validate_then_raise(ctx, state, decision):
        state = validate(ctx, state, decision)
        if ctx.case.reaction_id == "rx":
            raise KeyError("Te")
        return state

    monkeypatch.setitem(driver.HANDLERS, Action.VALIDATE_INTERMEDIATE, validate_then_raise)
    pes, again = fakes.triple_well(), ReactionInput(id="rx2", reactant="reactant",
                                                    product="product")
    system = SYSTEM.model_copy(update={"reactions": [*SYSTEM.reactions, again]})
    view = dft_view(tmp_run, pes)[0]
    with pytest.raises(KeyError):  # HFAUTO_STRICT=1 re-raises
        run_stage(fake_runtime, tmp_run, pes, view, system=system)
    monkeypatch.delenv("HFAUTO_STRICT")
    reactions, arts, _ = run_stage(fake_runtime, tmp_run, pes, view, system=system)
    rx = reactions["rx"]
    assert rx.outcome is O.UNRESOLVED and rx.reasons == ("error:KeyError",)
    log = (tmp_run / "stage" / rx.log).read_text().splitlines()
    assert json.loads(log[-1]) == {"action": "error", "reason": "error:KeyError", "detail": "'Te'"}
    out = {a.artifact_id: a.payload for a in arts}
    [well] = [m for m in out.values() if isinstance(m, MinimumRecord)]
    assert well.species_id.startswith("spc_rx_")  # registered by rx
    assert {species_artifact_id(well.species_id), well.opt_calc, well.freq_calc} <= set(out)
    assert reactions["rx2"].outcome is O.MULTI_STEP


@pytest.mark.parametrize("pes,points,outcome", [
    (fakes.symmetric_double_well, ENDS, O.DEGENERATE),
    (fakes.flat_uphill, ENDS, O.BARRIERLESS),  # the NEB between the DFT minima closes it
    (fakes.harmonic, {"reactant": "minimum", "product": "start"}, O.SAME_BASIN),
])
def test_outcomes_on_model_surfaces(tmp_run, fake_runtime, pes, points, outcome) -> None:
    view = dft_view(tmp_run, pes(), points)[0]
    reactions, _, engines = run_stage(fake_runtime, tmp_run, pes(), view)
    assert reactions["rx"].outcome is outcome
    assert outcome is not O.BARRIERLESS or reactions["rx"].reasons == ("screen:barrierless",)
    assert reactions["rx"].degenerate is (outcome is O.DEGENERATE)
    assert outcome is not O.SAME_BASIN or not any(e.calls for e in engines.values())  # no job


@pytest.mark.parametrize("screen", [True, False])  # the well on the NEB path or the string
def test_triple_well_splits_into_two_elementary_children(tmp_run, fake_runtime, screen) -> None:
    view, _ = dft_view(tmp_run, fakes.triple_well())
    reactions, arts, _ = run_stage(fake_runtime, tmp_run, fakes.triple_well(), view, screen=screen)
    assert reactions["rx"].outcome is O.MULTI_STEP
    assert reactions["rx"].barrier.verdict == "intermediate"
    children = [reactions[f"rx_split{i}"] for i in (1, 2)]
    assert [c.outcome for c in children] == [O.ELEMENTARY_STEP] * 2
    new_minima = {a.artifact_id for a in arts if a.type == T.MINIMUM}  # the intermediate
    assert children[0].minima[1] == children[1].minima[0] in new_minima


def test_a_split_child_queued_behind_a_sibling_that_used_its_budget_has_its_own(
        tmp_run, fake_runtime) -> None:
    """X7-1: the budget is counts only, per case. rx splits at the triple well's intermediate;
    rx_split1 (R -> I) runs first and uses up its saddle attempts (SCREEN's seed, then the
    string's); rx_split2 (I -> P), queued behind it, is still driven from 0 attempts and finds
    its TS."""
    pes = fakes.triple_well()
    view, _ = dft_view(tmp_run, pes)
    saddle = FirstStepSaddle(tmp_run, pes)
    reactions = run_stage(fake_runtime, tmp_run, pes, view, saddle=saddle)[0]
    first, second = (reactions[f"rx_split{i}"] for i in (1, 2))
    assert reactions["rx"].outcome is O.MULTI_STEP
    assert (first.outcome, first.reasons) == (O.UNRESOLVED, ("attempts_exhausted",))
    assert second.outcome is O.ELEMENTARY_STEP
    assert saddle.calls == ["refine:failed"] * 2 + ["refine"]
    log = [json.loads(line) for line in (tmp_run / "stage" / second.log).read_text().splitlines()]
    assert [e["action"] for e in log if "action" in e][:2] == ["screen", "refine_saddle"]


def test_a_saddle_search_that_falls_into_a_well_goes_on_to_the_screen_path(tmp_run,
                                                                          fake_runtime) -> None:
    """G3-P1: a saddle without an imaginary mode is a failed attempt like any rejected saddle;
    the intermediate comes from the screen profile's well (row 10), not from the saddle."""
    pes = fakes.triple_well()
    view, source = dft_view(tmp_run, pes)
    ts = fakes.write_geometry(tmp_run, "ts1.xyz", pes.symbols, pes.points["ts1"])
    found = DiscoveryRecord(discovery_id="d1", source_minimum=source, mechanism="nt2",
                            outcome="product", product_species="product", ts=ts)
    view.artifacts.append(Artifact(artifact_id="d1", type=T.DISCOVERY, payload=found))
    reactions, _, _ = run_stage(fake_runtime, tmp_run, pes, view,  # seeded by the shortcut
                                saddle=CollapsingSaddle(tmp_run, pes), max_split_depth=0)
    assert reactions["rx"].outcome is O.MULTI_STEP
    children = [reactions[f"rx_split{i}"] for i in (1, 2)]  # G8-P3: recorded, not driven
    assert [(c.outcome, c.reasons, c.log) for c in children] == [
        (O.UNRESOLVED, ("split_depth",), None)] * 2
    log = (tmp_run / "stage" / reactions["rx"].log).read_text()
    assert '"ts_rejected:no_imaginary_mode"' in log and '"reason": "path_intermediate"' in log
    assert '"saddle_collapsed"' not in log


def test_a_ts_joining_an_endpoint_to_a_new_basin_splits_and_its_child_replays_it(
        tmp_run, fake_runtime) -> None:
    """GEN-05: the low-level TS ts1 of R -> P joins R to the intermediate, so the case has two
    steps; the child R -> I validates the parent's TS first (the same freq and QRC jobs, JobStore
    hits) instead of searching again, and I -> P is searched as usual."""
    pes = fakes.triple_well()
    view, source = dft_view(tmp_run, pes)
    ts = fakes.write_geometry(tmp_run, "ts1.xyz", pes.symbols, pes.points["ts1"])
    found = DiscoveryRecord(discovery_id="d1", source_minimum=source, mechanism="nt2",
                            outcome="product", product_species="product", ts=ts)
    view.artifacts.append(Artifact(artifact_id="d1", type=T.DISCOVERY, payload=found))
    qm = KeyedQM(tmp_run, pes)
    reactions = run_stage(fake_runtime, tmp_run, pes, view, qm=qm)[0]
    parent, first, second = (reactions[r] for r in ("rx", "rx_split1", "rx_split2"))
    assert (parent.outcome, parent.reasons, parent.saddle) == (
        O.MULTI_STEP, ("intermediate_distinct",), None)
    assert (first.ts_calc, second.ts_calc) == (first.saddle.saddle_calc, None)
    assert [first.outcome, second.outcome] == [O.ELEMENTARY_STEP] * 2
    log = [json.loads(line) for line in (tmp_run / "stage" / first.log).read_text().splitlines()]
    assert [e["action"] for e in log if "action" in e] == ["validate_and_connect", "complete"]
    replayed = (first.saddle.freq_calc, *first.connection.side_calcs)
    assert [qm.jobs[c] for c in replayed] == [2, 2, 2]  # the parent's jobs, run again


def test_a_split_child_takes_the_result_of_a_queued_case_of_its_state_pair(
        tmp_run, fake_runtime) -> None:
    """G8-P7: rx splits at the triple well's intermediate; its child R -> I has the states of
    the declared ri, queued after rx, so it is not driven and takes ri's result (same_as:ri);
    I -> P is driven."""
    pes, names = fakes.triple_well(), {**ENDS, "intermediate": "intermediate"}
    view, _ = dft_view(tmp_run, pes, names)
    system = SystemConfig(
        system_id="t", species=[SpeciesInput(id=n, role="endpoint", xyz=Path(f"{n}.xyz"),
                                             multiplicity=1) for n in names],
        reactions=[*SYSTEM.reactions,
                   ReactionInput(id="ri", reactant="reactant", product="intermediate")])
    reactions = run_stage(fake_runtime, tmp_run, pes, view, system=system)[0]
    ri, (first, second) = reactions["ri"], (reactions[f"rx_split{i}"] for i in (1, 2))
    assert (reactions["rx"].outcome, ri.outcome) == (O.MULTI_STEP, O.ELEMENTARY_STEP)
    assert (first.outcome, first.reasons, first.log) == (ri.outcome, ("same_as:ri",), None)
    assert (first.saddle, first.connection) == (ri.saddle, ri.connection)
    assert second.outcome is O.ELEMENTARY_STEP and second.log is not None


def test_low_level_ts_shortcut_and_negative_discoveries_do_not_veto(
        tmp_run, fake_runtime) -> None:
    pes = fakes.double_well()
    view, source = dft_view(tmp_run, pes)
    ts = fakes.write_geometry(tmp_run, "ts.xyz", pes.symbols, pes.points["ts"])
    found = DiscoveryRecord(discovery_id="d1", source_minimum=source, mechanism="nt2",
                            outcome="product", product_species="product", ts=ts)
    trial = ReactionTrial(trial_id="t", source_minimum=source, kind="transfer",
                          associations=((1, 2),), dissociations=((0, 1),))
    negative = DiscoveryRecord(discovery_id="d2", source_minimum=source, mechanism="nt2",
                               outcome="negative", reason="no_nt2_maximum", trial=trial)
    for record in (found, negative):  # X1: a matching negative discovery changes nothing
        view.artifacts.append(Artifact(artifact_id=record.discovery_id, type=T.DISCOVERY,
                                       payload=record))
        reactions, _, engines = run_stage(fake_runtime, tmp_run, pes, view)
        assert reactions["rx"].outcome is O.ELEMENTARY_STEP
        assert not engines[Cap.PATH, "pysis_neb"].calls  # the low-level TS shortcut
