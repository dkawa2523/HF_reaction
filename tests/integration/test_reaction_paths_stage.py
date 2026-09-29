"""reaction-paths stage on fake engines (§7.3, §8.2), including the CH-05 regression."""

import json
from collections import Counter
from functools import partial
from pathlib import Path

import fakes
import pytest

from hfauto.backends.protocols import Capability as Cap
from hfauto.chemistry.xyz import composition_key
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType as T
from hfauto.core.records import CaseOutcome as O
from hfauto.core.records import DiscoveryRecord, ReactionRecord, ReactionTrial, SpeciesRecord
from hfauto.core.system import ReactionInput, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import Registry, calc_id, relax_to_minimum
from hfauto.stages import reaction_paths
from hfauto.stages.reaction_paths import ReactionPathsConfig, ReactionPathsStage

pytestmark = pytest.mark.integration
DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
XTB = MethodSpec(id="gfn2", kind="xtb", gfn=2)
ENDS = {"reactant": "reactant", "product": "product"}  # species id -> PES point
SYSTEM = SystemConfig(  # declared endpoints name an xyz, which this stage never reads
    system_id="t", reactions=[ReactionInput(id="rx", reactant="reactant", product="product")],
    species=[SpeciesInput(id=n, role="endpoint", xyz=Path(f"{n}.xyz")) for n in ENDS])
FAKES = {(Cap.QM, "nwchem"): fakes.FakeQM, (Cap.PATH, "nwchem_string"): fakes.FakePath,
         (Cap.QM, "xtb"): fakes.FakeQM, (Cap.PATH, "pysis_neb"): partial(fakes.FakePath, tsopt=True)}


class CollapsingSaddle(fakes.FakeSaddle):  # every saddle search falls into the intermediate
    def refine(self, seed, method, *, hessian, mode, deadline=None):
        key = self._key("collapse", seed.fingerprint())
        return self._evidence("saddle", seed, method, key, self._start(seed, key),
                              self.pes.points["intermediate"])


class KeyedQM(fakes.FakeQM):  # the calculation id of every result, as a JobStore sees it
    def __init__(self, root, pes):
        super().__init__(root, pes)
        self.jobs = Counter()

    def frequencies(self, mol, method, *, deadline=None):
        ev = super().frequencies(mol, method, deadline=deadline)
        self.jobs[calc_id(ev)] += 1
        return ev

    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        ev = super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)
        self.jobs[calc_id(ev)] += 1
        return ev


class SlowQM(fakes.FakeQM):  # an optimization uses up the rest of its hypothesis' walltime
    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        if deadline is not None:
            deadline.end = 0.0
        return super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)


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


def test_a_case_that_raises_is_unresolved_and_the_next_case_runs(
        tmp_run, fake_runtime, monkeypatch) -> None:
    real = reaction_paths.drive_case

    def drive_case(case, rt, rules, deadline):
        if case.reaction_id == "rx":
            raise KeyError("Te")
        return real(case, rt, rules, deadline)

    monkeypatch.setattr(reaction_paths, "drive_case", drive_case)
    pes, again = fakes.double_well(), ReactionInput(id="rx2", reactant="reactant", product="product")
    system = SYSTEM.model_copy(update={"reactions": [*SYSTEM.reactions, again]})
    view = dft_view(tmp_run, pes)[0]
    with pytest.raises(KeyError):  # HFAUTO_STRICT=1 re-raises
        run_stage(fake_runtime, tmp_run, pes, view, system=system)
    monkeypatch.delenv("HFAUTO_STRICT")
    reactions = run_stage(fake_runtime, tmp_run, pes, view, system=system)[0]
    rx = reactions["rx"]
    assert rx.outcome is O.UNRESOLVED and rx.reasons == ("error:KeyError",)
    assert json.loads((tmp_run / "stage" / rx.log).read_text()) == {
        "action": "error", "reason": "error:KeyError", "detail": "'Te'"}
    assert reactions["rx2"].outcome is O.ELEMENTARY_STEP


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


def test_split_children_run_on_the_rest_of_their_hypothesis_walltime(tmp_run,
                                                                    fake_runtime) -> None:
    """U9-P5: the well found at the parent's deadline still splits the case (evidence before the
    walltime row); its children inherit the expired deadline instead of 6 h each."""
    pes = fakes.triple_well()
    view, _ = dft_view(tmp_run, pes)
    reactions = run_stage(fake_runtime, tmp_run, pes, view, qm=SlowQM(tmp_run, pes))[0]
    assert reactions["rx"].outcome is O.MULTI_STEP
    children = [reactions[f"rx_split{i}"] for i in (1, 2)]
    assert [(c.outcome, c.reasons) for c in children] == [(O.UNRESOLVED, ("walltime",))] * 2


def test_collapsed_saddle_is_validated_as_an_intermediate(tmp_run, fake_runtime) -> None:
    pes = fakes.triple_well()
    view, source = dft_view(tmp_run, pes)
    ts = fakes.write_geometry(tmp_run, "ts1.xyz", pes.symbols, pes.points["ts1"])
    found = DiscoveryRecord(discovery_id="d1", source_minimum=source, mechanism="nt2",
                            outcome="product", product_species="product", ts=ts)
    view.artifacts.append(Artifact(artifact_id="d1", type=T.DISCOVERY, payload=found))
    reactions, _, _ = run_stage(fake_runtime, tmp_run, pes, view,  # seeded by the shortcut
                                saddle=CollapsingSaddle(tmp_run, pes), max_split_depth=0)
    assert reactions["rx"].outcome is O.MULTI_STEP and "rx_split1" not in reactions
    assert '"saddle_collapsed"' in (tmp_run / "stage" / reactions["rx"].log).read_text()


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
    assert [e["action"] for e in log if "action" in e] == ["validate_ts", "connect", "complete"]
    replayed = (first.saddle.freq_calc, *first.connection.side_calcs)
    assert [qm.jobs[c] for c in replayed] == [2, 2, 2]  # the parent's jobs, run again


def test_walltime_low_level_ts_shortcut_and_negative_discoveries_do_not_veto(
        tmp_run, fake_runtime) -> None:
    pes = fakes.double_well()
    view, source = dft_view(tmp_run, pes)
    rx = run_stage(fake_runtime, tmp_run, pes, view, walltime_h=0.0)[0]["rx"]
    assert rx.outcome is O.UNRESOLVED and rx.reasons == ("walltime",)
    ts = fakes.write_geometry(tmp_run, "ts.xyz", pes.symbols, pes.points["ts"])
    found = DiscoveryRecord(discovery_id="d1", source_minimum=source, mechanism="nt2",
                            outcome="product", product_species="product", ts=ts)
    trial = ReactionTrial(trial_id="t", source_minimum=source, kind="transfer", mechanism="afir",
                          associations=((1, 2),), dissociations=((0, 1),))
    negative = DiscoveryRecord(discovery_id="d2", source_minimum=source, mechanism="afir",
                               outcome="negative", reason="no_nt2_maximum", trial=trial)
    for record in (found, negative):  # X1: a matching negative discovery changes nothing
        view.artifacts.append(Artifact(artifact_id=record.discovery_id, type=T.DISCOVERY,
                                       payload=record))
        reactions, _, engines = run_stage(fake_runtime, tmp_run, pes, view)
        assert reactions["rx"].outcome is O.ELEMENTARY_STEP
        assert not engines[Cap.PATH, "pysis_neb"].calls  # the low-level TS shortcut
