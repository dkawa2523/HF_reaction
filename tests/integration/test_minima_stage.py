"""minima stage on fake surfaces (§4.1 #3, §7.2); K cases of test_minimum_mode_follow_stage."""

import fakes
import numpy as np
from scipy.spatial.distance import pdist

from hfauto.backends.protocols import Capability
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import composition_key
from hfauto.core.evidence import Evidence
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType as T
from hfauto.core.records import DiscoveryRecord, MinimumRecord, SpeciesRecord
from hfauto.core.system import SystemConfig
from hfauto.stages.minima import MinimaConfig, MinimaStage

METHODS = {"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2),
           "pbe0": MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")}
SCREEN = {"level": "screen", "engine": "xtb", "method": "gfn2"}
DFT = {"level": "dft", "engine": "nwchem", "method": "pbe0"}


class SoftQM(fakes.FakeQM):
    def frequencies(self, mol, method, *, deadline=None):  # a stuck −30 cm⁻¹ mode
        ev = super().frequencies(mol, method, deadline=deadline)
        return ev.model_copy(update={"frequencies_cm1": (-30.0, *ev.frequencies_cm1[1:]),
                                     "imaginary_modes": (tuple(np.eye(9)[3]),)})


class StuckQM(fakes.FakeQM):
    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        if not self.calls:  # the first input converges to the barrier top (a saddle endpoint)
            mol = self.pes.molecule("ts")
        return super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)


def species(root, pes, sid, point, shift=0.0):
    x = pes.points[point] + shift
    record = SpeciesRecord(species_id=sid, composition_id=composition_key(pes.symbols, 0, 1),
                           charge=0, multiplicity=1, source="input",
                           geometry=fakes.write_geometry(root, f"in/{sid}.xyz", pes.symbols, x),
                           state_label=state_label(pes.symbols, x))
    return Artifact(artifact_id=sid, type=T.SPECIES, payload=record)


def stage(fake_runtime, tmp_run, pes, low=fakes.FakeQM):  # run(id, view, **config) -> view+out
    xtb, dft = low(tmp_run, pes), fakes.FakeQM(tmp_run, pes)
    rt = fake_runtime(SystemConfig(system_id="t", species=[]),
                      {(Capability.QM, "xtb"): xtb, (Capability.QM, "nwchem"): dft}, METHODS)

    def run(stage_id, artifacts, **config):
        bound = rt.bind(stage_id, tmp_run / stage_id)
        bound.stage_dir.mkdir()
        view = Manifest(run_id="r", stage_id=stage_id, created_at="t", artifacts=artifacts)
        out = MinimaStage().run(view, MinimaConfig.model_validate(config), bound)
        produced = Manifest(run_id="r", stage_id=stage_id, created_at="t", artifacts=out)
        return Manifest.union([view, produced], run_id="r", stage_id=stage_id)

    return run, xtb, dft


def test_one_basin_gets_one_freq_job_per_tier_and_dft_starts_from_screen(fake_runtime, tmp_run):
    pes = fakes.harmonic()
    run, xtb, dft = stage(fake_runtime, tmp_run, pes)
    inputs = [species(tmp_run, pes, "a", "start"), species(tmp_run, pes, "b", "start", 0.01)]
    screen = run("screen", inputs, **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    assert low.tier == "screen" and low.members == ("a", "b")
    assert xtb.calls == ["optimize", "frequencies", "optimize"]  # b is known: no second freq
    tasks = [ev.task for ev in screen.records(T.CALCULATION, Evidence)]
    assert (tasks.count("opt"), tasks.count("freq")) == (2, 1)
    both = run("dft", screen.artifacts, **DFT, select={"include": "all"})
    (high,) = [m for m in both.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.members == ("a", "b") and dft.calls == ["optimize", "frequencies", "optimize"]
    dft.calls.clear()
    window = run("dft_window", screen.artifacts, **DFT)
    (high,) = [m for m in window.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.species_id == "a" and dft.calls == ["optimize", "frequencies"]
    start = window.evidence(high.opt_calc).start  # dft starts from the screen structure
    assert start.fingerprint == screen.evidence(low.opt_calc).final.fingerprint


def test_saddle_structure_gives_a_mode_follow_discovery(fake_runtime, tmp_run):
    pes = fakes.double_well()
    run = stage(fake_runtime, tmp_run, pes)[0]
    out = run("screen", [species(tmp_run, pes, "t", "ts")], **SCREEN)
    (found,) = out.records(T.DISCOVERY, DiscoveryRecord)
    sides = {m.species_id: m for m in out.records(T.MINIMUM, MinimumRecord)}
    assert set(sides) == {"t_mf1", "t_mf2"} and found.mechanism == "mode_follow"
    assert found.source_minimum == sides["t_mf1"].minimum_id and found.product_species == "t_mf2"
    assert found.ts is not None and found.ts_imag_cm1 < -50
    assert {s.source for s in out.records(T.SPECIES, SpeciesRecord)} == {"input", "mode_follow"}
    assert all(a.status == "success" for a in out.artifacts)  # no failed minimum for "t"


def test_an_endpoint_on_a_saddle_joins_the_side_nearer_its_input(fake_runtime, tmp_run):
    pes = fakes.double_well()  # declared as the product (H at O), above the reactant
    run = stage(fake_runtime, tmp_run, pes, low=StuckQM)[0]
    out = run("screen", [species(tmp_run, pes, "p", "product")], **SCREEN)
    minima = out.records(T.MINIMUM, MinimumRecord)
    (joined,) = [m for m in minima if "p" in m.members]
    assert joined.energy_hartree == max(m.energy_hartree for m in minima)  # nearer, not lower
    assert joined.notes == ("endpoint_was_saddle",) and len(minima) == 2
    assert len(out.records(T.DISCOVERY, DiscoveryRecord)) == 1


def test_soft_mode_known_minimum_and_init_hessian(fake_runtime, tmp_run):
    pes = fakes.double_well()  # the reactant is N-H + O: two fragments
    run, xtb, dft = stage(fake_runtime, tmp_run, pes, low=SoftQM)
    soft = run("screen", [species(tmp_run, pes, "r", "reactant")], **SCREEN)
    (record,) = soft.records(T.MINIMUM, MinimumRecord)
    assert record.notes == ("soft_imaginary_mode",)
    init = {"init_hessian": {"engine": "xtb", "method": "gfn2"}, "select": {"include": "all"}}
    first = run("dft", soft.artifacts, **DFT, **init)
    assert dft.calls == ["optimize+init_hessian", "frequencies"] and "frequencies" in xtb.calls
    dft.calls.clear()
    again = run("dft2", first.artifacts, **DFT, **init)
    assert dft.calls == ["optimize+init_hessian"]  # known basin: no freq job
    assert {m.members for m in again.records(T.MINIMUM, MinimumRecord)} == {("r",)}  # joined


def test_an_atom_is_relaxed_and_registered_like_any_species(fake_runtime, tmp_run):
    pes = fakes.PES(("Ar",), lambda x: -527.0, {"atom": np.zeros((1, 3))})
    run, xtb, dft = stage(fake_runtime, tmp_run, pes)
    screen = run("screen", [species(tmp_run, pes, "ar", "atom")], **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    freq = screen.evidence(low.freq_calc)
    assert xtb.calls == ["optimize", "frequencies"] and low.notes == ()
    assert (freq.frequencies_cm1, freq.n_external) == ((), 3)  # 3N - 3 = 0 modes
    both = run("dft", screen.artifacts, **DFT, select={"include": "all"})
    (high,) = [m for m in both.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.members == ("ar",) and dft.calls == ["optimize", "frequencies"]


def test_enantiomers_fall_into_one_chiral_basin(fake_runtime, tmp_run):
    ref = np.array([[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0],
                    [1.1, -1.1, -1.1]])  # CHFClBr held by pairwise springs
    pes = fakes.PES(("C", "H", "F", "Cl", "Br"),
                    lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - pdist(ref)) ** 2)),
                    {"r": ref, "s": ref * [-1.0, 1.0, 1.0]})
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    screen = run("screen", [species(tmp_run, pes, "r", "r"), species(tmp_run, pes, "s", "s")],
                 **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    assert low.chiral and low.members == ("r", "s")
    assert xtb.calls == ["optimize", "frequencies", "optimize"]  # s is known: no freq job
