"""minima stage on fake surfaces (§4.1 #3, §7.2); K cases of test_minimum_mode_follow_stage."""

import fakes
import numpy as np

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


class ShiftedQM(fakes.FakeQM):
    def optimize(self, mol, method, *, tight=False, **kw):  # "b" ends 3e-5 Eh high unless tight
        ev = super().optimize(mol, method, tight=tight, **kw)
        high = not tight and mol.xyz.coords[0, 0] > 0.005
        return ev.model_copy(update={"energy_hartree": ev.energy_hartree + 3e-5}) if high else ev


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


def test_ambiguous_duplicate_is_rejudged_into_one_basin_and_tiers_split(fake_runtime, tmp_run):
    pes = fakes.harmonic()
    run, _, dft = stage(fake_runtime, tmp_run, pes, low=ShiftedQM)
    inputs = [species(tmp_run, pes, "a", "start"), species(tmp_run, pes, "b", "start", 0.01)]
    screen = run("screen", inputs, **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    assert low.tier == "screen" and low.members == ("a", "b")
    assert [ev.task for ev in screen.records(T.CALCULATION, Evidence)].count("opt") == 3  # tight
    both = run("dft", screen.artifacts, **DFT)
    (high,) = [m for m in both.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.species_id == "a" and dft.calls == ["optimize", "frequencies"]
    start = both.evidence(high.opt_calc).start  # dft starts from the screen structure
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
