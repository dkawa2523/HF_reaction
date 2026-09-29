"""minima stage on fake surfaces (§4.1 #3, §7.2); K cases of test_minimum_mode_follow_stage;
R6 relaxation seeds."""

import json

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
from hfauto.stages.minima import MinimaConfig, MinimaStage, _image_groups, _Job

METHODS = {"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2),
           "pbe0": MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")}
SCREEN = {"level": "screen", "engine": "xtb", "method": "gfn2"}
DFT = {"level": "dft", "engine": "nwchem", "method": "pbe0"}


class SoftQM(fakes.FakeQM):
    def frequencies(self, mol, method, **kw):  # a stuck −30 cm⁻¹ mode
        ev = super().frequencies(mol, method, **kw)
        return ev.model_copy(update={"frequencies_cm1": (-30.0, *ev.frequencies_cm1[1:]),
                                     "imaginary_modes": (tuple(np.eye(9)[3]),)})


class CollapseQM(fakes.FakeQM):
    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        if np.allclose(mol.xyz.coords, self.pes.points["reactant"]):  # no barrier (GFN2 in S6)
            mol = self.pes.molecule("product")
        return super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)


class StuckQM(fakes.FakeQM):
    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        if not self.calls:  # the first input converges to the barrier top (a saddle endpoint)
            mol = self.pes.molecule("ts")
        return super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)


def species_at(root, sid, symbols, x):
    record = SpeciesRecord(species_id=sid, composition_id=composition_key(symbols, 0, 1),
                           charge=0, multiplicity=1, source="input",
                           geometry=fakes.write_geometry(root, f"in/{sid}.xyz", symbols, x),
                           state_label=state_label(symbols, x))
    return Artifact(artifact_id=sid, type=T.SPECIES, payload=record)


def species(root, pes, sid, point):
    return species_at(root, sid, pes.symbols, pes.points[point])


def stage(fake_runtime, tmp_run, pes, low=fakes.FakeQM, system=None, high=fakes.FakeQM):
    xtb, dft = low(tmp_run, pes), high(tmp_run, pes)  # run(id, view, **cfg) -> view + output
    rt = fake_runtime(system or SystemConfig(system_id="t", species=[]),
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
    pes = fakes.harmonic()  # "a" is a declared endpoint: its composition reacts
    system = SystemConfig(system_id="t", species=[{"id": "a", "xyz": "a.xyz", "role": "endpoint"}])
    run, xtb, dft = stage(fake_runtime, tmp_run, pes, system=system)
    inputs = [species(tmp_run, pes, "a", "start"), species(tmp_run, pes, "b", "minimum")]
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
    assert found.ts_calc is None  # an xTB saddle is a low-level TS
    assert {s.source for s in out.records(T.SPECIES, SpeciesRecord)} == {"input", "mode_follow"}
    assert all(a.status == "success" for a in out.artifacts)  # no failed minimum for "t"

    dft = run("dft", [species(tmp_run, pes, "t", "ts")], **DFT, select={"include": "all"})
    (verified,) = dft.records(T.DISCOVERY, DiscoveryRecord)
    saddle = dft.evidence(verified.ts_calc or "")  # the DFT saddle's opt, validated directly
    assert (saddle.task, saddle.level.method, saddle.final) == ("opt", "pbe0", verified.ts)


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
    nudge = np.zeros((5, 3))
    nudge[1, 0] = 0.03  # s starts near the mirror image, not at it: its own opt runs
    pes = fakes.PES(("C", "H", "F", "Cl", "Br"),
                    lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - pdist(ref)) ** 2)),
                    {"r": ref, "s": ref * [-1.0, 1.0, 1.0] + nudge})
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    screen = run("screen", [species(tmp_run, pes, "r", "r"), species(tmp_run, pes, "s", "s")],
                 **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    assert low.chiral and low.members == ("r", "s")
    assert xtb.calls == ["optimize", "frequencies", "optimize"]  # s is known: no freq job


def diagnostics(tmp_run, stage_id):
    return json.loads((tmp_run / stage_id / "diagnostics.json").read_text(encoding="utf-8"))


def test_mirror_sides_of_a_symmetric_saddle_give_one_basin_and_no_job(fake_runtime, tmp_run):
    pes = fakes.symmetric_double_well()  # F-H-F: the two sides are permutation images
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    out = run("screen", [species(tmp_run, pes, "t", "ts")], **SCREEN)
    (basin,) = out.records(T.MINIMUM, MinimumRecord)
    assert basin.members == ("t_mf1", "t_mf2", "t") and basin.notes == ("endpoint_was_saddle",)
    assert len(xtb.calls) == 6 and xtb.calls.count("frequencies") == 3  # the driver's jobs only
    (found,) = out.records(T.DISCOVERY, DiscoveryRecord)
    assert (found.source_minimum, found.product_species) == (basin.minimum_id, "t_mf2")


def test_an_exact_image_joins_a_minimum_with_no_job_and_settles_after_a_saddle(fake_runtime,
                                                                               tmp_run):
    pes = fakes.symmetric_double_well()  # reactant and product: F1 and F2 relabelled
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    inputs = [species(tmp_run, pes, "a", "reactant"), species(tmp_run, pes, "b", "product")]
    (basin,) = run("screen", inputs, **SCREEN).records(T.MINIMUM, MinimumRecord)
    assert basin.members == ("a", "b") and xtb.calls == ["optimize", "frequencies"]
    assert diagnostics(tmp_run, "screen")["b"] == ["image_of:a"]
    run, soft, _ = stage(fake_runtime, tmp_run, pes, low=SoftQM)  # a soft minimum: b settles
    run("soft", inputs, **SCREEN)
    assert diagnostics(tmp_run, "soft")["a"][0] == "soft_minimum"
    assert diagnostics(tmp_run, "soft")["b"][0] == "known" and soft.calls[-1] == "optimize"

    well = fakes.double_well()  # p's first opt stops on the saddle; q is p rotated by 180°
    run, stuck, _ = stage(fake_runtime, tmp_run, well, low=StuckQM)
    p = species(tmp_run, well, "p", "product")
    q = species_at(tmp_run, "q", well.symbols, well.points["product"] * [-1.0, 1.0, 1.0])
    jobs = [_Job(a.payload, a.payload.geometry) for a in (p, q)]
    assert [len(g) for g in _image_groups(jobs, fakes.xyz_loader(tmp_run))] == [2]
    out = run("saddle", [p, q], **SCREEN)
    (joined,) = [m for m in out.records(T.MINIMUM, MinimumRecord) if "p" in m.members]
    assert joined.notes == ("endpoint_was_saddle",) and "q" in joined.members
    assert diagnostics(tmp_run, "saddle")["q"][:2] == ["known", "opt"]  # its own opt
    assert stuck.calls[-1] == "optimize"


def test_only_reacting_compositions_and_their_monomers_are_refined(fake_runtime, tmp_run):
    tetra = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * 0.34
    springs = lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - 0.92) ** 2))
    pes = fakes.PES(("H", "F"), springs, {})  # every pair at 0.92 A: HF, H2 and tetrahedra
    system = SystemConfig(system_id="t", species=[{"id": "hf", "xyz": "hf.xyz"},
                                                  {"id": "h2", "xyz": "h2.xyz"}],
                          compositions=[{"id": "hf2", "components": {"hf": 2}},
                                        {"id": "hf_h2", "components": {"hf": 1, "h2": 1}}])
    run, _, dft = stage(fake_runtime, tmp_run, pes, system=system)
    inputs = [species_at(tmp_run, "hf", ("H", "F"), [[0, 0, 0], [0, 0, 0.95]]),
              species_at(tmp_run, "h2", ("H", "H"), [[0, 0, 0], [0, 0, 0.9]]),
              species_at(tmp_run, "hf2_c00", ("H", "F", "H", "F"), tetra + 0.02 * np.eye(4, 3)),
              species_at(tmp_run, "hf_h2_c00", ("H", "F", "H", "H"), tetra)]
    screen = run("screen", inputs, **SCREEN)
    (source,) = [m for m in screen.records(T.MINIMUM, MinimumRecord) if m.species_id == "hf2_c00"]
    found = DiscoveryRecord(discovery_id="disc", source_minimum=source.minimum_id,
                            mechanism="nt2", outcome="product", product_species="hf2_c00")
    disc = Artifact(artifact_id="disc", type=T.DISCOVERY, payload=found)
    out = run("dft", [*screen.artifacts, disc], **DFT, select={"rerank_sp": True})
    refined = {m.composition_id for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"}
    assert refined == {"F2H2_q0_m1", "FH_q0_m1"}  # (HF)2 reacts; HF is its monomer
    skipped = {k for k, v in diagnostics(tmp_run, "dft").items() if v == ["not_reacting"]}
    assert skipped == {"H2_q0_m1", "FH3_q0_m1"}
    assert "energy" not in dft.calls  # one candidate per group: no rerank single point


def test_a_relaxation_seed_is_asked_at_dft_once_from_its_own_geometry(fake_runtime, tmp_run):
    """R6 (S5, S6): the seed collapsed at screen and represents that basin, which DFT refines
    from the screen structure; the seed is asked once more, from its own geometry with no xTB
    start Hessian (no xTB stationary point), as its own species after the other jobs. A seed
    that collapses at DFT too is noted."""
    pes = fakes.double_well()  # the seed at the reactant (H at N); screen has no barrier to H-O
    seed = species(tmp_run, pes, "seed", "reactant")
    init = {"init_hessian": {"engine": "xtb", "method": "gfn2"}}  # both wells: two fragments
    for high in (fakes.FakeQM, CollapseQM):
        run, _, dft = stage(fake_runtime, tmp_run, pes, low=CollapseQM, high=high)
        screen = run(f"screen_{high.__name__}", [seed], **SCREEN)
        (lost,) = screen.records(T.MINIMUM, MinimumRecord)
        assert lost.species_id == "seed" and lost.state_label != seed.payload.state_label
        relax = DiscoveryRecord(discovery_id="relax_seed", source_minimum=lost.minimum_id,
                                mechanism="relaxation", outcome="product", product_species="seed")
        found = Artifact(artifact_id="relax_seed", type=T.DISCOVERY, payload=relax)
        out = run(f"dft_{high.__name__}", [*screen.artifacts, found], **DFT, **init)
        minima = {m.species_id: m for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"}
        history = diagnostics(tmp_run, f"dft_{high.__name__}")["spc_relax_seed"]
        own = next(s for s in out.records(T.SPECIES, SpeciesRecord)
                   if s.species_id == "spc_relax_seed")
        assert own.geometry == seed.payload.geometry
        if high is fakes.FakeQM:  # PBE0 keeps the seed state: two basins, two freq jobs
            assert set(minima) == {"seed", "spc_relax_seed"} and "collapsed" not in str(history)
            assert minima["spc_relax_seed"].state_label == seed.payload.state_label
            start = out.evidence(minima["spc_relax_seed"].opt_calc).start
            assert start.fingerprint == seed.payload.geometry.fingerprint
            assert dft.calls == ["optimize+init_hessian", "frequencies", "optimize", "frequencies"]
        else:  # the seed joins the collapse basin with no freq job
            assert minima["seed"].members == ("seed", "spc_relax_seed")
            assert history[-1] == "collapsed_at_dft_from_seed"
            assert dft.calls == ["optimize+init_hessian", "frequencies", "optimize"]
