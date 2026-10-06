"""minima stage on fake surfaces (§4.1 #3, §7.2); K cases of test_minimum_mode_follow_stage;
the DFT entry of the low-level edges (U4-P3); R6 relaxation seeds."""

import json
from dataclasses import replace

import fakes
import numpy as np
import pytest
from scipy.spatial.distance import pdist

from hfauto.backends.protocols import Capability
from hfauto.chemistry.gates import Policy
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import Evidence, Failure, FailureKind
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType as T
from hfauto.core.records import DiscoveryRecord, MinimumRecord, SpeciesRecord
from hfauto.core.system import SystemConfig
from hfauto.drivers.minimum import Registry, relax_to_minimum
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
    def optimize(self, mol, method, **kw):  # no barrier from the reactant (GFN2 in S6)
        ev = super().optimize(mol, method, **kw)
        if isinstance(ev, Evidence) and np.allclose(mol.xyz.coords, self.pes.points["reactant"]):
            x = self.pes.points["product"]
            final = fakes.write_geometry(self.root, f"fake/{ev.job_key[:16]}/product.xyz",
                                         self.pes.symbols, x)
            ev = ev.model_copy(update={"final": final, "energy_hartree": self.pes.energy(x)})
        return ev


class StalledQM(fakes.FakeQM):
    def optimize(self, mol, method, **kw):  # the first opt stops at its start; NWChem's gradient
        first = not self.calls
        ev = super().optimize(mol, method, **kw)
        if first:
            ev = ev.model_copy(update={"final": ev.start,
                                       "energy_hartree": self.pes.energy(mol.xyz.coords)})
        x = fakes.xyz_loader(self.root)(ev.final).coords
        return ev.model_copy(update={"gradient": tuple(self.pes.gradient(x) * BOHR_TO_ANGSTROM)})


class StuckQM(fakes.FakeQM):
    def optimize(self, mol, method, **kw):
        if not self.calls:  # the first input converges to the barrier top (a saddle endpoint)
            mol = self.pes.molecule("ts")
        return super().optimize(mol, method, **kw)


def species_at(root, sid, symbols, x):
    record = SpeciesRecord(species_id=sid, composition_id=composition_key(symbols, 0, 1),
                           charge=0, multiplicity=1, source="input",
                           geometry=fakes.write_geometry(root, f"in/{sid}.xyz", symbols, x),
                           state_label=state_label(symbols, x))
    return Artifact(artifact_id=sid, type=T.SPECIES, payload=record)


def species(root, pes, sid, point):
    return species_at(root, sid, pes.symbols, pes.points[point])


def stage(fake_runtime, tmp_run, pes, low=fakes.FakeQM, system=None, high=fakes.FakeQM,
          policy=None, cores=4):
    xtb, dft = low(tmp_run, pes), high(tmp_run, pes)  # run(id, view, **cfg) -> view + output
    rt = fake_runtime(system or SystemConfig(system_id="t", species=[]),
                      {(Capability.QM, "xtb"): xtb, (Capability.QM, "nwchem"): dft}, METHODS)
    rt = replace(rt, policy=policy) if policy else rt
    rt = replace(rt, site=rt.site.model_copy(update={"cores": cores}))

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
    system = SystemConfig(system_id="t", species=[{"id": "a", "xyz": "a.xyz", "role": "endpoint",
                                                   "multiplicity": 1}])
    run, xtb, dft = stage(fake_runtime, tmp_run, pes, system=system)
    inputs = [species(tmp_run, pes, "a", "start"), species(tmp_run, pes, "b", "minimum")]
    screen = run("screen", inputs, **SCREEN)
    (low,) = screen.records(T.MINIMUM, MinimumRecord)
    assert low.tier == "screen" and low.members == ("a", "b")
    assert xtb.calls == ["optimize", "optimize", "frequencies"]  # b is known: no second freq
    tasks = [ev.task for ev in screen.records(T.CALCULATION, Evidence)]
    assert (tasks.count("opt"), tasks.count("freq")) == (2, 1)
    both = run("dft", screen.artifacts, **DFT, select={"include": "all"})
    (high,) = [m for m in both.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.members == ("a", "b") and dft.calls == ["optimize", "optimize", "frequencies"]
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
    assert (found.source_species, found.product_species) == ("t_mf1", "t_mf2")  # its own ends
    assert found.ts is not None and found.dE_act_kcal > 0
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
    assert low.members == ("r", "s")
    assert xtb.calls == ["optimize", "optimize", "frequencies"]  # s is known: no freq job


def diagnostics(tmp_run, stage_id):
    return json.loads((tmp_run / stage_id / "diagnostics.json").read_text(encoding="utf-8"))


def test_mirror_sides_of_a_symmetric_saddle_give_one_basin_and_no_job(fake_runtime, tmp_run):
    pes = fakes.symmetric_double_well()  # F-H-F: the two sides are permutation images
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    out = run("screen", [species(tmp_run, pes, "t", "ts")], **SCREEN)
    (basin,) = out.records(T.MINIMUM, MinimumRecord)
    assert basin.members == ("t_mf1", "t_mf2", "t") and basin.notes == ("endpoint_was_saddle",)
    assert xtb.calls == ["optimize", "frequencies", "optimize+init_hessian",
                         "optimize+init_hessian", "frequencies"]  # the mirror side shares a freq
    (found,) = out.records(T.DISCOVERY, DiscoveryRecord)
    assert (found.source_species, found.product_species) == ("t_mf1", "t_mf2")


def test_an_exact_image_joins_a_minimum_with_no_job_and_settles_after_a_saddle(fake_runtime,
                                                                               tmp_run):
    pes = fakes.symmetric_double_well()  # reactant and product: F1 and F2 relabelled
    run, xtb, _ = stage(fake_runtime, tmp_run, pes)
    inputs = [species(tmp_run, pes, "a", "reactant"), species(tmp_run, pes, "b", "product")]
    (basin,) = run("screen", inputs, **SCREEN).records(T.MINIMUM, MinimumRecord)
    assert basin.members == ("a", "b") and xtb.calls == ["optimize", "frequencies"]
    assert diagnostics(tmp_run, "screen")["b"] == ["image_of:a"]
    run = stage(fake_runtime, tmp_run, pes, low=SoftQM)[0]  # a certified soft point
    run("soft", inputs, **SCREEN)  # is stationary and not pushed: its image joins with no job
    assert diagnostics(tmp_run, "soft")["a"] == ["minimum", "opt:reused", "freq:soft"]
    assert diagnostics(tmp_run, "soft")["b"] == ["image_of:a"]

    well = fakes.double_well()  # p's first opt stops on the saddle; q is p rotated by 180°
    run, stuck, _ = stage(fake_runtime, tmp_run, well, low=StuckQM)
    p = species(tmp_run, well, "p", "product")
    q = species_at(tmp_run, "q", well.symbols, well.points["product"] * [-1.0, 1.0, 1.0])
    jobs = [_Job(a.payload, a.payload.geometry) for a in (p, q)]
    assert [len(g) for g in _image_groups(jobs, fakes.xyz_loader(tmp_run))] == [2]
    out = run("saddle", [p, q], **SCREEN)
    (joined,) = [m for m in out.records(T.MINIMUM, MinimumRecord) if "p" in m.members]
    assert joined.notes == ("endpoint_was_saddle",) and "q" in joined.members
    assert diagnostics(tmp_run, "saddle")["q"][:2] == ["known", "opt:reused"]  # its own opt
    assert stuck.calls[-1] == "optimize"


def test_only_reacting_compositions_and_their_monomers_are_refined(fake_runtime, tmp_run):
    tetra = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * 0.34

    def springs(x):
        return float(np.sum((pdist(np.reshape(x, (-1, 3))) - 0.92) ** 2))

    pes = fakes.PES(("H", "F"), springs, {})  # every pair at 0.92 A: HF, H2 and tetrahedra
    system = SystemConfig(system_id="t", species=[{"id": "hf", "xyz": "hf.xyz", "multiplicity": 1},
                                                  {"id": "h2", "xyz": "h2.xyz", "multiplicity": 1}],
                          compositions=[{"id": "hf2", "components": {"hf": 2}},
                                        {"id": "hf_h2", "components": {"hf": 1, "h2": 1}}])
    run, _, dft = stage(fake_runtime, tmp_run, pes, system=system)
    inputs = [species_at(tmp_run, "hf", ("H", "F"), [[0, 0, 0], [0, 0, 0.95]]),
              species_at(tmp_run, "h2", ("H", "H"), [[0, 0, 0], [0, 0, 0.9]]),
              species_at(tmp_run, "hf2_c00", ("H", "F", "H", "F"), tetra + 0.02 * np.eye(4, 3)),
              species_at(tmp_run, "hf_h2_c00", ("H", "F", "H", "H"), tetra)]
    screen = run("screen", inputs, **SCREEN)
    found = DiscoveryRecord(discovery_id="disc", mechanism="nt2", outcome="product",
                            source_species="hf2_c00", product_species="hf2_c00")
    disc = Artifact(artifact_id="disc", type=T.DISCOVERY, payload=found)
    out = run("dft", [*screen.artifacts, disc], **DFT)
    refined = {m.composition_id for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"}
    assert refined == {"F2H2_q0_m1", "FH_q0_m1"}  # (HF)2 reacts; HF is its monomer
    skipped = {k for k, v in diagnostics(tmp_run, "dft").items() if v == ["not_reacting"]}
    assert skipped == {"H2_q0_m1", "FH3_q0_m1"}
    assert dft.calls.count("energy") == 1  # the edge's one structure; no crowded group


def test_a_crowded_state_is_reranked_by_single_points_at_its_screen_structures(fake_runtime,
                                                                                tmp_run):
    """Two bent conformers of one state (100° and 140°) and per_state 1: each gets a single
    point at its screen minimum's structure, and only the lower one is refined. The declared
    endpoint (three atoms apart, not screened) is always kept: it stands for itself and
    relaxes into the refined basin, after it."""
    cos = np.cos(np.radians([100.0, 140.0]))

    def bent(x):
        o, a, b = np.reshape(x, (3, 3))
        ra, rb = np.linalg.norm(a - o), np.linalg.norm(b - o)
        c = float((a - o) @ (b - o)) / (ra * rb)
        return float(0.5 * ((ra - 0.96) ** 2 + (rb - 0.96) ** 2)
                     + ((c - cos[0]) * (c - cos[1])) ** 2 + 1e-3 * c)  # 140° lies lower

    pes = fakes.PES(("O", "H", "H"), bent, {
        f"w{i}": np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0], [r * c, r * s, 0.0]])
        for i, (r, c, s) in enumerate(zip((0.96, 0.96, 2.5), (*cos, cos[1]),
                                           (*np.sin(np.arccos(cos)), np.sin(np.arccos(cos[1]))),
                                           strict=True), start=1)})
    endpoint = {"id": "w3", "xyz": "w3.xyz", "role": "endpoint", "multiplicity": 1}  # reacts
    system = SystemConfig(system_id="t", species=[endpoint])
    run, _, dft = stage(fake_runtime, tmp_run, pes, system=system)
    screen = run("screen", [species(tmp_run, pes, w, w) for w in ("w1", "w2")], **SCREEN)
    basins = screen.records(T.MINIMUM, MinimumRecord)
    assert len(basins) == 2 and len({m.state_label for m in basins}) == 1
    out = run("dft", [*screen.artifacts, species(tmp_run, pes, "w3", "w3")], **DFT,
              select={"per_state": 1})
    sps = [ev for ev in out.records(T.CALCULATION, Evidence) if ev.task == "sp"]
    starts = {ev.start.fingerprint for ev in sps}
    assert starts == {screen.evidence(m.opt_calc).final.fingerprint for m in basins}
    (high,) = [m for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"]
    assert high.members == ("w2", "w3")
    assert dft.calls == ["energy", "energy", "optimize", "optimize", "frequencies"]


def test_a_minimum_short_of_stationary_relaxes_once_from_its_trust_region_step(fake_runtime,
                                                                                tmp_run):
    """X1-2 through the stage: the runtime opens the freq Hessian, so a DFT opt stopped short
    of the minimum, with no imaginary mode, fails the trust-region certification and relaxes
    once from x + s into the certified minimum."""
    pes = fakes.harmonic()
    run = stage(fake_runtime, tmp_run, pes, high=StalledQM)[0]
    out = run("dft", [species(tmp_run, pes, "a", "start")], **DFT, select={"include": "all"})
    assert diagnostics(tmp_run, "dft")["a"] == ["minimum", "opt:reused", "freq:none",
                                                "tr_relax", "freq:none"]
    (record,) = out.records(T.MINIMUM, MinimumRecord)
    assert record.notes == () and record.energy_hartree == pytest.approx(0.0, abs=1e-8)


class ScfFailsAtTheTS(fakes.FakeQM):
    def energy(self, mol, method, **kw):
        if np.allclose(mol.xyz.coords, self.pes.points["ts"]):
            self.calls.append("energy")
            return Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf_unavailable:diverged")
        return super().energy(mol, method, **kw)


def test_the_dft_entry_admits_an_edge_by_single_points_at_its_ends_and_ts(fake_runtime,
                                                                          tmp_run):
    """U4-P3: one single point at each end's start and at the TS; the TS above the reactant
    (the start) is judged against the reaction window. Out of the window, or with a failed
    single point, the edge is recorded again as a negative and its product is not refined."""
    pes = fakes.double_well()  # barrier 0.01 Eh = 6.3 kcal/mol above the reactant
    found = DiscoveryRecord(discovery_id="edge", mechanism="nt2", outcome="product",
                            source_species="r", product_species="p",
                            ts=fakes.write_geometry(tmp_run, "in/ts.xyz", pes.symbols,
                                                    pes.points["ts"]))
    edge = Artifact(artifact_id="edge", type=T.DISCOVERY, payload=found)
    cases = (({}, None), ({"policy": Policy(reaction_window_kcal=5.0)}, "out_of_window"),
             ({"high": ScfFailsAtTheTS}, "not_evaluated:scf_not_converged"))
    for i, (kw, reason) in enumerate(cases):
        run, _, dft = stage(fake_runtime, tmp_run, pes, **kw)
        screen = run(f"screen{i}", [species(tmp_run, pes, "r", "reactant")], **SCREEN)
        out = run(f"dft{i}", [*screen.artifacts, species(tmp_run, pes, "p", "product"), edge],
                  **DFT)
        assert dft.calls[:3] == ["energy"] * 3 and dft.calls.count("energy") == 3
        refined = {m.species_id for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"}
        (record,) = out.records(T.DISCOVERY, DiscoveryRecord)
        verdict = diagnostics(tmp_run, f"dft{i}")["edge"]
        if reason is None:
            assert refined == {"r", "p"} and record == found
            assert verdict[0] == "admitted" and float(verdict[1]) == pytest.approx(6.275, abs=0.01)
        else:
            assert refined == set() and (record.outcome, record.reason) == ("negative", reason)
            assert verdict[0] == reason


def test_a_relaxation_seed_is_asked_at_dft_once_from_its_own_geometry(fake_runtime, tmp_run):
    """R6 (S5, S6): the seed collapsed at screen and represents that basin, which DFT refines
    from the screen structure; the seed copy, the source of the edge without a TS into that
    basin, is asked once more, from its own geometry with the xTB start Hessian of any 2+
    fragment start (G4-P2: the engine writes it as its positive-definite model), after the
    species a screen minimum holds."""
    pes = fakes.double_well()  # the seed at the reactant (H at N); screen has no barrier to H-O
    seed = species(tmp_run, pes, "seed", "reactant")
    copy = species(tmp_run, pes, "seed_copy", "reactant")  # as explore records the lost seed
    init = {"init_hessian": {"engine": "xtb", "method": "gfn2"}}  # both wells: two fragments
    for high in (fakes.FakeQM, CollapseQM):
        run, _, dft = stage(fake_runtime, tmp_run, pes, low=CollapseQM, high=high)
        screen = run(f"screen_{high.__name__}", [seed], **SCREEN)
        (lost,) = screen.records(T.MINIMUM, MinimumRecord)
        assert lost.species_id == "seed" and lost.state_label != seed.payload.state_label
        relax = DiscoveryRecord(discovery_id="relax_seed", mechanism="relaxation",
                                outcome="product", source_species="seed_copy",
                                product_species="seed")
        found = Artifact(artifact_id="relax_seed", type=T.DISCOVERY, payload=relax)
        out = run(f"dft_{high.__name__}", [*screen.artifacts, copy, found], **DFT, **init)
        minima = {m.species_id: m for m in out.records(T.MINIMUM, MinimumRecord) if m.tier == "dft"}
        assert diagnostics(tmp_run, f"dft_{high.__name__}")["relax_seed"][0] == "admitted"
        assert dft.calls[:2] == ["energy", "energy"]  # the edge's two ends
        if high is fakes.FakeQM:  # PBE0 keeps the seed state: two basins, two freq jobs
            assert set(minima) == {"seed", "seed_copy"}
            assert minima["seed_copy"].state_label == seed.payload.state_label
            start = out.evidence(minima["seed_copy"].opt_calc).start
            assert start.fingerprint == seed.payload.geometry.fingerprint
            assert dft.calls[2:] == ["optimize+init_hessian"] * 2 + ["frequencies"] * 2
        else:  # the seed copy joins the collapse basin with no freq job
            assert minima["seed"].members == ("seed", "seed_copy")
            assert dft.calls[2:] == ["optimize+init_hessian"] * 2 + ["frequencies"]


@pytest.mark.parametrize("pes,points", [
    (fakes.harmonic(), ("start", "minimum", "start", "minimum", "start")),
    (fakes.double_well(), ("reactant", "product", "reactant", "product", "reactant")),
    (fakes.triple_well(), ("intermediate", "reactant", "product", "intermediate", "product")),
])
def test_parallel_minima_give_the_registry_of_one_by_one(fake_runtime, tmp_run, pes, points):
    """U9-P2(a): the opts at once and the freqs of new basins at once give the basins, members
    and representatives of relaxing and registering one job after the other."""
    rng = np.random.default_rng(7)
    starts = {f"s{i}": pes.points[p] + rng.normal(scale=0.02, size=(3, 3))
              for i, p in enumerate(points)}
    inputs = [species_at(tmp_run, sid, pes.symbols, x) for sid, x in starts.items()]
    load, qm = fakes.xyz_loader(tmp_run), fakes.FakeQM(tmp_run, pes)
    serial = Registry([], load)
    for art in inputs:
        mol = Molecule(XYZ(list(pes.symbols), starts[art.artifact_id]), 0, 1)
        out = relax_to_minimum(mol, METHODS["gfn2"], qm, known=serial, load_xyz=load,
                               resolve=lambda ref: tmp_run / ref.path)
        serial.add(out, art.payload, tier="screen")
    expected = {i: (m.species_id, m.members) for i, (m, _) in serial.minima.items()}
    for cores in (4, 1):
        run, xtb, _ = stage(fake_runtime, tmp_run, pes, cores=cores)
        out = run(f"screen{cores}", inputs, **SCREEN)
        assert {m.minimum_id: (m.species_id, m.members)
                for m in out.records(T.MINIMUM, MinimumRecord)} == expected
        assert xtb.calls.count("frequencies") == len(expected)  # a known basin runs none
