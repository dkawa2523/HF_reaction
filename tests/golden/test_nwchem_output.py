"""NWChem parser on real outputs (design §10.2: G01, G03-G05, G07, G10, G13, G21-G34)."""

import json
import re
import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from hfauto.backends.nwchem import input as nw_in
from hfauto.backends.nwchem import output as nw
from hfauto.backends.nwchem.engine import NWChemEngine, NWChemSaddle
from hfauto.chemistry.gates import same_pes, spin_ok
from hfauto.chemistry.vibrations import external_basis, projected_frequencies
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.evidence import FailureKind as Kind
from hfauto.core.method import EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import STDOUT_NAME, CommandResult
from hfauto.pipeline.config import load_method

pytestmark = pytest.mark.golden
METHODS = Path(__file__).resolve().parents[2] / "configs" / "methods"
XFINE = MethodSpec(id="m", kind="dft", functional="pbe0", basis="def2-svpd",
                   dispersion="d3zero", grid="xfine", scf_energy_tol=1e-8)
SVPD = MethodSpec(id="pbe0-d3bj_def2-svpd", kind="dft", functional="pbe0", basis="def2-svpd",
                  dispersion="d3bj", grid="fine", scf_energy_tol=1e-7)  # G25-G28
CCSD_T = MethodSpec(id="ccsd-t_def2-tzvpd", kind="wft", wft_method="ccsd(t)", basis="def2-tzvpd")


def test_g01_two_vibrational_blocks(golden):
    assert nw.count_frequency_blocks(golden.text("nwchem/G01/nwchem.out")) == 2


@pytest.mark.parametrize(("stem", "energy"), [
    ("G03/hono_trans", -205.3455299945), ("G03/hono_cis", -205.3452826520),
    ("G21/water_reference", -76.2956762656), ("G21/water_distorted", -76.2956762657)])
def test_energy_level_and_final_structure(golden, stem, energy):
    text = golden.text(f"nwchem/{stem}.out")
    assert level_mismatches(XFINE, nw.observe_level(text), version_pin="7.2.3") == []
    assert nw.total_energy(text) == pytest.approx(energy, abs=1e-9)
    final = read_xyz(golden.path(f"nwchem/{stem}_final.xyz"))
    symbols, coords = nw.geometry_block(text, -1)
    assert symbols == tuple(final.symbols) and np.allclose(coords, final.coords, atol=1e-8)


@pytest.mark.parametrize("stem", ["G04/hnc", "G05/hcn"])
def test_linear_molecules_have_five_external_modes(golden, stem):
    text = golden.text(f"nwchem/{stem}.out")
    symbols, coords = nw.geometry_block(text, -1)
    assert external_basis(symbols, coords).shape[1] == 5
    printed = [float(v) for row in re.findall(r"^\s*P\.Frequency(.*)$", text, re.MULTILINE)
               for v in row.split()]  # NWChem's own projection, one block
    assert len([f for f in printed if abs(f) > 1.0]) == 3 * 3 - 5


def test_g07_hess_to_canonical_npy(golden, tmp_path):
    xyz = read_xyz(golden.path("nwchem/G07/final.xyz"))
    hess = golden.path("nwchem/G07/hfauto_job.hess")
    npy = nw.hess_to_npy(hess, 3, tmp_path / "hessian.npy")
    freqs, _, n_external = projected_frequencies(np.load(npy), xyz.symbols, xyz.coords)
    assert n_external == 6 and freqs[0] == pytest.approx(-1131.6, abs=1.0)
    with pytest.raises(ValueError, match="3N"):  # ported K case: 3N mismatch
        nw.read_hess(hess, 4)


def test_g10_grid_mismatch_and_frame_change(golden):
    text = golden.text("nwchem/G10/irc_000.001.qce_nwchem_stdout")
    level = nw.observe_level(text)
    assert (level.grid, level.dispersion, level.scf_tol, level.basis) == (
        "medium", "d3zero", 1e-6, "def2-svpd")
    assert any(m.startswith("grid") for m in level_mismatches(XFINE, level, version_pin="7.2.3"))
    rows = [r.split() for r in text.split("geometry units bohr\n")[1].split("\nend")[0].split("\n")
            if r.strip()]
    deck = np.array([[float(v) for v in r[1:]] for r in rows]) * BOHR_TO_ANGSTROM
    assert nw.frame_shift(text, [r[0] for r in rows], deck) > 1e-4  # rigidly moved: frame lost
    assert nw.frame_shift(text, *nw.geometry_block(text, 0)) == 0.0


def test_g13_final_bead_energies(golden):
    text = golden.text("nwchem/G13/nwchem_string.out")  # the '@zts Bead number' lines
    energies = nw.string_energies(text)
    blocks = re.findall(r"Path Energy #.*\n((?: string:[ \t]+\d+[ \t]+\S+[ \t]*\n)+)", text)
    last = [float(row.split()[2]) for row in blocks[-1].splitlines()]  # the last iteration
    assert len(energies) == 11 and energies == pytest.approx(last, abs=1e-9)
    assert (energies[0] - energies[-1]) * HARTREE_TO_KCAL_MOL == pytest.approx(0.93, abs=0.01)


def test_g22_g23_g24_and_other_failures(golden):
    timeout = json.loads(golden.text("nwchem/G22/command_result.json"))
    killed = json.loads(golden.text("nwchem/G24/command_result.json"))
    failures = [
        nw.classify_failure(golden.text("nwchem/G22/nwchem.out"),
                            returncode=timeout["returncode"], timed_out=timeout["timed_out"]),
        nw.classify_failure(golden.text("nwchem/G23/nwchem.out"), returncode=1, timed_out=False),
        nw.classify_failure("", returncode=killed["returncode"], timed_out=killed["timed_out"]),
        nw.classify_failure("Calculation failed to converge", returncode=1, timed_out=False),
        nw.classify_failure("", returncode=0, timed_out=False),
    ]
    assert [f.kind for f in failures] == [
        Kind.TIMEOUT, Kind.GEOMETRY_MAXITER, Kind.NONZERO_EXIT, Kind.SCF_NOT_CONVERGED,
        Kind.INCOMPLETE_OUTPUT]


@pytest.mark.parametrize(("stem", "returncode", "failure"), [
    ("G31/hcn_opt.out", 0, None),
    ("G32/saddle_maxiter.out", 255, (Kind.GEOMETRY_MAXITER, "maxiter")),
    ("G33/opt_autoz.out", 236, (Kind.INPUT_INVALID, "autoz"))])
def test_the_autoz_notice_is_no_failure_only_a_fatal_autoz_error_is(golden, stem, returncode,
                                                                     failure):
    """G4-P4: G31 (HCN opt) and G32 (r6 s17's saddle, rerun from its start as an autoz
    failure) print the notice and go on in Cartesian coordinates; G33 (VAL7 s10's QRC side)
    stops at frame 21 with insufficient internal variables."""
    text = golden.text(f"nwchem/{stem}")
    notice = "AUTOZ failed to generate good internal coordinates.\n Cartesian coordinates will be"
    assert (notice in text) is (stem != "G33/opt_autoz.out")
    result = nw.classify_failure(text, returncode=returncode, timed_out=False)
    assert (result and (result.kind, result.reason)) == failure


def _parse(golden, tmp_path, engine_type, kind, stem, files, *, returncode=0, **inputs):
    """(engine, task, job directory, CommandResult) of a recorded SVPD job: its stdout and
    ``files`` in the job directory, its input echo and observed state as the molecule."""
    text, work = golden.text(f"nwchem/{stem}"), tmp_path / "jobs" / kind
    work.mkdir(parents=True)
    (work / STDOUT_NAME).write_text(text, encoding="utf-8")
    for name, rel in files.items():
        shutil.copyfile(golden.path(f"nwchem/{rel}"), work / name)
    symbols, coords = nw.geometry_block(text, 0)
    level = nw.observe_level(text)
    mol = Molecule(XYZ(list(symbols), coords), level.charge, level.multiplicity)
    site = EngineSite(version="7.2.3")
    engine = engine_type(jobs=JobRunner(JobStore(tmp_path / "jobs"), cores=1), site=site)
    task = Task(engine=engine.name, version_pin="7.2.3", kind=kind, key_payload={},
                execution=site.execution,
                inputs={"mol": mol, "start": mol, "method": SVPD, **inputs})
    result = CommandResult(returncode, False, 1.0, work / STDOUT_NAME, work / "stderr.txt")
    return engine, task, work, result


def test_g31_an_opt_carries_its_final_gradient(golden, tmp_path):
    """G5-P1: the last DFT ENERGY GRADIENTS block (bohr, Eh/bohr) is at the final frame (as in
    487 NWChem opt / saddle Evidence of VAL7 and r9); its largest component is the gmax of the
    '@ 1' line. A final frame elsewhere fails the job."""
    final = read_xyz(golden.path("nwchem/G31/hcn_opt_final.xyz"))
    coords, gradient = nw.gradient(golden.text("nwchem/G31/hcn_opt.out"))
    assert np.allclose(coords * BOHR_TO_ANGSTROM, final.coords, atol=1e-6)
    assert np.abs(gradient).max() == pytest.approx(3.0e-5)
    engine, task, work, result = _parse(golden, tmp_path, NWChemEngine, "optimize",
                                        "G31/hcn_opt.out",
                                        {"final-001.xyz": "G31/hcn_opt_final.xyz"})
    ev = engine.parse(task, work, result)
    assert isinstance(ev, Evidence) and ev.task == "opt" and ev.s2 is None
    assert ev.gradient == pytest.approx(gradient.ravel()) and len(ev.gradient) == 9
    moved = (work / "final-001.xyz").read_text().replace("0.74709120", "0.74809120")
    (work / "final-002.xyz").write_text(moved)  # 1e-3 Å along x: beyond FRAME_TOL_A
    failure = engine.parse(task, work, result)
    assert (failure.kind, failure.reason) == (Kind.INCOMPLETE_OUTPUT, "gradient_not_at_final")


def test_g32_a_saddle_at_maxiter_fails_with_its_last_frame_and_its_energy(golden, tmp_path):
    """G1-P3: the energy of the last '@' step line is the last frame's (as in 20 saddles of
    VAL7 and r9 at maxiter)."""
    engine, task, work, result = _parse(golden, tmp_path, NWChemSaddle, "saddle",
                                        "G32/saddle_maxiter.out",
                                        {"final-050.xyz": "G32/final-050.xyz"}, returncode=255)
    failure = engine.parse(task, work, result)
    assert isinstance(failure, Failure) and failure.kind is Kind.GEOMETRY_MAXITER
    assert failure.energy_hartree == pytest.approx(-200.24793703, abs=1e-9)
    assert failure.final.file.path.endswith("final-050.xyz")


def test_g34_an_scf_rescue_is_cgmin_then_one_plain_scf_that_prints_s2(golden, tmp_path):
    """G2-P6, the r6 P3c H3 bead through NWChemEngine: attempt_00 failed its SCF; attempt_01
    (this deck; X3 now starts cgmin from the atomic guess, within 40 cycles) converged with
    cgmin, then two plain SCF iterations 3.7e-8 Eh away printed <S2> 0.9861. The rescue is
    judged on that plain SCF, never on cgmin (P0a): the Evidence carries its <S2> (spin_ok
    rejects the contaminated doublet) when it reaches its parent's spin state; a plain SCF that
    left cgmin's solution or the parent's state, or a rescue that did not converge, closes as
    scf_unavailable."""
    deck, text = golden.text("nwchem/G34/h3_rescue.nw"), golden.text("nwchem/G34/h3_rescue.out")
    engine, task, work, result = _parse(golden, tmp_path, NWChemEngine, "energy",
                                        "G34/h3_rescue.out", {}, scf_rescue=True)
    setup = nw_in.Setup(scratch_dir=re.search(r"scratch_dir (\S+)", deck)[1], memory_mb=2000,
                        scf_rescue=True)
    atomic = deck.replace("  iterations 100\n", "  iterations 40\n").replace(
        "  vectors input job.movecs\n  cgmin\n", "  cgmin\n")
    assert nw_in.render_energy(task.inputs["mol"], SVPD, setup) == atomic != deck
    energies = nw.dft_energies(text)
    assert energies == pytest.approx((-1.548056627840, -1.548056664703), abs=1e-12)
    ev = engine.parse(task, work, result)
    assert isinstance(ev, Evidence) and ev.energy_hartree == energies[-1]
    assert ev.s2 == pytest.approx(0.9861) and spin_ok(ev).reasons == ("spin_contaminated",)
    for parent_s2, reason in ((0.95, None), (0.7539, "scf_unavailable:spin_state")):
        parent = ev.model_copy(update={"s2": parent_s2})
        judged = engine.parse(replace(task, inputs={**task.inputs, "parent": parent}), work, result)
        assert (judged == ev) if reason is None else (judged.kind, judged.reason) == (
            Kind.SCF_NOT_CONVERGED, reason)
    off = text.replace("-1.548056664703", "-1.548156664703")  # 1e-4 Eh below cgmin's
    (work / STDOUT_NAME).write_text(off, encoding="utf-8")
    failure = engine.parse(task, work, result)
    assert (failure.kind, failure.reason) == (Kind.SCF_NOT_CONVERGED, "scf_unavailable:scf_noise")
    (work / STDOUT_NAME).write_text(text + "\n Calculation failed to converge\n", encoding="utf-8")
    stopped = engine.parse(task, work, CommandResult(1, False, 1.0, work / STDOUT_NAME,
                                                     work / "stderr.txt"))
    assert (stopped.kind, stopped.reason) == (Kind.SCF_NOT_CONVERGED, "scf_unavailable:scf")


def test_final_xyz_numbering_atom_order_and_missing_d3(golden, tmp_path):
    frame = "3\n geometry\n{}\n{}\n{}\n"
    atoms = ("C 0 0 0", "N 0 0 1.16", "H 0 0 -1.07")
    (tmp_path / "final-1000.xyz").write_text(frame.format(*atoms))
    for n in (9, 999):  # written later, but numbered lower
        (tmp_path / f"final-{n:03d}.xyz").write_text(frame.format(*atoms))
    assert nw.final_xyz(tmp_path, ["C", "N", "H"]).name == "final-1000.xyz"
    assert nw.final_xyz(tmp_path, ["N", "C", "H"]) is None  # ported K case: atom order
    text = golden.text("nwchem/G03/hono_trans.out").replace("DFT-D3 Model", "")
    mismatches = level_mismatches(XFINE, nw.observe_level(text), version_pin="7.2.3")
    assert [m.split(":")[0] for m in mismatches] == ["dispersion"]  # ported K case: no D3


@pytest.mark.parametrize(("stem", "charge", "multiplicity", "energy"), [
    ("G25/oh_opt", 0, 2, -75.601066898809), ("G25/oh_freq", 0, 2, -75.601067860254),
    ("G26/fhf_opt", -1, 1, -200.036913760096), ("G26/fhf_freq", -1, 1, -200.036913730359),
    ("G26/fhf_sp", -1, 1, -200.036913904367), ("G27/hi_sp", 0, 1, -298.309133268697),
    ("G28/cl_opt", -1, 1, -459.982266730126), ("G28/cl_freq", -1, 1, -459.982266730122)])
def test_doublet_anion_ecp_and_atom_level_s2_and_energy(golden, stem, charge, multiplicity,
                                                        energy):
    text = golden.text(f"nwchem/{stem}.out")
    level = nw.observe_level(text)
    assert level_mismatches(SVPD, level, version_pin="7.2.3") == []
    assert (level.charge, level.multiplicity) == (charge, multiplicity)
    assert nw.classify_failure(text, returncode=0, timed_out=False) is None
    assert nw.total_energy(text) == pytest.approx(energy, abs=1e-9)
    s2 = nw.s2(text)
    assert s2 is None if multiplicity == 1 else s2 == pytest.approx(0.753, abs=0.005)


@pytest.mark.parametrize(("stem", "charge", "multiplicity", "n_external", "n_modes"), [
    ("G25/oh_freq", 0, 2, 5, 1), ("G26/fhf_freq", -1, 1, 5, 4), ("G28/cl_freq", -1, 1, 3, 0)])
def test_doublet_anion_and_atom_freq_evidence(golden, tmp_path, stem, charge, multiplicity,
                                              n_external, n_modes):
    """The engine's own parse of a recorded job: 3N - n_external modes, spin_ok on <S^2>."""
    text, work = golden.text(f"nwchem/{stem}.out"), tmp_path / "jobs" / "freq"
    work.mkdir(parents=True)
    (work / STDOUT_NAME).write_text(text, encoding="utf-8")
    shutil.copyfile(golden.path(f"nwchem/{stem}.hess"), work / "job.hess")
    symbols, coords = nw.geometry_block(text, 0)
    mol = Molecule(XYZ(list(symbols), coords), charge, multiplicity)
    site = EngineSite(version="7.2.3")
    engine = NWChemEngine(jobs=JobRunner(JobStore(tmp_path / "jobs"), cores=1), site=site)
    task = Task(engine="nwchem", version_pin="7.2.3", kind="frequencies", key_payload={},
                execution=site.execution, inputs={"mol": mol, "start": mol, "method": SVPD})
    result = CommandResult(0, False, 1.0, work / STDOUT_NAME, work / "stderr.txt")
    ev = engine.parse(task, work, result)
    assert isinstance(ev, Evidence) and ev.n_external == n_external
    assert len(ev.frequencies_cm1) == n_modes and all(nu > 0 for nu in ev.frequencies_cm1)
    assert spin_ok(ev) and (ev.s2 is None) == (multiplicity == 1)
    if multiplicity > 1:  # cgmin prints no <S2>: no UKS Evidence without it
        (work / STDOUT_NAME).write_text(re.sub("<S2> =.*", "", text), encoding="utf-8")
        failure = engine.parse(task, work, result)
        assert (failure.kind, failure.reason) == (Kind.INCOMPLETE_OUTPUT, "s2_not_reported")


def test_g27_iodine_ecp_electrons_and_ccsd_t_frozen_core(golden):
    """def2-ECP on I only: 26 electrons, alpha 13; ``freeze atomic`` freezes no orbital of an
    ECP atom (its 4s4p4d are correlated), so render_wft writes the count (G6-P3)."""
    sp, cc = golden.text("nwchem/G27/hi_sp.out"), golden.text("nwchem/G27/hi_ccsdt.out")
    for text in (sp, cc):
        assert "\necp\n  I library def2-ecp\nend\n" in text and "* library def2-ecp" not in text
        assert re.search(r"I \(Iodine\) Replaces\s+28 electrons", text)
    assert re.search(r"Alpha electrons :\s+13\n", sp)
    assert re.findall(r"number of core\s+(\d+)", cc) == ["0"]


def test_g35_a_ccsd_t_deck_freezes_the_core_left_by_the_ecp_and_favours_global_arrays(golden):
    """G6-P3, I-···CH3I (VAL7 s2's PBE0 minimum) through NWChemEngine, as render_wft writes it:
    NWChem freezes 9 orbitals (C 1s, 2 x I 4s4p) and splits memory_mb_per_rank 2000 into heap
    100, stack 500 and global 1400 MB; 1024 s on 4 ranks."""
    deck, text = golden.text("nwchem/G35/i_ch3i_ccsdt.nw"), golden.text("nwchem/G35/i_ch3i_ccsdt.out")
    symbols, coords = nw.geometry_block(text, 0)
    setup = nw_in.Setup(scratch_dir=re.search(r"scratch_dir (\S+)", deck)[1], memory_mb=2000)
    assert nw_in.render_wft(Molecule(XYZ(list(symbols), coords), -1, 1), CCSD_T, setup) == deck
    assert re.findall(r"number of core\s+(\d+)", text) == [str(nw_in.frozen_core(symbols))]
    memory = re.findall(r"(heap|stack|global)\s+=\s+\d+ doubles =\s+(\S+) Mbytes", text)
    assert memory == [("heap", "100.0"), ("stack", "500.0"), ("global", "1400.0")]
    level = nw.observe_level(text)
    assert level_mismatches(CCSD_T, level, version_pin="7.2.3") == []
    assert (level.charge, level.multiplicity) == (-1, 1)
    assert nw.total_energy(text) == pytest.approx(-634.213649871623, abs=1e-9)


@pytest.mark.parametrize(("stem", "multiplicity", "energy"), [
    ("G27/hi_ccsdt", 1, -297.821560434727), ("G29/oh_ccsdt", 2, -75.640597871999560)])
def test_ccsd_t_level_and_energy_rhf_ccsd_module_and_rohf_tce(golden, stem, multiplicity,
                                                              energy):
    """G27: the ccsd module (RHF); G29: the TCE on a ROHF reference (``open shells = 1``)."""
    text = golden.text(f"nwchem/{stem}.out")
    level = nw.observe_level(text)
    assert level_mismatches(CCSD_T, level, version_pin="7.2.3") == []
    assert (level.charge, level.multiplicity) == (0, multiplicity)
    assert nw.total_energy(text) == pytest.approx(energy, abs=1e-9) and nw.s2(text) is None


# Lines of /home/user/hfauto_r10/probe/P0d/ch3o_uhf.out (sha256 2ab14bd3c9c9...), CH3O.
# UHF-CCSD(T)/def2-TZVPD through the TCE, the deck render_wft writes (P0d).
_UHF_TCE = """\
             Northwest Computational Chemistry Package (NWChem) 7.2.3
  ao basis        = "ao basis"
  functions       =   104
  atoms           =     5
  alpha electrons =     9
  beta  electrons =     8
  charge          =   0.00
  wavefunction    = UHF

 Summary of "ao basis" -> "ao basis" (spherical)
 ------------------------------------------------------------------------------
       Tag                 Description            Shells   Functions and Types
 ---------------- ------------------------------  ------  ---------------------
 C                         def2-tzvpd               13       37   6s3p3d1f
 O                         def2-tzvpd               14       40   6s4p3d1f
 H                         def2-tzvpd                5        9   3s2p

 CCSD[T] total energy / hartree       =      -114.874770459222532
 CCSD(T)  correction energy / hartree =        -0.012710846433504
 CCSD(T) correlation energy / hartree =        -0.404282987575066
 CCSD(T) total energy / hartree       =      -114.873944099508378
"""


def test_uhf_ccsd_t_level_and_energy():
    """P0d: a UHF reference states no ``open shells``; alpha - beta electrons give the
    multiplicity. Its total agrees with PySCF's UHF-CCSD(T) within 2.1e-7 Eh."""
    level = nw.observe_level(_UHF_TCE)
    assert level_mismatches(CCSD_T, level, version_pin="7.2.3") == []
    assert (level.charge, level.multiplicity) == (0, 2) and nw.s2(_UHF_TCE) is None
    assert nw.total_energy(_UHF_TCE) == pytest.approx(-114.873944099508, abs=1e-12)


def test_g30_energy_layer_matches_its_method_file_with_m06_2x_d3_on_another_pes(golden):
    """The default energy layer as NWChem 7.2.3 ran it (S6 TS, doublet): the observed Level
    matches configs/methods, the D3 term is M06-2X's zero damping, and the layer is not the
    PBE0 freq's PES (G25: PBE0-D3BJ/def2-SVPD freq of the same charge and multiplicity)."""
    text = golden.text("nwchem/G30/m062x_ts.out")
    level = nw.observe_level(text)
    layer = load_method(METHODS / "m06-2x-d3_def2-tzvpd.yaml")
    assert level_mismatches(layer, level, version_pin="7.2.3") == []
    assert (level.method, level.dispersion, level.basis, level.grid, level.scf_tol) == (
        "m06-2x", "d3zero", "def2-tzvpd", "fine", 1e-7)
    assert re.search(r"s6 scale factor\s+:\s+1\.0+\s+s8 scale factor\s+:\s+0\.0+\s+"
                     r"sr6 scale factor\s+:\s+1\.619", text)
    assert nw.total_energy(text) == pytest.approx(-116.224696434258, abs=1e-9)
    assert (level.charge, level.multiplicity, nw.s2(text)) == (0, 2, pytest.approx(0.7587))
    freq = nw.observe_level(golden.text("nwchem/G25/oh_freq.out"))
    assert same_pes(level, freq).reasons == (
        "pes_mismatch:method", "pes_mismatch:basis", "pes_mismatch:dispersion")
