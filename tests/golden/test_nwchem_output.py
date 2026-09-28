"""NWChem parser on real outputs (design §10.2: G01, G03-G05, G07, G10, G13, G21-G29)."""

import json
import re
import shutil

import numpy as np
import pytest

from hfauto.backends.nwchem import output as nw
from hfauto.backends.nwchem.engine import NWChemEngine
from hfauto.chemistry.gates import spin_ok
from hfauto.chemistry.vibrations import external_basis, projected_frequencies
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence
from hfauto.core.evidence import FailureKind as Kind
from hfauto.core.method import EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import STDOUT_NAME, CommandResult

pytestmark = pytest.mark.golden
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
        nw.classify_failure(" AUTOZ failed to generate", returncode=1, timed_out=False),
        nw.classify_failure("Calculation failed to converge", returncode=1, timed_out=False),
        nw.classify_failure("", returncode=0, timed_out=False),
    ]
    assert [f.kind for f in failures] == [
        Kind.TIMEOUT, Kind.GEOMETRY_MAXITER, Kind.NONZERO_EXIT, Kind.INPUT_INVALID,
        Kind.SCF_NOT_CONVERGED, Kind.INCOMPLETE_OUTPUT]
    ok = " AUTOZ failed to generate good internal coordinates.\n" + nw.NORMAL_END
    assert nw.classify_failure(ok, returncode=0, timed_out=False) is None  # NWChem fell back


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
    if multiplicity > 1:  # the SCF rescue (cgmin) prints no <S2>: no UKS Evidence without it
        (work / STDOUT_NAME).write_text(re.sub("<S2> =.*", "", text), encoding="utf-8")
        failure = engine.parse(task, work, result)
        assert (failure.kind, failure.reason) == (Kind.INCOMPLETE_OUTPUT, "s2_not_reported")


def test_g27_iodine_ecp_electrons_and_ccsd_t_frozen_core(golden):
    """def2-ECP on I only: 26 electrons, alpha 13; ``freeze atomic`` freezes no orbital of an
    ECP atom (its 4s4p4d are correlated)."""
    sp, cc = golden.text("nwchem/G27/hi_sp.out"), golden.text("nwchem/G27/hi_ccsdt.out")
    for text in (sp, cc):
        assert "\necp\n  I library def2-ecp\nend\n" in text and "* library def2-ecp" not in text
        assert re.search(r"I \(Iodine\) Replaces\s+28 electrons", text)
    assert re.search(r"Alpha electrons :\s+13\n", sp)
    assert re.findall(r"number of core\s+(\d+)", cc) == ["0"]


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
