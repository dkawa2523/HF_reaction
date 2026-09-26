"""NWChem engines through JobRunner with a stub executable that replays recorded outputs:
G07 (HCN TS: level lines, vibrational block and .hess), G01 (two blocks) and G13 (string)."""

import shutil
import sys

import numpy as np
import pytest

from hfauto.backends.nwchem.engine import NWChemEngine, NWChemString
from hfauto.backends.nwchem.output import geometry_block
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.evidence import Evidence, FailureKind, PathProfile
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore

pytestmark = pytest.mark.golden
FINE = MethodSpec(id="fine", kind="dft", functional="pbe0", basis="def2-svpd",
                  dispersion="d3zero", grid="fine", scf_energy_tol=1e-7)  # G07's level
XFINE = FINE.model_copy(update={"grid": "xfine", "scf_energy_tol": 1e-8})  # G01 / G13
STUB = r'''
import re, shutil, sys
from pathlib import Path
here, deck = Path(__file__).parent, Path(sys.argv[-1]).read_text()
atoms = re.search(r"^geometry[^\n]*\n(.*?)^end", deck, re.M | re.S).group(1).splitlines()
task = deck.rsplit("task ", 1)[1].strip()
frame = lambda rows: f"{len(rows)}\n geometry\n" + "\n".join(rows) + "\n"
if task == "dft string":
    Path("job.string_final.xyz").write_text(frame(atoms) * 11)
    sys.exit(print((here / "G13.out").read_text()))
if len(atoms) == 4:
    sys.exit(print((here / "G01.out").read_text()))
g07 = (here / "G07.out").read_text()
echo = "".join(f"{i:5d} {r.split()[0]} 1.0 {' '.join(r.split()[1:])}\n"
               for i, r in enumerate(atoms, 1))
print("".join(g07.splitlines(True)[:2]) + ' Geometry "geometry" -> ""\n'
      + " Output coordinates in angstroms\n" + echo + "".join(g07.splitlines(True)[13:50])
      + " Total DFT energy =  -93.166994283613\n")
if task == "dft frequencies":
    shutil.copyfile(here / "G07.hess", "job.hess")
    print(g07[g07.index("  Vibrational analysis"):g07.index(" Task  times")])
if task == "dft optimize":
    restarted = "vectors input" in deck
    moved = [f"{r.split()[0]} " + " ".join(f"{float(v) + 0.01:.8f}" for v in r.split()[1:])
             for r in atoms]
    for n, rows in enumerate([atoms] * (3 if restarted else 1) + ([] if restarted else [moved])):
        Path(f"final-{n:03d}.xyz").write_text(frame(rows))
    Path("job.movecs").write_text("v"), Path("job.drv.hess").write_text("h")
    print("@    0    -93.16699428\n@    1    -93.16699428\n")
    if not restarted:
        sys.exit(124)
    print("      Optimization converged\n")
print(" Total times  cpu:        1.0s")
'''


@pytest.fixture
def nwchem(golden, tmp_path):
    """(JobRunner, EngineSite) whose executables['nwchem'] runs STUB with sys.executable."""
    stub = tmp_path / "stub"
    stub.mkdir()
    for name, rel in (("G01.out", "G01/nwchem.out"), ("G13.out", "G13/nwchem_string.out"),
                      ("G07.out", "G07/nwchem.out"), ("G07.hess", "G07/hfauto_job.hess")):
        shutil.copyfile(golden.path(f"nwchem/{rel}"), stub / name)
    (stub / "nwchem.py").write_text(STUB)
    if sys.platform == "win32":
        exe = stub / "nwchem.cmd"
        exe.write_text(f'@"{sys.executable}" "{stub / "nwchem.py"}" %*\n')
    else:
        exe = stub / "nwchem"
        exe.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{stub / "nwchem.py"}" "$@"\n')
        exe.chmod(0o755)
    site = EngineSite(version="7.2.3", executables={"nwchem": str(exe)})
    return JobRunner(JobStore(tmp_path / "run" / "jobs"), cores=1), site


def _hcn_ts(golden) -> Molecule:
    return Molecule(read_xyz(golden.path("nwchem/G07/final.xyz")), 0, 1)


def test_frequencies_cache_and_input_checks(nwchem, golden):
    jobs, site = nwchem
    engine, mol = NWChemEngine(jobs=jobs, site=site), _hcn_ts(golden)
    ev = engine.frequencies(mol, FINE)
    assert isinstance(ev, Evidence) and ev.task == "freq" and ev.n_external == 6
    assert ev.frequencies_cm1[0] == pytest.approx(-1131.6, abs=1.0)
    assert len(ev.imaginary_modes) == 1
    assert np.load(jobs.store.run_dir / ev.hessian.path).shape == (9, 9)
    assert engine.frequencies(mol, FINE) == ev and jobs.stats().hits == 1
    near, far = (Molecule(XYZ(mol.xyz.symbols, mol.xyz.coords + [d, 0, 0]), 0, 1)
                 for d in (0.1, 0.6))
    side = engine.optimize(near, FINE, init_hessian=ev)  # a QRC side from its TS Hessian
    assert isinstance(side, Evidence) and side.task == "opt"
    deck = (jobs.store.attempt_dir(side.job_key, 0) / "job.nw").read_text()
    assert "inhess 2" in deck and "trust 0.1" in deck
    failure = engine.optimize(far, FINE, init_hessian=ev)  # beyond 0.5 Å per atom
    assert failure.kind is FailureKind.INPUT_INVALID
    assert failure.reason == "hessian_geometry_mismatch"
    triplet = engine.frequencies(Molecule(mol.xyz, 0, 3), FINE)  # ported K case: state conflict
    assert triplet.kind is FailureKind.METHOD_MISMATCH and "multiplicity" in triplet.reason
    grid = engine.frequencies(mol, XFINE)
    assert grid.kind is FailureKind.METHOD_MISMATCH and "grid" in grid.reason


def test_open_shell_wft_is_rejected_without_a_job(nwchem, golden):
    jobs, site = nwchem
    mp2 = MethodSpec(id="mp2", kind="wft", wft_method="mp2", basis="def2-svp")
    doublet = Molecule(_hcn_ts(golden).xyz, 1, 2)
    failure = NWChemEngine(jobs=jobs, site=site).energy(doublet, mp2)
    assert (failure.kind, failure.reason) == (FailureKind.INPUT_INVALID, "wft_closed_shell_only")
    assert (jobs.stats().hits, jobs.stats().misses) == (0, 0)


def test_g01_two_vibrational_blocks_are_not_a_freq_evidence(nwchem, golden):
    symbols, coords = geometry_block(golden.text("nwchem/G01/nwchem.out"))
    mol = Molecule(XYZ(list(symbols), coords), 0, 1)
    failure = NWChemEngine(jobs=nwchem[0], site=nwchem[1]).frequencies(mol, XFINE)
    assert (failure.kind, failure.reason) == (FailureKind.INCOMPLETE_OUTPUT, "frequency_blocks:2")


def test_timeout_continues_from_the_latest_frame(nwchem, golden):
    jobs, site = nwchem
    mol = _hcn_ts(golden)
    ev = NWChemEngine(jobs=jobs, site=site).optimize(mol, FINE)
    assert isinstance(ev, Evidence) and ev.task == "opt"
    assert len(ev.trajectory_energies_hartree) == 4
    second = jobs.store.attempt_dir(ev.job_key, 1)
    assert (second / "job.drv.hess").is_file() and (second / "job.movecs").is_file()
    assert "vectors input job.movecs" in (second / "job.nw").read_text()
    final = read_xyz(jobs.store.run_dir / ev.final.file.path)
    assert np.allclose(final.coords, mol.xyz.coords + 0.01, atol=1e-7)
    assert read_xyz(jobs.store.run_dir / ev.start.file.path).coords == pytest.approx(mol.xyz.coords)


def test_g13_string_profile(nwchem, golden):
    text = golden.text("nwchem/G13/nwchem_string.out")
    (s0, c0), (s1, c1) = geometry_block(text, 0), geometry_block(text, 1)
    start, end = Molecule(XYZ(list(s0), c0), 0, 1), Molecule(XYZ(list(s1), c1), 0, 1)
    engine = NWChemString(jobs=nwchem[0], site=nwchem[1])
    profile = engine.find_path(start, end, XFINE, images=11)
    assert isinstance(profile, PathProfile) and len(profile.energies_hartree) == 11
    assert len(profile.energy_history) == 3 and profile.gmax_history == () and profile.ts is None
