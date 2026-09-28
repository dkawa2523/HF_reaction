"""G16 / G17 / G28: xTB termination and parsing on real outputs (§10.2, BUG-07, CH-20)."""

import shutil

import numpy as np
import pytest

from hfauto.backends.xtb import XTBEngine
from hfauto.chemistry.xyz import Molecule, geometry_fingerprint, read_xyz
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def parse(tmp_path, golden, kind, stdout, stderr=None, version="6.7.1", case="G17", charge=0):
    wd = tmp_path / "run" / "jobs" / "ab" / "abc" / "attempt_00"
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(golden.path(stdout), wd / "stdout.txt")
    (wd / "stderr.txt").write_text(golden.text(stderr) if stderr else "")
    mol = Molecule(read_xyz(shutil.copy(golden.path(f"xtb/{case}/xtbopt.xyz"), wd)), charge, 1)
    mol.write(wd / "input.xyz")
    xtb = XTBEngine(jobs=JobRunner(JobStore(tmp_path / "run" / "jobs"), cores=1),
                    site=EngineSite(version=version))
    done = CommandResult(0, False, 1.0, wd / "stdout.txt", wd / "stderr.txt")
    return xtb.adapter.parse(xtb.task(kind, mol, GFN2), wd, done), mol


def test_g16_rc0_with_failed_to_converge_is_geometry_maxiter(tmp_path, golden):
    failure, _ = parse(tmp_path, golden, "optimize", "xtb/G16/stdout.txt", "xtb/G16/stderr.txt")
    assert failure.kind == "geometry_maxiter"


def test_g17_converged_optimization_level_and_hessian(tmp_path, golden):
    opt, mol = parse(tmp_path, golden, "optimize", "xtb/G17/xtb.out")
    assert opt.task == "opt" and opt.energy_hartree == pytest.approx(-9.6707491583, abs=1e-10)
    level = opt.level
    assert (level.version, level.method, level.electronic_temperature_K) == ("6.7.1", "gfn2", 300)
    assert opt.final.fingerprint == geometry_fingerprint(mol.xyz.symbols, mol.xyz.coords)
    hessian = tmp_path / "run" / "jobs" / "ab" / "abc" / "attempt_00" / "hessian"
    hessian.write_text("$hessian\n" + " ".join(map(str, (0.3 * np.eye(18)).ravel())) + "\n$end")
    freq, _ = parse(tmp_path, golden, "frequencies", "xtb/G17/xtb.out")
    assert freq.task == "freq" and freq.n_external == 6 and len(freq.frequencies_cm1) == 12
    assert freq.start == freq.final and freq.hessian.path.endswith("attempt_00/hessian.npy")
    hessian.write_text("$hessian\n 0.1 ************** 0.2\n$end\n")  # overflowing field
    assert parse(tmp_path, golden, "frequencies", "xtb/G17/xtb.out")[0].kind == "incomplete_output"
    pinned, _ = parse(tmp_path, golden, "optimize", "xtb/G17/xtb.out", version="6.6.1")
    assert pinned.kind == "method_mismatch" and "version" in pinned.reason


def test_g28_chloride_atom_opt_and_hessian(tmp_path, golden):
    opt, _ = parse(tmp_path, golden, "optimize", "xtb/G28/opt_stdout.txt", case="G28", charge=-1)
    assert opt.task == "opt" and opt.energy_hartree == pytest.approx(-4.785133953761, abs=1e-10)
    assert (opt.level.charge, opt.level.multiplicity) == (-1, 1)
    shutil.copy(golden.path("xtb/G28/hessian"), tmp_path / "run" / "jobs" / "ab" / "abc" /
                "attempt_00")
    freq, _ = parse(tmp_path, golden, "frequencies", "xtb/G28/hess_stdout.txt", case="G28",
                    charge=-1)
    assert (freq.frequencies_cm1, freq.n_external, freq.imaginary_modes) == ((), 3, ())
