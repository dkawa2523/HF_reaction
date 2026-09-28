"""xTB adapter: command line, environment and continuation (§6.3)."""

import numpy as np

from hfauto.backends.xtb import XTBEngine
from hfauto.chemistry.xyz import XYZ, Molecule, write_xyz
from hfauto.core.evidence import Failure
from hfauto.core.method import EngineSite, ExecutionSpec, MethodSpec
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore

WATER = XYZ(["O", "H", "H"], np.array([[0.0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]))
GFN1 = MethodSpec(id="gfn1", kind="xtb", gfn=1)


def engine(tmp_run) -> XTBEngine:
    site = EngineSite(version="6.7.1", execution=ExecutionSpec(threads=2))
    return XTBEngine(jobs=JobRunner(JobStore(tmp_run / "jobs"), cores=2), site=site)


def test_command_carries_charge_uhf_and_etemp(tmp_run):
    xtb = engine(tmp_run)
    method = MethodSpec(id="g", kind="xtb", gfn=2, electronic_temperature_K=1000.0)
    cmd = xtb.adapter.prepare(xtb.task("optimize", Molecule(WATER, -1, 2), method), tmp_run / "a")
    assert cmd.argv[1:] == ("input.xyz", "--opt", "vtight", "--gfn", "2", "--chrg", "-1", "--uhf",
                            "1", "--etemp", "1000.0")
    assert cmd.env["OMP_NUM_THREADS"] == "2,1" and cmd.env["OMP_STACKSIZE"] == "4G"
    hess = xtb.adapter.prepare(xtb.task("frequencies", Molecule(WATER, 0, 1), GFN1), tmp_run / "b")
    assert hess.argv[2:] == ("--hess", "--gfn", "1", "--chrg", "0", "--uhf", "0")
    assert xtb.supports(method) and not xtb.supports(MethodSpec(id="d", kind="dft"))


def test_optimization_continues_from_xtbopt(tmp_run):
    xtb, work = engine(tmp_run), tmp_run / "attempt_00"
    task = xtb.task("optimize", Molecule(WATER, 0, 1), GFN1)
    xtb.adapter.prepare(task, work)
    write_xyz(XYZ(WATER.symbols, WATER.coords * 1.01), work / "xtbopt.xyz")
    new = xtb.adapter.continuation(task, work, Failure(kind="geometry_maxiter", reason="")).inputs
    assert np.allclose(new["mol"].xyz.coords, WATER.coords * 1.01)
    assert new["start"].file.path == "attempt_00/input.xyz"
