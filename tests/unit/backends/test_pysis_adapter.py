"""pysis_neb: pysisyphus input, worker command / PATH and result parsing (§6.3)."""

import json
import os
import sys
from types import ModuleType
from types import SimpleNamespace as NS

import numpy as np
import pytest

from hfauto.backends.pysis.engine import NEB_MAX_CYCLES, PysisNEB, neb_input
from hfauto.backends.pysis.worker import run_neb
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz_trajectory, write_xyz_trajectory
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import FailureKind, FileRef
from hfauto.core.hashing import sha256_file
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def mol(symbols, *rows, charge=0):
    return Molecule(XYZ(list(symbols), np.array(rows, dtype=float)), charge, 1 + abs(charge))


HCN = mol("CNH", [0, 0, 0], [0, 0, 1.156], [0, 0, -1.066])
HNC = mol("CNH", [0, 0, 0], [0, 0, 1.17], [0, 0, 2.17])
MID = mol("CNH", [0, 0, 0], [0, 0, 1.16], [1.0, 0, 0.5])
HF_PAIR = mol("FHFH", [0, 0, 0], [0.92, 0, 0], [3.0, 0.4, 0], [3.9, 0.6, 0], charge=-1)


def trajectory(path, *mols):
    return write_xyz_trajectory([m.xyz for m in mols], path)


def test_one_cartesian_ci_neb_with_fixed_ends_then_tsopt():
    run = neb_input(HCN, GFN2, threads=4)
    assert run["geom"] == {"type": "cart", "fn": "initial.trj"} and "interpol" not in run
    assert run["cos"] == {"type": "neb", "climb": True}  # fix_first / fix_last by default
    assert run["opt"] == {"type": "lbfgs", "max_cycles": NEB_MAX_CYCLES}
    assert run["tsopt"] == {"type": "rsprfo"}
    assert run["calc"] == {"type": "xtb", "gfn": 2, "charge": 0, "mult": 1, "pal": 4}
    method = MethodSpec(id="s", kind="xtb", gfn=2, electronic_temperature_K=500.0)
    pair = neb_input(HF_PAIR, method, threads=1)
    assert pair["geom"]["type"] == "cart" and "geom" not in pair["tsopt"]
    assert [pair["calc"][k] for k in ("charge", "mult", "etemp")] == [-1, 2, 500]


def engine(tmp_run):  # any existing executable stands in for the site's xtb
    site = EngineSite(version="6.7.1", executables={"xtb": sys.executable}, python="worker-python")
    return PysisNEB(jobs=JobRunner(JobStore(tmp_run / "jobs"), cores=1), site=site)


def initial_path(tmp_run) -> FileRef:
    path = trajectory(tmp_run / "idpp.xyz", HCN, MID, HNC)
    return FileRef(path="idpp.xyz", sha256=sha256_file(path))


def test_worker_command_path_and_the_initial_path(tmp_run):
    neb, ref = engine(tmp_run), initial_path(tmp_run)
    task = neb.task(HCN, HNC, GFN2, images=3, initial_path=ref)
    (tmp_run / "a").mkdir()  # the attempt directory, made by JobRunner
    cmd = neb.adapter.prepare(task, tmp_run / "a")
    assert cmd.argv[:4] == ("worker-python", "-m", "hfauto.execution.worker",
                            "hfauto.backends.pysis.worker:run_neb")
    assert cmd.env["PATH"].split(os.pathsep)[0] == os.path.dirname(sys.executable)
    assert sha256_file(tmp_run / "a" / "initial.trj") == ref.sha256  # the IDPP, unchanged
    job = json.loads((tmp_run / "a" / "job.json").read_text())
    assert job["run_dict"] == neb_input(HCN, GFN2, threads=1)
    other = FileRef(path=ref.path, sha256="0" * 64)
    assert neb.task(HCN, HNC, GFN2, images=3, initial_path=other).key_payload != task.key_payload


def test_neb_images_become_a_path_profile(tmp_run):
    work, neb = tmp_run / "attempt_00", engine(tmp_run)
    (work / "qm_calcs").mkdir(parents=True)
    (work / "qm_calcs" / "calculator_000.000.xtb.out").write_text("   * xtb version 6.7.1 (x)")
    task = neb.task(HCN, HNC, GFN2, images=3, initial_path=initial_path(tmp_run))
    data = {"ts": MID.xyz.coords.tolist(), "xtb_version": "6.7.1"}
    (work / "result.json").write_text(json.dumps(data))
    (work / "stderr.txt").write_text("Traceback\nRuntimeError: boom\n")
    ok, crashed = (CommandResult(rc, False, 1.0, work / "stdout.txt", work / "stderr.txt")
                   for rc in (0, 1))
    assert neb.adapter.parse(task, work, ok).reason == "neb_images"  # no images at all
    trajectory(work / "current_geometries.trj", HCN, MID, HNC)  # an error stopped the NEB
    assert neb.adapter.parse(task, work, ok).images.path.endswith("current_geometries.trj")
    trajectory(work / "final_geometries.trj", HCN, MID, HNC)
    profile = neb.adapter.parse(task, work, ok)
    assert profile.images.path.endswith("final_geometries.trj") and profile.energies_hartree == ()
    assert np.allclose(read_xyz_trajectory(tmp_run / profile.images.path)[1].coords,
                       MID.xyz.coords)
    assert (tmp_run / profile.ts.file.path).is_file() and profile.level.method == "gfn2"
    trajectory(work / "final_geometries.trj", HCN, HNC)
    assert neb.adapter.parse(task, work, ok).kind is FailureKind.INCOMPLETE_OUTPUT
    assert neb.adapter.parse(task, work, crashed).reason.endswith("RuntimeError: boom")


def fake_pysisyphus(monkeypatch, run_from_dict):
    module = ModuleType("pysisyphus.run")
    module.run_from_dict = run_from_dict
    monkeypatch.setitem(sys.modules, "pysisyphus.run", module)


def test_worker_keeps_the_images_of_a_failed_tsopt_and_raises_before_any(tmp_run, monkeypatch):
    ts = NS(cart_coords=np.ravel(MID.xyz.coords) / BOHR_TO_ANGSTROM)

    def converged(run_dict, cwd):
        return NS(ts_geom=ts, ts_opt=NS(is_converged=True))

    fake_pysisyphus(monkeypatch, converged)
    assert np.allclose(run_neb({"run_dict": {}}, tmp_run)["ts"], MID.xyz.coords)
    fake_pysisyphus(monkeypatch, lambda run_dict, cwd: NS(ts_geom=ts, ts_opt=NS(is_converged=False)))
    assert run_neb({"run_dict": {}}, tmp_run)["ts"] is None

    def neb_then_tsopt_raises(run_dict, cwd):
        trajectory(cwd / "final_geometries.trj", HCN, MID, HNC)
        raise RuntimeError("TS optimization diverged")

    def raises_before_a_cycle(run_dict, cwd):
        raise ValueError("xtb failed")

    fake_pysisyphus(monkeypatch, raises_before_a_cycle)
    with pytest.raises(ValueError, match="xtb failed"):
        run_neb({"run_dict": {}}, tmp_run)
    fake_pysisyphus(monkeypatch, neb_then_tsopt_raises)
    data = run_neb({"run_dict": {}}, tmp_run)  # once, no rerun without the TS optimization
    assert data["ts"] is None and data["error"] == "RuntimeError: TS optimization diverged"
