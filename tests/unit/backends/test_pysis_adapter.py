"""pysis_gs: pysisyphus input, worker command / PATH and result parsing (§6.3)."""

import json
import os
import sys
from types import ModuleType
from types import SimpleNamespace as NS

import numpy as np
import pytest

from hfauto.backends.pysis.engine import PysisGrowingString, gs_input
from hfauto.backends.pysis.worker import run_growing_string, summarize
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def mol(symbols, *rows, charge=0):
    return Molecule(XYZ(list(symbols), np.array(rows, dtype=float)), charge, 1 + abs(charge))


HCN = mol("CNH", [0, 0, 0], [0, 0, 1.156], [0, 0, -1.066])
HNC = mol("CNH", [0, 0, 0], [0, 0, 1.17], [0, 0, 2.17])
BENT = mol("CNH", [0, 0, 0], [0, 0, 1.156], [0.185, 0, -1.05])  # 10° off linear: Cartesian
WATER = mol("OHH", [0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0])
HF_PAIR = mol("FHFH", [0, 0, 0], [0.92, 0, 0], [3.0, 0.4, 0], [3.9, 0.6, 0], charge=-1)


def test_coordinates_ts_optimization_and_calculator():
    linear = gs_input(BENT, HNC, GFN2, images=11, refine_ts=False, threads=4)
    assert linear["geom"]["type"] == "cart" and "tsopt" not in linear and "interpol" not in linear
    assert linear["cos"] == {"type": "gs", "max_nodes": 9, "climb": True, "climb_rms": 5e-3}
    assert linear["calc"] == {"type": "xtb", "gfn": 2, "charge": 0, "mult": 1, "pal": 4}
    bent = gs_input(WATER, WATER, GFN2, images=11, refine_ts=True, threads=1)
    assert bent["geom"]["type"] == "dlc" and bent["tsopt"] == {"type": "rsirfo", "thresh": "gau"}
    method = MethodSpec(id="s", kind="xtb", gfn=2, solvation="alpb:water",
                        electronic_temperature_K=500.0)
    pair = gs_input(HF_PAIR, HF_PAIR, method, images=7, refine_ts=True, threads=1)
    assert pair["geom"]["type"] == "cart" and pair["tsopt"]["geom"] == {"type": "tric"}
    assert [pair["calc"][k] for k in ("charge", "mult", "alpb", "etemp")] == [-1, 2, "water", 500]


def engine(tmp_run):  # any existing executable stands in for the site's xtb
    site = EngineSite(version="6.7.1", executables={"xtb": sys.executable}, python="worker-python")
    return PysisGrowingString(jobs=JobRunner(JobStore(tmp_run / "jobs"), cores=1), site=site)


def test_worker_command_path_and_linear_ends_moved_off_axis(tmp_run):
    gs = engine(tmp_run)
    cmd = gs.adapter.prepare(gs.task(HCN, HNC, GFN2, images=11, refine_ts=True), tmp_run / "a")
    assert cmd.argv[:4] == ("worker-python", "-m", "hfauto.execution.worker",
                            "hfauto.backends.pysis.worker:run_growing_string")
    assert cmd.env["PATH"].split(os.pathsep)[0] == os.path.dirname(sys.executable)
    moved = read_xyz(tmp_run / "a" / "start.xyz").coords  # ±0.01 Å across the axis
    assert np.allclose(moved - HCN.xyz.coords, [[0.01, 0, 0], [-0.01, 0, 0], [0.01, 0, 0]])


def image(coords, energy):  # a pysisyphus Geometry: Cartesian bohr and energy
    return NS(cart_coords=np.ravel(coords) / BOHR_TO_ANGSTROM, energy=energy)


def test_run_result_becomes_a_path_profile(tmp_run):
    work, mid = tmp_run / "attempt_00", (HCN.xyz.coords + BENT.xyz.coords) / 2
    (work / "qm_calcs").mkdir(parents=True)
    (work / "qm_calcs" / "calculator_000.000.xtb.out").write_text("   * xtb version 6.7.1 (x)")
    cos = NS(images=[image(HCN.xyz.coords, -5.50), image(mid, -5.38),
                     image(BENT.xyz.coords, -5.47)], started_climbing=True, get_hei_index=lambda: 1)
    result = NS(cos=cos, cos_opt=NS(max_forces=[1e-2, 1e-3], is_converged=True),
                ts_geom=image(mid, -5.387), ts_opt=NS(is_converged=True))
    data = summarize(result, work)
    assert np.allclose(data["images"][2], BENT.xyz.coords) and data["xtb_version"] == "6.7.1"
    assert summarize(NS(**{**vars(result), "ts_opt": NS(is_converged=False)}), work)["ts"] is None
    (work / "result.json").write_text(json.dumps(data))
    (work / "stderr.txt").write_text("Traceback\nRuntimeError: boom\n")
    ok, crashed = (CommandResult(rc, False, None, 1.0, work / "stdout.txt", work / "stderr.txt")
                   for rc in (0, 1))
    gs = engine(tmp_run)
    task = gs.task(HCN, BENT, GFN2, images=3, refine_ts=True)
    profile = gs.adapter.parse(task, work, ok)
    assert profile.energies_hartree == (-5.50, -5.38, -5.47) and profile.climbing_image == 1
    assert profile.ts_energy_hartree == -5.387 and (tmp_run / profile.ts.file.path).is_file()
    assert profile.program_converged and profile.gmax_history == (1e-2, 1e-3)
    assert gs.adapter.parse(task, work, crashed).reason.endswith("RuntimeError: boom")


def fake_pysisyphus(monkeypatch, run_from_dict):
    module = ModuleType("pysisyphus.run")
    module.run_from_dict = run_from_dict
    monkeypatch.setitem(sys.modules, "pysisyphus.run", module)


def test_tsopt_failure_keeps_the_growing_string_and_gs_failures_propagate(tmp_run, monkeypatch):
    gs_only = NS(cos=NS(images=[image(HCN.xyz.coords, -5.50), image(HNC.xyz.coords, -5.47)]),
                 cos_opt=NS(max_forces=[1e-3], is_converged=True), ts_geom=None, ts_opt=None)
    calls = []

    def tsopt_raises(run_dict, cwd):
        calls.append(sorted(run_dict))
        if "tsopt" in run_dict:
            raise RuntimeError("TS optimization diverged")
        return gs_only

    job = {"run_dict": {"cos": {"type": "gs"}, "tsopt": {"type": "rsirfo"}}}
    fake_pysisyphus(monkeypatch, tsopt_raises)
    data = run_growing_string(job, tmp_run)
    assert calls == [["cos", "tsopt"], ["cos"]] and data["ts"] is None
    assert data["tsopt_error"] == "RuntimeError: TS optimization diverged"
    assert len(data["images"]) == 2 and data["converged"]

    def always_raises(run_dict, cwd):
        raise ValueError("growing string failed")

    fake_pysisyphus(monkeypatch, always_raises)
    for failing in (job, {"run_dict": {"cos": {"type": "gs"}}}):
        with pytest.raises(ValueError, match="growing string failed"):
            run_growing_string(failing, tmp_run)
