"""GoodVibes adapter: QCData assembly in the worker, the engine's job and its cache key."""

import importlib.util
import json

import pytest
from fakes import FakeQM, harmonic

from hfauto.backends.goodvibes.engine import GoodVibesEngine
from hfauto.backends.goodvibes.worker import compute, qcdata_fields
from hfauto.chemistry.gates import zpe_hartree
from hfauto.chemistry.thermo import thermo_frequencies
from hfauto.core.method import EngineSite, MethodSpec, ThermoSettings
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore

GOODVIBES = importlib.util.find_spec("goodvibes") is not None
HNC = {  # G04 (NWChem, linear HNC): B = 2.170251 K, mol. weight 27.0109
    "version_pin": "4.3.0", "symbols": ["C", "N", "H"], "energy_hartree": -93.220271377032,
    "coords": [[0.74716204, 0.24207168, -0.79124657], [-0.03790437, -0.01228416, 0.04015256],
               [-0.70925766, -0.22978752, 0.75109401]], "charge": 0, "multiplicity": 1,
    "frequencies_cm1": [-8.0, 466.36, 2150.62, 3840.75], "n_external": 5,
    "temperatures_K": [298.15], "settings": [{"sha": "s", "qs": "grimme", "cutoff_cm1": 100.0,
                                              "vib_scale": 0.985, "symmetry": True}]}


def test_linear_molecule_qcdata_and_worker(tmp_path):  # port of test_goodvibes_nwchem
    f = qcdata_fields(HNC)
    assert f["linear_mol"] and f["rotemp"] == pytest.approx([2.17025] * 2, abs=1e-4)
    assert (f["frequency_wn"], f["im_frequency_wn"]) == ([466.36, 2150.62, 3840.75], [-8.0])
    assert f["molecular_mass"] == pytest.approx(27.0109, abs=1e-4)
    if GOODVIBES:  # the WSL environment; vib_scale scales the modes and the ZPE alike
        sent = list(thermo_frequencies(HNC["frequencies_cm1"], saddle=False))
        [r] = compute(HNC | {"frequencies_cm1": sent}, tmp_path)["results"]
        assert r["S_rot"] > 0 and r["n_real"] == 4
        assert r["zpe_hartree"] == pytest.approx(zpe_hartree(sent, scale=0.985), abs=1e-8)


def test_engine_job_key_and_fail_closed_worker(tmp_run):  # the key follows the modes sent
    method = MethodSpec(id="m", kind="dft", functional="pbe0", basis="def2-svpd")
    freq = FakeQM(tmp_run, pes := harmonic()).frequencies(pes.molecule("minimum"), method)
    store = JobStore(tmp_run / "jobs")
    engine = GoodVibesEngine(jobs=JobRunner(store, cores=1), site=EngineSite(version="4.3.0"))
    for saddle in (False, True):  # without goodvibes: a Failure, no fallback
        out = engine.thermo(freq, [ThermoSettings()], temperatures_K=(298.15,), saddle=saddle)
        assert out[0].job_key and out[0].S_rot > 0 if GOODVIBES else out.kind == "nonzero_exit"
    payloads = [json.loads(p.read_text())["key_payload"] for p in store.root.glob("*/*/job.json")]
    keys, n = {len(k["frequencies_cm1"]): k for k in payloads}, len(freq.frequencies_cm1)
    assert sorted(keys) == [n - 1, n]  # a saddle drops its lowest mode: another key
    key = keys[n]
    assert (key["freq"], key["hessian_sha"]) == (freq.job_key, freq.hessian.sha256)
    assert key["frequencies_cm1"] == sorted(abs(nu) for nu in freq.frequencies_cm1)
    assert key["temperatures_K"] == [298.15] and key["settings"][0]["vib_scale"] == 1.0
    assert set(key["settings"][0]) == {"sha", "qs", "cutoff_cm1", "vib_scale", "symmetry",
                                       "sensitivity"}
