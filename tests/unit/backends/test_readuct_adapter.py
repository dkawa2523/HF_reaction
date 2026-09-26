"""ReaDuct adapter without SCINE (projected counts: G20, test_vibrations): helpers, job, result."""

import json

import numpy as np
from fakes import write_geometry

from hfauto.backends import engines
from hfauto.backends.protocols import Capability, DiscoveryEngine, DiscoverySettings
from hfauto.backends.readuct import worker
from hfauto.backends.readuct.engine import _Adapter
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.constants import HARTREE_TO_KJ_MOL as KJ
from hfauto.core.evidence import FailureKind
from hfauto.core.method import EngineSite, ExecutionSpec, MethodSpec
from hfauto.core.records import ReactionTrial
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

WATER = (["O", "H", "H"], np.array([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]))
GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def test_worker_helpers():
    symbols, x = WATER
    swapped = x[[0, 2, 1]] + 0.01
    opened = np.array([[0, 0, 0], [0.96, 0, 0], [-0.6, 0.77, 0]])  # same graph, RMSD > 0.1
    moved = np.array([[0, 0, 0], [0.96, 0, 0], [1.7, 0, 0]])  # H–H bond: another graph
    assert worker.matches_source(symbols, x, swapped)
    assert worker.matches_source(symbols, x, opened)  # a conformer is the source (state label)
    assert not worker.matches_source(symbols, x, moved)
    assert worker.irc_product(symbols, x, [swapped, moved]) == (1, None)
    assert worker.irc_product(symbols, x, [moved, opened]) == (0, None)
    assert worker.irc_product(symbols, x, [moved, moved]) == (None, "irc_not_connected_to_source")
    assert worker.irc_product(symbols, x, [x, swapped]) == (None, "same_as_source")
    # after an SCC retry at 1000 K the window energies are the base-temperature ones
    found = worker.Found("product", structures={"source": 0, "ts": 0, "product": 0},
                         energies={"source": -5.0, "ts": -4.9, "product": -4.95})
    data = worker.result_dict(found, {"source": -5.5, "ts": -5.44, "product": -5.47},
                              temperature_K=1000.0, version="6.1.0")
    assert data["barrier_kj_mol"] == (-5.44 + 5.5) * KJ
    assert data["reaction_kj_mol"] == (-5.47 + 5.5) * KJ
    assert data["electronic_temperature_K"] == 1000.0
    assert worker.is_scc_failure("scf: Self consistent charge iterator did not converge")


def test_explore_writes_the_worker_job_in_the_attempt_directory(tmp_run):
    jobs = JobRunner(JobStore(tmp_run / "jobs"), cores=1)
    site = EngineSite(version="6.1.0", python=str(tmp_run / "missing-python"))
    engine = engines.create(Capability.DISCOVERY, "readuct", jobs=jobs, site=site)
    assert isinstance(engine, DiscoveryEngine) and engine.supports(GFN2)
    trial = ReactionTrial(trial_id="t", source_minimum="m", kind="h_shift", mechanism="nt2",
                          associations=((0, 2),), dissociations=((0, 1),))
    mol = Molecule(XYZ(["H", "C", "N"], np.eye(3)), 0, 1)
    assert engine.explore(mol, trial, GFN2, DiscoverySettings()).kind == "executable_missing"
    job = json.loads(next((tmp_run / "jobs").glob("*/*/attempt_00/job.json")).read_text())
    assert (job["method_family"], job["trial"]["associations"]) == ("GFN2", [[0, 2]])


def test_result_mapping(tmp_run):
    workdir = tmp_run / "jobs" / "a"
    write_geometry(tmp_run, "jobs/a/ts.xyz", ["H", "C", "N"], np.eye(3))
    (workdir / "stderr.txt").write_text("Traceback: boom")
    adapter = _Adapter(EngineSite(version="6.1.0"), JobRunner(JobStore(tmp_run / "jobs"), cores=1))
    task = Task(engine="readuct", version_pin="6.1.0", kind="discovery", key_payload={},
                execution=ExecutionSpec())

    def parse(data, returncode=0):
        (workdir / "result.json").write_text(json.dumps(data))
        result = CommandResult(returncode, False, 1.0, workdir / "o", workdir / "stderr.txt")
        return adapter.parse(task, workdir, result)

    found = worker.Found("negative", "irc_not_connected_to_source", {"ts": 0}, ts_imag_cm1=-500.0)
    data = worker.result_dict(found, {}, temperature_K=300.0, version="6.1.0")
    result = parse(data)
    assert result.ts.file.path == "jobs/a/ts.xyz" and result.product is None
    assert parse({**data, "version": "6.0.0"}).kind == FailureKind.METHOD_MISMATCH
    failure = {"version": "6.1.0", "failure": {"kind": "scf_not_converged", "reason": "scc"}}
    assert parse(failure).kind == FailureKind.SCF_NOT_CONVERGED
    assert "boom" in parse(data, returncode=1).reason
