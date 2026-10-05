"""ReaDuct adapter without SCINE (projected counts: G20, test_vibrations): helpers, job, result,
and the minimum optimization on a scripted surface."""

import json

import numpy as np
import pytest
from fakes import NH3_HF, NH3_HF_EXCHANGED, NH3_HF_SYMBOLS, write_geometry

from hfauto.backends import engines
from hfauto.backends.protocols import Capability, DiscoveryEngine, DiscoverySettings
from hfauto.backends.readuct import worker
from hfauto.backends.readuct.engine import _Adapter
from hfauto.chemistry.modes import amplitude
from hfauto.chemistry.topology import bond_changes, same_bonding, state_label
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as KCAL
from hfauto.core.evidence import FailureKind
from hfauto.core.method import EngineSite, ExecutionSpec, MethodSpec
from hfauto.core.records import ReactionTrial
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

WATER = (["O", "H", "H"], np.array([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]))
GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def test_worker_helpers(monkeypatch):
    symbols, x = WATER
    swapped = x[[0, 2, 1]] + 0.01  # the same bonds atom by atom
    opened = np.array([[0, 0, 0], [0.96, 0, 0], [-0.6, 0.77, 0]])  # same graph, RMSD > 0.1
    moved = np.array([[0, 0, 0], [0.96, 0, 0], [1.7, 0, 0]])  # H–H bond: another graph
    cut = np.array([[0, 0, 0], [0.96, 0, 0], [-3.0, 0, 0]])  # O–H cut: a third graph
    assert worker.irc_ends(symbols, x, [swapped, moved]) == ((0, 1), True)
    assert worker.irc_ends(symbols, x, [moved, opened]) == ((1, 0), True)  # source end first
    assert worker.irc_ends(symbols, x, [moved, cut]) == ((0, 1), False)  # two other states
    assert worker.irc_ends(symbols, x, [x, swapped]) is None  # no bond change: no edge
    assert worker.irc_ends(symbols, x, [moved, moved + 0.01]) is None
    # energies from the first end; after an SCC retry at 1000 K they are the base-temperature ones
    found = worker.Found("product", structures={"source": 0, "ts": 0, "end0": 0, "end1": 0})
    data = worker.result_dict(found, {"source": -5.5, "ts": -5.44, "end0": -5.49, "end1": -5.47},
                              version="v")
    assert data["ends"] == list(worker.END_FILES) and data["ts"] == worker.TS_FILE
    assert data["dE_act_kcal"] == (-5.44 + 5.49) * KCAL
    assert data["dE_rxn_kcal"] == (-5.47 + 5.49) * KCAL
    negative = worker.Found("negative", "no_bond_change", {"source": 0, "ts": 0},
                            {"source": -5.0, "ts": -4.9})  # a negative with a TS keeps its barrier
    data = worker.result_dict(negative, negative.energies, version="v")
    assert data["dE_act_kcal"] == pytest.approx(0.1 * KCAL) and data["dE_rxn_kcal"] is None
    assert data["ends"] is None
    relaxed = worker.Found("product", structures={"end0": 0, "end1": 0})  # no TS, no energies
    data = worker.result_dict(relaxed, {}, version="v")
    assert (data["ts"], data["dE_act_kcal"], data["dE_rxn_kcal"]) == (None, None, None)
    assert worker.is_scc_failure("scf: Self consistent charge iterator did not converge")
    monkeypatch.setattr(worker.importlib.metadata, "version", lambda name: "1.0")
    assert worker.version() == (
        "scine-readuct==1.0,scine-utilities==1.0,scine-xtb-wrapper==1.0")


def test_a_degenerate_rearrangement_is_an_edge():
    """U4-P5: the amine_pilot2 double H exchange (−1265i) has the source's state label but other
    bonds atom by atom: an edge from the source. An end that is the source relabelled (R5a: S19,
    S6) is another structure too: the ends stay in IRC order for the caller to place as states."""
    symbols, x, exchanged = NH3_HF_SYMBOLS, NH3_HF, NH3_HF_EXCHANGED
    assert state_label(symbols, exchanged) == state_label(symbols, x)
    assert bond_changes(symbols, x, exchanged) == ({(0, 4), (3, 5)}, {(0, 3), (4, 5)})
    assert worker.irc_ends(symbols, x, [x + 0.01, exchanged]) == ((0, 1), True)
    product = x.copy()
    product[4] = [0.31, -0.98, 0.0]  # the HF proton on N: NH4+ ... F-
    relabelled = exchanged + [1.0, -2.0, 0.5]
    assert not same_bonding(symbols, x, relabelled)
    assert worker.irc_ends(symbols, x, [product, relabelled]) == ((0, 1), False)


SADDLE = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])
MODE = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])


class _Surface(worker._Scine):
    """_Scine without SCINE: an optimization from SADDLE stops there (ν −170 cm⁻¹ along MODE);
    one from a displaced start reaches the minimum 0.5 Å further on its side (energy = x)."""

    def __init__(self, start: np.ndarray) -> None:
        self.symbols, self.job = ["H", "H"], {"settings": {"imag_cutoff_cm1": 50.0}}
        self.systems, self.calls, self.loaded = {"end": start}, [], []

    def load(self, name, coords):
        self.loaded.append(coords)
        self.systems[name] = coords

    def task(self, run, name, output=(), *, strict=True, **settings):
        self.calls.append((run, settings))
        side = np.sign(self.systems[name][1, 0] - SADDLE[1, 0])
        self.systems[output[0]] = SADDLE + 0.5 * side * MODE.reshape(2, 3)
        return True

    def imaginary(self, name):
        stuck = np.array_equal(self.systems[name], SADDLE)
        return (1, -170.0, MODE) if stuck else (0, 100.0, MODE)

    def coords(self, name):
        return self.systems[name]

    def energy(self, name):
        return float(self.systems[name][1, 0])


def test_a_minimum_left_on_a_saddle_is_displaced_along_its_lowest_mode():
    """U4-P4: an end that stops with ν < −cutoff is displaced ± by modes.amplitude, each side
    optimized once and the lower minimum taken; every optimization may take 500 iterations."""
    run = _Surface(SADDLE)
    assert run.minimum("end", "end_opt")
    step = amplitude(-170.0, MODE, run.symbols)
    assert [x[1, 0] - SADDLE[1, 0] for x in run.loaded] == pytest.approx([step, -step])
    assert run.coords("end_opt")[1, 0] == pytest.approx(SADDLE[1, 0] - 0.5)
    assert [s for _, s in run.calls] == [worker.MAX_ITERATIONS] * 3
    clean = _Surface(SADDLE + 0.1 * MODE.reshape(2, 3))  # already a minimum: no displacement
    assert clean.minimum("end", "end_opt") and not clean.loaded


def test_explore_writes_the_worker_job_in_the_attempt_directory(tmp_run):
    jobs = JobRunner(JobStore(tmp_run / "jobs"), cores=1)
    site = EngineSite(version="6.1.0", python=str(tmp_run / "missing-python"))
    engine = engines.create(Capability.DISCOVERY, "readuct", jobs=jobs, site=site)
    assert isinstance(engine, DiscoveryEngine) and engine.supports(GFN2)
    trial = ReactionTrial(trial_id="t", kind="f1b1", associations=((0, 2),),
                          dissociations=((0, 1),))
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

    found = worker.Found("negative", "no_bond_change", {"ts": 0})
    data = worker.result_dict(found, {}, version="6.1.0")
    result = parse(data)
    assert result.ts.file.path == "jobs/a/ts.xyz" and result.ends is None
    for name in worker.END_FILES:
        write_geometry(tmp_run, f"jobs/a/{name}", ["H", "C", "N"], np.eye(3))
    edge = worker.Found("product", None, {"ts": 0, "end0": 0, "end1": 0}, irc_connected=True)
    result = parse(worker.result_dict(edge, {}, version="6.1.0"))
    assert [g.file.path for g in result.ends] == ["jobs/a/end0.xyz", "jobs/a/end1.xyz"]
    assert result.irc_connected_to_source
    assert parse({**data, "version": "6.0.0"}).kind == FailureKind.METHOD_MISMATCH
    failure = {"version": "6.1.0", "failure": {"kind": "scf_not_converged", "reason": "scc"}}
    assert parse(failure).kind == FailureKind.SCF_NOT_CONVERGED
    assert "boom" in parse(data, returncode=1).reason
