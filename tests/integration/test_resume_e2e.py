"""Resume after an interruption or an external stop (SIGTERM), --from and the SiteLock on a
whole run (design §7.1, §7.4)."""

import signal
import sys
from functools import partialmethod
from types import SimpleNamespace

import fakes
import pytest
import yaml
from typer.testing import CliRunner

from hfauto.backends import engines
from hfauto.backends.protocols import Capability
from hfauto.chemistry.xyz import write_xyz
from hfauto.cli.main import app
from hfauto.core.evidence import Evidence
from hfauto.core.manifest import load_manifest
from hfauto.core.method import ExecutionSpec
from hfauto.execution import process
from hfauto.execution.jobs import Task
from hfauto.execution.lock import SiteLock
from hfauto.execution.process import Command
from hfauto.pipeline import preflight
from hfauto.pipeline.config import ResolvedConfig
from hfauto.pipeline.layout import RunLayout
from hfauto.pipeline.runner import run_pipeline

pytestmark = pytest.mark.integration
STAGES = [{"id": "structures", "stage": "structures"},
          {"id": "screen", "stage": "minima", "level": "screen", "engine": "nwchem",
           "method": "tzvp"},
          {"id": "dft", "stage": "minima", "level": "dft", "engine": "nwchem", "method": "svp",
           "select": {"include": "all"}},
          {"id": "report", "stage": "report"}]
ADAPTER = SimpleNamespace(  # one process per job (``code``); the result is computed in parse
    result_type=Evidence, continuation=lambda *_: None,
    prepare=lambda task, workdir: Command(argv=(sys.executable, "-c", task.inputs["code"]),
                                          cwd=workdir),
    parse=lambda task, workdir, result: task.inputs["compute"]())
SIGTERM_PARENT = "import os, signal, time; os.kill(os.getppid(), signal.SIGTERM); time.sleep(60)"


class JobQM(fakes.FakeQM):
    """FakeQM whose calls are JobRunner jobs; ``crash`` interrupts the seventh call (the dft
    stage's second optimization); the process of the call numbered ``term`` (from 0) sends
    SIGTERM to this one, as a batch system would."""

    def __init__(self, root, jobs, crash=False, term=None):
        super().__init__(root, fakes.double_well())
        self.jobs, self.crash, self.term = jobs, crash, term

    def _job(self, kind, mol, method, **kw):
        def compute():
            if self.crash and len(self.calls) == 6:
                raise RuntimeError("interrupted")
            return getattr(fakes.FakeQM, kind)(self, mol, method, **kw)

        key = {"mol": mol.fingerprint(), "method": method.signature(), "kw": sorted(kw.items())}
        code = SIGTERM_PARENT if len(self.calls) == self.term else ""
        task = Task("fake", "0", kind, key, ExecutionSpec(), {"compute": compute, "code": code})
        return self.jobs.run(task, ADAPTER)

    energy = partialmethod(_job, "energy")
    optimize = partialmethod(_job, "optimize")
    frequencies = partialmethod(_job, "frequencies")


def test_interrupted_run_resumes_and_from_reruns_the_downstream_stages(tmp_path, tmp_run,
                                                                       monkeypatch):
    ends = [{"id": n, "role": "endpoint", "multiplicity": 1,
             "xyz": write_xyz(fakes.double_well().molecule(n).xyz, tmp_path / f"{n}.xyz")}
            for n in ("reactant", "product")]
    resolved = ResolvedConfig.model_validate({
        "pipeline": {"pipeline_id": "resume", "stages": STAGES}, "code_version": "test",
        "system": {"system_id": "dw", "species": ends},
        "site": {"site": "fake", "scratch_root": str(tmp_path / "scratch"), "cores": 1,
                 "engines": {"nwchem": {"version": "0"}}},
        "methods": {b: {"id": b, "kind": "dft", "functional": "pbe0", "basis": b}
                    for b in ("svp", "tzvp")}})
    made: list[JobQM] = []  # one engine (and JobRunner) per run; the first one crashes

    def create(capability, name, *, jobs, site):  # engines.create for the only engine
        assert (capability, name) == (Capability.QM, "nwchem")
        made.append(JobQM(tmp_run, jobs, crash=not made))
        return made[-1]

    def run(**kw):  # -> (status per stage, (hits, misses) of this run's JobRunner)
        states, stats = run_pipeline(resolved, tmp_run, **kw).read_state(), made[-1].jobs.stats()
        return {s.stage_id: s.status for s in states}, (stats.hits, stats.misses)

    monkeypatch.setattr(engines, "create", create)
    with pytest.raises(RuntimeError, match="interrupted"):
        run()
    assert len(made[0].calls) == 6  # opt + freq per endpoint in screen, one in dft
    done = {s["id"]: "done" for s in STAGES}
    assert run() == (done, (2, 2))  # screen is skipped; dft's finished opt + freq are hits
    assert made[1].calls == ["optimize", "frequencies"]
    assert run(start="screen", stop="screen") == (
        {**done, "dft": "stale", "report": "stale"}, (4, 0))
    assert not made[2].calls and run() == (done, (4, 0))  # dft and report run again
    with SiteLock(tmp_path / "scratch", tmp_path / "other"), pytest.raises(RuntimeError,
                                                                           match="site lock"):
        run_pipeline(resolved, tmp_run)  # another run holds the site


@pytest.fixture
def stop_run(tmp_path, monkeypatch):
    """``stop_run(run_dir, term=None)``: ``hfauto run`` of structures → dft → reaction-paths on
    the double well, its QM calls JobRunner jobs (JobQM, the one numbered ``term`` sending
    SIGTERM); -> (exit code, {stage id: StageState}). Every run starts as a new process would."""
    pes = fakes.double_well()
    (tmp_path / "pipelines" / "methods").mkdir(parents=True)
    files = {
        "pipelines/stop.yaml": {"pipeline_id": "stop", "stages": [
            {"id": "structures", "stage": "structures"},
            {"id": "dft", "stage": "minima", "level": "dft", "engine": "nwchem", "method": "svp",
             "select": {"include": "all"}},
            {"id": "paths", "stage": "reaction-paths", "method": "svp", "engines": {
                "qm": "nwchem", "saddle": "nwchem_saddle", "path": "nwchem_string"}}]},
        "pipelines/methods/svp.yaml": {"id": "svp", "kind": "dft", "functional": "pbe0",
                                       "basis": "def2-svp"},
        "system.yaml": {"system_id": "dw", "species": [
            {"id": n, "role": "endpoint", "multiplicity": 1,
             "xyz": str(write_xyz(pes.molecule(n).xyz, tmp_path / f"{n}.xyz"))}
            for n in ("reactant", "product")],
            "reactions": [{"id": "rx", "reactant": "reactant", "product": "product"}]},
        "site.yaml": {"site": "fake", "scratch_root": str(tmp_path / "scratch"), "cores": 1,
                      "engines": {n: {"version": "0"}
                                  for n in ("nwchem", "nwchem_saddle", "nwchem_string")}},
    }
    for name, data in files.items():
        (tmp_path / name).write_text(yaml.safe_dump(data), encoding="utf-8")
    saved = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
    monkeypatch.setattr(preflight, "preflight", lambda resolved, dry_run=False: [])

    def run(run_dir, term=None):
        def create(capability, name, *, jobs, site):
            root = jobs.store.run_dir
            if capability is Capability.QM:
                return JobQM(root, jobs, term=term)
            return (fakes.FakeSaddle if capability is Capability.SADDLE else fakes.FakePath)(
                root, pes)

        monkeypatch.setattr(engines, "create", create)
        monkeypatch.setattr(process, "_stopping", False)  # a new process
        result = CliRunner().invoke(app, [
            "run", str(tmp_path / "pipelines" / "stop.yaml"), "--system",
            str(tmp_path / "system.yaml"), "--site", str(tmp_path / "site.yaml"),
            "--run-dir", str(run_dir)])
        return result.exit_code, {s.stage_id: s for s in RunLayout(run_dir).read_state()}

    yield run
    for signum, handler in saved.items():
        signal.signal(signum, handler)


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals")
def test_sigterm_leaves_the_stage_incomplete_and_a_rerun_resumes_it(tmp_path, stop_run):
    """X7-1: SIGTERM during reaction-paths (from the process of its first QM job, the seed's
    freq) kills that job, stores nothing for it and leaves the stage incomplete with exit code 1;
    the rerun skips the done stages and gives the records of an uninterrupted run."""
    run, other = tmp_path / "run", tmp_path / "other"
    code, states = stop_run(run, term=4)  # calls 0-3: opt + freq per endpoint in dft
    assert code == 1 and {k: s.status for k, s in states.items()} == {
        "structures": "done", "dft": "done", "paths": "incomplete"}
    assert len(list((run / "jobs").glob("*/*/result.json"))) == 4  # the killed job left none
    assert CliRunner().invoke(app, ["status", str(run)]).exit_code == 1
    code, resumed = stop_run(run)
    assert code == 0 and resumed["paths"].status == "done"
    assert resumed["dft"] == states["dft"]  # skipped: not executed again
    assert resumed["paths"].jobs.misses > 0
    assert stop_run(other)[0] == 0

    def records(run_dir):
        manifest = load_manifest(RunLayout(run_dir).manifest_path("paths"))
        return [(a.artifact_id, a.type, a.payload) for a in manifest.artifacts]

    assert records(run) == records(other) and records(run)
