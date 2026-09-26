"""Resume after an interruption, --from and the SiteLock on a whole run (design §7.1, §7.4)."""

import sys
from functools import partialmethod
from types import SimpleNamespace

import fakes
import pytest

from hfauto.backends import engines
from hfauto.backends.protocols import Capability
from hfauto.chemistry.xyz import write_xyz
from hfauto.core.evidence import Evidence
from hfauto.core.method import ExecutionSpec
from hfauto.execution.jobs import Task
from hfauto.execution.lock import SiteLock
from hfauto.execution.process import Command
from hfauto.pipeline.config import ResolvedConfig
from hfauto.pipeline.runner import run_pipeline

pytestmark = pytest.mark.integration
STAGES = [{"id": "structures", "stage": "structures"},
          {"id": "dft", "stage": "minima", "level": "dft", "engine": "nwchem", "method": "svp",
           "select": {"include": "all"}},
          {"id": "sp", "stage": "sp", "engine": "nwchem", "methods": ["tzvp"],
           "targets": "all_minima"},
          {"id": "report", "stage": "report"}]
ADAPTER = SimpleNamespace(  # one trivial process per job; the result is computed in parse
    result_type=Evidence, continuation=lambda *_: None,
    prepare=lambda task, workdir: Command(argv=(sys.executable, "-c", ""), cwd=workdir),
    parse=lambda task, workdir, result: task.inputs["compute"]())


class JobQM(fakes.FakeQM):
    """FakeQM whose calls are JobRunner jobs; ``crash`` interrupts the second single point."""

    def __init__(self, root, jobs, crash):
        super().__init__(root, fakes.double_well())
        self.jobs, self.crash = jobs, crash

    def _job(self, kind, mol, method, *, deadline=None, **kw):
        def compute():
            if self.crash and kind == "energy" and "energy" in self.calls:
                raise RuntimeError("interrupted")
            return getattr(fakes.FakeQM, kind)(self, mol, method, **kw)

        key = {"mol": mol.fingerprint(), "method": method.signature(), "kw": sorted(kw.items())}
        task = Task("fake", "0", kind, key, ExecutionSpec(), {"compute": compute})
        return self.jobs.run(task, ADAPTER, deadline=deadline)

    energy = partialmethod(_job, "energy")
    optimize = partialmethod(_job, "optimize")
    frequencies = partialmethod(_job, "frequencies")


def test_interrupted_run_resumes_and_from_reruns_the_downstream_stages(tmp_path, tmp_run):
    ends = [{"id": n, "role": "endpoint", "xyz": write_xyz(fakes.double_well().molecule(n).xyz,
                                                           tmp_path / f"{n}.xyz")}
            for n in ("reactant", "product")]
    resolved = ResolvedConfig.model_validate({
        "pipeline": {"pipeline_id": "resume", "stages": STAGES}, "code_version": "test",
        "system": {"system_id": "dw", "species": ends},
        "site": {"site": "fake", "scratch_root": str(tmp_path / "scratch"), "cores": 1,
                 "engines": {"nwchem": {"version": "0"}}},
        "methods": {b: {"id": b, "kind": "dft", "functional": "pbe0", "basis": b}
                    for b in ("svp", "tzvp")}})
    made: list[JobQM] = []  # one engine (and JobRunner) per run; the first one crashes

    def factory(*, jobs, site):
        made.append(JobQM(tmp_run, jobs, crash=not made))
        return made[-1]

    def run(**kw):  # -> (status per stage, (hits, misses) of this run's JobRunner)
        states, stats = run_pipeline(resolved, tmp_run, **kw).read_state(), made[-1].jobs.stats()
        return {s.stage_id: s.status for s in states}, (stats.hits, stats.misses)

    with engines.override(Capability.QM, "nwchem", factory):
        with pytest.raises(RuntimeError, match="interrupted"):
            run()
        assert len(made[0].calls) == 5  # opt + freq per endpoint and one single point
        done = {s["id"]: "done" for s in STAGES}
        assert run() == (done, (1, 1))  # dft is skipped; the finished single point is a hit
        assert made[1].calls == ["energy"]
        assert run(start="dft", stop="dft") == ({**done, "sp": "stale", "report": "stale"}, (4, 0))
        assert not made[2].calls and run() == (done, (2, 0))  # sp and report run again
    with SiteLock(tmp_path / "scratch", tmp_path / "other"), pytest.raises(RuntimeError,
                                                                           match="site lock"):
        run_pipeline(resolved, tmp_run)  # another run holds the site
