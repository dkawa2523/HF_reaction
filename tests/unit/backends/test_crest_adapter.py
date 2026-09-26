"""CREST adapter: settings → CLI, run classification, registry construction (§6.3)."""

import numpy as np
import pytest
from pydantic import ValidationError

from hfauto.backends import engines
from hfauto.backends.crest import command, topology_removed
from hfauto.backends.protocols import Capability
from hfauto.backends.protocols import ConformerSettings as CS
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.evidence import FailureKind as K
from hfauto.core.method import EngineSite, ExecutionSpec, MethodSpec
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)
ION = Molecule(XYZ(["N", "H", "H", "H", "H", "F"], np.arange(18.0).reshape(6, 3)), -1, 2)


def test_command_maps_every_setting():
    water = GFN2.model_copy(update={"solvation": "alpb:water"})
    argv = command(ION, water, CS(nci=True, ewin_kcal=4.0, notopo_atoms=(0, 4, 5)), threads=2,
                   executable="/opt/crest")  # a composition: --noopt
    assert " ".join(argv) == ("/opt/crest input.xyz --gfn2 --nci --quick -T 2 --ewin 4 --chrg -1"
                              " --uhf 1 --notopo 1,5,6 --noopt --alpb water")
    plain = command(ION, GFN2, CS(quick=False, topology="noref"), threads=1)  # a monomer retry
    assert not {"--nci", "--noopt", "--quick", "--notopo"} & set(plain)
    assert plain[-1] == "--noreftopo"
    with pytest.raises(ValidationError):
        CS(topology="off")
    assert topology_removed("CREGEN> number of topology-based structure removals: 23\n"
                            "CREGEN> number of topology-based structure removals: 22\n") == 45


def test_search_keys_noopt_and_runs_on_site_threads(tmp_path, monkeypatch):
    site = EngineSite(version="3.0.2", execution=ExecutionSpec(threads=3))
    crest = engines.create(Capability.CONFORMERS, "crest", site=site,
                           jobs=JobRunner(JobStore(tmp_path / "jobs"), cores=4))
    tasks = []
    monkeypatch.setattr(crest.jobs, "run", lambda task, adapter, deadline=None: tasks.append(task))
    for nci in (False, True):
        crest.search(ION, GFN2, CS(nci=nci))
    assert [t.key_payload["noopt"] for t in tasks] == [False, True]
    assert all(t.execution.threads == 3 for t in tasks)
    cmd = crest._adapter.prepare(tasks[1], tmp_path / "w")
    assert cmd.argv[cmd.argv.index("-T") + 1] == "3" and "--noopt" in cmd.argv
    assert cmd.env["OMP_NUM_THREADS"] == "3,1"


@pytest.mark.parametrize("rc,timed_out,pin,kind", [  # the last: no crest_conformers.xyz
    (0, True, "3.0.2", K.TIMEOUT), (0, False, "3.1.0", K.METHOD_MISMATCH),
    (1, False, "3.0.2", K.NONZERO_EXIT), (0, False, "3.0.2", K.INCOMPLETE_OUTPUT)])
def test_run_classification(tmp_path, rc, timed_out, pin, kind):
    workdir = tmp_path / "jobs" / "ab" / "abcd" / "attempt_00"
    site = EngineSite(version=pin, scratch_dir=str(tmp_path / "scratch"))  # never used
    crest = engines.create(Capability.CONFORMERS, "crest", site=site,  # checks the Protocol
                           jobs=JobRunner(JobStore(tmp_path / "jobs"), cores=4))
    assert crest.supports(GFN2) and not crest.supports(GFN2.model_copy(update={"gfn": None}))
    task = Task(engine="crest", version_pin=pin, kind="conformers", key_payload={},
                execution=crest.site.execution,
                inputs={"molecule": ION, "method": GFN2, "settings": CS(nci=True)})
    cmd = crest._adapter.prepare(task, workdir)
    assert (workdir / "input.xyz").is_file() and cmd.env["OMP_NUM_THREADS"] == "1,1"
    assert "--scratch" not in cmd.argv and not (tmp_path / "scratch").exists()
    (workdir / "stdout.txt").write_text("Version 3.0.2, Sun\n", encoding="utf-8")
    result = CommandResult(rc, timed_out, None, 1.0, workdir / "stdout.txt", workdir / "e.txt")
    assert crest._adapter.parse(task, workdir, result).kind is kind
