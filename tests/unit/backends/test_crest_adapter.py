"""CREST adapter: settings → CLI, run classification, registry construction (§6.3)."""

import numpy as np
import pytest

from hfauto.backends import engines
from hfauto.backends.crest import command, topology_removed
from hfauto.backends.protocols import Capability
from hfauto.backends.protocols import ConformerSettings as CS
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.evidence import FailureKind as K
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import CommandResult

GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)
ION = Molecule(XYZ(["N", "H", "H", "H", "H", "F"], np.arange(18.0).reshape(6, 3)), -1, 2)


def test_command_maps_every_setting():
    water = GFN2.model_copy(update={"solvation": "alpb:water"})
    argv = command(ION, water, CS(nci=True, threads=2, ewin_kcal=4.0, notopo_atoms=(0, 4, 5)),
                   executable="/opt/crest")
    assert " ".join(argv) == ("/opt/crest input.xyz --gfn2 --nci --quick -T 2 --ewin 4 --chrg -1"
                              " --uhf 1 --notopo 1,5,6 --alpb water")
    plain = command(ION, GFN2, CS(nci=False, quick=False, topology="noref"))
    assert "--nci" not in plain and "--quick" not in plain and plain[-1] == "--noreftopo"
    assert command(ION, GFN2, CS(nci=False, topology="off"))[-1] == "--notopo"
    assert topology_removed("CREGEN> number of topology-based structure removals: 23\n"
                            "CREGEN> number of topology-based structure removals: 22\n") == 45


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
                inputs={"molecule": ION, "method": GFN2, "settings": CS(nci=True, threads=3)})
    cmd = crest._adapter.prepare(task, workdir)
    assert (workdir / "input.xyz").is_file() and cmd.env["OMP_NUM_THREADS"] == "3,1"
    assert "--scratch" not in cmd.argv and not (tmp_path / "scratch").exists()
    (workdir / "stdout.txt").write_text("Version 3.0.2, Sun\n", encoding="utf-8")
    result = CommandResult(rc, timed_out, None, 1.0, workdir / "stdout.txt", workdir / "e.txt")
    assert crest._adapter.parse(task, workdir, result).kind is kind
