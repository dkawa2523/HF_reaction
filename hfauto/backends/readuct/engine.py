"""SCINE ReaDuct 6.1 discovery engine (design §6.3, DISCOVERY).

Each ``explore`` call is one JobStore job: an attempt runs ``worker.run_attempt`` in a worker
subprocess whose cwd is the attempt directory (CH-31), with the timeout of
``DiscoverySettings.timeout_s``. The worker reports the observed scine-readuct version, which
must equal the site pin. There is no continuation: a failed attempt is final.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, ClassVar

from hfauto.backends.protocols import DiscoveryResult, DiscoverySettings, Requirements
from hfauto.chemistry.xyz import Molecule, geometry_fingerprint, read_xyz
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.method import Deadline, EngineSite, MethodSpec
from hfauto.core.records import ReactionTrial
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.process import Command, CommandResult
from hfauto.execution.worker import RESULT_NAME

WORKER_TARGET = "hfauto.backends.readuct.worker:run_attempt"
JOB_NAME = "job.json"
_STDERR_TAIL = 300


class ReaDuctEngine:
    name: ClassVar[str] = "readuct"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self.jobs, self.site = jobs, site

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(
            python_modules=("scine_readuct", "scine_utilities", "scine_xtb_wrapper"))

    def supports(self, method: MethodSpec) -> bool:
        return method.kind == "xtb" and method.gfn is not None and method.solvation is None

    def explore(self, source: Molecule, trial: ReactionTrial, method: MethodSpec,
                settings: DiscoverySettings, *,
                deadline: Deadline | None = None) -> DiscoveryResult | Failure:
        drive = trial.model_dump(mode="json",
                                 include={"mechanism", "associations", "dissociations"})
        task = Task(
            engine=self.name, version_pin=self.site.version, kind="discovery",
            key_payload={"method": method.signature(), "source": source.fingerprint(),
                         "trial": drive,
                         "settings": settings.model_dump(mode="json", exclude={"timeout_s"})},
            execution=self.site.execution.model_copy(update={"timeout_s": settings.timeout_s}),
            inputs={"job": {"symbols": list(source.xyz.symbols),
                            "coords": source.xyz.coords.tolist(), "charge": source.charge,
                            "multiplicity": source.multiplicity,
                            "method_family": f"GFN{method.gfn}", "trial": drive,
                            "settings": settings.model_dump(mode="json")}},
        )
        return self.jobs.run(task, _Adapter(self.site, self.jobs), deadline=deadline)


class _Adapter:
    result_type = DiscoveryResult

    def __init__(self, site: EngineSite, jobs: JobRunner) -> None:
        self.site, self.jobs = site, jobs

    def prepare(self, task: Task, workdir: Path) -> Command:
        job_path = workdir / JOB_NAME
        job_path.write_text(json.dumps(task.inputs["job"], indent=1), encoding="utf-8")
        threads = task.execution.threads
        env = {"OMP_NUM_THREADS": f"{threads},1", "OMP_STACKSIZE": "4G", **task.execution.env}
        python = self.site.python or sys.executable
        argv = (python, "-m", "hfauto.execution.worker", WORKER_TARGET, str(job_path))
        return Command(argv=argv, cwd=workdir, env=env)

    def parse(self, task: Task, workdir: Path, result: CommandResult
              ) -> DiscoveryResult | Failure:
        path = workdir / RESULT_NAME
        if result.timed_out:
            return Failure(kind=FailureKind.TIMEOUT, reason="readuct attempt timed out")
        if result.returncode != 0 or not path.is_file():
            kind = (FailureKind.NONZERO_EXIT if result.returncode != 0
                    else FailureKind.INCOMPLETE_OUTPUT)
            tail = result.stderr.read_text(encoding="utf-8", errors="replace")[-_STDERR_TAIL:]
            return Failure(kind=kind, reason=f"worker rc={result.returncode}: {tail.strip()}")
        data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
        if data["version"] != task.version_pin:
            return Failure(kind=FailureKind.METHOD_MISMATCH,
                           reason=f"scine-readuct {data['version']} != pin {task.version_pin}")
        if "failure" in data:
            return Failure(kind=FailureKind(data["failure"]["kind"]),
                           reason=data["failure"]["reason"])
        return DiscoveryResult(
            outcome=data["outcome"], reason=data["reason"],
            product=self._geometry(workdir, data["product"]),
            ts=self._geometry(workdir, data["ts"]),
            ts_imag_cm1=data["ts_imag_cm1"], dE_act_kcal=data["dE_act_kcal"],
            dE_rxn_kcal=data["dE_rxn_kcal"],
            irc_connected_to_source=data["irc_connected_to_source"],
            electronic_temperature_K=data["electronic_temperature_K"], job_key="",
        )

    def _geometry(self, workdir: Path, name: str | None) -> Geometry | None:
        if name is None:
            return None
        xyz = read_xyz(workdir / name)
        return Geometry(file=self.jobs.store.file_ref(workdir / name), symbols=tuple(xyz.symbols),
                        fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords))

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None:
        return None
