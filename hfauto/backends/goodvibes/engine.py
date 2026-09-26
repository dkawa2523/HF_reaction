"""GoodVibes thermochemistry engine (design §6.3): the Python API run in a worker process.

The engine hands the worker its own projected frequencies (as ``thermo_frequencies``: no
negative mode reaches GoodVibes), isotopic masses and rotational constants from a freq
Evidence; it never reads GoodVibes text output. The cache key is the freq job_key, the Hessian
sha, the frequencies sent, the temperatures and the full settings.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import ClassVar, Literal

from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import Requirements, ThermoResult
from hfauto.chemistry.thermo import settings_sha, thermo_frequencies
from hfauto.chemistry.xyz import geometry_fingerprint, read_xyz
from hfauto.core.evidence import Evidence, Failure, FailureKind
from hfauto.core.method import Deadline, EngineSite, MethodSpec, ThermoSettings
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.process import Command, CommandResult
from hfauto.execution.worker import RESULT_NAME

_WORKER = "hfauto.backends.goodvibes.worker:compute"
_JOB_NAME = "job.json"


class ThermoBatch(BaseModel):
    """What the JobStore keeps for one worker run: every (settings, temperature) result."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    kind: Literal["thermo_batch"] = "thermo_batch"
    results: tuple[ThermoResult, ...]
    job_key: str = ""


def _tail(path: Path, n: int = 300) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[-n:] if path.is_file() else ""


class _Adapter:
    result_type = ThermoBatch

    def __init__(self, python: str) -> None:
        self._python = python

    def prepare(self, task: Task, workdir: Path) -> Command:
        job = workdir / _JOB_NAME
        job.write_text(json.dumps(task.inputs["job"]), encoding="utf-8")
        argv = (self._python, "-m", "hfauto.execution.worker", _WORKER, str(job))
        return Command(argv=argv, cwd=workdir)

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> ThermoBatch | Failure:
        if result.timed_out:
            return Failure(kind=FailureKind.TIMEOUT, reason="goodvibes worker timed out")
        path = workdir / RESULT_NAME
        if result.returncode != 0 or not path.is_file():
            reason = f"worker rc={result.returncode}: {_tail(result.stderr)}"
            return Failure(kind=FailureKind.NONZERO_EXIT, reason=reason)
        data = json.loads(path.read_text(encoding="utf-8"))
        if "error" in data:
            return Failure(kind=FailureKind(data["error"]["kind"]), reason=data["error"]["reason"])
        return ThermoBatch(results=tuple(
            ThermoResult(**r, notes=(), job_key="") for r in data["results"]))

    def monitor(self, task: Task) -> Callable[[Path], str | None] | None:
        return None

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None:
        return None  # a failed thermochemistry is final (G = None downstream)


class GoodVibesEngine:
    """ThermoEngine; ``site.python`` is the worker interpreter (goodvibes and pymsym)."""

    name: ClassVar[str] = "goodvibes"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self._jobs = jobs
        self._site = site

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(
            python_modules=("goodvibes", "pymsym"),
            version_command=("python", "-c", "print(__import__('goodvibes').__version__)"),
        )

    def supports(self, method: MethodSpec) -> bool:
        return True  # thermochemistry takes no method; the level comes from the Evidence

    def thermo(
        self,
        freq: Evidence,
        settings: Sequence[ThermoSettings],
        *,
        temperatures_K: Sequence[float],
        saddle: bool = False,
        deadline: Deadline | None = None,
    ) -> list[ThermoResult] | Failure:
        if freq.task != "freq" or freq.frequencies_cm1 is None or freq.hessian is None:
            return Failure(kind=FailureKind.INPUT_INVALID, reason=f"not_freq_task:{freq.task}")
        xyz = read_xyz(self._jobs.store.run_dir / freq.final.file.path)
        if geometry_fingerprint(xyz.symbols, xyz.coords) != freq.final.fingerprint:
            return Failure(kind=FailureKind.INPUT_INVALID, reason="geometry_mismatch")
        job = {
            "version_pin": self._site.version,
            "symbols": xyz.symbols, "coords": xyz.coords.tolist(),
            "energy_hartree": freq.energy_hartree,
            "frequencies_cm1": list(thermo_frequencies(freq.frequencies_cm1, saddle=saddle)),
            "n_external": freq.n_external, "charge": freq.level.charge,
            "multiplicity": freq.level.multiplicity, "temperatures_K": list(temperatures_K),
            "settings": [s.model_dump(mode="json") | {"sha": settings_sha(s)} for s in settings],
        }
        task = Task(
            engine=self.name, version_pin=self._site.version, kind="thermo",
            key_payload={"freq": freq.job_key, "hessian_sha": freq.hessian.sha256} | {
                k: job[k] for k in ("frequencies_cm1", "temperatures_K", "settings")},
            execution=self._site.execution, inputs={"job": job},
        )
        out = self._jobs.run(task, _Adapter(self._site.python or sys.executable),
                             deadline=deadline)
        if isinstance(out, Failure):
            return out
        return [r.model_copy(update={"job_key": out.job_key}) for r in out.results]
