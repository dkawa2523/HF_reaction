"""``pysis_neb``: xTB climbing-image NEB between fixed ends with pysisyphus (design §6.3).

xTB only, never DFT (CH-01, CH-02). The caller's initial path (hfauto's IDPP between the DFT
minima) is relaxed in Cartesian coordinates by a climbing-image NEB whose end images stay
fixed (the COS default), with LBFGS for at most ``NEB_MAX_CYCLES``; an unconverged NEB is still
a path. The same input then optimizes the TS from the climbing image (rsprfo on the xTB
Hessian, following the imaginary mode with the largest overlap with the HEI tangent). The run
happens in a worker subprocess whose PATH starts with the site's xtb directory, because the
pysisyphus calculator starts ``xtb`` from PATH. The version pin is the xTB version (the PES);
pysisyphus itself is pinned by the production extra.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from hfauto.backends.protocols import Requirements
from hfauto.backends.pysis.worker import neb_images
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz_trajectory, write_xyz, written_geometry
from hfauto.core.evidence import Failure, FailureKind, FileRef, Geometry, Level, PathProfile
from hfauto.core.method import Deadline, EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import Task
from hfauto.execution.process import Command, CommandResult, resolve_executable
from hfauto.execution.worker import RESULT_NAME

if TYPE_CHECKING:
    from hfauto.execution.jobs import JobRunner

INITIAL_NAME, JOB_NAME = "initial.trj", "job.json"
WORKER = "hfauto.backends.pysis.worker:run_neb"
NEB_MAX_CYCLES = 100
_XTB_DEFAULT_ETEMP_K = 300.0


def neb_input(mol: Molecule, method: MethodSpec, *, threads: int) -> dict[str, Any]:
    """The pysisyphus run dictionary: CI-NEB on the initial path, then TSOpt from its CI."""
    calc: dict[str, Any] = {"type": "xtb", "gfn": method.gfn, "charge": mol.charge,
                            "mult": mol.multiplicity, "pal": threads}
    if method.electronic_temperature_K is not None:
        calc["etemp"] = method.electronic_temperature_K
    return {"geom": {"type": "cart", "fn": INITIAL_NAME}, "calc": calc,
            "cos": {"type": "neb", "climb": True},
            "opt": {"type": "lbfgs", "max_cycles": NEB_MAX_CYCLES},
            "tsopt": {"type": "rsprfo"}}


def _fail(kind: FailureKind, reason: str) -> Failure:
    return Failure(kind=kind, reason=reason)


class NEBAdapter:
    """render → worker → parse for one ``pysis_neb`` task (``jobs.Adapter``)."""

    result_type = PathProfile

    def __init__(self, site: EngineSite, file_ref: Callable[[Path], FileRef]) -> None:
        self.site = site
        self.file_ref = file_ref

    def prepare(self, task: Task, workdir: Path) -> Command:
        shutil.copyfile(task.inputs["initial_path"], workdir / INITIAL_NAME)
        job = workdir / JOB_NAME
        job.write_text(json.dumps({"run_dict": task.inputs["run_dict"]}), encoding="utf-8")
        python = self.site.python or sys.executable
        argv = (python, "-m", "hfauto.execution.worker", WORKER, str(job))
        return Command(argv=argv, cwd=workdir, env=self.env(task))

    def env(self, task: Task) -> dict[str, str]:
        env = dict(task.execution.env)
        xtb = resolve_executable("xtb", self.site.executables.get("xtb"))
        if xtb is not None:
            path = env.get("PATH", os.environ.get("PATH", ""))
            env["PATH"] = os.pathsep.join((str(Path(xtb).parent), path))
        return env

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> None:
        return None

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> PathProfile | Failure:
        if result.timed_out:
            return _fail(FailureKind.TIMEOUT, "pysisyphus timed out")
        if result.returncode != 0:
            lines = result.stderr.read_text(encoding="utf-8", errors="replace").splitlines()
            last = lines[-1] if lines else ""
            return _fail(FailureKind.NONZERO_EXIT, f"worker exited {result.returncode}: {last}")
        path = workdir / RESULT_NAME
        if not path.is_file():
            return _fail(FailureKind.INCOMPLETE_OUTPUT, f"{RESULT_NAME} missing")
        return self.profile(task, workdir, json.loads(path.read_text(encoding="utf-8")))

    def profile(self, task: Task, workdir: Path, data: dict[str, Any]) -> PathProfile | Failure:
        """PathProfile of the NEB's final images; the TS becomes an xyz file."""
        start: Molecule = task.inputs["start"]
        method: MethodSpec = task.inputs["method"]
        if not data.get("xtb_version"):
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "no xtb version in the calculator output")
        symbols, images = list(start.xyz.symbols), neb_images(workdir)
        try:
            frames = [] if images is None else read_xyz_trajectory(images)
        except (OSError, ValueError):
            frames = []
        wrong = len(frames) != task.inputs["images"] or any(f.symbols != symbols for f in frames)
        if images is None or wrong:
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "neb_images")
        etemp = method.electronic_temperature_K
        level = Level(program="xtb", version=data["xtb_version"], method=f"gfn{method.gfn}",
                      charge=start.charge, multiplicity=start.multiplicity,
                      electronic_temperature_K=_XTB_DEFAULT_ETEMP_K if etemp is None else etemp)
        problems = level_mismatches(method, level, version_pin=task.version_pin)
        if problems:
            return _fail(FailureKind.METHOD_MISMATCH, "; ".join(problems))
        ts = data.get("ts")
        return PathProfile(engine=PysisNEB.name, level=level, images=self.file_ref(images),
                           ts=self._geometry(workdir / "ts.xyz", symbols, ts) if ts else None)

    def _geometry(self, path: Path, symbols: list[str], coords: Any) -> Geometry:
        return written_geometry(write_xyz(XYZ(symbols, np.asarray(coords, dtype=float)), path),
                                self.file_ref)


class PysisNEB:
    """PathEngine: xTB CI-NEB from ``initial_path`` with fixed ends, TS from the CI."""

    name: ClassVar[str] = "pysis_neb"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self.jobs = jobs
        self.site = site
        self.adapter = NEBAdapter(site, jobs.store.file_ref)

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(executables=("xtb",), python_modules=("pysisyphus",),
                            version_command=("xtb", "--version"))

    def supports(self, method: MethodSpec) -> bool:
        return method.kind == "xtb" and method.gfn is not None

    def task(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
             initial_path: FileRef) -> Task:
        return Task(
            engine=self.name, version_pin=self.site.version, kind="neb",
            key_payload={"method": method.signature(), "start": start.fingerprint(),
                         "end": end.fingerprint(), "images": images,
                         "initial_path": initial_path.sha256},
            execution=self.site.execution,
            inputs={"start": start, "method": method, "images": images,
                    "initial_path": self.jobs.store.resolve(initial_path),
                    "run_dict": neb_input(start, method, threads=self.site.execution.threads)},
        )

    def find_path(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
                  initial_path: FileRef, deadline: Deadline | None = None
                  ) -> PathProfile | Failure:
        task = self.task(start, end, method, images=images, initial_path=initial_path)
        return self.jobs.run(task, self.adapter, deadline=deadline)
