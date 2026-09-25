"""``pysis_gs``: pysisyphus growing string (GS) on the native xTB calculator (design §6.3).

xTB only, never DFT (CH-01, CH-02). One pysisyphus input holds the GS and, with
``refine_ts=True``, the TS optimization (``run_tsopt_from_cos`` picks the root from the
overlap with the splined HEI tangent, chem 1). The run happens in a worker subprocess whose
PATH starts with the site's xtb directory, because the pysisyphus calculator starts ``xtb``
from PATH.

- GS: ``max_nodes = images - 2``, ``climb``, ``climb_rms 5e-3``, the ``string`` optimizer.
  Coordinates are Cartesian for a linear endpoint or two or more fragments, else DLC.
  "Linear" includes near-linear ends (every atom within 0.3 Å of the principal axis): DLC
  crashes there (HCN bent by 10° on the way to HNC, while 30° works). A Cartesian string
  between exactly linear ends stays on the axis by symmetry and runs atoms into each other,
  so such an end is first offset by ±0.01 Å (alternating atoms) across its axis.
- No ``interpol`` section: pysisyphus 1.0's GrowingString grows new nodes from the two
  endpoints itself and takes exactly two input geometries (checked in its source and on a
  real HCN/HNC run), so the default is kept.
- The version pin is the xTB version (the PES); pysisyphus itself is pinned by the
  production extra.
"""

from __future__ import annotations

import json
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from hfauto.backends.protocols import Requirements
from hfauto.chemistry.topology import fragments
from hfauto.chemistry.vibrations import external_basis
from hfauto.chemistry.xyz import XYZ, Molecule, geometry_fingerprint, read_xyz, write_xyz
from hfauto.chemistry.xyz_trajectory import write_xyz_trajectory
from hfauto.core.evidence import Failure, FailureKind, FileRef, Geometry, Level, PathProfile
from hfauto.core.method import Deadline, EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import Task
from hfauto.execution.process import Command, CommandResult, resolve_executable
from hfauto.execution.worker import RESULT_NAME

if TYPE_CHECKING:
    from hfauto.execution.jobs import JobRunner

START_NAME, END_NAME, JOB_NAME = "start.xyz", "end.xyz", "job.json"
WORKER = "hfauto.backends.pysis.worker:run_growing_string"
_XTB_DEFAULT_ETEMP_K = 300.0
_NEAR_LINEAR_A = 0.3
_OFF_AXIS_A = 0.01


def _axis(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unit principal axis and each atom's offset from the line through the centroid."""
    centered = coords - coords.mean(axis=0)
    axis = np.linalg.svd(centered)[2][0]
    return axis, centered - np.outer(centered @ axis, axis)


def _linear(mol: Molecule) -> bool:
    offsets = _axis(mol.xyz.coords)[1]
    return float(np.linalg.norm(offsets, axis=1).max()) <= _NEAR_LINEAR_A


def off_axis(mol: Molecule) -> Molecule:
    """An exactly linear end (5 external modes) offset by ±0.01 Å across its axis."""
    if external_basis(mol.xyz.symbols, mol.xyz.coords).shape[1] != 5:
        return mol
    axis = _axis(mol.xyz.coords)[0]
    ref = np.eye(3)[int(np.argmin(np.abs(axis)))]
    across = ref - (ref @ axis) * axis
    signs = np.where(np.arange(len(mol.xyz.symbols)) % 2 == 0, 1.0, -1.0)
    shift = _OFF_AXIS_A * signs[:, None] * across / np.linalg.norm(across)
    return Molecule(XYZ(list(mol.xyz.symbols), mol.xyz.coords + shift), mol.charge,
                    mol.multiplicity)


def _split(mol: Molecule) -> bool:
    return len(fragments(mol.xyz.symbols, mol.xyz.coords)) >= 2


def gs_input(start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
             refine_ts: bool, threads: int) -> dict[str, Any]:
    """The pysisyphus run dictionary for one GS (+ TS optimization)."""
    split = _split(start) or _split(end)
    cartesian = split or _linear(start) or _linear(end)
    calc: dict[str, Any] = {"type": "xtb", "gfn": method.gfn, "charge": start.charge,
                            "mult": start.multiplicity, "pal": threads}
    if method.electronic_temperature_K is not None:
        calc["etemp"] = method.electronic_temperature_K
    if method.solvation:
        calc["alpb"] = method.solvation.split(":", 1)[1]
    run: dict[str, Any] = {
        "geom": {"type": "cart" if cartesian else "dlc", "fn": [START_NAME, END_NAME]},
        "calc": calc,
        "cos": {"type": "gs", "max_nodes": images - 2, "climb": True, "climb_rms": 5e-3},
        "opt": {"type": "string"},
    }
    if refine_ts:
        tsopt: dict[str, Any] = {"type": "rsirfo", "thresh": "gau"}
        if split:
            tsopt["geom"] = {"type": "tric"}
        run["tsopt"] = tsopt
    return run


def _fail(kind: FailureKind, reason: str) -> Failure:
    return Failure(kind=kind, reason=reason)


class GrowingStringAdapter:
    """render → worker → parse for one ``pysis_gs`` task (``jobs.Adapter``)."""

    result_type = PathProfile

    def __init__(self, site: EngineSite, file_ref: Callable[[Path], FileRef]) -> None:
        self.site = site
        self.file_ref = file_ref

    def prepare(self, task: Task, workdir: Path) -> Command:
        task.inputs["start"].write(workdir / START_NAME)
        task.inputs["end"].write(workdir / END_NAME)
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

    def monitor(self, task: Task) -> None:
        return None

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> None:
        return None

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> PathProfile | Failure:
        if result.timed_out:
            return _fail(FailureKind.TIMEOUT, "pysisyphus timed out")
        if result.stopped is not None:
            return _fail(FailureKind.STAGNATED, result.stopped)
        if result.returncode != 0:
            lines = result.stderr.read_text(encoding="utf-8", errors="replace").splitlines()
            last = lines[-1] if lines else ""
            return _fail(FailureKind.NONZERO_EXIT, f"worker exited {result.returncode}: {last}")
        path = workdir / RESULT_NAME
        if not path.is_file():
            return _fail(FailureKind.INCOMPLETE_OUTPUT, f"{RESULT_NAME} missing")
        return self.profile(task, workdir, json.loads(path.read_text(encoding="utf-8")))

    def profile(self, task: Task, workdir: Path, data: dict[str, Any]) -> PathProfile | Failure:
        """PathProfile from the worker's summary; images and TS become xyz files."""
        start: Molecule = task.inputs["start"]
        method: MethodSpec = task.inputs["method"]
        if not data.get("xtb_version"):
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "no xtb version in the calculator output")
        if len(data["images"]) < 2 or len(data["images"]) != len(data["energies"]):
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "images and energies do not match")
        etemp = method.electronic_temperature_K
        level = Level(program="xtb", version=data["xtb_version"], method=f"gfn{method.gfn}",
                      solvation=method.solvation, charge=start.charge,
                      multiplicity=start.multiplicity,
                      electronic_temperature_K=_XTB_DEFAULT_ETEMP_K if etemp is None else etemp)
        problems = level_mismatches(method, level, version_pin=task.version_pin)
        if problems:
            return _fail(FailureKind.METHOD_MISMATCH, "; ".join(problems))
        symbols = list(start.xyz.symbols)
        frames = [XYZ(symbols, np.asarray(c, dtype=float)) for c in data["images"]]
        images = write_xyz_trajectory(frames, workdir / "images.xyz")
        ts = data.get("ts")
        return PathProfile(
            engine=PysisGrowingString.name, level=level, images=self.file_ref(images),
            energies_hartree=tuple(data["energies"]), gmax_history=tuple(data["max_forces"]),
            program_converged=bool(data["converged"]), climbing_image=data["climbing_image"],
            ts=self._geometry(workdir / "ts.xyz", symbols, ts["coords"]) if ts else None,
            ts_energy_hartree=ts["energy"] if ts else None, job_key="",
        )

    def _geometry(self, path: Path, symbols: list[str], coords: Any) -> Geometry:
        xyz = read_xyz(write_xyz(XYZ(symbols, np.asarray(coords, dtype=float)), path))
        return Geometry(file=self.file_ref(path), symbols=tuple(xyz.symbols),
                        fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords))


class PysisGrowingString:
    """PathEngine; ``initial_path`` is ignored because the string grows from the endpoints."""

    name: ClassVar[str] = "pysis_gs"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self.jobs = jobs
        self.site = site
        self.adapter = GrowingStringAdapter(site, jobs.store.file_ref)

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(executables=("xtb",), python_modules=("pysisyphus",),
                            version_command=("xtb", "--version"))

    def supports(self, method: MethodSpec) -> bool:
        solvation_ok = method.solvation is None or method.solvation.startswith("alpb:")
        return method.kind == "xtb" and method.gfn is not None and solvation_ok

    def task(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
             refine_ts: bool) -> Task:
        run_dict = gs_input(start, end, method, images=images, refine_ts=refine_ts,
                            threads=self.site.execution.threads)
        return Task(
            engine=self.name, version_pin=self.site.version, kind="growing_string",
            key_payload={"method": method.signature(), "start": start.fingerprint(),
                         "end": end.fingerprint(), "images": images, "refine_ts": refine_ts},
            execution=self.site.execution,
            inputs={"start": off_axis(start), "end": off_axis(end), "method": method,
                    "run_dict": run_dict},
        )

    def find_path(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
                  initial_path: FileRef | None = None, refine_ts: bool = False,
                  deadline: Deadline | None = None) -> PathProfile | Failure:
        task = self.task(start, end, method, images=images, refine_ts=refine_ts)
        return self.jobs.run(task, self.adapter, deadline=deadline)
