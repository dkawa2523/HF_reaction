"""xTB as a QMEngine (design §6.3): single points, optimizations and Hessians.

Runs ``xtb input.xyz`` with ``--sp`` / ``--opt vtight`` / ``--hess`` plus ``--gfn``,
``--chrg``, ``--uhf`` (multiplicity - 1), ``--alpb <solvent>`` and ``--etemp`` when the
method asks for them. An optimization counts as converged only with return code 0, no
``FAILED TO CONVERGE`` in the output and no ``NOT_CONVERGED`` file (BUG-07, CH-20); a
continuation restarts from ``xtbopt.xyz``. The ``hessian`` file (Eh/bohr², input frame)
becomes the canonical ``.npy`` and the frequencies come from ``chemistry.vibrations``.
``init_hessian`` is accepted and ignored: every optimization is vtight and xTB builds its
own model Hessian.
"""

from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from hfauto.backends.protocols import Requirements
from hfauto.chemistry.vibrations import projected_frequencies, to_canonical_npy
from hfauto.chemistry.xyz import XYZ, Molecule, geometry_fingerprint, read_xyz
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry, Level
from hfauto.core.method import Deadline, EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import Task
from hfauto.execution.process import STDOUT_NAME, Command, CommandResult, resolve_executable

if TYPE_CHECKING:
    from collections.abc import Callable

    from hfauto.execution.jobs import JobRunner

INPUT_NAME = "input.xyz"
_RUN_FLAGS = {"energy": ("--sp",), "optimize": ("--opt", "vtight"), "frequencies": ("--hess",)}
_TASKS = {"energy": "sp", "optimize": "opt", "frequencies": "freq"}
_TOTAL_ENERGY = re.compile(r"\|\s*TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh")
_CYCLE_ENERGY = re.compile(r"^\s*\* total energy\s*:\s*(-?\d+\.\d+)\s+Eh", re.MULTILINE)
_VERSION = re.compile(r"xtb version (\S+)")
_HAMILTONIAN = re.compile(r"Hamiltonian\s+(GFN\d)-xTB", re.IGNORECASE)
_ETEMP = re.compile(r"electronic temp\.\s+(\d+(?:\.\d*)?)\s+K")
_CHARGE = re.compile(r"net charge\s+(-?\d+)")
_UNPAIRED = re.compile(r"unpaired electrons\s+(\d+)")
_SOLVATION = re.compile(r"Solvation model:\s+(\S+)\s*\n\s*Solvent\s+(\S+)")
_SCF_FAILED = re.compile(r"SCF not converged|Self consistent charge iterator did not converge")


def _fail(kind: FailureKind, reason: str) -> Failure:
    return Failure(kind=kind, reason=reason)


def _termination(workdir: Path, result: CommandResult, text: str) -> Failure | None:
    if result.timed_out:
        return _fail(FailureKind.TIMEOUT, "xtb timed out")
    if "FAILED TO CONVERGE" in text or (workdir / "NOT_CONVERGED").exists():
        return _fail(FailureKind.GEOMETRY_MAXITER, "xtb optimization did not converge")
    if result.returncode != 0:
        scf = _SCF_FAILED.search(text) is not None
        kind = FailureKind.SCF_NOT_CONVERGED if scf else FailureKind.NONZERO_EXIT
        return _fail(kind, f"xtb exited with {result.returncode}")
    return None


def observed_level(text: str) -> Level | None:
    """The level xTB reports in its setup block; None when a field is missing."""
    found = [p.search(text) for p in (_VERSION, _HAMILTONIAN, _ETEMP, _CHARGE, _UNPAIRED)]
    if any(match is None for match in found):
        return None
    version, hamiltonian, etemp, charge, unpaired = (m.group(1) for m in found if m)
    solvation = _SOLVATION.search(text)
    return Level(
        program="xtb", version=version, method=hamiltonian, charge=int(charge),
        multiplicity=int(unpaired) + 1, electronic_temperature_K=float(etemp),
        solvation=f"{solvation.group(1)}:{solvation.group(2)}" if solvation else None,
    )


def read_hessian(path: Path, n_atoms: int) -> np.ndarray | None:
    """The ``$hessian`` block as a (3N, 3N) array; None when missing or malformed."""
    if not path.is_file():
        return None
    body = path.read_text(encoding="utf-8", errors="replace").split("$hessian", 1)[-1]
    try:
        values = np.array(body.split("$end", 1)[0].split(), dtype=float)
    except ValueError:  # overflowing fields are printed as asterisks
        return None
    size = 3 * n_atoms
    return values.reshape(size, size) if values.size == size * size else None


class _Adapter:
    """render → run → parse for one xTB task (``jobs.Adapter``)."""

    result_type = Evidence

    def __init__(self, site: EngineSite, file_ref: Callable[[Path], FileRef]) -> None:
        self.site = site
        self.file_ref = file_ref

    def prepare(self, task: Task, workdir: Path) -> Command:
        mol: Molecule = task.inputs["mol"]
        method: MethodSpec = task.inputs["method"]
        mol.write(workdir / INPUT_NAME)
        explicit = self.site.executables.get("xtb")
        exe = resolve_executable("xtb", explicit) or explicit or "xtb"
        argv = [exe, INPUT_NAME, *_RUN_FLAGS[task.kind], "--gfn", str(method.gfn),
                "--chrg", str(mol.charge), "--uhf", str(mol.multiplicity - 1)]
        if method.solvation:
            argv += ["--alpb", method.solvation.split(":", 1)[1]]
        if method.electronic_temperature_K is not None:
            argv += ["--etemp", str(method.electronic_temperature_K)]
        if task.kind == "optimize" and task.execution.maxiter is not None:
            argv += ["--cycles", str(task.execution.maxiter)]
        threads = task.execution.threads
        env = {"OMP_NUM_THREADS": f"{threads},1", "OMP_STACKSIZE": "4G", **task.execution.env}
        return Command(argv=tuple(argv), cwd=workdir, env=env)

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None:
        """Restart an optimization from ``xtbopt.xyz`` after a timeout or maxiter."""
        restart = workdir / "xtbopt.xyz"
        if task.kind != "optimize" or not restart.is_file():
            return None
        mol: Molecule = task.inputs["mol"]
        text = (workdir / STDOUT_NAME).read_text(encoding="utf-8", errors="replace")
        inputs = {
            **task.inputs,
            "mol": Molecule(read_xyz(restart), mol.charge, mol.multiplicity),
            "start": task.inputs.get("start") or self._geometry(workdir / INPUT_NAME, mol.xyz),
            "trajectory": (*task.inputs.get("trajectory", ()), *_cycle_energies(text)),
        }
        return replace(task, inputs=inputs)

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> Evidence | Failure:
        text = result.stdout.read_text(encoding="utf-8", errors="replace")
        text += result.stderr.read_text(encoding="utf-8", errors="replace")
        failure = _termination(workdir, result, text)
        if failure is not None:
            return failure
        mol: Molecule = task.inputs["mol"]
        level = observed_level(text)
        if level is None:
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "xtb setup block incomplete")
        problems = level_mismatches(task.inputs["method"], level, version_pin=task.version_pin)
        if (level.charge, level.multiplicity) != (mol.charge, mol.multiplicity):
            problems.append(f"charge/multiplicity: observed {level.charge}/{level.multiplicity}")
        if problems:
            return _fail(FailureKind.METHOD_MISMATCH, "; ".join(problems))
        energies = _TOTAL_ENERGY.findall(text)
        if not energies:
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "no TOTAL ENERGY in the xtb output")
        start = task.inputs.get("start") or self._geometry(workdir / INPUT_NAME, mol.xyz)
        base: dict[str, Any] = {
            "engine": XTBEngine.name, "task": _TASKS[task.kind], "level": level,
            "start": start, "energy_hartree": float(energies[-1]),
            "output": self.file_ref(result.stdout), "job_key": "",
        }
        if task.kind == "optimize":
            return self._optimized(task, workdir, text, base)
        if task.kind == "frequencies":
            return self._frequencies(workdir, mol, base)
        return Evidence(final=start, **base)

    def _optimized(self, task: Task, workdir: Path, text: str, base: dict[str, Any]
                   ) -> Evidence | Failure:
        path = workdir / "xtbopt.xyz"
        if not path.is_file():
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "xtbopt.xyz missing")
        final = read_xyz(path)
        if list(final.symbols) != list(task.inputs["mol"].xyz.symbols):
            return _fail(FailureKind.METHOD_MISMATCH, "atom order changed")
        trajectory = (*task.inputs.get("trajectory", ()), *_cycle_energies(text),
                      base["energy_hartree"])
        return Evidence(final=self._geometry(path, final), trajectory_energies_hartree=trajectory,
                        **base)

    def _frequencies(self, workdir: Path, mol: Molecule, base: dict[str, Any]
                     ) -> Evidence | Failure:
        hessian = read_hessian(workdir / "hessian", len(mol.xyz.symbols))
        if hessian is None:
            return _fail(FailureKind.INCOMPLETE_OUTPUT, "hessian file missing or malformed")
        freqs, modes, n_external = projected_frequencies(hessian, mol.xyz.symbols, mol.xyz.coords)
        npy = to_canonical_npy(hessian, workdir / "hessian.npy")
        imaginary = tuple(tuple(float(x) for x in mode)
                          for f, mode in zip(freqs, modes, strict=True) if f < 0)
        base |= {"n_external": n_external}  # 3 (an atom: no mode), 5 (linear) or 6
        return Evidence(final=base["start"], frequencies_cm1=tuple(float(f) for f in freqs),
                        imaginary_modes=imaginary, hessian=self.file_ref(npy), **base)

    def _geometry(self, path: Path, xyz: XYZ) -> Geometry:
        """Fingerprint of the coordinates in memory, so an opt's final geometry read back
        from its file and handed to a freq job keeps the same fingerprint."""
        return Geometry(file=self.file_ref(path), symbols=tuple(xyz.symbols),
                        fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords))


def _cycle_energies(text: str) -> tuple[float, ...]:
    """Energies of the optimization cycles; the first one is the starting structure's."""
    return tuple(float(e) for e in _CYCLE_ENERGY.findall(text))


class XTBEngine:
    """QMEngine on the xtb executable; every job goes through the JobRunner."""

    name: ClassVar[str] = "xtb"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self.jobs = jobs
        self.site = site
        self.adapter = _Adapter(site, jobs.store.file_ref)

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(executables=("xtb",), version_command=("xtb", "--version"))

    def supports(self, method: MethodSpec) -> bool:
        solvation_ok = method.solvation is None or method.solvation.startswith("alpb:")
        return method.kind == "xtb" and method.gfn is not None and solvation_ok

    def task(self, kind: str, mol: Molecule, method: MethodSpec) -> Task:
        return Task(
            engine=self.name, version_pin=self.site.version, kind=kind,
            key_payload={"method": method.signature(), "molecule": mol.fingerprint()},
            execution=self.site.execution, inputs={"mol": mol, "method": method},
        )

    def energy(self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None
               ) -> Evidence | Failure:
        return self.jobs.run(self.task("energy", mol, method), self.adapter, deadline=deadline)

    def optimize(self, mol: Molecule, method: MethodSpec, *,
                 init_hessian: Evidence | None = None, deadline: Deadline | None = None
                 ) -> Evidence | Failure:
        return self.jobs.run(self.task("optimize", mol, method), self.adapter, deadline=deadline)

    def frequencies(self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None
                    ) -> Evidence | Failure:
        task = self.task("frequencies", mol, method)
        return self.jobs.run(task, self.adapter, deadline=deadline)
