"""CREST conformer engine (design §6.3): ConformerSettings and Molecule → CLI → ConformerEnsemble.

One JobRunner task per search on the site's execution threads (``-T``, ``OMP_NUM_THREADS``);
``--nci`` searches add ``--noopt``. The attempt directory receives ``input.xyz``; the parser
reads ``crest_conformers.xyz`` (energies from the comment lines, a missing one stays None,
BUG-08), sums the CREGEN topology-based removals, and turns a stop on
'Change in topology detected' into ``topology_stops`` (the last ``crestopt.log`` frame).
The observed version (banner of every run, the same text as ``crest --version``) must match
the site pin. There is no RDKit fallback.

CREST runs in the attempt directory (under ``<run>/jobs``) without ``--scratch``: CREST 3.0.2
copies its working directory, captured stdout included, to the scratch directory and back,
which overwrites the log holding the CREGEN counts and the topology-stop message.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

from hfauto.backends.protocols import ConformerEnsemble, ConformerSettings, Requirements
from hfauto.chemistry.xyz import XYZ, Molecule, geometry_fingerprint, read_xyz, write_xyz
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.evidence import Failure, FailureKind, FileRef, Geometry
from hfauto.core.method import Deadline, EngineSite, MethodSpec
from hfauto.execution.jobs import Task
from hfauto.execution.process import Command, CommandResult, resolve_executable

if TYPE_CHECKING:
    from hfauto.execution.jobs import JobRunner
    from hfauto.execution.jobstore import JobStore

INPUT = "input.xyz"
CONFORMERS = "crest_conformers.xyz"
OPT_LOG = "crestopt.log"
TOPOLOGY_STOP = "Change in topology detected"
_VERSION = re.compile(r"Version\s+(\d[\w.+-]*)")
_REMOVALS = re.compile(r"CREGEN> number of topology-based structure removals:\s*(\d+)")
_FLOAT = re.compile(r"[-+]?\d+\.\d*(?:[Ee][-+]?\d+)?")


def command(mol: Molecule, method: MethodSpec, settings: ConformerSettings, *, threads: int,
            executable: str = "crest") -> tuple[str, ...]:
    """``--gfnN [--nci] [--quick] -T n --ewin e --chrg q --uhf m-1 [--notopo atoms] [--noopt]
    [--alpb solvent]``; atoms are 1-based on the CLI."""
    argv = [executable, INPUT, f"--gfn{method.gfn}"]
    argv += ["--nci"] if settings.nci else []
    argv += ["--quick"] if settings.quick else []
    argv += ["-T", str(threads), "--ewin", f"{settings.ewin_kcal:g}",
             "--chrg", str(mol.charge), "--uhf", str(mol.multiplicity - 1)]
    if settings.notopo_atoms:
        argv += ["--notopo", ",".join(str(i + 1) for i in settings.notopo_atoms)]
    # CREST 3.0.2's initial topology check after the pre-optimization ignores --notopo
    # (setuptest.f90), so acid-base and anion complexes stopped there; --noopt skips it.
    argv += ["--noopt"] if settings.nci else []
    if method.solvation:
        argv += ["--alpb", method.solvation.split(":", 1)[1]]
    return tuple(argv)


def observed_version(stdout: str) -> str | None:
    match = _VERSION.search(stdout)
    return match.group(1) if match else None


def topology_removed(stdout: str) -> int:
    """Sum over every CREGEN call of the topology-based structure removals."""
    return sum(int(n) for n in _REMOVALS.findall(stdout))


def _energy(comment: str) -> float | None:
    match = _FLOAT.search(comment)
    return float(match.group(0)) if match else None


def _geometry(xyz: XYZ, path: Path, file_ref: Callable[[Path], FileRef]) -> Geometry:
    written = read_xyz(write_xyz(XYZ(xyz.symbols, xyz.coords, "generated_by=crest"), path))
    return Geometry(file=file_ref(path), symbols=tuple(written.symbols),
                    fingerprint=geometry_fingerprint(written.symbols, written.coords))


def parse_outputs(workdir: Path, stdout: str, symbols: Sequence[str],
                  file_ref: Callable[[Path], FileRef]) -> ConformerEnsemble | Failure:
    """Ensemble of a finished CREST directory (return code already checked)."""
    removed, version = topology_removed(stdout), observed_version(stdout) or "unknown"
    if TOPOLOGY_STOP in stdout:
        if not (workdir / OPT_LOG).is_file():
            reason = f"topology stop without {OPT_LOG}"
            return Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason=reason)
        last = read_xyz_trajectory(workdir / OPT_LOG)[-1]
        stop = _geometry(last, workdir / "topology_stop.xyz", file_ref)
        return ConformerEnsemble(members=(), topology_removed=removed, topology_stops=(stop,),
                                 version=version, job_key="")
    if not (workdir / CONFORMERS).is_file():
        return Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason=f"{CONFORMERS} missing")
    frames = read_xyz_trajectory(workdir / CONFORMERS)
    if any(list(f.symbols) != list(symbols) for f in frames):
        return Failure(kind=FailureKind.METHOD_MISMATCH, reason="atom_order_changed")
    members = tuple((_geometry(f, workdir / f"conformer_{k:03d}.xyz", file_ref), _energy(f.comment))
                    for k, f in enumerate(frames))
    return ConformerEnsemble(members=members, topology_removed=removed, topology_stops=(),
                             version=version, job_key="")


class _Adapter:
    result_type = ConformerEnsemble

    def __init__(self, store: JobStore, site: EngineSite) -> None:
        self.store, self.site = store, site

    def prepare(self, task: Task, workdir: Path) -> Command:
        mol: Molecule = task.inputs["molecule"]
        mol.write(workdir / INPUT)
        exe = resolve_executable("crest", self.site.executables.get("crest")) or "crest"
        threads = task.execution.threads
        env = {"OMP_NUM_THREADS": f"{threads},1", "OMP_STACKSIZE": "4G", **task.execution.env}
        argv = command(mol, task.inputs["method"], task.inputs["settings"], threads=threads,
                       executable=exe)
        return Command(argv=argv, cwd=workdir, env=env)

    def parse(self, task: Task, workdir: Path, result: CommandResult
              ) -> ConformerEnsemble | Failure:
        if result.timed_out:
            return Failure(kind=FailureKind.TIMEOUT, reason=f"{result.duration_s:.0f} s")
        stdout = result.stdout.read_text(encoding="utf-8", errors="replace")
        version = observed_version(stdout)
        if (version or "").lower() != self.site.version.lower():
            return Failure(kind=FailureKind.METHOD_MISMATCH,
                           reason=f"version: requested {self.site.version!r}, observed {version!r}")
        if result.returncode != 0 and TOPOLOGY_STOP not in stdout:
            return Failure(kind=FailureKind.NONZERO_EXIT, reason=f"returncode {result.returncode}")
        mol: Molecule = task.inputs["molecule"]
        return parse_outputs(workdir, stdout, mol.xyz.symbols, self.store.file_ref)

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> None:
        return None


class CRESTEngine:
    """ConformerEngine on the ``crest`` executable (GFN-xTB methods, ALPB solvation only)."""

    name: ClassVar[str] = "crest"

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self.jobs, self.site = jobs, site

    @property
    def _adapter(self) -> _Adapter:
        return _Adapter(self.jobs.store, self.site)

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(executables=("crest",), version_command=("crest", "--version"))

    def supports(self, method: MethodSpec) -> bool:
        solvation_ok = method.solvation is None or method.solvation.startswith("alpb:")
        etemp_ok = method.electronic_temperature_K in (None, 300.0)
        return method.kind == "xtb" and method.gfn is not None and solvation_ok and etemp_ok

    def search(self, mol: Molecule, method: MethodSpec, settings: ConformerSettings, *,
               deadline: Deadline | None = None) -> ConformerEnsemble | Failure:
        if not self.supports(method):
            return Failure(kind=FailureKind.INPUT_INVALID, reason=f"crest cannot run {method.id!r}")
        task = Task(
            engine=self.name, version_pin=self.site.version, kind="conformers",
            key_payload={"molecule": mol.fingerprint(), "method": method.signature(),
                         "settings": settings.model_dump(mode="json"), "noopt": settings.nci},
            execution=self.site.execution,
            inputs={"molecule": mol, "method": method, "settings": settings},
        )
        return self.jobs.run(task, self._adapter, deadline=deadline)
