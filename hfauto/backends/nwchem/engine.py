"""NWChem engines (design §6.3): QM (energy / optimize / frequencies, plus CCSD(T) energies),
the ZTS string (PATH) and the saddle refiner (SADDLE); gas phase only.

Each engine is also the JobRunner Adapter of its own jobs (render -> run -> parse). An
Evidence or PathProfile is returned only when the job terminated normally, the observed
Level matches the request (version pin, charge and multiplicity included), the echoed
input geometry equals the input within 1e-4 Å and a driver job's last DFT gradient is at
its final frame; anything else is a Failure.
"""

from __future__ import annotations

import json
import math
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any, ClassVar, Literal

import numpy as np
from pydantic import BaseModel

from hfauto.backends.nwchem import input as nw_in
from hfauto.backends.nwchem import output as nw_out
from hfauto.backends.protocols import Requirements
from hfauto.chemistry.elements import atomic_number
from hfauto.chemistry.gates import Policy, qrc_drop
from hfauto.chemistry.vibrations import projected_frequencies, shape_hessian
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    read_xyz,
    read_xyz_trajectory,
    write_xyz,
    written_geometry,
)
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import (
    Evidence,
    Failure,
    FailureKind,
    FileRef,
    Geometry,
    Level,
    PathProfile,
)
from hfauto.core.hashing import sha256_text
from hfauto.core.method import Deadline, EngineSite, MethodSpec, level_mismatches
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.process import STDOUT_NAME, Command, CommandResult, resolve_executable

NAME = "job"  # NWChem file prefix: job.nw, job.movecs, job.hess, job.drv.hess
FRAME_TOL_A = 1.0e-4
HESSIAN_NEAR_A = 0.5  # largest per-atom distance from a start to its init Hessian's structure
_SADDLE_CM1 = Policy().saddle_cm1  # a saddle direction, as the gates count them
_TASK: dict[str, Literal["sp", "opt", "freq", "saddle"]] = {
    "energy": "sp", "optimize": "opt", "frequencies": "freq", "saddle": "saddle"}
_DRIVER_JOBS = frozenset({"optimize", "saddle"})
# Failures a driver job continues from its latest frame, besides an opt's autoz failure; a
# saddle at maxiter is not continued (its updated Hessian stalls it): it returns its last frame
# for a restart with a fresh Hessian.
_CONTINUED = {"optimize": frozenset({FailureKind.TIMEOUT, FailureKind.GEOMETRY_MAXITER}),
              "saddle": frozenset({FailureKind.TIMEOUT})}
# double hybrids: a DFT task of NWChem silently leaves out their PT2 part
_DOUBLE_HYBRIDS = ("b2plyp", "b2gpplyp", "dsd-", "pwpb95")


def _invalid(reason: str) -> Failure:
    return Failure(kind=FailureKind.INPUT_INVALID, reason=reason)


def _incomplete(reason: str) -> Failure:
    return Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason=reason)


def _dft_supported(method: MethodSpec) -> bool:
    xc = (method.functional or "").lower()
    return method.kind == "dft" and bool(xc and method.basis) and not xc.startswith(_DOUBLE_HYBRIDS)


def _wft_supported(method: MethodSpec) -> bool:
    return method.kind == "wft" and method.wft_method is not None and bool(method.basis)


def _max_shift_A(xyz: XYZ, mol: Molecule) -> float:
    """Largest per-atom distance between two structures in one frame; inf for other atoms."""
    if list(xyz.symbols) != list(mol.xyz.symbols):
        return math.inf
    delta = np.asarray(xyz.coords, dtype=float) - np.asarray(mol.xyz.coords, dtype=float)
    return float(np.linalg.norm(delta.reshape(-1, 3), axis=1).max())


def _fixed_bond(mol: Molecule, bond: tuple[int, int, float] | None
                ) -> dict[str, tuple[int, int, float]] | Failure:
    """{"fixed_bond": (i < j, r rounded to 1e-4 Å as in the job key)} ({} without one), or a
    Failure unless mol has atoms i and j r apart within FRAME_TOL_A: NWChem moves an input to
    its zcoord value (the frame check)."""
    if bond is None:
        return {}
    coords = np.asarray(mol.xyz.coords, dtype=float).reshape(-1, 3)
    (i, j), r = sorted(map(int, bond[:2])), round(float(bond[2]), 4)
    if not 0 <= i < j < len(coords) or abs(np.linalg.norm(coords[j] - coords[i]) - r) > FRAME_TOL_A:
        return _invalid("fixed_bond_mismatch")
    return {"fixed_bond": (i, j, r)}


def _rescue_changed(text: str, level: Level) -> bool:
    """The plain SCF after cgmin (input._task) ended more than the SCF noise away from cgmin's
    energy: another solution than the one the job converged on."""
    energies = nw_out.dft_energies(text)
    return len(energies) < 2 or abs(energies[-1] - energies[-2]) > qrc_drop(level)


def _final_gradient(text: str, final: XYZ) -> tuple[float, ...] | None:
    """The last DFT gradient (3N, Eh/bohr), or None unless its block is at ``final`` within
    FRAME_TOL_A per atom."""
    block = nw_out.gradient(text)
    coords = np.asarray(final.coords, dtype=float).reshape(-1, 3)
    if block is None or block[0].shape != coords.shape:
        return None
    shift = np.linalg.norm(block[0] * BOHR_TO_ANGSTROM - coords, axis=1).max()
    return None if shift > FRAME_TOL_A else tuple(map(float, block[1].ravel()))


def _resumed(task: Task, workdir: Path, latest: Path | None, *, autoz: bool, rescue: bool
             ) -> Task:
    """The next attempt of ``task``: from ``latest`` (if any) without an initial Hessian, with
    the old vectors and the driver Hessian (not after autoz: it is in internal coordinates);
    Cartesian after autoz; cgmin after a DFT SCF failure."""
    kept = (".movecs",) if autoz else (".movecs", ".drv.hess")
    restart = [workdir / f"{NAME}{suffix}" for suffix in kept]
    inputs = {**task.inputs, "restart": [p for p in restart if p.is_file()],
              "cartesian": autoz or bool(task.inputs.get("cartesian")),
              "scf_rescue": bool(task.inputs.get("scf_rescue"))
              or (rescue and task.inputs["method"].kind == "dft")}
    if latest is not None:
        mol = task.inputs["mol"]
        inputs |= {"mol": Molecule(read_xyz(latest), mol.charge, mol.multiplicity),
                   "hessian": None}
    return replace(task, inputs=inputs)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""


def _render(task: Task, setup: nw_in.Setup) -> str:
    i = task.inputs
    mol, method, hessian = i["mol"], i["method"], i.get("hessian") is not None
    if task.kind == "energy":
        render = nw_in.render_wft if method.kind == "wft" else nw_in.render_energy
        return render(mol, method, setup)
    if task.kind == "optimize":
        return nw_in.render_optimize(mol, method, setup, init_hessian=hessian,
                                     fixed_bond=i.get("fixed_bond"))
    if task.kind == "frequencies":
        return nw_in.render_frequencies(mol, method, setup)
    if task.kind == "saddle":
        return nw_in.render_saddle(mol, method, setup, init_hessian=hessian)
    return nw_in.render_string(mol, i["end"], method, setup, nbeads=i["images"],
                               initial_path=i.get("initial_path") is not None)


class _NWChem:
    """Shared engine and Adapter behaviour; ``result_type`` is Evidence or PathProfile."""

    name: ClassVar[str]
    result_type: type[BaseModel] = Evidence  # the Adapter's result model

    def __init__(self, *, jobs: JobRunner, site: EngineSite) -> None:
        self._jobs, self._site = jobs, site

    @classmethod
    def requirements(cls) -> Requirements:
        return Requirements(executables=("nwchem",))  # NWChem has no version flag

    def supports(self, method: MethodSpec) -> bool:
        return _dft_supported(method)

    # --- engine side -------------------------------------------------------------------

    def _run(self, kind: str, payload: Mapping[str, Any], inputs: Mapping[str, Any],
             deadline: Deadline | None) -> Any:
        try:  # no all-electron deck for an element that needs an ECP
            nw_in.ecp(inputs["mol"].xyz.symbols, inputs["method"].basis)
        except ValueError as exc:
            return _invalid(str(exc))
        task = Task(engine=self.name, version_pin=self._site.version, kind=kind,
                    key_payload={"method": inputs["method"].signature(), **payload},
                    execution=self._site.execution, inputs=inputs)
        return self._jobs.run(task, self, deadline=deadline)

    def _hessian_file(self, freq: Evidence, mol: Molecule) -> Path | Failure:
        """The canonical .npy of a freq Evidence (any Level) computed at the same atoms in the
        same frame within HESSIAN_NEAR_A per atom of ``mol``."""
        store = self._jobs.store
        if (freq.task != "freq" or freq.hessian is None
                or _max_shift_A(read_xyz(store.resolve(freq.final.file)), mol) > HESSIAN_NEAR_A):
            return _invalid("hessian_geometry_mismatch")
        return store.resolve(freq.hessian)

    # --- Adapter side ------------------------------------------------------------------

    def prepare(self, task: Task, workdir: Path) -> Command:
        for source in task.inputs.get("restart", ()):
            shutil.copyfile(source, workdir / Path(source).name)
        if task.inputs.get("hessian") is not None:
            h = np.load(task.inputs["hessian"])
            if task.kind == "saddle":
                h = shape_hessian(h, task.inputs["mol"].xyz.coords, task.inputs["mode"])
            elif task.inputs.get("hessian_model") == "positive":
                h = shape_hessian(h, task.inputs["mol"].xyz.coords)
            (workdir / f"{NAME}.hess").write_text(nw_in.hess_text(h), encoding="ascii")
        if task.inputs.get("initial_path") is not None:
            shutil.copyfile(task.inputs["initial_path"], workdir / nw_in.INITIAL_PATH)
        (workdir / f"{NAME}.nw").write_text(_render(task, self._setup(task, workdir)),
                                            encoding="utf-8")
        env = {"OMP_NUM_THREADS": "1", **task.execution.env}
        return Command(argv=self._argv(task), cwd=workdir, env=env)

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> Any:
        if self._site.scratch_dir:
            shutil.rmtree(self._scratch(workdir), ignore_errors=True)
        text = _read(workdir / STDOUT_NAME)
        failure = nw_out.classify_failure(text, returncode=result.returncode,
                                          timed_out=result.timed_out)
        if failure is not None:
            stalled = task.kind == "saddle" and failure.kind is FailureKind.GEOMETRY_MAXITER
            latest = nw_out.final_xyz(workdir, task.inputs["mol"].xyz.symbols) if stalled else None
            return failure if latest is None else failure.model_copy(update={
                "final": self._geometry(latest), "energy_hartree": nw_out.step_energy(text)})
        level = self._observed(task, text)
        if isinstance(level, Failure):
            return level
        if task.kind == "string":
            return self._path(task, workdir, text, level)
        return self._evidence(task, workdir, text, level)

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None:
        """A driver job continues from its latest frame with its old vectors (design §7.1):
        after a timeout (and an opt's maxiter) with its driver Hessian; after an opt's autoz
        failure in Cartesian coordinates without it (a fixed bond becomes a spring restraint,
        input._fixed). An autoz failure before the first frame, or of a saddle, restarts in
        Cartesian coordinates from the same start. An SCF failure continues with the old
        vectors (DFT: cgmin, then one plain SCF for <S2>, input._task)."""
        autoz = failure.kind is FailureKind.INPUT_INVALID and failure.reason == "autoz"
        rescue = failure.kind is FailureKind.SCF_NOT_CONVERGED
        symbols = task.inputs["mol"].xyz.symbols
        latest = nw_out.final_xyz(workdir, symbols) if task.kind in _DRIVER_JOBS else None
        if autoz and (latest is None or task.kind == "saddle"):
            return replace(task, inputs={**task.inputs, "cartesian": True})
        continued = autoz or failure.kind in _CONTINUED.get(task.kind, ())
        if not (rescue or (continued and latest is not None)):
            return None
        return _resumed(task, workdir, latest, autoz=autoz, rescue=rescue)

    def _scratch(self, workdir: Path) -> Path:
        return Path(self._site.scratch_dir or ".") / f"{workdir.parent.name[:16]}_{workdir.name}"

    def _setup(self, task: Task, workdir: Path) -> nw_in.Setup:
        scratch = None
        if self._site.scratch_dir:
            self._scratch(workdir).mkdir(parents=True, exist_ok=True)
            scratch = self._scratch(workdir).as_posix()
        restart = task.inputs.get("restart", ())
        return nw_in.Setup(
            name=NAME, scratch_dir=scratch, memory_mb=task.execution.memory_mb_per_rank,
            cartesian=bool(task.inputs.get("cartesian")),
            restart_vectors=any(Path(p).suffix == ".movecs" for p in restart),
            scf_rescue=bool(task.inputs.get("scf_rescue")),
        )

    def _argv(self, task: Task) -> tuple[str, ...]:
        """nwchem <deck>, under ``mpirun -np <ranks> --bind-to none`` only when the site names
        mpirun. Never bound: every mpirun binds its ranks from core 0, so jobs run at once (a
        rank share of jobs.thread_map, or full-rank jobs of a site with ranks < cores) would
        share cores."""
        exe = self._site.executables
        nwchem = resolve_executable("nwchem", exe.get("nwchem")) or exe.get("nwchem", "nwchem")
        if "mpirun" not in exe:
            return (nwchem, f"{NAME}.nw")
        mpirun = resolve_executable("mpirun", exe["mpirun"]) or exe["mpirun"]
        return (mpirun, "-np", str(task.execution.ranks), "--bind-to", "none", nwchem,
                f"{NAME}.nw")

    def _observed(self, task: Task, text: str) -> Level | Failure:
        level = nw_out.observe_level(text)
        if level is None:
            return _incomplete("level_not_observed")
        mol, method = task.inputs["mol"], task.inputs["method"]
        problems = level_mismatches(method, level, version_pin=task.version_pin)
        problems += [f"{field}: requested {want}, observed {got}" for field, want, got in (
            ("charge", mol.charge, level.charge),
            ("multiplicity", mol.multiplicity, level.multiplicity)) if want != got]
        if problems:
            return Failure(kind=FailureKind.METHOD_MISMATCH, reason="; ".join(problems))
        shift = nw_out.frame_shift(text, mol.xyz.symbols, mol.xyz.coords)
        if shift is None:
            return _incomplete("input_geometry_not_echoed")
        if shift > FRAME_TOL_A:
            return Failure(kind=FailureKind.METHOD_MISMATCH, reason="frame_changed")
        return level

    def _geometry(self, path: Path) -> Geometry:
        return written_geometry(path, self._jobs.store.file_ref)

    def _evidence(self, task: Task, workdir: Path, text: str, level: Level) -> Evidence | Failure:
        energy, start_mol, s2 = nw_out.total_energy(text), task.inputs["start"], nw_out.s2(text)
        if energy is None:
            return _incomplete("no_energy")
        if s2 is None and level.multiplicity > 1 and task.inputs["method"].kind == "dft":
            return _incomplete("s2_not_reported")  # spin_ok would pass a missing <S2>
        if task.inputs.get("scf_rescue") and _rescue_changed(text, level):
            return _incomplete("rescue_solution_changed")
        start = self._geometry(write_xyz(start_mol.xyz, workdir / "start.xyz"))
        final, extra = start, dict[str, Any]()
        if task.kind in _DRIVER_JOBS:
            path = nw_out.final_xyz(workdir, start.symbols)
            if path is None:
                return _incomplete("final_xyz_missing_or_atom_order_changed")
            gradient = _final_gradient(text, read_xyz(path))
            if gradient is None:
                return _incomplete("gradient_not_at_final")
            final, extra["gradient"] = self._geometry(path), gradient
        elif task.kind == "frequencies":
            vibrations = self._frequencies(workdir, text, start_mol)
            if isinstance(vibrations, Failure):
                return vibrations
            extra.update(vibrations)
        return Evidence(engine=task.engine, task=_TASK[task.kind], level=level, start=start,
                        final=final, energy_hartree=energy, s2=s2,
                        output=self._jobs.store.file_ref(workdir / STDOUT_NAME), job_key="",
                        **extra)

    def _frequencies(self, workdir: Path, text: str, mol: Molecule) -> dict[str, Any] | Failure:
        """Projected frequencies from the job's own .hess (one vibrational block only)."""
        blocks = nw_out.count_frequency_blocks(text)
        if blocks != 1:
            return _incomplete(f"frequency_blocks:{blocks}")
        symbols = mol.xyz.symbols
        try:
            npy = nw_out.hess_to_npy(workdir / f"{NAME}.hess", len(symbols),
                                     workdir / "hessian.npy")
        except (OSError, ValueError) as exc:
            return _incomplete(f"hessian_unreadable:{type(exc).__name__}")
        freqs, modes, n_external = projected_frequencies(np.load(npy), symbols, mol.xyz.coords)
        imaginary = tuple(tuple(map(float, m)) for f, m in zip(freqs, modes, strict=True) if f < 0)
        return {"frequencies_cm1": tuple(map(float, freqs)), "n_external": n_external,
                "imaginary_modes": imaginary, "hessian": self._jobs.store.file_ref(npy)}

    def _path(self, task: Task, workdir: Path, text: str, level: Level) -> PathProfile | Failure:
        images, symbols = workdir / f"{NAME}.string_final.xyz", task.inputs["mol"].xyz.symbols
        energies = nw_out.string_energies(text)
        try:
            frames = read_xyz_trajectory(images)
        except (OSError, ValueError):
            return _incomplete("string_path_missing")
        n = task.inputs["images"]
        if len(frames) != n or len(energies) != n or any(f.symbols != symbols for f in frames):
            return _incomplete("string_path_images")
        return PathProfile(engine=task.engine, level=level,
                           images=self._jobs.store.file_ref(images), energies_hartree=energies)


class NWChemEngine(_NWChem):
    """QM: DFT energy / optimize (also with a fixed bond) / frequencies; CCSD(T) energies
    (ROHF-CCSD(T) for open shells)."""

    name: ClassVar[str] = "nwchem"

    def supports(self, method: MethodSpec) -> bool:
        return _dft_supported(method) or _wft_supported(method)

    def _qm(self, kind: str, mol: Molecule, method: MethodSpec, deadline: Deadline | None, *,
            payload: Mapping[str, Any] | None = None, **inputs: Any) -> Evidence | Failure:
        if not (_dft_supported(method) or (kind == "energy" and _wft_supported(method))):
            return _invalid(f"unsupported_method:{method.id}")
        return self._run(kind, {"molecule": mol.fingerprint(), **(payload or {})},
                         {"mol": mol, "start": mol, "method": method, **inputs}, deadline)

    def energy(self, mol: Molecule, method: MethodSpec, *,
               deadline: Deadline | None = None) -> Evidence | Failure:
        """A CCSD(T) key names the frozen core when an atom beyond Kr makes it differ from the
        old ``freeze atomic`` (input.frozen_core): no result of that deck is reused."""
        symbols = mol.xyz.symbols
        ecp = method.kind == "wft" and any(atomic_number(s) > 36 for s in symbols)
        payload = {"frozen_core": nw_in.frozen_core(symbols)} if ecp else {}
        return self._qm("energy", mol, method, deadline, payload=payload)

    def _guess(self, scf_guess: Evidence | None) -> tuple[dict[str, Any], dict[str, Any]]:
        """(key payload, inputs) starting the SCF from the converged vectors of ``scf_guess``'s
        job; empty when there is none (no job.movecs next to its output: another engine)."""
        if scf_guess is None:
            return {}, {}
        movecs = self._jobs.store.resolve(scf_guess.output).with_name(f"{NAME}.movecs")
        if not movecs.is_file():
            return {}, {}
        return {"scf_guess": scf_guess.job_key}, {"restart": [movecs]}

    def optimize(self, mol: Molecule, method: MethodSpec, *,
                 init_hessian: Evidence | None = None,
                 fixed_bond: tuple[int, int, float] | None = None,
                 scf_guess: Evidence | None = None,
                 deadline: Deadline | None = None) -> Evidence | Failure:
        """``init_hessian`` may come from a nearby structure (a QRC or mode-follow side from its
        saddle). A first-order saddle's (one mode below -saddle_cm1: the side's only negative
        direction is its displacement) is written as its positive-definite model
        (vibrations.shape_hessian; why in input.render_optimize). Any other is written as it
        is: from a higher-order saddle the side must stay free to leave its other saddle
        directions (DME C2v seed: 29 steps as it is, unconverged after 207 as the model).
        ``fixed_bond`` and ``scf_guess`` (a relaxed scan's point and its predecessor) enter the
        job key only when given, so every other key is unchanged."""
        hessian = None if init_hessian is None else self._hessian_file(init_hessian, mol)
        if isinstance(hessian, Failure):
            return hessian
        fixed = _fixed_bond(mol, fixed_bond)
        if isinstance(fixed, Failure):
            return fixed
        sha = init_hessian and init_hessian.hessian and init_hessian.hessian.sha256
        first_order = init_hessian is not None and hessian is not None and sum(
            f < -_SADDLE_CM1 for f in init_hessian.frequencies_cm1 or ()) == 1
        model = {"hessian_model": "positive"} if first_order else {}
        guess, restart = self._guess(scf_guess)
        return self._qm("optimize", mol, method, deadline,
                        payload={"hessian": sha, **model, **fixed, **guess},
                        hessian=hessian, **model, **fixed, **restart)

    def frequencies(self, mol: Molecule, method: MethodSpec, *, scf_guess: Evidence | None = None,
                    deadline: Deadline | None = None) -> Evidence | Failure:
        """``scf_guess`` (the opt or saddle at mol) starts the SCF from its converged vectors, so
        the freq stays on its electronic state: from scratch, the UKS OH···CH4 complex found the
        other OH π component, 7.2e-5 Eh above its opt."""
        guess, restart = self._guess(scf_guess)
        return self._qm("frequencies", mol, method, deadline, payload=guess, **restart)


class NWChemSaddle(_NWChem):
    """SADDLE: eigenvector following from a freq Hessian (any Level) computed at the seed or
    within HESSIAN_NEAR_A of it, written as its model with ``mode`` (the reaction direction,
    3N) its only negative curvature (vibrations.shape_hessian; the job key names the model). A
    search stopped at maxiter fails with its last frame and that frame's energy."""

    name: ClassVar[str] = "nwchem_saddle"

    def refine(self, seed: Molecule, method: MethodSpec, *, hessian: Evidence,
               mode: Sequence[float], deadline: Deadline | None = None) -> Evidence | Failure:
        if not _dft_supported(method):
            return _invalid(f"unsupported_method:{method.id}")
        path = self._hessian_file(hessian, seed)
        if isinstance(path, Failure):
            return path
        unit = np.ravel(mode) / np.linalg.norm(mode)
        payload = {"molecule": seed.fingerprint(), "hessian": hessian.hessian and
                   hessian.hessian.sha256, "hessian_model": "negative_along_mode",
                   "mode": sha256_text(json.dumps((np.round(unit, 6) + 0.0).tolist()))}
        inputs = {"mol": seed, "start": seed, "method": method, "hessian": path, "mode": unit}
        return self._run("saddle", payload, inputs, deadline)


class NWChemString(_NWChem):
    """PATH: one zero-temperature string chunk with frozen ends from the initial path;
    ``PathProfile.ts`` stays None (NWChem does not optimize a TS inside the string job)."""

    name: ClassVar[str] = "nwchem_string"
    result_type: type[BaseModel] = PathProfile

    def find_path(self, start: Molecule, end: Molecule, method: MethodSpec, *, images: int,
                  initial_path: FileRef, deadline: Deadline | None = None
                  ) -> PathProfile | Failure:
        if not _dft_supported(method):
            return _invalid(f"unsupported_method:{method.id}")
        if start.xyz.symbols != end.xyz.symbols or images < 3:
            return _invalid("string_endpoints_or_images")
        payload = {"start": start.fingerprint(), "end": end.fingerprint(), "images": images,
                   "initial_path": initial_path.sha256}
        inputs = {"mol": start, "start": start, "end": end, "method": method, "images": images,
                  "initial_path": self._jobs.store.resolve(initial_path)}
        return self._run("string", payload, inputs, deadline)
