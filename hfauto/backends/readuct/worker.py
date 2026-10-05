"""One ReaDuct NT2 attempt, run by ``hfauto.execution.worker`` in the attempt directory.

SCINE is imported only inside ``_Scine``. The decisions (imaginary-mode count on projected
frequencies, the order of the IRC ends, result assembly) are pure helpers.

Flow (design §6.3): the start is relaxed first (reference energy and structure; linear starts
arrive bent); a start that relaxed into other bonds is an edge without a TS, from the start to
that minimum. NT2 → Bofill TS optimization with
``automatic_mode_selection = sorted(associations ∪ dissociations)`` → projected frequencies
(exactly one ν < −cutoff) → IRC → both ends optimized to minima. The two ends are an edge unless
they have the same bonds (topology.same_bonding, atom by atom: a relabelled image of the source,
a degenerate rearrangement's product, is another structure); the end with the source's bonds
comes first, and the caller identifies both as states. The IRC and minimum optimizations take up
to 500 iterations. An SCC failure reruns the attempt once at the retry electronic temperature;
the energies are then recomputed at the base temperature.
"""

from __future__ import annotations

import importlib.metadata
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.modes import amplitude, displace
from hfauto.chemistry.topology import same_bonding
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL

TS_FILE, END_FILES = "ts.xyz", ("end0.xyz", "end1.xyz")
# the distributions whose versions are the engine's pin (EngineSite.version), in this order
DISTRIBUTIONS = ("scine-readuct", "scine-utilities", "scine-xtb-wrapper")
# scine-xtb-wrapper 3.0.2 accepts only "any" and "restricted_open_shell"; the latter is xTB's
# own treatment (unpaired electrons from the multiplicity), so it is set explicitly.
SPIN_MODE = "restricted_open_shell"
MAX_ITERATIONS = {"convergence_max_iterations": 500}  # IRC and minimum optimizations


class SccFailure(Exception):
    """The xTB SCC did not converge."""


def is_scc_failure(message: str) -> bool:
    return "self consistent charge" in message.lower()


def imaginary_count(freqs_cm1: Sequence[float] | np.ndarray, cutoff_cm1: float) -> int:
    return int(np.count_nonzero(np.asarray(freqs_cm1) < -cutoff_cm1))


def irc_ends(symbols: Sequence[str], source: np.ndarray, ends: Sequence[np.ndarray]
             ) -> tuple[tuple[int, int], bool] | None:
    """(the order of the two ends, the one with the source's bonds first; whether one has
    them), or None when the ends have the same bonds (no edge: a conformer change)."""
    a, b = ends
    if same_bonding(symbols, a, b):
        return None
    if same_bonding(symbols, source, b):
        return (1, 0), True
    return (0, 1), same_bonding(symbols, source, a)


@dataclass(frozen=True)
class Found:
    outcome: str  # "product" | "negative"
    reason: str | None = None
    structures: dict[str, np.ndarray] = field(default_factory=dict)  # source/ts/end0/end1, Å
    energies: dict[str, float] = field(default_factory=dict)  # same keys, Eh
    irc_connected: bool = False  # end0 has the source's bonds


def version() -> str:
    """The engine's version pin as installed: each of DISTRIBUTIONS as ``name==version``."""
    return ",".join(f"{d}=={importlib.metadata.version(d)}" for d in DISTRIBUTIONS)


def result_dict(found: Found, energies: Mapping[str, float], *, version: str) -> dict[str, Any]:
    """result.json of the attempt: energies (kcal/mol) from the first end, or without ends
    from the relaxed start (a negative with a TS keeps its barrier)."""
    zero = energies.get("end0", energies.get("source"))

    def relative(name: str) -> float | None:
        known = name in energies and zero is not None
        return (energies[name] - zero) * HARTREE_TO_KCAL_MOL if known else None

    has_ends = found.outcome == "product"
    return {
        "version": version, "outcome": found.outcome, "reason": found.reason,
        "ends": list(END_FILES) if has_ends else None,
        "ts": TS_FILE if "ts" in found.structures else None,
        "dE_act_kcal": relative("ts"), "dE_rxn_kcal": relative("end1"),
        "irc_connected_to_source": found.irc_connected,
    }


class _Scine:
    """SCINE calls of one attempt at one electronic temperature (systems keyed by name)."""

    def __init__(self, job: Mapping[str, Any], workdir: Path, temperature_K: float) -> None:
        import scine_readuct
        import scine_utilities
        import scine_xtb_wrapper  # noqa: F401  (registers the XTB calculator)

        self.readuct = scine_readuct
        self.utils: Any = scine_utilities  # its type stub lacks the core module
        self.job, self.workdir, self.T = job, workdir, temperature_K
        self.symbols: list[str] = list(job["symbols"])
        self.systems: dict[str, Any] = {}

    def load(self, name: str, coords: np.ndarray) -> None:
        path = write_xyz(XYZ(self.symbols, np.asarray(coords)),
                         self.workdir / f"{name}_{self.T:.0f}K.xyz")
        self.systems[name] = self.utils.core.load_system_into_calculator(
            str(path), self.job["method_family"], program="XTB",
            molecular_charge=self.job["charge"], spin_multiplicity=self.job["multiplicity"],
            spin_mode=SPIN_MODE, max_scf_iterations=self.job["settings"]["max_scf_iterations"],
            electronic_temperature=self.T)

    def task(self, run: str, name: str, output: Sequence[str] = (), *, strict: bool = True,
             **settings: Any) -> bool:
        """Whether the task succeeded; SCC failures raise SccFailure.

        ``strict=False`` keeps what an unconverged task produced (the IRC ends are optimized
        afterwards anyway).
        """
        if output:
            settings["output"] = list(output)
        try:
            self.systems, ok = getattr(self.readuct, run)(
                self.systems, [name], stop_on_error=strict, **settings)
        except RuntimeError as exc:  # SCINE reports every failed task as RuntimeError
            if is_scc_failure(str(exc)):
                raise SccFailure(str(exc)) from exc
            print(f"{run}: {exc}", file=sys.stderr)  # the cause of the failure, for diagnosis
            return False
        return (bool(ok) or not strict) and all(o in self.systems for o in output)

    def coords(self, name: str) -> np.ndarray:
        return np.asarray(self.systems[name].positions, dtype=float) * BOHR_TO_ANGSTROM

    def energy(self, name: str) -> float:
        return float(self.systems[name].get_results().energy)

    def imaginary(self, name: str) -> tuple[int, float, np.ndarray]:
        """(ν < −cutoff count, lowest ν, its mode) from projected frequencies; −1 if it failed."""
        if not self.task("run_hessian_task", name):
            return -1, float("nan"), np.zeros(0)
        freqs, modes, _ = projected_frequencies(
            self.systems[name].get_results().hessian, self.symbols, self.coords(name))
        cutoff = self.job["settings"]["imag_cutoff_cm1"]
        return imaginary_count(freqs, cutoff), float(freqs[0]), modes[0]

    def minimum(self, name: str, output: str) -> bool:
        """Optimize ``name`` into the minimum ``output``; one that ends with ν < −cutoff is
        displaced ± along its lowest mode (modes.amplitude) and the lower side that optimizes
        to n_imag = 0 is taken."""
        if not self.task("run_opt_task", name, [output], **MAX_ITERATIONS):
            return False
        n_imag, nu, mode = self.imaginary(output)
        if n_imag <= 0:
            return n_imag == 0
        sides = []
        step = amplitude(nu, mode, self.symbols)
        for tag, coords in zip("pm", displace(self.coords(output), mode, step), strict=True):
            start, end = f"{output}_{tag}", f"{output}_{tag}_opt"
            self.load(start, coords)
            if (self.task("run_opt_task", start, [end], **MAX_ITERATIONS)
                    and self.imaginary(end)[0] == 0):
                sides.append(end)
        if sides:
            self.systems[output] = self.systems[min(sides, key=self.energy)]
        return bool(sides)


def _flat(pairs: Sequence[Sequence[int]]) -> list[int]:
    return [int(i) for pair in pairs for i in pair]


def _nt2(run: _Scine, source: np.ndarray, trial: Mapping[str, Any]) -> Found:
    atoms = sorted(set(_flat(trial["associations"]) + _flat(trial["dissociations"])))
    if not run.task("run_nt2_task", "start", ["guess"],
                    nt_associations=_flat(trial["associations"]),
                    nt_dissociations=_flat(trial["dissociations"])):
        return Found("negative", "no_nt2_maximum")
    if not run.task("run_tsopt_task", "guess", ["ts"], optimizer="bofill",
                    automatic_mode_selection=atoms):
        return Found("negative", "ts_not_converged")
    n_imag = run.imaginary("ts")[0]
    if n_imag != 1:
        return Found("negative", f"ts_imaginary_modes:{n_imag}")
    ts, energy = {"ts": run.coords("ts")}, {"ts": run.energy("ts")}
    irc = run.task("run_irc_task", "ts", ["irc_f", "irc_b"], strict=False, **MAX_ITERATIONS)
    if not irc or not all(run.minimum(end, f"{end}_opt") for end in ("irc_f", "irc_b")):
        return Found("negative", "irc_end_not_minimum", ts, energy)
    names = ("irc_f_opt", "irc_b_opt")
    order = irc_ends(run.symbols, source, [run.coords(name) for name in names])
    if order is None:
        return Found("negative", "no_bond_change", ts, energy)
    (i, j), connected = order
    ends = {"end0": names[i], "end1": names[j]}
    return Found("product", None, ts | {k: run.coords(n) for k, n in ends.items()},
                 energy | {k: run.energy(n) for k, n in ends.items()}, irc_connected=connected)


def _explore(job: Mapping[str, Any], workdir: Path, temperature_K: float) -> Found:
    run = _Scine(job, workdir, temperature_K)
    start = np.asarray(job["coords"], dtype=float)
    run.load("start", start)
    if not run.task("run_opt_task", "start", ["source"]):
        return Found("negative", "source_not_converged")
    source = run.coords("source")
    if not same_bonding(run.symbols, start, source):  # no basin of the start's bonds there
        return Found("product", None, {"end0": start, "end1": source})
    found = _nt2(run, source, job["trial"])
    return replace(found, structures={"source": source, **found.structures},
                   energies={"source": run.energy("source"), **found.energies})


def _single_points(job: Mapping[str, Any], workdir: Path, found: Found,
                   temperature_K: float) -> dict[str, float]:
    run = _Scine(job, workdir, temperature_K)
    energies: dict[str, float] = {}
    for name, coords in found.structures.items():
        run.load(f"{name}_sp", coords)
        if not run.task("run_sp_task", f"{name}_sp"):
            raise RuntimeError(f"{name}: single point at {temperature_K} K failed")
        energies[name] = run.energy(f"{name}_sp")
    return energies


def _attempt(job: Mapping[str, Any], workdir: Path, temperature_K: float,
             version: str) -> dict[str, Any]:
    found, base_K = _explore(job, workdir, temperature_K), job["settings"]["electronic_temperature_K"]
    energies = found.energies
    if temperature_K != base_K and len(energies) > 1:  # a TS or ends besides the source
        energies = _single_points(job, workdir, found, base_K)
    for name, file in (("ts", TS_FILE), *zip(("end0", "end1"), END_FILES, strict=True)):
        if name in found.structures:
            write_xyz(XYZ(list(job["symbols"]), found.structures[name]), workdir / file)
    return result_dict(found, energies, version=version)


def run_attempt(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """Worker entry: the attempt at the base temperature, once more at the retry one on SCC
    failure; an SCC failure there too is a ``scf_not_converged`` failure."""
    pin, settings = version(), job["settings"]
    try:
        return _attempt(job, workdir, settings["electronic_temperature_K"], pin)
    except SccFailure:
        pass
    try:
        return _attempt(job, workdir, settings["scc_retry_temperature_K"], pin)
    except SccFailure as exc:
        failure = {"kind": "scf_not_converged", "reason": str(exc)[-300:]}
        return {"version": pin, "failure": failure}
