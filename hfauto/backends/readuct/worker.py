"""One ReaDuct NT2 or AFIR attempt, run by ``hfauto.execution.worker`` in the attempt directory.

SCINE is imported only inside ``_Scine``. The decisions (imaginary-mode count on projected
frequencies, match with the source, choice of the IRC end, result assembly) are pure helpers.

Flow (design §6.3): the source is relaxed first (reference energy and structure; linear
sources arrive bent). NT2 → Bofill TS optimization with
``automatic_mode_selection = sorted(associations ∪ dissociations)`` → projected frequencies
(exactly one ν < −cutoff) → IRC → both ends optimized to minima → one end must have the
source's bond set, the other is the product. Bond sets are compared atom by atom in the source's
atom order, so a relabelled image of the source (a degenerate rearrangement such as a double H
exchange) is a product and a conformer or stereoisomer is not; the product keeps the source's
atom order. AFIR (γ, then the retry γ once) → unbiased optimization to a minimum. The IRC and
minimum optimizations take up to 500 iterations (C13). An SCC failure reruns the attempt once
at the retry electronic temperature; the energies of the source, TS and product are then
recomputed at the base temperature (chem 11).
"""

from __future__ import annotations

import importlib.metadata
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.backends.protocols import NO_NT2_MAXIMUM
from hfauto.chemistry.modes import amplitude, displace
from hfauto.chemistry.topology import bonds
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL

PRODUCT_FILE, TS_FILE = "product.xyz", "ts.xyz"
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


def matches_source(symbols: Sequence[str], source: np.ndarray, end: np.ndarray) -> bool:
    """Same atom-indexed bond set: conformers and stereoisomers of the source match, a
    relabelled image of it (degenerate rearrangement) does not."""
    return bonds(symbols, end) == bonds(symbols, source)


def irc_product(symbols: Sequence[str], source: np.ndarray,
                ends: Sequence[np.ndarray]) -> tuple[int | None, str | None]:
    """(index of the product end, None) or (None, negative reason)."""
    hits = [matches_source(symbols, source, end) for end in ends]
    if not any(hits):
        return None, "irc_not_connected_to_source"
    if all(hits):
        return None, "same_as_source"
    return hits.index(False), None


def afir_pair(trial: Mapping[str, Any]) -> tuple[tuple[int, int], bool] | None:
    """The biased pair: the first association (attractive), else the first dissociation."""
    for key, attractive in (("associations", True), ("dissociations", False)):
        if trial[key]:
            i, j = trial[key][0]
            return (int(i), int(j)), attractive
    return None


@dataclass(frozen=True)
class Found:
    outcome: str  # "product" | "negative"
    reason: str | None = None
    structures: dict[str, np.ndarray] = field(default_factory=dict)  # source / ts / product, Å
    energies: dict[str, float] = field(default_factory=dict)  # same keys, Eh
    ts_imag_cm1: float | None = None
    irc_connected: bool = False


def result_dict(found: Found, energies: Mapping[str, float], *, temperature_K: float,
                version: str) -> dict[str, Any]:
    """result.json of the attempt; barrier and reaction energy come from ``energies`` (a
    negative with a TS keeps its barrier)."""

    def relative(name: str) -> float | None:
        if name not in energies or "source" not in energies:
            return None
        return (energies[name] - energies["source"]) * HARTREE_TO_KCAL_MOL

    return {
        "version": version, "outcome": found.outcome, "reason": found.reason,
        "product": PRODUCT_FILE if found.outcome == "product" else None,
        "ts": TS_FILE if "ts" in found.structures else None,
        "ts_imag_cm1": found.ts_imag_cm1, "dE_act_kcal": relative("ts"),
        "dE_rxn_kcal": relative("product"), "irc_connected_to_source": found.irc_connected,
        "electronic_temperature_K": temperature_K,
    }


class _Scine:
    """SCINE calls of one attempt at one electronic temperature (systems keyed by name)."""

    def __init__(self, job: Mapping[str, Any], workdir: Path, temperature_K: float) -> None:
        import scine_readuct
        import scine_utilities
        import scine_xtb_wrapper  # noqa: F401  (registers the XTB calculator)

        self.readuct, self.utils = scine_readuct, scine_utilities
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
        return Found("negative", NO_NT2_MAXIMUM)
    if not run.task("run_tsopt_task", "guess", ["ts"], optimizer="bofill",
                    automatic_mode_selection=atoms):
        return Found("negative", "ts_not_converged")
    n_imag, imag, _ = run.imaginary("ts")
    if n_imag != 1:
        return Found("negative", f"ts_imaginary_modes:{n_imag}")
    ts, energy = {"ts": run.coords("ts")}, {"ts": run.energy("ts")}
    irc = run.task("run_irc_task", "ts", ["irc_f", "irc_b"], strict=False, **MAX_ITERATIONS)
    if not irc or not all(run.minimum(end, f"{end}_opt") for end in ("irc_f", "irc_b")):
        return Found("negative", "irc_end_not_minimum", ts, energy, ts_imag_cm1=imag)
    ends = [run.coords("irc_f_opt"), run.coords("irc_b_opt")]
    index, reason = irc_product(run.symbols, source, ends)
    if index is None:
        return Found("negative", reason, ts, energy, ts_imag_cm1=imag,
                     irc_connected=reason == "same_as_source")
    return Found("product", None, {**ts, "product": ends[index]},
                 {**energy, "product": run.energy(("irc_f_opt", "irc_b_opt")[index])},
                 ts_imag_cm1=imag, irc_connected=True)


def _afir(run: _Scine, source: np.ndarray, trial: Mapping[str, Any]) -> Found:
    chosen = afir_pair(trial)
    if chosen is None:
        return Found("negative", "afir_without_pair")
    (lhs, rhs), attractive = chosen
    reason = "afir_not_converged"
    settings = run.job["settings"]
    for gamma in (settings["afir_gamma_kj_mol"], settings["afir_gamma_retry_kj_mol"]):
        biased, relaxed = f"afir_{gamma:.0f}", f"afir_{gamma:.0f}_opt"
        if not run.task("run_afir_task", "start", [biased], afir_lhs_list=[lhs],
                        afir_rhs_list=[rhs], afir_attractive=attractive,
                        afir_energy_allowance=float(gamma)):
            reason = "afir_not_converged"
        elif not run.minimum(biased, relaxed):
            reason = "product_not_converged"
        elif matches_source(run.symbols, source, run.coords(relaxed)):
            reason = "same_as_source"
        else:
            return Found("product", None, {"product": run.coords(relaxed)},
                         {"product": run.energy(relaxed)})
    return Found("negative", reason)


def _explore(job: Mapping[str, Any], workdir: Path, temperature_K: float) -> Found:
    run = _Scine(job, workdir, temperature_K)
    run.load("start", np.asarray(job["coords"], dtype=float))
    if not run.task("run_opt_task", "start", ["source"]):
        return Found("negative", "source_not_converged")
    source = run.coords("source")
    explore = _nt2 if job["trial"]["mechanism"] == "nt2" else _afir
    found = explore(run, source, job["trial"])
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
    found = _explore(job, workdir, temperature_K)
    base_K = job["settings"]["electronic_temperature_K"]
    energies = found.energies
    if temperature_K != base_K and len(energies) > 1:  # a TS or a product besides the source
        energies = _single_points(job, workdir, found, base_K)
    for name, file in (("ts", TS_FILE), ("product", PRODUCT_FILE)):
        if name in found.structures:
            write_xyz(XYZ(list(job["symbols"]), found.structures[name]), workdir / file)
    return result_dict(found, energies, temperature_K=temperature_K, version=version)


def run_attempt(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """Worker entry: the attempt at the base temperature, once more at the retry one on SCC
    failure; an SCC failure there too is a ``scf_not_converged`` failure."""
    version = importlib.metadata.version("scine-readuct")
    settings = job["settings"]
    try:
        return _attempt(job, workdir, settings["electronic_temperature_K"], version)
    except SccFailure:
        pass
    try:
        return _attempt(job, workdir, settings["scc_retry_temperature_K"], version)
    except SccFailure as exc:
        return {"version": version,
                "failure": {"kind": "scf_not_converged", "reason": str(exc)[-300:]}}
