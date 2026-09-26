"""NWChem input decks (design §6.3, §8.2): pure functions from typed requests to text.

Every deck writes the top-level ``charge``, the geometry in the input frame
(``units angstrom nocenter noautosym``, plus ``noautoz`` for Cartesian coordinates), a
spherical basis and, for DFT, xc / mult (``odft`` when open shell) / grid / energy
convergence / dispersion. Optimizations and saddles never compute a Hessian: frequencies
are a job of their own, and an initial Hessian is read from ``<name>.hess`` (``inhess 2``).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec

OPT_MAXITER = 100
SADDLE_MAXITER = 50
STRING_MAXITER = 20
INITIAL_PATH = "initial_path.xyz"  # xyz_path is read from permanent_dir (the job's cwd)
_VDW = {"d3zero": 3, "d3bj": 4}
# NWChem tightens frequency jobs to 1e-7 on its own (opt default 1e-6); writing one value for
# every task keeps opt, saddle, freq and sp on the same numerics layer (same_pes).
DEFAULT_SCF_ENERGY_TOL = 1e-7


@dataclass(frozen=True)
class Setup:
    """Per-attempt settings that are not chemistry (never part of the job key)."""

    name: str = "job"  # file prefix: <name>.movecs, <name>.hess, <name>.drv.hess
    scratch_dir: str | None = None  # EngineSite.scratch_dir; permanent_dir is the cwd
    memory_mb: int = 1200  # per rank
    cartesian: bool = False  # ExecutionSpec.coordinates == "cartesian" (noautoz)
    maxiter: int | None = None  # ExecutionSpec.maxiter overrides the driver / string default
    restart_vectors: bool = False  # start from <name>.movecs of the previous attempt
    scf_rescue: bool = False  # damping and level shift after SCF_NOT_CONVERGED


_DEFAULT = Setup()


def cosmo_dielectric(method: MethodSpec) -> float | None:
    """The dielectric constant of ``solvation: cosmo:<eps>``; other solvation is unsupported."""
    if method.solvation is None:
        return None
    model, _, value = method.solvation.partition(":")
    if model.lower() != "cosmo":
        raise ValueError(f"NWChem supports only cosmo:<eps> solvation, not {method.solvation!r}")
    return float(value)


def _deck(setup: Setup, *blocks: Sequence[str]) -> str:
    header = ["echo", f"start {setup.name}", "permanent_dir ."]
    if setup.scratch_dir:
        header.append(f"scratch_dir {setup.scratch_dir}")
    header.append(f"memory total {setup.memory_mb} mb")
    return "\n\n".join("\n".join(block) for block in (header, *blocks) if block) + "\n"


def _geometry(xyz: XYZ, *, cartesian: bool, label: str = "") -> list[str]:
    flags = "units angstrom nocenter noautosym" + (" noautoz" if cartesian else "")
    rows = np.asarray(xyz.coords, dtype=float).reshape(-1, 3)
    lines = [f"geometry {label} {flags}" if label else f"geometry {flags}"]
    lines += [f"  {s:2s} {x: .10f} {y: .10f} {z: .10f}" for s, (x, y, z) in
              zip(xyz.symbols, rows, strict=True)]
    return [*lines, "end"]


def _system(mol: Molecule, method: MethodSpec, setup: Setup, *, cartesian: bool = False,
            end: Molecule | None = None) -> list[str]:
    """Geometry (and the string's end geometry), charge, basis and COSMO."""
    if not method.basis:
        raise ValueError(f"method {method.id!r} names no basis set")
    noautoz = setup.cartesian or cartesian
    lines = _geometry(mol.xyz, cartesian=noautoz)
    if end is not None:
        lines += _geometry(end.xyz, cartesian=noautoz, label="endgeom")
    lines += [f"charge {mol.charge}", "basis spherical", f"  * library {method.basis}", "end"]
    eps = cosmo_dielectric(method)
    return lines + ([] if eps is None else ["cosmo", f"  dielec {eps:g}", "end"])


def _dft(mol: Molecule, method: MethodSpec, setup: Setup) -> list[str]:
    if method.kind != "dft" or not method.functional:
        raise ValueError(f"method {method.id!r} is not a DFT method")
    if method.dispersion is not None and method.dispersion not in _VDW:
        raise ValueError(f"NWChem has no {method.dispersion} dispersion")
    lines = ["dft", f"  xc {method.functional.lower()}", f"  mult {mol.multiplicity}"]
    if mol.multiplicity > 1:
        lines.append("  odft")
    if method.grid:
        lines.append(f"  grid {method.grid}")
    lines.append(f"  convergence energy {method.scf_energy_tol or DEFAULT_SCF_ENERGY_TOL:.1e}")
    if method.dispersion:
        lines.append(f"  disp vdw {_VDW[method.dispersion]}")
    if setup.restart_vectors:
        lines.append(f"  vectors input {setup.name}.movecs")
    if setup.scf_rescue:
        lines.append("  convergence damp 40 ncydp 30 lshift 0.5")
    return [*lines, "end"]


def _driver(setup: Setup, maxiter: int, options: Sequence[str]) -> list[str]:
    return ["driver", f"  maxiter {setup.maxiter or maxiter}", *options, "  xyz final", "end"]


def render_energy(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 ["task dft energy"])


def render_optimize(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT, *,
                    init_hessian: bool = False) -> str:
    """Default driver thresholds; ``init_hessian`` reads <name>.hess.

    ``trust 0.3`` (the NWChem default) with an initial Hessian (a QRC side from its TS, a
    complex from xTB), ``trust 0.1`` without: from a mode-follow displacement the diagonal
    guess overshoots back above the TS energy (HCN->HNC in NWChem 7.2.3).
    """
    options = ["  trust 0.3", "  inhess 2"] if init_hessian else ["  trust 0.1"]
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 _driver(setup, OPT_MAXITER, options), ["task dft optimize"])


def render_frequencies(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """The only deck with a Hessian task; <name>.hess is written to permanent_dir."""
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 ["task dft frequencies"])


def render_saddle(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT, *,
                  moddir: int = 0, cartesian: bool = False, init_hessian: bool = True) -> str:
    """Eigenvector following from <name>.hess along the driver's mode ``moddir`` (1-based;
    0 leaves the choice to NWChem); ``cartesian`` when that number counts Cartesian modes."""
    options = ["  trust 0.1", "  sadstp 0.1"] + (["  inhess 2"] if init_hessian else [])
    options += [f"  moddir {moddir}"] if moddir > 0 else []
    return _deck(setup, _system(mol, method, setup, cartesian=cartesian), _dft(mol, method, setup),
                 _driver(setup, SADDLE_MAXITER, options), ["task dft saddle"])


def render_string(start: Molecule, end: Molecule, method: MethodSpec, setup: Setup = _DEFAULT,
                  *, nbeads: int, initial_path: bool = False) -> str:
    """Zero-temperature string with frozen ends; ``initial_path`` reads INITIAL_PATH."""
    string = ["string", f"  nbeads {nbeads}", f"  maxiter {setup.maxiter or STRING_MAXITER}",
              "  stepsize 0.05", "  interpol 3", "  tol 1e-5", "  freeze1 .true.",
              "  freezeN .true.", "  impose"]
    string += [f"  xyz_path {INITIAL_PATH}"] if initial_path else []
    return _deck(setup, _system(start, method, setup, end=end), _dft(start, method, setup),
                 [*string, "end"], ["task dft string"])


def render_wft(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """MP2 or CCSD(T) single point with frozen atomic cores; closed shell only (RHF reference).

    CCSD may take 50 iterations (default 20); the scf block only restarts from <name>.movecs.
    """
    if method.kind != "wft" or method.wft_method is None:
        raise ValueError(f"method {method.id!r} is not a wave-function method")
    scf = ["scf", f"  vectors input {setup.name}.movecs", "end"] if setup.restart_vectors else []
    module = "mp2" if method.wft_method == "mp2" else "ccsd"
    return _deck(setup, _system(mol, method, setup), scf,
                 [module, "  freeze atomic", *(["  maxiter 50"] if module == "ccsd" else []),
                  "end"], [f"task {method.wft_method} energy"])


def hess_text(hessian: np.ndarray) -> str:
    """NWChem ``.hess``: the lower triangle row by row in Eh/bohr², Fortran D exponents."""
    h = np.asarray(hessian, dtype=float)
    return "".join(f"{v:20.10E}\n".replace("E", "D") for v in h[np.tril_indices(len(h))])
