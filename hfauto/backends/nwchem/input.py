"""NWChem input decks (design §6.3, §8.2): pure functions from typed requests to text.

Every deck writes the top-level ``charge``, the geometry in the input frame
(``units angstrom nocenter noautosym``, plus ``noautoz`` for Cartesian coordinates after an
autoz failure), a spherical basis with the def2-ECP of every element beyond Kr, at most
SCF_MAXITER SCF cycles and, for DFT, xc / mult (``odft`` when open shell) / grid / energy
convergence / dispersion. Gas phase only.
Optimizations and saddles never compute a Hessian: frequencies are a job of their own, and an
initial Hessian is read from ``<name>.hess`` (``inhess 2``).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from hfauto.chemistry.elements import atomic_number
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec

OPT_MAXITER = 100
SADDLE_MAXITER = 50
STRING_MAXITER = 20
SCF_MAXITER = 100  # NWChem 7.2.3 defaults: DFT 50, SCF 20
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
    cartesian: bool = False  # noautoz: the continuation of an autoz failure
    restart_vectors: bool = False  # start from <name>.movecs of the previous attempt
    scf_rescue: bool = False  # quadratic SCF (cgmin) after SCF_NOT_CONVERGED


_DEFAULT = Setup()


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


def ecp(symbols: Sequence[str], basis: str) -> list[str]:
    """The ``ecp`` block: def2 defines Rb-Rn together with the def2-ECP (Weigend & Ahlrichs,
    PCCP 7, 3297 (2005)), one line per element (``*`` aborts NWChem 7.2.3 on H-Kr). Any other
    basis with such an element raises ValueError: an all-electron deck would be silently wrong.
    """
    heavy = sorted({s for s in symbols if atomic_number(s) > 36}, key=atomic_number)
    if not heavy:
        return []
    if not basis.lower().startswith("def2"):
        raise ValueError(f"no_ecp_for_basis:{basis}:{','.join(heavy)}")
    return ["ecp", *(f"  {s} library def2-ecp" for s in heavy), "end"]


def _system(mol: Molecule, method: MethodSpec, setup: Setup, *, end: Molecule | None = None
            ) -> list[str]:
    """Geometry (and the string's end geometry), charge, basis and ECP."""
    if not method.basis:
        raise ValueError(f"method {method.id!r} names no basis set")
    lines = _geometry(mol.xyz, cartesian=setup.cartesian)
    if end is not None:
        lines += _geometry(end.xyz, cartesian=setup.cartesian, label="endgeom")
    lines += [f"charge {mol.charge}", "basis spherical", f"  * library {method.basis}", "end"]
    return lines + ecp(mol.xyz.symbols, method.basis)


def _dft(mol: Molecule, method: MethodSpec, setup: Setup) -> list[str]:
    if method.kind != "dft" or not method.functional:
        raise ValueError(f"method {method.id!r} is not a DFT method")
    lines = ["dft", f"  xc {method.functional.lower()}", f"  mult {mol.multiplicity}"]
    if mol.multiplicity > 1:
        lines.append("  odft")
    if method.grid:
        lines.append(f"  grid {method.grid}")
    lines.append(f"  convergence energy {method.scf_energy_tol or DEFAULT_SCF_ENERGY_TOL:.1e}")
    lines.append(f"  iterations {SCF_MAXITER}")
    if method.dispersion:
        lines.append(f"  disp vdw {_VDW[method.dispersion]}")
    if setup.restart_vectors:
        lines.append(f"  vectors input {setup.name}.movecs")
    if setup.scf_rescue:  # converged the P3c H3 bead where DIIS, damping and rabuck oscillated
        lines.append("  cgmin")  # (docs/validation.md); it prints no <S2>
    return [*lines, "end"]


def _driver(maxiter: int, options: Sequence[str]) -> list[str]:
    return ["driver", f"  maxiter {maxiter}", *options, "  xyz final", "end"]


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
                 _driver(OPT_MAXITER, options), ["task dft optimize"])


def render_frequencies(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """The only deck with a Hessian task; <name>.hess is written to permanent_dir."""
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 ["task dft frequencies"])


def render_saddle(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT, *,
                  init_hessian: bool = True) -> str:
    """Eigenvector following of mode 1 (``moddir 1``): the only negative mode of <name>.hess,
    which the adapter shapes along the reaction direction (vibrations.shape_hessian)."""
    options = ["  trust 0.1", "  sadstp 0.1", *(["  inhess 2"] if init_hessian else []),
               "  moddir 1"]
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 _driver(SADDLE_MAXITER, options), ["task dft saddle"])


def render_string(start: Molecule, end: Molecule, method: MethodSpec, setup: Setup = _DEFAULT,
                  *, nbeads: int, initial_path: bool = False) -> str:
    """Zero-temperature string with frozen ends; ``initial_path`` reads INITIAL_PATH."""
    string = ["string", f"  nbeads {nbeads}", f"  maxiter {STRING_MAXITER}",
              "  stepsize 0.05", "  interpol 3", "  tol 1e-5", "  freeze1 .true.",
              "  freezeN .true.", "  impose"]
    string += [f"  xyz_path {INITIAL_PATH}"] if initial_path else []
    return _deck(setup, _system(start, method, setup, end=end), _dft(start, method, setup),
                 [*string, "end"], ["task dft string"])


def render_wft(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """CCSD(T) single point with frozen atomic cores. Closed shells use the RHF ``ccsd``
    module; open shells a high-spin ROHF reference (``nopen`` = m - 1) and the TCE, with
    ``2eorb 2emet 13`` (without them CH3O./def2-TZVPD exceeds 1200 MB/rank).

    ``freeze atomic`` freezes no orbital of an ECP atom (I keeps 4s4p4d correlated, G27).
    The ccsd module may take 50 iterations (default 20).
    """
    if method.kind != "wft" or method.wft_method is None:
        raise ValueError(f"method {method.id!r} is not a wave-function method")
    open_shell = mol.multiplicity > 1
    scf = ["scf", *(["  rohf", f"  nopen {mol.multiplicity - 1}"] if open_shell else []),
           f"  maxiter {SCF_MAXITER}",
           *([f"  vectors input {setup.name}.movecs"] if setup.restart_vectors else []), "end"]
    if open_shell:
        body, task = ["tce", "  2eorb", "  2emet 13", "  ccsd(t)", "  freeze atomic", "end"], "tce"
    else:
        body, task = ["ccsd", "  freeze atomic", "  maxiter 50", "end"], "ccsd(t)"
    return _deck(setup, _system(mol, method, setup), scf, body, [f"task {task} energy"])


def hess_text(hessian: np.ndarray) -> str:
    """NWChem ``.hess``: the lower triangle row by row in Eh/bohr², Fortran D exponents."""
    h = np.asarray(hessian, dtype=float)
    return "".join(f"{v:20.10E}\n".replace("E", "D") for v in h[np.tril_indices(len(h))])
