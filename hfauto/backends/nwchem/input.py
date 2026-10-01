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
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.method import MethodSpec

OPT_MAXITER = 100
SADDLE_MAXITER = 50
STRING_MAXITER = 20
SCF_MAXITER = 100  # NWChem 7.2.3 defaults: DFT 50, SCF 20
INITIAL_PATH = "initial_path.xyz"  # xyz_path is read from permanent_dir (the job's cwd)
# Eh/bohr²: NWChem 7.2.3 adds k (r - r0)² (r, r0 in bohr; not in "Total DFT energy"), so r
# settles F/2k short of r0: 7e-4 Å at 0.05 Eh/bohr, the steepest CH3 + O2 scan force (a
# water O-H held at 1.2 Å: k 5 / 20 / 50 -> 1.1948 / 1.1987 / 1.1995 Å in 6 / 8 / 10 steps).
SPRING_K = 20.0
_VDW = {"d3zero": 3, "d3bj": 4}
_NOBLE_GAS_Z = (2, 10, 18, 36, 54)  # elements.py ends at Rn
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
    scf_rescue: bool = False  # after SCF_NOT_CONVERGED: cgmin, then one plain SCF (_task)


_DEFAULT = Setup()


def _deck(setup: Setup, *blocks: Sequence[str], memory: str = "") -> str:
    header = ["echo", f"start {setup.name}", "permanent_dir ."]
    if setup.scratch_dir:
        header.append(f"scratch_dir {setup.scratch_dir}")
    header.append(memory or f"memory total {setup.memory_mb} mb")
    return "\n\n".join("\n".join(block) for block in (header, *blocks) if block) + "\n"


def _task(setup: Setup, operation: str) -> tuple[list[str], ...]:
    """``task dft <operation>``; after a cgmin rescue, one plain SCF from its vectors at the
    same structure prints the <S2> cgmin does not (CH3: 2 iterations, 6e-9 Eh apart)."""
    if not setup.scf_rescue:
        return ([f"task dft {operation}"],)
    return ([f"task dft {operation}"], ["unset dft:cgmin"],
            ["dft", f"  vectors input {setup.name}.movecs", "end"], ["task dft energy"])


def _geometry(xyz: XYZ, *, cartesian: bool, label: str = "", zcoord: Sequence[str] = ()
              ) -> list[str]:
    flags = "units angstrom nocenter noautosym" + (" noautoz" if cartesian else "")
    rows = np.asarray(xyz.coords, dtype=float).reshape(-1, 3)
    lines = [f"geometry {label} {flags}" if label else f"geometry {flags}"]
    lines += [f"  {s:2s} {x: .10f} {y: .10f} {z: .10f}" for s, (x, y, z) in
              zip(xyz.symbols, rows, strict=True)]
    return [*lines, *zcoord, "end"]


def _fixed(bond: tuple[int, int, float] | None, cartesian: bool) -> tuple[list[str], list[str]]:
    """(zcoord lines inside the geometry, constraints block) holding atoms i, j (0-based) r Å
    apart: a frozen ``zcoord`` bond (NWChem moves the input to its value, so the start must
    have it), or in Cartesian coordinates (noautoz) a SPRING_K ``spring bond`` restraint."""
    if bond is None:
        return [], []
    i, j, r = bond
    if cartesian:
        r0 = r / BOHR_TO_ANGSTROM
        return [], ["constraints", f"  spring bond {i + 1} {j + 1} {SPRING_K} {r0:.6f}", "end"]
    return ["  zcoord", f"    bond {i + 1} {j + 1} {r:.4f} rc constant", "  end"], []


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


def frozen_core(symbols: Sequence[str]) -> int:
    """Frozen orbitals: per atom, NWChem's ``freeze atomic`` core (the noble gas before its
    row, src/geom/geom_numcore.F) less its def2-ECP electrons (NWChem's library def2-ecp: 28
    for Rb-Xe, 46 for Cs-La, 60 for Hf-Rn), never below zero. ``freeze atomic`` itself
    freezes nothing on an ECP atom (G27: I keeps 4s4p correlated)."""
    total = 0
    for z in map(atomic_number, symbols):
        core = max((g for g in _NOBLE_GAS_Z if g < z), default=0)
        ecp_electrons = 0 if z <= 36 else 28 if z <= 54 else 46 if z <= 57 else 60
        total += max(0, core - ecp_electrons) // 2
    return total


def _system(mol: Molecule, method: MethodSpec, setup: Setup, *, end: Molecule | None = None,
            zcoord: Sequence[str] = ()) -> list[str]:
    """Geometry (and the string's end geometry), charge, basis and ECP."""
    if not method.basis:
        raise ValueError(f"method {method.id!r} names no basis set")
    lines = _geometry(mol.xyz, cartesian=setup.cartesian, zcoord=zcoord)
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
    if setup.scf_rescue:
        lines.append("  cgmin")
    return [*lines, "end"]


def _driver(maxiter: int, options: Sequence[str]) -> list[str]:
    return ["driver", f"  maxiter {maxiter}", *options, "  xyz final", "end"]


def render_energy(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 *_task(setup, "energy"))


def render_optimize(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT, *,
                    init_hessian: bool = False,
                    fixed_bond: tuple[int, int, float] | None = None) -> str:
    """Default driver thresholds; ``init_hessian`` reads <name>.hess; ``fixed_bond`` (i, j, r
    Å; 0-based) holds a bond (_fixed).

    ``trust 0.3`` (the NWChem default) with an initial Hessian (a QRC or mode-follow side from
    its saddle, a complex from xTB), ``trust 0.1`` without: from a displaced start the diagonal
    guess overshoots back above the TS energy (HCN->HNC in NWChem 7.2.3). NWChem clamps a
    minimization step along a negative eigenvalue to 0.03-0.3 x trust and its BFGS update keeps
    the sign (opt_drv.F 7.2.3): the adapter writes every initial Hessian but a higher-order
    saddle's as its positive-definite model.
    """
    options = ["  trust 0.3", "  inhess 2"] if init_hessian else ["  trust 0.1"]
    zcoord, constraints = _fixed(fixed_bond, setup.cartesian)
    return _deck(setup, _system(mol, method, setup, zcoord=zcoord), constraints,
                 _dft(mol, method, setup), _driver(OPT_MAXITER, options),
                 *_task(setup, "optimize"))


def render_frequencies(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """The only deck with a Hessian task; <name>.hess is written to permanent_dir."""
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 *_task(setup, "frequencies"))


def render_saddle(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT, *,
                  init_hessian: bool = True) -> str:
    """Eigenvector following of mode 1 (``moddir 1``): the only negative mode of <name>.hess,
    which the adapter shapes along the reaction direction (vibrations.shape_hessian)."""
    options = ["  trust 0.1", "  sadstp 0.1", *(["  inhess 2"] if init_hessian else []),
               "  moddir 1"]
    return _deck(setup, _system(mol, method, setup), _dft(mol, method, setup),
                 _driver(SADDLE_MAXITER, options), *_task(setup, "saddle"))


def render_string(start: Molecule, end: Molecule, method: MethodSpec, setup: Setup = _DEFAULT,
                  *, nbeads: int, initial_path: bool = False) -> str:
    """Zero-temperature string with frozen ends; ``initial_path`` reads INITIAL_PATH. No SCF
    follows a rescue: a PathProfile carries no <S2>."""
    string = ["string", f"  nbeads {nbeads}", f"  maxiter {STRING_MAXITER}",
              "  stepsize 0.05", "  interpol 3", "  tol 1e-5", "  freeze1 .true.",
              "  freezeN .true.", "  impose"]
    string += [f"  xyz_path {INITIAL_PATH}"] if initial_path else []
    return _deck(setup, _system(start, method, setup, end=end), _dft(start, method, setup),
                 [*string, "end"], ["task dft string"])


def render_wft(mol: Molecule, method: MethodSpec, setup: Setup = _DEFAULT) -> str:
    """CCSD(T) single point with the frozen core of frozen_core. Closed shells use the RHF
    ``ccsd`` module (up to 50 iterations; default 20); open shells a high-spin ROHF reference
    (``nopen`` = m - 1) and the TCE, with ``2eorb 2emet 13`` (without them CH3O./def2-TZVPD
    exceeds 1200 MB/rank).

    Of the memory per rank, Global Arrays (the (T) amplitudes) get 70 % instead of the 50 % of
    ``memory total``: malonaldehyde/def2-TZVPD (227 functions) failed to allocate at ``total
    1200`` and ran at ``heap 100 mb stack 500 mb global 1300 mb`` (a unit after every size).
    """
    if method.kind != "wft" or method.wft_method is None:
        raise ValueError(f"method {method.id!r} is not a wave-function method")
    open_shell = mol.multiplicity > 1
    scf = ["scf", *(["  rohf", f"  nopen {mol.multiplicity - 1}"] if open_shell else []),
           f"  maxiter {SCF_MAXITER}",
           *([f"  vectors input {setup.name}.movecs"] if setup.restart_vectors else []), "end"]
    freeze = f"  freeze {frozen_core(mol.xyz.symbols)}"
    if open_shell:
        body, task = ["tce", "  2eorb", "  2emet 13", "  ccsd(t)", freeze, "end"], "tce"
    else:
        body, task = ["ccsd", freeze, "  maxiter 50", "end"], "ccsd(t)"
    mb = setup.memory_mb
    heap, stack = mb * 5 // 100, mb * 25 // 100
    return _deck(setup, _system(mol, method, setup), scf, body, [f"task {task} energy"],
                 memory=f"memory heap {heap} mb stack {stack} mb global {mb - heap - stack} mb")


def hess_text(hessian: np.ndarray) -> str:
    """NWChem ``.hess``: the lower triangle row by row in Eh/bohr², Fortran D exponents."""
    h = np.asarray(hessian, dtype=float)
    return "".join(f"{v:20.10E}\n".replace("E", "D") for v in h[np.tril_indices(len(h))])
