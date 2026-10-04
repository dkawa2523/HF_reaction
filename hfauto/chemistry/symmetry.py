"""The point group of a stationary structure, one judgement for the thermochemistry (design §8.2
thermo).

libmsym (pymsym 0.3.5, the detector GoodVibes uses) proposes a group at three threshold
settings; an exception is no proposal. A linear proposal is symmetrized here: projected onto
the principal axis (C∞v) and, averaged with its relabelled inversion, D∞h. A proposal counts
only when its symmetrized structure stays in the basin of the structure by the criterion that
makes two minima one (identity): the harmonic rise ½ δxᵀ|H|δx <= BASIN_DE_HARTREE, |H| from the
absolute eigenvalues of the freq Hessian and δx after a proper Kabsch fit with the atom order
kept, and RMSD <= BASIN_A. The accepted group of highest order wins; C1, the structure itself,
is always accepted. The thermochemistry evaluates the modes and the moments at the symmetrized
structure.

The group's operations give sigma and m (Fernández-Ramos et al., Theor. Chem. Acc. 118, 813
(2007)): sigma counts the identity and the proper rotations (C∞v 1, D∞h 2), and m = 2 when
there is no improper operation (S_n, reflection, inversion), else 1.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

import numpy as np

from hfauto.chemistry.elements import mass
from hfauto.chemistry.geometry import kabsch
from hfauto.chemistry.identity import BASIN_A, BASIN_DE_HARTREE
from hfauto.core.constants import BOHR_TO_ANGSTROM

# libmsym thresholds: equivalence 2e-3 (its default 5e-4 split a C3v complex left 0.03° off the
# axis into Cs), then all seven at 1e-2 and at 5e-2. All seven at 1e-2 find the group of every
# validation minimum and TS under 1e-3 Å noise (tools/symcheck.py); a default 'zero' of 1e-3
# fails there. 5e-2 alone errs (I-...CH3I as Cs): the basin test and the highest order decide.
_THRESHOLDS = ("zero", "geometry", "angle", "equivalence", "eigfact", "permutation",
               "orthogonalization")
_SETTINGS = ({"equivalence": 2.0e-3}, dict.fromkeys(_THRESHOLDS, 1.0e-2),
             dict.fromkeys(_THRESHOLDS, 5.0e-2))
_IMPROPER = frozenset({2, 3, 4})  # libmsym operation types: S_n, reflection, inversion
Rank = tuple[bool, int]  # (linear, order); D∞h ranks 2 and C∞v 1 among the linear groups


@dataclass(frozen=True)
class Symmetry:
    point_group: str  # Schoenflies as libmsym names it; Cinfv or Dinfh if linear, Kh an atom
    sigma: int  # external rotational symmetry number: identity + proper rotations of the group
    linear: bool
    m: int  # optical isomer number: 2 for a group without an improper operation, else 1
    coords: np.ndarray  # the symmetrized structure (Å, in the frame of the input)


def analyze(symbols: Sequence[str], coords: np.ndarray, hessian: np.ndarray) -> Symmetry:
    """The accepted point group of highest order of a structure with its freq Hessian
    (Eh/bohr², in the frame of ``coords``)."""
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) == 1:
        return Symmetry("Kh", 1, False, 1, x)
    h = np.asarray(hessian, dtype=float)
    if h.shape != (x.size, x.size):
        raise ValueError(f"the Hessian must be ({x.size}, {x.size})")
    curvature = np.linalg.eigh(0.5 * (h + h.T))
    best, top = Symmetry("C1", 1, False, 2, x), (False, 1)
    for rank, proposal in _proposals(symbols, x).values():
        fitted = _in_basin(x, proposal.coords, curvature)
        if fitted is not None and rank > top:
            best, top = replace(proposal, coords=fitted), rank
    return best


def _proposals(symbols: Sequence[str], x: np.ndarray) -> dict[str, tuple[Rank, Symmetry]]:
    """Each setting's group with its rank and symmetrized structure (the tightest first)."""
    import pymsym

    out: dict[str, tuple[Rank, Symmetry]] = {}
    for thresholds in _SETTINGS:
        elements = [pymsym.Element(name=s, coordinates=[float(v) for v in c])
                    for s, c in zip(symbols, x, strict=True)]
        try:
            with pymsym.Context(elements=elements) as ctx:
                ctx.set_thresholds(**thresholds)
                group = ctx.find_symmetry()
                types = [op.type for op in ctx.symmetry_operations]
                y = np.array([e.coordinates for e in ctx.symmetrize_elements()], dtype=float)
        except Exception:  # libmsym raises when no group fits
            continue
        if group is None:
            continue
        if group in ("C0v", "D0h"):
            found = _linear(symbols, x)
        else:
            proper = sum(t not in _IMPROPER for t in types)
            m = 2 if proper == len(types) else 1
            found = {group: ((False, len(types)), Symmetry(group, proper, False, m, y))}
        for name, proposal in found.items():
            out.setdefault(name, proposal)
    return out


def _linear(symbols: Sequence[str], x: np.ndarray) -> dict[str, tuple[Rank, Symmetry]]:
    """C∞v: x projected onto its principal axis through the centre of mass; D∞h: that projection
    averaged with its inversion, the atoms of each element paired in reverse order along it."""
    masses = np.array([mass(s) for s in symbols])
    com = masses @ x / masses.sum()
    c = x - com
    axis = np.linalg.eigh((masses[:, None] * c).T @ c)[1][:, -1]  # the least moment of inertia
    z = c @ axis
    pair = np.arange(len(z))
    for element in set(symbols):
        same = [i for i in np.argsort(z) if symbols[i] == element]
        pair[same] = same[::-1]
    return {"Cinfv": ((True, 1), Symmetry("Cinfv", 1, True, 1, com + np.outer(z, axis))),
            "Dinfh": ((True, 2), Symmetry("Dinfh", 2, True, 1,
                                          com + np.outer(0.5 * (z - z[pair]), axis)))}


def _in_basin(x: np.ndarray, y: np.ndarray, curvature: tuple[np.ndarray, np.ndarray]
              ) -> np.ndarray | None:
    """y fitted onto x (proper rotation, atom order kept) when it lies in the basin of x."""
    yc, centre = y - y.mean(axis=0), x.mean(axis=0)
    fitted = yc @ kabsch(yc, x - centre) + centre
    d = fitted - x
    values, vectors = curvature
    rise = 0.5 * float(np.abs(values) @ (vectors.T @ d.ravel() / BOHR_TO_ANGSTROM) ** 2)
    rmsd = float(np.sqrt(np.mean(np.sum(d**2, axis=1))))
    return fitted if rise <= BASIN_DE_HARTREE and rmsd <= BASIN_A else None
