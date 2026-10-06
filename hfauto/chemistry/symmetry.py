"""The point group a stationary structure is counted in, one judgement for the thermochemistry
(design §7.1 thermo).

Candidates. A non-linear group is a set of operations, each a permutation of like atoms that
keeps every interatomic distance within CANDIDATE_A (a backtracking search; the tolerance only
bounds the search) with the proper or improper orthogonal map that fits it; the identity with an
improper map is the plane of a planar structure. A linear group is the projection onto the
principal axis (C∞v) and, averaged with its relabelled inversion, D∞h.

Acceptance, the zero-point criterion. A group's symmetrised structure y is the average of the
structure's fitted images under its operations. Along d = y − x, mass-weighted with the external
motion projected out, the group counts when the harmonic rise stays within the zero-point energy
of d's own mode: ½ dᵀ|H|d <= ½ħω_d with ω_d² = dᵀ|H|d / dᵀMd, |H| the freq Hessian with its
eigenvalues by magnitude (a saddle counts its reaction mode as a curvature). It depends on no
temperature and no constant, and the harmonic rise overrates a double well's barrier, so it errs
towards the lower group; the boundary moves away, it does not vanish. Exception: an open-shell
structure enters a group with degenerate representations (an operation of order >= 3, or a
linear group) only when y is x within identity's noise (BASIN_DE_HARTREE, BASIN_A), as the
symmetrisation may cross an electronic degeneracy (Jahn–Teller; CH3O in C3v is an E state),
where the harmonic estimate fails.

The cyclic groups of the operations that pass are joined in the order of their rise while the
closure still passes; a linear group ranks above any other. sigma counts the identity and the
proper rotations (C∞v 1, D∞h 2) and m = 2 without an improper operation, else 1 (Fernández-Ramos
et al., Theor. Chem. Acc. 118, 813 (2007)). The thermochemistry takes the energy and the modes
from the structure itself and only the moments of inertia from y.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from hfauto.chemistry.elements import mass
from hfauto.chemistry.geometry import kabsch
from hfauto.chemistry.identity import BASIN_A, BASIN_DE_HARTREE
from hfauto.chemistry.vibrations import external_basis
from hfauto.core.constants import AMU_TO_ME, BOHR_TO_ANGSTROM

# Å; an accepted operation changed a pair distance by at most 0.09 Å over the B1, VAL9, W3 and
# W4 subjects, and 0.3 or 0.5 Å accept the same groups there
CANDIDATE_A = 0.15
Op = tuple[tuple[int, ...], bool]  # (perm, improper): the image has atom i at R·x[perm[i]]
_MIRROR = np.diag([-1.0, 1.0, 1.0])


@dataclass(frozen=True)
class Symmetry:
    point_group: str  # Schoenflies; Cinfv or Dinfh if linear, Kh an atom
    sigma: int  # external rotational symmetry number: identity + proper rotations of the group
    linear: bool
    m: int  # optical isomer number: 2 for a group without an improper operation, else 1
    coords: np.ndarray  # the symmetrised structure (Å, in the frame of the input)


def analyze(symbols: Sequence[str], coords: np.ndarray, hessian: np.ndarray, *,
            open_shell: bool = False) -> Symmetry:
    """The accepted point group of a structure with its freq Hessian (Eh/bohr², in the frame of
    ``coords``); ``open_shell``: the declared multiplicity is above 1."""
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) == 1:
        return Symmetry("Kh", 1, False, 1, x)
    h = np.asarray(hessian, dtype=float)
    if h.shape != (x.size, x.size):
        raise ValueError(f"the Hessian must be ({x.size}, {x.size})")
    accepts = _criterion(symbols, x, 0.5 * (h + h.T), open_shell)
    for name, sigma, y in _linear(symbols, x):
        if accepts(y, True) is not None:
            return Symmetry(name, sigma, True, 1, y)
    centre = x.mean(axis=0)
    group = _group(symbols, x - centre, accepts)
    proper = sum(not improper for _, improper in group)
    y = _average(x - centre, group) + centre
    return Symmetry(_name(symbols, y, group), proper, False, 2 if proper == len(group) else 1, y)


def _criterion(symbols: Sequence[str], x: np.ndarray, h: np.ndarray, open_shell: bool
               ) -> Callable[[np.ndarray, bool], float | None]:
    """accepts(y, degenerate): the rise ½ dᵀ|H|d (Eh) of y fitted onto x when it passes, else
    None (atomic units: ħ = 1, masses in electron masses)."""
    m = np.repeat([mass(s) for s in symbols], 3) * AMU_TO_ME
    b = external_basis(symbols, x)
    project = np.eye(m.size) - b @ b.T
    values, vectors = np.linalg.eigh(project @ (h / np.sqrt(np.outer(m, m))) @ project)

    def accepts(y: np.ndarray, degenerate: bool) -> float | None:
        d = _fit(y, x) - x
        q = project @ (np.sqrt(m) * d.ravel() / BOHR_TO_ANGSTROM)
        curvature, norm = float(np.abs(values) @ (vectors.T @ q) ** 2), float(q @ q)
        rise = 0.5 * curvature
        if degenerate and open_shell:
            ok = rise <= BASIN_DE_HARTREE and math.sqrt(np.mean(np.sum(d**2, axis=1))) <= BASIN_A
        else:
            ok = norm == 0.0 or rise <= 0.5 * math.sqrt(curvature / norm)  # ½ħω_d
        return rise if ok else None

    return accepts


def _fit(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """y on x: proper rotation about the centroids, atom order kept."""
    yc, centre = y - y.mean(axis=0), x.mean(axis=0)
    return yc @ kabsch(yc, x - centre) + centre


def _linear(symbols: Sequence[str], x: np.ndarray) -> list[tuple[str, int, np.ndarray]]:
    """D∞h then C∞v: x projected onto its principal axis through the centre of mass, D∞h that
    projection averaged with its inversion, the atoms of each element paired in reverse order
    along it."""
    masses = np.array([mass(s) for s in symbols])
    com = masses @ x / masses.sum()
    c = x - com
    axis = np.linalg.eigh((masses[:, None] * c).T @ c)[1][:, -1]  # the least moment of inertia
    z = c @ axis
    pair = np.arange(len(z))
    for element in set(symbols):
        same = [i for i in np.argsort(z) if symbols[i] == element]
        pair[same] = same[::-1]
    return [("Dinfh", 2, com + np.outer(0.5 * (z - z[pair]), axis)),
            ("Cinfv", 1, com + np.outer(z, axis))]


def _permutations(symbols: Sequence[str], x: np.ndarray) -> list[tuple[int, ...]]:
    """Every permutation of like atoms that changes no pair distance by more than CANDIDATE_A,
    placing the atoms farthest from the rest first."""
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    n, order = len(x), np.argsort(-d.sum(axis=1), kind="stable")
    perm, used, out = np.full(n, -1), np.zeros(n, dtype=bool), []

    def place(k: int) -> None:
        if k == n:
            out.append(tuple(int(j) for j in perm))
            return
        i, placed = order[k], order[:k]
        for j in range(n):
            if not used[j] and symbols[j] == symbols[i] and np.all(
                    np.abs(d[i, placed] - d[j, perm[placed]]) <= CANDIDATE_A):
                perm[i], used[j] = j, True
                place(k + 1)
                perm[i], used[j] = -1, False

    place(0)
    return out


def _group(symbols: Sequence[str], x: np.ndarray,
           accepts: Callable[[np.ndarray, bool], float | None]) -> list[Op]:
    """The operations of the accepted group of the centred structure x."""
    identity: Op = (tuple(range(len(x))), False)
    passed: list[tuple[float, Op]] = []
    for perm in _permutations(symbols, x):
        for op in ((perm, False), (perm, True)):
            if op != identity:
                cyclic = _closure([op])
                rise = accepts(_average(x, cyclic), _degenerate(cyclic))
                if rise is not None:
                    passed.append((rise, op))
    group, generators = [identity], []
    for _, op in sorted(passed):
        if op not in group:
            trial = _closure([*generators, op])
            if accepts(_average(x, trial), _degenerate(trial)) is not None:
                group, generators = trial, [*generators, op]
    return group


def _closure(generators: list[Op]) -> list[Op]:
    """The group the operations generate."""
    n = len(generators[0][0])
    group: set[Op] = {(tuple(range(n)), False)}
    frontier = list(group)
    while frontier:
        new = {(tuple(p[i] for i in g), gi != pi) for g, gi in frontier for p, pi in generators}
        frontier = list(new - group)
        group |= new
    return sorted(group)


def _degenerate(group: list[Op]) -> bool:
    """An operation of order >= 3, i.e. whose permutation squared moves an atom (in a
    non-linear structure the map follows from the permutation up to its handedness): the group
    has a degenerate representation."""
    return any(perm[perm[i]] != i for perm, _ in group for i in range(len(perm)))


def _image(x: np.ndarray, op: Op) -> np.ndarray:
    """The centred x under op, fitted onto x."""
    perm, improper = op
    y = x[list(perm)] @ _MIRROR if improper else x[list(perm)]
    return y @ kabsch(y, x)


def _average(x: np.ndarray, group: list[Op]) -> np.ndarray:
    return sum((_image(x, op) for op in group), np.zeros_like(x)) / len(group)


def _name(symbols: Sequence[str], y: np.ndarray, group: list[Op]) -> str:
    """The Schoenflies name libmsym gives the symmetrised structure; its order is the group's."""
    if len(group) == 1:
        return "C1"
    import pymsym

    elements = [pymsym.Element(name=s, coordinates=[float(v) for v in c])
                for s, c in zip(symbols, y, strict=True)]
    with pymsym.Context(elements=elements) as ctx:
        return str(ctx.find_symmetry())
