"""Permutation-invariant structure identity (§5.5).

Same-element atoms are matched with scipy's Hungarian solver, alternating with a proper Kabsch
fit from several starting orientations (the alternation only finds a local optimum). A basin
(``assign``) admits proper and improper rotations: mirror images are one minimum, and
``basin_coords`` keeps the handedness of a chiral one (``is_chiral``). Labelled comparisons
(``mapped_rmsd``, ``same_as_labelled``) stay proper, so that a degenerate rearrangement such as
the NH3 inversion differs from the identity. ``carry`` moves a structure along the relabelling,
mirror and rotation that match two others; ``is_image`` judges an exact image (IMAGE_A), and
``member_coords`` gives a basin's structure in a member species' atom order.
"""

from __future__ import annotations

import math
from collections.abc import Hashable, Mapping, Sequence

import numpy as np
from scipy.optimize import linear_sum_assignment

from hfauto.chemistry import topology
from hfauto.chemistry.geometry import kabsch

IMAGE_A = 0.005  # an exact image: QRC ± starts at a symmetric TS <= 0.0009 A, others >= 0.025 A
BASIN_A = 0.05  # one basin: permutation-invariant RMSD (mirror image included) ...
BASIN_DE_HARTREE = 5.0e-5  # ... and |dE|, with the best match clearly ahead of the runner-up
_RUNNER_UP_RATIO, _RUNNER_UP_GAP_A = 3.0, 0.1
_MAX_REFINE = 5
_DEGENERATE_REL = 0.05  # principal moments this close (relative) count as degenerate
_IN_PLANE_STEP_DEG = 30
Labels = tuple[Sequence[Hashable], Sequence[Hashable]]  # per atom of a and of b
# Proper sign flips of the principal axes.
_SIGN_FLIPS = tuple(np.diag(s) for s in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))


def _centered(coords: np.ndarray) -> np.ndarray:
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    return x - x.mean(axis=0)


def _mirror(coords: np.ndarray) -> np.ndarray:
    return np.asarray(coords, dtype=float).reshape(-1, 3) * [-1.0, 1.0, 1.0]


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def _principal_frame(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ascending moments of the unit-mass inertia tensor and a right-handed axis matrix."""

    moments, axes = np.linalg.eigh(np.trace(x.T @ x) * np.eye(3) - x.T @ x)
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1.0
    return moments, axes


def _spins(moments: np.ndarray) -> list[np.ndarray]:
    """In-plane rotations (principal frame) for each degenerate pair of moments."""

    spins = [np.eye(3)]
    scale = max(float(moments.max()), 1e-12)
    for i, j in ((1, 2), (0, 2), (0, 1)):
        if abs(moments[i] - moments[j]) > _DEGENERATE_REL * scale:
            continue
        for step in range(1, 360 // _IN_PLANE_STEP_DEG):
            angle = math.radians(step * _IN_PLANE_STEP_DEG)
            spin = np.eye(3)
            spin[i, i] = spin[j, j] = math.cos(angle)
            spin[i, j], spin[j, i] = -math.sin(angle), math.sin(angle)
            spins.append(spin)
    return spins


def _assign(labels: Labels, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """perm with b[perm[i]] matched to a[i], atoms of one label (element or WL class) only."""

    perm = np.arange(len(a))
    la, lb = (np.asarray(x) for x in labels)
    for label in set(labels[0]):
        ia, ib = np.flatnonzero(la == label), np.flatnonzero(lb == label)
        cost = np.sum((a[ia, None, :] - b[None, ib, :]) ** 2, axis=-1)
        rows, cols = linear_sum_assignment(cost)
        perm[ia[rows]] = ib[cols]
    return perm


def _check(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xa, xb = _centered(a), _centered(b)
    if xa.shape != xb.shape or len(xa) != len(symbols):
        raise ValueError("structures must have the same atoms as symbols")
    return xa, xb


def permutation_invariant_rmsd(
    symbols: Sequence[str], a: np.ndarray, b: np.ndarray, labels: Labels | None = None
) -> tuple[float, np.ndarray]:
    """(rmsd, perm) such that b[perm] is the best proper-rotation match of a; atoms match
    within their ``labels`` of a and b (default: the elements)."""

    xa, xb = _check(symbols, a, b)
    labels = labels or (symbols, symbols)
    moments, frame_a = _principal_frame(xa)
    _, frame_b = _principal_frame(xb)
    starts = [kabsch(xb, xa)]  # identity mapping first
    spins = _spins(moments)
    starts += [frame_b @ flip @ spin @ frame_a.T for flip in _SIGN_FLIPS for spin in spins]
    best = (math.inf, np.arange(len(xa)))
    for start in starts:
        perm = _assign(labels, xa, xb @ start)
        for _ in range(_MAX_REFINE):
            new_perm = _assign(labels, xa, xb @ kabsch(xb[perm], xa))
            if np.array_equal(new_perm, perm):
                break
            perm = new_perm
        value = _rmsd(xb[perm] @ kabsch(xb[perm], xa), xa)
        if value < best[0]:
            best = (value, perm)
    return best


def mapped_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Kabsch RMSD with the identity atom mapping (proper rotations only)."""

    xa, xb = _centered(a), _centered(b)
    if xa.shape != xb.shape:
        raise ValueError("structures differ in atom count")
    return _rmsd(xb @ kabsch(xb, xa), xa)


def _basin_match(symbols: Sequence[str], a: np.ndarray, b: np.ndarray,
                 labels: Labels | None = None) -> tuple[float, np.ndarray, bool]:
    """(rmsd, perm, mirrored): the better of b and its mirror image matched onto a."""

    rmsd, perm = permutation_invariant_rmsd(symbols, a, b, labels)
    m_rmsd, m_perm = permutation_invariant_rmsd(symbols, a, _mirror(b), labels)
    return (m_rmsd, m_perm, True) if m_rmsd < rmsd else (rmsd, perm, False)


def carry(symbols: Sequence[str], ref: np.ndarray, x: np.ndarray, other: np.ndarray
          ) -> tuple[float, np.ndarray]:
    """(rmsd, carried): the basin RMSD of x matched onto ref (relabelling and mirror image
    included), and ``other``, a structure in the atom order and frame of x, carried by that
    relabelling, mirror, Kabsch rotation and centroid shift into the atom order and frame of ref.
    The caller judges the RMSD (IMAGE_A for an exact image, BASIN_A for one basin)."""

    rmsd, perm, mirrored = _basin_match(symbols, ref, x)
    sign = np.array([-1.0 if mirrored else 1.0, 1.0, 1.0])
    xs, ys = (np.asarray(y, dtype=float).reshape(-1, 3)[perm] * sign for y in (x, other))
    centre, target = xs.mean(axis=0), np.asarray(ref, dtype=float).reshape(-1, 3).mean(axis=0)
    return rmsd, (ys - centre) @ kabsch(xs - centre, _centered(ref)) + target


def is_image(symbols: Sequence[str], ref: np.ndarray, x: np.ndarray) -> bool:
    """x is an exact image of ref: ``carry``'s basin RMSD within IMAGE_A. On a PES invariant
    under the exchange of like nuclei and inversion, what x relaxes to is then the image of
    what ref relaxes to."""

    return _basin_match(symbols, ref, x)[0] <= IMAGE_A


def same_as_labelled(xa: np.ndarray, xb: np.ndarray, ea: float, eb: float) -> bool:
    """One structure as labelled: |ea - eb| <= 5e-5 Eh and the identity-mapped RMSD (proper
    rotations) <= 0.05 A. The two structures of a degenerate rearrangement (the NH3 inversion,
    an enantiomerization) are one basin but two structures as labelled."""

    return abs(ea - eb) <= BASIN_DE_HARTREE and mapped_rmsd(xa, xb) <= BASIN_A


def is_chiral(symbols: Sequence[str], x: np.ndarray) -> bool:
    """True when x and its mirror image are distinct under proper rotations and relabelling."""

    return permutation_invariant_rmsd(symbols, x, _mirror(x))[0] > BASIN_A


def assign(
    symbols: Sequence[str],
    coords: np.ndarray,
    energy: float,
    candidates: Mapping[str, tuple[np.ndarray, float]],
) -> str | None:
    """Id of the candidate (coords, energy) whose basin holds the structure, else None.

    Candidates within 5e-5 Eh come first; the closest of them (mirror image included) must lie
    within 0.05 A and clearly ahead of the runner-up."""

    scored = sorted((_basin_match(symbols, coords, xyz)[0], key)
                    for key, (xyz, e) in candidates.items() if abs(energy - e) <= BASIN_DE_HARTREE)
    if not scored or scored[0][0] > BASIN_A:
        return None
    if len(scored) > 1:
        best, runner = scored[0][0], scored[1][0]
        separated = runner - best >= _RUNNER_UP_GAP_A or (
            runner >= _RUNNER_UP_RATIO * best and runner > best)
        if not separated:
            return None
    return scored[0][1]


def basin_coords(symbols: Sequence[str], basin: np.ndarray, own: np.ndarray) -> np.ndarray:
    """The basin's structure in the atom order of ``own`` (a structure of that basin), mirrored
    when the basin is chiral and ``own`` has the other handedness. When ``own`` has the basin's
    bond graph, an atom maps only onto one of its WL class, so the relabelling keeps own's
    atom-indexed bonds however far own lies from the basin (an xTB product of a DFT basin);
    else (own changed state in the basin's relaxation) same elements match."""

    x = np.asarray(basin, dtype=float).reshape(-1, 3)
    y = np.asarray(own, dtype=float).reshape(-1, 3)
    labels = None
    if topology.state_label(symbols, y) == topology.state_label(symbols, x):
        la, lb = (topology.wl_classes(symbols, topology.bonds(symbols, c)) for c in (y, x))
        labels = la, lb
    if not is_chiral(symbols, x):
        return x[permutation_invariant_rmsd(symbols, y, x, labels)[1]]
    _, perm, mirrored = _basin_match(symbols, y, x, labels)
    return (_mirror(x) if mirrored else x)[perm]


def member_coords(symbols: Sequence[str], basin: np.ndarray, rep_species_id: str,
                  species_id: str, own: np.ndarray) -> np.ndarray:
    """A basin's optimized structure (its representative's, ``rep_species_id``) for a member
    species ``species_id`` whose own structure is ``own``: the input coordinates unchanged,
    bit for bit, for the representative itself (a numerical alignment could move the job keys
    downstream), else ``basin_coords``."""

    same = species_id == rep_species_id
    return np.asarray(basin, dtype=float) if same else basin_coords(symbols, basin, own)


def mapped_equivalent(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> bool:
    """True when a and b differ as labelled (proper rotations) but are one basin structure
    (mirror image included): a degenerate rearrangement such as the NH3 inversion or the
    enantiomerization of a chiral minimum."""

    return mapped_rmsd(a, b) > BASIN_A and _basin_match(symbols, a, b)[0] <= BASIN_A


def periodic_nearest(value_deg: float, targets_deg: Sequence[float]) -> int:
    """Index of the target closest to value on a 360° circle (first one on ties)."""

    if not targets_deg:
        raise ValueError("no targets")
    distances = [abs((value_deg - t + 180.0) % 360.0 - 180.0) for t in targets_deg]
    return int(np.argmin(distances))
