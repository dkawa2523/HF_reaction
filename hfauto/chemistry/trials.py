"""Reaction trials for single-ended discovery (§8.2 explore): one bond-graph edit enumerator.

A reaction is a local edit of the bond graph G (ZStruct, YARP): 1–2 bonds formed, 0–2 broken.
Mechanism types (1,2-elimination, [3,3], [4+2], 1,n-shifts, SN2, abstraction, adducts) are
results of the enumeration, never inputs.

- A formed pair is any pair at graph distance ≥ 2, inside or across fragments.
- Locality: the edited bonds form one simple path (chain), or every changed atom lies on one
  simple cycle of ≤ 6 atoms of G + formed (ring).
- Fit: an atom that gains bonds ends with at most ``elements.max_coordination`` neighbours.
- The product graph has a Lewis structure (the source graph is not tested, so TMA·(HF)2 with an
  H bonded to N and F is still enumerated): bond orders 1–3, formal charges in {−1, 0, +1},
  unpaired electrons = multiplicity − 1, and every atom's bond orders plus radicals in its
  valence set: min(ve', s − ve') of its ve' = ve − charge valence electrons (s = 2 for H and He,
  else 8); groups 15–16 from period 3 on also take ve' − 2k up to 6 (P {3, 5}, S {2, 4, 6}).
  d-block atoms have no such valence, so a source with a bonded one gives no class (out of
  scope).
  Charge separation is not limited: probe P0f's limit (ions, ion pairs and 1,2-dipoles only)
  lost a product state that P0c found (the ethyl formate 1,5-H shift, a 1,5-zwitterion).

A class is an orbit of edits under the automorphisms of the element-labelled source graph: the
canonical SMILES of G with kept, formed and broken bonds as three bond types. So classes do not
depend on atom numbering or geometry. State labels distinguish stereoisomers, but edit
classes do not: paths that only geometry tells apart (syn/anti, diastereotopic H, E/Z) share
one class. Only one realization per class and state is tried, so distinct stereochemical
paths within that state are not exhaustively sampled.

A class runs once per state, realised on the conformer and member edit of least drive value
Σ r/Σr_cov over its formed pairs. When no intermolecular formed pair is in contact (r ≤ Σr_vdW)
there, the closest one is placed rigidly at PLACE_RATIO·Σr_cov: the smaller fragment turns
about its atom to face the partner along the partner's outward direction (from its fragment's
centroid) and moves there. A linear start is bent by 10° so that NT2 does not start on a
symmetry line. Trials come in budget order: fewest bond changes, then the drive value in the
conformer, then the class. The budget is the caller's.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.spatial.transform import Rotation

from hfauto.chemistry import topology
from hfauto.chemistry.elements import covalent_radius, max_coordination, vdw_radius
from hfauto.chemistry.geometry import neighbours
from hfauto.chemistry.vibrations import external_basis

Pair = tuple[int, int]
Edit = tuple[frozenset[Pair], frozenset[Pair]]  # (formed, broken)

MAX_RING = 6
PLACE_RATIO = 1.5  # r/Σr_cov of a placed pair: the lower end of review U4-P2's 1.5–2, not tuned
LINEAR_BEND_DEG = 10.0
_EXPANDED = frozenset({"P", "As", "Sb", "Bi", "S", "Se", "Te", "Po"})  # groups 15–16, period ≥ 3


@dataclass(frozen=True, eq=False)
class Trial:
    """One edit class of a state, realised on one of its conformers."""

    key: str  # the class
    formed: tuple[Pair, ...]
    broken: tuple[Pair, ...]
    conformer: int  # index into the conformers given to ``trials``
    drive: float  # Σ r/Σr_cov of the formed pairs in that conformer
    start: np.ndarray  # NT2 start structure

    @property
    def order(self) -> tuple[int, float, str]:
        return len(self.formed) + len(self.broken), self.drive, self.key


def _pair(i: int, j: int) -> Pair:
    return (i, j) if i < j else (j, i)


def _subsets(items: Sequence[Pair], sizes: Sequence[int]) -> Iterator[frozenset[Pair]]:
    return (frozenset(c) for k in sizes for c in itertools.combinations(items, k))


def _chains(nbr: Sequence[set[int]]) -> set[Edit]:
    """Edits whose edited bonds form one simple path."""
    out: set[Edit] = set()

    def grow(path: list[int], formed: tuple[Pair, ...], broken: tuple[Pair, ...]) -> None:
        if formed:
            out.add((frozenset(formed), frozenset(broken)))
        for w in set(range(len(nbr))).difference(path):
            step = (_pair(path[-1], w),)
            if w in nbr[path[-1]] and len(broken) < 2:
                grow([*path, w], formed, broken + step)
            elif w not in nbr[path[-1]] and len(formed) < 2:
                grow([*path, w], formed + step, broken)

    for v in range(len(nbr)):
        grow([v], (), ())
    return out


def _rings(nbr: Sequence[set[int]], bonds: frozenset[Pair]) -> set[Edit]:
    """Edits whose changed atoms lie on one simple cycle of ≤ MAX_RING atoms of G + formed."""
    cycles: defaultdict[frozenset[int], set[frozenset[Pair]]] = defaultdict(set)

    def walk(path: list[int], new: frozenset[Pair]) -> None:  # new: the formed edges it uses
        for w in (range(len(nbr)) if len(new) < 2 else nbr[path[-1]]):
            used = new if w in nbr[path[-1]] else new | {_pair(path[-1], w)}
            if w == path[0] and len(path) >= 3:
                cycles[frozenset(path)].add(used)
            elif w not in path and len(path) < MAX_RING:
                walk([*path, w], used)

    for v in range(len(nbr)):
        walk([v], frozenset())
    out: set[Edit] = set()
    for atoms, uses in cycles.items():
        pairs = list(itertools.combinations(sorted(atoms), 2))
        inner = [p for p in pairs if p in bonds]
        out.update((formed, broken)
                   for formed in _subsets([p for p in pairs if p not in bonds], (1, 2))
                   if any(u <= formed for u in uses) for broken in _subsets(inner, (0, 1, 2)))
    return out


def _class(symbols: Sequence[str], bonds: frozenset[Pair], edit: Edit) -> str:
    """Canonical SMILES with kept, formed and broken bonds as single, double and triple."""
    from rdkit import Chem
    formed, broken = edit
    mol = Chem.RWMol(Chem.MolFromSmiles(".".join(f"[{s}]" for s in symbols), sanitize=False))
    for i, j in bonds | formed:
        mol.AddBond(i, j, Chem.BondType.DOUBLE if (i, j) in formed else
                    Chem.BondType.TRIPLE if (i, j) in broken else Chem.BondType.SINGLE)
    return Chem.MolToSmiles(mol, allBondsExplicit=True)


def _valences(symbol: str, ve: int) -> range:
    shell = 2 if symbol in ("H", "He") else 8
    base = min(ve, shell - ve)  # negative: no valence
    top = min(ve, 6) if symbol in _EXPANDED else base
    return range(base, top + 1, 2) if base >= 0 else range(0)


def lewis(symbols: Sequence[str], bonds: frozenset[Pair], charge: int, multiplicity: int) -> bool:
    """Whether the bond graph has a Lewis structure under the module's rules: a MILP over one
    option (charge, radicals, valence) per atom and a π order 0–2 per bond."""
    from rdkit import Chem
    outer, n = Chem.GetPeriodicTable().GetNOuterElecs, len(symbols)
    edges = np.array(sorted(bonds), dtype=int).reshape(-1, 2)
    deg = np.bincount(edges.ravel(), minlength=n)
    opts = np.array([(i, c, r, v - r - deg[i]) for i, s in enumerate(symbols) for c in (-1, 0, 1)
                     for v in _valences(s, outer(s) - c) for r in range(v - deg[i] + 1)],
                    dtype=int).reshape(-1, 4)
    if len(set(opts[:, 0])) < n:
        return False
    atom, c, r, spare = opts.T
    pick = np.equal.outer(np.arange(n), atom).astype(float)
    incidence = np.eye(n)[:, edges[:, 0]] + np.eye(n)[:, edges[:, 1]]
    a = np.vstack([np.hstack([pick, 0 * incidence]),  # one option per atom
                   np.hstack([-pick * spare, incidence]),  # π bonds use the option's spare valence
                   np.r_[c, 0 * edges[:, 0]], np.r_[r, 0 * edges[:, 0]]])
    target = np.r_[np.ones(n), np.zeros(n), charge, multiplicity - 1]
    upper = np.r_[np.ones(len(opts)), np.full(len(edges), 2.0)]
    return milp(np.zeros(len(upper)), integrality=np.ones(len(upper)), bounds=Bounds(0, upper),
                constraints=LinearConstraint(a, target, target)).status == 0


def edits(symbols: Sequence[str], coords: np.ndarray, charge: int,
          multiplicity: int) -> dict[str, list[Edit]]:
    """The structure's edit classes whose product has a Lewis structure: class -> its member
    edits in this atom numbering."""
    bonds, members = topology.bonds(symbols, coords), defaultdict(list)
    partners = neighbours(bonds)
    nbr = [partners.get(i, set()) for i in range(len(symbols))]
    room = [max_coordination(s) - len(nbr[i]) for i, s in enumerate(symbols)]
    for formed, broken in _chains(nbr) | _rings(nbr, bonds):
        net = Counter(i for pair in formed for i in pair)
        net.subtract(i for pair in broken for i in pair)
        if all(gain <= room[i] for i, gain in net.items() if gain > 0):  # it fits
            members[_class(symbols, bonds, (formed, broken))].append((formed, broken))
    return {key: es for key, es in members.items()
            if lewis(symbols, (bonds - es[0][1]) | es[0][0], charge, multiplicity)}


def trials(symbols: Sequence[str], conformers: Sequence[np.ndarray], charge: int,
           multiplicity: int) -> list[Trial]:
    """Every edit class of one state (its conformers in one atom order), realised on the
    conformer and member edit of least drive, in budget order."""
    cov = np.array([covalent_radius(s) for s in symbols])
    xs = [np.asarray(c, dtype=float).reshape(-1, 3) for c in conformers]
    found: dict[frozenset[Pair], dict[str, list[Edit]]] = {}
    best: dict[str, tuple[float, int, tuple[Pair, ...], tuple[Pair, ...]]] = {}
    for k, x in enumerate(xs):
        bonds = topology.bonds(symbols, x)
        if bonds not in found:
            found[bonds] = edits(symbols, x, charge, multiplicity)
        ratio = np.linalg.norm(x[:, None] - x[None], axis=-1) / np.add.outer(cov, cov)
        for key, members in found[bonds].items():
            drive, formed, broken = min((float(sum(ratio[p] for p in f)), tuple(sorted(f)),
                                         tuple(sorted(b))) for f, b in members)
            if key not in best or drive < best[key][0]:
                best[key] = (drive, k, formed, broken)
    out, starts = [], [_bend_linear(symbols, x) for x in xs]
    for key, (drive, k, formed, broken) in best.items():
        x = _place(symbols, xs[k], formed)
        start = starts[k] if x is xs[k] else _bend_linear(symbols, x)
        out.append(Trial(key, formed, broken, k, drive, start))
    return sorted(out, key=lambda t: t.order)


def _place(symbols: Sequence[str], x: np.ndarray, formed: Sequence[Pair]) -> np.ndarray:
    """x itself, or a copy with one fragment placed (see the module docstring)."""
    group = {i: g for g in topology.fragments(symbols, x) for i in g}
    cov, vdw = ({p: f(symbols[p[0]]) + f(symbols[p[1]]) for p in formed}
                for f in (covalent_radius, vdw_radius))
    r = {p: float(np.linalg.norm(x[p[0]] - x[p[1]])) for p in formed}
    inter = [p for p in formed if group[p[0]] != group[p[1]]]
    if not inter or any(r[p] <= vdw[p] for p in inter):
        return x
    a, b = min(inter, key=lambda p: (r[p] / cov[p], p))
    distance = PLACE_RATIO * cov[a, b]
    if len(group[a]) < len(group[b]):
        a, b = b, a  # b's fragment, the smaller one, moves

    def outward(i: int, default: np.ndarray) -> np.ndarray:
        v = x[i] - x[list(group[i])].mean(axis=0)
        return v / np.linalg.norm(v) if np.linalg.norm(v) > 0.1 else default

    axis = outward(a, (x[b] - x[a]) / np.linalg.norm(x[b] - x[a]))
    turn = Rotation.align_vectors(-axis, outward(b, -axis))[0]
    moved, guest = x.copy(), list(group[b])
    moved[guest] = turn.apply(x[guest] - x[b]) + x[a] + distance * axis
    return moved


def _bend_linear(symbols: Sequence[str], x: np.ndarray) -> np.ndarray:
    """A linear structure with its atoms on one side of the middle atom turned by 10° about a
    fixed normal of the axis; any other structure unchanged."""
    if len(symbols) < 3 or external_basis(symbols, x).shape[1] != 5:
        return x
    axis = np.linalg.svd(x - x.mean(axis=0))[2][0]
    along = (x - x.mean(axis=0)) @ axis
    pivot = int(np.argsort(along)[len(along) // 2])
    normal = np.cross(axis, np.eye(3)[int(np.argmin(np.abs(axis)))])
    turn = Rotation.from_rotvec(np.radians(LINEAR_BEND_DEG) * normal / np.linalg.norm(normal))
    bent, side = x.copy(), along < along[pivot]
    bent[side] = turn.apply(x[side] - x[pivot]) + x[pivot]
    return bent
