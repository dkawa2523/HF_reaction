"""Reaction hypotheses: pairs of minima that become reaction cases (design §8.2).

Priority: declared reactions, discovery products (low-level TS first), mode-follow TS candidates,
conformer pairs.  Declared reactions are always kept so that decide() classifies them (blocked,
same basin, out of window); the others must join two DFT basins of one level inside the window
and show a change (CH-07).  ``minima`` includes the screen minima that discoveries start from.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from hfauto.chemistry import identity, topology
from hfauto.chemistry.gates import Policy
from hfauto.chemistry.xyz import XYZ, composition_key
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Geometry
from hfauto.core.hashing import sha256_text
from hfauto.core.ids import reaction_id
from hfauto.core.records import (
    CoordinateTerm,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
    StoichTerm,
)
from hfauto.core.system import ReactionInput

LoadXYZ = Callable[[Geometry], XYZ]
Source = Literal["declared", "discovery", "mode_follow", "conformer"]
Candidate = tuple[Source, MinimumRecord | None, MinimumRecord | None, Geometry | None]
Bonds = frozenset[topology.Bond]


def pick_endpoints(
    members_a: Sequence[SpeciesRecord], members_b: Sequence[SpeciesRecord], load_xyz: LoadXYZ
) -> tuple[str, str] | None:
    """Member pair with the smallest identity-mapped RMSD, so that a relabelled representative
    is never used (chem 13).  None when no pair shares the element order."""
    scored = [
        (identity.mapped_rmsd(load_xyz(a.geometry).coords, load_xyz(b.geometry).coords), i, j)
        for i, a in enumerate(members_a)
        for j, b in enumerate(members_b)
        if a.geometry.symbols == b.geometry.symbols
    ]
    if not scored:
        return None
    _, i, j = min(scored)
    return members_a[i].species_id, members_b[j].species_id


def _dihedrals_deg(x: np.ndarray, quads: list[tuple[int, int, int, int]]) -> np.ndarray:
    p = x[np.asarray(quads)]
    b0, b1, b2 = p[:, 0] - p[:, 1], p[:, 2] - p[:, 1], p[:, 3] - p[:, 2]
    b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - np.sum(b0 * b1, axis=1)[:, None] * b1
    w = b2 - np.sum(b2 * b1, axis=1)[:, None] * b1
    return np.degrees(np.arctan2(np.sum(np.cross(b1, v) * w, axis=1), np.sum(v * w, axis=1)))


def _max_torsion_change(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> float:
    """Largest periodic change of a dihedral i-j-k-l along bonds (bonds taken from a)."""
    bonded = topology.bonds(symbols, a)
    neighbours: defaultdict[int, set[int]] = defaultdict(set)
    for i, j in bonded:
        neighbours[i].add(j)
        neighbours[j].add(i)
    quads = [(i, j, k, m) for j, k in bonded
             for i in neighbours[j] - {k} for m in neighbours[k] - {j} if i != m]
    if not quads:
        return 0.0
    delta = _dihedrals_deg(b, quads) - _dihedrals_deg(a, quads)
    return float(np.max(np.abs((delta + 180.0) % 360.0 - 180.0)))


def _max_distance_change(a: np.ndarray, b: np.ndarray) -> float:
    da = np.linalg.norm(a[:, None] - a[None], axis=-1)
    db = np.linalg.norm(b[:, None] - b[None], axis=-1)
    return float(np.max(np.abs(da - db)))


def _stoich(species: SpeciesRecord | None) -> tuple[StoichTerm, ...]:
    if species is None:
        return ()
    key = composition_key(species.geometry.symbols, species.charge, species.multiplicity)
    return (StoichTerm(composition_id=key, coefficient=1),)


def _pairs(atoms: Iterable[tuple[int, int]]) -> Bonds:
    return frozenset((min(p), max(p)) for p in atoms)


@dataclass(frozen=True)
class _Pool:
    species: dict[str, SpeciesRecord]
    minima: dict[str, MinimumRecord]  # every tier, by minimum_id
    basin_of: dict[str, MinimumRecord]  # species id -> minimum holding it (DFT first)
    negatives: tuple[DiscoveryRecord, ...]
    load: LoadXYZ
    window_kcal: float
    min_distance_A: float
    min_angle_deg: float

    def coords(self, species: SpeciesRecord) -> np.ndarray:
        return np.asarray(self.load(species.geometry).coords, dtype=float)

    def members(self, minimum: MinimumRecord) -> list[SpeciesRecord]:
        ids = dict.fromkeys((minimum.species_id, *minimum.members))
        return [self.species[s] for s in ids if s in self.species]

    def dft_basin(self, minimum_id: str) -> MinimumRecord | None:
        """DFT minimum holding the representative or a member of ``minimum_id`` (any tier)."""
        source = self.minima.get(minimum_id)
        held = (source.species_id, *source.members) if source is not None else ()
        found = (self.basin_of.get(s) for s in held)
        return next((m for m in found if m is not None and m.tier == "dft"), None)

    def endpoint_basin(self, species_id: str) -> MinimumRecord | None:
        """Basin of a declared endpoint: a seed that collapsed into a screen basin is not
        refined itself (§8.2), so it takes the DFT minimum of that basin when there is one."""
        own = self.basin_of.get(species_id)
        if own is None or own.tier == "dft":
            return own
        return self.dft_basin(own.minimum_id) or own

    def negative_evidence(self, composition: str, formed: Bonds, broken: Bonds) -> tuple[str, ...]:
        """Negative discoveries of the same composition and bond change (either direction)."""
        wanted = {(formed, broken), (broken, formed)} if formed or broken else set()
        found = []
        for d in self.negatives:
            source = self.minima.get(d.source_minimum)
            if d.trial is None or source is None or source.composition_id != composition:
                continue
            if (_pairs(d.trial.associations), _pairs(d.trial.dissociations)) in wanted:
                found.append(f"{d.discovery_id}:{d.reason or 'negative'}")
        return tuple(found)


def _record(pool: _Pool, rid: str, source: Source, minima: tuple[MinimumRecord, MinimumRecord],
            ends: tuple[SpeciesRecord, SpeciesRecord], *,
            coordinate: tuple[CoordinateTerm, ...] = (), torsional: bool | None = None,
            low_level_ts: Geometry | None = None) -> ReactionRecord:
    (ma, mb), (sa, sb) = minima, ends
    symbols = sa.geometry.symbols
    xa, xb = pool.coords(sa), pool.coords(sb)
    formed, broken = topology.bond_changes(symbols, xa, xb)
    return ReactionRecord(
        reaction_id=rid, source=source, reactants=_stoich(sa), products=_stoich(sb),
        minima=(ma.minimum_id, mb.minimum_id), endpoints=(sa.species_id, sb.species_id),
        # CH-35: one basin reached through a relabelling (NH3 inversion) stays a reaction.
        degenerate=ma.basin_id == mb.basin_id and identity.mapped_equivalent(symbols, xa, xb),
        coordinate=coordinate, torsional=not (formed or broken) if torsional is None else torsional,
        n_h_transferred=topology.transferred_hydrogens(symbols, xa, xb), low_level_ts=low_level_ts,
        negative_evidence=pool.negative_evidence(ma.composition_id, formed, broken),
    )


def _declared(pool: _Pool, reaction: ReactionInput) -> ReactionRecord:
    sa, sb = pool.species.get(reaction.reactant), pool.species.get(reaction.product)
    ma, mb = pool.endpoint_basin(reaction.reactant), pool.endpoint_basin(reaction.product)
    coordinate = tuple(reaction.coordinate)
    if sa is None or sb is None or ma is None or mb is None:
        # An endpoint without a minimum is kept so that decide() blocks it (row 1).
        minima = (ma.minimum_id if ma else "", mb.minimum_id if mb else "")
        return ReactionRecord(reaction_id=reaction.id, source="declared", minima=minima,
                              reactants=_stoich(sa), products=_stoich(sb), coordinate=coordinate,
                              endpoints=(reaction.reactant, reaction.product),
                              torsional=bool(reaction.torsional))
    if sa.geometry.symbols != sb.geometry.symbols:
        raise ValueError(f"reaction {reaction.id}: endpoints differ in atom order")
    return _record(pool, reaction.id, "declared", (ma, mb), (sa, sb),
                   coordinate=coordinate, torsional=reaction.torsional)


def _candidates(pool: _Pool, discoveries: Sequence[DiscoveryRecord]) -> Iterator[Candidate]:
    products = [d for d in discoveries if d.outcome == "product" and d.product_species]
    for d in sorted(products, key=lambda d: (d.mechanism == "mode_follow", d.ts is None)):
        source: Source = "mode_follow" if d.mechanism == "mode_follow" else "discovery"
        product = pool.basin_of.get(d.product_species or "")
        yield source, pool.dft_basin(d.source_minimum), product, d.ts
    groups: defaultdict[tuple[str, str, str], list[MinimumRecord]] = defaultdict(list)
    for m in sorted(pool.minima.values(), key=lambda m: (m.energy_hartree, m.minimum_id)):
        if m.tier == "dft":
            groups[(m.composition_id, m.state_label, m.level_key)].append(m)
    for group in groups.values():
        for ma, mb in itertools.combinations(group, 2):
            yield "conformer", ma, mb, None


def _auto(pool: _Pool, source: Source, ma: MinimumRecord, mb: MinimumRecord,
          ts: Geometry | None) -> ReactionRecord | None:
    same_level = ma.tier == mb.tier == "dft" and ma.level_key == mb.level_key
    if not same_level or ma.composition_id != mb.composition_id or ma.basin_id == mb.basin_id:
        return None
    if (mb.energy_hartree - ma.energy_hartree) * HARTREE_TO_KCAL_MOL > pool.window_kcal:
        return None
    ends = pick_endpoints(pool.members(ma), pool.members(mb), pool.load)
    if ends is None:
        return None
    sa, sb = pool.species[ends[0]], pool.species[ends[1]]
    xa, xb = pool.coords(sa), pool.coords(sb)
    changed = any(topology.bond_changes(sa.geometry.symbols, xa, xb))
    twisted = _max_torsion_change(sa.geometry.symbols, xa, xb) >= pool.min_angle_deg
    # A conformer pair twists without a bond change; others need a large enough change (CH-07).
    moved = changed or twisted or _max_distance_change(xa, xb) >= pool.min_distance_A
    if not ((twisted and not changed) if source == "conformer" else moved):
        return None
    rid = reaction_id(source, sha256_text(f"{ma.minimum_id}|{mb.minimum_id}")[:10])
    return _record(pool, rid, source, (ma, mb), (sa, sb), low_level_ts=ts)


def _pool(minima: Iterable[MinimumRecord], species: Iterable[SpeciesRecord],
          discoveries: Sequence[DiscoveryRecord], load_xyz: LoadXYZ, window_kcal: float,
          min_distance_A: float, min_angle_deg: float) -> _Pool:
    by_id = {m.minimum_id: m for m in minima}
    basin_of: dict[str, MinimumRecord] = {}
    for m in sorted(by_id.values(), key=lambda m: m.tier != "dft"):  # DFT minima first
        for s in (m.species_id, *m.members):
            basin_of.setdefault(s, m)
    return _Pool(species={s.species_id: s for s in species}, minima=by_id, basin_of=basin_of,
                 negatives=tuple(d for d in discoveries if d.outcome == "negative"),
                 load=load_xyz, window_kcal=window_kcal, min_distance_A=min_distance_A,
                 min_angle_deg=min_angle_deg)


def _lend_ts(record: ReactionRecord, ts: Geometry | None) -> ReactionRecord:
    """A repeated minima pair keeps its first hypothesis, which borrows a low-level TS it lacks."""
    if record.low_level_ts is None and ts is not None:
        return record.model_copy(update={"low_level_ts": ts})
    return record


def select(minima: Iterable[MinimumRecord], species: Iterable[SpeciesRecord],
           discoveries: Iterable[DiscoveryRecord], declared: Sequence[ReactionInput],
           load_xyz: LoadXYZ, *, window_kcal: float = Policy.reaction_window_kcal,
           min_distance_A: float = 0.2, min_angle_deg: float = 30.0,
           max_per_composition: int = 6) -> list[ReactionRecord]:
    """One ReactionRecord per hypothesis; a repeated minima pair keeps the first hypothesis."""
    found = list(discoveries)
    pool = _pool(minima, species, found, load_xyz, window_kcal, min_distance_A, min_angle_deg)
    records = [_declared(pool, r) for r in declared]
    index = {frozenset(r.minima): i for i, r in enumerate(records)}
    per_composition = Counter(r.reactants[0].composition_id for r in records if r.reactants)
    for source, ma, mb, ts in _candidates(pool, found):
        if ma is None or mb is None:
            continue
        pair = frozenset((ma.minimum_id, mb.minimum_id))
        if pair in index:
            records[index[pair]] = _lend_ts(records[index[pair]], ts)
            continue
        if per_composition[ma.composition_id] >= max_per_composition:
            continue
        record = _auto(pool, source, ma, mb, ts)
        if record is not None:
            index[pair] = len(records)
            records.append(record)
            per_composition[ma.composition_id] += 1
    return records
