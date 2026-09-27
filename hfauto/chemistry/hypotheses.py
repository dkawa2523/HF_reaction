"""Reaction hypotheses: pairs of minima that become reaction cases (design §8.2).

Priority: declared reactions, discovery products (low-level TS first), mode-follow TS candidates.
Declared reactions are always kept so that decide() classifies them (blocked, same basin, out of
window); the others must join two DFT basins of one level inside the window and show a change
(CH-07).  Two conformers of one state are no hypothesis: a fast pre-equilibrium does not enter
the ranking (Curtin-Hammett), and a torsion is studied only when declared.  A discovery whose
two ends fall into one DFT basin is kept when its ends map onto each other by a non-identity
permutation: a degenerate rearrangement, evaluated like a declared one.  ``minima`` includes
the screen minima that discoveries start from.
Degeneracy and bond changes are judged on the basins' optimized structures in the endpoints'
atom order and handedness (``identity.basin_coords``), never on input coordinates.
Negative discoveries never veto a hypothesis (review X1); only the summary reports them.
"""

from __future__ import annotations

from collections import Counter
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
Source = Literal["declared", "discovery", "mode_follow"]
Ends = tuple[str, str] | None  # a discovery's source and product species
Candidate = tuple[Source, MinimumRecord | None, MinimumRecord | None, Geometry | None, Ends]


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


def _max_distance_change(a: np.ndarray, b: np.ndarray) -> float:
    da = np.linalg.norm(a[:, None] - a[None], axis=-1)
    db = np.linalg.norm(b[:, None] - b[None], axis=-1)
    return float(np.max(np.abs(da - db)))


def _stoich(species: SpeciesRecord | None) -> tuple[StoichTerm, ...]:
    if species is None:
        return ()
    key = composition_key(species.geometry.symbols, species.charge, species.multiplicity)
    return (StoichTerm(composition_id=key, coefficient=1),)


@dataclass(frozen=True)
class _Pool:
    species: dict[str, SpeciesRecord]
    minima: dict[str, MinimumRecord]  # every tier, by minimum_id
    structure: dict[str, Geometry]  # minimum_id -> its optimized structure
    basin_of: dict[str, MinimumRecord]  # species id -> minimum holding it (DFT first)
    load: LoadXYZ
    window_kcal: float
    min_distance_A: float

    def coords(self, minimum: MinimumRecord, species: SpeciesRecord) -> np.ndarray:
        """The basin's optimized structure in the atom order and handedness of ``species``."""
        basin = self.load(self.structure[minimum.minimum_id])
        if species.species_id == minimum.species_id:  # the representative itself
            return np.asarray(basin.coords, dtype=float)
        own = self.load(species.geometry).coords
        return identity.basin_coords(basin.symbols, basin.coords, own)

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


def _record(rid: str, source: Source, minima: tuple[MinimumRecord, MinimumRecord],
            ends: tuple[SpeciesRecord, SpeciesRecord], coords: tuple[np.ndarray, np.ndarray], *,
            coordinate: tuple[CoordinateTerm, ...] = (), torsional: bool | None = None,
            low_level_ts: Geometry | None = None) -> ReactionRecord:
    (ma, mb), (sa, sb), (xa, xb) = minima, ends, coords
    symbols = sa.geometry.symbols
    formed, broken = topology.bond_changes(symbols, xa, xb)
    return ReactionRecord(
        reaction_id=rid, source=source, reactants=_stoich(sa), products=_stoich(sb),
        minima=(ma.minimum_id, mb.minimum_id), endpoints=(sa.species_id, sb.species_id),
        # One basin reached through a relabelling (NH3 inversion) or as its mirror image
        # (enantiomerization) stays a reaction (CH-35).
        degenerate=ma.basin_id == mb.basin_id and identity.mapped_equivalent(symbols, xa, xb),
        coordinate=coordinate, torsional=not (formed or broken) if torsional is None else torsional,
        low_level_ts=low_level_ts,
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
    coords = (pool.coords(ma, sa), pool.coords(mb, sb))
    return _record(reaction.id, "declared", (ma, mb), (sa, sb), coords,
                   coordinate=coordinate, torsional=reaction.torsional)


def _candidates(pool: _Pool, discoveries: Sequence[DiscoveryRecord]) -> Iterator[Candidate]:
    products = [d for d in discoveries if d.outcome == "product" and d.product_species]
    for d in sorted(products, key=lambda d: (d.mechanism == "mode_follow", d.ts is None)):
        source: Source = "mode_follow" if d.mechanism == "mode_follow" else "discovery"
        product = pool.basin_of.get(d.product_species or "")
        start = pool.minima.get(d.source_minimum)
        ends = (start.species_id, d.product_species or "") if start else None
        yield source, pool.dft_basin(d.source_minimum), product, d.ts, ends


def _degenerate(pool: _Pool, rid: str, source: Source, basin: MinimumRecord, ends: Ends,
                ts: Geometry | None) -> ReactionRecord | None:
    """Both ends of a discovery in one basin: a degenerate rearrangement when the basin in the
    atom orders of the two ends differs as labelled (``_record``), else no reaction."""
    sa, sb = (pool.species.get(e) for e in ends) if ends else (None, None)
    if sa is None or sb is None:
        return None
    coords = (pool.coords(basin, sa), pool.coords(basin, sb))
    record = _record(rid, source, (basin, basin), (sa, sb), coords, low_level_ts=ts)
    return record if record.degenerate else None


def _auto(pool: _Pool, source: Source, ma: MinimumRecord, mb: MinimumRecord,
          ts: Geometry | None, ends: Ends) -> ReactionRecord | None:
    same_level = ma.tier == mb.tier == "dft" and ma.level_key == mb.level_key
    if not same_level or ma.composition_id != mb.composition_id:
        return None
    rid = reaction_id(source, sha256_text(f"{ma.minimum_id}|{mb.minimum_id}")[:10])
    if ma.basin_id == mb.basin_id:
        return _degenerate(pool, rid, source, ma, ends, ts)
    if (mb.energy_hartree - ma.energy_hartree) * HARTREE_TO_KCAL_MOL > pool.window_kcal:
        return None
    picked = pick_endpoints(pool.members(ma), pool.members(mb), pool.load)
    if picked is None:
        return None
    sa, sb = pool.species[picked[0]], pool.species[picked[1]]
    xa, xb = pool.coords(ma, sa), pool.coords(mb, sb)
    changed = any(topology.bond_changes(sa.geometry.symbols, xa, xb))
    if not (changed or _max_distance_change(xa, xb) >= pool.min_distance_A):  # CH-07
        return None
    return _record(rid, source, (ma, mb), (sa, sb), (xa, xb), low_level_ts=ts)


def _pool(minima: Iterable[tuple[MinimumRecord, Geometry]], species: Iterable[SpeciesRecord],
          load_xyz: LoadXYZ, window_kcal: float, min_distance_A: float) -> _Pool:
    pairs = list(minima)
    by_id = {m.minimum_id: m for m, _ in pairs}
    basin_of: dict[str, MinimumRecord] = {}
    for m in sorted(by_id.values(), key=lambda m: m.tier != "dft"):  # DFT minima first
        for s in (m.species_id, *m.members):
            basin_of.setdefault(s, m)
    return _Pool(species={s.species_id: s for s in species}, minima=by_id,
                 structure={m.minimum_id: g for m, g in pairs}, basin_of=basin_of,
                 load=load_xyz, window_kcal=window_kcal, min_distance_A=min_distance_A)


def _lend_ts(record: ReactionRecord, ts: Geometry | None) -> ReactionRecord:
    """A repeated minima pair keeps its first hypothesis, which borrows a low-level TS it lacks."""
    if record.low_level_ts is None and ts is not None:
        return record.model_copy(update={"low_level_ts": ts})
    return record


def select(minima: Iterable[tuple[MinimumRecord, Geometry]], species: Iterable[SpeciesRecord],
           discoveries: Iterable[DiscoveryRecord], declared: Sequence[ReactionInput],
           load_xyz: LoadXYZ, *, window_kcal: float = Policy.reaction_window_kcal,
           min_distance_A: float = 0.2, max_per_composition: int = 6) -> list[ReactionRecord]:
    """One ReactionRecord per hypothesis; a repeated minima pair keeps the first hypothesis.
    ``minima`` pairs every minimum (any tier) with its optimized structure."""
    found = list(discoveries)
    pool = _pool(minima, species, load_xyz, window_kcal, min_distance_A)
    records = [_declared(pool, r) for r in declared]
    index = {frozenset(r.minima): i for i, r in enumerate(records)}
    per_composition = Counter(r.reactants[0].composition_id for r in records if r.reactants)
    for source, ma, mb, ts, ends in _candidates(pool, found):
        if ma is None or mb is None:
            continue
        pair = frozenset((ma.minimum_id, mb.minimum_id))
        if pair in index:
            records[index[pair]] = _lend_ts(records[index[pair]], ts)
            continue
        if per_composition[ma.composition_id] >= max_per_composition:
            continue
        record = _auto(pool, source, ma, mb, ts, ends)
        if record is not None:
            index[pair] = len(records)
            records.append(record)
            per_composition[ma.composition_id] += 1
    return records
