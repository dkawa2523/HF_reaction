"""Reaction hypotheses: pairs of minima that become reaction cases (design §8.2).

Priority: declared reactions, discovery products (low-level TS first), mode-follow TS candidates.
Declared reactions, torsions and inversions included, are always kept so that decide() classifies
them. An undeclared pair joins two DFT minima of one level inside the window, or one basin, and
changes bonds between its ends (CH-07): a conformer change, torsion or enantiomerization is studied
only when declared (Curtin-Hammett), so one basin gives a degenerate rearrangement that exchanges
bonded partners (S10). Bonds and degeneracy are judged on the basins' optimized structures in the
endpoints' atom order and handedness (``identity.basin_coords``), never on input coordinates.
A discovery's ends are its own ``source_species`` and ``product_species``, labelled along the
structures it followed, so a basin representative's arbitrary labelling never makes or hides a
bond change (analysis X2); a discovery without both ends gives no hypothesis.
A mode-follow saddle of the DFT tier is a verified saddle (``ts_calc``); any other discovery TS
is a low-level TS. Negative discoveries never veto a hypothesis (review X1). A relaxation product
(a seed state lost at screen, R6) runs from the DFT basin of the seed's own job, while it keeps
the seed's state, to the DFT basin of its collapse basin; it offers no TS.

A hypothesis that forms exactly one bond, between two fragments of its reactant, and breaks
none, in a composition with monomer states (``thermo.monomer_states``) whose DFT minima lie on
its level, is an association (design X4): its reactant side is the separated monomers
(``ReactionRecord.monomers``), each the lowest DFT minimum of its state on the complex's level
(charge and multiplicity aside), as a barrierless association has an asymptote, not a minimum,
there. The complex stays ``minima[0]``: it identifies the hypothesis and gives the formed bond,
but ends no path. A complex without a DFT minimum of its own (it relaxed into the adduct, as
BH3 + NH3 does) is judged on its own input structure, and ``minima[0]`` is then the adduct's
basin (the driver takes that structure as the reactant end). A formation of two or more bonds (a cycloaddition)
stays an ordinary hypothesis: a one-bond relaxed scan cannot follow it.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from hfauto.chemistry import identity, topology
from hfauto.chemistry.gates import Policy, same_pes
from hfauto.chemistry.thermo import Monomers, State
from hfauto.chemistry.xyz import XYZ, composition_key, hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Geometry, Level
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
Ends = tuple[SpeciesRecord, SpeciesRecord]


def seed_species_id(relaxation: DiscoveryRecord) -> str:
    """The species a relaxation seed is refined as at DFT (R6), named as explore names a
    discovery's new structure: the seed itself may represent its collapse basin."""
    return f"spc_{relaxation.discovery_id}"


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
    monomers: Monomers  # thermo.monomer_states
    levels: Mapping[str, Level]  # DFT minimum id -> the Level of its optimization

    def lowest(self, state: State, level: Level) -> MinimumRecord | None:
        """The lowest DFT minimum of ``state`` on ``level``, charge and multiplicity aside."""
        found = [m for m in self.minima.values()
                 if m.tier == "dft" and (m.composition_id, m.state_label) == state
                 and (own := self.levels.get(m.minimum_id)) is not None
                 and same_pes(level, own, state=False)]
        return min(found, key=lambda m: m.energy_hartree) if found else None

    def coords(self, minimum: MinimumRecord, species: SpeciesRecord) -> np.ndarray:
        """The basin's optimized structure in the atom order and handedness of ``species``."""
        basin = self.load(self.structure[minimum.minimum_id])
        if species.species_id == minimum.species_id:  # the representative itself
            return np.asarray(basin.coords, dtype=float)
        own = self.load(species.geometry).coords
        return identity.basin_coords(basin.symbols, basin.coords, own)

    def ends(self, d: DiscoveryRecord) -> Ends | None:
        """The discovery's own ends (source_species, product_species); None unless both are
        species."""
        sa, sb = (self.species.get(s or "") for s in (d.source_species, d.product_species))
        return None if sa is None or sb is None else (sa, sb)

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


@dataclass(frozen=True)
class _Candidate:
    """A discovery product as a pair of basins, with the TS the discovery offers."""

    source: Source
    start: MinimumRecord | None  # the DFT basin of the discovery's source (a relaxation's seed)
    end: MinimumRecord | None  # the basin of its product
    ends: Ends | None  # the discovery's own source and product species (_Pool.ends)
    low_level_ts: Geometry | None
    ts_calc: str | None  # a verified DFT saddle: a mode-follow saddle of the DFT tier


def _monomers(pool: _Pool, complex_: MinimumRecord, species: SpeciesRecord, x: np.ndarray,
              change: tuple[frozenset[topology.Bond], frozenset[topology.Bond]]
              ) -> tuple[MinimumRecord, ...]:
    """The separated monomers of an association, each repeated by its count (module doc); ()
    for any other hypothesis. ``change``: (formed, broken) from the reactant ``x``."""
    (formed, broken), symbols = change, species.geometry.symbols
    parts = pool.monomers.get((hill_formula(symbols), species.charge), [])
    level = pool.levels.get(complex_.minimum_id)
    if len(formed) != 1 or broken or not parts or level is None:
        return ()
    [(i, j)] = formed
    if any(i in f and j in f for f in topology.fragments(symbols, x)):
        return ()  # a ring closed within one fragment
    found = [(pool.lowest(state, level), count) for state, count in parts]
    monomers = tuple(m for m, count in found if m is not None for _ in range(count))
    return monomers if all(m is not None for m, _ in found) else ()


def _separated(monomers: Sequence[MinimumRecord]) -> tuple[StoichTerm, ...]:
    counts = Counter(m.composition_id for m in monomers)
    return tuple(StoichTerm(composition_id=c, coefficient=n) for c, n in counts.items())


def _record(pool: _Pool, rid: str, source: Source, minima: tuple[MinimumRecord, MinimumRecord],
            ends: Ends, coords: tuple[np.ndarray, np.ndarray], *,
            coordinate: tuple[CoordinateTerm, ...] = (), torsional: bool | None = None,
            low_level_ts: Geometry | None = None, ts_calc: str | None = None) -> ReactionRecord:
    (ma, mb), (sa, sb), (xa, xb) = minima, ends, coords
    symbols = sa.geometry.symbols
    if ma.minimum_id == mb.minimum_id:  # a complex relaxed into the adduct: its own structure
        own = np.asarray(pool.load(sa.geometry).coords, dtype=float)
        xa = own if _monomers(pool, ma, sa, own, topology.bond_changes(symbols, own, xb)) else xa
    formed, broken = topology.bond_changes(symbols, xa, xb)
    monomers = _monomers(pool, ma, sa, xa, (formed, broken))
    return ReactionRecord(
        reaction_id=rid, source=source, products=_stoich(sb),
        reactants=_separated(monomers) if monomers else _stoich(sa),
        monomers=tuple(m.minimum_id for m in monomers),
        minima=(ma.minimum_id, mb.minimum_id), endpoints=(sa.species_id, sb.species_id),
        # One basin reached through a relabelling (NH3 inversion) or as its mirror image
        # (enantiomerization) stays a reaction (CH-35).
        degenerate=ma.basin_id == mb.basin_id and identity.mapped_equivalent(symbols, xa, xb),
        coordinate=coordinate, torsional=not (formed or broken) if torsional is None else torsional,
        low_level_ts=low_level_ts, ts_calc=ts_calc,
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
    return _record(pool, reaction.id, "declared", (ma, mb), (sa, sb), coords,
                   coordinate=coordinate, torsional=reaction.torsional)


def _relaxation(pool: _Pool, d: DiscoveryRecord) -> _Candidate:
    """Seed -> collapse (R6): the DFT basin of the seed's own job (``source_species``), only
    while it keeps the seed's state (a seed that collapsed at DFT too gives none), and its
    collapse basin's. Both ends carry the seed's labelling."""
    ends = pool.ends(d)
    basin = pool.basin_of.get(ends[0].species_id) if ends else None
    kept = ends is not None and basin is not None and basin.state_label == ends[0].state_label
    return _Candidate(source="discovery", start=basin if kept else None,
                      end=pool.dft_basin(d.source_minimum), ends=ends, low_level_ts=None,
                      ts_calc=None)


def _candidates(pool: _Pool, discoveries: Iterable[DiscoveryRecord]) -> Iterator[_Candidate]:
    products = [d for d in discoveries if d.outcome == "product" and d.product_species]
    for d in sorted(products, key=lambda d: (d.mechanism == "mode_follow", d.ts is None)):
        if d.mechanism == "relaxation":
            yield _relaxation(pool, d)
            continue
        yield _Candidate(
            source="mode_follow" if d.mechanism == "mode_follow" else "discovery",
            start=pool.dft_basin(d.source_minimum), end=pool.basin_of.get(d.product_species or ""),
            ends=pool.ends(d), low_level_ts=None if d.ts_calc else d.ts, ts_calc=d.ts_calc)


def _auto(pool: _Pool, c: _Candidate, ma: MinimumRecord, mb: MinimumRecord
          ) -> ReactionRecord | None:
    """An undeclared hypothesis: DFT minima of one level and composition inside the window whose
    ends, the discovery's own, differ in bonds (CH-07); in one basin, a degenerate
    rearrangement."""
    same_level = ma.tier == mb.tier == "dft" and ma.level_key == mb.level_key
    if not same_level or ma.composition_id != mb.composition_id or c.ends is None:
        return None
    if (mb.energy_hartree - ma.energy_hartree) * HARTREE_TO_KCAL_MOL > pool.window_kcal:
        return None
    sa, sb = c.ends
    rid = reaction_id(c.source, sha256_text(f"{ma.minimum_id}|{mb.minimum_id}")[:10])
    record = _record(pool, rid, c.source, (ma, mb), c.ends,
                     (pool.coords(ma, sa), pool.coords(mb, sb)), low_level_ts=c.low_level_ts,
                     ts_calc=c.ts_calc)
    return None if record.torsional else record


def _pool(minima: Iterable[tuple[MinimumRecord, Geometry]], species: Iterable[SpeciesRecord],
          load_xyz: LoadXYZ, window_kcal: float, monomers: Monomers | None = None,
          levels: Mapping[str, Level] | None = None) -> _Pool:
    pairs = list(minima)
    by_id = {m.minimum_id: m for m, _ in pairs}
    basin_of: dict[str, MinimumRecord] = {}
    for m in sorted(by_id.values(), key=lambda m: m.tier != "dft"):  # DFT minima first
        for s in (m.species_id, *m.members):
            basin_of.setdefault(s, m)
    return _Pool(species={s.species_id: s for s in species}, minima=by_id,
                 structure={m.minimum_id: g for m, g in pairs}, basin_of=basin_of,
                 load=load_xyz, window_kcal=window_kcal, monomers=monomers or {},
                 levels=levels or {})


def _lend(record: ReactionRecord, c: _Candidate) -> ReactionRecord:
    """The TS a repeated pair's first hypothesis lacks: a low-level TS, and a verified DFT saddle
    even when a low-level TS was lent first (the case validates the saddle first)."""
    return record.model_copy(update={"low_level_ts": record.low_level_ts or c.low_level_ts,
                                     "ts_calc": record.ts_calc or c.ts_calc})


def select(minima: Iterable[tuple[MinimumRecord, Geometry]], species: Iterable[SpeciesRecord],
           discoveries: Iterable[DiscoveryRecord], declared: Sequence[ReactionInput],
           load_xyz: LoadXYZ, *, window_kcal: float = Policy.reaction_window_kcal,
           max_per_composition: int = 6, monomers: Monomers | None = None,
           levels: Mapping[str, Level] | None = None) -> list[ReactionRecord]:
    """One ReactionRecord per hypothesis. A repeated minima pair keeps the first one, which
    borrows the TS it lacks, also from a discovery that is no hypothesis (an inversion's saddle).
    ``minima`` pairs every minimum (any tier) with its optimized structure; ``monomers``
    (thermo.monomer_states) and ``levels`` (DFT minimum id -> its opt Level) find an
    association's separated monomers."""
    pool = _pool(minima, species, load_xyz, window_kcal, monomers, levels)
    records = [_declared(pool, r) for r in declared]
    index = {frozenset(r.minima): i for i, r in enumerate(records)}
    # by the composition of the case's ends (an association's reactants are its monomers)
    per_composition = Counter(t.composition_id for r in records
                              for t in (r.products or r.reactants)[:1])
    for c in _candidates(pool, discoveries):
        ma, mb = c.start, c.end
        if ma is None or mb is None:
            continue
        pair = frozenset((ma.minimum_id, mb.minimum_id))
        if pair in index:
            records[index[pair]] = _lend(records[index[pair]], c)
            continue
        if per_composition[ma.composition_id] >= max_per_composition:
            continue
        record = _auto(pool, c, ma, mb)
        if record is not None:
            index[pair] = len(records)
            records.append(record)
            per_composition[ma.composition_id] += 1
    return records
