"""Reaction hypotheses: pairs of minima that become reaction cases (design §8.2).

Priority: declared reactions, discovery products (low-level TS first), mode-follow TS candidates.
Declared reactions, torsions and inversions included, are always kept. A hypothesis stands for
its case key (``pair_key``): a pair that differs in chemical state by its two states, else by its
two minima. An undeclared pair whose key a hypothesis already has makes none; every hypothesis of
a key holds each distinct low-level TS of that key (its own first, then in priority order), the
seeds its case tries in turn, and the first verified DFT saddle. An undeclared pair joins two
DFT minima of one level, or one basin, and changes bonds between its ends: a conformer change,
torsion or enantiomerization is studied only when declared (Curtin-Hammett), so one basin gives
a degenerate rearrangement that exchanges bonded partners. A hypothesis that a static check
decides is closed here, a record without a job (``_closed``): an end without a DFT minimum
(BLOCKED), both ends in one basin when the case is neither degenerate nor an association
(SAME_BASIN), or the product more than the window above the reactant asymptote (``uphill``:
OUT_OF_WINDOW). Bonds and degeneracy are judged on the basins' optimized structures in the
endpoints' atom order and handedness (``identity.member_coords``), never on input coordinates.
A discovery's ends are its own ``source_species`` and ``product_species`` (an edge's two
species, labelled along the structures it followed), so a basin representative's arbitrary
labelling never makes or hides a bond change; a discovery without both ends, or not an edge
reached from the start states (``product``), gives no hypothesis. Each end takes the DFT
minimum holding its species, or holding the screen basin it fell into. A mode-follow saddle of
the DFT tier is a verified saddle (``ts_calc``); any other discovery TS is a low-level TS.
Negative discoveries never veto a hypothesis. An edge without a TS (a structure that relaxed
into another state at the low level) gives one only while its source keeps its own state at
DFT; it offers no TS.

A hypothesis that forms exactly one bond, between two fragments of its reactant, and breaks
none, whose reactant fragments are a declared composition's monomers
(``thermo.separated_states``) with DFT minima on its level, is an association: its reactant side is the separated monomers
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
from hfauto.chemistry.thermo import Monomers, State, separated_states
from hfauto.chemistry.xyz import XYZ, composition_key
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Geometry, Level
from hfauto.core.hashing import sha256_text
from hfauto.core.ids import reaction_id
from hfauto.core.records import (
    CaseOutcome,
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
CaseKey = frozenset[tuple[str, str]] | frozenset[str]


def pair_key(a: MinimumRecord, b: MinimumRecord) -> CaseKey:
    """What a case between ``a`` and ``b`` answers, at the granularity it is judged
    (drivers.reaction_case.connection._key): the two chemical states (composition, state
    label) when they differ, as the basins of a state interconvert faster than it reacts
    (Curtin-Hammett); else (a torsion, a degenerate rearrangement) the two minima."""
    sa, sb = (a.composition_id, a.state_label), (b.composition_id, b.state_label)
    return frozenset((sa, sb)) if sa != sb else frozenset((a.minimum_id, b.minimum_id))


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
    monomers: Monomers  # thermo.declared_monomers
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
        basin, own = self.load(self.structure[minimum.minimum_id]), self.load(species.geometry)
        return identity.member_coords(basin.symbols, basin.coords, minimum.species_id,
                                      species.species_id, own.coords)

    def ends(self, d: DiscoveryRecord) -> Ends | None:
        """The discovery's own ends (source_species, product_species); None unless both are
        species."""
        sa, sb = (self.species.get(s or "") for s in (d.source_species, d.product_species))
        return None if sa is None or sb is None else (sa, sb)

    def basin(self, species_id: str) -> MinimumRecord | None:
        """The minimum holding a species: a structure that collapsed into a screen basin is not
        refined itself (§8.2), so it takes the DFT minimum of that basin when there is one."""
        own = self.basin_of.get(species_id)
        if own is None or own.tier == "dft":
            return own
        held = (self.basin_of.get(s) for s in (own.species_id, *own.members))
        return next((m for m in held if m is not None and m.tier == "dft"), own)


@dataclass(frozen=True)
class _Candidate:
    """A discovery product as a pair of basins, with the TS the discovery offers."""

    source: Source
    start: MinimumRecord | None  # the basin of the discovery's source
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
    states = separated_states(XYZ(list(symbols), x), species.charge, pool.monomers)
    level = pool.levels.get(complex_.minimum_id)
    if len(formed) != 1 or broken or not states or level is None:
        return ()
    [(i, j)] = formed
    if any(i in f and j in f for f in topology.fragments(symbols, x)):
        return ()  # a ring closed within one fragment
    found = [pool.lowest(state, level) for state in states]
    return tuple(m for m in found if m is not None) if None not in found else ()


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
        # (enantiomerization) stays a reaction (design §8.2).
        degenerate=ma.basin_id == mb.basin_id and identity.mapped_equivalent(symbols, xa, xb),
        coordinate=coordinate, torsional=not (formed or broken) if torsional is None else torsional,
        low_level_ts=() if low_level_ts is None else (low_level_ts,), ts_calc=ts_calc,
    )


def _declared(pool: _Pool, reaction: ReactionInput) -> ReactionRecord:
    sa, sb = pool.species.get(reaction.reactant), pool.species.get(reaction.product)
    ma, mb = pool.basin(reaction.reactant), pool.basin(reaction.product)
    coordinate = tuple(reaction.coordinate)
    if sa is None or sb is None or ma is None or mb is None:
        # An endpoint without a minimum is kept, closed BLOCKED (``_closed``).
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


def _candidates(pool: _Pool, discoveries: Iterable[DiscoveryRecord]) -> Iterator[_Candidate]:
    products = [d for d in discoveries if d.outcome == "product" and d.product_species]
    for d in sorted(products, key=lambda d: (d.mechanism == "mode_follow", d.ts is None)):
        ends, start = pool.ends(d), pool.basin(d.source_species or "")
        relaxed = d.ts is None and d.ts_calc is None
        if relaxed and (ends is None or start is None or start.state_label != ends[0].state_label):
            start = None  # its source collapsed at DFT as well: nothing to ask
        yield _Candidate(
            source="mode_follow" if d.mechanism == "mode_follow" else "discovery", start=start,
            end=pool.basin(d.product_species or ""), ends=ends,
            low_level_ts=None if d.ts_calc else d.ts, ts_calc=d.ts_calc)


def _auto(pool: _Pool, c: _Candidate, ma: MinimumRecord, mb: MinimumRecord
          ) -> ReactionRecord | None:
    """An undeclared hypothesis: DFT minima of one level and composition whose ends, the
    discovery's own, differ in bonds (design §8.2); in one basin, a degenerate rearrangement."""
    same_level = ma.tier == mb.tier == "dft" and ma.level_key == mb.level_key
    if not same_level or ma.composition_id != mb.composition_id or c.ends is None:
        return None
    sa, sb = c.ends
    rid = reaction_id(c.source, sha256_text(f"{ma.minimum_id}|{mb.minimum_id}")[:10])
    record = _record(pool, rid, c.source, (ma, mb), c.ends,
                     (pool.coords(ma, sa), pool.coords(mb, sb)), low_level_ts=c.low_level_ts,
                     ts_calc=c.ts_calc)
    return None if record.torsional else record


def uphill(minima: Mapping[str, MinimumRecord], case: ReactionRecord) -> float:
    """ΔE (kcal/mol) from the reactant asymptote to the product: from the separated monomers'
    energy sum for an association, else from the reactant minimum."""
    reactant = sum(minima[m].energy_hartree for m in case.monomers or case.minima[:1])
    return (minima[case.minima[1]].energy_hartree - reactant) * HARTREE_TO_KCAL_MOL


def _closed(pool: _Pool, case: ReactionRecord) -> ReactionRecord:
    """The hypothesis, or its closing record when a static check decides it (module doc)."""
    a, b = (pool.minima.get(m) for m in case.minima)
    if a is None or b is None or a.tier != "dft" or b.tier != "dft":
        outcome, reason = CaseOutcome.BLOCKED, "endpoint_without_dft_minimum"
    elif a.basin_id == b.basin_id and not case.degenerate and not case.monomers:
        outcome, reason = CaseOutcome.SAME_BASIN, "same_basin"
    elif uphill(pool.minima, case) > pool.window_kcal:
        outcome, reason = CaseOutcome.OUT_OF_WINDOW, "out_of_window"
    else:
        return case
    return case.model_copy(update={"outcome": outcome, "reasons": (reason,)})


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


def _lend(pool: _Pool, record: ReactionRecord, c: _Candidate) -> ReactionRecord:
    """The TSs of a candidate of the hypothesis's key: its low-level TS is appended unless it
    lies in one basin with a TS held (identity.assign on structure alone: a low-level TS has no
    energy on the case's PES); its verified DFT saddle when the hypothesis has none (the case
    validates it first)."""
    held, ts = record.low_level_ts, c.low_level_ts
    if ts is not None:
        x = pool.load(ts)
        known = {str(k): (pool.load(g).coords, 0.0) for k, g in enumerate(held)}
        if identity.assign(x.symbols, x.coords, 0.0, known) is None:
            held = (*held, ts)
    return record.model_copy(update={"low_level_ts": held, "ts_calc": record.ts_calc or c.ts_calc})


def select(minima: Iterable[tuple[MinimumRecord, Geometry]], species: Iterable[SpeciesRecord],
           discoveries: Iterable[DiscoveryRecord], declared: Sequence[ReactionInput],
           load_xyz: LoadXYZ, *, window_kcal: float = Policy.reaction_window_kcal,
           monomers: Monomers | None = None,
           levels: Mapping[str, Level] | None = None) -> list[ReactionRecord]:
    """One ReactionRecord per hypothesis: each declared reaction, and the first discovery
    candidate of a case key (``pair_key``) no hypothesis has, open or closed (``_closed``). Every
    open hypothesis of a key then takes the TSs of all candidates of that key in priority order
    (``_lend``), also of one that is no hypothesis itself (an inversion's saddle). ``minima``
    pairs every minimum (any tier) with its optimized structure; ``monomers``
    (thermo.declared_monomers) and ``levels`` (DFT minimum id -> its opt Level) find an
    association's separated monomers."""
    pool = _pool(minima, species, load_xyz, window_kcal, monomers, levels)
    records = [_closed(pool, _declared(pool, r)) for r in declared]
    index: dict[CaseKey, list[int]] = {}
    for i, r in enumerate(records):
        if r.outcome is not CaseOutcome.BLOCKED:
            a, b = (pool.minima[m] for m in r.minima)
            index.setdefault(pair_key(a, b), []).append(i)
    candidates = list(_candidates(pool, discoveries))
    for c in candidates:
        ma, mb = c.start, c.end
        if ma is None or mb is None or (key := pair_key(ma, mb)) in index:
            continue
        record = _auto(pool, c, ma, mb)
        if record is not None:
            index[key] = [len(records)]
            records.append(_closed(pool, record))
    for c in candidates:
        if c.start is not None and c.end is not None:
            for i in index.get(pair_key(c.start, c.end), []):
                if records[i].outcome is None:
                    records[i] = _lend(pool, records[i], c)
    return records
