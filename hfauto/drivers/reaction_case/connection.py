"""The connection and intermediate actions of a reaction case (design §7.3): the TS checks and
QRC from the TS, assignment of its sides to DFT basins, and the wells that split a case. A well
that becomes an endpoint is the structure the case reached, in its atom order; a case is
judged at the granularity its ends differ in (``_key``). A connection that gives no answer is one
failure token: it leaves neither claim nor connection, and the case searches on."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from typing import Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry.gates import connection
from hfauto.chemistry.geometry import declared_coordinate
from hfauto.chemistry.identity import mapped_equivalent, periodic_nearest
from hfauto.chemistry.modes import qrc_step
from hfauto.chemistry.profile import Profile
from hfauto.chemistry.xyz import composition_key
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.ids import species_id
from hfauto.core.records import ConnectionClaim, MinimumRecord, SpeciesRecord
from hfauto.drivers.minimum import MinimumOutcome, Relaxer, Settled, calc_id, relax_to_minimum
from hfauto.drivers.reaction_case.actions import Ctx, validate_ts
from hfauto.drivers.reaction_case.paths import judge, peak_seeds
from hfauto.drivers.reaction_case.state import CaseState


class Reached(NamedTuple):
    """A structure the case reached (a QRC side or a relaxed well, in the case atom order): its
    DFT basin, the optimization that reached it and, when the case registered it, the case's
    species there (None: the basin was found in the Registry)."""

    record: MinimumRecord
    opt: Evidence
    species: SpeciesRecord | None = None


def _state(record: MinimumRecord) -> str:
    return f"{record.composition_id}/{record.state_label}"


def _key(ctx: Ctx, record: MinimumRecord) -> str:
    """What a minimum counts as in this case. A case whose ends differ in chemical state
    (composition, state label) is judged by state, as thermo prices it (Curtin–Hammett:
    the basins of a state interconvert faster than it reacts); a same-state case (a declared
    torsion, a conformer change, a degenerate rearrangement) by basin."""
    a, b = ctx.end_records()
    return _state(record) if _state(a) != _state(b) else record.basin_id


def _as_end(ctx: Ctx, well: MinimumRecord, end: MinimumRecord, path: Sequence[float]) -> bool:
    """Another basin of the end's state, within a resolution of it and with no hill of a
    resolution between them on ``path`` (the energies joining them), is that end at the
    resolution: a basin split finer than the resolution is no intermediate."""
    e_well, e_end, res = well.energy_hartree, end.energy_hartree, ctx.resolution
    return (_state(well) == _state(end) and abs(e_well - e_end) < res
            and max(path) - max(e_well, e_end) < res)


def _species(ctx: Ctx, opt: Evidence, name: str,
             source: Literal["connection", "intermediate"]) -> SpeciesRecord:
    """The optimized structure of ``opt`` as a species of this case, in the case atom order."""
    return SpeciesRecord(
        species_id=species_id(ctx.case.reaction_id, name),
        composition_id=composition_key(ctx.symbols, ctx.charge, ctx.multiplicity),
        charge=ctx.charge, multiplicity=ctx.multiplicity, geometry=opt.final, source=source,
        state_label=topology.state_label(ctx.symbols, ctx.coords(opt.final)),
        energy_hartree=opt.energy_hartree, level_key=opt.level.full_key(),
    )


def _enter(ctx: Ctx, record: MinimumRecord, species: SpeciesRecord) -> MinimumRecord:
    """``species`` and ``record`` (the Registry's, ``species`` among its members) join the
    runtime and the case's artifacts."""
    ctx.work.species[species.species_id] = ctx.rt.species[species.species_id] = species
    ctx.work.minima[record.minimum_id] = record
    return record


def _registered(ctx: Ctx, out: MinimumOutcome, name: str,
                source: Literal["connection", "intermediate"]) -> Reached | None:
    """A known basin as it is; a new minimum (its opt and freq kept) registered with the case's
    species as its member; anything else None, noted."""
    rt = ctx.rt
    if out.status == "known" and out.known_basin is not None and out.opt is not None:
        return Reached(rt.registry.basin(out.known_basin), out.opt)
    if out.status != "minimum" or out.opt is None or out.freq is None:
        ctx.note(f"{name}:{out.status}:{out.failure.reason if out.failure else ''}")
        return None
    species = _species(ctx, out.opt, name, source)
    record = _enter(ctx, rt.registry.add(out, species, tier="dft"), species)
    ctx.keep(out.opt)
    ctx.keep(out.freq)
    return Reached(record, out.opt, species)


def _member(ctx: Ctx, well: Reached) -> tuple[MinimumRecord, SpeciesRecord]:
    """The intermediate as the case reached it: its own species when the case registered the
    basin, else the structure reached, joined to the known basin as a member without a job. A
    split child's endpoint (driver._endpoint → identity.basin_coords) then continues the case's
    atom labelling, never the representative's arbitrary one. A case has one intermediate (it
    completes the case, row 2), so the member's id is unique."""
    if well.species is not None:
        return ctx.rt.registry.minima[well.record.minimum_id][0], well.species
    species = _species(ctx, well.opt, "intermediate", "intermediate")
    ctx.keep(well.opt)
    record = ctx.rt.registry.join(well.record.basin_id, species.species_id)
    return _enter(ctx, record, species), species


def _known(ctx: Ctx, opt: Evidence, x: np.ndarray) -> str | None:
    """The basin of a QRC side at ``x`` (its own structure, or its image's): Registry.find; a
    torsional case falls back to the end of the nearest declared dihedral (design §7.3)."""
    registry, terms = ctx.rt.registry, ctx.case.coordinate
    if (basin := registry.find(opt, coords=x)) is not None:
        return basin
    if ctx.case.torsional and terms and all(t.kind == "dihedral" for t in terms):
        values = [declared_coordinate(terms, end) for end in ctx.raw]
        end = ctx.case.minima[periodic_nearest(declared_coordinate(terms, x), values)]
        return registry.minima[end][0].basin_id
    return None


def _assign(ctx: Ctx, side: Settled, name: str) -> Reached | None:
    """A settled QRC side's basin: known (found again: the side before it may have registered
    it), else registered as a new minimum; a side still on a saddle is unassigned."""
    basin = side.basin or _known(ctx, side.opt, side.x)
    if basin is not None:
        return Reached(ctx.rt.registry.basin(basin), side.opt)
    out = MinimumOutcome("saddle" if side.tier == "saddle" else "minimum", side.opt, side.freq,
                         side.history, notes=side.notes)
    return _registered(ctx, out, name, "connection")


def validate_and_connect(ctx: Ctx, state: CaseState) -> CaseState:
    """The saddle's TS checks (actions.validate_ts), then QRC from an accepted TS."""
    state = validate_ts(ctx, state)
    if state.claim is None:
        return state
    return connect(ctx, state, ctx.work.calcs[state.claim.freq_calc])


def connect(ctx: Ctx, state: CaseState, freq: Evidence) -> CaseState:
    """QRC from the TS ``freq``: ± along its mode (modes.qrc_step), down to minima
    (minimum.Relaxer.descend: an image side carried, a side on a saddle taken down once more),
    assigned and gated. No answer is a failure token: the claim goes with it."""
    rt, x = ctx.rt, ctx.coords(freq.final)
    relaxer = Relaxer(ctx.mol(x), rt.method, rt.qm, rt.load_xyz, ctx.rules.gates, rt.resolve,
                      lambda opt, y: _known(ctx, opt, y), rt.map)
    sides = relaxer.descend(x, freq, qrc_step(freq, 0, ctx.symbols), both=True, name="qrc",
                            carried=True, past_saddles=True)
    if any(isinstance(s, Failure) for s in sides):
        ctx.note("qrc:side_failed")
        return replace(state, claim=None, connection=None)
    plus, minus = (s for s in sides if isinstance(s, Settled))
    if minus.opt is plus.opt:
        ctx.note("qrc:minus_is_image")
    a, b = (_assign(ctx, s, f"qrc1_{i}") for i, s in enumerate((plus, minus)))
    first, second = _keys(ctx, (a, b), freq.energy_hartree)
    distinct = mapped_equivalent(ctx.symbols, plus.x, minus.x) if ctx.case.degenerate else True
    exchange = ctx.change(), topology.bond_changes(ctx.symbols, plus.x, minus.x)
    ends = frozenset(_key(ctx, r) for r in ctx.end_records())
    gate, label = connection(freq, (plus.opt, minus.opt), (first, second), ends,
                             degenerate=ctx.case.degenerate, sides_distinct=distinct,
                             exchange=exchange)
    if label == "failed" or a is None or b is None:
        ctx.note(f"qrc:failed:{','.join(gate.reasons)}")
        return replace(state, claim=None, connection=None)
    side_calcs = (calc_id(ctx.keep(plus.opt)), calc_id(ctx.keep(minus.opt)))
    if (two := _two_steps(ctx, state, label, (a, b), (first, second))) is not None:
        return two
    ctx.work.connection = ConnectionClaim(side_calcs=side_calcs,
                                          minima=(a.record.minimum_id, b.record.minimum_id))
    return replace(state, connection=label)


def _keys(ctx: Ctx, sides: tuple[Reached | None, Reached | None], ts: float
          ) -> tuple[str | None, str | None]:
    """The sides' case keys. A side in a new basin of the state of the end the other side
    reached is that end when it and the TS, the one hill between them, lie within a resolution
    of that end (``_as_end`` on the QRC path)."""
    ends = {_key(ctx, e): e for e in ctx.end_records()}
    keys = [None if r is None else _key(ctx, r.record) for r in sides]
    for well, mine, other in ((sides[0], keys[0], keys[1]), (sides[1], keys[1], keys[0])):
        end = None if other is None else ends.get(other)
        if (well is not None and end is not None and mine not in ends
                and _as_end(ctx, well.record, end, (ts,))):
            return other, other
    return keys[0], keys[1]


def _two_steps(ctx: Ctx, state: CaseState, label: str, sides: tuple[Reached, Reached],
               keys: tuple[str | None, str | None]) -> CaseState | None:
    """A reassigned TS with exactly one side at an end (by its case key) makes the case
    two steps (row 2), the other side being the intermediate: a new chemical state for a
    bond-changing case, a new basin for a same-state one; the split child between those two
    validates this TS (a degenerate case: child 1). None for any other connection."""
    ends = [_key(ctx, r) for r in ctx.end_records()]
    inside = [ends.index(k) for k in keys if k in ends]
    if label != "reassigned" or len(inside) != 1 or state.claim is None:
        return None
    well = next(r for r, k in zip(sides, keys, strict=True) if k not in ends)
    ctx.note(f"qrc:end{inside[0]}_to_new_basin:{well.record.minimum_id}")
    ctx.work.intermediate = _member(ctx, well)
    ctx.work.split_ts = (inside[0] + 1, state.claim.saddle_calc)
    return replace(state, connection=None, claim=None, intermediate="distinct")


def _past_the_well(ctx: Ctx, state: CaseState, path: Profile, w: int, end: int, name: str
                   ) -> CaseState:
    """The well (node ``w``) relaxed into end ``end``: the case goes on with the profile from
    the other end to the well, judged like any profile (``paths.judge`` with its energy
    function, densified when barrierless). Barrierless closes the case; otherwise its highest
    peak seeds the saddle search."""
    sample = ctx.work.sample or ctx.points
    start, stop = (w, len(path.energies)) if end == 0 else (0, w + 1)
    part = Profile(path.frames[start:stop], path.energies[start:stop], path.source,
                   path.s2[start:stop])
    verdict, part, _ = judge(ctx, part, ctx.work.points[start:stop],
                             lambda new: sample([(i + start, x) for i, x in new]))
    ctx.note(f"{name}:end{end}:{verdict.verdict}")
    if verdict.verdict == "barrierless":
        return replace(state, screen=verdict)
    return replace(state, seeds=(*state.seeds, *peak_seeds(ctx, f"{name}_hei", part)))


def validate_intermediate(ctx: Ctx, state: CaseState) -> CaseState:
    """relax_to_minimum → Registry on the latest profile's lowest well (row 6); the claim of an
    earlier TS is dropped with its connection. The well is an end when it has an end's case key,
    or (``_as_end``) lies within a resolution of an end with no hill of a resolution between
    them on the profile. A failed relaxation is no chemical result: the profile goes on from its
    highest peak."""
    path, state = ctx.work.path, replace(state, claim=None, connection=None)
    if path is None:  # row 6 follows a DFT profile (CaseState.screen)
        return replace(state, intermediate="relax_failed")
    e, rt = path.energies, ctx.rt
    w = min(profile.interior_maxima([-v for v in e], ctx.resolution), key=e.__getitem__)
    name = f"int{state.attempts}_{state.path_runs}"
    out = relax_to_minimum(ctx.mol(path.frames[w]), rt.method, rt.qm, known=rt.registry,
                           gates=ctx.rules.gates, resolve=rt.resolve)
    well = _registered(ctx, out, name, "intermediate")
    if well is None:
        ctx.note("int:relax_failed")  # the latest profile's highest peak seeds
        return replace(state, intermediate="relax_failed",
                       seeds=(*state.seeds, *peak_seeds(ctx, f"{name}_hei")))
    ends, legs = ctx.end_records(), (e[:w + 1], e[w:])  # the profile from the well to each end
    for k, (end, leg) in enumerate(zip(ends, legs, strict=True)):
        if _key(ctx, well.record) == _key(ctx, end) or _as_end(ctx, well.record, end, leg):
            return _past_the_well(ctx, replace(state, intermediate="same_as_endpoint"), path,
                                  w, k, name)
    ctx.work.intermediate = _member(ctx, well)
    return replace(state, intermediate="distinct")
