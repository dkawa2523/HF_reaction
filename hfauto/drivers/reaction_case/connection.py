"""The connection and intermediate actions of a reaction case (design §7.3): QRC from a
validated TS, assignment of its sides to DFT basins, and the wells that split a case. A well
that becomes an endpoint is the structure the case reached, in its atom order; a case is
judged at the granularity its ends differ in (``_key``)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from typing import Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry.gates import connection
from hfauto.chemistry.geometry import declared_coordinate
from hfauto.chemistry.identity import IMAGE_A, carry, mapped_equivalent, periodic_nearest
from hfauto.chemistry.modes import BOUNDS_A, displace
from hfauto.chemistry.xyz import composition_key
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.ids import species_id
from hfauto.core.records import ConnectionClaim, MinimumRecord, SpeciesRecord
from hfauto.drivers.minimum import calc_id, relax_to_minimum
from hfauto.drivers.reaction_case.actions import Ctx, Profile, mode_amplitude, peak_seed
from hfauto.drivers.reaction_case.state import CaseState, Decision

QRC_RETRY_FACTOR = 2.0  # the second QRC amplitude, capped at BOUNDS_A[1]


class Reached(NamedTuple):
    """A structure the case reached (a QRC side or a relaxed well, in the case atom order): its
    DFT basin, the optimization that reached it and, when the case registered it, the case's
    species there (None: the basin was found in the Registry)."""

    record: MinimumRecord
    opt: Evidence
    species: SpeciesRecord | None = None


def _state(record: MinimumRecord) -> str:
    return f"{record.composition_id}/{record.state_label}"


def _end_records(ctx: Ctx) -> tuple[MinimumRecord, MinimumRecord]:
    return ctx.rt.minima[ctx.case.minima[0]][0], ctx.rt.minima[ctx.case.minima[1]][0]


def _key(ctx: Ctx, record: MinimumRecord) -> str:
    """What a minimum counts as in this case. A case whose ends differ in chemical state
    (composition, state label) is judged by state, as thermo prices it (Curtin–Hammett:
    the basins of a state interconvert faster than it reacts); a same-state case (a declared
    torsion, a conformer change, a degenerate rearrangement) by basin."""
    a, b = _end_records(ctx)
    return _state(record) if _state(a) != _state(b) else record.basin_id


def _as_end(ctx: Ctx, well: MinimumRecord, end: MinimumRecord, path: Sequence[float]) -> bool:
    """G5-P2: another basin of the end's state, within a resolution of it and with no hill of a
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
    runtime and the case's artifacts; a basin already known keeps its representative geometry."""
    rt = ctx.rt
    ctx.work.species[species.species_id] = rt.species[species.species_id] = species
    ctx.work.minima[record.minimum_id] = record
    known = rt.minima.get(record.minimum_id)
    rt.minima[record.minimum_id] = (record, species.geometry if known is None else known[1])
    return record


def _register(ctx: Ctx, coords: np.ndarray, name: str,
              source: Literal["connection", "intermediate"], opt: Evidence | None = None
              ) -> Reached | None:
    """relax_to_minimum (from ``opt`` when it has converged already, not optimized again) →
    Registry: a known basin, a new one (or the one a push reached) with the case's species as a
    member, or None (no minimum)."""
    rt = ctx.rt
    out = relax_to_minimum(ctx.mol(coords), rt.method, rt.qm, known=rt.registry, opt=opt,
                           deadline=ctx.deadline, gates=ctx.rules.gates)
    if out.status == "known" and out.known_basin is not None and out.opt is not None:
        return Reached(ctx.record(out.known_basin), out.opt)
    if out.status not in ("minimum", "soft_minimum") or out.opt is None or out.freq is None:
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
    completes the case, row 5), so the member's id is unique."""
    if well.species is not None:
        return ctx.rt.minima[well.record.minimum_id][0], well.species
    species = _species(ctx, well.opt, "intermediate", "intermediate")
    ctx.keep(well.opt)
    record = ctx.rt.registry.join(well.record.basin_id, species.species_id)
    return _enter(ctx, record, species), species


def _assign(ctx: Ctx, side: Evidence, x: np.ndarray, name: str) -> Reached | None:
    """Registry.find of the side at ``x`` (its own structure, or its image's); a torsional case
    falls back to the nearest declared dihedral (CH-04); else the converged side is registered."""
    basin = ctx.rt.registry.find(side, coords=x)
    if basin is not None:
        return Reached(ctx.record(basin), side)
    terms = ctx.case.coordinate
    if ctx.case.torsional and terms and all(t.kind == "dihedral" for t in terms):
        values = [declared_coordinate(terms, end) for end in ctx.raw]
        end = ctx.case.minima[periodic_nearest(declared_coordinate(terms, x), values)]
        return Reached(ctx.rt.minima[end][0], side)
    return _register(ctx, x, name, "connection", opt=side)


def _sides(ctx: Ctx, freq: Evidence, starts: tuple[np.ndarray, np.ndarray], attempt: int
           ) -> tuple[tuple[Evidence, Evidence], tuple[np.ndarray, np.ndarray]] | None:
    """The QRC sides optimized from the TS Hessian and their final structures; None when one
    fails. At a symmetric TS the minus start is an exact image of the plus start (identity.carry
    within IMAGE_A), so on the invariant PES the minus optimum is the plus one's image: only the
    plus side runs and stands for both, the minus structure carried by that image. Two sides
    run at once (CaseRuntime.map)."""
    rt = ctx.rt
    image = carry(ctx.symbols, starts[1], starts[0], starts[0])[0] <= IMAGE_A
    if image:
        ctx.note(f"qrc{attempt}:minus_is_image")
    runs = rt.map(lambda y: rt.qm.optimize(ctx.mol(y), rt.method, init_hessian=freq,
                                           deadline=ctx.deadline), starts[:1 if image else 2])
    plus, minus = runs[0], runs[-1]
    if isinstance(plus, Failure) or isinstance(minus, Failure):
        ctx.note(f"qrc{attempt}:side_failed")
        return None
    x = ctx.coords(plus.final)
    y = carry(ctx.symbols, starts[1], starts[0], x)[1] if image else ctx.coords(minus.final)
    return (plus, minus), (x, y)


def connect(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """QRC: displace ± along the TS mode by an energy target, optimize, assign, gate."""
    freq = ctx.work.ts_freq
    attempt = state.connection_attempts + 1
    state = replace(state, connection_attempts=attempt)
    if freq is None or not freq.imaginary_modes:
        return replace(state, connection="failed")
    step = min(mode_amplitude(ctx, freq, 0) * QRC_RETRY_FACTOR ** (attempt - 1), BOUNDS_A[1])
    optimized = _sides(ctx, freq, displace(ctx.coords(freq.final),
                                           np.asarray(freq.imaginary_modes[0]), step), attempt)
    if optimized is None:
        return replace(state, connection="failed")
    sides, finals = optimized
    a, b = (_assign(ctx, s, x, f"qrc{attempt}_{i}")
            for i, (s, x) in enumerate(zip(sides, finals, strict=True)))
    first, second = _keys(ctx, (a, b), freq.energy_hartree)
    distinct = mapped_equivalent(ctx.symbols, *finals) if ctx.case.degenerate else True
    bonds = tuple(topology.bonds(ctx.symbols, x) for x in (*ctx.ends, finals[1], finals[0]))
    ends = frozenset(_key(ctx, r) for r in _end_records(ctx))
    gate, label = connection(freq, sides, (first, second), ends,
                             degenerate=ctx.case.degenerate, sides_distinct=distinct,
                             bond_sets=bonds)
    if label == "failed" or a is None or b is None:
        verdict = _rejected(gate.reasons, first == second, a, b)
        ctx.note(f"qrc{attempt}:{verdict}:{','.join(gate.reasons)}")
        return replace(state, connection=verdict)
    side_calcs = (calc_id(ctx.keep(sides[0])), calc_id(ctx.keep(sides[1])))
    if (two := _two_steps(ctx, state, label, (a, b), (first, second))) is not None:
        return two
    ctx.work.connection = ConnectionClaim(side_calcs=side_calcs,
                                          minima=(a.record.minimum_id, b.record.minimum_id))
    return replace(state, connection=label)


def _keys(ctx: Ctx, sides: tuple[Reached | None, Reached | None], ts: float
          ) -> tuple[str | None, str | None]:
    """The sides' case keys. A side in a new basin of the state of the end the other side
    reached is that end when it and the TS, the one hill between them, lie within a resolution
    of that end (G5-P2 on the QRC path)."""
    ends = {_key(ctx, e): e for e in _end_records(ctx)}
    keys = [None if r is None else _key(ctx, r.record) for r in sides]
    for well, mine, other in ((sides[0], keys[0], keys[1]), (sides[1], keys[1], keys[0])):
        end = None if other is None else ends.get(other)
        if (well is not None and end is not None and mine not in ends
                and _as_end(ctx, well.record, end, (ts,))):
            return other, other
    return keys[0], keys[1]


def _rejected(reasons: tuple[str, ...], one_key: bool, a: Reached | None, b: Reached | None
              ) -> Literal["same_basin", "same_state", "failed"]:
    """A connection that failed the gate. Both sides in one basin may be a displacement too
    small to leave it: worth a wider one (row 8). Both in one key but two basins (one state of a
    bond-changing case, or a basin split finer than the resolution) make a saddle of another
    process, never this case's (row 8 searches on). Anything else failed."""
    if a is None or b is None or not one_key:
        return "failed"
    if a.record.basin_id != b.record.basin_id:
        return "same_state"
    return "same_basin" if "sides_same_basin" in reasons else "failed"


def _two_steps(ctx: Ctx, state: CaseState, label: str, sides: tuple[Reached, Reached],
               keys: tuple[str | None, str | None]) -> CaseState | None:
    """GEN-05: a reassigned TS with exactly one side at an end (by its case key) makes the case
    two steps (row 5), the other side being the intermediate: a new chemical state for a
    bond-changing case, a new basin for a same-state one; the split child between those two
    validates this TS (a degenerate case: child 1). None for any other connection."""
    ends = [_key(ctx, r) for r in _end_records(ctx)]
    inside = [ends.index(k) for k in keys if k in ends]
    if label != "reassigned" or len(inside) != 1 or state.claim is None:
        return None
    well = next(r for r, k in zip(sides, keys, strict=True) if k not in ends)
    ctx.note(f"qrc{state.connection_attempts}:end{inside[0]}_to_new_basin:"
             f"{well.record.minimum_id}")
    ctx.work.intermediate = _member(ctx, well)
    ctx.work.split_ts = (inside[0] + 1, state.claim.saddle_calc)
    return replace(state, connection=None, claim=None, intermediate="distinct")


def _with_peak(ctx: Ctx, state: CaseState, name: str) -> CaseState:
    """The latest profile's highest peak joins the seeds."""
    v = state.screen
    source = "screen_hei" if v is not None and v.source == "screen" else "path_hei"
    seed = peak_seed(ctx, f"{name}_hei", source)
    return replace(state, seeds=(*state.seeds, seed)) if seed else state


def _past_the_well(ctx: Ctx, state: CaseState, path: Profile, name: str) -> CaseState:
    """The well relaxed into an endpoint: barrierless when no interior point rises a resolution
    above the higher end, else the highest peak seeds the saddle search."""
    v, e = state.screen, path.energies
    if v is not None and max(e[1:-1]) - max(e[0], e[-1]) < ctx.resolution:
        return replace(state, screen=v.model_copy(update={"verdict": "barrierless",
                                                          "reasons": ("well_is_endpoint",)}))
    return _with_peak(ctx, state, name)


def validate_intermediate(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """relax_to_minimum → Registry on the collapsed saddle (or soft TS whose QRC failed: its
    claim is withdrawn) or on the latest profile's lowest well; consumes its trigger. The well
    is an end when it has an end's case key, or (G5-P2) lies within a resolution of an end with
    no hill of a resolution between them on its profile; a collapsed saddle has no profile to
    show that. A failed relaxation is no chemical result: a well's profile goes on from its
    highest peak."""
    path = ctx.work.path if decision.reason == "path_intermediate" else None
    state = replace(state, ts_check=None, last_saddle=None, claim=None, connection=None)
    legs: tuple[Sequence[float], ...] = ()  # the profile from the well to end 0 and to end 1
    if path is not None:
        e = path.energies
        w = min(profile.interior_maxima([-v for v in e], ctx.resolution), key=e.__getitem__)
        x, legs = path.frames[w], (e[:w + 1], e[w:])
    elif ctx.work.saddle is not None:
        x = ctx.coords(ctx.work.saddle.final)
    else:
        return replace(state, intermediate="same_as_endpoint")
    name = f"int{state.saddle_attempts}_{state.path_runs}"
    well = _register(ctx, x, name, "intermediate")
    if well is None:
        ctx.note("int:relax_failed")
        state = replace(state, intermediate="relax_failed")
        return state if path is None else _with_peak(ctx, state, name)
    ends = _end_records(ctx)
    if _key(ctx, well.record) in {_key(ctx, e) for e in ends} or any(
            _as_end(ctx, well.record, e, leg) for e, leg in zip(ends, legs, strict=False)):
        state = replace(state, intermediate="same_as_endpoint")
        return state if path is None else _past_the_well(ctx, state, path, name)
    ctx.work.intermediate = _member(ctx, well)
    return replace(state, intermediate="distinct")
