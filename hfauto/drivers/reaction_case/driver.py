"""ReactionCaseDriver (design §7.3): decide → action until a terminal decision. Backends are
seen only through ``hfauto.backends.protocols``; the stage builds the CaseRuntime."""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from hfauto.chemistry.classification import finalize, split
from hfauto.chemistry.identity import basin_coords
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import Evidence, FileRef, Geometry
from hfauto.core.ids import path_token, species_artifact_id
from hfauto.core.manifest import Artifact
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case import actions, connection, paths
from hfauto.drivers.reaction_case.state import Action, CaseRules, CaseState, Decision, decide

if TYPE_CHECKING:
    from hfauto.backends.protocols import PathEngine, QMEngine, SaddleRefiner


@dataclass(frozen=True)
class CaseRuntime:
    """``minima`` (DFT minimum id → record, optimized geometry) and ``species`` grow as cases
    register basins, ``calcs`` (calculation id → Evidence) as cases finish; ``resolve`` opens a
    FileRef of the run (trajectories, Hessians); ``map`` (StageRuntime.thread_map) runs
    independent jobs at once, in input order."""

    qm: QMEngine
    saddle: SaddleRefiner
    path: PathEngine
    screen_qm: QMEngine | None
    screen_path: PathEngine | None
    method: MethodSpec
    screen_method: MethodSpec | None
    registry: Registry
    load_xyz: Callable[[Geometry], XYZ]
    file_ref: Callable[[Path], FileRef]
    case_dir: Path  # cases/; each case writes below case_dir/<reaction_id>
    resolve: Callable[[FileRef], Path]
    map: Callable[[Callable[[Any], Any], Sequence[Any]], list[Any]]
    minima: dict[str, tuple[MinimumRecord, Geometry]]
    species: dict[str, SpeciesRecord]
    calcs: dict[str, Evidence]


@dataclass(frozen=True)
class CaseResult:
    reaction: ReactionRecord
    children: tuple[ReactionRecord, ...]  # split cases still to drive
    artifacts: tuple[Artifact, ...]  # reaction, calculation, new species and minima


def _endpoint(rt: CaseRuntime, minimum_id: str, species_id: str) -> np.ndarray:
    """The basin's optimized structure in a member species' atom order and handedness."""
    record, geometry = rt.minima[minimum_id]
    xyz = rt.load_xyz(geometry)
    if species_id == record.species_id:  # the representative itself
        return np.asarray(xyz.coords, dtype=float)
    own = rt.load_xyz(rt.species[species_id].geometry)
    return basin_coords(xyz.symbols, xyz.coords, own.coords)


def open_case(case: ReactionRecord, rt: CaseRuntime, rules: CaseRules, deadline: Deadline,
              folder: Path, log: Callable[[dict[str, object]], None]) -> actions.Ctx:
    """The case context; a given DFT stationary point (``ts_calc``) is its saddle. An
    association's reactant energy is its separated monomers' sum."""
    first, reactant = rt.species[case.endpoints[0]], case.monomers or case.minima[:1]
    # an association's complex that relaxed into the adduct ends at its own structure (hypotheses)
    collapsed = bool(case.monomers) and case.minima[0] == case.minima[1]
    raw = (np.asarray(rt.load_xyz(first.geometry).coords, dtype=float) if collapsed
           else _endpoint(rt, case.minima[0], case.endpoints[0]),
           _endpoint(rt, case.minima[1], case.endpoints[1]))
    return actions.Ctx(
        case=case, rt=rt, rules=rules, deadline=deadline, folder=folder,
        symbols=list(first.geometry.symbols), charge=first.charge,
        multiplicity=first.multiplicity, raw=raw, ends=(raw[0], align_mapped(raw[0], raw[1])),
        energies=(sum(rt.minima[m][0].energy_hartree for m in reactant),
                  rt.minima[case.minima[1]][0].energy_hartree),
        log=log, work=actions.Work(saddle=rt.calcs[case.ts_calc] if case.ts_calc else None),
    )


def _artifacts(record: ReactionRecord, work: actions.Work) -> tuple[Artifact, ...]:
    out = [Artifact(artifact_id=k, type=ArtifactType.CALCULATION, payload=ev)
           for k, ev in work.calcs.items()]
    out += [Artifact(artifact_id=species_artifact_id(s.species_id), type=ArtifactType.SPECIES,
                     payload=s) for s in work.species.values()]
    out += [Artifact(artifact_id=m.minimum_id, type=ArtifactType.MINIMUM, payload=m,
                     parents=(m.opt_calc, m.freq_calc)) for m in work.minima.values()]
    out.append(Artifact(artifact_id=record.reaction_id, type=ArtifactType.REACTION,
                        payload=record, parents=tuple(m for m in record.minima if m)))
    return tuple(out)


Handler = Callable[[actions.Ctx, CaseState, Decision], CaseState]
HANDLERS: dict[Action, Handler] = {
    Action.SCREEN: paths.screen,
    Action.REFINE_SADDLE: actions.refine_saddle,
    Action.VALIDATE_TS: actions.validate_ts,
    Action.FIND_PATH: paths.find_path,
    Action.CONNECT: connection.connect,
    Action.VALIDATE_INTERMEDIATE: connection.validate_intermediate,
}


def drive_case(case: ReactionRecord, rt: CaseRuntime, rules: CaseRules, deadline: Deadline
               ) -> CaseResult:
    """Loop decide → action until ``deadline`` (the hypothesis', shared by its split children);
    a re-run replays finished jobs from the JobStore."""
    folder = rt.case_dir / path_token(case.reaction_id)
    folder.mkdir(parents=True, exist_ok=True)
    log_path = folder / "log.jsonl"
    log_path.write_text("", encoding="utf-8")  # written only, never read back

    def log(entry: dict[str, object]) -> None:
        with log_path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(entry, sort_keys=True) + "\n")

    minima = tuple(rt.minima[m][0] if m in rt.minima else None for m in case.minima)
    state = CaseState(minima=(minima[0], minima[1]),
                      last_saddle="converged" if case.ts_calc else None)
    if case.ts_calc:  # validated first (row 10); a failed gate goes on to SCREEN
        log({"note": f"ts_calc:{case.ts_calc}"})
    ctx: actions.Ctx | None = None
    while True:
        state = replace(state, expired=deadline.expired())
        decision = decide(case, state, rules)
        log({"action": decision.action.value, "reason": decision.reason,
             "outcome": decision.outcome.value if decision.outcome else None})
        if decision.outcome is not None:
            break
        ctx = ctx or open_case(case, rt, rules, deadline, folder, log)
        state = HANDLERS[decision.action](ctx, state, decision)
    work = ctx.work if ctx is not None else actions.Work()
    record = finalize(case, decision, barrier=state.screen, claim=state.claim,
                      connection=work.connection)
    record = record.model_copy(update={"log": f"cases/{folder.name}/log.jsonl"})
    return CaseResult(record, _children(record, rt, ctx), _artifacts(record, work))


def _children(record: ReactionRecord, rt: CaseRuntime, ctx: actions.Ctx | None
              ) -> tuple[ReactionRecord, ...]:
    """The steps R→I and I→P of a multi-step case, I in the atom order of the structure this
    case reached; the one whose ends a validated TS of this case joins validates it again (a
    JobStore replay) instead of searching."""
    if record.outcome is not CaseOutcome.MULTI_STEP or ctx is None or not ctx.work.intermediate:
        return ()
    well, well_species = ctx.work.intermediate
    middle = _endpoint(rt, well.minimum_id, well_species.species_id)
    return split(record, well, well_species, (ctx.raw[0], middle, ctx.raw[1]),
                 ts=ctx.work.split_ts)
