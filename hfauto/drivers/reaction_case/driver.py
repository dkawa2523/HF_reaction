"""ReactionCaseDriver (design §7.3): decide → action until a terminal decision. Backends are
seen only through ``hfauto.backends.protocols``; the stage builds the CaseRuntime."""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from hfauto.chemistry.classification import finalize, split
from hfauto.chemistry.identity import permutation_invariant_rmsd
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.ids import path_token
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
from hfauto.drivers.reaction_case import actions
from hfauto.drivers.reaction_case.state import CasePolicy, CaseState, decide

if TYPE_CHECKING:
    from hfauto.backends.protocols import PathEngine, QMEngine, SaddleRefiner


@dataclass(frozen=True)
class CaseRuntime:
    """``minima`` (DFT minimum id → record, optimized geometry) and ``species`` grow as cases
    register basins; ``resolve`` opens a FileRef of the run (trajectories, Hessians)."""

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
    deadline: Callable[[float], Deadline]
    resolve: Callable[[FileRef], Path]
    minima: dict[str, tuple[MinimumRecord, Geometry]]
    species: dict[str, SpeciesRecord]


@dataclass(frozen=True)
class CaseResult:
    reaction: ReactionRecord
    children: tuple[ReactionRecord, ...]  # split cases still to drive
    artifacts: tuple[Artifact, ...]  # reaction, calculation, new species and minima


def _endpoint(rt: CaseRuntime, minimum_id: str, species_id: str) -> np.ndarray:
    """The basin's optimized structure, relabelled into a member species' atom order."""
    record, geometry = rt.minima[minimum_id]
    xyz = rt.load_xyz(geometry)
    if species_id == record.species_id:  # the representative itself
        return np.asarray(xyz.coords, dtype=float)
    own = rt.load_xyz(rt.species[species_id].geometry)
    _, perm = permutation_invariant_rmsd(xyz.symbols, own.coords, xyz.coords)
    return np.asarray(xyz.coords, dtype=float)[perm]


def open_case(case: ReactionRecord, rt: CaseRuntime, policy: CasePolicy, deadline: Deadline,
              folder: Path, log: Callable[[dict[str, object]], None]) -> actions.Ctx:
    raw = (_endpoint(rt, case.minima[0], case.endpoints[0]),
           _endpoint(rt, case.minima[1], case.endpoints[1]))
    first = rt.species[case.endpoints[0]]
    return actions.Ctx(
        case=case, rt=rt, policy=policy, deadline=deadline, folder=folder,
        symbols=list(first.geometry.symbols), charge=first.charge,
        multiplicity=first.multiplicity, raw=raw, ends=(raw[0], align_mapped(raw[0], raw[1])),
        energies=(rt.minima[case.minima[0]][0].energy_hartree,
                  rt.minima[case.minima[1]][0].energy_hartree),
        log=log,
    )


def _artifacts(record: ReactionRecord, work: actions.Work) -> tuple[Artifact, ...]:
    out = [Artifact(artifact_id=k, type=ArtifactType.CALCULATION, payload=ev)
           for k, ev in work.calcs.items()]
    out += [Artifact(artifact_id=f"species_{s.species_id}", type=ArtifactType.SPECIES, payload=s)
            for s in work.species.values()]
    out += [Artifact(artifact_id=m.minimum_id, type=ArtifactType.MINIMUM, payload=m,
                     parents=(m.opt_calc, m.freq_calc)) for m in work.minima.values()]
    out.append(Artifact(artifact_id=record.reaction_id, type=ArtifactType.REACTION,
                        payload=record, parents=tuple(m for m in record.minima if m)))
    return tuple(out)


def drive_case(case: ReactionRecord, rt: CaseRuntime, policy: CasePolicy) -> CaseResult:
    """Loop decide → action; a re-run replays finished jobs from the JobStore."""
    deadline = rt.deadline(policy.walltime_s)
    folder = rt.case_dir / path_token(case.reaction_id)
    folder.mkdir(parents=True, exist_ok=True)
    log_path = folder / "log.jsonl"
    log_path.write_text("", encoding="utf-8")  # written only, never read back

    def log(entry: dict[str, object]) -> None:
        with log_path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(entry, sort_keys=True) + "\n")

    minima = tuple(rt.minima[m][0] if m in rt.minima else None for m in case.minima)
    state = CaseState(minima=(minima[0], minima[1]))
    ctx: actions.Ctx | None = None
    while True:
        state = replace(state, expired=deadline.expired())
        decision = decide(case, state, policy)
        log({"action": decision.action.value, "reason": decision.reason,
             "outcome": decision.outcome.value if decision.outcome else None})
        if decision.outcome is not None:
            break
        ctx = ctx or open_case(case, rt, policy, deadline, folder, log)
        state = actions.HANDLERS[decision.action](ctx, state, decision)
    work = ctx.work if ctx is not None else actions.Work()
    record = finalize(case, decision, barrier=state.screen, claim=state.claim,
                      connection=work.connection)
    record = record.model_copy(update={"log": f"cases/{folder.name}/log.jsonl"})
    children: tuple[ReactionRecord, ...] = ()
    if record.outcome is CaseOutcome.MULTI_STEP and work.intermediate is not None:
        children = split(record, *work.intermediate)
    return CaseResult(record, children, _artifacts(record, work))
