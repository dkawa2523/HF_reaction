"""minima stage (design §4.1 #3, §7.2, §8.2): selection → relax_to_minimum → Registry.

``level: screen`` relaxes every input species; ``level: dft`` refines a selection
(``chemistry.selection``) of the reacting compositions, started from the screen minima's
optimized structures. Its entry takes DFT single points at low-level structures, once per
structure: at both ends and the TS of every low-level edge (a product discovery), which
``selection.admit`` admits or records again as a negative with its reason, and at the
conformers of a crowded state (``selection.rerank``). The jobs are ordered with the species a
screen minimum holds first, each group in species-id order. Their opts run at once
(StageRuntime.thread_map), then at once the freq (and mode following) of each opt that neither a
known basin nor an earlier opt holds (one identity criterion: identity.assign); the jobs are
then registered one by one in that order, and a job that falls into a basin registered before
it (mirror image included) is ``known`` and runs no freq job, as when they run one by one. An
exception in a species' jobs fails only its minimum (StageRuntime.contain). An exact permutation
or mirror image of an earlier start (identity.is_image) runs no job: the PES is invariant under
both, so it joins that start's basin when that start relaxed straight into a minimum. A saddle whose ±
displacements reach two distinct minima gives two ``mode_follow`` species at the driver's side
minima, registered by identity with no new job, and a ``mode_follow`` discovery between them
(its ends, in the saddle's atom order); the saddle's own species joins the side basin nearer
its input structure (note ``endpoint_was_saddle``).
"""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import ClassVar, Literal, cast

import numpy as np
from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.identity import assign, is_image, permutation_invariant_rmsd
from hfauto.chemistry.selection import (
    Edge,
    admit,
    reacting_candidates,
    rerank,
    select_for_refinement,
)
from hfauto.chemistry.topology import fragments, state_label
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FailureKind, Geometry
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, level_mismatches
from hfauto.core.records import ArtifactType, DiscoveryRecord, MinimumRecord, Payload, SpeciesRecord
from hfauto.drivers import minimum as driver
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

T = ArtifactType
_TIE_A = 0.05  # side RMSDs this close to the input count as equally near: the lower minimum wins
_MOVED = ("follow",)  # steps along a mode, whose sign an image's relaxation need not share


def _artifact(artifact_id: str, payload: Payload) -> Artifact:
    return Artifact(artifact_id=artifact_id, type=ArtifactType(payload.kind), payload=payload)


class SelectConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    include: Literal["window", "all"] = "window"
    per_state: int = 3
    window_kcal: float = 6.0
    rerank_top: int = 8  # a crowded state's lowest by screen energy, scored by single points
    max_edges: int = 6  # admitted low-level edges (pairs of ends) per composition


class EngineMethod(BaseModel):
    model_config = ConfigDict(extra="forbid")
    engine: str
    method: str


class MinimaConfig(StageConfig):
    level: Literal["screen", "dft"]
    engine: str
    method: str
    select: SelectConfig = SelectConfig()
    init_hessian: EngineMethod | None = None  # low-level freq at the start of 2+ fragments
    mode_follow: int = 2


@dataclass(frozen=True)
class _Job:
    species: SpeciesRecord
    start: Geometry


_Relaxed = tuple[_Job, driver.MinimumOutcome]
_Side = tuple[SpeciesRecord, MinimumRecord]


def _image_groups(jobs: Sequence[_Job], load_xyz: Callable[[Geometry], XYZ]
                  ) -> list[list[_Job]]:
    """Jobs of one composition (charge and multiplicity included) and atom order whose starts
    are exact permutation or mirror images (identity.is_image), in input order."""
    groups: list[tuple[XYZ, list[_Job]]] = []
    for job in jobs:
        x = load_xyz(job.start)
        group = next((members for ref, members in groups
                      if members[0].species.composition_id == job.species.composition_id
                      and ref.symbols == x.symbols
                      and is_image(x.symbols, ref.coords, x.coords)), None)
        if group is None:
            groups.append((x, [job]))
        else:
            group.append(job)
    return [members for _, members in groups]


@dataclass
class _Run:
    rt: StageRuntime
    cfg: MinimaConfig
    qm: QMEngine
    method: MethodSpec
    registry: driver.Registry
    calcs: dict[str, Evidence] = field(default_factory=dict)
    minima: dict[str, MinimumRecord] = field(default_factory=dict)  # latest record per basin
    species: list[SpeciesRecord] = field(default_factory=list)  # mode_follow sides
    extra: list[Artifact] = field(default_factory=list)  # discoveries and failed minima
    # single points by (structure fingerprint, charge, multiplicity)
    sps: dict[tuple[str, int, int], Evidence | Failure] = field(default_factory=dict)
    records: dict[str, MinimumRecord | None] = field(default_factory=dict)  # per species id
    history: dict[str, list[str]] = field(default_factory=dict)  # diagnostics.json
    early: dict[str, _Relaxed] = field(default_factory=dict)  # new basins relaxed at once

    def molecule(self, species: SpeciesRecord, geometry: Geometry) -> Molecule:
        return Molecule(self.rt.load_xyz(geometry), species.charge, species.multiplicity)

    def keep(self, *evidence: Evidence | None) -> None:
        self.calcs.update((driver.calc_id(ev), ev) for ev in evidence if ev is not None)

    def relax_all(self, jobs: Sequence[_Job]) -> list[Evidence | Failure]:
        """The opts of ``jobs`` at once; then at once the relaxation of each new basin
        (``new_basins``), kept in ``early`` for ``settle``."""
        opts = self.rt.thread_map(self.optimize, jobs)
        new = self.new_basins(opts)
        relaxed = self.rt.thread_map(lambda i: self.relax(jobs[i], opts[i]), new)
        self.early.update((jobs[i].species.species_id, done)
                          for i, done in zip(new, relaxed, strict=True))
        return opts

    def new_basins(self, opts: Sequence[Evidence | Failure]) -> list[int]:
        """The opts that neither a known basin nor an earlier one of them holds: their freq jobs
        run at once. A later opt in one of their basins waits for the registration (known)."""
        seen: dict[str, tuple[Evidence, np.ndarray]] = {}
        for i, opt in enumerate(opts):
            if isinstance(opt, Failure) or self.registry.find(opt) is not None:
                continue
            x = np.asarray(self.rt.load_xyz(opt.final).coords, dtype=float)
            alike = {k: (y, ev.energy_hartree) for k, (ev, y) in seen.items()
                     if (ev.final.symbols, ev.level) == (opt.final.symbols, opt.level)}
            if assign(opt.final.symbols, x, opt.energy_hartree, alike) is None:
                seen[str(i)] = (opt, x)
        return [int(i) for i in seen]

    def optimize(self, job: _Job) -> Evidence | Failure:
        """The job's opt, from the low-level Hessian of a start of 2+ fragments; contained."""
        def opt() -> Evidence | Failure:
            mol = self.molecule(job.species, job.start)
            return self.qm.optimize(mol, self.method, init_hessian=self.init_hessian(mol))

        return self.rt.contain(job.species.species_id, opt)

    def relax(self, job: _Job, opt: Evidence | Failure) -> _Relaxed:
        """relax_to_minimum from the job's opt (``known``: no freq job); contained."""
        def relaxed(opt: Evidence) -> driver.MinimumOutcome:
            return driver.relax_to_minimum(
                self.molecule(job.species, job.start), self.method, self.qm, known=self.registry,
                opt=opt, max_mode_follow=self.cfg.mode_follow, gates=self.rt.policy,
                load_xyz=self.rt.load_xyz, resolve=self.rt.resolve)

        out = opt if isinstance(opt, Failure) else self.rt.contain(job.species.species_id,
                                                                   lambda: relaxed(opt))
        if isinstance(out, Failure):
            history = ("opt", f"failed:{out.kind.value}")
            return job, driver.MinimumOutcome("failed", None, None, history, failure=out)
        return job, out

    def settle(self, job: _Job, opt: Evidence | Failure) -> _Relaxed:
        """Register the job's relaxation (``early``, or now: a basin registered before it is
        known); at a TS candidate, its two side minima, the join of the saddle species to the
        nearer side and a mode_follow discovery. The Registry is shared: not contained."""
        done = self.early.pop(job.species.species_id, None) or self.relax(job, opt)
        self.register(done)
        sides = self.sides(done)
        if len(sides) == 2:
            self.join_nearer_side(job, sides)
            found = self.discovery(done, sides)
            if found is not None:
                self.extra.append(_artifact(found.discovery_id, found))
        return done

    def process(self, job: _Job) -> _Relaxed:
        """opt, relax and register one job."""
        return self.settle(job, self.optimize(job))

    def init_hessian(self, mol: Molecule) -> Evidence | None:
        """A low-level freq at the start of 2+ fragments, an unrelaxed seed's included: the
        engine writes it as its positive-definite model (the declared default)."""
        spec = self.cfg.init_hessian
        if spec is None or len(fragments(mol.xyz.symbols, mol.xyz.coords)) < 2:
            return None
        low = cast(QMEngine, self.rt.engine(Capability.QM, spec.engine))
        freq = low.frequencies(mol, self.rt.method(spec.method))
        return None if isinstance(freq, Failure) else freq

    def energy(self, job: _Job) -> Evidence | Failure:
        """A single point at ``job.start``, once per structure and electronic state."""
        key = (job.start.fingerprint, job.species.charge, job.species.multiplicity)
        if key not in self.sps:
            self.sps[key] = ev = self.rt.contain(job.species.species_id, lambda: self.qm.energy(
                self.molecule(job.species, job.start), self.method))
            self.keep(None if isinstance(ev, Failure) else ev)
        return self.sps[key]

    def single_points(self, jobs: Sequence[_Job]) -> dict[str, float]:
        done = {j.species.species_id: self.energy(j) for j in jobs}
        return {sid: ev.energy_hartree for sid, ev in done.items() if isinstance(ev, Evidence)}

    def admitted(self, edges: Sequence[DiscoveryRecord], jobs: Mapping[str, _Job]
                 ) -> list[DiscoveryRecord]:
        """The edges ``selection.admit`` admits, on single points at their ends' starts and at
        their TS, taken only for the edges it asks; each other one is recorded again as a
        negative with its reason. Verdicts and heights go to diagnostics.json."""
        def point(job: _Job) -> float | str:
            ev = self.energy(job)
            return ev.kind.value if isinstance(ev, Failure) else ev.energy_hartree

        def state(job: _Job) -> str:  # of the structure the single point is taken at
            x = self.rt.load_xyz(job.start)
            return f"{job.species.composition_id}|{state_label(x.symbols, x.coords)}"

        def edge(d: DiscoveryRecord) -> Edge:
            source, product = jobs[d.source_species or ""], jobs[d.product_species or ""]
            ts = () if d.ts is None else (_Job(source.species, d.ts),)
            return Edge(d.discovery_id, source.species.composition_id, state(source),
                        state(product), lambda: tuple(point(j) for j in (source, product, *ts)),
                        first=(d.generation or 1) == 1)

        verdicts = admit([edge(d) for d in edges],
                         window_kcal=self.rt.policy.reaction_window_kcal,
                         per_composition=self.cfg.select.max_edges)
        for d in edges:
            reason, height = verdicts[d.discovery_id]
            self.history[d.discovery_id] = [reason or "admitted",
                                            *([] if height is None else [f"{height:.2f}"])]
            if reason is not None:
                self.extra.append(_artifact(d.discovery_id, d.model_copy(
                    update={"outcome": "negative", "reason": reason})))
        return [d for d in edges if verdicts[d.discovery_id][0] is None]

    def register(self, done: _Relaxed) -> MinimumRecord | None:
        """Registry.add (a known or identical basin is joined), at once: a later job in this
        basin is known (no freq job)."""
        job, outcome = done
        sid, record = job.species.species_id, None
        self.keep(outcome.opt, outcome.freq)
        self.history[sid] = [outcome.status, *outcome.history]
        if outcome.status in ("minimum", "known"):
            record = self.registry.add(outcome, job.species, tier=self.cfg.level)
            self.minima[record.basin_id] = record
        elif outcome.ts_candidate is None:  # a TS candidate is reported as a discovery instead
            failure = outcome.failure or Failure(kind=FailureKind.GATE_REJECTED,
                                                 reason=f"relaxed_to_{outcome.status}")
            self.extra.append(Artifact(
                artifact_id=f"min_{sid}_{self.rt.stage_id}", type=T.MINIMUM, status="failed",
                failure=failure, parents=(species_artifact_id(sid),)))
        self.records[sid] = record
        return record

    def join_image(self, job: _Job, first: _Relaxed) -> bool:
        """``job`` starts at an exact image of ``first``'s start and joins its basin with no job
        when ``first`` relaxed straight into a minimum or a known basin; False otherwise (a
        saddle endpoint, a failure, or a step along a mode)."""
        rep, outcome = first
        record = self.records.get(rep.species.species_id)
        moved = any(step.startswith(_MOVED) for step in outcome.history)
        if record is None or outcome.status not in ("minimum", "known") or moved:
            return False
        sid = job.species.species_id
        self.minima[record.basin_id] = self.records[sid] = self.registry.join(record.basin_id, sid)
        self.history[sid] = [f"image_of:{rep.species.species_id}"]
        return True

    def sides(self, done: _Relaxed) -> list[_Side]:
        """New species at the two minima the driver reached from a TS candidate, registered by
        identity with no job (the mirror side of a symmetric saddle joins the other's basin)."""
        parent, out = done[0].species, []
        for i, side in enumerate(done[1].ts_candidate or (), start=1):
            if side.opt is None:
                continue
            xyz = self.rt.load_xyz(side.opt.final)
            species = parent.model_copy(update={
                "species_id": f"{parent.species_id}_mf{i}", "geometry": side.opt.final,
                "source": "mode_follow", "state_label": state_label(xyz.symbols, xyz.coords),
                "energy_hartree": side.opt.energy_hartree, "level_key": side.opt.level.full_key()})
            self.species.append(species)
            record = self.register((_Job(species, side.opt.final), side))
            if record is not None:
                out.append((species, record))
        return out

    def join_nearer_side(self, job: _Job, sides: list[_Side]) -> None:
        """The species of a TS candidate joins the side basin nearer its input structure
        (permutation-invariant RMSD; the lower minimum when both are within _TIE_A)."""
        x = self.rt.load_xyz(job.start)
        near = [(permutation_invariant_rmsd(x.symbols, x.coords,
                                            self.rt.load_xyz(s.geometry).coords)[0], m)
                for s, m in sides]
        nearest = min(rmsd for rmsd, _ in near)
        basin = min((m for rmsd, m in near if rmsd - nearest <= _TIE_A),
                    key=lambda m: m.energy_hartree)
        sid = job.species.species_id
        record = self.registry.join(basin.basin_id, sid, "endpoint_was_saddle")
        self.minima[record.basin_id] = self.records[sid] = record
        self.history[sid].append(f"joined:{record.basin_id}")

    def discovery(self, done: _Relaxed, sides: list[_Side]) -> DiscoveryRecord | None:
        """An edge between the two side species (both in the saddle's atom order) with ts = the
        saddle; at the DFT tier ts_calc is its opt: a verified DFT saddle, validated directly by
        the case."""
        (source, a), (product, b) = sides
        parent, saddle = done[0].species.species_id, done[1].opt
        if saddle is None:
            return None
        return DiscoveryRecord(
            discovery_id=f"disc_mode_follow_{parent}", mechanism="mode_follow",
            outcome="product", source_species=source.species_id,
            product_species=product.species_id,
            ts=saddle.final, ts_calc=driver.calc_id(saddle) if self.cfg.level == "dft" else None,
            dE_act_kcal=(saddle.energy_hartree - a.energy_hartree) * HARTREE_TO_KCAL_MOL,
            dE_rxn_kcal=(b.energy_hartree - a.energy_hartree) * HARTREE_TO_KCAL_MOL,
        )


def _known(inputs: Manifest, cfg: MinimaConfig, rt: StageRuntime, qm: QMEngine,
           method: MethodSpec) -> driver.Registry:
    """The view's minima of this tier optimized with this engine and method."""
    pairs = []
    for m in inputs.records(T.MINIMUM, MinimumRecord):
        opt = inputs.evidence(m.opt_calc)
        same = not level_mismatches(method, opt.level, version_pin=opt.level.version)
        if m.tier == cfg.level and opt.engine == qm.name and same:
            pairs.append((m, opt.final))
    return driver.Registry(pairs, rt.load_xyz)


def _jobs(inputs: Manifest, cfg: MinimaConfig, run: _Run) -> list[_Job]:
    """Every species (screen, all), or the selection (chemistry.selection) of the reacting
    compositions (the others are noted not_reacting), each started from the optimized structure
    of the screen minimum it represents or else from its own geometry; in species-id order, the
    species a screen minimum holds first."""
    species = {s.species_id: s for s in inputs.records(T.SPECIES, SpeciesRecord)}
    if cfg.level == "screen" or cfg.select.include == "all":
        return [_Job(species[sid], species[sid].geometry) for sid in sorted(species)]
    screen = [m for m in inputs.records(T.MINIMUM, MinimumRecord)
              if m.tier == "screen" and m.species_id in species]
    starts = {m.species_id: inputs.evidence(m.opt_calc).final for m in screen}
    jobs = {sid: _Job(s, starts.get(sid, s.geometry)) for sid, s in species.items()}
    held = {s for m in screen for s in (m.species_id, *m.members)}
    edges = [d for d in inputs.records(T.DISCOVERY, DiscoveryRecord) if d.outcome == "product"
             and d.source_species in species and d.product_species in species]
    pool, idle = reacting_candidates(species, run.admitted(edges, jobs), screen,
                                     run.rt.system)
    run.history.update((composition, ["not_reacting"]) for composition in idle)
    sel = cfg.select
    chosen = rerank(select_for_refinement(pool, per_state=sel.rerank_top,
                                          window_kcal=sel.window_kcal),
                    lambda crowd: run.single_points([jobs[c.species_id] for c in crowd]),
                    sel.per_state, sel.window_kcal)
    return sorted((jobs[c.species_id] for c in chosen),
                  key=lambda j: (j.species.species_id not in held, j.species.species_id))


class MinimaStage:
    spec: ClassVar[StageSpec] = StageSpec(  # discovery is optional input (none before explore)
        "minima", MinimaConfig, (T.SPECIES,), (T.CALCULATION, T.MINIMUM, T.SPECIES, T.DISCOVERY))

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(MinimaConfig, config)
        method = rt.method(cfg.method)
        qm = cast(QMEngine, rt.engine(Capability.QM, cfg.engine))
        run = _Run(rt, cfg, qm, method, _known(inputs, cfg, rt, qm, method))
        groups = _image_groups(_jobs(inputs, cfg, run), rt.load_xyz)
        firsts = [first for first, *_ in groups]
        for (first, *images), opt in zip(groups, run.relax_all(firsts), strict=True):
            done = run.settle(first, opt)
            for image in images:
                if not run.join_image(image, done):
                    run.process(image)
        (rt.stage_dir / "diagnostics.json").write_text(json.dumps(run.history, indent=1))
        return [*(_artifact(k, ev) for k, ev in run.calcs.items()),
                *(_artifact(species_artifact_id(s.species_id), s) for s in run.species),
                *(_artifact(m.minimum_id, m) for m in run.minima.values()), *run.extra]
