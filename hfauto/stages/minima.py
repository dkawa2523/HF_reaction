"""minima stage (design §4.1 #3, §7.2, §8.2): selection → relax_to_minimum → Registry.

``level: screen`` relaxes every input species; ``level: dft`` refines a selection
(``chemistry.selection``) started from the screen minima's optimized structures. Jobs run
serially in species-id order and each is registered as soon as it is relaxed, so a later job
that falls into a registered basin (its mirror image included) is ``known`` and skips its freq
job (one identity criterion: identity.assign). A saddle whose ± displacements reach two
distinct minima gives two ``mode_follow`` species, relaxed and registered right after it, their
minima and a ``mode_follow`` discovery; the saddle's own species joins the side basin nearer its
input structure (note ``endpoint_was_saddle``).
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import ClassVar, Literal, cast

from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.identity import permutation_invariant_rmsd
from hfauto.chemistry.selection import Candidate, rerank, select_for_refinement
from hfauto.chemistry.topology import fragments, state_label
from hfauto.chemistry.xyz import Molecule
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FailureKind, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, level_mismatches
from hfauto.core.records import ArtifactType, DiscoveryRecord, MinimumRecord, Payload, SpeciesRecord
from hfauto.drivers import minimum as driver
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

T = ArtifactType
_TIE_A = 0.05  # side RMSDs this close to the input count as equally near: the lower minimum wins


def _artifact(artifact_id: str, payload: Payload) -> Artifact:
    return Artifact(artifact_id=artifact_id, type=ArtifactType(payload.kind), payload=payload)


class SelectConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    include: Literal["window", "all"] = "window"
    per_state: int = 3
    window_kcal: float = 6.0
    rerank_sp: bool = False
    rerank_top: int = 8


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


@dataclass
class _Run:
    rt: StageRuntime
    cfg: MinimaConfig
    qm: QMEngine
    method: MethodSpec
    registry: driver.Registry
    calcs: dict[str, Evidence] = field(default_factory=dict)
    minima: dict[str, MinimumRecord] = field(default_factory=dict)  # latest record per basin
    extra: list[Artifact] = field(default_factory=list)  # discoveries and failed minima
    records: dict[str, MinimumRecord | None] = field(default_factory=dict)  # per species id
    history: dict[str, list[str]] = field(default_factory=dict)  # diagnostics.json

    def molecule(self, species: SpeciesRecord, geometry: Geometry) -> Molecule:
        return Molecule(self.rt.load_xyz(geometry), species.charge, species.multiplicity)

    def keep(self, *evidence: Evidence | None) -> None:
        self.calcs.update((driver.calc_id(ev), ev) for ev in evidence if ev is not None)

    def settle(self, job: _Job) -> _Relaxed:
        """relax, then register at once: a later job in this basin is known (no freq job)."""
        done = self.relax(job)
        sid = job.species.species_id
        self.records[sid] = self.register(done)
        self.history[sid] = [done[1].status, *done[1].history]
        return done

    def relax(self, job: _Job) -> _Relaxed:
        mol = self.molecule(job.species, job.start)
        return job, driver.relax_to_minimum(
            mol, self.method, self.qm, known=self.registry, init_hessian=self.init_hessian(mol),
            max_mode_follow=self.cfg.mode_follow, gates=self.rt.policy, load_xyz=self.rt.load_xyz)

    def init_hessian(self, mol: Molecule) -> Evidence | None:
        spec = self.cfg.init_hessian
        if spec is None or len(fragments(mol.xyz.symbols, mol.xyz.coords)) < 2:
            return None
        low = cast(QMEngine, self.rt.engine(Capability.QM, spec.engine))
        freq = low.frequencies(mol, self.rt.method(spec.method))
        return None if isinstance(freq, Failure) else freq

    def single_points(self, jobs: Sequence[_Job]) -> dict[str, float]:
        runs = {j.species.species_id: self.qm.energy(self.molecule(j.species, j.start), self.method)
                for j in jobs}
        done = {sid: ev for sid, ev in runs.items() if not isinstance(ev, Failure)}
        self.keep(*done.values())
        return {sid: ev.energy_hartree for sid, ev in done.items()}

    def register(self, done: _Relaxed) -> MinimumRecord | None:
        job, outcome = done
        self.keep(outcome.opt, outcome.freq, *(outcome.ts_candidate or ()))
        if outcome.status in ("minimum", "soft_minimum", "known"):
            record = self.registry.add(outcome, job.species, tier=self.cfg.level)
            self.minima[record.basin_id] = record
            return record
        if outcome.ts_candidate is None:  # a TS candidate is reported as a discovery instead
            failure = outcome.failure or Failure(kind=FailureKind.GATE_REJECTED,
                                                 reason=f"relaxed_to_{outcome.status}")
            self.extra.append(Artifact(
                artifact_id=f"min_{job.species.species_id}_{self.rt.stage_id}", type=T.MINIMUM,
                status="failed", failure=failure, parents=(f"species_{job.species.species_id}",)))
        return None

    def sides(self, done: _Relaxed) -> list[_Job]:
        """New species at the two minima reached from a TS candidate."""
        parent, jobs = done[0].species, []
        for i, freq in enumerate(done[1].ts_candidate or (), start=1):
            xyz = self.rt.load_xyz(freq.start)
            species = parent.model_copy(update={
                "species_id": f"{parent.species_id}_mf{i}", "geometry": freq.start,
                "source": "mode_follow", "state_label": state_label(xyz.symbols, xyz.coords),
                "energy_hartree": freq.energy_hartree, "level_key": freq.level.full_key()})
            jobs.append(_Job(species, freq.start))
        return jobs

    def join_nearer_side(self, done: _Relaxed) -> None:
        """The species of a TS candidate joins the side basin nearer its input structure
        (permutation-invariant RMSD; the lower minimum when both are within _TIE_A)."""
        job, outcome = done
        if outcome.ts_candidate is None:
            return
        sid, x = job.species.species_id, self.rt.load_xyz(job.start)
        sides = [(permutation_invariant_rmsd(x.symbols, x.coords,
                                             self.rt.load_xyz(freq.start).coords)[0], record)
                 for i, freq in enumerate(outcome.ts_candidate, start=1)
                 if (record := self.records.get(f"{sid}_mf{i}")) is not None]
        if len(sides) != 2:
            return
        nearest = min(rmsd for rmsd, _ in sides)
        basin = min((m for rmsd, m in sides if rmsd - nearest <= _TIE_A),
                    key=lambda m: m.energy_hartree)
        record = self.registry.join(basin.basin_id, sid, "endpoint_was_saddle")
        self.minima[record.basin_id] = self.records[sid] = record
        self.history[sid].append(f"joined:{record.basin_id}")

    def discovery(self, done: _Relaxed) -> DiscoveryRecord | None:
        """source_minimum = side 1's minimum, product_species = side 2, ts = the saddle."""
        parent, saddle, freq = done[0].species.species_id, done[1].opt, done[1].freq
        a, b = self.records.get(f"{parent}_mf1"), self.records.get(f"{parent}_mf2")
        if a is None or b is None or saddle is None or freq is None:
            return None
        return DiscoveryRecord(
            discovery_id=f"disc_mode_follow_{parent}", source_minimum=a.minimum_id,
            mechanism="mode_follow", outcome="product", product_species=f"{parent}_mf2",
            ts=saddle.final, ts_imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
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


def _pool(inputs: Manifest, species: dict[str, SpeciesRecord]
          ) -> dict[str, tuple[Candidate, _Job]]:
    """Screen minima (from their optimized structure) plus discovery products and sources."""
    found = [d for d in inputs.records(T.DISCOVERY, DiscoveryRecord) if d.outcome == "product"]
    always = {d.source_minimum for d in found}
    always |= {d.product_species for d in found if d.product_species is not None}
    screen = [m for m in inputs.records(T.MINIMUM, MinimumRecord)
              if m.tier == "screen" and m.species_id in species]
    pool: dict[str, tuple[Candidate, _Job]] = {}
    for m in screen:
        kept = bool(always & {m.minimum_id, m.species_id, *m.members})
        candidate = Candidate(m.species_id, m.composition_id, m.state_label, m.energy_hartree,
                              always=kept)
        pool[m.species_id] = (candidate, _Job(species[m.species_id],
                                              inputs.evidence(m.opt_calc).final))
    held = {s for m in screen for s in (m.species_id, *m.members)}
    for sid in sorted((always & species.keys()) - held):  # products that skipped screen
        s = species[sid]
        pool[sid] = (Candidate(sid, s.composition_id, s.state_label, s.energy_hartree,
                               always=True), _Job(s, s.geometry))
    return pool


def _jobs(inputs: Manifest, cfg: MinimaConfig, run: _Run) -> list[_Job]:
    species = {s.species_id: s for s in inputs.records(T.SPECIES, SpeciesRecord)}
    if cfg.level == "screen" or cfg.select.include == "all":
        return [_Job(s, s.geometry) for s in species.values()]
    pool, sel = _pool(inputs, species), cfg.select
    chosen = select_for_refinement(
        [candidate for candidate, _ in pool.values()],
        per_state=sel.rerank_top if sel.rerank_sp else sel.per_state, window_kcal=sel.window_kcal)
    if sel.rerank_sp:
        energies = run.single_points([pool[c.species_id][1] for c in chosen if not c.always])
        chosen = rerank(chosen, energies, sel.per_state, sel.window_kcal)
    return [pool[c.species_id][1] for c in chosen]


class MinimaStage:
    spec: ClassVar[StageSpec] = StageSpec(  # discovery is optional input (none before explore)
        "minima", MinimaConfig, (T.SPECIES,), (T.CALCULATION, T.MINIMUM, T.SPECIES, T.DISCOVERY))

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(MinimaConfig, config)
        method = rt.method(cfg.method)
        qm = cast(QMEngine, rt.engine(Capability.QM, cfg.engine))
        run = _Run(rt, cfg, qm, method, _known(inputs, cfg, rt, qm, method))
        side_jobs: list[_Job] = []
        for job in sorted(_jobs(inputs, cfg, run), key=lambda j: j.species.species_id):
            done = run.settle(job)
            sides = run.sides(done)  # mode-follow sides right after their parent
            for side in sides:
                run.settle(side)
            side_jobs += sides
            run.join_nearer_side(done)
            found = run.discovery(done)
            if found is not None:
                run.extra.append(_artifact(found.discovery_id, found))
        (rt.stage_dir / "diagnostics.json").write_text(json.dumps(run.history, indent=1))
        return [*(_artifact(k, ev) for k, ev in run.calcs.items()),
                *(_artifact(f"species_{j.species.species_id}", j.species) for j in side_jobs),
                *(_artifact(m.minimum_id, m) for m in run.minima.values()), *run.extra]
