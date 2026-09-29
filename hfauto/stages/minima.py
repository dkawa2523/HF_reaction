"""minima stage (design §4.1 #3, §7.2, §8.2): selection → relax_to_minimum → Registry.

``level: screen`` relaxes every input species; ``level: dft`` refines a selection
(``chemistry.selection``) of the reacting compositions, started from the screen minima's
optimized structures. Jobs run serially in species-id order and each is registered at once, so a
later job that falls into a registered basin (mirror image included) is ``known`` and skips its
freq job (one identity criterion: identity.assign). An exact permutation or mirror image of an
earlier start (identity.carry within IMAGE_A) runs no job: the PES is invariant under both, so it
joins that start's basin when that start relaxed straight into a minimum. A saddle whose ±
displacements reach two distinct minima gives two ``mode_follow`` species at the driver's side
minima, registered by identity with no new job, and a ``mode_follow`` discovery; the saddle's
own species joins the side basin nearer its input structure (note ``endpoint_was_saddle``).
A relaxation product's seed (R6) is asked at DFT once, last, from its unrelaxed geometry as its
own species (``hypotheses.seed_species_id``; the seed may carry its collapse basin's job); one that
lands in another state than its own is noted ``collapsed_at_dft_from_seed``.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from typing import ClassVar, Literal, cast

from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.hypotheses import seed_species_id
from hfauto.chemistry.identity import IMAGE_A, carry, permutation_invariant_rmsd
from hfauto.chemistry.selection import Candidate, crowded, rerank, select_for_refinement
from hfauto.chemistry.thermo import monomer_states
from hfauto.chemistry.topology import fragments, state_label
from hfauto.chemistry.xyz import XYZ, Molecule, hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FailureKind, Geometry
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, level_mismatches
from hfauto.core.records import ArtifactType, DiscoveryRecord, MinimumRecord, Payload, SpeciesRecord
from hfauto.core.system import SystemConfig
from hfauto.drivers import minimum as driver
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

T = ArtifactType
_TIE_A = 0.05  # side RMSDs this close to the input count as equally near: the lower minimum wins
_MOVED = ("follow", "soft")  # steps along a mode, whose sign an image's relaxation need not share


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
_Side = tuple[SpeciesRecord, MinimumRecord]


def _image_groups(jobs: Sequence[_Job], load_xyz: Callable[[Geometry], XYZ]
                  ) -> list[list[_Job]]:
    """Jobs of one composition (charge and multiplicity included) and atom order whose starts
    are exact permutation or mirror images (identity.carry within IMAGE_A), in input order."""
    groups: list[tuple[XYZ, list[_Job]]] = []
    for job in jobs:
        x = load_xyz(job.start)
        group = next((members for ref, members in groups
                      if members[0].species.composition_id == job.species.composition_id
                      and ref.symbols == x.symbols
                      and carry(x.symbols, ref.coords, x.coords, x.coords)[0] <= IMAGE_A), None)
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
    seeds: list[SpeciesRecord] = field(default_factory=list)  # relaxation seeds (R6)
    extra: list[Artifact] = field(default_factory=list)  # discoveries and failed minima
    records: dict[str, MinimumRecord | None] = field(default_factory=dict)  # per species id
    history: dict[str, list[str]] = field(default_factory=dict)  # diagnostics.json

    def molecule(self, species: SpeciesRecord, geometry: Geometry) -> Molecule:
        return Molecule(self.rt.load_xyz(geometry), species.charge, species.multiplicity)

    def keep(self, *evidence: Evidence | None) -> None:
        self.calcs.update((driver.calc_id(ev), ev) for ev in evidence if ev is not None)

    def process(self, job: _Job) -> _Relaxed:
        """relax and register; at a TS candidate, its two side minima, the join of the saddle
        species to the nearer side and a mode_follow discovery."""
        mol = self.molecule(job.species, job.start)
        done = job, driver.relax_to_minimum(
            mol, self.method, self.qm, known=self.registry, init_hessian=self.init_hessian(job, mol),
            max_mode_follow=self.cfg.mode_follow, gates=self.rt.policy, load_xyz=self.rt.load_xyz)
        self.register(done)
        sides = self.sides(done)
        if len(sides) == 2:
            self.join_nearer_side(job, sides)
            found = self.discovery(done, sides)
            if found is not None:
                self.extra.append(_artifact(found.discovery_id, found))
        return done

    def init_hessian(self, job: _Job, mol: Molecule) -> Evidence | None:
        """A low-level freq at the start of 2+ fragments; none at a relaxation seed, which is no
        low-level stationary point (a negative xTB eigenvalue at each seed of S5, S6 and S19)."""
        spec, single = self.cfg.init_hessian, len(fragments(mol.xyz.symbols, mol.xyz.coords)) < 2
        if spec is None or single or job.species in self.seeds:
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
        """Registry.add (a known or identical basin is joined), at once: a later job in this
        basin is known (no freq job)."""
        job, outcome = done
        sid, record = job.species.species_id, None
        self.keep(outcome.opt, outcome.freq)
        self.history[sid] = [outcome.status, *outcome.history]
        if outcome.status in ("minimum", "soft_minimum", "known"):
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
        saddle endpoint, a soft minimum, a failure or a step along a mode)."""
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

    def note_collapses(self) -> None:
        """R6: a seed whose DFT minimum lies in another state than its own collapsed again."""
        for seed in self.seeds:
            record = self.records.get(seed.species_id)
            if record is not None and record.state_label != seed.state_label:
                self.history[seed.species_id].append("collapsed_at_dft_from_seed")

    def discovery(self, done: _Relaxed, sides: list[_Side]) -> DiscoveryRecord | None:
        """source_minimum = side 1's minimum, product_species = side 2, ts = the saddle; at the
        DFT tier ts_calc is its opt: a verified DFT saddle, validated directly by the case."""
        (_, a), (product, b) = sides
        parent, saddle, freq = done[0].species.species_id, done[1].opt, done[1].freq
        if saddle is None or freq is None:
            return None
        return DiscoveryRecord(
            discovery_id=f"disc_mode_follow_{parent}", source_minimum=a.minimum_id,
            mechanism="mode_follow", outcome="product", product_species=product.species_id,
            ts=saddle.final, ts_calc=driver.calc_id(saddle) if self.cfg.level == "dft" else None,
            ts_imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
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


def _seeds(inputs: Manifest, species: dict[str, SpeciesRecord]) -> list[SpeciesRecord]:
    """R6: the seed of each relaxation product, unrelaxed, as its own species
    (hypotheses.seed_species_id): the seed itself may carry its collapse basin's job."""
    return [species[d.product_species].model_copy(update={"species_id": seed_species_id(d)})
            for d in inputs.records(T.DISCOVERY, DiscoveryRecord)
            if d.mechanism == "relaxation" and d.product_species in species]


def _pool(inputs: Manifest, species: dict[str, SpeciesRecord], seeds: list[SpeciesRecord]
          ) -> dict[str, tuple[Candidate, _Job]]:
    """Screen minima (from their optimized structure) plus discovery products and sources;
    products that skipped screen and the relaxation seeds start from their own geometry."""
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
    for s in [*(species[sid] for sid in sorted((always & species.keys()) - held)), *seeds]:
        pool[s.species_id] = (Candidate(s.species_id, s.composition_id, s.state_label,
                                        s.energy_hartree, always=True), _Job(s, s.geometry))
    return pool


def _reacting(candidates: Iterable[Candidate], species: dict[str, SpeciesRecord],
              system: SystemConfig) -> set[str]:
    """Compositions of the declared endpoints and of the discovery sources and products (the
    always-kept candidates), plus the monomers of those complexes: the association and
    separated references of the thermo stage (thermo.monomer_states)."""
    ends = {species[s.id].composition_id for s in system.species
            if s.role == "endpoint" and s.id in species}
    reacting = ends | {c.composition_id for c in candidates if c.always}
    monomers = monomer_states(species.values(), system.compositions)
    formulas = {(hill_formula(s.geometry.symbols), s.charge) for s in species.values()
                if s.composition_id in reacting}
    return reacting | {state[0] for f in formulas for state, _ in monomers.get(f, ())}


def _jobs(inputs: Manifest, cfg: MinimaConfig, run: _Run) -> list[_Job]:
    """Every species (screen, all), or the window selection of the reacting compositions (the
    others are noted not_reacting) and the relaxation seeds (run.seeds, always kept); single
    points rerank only the crowded groups."""
    species = {s.species_id: s for s in inputs.records(T.SPECIES, SpeciesRecord)}
    if cfg.level == "screen" or cfg.select.include == "all":
        return [_Job(s, s.geometry) for s in species.values()]
    run.seeds = _seeds(inputs, species)
    pool, sel = _pool(inputs, species, run.seeds), cfg.select
    candidates = [candidate for candidate, _ in pool.values()]
    reacting = _reacting(candidates, species, run.rt.system)
    for composition in sorted({c.composition_id for c in candidates} - reacting):
        run.history[composition] = ["not_reacting"]
    chosen = select_for_refinement(
        [c for c in candidates if c.composition_id in reacting],
        per_state=sel.rerank_top if sel.rerank_sp else sel.per_state, window_kcal=sel.window_kcal)
    if sel.rerank_sp:
        crowd = crowded(chosen, sel.per_state)
        chosen = rerank(chosen, run.single_points([pool[c.species_id][1] for c in crowd]),
                        sel.per_state, sel.window_kcal)
    return [pool[c.species_id][1] for c in chosen]


class MinimaStage:
    spec: ClassVar[StageSpec] = StageSpec(  # discovery is optional input (none before explore)
        "minima", MinimaConfig, (T.SPECIES,), (T.CALCULATION, T.MINIMUM, T.SPECIES, T.DISCOVERY))

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(MinimaConfig, config)
        method = rt.method(cfg.method)
        qm = cast(QMEngine, rt.engine(Capability.QM, cfg.engine))
        run = _Run(rt, cfg, qm, method, _known(inputs, cfg, rt, qm, method))
        jobs = sorted(_jobs(inputs, cfg, run),  # relaxation seeds last: basins keep their species
                      key=lambda j: (j.species in run.seeds, j.species.species_id))
        for first, *images in _image_groups(jobs, rt.load_xyz):
            done = run.process(first)
            for image in images:
                if not run.join_image(image, done):
                    run.process(image)
        run.note_collapses()
        (rt.stage_dir / "diagnostics.json").write_text(json.dumps(run.history, indent=1))
        return [*(_artifact(k, ev) for k, ev in run.calcs.items()),
                *(_artifact(species_artifact_id(s.species_id), s)
                  for s in (*run.seeds, *run.species)),
                *(_artifact(m.minimum_id, m) for m in run.minima.values()), *run.extra]
