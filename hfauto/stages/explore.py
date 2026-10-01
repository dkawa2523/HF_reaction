"""explore stage (design §4.1 #4, §8.2): reaction trials on screen minima with a discovery engine.

Sources are the ``sources_per_state`` lowest screen minima of each composition × state label,
started from their screen-optimized structures. Each (source, trial) unit runs NT2 once (through
``thread_map``); every attempt is recorded in input order: product, negative or failed, with its
source species (the source's representative, in whose atom order the trial ran). A kept product
joins a known basin (a screen minimum or an earlier product of its composition) only as
labelled: the same atom-indexed bonds and permutation-invariant RMSD (ReaDuct's and the screen's
xTB energies are not compared); else it becomes a species. So a discovery's ends keep the
labelling it followed, and a degenerate product never joins its source. A lost
seed state (a bond change from the basin it collapsed into) becomes one ``relaxation`` product,
from the seed refined as its own species.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import ClassVar, cast

from hfauto.backends.protocols import (
    Capability,
    DiscoveryEngine,
    DiscoveryResult,
    DiscoverySettings,
)
from hfauto.chemistry import gates, identity, topology, trials
from hfauto.chemistry.hypotheses import seed_species_id
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.evidence import Failure, Geometry
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionTrial,
    SpeciesRecord,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class ExploreConfig(StageConfig):
    engine: str
    method: str
    sources_per_state: int = 2
    max_trials_per_source: int = 10
    settings: DiscoverySettings = DiscoverySettings()


@dataclass(frozen=True)
class _Unit:
    minimum: MinimumRecord
    start: Molecule
    trial: ReactionTrial


@dataclass(frozen=True)
class _Basin:
    """A known structure of one composition, as labelled: a screen minimum's or a kept product's."""

    species_id: str
    xyz: XYZ


def _basin_of(product: _Basin, known: Sequence[_Basin]) -> str | None:
    """The species id of the known basin that holds ``product`` as labelled (identity.assign
    without the energy criterion, among the known structures with its atom-indexed bonds), else
    None: a relabelled copy (a degenerate one of its source included) is a species of its own."""
    symbols, coords = product.xyz.symbols, product.xyz.coords
    bonded = topology.bonds(symbols, coords)
    candidates = {k.species_id: (k.xyz.coords, 0.0) for k in known
                  if list(k.xyz.symbols) == list(symbols)
                  and topology.bonds(k.xyz.symbols, k.xyz.coords) == bonded}
    return identity.assign(symbols, coords, 0.0, candidates)


def _sources(minima: list[MinimumRecord], per_state: int) -> list[MinimumRecord]:
    groups: defaultdict[tuple[str, str], list[MinimumRecord]] = defaultdict(list)
    for m in sorted(minima, key=lambda m: (m.energy_hartree, m.minimum_id)):
        groups[(m.composition_id, m.state_label)].append(m)
    return [m for group in groups.values() for m in group[:per_state]]


def _left_state(seed: XYZ, basin: XYZ) -> bool:
    """True when a bond change (topology.bond_changes, resolved by its band) separates the seed
    from its basin's structure carried into the seed's labelling (identity.basin_coords)."""
    symbols = seed.symbols
    return any(topology.bond_changes(
        symbols, seed.coords, identity.basin_coords(symbols, basin.coords, seed.coords)))


def _relaxations(species: Mapping[str, SpeciesRecord], minima: Sequence[MinimumRecord],
                 final: Mapping[str, XYZ], load: Callable[[Geometry], XYZ]
                 ) -> list[DiscoveryRecord]:
    """One ``relaxation`` product per seed state the screen lost: no screen minimum of the
    composition keeps its label, and a bond change separates a seed of that label from
    the basin it collapsed into (``final``: screen-optimized structures by minimum id). A
    threshold crossing inside one basin loses no state (TMA·(HF)2: N···H 1.426 Å in a seed,
    1.414 Å in its basin, r_thr 1.42 Å); one resolved change suffices (OH + CH4: H2–O5 +2.00 /
    −0.35 Å about r_thr). The product is the lost seed of lowest low-level energy (then species id),
    unrelaxed, from its basin, for the DFT tier to ask once. It has no low-level stationary
    step, so discovery_verdict does not apply."""
    kept = {(m.composition_id, m.state_label) for m in minima}
    seeds = [(species[sid], m) for m in minima
             for sid in dict.fromkeys((m.species_id, *m.members)) if sid in species]
    lost: dict[tuple[str, str], tuple[SpeciesRecord, MinimumRecord]] = {}
    for seed, m in sorted(seeds, key=lambda p: (p[0].energy_hartree is None,
                                                p[0].energy_hartree or 0.0, p[0].species_id)):
        state = (m.composition_id, seed.state_label)
        if state not in kept and state not in lost and _left_state(
                load(seed.geometry), final[m.minimum_id]):
            lost[state] = (seed, m)
    found = [DiscoveryRecord(discovery_id=f"relax_{seed.species_id}", source_minimum=m.minimum_id,
                             mechanism="relaxation", outcome="product",
                             product_species=seed.species_id)
             for seed, m in lost.values()]
    return [d.model_copy(update={"source_species": seed_species_id(d)}) for d in found]


def _units(minimum: MinimumRecord, species: SpeciesRecord, xyz: XYZ,
           max_trials: int) -> list[_Unit]:
    start, drives = trials.generate(minimum.minimum_id, xyz.symbols, xyz.coords,
                                    charge=species.charge, multiplicity=species.multiplicity,
                                    max_trials=max_trials)
    mol = Molecule(XYZ(list(xyz.symbols), start), species.charge, species.multiplicity)
    return [_Unit(minimum, mol, trial) for trial in drives]


class _Recorder:
    """Discovery artifacts; a kept product is identified against ``known`` (in call order)."""

    def __init__(self, rt: StageRuntime, known: dict[str, list[_Basin]]) -> None:
        self.rt, self.known = rt, known

    def record(self, unit: _Unit, result: DiscoveryResult | Failure) -> list[Artifact]:
        discovery_id = f"disc_{unit.trial.trial_id}_nt2"
        base = DiscoveryRecord(discovery_id=discovery_id, source_minimum=unit.minimum.minimum_id,
                               source_species=unit.minimum.species_id, mechanism="nt2",
                               trial=unit.trial, outcome="failed")
        parents = (unit.minimum.minimum_id,)
        if isinstance(result, Failure):
            payload = base.model_copy(update={"reason": f"{result.kind}:{result.reason}"})
            return [Artifact(artifact_id=discovery_id, type=ArtifactType.DISCOVERY,
                             parents=parents, status="failed", payload=payload,
                             failure=result)]
        reason, product, new = result.reason, None, None
        if result.outcome == "product" and result.product is None:
            reason = "no_product_structure"
        elif result.outcome == "product" and result.product is not None:
            reason = gates.discovery_verdict(
                ts_validated=result.ts is not None and result.irc_connected_to_source,
                dE_act_kcal=result.dE_act_kcal, dE_rxn_kcal=result.dE_rxn_kcal,
                policy=self.rt.policy)
            if reason is None:
                product, new = self.identify(discovery_id, result.product, unit)
        record = base.model_copy(update={
            "outcome": "negative" if product is None else "product", "reason": reason,
            "product_species": product, "ts": result.ts, "ts_imag_cm1": result.ts_imag_cm1,
            "dE_act_kcal": result.dE_act_kcal, "dE_rxn_kcal": result.dE_rxn_kcal,
            "electronic_temperature_K": result.electronic_temperature_K})
        out = [Artifact(artifact_id=discovery_id, type=ArtifactType.DISCOVERY, parents=parents,
                        payload=record)]
        if new is not None:
            out.append(Artifact(artifact_id=species_artifact_id(new.species_id),
                                type=ArtifactType.SPECIES, parents=(discovery_id,), payload=new))
        return out

    def identify(self, discovery_id: str, geometry: Geometry,
                 unit: _Unit) -> tuple[str, SpeciesRecord | None]:
        """(product species id, the new species or None when a known basin holds it)."""
        xyz, source, mol = self.rt.load_xyz(geometry), unit.minimum, unit.start
        product = _Basin(f"spc_{discovery_id}", xyz)
        known = self.known[source.composition_id]
        found = _basin_of(product, known)
        if found is not None:
            return found, None
        known.append(product)
        return product.species_id, SpeciesRecord(
            species_id=product.species_id, composition_id=source.composition_id,
            charge=mol.charge, multiplicity=mol.multiplicity, geometry=geometry,
            source="discovery", state_label=topology.state_label(xyz.symbols, xyz.coords))


class ExploreStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="explore",
        config=ExploreConfig,
        consumes=(ArtifactType.MINIMUM, ArtifactType.SPECIES),
        produces=(ArtifactType.DISCOVERY, ArtifactType.SPECIES),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(ExploreConfig, config)
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)}
        screen = [m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)
                  if m.tier == "screen"]
        final = {m.minimum_id: rt.load_xyz(inputs.evidence(m.opt_calc).final) for m in screen}
        known: dict[str, list[_Basin]] = defaultdict(list)
        for m in screen:
            known[m.composition_id].append(_Basin(m.species_id, final[m.minimum_id]))
        units = [unit for m in _sources(screen, cfg.sources_per_state)
                 for unit in _units(m, species[m.species_id], final[m.minimum_id],
                                    cfg.max_trials_per_source)]
        engine = cast(DiscoveryEngine, rt.engine(Capability.DISCOVERY, cfg.engine))
        method = rt.method(cfg.method)

        def attempt(unit: _Unit) -> DiscoveryResult | Failure:
            return engine.explore(unit.start, unit.trial, method, cfg.settings)

        out = [Artifact(artifact_id=d.discovery_id, type=ArtifactType.DISCOVERY,
                        parents=(d.source_minimum,), payload=d)
               for d in _relaxations(species, screen, final, rt.load_xyz)]
        recorder = _Recorder(rt, known)
        for unit, result in zip(units, rt.thread_map(attempt, units), strict=True):
            out += recorder.record(unit, result)
        return out
