"""explore stage (design §4.1 #4, §8.2): reaction trials on screen minima with a discovery engine.

Sources are the ``sources_per_state`` lowest screen minima of each composition × state label,
started from their screen-optimized structures. Each (source, trial) unit runs NT2, and AFIR on
the same drive only when NT2 found no maximum; the units run through ``thread_map`` and every
attempt is recorded in input order: product, negative or failed. A kept product joins a known
basin (a screen minimum or an earlier product of its composition: state label and permutation-
invariant RMSD; ReaDuct's and the screen's xTB energies are not compared) or becomes a species;
a degenerate one (the source's label) also needs the same atom-indexed bonds, so it never joins
its source. Seeds that relaxed into another state become ``relaxation`` discoveries.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import ClassVar, cast

from hfauto.backends.protocols import (
    NO_NT2_MAXIMUM,
    Capability,
    DiscoveryEngine,
    DiscoveryResult,
    DiscoverySettings,
)
from hfauto.chemistry import gates, identity, topology, trials
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.evidence import Failure, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionTrial,
    SpeciesRecord,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

Attempt = tuple[ReactionTrial, DiscoveryResult | Failure]


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
    """A known structure of one composition: a screen minimum's or a kept product's."""

    species_id: str
    state_label: str
    xyz: XYZ


def _basin_of(product: _Basin, known: Sequence[_Basin], *, degenerate: bool) -> str | None:
    """The species id of the known basin that holds ``product`` (identity.assign without the
    energy criterion), else None. A degenerate product also needs the atom-indexed bonds."""
    bonded = topology.bonds(product.xyz.symbols, product.xyz.coords)
    candidates = {
        k.species_id: (k.xyz.coords, 0.0) for k in known
        if k.state_label == product.state_label and list(k.xyz.symbols) == list(product.xyz.symbols)
        and not (degenerate and topology.bonds(k.xyz.symbols, k.xyz.coords) != bonded)
    }
    return identity.assign(product.xyz.symbols, product.xyz.coords, 0.0, candidates)


def _sources(minima: list[MinimumRecord], per_state: int) -> list[MinimumRecord]:
    groups: defaultdict[tuple[str, str], list[MinimumRecord]] = defaultdict(list)
    for m in sorted(minima, key=lambda m: (m.energy_hartree, m.minimum_id)):
        groups[(m.composition_id, m.state_label)].append(m)
    return [m for group in groups.values() for m in group[:per_state]]


def _relaxations(species: Iterable[SpeciesRecord],
                 minima: Iterable[MinimumRecord]) -> list[DiscoveryRecord]:
    """A seed whose state label differs from its basin's collapsed without a barrier."""
    labels = {s.species_id: s.state_label for s in species}
    return [
        DiscoveryRecord(discovery_id=f"relax_{m.minimum_id}_{sid}", source_minimum=m.minimum_id,
                        mechanism="relaxation", outcome="negative",
                        reason=f"collapsed_to:{m.state_label}")
        for m in minima for sid in dict.fromkeys((m.species_id, *m.members))
        if labels.get(sid, m.state_label) != m.state_label
    ]


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

    def record(self, unit: _Unit, trial: ReactionTrial,
               result: DiscoveryResult | Failure) -> list[Artifact]:
        discovery_id = f"disc_{trial.trial_id}_{trial.mechanism}"
        base = DiscoveryRecord(discovery_id=discovery_id, source_minimum=unit.minimum.minimum_id,
                               mechanism=trial.mechanism, trial=trial, outcome="failed")
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
                trial.mechanism, ts_validated=result.ts is not None and result.irc_connected_to_source,
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
            out.append(Artifact(artifact_id=f"species_{new.species_id}",
                                type=ArtifactType.SPECIES, parents=(discovery_id,), payload=new))
        return out

    def identify(self, discovery_id: str, geometry: Geometry,
                 unit: _Unit) -> tuple[str, SpeciesRecord | None]:
        """(product species id, the new species or None when a known basin holds it)."""
        xyz, source, mol = self.rt.load_xyz(geometry), unit.minimum, unit.start
        product = _Basin(f"spc_{discovery_id}", topology.state_label(xyz.symbols, xyz.coords), xyz)
        known = self.known[source.composition_id]
        found = _basin_of(product, known, degenerate=product.state_label == source.state_label)
        if found is not None:
            return found, None
        known.append(product)
        return product.species_id, SpeciesRecord(
            species_id=product.species_id, composition_id=source.composition_id,
            charge=mol.charge, multiplicity=mol.multiplicity, geometry=geometry,
            source="discovery", state_label=product.state_label)


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
            known[m.composition_id].append(_Basin(m.species_id, m.state_label, final[m.minimum_id]))
        units = [unit for m in _sources(screen, cfg.sources_per_state)
                 for unit in _units(m, species[m.species_id], final[m.minimum_id],
                                    cfg.max_trials_per_source)]
        engine = cast(DiscoveryEngine, rt.engine(Capability.DISCOVERY, cfg.engine))
        method = rt.method(cfg.method)

        def attempts(unit: _Unit) -> list[Attempt]:
            result = engine.explore(unit.start, unit.trial, method, cfg.settings)
            done: list[Attempt] = [(unit.trial, result)]
            if isinstance(result, DiscoveryResult) and result.reason == NO_NT2_MAXIMUM:
                afir = unit.trial.model_copy(update={"mechanism": "afir"})
                done.append((afir, engine.explore(unit.start, afir, method, cfg.settings)))
            return done

        out = [Artifact(artifact_id=d.discovery_id, type=ArtifactType.DISCOVERY,
                        parents=(d.source_minimum,), payload=d)
               for d in _relaxations(species.values(), screen)]
        recorder = _Recorder(rt, known)
        for unit, done in zip(units, rt.thread_map(attempts, units), strict=True):
            for trial, result in done:
                out += recorder.record(unit, trial, result)
        return out
