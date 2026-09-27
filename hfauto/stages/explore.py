"""explore stage (design §4.1 #4, §8.2): reaction trials on screen minima with a discovery engine.

Sources are the ``sources_per_state`` lowest screen minima of each composition × state label,
each started from its screen-optimized structure (the final geometry of ``opt_calc``). Every
trial runs NT2 first and falls back to AFIR on the same drive when NT2 gives no product (an
out-of-window product counts as found). Seeds that relaxed into another state become
``relaxation`` discoveries. Every attempt is recorded: product, negative or failed.
"""

from __future__ import annotations

from collections import defaultdict
from typing import ClassVar, cast

from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import (
    Capability,
    DiscoveryEngine,
    DiscoveryResult,
    DiscoverySettings,
)
from hfauto.chemistry import topology, trials
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key
from hfauto.core.evidence import Failure, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionTrial,
    SpeciesRecord,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class Window(BaseModel):
    model_config = ConfigDict(extra="forbid")
    barrier_kj: float = 150.0  # ΔE‡ at 300 K (low level)
    reaction_kj: float = 100.0


class ExploreConfig(StageConfig):
    engine: str
    method: str
    sources_per_state: int = 2
    max_trials_per_source: int = 10
    settings: DiscoverySettings = DiscoverySettings()
    window: Window = Window()


def _sources(minima: list[MinimumRecord], per_state: int) -> list[MinimumRecord]:
    groups: defaultdict[tuple[str, str], list[MinimumRecord]] = defaultdict(list)
    for m in sorted(minima, key=lambda m: (m.energy_hartree, m.minimum_id)):
        groups[(m.composition_id, m.state_label)].append(m)
    return [m for group in groups.values() for m in group[:per_state]]


class _Explorer:
    def __init__(self, config: ExploreConfig, rt: StageRuntime) -> None:
        self.config, self.rt = config, rt
        self.engine = cast(DiscoveryEngine, rt.engine(Capability.DISCOVERY, config.engine))
        self.method: MethodSpec = rt.method(config.method)

    def source(self, minimum: MinimumRecord, species: SpeciesRecord,
               optimized: Geometry) -> list[Artifact]:
        xyz = self.rt.load_xyz(optimized)
        start, drives = trials.generate(minimum.minimum_id, xyz.symbols, xyz.coords,
                                        charge=species.charge, multiplicity=species.multiplicity,
                                        max_trials=self.config.max_trials_per_source)
        mol = Molecule(XYZ(list(xyz.symbols), start), species.charge, species.multiplicity)
        out: list[Artifact] = []
        for trial in drives:
            for mechanism in ("nt2", "afir"):
                attempt = trial.model_copy(update={"mechanism": mechanism})
                result = self.engine.explore(mol, attempt, self.method, self.config.settings)
                artifacts, found = self.record(minimum, attempt, result, species)
                out += artifacts
                if found:
                    break
        return out

    def record(self, minimum: MinimumRecord, trial: ReactionTrial,
               result: DiscoveryResult | Failure,
               source: SpeciesRecord) -> tuple[list[Artifact], bool]:
        """(artifacts, whether a product was found, kept or out of window)."""
        discovery_id = f"disc_{trial.trial_id}_{trial.mechanism}"
        base = DiscoveryRecord(discovery_id=discovery_id, source_minimum=minimum.minimum_id,
                               mechanism=trial.mechanism, trial=trial, outcome="failed")
        parents = (minimum.minimum_id,)
        if isinstance(result, Failure):
            payload = base.model_copy(update={"reason": f"{result.kind}:{result.reason}"})
            return [Artifact(artifact_id=discovery_id, type=ArtifactType.DISCOVERY,
                             parents=parents, status="failed", payload=payload,
                             failure=result)], False
        reason, window = result.reason, self.config.window
        validated = result.ts is not None and result.irc_connected_to_source
        if result.outcome == "product":
            reason = "no_product_structure" if result.product is None else trials.product_verdict(
                trial.mechanism, ts_validated=validated, barrier_kj=result.barrier_kj_mol,
                reaction_kj=result.reaction_kj_mol, barrier_max_kj=window.barrier_kj,
                reaction_max_kj=window.reaction_kj)
        kept = result.outcome == "product" and reason is None
        product = self.product(discovery_id, result.product, source) if kept else None
        record = base.model_copy(update={
            "outcome": "product" if kept else "negative", "reason": reason,
            "product_species": product.species_id if product else None, "ts": result.ts,
            "ts_imag_cm1": result.ts_imag_cm1, "barrier_kj_mol": result.barrier_kj_mol,
            "reaction_kj_mol": result.reaction_kj_mol,
            "electronic_temperature_K": result.electronic_temperature_K})
        out = [Artifact(artifact_id=discovery_id, type=ArtifactType.DISCOVERY, parents=parents,
                        payload=record)]
        if product is not None:
            out.append(Artifact(artifact_id=f"species_{product.species_id}",
                                type=ArtifactType.SPECIES, parents=(discovery_id,),
                                payload=product))
        return out, kept or reason == "out_of_window"

    def product(self, discovery_id: str, geometry: Geometry | None,
                source: SpeciesRecord) -> SpeciesRecord | None:
        if geometry is None:
            return None
        xyz = self.rt.load_xyz(geometry)
        return SpeciesRecord(
            species_id=f"spc_{discovery_id}",
            composition_id=composition_key(xyz.symbols, source.charge, source.multiplicity),
            charge=source.charge, multiplicity=source.multiplicity, geometry=geometry,
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
        out = [Artifact(artifact_id=d.discovery_id, type=ArtifactType.DISCOVERY,
                        parents=(d.source_minimum,), payload=d)
               for d in trials.relaxation_discoveries(species.values(), screen)]
        explorer = _Explorer(cfg, rt)
        for minimum in _sources(screen, cfg.sources_per_state):
            out += explorer.source(minimum, species[minimum.species_id],
                                   inputs.evidence(minimum.opt_calc).final)
        return out
