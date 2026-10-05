"""explore stage (design §4.1 #4, §8.2): low-level edges between states, over generations.

A state is (composition, state label); the start states are those of the screen minima and of
the input species. Generation 1 runs from the screen states, on their minima within
SCREEN_WINDOW_KCAL of the state's lowest; each later one from the reached states not yet run,
on their IRC-end minima, lowest GFN2 bottleneck from a screen state first (an order, not a
threshold). A state's edit classes (chemistry.trials) run once, each on one structure;
``max_trials`` attempts run over all generations in that order, then the classes' own; the
classes beyond it, and the states left once it is spent, are ``not_attempted``.

An edge (``product``) joins two states by an NT2 TS and its two IRC-end minima, or without a TS
by a structure and the minimum it relaxed into: an attempt's start, or the lowest seed of a state
the screen lost that a resolved bond change separates from its basin (``relaxation``). Its ends
are species in the attempt's atom order: the attempt's structure for an end with its bonds, else
a known minimum holding the end as labelled (same bonds, one basin), else a new species; a
structure that relaxed is its own species. A TS in one basin with an earlier edge's adds no edge
(``same_edge``). An edge runs from its end nearer the start states; one with no reached end is
``unconnected``. ``generation``: 1 + the least number of edges from a start state to the source
end.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, ClassVar, cast

from hfauto.backends.protocols import (
    Capability,
    DiscoveryEngine,
    DiscoveryResult,
    DiscoverySettings,
)
from hfauto.chemistry import identity, topology, trials
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, Geometry
from hfauto.core.hashing import sha256_text
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

State = tuple[str, str]  # (composition id, state label)
SCREEN_WINDOW_KCAL = 6.0  # the conformer search's energy window
_FAR = (float("inf"), float("inf"))


class ExploreConfig(StageConfig):
    engine: str
    method: str
    max_trials: int = 3000  # NT2 attempts over all generations
    settings: DiscoverySettings = DiscoverySettings()


@dataclass(frozen=True)
class _Unit:
    state: State
    source: SpeciesRecord  # the structure the class is realised on
    start: Molecule
    trial: ReactionTrial
    order: tuple[Any, ...]  # trials.Trial.order


def _reversed(d: DiscoveryRecord) -> DiscoveryRecord:
    act, rxn = d.dE_act_kcal, d.dE_rxn_kcal
    return d.model_copy(update={
        "source_species": d.product_species, "product_species": d.source_species,
        "dE_act_kcal": None if act is None or rxn is None else act - rxn,
        "dE_rxn_kcal": None if rxn is None else -rxn})


@dataclass
class _Network:
    """The states, the minima that represent them and the records so far, in order."""

    load: Callable[[Geometry], XYZ]
    starts: set[State]
    state: dict[str, State] = field(default_factory=dict)  # species id -> its state
    minima: defaultdict[State, list[tuple[SpeciesRecord, XYZ]]] = field(
        default_factory=lambda: defaultdict(list))
    saddles: list[tuple[str, str, XYZ]] = field(default_factory=list)  # composition, edge, TS
    level: dict[State, tuple[float, float]] = field(default_factory=dict)  # bottleneck, energy
    records: list[tuple[DiscoveryRecord, Failure | None]] = field(default_factory=list)
    species: list[SpeciesRecord] = field(default_factory=list)  # the new ones

    def new(self, sid: str, like: SpeciesRecord, geometry: Geometry, *, minimum: bool) -> str:
        xyz = self.load(geometry)
        state = (like.composition_id, topology.state_label(xyz.symbols, xyz.coords))
        record = like.model_copy(update={
            "species_id": sid, "geometry": geometry, "source": "discovery",
            "state_label": state[1], "energy_hartree": None, "level_key": None})
        self.species.append(record)
        self.state[sid] = state
        if minimum:
            self.minima[state].append((record, xyz))
        return sid

    def end(self, like: SpeciesRecord, geometry: Geometry, sid: str) -> str:
        """An IRC-end minimum's species: a known minimum holding it as labelled, else new."""
        x = self.load(geometry)
        held = {s.species_id: (y.coords, 0.0) for (c, _), items in self.minima.items()
                if c == like.composition_id for s, y in items
                if y.symbols == x.symbols and topology.same_bonding(x.symbols, x.coords, y.coords)}
        return (identity.assign(x.symbols, x.coords, 0.0, held)
                or self.new(sid, like, geometry, minimum=True))

    def same_edge(self, composition: str, ts: Geometry) -> str | None:
        x = self.load(ts)
        held = {k: (y.coords, 0.0) for c, k, y in self.saddles
                if c == composition and y.symbols == x.symbols}
        return identity.assign(x.symbols, x.coords, 0.0, held)

    def attempt(self, unit: _Unit, result: DiscoveryResult | Failure) -> None:
        did = f"disc_{unit.trial.trial_id}"
        record = DiscoveryRecord(discovery_id=did, mechanism="nt2", trial=unit.trial,
                                 outcome="failed", source_species=unit.source.species_id)
        if isinstance(result, Failure):
            reason = f"{result.kind}:{result.reason}"
            self.records.append((record.model_copy(update={"reason": reason}), result))
            return
        fields: dict[str, Any] = {"outcome": "negative", "reason": result.reason, "ts": result.ts,
                                  "dE_act_kcal": result.dE_act_kcal,
                                  "dE_rxn_kcal": result.dE_rxn_kcal}
        composition, ts = unit.state[0], result.ts
        same = None if result.ends is None or ts is None else self.same_edge(composition, ts)
        if same is not None:
            fields["reason"] = f"same_edge:{same}"
        elif result.ends is not None:
            first, second = result.ends
            if ts is None:  # the start relaxed: a species of its own
                source = self.new(f"spc_{did}_0", unit.source, first, minimum=False)
            elif result.irc_connected_to_source:
                source = unit.source.species_id
            else:
                source = self.end(unit.source, first, f"spc_{did}_0")
            fields |= {"outcome": "product", "source_species": source,
                       "product_species": self.end(unit.source, second, f"spc_{did}_1")}
            if ts is not None:
                self.saddles.append((composition, did, self.load(ts)))
        record = record.model_copy(update=fields)
        self.records.append((record, None))
        if record.outcome == "product":
            self.reach(record)

    def reach(self, d: DiscoveryRecord) -> None:
        """The product's bottleneck from a screen state (kcal/mol): the least, over its edges,
        of the highest point on the way (the TS, else the higher end)."""
        source = self.state[d.source_species or ""]
        if source in self.level:
            top, energy = self.level[source]
            after = energy + (d.dE_rxn_kcal or 0.0)
            point = (max(top, after, energy + (d.dE_act_kcal or 0.0)), after)
            product = self.state[d.product_species or ""]
            self.level[product] = min(self.level.get(product, _FAR), point)

    def units(self, state: State) -> list[_Unit]:
        """One unit per class of the state. Its structures are enumerated per atom order (two
        declared compositions of the same atoms, S6, reach one state in two orders); a class
        found in both runs once, on the realisation of least order."""
        orders: defaultdict[tuple[str, ...], list[tuple[SpeciesRecord, XYZ]]] = defaultdict(list)
        for source, xyz in self.minima[state]:
            orders[tuple(xyz.symbols)].append((source, xyz))
        best: dict[str, _Unit] = {}
        for symbols, members in orders.items():
            first = members[0][0]
            for t in trials.trials(symbols, [x.coords for _, x in members], first.charge,
                                   first.multiplicity):
                trial = ReactionTrial(trial_id="trial_" + sha256_text(f"{state[0]}|{t.key}"),
                                      kind=f"f{len(t.formed)}b{len(t.broken)}",
                                      associations=t.formed, dissociations=t.broken)
                source = members[t.conformer][0]
                start = Molecule(XYZ(list(symbols), t.start), source.charge, source.multiplicity)
                if t.key not in best or t.order < best[t.key].order:
                    best[t.key] = _Unit(state, source, start, trial, t.order)
        return list(best.values())

    def depths(self) -> dict[State, int]:
        """The least number of edges from a start state, per reached state."""
        depth = dict.fromkeys(self.starts, 0)
        edges = [(self.state[d.source_species or ""], self.state[d.product_species or ""])
                 for d, _ in self.records if d.outcome == "product"]
        changed = True
        while changed:
            changed = False
            for a, b in edges + [(b, a) for a, b in edges]:
                if a in depth and (b not in depth or depth[a] + 1 < depth[b]):
                    depth[b], changed = depth[a] + 1, True
        return depth

    def finish(self, d: DiscoveryRecord, depth: Mapping[State, int]) -> DiscoveryRecord:
        """The generation; an edge runs from its end nearer the start states."""
        a, b = (depth.get(self.state.get(s or "", ("", "")))
                for s in (d.source_species, d.product_species))
        if d.outcome == "product":
            if a is None and b is None:
                return d.model_copy(update={"outcome": "unconnected"})
            if a is None or (b is not None and b < a):
                d, a = _reversed(d), b
        return d if a is None else d.model_copy(update={"generation": a + 1})

    def expand(self, attempt: Callable[[list[_Unit]], list[DiscoveryResult | Failure]],
               budget: int) -> None:
        """Generations while a reached state has not been run from; ``budget`` attempts in all.
        A generation enumerates its states in bottleneck order, one bottleneck at a time, only
        until the budget is covered: the classes beyond it are ``not_attempted``, and so is each
        state left (no trial: its classes are never enumerated; S19 and S10 reach hundreds)."""
        done: set[State] = set()
        frontier = sorted(self.minima)
        while frontier and budget:
            units: list[_Unit] = []
            for _, states in itertools.groupby(
                    sorted(frontier, key=lambda s: (self.level.get(s, _FAR)[0], s)),
                    key=lambda s: self.level.get(s, _FAR)[0]):
                if len(units) >= budget:
                    break
                group = list(states)
                done.update(group)
                units += sorted((u for s in group for u in self.units(s)), key=lambda u: u.order)
            run, cut, budget = units[:budget], units[budget:], max(budget - len(units), 0)
            for unit, result in zip(run, attempt(run), strict=True):
                self.attempt(unit, result)
            self.records += [(DiscoveryRecord(
                discovery_id=f"disc_{u.trial.trial_id}", mechanism="nt2", trial=u.trial,
                outcome="not_attempted", source_species=u.source.species_id), None) for u in cut]
            reached = self.depths()
            frontier = sorted(s for s in self.minima if s in reached and s not in done)
        self.records += [(DiscoveryRecord(
            discovery_id="disc_state_" + sha256_text("|".join(s)), mechanism="nt2",
            outcome="not_attempted", source_species=self.minima[s][0][0].species_id), None)
            for s in frontier]

    def artifacts(self) -> list[Artifact]:
        depth = self.depths()
        records = [(self.finish(d, depth), failure) for d, failure in self.records]
        return [*(Artifact(artifact_id=d.discovery_id, type=ArtifactType.DISCOVERY, payload=d,
                           status="success" if f is None else "failed", failure=f)
                  for d, f in records),
                *(Artifact(artifact_id=species_artifact_id(s.species_id),
                           type=ArtifactType.SPECIES, payload=s) for s in self.species)]


def _relaxations(net: _Network, species: Mapping[str, SpeciesRecord],
                 screen: Sequence[MinimumRecord], final: Mapping[str, XYZ]) -> None:
    """One edge without a TS per seed state the screen lost (no screen minimum of the
    composition has its label): from the seed of lowest low-level energy (then species id) that
    a resolved bond change separates from the basin it collapsed into (``final``, carried into
    the seed's labelling by identity.basin_coords), unrelaxed as a species of its own, to that
    basin. A threshold crossing inside one basin loses no state (TMA·(HF)2: N···H 1.426 Å in a
    seed, 1.414 Å in its basin, r_thr 1.42 Å)."""
    kept = {(m.composition_id, m.state_label) for m in screen}
    seeds = sorted(((species[sid], m) for m in screen
                    for sid in dict.fromkeys((m.species_id, *m.members)) if sid in species),
                   key=lambda p: (p[0].energy_hartree is None, p[0].energy_hartree or 0.0,
                                  p[0].species_id))
    lost: dict[State, tuple[SpeciesRecord, MinimumRecord]] = {}
    for seed, m in seeds:
        state, x = (m.composition_id, seed.state_label), net.load(seed.geometry)
        basin = identity.basin_coords(x.symbols, final[m.minimum_id].coords, x.coords)
        if state not in kept | lost.keys() and not topology.same_bonding(
                x.symbols, x.coords, basin):
            lost[state] = (seed, m)
    for seed, m in lost.values():
        did = f"relax_{seed.species_id}"
        source = net.new(f"spc_{did}_0", seed, seed.geometry, minimum=False)
        net.records.append((DiscoveryRecord(
            discovery_id=did, mechanism="relaxation", outcome="product", source_species=source,
            product_species=m.species_id), None))


def _network(species: Mapping[str, SpeciesRecord], screen: Sequence[MinimumRecord],
             final: Mapping[str, XYZ], load: Callable[[Geometry], XYZ]) -> _Network:
    """The start states, each screen state's minima within SCREEN_WINDOW_KCAL of its lowest
    (lowest first) at bottleneck 0, and the relaxations."""
    net = _Network(load, {(m.composition_id, m.state_label) for m in screen}
                   | {(s.composition_id, s.state_label) for s in species.values()})
    lowest: dict[State, float] = {}
    for m in sorted(screen, key=lambda m: (m.energy_hartree, m.minimum_id)):
        state = net.state[m.species_id] = (m.composition_id, m.state_label)
        lowest.setdefault(state, m.energy_hartree)
        if HARTREE_TO_KCAL_MOL * (m.energy_hartree - lowest[state]) <= SCREEN_WINDOW_KCAL:
            net.minima[state].append((species[m.species_id], final[m.minimum_id]))
        net.level[state] = (0.0, 0.0)
    _relaxations(net, species, screen, final)
    return net


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
        net = _network(species, screen, final, rt.load_xyz)
        engine = cast(DiscoveryEngine, rt.engine(Capability.DISCOVERY, cfg.engine))
        method = rt.method(cfg.method)

        def attempt(unit: _Unit) -> DiscoveryResult | Failure:
            return engine.explore(unit.start, unit.trial, method, cfg.settings)

        net.expand(lambda units: rt.thread_map(attempt, units), cfg.max_trials)
        return net.artifacts()
