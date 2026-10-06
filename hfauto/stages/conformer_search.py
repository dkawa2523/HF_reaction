"""conformers stage (design §4.1 #2, §8.2): CREST conformer search of monomers, placement seeds
and ``--nci`` search of compositions from seed00, state labels, and per state label the
``keep_per_state`` lowest structures (no energy window and no own duplicate test: CREGEN's output
is unique, and the screen stage selects).

Small rigid monomers skip the search and pass their input through. A topology-change stop is kept
as a ``crest_topology`` species and the search reruns once from that structure with the same
settings (a second stop fails). Compositions take the summed charge of their components and the
declared multiplicity (or the only one spin coupling allows; otherwise INPUT_INVALID; a low-spin
coupling such as CH3·O2 doublet is accepted, its singlet is not: ``composition_multiplicity``),
and get ``--notopo`` on every atom: hfauto's state label decides the state, CREST only samples.

When no candidate carries the state label of the input (monomer) or of the seeds
(composition), a failed search included, those structures are output as well and the reason is
the item's Failure. Observed in VAL9 R2: CREST 3.0.2 failed on CH3·O2 (rc -11) and OH·CH4 (rc 1),
whose seeds are therefore the output; the other five composition systems returned rc 0. Users
set quick, ewin_kcal, seeds_per_composition and keep_per_state; CREST threads come from the site.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
from typing import ClassVar, Literal, cast

from pydantic import PositiveInt

from hfauto.backends import protocols as bp
from hfauto.chemistry import placement
from hfauto.chemistry.electronic_state import composition_multiplicity
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key, write_xyz, written_geometry
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, SpeciesRecord
from hfauto.core.system import CompositionInput
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

Source = Literal["conformer", "placement", "crest_topology"]
_TAG = {"conformer": "c", "placement": "p", "crest_topology": "t"}


class ConformersConfig(StageConfig):
    engine: str
    method: str
    quick: bool = True
    ewin_kcal: float = 6.0
    seeds_per_composition: PositiveInt = 6  # seed00 starts CREST
    keep_per_state: int = 6


@dataclass(frozen=True)
class _Cand:
    geometry: Geometry
    energy: float | None
    source: Source
    xyz: XYZ
    label: str


def _cand(geometry: Geometry, energy: float | None, source: Source, xyz: XYZ) -> _Cand:
    return _Cand(geometry, energy, source, xyz, state_label(xyz.symbols, xyz.coords))


@dataclass
class _Item:
    base: str  # species or composition id
    parents: tuple[str, ...]
    mol: Molecule
    settings: bp.ConformerSettings
    fallback: tuple[_Cand, ...]  # the input or the seeds: the declared state
    skip: bool = False


@dataclass
class _Found:
    item: _Item
    candidates: list[_Cand]
    failure: Failure | None = None


def _failed(base: str, parents: tuple[str, ...], failure: Failure) -> Artifact:
    return Artifact(artifact_id=f"conformers_{base}", type=ArtifactType.SPECIES, parents=parents,
                    status="failed", failure=failure)


def _crest(engine: bp.ConformerEngine, method: MethodSpec, load_xyz: Callable[[Geometry], XYZ],
           item: _Item) -> tuple[list[_Cand], Failure | None]:
    """The topology stops, then the members of the search (rerun once from the first stop)."""
    result = engine.search(item.mol, method, item.settings)
    cands: list[_Cand] = []
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:
        cands += [_cand(g, None, "crest_topology", load_xyz(g)) for g in result.topology_stops]
        stop = Molecule(cands[0].xyz, item.mol.charge, item.mol.multiplicity)
        result = engine.search(stop, method, item.settings)
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:  # no second retry
        result = Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason="topology_stop_repeated")
    if isinstance(result, Failure):
        return cands, result
    return cands + [_cand(g, e, "conformer", load_xyz(g)) for g, e in result.members], None


def _search(engine: bp.ConformerEngine, method: MethodSpec, load_xyz: Callable[[Geometry], XYZ],
            item: _Item) -> _Found:
    if item.skip:
        return _Found(item, list(item.fallback))
    cands, failure = _crest(engine, method, load_xyz, item)
    declared = {c.label for c in item.fallback}
    if all(c.label not in declared for c in cands):
        cands += item.fallback
        failure = failure or Failure(kind=FailureKind.GATE_REJECTED,
                                     reason=f"input_state_lost:{','.join(sorted(declared))}")
    return _Found(item, cands, failure)


def _select(cands: list[_Cand], keep_per_state: int) -> list[_Cand]:
    """Per state label the lowest ``keep_per_state``. A missing energy is dropped when the state
    has energies; a state without energies keeps input order."""
    ranked = sorted(cands, key=lambda c: (c.energy is None, c.energy or 0.0))
    kept: list[_Cand] = []
    for label in dict.fromkeys(c.label for c in ranked):
        group = [c for c in ranked if c.label == label]
        lowest = group[0].energy
        kept += [c for c in group if lowest is None or c.energy is not None][:keep_per_state]
    return kept


def _artifacts(found: _Found, cfg: ConformersConfig) -> list[Artifact]:
    item, charge, mult = found.item, found.item.mol.charge, found.item.mol.multiplicity
    out = [] if found.failure is None else [_failed(item.base, item.parents, found.failure)]
    for k, c in enumerate(_select(found.candidates, cfg.keep_per_state)):
        record = SpeciesRecord(
            species_id=f"{item.base}_{_TAG[c.source]}{k:02d}",
            composition_id=composition_key(c.xyz.symbols, charge, mult), charge=charge,
            multiplicity=mult, geometry=c.geometry, source=c.source, state_label=c.label,
            energy_hartree=c.energy,
        )
        out.append(Artifact(artifact_id=species_artifact_id(record.species_id),
                            type=ArtifactType.SPECIES, parents=item.parents, payload=record))
    return out


def _monomer(sp: SpeciesRecord, rt: StageRuntime, cfg: ConformersConfig) -> _Item:
    xyz = rt.load_xyz(sp.geometry)
    settings = bp.ConformerSettings(quick=cfg.quick, ewin_kcal=cfg.ewin_kcal)
    return _Item(sp.species_id, (species_artifact_id(sp.species_id),),
                 Molecule(xyz, sp.charge, sp.multiplicity), settings,
                 (_cand(sp.geometry, None, "conformer", xyz),),
                 skip=placement.is_small_rigid(xyz.symbols, xyz.coords))


def _monomers(species: dict[str, SpeciesRecord], rt: StageRuntime, cfg: ConformersConfig
              ) -> list[_Item]:
    return [_monomer(species[s.id], rt, cfg) for s in rt.system.species
            if s.role == "monomer" and s.id in species]


def _lowest(found: _Found) -> XYZ:
    """The lowest candidate; without energies a relaxed topology stop beats the input."""
    return min(found.candidates,
               key=lambda c: (c.energy is None, c.source != "crest_topology", c.energy or 0)).xyz


def _composition(comp: CompositionInput, species: dict[str, SpeciesRecord],
                 best: dict[str, XYZ], rt: StageRuntime, cfg: ConformersConfig
                 ) -> _Item | Artifact:
    def reject(reason: str) -> Artifact:
        return _failed(comp.id, (), Failure(kind=FailureKind.INPUT_INVALID, reason=reason))

    missing = sorted(set(comp.components) - set(species))
    if missing:
        return reject(f"missing_component:{','.join(missing)}")
    parts = sorted(((sid, best.get(sid) or rt.load_xyz(species[sid].geometry))
                    for sid, n in comp.components.items() for _ in range(n)),
                   key=lambda p: (-len(p[1].symbols), p[0]))
    charge = sum(species[sid].charge for sid, _ in parts)
    try:
        mult = composition_multiplicity([species[sid].multiplicity for sid, _ in parts],
                                        comp.multiplicity)
    except ValueError as exc:
        return reject(str(exc))
    folder = rt.stage_dir / "placement"
    seeds = []
    for k, seed in enumerate(placement.seeds(parts[0][1], [xyz for _, xyz in parts[1:]],
                                             n_seeds=cfg.seeds_per_composition)):
        geometry = written_geometry(write_xyz(seed, folder / f"{comp.id}_seed{k:02d}.xyz"),
                                    rt.file_ref)  # CREST starts from the coordinates as written
        seeds.append(_cand(geometry, None, "placement", rt.load_xyz(geometry)))
    settings = bp.ConformerSettings(nci=True, quick=cfg.quick, ewin_kcal=cfg.ewin_kcal)
    parents = tuple(dict.fromkeys(species_artifact_id(sid) for sid, _ in parts))
    return _Item(comp.id, parents, Molecule(seeds[0].xyz, charge, mult), settings, tuple(seeds))


class ConformersStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="conformers", config=ConformersConfig, consumes=(ArtifactType.SPECIES,),
        produces=(ArtifactType.SPECIES,),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(ConformersConfig, config)
        engine = cast(bp.ConformerEngine, rt.engine(bp.Capability.CONFORMERS, cfg.engine))
        run = partial(_search, engine, rt.method(cfg.method), rt.load_xyz)

        def search(item: _Item) -> _Found:  # a raising search fails its item alone
            found = rt.contain(item.base, lambda: run(item))
            return _Found(item, list(item.fallback), found) if isinstance(found, Failure) else found
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)
                   if s.source == "input"}
        # the JobRunner core semaphore limits concurrent CREST jobs by the site threads
        found = rt.thread_map(search, _monomers(species, rt, cfg))
        best = {f.item.base: _lowest(f) for f in found}
        built = [_composition(c, species, best, rt, cfg) for c in rt.system.compositions]
        items = [b for b in built if isinstance(b, _Item)]
        found += rt.thread_map(search, items)
        failed = [b for b in built if isinstance(b, Artifact)]
        return failed + [a for f in found for a in _artifacts(f, cfg)]
