"""conformers stage (design §4.1 #2, §8.2): conformer search of monomers, placement seeds and
``--nci`` search of compositions, state labels and per state label the ``keep_per_state``
lowest structures (no energy window and no own duplicate test: CREGEN's output is unique, and
the screen stage selects).

Small rigid monomers skip the search and pass their input through. A topology-change stop is kept
as a ``crest_topology`` species and the search reruns once from that structure with the same
settings (a second stop fails). Compositions take the summed charge of their components and the
declared multiplicity (or the only one spin coupling allows; otherwise INPUT_INVALID), and get
``--notopo`` on every atom: hfauto's state label decides the state, CREST only samples. When a
search fails its input (monomer) or placement seeds (composition) are output instead (CREST 3.0.2
fails on open-shell and some small ionic compositions). Users set quick, ewin_kcal,
seeds_per_composition and keep_per_state; CREST threads come from the site.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
from typing import ClassVar, Literal, cast

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
Candidate = tuple[Geometry, float | None, Source]
_TAG = {"conformer": "c", "placement": "p", "crest_topology": "t"}


class ConformersConfig(StageConfig):
    engine: str
    method: str
    quick: bool = True
    ewin_kcal: float = 6.0
    seeds_per_composition: int = 6
    keep_per_state: int = 6


@dataclass
class _Item:
    base: str  # species or composition id
    parents: tuple[str, ...]
    mol: Molecule
    settings: bp.ConformerSettings
    fallback: tuple[Candidate, ...]  # output when the search is skipped or fails
    skip: bool = False


@dataclass
class _Found:
    item: _Item
    candidates: list[Candidate]
    failure: Failure | None = None


@dataclass(frozen=True)
class _Cand:
    geometry: Geometry
    energy: float | None
    source: Source
    xyz: XYZ
    label: str


def _failed(base: str, parents: tuple[str, ...], failure: Failure) -> Artifact:
    return Artifact(artifact_id=f"conformers_{base}", type=ArtifactType.SPECIES, parents=parents,
                    status="failed", failure=failure)


def _search(engine: bp.ConformerEngine, method: MethodSpec, load_xyz: Callable[[Geometry], XYZ],
            item: _Item) -> _Found:
    if item.skip:
        return _Found(item, list(item.fallback))
    result, stops = engine.search(item.mol, method, item.settings), []
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:
        stops = [(g, None, cast(Source, "crest_topology")) for g in result.topology_stops]
        stop = Molecule(load_xyz(result.topology_stops[0]), item.mol.charge, item.mol.multiplicity)
        result = engine.search(stop, method, item.settings)
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:  # no second retry
        result = Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason="topology_stop_repeated")
    if isinstance(result, Failure):
        return _Found(item, stops + list(item.fallback), result)
    return _Found(item, stops + [(g, e, cast(Source, "conformer")) for g, e in result.members])


def _select(cands: list[_Cand], keep_per_state: int) -> list[_Cand]:
    """Per state label the lowest ``keep_per_state``. A missing energy is dropped when the state
    has energies (BUG-08); a state without energies keeps input order."""
    ranked = sorted(cands, key=lambda c: (c.energy is None, c.energy or 0.0))
    kept: list[_Cand] = []
    for label in dict.fromkeys(c.label for c in ranked):
        group = [c for c in ranked if c.label == label]
        lowest = group[0].energy
        kept += [c for c in group if lowest is None or c.energy is not None][:keep_per_state]
    return kept


def _artifacts(found: _Found, rt: StageRuntime, cfg: ConformersConfig) -> list[Artifact]:
    item, charge, mult = found.item, found.item.mol.charge, found.item.mol.multiplicity
    out = [] if found.failure is None else [_failed(item.base, item.parents, found.failure)]
    cands = []
    for geometry, energy, source in found.candidates:
        xyz = rt.load_xyz(geometry)
        cands.append(_Cand(geometry, energy, source, xyz, state_label(xyz.symbols, xyz.coords)))
    for k, c in enumerate(_select(cands, cfg.keep_per_state)):
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
                 ((sp.geometry, None, "conformer"),),
                 skip=placement.is_small_rigid(xyz.symbols, xyz.coords))


def _monomers(species: dict[str, SpeciesRecord], rt: StageRuntime, cfg: ConformersConfig
              ) -> list[_Item]:
    return [_monomer(species[s.id], rt, cfg) for s in rt.system.species
            if s.role == "monomer" and s.id in species]


def _lowest(found: _Found, rt: StageRuntime) -> XYZ:
    ranked = sorted(found.candidates,  # without energies a relaxed topology stop beats the input
                    key=lambda c: (c[1] is None, c[2] != "crest_topology", c[1] or 0))
    return rt.load_xyz((ranked or list(found.item.fallback))[0][0])


def _seed_geometries(comp_id: str, seeds: list[XYZ], rt: StageRuntime) -> list[Geometry]:
    folder = rt.stage_dir / "placement"
    return [written_geometry(write_xyz(seed, folder / f"{comp_id}_seed{k:02d}.xyz"), rt.file_ref)
            for k, seed in enumerate(seeds)]


def _composition(comp: CompositionInput, species: dict[str, SpeciesRecord],
                 best: dict[str, XYZ], rt: StageRuntime, cfg: ConformersConfig
                 ) -> _Item | Artifact:
    def reject(kind: FailureKind, reason: str) -> Artifact:
        return _failed(comp.id, (), Failure(kind=kind, reason=reason))

    missing = sorted(set(comp.components) - set(species))
    if missing:
        return reject(FailureKind.INPUT_INVALID, f"missing_component:{','.join(missing)}")
    parts = sorted(((sid, best.get(sid) or rt.load_xyz(species[sid].geometry))
                    for sid, n in comp.components.items() for _ in range(n)),
                   key=lambda p: (-len(p[1].symbols), p[0]))
    charge = sum(species[sid].charge for sid, _ in parts)
    try:
        mult = composition_multiplicity([species[sid].multiplicity for sid, _ in parts],
                                        comp.multiplicity)
    except ValueError as exc:
        return reject(FailureKind.INPUT_INVALID, str(exc))
    seeds = placement.seeds(parts[0][1], [xyz for _, xyz in parts[1:]],
                            n_seeds=cfg.seeds_per_composition)
    if not seeds:
        return reject(FailureKind.GATE_REJECTED, "no_collision_free_seed")
    geoms = _seed_geometries(comp.id, seeds, rt)
    first = rt.load_xyz(geoms[0])
    settings = bp.ConformerSettings(nci=True, quick=cfg.quick, ewin_kcal=cfg.ewin_kcal)
    parents = tuple(dict.fromkeys(species_artifact_id(sid) for sid, _ in parts))
    return _Item(comp.id, parents, Molecule(first, charge, mult), settings,
                 tuple((g, None, "placement") for g in geoms))


class ConformersStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="conformers", config=ConformersConfig, consumes=(ArtifactType.SPECIES,),
        produces=(ArtifactType.SPECIES,),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(ConformersConfig, config)
        engine = cast(bp.ConformerEngine, rt.engine(bp.Capability.CONFORMERS, cfg.engine))
        search = partial(_search, engine, rt.method(cfg.method), rt.load_xyz)
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)
                   if s.source == "input"}
        # the JobRunner core semaphore limits concurrent CREST jobs by the site threads
        found = rt.thread_map(search, _monomers(species, rt, cfg))
        best = {f.item.base: _lowest(f, rt) for f in found}
        built = [_composition(c, species, best, rt, cfg) for c in rt.system.compositions]
        items = [b for b in built if isinstance(b, _Item)]
        found += rt.thread_map(search, items)
        failed = [b for b in built if isinstance(b, Artifact)]
        return failed + [a for f in found for a in _artifacts(f, rt, cfg)]
