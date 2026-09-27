"""conformers stage (design §4.1 #2, §8.2): conformer search of monomers, placement seeds and
``--nci`` search of compositions, de-duplication, state labels and per state label the
``keep_per_state`` lowest structures (no energy window: the screen stage selects).

Small rigid monomers skip the search and pass their input through. A topology-change stop is kept
as a ``crest_topology`` species and the search reruns once from that structure with the same
settings (a second stop fails). Compositions take the summed charge of their components and the
declared multiplicity (or the only one spin coupling allows; otherwise INPUT_INVALID), and get
``--notopo`` on labile H and acceptor atoms when they have both. When a search fails its input
(monomer) or placement seeds (composition) are output instead. Users set quick, ewin_kcal,
seeds_per_composition and keep_per_state; CREST threads come from the site.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import partial
from typing import ClassVar, Literal, cast

from hfauto.backends import protocols as bp
from hfauto.chemistry import placement
from hfauto.chemistry.electronic_state import composition_multiplicity
from hfauto.chemistry.identity import permutation_invariant_rmsd
from hfauto.chemistry.topology import acceptor_atoms, labile_hydrogens, state_label
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    composition_key,
    geometry_fingerprint,
    read_xyz,
    write_xyz,
)
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, SpeciesRecord
from hfauto.core.system import CompositionInput
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

Source = Literal["conformer", "placement", "crest_topology"]
Candidate = tuple[Geometry, float | None, Source]
_TAG = {"conformer": "c", "placement": "p", "crest_topology": "t"}
DUPLICATE_RMSD_A = 0.1
DUPLICATE_DE_KCAL = 0.1


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
    diagnostics: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class _Cand:
    geometry: Geometry
    energy: float | None
    source: Source
    xyz: XYZ
    label: str


def artifact_id(species_id: str) -> str:
    return f"species_{species_id}"


def _failed(base: str, parents: tuple[str, ...], failure: Failure) -> Artifact:
    return Artifact(artifact_id=f"conformers_{base}", type=ArtifactType.SPECIES, parents=parents,
                    status="failed", failure=failure)


def _search(engine: bp.ConformerEngine, method: MethodSpec, load_xyz: Callable[[Geometry], XYZ],
            item: _Item) -> _Found:
    if item.skip:
        return _Found(item, list(item.fallback), diagnostics={"skipped": True})
    result, stops = engine.search(item.mol, method, item.settings), []
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:
        stops = [(g, None, cast(Source, "crest_topology")) for g in result.topology_stops]
        stop = Molecule(load_xyz(result.topology_stops[0]), item.mol.charge, item.mol.multiplicity)
        result = engine.search(stop, method, item.settings)
    if isinstance(result, bp.ConformerEnsemble) and result.topology_stops:  # no second retry
        result = Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason="topology_stop_repeated")
    diagnostics: dict[str, object] = {"topology_retry": bool(stops)}
    if isinstance(result, Failure):
        diagnostics["failure"] = f"{result.kind}: {result.reason}"
        return _Found(item, stops + list(item.fallback), result, diagnostics)
    diagnostics |= {"crest_version": result.version, "topology_removed": result.topology_removed}
    members = [(g, e, cast(Source, "conformer")) for g, e in result.members]
    return _Found(item, stops + members, None, diagnostics)


def _same(a: _Cand, b: _Cand) -> bool:
    de = None if a.energy is None or b.energy is None else abs(a.energy - b.energy)
    if a.xyz.symbols != b.xyz.symbols or (de or 0.0) * HARTREE_TO_KCAL_MOL >= DUPLICATE_DE_KCAL:
        return False
    rmsd, _ = permutation_invariant_rmsd(a.xyz.symbols, a.xyz.coords, b.xyz.coords)
    return rmsd < DUPLICATE_RMSD_A


def _select(cands: list[_Cand], keep_per_state: int) -> list[_Cand]:
    """Unique structures; per state label the lowest ``keep_per_state``. A missing energy is
    dropped when the state has energies (BUG-08); a state without energies keeps input order."""
    unique: list[_Cand] = []
    for c in sorted(cands, key=lambda c: (c.energy is None, c.energy or 0.0)):
        if not any(_same(c, u) for u in unique):
            unique.append(c)
    kept: list[_Cand] = []
    for label in dict.fromkeys(c.label for c in unique):
        group = [c for c in unique if c.label == label]
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
        out.append(Artifact(artifact_id=artifact_id(record.species_id), type=ArtifactType.SPECIES,
                            parents=item.parents, payload=record))
    return out


def _monomer(sp: SpeciesRecord, rt: StageRuntime, cfg: ConformersConfig) -> _Item:
    xyz = rt.load_xyz(sp.geometry)
    settings = bp.ConformerSettings(quick=cfg.quick, ewin_kcal=cfg.ewin_kcal)
    return _Item(sp.species_id, (artifact_id(sp.species_id),),
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
    geoms = []
    for k, seed in enumerate(seeds):
        path = write_xyz(seed, rt.stage_dir / "placement" / f"{comp_id}_seed{k:02d}.xyz")
        xyz = read_xyz(path)  # fingerprint the file as written
        geoms.append(Geometry(file=rt.file_ref(path), symbols=tuple(xyz.symbols),
                              fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords)))
    return geoms


def _notopo(xyz: XYZ) -> tuple[int, ...]:
    """Labile H and acceptor atoms when the structure has both, else none."""
    labile = labile_hydrogens(xyz.symbols, xyz.coords)
    acceptors = acceptor_atoms(xyz.symbols, xyz.coords)
    return tuple(sorted({*labile, *acceptors})) if labile and acceptors else ()


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
    settings = bp.ConformerSettings(nci=True, quick=cfg.quick, ewin_kcal=cfg.ewin_kcal,
                                    notopo_atoms=_notopo(first))
    parents = tuple(dict.fromkeys(artifact_id(sid) for sid, _ in parts))
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
        found = rt.thread_map(search, _monomers(species, rt, cfg), threads_per_item=1)
        best = {f.item.base: _lowest(f, rt) for f in found}
        built = [_composition(c, species, best, rt, cfg) for c in rt.system.compositions]
        items = [b for b in built if isinstance(b, _Item)]
        found += rt.thread_map(search, items, threads_per_item=1)
        diagnostics = json.dumps({f.item.base: f.diagnostics for f in found}, indent=1)
        (rt.stage_dir / "diagnostics.json").write_text(diagnostics, encoding="utf-8")
        failed = [b for b in built if isinstance(b, Artifact)]
        return failed + [a for f in found for a in _artifacts(f, rt, cfg)]
