"""structures stage (design §4.1 #1, §8.2): system species → ``species`` artifacts.

xyz files are copied into the stage directory and SMILES are embedded once (RDKit ETKDG).
Charge and multiplicity are checked with ``check_electronic_state``; both ends of a
declared reaction must list the same elements in the same order. Every problem becomes a
failed artifact with ``INPUT_INVALID``.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import ClassVar

from hfauto.chemistry.electronic_state import check_electronic_state
from hfauto.chemistry.smiles import smiles_to_molecule
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import XYZ, composition_key, geometry_fingerprint, read_xyz
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import ArtifactType, SpeciesRecord
from hfauto.core.system import SpeciesInput
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class StructuresConfig(StageConfig):
    pass


def artifact_id(species_id: str) -> str:
    return f"species_{species_id}"


def geometry_of(path: Path, rt: StageRuntime) -> tuple[Geometry, XYZ]:
    """Geometry of an xyz file in the run directory, fingerprinted as written."""
    xyz = read_xyz(path)
    fingerprint = geometry_fingerprint(xyz.symbols, xyz.coords)
    geometry = Geometry(file=rt.file_ref(path), fingerprint=fingerprint, symbols=tuple(xyz.symbols))
    return geometry, xyz


def _failed(species_id: str, kind: FailureKind, reason: str) -> Artifact:
    return Artifact(artifact_id=artifact_id(species_id), type=ArtifactType.SPECIES,
                    status="failed", failure=Failure(kind=kind, reason=reason))


def _write(species: SpeciesInput, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    if species.xyz is not None and species.smiles is None:
        shutil.copyfile(species.xyz, target)
        return target
    if species.smiles is not None and species.xyz is None:
        mol = smiles_to_molecule(species.smiles, species.charge, species.multiplicity)
        return mol.write(target)
    raise ValueError("give exactly one of xyz or smiles")


def _species(species: SpeciesInput, rt: StageRuntime) -> Artifact:
    try:
        geometry, xyz = geometry_of(_write(species, rt.stage_dir / "xyz" / f"{species.id}.xyz"), rt)
        check_electronic_state(xyz.symbols, species.charge, species.multiplicity)
        record = SpeciesRecord(
            species_id=species.id,
            composition_id=composition_key(xyz.symbols, species.charge, species.multiplicity),
            charge=species.charge, multiplicity=species.multiplicity, geometry=geometry,
            source="input", state_label=state_label(xyz.symbols, xyz.coords),
        )
    except ImportError as exc:
        return _failed(species.id, FailureKind.EXECUTABLE_MISSING, f"rdkit: {exc}")
    except (OSError, ValueError) as exc:
        return _failed(species.id, FailureKind.INPUT_INVALID, str(exc))
    return Artifact(artifact_id=artifact_id(species.id), type=ArtifactType.SPECIES, payload=record)


class StructuresStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="structures", config=StructuresConfig, consumes=(), produces=(ArtifactType.SPECIES,)
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        out = {s.id: _species(s, rt) for s in rt.system.species}
        for reaction in rt.system.reactions:
            a, b = out[reaction.reactant].payload, out[reaction.product].payload
            if not (isinstance(a, SpeciesRecord) and isinstance(b, SpeciesRecord)):
                continue  # a failed endpoint already carries its reason
            if a.geometry.symbols != b.geometry.symbols:
                reason = f"reaction {reaction.id}: endpoints differ in atom order"
                for end in (reaction.reactant, reaction.product):
                    out[end] = _failed(end, FailureKind.INPUT_INVALID, reason)
        return list(out.values())
