"""structures stage (design §4.1 #1, §8.2): system species → ``species`` artifacts.

xyz files are copied into the stage directory and SMILES are embedded once (RDKit ETKDG).
An undeclared multiplicity is the SMILES radical electrons + 1, or 1 for xyz. Elements,
charge and multiplicity are checked once here with ``check_electronic_state``; both ends of a
declared reaction must share composition id (formula, charge, multiplicity) and atom order,
and its coordinate indices must fall inside the molecule. Every problem becomes a failed
artifact with ``INPUT_INVALID``.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import ClassVar, cast

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


def _write(species: SpeciesInput, target: Path) -> int:
    """Write the input structure to ``target`` and return its multiplicity."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if species.xyz is not None:  # SpeciesInput guarantees exactly one of xyz and smiles
        shutil.copyfile(species.xyz, target)
        return 1 if species.multiplicity is None else species.multiplicity
    mol = smiles_to_molecule(cast(str, species.smiles), species.charge, species.multiplicity)
    mol.write(target)
    return mol.multiplicity


def _species(species: SpeciesInput, rt: StageRuntime) -> Artifact:
    path = rt.stage_dir / "xyz" / f"{species.id}.xyz"
    try:
        mult = _write(species, path)
        geometry, xyz = geometry_of(path, rt)
        check_electronic_state(xyz.symbols, species.charge, mult,
                               declared=species.multiplicity is not None)
        record = SpeciesRecord(
            species_id=species.id,
            composition_id=composition_key(xyz.symbols, species.charge, mult),
            charge=species.charge, multiplicity=mult, geometry=geometry,
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
            if a.composition_id != b.composition_id or a.geometry.symbols != b.geometry.symbols:
                reason = (f"reaction {reaction.id}: endpoints differ in charge, multiplicity "
                          "or atom order")
            elif any(i >= len(a.geometry.symbols) for t in reaction.coordinate for i in t.atoms):
                reason = f"reaction {reaction.id}: coordinate atom index out of range"
            else:
                continue
            for end in (reaction.reactant, reaction.product):
                out[end] = _failed(end, FailureKind.INPUT_INVALID, reason)
        return list(out.values())
