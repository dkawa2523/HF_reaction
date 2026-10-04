"""structures stage (design §4.1 #1, §8.2): system species → ``species`` artifacts.

xyz files are copied into the stage directory and SMILES are embedded once (RDKit ETKDG).
The multiplicity is always the declared one. Elements and the parity of the declared charge
and multiplicity are checked once here with ``check_electronic_state``; both ends of a
declared reaction must share composition id (formula, charge, multiplicity) and atom order,
and its coordinate indices must fall inside the molecule. Ends that differ only in
multiplicity are a spin crossing (``spin_crossing_reaction_unsupported``: no MECP search; each
surface can be declared as its own reaction). Every problem becomes a failed artifact with
``INPUT_INVALID``.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import ClassVar, cast

from hfauto.chemistry.electronic_state import check_electronic_state
from hfauto.chemistry.smiles import smiles_to_xyz
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import composition_key, read_xyz, write_xyz, written_geometry
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.ids import species_artifact_id
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import ArtifactType, SpeciesRecord
from hfauto.core.system import SpeciesInput
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class StructuresConfig(StageConfig):
    pass


def _failed(species_id: str, kind: FailureKind, reason: str) -> Artifact:
    return Artifact(artifact_id=species_artifact_id(species_id), type=ArtifactType.SPECIES,
                    status="failed", failure=Failure(kind=kind, reason=reason))


def _write(species: SpeciesInput, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if species.xyz is not None:  # SpeciesInput guarantees exactly one of xyz and smiles
        shutil.copyfile(species.xyz, target)
    else:
        write_xyz(smiles_to_xyz(cast(str, species.smiles), species.charge), target)


def _species(species: SpeciesInput, rt: StageRuntime) -> Artifact:
    path = rt.stage_dir / "xyz" / f"{species.id}.xyz"
    mult = species.multiplicity
    try:
        _write(species, path)
        geometry, xyz = written_geometry(path, rt.file_ref), read_xyz(path)
        check_electronic_state(xyz.symbols, species.charge, mult)
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
    return Artifact(artifact_id=species_artifact_id(species.id), type=ArtifactType.SPECIES,
                    payload=record)


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
            if a.charge != b.charge or a.geometry.symbols != b.geometry.symbols:
                reason = f"reaction {reaction.id}: endpoints differ in charge or atom order"
            elif a.multiplicity != b.multiplicity:
                reason = (f"spin_crossing_reaction_unsupported: reaction {reaction.id} joins "
                          f"multiplicities {a.multiplicity} and {b.multiplicity}")
            elif any(i >= len(a.geometry.symbols) for t in reaction.coordinate for i in t.atoms):
                reason = f"reaction {reaction.id}: coordinate atom index out of range"
            else:
                continue
            for end in (reaction.reactant, reaction.product):
                out[end] = _failed(end, FailureKind.INPUT_INVALID, reason)
        return list(out.values())
