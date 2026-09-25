"""Enumerate bounded molecular states before conformer and reaction discovery."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.microstates import enumerate_molecular_states
from hfauto.core.ids import slug
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class EnumerateStatesStage(Stage):
    """Create unique tautomer/protomer records; do not predict populations."""

    name = "enumerate-states"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("enumerate-states requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        records: list[dict[str, Any]] = []
        include_tautomers = bool(config.get("include_tautomers", True))
        include_protomers = bool(config.get("include_protomers", True))
        max_states = int(config.get("max_states_per_molecule", 12))

        for molecule in manifest.latest_artifacts("molecule"):
            if molecule.status.status != "success":
                continue
            smiles = molecule.data.get("isomeric_smiles") or molecule.data.get(
                "canonical_smiles"
            )
            if not smiles:
                out.add_artifact(
                    Artifact.failure(
                        f"state_failed_{molecule.artifact_id}",
                        "molecule_state",
                        "molecule has no usable SMILES",
                        category="missing_input",
                        parents=[molecule.artifact_id],
                    )
                )
                continue
            try:
                candidates = enumerate_molecular_states(
                    str(smiles),
                    include_tautomers=include_tautomers,
                    include_protomers=include_protomers,
                    max_tautomers=int(config.get("max_tautomers", 8)),
                    max_protomers=int(config.get("max_protomers", 8)),
                )[:max_states]
            except (ImportError, RuntimeError, ValueError) as exc:
                out.add_artifact(
                    Artifact.failure(
                        f"state_failed_{molecule.artifact_id}",
                        "molecule_state",
                        str(exc),
                        category="state_enumeration_failed",
                        parents=[molecule.artifact_id],
                    )
                )
                continue

            from rdkit import Chem

            for index, candidate in enumerate(candidates):
                state_mol = Chem.AddHs(Chem.MolFromSmiles(candidate.smiles))
                charge = int(Chem.GetFormalCharge(state_mol))
                multiplicity = int(molecule.data.get("multiplicity", 1) or 1)
                symbols = [atom.GetSymbol() for atom in state_mol.GetAtoms()]
                try:
                    resolve_electronic_state(
                        {"charge": charge, "multiplicity": multiplicity}, symbols=symbols
                    )
                except ValueError:
                    # Do not guess an alternative spin state.  Open-shell state
                    # generation is an explicit extension point.
                    continue
                state_id = f"{slug(str(molecule.data['mol_id']))}_state{index:03d}"
                data = {
                    **molecule.data,
                    "state_id": state_id,
                    "mol_id": state_id,
                    "source_mol_id": molecule.data["mol_id"],
                    "source_molecule_artifact_id": molecule.artifact_id,
                    "canonical_smiles": candidate.smiles,
                    "isomeric_smiles": candidate.smiles,
                    "formal_charge": charge,
                    "multiplicity": multiplicity,
                    "state_kind": candidate.kind,
                    "source_atom_index": candidate.source_atom_index,
                    "state_population_claimed": False,
                }
                artifact = Artifact(
                    artifact_id=state_id,
                    artifact_type="molecule_state",
                    parents=[molecule.artifact_id],
                    data=data,
                    qc={
                        "chemically_sanitized": True,
                        "electron_state_validated": True,
                        "population_model": "not_evaluated",
                    },
                )
                out.add_artifact(artifact)
                records.append(data)
        write_jsonl(records, out_dir / "molecular_state_records.jsonl")
        return out
