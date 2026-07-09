from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.rdkit_utils import embed_conformers_from_smiles, mol_to_xyz_block
from hfauto.core.schemas.artifact import Artifact


class RDKitConformerBackend:
    """Local RDKit ETKDGv3 + MMFF/UFF conformer backend."""

    name = "rdkit"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def generate(self, molecule: Artifact, config: dict, workdir: str) -> list[Artifact]:
        smiles = molecule.data.get("canonical_smiles")
        if not smiles:
            return [
                Artifact.failure(
                    f"conf_failed_{molecule.artifact_id}",
                    "conformer",
                    "missing canonical_smiles",
                    category="missing_input",
                    parents=[molecule.artifact_id],
                )
            ]
        cfg = {**self.config, **(config or {})}
        try:
            confs = embed_conformers_from_smiles(
                smiles,
                max_conformers=int(cfg.get("n_embed", cfg.get("max_conformers", 20))),
                seed=int(cfg.get("seed", 20260708)),
                prune_rms_thresh=float(cfg.get("prune_rms_thresh", 0.5)),
                optimize=bool(cfg.get("optimize", True)),
            )
        except Exception as exc:
            return [
                Artifact.failure(
                    f"conf_failed_{molecule.artifact_id}",
                    "conformer",
                    str(exc),
                    category="conformer_generation_failed",
                    parents=[molecule.artifact_id],
                )
            ]
        energy_window = float(cfg.get("energy_window_kcal_mol", 8.0))
        selected = [c for c in confs if c.relative_energy_kcal_mol <= energy_window]
        if not selected and confs:
            selected = confs[:1]
        selected_max = int(cfg.get("selected_max", cfg.get("max_conformers", len(selected))))
        selected = selected[:selected_max]

        out: list[Artifact] = []
        base = Path(workdir) / molecule.data["mol_id"]
        base.mkdir(parents=True, exist_ok=True)
        for idx, conf in enumerate(selected):
            conformer_id = f"{molecule.data['mol_id']}_conf{idx:04d}"
            xyz_path = base / f"conf{idx:04d}.xyz"
            comment = (
                f"generated_by=hfauto source=rdkit conf_id={conf.conf_id} "
                f"relE_kcal_mol={conf.relative_energy_kcal_mol:.6f}"
            )
            xyz_path.write_text(mol_to_xyz_block(conf.mol, conf.conf_id, comment=comment), encoding="utf-8")
            rec = {
                "conformer_id": conformer_id,
                "mol_id": molecule.data["mol_id"],
                "source": self.name,
                "rdkit_conf_id": int(conf.conf_id),
                "forcefield_energy_kcal_mol": float(conf.energy_kcal_mol),
                "relative_energy_kcal_mol": float(conf.relative_energy_kcal_mol),
                "boltzmann_weight_298K": float(conf.boltzmann_weight_298K),
                "xyz_path": str(xyz_path),
                "selected_for_hf_build": True,
                "rmsd_cluster_id": f"rdkit_rank_{idx:04d}",
                "formal_charge": int(molecule.data.get("formal_charge", 0) or 0),
                "multiplicity": int(molecule.data.get("multiplicity", 1) or 1),
                "extras": {"n_generated": len(confs), "n_selected": len(selected)},
            }
            out.append(
                Artifact(
                    artifact_id=conformer_id,
                    artifact_type="conformer",
                    parents=[molecule.artifact_id],
                    paths={"xyz": str(xyz_path)},
                    data=rec,
                    method={"engine": self.name, "algorithm": "ETKDGv3+MMFF/UFF"},
                    qc={
                        "generated": True,
                        "relative_energy_kcal_mol": float(conf.relative_energy_kcal_mol),
                        "boltzmann_weight_298K": float(conf.boltzmann_weight_298K),
                        "conformer_backend": self.name,
                    },
                )
            )
        return out
