from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from hfauto.core.hashing import sha256_file
from hfauto.core.ids import mol_id_from_index
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.schemas.records import MoleculeRecord
from hfauto.stages.base import Stage, StageContext


def _has_usable_3d(mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    if not conf.Is3D():
        return False
    coords = []
    for idx in range(min(mol.GetNumAtoms(), 10)):
        p = conf.GetAtomPosition(idx)
        coords.append((round(p.x, 6), round(p.y, 6), round(p.z, 6)))
    return len(set(coords)) > 1


class IngestStage(Stage):
    """Read candidate SDF and create canonical molecule records."""

    name = "ingest"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        sdf_path = Path(config["sdf"])
        out_dir = ensure_dir(context.out_dir)
        records: list[dict] = []
        artifacts: list[Artifact] = []
        sdf_hash = sha256_file(sdf_path)

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
        except Exception:
            Chem = None

        if Chem is None:
            raw = sdf_path.read_text(encoding="utf-8", errors="ignore")
            n = raw.count("$$$$")
            for idx in range(n):
                mol_id = mol_id_from_index(idx)
                rec = MoleculeRecord(
                    mol_id=mol_id,
                    source_sdf_index=idx,
                    name=mol_id,
                    canonical_smiles=None,
                    inchikey=None,
                    formal_charge=0,
                    multiplicity=1,
                    has_3d=False,
                    ingest_qc={"sanitize_ok": False, "reason": "rdkit_not_available"},
                ).model_dump()
                records.append(rec)
                artifacts.append(
                    Artifact.failure(
                        mol_id,
                        "molecule",
                        "RDKit is not available; only record counting was performed.",
                        category="rdkit_not_available",
                        recoverable=True,
                        **rec,
                    )
                )
        else:
            suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
            writer = Chem.SDWriter(str(out_dir / "normalized.sdf"))
            rejects = Chem.SDWriter(str(out_dir / "rejects.sdf"))
            tmp_records: list[dict] = []
            tmp_artifacts: list[Artifact] = []
            for idx, mol in enumerate(suppl):
                mol_id = mol_id_from_index(idx)
                if mol is None:
                    rec = MoleculeRecord(
                        mol_id=mol_id,
                        source_sdf_index=idx,
                        name=mol_id,
                        formal_charge=None,
                        multiplicity=1,
                        has_3d=False,
                        ingest_qc={"sanitize_ok": False, "reject_reason": "rdkit_parse_failed"},
                    ).model_dump()
                    tmp_records.append(rec)
                    tmp_artifacts.append(
                        Artifact.failure(
                            mol_id,
                            "molecule",
                            "rdkit_parse_failed",
                            category="sanitize_failed",
                            recoverable=False,
                            **rec,
                        )
                    )
                    continue

                mol_no_h = Chem.RemoveHs(mol, sanitize=False)
                try:
                    Chem.SanitizeMol(mol_no_h)
                    sanitize_ok = True
                    reject_reason = None
                except Exception as exc:
                    sanitize_ok = False
                    reject_reason = f"sanitize_after_removeHs_failed: {exc}"

                smiles = Chem.MolToSmiles(mol_no_h, canonical=True) if sanitize_ok else None
                isomeric = Chem.MolToSmiles(mol_no_h, canonical=True, isomericSmiles=True) if sanitize_ok else None
                try:
                    inchi = Chem.MolToInchi(mol_no_h) if sanitize_ok else None
                    inchikey = Chem.InchiToInchiKey(inchi) if inchi else None
                except Exception:
                    inchi = None
                    inchikey = None

                frags = Chem.GetMolFrags(mol_no_h, asMols=False, sanitizeFrags=False) if sanitize_ok else ()
                n_frags = len(frags) if frags else 1
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else mol_id
                props = {p: mol.GetProp(p) for p in mol.GetPropNames()}
                mol.SetProp("hfauto_mol_id", mol_id)
                mol.SetProp("hfauto_source_sdf_sha256", sdf_hash)
                if sanitize_ok:
                    writer.write(mol)
                else:
                    rejects.write(mol)

                formula = rdMolDescriptors.CalcMolFormula(mol_no_h) if sanitize_ok else None
                exact_mw = float(Descriptors.ExactMolWt(mol_no_h)) if sanitize_ok else None
                qc = {
                    "sanitize_ok": sanitize_ok,
                    "valence_ok": sanitize_ok,
                    "reject_reason": reject_reason,
                    "mixture_or_salt_flag": n_frags > 1,
                    "needs_manual_charge_review": abs(int(Chem.GetFormalCharge(mol_no_h))) > 1 if sanitize_ok else True,
                }
                rec = MoleculeRecord(
                    mol_id=mol_id,
                    source_sdf_index=idx,
                    name=name,
                    canonical_smiles=smiles,
                    isomeric_smiles=isomeric,
                    inchi=inchi,
                    inchikey=inchikey,
                    formula=formula,
                    exact_mw=exact_mw,
                    formal_charge=int(Chem.GetFormalCharge(mol_no_h)) if sanitize_ok else None,
                    multiplicity=int(config.get("default_multiplicity", 1)),
                    num_atoms=int(mol.GetNumAtoms()),
                    num_fragments=n_frags,
                    has_3d=_has_usable_3d(mol),
                    sdf_props=props,
                    ingest_qc=qc,
                ).model_dump()
                tmp_records.append(rec)
                if sanitize_ok:
                    tmp_artifacts.append(
                        Artifact(
                            artifact_id=mol_id,
                            artifact_type="molecule",
                            paths={"source_sdf": str(sdf_path), "normalized_sdf": str(out_dir / "normalized.sdf")},
                            data=rec,
                            provenance={"source_sdf_sha256": sdf_hash},
                            qc=qc,
                        )
                    )
                else:
                    tmp_artifacts.append(
                        Artifact.failure(
                            mol_id,
                            "molecule",
                            reject_reason or "sanitize_failed",
                            category="sanitize_failed",
                            recoverable=False,
                            **rec,
                        )
                    )

            writer.close()
            rejects.close()

            inchikey_counts = Counter(r.get("inchikey") for r in tmp_records if r.get("inchikey"))
            duplicate_groups = {k: f"dup_{i:04d}" for i, (k, c) in enumerate(inchikey_counts.items()) if c > 1}
            for rec in tmp_records:
                dk = duplicate_groups.get(rec.get("inchikey"))
                rec["duplicate_group_id"] = dk
                rec["ingest_qc"]["duplicate_flag"] = dk is not None
            for art in tmp_artifacts:
                if art.data.get("mol_id"):
                    match = next((r for r in tmp_records if r["mol_id"] == art.data["mol_id"]), None)
                    if match:
                        art.data = match
                        art.qc.update(match.get("ingest_qc", {}))

            records.extend(tmp_records)
            artifacts.extend(tmp_artifacts)

        write_jsonl(records, out_dir / "molecule_records.jsonl")
        summary = {
            "n_records": len(records),
            "n_success": sum(1 for r in records if r.get("ingest_qc", {}).get("sanitize_ok")),
            "n_failed": sum(1 for r in records if not r.get("ingest_qc", {}).get("sanitize_ok")),
            "n_mixture_or_salt": sum(1 for r in records if r.get("ingest_qc", {}).get("mixture_or_salt_flag")),
            "n_duplicates": sum(1 for r in records if r.get("ingest_qc", {}).get("duplicate_flag")),
        }
        import json
        (out_dir / "ingest_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        out = Manifest.new(run_id=context.run_id, stage=self.name, metadata={"input_sdf": str(sdf_path), "summary": summary})
        out.extend(artifacts)
        return out
