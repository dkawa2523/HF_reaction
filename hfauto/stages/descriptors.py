from __future__ import annotations

from typing import Any

from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.core.artifacts import canonical_species_id, preferred_species
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _preferred_reactant_species(manifest: Manifest) -> list[Artifact]:
    return preferred_species(manifest, states={"reactant_complex"})


class DescriptorsStage(Stage):
    name = "descriptors"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        records: list[dict] = []
        calcs_by_species: dict[str, list[Artifact]] = {}
        for calc in manifest.iter_artifacts("calculation"):
            sid = calc.data.get("species_id")
            if sid:
                calcs_by_species.setdefault(str(sid), []).append(calc)
        for species in _preferred_reactant_species(manifest):
            canonical_id = canonical_species_id(species)
            desc = hf_descriptors_from_species(species.data)
            hf_stretches = [c.data.get("hf_stretch_cm1") for c in calcs_by_species.get(canonical_id, []) if c.data.get("hf_stretch_cm1")]
            if hf_stretches:
                desc["nu_HF_cm1"] = float(hf_stretches[-1])
                desc["delta_nu_HF_cm1"] = float(hf_stretches[-1]) - 4100.0
            real_qm = any(c.qc.get("real_qm_executed") or c.qc.get("real_orca_executed") for c in calcs_by_species.get(canonical_id, []))
            rec = {
                "species_id": canonical_id,
                "geometry_artifact_id": species.artifact_id,
                "geometry_artifact_type": species.artifact_type,
                "mol_id": species.data.get("mol_id"),
                "site_id": species.data.get("site_id"),
                "site_type": species.data.get("site_type"),
                "hf_n": species.data.get("hf_n"),
                "descriptor_level": "dft_optimized" if real_qm and species.artifact_type == "species_optimized" else species.artifact_type,
                **desc,
            }
            records.append(rec)
            out.add_artifact(
                Artifact(
                    artifact_id="desc_" + canonical_id,
                    artifact_type="descriptor",
                    parents=[species.artifact_id],
                    data=rec,
                    qc={"descriptor_status": "success" if desc else "partial", "geometry_source": species.artifact_type, "real_qm_support": real_qm},
                )
            )
        write_jsonl(records, out_dir / "descriptor_records.jsonl")
        return out
