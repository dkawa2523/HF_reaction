from __future__ import annotations

from typing import Any

from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _preferred_reactant_species(manifest: Manifest) -> list[Artifact]:
    optimized: dict[str, Artifact] = {}
    for art in manifest.iter_artifacts("species_optimized"):
        if art.data.get("state") == "reactant_complex" and art.status.status == "success":
            optimized[_canonical_species_id(art)] = art
    preopt: dict[str, Artifact] = {}
    for art in manifest.iter_artifacts("species_preopt"):
        if art.data.get("state") == "reactant_complex" and art.status.status == "success":
            preopt[_canonical_species_id(art)] = art
    out: list[Artifact] = []
    seen: set[str] = set()
    for species in manifest.latest_artifacts("species"):
        if species.data.get("state") != "reactant_complex":
            continue
        cid = _canonical_species_id(species)
        if cid in seen:
            continue
        seen.add(cid)
        out.append(optimized.get(cid) or preopt.get(cid) or species)
    return out


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
            canonical_id = _canonical_species_id(species)
            desc = hf_descriptors_from_species(species.data)
            hf_stretches = [c.data.get("hf_stretch_cm1") for c in calcs_by_species.get(canonical_id, []) if c.data.get("hf_stretch_cm1")]
            if hf_stretches:
                desc["nu_HF_cm1"] = float(hf_stretches[-1])
                desc["delta_nu_HF_cm1"] = float(hf_stretches[-1]) - 4100.0
            real_orca = any(c.qc.get("real_orca_executed") for c in calcs_by_species.get(canonical_id, []))
            rec = {
                "species_id": canonical_id,
                "geometry_artifact_id": species.artifact_id,
                "geometry_artifact_type": species.artifact_type,
                "mol_id": species.data.get("mol_id"),
                "site_id": species.data.get("site_id"),
                "site_type": species.data.get("site_type"),
                "hf_n": species.data.get("hf_n"),
                "descriptor_level": "orca_optimized" if real_orca and species.artifact_type == "species_optimized" else species.artifact_type,
                **desc,
            }
            records.append(rec)
            out.add_artifact(
                Artifact(
                    artifact_id="desc_" + canonical_id,
                    artifact_type="descriptor",
                    parents=[species.artifact_id],
                    data=rec,
                    qc={"descriptor_status": "success" if desc else "partial", "geometry_source": species.artifact_type, "real_orca_support": real_orca},
                )
            )
        write_jsonl(records, out_dir / "descriptor_records.jsonl")
        return out
