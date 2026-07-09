from __future__ import annotations

from typing import Any

from hfauto.chemistry.hf_builder import build_hf_complexes, copy_candidate_species, write_hf_cluster
from hfauto.core.ids import reaction_id as make_reaction_id
from hfauto.core.ids import species_id as make_species_id
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class BuildHFStage(Stage):
    """Build isolated candidate, HF cluster, reactant-complex and ion-pair endpoints."""

    name = "build-hf"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        hf_ns = [int(x) for x in config.get("hf_n", [1, 2, 3])]
        species_records: list[dict] = []
        reaction_records: list[dict] = []

        conformers = [c for c in manifest.iter_artifacts("conformer") if c.status.status == "success"]
        sites = [s for s in manifest.iter_artifacts("site") if s.status.status == "success"]
        sites_by_mol: dict[str, list[Artifact]] = {}
        for site in sites:
            if site.data.get("excluded") is True and not config.get("include_excluded_sites", False):
                continue
            sites_by_mol.setdefault(site.data["mol_id"], []).append(site)

        hf_cluster_ids: dict[int, str] = {}
        for n in hf_ns:
            hf_id = make_species_id("hf_cluster", f"hf{n}")
            hf_cluster_ids[n] = hf_id
            hf_xyz = out_dir / "hf_clusters" / f"hf{n}.xyz"
            hf_qc = write_hf_cluster(n, hf_xyz)
            rec = {
                "species_id": hf_id,
                "mol_id": None,
                "site_id": None,
                "conformer_id": None,
                "state": "hf_cluster",
                "hf_n": n,
                "charge": 0,
                "multiplicity": 1,
                "xyz_path": str(hf_xyz),
                "atom_order_key": f"hf{n}_cluster_v1",
                "components": [{"name": f"HF_{i+1}", "atom_indices": [2 * i, 2 * i + 1]} for i in range(n)],
            }
            species_records.append(rec)
            out.add_artifact(Artifact(artifact_id=hf_id, artifact_type="species", paths={"xyz": str(hf_xyz)}, data=rec, qc={"geometry_qc": hf_qc}))

        for conf in conformers:
            mol_id = conf.data["mol_id"]
            cand_id = make_species_id(mol_id, conf.artifact_id, "candidate")
            cand_xyz = out_dir / mol_id / "candidate" / f"{conf.artifact_id}.xyz"
            cand_qc = copy_candidate_species(conf.paths["xyz"], cand_xyz)
            cand_rec = {
                "species_id": cand_id,
                "mol_id": mol_id,
                "site_id": None,
                "conformer_id": conf.artifact_id,
                "state": "bare_candidate",
                "hf_n": 0,
                "charge": int(conf.data.get("formal_charge", 0) or 0),
                "multiplicity": int(conf.data.get("multiplicity", 1) or 1),
                "xyz_path": str(cand_xyz),
                "atom_order_key": f"{mol_id}_{conf.artifact_id}_candidate_v1",
                "components": [{"name": "B", "atom_indices": list(range(int(cand_qc["n_atoms"]))) }],
            }
            species_records.append(cand_rec)
            out.add_artifact(Artifact(artifact_id=cand_id, artifact_type="species", parents=[conf.artifact_id], paths={"xyz": str(cand_xyz)}, data=cand_rec, qc={"geometry_qc": cand_qc}))

            for site in sites_by_mol.get(mol_id, []):
                for n in hf_ns:
                    stem = f"{mol_id}_{site.artifact_id}_hf{n}_{conf.artifact_id}"
                    atom_order_key = stem + "_map_v1"
                    rc_id = make_species_id(stem, "rc")
                    ip_id = make_species_id(stem, "ip")
                    rc_xyz = out_dir / mol_id / site.artifact_id / f"hf{n}" / f"{conf.artifact_id}_rc.xyz"
                    ip_xyz = out_dir / mol_id / site.artifact_id / f"hf{n}" / f"{conf.artifact_id}_ip.xyz"
                    try:
                        rc_qc, ip_qc = build_hf_complexes(conf.paths["xyz"], int(site.data["atom_index"]), n, rc_xyz, ip_xyz)
                    except Exception as exc:
                        out.add_artifact(Artifact.failure(f"build_hf_failed_{stem}", "species", str(exc), category="hf_build_failed", parents=[conf.artifact_id, site.artifact_id], recommended_fallback="review_site_atom_index_or_conformer"))
                        continue
                    base_atom_count = int(rc_qc["transfer_h"])
                    hf_components = [{"name": f"HF_{i+1}", "atom_indices": [base_atom_count + 2 * i, base_atom_count + 2 * i + 1]} for i in range(n)]
                    reaction_coordinate = {
                        "type": "distance_difference",
                        "definition": "r(B-H)-r(H-F)",
                        "atoms": {"base_atom": rc_qc["base_atom"], "transfer_h": rc_qc["transfer_h"], "leaving_f": rc_qc["leaving_f"]},
                    }
                    for state, sid, xyz_path, qc in [("reactant_complex", rc_id, rc_xyz, rc_qc), ("ion_pair", ip_id, ip_xyz, ip_qc)]:
                        rec = {
                            "species_id": sid,
                            "mol_id": mol_id,
                            "site_id": site.artifact_id,
                            "site_type": site.data.get("site_type"),
                            "conformer_id": conf.artifact_id,
                            "state": state,
                            "hf_n": n,
                            "charge": int(conf.data.get("formal_charge", 0) or 0),
                            "multiplicity": int(conf.data.get("multiplicity", 1) or 1),
                            "xyz_path": str(xyz_path),
                            "atom_order_key": atom_order_key,
                            "components": [{"name": "B", "atom_indices": list(range(base_atom_count))}, *hf_components],
                            "reaction_coordinate": reaction_coordinate,
                        }
                        species_records.append(rec)
                        out.add_artifact(Artifact(artifact_id=sid, artifact_type="species", parents=[conf.artifact_id, site.artifact_id], paths={"xyz": str(xyz_path)}, data=rec, qc={"geometry_qc": qc}))

                    rxn_id = make_reaction_id(mol_id, site.artifact_id, f"hf{n}", conf.artifact_id, "PT")
                    rxn = {
                        "reaction_id": rxn_id,
                        "reaction_type": "proton_transfer",
                        "mol_id": mol_id,
                        "site_id": site.artifact_id,
                        "site_type": site.data.get("site_type"),
                        "conformer_id": conf.artifact_id,
                        "hf_n": n,
                        "candidate_species_id": cand_id,
                        "hf_cluster_species_id": hf_cluster_ids[n],
                        "reactant_species_id": rc_id,
                        "product_species_id": ip_id,
                        "ts_species_id": None,
                        "atom_order_key": atom_order_key,
                        "reaction_coordinate": reaction_coordinate,
                    }
                    reaction_records.append(rxn)
                    out.add_artifact(Artifact(artifact_id=rxn_id, artifact_type="reaction", parents=[cand_id, hf_cluster_ids[n], rc_id, ip_id], data=rxn))

        write_jsonl(species_records, out_dir / "species_records.jsonl")
        write_jsonl(reaction_records, out_dir / "reaction_seed_records.jsonl")
        return out
