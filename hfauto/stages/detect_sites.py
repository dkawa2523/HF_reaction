from __future__ import annotations

from typing import Any

from hfauto.chemistry.site_detection import detect_basic_sites_from_smiles
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class DetectSitesStage(Stage):
    name = "detect-sites"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        records: list[dict] = []
        include_excluded = bool(config.get("include_excluded", False))
        max_sites = int(config.get("max_sites_per_molecule", 9999))

        mol_artifacts = list(manifest.iter_artifacts("molecule_enriched")) or list(manifest.iter_artifacts("molecule"))
        for mol in mol_artifacts:
            smiles = mol.data.get("canonical_smiles")
            if not smiles:
                out.add_artifact(Artifact.failure(f"site_failed_{mol.artifact_id}", "site", "missing_smiles", parents=[mol.artifact_id]))
                continue
            sites = detect_basic_sites_from_smiles(smiles)
            sites = [s for s in sites if include_excluded or not s.get("excluded")]
            sites = sites[:max_sites]
            if not sites:
                out.add_artifact(
                    Artifact.failure(
                        f"site_none_{mol.artifact_id}",
                        "site",
                        "no_basic_site",
                        category="no_basic_site",
                        parents=[mol.artifact_id],
                        recoverable=False,
                    )
                )
                continue
            for site in sites:
                site_id = f"{mol.data['mol_id']}_{site['site_type']}_{site['atom_index']}"
                rec = {"site_id": site_id, "mol_id": mol.data["mol_id"], **site}
                records.append(rec)
                out.add_artifact(
                    Artifact(
                        artifact_id=site_id,
                        artifact_type="site",
                        parents=[mol.artifact_id],
                        data=rec,
                        qc={
                            "site_detected": True,
                            "site_confidence": site["site_confidence"],
                            "priority": site.get("priority"),
                            "excluded": site.get("excluded", False),
                        },
                    )
                )
        write_jsonl(records, out_dir / "site_records.jsonl")
        return out
