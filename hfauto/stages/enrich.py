from __future__ import annotations

from typing import Any

import pandas as pd

from hfauto.backends.registry import get_db_provider
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class EnrichStage(Stage):
    name = "enrich"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        providers = config.get("providers", ["dummy"])
        provider_config = config.get("provider_config", {})
        records = []
        summary_rows = []
        conflicts = []

        for mol in manifest.iter_artifacts("molecule"):
            if mol.status.status == "failed":
                continue
            data = mol.data.copy()
            provider_status: dict[str, Any] = {}
            for provider_spec in providers:
                if isinstance(provider_spec, dict):
                    provider_name = provider_spec.get("name") or provider_spec.get("provider")
                    kwargs = {k: v for k, v in provider_spec.items() if k not in {"name", "provider"}}
                else:
                    provider_name = provider_spec
                    kwargs = provider_config.get(provider_name, {}) if isinstance(provider_config, dict) else {}
                try:
                    provider = get_db_provider(provider_spec, **kwargs)
                    data = provider.enrich(Artifact(**{**mol.model_dump(), "data": data}))
                    status = data.get("db_provider_status", {}).get(provider_name, {"status": "success"})
                    provider_status[provider_name] = status
                except Exception as exc:
                    provider_status[provider_name or "unknown"] = {"status": "failed", "reason": str(exc), "matched": False}
                    data.setdefault("public_data", {}).setdefault(provider_name or "unknown", {"status": "failed", "reason": str(exc)})
            identity = data.get("identity", {}) or {}
            if identity.get("identity_conflict"):
                conflicts.append({
                    "mol_id": data.get("mol_id"),
                    "name": data.get("name"),
                    "inchikey": data.get("inchikey"),
                    "pubchem_cid": identity.get("pubchem_cid"),
                    "reason": "identity_conflict",
                })
            enriched_id = f"{mol.artifact_id}_enriched"
            enriched = Artifact(
                artifact_id=enriched_id,
                artifact_type="molecule_enriched",
                parents=[mol.artifact_id],
                data=data,
                qc={"db_enrichment_status": "success", "providers": providers, "provider_status": provider_status},
            )
            out.add_artifact(enriched)
            records.append(data)
            for provider_name, status in provider_status.items():
                summary_rows.append({
                    "mol_id": data.get("mol_id"),
                    "name": data.get("name"),
                    "provider": provider_name,
                    "status": status.get("status"),
                    "matched": bool(status.get("matched")),
                    "reason": status.get("reason"),
                })

        write_jsonl(records, out_dir / "enriched_molecule_records.jsonl")
        summary_path = out_dir / "public_db_enrichment_summary.csv"
        conflicts_path = out_dir / "identity_conflicts.csv"
        pd.DataFrame(summary_rows, columns=["mol_id", "name", "provider", "status", "matched", "reason"]).to_csv(summary_path, index=False)
        pd.DataFrame(conflicts, columns=["mol_id", "name", "inchikey", "pubchem_cid", "reason"]).to_csv(conflicts_path, index=False)
        out.add_artifact(Artifact(
            artifact_id="public_db_enrichment_summary",
            artifact_type="table",
            paths={"csv": str(summary_path)},
            data={"n_rows": len(summary_rows), "providers": [p.get("name") if isinstance(p, dict) else p for p in providers]},
        ))
        out.add_artifact(Artifact(
            artifact_id="identity_conflicts",
            artifact_type="table",
            paths={"csv": str(conflicts_path)},
            data={"n_rows": len(conflicts)},
            qc={"manual_review_required": bool(conflicts)},
        ))
        return out
