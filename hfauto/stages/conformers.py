from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_conformer_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class ConformersStage(Stage):
    """Generate molecule conformer ensembles through a pluggable backend.

    Backends share the same ConformerRecord/Artifact contract:
      - rdkit: ETKDGv3 + MMFF/UFF local generation
      - crest: CREST/GFN-xTB search, with RDKit fallback when unavailable
    """

    name = "conformers"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine_spec = config.get("engine", "rdkit")
        engine_kwargs = config.get("engine_settings", {}) or {}
        engine = get_conformer_engine(engine_spec, **engine_kwargs)
        records: list[dict] = []
        mol_artifacts = list(manifest.iter_artifacts("molecule_enriched")) or list(manifest.iter_artifacts("molecule"))

        for mol_art in mol_artifacts:
            if mol_art.status.status == "failed":
                continue
            artifacts = engine.generate(mol_art, config, str(out_dir))
            for artifact in artifacts:
                out.add_artifact(artifact)
                if artifact.status.status == "success":
                    records.append(artifact.data)
                else:
                    records.append({"artifact_id": artifact.artifact_id, "status": "failed", "reason": artifact.status.reason})
        write_jsonl(records, out_dir / "conformer_records.jsonl")
        return out
