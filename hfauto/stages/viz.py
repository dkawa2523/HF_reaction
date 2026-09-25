from __future__ import annotations

from typing import Any

from hfauto.core.io import ensure_dir
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class VizStage(Stage):
    """Generate visualization artifacts as an optional post-processing stage.

    The visualization code lives in the separate ``hfauto_viz`` package and is a
    read-only consumer of the run directory. This stage only bridges the normal
    ``hfauto pipeline`` contract to that companion package.
    """

    name = "viz"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        try:
            from hfauto_viz.reports.run_report import generate_visualizations
        except Exception as exc:  # pragma: no cover - dependency/import failure branch
            out.add_artifact(
                Artifact.failure(
                    "viz_import_failed",
                    "visualization_bundle",
                    f"hfauto_viz import failed: {exc}",
                    category="visualization_import_failed",
                    recoverable=True,
                    recommended_fallback="install hfauto_viz optional dependencies",
                )
            )
            return out

        run_dir = context.out_dir.parent
        settings = dict(config.get("settings") or {})
        # Allow concise top-level keys too.
        for key in [
            "max_molecule_dossiers",
            "max_reaction_dossiers",
            "top_energy_profiles",
            "max_network_reactions",
            "library_mode",
            "asset_mode",
            "ranking_thumbnail_rows",
        ]:
            if key in config and key not in settings:
                settings[key] = config[key]

        try:
            report_path, viz_manifest = generate_visualizations(run_dir, out_dir=context.out_dir, config=settings)
            out.add_artifact(
                Artifact(
                    artifact_id="visualization_bundle",
                    artifact_type="visualization_bundle",
                    parents=[manifest.manifest_id],
                    paths={
                        "html": str(report_path),
                        "viz_manifest": str(context.out_dir / "viz_manifest.json"),
                        "directory": str(context.out_dir),
                    },
                    data={
                        "visualization_package": "hfauto_viz",
                        "n_visualizations": len(getattr(viz_manifest, "records", [])),
                        "warnings": list(getattr(viz_manifest, "warnings", [])),
                    },
                    provenance={"stage": self.name, "created_at": Artifact.now_iso()},
                    qc={
                        "visualization_generated": True,
                        "read_only_consumer": True,
                        "contains_dummy_or_fallback_warning": bool(getattr(viz_manifest, "warnings", [])),
                    },
                )
            )
        except Exception as exc:
            out.add_artifact(
                Artifact.failure(
                    "viz_generation_failed",
                    "visualization_bundle",
                    f"visualization generation failed: {exc}",
                    category="visualization_failed",
                    recoverable=True,
                    recommended_fallback="run hfauto-viz report manually with --library-mode none",
                )
            )
        return out
