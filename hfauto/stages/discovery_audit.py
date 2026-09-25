"""Publish bounded-search coverage without claiming chemical exhaustiveness."""

from __future__ import annotations

import csv
from typing import Any

from hfauto.core.io import (
    ensure_dir,
    read_manifest,
    write_json,
    write_jsonl,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.discovery_coverage import (
    build_discovery_coverage,
    build_discovery_saturation,
)


class DiscoveryAuditStage(Stage):
    """Summarize evidence attrition from generated trials through validated IRC."""

    name = "discovery-audit"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("discovery-audit requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        rows, summary = build_discovery_coverage(manifest)
        jsonl_path = write_jsonl(rows, out_dir / "discovery_coverage.jsonl")
        summary_path = write_json(out_dir / "discovery_coverage_summary.json", summary)
        csv_path = out_dir / "discovery_coverage.csv"
        fields = sorted({key for row in rows for key in row})
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if fields:
                writer.writeheader()
                writer.writerows(rows)
        out.add_artifact(
            Artifact(
                artifact_id="discovery_coverage",
                artifact_type="discovery_coverage",
                parents=[
                    item.artifact_id
                    for artifact_type in (
                        "reaction_trial",
                        "reaction_discovery_attempt",
                        "reaction_candidate",
                        "minimum_registry",
                        "reaction",
                        "reaction_validated",
                        "reaction_path_validated",
                    )
                    for item in manifest.latest_artifacts(artifact_type)
                ],
                paths={
                    "csv": str(csv_path),
                    "jsonl": str(jsonl_path),
                    "summary_json": str(summary_path),
                },
                data=summary,
                qc={
                    "finite_search_only": True,
                    "exhaustive_claim_allowed": False,
                    "source_rows_present": bool(rows),
                },
                provenance={"created_by": self.name},
            )
        )
        comparison_paths = [
            str(path) for path in config.get("comparison_manifests", [])
        ]
        if comparison_paths:
            comparison_manifests = [
                read_manifest(path) for path in comparison_paths
            ]
            campaigns = [
                (item.run_id, item) for item in comparison_manifests
            ]
            if not any(item.run_id == manifest.run_id for _label, item in campaigns):
                campaigns.append((manifest.run_id, manifest))
            saturation_rows, saturation_summary = build_discovery_saturation(
                campaigns
            )
            saturation_jsonl = write_jsonl(
                saturation_rows, out_dir / "discovery_saturation.jsonl"
            )
            saturation_summary_path = write_json(
                out_dir / "discovery_saturation_summary.json",
                saturation_summary,
            )
            saturation_csv = out_dir / "discovery_saturation.csv"
            saturation_fields = sorted(
                {key for row in saturation_rows for key in row}
            )
            with saturation_csv.open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=saturation_fields
                )
                if saturation_fields:
                    writer.writeheader()
                    writer.writerows(saturation_rows)
            out.add_artifact(
                Artifact(
                    artifact_id="discovery_saturation",
                    artifact_type="discovery_saturation",
                    parents=["discovery_coverage"],
                    paths={
                        "csv": str(saturation_csv),
                        "jsonl": str(saturation_jsonl),
                        "summary_json": str(saturation_summary_path),
                    },
                    data=saturation_summary,
                    qc={
                        "finite_campaign_comparison_only": True,
                        "exhaustive_claim_allowed": False,
                    },
                    provenance={
                        "created_by": self.name,
                        "comparison_manifests": comparison_paths,
                    },
                )
            )
        return out
