from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.core.io import ensure_dir, write_json
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto_ops.cache.registry import annotate_job_cache_status, duplicate_report, index_manifest_calculations
from hfauto_ops.core.run_index import run_summary
from hfauto_ops.reports.ops_report import render_ops_report
from hfauto_ops.retry.policy import build_retry_plan, write_retry_shell
from hfauto_ops.scheduler.resources import jobs_to_dataframe, make_stage_job_plan
from hfauto_ops.scheduler.templates import render_snakemake_scaffold


class OpsPlanStage(Stage):
    """Generate HPC/retry/cache planning artifacts without launching jobs."""

    name = "ops-plan"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        run_dir = Path(config.get("run_dir") or context.out_dir.parent)
        stages = list(config.get("plan_stages") or ["conformers", "preopt", "dft-minima", "ts-search", "irc", "sp", "thermo", "kinetics", "calibrate", "viz"])

        job_plan = jobs_to_dataframe(make_stage_job_plan(manifest, run_dir, stages, config=config))
        registry_db = out_dir / "calculation_registry.sqlite"
        cache_rows = index_manifest_calculations(manifest, registry_db, run_id=context.run_id)
        job_plan = annotate_job_cache_status(job_plan, registry_db)
        dup = duplicate_report(cache_rows, registry_db)
        retry = build_retry_plan(manifest, run_dir, config=config)
        retry_sh = write_retry_shell(retry, out_dir / "retry_commands.sh")
        scaffold_paths = render_snakemake_scaffold(job_plan, out_dir / "scheduler", config=config.get("scheduler", {}) or {})

        job_plan_path = out_dir / "hpc_job_plan.csv"
        retry_path = out_dir / "retry_plan.csv"
        dup_path = out_dir / "duplicate_report.csv"
        cache_path = out_dir / "calculation_cache_index.csv"
        summary_path = out_dir / "ops_summary.json"
        job_plan.to_csv(job_plan_path, index=False)
        retry.to_csv(retry_path, index=False)
        dup.to_csv(dup_path, index=False)
        cache_rows.to_csv(cache_path, index=False)
        summary = run_summary(manifest)
        summary.update({
            "n_planned_jobs": int(len(job_plan)),
            "n_retry_items": int(len(retry)),
            "n_cache_index_rows": int(len(cache_rows)),
            "n_duplicate_rows": int(len(dup)),
            "planned_stages": stages,
        })
        write_json(summary_path, summary)
        report_path = render_ops_report(out_dir / "ops_report.html", summary, job_plan, retry, dup)

        out.add_artifact(
            Artifact(
                artifact_id="hpc_job_plan",
                artifact_type="ops_job_plan",
                paths={"csv": str(job_plan_path)},
                data={
                    "n_jobs": int(len(job_plan)),
                    "planned_stages": stages,
                    "cache_hits": int((job_plan.get("cache_status") == "hit").sum()) if not job_plan.empty and "cache_status" in job_plan else 0,
                    "cache_misses": int((job_plan.get("cache_status") == "miss").sum()) if not job_plan.empty and "cache_status" in job_plan else 0,
                },
                qc={"review_required_before_submission": True},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="retry_plan",
                artifact_type="ops_retry_plan",
                paths={"csv": str(retry_path), "shell": str(retry_sh)},
                data={"n_retry_items": int(len(retry))},
                qc={"manual_review_required": bool(len(retry))},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="calculation_cache_index",
                artifact_type="ops_cache_index",
                paths={"csv": str(cache_path), "sqlite": str(registry_db)},
                data={"n_indexed": int(len(cache_rows)), "n_duplicates": int(len(dup))},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="duplicate_report",
                artifact_type="ops_duplicate_report",
                paths={"csv": str(dup_path)},
                data={"n_duplicates": int(len(dup))},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="scheduler_scaffold",
                artifact_type="ops_scheduler_scaffold",
                paths=scaffold_paths,
                data={"scheduler": "snakemake_slurm_pbs_scaffold", "n_jobs": int(len(job_plan))},
                qc={"submission_ready": False, "reason": "review_account_partition_modules_before_use"},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="ops_report",
                artifact_type="ops_report",
                paths={"html": str(report_path), "json": str(summary_path)},
                data=summary,
            )
        )
        return out
