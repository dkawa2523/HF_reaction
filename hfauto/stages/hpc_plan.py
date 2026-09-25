from __future__ import annotations

from typing import Any

from hfauto.core.io import ensure_dir
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.hpc.cache import build_cache_index, duplicate_cache_keys, write_cache_index
from hfauto.hpc.dashboard import render_hpc_dashboard
from hfauto.hpc.job_plan import (
    job_artifact,
    planned_stage_jobs,
    retry_jobs_from_manifest,
    write_job_plan,
)
from hfauto.hpc.schedulers import (
    write_snakefile_from_job_plan,
    write_snakemake_profile,
    write_submit_scripts,
)
from hfauto.stages.base import Stage, StageContext


class HPCPlanStage(Stage):
    """Generate reviewable large-scale execution artifacts without running jobs."""

    name = "hpc-plan"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        scheduler = str(config.get("scheduler", config.get("profile", "slurm")))
        planned_stages = config.get("planned_stages") or config.get("stages") or []
        run_dir = config.get("run_dir") or str(out_dir.parent)
        input_manifest = config.get("input_manifest")
        jobs = []
        if planned_stages:
            jobs = planned_stage_jobs(
                run_id=context.run_id,
                stages=planned_stages,
                run_dir=run_dir,
                input_manifest=input_manifest,
                global_config=context.global_config,
            )
        job_paths = write_job_plan(jobs, out_dir, prefix="job_plan")
        out.add_artifact(job_artifact(job_paths, len(jobs), kind="job_plan"))

        retry_root = out_dir / "retries"
        retry_jobs = retry_jobs_from_manifest(manifest, run_dir=run_dir, retry_root=retry_root, config=config)
        retry_paths = write_job_plan(retry_jobs, out_dir, prefix="retry_plan")
        out.add_artifact(job_artifact(retry_paths, len(retry_jobs), kind="retry_plan"))

        cache_paths = write_cache_index(manifest, out_dir, run_dir=run_dir)
        cache_df = build_cache_index(manifest, run_dir=run_dir)
        duplicate_df = duplicate_cache_keys(cache_df)
        duplicate_path = out_dir / "duplicate_cache_keys.csv"
        duplicate_df.to_csv(duplicate_path, index=False)
        out.add_artifact(Artifact(
            artifact_id="calculation_cache_index",
            artifact_type="cache_index",
            paths={"csv": str(cache_paths["csv"]), "json": str(cache_paths["json"]), "duplicates_csv": str(duplicate_path)},
            data={"n_rows": len(cache_df), "n_duplicate_rows": len(duplicate_df)},
        ))

        snakefile = write_snakefile_from_job_plan(job_paths["csv"], out_dir / "snakemake")
        profile_paths = write_snakemake_profile(out_dir / "snakemake_profile", scheduler=scheduler)
        submit_scripts = write_submit_scripts(job_paths["csv"], out_dir / "submit", scheduler=scheduler)
        retry_submit_scripts = write_submit_scripts(retry_paths["csv"], out_dir / "retry_submit", scheduler=scheduler)
        dashboard = render_hpc_dashboard(out_dir, job_plan_csv=job_paths["csv"], retry_plan_csv=retry_paths["csv"], cache_index_csv=cache_paths["csv"])
        out.add_artifact(Artifact(
            artifact_id="hpc_execution_bundle",
            artifact_type="hpc_bundle",
            paths={
                "snakefile": str(snakefile),
                "profile_config": str(profile_paths["config"]),
                "profile_manifest": str(profile_paths["manifest"]),
                "submit_dir": str(out_dir / "submit"),
                "retry_submit_dir": str(out_dir / "retry_submit"),
                "dashboard": str(dashboard),
            },
            data={
                "scheduler": scheduler,
                "n_submit_scripts": len(submit_scripts),
                "n_retry_submit_scripts": len(retry_submit_scripts),
                "n_jobs": len(jobs),
                "n_retry_jobs": len(retry_jobs),
                "scientific_note": "This stage plans execution only; it does not alter scientific artifacts.",
            },
        ))
        return out
