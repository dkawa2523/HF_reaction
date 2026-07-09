from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.hpc.cache import build_cache_index, write_cache_index
from hfauto.hpc.compare import compare_runs, compare_backends
from hfauto.hpc.job_plan import planned_stage_jobs, retry_jobs_from_manifest, write_job_plan
from hfauto.hpc.schedulers import render_submit_script, write_snakemake_profile, write_snakefile_from_job_plan, write_submit_scripts
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_phase9_resource_job_plan_contract(tmp_path: Path):
    stages = [
        {"name": "ingest", "enabled": True},
        {"name": "dft-minima", "enabled": True},
        {"name": "ts-search", "enabled": True},
    ]
    jobs = planned_stage_jobs("rtest", stages, tmp_path / "runs" / "rtest")
    assert len(jobs) == 3
    assert jobs[1].ncores >= 8
    assert jobs[2].memory_gb >= jobs[1].memory_gb
    paths = write_job_plan(jobs, tmp_path / "ops")
    assert paths["csv"].exists()
    df = pd.read_csv(paths["csv"])
    assert {"job_id", "stage", "ncores", "memory_gb", "walltime", "command"} <= set(df.columns)


def test_phase9_scheduler_outputs(tmp_path: Path):
    jobs = planned_stage_jobs("rtest", [{"name": "ingest", "enabled": True}], tmp_path / "runs" / "rtest")
    paths = write_job_plan(jobs, tmp_path / "ops")
    scripts = write_submit_scripts(paths["csv"], tmp_path / "ops" / "submit", scheduler="slurm")
    assert any(p.name.startswith("submit_all") for p in scripts)
    text = scripts[0].read_text(encoding="utf-8")
    assert "#SBATCH" in text
    snakefile = write_snakefile_from_job_plan(paths["csv"], tmp_path / "ops" / "snakemake")
    assert "rule all" in snakefile.read_text(encoding="utf-8")
    profile = write_snakemake_profile(tmp_path / "ops" / "profile", scheduler="slurm")
    assert profile["config"].exists()
    assert "executor" in profile["config"].read_text(encoding="utf-8")
    assert "#PBS" in render_submit_script(jobs[0].to_dict(), scheduler="pbs")


def test_phase9_pipeline_hpc_plan_outputs(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase9_hpc_ops_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = Path(run_pipeline(cfg, "phase9_ops"))
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage in {"hpc-plan", "ops"}
    ops_dir = run_dir / ("15_hpc-plan" if (run_dir / "15_hpc-plan").exists() else "16_ops")
    assert ops_dir.exists()
    assert (ops_dir / "retry_plan.csv").exists()
    assert (ops_dir / "operations_report.html").exists() or (ops_dir / "hpc_dashboard.html").exists()
    assert list(manifest.iter_artifacts("hpc_bundle")) or list(manifest.iter_artifacts("operations_bundle"))


def test_phase9_cache_and_compare(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase7_full_publicdb_rank_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_a = Path(run_pipeline(cfg, "run_a"))
    run_b = Path(run_pipeline(cfg, "run_b"))
    ma = read_manifest(latest_manifest_path(run_a))
    df = build_cache_index(ma, run_dir=run_a)
    assert not df.empty
    paths = write_cache_index(ma, tmp_path / "cache", run_dir=run_a)
    assert paths["csv"].exists()
    comp = compare_runs(run_a, run_b, tmp_path / "compare_runs")
    assert comp["csv"].exists()
    backend = compare_backends(run_a, tmp_path / "compare_backends")
    assert backend["csv"].exists()


def test_phase9_retry_plan_uses_failed_artifacts(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase7_full_publicdb_rank_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = Path(run_pipeline(cfg, "retry_base"))
    manifest = read_manifest(latest_manifest_path(run_dir))
    # Inject a synthetic failure to test retry planning without needing a failed scientific run.
    from hfauto.core.schemas.artifact import Artifact
    manifest.add_artifact(Artifact.failure("synthetic_failed_ts", "calculation", "mock ts failure", category="ts_failed", recommended_fallback="scan_optts", data={"stage": "ts-search"}))
    jobs = retry_jobs_from_manifest(manifest, run_dir=run_dir, retry_root=tmp_path / "retry")
    assert len(jobs) >= 1
    assert jobs[0].retry_of == "synthetic_failed_ts"
