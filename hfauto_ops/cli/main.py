from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from rich.console import Console

from hfauto.core.config import load_yaml
from hfauto_ops.bundle import build_operations_bundle
from hfauto_ops.core.artifact_jobs import write_artifact_job_plan
from hfauto_ops.core.autotune import build_resource_autotune_suggestions, write_resource_autotune
from hfauto_ops.core.backend_compare import write_backend_comparison
from hfauto_ops.core.compare import compare_runs
from hfauto_ops.core.qcarchive import write_qcarchive_bundle
from hfauto_ops.core.resource import build_resource_plan, write_resource_plan
from hfauto_ops.core.retry import write_retry_plan
from hfauto_ops.core.reuse import write_reuse_plan
from hfauto_ops.core.run_index import write_index_bundle
from hfauto_ops.schedulers.controller import scheduler_cancel, scheduler_status, submit_scripts
from hfauto_ops.schedulers.portable import render_scheduler_bundle
from hfauto_ops.schedulers.slurm import render_slurm_artifact_array, render_slurm_bundle
from hfauto_ops.status.snapshot import write_status_snapshot

app = typer.Typer(help="HF Auto operations, HPC planning, retry, reuse, scheduler, archive, and run comparison tools")
console = Console()


def _print_paths(paths: dict[str, Path]) -> None:
    console.print_json(data={k: str(v) for k, v in paths.items()})


@app.command("index")
def index_cmd(run_dir: Path, out: Path = typer.Option(..., help="Output ops directory")):
    _print_paths(write_index_bundle(run_dir, out))


@app.command("plan")
def plan_cmd(run_dir: Path, pipeline_config: Path = typer.Option(..., "--pipeline-config"), out: Path = typer.Option(...), ops_config: Optional[Path] = None):
    cfg = load_yaml(ops_config)
    _print_paths(write_resource_plan(build_resource_plan(pipeline_config, run_dir=run_dir, ops_config=cfg), out))


@app.command("artifact-plan")
def artifact_plan_cmd(run_dir: Path, out: Path = typer.Option(...), pipeline_config: Optional[Path] = typer.Option(None, "--pipeline-config"), ops_config: Optional[Path] = None):
    cfg = load_yaml(ops_config)
    _print_paths(write_artifact_job_plan(run_dir, out, pipeline_config=pipeline_config, ops_config=cfg))


@app.command("scheduler-bundle")
def scheduler_bundle_cmd(run_dir: Path, pipeline_config: Path = typer.Option(..., "--pipeline-config"), out: Path = typer.Option(...), scheduler: str = "slurm", run_id: Optional[str] = None, ops_config: Optional[Path] = None, artifact_level: bool = False):
    cfg = load_yaml(ops_config)
    if artifact_level:
        ap = write_artifact_job_plan(run_dir, out, pipeline_config=pipeline_config, ops_config=cfg)
        df = pd.read_csv(ap["artifact_job_csv"] if "artifact_job_csv" in ap else ap["artifact_jobs_csv"])
        if not df.empty and "stage" not in df.columns:
            df = df.rename(columns={"target_stage": "stage", "artifact_job_id": "job_id"})
        paths = render_scheduler_bundle(df, out, scheduler=scheduler, run_id=run_id or run_dir.name, project_root=cfg.get("project_root", "$PWD"), prefix="artifact")
    else:
        plan = build_resource_plan(pipeline_config, run_dir=run_dir, ops_config=cfg)
        write_resource_plan(plan, out)
        paths = render_scheduler_bundle(plan, out, scheduler=scheduler, run_id=run_id or run_dir.name, project_root=cfg.get("project_root", "$PWD"), prefix="stage")
    _print_paths(paths)


@app.command("slurm")
def slurm_cmd(run_dir: Path, pipeline_config: Path = typer.Option(..., "--pipeline-config"), out: Path = typer.Option(...), run_id: Optional[str] = None, ops_config: Optional[Path] = None):
    cfg = load_yaml(ops_config)
    plan = build_resource_plan(pipeline_config, run_dir=run_dir, ops_config=cfg)
    write_resource_plan(plan, out)
    _print_paths(render_slurm_bundle(plan, out, run_id=run_id or run_dir.name, project_root=cfg.get("project_root", "$PWD")))


@app.command("slurm-artifact-array")
def slurm_artifact_array_cmd(artifact_job_plan_csv: Path, out: Path = typer.Option(...), run_id: str = "run", project_root: str = "$PWD", chunk_size: int = 1):
    _print_paths(render_slurm_artifact_array(artifact_job_plan_csv, out, run_id=run_id, project_root=project_root, chunk_size=chunk_size))


@app.command("submit")
def submit_cmd(ops_dir: Path, scheduler: str = "slurm", dry_run: bool = True, allow_execute: bool = False, limit: Optional[int] = None):
    _print_paths(submit_scripts(ops_dir, scheduler=scheduler, dry_run=dry_run, allow_execute=allow_execute, limit=limit))


@app.command("scheduler-status")
def scheduler_status_cmd(ops_dir: Path, scheduler: str = "slurm", job_id: Optional[str] = None, user: Optional[str] = None, dry_run: bool = True, allow_execute: bool = False):
    _print_paths(scheduler_status(ops_dir, scheduler=scheduler, job_id=job_id, user=user, dry_run=dry_run, allow_execute=allow_execute))


@app.command("scheduler-cancel")
def scheduler_cancel_cmd(ops_dir: Path, scheduler: str = "slurm", job_id: Optional[str] = None, job_file: Optional[Path] = None, dry_run: bool = True, allow_execute: bool = False):
    _print_paths(scheduler_cancel(ops_dir, scheduler=scheduler, job_id=job_id, job_file=job_file, dry_run=dry_run, allow_execute=allow_execute))


@app.command("status-snapshot")
def status_snapshot_cmd(plan_csv: Path = typer.Option(...), out: Path = typer.Option(...), scheduler: str = "slurm"):
    _print_paths(write_status_snapshot(plan_csv, out, scheduler=scheduler))


@app.command("retry")
def retry_cmd(run_dir: Path, out: Path = typer.Option(...), pipeline_config: Optional[Path] = typer.Option(None, "--pipeline-config")):
    _print_paths(write_retry_plan(run_dir, out, pipeline_config))


@app.command("reuse-plan")
def reuse_plan_cmd(run_dir: Path, out: Path = typer.Option(...), cache_run: list[Path] = typer.Option([], "--cache-run")):
    _print_paths(write_reuse_plan(run_dir, out, cache_runs=cache_run))


@app.command("qcarchive-plan")
def qcarchive_plan_cmd(artifact_job_plan_csv: Path, out: Path = typer.Option(...), allow_submit: bool = False):
    _print_paths(write_qcarchive_bundle(artifact_job_plan_csv, out, allow_submit=allow_submit))


@app.command("backend-compare")
def backend_compare_cmd(run_dir: Path, out: Path = typer.Option(...)):
    _print_paths(write_backend_comparison(run_dir, out))


@app.command("autotune")
def autotune_cmd(run_dir: Path, out: Path = typer.Option(...), pipeline_config: Optional[Path] = typer.Option(None, "--pipeline-config"), ops_config: Optional[Path] = None):
    cfg = load_yaml(ops_config)
    df = build_resource_autotune_suggestions(run_dir, pipeline_config, cfg)
    _print_paths(write_resource_autotune(df, out))


@app.command("compare-runs")
def compare_runs_cmd(run_a: Path, run_b: Path, out: Path = typer.Option(...)):
    _print_paths(compare_runs(run_a, run_b, out))


@app.command("bundle")
def bundle_cmd(run_dir: Path, out: Path = typer.Option(...), pipeline_config: Optional[Path] = typer.Option(None, "--pipeline-config"), ops_config: Optional[Path] = None):
    cfg = load_yaml(ops_config)
    _print_paths(build_operations_bundle(run_dir, out, pipeline_config=pipeline_config, ops_config=cfg))


if __name__ == "__main__":
    app()
