from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest, write_manifest
from hfauto.reporting.html_report import latest_manifest_path, render_report
from hfauto.stages.base import StageContext
from hfauto.stages.registry import get_stage
from hfauto.workflow.runner import run_pipeline
from hfauto.hpc.job_plan import planned_stage_jobs, retry_jobs_from_manifest, write_job_plan
from hfauto.hpc.cache import write_cache_index
from hfauto.hpc.compare import compare_runs as hpc_compare_runs, compare_backends as hpc_compare_backends
from hfauto.hpc.schedulers import write_submit_scripts, write_snakemake_profile, write_snakefile_from_job_plan
from hfauto.hpc.dashboard import render_hpc_dashboard

app = typer.Typer(help="HF gas-phase reactivity workflow")
console = Console()


@app.command()
def ingest(
    sdf: Path = typer.Option(..., help="Input candidates.sdf"),
    out: Path = typer.Option(..., help="Output stage directory"),
    run_id: str = typer.Option("manual", help="Run id"),
):
    stage = get_stage("ingest")
    context = StageContext(out_dir=out, run_id=run_id, global_config={})
    manifest = stage.run(None, {"sdf": str(sdf)}, context)
    path = write_manifest(manifest, out)
    console.print(f"[green]wrote[/green] {path}")


@app.command("run-stage")
def run_stage(
    stage_name: str = typer.Argument(..., help="Stage name, e.g. enrich, detect-sites, dft-minima"),
    in_manifest: Path = typer.Option(..., "--in", help="Input manifest.json"),
    out: Path = typer.Option(..., help="Output stage directory"),
    config: Optional[Path] = typer.Option(None, help="YAML config for this stage"),
    global_config: Optional[Path] = typer.Option(None, help="Pipeline/global YAML to inherit temperature, pressure, and mode settings"),
    run_id: Optional[str] = typer.Option(None, help="Override run id"),
):
    manifest = read_manifest(in_manifest)
    cfg = load_yaml(config)
    inherited = dict((manifest.metadata or {}).get("global_config") or {})
    if global_config:
        loaded = load_yaml(global_config)
        inherited.update(loaded.get("global", loaded))
    if isinstance(cfg.get("global"), dict):
        inherited.update(cfg.pop("global"))
    stage = get_stage(stage_name)
    context = StageContext(out_dir=out, run_id=run_id or manifest.run_id, global_config=inherited)
    out_manifest = stage.run(manifest, cfg, context)
    out_manifest.metadata.setdefault("global_config", inherited)
    path = write_manifest(out_manifest, out)
    console.print(f"[green]wrote[/green] {path}")


@app.command()
def pipeline(
    config: Path = typer.Option(..., help="Pipeline YAML"),
    run_id: str = typer.Option(..., help="Run id"),
    from_stage: Optional[str] = typer.Option(None, "--from", help="Start stage name for partial run"),
    to_stage: Optional[str] = typer.Option(None, "--to", help="End stage name for partial run"),
    start_manifest: Optional[Path] = typer.Option(None, help="Input manifest for partial run"),
):
    cfg = load_yaml(config)
    cfg["__config_path"] = str(config)
    run_dir = run_pipeline(cfg, run_id, from_stage=from_stage, to_stage=to_stage, start_manifest=start_manifest)
    console.print(f"[green]pipeline completed[/green] {run_dir}")


@app.command()
def status(run_dir: Path):
    path = latest_manifest_path(run_dir)
    manifest = read_manifest(path)
    console.print(f"run_id: {manifest.run_id}")
    console.print(f"latest_manifest: {path}")
    console.print(f"latest_stage: {manifest.stage}")
    counts = {}
    failures = 0
    for a in manifest.artifacts:
        counts[a.artifact_type] = counts.get(a.artifact_type, 0) + 1
        if a.status.status == "failed":
            failures += 1
    for k in sorted(counts):
        console.print(f"{k:20s} {counts[k]}")
    console.print(f"failures             {failures}")


@app.command()
def failures(run_dir: Path, limit: int = typer.Option(50, help="Maximum rows")):
    manifest = read_manifest(latest_manifest_path(run_dir))
    rows = [a for a in manifest.artifacts if a.status.status == "failed"]
    for a in rows[:limit]:
        console.print(f"{a.artifact_id}\t{a.artifact_type}\t{a.status.category}\t{a.status.reason}\t{a.status.recommended_fallback}")
    console.print(f"total_failures: {len(rows)}")


@app.command()
def inspect(run_dir: Path, artifact_id: str):
    manifest = read_manifest(latest_manifest_path(run_dir))
    artifact = manifest.find(artifact_id)
    if artifact is None:
        raise typer.BadParameter(f"Artifact not found: {artifact_id}")
    console.print_json(artifact.model_dump_json(indent=2))


@app.command()
def report(
    run_dir: Path,
    out: Path = typer.Option(..., help="Output HTML report path"),
):
    path = render_report(run_dir, out)
    console.print(f"[green]wrote[/green] {path}")


@app.command("hpc-plan")
def hpc_plan(
    run_dir: Path = typer.Argument(..., help="Existing run directory"),
    out: Path = typer.Option(..., help="Output operations directory"),
    config: Optional[Path] = typer.Option(None, help="Pipeline YAML containing stages to plan"),
    scheduler: str = typer.Option("slurm", help="slurm | pbs | local"),
):
    """Create job plans, submit script templates, Snakemake files, retry plans, and cache index."""
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest(latest_manifest_path(run_dir))
    cfg = load_yaml(config) if config else {}
    stages = cfg.get("stages", [])
    jobs = planned_stage_jobs(manifest.run_id, stages, run_dir=run_dir, global_config=cfg.get("global", {})) if stages else []
    job_paths = write_job_plan(jobs, out, prefix="job_plan")
    retry_jobs = retry_jobs_from_manifest(manifest, run_dir=run_dir, retry_root=out / "retries", config=cfg)
    retry_paths = write_job_plan(retry_jobs, out, prefix="retry_plan")
    cache_paths = write_cache_index(manifest, out, run_dir=run_dir)
    snakefile = write_snakefile_from_job_plan(job_paths["csv"], out / "snakemake")
    profile_paths = write_snakemake_profile(out / "snakemake_profile", scheduler=scheduler)
    write_submit_scripts(job_paths["csv"], out / "submit", scheduler=scheduler)
    write_submit_scripts(retry_paths["csv"], out / "retry_submit", scheduler=scheduler)
    dashboard = render_hpc_dashboard(out, job_plan_csv=job_paths["csv"], retry_plan_csv=retry_paths["csv"], cache_index_csv=cache_paths["csv"])
    console.print(f"[green]wrote[/green] {out}")
    console.print(f"job_plan: {job_paths['csv']}")
    console.print(f"retry_plan: {retry_paths['csv']}")
    console.print(f"cache_index: {cache_paths['csv']}")
    console.print(f"snakefile: {snakefile}")
    console.print(f"profile: {profile_paths['config']}")
    console.print(f"dashboard: {dashboard}")


@app.command("cache-index")
def cache_index(
    run_dir: Path = typer.Argument(..., help="Existing run directory"),
    out: Path = typer.Option(..., help="Output directory"),
):
    manifest = read_manifest(latest_manifest_path(run_dir))
    paths = write_cache_index(manifest, out, run_dir=run_dir)
    console.print(f"[green]wrote[/green] {paths['csv']}")


@app.command("retry-plan")
def retry_plan(
    run_dir: Path = typer.Argument(..., help="Existing run directory"),
    out: Path = typer.Option(..., help="Output directory"),
    scheduler: str = typer.Option("slurm", help="slurm | pbs | local"),
):
    manifest = read_manifest(latest_manifest_path(run_dir))
    jobs = retry_jobs_from_manifest(manifest, run_dir=run_dir, retry_root=out / "retries")
    paths = write_job_plan(jobs, out, prefix="retry_plan")
    write_submit_scripts(paths["csv"], out / "retry_submit", scheduler=scheduler)
    console.print(f"[green]wrote[/green] {paths['csv']} ({len(jobs)} retry jobs)")


@app.command("compare-runs")
def compare_runs_cmd(
    run_a: Path = typer.Argument(...),
    run_b: Path = typer.Argument(...),
    out: Path = typer.Option(..., help="Output comparison directory"),
):
    paths = hpc_compare_runs(run_a, run_b, out)
    console.print(f"[green]wrote[/green] {paths['csv']}")


@app.command("science-status")
def science_status(run_dir: Path):
    """Show simple scientific-readiness counts for the latest run."""
    manifest = read_manifest(latest_manifest_path(run_dir))
    console.print(f"run_id: {manifest.run_id}")
    rank_files = []
    for art in manifest.iter_artifacts("ranking"):
        if art.paths.get("csv"):
            rank_files.append(Path(art.paths["csv"]))
    summary_art = manifest.find("candidate_summary")
    if summary_art and summary_art.paths.get("csv"):
        df = pd.read_csv(summary_art.paths["csv"])
        console.print(f"candidates: {len(df)}")
        if "quality_tier" in df:
            console.print("quality tiers:")
            for tier, count in df["quality_tier"].fillna("NA").value_counts().sort_index().items():
                console.print(f"  {tier}: {count}")
        for col in ["scientific_rank_eligible", "production_rank_eligible", "main_values_are_dummy", "main_values_are_fallback"]:
            if col in df:
                console.print(f"{col}: {int(df[col].fillna(False).astype(bool).sum())}/{len(df)}")
        if "science_next_action" in df:
            console.print("science next actions:")
            for action, count in df["science_next_action"].fillna("unknown").value_counts().items():
                console.print(f"  {action}: {count}")
    else:
        console.print("candidate_summary.csv not found; run rank stage first")


@app.command("next-actions")
def next_actions(run_dir: Path, limit: int = typer.Option(20, help="Maximum candidate rows")):
    """Print concise next actions from candidate_summary.csv."""
    manifest = read_manifest(latest_manifest_path(run_dir))
    summary_art = manifest.find("candidate_summary")
    if not summary_art or not summary_art.paths.get("csv"):
        console.print("candidate_summary.csv not found; run rank stage first")
        return
    df = pd.read_csv(summary_art.paths["csv"]).head(limit)
    cols = [c for c in ["mol_id", "name", "quality_tier", "confidence_score", "production_rank_eligible", "science_next_action", "process_next_action", "ops_next_action"] if c in df.columns]
    console.print(df[cols].to_string(index=False))


@app.command("outputs")
def outputs(run_dir: Path):
    """Print the small set of files a normal user should open first."""
    manifest_path = latest_manifest_path(run_dir)
    manifest = read_manifest(manifest_path)
    candidates = [
        ("Latest manifest", manifest_path),
        ("Main HTML report", run_dir / "15_viz" / "report.html"),
        ("Candidate summary", run_dir / "14_rank" / "candidate_summary.csv"),
        ("Screening ranking", run_dir / "14_rank" / "rank_screening.csv"),
        ("Scientific ranking", run_dir / "14_rank" / "rank_scientific.csv"),
        ("Production ranking", run_dir / "14_rank" / "rank_production.csv"),
        ("Ranking summary", run_dir / "14_rank" / "ranking_summary.json"),
        ("Connector audit", run_dir / "13_connector-audit" / "production_connector_audit.html"),
        ("Operations report", run_dir / "16_ops" / "operations_report.html"),
    ]
    console.print(f"run_id: {manifest.run_id}")
    console.print(f"latest_stage: {manifest.stage}")
    for label, path in candidates:
        mark = "✓" if Path(path).exists() else "-"
        console.print(f"{mark} {label:22s} {path}")


@app.command("compare-backends")
def compare_backends_cmd(
    run_dir: Path = typer.Argument(...),
    out: Path = typer.Option(..., help="Output comparison directory"),
):
    paths = hpc_compare_backends(run_dir, out)
    console.print(f"[green]wrote[/green] {paths['csv']}")


if __name__ == "__main__":
    app()
