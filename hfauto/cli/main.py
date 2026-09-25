from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from hfauto.core.config import load_yaml
from hfauto.core.environment import environment_report
from hfauto.core.io import read_manifest, write_manifest
from hfauto.reporting.html_report import latest_manifest_path, render_report
from hfauto.stages.base import StageContext
from hfauto.stages.registry import get_stage
from hfauto.workflow.runner import run_pipeline

app = typer.Typer(help="Molecular complex, reaction-path, and energy workflow")
console = Console()


@app.command("doctor")
def doctor(
    config: Path | None = typer.Option(None, help="Production pipeline YAML to validate"),
    strict: bool = typer.Option(False, help="Exit non-zero when required tools are missing"),
):
    """Check Python packages and external chemistry/HPC executables."""
    report = environment_report(config)
    console.print_json(data=report)
    if strict and not report["ready"]:
        raise typer.Exit(code=1)


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
    config: Path | None = typer.Option(None, help="YAML config for this stage"),
    global_config: Path | None = typer.Option(None, help="Pipeline/global YAML to inherit temperature, pressure, and mode settings"),
    run_id: str | None = typer.Option(None, help="Override run id"),
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
    from_stage: str | None = typer.Option(None, "--from", help="Start stage name for partial run"),
    to_stage: str | None = typer.Option(None, "--to", help="End stage name for partial run"),
    start_manifest: Path | None = typer.Option(None, help="Input manifest for partial run"),
    skip_preflight: bool = typer.Option(False, help="Skip production dependency preflight (advanced use only)"),
):
    cfg = load_yaml(config)
    cfg["__config_path"] = str(config)
    if str((cfg.get("global") or {}).get("mode", "")).lower() == "production" and not skip_preflight:
        report = environment_report(cfg)
        if not report["ready"]:
            console.print_json(data=report)
            missing = ", ".join(report["missing_python"] + report["missing_external"])
            console.print(f"[red]production preflight failed[/red]: missing or unusable: {missing}")
            raise typer.Exit(code=2)
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


if __name__ == "__main__":
    app()
