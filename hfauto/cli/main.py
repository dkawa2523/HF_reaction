"""hfauto command line (design §7.4): doctor, run, status, report and case.

PIPELINE, --system and --site take a file path or a name under ``configs/<kind>/`` of the
working directory.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import NoReturn

import typer
from rich.console import Console
from rich.markup import escape
from rich.table import Table

from hfauto.chemistry.xyz import read_xyz
from hfauto.core.evidence import FailureKind
from hfauto.core.ids import path_token
from hfauto.core.system import SystemConfig
from hfauto.pipeline import config, preflight, runner
from hfauto.pipeline.layout import RunLayout
from hfauto.reporting import html

app = typer.Typer(
    help="Reaction discovery for gas-phase molecules and complexes.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()
CONFIGS = Path("configs")
RUNS = Path("runs")


def _fail(message: object) -> NoReturn:
    console.print(f"[red]error[/red] {escape(str(message))}")
    raise typer.Exit(2)


def config_path(kind: str, value: str) -> Path:
    """``value`` when it names a .yaml file, else ``configs/<kind>/<value>.yaml``."""
    path = Path(value)
    if path.suffix not in (".yaml", ".yml"):
        path = CONFIGS / kind / f"{value}.yaml"
    if not path.is_file():
        _fail(f"{kind} file not found: {path}")
    return path


def _layout(run_dir: Path) -> RunLayout:
    layout = RunLayout(run_dir)
    if not layout.run_state_path.is_file():
        _fail(f"{run_dir} is not a run directory (no run_state.json)")
    return layout


def _failure_kinds(value: str) -> tuple[FailureKind, ...]:
    try:
        return tuple(FailureKind(kind.strip()) for kind in value.split(",") if kind.strip())
    except ValueError:
        _fail(f"--retry-failed takes a comma-separated list of {', '.join(FailureKind)}")


def _site_probe(site: config.SiteConfig) -> config.ResolvedConfig:
    """A pipeline naming every engine of the site, so that preflight checks all of them."""
    stages = [config.StageEntry(id=name, stage="doctor", engine=name) for name in site.engines]
    return config.ResolvedConfig(
        pipeline=config.PipelineConfig(pipeline_id="doctor", stages=stages),
        system=SystemConfig(system_id="doctor", species=[]),
        site=site,
        methods={},
        code_version=config.code_version(),
    )


def _failed_stages(layout: RunLayout) -> list[str]:
    """Stages that failed or saved an artifact carrying a failure."""
    return [s.stage_id for s in layout.read_state()
            if s.status == "failed" or (s.status == "done" and s.n_failed > 0)]


def _print_problems(problems: list[str]) -> None:
    for problem in problems:
        console.print(f"[red]problem[/red] {escape(problem)}")


@app.command()
def doctor(
    site: str | None = typer.Option(None, help="Site name or file (default: configs/sites/*)"),
) -> None:
    """Preflight every engine of a site: executables, version pins, worker modules, scratch."""
    paths = [config_path("sites", site)] if site else sorted((CONFIGS / "sites").glob("*.yaml"))
    if not paths:
        _fail(f"no site files under {CONFIGS / 'sites'}")
    console.print(f"hfauto {config.code_version()}")
    ready = True
    for path in paths:
        probe = _site_probe(config.load_site(path))
        problems = preflight.preflight(probe)
        table = Table("engine", "version pin", "check", title=f"site {probe.site.site} ({path})")
        for name, engine in probe.site.engines.items():
            own = [p.split(": ", 1)[1] for p in problems if p.startswith(f"{name}: ")]
            table.add_row(name, engine.version, escape("\n".join(own)) or "[green]ok[/green]")
        console.print(table)
        _print_problems([p for p in problems if p.split(": ", 1)[0] not in probe.site.engines])
        ready = ready and not problems
    if not ready:
        raise typer.Exit(1)


@app.command()
def run(
    pipeline: str = typer.Argument(..., help="Pipeline name (configs/pipelines/<name>.yaml) or file"),
    system: str = typer.Option(..., help="System name or file"),
    site: str = typer.Option(..., help="Site name or file"),
    run_dir: Path | None = typer.Option(None, help="Default: runs/<system_id>_<pipeline_id>"),
    start: str | None = typer.Option(None, "--from", help="Re-run from this stage id"),
    stop: str | None = typer.Option(None, "--to", help="Stop after this stage id"),
    dry_run: bool = typer.Option(False, help="Check the configs and print the plan only"),
    retry_failed: str = typer.Option("", help="FailureKinds whose cached failures are retried"),
) -> None:
    """Run a pipeline for a system on a site (resumes a run directory).

    Exit code 1 when a stage failed or a saved artifact carries a failure."""
    kinds = _failure_kinds(retry_failed)
    paths = [config_path(k, v) for k, v in (("pipelines", pipeline), ("systems", system),
                                            ("sites", site))]
    try:
        resolved = config.load(*paths)
        pipeline_id = resolved.pipeline.pipeline_id
        target = run_dir or RUNS / f"{resolved.system.system_id}_{pipeline_id}"
        steps = runner.plan(resolved, RunLayout(target), start=start, stop=stop, retry_failed=kinds)
    except (OSError, ValueError, KeyError) as exc:
        _fail(exc)
    problems = preflight.preflight(resolved, dry_run=dry_run)
    if problems:
        _print_problems(problems)
        raise typer.Exit(2)
    if dry_run:
        table = Table("stage id", "stage", "action")
        for step in steps:
            table.add_row(step.stage_id, step.stage, step.action)
        console.print(f"dry run of {pipeline_id} in {escape(str(target))}", soft_wrap=True)
        console.print(table)
        return
    runner.run_pipeline(resolved, target, start=start, stop=stop, retry_failed=kinds)
    failed = _failed_stages(RunLayout(target))
    if failed:
        console.print(f"[yellow]done with failures[/yellow] in {escape(', '.join(failed))}:"
                      f" {escape(str(target))}", soft_wrap=True)
        raise typer.Exit(1)
    console.print(f"[green]done[/green] {escape(str(target))}", soft_wrap=True)


@app.command()
def status(run_dir: Path) -> None:
    """Stages in execution order, job failures by FailureKind and the JobStore reuse rate."""
    states = _layout(run_dir).read_state()
    table = Table("stage", "pipeline", "status", "ok", "failed", "job hits", "job misses")
    failures: Counter[str] = Counter()
    for s in states:
        jobs = s.jobs
        table.add_row(s.stage_id, s.pipeline_id, s.status, str(s.n_ok), str(s.n_failed),
                      str(jobs.hits), str(jobs.misses))
        failures.update(jobs.failures_by_kind)
    console.print(table)
    counts = ", ".join(f"{kind}={n}" for kind, n in sorted(failures.items()))
    console.print(f"job failures: {counts or 'none'}")
    hits = sum(s.jobs.hits for s in states)
    total = hits + sum(s.jobs.misses for s in states)
    console.print(f"job reuse: {hits}/{total}" + (f" ({hits / total:.0%})" if total else ""))


@app.command()
def report(run_dir: Path) -> None:
    """Write RUN_DIR/report.html from the view of every done stage."""
    layout = _layout(run_dir)
    path = html.render(layout.view(), layout.run_dir,
                       load_xyz=lambda geometry: read_xyz(layout.run_dir / geometry.file.path))
    console.print(f"[green]wrote[/green] {escape(str(path))}", soft_wrap=True)


@app.command()
def case(run_dir: Path, reaction_id: str) -> None:
    """Print cases/<REACTION_ID>/log.jsonl of every stage that drove the reaction."""
    layout = _layout(run_dir)
    logs = [layout.cases_dir(s.stage_id) / path_token(reaction_id) / "log.jsonl"
            for s in layout.read_state()]
    found = [path for path in logs if path.is_file()]
    if not found:
        _fail(f"no case log for {reaction_id!r} in {run_dir}")
    for path in found:
        console.rule(escape(path.relative_to(layout.run_dir).as_posix()))
        typer.echo(path.read_text(encoding="utf-8"), nl=False)


if __name__ == "__main__":
    app()
