"""Checks before a run (design §7.4): executables, version pins, worker modules, method
support and the scratch location. Problems are returned as messages, never raised.
"""

from __future__ import annotations

import re
import sys
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from hfauto.backends.protocols import Capability, Requirements
from hfauto.core.method import EngineSite
from hfauto.execution.process import Command, resolve_executable, run_command
from hfauto.pipeline.config import PipelineConfig, ResolvedConfig, SiteConfig, method_ids

_CAPABILITIES = frozenset(c.value for c in Capability)
_DRIVE_MOUNT = re.compile(r"^/mnt/[A-Za-z](/|$)")  # Windows drives seen from WSL (drvfs)
_COMMAND_TIMEOUT_S = 60.0
_FIND_MODULES = (
    "import importlib.util, sys; "
    "print(' '.join(m for m in sys.argv[1:] if importlib.util.find_spec(m) is None))"
)


@dataclass(frozen=True)
class EngineUse:
    capability: Capability | None  # None for a bare ``engine:`` key (found in the registry)
    name: str
    method: str | None


def _named_engines(node: dict) -> list[tuple[Capability | None, str]]:
    named: list[tuple[Capability | None, str]] = []
    if isinstance(node.get("engine"), str):
        named.append((None, node["engine"]))
    engines = node.get("engines")
    for mapping in (node, engines if isinstance(engines, dict) else {}):
        named += [
            (Capability(key), value)
            for key, value in mapping.items()
            if key in _CAPABILITIES and isinstance(value, str)
        ]
    return named


def _uses(node: object) -> list[EngineUse]:
    if isinstance(node, list):
        return [use for item in node for use in _uses(item)]
    if not isinstance(node, dict):
        return []
    own = {k: v for k, v in node.items() if k in ("method", "methods")}
    methods: list[str | None] = [*sorted(method_ids(own))] or [None]
    uses = [EngineUse(cap, name, m) for cap, name in _named_engines(node) for m in methods]
    for key, value in node.items():
        if key != "engines":
            uses += _uses(value)
    return uses


def engine_uses(pipeline: PipelineConfig) -> list[EngineUse]:
    """(capability, engine, method) triples referenced by the pipeline's stage settings.

    ``engine:`` / ``engines: {<capability>: name}`` / ``<capability>: name`` keys pair with
    the ``method:`` / ``methods:`` keys of the same mapping.
    """
    uses = _uses([entry.settings() for entry in pipeline.stages])
    return list(dict.fromkeys(uses))


def _run(tmp: Path, argv: tuple[str, ...]) -> tuple[int | None, str]:
    try:
        result = run_command(Command(argv=argv, cwd=tmp), timeout_s=_COMMAND_TIMEOUT_S)
    except OSError as exc:
        return None, str(exc)
    output = result.stdout.read_text(encoding="utf-8", errors="replace")
    output += result.stderr.read_text(encoding="utf-8", errors="replace")
    return result.returncode, output


def _executable_problems(name: str, req: Requirements, site: EngineSite) -> list[str]:
    return [
        f"{name}: executable {key!r} not found ({site.executables.get(key) or 'PATH'})"
        for key in req.executables
        if resolve_executable(key, site.executables.get(key)) is None
    ]


def _version_problems(tmp: Path, name: str, req: Requirements, site: EngineSite) -> list[str]:
    if not req.version_command:
        return []
    head, *rest = req.version_command
    executable = resolve_executable(head, site.executables.get(head))
    if executable is None:
        return [f"{name}: version command {head!r} not found"]
    _, output = _run(tmp, (executable, *rest))
    if re.search(rf"(?<![\w.]){re.escape(site.version)}(?![\w.])", output) is None:
        return [f"{name}: version pin {site.version!r} not in the output of {req.version_command}"]
    return []


def _module_problems(tmp: Path, name: str, req: Requirements, site: EngineSite) -> list[str]:
    if not req.python_modules:
        return []
    python = site.python or sys.executable
    code, output = _run(tmp, (python, "-c", _FIND_MODULES, *req.python_modules))
    if code != 0:
        return [f"{name}: worker python {python!r} failed: {output.strip()[-200:]}"]
    return [f"{name}: python module {m!r} missing in {python}" for m in output.split()]


def scratch_problems(site: SiteConfig) -> list[str]:
    paths = [site.scratch_root] + [e.scratch_dir for e in site.engines.values() if e.scratch_dir]
    return [
        f"scratch {path} is on a Windows drive mount; use an ext4 path"
        for path in paths
        if _DRIVE_MOUNT.match(Path(path).as_posix())
    ]


def check_site(
    site: SiteConfig, requirements: dict[str, Requirements], *, dry_run: bool = False
) -> list[str]:
    """Scratch location, and unless ``dry_run`` executables, version pins and worker modules."""
    problems = scratch_problems(site)
    missing = sorted(set(requirements) - set(site.engines))
    problems += [f"{name}: not configured in site {site.site!r}" for name in missing]
    if dry_run:
        return problems
    with tempfile.TemporaryDirectory(prefix="hfauto_preflight_") as tmp:
        for name, req in sorted(requirements.items()):
            if name in site.engines:
                engine_site = site.engines[name]
                problems += _executable_problems(name, req, engine_site)
                problems += _version_problems(Path(tmp), name, req, engine_site)
                problems += _module_problems(Path(tmp), name, req, engine_site)
    return problems


Use = tuple[Capability, str, str | None]  # registered (capability, engine, method)


def _registered(use: EngineUse) -> Capability | None:
    """The capability under which ``use.name`` is in the registry table, if any."""
    from hfauto.backends.engines import _TABLE

    candidates = [use.capability] if use.capability is not None else list(_TABLE)
    return next((c for c in candidates if use.name in _TABLE[c]), None)


def _support_problems(resolved: ResolvedConfig, uses: Iterable[Use]) -> list[str]:
    from hfauto.backends.engines import create
    from hfauto.execution.jobs import JobRunner
    from hfauto.execution.jobstore import JobStore

    problems: list[str] = []
    with tempfile.TemporaryDirectory(prefix="hfauto_preflight_") as tmp:
        jobs = JobRunner(JobStore(Path(tmp) / "jobs"), cores=resolved.site.cores)
        for capability, name, method in uses:
            if method is None or name not in resolved.site.engines:
                continue
            try:
                engine = create(capability, name, jobs=jobs, site=resolved.site.engines[name])
            except Exception as exc:  # an engine that cannot even be built is a problem
                problems.append(f"{name}: cannot create engine: {exc}")
                continue
            if not engine.supports(resolved.methods[method]):
                problems.append(f"{name}: does not support method {method!r}")
    return problems


def preflight(resolved: ResolvedConfig, *, dry_run: bool = False) -> list[str]:
    """Every problem found for the engines the pipeline uses; empty means ready to run.

    ``dry_run`` skips everything that runs an executable (executables, version pins and
    the worker python's modules).
    """
    from hfauto.backends.engines import requirements_for

    uses: list[Use] = []
    problems: list[str] = []
    for use in engine_uses(resolved.pipeline):
        capability = _registered(use)
        if capability is None:
            problems.append(f"{use.name}: not a registered {use.capability or 'engine'}")
        else:
            uses.append((capability, use.name, use.method))
    selection: dict[Capability, list[str]] = {}
    for capability, name, _ in uses:
        selection.setdefault(capability, []).append(name)
    problems += check_site(resolved.site, requirements_for(selection), dry_run=dry_run)
    return problems + _support_problems(resolved, dict.fromkeys(uses))
