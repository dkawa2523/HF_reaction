from __future__ import annotations

"""Resource estimation for independent workflow stages and backend jobs.

The estimates are intentionally conservative defaults.  They are not used to
change scientific results; they only help create reviewable job plans and HPC
submit templates.
"""

from dataclasses import dataclass, asdict
from typing import Any


@dataclass(frozen=True)
class ResourceRequest:
    stage: str
    backend: str = "local"
    ncores: int = 1
    memory_gb: float = 2.0
    walltime: str = "00:30:00"
    partition: str | None = None
    queue: str | None = None
    gpus: int = 0
    priority: int = 50
    reason: str = "default"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


DEFAULT_STAGE_RESOURCES: dict[str, ResourceRequest] = {
    "ingest": ResourceRequest("ingest", ncores=1, memory_gb=1, walltime="00:10:00", priority=10, reason="SDF parsing"),
    "enrich": ResourceRequest("enrich", ncores=1, memory_gb=1, walltime="00:30:00", priority=20, reason="cache-first DB enrichment"),
    "detect-sites": ResourceRequest("detect-sites", ncores=1, memory_gb=1, walltime="00:10:00", priority=15, reason="SMARTS matching"),
    "conformers": ResourceRequest("conformers", backend="crest_or_rdkit", ncores=4, memory_gb=8, walltime="04:00:00", priority=50, reason="conformer ensemble generation"),
    "build-hf": ResourceRequest("build-hf", ncores=1, memory_gb=2, walltime="00:30:00", priority=30, reason="HF endpoint construction"),
    "preopt": ResourceRequest("preopt", backend="xtb", ncores=4, memory_gb=8, walltime="02:00:00", priority=55, reason="xTB preoptimization"),
    "dft-minima": ResourceRequest("dft-minima", backend="orca", ncores=16, memory_gb=32, walltime="24:00:00", priority=70, reason="DFT opt/freq"),
    "ts-search": ResourceRequest("ts-search", backend="orca_nebts", ncores=32, memory_gb=64, walltime="48:00:00", priority=85, reason="NEB-TS / OptTS"),
    "irc": ResourceRequest("irc", backend="orca", ncores=16, memory_gb=32, walltime="24:00:00", priority=80, reason="IRC validation"),
    "sp": ResourceRequest("sp", backend="orca", ncores=16, memory_gb=32, walltime="12:00:00", priority=75, reason="high-level single point"),
    "thermo": ResourceRequest("thermo", backend="goodvibes_internal", ncores=1, memory_gb=2, walltime="00:30:00", priority=35, reason="thermal corrections"),
    "kinetics": ResourceRequest("kinetics", backend="tst_cantera", ncores=1, memory_gb=2, walltime="00:30:00", priority=35, reason="TST / mechanism export"),
    "calibrate": ResourceRequest("calibrate", ncores=1, memory_gb=2, walltime="00:30:00", priority=25, reason="public reference audit"),
    "rank": ResourceRequest("rank", ncores=1, memory_gb=2, walltime="00:10:00", priority=20, reason="scoring and tables"),
    "viz": ResourceRequest("viz", ncores=1, memory_gb=4, walltime="00:30:00", priority=20, reason="HTML/report rendering"),
}

STAGE_ALIASES = {
    "detect_sites": "detect-sites",
    "build_hf": "build-hf",
    "dft_minima": "dft-minima",
    "ts_search": "ts-search",
}


def canonical_stage_name(stage: str) -> str:
    return STAGE_ALIASES.get(stage, stage)


def parse_walltime_to_seconds(walltime: str) -> int:
    parts = [int(p) for p in walltime.split(":")]
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h, m, s = 0, parts[0], parts[1]
    else:
        h, m, s = 0, 0, parts[0]
    return h * 3600 + m * 60 + s


def merge_resource(base: ResourceRequest, override: dict[str, Any] | None = None) -> ResourceRequest:
    if not override:
        return base
    data = base.to_dict()
    data.update({k: v for k, v in override.items() if v is not None})
    return ResourceRequest(**data)


def estimate_stage_resources(stage: str, config: dict[str, Any] | None = None) -> ResourceRequest:
    name = canonical_stage_name(stage)
    base = DEFAULT_STAGE_RESOURCES.get(name, ResourceRequest(name, reason="generic stage"))
    cfg = config or {}
    stage_overrides = (cfg.get("resources") or {}).get(name) or cfg.get("resources") or {}
    # If this is a per-stage config with a nested resources block, use it.
    if any(k in stage_overrides for k in ["ncores", "memory_gb", "walltime", "partition", "queue", "gpus"]):
        return merge_resource(base, stage_overrides)
    return base


def estimate_artifact_resources(artifact: Any, stage: str | None = None, config: dict[str, Any] | None = None) -> ResourceRequest:
    """Return a resource estimate using artifact metadata when possible."""
    data = getattr(artifact, "data", {}) or {}
    method = getattr(artifact, "method", None) or {}
    qc = getattr(artifact, "qc", {}) or {}
    inferred_stage = stage or method.get("stage") or method.get("task") or data.get("stage") or "unknown"
    req = estimate_stage_resources(str(inferred_stage), config)
    # Heuristic: larger HF clusters and TS jobs are a little heavier.
    try:
        hf_n = int(data.get("hf_n") or 0)
    except Exception:
        hf_n = 0
    if hf_n >= 3 and req.memory_gb < 16:
        req = merge_resource(req, {"memory_gb": max(req.memory_gb, 12), "reason": req.reason + "; hf_n>=3"})
    if qc.get("retry_recommended") or getattr(getattr(artifact, "status", None), "status", None) == "failed":
        req = merge_resource(req, {"walltime": "08:00:00" if parse_walltime_to_seconds(req.walltime) < 8 * 3600 else req.walltime, "reason": req.reason + "; retry"})
    return req
