from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.artifact_types import STAGE_CONTRACTS
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest

DEFAULT_RESOURCE_PROFILES: dict[str, dict[str, Any]] = {
    "conformers": {
        "backend": "crest",
        "threads": 8,
        "memory_gb": 8,
        "walltime": "04:00:00",
        "queue": "short",
    },
    "build-complexes": {
        "backend": "crest",
        "threads": 8,
        "memory_gb": 8,
        "walltime": "04:00:00",
        "queue": "short",
    },
    "explore-reactions": {
        "backend": "readuct",
        "threads": 8,
        "memory_gb": 16,
        "walltime": "12:00:00",
        "queue": "medium",
    },
    "preopt": {
        "backend": "xtb",
        "threads": 4,
        "memory_gb": 4,
        "walltime": "01:00:00",
        "queue": "short",
    },
    "relaxation-discovery": {
        "backend": "internal-connectivity",
        "threads": 1,
        "memory_gb": 2,
        "walltime": "00:30:00",
        "queue": "short",
    },
    "dft-minima": {
        "backend": "configured-qm",
        "threads": 16,
        "memory_gb": 32,
        "walltime": "12:00:00",
        "queue": "medium",
    },
    "sp": {
        "backend": "configured-qm",
        "threads": 24,
        "memory_gb": 48,
        "walltime": "08:00:00",
        "queue": "medium",
    },
    "ts-search": {
        "backend": "configured-path",
        "threads": 24,
        "memory_gb": 64,
        "walltime": "24:00:00",
        "queue": "long",
    },
    "irc": {
        "backend": "configured-irc",
        "threads": 16,
        "memory_gb": 32,
        "walltime": "12:00:00",
        "queue": "medium",
    },
    "thermo": {
        "backend": "goodvibes",
        "threads": 1,
        "memory_gb": 2,
        "walltime": "00:30:00",
        "queue": "short",
    },
}

STAGE_INPUT_TYPES: dict[str, list[str]] = {
    name: list(contract["in"]) for name, contract in STAGE_CONTRACTS.items()
}

STAGE_EXPECTED_OUTPUT_TYPES: dict[str, list[str]] = {
    name: list(contract["out"]) for name, contract in STAGE_CONTRACTS.items()
}


@dataclass
class JobSpec:
    job_id: str
    stage: str
    array_index: int
    input_artifact_id: str
    input_artifact_type: str
    expected_outputs: str
    command: str
    threads: int
    memory_gb: float
    walltime: str
    queue: str
    backend: str
    priority: int
    fingerprint: str
    cache_status: str = "unknown"
    retry_of: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _profile(stage: str, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    prof = dict(DEFAULT_RESOURCE_PROFILES.get(stage, DEFAULT_RESOURCE_PROFILES["preopt"]))
    if overrides:
        prof.update(overrides)
    return prof


def _representative_inputs(manifest: Manifest, stage: str, limit: int | None = None) -> list[Artifact]:
    types = STAGE_INPUT_TYPES.get(stage, [])
    seen: set[str] = set()
    selected: list[Artifact] = []
    for typ in types:
        for art in manifest.latest_artifacts(typ):
            key = art.data.get("species_id") or art.data.get("reaction_id") or art.data.get("mol_id") or art.artifact_id
            key = str(key)
            if key in seen or art.status.status == "failed":
                continue
            seen.add(key)
            selected.append(art)
            if limit and len(selected) >= limit:
                return selected
    return selected


def make_stage_job_plan(
    manifest: Manifest,
    run_dir: str | Path,
    stages: list[str],
    config: dict[str, Any] | None = None,
) -> list[JobSpec]:
    """Create a scheduler-neutral job plan from current artifacts.

    The plan intentionally does not launch calculations. It gives the HPC layer
    a machine-readable description of candidate jobs, expected resources, and
    cache fingerprints.
    """
    cfg = config or {}
    profiles = cfg.get("resource_profiles", {}) or {}
    pipeline_config = str(
        cfg.get("pipeline_config") or "$HFAUTO_PIPELINE_CONFIG"
    )
    max_jobs_per_stage = cfg.get("max_jobs_per_stage")
    run_dir = Path(run_dir)
    jobs: list[JobSpec] = []
    priority_base = int(cfg.get("priority_base", 100))
    for stage_idx, stage in enumerate(stages):
        inputs = _representative_inputs(manifest, stage, limit=max_jobs_per_stage)
        prof = _profile(stage, profiles.get(stage))
        for i, art in enumerate(inputs):
            payload = {
                "run_id": manifest.run_id,
                "stage": stage,
                "input_artifact_id": art.artifact_id,
                "input_artifact_type": art.artifact_type,
                "input_data_key": art.data.get("species_id") or art.data.get("reaction_id") or art.data.get("mol_id"),
                "backend": prof.get("backend"),
            }
            fp = fingerprint_dict(payload)
            job_id = f"job_{stage.replace('-', '_')}_{fp}"
            command = (
                f"hfauto pipeline --config {pipeline_config} "
                f"--run-id {manifest.run_id} --from {stage} --to {stage} "
                f"--start-manifest $(cat {run_dir / 'manifest.path'})"
            )
            jobs.append(
                JobSpec(
                    job_id=job_id,
                    stage=stage,
                    array_index=i,
                    input_artifact_id=art.artifact_id,
                    input_artifact_type=art.artifact_type,
                    expected_outputs=",".join(STAGE_EXPECTED_OUTPUT_TYPES.get(stage, [])),
                    command=command,
                    threads=int(prof.get("threads", 1)),
                    memory_gb=float(prof.get("memory_gb", 2)),
                    walltime=str(prof.get("walltime", "01:00:00")),
                    queue=str(prof.get("queue", "short")),
                    backend=str(prof.get("backend", "unknown")),
                    priority=priority_base - stage_idx,
                    fingerprint=fp,
                    notes=str(prof.get("notes", "")),
                )
            )
    return jobs


def jobs_to_dataframe(jobs: list[JobSpec]) -> pd.DataFrame:
    return pd.DataFrame([j.to_dict() for j in jobs])
