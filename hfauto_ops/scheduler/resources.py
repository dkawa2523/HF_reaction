from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


DEFAULT_RESOURCE_PROFILES: dict[str, dict[str, Any]] = {
    "conformers": {"backend": "crest", "threads": 8, "memory_gb": 8, "walltime": "04:00:00", "queue": "short"},
    "preopt": {"backend": "xtb", "threads": 4, "memory_gb": 4, "walltime": "01:00:00", "queue": "short"},
    "dft-minima": {"backend": "orca", "threads": 16, "memory_gb": 32, "walltime": "12:00:00", "queue": "medium"},
    "sp": {"backend": "orca", "threads": 24, "memory_gb": 48, "walltime": "08:00:00", "queue": "medium"},
    "ts-search": {"backend": "orca_nebts", "threads": 24, "memory_gb": 64, "walltime": "24:00:00", "queue": "long"},
    "irc": {"backend": "orca", "threads": 16, "memory_gb": 32, "walltime": "12:00:00", "queue": "medium"},
    "thermo": {"backend": "goodvibes", "threads": 1, "memory_gb": 2, "walltime": "00:30:00", "queue": "short"},
    "kinetics": {"backend": "arkane_cantera", "threads": 2, "memory_gb": 4, "walltime": "01:00:00", "queue": "short"},
    "calibrate": {"backend": "publicdb", "threads": 1, "memory_gb": 2, "walltime": "01:00:00", "queue": "short"},
    "viz": {"backend": "hfauto_viz", "threads": 1, "memory_gb": 4, "walltime": "00:30:00", "queue": "short"},
}

STAGE_INPUT_TYPES: dict[str, list[str]] = {
    "conformers": ["molecule", "molecule_enriched"],
    "preopt": ["species"],
    "dft-minima": ["species_preopt", "species"],
    "sp": ["species_optimized", "species_preopt", "species"],
    "ts-search": ["reaction"],
    "irc": ["reaction_validated"],
    "thermo": ["calculation", "reaction_path_validated"],
    "kinetics": ["thermo"],
    "calibrate": ["molecule_enriched", "molecule", "thermo"],
    "viz": ["ranking", "table", "thermo", "kinetics"],
}

STAGE_EXPECTED_OUTPUT_TYPES: dict[str, list[str]] = {
    "conformers": ["conformer"],
    "preopt": ["species_preopt", "preopt_geometry"],
    "dft-minima": ["species_optimized", "calculation"],
    "sp": ["calculation"],
    "ts-search": ["ts_path", "reaction_validated"],
    "irc": ["irc", "reaction_path_validated"],
    "thermo": ["thermo"],
    "kinetics": ["kinetics"],
    "calibrate": ["calibration", "method_validation"],
    "viz": ["visualization_bundle"],
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
                f"hfauto pipeline --config configs/pipelines/phase9_hpc_template.yaml "
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
