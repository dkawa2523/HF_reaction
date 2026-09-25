"""Normalize backend configuration into a small scientific method record."""

from __future__ import annotations

from typing import Any

from hfauto.core.schemas.method import ElectronicStructureMethodRecord

_NON_METHOD_KEYS = {
    "allow_subprocess",
    "basin_assessment",
    "dry_run",
    "env",
    "executable",
    "fallback_to_dummy",
    "memory_mb",
    "ncores",
    "nprocs",
    "saddle_seed_hessian_evidence",
    "saddle_seed_candidate",
    "timeout_s",
}


def electronic_structure_method(
    *,
    engine: str,
    task: str,
    config: dict[str, Any],
    backend: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
) -> dict[str, Any]:
    """Return one backend-neutral method payload without execution controls."""

    known = {
        "method_id",
        "functional",
        "xc",
        "basis",
        "disp_vdw",
        "dispersion",
        "solvation_model",
        "dielectric",
        "required_program_version",
        "settings",
    }
    settings = {
        **dict(config.get("settings", {}) or {}),
        **{
        key: value
        for key, value in config.items()
        if key not in known | _NON_METHOD_KEYS and value is not None
        },
    }
    record = ElectronicStructureMethodRecord(
        engine=engine,
        backend=backend,
        task=task,
        method_id=config.get("method_id"),
        functional=config.get("functional") or config.get("xc"),
        basis=config.get("basis"),
        dispersion=config.get("dispersion", config.get("disp_vdw")),
        solvation_model=config.get("solvation_model"),
        dielectric=config.get("dielectric"),
        required_program_version=config.get("required_program_version"),
        charge=charge,
        multiplicity=multiplicity,
        settings=settings,
    )
    return record.model_dump(exclude_none=True, exclude_defaults=True)
