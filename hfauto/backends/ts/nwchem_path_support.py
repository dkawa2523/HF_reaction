"""Shared, fail-closed support for NWChem double-ended path engines.

This module validates endpoint identity and raw execution evidence.  It does
not render, run, or classify a reaction path; NEB and string engines retain
those method-specific responsibilities.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from hfauto.backends.qm.nwchem import species_xyz_path
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.geometry import kabsch_rmsd
from hfauto.chemistry.path_initialization import initialize_reaction_path
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.chemistry.xyz_trajectory import (
    read_xyz_trajectory,
    resample_xyz_trajectory,
    write_xyz_trajectory,
)
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact

_PATH_PROVENANCE_KEYS = {
    "initial_path_xyz",
    "path_initial_trajectory",
    "path_initialization",
}


def existing_file_hash(path: str | Path | None) -> str | None:
    """Hash an existing evidence file; return ``None`` for absent evidence."""

    if path is None or not Path(path).is_file():
        return None
    return sha256_file(path)


def electronic_method_config(method: dict[str, Any]) -> dict[str, Any]:
    """Exclude geometry-file provenance from electronic method identity."""

    return {
        key: value
        for key, value in method.items()
        if key not in _PATH_PROVENANCE_KEYS
    }


def resolve_path_endpoints(
    reactant: Artifact,
    product: Artifact,
    method: dict[str, Any],
) -> tuple[XYZ, XYZ, dict[str, Any]]:
    """Load compatible endpoints and resolve their common electronic state."""

    start = read_xyz(species_xyz_path(reactant))
    end = read_xyz(species_xyz_path(product))
    if start.symbols != end.symbols:
        raise ValueError("reactant/product atom order mismatch")

    state = resolve_electronic_state(reactant.data, method, start.symbols)
    end_state = resolve_electronic_state(product.data, method, end.symbols)
    state_keys = ("charge", "multiplicity", "electron_count")
    if any(state.get(key) != end_state.get(key) for key in state_keys):
        raise ValueError("reactant/product electronic state mismatch")
    return start, end, state


def validate_path_trajectory(
    images: list[XYZ],
    start: XYZ,
    end: XYZ,
    *,
    endpoint_tolerance_A: float,
) -> None:
    """Require a multi-image path with ordered atoms and matching endpoints."""

    if len(images) < 3:
        raise ValueError("reaction path requires at least three images")
    if any(image.symbols != start.symbols for image in images):
        raise ValueError("reaction path atom symbols/order mismatch")
    if kabsch_rmsd(images[0].coords, start.coords) > endpoint_tolerance_A:
        raise ValueError("reaction path does not start at the reactant")
    if kabsch_rmsd(images[-1].coords, end.coords) > endpoint_tolerance_A:
        raise ValueError("reaction path does not end at the product")


def validated_initial_path(
    source: str | Path,
    start: XYZ,
    end: XYZ,
    destination: str | Path,
    *,
    endpoint_tolerance_A: float,
    image_count: int | None = None,
) -> tuple[Path, int]:
    """Validate, optionally rediscretize, and localize an input trajectory."""

    images = read_xyz_trajectory(source)
    validate_path_trajectory(
        images,
        start,
        end,
        endpoint_tolerance_A=endpoint_tolerance_A,
    )
    if image_count is not None and len(images) != int(image_count):
        images = resample_xyz_trajectory(images, int(image_count))
    images[0] = XYZ(list(start.symbols), start.coords.copy(), images[0].comment)
    images[-1] = XYZ(list(end.symbols), end.coords.copy(), images[-1].comment)
    return write_xyz_trajectory(images, destination), len(images)


def prepare_initial_path(
    reaction: Artifact,
    start: XYZ,
    end: XYZ,
    method: dict[str, Any],
    destination: str | Path,
    *,
    image_count: int,
) -> tuple[Path | None, dict[str, Any]]:
    """Resolve a supplied path or generate one from a supported coordinate.

    No generic Cartesian fallback is created here.  Returning ``None`` makes
    the backend's documented default explicit and prevents unsupported
    internal coordinates from being mislabeled as chemistry-aware paths.
    """

    requested = method.get("path_initial_trajectory") or method.get(
        "initial_path_xyz"
    )
    if requested:
        source = Path(str(requested))
        if not source.is_file():
            raise ValueError(f"initial_path_xyz is not readable: {requested}")
        return source, {
            "strategy": "provided_xyz_path",
            "source_path": str(source),
            "source_sha256": existing_file_hash(source),
        }
    if not bool(method.get("auto_initialize_reaction_path", True)):
        return None, {"strategy": "backend_default_cartesian"}

    initialized = initialize_reaction_path(
        reaction,
        start,
        end,
        image_count=image_count,
    )
    if initialized is None:
        return None, {
            "strategy": "backend_default_cartesian",
            "reason": "no_supported_explicit_internal_coordinate",
        }
    path = write_xyz_trajectory(initialized.images, destination)
    return path, {
        **initialized.metadata,
        "generated_path": str(path),
        "generated_path_sha256": existing_file_hash(path),
    }


def assess_nwchem_rendered_method(
    method: dict[str, Any], rendered_input: str
) -> dict[str, Any]:
    """Verify the scientifically relevant NWChem method in a rendered deck."""

    functional = str(method.get("functional", method.get("xc", "pbe0")))
    basis = str(method.get("basis", "def2-svp"))
    grid = str(method.get("grid", "fine"))
    scf_tolerance = f"{float(method.get('scf_energy_tolerance', 1.0e-7)):.1e}"
    expected = {
        "functional": rf"^\s*xc\s+{re.escape(functional)}\s*$",
        "basis": rf"^\s*\*\s+library\s+{re.escape(basis)}\s*$",
        "grid": rf"^\s*grid\s+{re.escape(grid)}\s*$",
        "scf_energy_tolerance": (
            rf"^\s*convergence\s+energy\s+{re.escape(scf_tolerance)}\s*$"
        ),
        "charge": rf"^\s*charge\s+{int(method.get('charge', 0))}\s*$",
        "multiplicity": (
            rf"^\s*mult\s+{int(method.get('multiplicity', 1))}\s*$"
        ),
    }
    optional_values: dict[str, tuple[str, Any]] = {
        "disp_vdw": ("disp\\s+vdw", method.get("disp_vdw")),
        "path_image_count": ("nbeads", method.get("path_image_count")),
    }
    backend = str(method.get("backend") or "").lower()
    if "string" in backend:
        optional_values.update(
            {
                "string_maxiter": ("maxiter", method.get("string_maxiter")),
                "string_stepsize": ("stepsize", method.get("string_stepsize")),
                "string_history": ("nhist", method.get("string_nhist")),
                "string_interpolation": (
                    "interpol",
                    method.get("string_interpol"),
                ),
                "string_tolerance": ("tol", method.get("string_tolerance")),
            }
        )
    elif "neb" in backend:
        optional_values.update(
            {
                "neb_maxiter": ("maxiter", method.get("neb_maxiter")),
                "neb_stepsize": ("stepsize", method.get("neb_stepsize")),
            }
        )
    expected.update(
        {
            key: rf"^\s*{keyword}\s+{re.escape(str(value))}\s*$"
            for key, (keyword, value) in optional_values.items()
            if value is not None
        }
    )
    checks = {
        key: bool(
            re.search(pattern, rendered_input, flags=re.IGNORECASE | re.MULTILINE)
        )
        for key, pattern in expected.items()
    }
    return {"accepted": all(checks.values()), "checks": checks}


def assess_nwchem_path_evidence(
    parsed_output: dict[str, Any],
    method: dict[str, Any],
    *,
    dispersion_in_input: bool,
    rendered_input: str | None = None,
) -> dict[str, Any]:
    """Assess method lineage from the rendered input and raw NWChem output."""

    observed_version = parsed_output.get("program_version")
    required_version = method.get("required_program_version")
    program_version_ok = bool(
        observed_version is not None
        and (
            required_version is None
            or str(observed_version) == str(required_version)
        )
    )
    dispersion_requested = method.get("disp_vdw") is not None
    dispersion_evidence_ok = bool(
        not dispersion_requested
        or (
            dispersion_in_input
            and parsed_output.get("dft_d3_applied") is True
        )
    )
    input_checks: dict[str, bool] = {}
    if rendered_input is not None:
        input_checks = assess_nwchem_rendered_method(
            method, rendered_input
        )["checks"]
    input_method_evidence_ok = bool(
        rendered_input is None or all(input_checks.values())
    )
    method_identity_validated = bool(
        program_version_ok
        and dispersion_evidence_ok
        and input_method_evidence_ok
    )
    normal_termination = parsed_output.get("normal_termination") is True
    validated = bool(normal_termination and method_identity_validated)
    return {
        "validated": validated,
        "method_identity_validated": method_identity_validated,
        "normal_termination_validated": normal_termination,
        "program_version_ok": program_version_ok,
        "dispersion_evidence_ok": dispersion_evidence_ok,
        "input_method_evidence_ok": input_method_evidence_ok,
        "input_method_checks": input_checks,
        "observed_method": {
            "program_version": observed_version,
            "disp_vdw": (
                method.get("disp_vdw")
                if parsed_output.get("dft_d3_applied") is True
                else None
            ),
            "functional": method.get("functional", method.get("xc", "pbe0")),
            "basis": method.get("basis", "def2-svp"),
            "grid": method.get("grid", "fine"),
            "scf_energy_tolerance": float(
                method.get("scf_energy_tolerance", 1.0e-7)
            ),
        },
    }
