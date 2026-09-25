"""Small fail-closed helpers for keeping endpoint, TS, and IRC on one PES."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.nwchem_evidence import exact_int
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact

METHOD_KEYS = ("engine", "functional", "basis", "disp_vdw", "program_version")
PATH_NUMERICAL_KEYS = ("grid", "scf_energy_tolerance")
ENDPOINT_SOURCE_KEYS = ("input", "output", "stderr", "final_xyz")


def _lower_text(value: Any) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    return value.strip().lower()


def evaluate_endpoint_source_integrity(endpoint: Artifact) -> dict[str, Any]:
    """Rehash the raw files previously bound to an optimized minimum."""

    evidence = (endpoint.provenance or {}).get(
        "validated_source_integrity"
    )
    reasons: list[str] = []
    if not isinstance(evidence, dict):
        evidence = {}
        reasons.append("validated_source_integrity_missing")
    if evidence.get("accepted") is not True:
        reasons.append("validated_source_integrity_not_accepted")
    files = evidence.get("files")
    if not isinstance(files, dict):
        files = {}
    for key in ENDPOINT_SOURCE_KEYS:
        item = files.get(key)
        if not isinstance(item, dict):
            reasons.append(f"validated_source_{key}_missing")
            continue
        path = item.get("path")
        expected_hash = item.get("sha256")
        try:
            observed_hash = (
                sha256_file(path) if path and Path(path).is_file() else None
            )
        except OSError:
            observed_hash = None
        if not expected_hash or observed_hash != expected_hash:
            reasons.append(f"validated_source_{key}_hash_mismatch")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "files": files,
        "endpoint_artifact_id": endpoint.artifact_id,
    }


def evaluate_endpoint_pair_lineage(
    reactant: Artifact,
    product: Artifact,
) -> dict[str, Any]:
    """Require two stationary endpoints to be hash-bound points on one PES."""

    endpoint_methods: dict[str, dict[str, Any]] = {}
    endpoint_states: dict[str, dict[str, int | None]] = {}
    endpoint_numerical: dict[str, dict[str, Any]] = {}
    source_integrity: dict[str, dict[str, Any]] = {}
    reasons: list[str] = []
    for role, endpoint in (("reactant", reactant), ("product", product)):
        method, method_reasons = _endpoint_method(endpoint)
        state, state_reasons = _endpoint_state(endpoint)
        numerical = (endpoint.provenance or {}).get(
            "validated_numerical_settings"
        )
        if not isinstance(numerical, dict) or any(
            numerical.get(key) is None for key in PATH_NUMERICAL_KEYS
        ):
            numerical = {}
            reasons.append(f"{role}_validated_numerical_settings_missing")
        if endpoint.qc.get("numerical_method_evidence_validated") is not True:
            reasons.append(
                f"{role}_numerical_method_evidence_validated_is_not_true"
            )
        source = evaluate_endpoint_source_integrity(endpoint)
        endpoint_methods[role] = method
        endpoint_states[role] = state
        observed_numerical = path_numerical_settings(endpoint.method or {})
        if any(
            observed_numerical.get(key) != numerical.get(key)
            for key in PATH_NUMERICAL_KEYS
        ):
            reasons.append(
                f"{role}_numerical_settings_do_not_match_validated_provenance"
            )
        endpoint_numerical[role] = numerical
        source_integrity[role] = source
        reasons.extend(f"{role}_{reason}" for reason in method_reasons)
        reasons.extend(f"{role}_{reason}" for reason in state_reasons)
        reasons.extend(
            f"{role}_{reason}" for reason in source["reasons"]
        )
    if endpoint_methods["reactant"] != endpoint_methods["product"]:
        reasons.append("endpoint_methods_do_not_match")
    if endpoint_states["reactant"] != endpoint_states["product"]:
        reasons.append("endpoint_electronic_states_do_not_match")
    if endpoint_numerical["reactant"] != endpoint_numerical["product"]:
        reasons.append("endpoint_numerical_settings_do_not_match")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "endpoint_methods": endpoint_methods,
        "endpoint_electronic_states": endpoint_states,
        "endpoint_numerical_settings": endpoint_numerical,
        "source_integrity": source_integrity,
    }


def configured_qm_method(method: dict[str, Any], backend_name: str | None) -> dict[str, Any]:
    """Return the calculator method hidden behind a path-search wrapper."""

    program = _lower_text(method.get("program") or method.get("calculator_program"))
    backend = _lower_text(backend_name)
    if program is None and backend is not None:
        if "nwchem" in backend:
            program = "nwchem"
        elif "orca" in backend:
            program = "orca"
    return {
        "engine": program,
        "functional": _lower_text(method.get("functional") or method.get("xc")),
        "basis": _lower_text(method.get("basis")),
        "disp_vdw": exact_int(method.get("disp_vdw")),
        "program_version": _lower_text(method.get("required_program_version")),
    }


def path_numerical_settings(method: dict[str, Any]) -> dict[str, Any]:
    """Return numerical settings that matter for sub-kcal path comparison."""

    settings = {**method, **(method.get("settings") or {})}
    engine = _lower_text(settings.get("engine") or settings.get("program"))
    if engine is None:
        backend = _lower_text(settings.get("backend")) or ""
        engine = "nwchem" if "nwchem" in backend else None
    return {
        "grid": _lower_text(settings.get("grid") or ("fine" if engine == "nwchem" else None)),
        "scf_energy_tolerance": (
            float(settings.get("scf_energy_tolerance", 1.0e-7))
            if engine == "nwchem"
            or settings.get("scf_energy_tolerance") is not None
            else None
        ),
    }


def evaluate_endpoint_path_numerical_lineage(
    reactant: Artifact,
    product: Artifact,
    configured_method: dict[str, Any],
) -> dict[str, Any]:
    """Require stationary endpoints at the path's declared numerical resolution."""

    path_settings = path_numerical_settings(configured_method)
    endpoint_settings: dict[str, dict[str, Any]] = {}
    reasons: list[str] = []
    for role, endpoint in (("reactant", reactant), ("product", product)):
        if endpoint.qc.get("numerical_method_evidence_validated") is not True:
            reasons.append(
                f"{role}_numerical_method_evidence_validated_is_not_true"
            )
        validated = (endpoint.provenance or {}).get(
            "validated_numerical_settings"
        )
        if not isinstance(validated, dict):
            reasons.append(f"{role}_validated_numerical_settings_missing")
            validated = {}
        observed = path_numerical_settings(endpoint.method or {})
        endpoint_settings[role] = observed
        for key in PATH_NUMERICAL_KEYS:
            if observed.get(key) != validated.get(key):
                reasons.append(
                    f"{role}_{key}_does_not_match_validated_provenance"
                )
        if observed.get("grid") != path_settings.get("grid"):
            reasons.append(f"{role}_grid_does_not_match_path_setting")
        endpoint_scf = observed.get("scf_energy_tolerance")
        path_scf = path_settings.get("scf_energy_tolerance")
        if (
            endpoint_scf is None
            or path_scf is None
            or float(endpoint_scf) > float(path_scf)
        ):
            reasons.append(
                f"{role}_scf_energy_tolerance_is_coarser_than_path_setting"
            )
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "path_numerical_settings": path_settings,
        "endpoint_numerical_settings": endpoint_settings,
    }


def _endpoint_method(endpoint: Artifact) -> tuple[dict[str, Any], list[str]]:
    reasons: list[str] = []
    method = endpoint.method or {}
    provenance = endpoint.provenance or {}
    validated = provenance.get("validated_qm_method") or provenance.get("validated_nwchem_method")
    if not isinstance(validated, dict):
        validated = {}
        reasons.append("validated_method_provenance_missing")

    values = {
        "engine": _lower_text(method.get("engine")),
        "functional": _lower_text(method.get("functional")),
        "basis": _lower_text(method.get("basis")),
        "disp_vdw": exact_int(method.get("disp_vdw")),
        "program_version": _lower_text(validated.get("program_version")),
    }
    validated_values = {
        "engine": _lower_text(validated.get("engine")),
        "functional": _lower_text(validated.get("functional")),
        "basis": _lower_text(validated.get("basis")),
        "disp_vdw": exact_int(validated.get("disp_vdw")),
        "program_version": _lower_text(validated.get("program_version")),
    }
    for key in METHOD_KEYS:
        if values[key] is None:
            reasons.append(f"{key}_missing")
        if validated_values[key] is None:
            reasons.append(f"validated_{key}_missing")
        if values[key] != validated_values[key]:
            reasons.append(f"{key}_does_not_match_validated_provenance")
    if endpoint.qc.get("state_method_evidence_validated") is not True:
        reasons.append("state_method_evidence_validated_is_not_true")
    return values, reasons


def _endpoint_state(endpoint: Artifact) -> tuple[dict[str, int | None], list[str]]:
    reasons: list[str] = []
    data = endpoint.data or {}
    state = {
        "charge": exact_int(data.get("resolved_charge")),
        "multiplicity": exact_int(data.get("resolved_multiplicity")),
        "electron_count": exact_int(data.get("electron_count")),
    }
    validated = (endpoint.provenance or {}).get("validated_electronic_state")
    if not isinstance(validated, dict):
        validated = {}
        reasons.append("validated_electronic_state_provenance_missing")
    for key, value in state.items():
        if value is None:
            reasons.append(f"{key}_missing")
        if exact_int(validated.get(key)) != value:
            reasons.append(f"{key}_does_not_match_validated_provenance")
    return state, reasons


def evaluate_endpoint_method_lineage(
    reactant: Artifact,
    product: Artifact,
    configured_method: dict[str, Any],
    backend_name: str | None,
) -> dict[str, Any]:
    """Require both validated endpoints and the next calculation to share a PES."""

    configured = configured_qm_method(configured_method, backend_name)
    reasons: list[str] = []
    for key, value in configured.items():
        if value is None:
            reasons.append(f"configured_{key}_missing")

    endpoint_methods: dict[str, dict[str, Any]] = {}
    endpoint_states: dict[str, dict[str, int | None]] = {}
    for role, endpoint in (("reactant", reactant), ("product", product)):
        method, method_reasons = _endpoint_method(endpoint)
        state, state_reasons = _endpoint_state(endpoint)
        endpoint_methods[role] = method
        endpoint_states[role] = state
        reasons.extend(f"{role}_{reason}" for reason in method_reasons)
        reasons.extend(f"{role}_{reason}" for reason in state_reasons)
        for key in METHOD_KEYS:
            if method[key] != configured[key]:
                reasons.append(f"{role}_{key}_does_not_match_configured_method")

    if endpoint_methods["reactant"] != endpoint_methods["product"]:
        reasons.append("endpoint_methods_do_not_match")
    if endpoint_states["reactant"] != endpoint_states["product"]:
        reasons.append("endpoint_electronic_states_do_not_match")

    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "configured_method": configured,
        "endpoint_methods": endpoint_methods,
        "electronic_state": endpoint_states["reactant"],
        "reactant_endpoint_artifact_id": reactant.artifact_id,
        "product_endpoint_artifact_id": product.artifact_id,
    }


def evaluate_artifact_method_lineage(
    artifact: Artifact, endpoint_lineage: dict[str, Any]
) -> dict[str, Any]:
    """Check that a successful TS/path artifact carries the same audited lineage."""

    expected_method = endpoint_lineage.get("configured_method") or {}
    expected_state = endpoint_lineage.get("electronic_state") or {}
    method = artifact.method or {}
    data = artifact.data or {}
    observed_method = {
        "engine": _lower_text(method.get("engine") or method.get("program")),
        "functional": _lower_text(method.get("functional")),
        "basis": _lower_text(method.get("basis")),
        "disp_vdw": exact_int(method.get("disp_vdw")),
        "program_version": _lower_text(data.get("program_version")),
    }
    reasons: list[str] = []
    if _lower_text(method.get("required_program_version")) != observed_method["program_version"]:
        reasons.append("required_program_version_does_not_match_observed_version")
    for key in METHOD_KEYS:
        if observed_method[key] is None:
            reasons.append(f"{key}_missing")
        if observed_method[key] != expected_method.get(key):
            reasons.append(f"{key}_does_not_match_endpoint_method")

    observed_state = {
        "charge": exact_int(data.get("resolved_charge")),
        "multiplicity": exact_int(data.get("resolved_multiplicity")),
        "electron_count": exact_int(data.get("electron_count")),
    }
    provenance_state = (artifact.provenance or {}).get("electronic_state")
    if not isinstance(provenance_state, dict):
        provenance_state = {}
        reasons.append("electronic_state_provenance_missing")
    for key, expected in expected_state.items():
        if observed_state.get(key) != expected:
            reasons.append(f"{key}_does_not_match_endpoint_state")
        if exact_int(provenance_state.get(key)) != expected:
            reasons.append(f"provenance_{key}_does_not_match_endpoint_state")

    expected_endpoint_ids = {
        "reactant": endpoint_lineage.get("reactant_endpoint_artifact_id"),
        "product": endpoint_lineage.get("product_endpoint_artifact_id"),
    }
    if (artifact.provenance or {}).get("endpoint_artifact_ids") != expected_endpoint_ids:
        reasons.append("endpoint_artifact_ids_provenance_mismatch")
    if artifact.qc.get("method_evidence_validated") is not True:
        reasons.append("method_evidence_validated_is_not_true")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "artifact_id": artifact.artifact_id,
        "observed_method": observed_method,
        "observed_electronic_state": observed_state,
        "expected_method": expected_method,
        "expected_electronic_state": expected_state,
    }
