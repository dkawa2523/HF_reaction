from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from typing import Any

from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.reaction_path_qc import (
    estimate_reaction_mode_overlap,
    make_midpoint_ts_xyz,
    reaction_coordinate_value,
    select_reaction_coordinate_imaginary_mode,
)
from hfauto.chemistry.xyz import xyz_files_have_same_atom_order
from hfauto.core.artifacts import canonical_species_id
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.frequency_qc import (
    imaginary_frequency_cutoff_from_method,
    is_significant_imaginary_frequency,
)
from hfauto.core.qc import ts_qc
from hfauto.core.schemas.artifact import Artifact


@dataclass
class TSSearchResult:
    """One explicit result shape shared by all TS-search backends."""

    artifacts: list[Artifact]
    record: Any = None
    success: bool | None = None
    ts_species_id: str | None = None
    ts_calc_id: str | None = None
    ts_species: Artifact | None = field(init=False, default=None)
    calculation: Artifact | None = field(init=False, default=None)

    def __post_init__(self) -> None:
        self.ts_species = next(
            (
                artifact
                for artifact in self.artifacts
                if artifact.artifact_type == "species"
                and artifact.data.get("state") == "transition_state"
            ),
            None,
        )
        self.calculation = next(
            (
                artifact
                for artifact in self.artifacts
                if artifact.artifact_type == "calculation"
            ),
            None,
        )
        if self.ts_species is not None and self.ts_species_id is None:
            self.ts_species_id = self.ts_species.artifact_id
        if self.calculation is not None and self.ts_calc_id is None:
            self.ts_calc_id = self.calculation.artifact_id
        if self.success is None:
            self.success = not any(
                artifact.status.status == "failed" for artifact in self.artifacts
            )
        self.success = bool(self.success)


@dataclass
class IRCResult:
    """One explicit result shape shared by all IRC backends."""

    artifacts: list[Artifact]
    record: Any = None
    success: bool | None = None
    artifact: Artifact | None = field(init=False, default=None)

    def __post_init__(self) -> None:
        self.artifact = next(
            (
                artifact
                for artifact in self.artifacts
                if artifact.artifact_type == "irc"
            ),
            self.artifacts[0] if self.artifacts else None,
        )
        if self.success is None:
            self.success = not any(
                artifact.status.status == "failed" for artifact in self.artifacts
            )
        self.success = bool(self.success)


def midpoint_ts_xyz(reactant: Artifact, product: Artifact, out_xyz: str | Path) -> Path:
    return make_midpoint_ts_xyz(
        reactant.data.get("xyz_path") or reactant.paths.get("xyz"),
        product.data.get("xyz_path") or product.paths.get("xyz"),
        out_xyz,
    )


def make_ts_species_artifact(
    reaction: Artifact,
    reactant: Artifact,
    product: Artifact,
    ts_xyz: str | Path,
    ts_id: str,
    source: str,
    extra_qc: dict[str, Any] | None = None,
) -> Artifact:
    data = {
        **reactant.data,
        "species_id": ts_id,
        "source_species_id": ts_id,
        "state": "transition_state",
        "xyz_path": str(ts_xyz),
        "reactant_species_id": canonical_species_id(reactant),
        "product_species_id": canonical_species_id(product),
        "reaction_id": reaction.data.get("reaction_id", reaction.artifact_id),
        "reaction_type": reaction.data.get("reaction_type"),
        "mechanism_family": reaction.data.get("mechanism_family"),
        "bond_changes": list(reaction.data.get("bond_changes", []) or []),
        "reaction_coordinate": reaction.data.get("reaction_coordinate")
        or reactant.data.get("reaction_coordinate"),
        "ts_builder": source,
    }
    geom_qc = geometry_qc_from_xyz(data, ts_xyz)
    progress = reaction_coordinate_progress_score(
        reactant,
        product,
        ts_xyz,
        {
            "reaction_coordinate": data.get("reaction_coordinate") or {},
            "bond_changes": data.get("bond_changes") or [],
        },
    )
    return Artifact(
        artifact_id=ts_id,
        artifact_type="species",
        parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
        paths={
            "xyz": str(ts_xyz),
            "reactant_xyz": reactant.data.get("xyz_path", "")
            or reactant.paths.get("xyz", ""),
            "product_xyz": product.data.get("xyz_path", "")
            or product.paths.get("xyz", ""),
        },
        data=data,
        method={"engine": source, "task": "ts_guess"},
        qc={
            "ts_guess_source": source,
            "geometry_qc": geom_qc,
            **progress,
            **(extra_qc or {}),
        },
    )


def endpoint_electronic_energies(
    basin_assessment: dict[str, Any] | None,
) -> dict[str, float] | None:
    """Return a complete finite same-PES endpoint energy pair, or no claim."""

    energies: dict[str, float] = {}
    basin = dict(basin_assessment or {})
    for label in ("reactant", "product"):
        evidence = dict(basin.get(f"{label}_evidence") or {})
        try:
            energy = float(evidence["electronic_energy_hartree"])
        except (KeyError, TypeError, ValueError):
            return None
        if not isfinite(energy):
            return None
        energies[label] = energy
    return energies


def validate_ts_frequency_calculation(
    ts_species: Artifact,
    calculation: Artifact,
    method: dict[str, Any],
    *,
    backend: str,
    search_method_evidence_validated: bool,
    overlap_evaluator=estimate_reaction_mode_overlap,
) -> dict[str, Any]:
    """Apply the one TS frequency/mode gate shared by path-search backends."""

    projection = _reaction_mode_projection(
        ts_species,
        calculation,
        method,
        overlap_evaluator=overlap_evaluator,
    )
    geometry_qc = projection["geometry_qc"]
    overlap = projection["mode_overlap_score"]
    overlap_method = projection["mode_overlap_method"]
    imaginary_cutoff = imaginary_frequency_cutoff_from_method(method)
    calculation.qc.update(
        ts_qc(
            calculation.data.get("n_imag"),
            calculation.data.get("imag_freq_cm1"),
            overlap,
            mode_overlap_threshold=float(
                method.get("ts_mode_overlap_threshold", 0.50)
            ),
            imaginary_frequency_cutoff_cm1=imaginary_cutoff,
        )
    )
    calculation.qc.update(
        {
            "mode_overlap_score": overlap,
            "mode_overlap_method": overlap_method,
            "geometry_qc": geometry_qc,
            "real_ts_search_executed": True,
        }
    )
    calculation.data.update(
        {
            "mode_overlap_score": overlap,
            "mode_overlap_method": overlap_method,
            "ts_backend": backend,
        }
    )
    endpoint_energies = endpoint_electronic_energies(
        method.get("basin_assessment")
    )
    ts_energy = calculation.data.get("electronic_energy_hartree")
    energy_order_validated = None
    activation_energy_hartree = None
    activation_energy_resolved = None
    if endpoint_energies is not None and ts_energy is not None:
        tolerance = float(
            method.get("ts_endpoint_energy_tolerance_hartree", 1.0e-5)
        )
        activation_energy_hartree = float(ts_energy) - max(
            endpoint_energies.values()
        )
        resolution = float(
            method.get("ts_barrier_resolution_hartree", tolerance)
        )
        if resolution <= 0.0:
            raise ValueError("ts_barrier_resolution_hartree must be positive")
        activation_energy_resolved = bool(
            activation_energy_hartree >= resolution
        )
        energy_order_validated = bool(
            float(ts_energy) + tolerance >= max(endpoint_energies.values())
        )
        calculation.qc.update(
            {
                "ts_energy_above_endpoints": energy_order_validated,
                "ts_endpoint_energy_tolerance_hartree": tolerance,
                "activation_energy_resolved": activation_energy_resolved,
                "activation_energy_resolution_hartree": resolution,
            }
        )
        calculation.data.update(
            {
                "electronic_activation_energy_hartree": (
                    activation_energy_hartree
                ),
                "electronic_activation_energy_kcal_mol": (
                    activation_energy_hartree * HARTREE_TO_KCAL_MOL
                ),
                "activation_energy_resolved": activation_energy_resolved,
                "activation_energy_resolution_hartree": resolution,
            }
        )
    required_version = method.get("required_program_version")
    validated = bool(
        calculation.status.status == "success"
        and calculation.qc.get("ts_validated_by_frequency") is True
        and calculation.qc.get("real_qm_executed") is True
        and calculation.qc.get("fallback_dummy") is False
        and search_method_evidence_validated
        and energy_order_validated is not False
        and (
            required_version is None
            or str(calculation.data.get("program_version"))
            == str(required_version)
        )
        and (
            method.get("disp_vdw") is None
            or calculation.qc.get("dispersion_applied") is True
        )
    )
    ts_species.data["program_version"] = calculation.data.get(
        "program_version"
    )
    ts_species.qc.update(
        {
            "method_evidence_validated": validated,
            "dispersion_applied": calculation.qc.get("dispersion_applied"),
            "real_qm_executed": calculation.qc.get("real_qm_executed") is True,
            "fallback_dummy": bool(
                calculation.qc.get("fallback_dummy", False)
            ),
        }
    )
    return {
        **projection,
        "validated": validated,
        "ts_energy_above_endpoints": energy_order_validated,
        "activation_energy_resolved": activation_energy_resolved,
        "electronic_activation_energy_hartree": activation_energy_hartree,
    }


def _reaction_mode_projection(
    species: Artifact,
    calculation: Artifact,
    method: dict[str, Any],
    *,
    overlap_evaluator=estimate_reaction_mode_overlap,
) -> dict[str, Any]:
    """Project the reported imaginary mode onto the declared reaction coordinate."""

    final_xyz = calculation.paths.get("final_xyz") or species.data.get(
        "xyz_path"
    )
    species.data["xyz_path"] = str(final_xyz)
    species.paths["xyz"] = str(final_xyz)
    geometry_qc = geometry_qc_from_xyz(species.data, final_xyz)
    output_text = None
    output_path = calculation.paths.get("output")
    if output_path and Path(output_path).is_file():
        output_text = Path(output_path).read_text(
            encoding="utf-8", errors="ignore"
        )
    imaginary_cutoff = imaginary_frequency_cutoff_from_method(method)
    overlap, overlap_method = overlap_evaluator(
        species.data,
        final_xyz,
        calculation.data.get("n_imag"),
        calculation.data.get("imag_freq_cm1"),
        output_text,
        mode_displacements=calculation.data.get(
            "imaginary_mode_displacements"
        ),
        mode_component_units=(
            (calculation.data.get("projected_imaginary_mode") or {}).get(
                "component_units"
            )
        ),
        allow_legacy_geometry_heuristic=bool(
            method.get("allow_legacy_mode_overlap_heuristic", False)
        ),
        imaginary_frequency_cutoff_cm1=imaginary_cutoff,
    )
    return {
        "final_xyz": str(final_xyz),
        "geometry_qc": geometry_qc,
        "mode_overlap_score": overlap,
        "mode_overlap_method": overlap_method,
    }


def validate_saddle_seed_hessian(
    seed_species: Artifact,
    calculation: Artifact,
    method: dict[str, Any],
    *,
    backend: str,
    overlap_evaluator=estimate_reaction_mode_overlap,
) -> dict[str, Any]:
    """Accept a fixed-geometry Hessian only as a mode-following seed.

    A non-stationary structure can have one imaginary mode, so this gate never
    validates a transition state.  It only decides whether saddle optimization
    may follow the mode.
    """

    modes = calculation.data.get("projected_imaginary_modes")
    selected_mode = select_reaction_coordinate_imaginary_mode(
        seed_species.data,
        calculation.paths.get("final_xyz") or seed_species.data.get("xyz_path"),
        modes if isinstance(modes, list) else None,
    )
    if selected_mode is not None:
        final_xyz = calculation.paths.get("final_xyz") or seed_species.data.get(
            "xyz_path"
        )
        projection = {
            "final_xyz": str(final_xyz),
            "geometry_qc": geometry_qc_from_xyz(seed_species.data, final_xyz),
            "mode_overlap_score": selected_mode["mode_overlap_score"],
            "mode_overlap_method": selected_mode["mode_overlap_method"],
        }
    else:
        # Compatibility for legacy/test artifacts that expose only one mode.
        projection = _reaction_mode_projection(
            seed_species,
            calculation,
            method,
            overlap_evaluator=overlap_evaluator,
        )
        if calculation.data.get("n_imag") == 1:
            selected_mode = {
                "mode_number": (
                    (calculation.data.get("projected_imaginary_mode") or {}).get(
                        "mode_number", 1
                    )
                ),
                "frequency_cm1": calculation.data.get("imag_freq_cm1"),
                "mode_overlap_score": projection["mode_overlap_score"],
                "mode_overlap_method": projection["mode_overlap_method"],
            }
    cutoff = imaginary_frequency_cutoff_from_method(method)
    overlap_threshold = float(
        method.get("saddle_seed_mode_overlap_threshold", 0.50)
    )
    overlap_margin_threshold = float(
        method.get("saddle_seed_mode_overlap_margin", 0.10)
    )
    overlap_margin = float(
        (selected_mode or {}).get("overlap_margin", 1.0)
    )
    target_mode_present = bool(
        selected_mode is not None
        and is_significant_imaginary_frequency(
            selected_mode.get("frequency_cm1"), cutoff
        )
        and float(selected_mode.get("mode_overlap_score") or 0.0)
        >= overlap_threshold
        and (
            calculation.data.get("n_imag") == 1
            or overlap_margin >= overlap_margin_threshold
        )
    )
    spectral_qc = {
        "seed_negative_mode_count": calculation.data.get("n_imag"),
        "seed_target_mode_present": target_mode_present,
        "seed_mode_overlap_threshold": overlap_threshold,
        "seed_mode_overlap_margin_threshold": overlap_margin_threshold,
        "imaginary_frequency_cutoff_cm1": cutoff,
        "selected_mode": selected_mode,
    }
    required_version = method.get("required_program_version")
    accepted = bool(
        calculation.status.status == "success"
        and calculation.data.get("frequency_count_complete") is True
        and target_mode_present
        and calculation.qc.get("real_qm_executed") is True
        and calculation.qc.get("fallback_dummy") is False
        and projection["geometry_qc"].get("geometry_sane") is True
        and (
            required_version is None
            or str(calculation.data.get("program_version"))
            == str(required_version)
        )
        and (
            method.get("disp_vdw") is None
            or calculation.qc.get("dispersion_applied") is True
        )
    )
    calculation.qc.update(
        {
            "saddle_seed_hessian_accepted": accepted,
            "seed_has_one_target_imaginary_mode": bool(
                calculation.data.get("n_imag") == 1 and target_mode_present
            ),
            "seed_has_target_imaginary_mode": target_mode_present,
            "seed_selected_mode_number": (
                selected_mode.get("mode_number") if selected_mode else None
            ),
            "seed_selected_mode_frequency_cm1": (
                selected_mode.get("frequency_cm1") if selected_mode else None
            ),
            "stationary_point_validated": False,
            "ts_validated_by_frequency": False,
            "mode_overlap_score": projection["mode_overlap_score"],
            "mode_overlap_method": projection["mode_overlap_method"],
            "geometry_qc": projection["geometry_qc"],
        }
    )
    calculation.data.update(
        {
            "mode_overlap_score": projection["mode_overlap_score"],
            "mode_overlap_method": projection["mode_overlap_method"],
            "saddle_seed_backend": backend,
            "scientific_role": "mode_following_seed_only",
            "selected_mode": selected_mode,
        }
    )
    return {
        **projection,
        "accepted": accepted,
        "stationary_point_validated": False,
        "spectral_qc": spectral_qc,
    }


def plan_saddle_recovery(
    calculation: Artifact,
    assessment: dict[str, Any],
    method: dict[str, Any],
) -> dict[str, Any]:
    """Choose one bounded recovery from evidence, not from backend errors alone."""

    final_xyz = str(calculation.paths.get("final_xyz") or "")
    final_geometry_available = bool(final_xyz and Path(final_xyz).is_file())
    restart_count = int(method.get("saddle_hessian_restart_count", 0))
    maximum_restarts = int(method.get("saddle_max_hessian_restarts", 1))
    if (
        calculation.status.status != "success"
        and calculation.data.get("geometry_converged") is False
    ):
        restart_allowed = bool(
            final_geometry_available and restart_count < maximum_restarts
        )
        if restart_allowed:
            trust = float(method.get("driver_trust", 0.15))
            saddle_step = float(method.get("driver_saddle_step", 0.05))
            return {
                "diagnosis": "saddle_optimizer_failed",
                "next_action": "refresh_hessian_and_restart",
                "restart_allowed": True,
                "restart_seed_xyz": final_xyz,
                "method_updates": {
                    "saddle_hessian_restart_count": restart_count + 1,
                    "driver_trust": max(0.03, trust * 0.5),
                    "driver_saddle_step": max(0.01, saddle_step * 0.4),
                    "saddle_initial_hessian_only": True,
                    "saddle_follow_first_negative": False,
                },
            }
        return {
            "diagnosis": "saddle_optimizer_failed",
            "next_action": "adaptive_double_ended_path",
            "restart_allowed": False,
            "restart_seed_xyz": final_xyz or None,
            "method_updates": {},
        }
    if calculation.data.get("geometry_converged") is True and not bool(
        calculation.data.get("frequency_count_complete")
    ):
        return {
            "diagnosis": "final_frequency_missing",
            "next_action": "run_fixed_geometry_frequency",
            "restart_allowed": False,
            "restart_seed_xyz": final_xyz or None,
            "method_updates": {},
        }
    if assessment.get("ts_energy_above_endpoints") is False:
        return {
            "diagnosis": "saddle_energy_order_invalid",
            "next_action": "adaptive_double_ended_path",
            "restart_allowed": False,
            "restart_seed_xyz": final_xyz or None,
            "method_updates": {},
        }
    if calculation.data.get("n_imag") != 1:
        return {
            "diagnosis": "saddle_order_invalid",
            "next_action": "adaptive_double_ended_path",
            "restart_allowed": False,
            "restart_seed_xyz": final_xyz or None,
            "method_updates": {},
        }
    if assessment.get("mode_overlap_score") is None:
        diagnosis = "saddle_mode_evidence_missing"
    else:
        diagnosis = "saddle_mode_mismatch"
    return {
        "diagnosis": diagnosis,
        "next_action": "adaptive_double_ended_path",
        "restart_allowed": False,
        "restart_seed_xyz": final_xyz or None,
        "method_updates": {},
    }


def endpoints_compatible(reactant: Artifact, product: Artifact) -> tuple[bool, dict[str, Any]]:
    rc_xyz = reactant.data.get("xyz_path") or reactant.paths.get("xyz")
    ip_xyz = product.data.get("xyz_path") or product.paths.get("xyz")
    if not rc_xyz or not ip_xyz:
        return False, {"reason": "missing_xyz"}
    try:
        atom_order_ok = xyz_files_have_same_atom_order(rc_xyz, ip_xyz)
    except (OSError, TypeError, ValueError) as exc:
        return False, {"reason": str(exc), "atom_order_ok": False}
    map_ok = reactant.data.get("atom_order_key") == product.data.get("atom_order_key")
    return bool(atom_order_ok and map_ok), {
        "atom_order_ok": atom_order_ok,
        "atom_order_key_ok": map_ok,
    }


def reaction_coordinate_progress_score(
    reactant: Artifact, product: Artifact, ts_xyz: str | Path, reaction_definition: dict[str, Any]
) -> dict[str, Any]:
    rc_data = {**reactant.data, **reaction_definition}
    ip_data = {**product.data, **reaction_definition}
    ts_data = {**reactant.data, **reaction_definition, "state": "transition_state"}
    q_rc = reaction_coordinate_value(
        rc_data, reactant.data.get("xyz_path") or reactant.paths.get("xyz")
    )
    q_ip = reaction_coordinate_value(
        ip_data, product.data.get("xyz_path") or product.paths.get("xyz")
    )
    q_ts = reaction_coordinate_value(ts_data, ts_xyz)
    if q_rc is None or q_ip is None or q_ts is None or abs(q_ip - q_rc) < 1e-12:
        return {
            "reaction_coordinate_progress_score": 0.0,
            "reaction_coordinate_between_endpoints": False,
            "q_reactant_A": q_rc,
            "q_product_A": q_ip,
            "q_ts_A": q_ts,
        }
    lo, hi = sorted([q_rc, q_ip])
    between = lo <= q_ts <= hi
    progress = (q_ts - q_rc) / (q_ip - q_rc)
    score = max(0.0, 1.0 - 2.0 * abs(float(progress) - 0.5))
    return {
        "reaction_coordinate_progress_score": float(score),
        "reaction_coordinate_between_endpoints": bool(between),
        "reaction_coordinate_progress_fraction": float(progress),
        "q_reactant_A": q_rc,
        "q_product_A": q_ip,
        "q_ts_A": q_ts,
    }
