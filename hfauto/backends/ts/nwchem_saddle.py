"""NWChem interior-bracket saddle search for narrow molecular barriers."""

from __future__ import annotations

from math import isfinite
from pathlib import Path
from typing import Any

from hfauto.backends.qm.nwchem import NWChemEngine, species_xyz_path
from hfauto.backends.ts.base import (
    TSSearchResult,
    endpoint_electronic_energies,
    make_ts_species_artifact,
    plan_saddle_recovery,
    validate_saddle_seed_hessian,
    validate_ts_frequency_calculation,
)
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.methods import electronic_structure_method
from hfauto.chemistry.saddle_seeds import (
    aligned_interpolation_xyz,
    bracket_energy_maximum,
    endpoint_biased_fractions,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.chemistry.xyz_trajectory import endpoint_mode_reference
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


def _scan_calculation_accepted(
    calculation: Artifact, method: dict[str, Any]
) -> bool:
    energy = calculation.data.get("electronic_energy_hartree")
    try:
        finite_energy = isfinite(float(energy))
    except (TypeError, ValueError):
        finite_energy = False
    required_version = method.get("required_program_version")
    return bool(
        calculation.status.status == "success"
        and finite_energy
        and calculation.qc.get("scf_converged") is True
        and calculation.qc.get("normal_termination") is True
        and calculation.qc.get("real_qm_executed") is True
        and calculation.qc.get("fallback_dummy") is False
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


class NWChemSaddleEngine:
    """Bracket a narrow barrier, then refine and frequency-check one saddle."""

    name = "nwchem_saddle"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _refine_candidate(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict[str, Any],
        workdir: Path,
        *,
        candidate: dict[str, Any],
        state: dict[str, Any],
        search_evidence_validated: bool,
        artifacts: list[Artifact],
        attempt: Artifact,
    ) -> TSSearchResult:
        """Refine one evidenced seed; frequency/mode validation remains final."""

        reaction_id = _reaction_id(reaction)
        strategy = str(method.get("path_strategy") or "saddle_refinement")
        seed_path = Path(str(candidate.get("xyz_path") or ""))
        if not seed_path.is_file() or not search_evidence_validated:
            fail = Artifact.failure(
                f"saddle_refinement_failed_{reaction_id}",
                "ts_result",
                "Saddle refinement requires a readable, evidenced seed",
                category="saddle_refinement_seed_invalid",
                parents=[attempt.artifact_id],
            )
            return TSSearchResult(
                [*artifacts, fail], record=fail.model_dump(), success=False
            )
        ts_id = (
            reaction_id.replace("rxn_", "ts_", 1)
            if reaction_id.startswith("rxn_")
            else f"ts_{reaction_id}"
        )
        ts_species = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            seed_path,
            ts_id,
            self.name,
            {
                "ts_guess_strategy": strategy,
                "directional_curvature_hartree": candidate.get(
                    "directional_curvature_hartree"
                ),
                "real_ts_search_executed": True,
            },
        )
        mode_reference = candidate.get("reaction_mode_reference")
        if not isinstance(mode_reference, dict):
            try:
                mode_reference = endpoint_mode_reference(
                    read_xyz(species_xyz_path(reactant)),
                    read_xyz(species_xyz_path(product)),
                )
            except ValueError:
                mode_reference = None
        if isinstance(mode_reference, dict):
            ts_species.data["reaction_mode_reference"] = mode_reference
        refinement_method = {
            **method,
            "method_id": method.get(
                "saddle_method_id", method.get("method_id", "nwchem_saddle")
            ),
        }
        if method.get("saddle_seed_hessian_precheck", True):
            supplied_hessian = method.get("saddle_seed_hessian_evidence")
            seed_hessian = (
                supplied_hessian
                if isinstance(supplied_hessian, Artifact)
                else Artifact.model_validate(supplied_hessian)
                if isinstance(supplied_hessian, dict)
                else NWChemEngine().frequency(
                    ts_species,
                    refinement_method,
                    str(workdir / "seed_hessian"),
                )
            )
            seed_assessment = validate_saddle_seed_hessian(
                ts_species,
                seed_hessian,
                method,
                backend=self.name,
            )
            artifacts.extend([ts_species, seed_hessian])
            attempt.data["seed_hessian_calculation_id"] = (
                seed_hessian.artifact_id
            )
            attempt.qc.update(
                {
                    "saddle_seed_hessian_accepted": seed_assessment[
                        "accepted"
                    ],
                    "seed_n_imag": seed_hessian.data.get("n_imag"),
                    "seed_mode_overlap_score": seed_assessment[
                        "mode_overlap_score"
                    ],
                    "seed_stationary_point_validated": False,
                    "seed_hessian_reused": supplied_hessian is not None,
                    "seed_selected_mode": seed_assessment.get("spectral_qc", {}).get(
                        "selected_mode"
                    ),
                }
            )
            if not seed_assessment["accepted"]:
                attempt.data.update(
                    {
                        "diagnosis": "saddle_seed_hessian_rejected",
                        "next_action": "adaptive_path_refinement",
                    }
                )
                fail = Artifact.failure(
                    f"saddle_seed_hessian_failed_{reaction_id}",
                    "ts_result",
                    "Saddle seed Hessian did not contain a negative mode aligned "
                    "with the evidenced reaction direction",
                    category="saddle_seed_hessian_rejected",
                    parents=[
                        attempt.artifact_id,
                        ts_species.artifact_id,
                        seed_hessian.artifact_id,
                    ],
                    n_imag=seed_hessian.data.get("n_imag"),
                    mode_overlap_score=seed_assessment[
                        "mode_overlap_score"
                    ],
                )
                return TSSearchResult(
                    [*artifacts, fail], record=fail.model_dump(), success=False
                )
            selected_mode = seed_assessment.get("spectral_qc", {}).get(
                "selected_mode"
            )
            negative_mode_count = seed_hessian.data.get("n_imag")
            driver_mode_mapping_validated = negative_mode_count == 1
            attempt.qc["seed_driver_mode_mapping_validated"] = (
                driver_mode_mapping_validated
            )
            if not driver_mode_mapping_validated:
                attempt.data.update(
                    {
                        "diagnosis": "saddle_seed_mode_mapping_ambiguous",
                        "next_action": "adaptive_path_refinement",
                    }
                )
                fail = Artifact.failure(
                    f"saddle_seed_mode_mapping_failed_{reaction_id}",
                    "ts_result",
                    "NWChem DRIVER mode numbering is not uniquely mappable from "
                    "a seed Hessian with multiple negative modes",
                    category="saddle_seed_mode_mapping_ambiguous",
                    parents=[
                        attempt.artifact_id,
                        ts_species.artifact_id,
                        seed_hessian.artifact_id,
                    ],
                    n_imag=negative_mode_count,
                    selected_frequency_mode=(selected_mode or {}).get(
                        "mode_number"
                    ),
                )
                return TSSearchResult(
                    [*artifacts, fail], record=fail.model_dump(), success=False
                )
            if isinstance(selected_mode, dict) and selected_mode.get(
                "mode_number"
            ) is not None:
                refinement_method.update(
                    {
                        # One negative direction makes the DRIVER mapping unique.
                        # Frequency-output mode indices are otherwise not DRIVER
                        # mode indices and must not be copied to MODDIR.
                        "saddle_mode_number": 1,
                        "saddle_follow_first_negative": False,
                        "saddle_initial_hessian_only": True,
                    }
                )
        else:
            artifacts.append(ts_species)
        calculation = NWChemEngine().saddle_frequency(
            ts_species,
            refinement_method,
            str(workdir / "saddle_frequency"),
        )
        assessment = validate_ts_frequency_calculation(
            ts_species,
            calculation,
            method,
            backend=self.name,
            search_method_evidence_validated=search_evidence_validated,
        )
        validated = bool(assessment["validated"])
        recovery = (
            {
                "diagnosis": "saddle_validated",
                "next_action": "validate_connectivity",
                "restart_allowed": False,
                "method_updates": {},
            }
            if validated
            else plan_saddle_recovery(calculation, assessment, method)
        )
        attempt.data.update(
            {
                "diagnosis": recovery["diagnosis"],
                "next_action": recovery["next_action"],
                "saddle_recovery": recovery,
                "saddle_frequency_calculation_id": calculation.artifact_id,
            }
        )
        attempt.qc.update(
            {
                "ts_validated_by_frequency": validated,
                "n_imag": calculation.data.get("n_imag"),
                "mode_overlap_score": assessment["mode_overlap_score"],
                "ts_energy_above_endpoints": assessment.get(
                    "ts_energy_above_endpoints"
                ),
            }
        )
        artifacts.append(calculation)
        if not validated:
            fail = Artifact.failure(
                f"saddle_frequency_failed_{reaction_id}",
                "ts_result",
                "Saddle refinement did not produce a validated first-order saddle",
                category=str(recovery["diagnosis"]),
                parents=[
                    attempt.artifact_id,
                    ts_species.artifact_id,
                    calculation.artifact_id,
                ],
                n_imag=calculation.data.get("n_imag"),
                mode_overlap_score=assessment["mode_overlap_score"],
            )
            return TSSearchResult(
                [*artifacts, fail], record=fail.model_dump(), success=False
            )

        method_record = electronic_structure_method(
            engine="nwchem",
            backend=self.name,
            task="bracketed_saddle_frequency",
            config=method,
            charge=state["charge"],
            multiplicity=state["multiplicity"],
        )
        method_record["disp_vdw"] = method_record.pop("dispersion", None)
        reaction_validated = Artifact(
            artifact_id=reaction.artifact_id + "_with_ts",
            artifact_type="reaction_validated",
            parents=[
                reaction.artifact_id,
                attempt.artifact_id,
                ts_species.artifact_id,
                calculation.artifact_id,
            ],
            data={
                **reaction.data,
                "ts_species_id": ts_species.artifact_id,
                "ts_calc_id": calculation.artifact_id,
                "ts_backend": self.name,
                "resolved_charge": state["charge"],
                "resolved_multiplicity": state["multiplicity"],
                "electron_count": state["electron_count"],
                "program_version": calculation.data.get("program_version"),
            },
            qc={
                "ts_engine": self.name,
                "ts_found": True,
                "ts_validated_by_frequency": True,
                "n_imag": calculation.data.get("n_imag"),
                "imag_freq_cm1": calculation.data.get("imag_freq_cm1"),
                "mode_overlap_score": assessment["mode_overlap_score"],
                "real_ts_search_executed": True,
                "real_qm_executed": True,
                "method_evidence_validated": True,
                "dispersion_applied": calculation.qc.get("dispersion_applied"),
                "fallback_dummy": False,
            },
            method=method_record,
            provenance={
                "created_by": "NWChemSaddleEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "saddle_attempt_artifact_id": attempt.artifact_id,
                "saddle_frequency_calculation_id": calculation.artifact_id,
            },
        )
        artifacts.append(reaction_validated)
        return TSSearchResult(
            artifacts,
            record=reaction_validated.data,
            success=True,
            ts_species_id=ts_species.artifact_id,
            ts_calc_id=calculation.artifact_id,
        )

    def search_ts(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method_in: dict[str, Any],
        workdir: str | Path,
    ) -> TSSearchResult:
        method = {**self.config, **(method_in or {})}
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        reaction_id = _reaction_id(reaction)
        strategy = str(method.get("path_strategy") or "")
        if strategy not in {"bracketed_saddle_search", "saddle_refinement"}:
            fail = Artifact.failure(
                f"saddle_failed_{reaction_id}",
                "ts_result",
                "NWChem saddle bracketing was not authorized",
                category="saddle_strategy_not_authorized",
                parents=[reaction.artifact_id],
            )
            return TSSearchResult([fail], record=fail.model_dump(), success=False)

        if strategy == "saddle_refinement":
            candidate = dict(method.get("saddle_seed_candidate") or {})
            structure = read_xyz(species_xyz_path(reactant))
            state = resolve_electronic_state(
                reactant.data, method, structure.symbols
            )
            attempt = Artifact(
                artifact_id="saddle_attempt_"
                + fingerprint_dict(
                    {
                        "reaction": reaction_id,
                        "strategy": strategy,
                        "seed": candidate.get("xyz_path"),
                    }
                ),
                artifact_type="saddle_attempt",
                parents=[reaction.artifact_id],
                paths={"seed_xyz": str(candidate.get("xyz_path") or "")},
                data={
                    "reaction_id": reaction_id,
                    "strategy": strategy,
                    "engine": self.name,
                    "diagnosis": "resolved_saddle_candidate",
                    "next_action": "refine_saddle",
                    "candidate": candidate,
                    "all_scan_points_valid": bool(
                        method.get("saddle_seed_evidence_validated")
                    ),
                },
                qc={
                    "saddle_seed_resolved": bool(candidate),
                    "all_scan_points_valid": bool(
                        method.get("saddle_seed_evidence_validated")
                    ),
                },
                provenance={"created_by": "NWChemSaddleEngine"},
            )
            return self._refine_candidate(
                reaction,
                reactant,
                product,
                method,
                workdir,
                candidate=candidate,
                state=state,
                search_evidence_validated=bool(
                    method.get("saddle_seed_evidence_validated")
                ),
                artifacts=[attempt],
                attempt=attempt,
            )

        basin = dict(method.get("basin_assessment") or {})
        segment = dict(method.get("saddle_bracket_segment") or {})
        endpoint_energies = endpoint_electronic_energies(basin)
        if segment:
            try:
                segment_start_xyz = Path(str(segment["start_xyz"])).resolve()
                segment_target_xyz = Path(str(segment["target_xyz"])).resolve()
                segment_start_energy = float(segment["start_energy_hartree"])
                segment_target_energy = float(segment["target_energy_hartree"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("invalid saddle_bracket_segment") from exc
            if (
                segment.get("accepted") is not True
                or not segment_start_xyz.is_file()
                or not segment_target_xyz.is_file()
                or not isfinite(segment_start_energy)
                or not isfinite(segment_target_energy)
            ):
                raise ValueError("saddle_bracket_segment evidence is incomplete")
            endpoint_energies = {
                "reactant": segment_start_energy,
                "product": segment_target_energy,
            }
        if endpoint_energies is None:
            fail = Artifact.failure(
                f"saddle_failed_{reaction_id}",
                "ts_result",
                "Validated endpoint energies are required for saddle bracketing",
                category="saddle_bracket_endpoint_energy_missing",
                parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            )
            return TSSearchResult([fail], record=fail.model_dump(), success=False)
        reactant_energy = endpoint_energies["reactant"]
        product_energy = endpoint_energies["product"]

        if segment:
            start, target = reactant, product
            start_label = str(segment.get("start_label", "path_segment_start"))
            start_energy, target_energy = reactant_energy, product_energy
            start_xyz, target_xyz = segment_start_xyz, segment_target_xyz
            if reactant_energy < product_energy:
                start, target = target, start
                start_label = str(segment.get("target_label", "path_segment_end"))
                start_energy, target_energy = product_energy, reactant_energy
                start_xyz, target_xyz = target_xyz, start_xyz
        elif reactant_energy >= product_energy:
            start, target = reactant, product
            start_label = "reactant"
            start_energy, target_energy = reactant_energy, product_energy
        else:
            start, target = product, reactant
            start_label = "product"
            start_energy, target_energy = product_energy, reactant_energy
        if not segment:
            start_xyz = species_xyz_path(start)
            target_xyz = species_xyz_path(target)
        structure = read_xyz(start_xyz)
        state = resolve_electronic_state(start.data, method, structure.symbols)
        fractions = endpoint_biased_fractions(
            method.get("saddle_bracket_fractions")
        )
        endpoint_source = (
            "validated_reaction_path_image"
            if segment
            else "frequency_validated_minimum"
        )
        points: list[dict[str, Any]] = [
            {
                "fraction": 0.0,
                "energy_hartree": start_energy,
                "xyz_path": str(start_xyz),
                "source": endpoint_source,
            },
            {
                "fraction": 1.0,
                "energy_hartree": target_energy,
                "xyz_path": str(target_xyz),
                "source": endpoint_source,
            },
        ]
        artifacts: list[Artifact] = []
        scan_ok = True
        for index, fraction in enumerate(fractions[1:-1], start=1):
            point_dir = workdir / f"point_{index:02d}"
            seed_path = aligned_interpolation_xyz(
                start_xyz,
                target_xyz,
                point_dir / "seed.xyz",
                fraction,
            )
            token = fingerprint_dict(
                {"reaction": reaction_id, "fraction": fraction, "start": start_label}
            )
            seed = Artifact(
                artifact_id=f"saddle_seed_{token}",
                artifact_type="species",
                parents=[reaction.artifact_id, start.artifact_id, target.artifact_id],
                paths={"xyz": str(seed_path)},
                data={
                    **start.data,
                    "species_id": f"saddle_seed_{token}",
                    "state": "saddle_seed",
                    "xyz_path": str(seed_path),
                    "reaction_id": reaction_id,
                    "interpolation_fraction": fraction,
                    "interpolation_start": start_label,
                    "reaction_coordinate": reaction.data.get("reaction_coordinate"),
                    "bond_changes": list(reaction.data.get("bond_changes") or []),
                },
                qc={"scientific_role": "saddle_seed_only"},
                provenance={"created_by": "NWChemSaddleEngine"},
            )
            calculation = NWChemEngine().single_point(
                seed,
                {
                    **method,
                    "method_id": method.get(
                        "saddle_scan_method_id",
                        method.get("method_id", "nwchem_saddle_scan"),
                    ),
                },
                str(point_dir / "single_point"),
            )
            accepted = _scan_calculation_accepted(calculation, method)
            calculation.qc["saddle_seed_energy_accepted"] = accepted
            scan_ok = scan_ok and accepted
            artifacts.extend([seed, calculation])
            if accepted:
                points.append(
                    {
                        "fraction": fraction,
                        "energy_hartree": float(
                            calculation.data["electronic_energy_hartree"]
                        ),
                        "xyz_path": str(seed_path),
                        "species_artifact_id": seed.artifact_id,
                        "calculation_artifact_id": calculation.artifact_id,
                        "source": "nwchem_single_point",
                    }
                )

        profile = bracket_energy_maximum(
            points,
            minimum_prominence_hartree=float(
                method.get(
                    "saddle_bracket_minimum_prominence_hartree",
                    method.get(
                        "ts_barrier_resolution_hartree",
                        method.get(
                            "ts_endpoint_energy_tolerance_hartree", 1.0e-5
                        ),
                    ),
                )
            ),
        )
        scan_id = "saddle_seed_scan_" + fingerprint_dict(
            {"reaction": reaction_id, "fractions": fractions, "start": start_label}
        )
        scan = Artifact(
            artifact_id=scan_id,
            artifact_type="saddle_seed_scan",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            data={
                "reaction_id": reaction_id,
                "strategy": strategy,
                "engine": self.name,
                "start_endpoint": start_label,
                "bracket_segment": segment or None,
                "points": sorted(points, key=lambda point: point["fraction"]),
                **profile,
            },
            method=electronic_structure_method(
                engine="nwchem",
                backend=self.name,
                task="saddle_seed_scan",
                config=method,
                charge=state["charge"],
                multiplicity=state["multiplicity"],
            ),
            qc={
                "all_scan_points_valid": scan_ok,
                "internal_maximum_resolved": profile["candidate"] is not None,
                "real_qm_executed": bool(artifacts),
                "fallback_dummy": False,
            },
            provenance={"created_by": "NWChemSaddleEngine"},
        )
        artifacts.append(scan)
        candidate = profile["candidate"] if scan_ok else None
        if candidate is not None:
            candidate = {
                **candidate,
                "reaction_mode_reference": endpoint_mode_reference(
                    read_xyz(start_xyz),
                    read_xyz(target_xyz),
                    source=(
                        f"validated_path_segment:{segment.get('start_index')}-"
                        f"{segment.get('end_index')}"
                        if segment
                        else "validated_reaction_endpoints"
                    ),
                ),
            }
        attempt = Artifact(
            artifact_id="saddle_attempt_"
            + fingerprint_dict({"reaction": reaction_id, "scan": scan_id}),
            artifact_type="saddle_attempt",
            parents=[scan.artifact_id],
            paths={"seed_xyz": str(candidate["xyz_path"]) if candidate else ""},
            data={
                "reaction_id": reaction_id,
                "strategy": strategy,
                "engine": self.name,
                "diagnosis": (
                    "resolved_saddle_candidate"
                    if candidate
                    else "multiple_step_candidate"
                    if profile["classification"] == "multiple_internal_maxima"
                    else "saddle_seed_bracket_failed"
                ),
                "next_action": (
                    "refine_saddle"
                    if candidate
                    else "validate_intermediate_basins"
                    if profile["classification"] == "multiple_internal_maxima"
                    else "adaptive_path_refinement"
                ),
                "candidate": candidate,
                "all_scan_points_valid": scan_ok,
            },
            qc={
                "all_scan_points_valid": scan_ok,
                "saddle_seed_resolved": candidate is not None,
                "method_evidence_validated": scan_ok,
            },
            provenance={"created_by": "NWChemSaddleEngine"},
        )
        artifacts.append(attempt)
        if candidate is None or method.get("saddle_bracket_only", False):
            fail = Artifact.failure(
                f"saddle_failed_{reaction_id}",
                "ts_result",
                (
                    "Multiple internal maxima require intermediate-basin validation"
                    if profile["classification"] == "multiple_internal_maxima"
                    else "No internal maximum was resolved by the saddle seed scan"
                    if candidate is None
                    else "Saddle bracket-only execution completed without refinement"
                ),
                category=(
                    "multiple_step_candidate"
                    if profile["classification"] == "multiple_internal_maxima"
                    else "saddle_seed_bracket_failed"
                    if candidate is None
                    else "saddle_bracket_only_completed"
                ),
                parents=[attempt.artifact_id],
                profile=profile,
            )
            return TSSearchResult(
                [*artifacts, fail], record=fail.model_dump(), success=False
            )

        return self._refine_candidate(
            reaction,
            reactant,
            product,
            method,
            workdir,
            candidate=candidate,
            state=state,
            search_evidence_validated=scan_ok,
            artifacts=artifacts,
            attempt=attempt,
        )
