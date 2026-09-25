"""NWChem NEB and saddle-point backend.

NEB supplies a reaction-path TS guess.  A separate saddle optimization and
frequency calculation supplies the scientific TS check.  NWChem NEB is never
reported as IRC: callers receive an explicit unsupported artifact instead.
"""

from __future__ import annotations

import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

from hfauto.backends.qm.nwchem import (
    NWChemEngine,
    NWChemInputRenderer,
    allow_nwchem_subprocess,
    parse_nwchem_output,
    run_nwchem_input,
)
from hfauto.backends.ts.base import (
    IRCResult,
    TSSearchResult,
    canonical_species_id,
    make_ts_species_artifact,
    validate_ts_frequency_calculation,
)
from hfauto.backends.ts.nwchem_path_support import (
    assess_nwchem_path_evidence,
    electronic_method_config,
    existing_file_hash,
    prepare_initial_path,
    resolve_path_endpoints,
    validate_path_trajectory,
    validated_initial_path,
)
from hfauto.chemistry.methods import electronic_structure_method
from hfauto.chemistry.path_diagnostics import (
    make_path_attempt_record,
    parse_neb_optimization_history,
)
from hfauto.chemistry.reaction_path_qc import estimate_reaction_mode_overlap
from hfauto.chemistry.reaction_profile import (
    analyze_reaction_path,
    select_path_ts_guess,
)
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.chemistry.xyz_trajectory import (
    reaction_mode_reference,
    read_xyz_trajectory,
)
from hfauto.core.executables import resolve_executable
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


def nwchem_neb_converged(output: str) -> bool:
    """Recognize only NWChem's affirmative NEB termination record."""

    return bool(
        re.search(
            r"^@neb\s+NEB\s+calculation\s+converged\s*$",
            output,
            flags=re.MULTILINE,
        )
    )


def _ts_id(reaction: Artifact) -> str:
    reaction_id = _reaction_id(reaction)
    return (
        reaction_id.replace("rxn_", "ts_", 1)
        if reaction_id.startswith("rxn_")
        else f"ts_{reaction_id}"
    )


def select_neb_ts_guess(path: str | Path) -> tuple[XYZ, int, str]:
    """Select the highest-energy internal bead, or the middle bead if unlabelled."""
    xyz, index, reason = select_path_ts_guess(str(path))
    compatibility_reason = {
        "highest_energy_internal_image": "highest_energy_internal_bead",
        "middle_internal_image_no_energy_labels": "middle_internal_bead_no_energy_labels",
    }[reason]
    return xyz, index, compatibility_reason


def classify_neb_energy_profile(
    path: str | Path,
    barrier_threshold_kcal_mol: float = 0.5,
    *,
    reaction_id: str = "unknown",
    engine: str = "nwchem_neb",
    converged: bool = True,
) -> dict[str, Any]:
    """Compatibility view backed by the normalized reaction-path record."""

    images = read_xyz_trajectory(path)
    record = analyze_reaction_path(
        reaction_id=reaction_id,
        engine=engine,
        comments=[image.comment for image in images],
        converged=converged,
        barrier_threshold_kcal_mol=barrier_threshold_kcal_mol,
    )
    classification = (
        "unresolved_missing_energies"
        if record.classification == "missing_profile"
        else record.classification
    )
    return {
        "reaction_path": record.model_dump(),
        "neb_profile_classification": classification,
        "neb_energies_hartree": [image.energy_hartree for image in record.images],
        "neb_highest_internal_bead": record.highest_internal_image_index,
        "neb_barrier_from_reactant_kcal_mol": record.barrier_from_reactant_kcal_mol,
        "neb_reaction_energy_kcal_mol": record.reaction_energy_kcal_mol,
        "neb_barrier_threshold_kcal_mol": record.barrier_threshold_kcal_mol,
    }


class NWChemNEBEngine:
    name = "nwchem_neb"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        return {**self.config, **(method or {})}

    def render_neb_input(
        self, reactant: Artifact, product: Artifact, method_in: dict[str, Any], path: str | Path
    ) -> Path:
        method = self._method(method_in)
        start, end, state = resolve_path_endpoints(reactant, product, method)
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        nbeads = max(
            3,
            int(method.get("path_image_count", method.get("nbeads", 9))),
        )
        initial_path: Path | None = None
        initial_path_source = method.get("path_initial_trajectory") or method.get(
            "initial_path_xyz"
        )
        if initial_path_source:
            if not Path(str(initial_path_source)).is_file():
                raise ValueError(
                    f"initial_path_xyz is not readable: {initial_path_source}"
                )
            initial_path, nbeads = validated_initial_path(
                str(initial_path_source),
                start,
                end,
                output.parent / "initial_path.xyz",
                endpoint_tolerance_A=float(
                    method.get("path_endpoint_rmsd_tolerance_A", 0.10)
                ),
                image_count=(
                    nbeads
                    if method.get(
                        "resample_path_initial_trajectory",
                        method.get("resample_initial_path", False),
                    )
                    else None
                ),
            )
        basis = str(method.get("basis", "def2-svp"))
        lines = [
            "start hfauto_neb",
            f'title "hfauto NEB {canonical_species_id(reactant)} to {canonical_species_id(product)}"',
            f"memory total {int(method.get('memory_mb', 2000))} mb",
            "",
            *NWChemInputRenderer._geometry("", start),
            *NWChemInputRenderer._geometry("endgeom", end),
            f"charge {state['charge']}",
            "",
            "basis spherical",
            f"  * library {basis}",
            "end",
            "",
            *NWChemInputRenderer.method_block(method, int(state["multiplicity"])),
        ]
        for block in method.get("extra_blocks", []) or []:
            lines += ["", str(block).rstrip()]
        lines += [
            "",
            "neb",
            f"  nbeads {nbeads}",
            f"  maxiter {int(method.get('neb_maxiter', 100))}",
            f"  stepsize {float(method.get('neb_stepsize', 0.1))}",
        ]
        if initial_path is not None:
            lines.append(f"  xyz_path {initial_path.name}")
        lines += [
            f"  print_shift {int(method.get('neb_print_shift', 1))}",
            "end",
            "",
            "task dft neb ignore",
        ]
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return output

    @staticmethod
    def _path_file(workdir: Path) -> Path | None:
        candidates = [*workdir.glob("*.nebpath*.xyz"), *workdir.glob("*neb*path*.xyz")]
        return max(candidates, key=lambda item: item.stat().st_mtime) if candidates else None

    def search_ts(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method_in: dict[str, Any],
        workdir: str | Path,
    ) -> TSSearchResult:
        method = self._method(method_in)
        wd = Path(workdir)
        try:
            start, end, state = resolve_path_endpoints(reactant, product, method)
            nbeads = max(
                3,
                int(method.get("path_image_count", method.get("nbeads", 9))),
            )
            initial_path, initialization = prepare_initial_path(
                reaction,
                start,
                end,
                method,
                wd / "neb" / "generated_initial_path.xyz",
                image_count=nbeads,
            )
            method["path_initialization"] = initialization
            if initial_path is not None:
                method["path_initial_trajectory"] = str(initial_path)
            input_path = self.render_neb_input(
                reactant, product, method, wd / "neb" / "nwchem_neb.nw"
            )
        except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
            fail = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                str(exc),
                category="endpoint_mismatch",
                parents=[reaction.artifact_id],
            )
            return TSSearchResult([fail], record=fail.model_dump(), success=False)
        method_record = electronic_structure_method(
            engine="nwchem",
            backend=self.name,
            task="neb_saddle_frequency",
            config={
                "functional": method.get("functional", method.get("xc", "pbe0")),
                "basis": method.get("basis", "def2-svp"),
                **electronic_method_config(method),
            },
            charge=state["charge"],
            multiplicity=state["multiplicity"],
        )
        # Compatibility spelling retained at the artifact boundary.
        method_record["disp_vdw"] = method_record.pop("dispersion", None)
        input_text = input_path.read_text(encoding="utf-8", errors="ignore")
        expected_dispersion_line = (
            f"disp vdw {int(method['disp_vdw'])}" if method.get("disp_vdw") is not None else None
        )
        dispersion_in_input = bool(
            expected_dispersion_line and expected_dispersion_line.lower() in input_text.lower()
        )
        run_token = fingerprint_dict(
            {
                "reaction": reaction.artifact_id,
                "engine": self.name,
                "input": existing_file_hash(input_path),
            }
        )
        path_art = Artifact(
            artifact_id="reaction_path_" + run_token,
            artifact_type="reaction_path",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            paths={"input": str(input_path)},
            method=method_record,
            data={
                "reaction_id": _reaction_id(reaction),
                "reaction_path": analyze_reaction_path(
                    reaction_id=_reaction_id(reaction),
                    engine=self.name,
                    comments=[],
                    converged=False,
                    energy_source="input_rendered_path_not_executed",
                ).model_dump(),
                "resolved_charge": state["charge"],
                "resolved_multiplicity": state["multiplicity"],
                "electron_count": state["electron_count"],
                "path_ensemble": method.get("path_ensemble"),
                "path_initialization": method.get("path_initialization"),
            },
            qc={
                "input_rendered": True,
                "neb_profile_classification": "missing_profile",
                "dispersion_in_input": dispersion_in_input,
                "real_neb_executed": False,
                "fallback_dummy": False,
            },
            provenance={
                "created_by": "NWChemNEBEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "input_sha256": existing_file_hash(input_path),
                "path_ensemble": method.get("path_ensemble"),
                "path_initialization": method.get("path_initialization"),
            },
        )
        executable = resolve_executable("nwchem", method.get("executable"))
        if (
            method.get("dry_run", False)
            or not allow_nwchem_subprocess(method)
            or executable is None
        ):
            reason = (
                "dry_run=True"
                if method.get("dry_run", False)
                else (
                    "NWChem execution disabled"
                    if not allow_nwchem_subprocess(method)
                    else "NWChem executable not found"
                )
            )
            fail = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                reason,
                category="nwchem_neb_not_run",
                parents=[path_art.artifact_id],
                input=str(input_path),
                recommended_fallback="configure NWChem and set allow_subprocess=true",
            )
            return TSSearchResult([path_art, fail], record=fail.model_dump(), success=False)
        try:
            result, output = run_nwchem_input(input_path, method, "nwchem_neb.out")
        except (OSError, TypeError, ValueError) as exc:
            fail = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                str(exc),
                category="nwchem_neb_failed",
                parents=[path_art.artifact_id],
            )
            return TSSearchResult([path_art, fail], record=fail.model_dump(), success=False)
        neb_path = self._path_file(input_path.parent)
        converged = nwchem_neb_converged(output)
        parsed_neb = parse_nwchem_output(output)
        method_evidence = assess_nwchem_path_evidence(
            parsed_neb,
            method,
            dispersion_in_input=dispersion_in_input,
            rendered_input=input_text,
        )
        required_version = method.get("required_program_version")
        program_version_ok = method_evidence["program_version_ok"]
        dispersion_evidence_ok = method_evidence["dispersion_evidence_ok"]
        raw_method_evidence_ok = method_evidence["validated"]
        method_identity_ok = method_evidence["method_identity_validated"]
        path_art.paths.update(
            {
                "output": result.stdout_path,
                "stderr": result.stderr_path,
                "neb_path_xyz": str(neb_path) if neb_path else "",
            }
        )
        path_art.data.update(
            {
                "program_version": parsed_neb.get("program_version"),
                "dft_d3_applied": parsed_neb.get("dft_d3_applied"),
                "dispersion_correction_hartree": parsed_neb.get("dispersion_correction_hartree"),
            }
        )
        path_art.qc.update(
            {
                "real_neb_executed": True,
                "neb_converged": converged,
                "command_ok": result.ok,
                "normal_termination": parsed_neb.get("normal_termination"),
                "program_version_ok": program_version_ok,
                "dispersion_applied": parsed_neb.get("dft_d3_applied"),
                "method_evidence_validated": raw_method_evidence_ok,
                "method_identity_validated": method_identity_ok,
                "input_method_evidence_validated": method_evidence[
                    "input_method_evidence_ok"
                ],
            }
        )
        path_art.provenance.update(
            {
                "command": asdict(result),
                "output_sha256": existing_file_hash(result.stdout_path),
                "stderr_sha256": existing_file_hash(result.stderr_path),
                "neb_path_sha256": existing_file_hash(neb_path),
                "observed_nwchem_method": method_evidence["observed_method"],
            }
        )
        path_geometry_error: str | None = None
        if neb_path is not None:
            try:
                validate_path_trajectory(
                    read_xyz_trajectory(neb_path),
                    start,
                    end,
                    endpoint_tolerance_A=float(
                        method.get("path_endpoint_rmsd_tolerance_A", 0.10)
                    ),
                )
            except (OSError, TypeError, ValueError) as exc:
                path_geometry_error = str(exc)
        path_art.qc["path_geometry_validated"] = bool(
            neb_path is not None and path_geometry_error is None
        )
        if path_geometry_error:
            path_art.qc["path_geometry_error"] = path_geometry_error
        if neb_path is not None and path_geometry_error is None:
            profile = classify_neb_energy_profile(
                neb_path,
                float(
                    method.get(
                        "path_barrier_threshold_kcal_mol",
                        method.get("neb_barrier_threshold_kcal_mol", 0.5),
                    )
                ),
                reaction_id=_reaction_id(reaction),
                engine=self.name,
                converged=converged,
            )
            path_art.data.update(profile)
            path_art.qc["neb_profile_classification"] = profile[
                "neb_profile_classification"
            ]
            path_record = profile["reaction_path"]
        else:
            path_record = path_art.data["reaction_path"]
        optimization_files = list(input_path.parent.glob("*.neb_epath"))
        optimization_path = (
            max(optimization_files, key=lambda item: item.stat().st_mtime)
            if optimization_files
            else None
        )
        optimization = parse_neb_optimization_history(
            optimization_path,
            converged=converged,
            stagnation_window=int(method.get("path_stagnation_window", 3)),
            stagnation_relative_range=float(
                method.get("path_stagnation_relative_range", 0.10)
            ),
            stagnation_displacement_A=float(
                method.get("path_stagnation_displacement_A", 1.0e-4)
            ),
        )
        attempt_record = make_path_attempt_record(
            attempt_id="path_attempt_"
            + fingerprint_dict(
                {
                    "reaction": reaction.artifact_id,
                    "engine": self.name,
                    "strategy": method.get("path_strategy", "double_ended_path"),
                    "input": existing_file_hash(input_path),
                }
            ),
            reaction_id=_reaction_id(reaction),
            strategy=str(method.get("path_strategy", "double_ended_path")),
            engine=self.name,
            endpoint_basin_status=str(
                method.get("endpoint_basin_status", "unresolved")
            ),
            path=path_record,
            optimization=optimization,
            evidence={
                "reaction_case_artifact_id": method.get(
                    "reaction_case_artifact_id"
                ),
                "path_xyz": str(neb_path) if neb_path else None,
                "optimization_history_path": (
                    str(optimization_path) if optimization_path else None
                ),
                "path_ensemble": method.get("path_ensemble"),
                "path_initialization": method.get("path_initialization"),
            },
        )
        attempt_artifact = Artifact(
            artifact_id=attempt_record.attempt_id,
            artifact_type="path_attempt",
            parents=[path_art.artifact_id],
            paths={
                "path_xyz": str(neb_path) if neb_path else "",
                "optimization_history": (
                    str(optimization_path) if optimization_path else ""
                ),
            },
            data=attempt_record.model_dump(mode="json"),
            method={**method_record, "task": "reaction_path"},
            qc={
                "path_converged": converged,
                "path_stagnant": optimization.stagnant,
                "saddle_refinement_allowed": (
                    attempt_record.next_action == "refine_saddle"
                ),
            },
            provenance={
                "created_by": "NWChemNEBEngine",
                "reaction_path_artifact_id": path_art.artifact_id,
            },
        )
        seed_only = bool(
            not converged
            and method.get("allow_unconverged_path_seed", True)
            and (
                result.ok
                or method.get("allow_interrupted_path_seed", False)
            )
            and neb_path is not None
            and path_geometry_error is None
            and method_identity_ok
            and attempt_record.next_action == "refine_saddle"
        )
        seed_method_evidence_sufficient = bool(
            raw_method_evidence_ok or seed_only
        )
        path_art.qc.update(
            {
                "path_energy_publishable": converged,
                "unconverged_path_used_as_seed_only": seed_only,
                "saddle_seed_method_evidence_sufficient": (
                    seed_method_evidence_sufficient
                ),
            }
        )
        attempt_artifact.qc.update(
            {
                "saddle_refinement_allowed": (
                    attempt_record.next_action == "refine_saddle"
                    and (converged or seed_only)
                ),
                "path_energy_publishable": converged,
                "unconverged_path_used_as_seed_only": seed_only,
            }
        )
        path_execution_usable = bool(
            (result.ok and converged and raw_method_evidence_ok) or seed_only
        )
        if (
            not path_execution_usable
            or neb_path is None
            or path_geometry_error is not None
            or not seed_method_evidence_sufficient
        ):
            category = (
                "nwchem_neb_path_geometry_invalid"
                if path_geometry_error is not None
                else (
                    "nwchem_neb_method_evidence_invalid"
                    if result.ok and neb_path is not None and not raw_method_evidence_ok
                    else "nwchem_neb_failed"
                )
            )
            fail = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                "NWChem NEB failed validation: "
                f"returncode={result.returncode}, converged={converged}, "
                f"path_found={neb_path is not None}, version_ok={program_version_ok}, "
                f"d3_evidence_ok={dispersion_evidence_ok}",
                category=category,
                parents=[path_art.artifact_id, attempt_artifact.artifact_id],
                output=result.stdout_path,
                path_diagnosis=attempt_record.diagnosis,
                next_action=attempt_record.next_action,
            )
            return TSSearchResult(
                [path_art, attempt_artifact, fail],
                record=fail.model_dump(),
                success=False,
            )
        profile_classification = profile["neb_profile_classification"]
        if attempt_record.next_action != "refine_saddle":
            monotonic = profile_classification == "monotonic_no_internal_maximum"
            assessment = Artifact.failure(
                f"ts_unresolved_{_reaction_id(reaction)}",
                "ts_result",
                "NWChem path is not authorized for saddle refinement",
                category=(
                    "endpoint_path_inconsistency"
                    if attempt_record.diagnosis == "endpoint_path_inconsistency"
                    else "same_basin"
                    if attempt_record.diagnosis == "same_basin"
                    else "neb_no_internal_maximum"
                    if monotonic
                    else "neb_internal_maximum_unresolved"
                ),
                parents=[path_art.artifact_id, attempt_artifact.artifact_id],
                recoverable=True,
                recommended_fallback=(
                    "follow the path-attempt next_action; do not repeat an identical "
                    "fixed-image path or force a saddle"
                ),
                data={
                    "reaction_id": _reaction_id(reaction),
                    "path_diagnosis": attempt_record.diagnosis,
                    "next_action": attempt_record.next_action,
                    **profile,
                },
            )
            assessment.qc = {
                "possible_barrierless_path": bool(
                    monotonic
                    and attempt_record.endpoint_basin_status != "distinct_basin"
                ),
                "ts_validated_by_frequency": False,
                "fallback_dummy": False,
            }
            return TSSearchResult(
                [path_art, attempt_artifact, assessment],
                record=assessment.model_dump(),
                success=False,
            )
        try:
            guess, bead_index, selection = select_neb_ts_guess(neb_path)
        except (IndexError, OSError, TypeError, ValueError) as exc:
            fail = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                str(exc),
                category="nwchem_neb_path_parse_failed",
                parents=[path_art.artifact_id, attempt_artifact.artifact_id],
            )
            return TSSearchResult(
                [path_art, attempt_artifact, fail],
                record=fail.model_dump(),
                success=False,
            )
        ts_id = _ts_id(reaction)
        guess_path = wd / f"{ts_id}_neb_guess.xyz"
        write_xyz(guess, guess_path)
        ts_species = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            guess_path,
            ts_id,
            self.name,
            {
                "real_neb_executed": True,
                "neb_converged": converged,
                "unconverged_path_used_as_seed_only": seed_only,
                "neb_bead_index": bead_index,
                "neb_selection": selection,
            },
        )
        try:
            ts_species.data["reaction_mode_reference"] = reaction_mode_reference(
                read_xyz_trajectory(neb_path),
                bead_index,
                source="converged_nwchem_neb",
            )
        except ValueError:
            # A degenerate test/path tangent supplies no scientific evidence;
            # the final mode gate must then rely on an explicit coordinate.
            pass
        ts_species.method = {**method_record, "task": "ts_guess"}
        ts_species.data.update(
            {
                "resolved_charge": state["charge"],
                "resolved_multiplicity": state["multiplicity"],
                "electron_count": state["electron_count"],
                "program_version": parsed_neb.get("program_version"),
            }
        )
        ts_species.provenance.update(
            {
                "created_by": "NWChemNEBEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "neb_path_artifact_id": path_art.artifact_id,
                "neb_path_sha256": existing_file_hash(neb_path),
            }
        )
        calc = NWChemEngine().saddle_frequency(
            ts_species,
            {
                **method,
                "method_id": method.get(
                    "saddle_method_id", method.get("method_id", "nwchem_saddle")
                ),
            },
            str(wd / "saddle_freq"),
        )
        calc.method = {
            **method_record,
            **(calc.method or {}),
            "engine": "nwchem",
            "backend": self.name,
            "stage": "ts-search",
            "backend_step": "saddle_freq",
        }
        calc.provenance.update(
            {
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "neb_path_artifact_id": path_art.artifact_id,
                "neb_method_evidence_validated": raw_method_evidence_ok,
                "neb_method_identity_validated": method_identity_ok,
                "neb_seed_method_evidence_sufficient": (
                    seed_method_evidence_sufficient
                ),
            }
        )
        ts_assessment = validate_ts_frequency_calculation(
            ts_species,
            calc,
            method,
            backend=self.name,
            search_method_evidence_validated=seed_method_evidence_sufficient,
            overlap_evaluator=estimate_reaction_mode_overlap,
        )
        overlap = ts_assessment["mode_overlap_score"]
        validated = ts_assessment["validated"]
        ts_species.provenance["saddle_frequency_calculation_id"] = calc.artifact_id
        ts_species.provenance["saddle_frequency_output_sha256"] = existing_file_hash(
            calc.paths.get("output")
        )
        reaction_validated = Artifact(
            artifact_id=reaction.artifact_id + "_with_ts",
            artifact_type="reaction_validated",
            parents=[
                reaction.artifact_id,
                path_art.artifact_id,
                ts_species.artifact_id,
                calc.artifact_id,
            ],
            data={
                **reaction.data,
                "ts_species_id": ts_species.artifact_id,
                "ts_calc_id": calc.artifact_id,
                "ts_backend": self.name,
                "resolved_charge": state["charge"],
                "resolved_multiplicity": state["multiplicity"],
                "electron_count": state["electron_count"],
                "program_version": calc.data.get("program_version"),
            },
            qc={
                "ts_engine": self.name,
                "ts_found": calc.status.status == "success",
                "ts_validated_by_frequency": validated,
                "n_imag": calc.data.get("n_imag"),
                "imag_freq_cm1": calc.data.get("imag_freq_cm1"),
                "mode_overlap_score": overlap,
                "real_neb_executed": True,
                "neb_converged": converged,
                "unconverged_path_used_as_seed_only": seed_only,
                "real_ts_search_executed": True,
                "real_qm_executed": calc.qc.get("real_qm_executed") is True,
                "method_evidence_validated": validated,
                "dispersion_applied": calc.qc.get("dispersion_applied"),
                "fallback_dummy": bool(calc.qc.get("fallback_dummy", False)),
            },
            method={**method_record, "task": "ts_search"},
            provenance={
                "created_by": "NWChemNEBEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "neb_path_artifact_id": path_art.artifact_id,
                "saddle_frequency_calculation_id": calc.artifact_id,
                "method_evidence": {
                    "neb": raw_method_evidence_ok,
                    "neb_method_identity": method_identity_ok,
                    "neb_seed_sufficient": seed_method_evidence_sufficient,
                    "saddle_frequency": bool(
                        calc.status.status == "success"
                        and str(calc.data.get("program_version")) == str(required_version)
                        and calc.qc.get("dispersion_applied") is True
                    ),
                },
            },
        )
        return TSSearchResult(
            [path_art, attempt_artifact, ts_species, calc, reaction_validated],
            record=reaction_validated.data,
            success=validated,
            ts_species_id=ts_id,
            ts_calc_id=calc.artifact_id,
        )

    def run_irc(
        self,
        reaction: Artifact,
        ts_species: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict[str, Any],
        workdir: str | Path,
    ) -> IRCResult:
        """Fail explicitly: NWChem NEB validates a path but is not an IRC."""
        Path(workdir).mkdir(parents=True, exist_ok=True)
        fail = Artifact.failure(
            f"irc_failed_{reaction.artifact_id}",
            "irc",
            "NWChem NEB is not an intrinsic reaction coordinate calculation",
            category="nwchem_irc_not_supported",
            parents=[reaction.artifact_id, ts_species.artifact_id],
            recoverable=True,
            recommended_fallback="use a gradient-capable IRC backend and keep the NWChem NEB result as path evidence",
            reaction_id=_reaction_id(reaction),
            real_irc_executed=False,
        )
        fail.qc = {"real_irc_executed": False, "irc_validated": False, "fallback_dummy": False}
        fail.method = {"engine": self.name, "task": "irc", "supported": False}
        return IRCResult([fail], record=fail.model_dump(), success=False)
