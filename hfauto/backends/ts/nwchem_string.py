"""NWChem zero-temperature-string backend for adaptive molecular paths.

The backend owns path optimization only.  A converged internal maximum is
published as an evidenced saddle seed and is refined by ``NWChemSaddleEngine``;
string convergence is never promoted directly to a transition state.
"""

from __future__ import annotations

import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

from hfauto.backends.qm.nwchem import (
    NWChemInputRenderer,
    allow_nwchem_subprocess,
    parse_nwchem_output,
    run_nwchem_input,
)
from hfauto.backends.ts.base import TSSearchResult
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
    parse_string_optimization_history,
)
from hfauto.chemistry.reaction_profile import (
    analyze_reaction_path,
)
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.chemistry.xyz_trajectory import (
    reaction_mode_reference,
    read_xyz_trajectory,
)
from hfauto.core.artifacts import canonical_species_id
from hfauto.core.executables import resolve_executable
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact

_SUPPORTED_STRATEGIES = {
    "adaptive_double_ended_path",
    "reparameterized_double_ended_path",
}

_STRING_ENERGY_ROW = re.compile(
    r"^\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$", re.MULTILINE
)


def nwchem_string_converged(output: str) -> bool:
    """Require an affirmative NWChem ZTS convergence record."""

    return bool(
        re.search(
            r"^\s*@zts\s+(?:the\s+)?string calculation converged\s*$",
            output,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    ) and not bool(
        re.search(
            r"^\s*@zts\s+.*failed to converge\s*$",
            output,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    )


def parse_string_path_energies(
    source: str | Path | None, expected_images: int
) -> list[float]:
    """Read the final NWChem string energy block without guessing XYZ comments."""

    if source is None or not Path(source).is_file():
        return []
    text = Path(source).read_text(encoding="utf-8", errors="ignore")
    rows = [float(energy) for _fraction, energy in _STRING_ENERGY_ROW.findall(text)]
    count = int(expected_images)
    return rows[-count:] if count > 0 and len(rows) >= count else []


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


class NWChemStringEngine:
    """Optimize a reparameterized minimum-energy path and publish one seed."""

    name = "nwchem_string"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        return {**self.config, **(method or {})}

    @staticmethod
    def _path_file(workdir: Path) -> Path | None:
        candidates = [
            *workdir.glob("*.string_final.xyz"),
            *workdir.glob("*.stringpath*.xyz"),
            *workdir.glob("*string*path*.xyz"),
        ]
        finals = [path for path in candidates if "final" in path.name.lower()]
        pool = finals or candidates
        return max(pool, key=lambda path: path.stat().st_mtime) if pool else None

    @staticmethod
    def _energy_file(workdir: Path) -> Path | None:
        finals = list(workdir.glob("*.string_final_epath"))
        candidates = finals or list(workdir.glob("*.string_epath"))
        return (
            max(candidates, key=lambda path: path.stat().st_mtime)
            if candidates
            else None
        )

    def render_string_input(
        self,
        reactant: Artifact,
        product: Artifact,
        method_in: dict[str, Any],
        path: str | Path,
    ) -> Path:
        method = self._method(method_in)
        start, end, state = resolve_path_endpoints(reactant, product, method)
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        nbeads = max(
            3,
            int(
                method.get(
                    "path_image_count",
                    method.get("string_nbeads", method.get("nbeads", 9)),
                )
            ),
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

        seed = dict(method.get("adaptive_path_seed_candidate") or {})
        seed_path = Path(str(seed.get("xyz_path") or ""))
        middle = read_xyz(seed_path) if not initial_path and seed_path.is_file() else None
        if middle is not None and middle.symbols != start.symbols:
            raise ValueError("adaptive string seed atom symbols/order mismatch")

        basis = str(method.get("basis", "def2-svp"))
        lines = [
            "start hfauto_string",
            (
                f'title "hfauto string {canonical_species_id(reactant)} '
                f'to {canonical_species_id(product)}"'
            ),
            f"memory total {int(method.get('memory_mb', 2000))} mb",
            "",
            *NWChemInputRenderer._geometry("", start),
        ]
        if middle is not None:
            lines += [*NWChemInputRenderer._geometry("midgeom", middle)]
        lines += [
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
            "string",
            f"  nbeads {nbeads}",
            f"  maxiter {int(method.get('string_maxiter', 100))}",
            f"  stepsize {float(method.get('string_stepsize', 0.1))}",
            f"  nhist {int(method.get('string_nhist', 10))}",
            f"  interpol {int(method.get('string_interpol', 3))}",
            f"  tol {float(method.get('string_tolerance', 0.00045))}",
            "  freeze1 .true.",
            "  freezeN .true.",
            "  impose",
        ]
        if initial_path is not None:
            lines.append(f"  xyz_path {initial_path.name}")
        elif middle is not None:
            lines.append("  hasmiddle")
        lines += [
            f"  print_shift {int(method.get('string_print_shift', 1))}",
            "end",
            "",
            "task dft string ignore",
        ]
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return output

    def search_ts(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method_in: dict[str, Any],
        workdir: str | Path,
    ) -> TSSearchResult:
        method = self._method(method_in)
        strategy = str(method.get("path_strategy") or "")
        if strategy not in _SUPPORTED_STRATEGIES:
            failure = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                "NWChem string was not authorized for the reaction case",
                category="string_strategy_not_authorized",
                parents=[reaction.artifact_id],
            )
            return TSSearchResult([failure], record=failure.model_dump(), success=False)

        wd = Path(workdir)
        try:
            start, end, state = resolve_path_endpoints(reactant, product, method)
            nbeads = max(
                3,
                int(
                    method.get(
                        "path_image_count",
                        method.get(
                            "string_nbeads", method.get("nbeads", 9)
                        ),
                    )
                ),
            )
            initial_path, initialization = prepare_initial_path(
                reaction,
                start,
                end,
                method,
                wd / "string" / "generated_initial_path.xyz",
                image_count=nbeads,
            )
            method["path_initialization"] = initialization
            if initial_path is not None:
                method["path_initial_trajectory"] = str(initial_path)
            input_path = self.render_string_input(
                reactant,
                product,
                method,
                wd / "string" / "nwchem_string.nw",
            )
        except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
            failure = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}",
                "ts_result",
                str(exc),
                category="invalid_string_path_input",
                parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            )
            return TSSearchResult([failure], record=failure.model_dump(), success=False)

        method_record = electronic_structure_method(
            engine="nwchem",
            backend=self.name,
            task="zero_temperature_string",
            config=electronic_method_config(method),
            charge=state["charge"],
            multiplicity=state["multiplicity"],
        )
        method_record["disp_vdw"] = method_record.pop("dispersion", None)
        input_text = input_path.read_text(encoding="utf-8", errors="ignore")
        dispersion_line = (
            f"disp vdw {int(method['disp_vdw'])}"
            if method.get("disp_vdw") is not None
            else None
        )
        dispersion_in_input = bool(
            dispersion_line and dispersion_line.lower() in input_text.lower()
        )
        run_token = fingerprint_dict(
            {
                "reaction": reaction.artifact_id,
                "engine": self.name,
                "strategy": strategy,
                "input": existing_file_hash(input_path),
            }
        )
        path_artifact = Artifact(
            artifact_id=f"reaction_path_{run_token}",
            artifact_type="reaction_path",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            paths={"input": str(input_path)},
            method=method_record,
            data={
                "reaction_id": _reaction_id(reaction),
                "path_method": "zero_temperature_string",
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
                "real_path_executed": False,
                "path_converged": False,
                "dispersion_in_input": dispersion_in_input,
                "fallback_dummy": False,
            },
            provenance={
                "created_by": "NWChemStringEngine",
                "input_sha256": existing_file_hash(input_path),
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
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
                else "NWChem execution disabled"
                if not allow_nwchem_subprocess(method)
                else "NWChem executable not found"
            )
            failure = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}_{run_token}",
                "ts_result",
                reason,
                category="nwchem_string_not_run",
                parents=[path_artifact.artifact_id],
                recommended_fallback="configure NWChem and set allow_subprocess=true",
            )
            return TSSearchResult(
                [path_artifact, failure], record=failure.model_dump(), success=False
            )

        try:
            result, output = run_nwchem_input(
                input_path, method, "nwchem_string.out"
            )
        except (OSError, TypeError, ValueError) as exc:
            failure = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}_{run_token}",
                "ts_result",
                str(exc),
                category="nwchem_string_failed",
                parents=[path_artifact.artifact_id],
            )
            return TSSearchResult(
                [path_artifact, failure], record=failure.model_dump(), success=False
            )

        path_file = self._path_file(input_path.parent)
        energy_file = self._energy_file(input_path.parent)
        converged = nwchem_string_converged(output)
        parsed = parse_nwchem_output(output)
        method_evidence = assess_nwchem_path_evidence(
            parsed,
            method,
            dispersion_in_input=dispersion_in_input,
            rendered_input=input_text,
        )
        version_ok = method_evidence["program_version_ok"]
        dispersion_ok = method_evidence["dispersion_evidence_ok"]
        method_evidence_ok = method_evidence["validated"]
        method_identity_ok = method_evidence["method_identity_validated"]
        path_artifact.paths.update(
            {
                "output": result.stdout_path,
                "stderr": result.stderr_path,
                "path_xyz": str(path_file) if path_file else "",
                "path_energies": str(energy_file) if energy_file else "",
            }
        )
        path_artifact.data.update(
            {
                "program_version": parsed.get("program_version"),
                "dft_d3_applied": parsed.get("dft_d3_applied"),
                "dispersion_correction_hartree": parsed.get(
                    "dispersion_correction_hartree"
                ),
            }
        )
        path_artifact.qc.update(
            {
                "real_path_executed": True,
                "path_converged": converged,
                "command_ok": result.ok,
                "normal_termination": parsed.get("normal_termination"),
                "program_version_ok": version_ok,
                "dispersion_applied": parsed.get("dft_d3_applied"),
                "method_evidence_validated": method_evidence_ok,
                "method_identity_validated": method_identity_ok,
                "input_method_evidence_validated": method_evidence[
                    "input_method_evidence_ok"
                ],
            }
        )
        path_artifact.provenance.update(
            {
                "command": asdict(result),
                "output_sha256": existing_file_hash(result.stdout_path),
                "stderr_sha256": existing_file_hash(result.stderr_path),
                "path_sha256": existing_file_hash(path_file),
                "path_energies_sha256": existing_file_hash(energy_file),
                "observed_nwchem_method": method_evidence["observed_method"],
            }
        )
        images: list[XYZ] = []
        path_geometry_error: str | None = None
        if path_file:
            try:
                images = read_xyz_trajectory(path_file)
                validate_path_trajectory(
                    images,
                    start,
                    end,
                    endpoint_tolerance_A=float(
                        method.get("path_endpoint_rmsd_tolerance_A", 0.10)
                    ),
                )
            except (OSError, TypeError, ValueError) as exc:
                path_geometry_error = str(exc)
        path_artifact.qc["path_geometry_validated"] = bool(
            path_file is not None and path_geometry_error is None
        )
        if path_geometry_error:
            path_artifact.qc["path_geometry_error"] = path_geometry_error
        energies = parse_string_path_energies(energy_file, len(images))
        comments = (
            [f"energy_hartree={energy:.15f}" for energy in energies]
            if energies
            else [image.comment for image in images]
        )
        path_record = analyze_reaction_path(
            reaction_id=_reaction_id(reaction),
            engine=self.name,
            comments=comments,
            converged=converged,
            barrier_threshold_kcal_mol=float(
                method.get("path_barrier_threshold_kcal_mol", 0.5)
            ),
            energy_source=(
                "nwchem_string_final_epath"
                if energies
                else "xyz_comment"
            ),
        )
        path_artifact.data["reaction_path"] = path_record.model_dump()
        optimization = parse_string_optimization_history(
            output, converged=converged
        )
        attempt_record = make_path_attempt_record(
            attempt_id="path_attempt_"
            + fingerprint_dict(
                {
                    "reaction": reaction.artifact_id,
                    "engine": self.name,
                    "strategy": strategy,
                    "input": existing_file_hash(input_path),
                }
            ),
            reaction_id=_reaction_id(reaction),
            strategy=strategy,
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
                "path_xyz": str(path_file) if path_file else None,
                "initial_path_xyz": method.get("initial_path_xyz"),
                "adaptive_path_seed_candidate": method.get(
                    "adaptive_path_seed_candidate"
                ),
                "path_ensemble": method.get("path_ensemble"),
                "path_initialization": method.get("path_initialization"),
            },
        )
        attempt = Artifact(
            artifact_id=attempt_record.attempt_id,
            artifact_type="path_attempt",
            parents=[path_artifact.artifact_id],
            paths={"path_xyz": str(path_file) if path_file else ""},
            data=attempt_record.model_dump(mode="json"),
            method={**method_record, "task": "reaction_path"},
            qc={
                "path_converged": converged,
                "saddle_refinement_allowed": (
                    attempt_record.next_action == "refine_saddle"
                ),
            },
            provenance={
                "created_by": "NWChemStringEngine",
                "path_ensemble": method.get("path_ensemble"),
            },
        )
        seed_only = bool(
            not converged
            and method.get("allow_unconverged_path_seed", True)
            and (
                result.ok
                or method.get("allow_interrupted_path_seed", False)
            )
            and path_file is not None
            and path_geometry_error is None
            and method_identity_ok
            and attempt_record.next_action == "refine_saddle"
        )
        path_artifact.qc.update(
            {
                "path_energy_publishable": converged,
                "unconverged_path_used_as_seed_only": seed_only,
            }
        )
        attempt.qc.update(
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
            (result.ok and converged and method_evidence_ok) or seed_only
        )
        if (
            not path_execution_usable
            or path_file is None
            or path_geometry_error is not None
            or not (method_evidence_ok or seed_only)
        ):
            failure = Artifact.failure(
                f"ts_failed_{_reaction_id(reaction)}_{run_token}",
                "ts_result",
                "NWChem string failed validation: "
                f"returncode={result.returncode}, converged={converged}, "
                f"path_found={path_file is not None}, version_ok={version_ok}, "
                f"d3_evidence_ok={dispersion_ok}",
                category=(
                    "nwchem_string_path_geometry_invalid"
                    if path_geometry_error is not None
                    else (
                        "nwchem_string_method_evidence_invalid"
                        if result.ok and path_file and not method_evidence_ok
                        else "nwchem_string_failed"
                    )
                ),
                parents=[path_artifact.artifact_id, attempt.artifact_id],
                path_diagnosis=attempt_record.diagnosis,
                next_action=attempt_record.next_action,
            )
            return TSSearchResult(
                [path_artifact, attempt, failure],
                record=failure.model_dump(),
                success=False,
            )
        if attempt_record.next_action != "refine_saddle":
            failure = Artifact.failure(
                f"ts_unresolved_{_reaction_id(reaction)}_{run_token}",
                "ts_result",
                "Converged string did not resolve exactly one internal saddle candidate",
                category=(
                    "multiple_step_candidate"
                    if attempt_record.diagnosis == "multiple_step_candidate"
                    else "string_no_internal_maximum"
                ),
                parents=[path_artifact.artifact_id, attempt.artifact_id],
                path_diagnosis=attempt_record.diagnosis,
                next_action=attempt_record.next_action,
            )
            return TSSearchResult(
                [path_artifact, attempt, failure],
                record=failure.model_dump(),
                success=False,
            )

        image_index = path_record.highest_internal_image_index
        if image_index is None or not 0 < image_index < len(images) - 1:
            failure = Artifact.failure(
                f"ts_unresolved_{_reaction_id(reaction)}_{run_token}",
                "ts_result",
                "Resolved path record has no valid internal maximum index",
                category="string_internal_maximum_index_invalid",
                parents=[path_artifact.artifact_id, attempt.artifact_id],
            )
            return TSSearchResult(
                [path_artifact, attempt, failure],
                record=failure.model_dump(),
                success=False,
            )
        guess = images[image_index]
        selection = "highest_energy_internal_image"
        seed_id = "saddle_seed_" + fingerprint_dict(
            {
                "reaction": reaction.artifact_id,
                "path": existing_file_hash(path_file),
                "image": image_index,
            }
        )
        seed_path = wd / f"{seed_id}.xyz"
        write_xyz(guess, seed_path)
        seed = Artifact(
            artifact_id=seed_id,
            artifact_type="species",
            parents=[
                reaction.artifact_id,
                reactant.artifact_id,
                product.artifact_id,
                path_artifact.artifact_id,
            ],
            paths={"xyz": str(seed_path)},
            data={
                **reactant.data,
                "species_id": seed_id,
                "state": "saddle_seed",
                "xyz_path": str(seed_path),
                "reaction_id": _reaction_id(reaction),
                "reaction_coordinate": reaction.data.get("reaction_coordinate"),
                "bond_changes": list(reaction.data.get("bond_changes") or []),
            },
            method={**method_record, "task": "saddle_seed"},
            qc={
                "scientific_role": "saddle_seed_only",
                "path_converged": converged,
                "unconverged_path_used_as_seed_only": seed_only,
                "path_image_index": image_index,
                "selection": selection,
                "method_evidence_validated": method_evidence_ok,
                "method_identity_validated": method_identity_ok,
                "saddle_seed_method_evidence_sufficient": bool(
                    method_evidence_ok or seed_only
                ),
            },
            provenance={"created_by": "NWChemStringEngine"},
        )
        seed_source = (
            "converged_zero_temperature_string"
            if converged
            else "unconverged_zero_temperature_string_seed_only"
        )
        try:
            mode_reference = reaction_mode_reference(
                images,
                image_index,
                source=seed_source,
            )
        except ValueError:
            mode_reference = None
        candidate = {
            "xyz_path": str(seed_path),
            "species_artifact_id": seed.artifact_id,
            "path_artifact_id": path_artifact.artifact_id,
            "path_image_index": image_index,
            "source": seed_source,
            "reaction_mode_reference": mode_reference,
        }
        saddle_attempt = Artifact(
            artifact_id="saddle_attempt_"
            + fingerprint_dict(
                {"reaction": reaction.artifact_id, "path": path_artifact.artifact_id}
            ),
            artifact_type="saddle_attempt",
            parents=[attempt.artifact_id, seed.artifact_id],
            paths={"seed_xyz": str(seed_path)},
            data={
                "reaction_id": _reaction_id(reaction),
                "strategy": strategy,
                "engine": self.name,
                "diagnosis": "resolved_saddle_candidate",
                "next_action": "refine_saddle",
                "candidate": candidate,
                "search_evidence_validated": True,
            },
            qc={
                "path_converged": converged,
                "unconverged_path_used_as_seed_only": seed_only,
                "saddle_seed_resolved": True,
                "method_evidence_validated": method_evidence_ok,
                "method_identity_validated": method_identity_ok,
                "saddle_seed_method_evidence_sufficient": bool(
                    method_evidence_ok or seed_only
                ),
            },
            provenance={"created_by": "NWChemStringEngine"},
        )
        pending = Artifact.failure(
            f"ts_pending_{_reaction_id(reaction)}_{run_token}",
            "ts_result",
            "String path resolved a seed; independent Hessian/saddle refinement is required",
            category="saddle_refinement_required",
            parents=[saddle_attempt.artifact_id],
            recoverable=True,
            candidate=candidate,
        )
        return TSSearchResult(
            [path_artifact, attempt, seed, saddle_attempt, pending],
            record=pending.model_dump(),
            success=False,
        )
