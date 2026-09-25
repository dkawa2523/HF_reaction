"""pysisyphus path adapter, currently used for NWChem-backed IRC.

The previous compatibility skeleton silently returned dummy IRC endpoints.  This
adapter either executes a real gradient/Hessian workflow or returns a structured
failure; it never upgrades synthetic endpoints to scientific evidence.
"""

from __future__ import annotations

import os
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from hfauto.backends.qm.nwchem import (
    NWChemEngine,
    nwchem_environment,
    parse_nwchem_output,
    species_xyz_path,
)
from hfauto.backends.ts.base import (
    IRCResult,
    TSSearchResult,
    make_ts_species_artifact,
    validate_ts_frequency_calculation,
)
from hfauto.backends.ts.nwchem_neb import NWChemNEBEngine
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.methods import electronic_structure_method
from hfauto.chemistry.reaction_path_qc import (
    endpoint_pair_match_qc,
    make_displaced_reaction_coordinate_seed,
    reaction_coordinate_vector,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.schemas.artifact import Artifact


def allow_pysisyphus_subprocess(method: dict[str, Any]) -> bool:
    return bool(
        method.get("allow_subprocess", False)
        or os.environ.get("HFAUTO_ALLOW_PYSISYPHUS") == "1"
        or os.environ.get("HFAUTO_ALLOW_SUBPROCESS") == "1"
    )


def _pysisyphus_environment(
    pysis: str, nwchem: str, method: dict[str, Any]
) -> dict[str, str]:
    """Build one resource-bounded QCEngine environment.

    QCEngine does not launch MPI on a single node by default.  When an NWChem
    installation supplies its matching ``prterun``, opt in explicitly so that
    ``ncores`` is effective.  Serial installations remain serial.
    """

    env = nwchem_environment(nwchem, method.get("env"))
    env["PATH"] = str(Path(pysis).parent) + os.pathsep + env.get(
        "PATH", os.environ.get("PATH", "")
    )
    ncores = max(1, int(method.get("ncores", 1)))
    env["QCENGINE_NCORES"] = str(ncores)
    env["QCENGINE_MEMORY"] = str(
        max(0.25, float(method.get("memory_mb", 2000)) / 1024.0)
    )
    mpi_launcher = Path(nwchem).resolve().with_name("prterun")
    if ncores > 1 and mpi_launcher.is_file():
        env["QCENGINE_USE_MPIEXEC"] = "true"
        env["QCENGINE_MPIEXEC_COMMAND"] = (
            f"{mpi_launcher} -np {{total_ranks}}"
        )
    return env


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


def _file_hash(path: str | Path | None) -> str | None:
    if path is None or not Path(path).is_file():
        return None
    return sha256_file(path)


def _method_record(
    method: dict[str, Any],
    state: dict[str, Any],
    *,
    task: str = "irc",
) -> dict[str, Any]:
    return {
        "engine": "nwchem",
        "backend": "pysisyphus",
        "program": "nwchem",
        "task": task,
        "functional": method.get("functional", method.get("xc", "pbe0")),
        "basis": method.get("basis", "def2-svp"),
        "disp_vdw": method.get("disp_vdw"),
        "required_program_version": method.get("required_program_version"),
        "charge": state["charge"],
        "multiplicity": state["multiplicity"],
    }


def _nwchem_qcengine_evidence(
    workdir: Path,
    required_version: str,
    *,
    expected_basis: str | None = None,
    require_spherical: bool = False,
) -> dict[str, Any]:
    """Audit every QCEngine/NWChem output against the requested PES.

    ``basis__spherical`` is an input contract, but scientific evidence must come
    from the program output.  Keeping the stricter checks opt-in preserves the
    helper's use for historical manifests that predate basis metadata.
    """

    output_paths = sorted(
        path for path in workdir.rglob("*qce_nwchem_stdout*") if path.is_file()
    )
    records: list[dict[str, Any]] = []
    for path in output_paths:
        text = path.read_text(encoding="utf-8", errors="ignore")
        parsed = parse_nwchem_output(text)
        basis_representations = sorted(
            {
                value.lower()
                for value in re.findall(
                    r'Basis\s+"ao basis".*?\((spherical|cartesian)\)',
                    text,
                    flags=re.IGNORECASE,
                )
            }
        )
        basis_name_matched = bool(
            expected_basis is None
            or str(expected_basis).strip().lower() in text.lower()
        )
        spherical_matched = bool(
            not require_spherical or basis_representations == ["spherical"]
        )
        records.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "program_version": parsed.get("program_version"),
                "normal_termination": parsed.get("normal_termination"),
                "dft_d3_applied": parsed.get("dft_d3_applied"),
                "dispersion_correction_hartree": parsed.get(
                    "dispersion_correction_hartree"
                ),
                "basis_representations": basis_representations,
                "basis_name_matched": basis_name_matched,
                "spherical_basis_matched": spherical_matched,
            }
        )
    accepted = bool(
        records
        and all(
            record["program_version"] == str(required_version)
            and record["normal_termination"] is True
            and record["dft_d3_applied"] is True
            and record["basis_name_matched"] is True
            and record["spherical_basis_matched"] is True
            for record in records
        )
    )
    return {
        "accepted": accepted,
        "required_program_version": str(required_version),
        "expected_basis": expected_basis,
        "require_spherical": require_spherical,
        "raw_output_count": len(records),
        "raw_outputs": records,
    }


def _qcengine_keywords(method: dict[str, Any]) -> dict[str, Any]:
    """Translate the shared NWChem PES definition to QCEngine keywords."""

    if method.get("disp_vdw") != 3:
        raise ValueError("pysisyphus/NWChem IRC requires disp_vdw=3 (D3(0))")
    keywords = dict(method.get("qcengine_keywords", {}) or {})
    requested = keywords.get("dft__disp")
    if requested is None:
        keywords["dft__disp"] = "vdw 3"
    elif str(requested).strip().lower() != "vdw 3":
        raise ValueError("qcengine_keywords.dft__disp conflicts with disp_vdw=3")
    spherical = keywords.get("basis__spherical")
    if spherical is None:
        keywords["basis__spherical"] = True
    elif spherical is not True:
        raise ValueError(
            "qcengine_keywords.basis__spherical must remain true to match endpoints"
        )
    return keywords


def _endpoint_optimization_evidence(workdir: Path) -> dict[str, Any]:
    """Read branch optimizer conclusions; file creation alone is not convergence."""

    branches: dict[str, dict[str, Any]] = {}
    for direction in ("forward", "backward"):
        log_path = workdir / f"{direction}_end_optimizer.log"
        geometry_path = workdir / f"{direction}_end_opt.xyz"
        converged = _optimizer_converged(log_path, geometry_path)
        branches[direction] = {
            "converged": converged,
            "geometry_path": str(geometry_path) if geometry_path.is_file() else None,
            "optimizer_log_path": str(log_path) if log_path.is_file() else None,
            "optimizer_log_sha256": _file_hash(log_path),
        }
    return {
        "accepted": all(branch["converged"] for branch in branches.values()),
        "backend": "pysisyphus",
        "branches": branches,
    }


def _optimize_irc_endpoints_with_nwchem(
    raw_endpoints: dict[str, Path],
    *,
    reaction_id: str,
    state: dict[str, Any],
    method: dict[str, Any],
    workdir: Path,
) -> tuple[dict[str, Any], list[Artifact]]:
    """Relax raw IRC termini with the shared, unconstrained QM optimizer."""

    calculations: list[Artifact] = []
    branches: dict[str, dict[str, Any]] = {}
    engine = NWChemEngine()
    for direction, raw_path in raw_endpoints.items():
        endpoint = Artifact(
            artifact_id=f"{reaction_id}_{direction}_irc_terminus",
            artifact_type="species",
            paths={"xyz": str(raw_path)},
            data={
                "species_id": f"{reaction_id}_{direction}_irc_terminus",
                "state": "irc_endpoint",
                "xyz_path": str(raw_path),
                "charge": state["charge"],
                "multiplicity": state["multiplicity"],
            },
            provenance={"created_by": "PysisyphusEngine"},
        )
        calculation = engine.optimize(
            endpoint,
            {
                **method,
                "method_id": method.get(
                    "method_id", "irc_endpoint_unconstrained_optimization"
                ),
                "allow_subprocess": True,
                "fallback_to_dummy": False,
            },
            str(workdir / f"{direction}_endpoint_optimization"),
        )
        calculations.append(calculation)
        final_path = calculation.paths.get("final_xyz")
        converged = bool(
            calculation.status.status == "success"
            and calculation.qc.get("real_qm_executed") is True
            and calculation.qc.get("geometry_converged") is True
            and calculation.qc.get("normal_termination") is True
            and final_path
            and Path(final_path).is_file()
        )
        branches[direction] = {
            "converged": converged,
            "raw_geometry_path": str(raw_path),
            "geometry_path": str(final_path) if converged else None,
            "calculation_id": calculation.artifact_id,
            "status": calculation.status.model_dump(),
        }
    return (
        {
            "accepted": all(branch["converged"] for branch in branches.values()),
            "backend": "nwchem",
            "branches": branches,
        },
        calculations,
    )


def _optimizer_converged(log_path: Path, geometry_path: Path) -> bool:
    try:
        log_text = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    return bool(
        geometry_path.is_file()
        and "Converged!" in log_text
        and "Number of cycles exceeded!" not in log_text
    )


def _dimer_convergence_evidence(workdir: Path) -> dict[str, Any]:
    optimizer_log = workdir / "ts_optimizer.log"
    geometry_path = workdir / "ts_opt.xyz"
    dimer_log = workdir / "dimer.log"
    try:
        text = dimer_log.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        text = ""
    final_rotation = text.rsplit("Doing dimer rotations", 1)[-1]
    curvatures = [
        float(value)
        for value in re.findall(r"\bC=\s*([-+0-9.Ee]+)", final_rotation)
    ]
    final_curvature = curvatures[-1] if curvatures else None
    optimizer_converged = _optimizer_converged(
        optimizer_log, geometry_path
    )
    rotation_converged = "Dimer rotation converged in" in final_rotation
    negative_curvature = (
        final_curvature is not None
        and final_curvature < 0.0
        and "Rotation did not yield a negative curvature" not in final_rotation
    )
    return {
        "accepted": bool(
            optimizer_converged and rotation_converged and negative_curvature
        ),
        "optimizer_converged": optimizer_converged,
        "rotation_converged": rotation_converged,
        "negative_curvature": negative_curvature,
        "final_curvature_hartree_per_bohr2": final_curvature,
        "optimizer_log_path": (
            str(optimizer_log) if optimizer_log.is_file() else None
        ),
        "optimizer_log_sha256": _file_hash(optimizer_log),
        "dimer_log_path": str(dimer_log) if dimer_log.is_file() else None,
        "dimer_log_sha256": _file_hash(dimer_log),
        "ts_guess_path": str(geometry_path) if geometry_path.is_file() else None,
    }


class PysisyphusEngine:
    """Run pysisyphus with QCEngine and NWChem gradients/Hessians."""

    name = "pysisyphus"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        return {**self.config, **(method or {})}

    def search_ts(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict[str, Any],
        workdir: str | Path,
    ) -> TSSearchResult:
        merged = self._method(method)
        if merged.get("path_strategy") == "endpoint_local_saddle_search":
            return self._search_ts_with_endpoint_dimers(
                reaction, reactant, product, merged, Path(workdir)
            )
        program = str(merged.get("program", merged.get("calculator_program", "nwchem"))).lower()
        if program == "nwchem":
            return NWChemNEBEngine(**self.config).search_ts(
                reaction, reactant, product, merged, workdir
            )
        fail = Artifact.failure(
            "ts_failed_" + _reaction_id(reaction),
            "ts_result",
            f"pysisyphus TS adapter does not support calculator program {program!r}",
            category="unsupported_calculator",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            recommended_fallback="use nwchem_neb or add a tested calculator adapter",
        )
        return TSSearchResult([fail], record=fail.model_dump(), success=False)

    def render_dimer_input(
        self,
        reaction: Artifact,
        endpoint: Artifact,
        method_in: dict[str, Any],
        path: str | Path,
        *,
        geometry_path: str | Path | None = None,
    ) -> Path:
        """Render one endpoint-local dimer search along declared chemistry."""

        method = self._method(method_in)
        program = str(
            method.get("program", method.get("calculator_program", "nwchem"))
        ).lower()
        if program != "nwchem":
            raise ValueError(
                "pysisyphus endpoint dimer currently supports only NWChem QCEngine"
            )
        if str(method.get("required_program_version") or "") != "7.2.3":
            raise ValueError(
                "pysisyphus/NWChem dimer requires required_program_version=7.2.3"
            )
        endpoint_xyz = species_xyz_path(endpoint).resolve()
        search_xyz = Path(geometry_path or endpoint_xyz).resolve()
        structure = read_xyz(search_xyz)
        state = resolve_electronic_state(
            endpoint.data, method, structure.symbols
        )
        orientation = reaction_coordinate_vector(reaction.data, search_xyz)
        if (
            orientation is None
            or orientation.shape != structure.coords.shape
            or not np.all(np.isfinite(orientation))
            or float(np.linalg.norm(orientation)) <= 1.0e-12
        ):
            raise ValueError(
                "declared reaction-coordinate vector is unavailable for endpoint dimer"
            )
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        orientation_path = output.parent / "reaction_coordinate_orientation.txt"
        np.savetxt(orientation_path, orientation.reshape(-1))
        config = {
            "geom": {"type": "cart", "fn": str(search_xyz)},
            "calc": {
                "type": "dimer",
                "N_raw": str(orientation_path.resolve()),
                "rotation_max_cycles": int(
                    method.get("dimer_rotation_maxiter", 15)
                ),
                "rotation_thresh": float(
                    method.get("dimer_rotation_threshold", 1.0e-4)
                ),
                "forward_hessian": False,
                "calc": {
                    "type": "qcengine",
                    "program": program,
                    "model": {
                        "method": str(
                            method.get("functional", method.get("xc", "pbe0"))
                        ),
                        "basis": str(method.get("basis", "def2-svp")),
                    },
                    "keywords": _qcengine_keywords(method),
                    "charge": state["charge"],
                    "mult": state["multiplicity"],
                    "pal": max(1, int(method.get("ncores", 1))),
                },
            },
            "tsopt": {
                "type": "plbfgs",
                "max_cycles": int(method.get("dimer_maxiter", 150)),
                "thresh": method.get("dimer_opt_threshold", "gau_loose"),
                "dump": True,
            },
        }
        output.write_text(
            yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
        )
        return output

    def _run_endpoint_dimer(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        endpoint: Artifact,
        endpoint_label: str,
        method: dict[str, Any],
        workdir: Path,
    ) -> TSSearchResult:
        workdir.mkdir(parents=True, exist_ok=True)
        try:
            structure = read_xyz(species_xyz_path(endpoint))
            state = resolve_electronic_state(
                endpoint.data, method, structure.symbols
            )
            target = product if endpoint_label == "reactant" else reactant
            seed = make_displaced_reaction_coordinate_seed(
                reaction.data,
                species_xyz_path(endpoint),
                species_xyz_path(target),
                workdir / "dimer_seed.xyz",
                displacement_A=float(
                    method.get("dimer_seed_displacement_A", 0.15)
                ),
            )
            input_path = self.render_dimer_input(
                reaction,
                endpoint,
                method,
                workdir / "pysis_dimer.yaml",
                geometry_path=seed["path"],
            )
        except (KeyError, OSError, TypeError, ValueError, yaml.YAMLError) as exc:
            fail = Artifact.failure(
                f"dimer_failed_{_reaction_id(reaction)}_{endpoint_label}",
                "ts_result",
                str(exc),
                category="pysisyphus_dimer_input_invalid",
                parents=[reaction.artifact_id, endpoint.artifact_id],
            )
            return TSSearchResult([fail], record=fail.model_dump(), success=False)

        pysis = resolve_executable("pysis", method.get("pysis_executable"))
        nwchem = resolve_executable("nwchem", method.get("executable"))
        runnable = allow_pysisyphus_subprocess(method) and not method.get(
            "dry_run", False
        )
        if not runnable or pysis is None or nwchem is None:
            reason = (
                "dry_run=True"
                if method.get("dry_run", False)
                else "pysisyphus execution disabled"
                if not allow_pysisyphus_subprocess(method)
                else "pysis executable not found"
                if pysis is None
                else "NWChem executable not found"
            )
            fail = Artifact.failure(
                f"dimer_failed_{_reaction_id(reaction)}_{endpoint_label}",
                "ts_result",
                reason,
                category="pysisyphus_dimer_not_run",
                parents=[reaction.artifact_id, endpoint.artifact_id],
                input=str(input_path),
            )
            return TSSearchResult([fail], record=fail.model_dump(), success=False)

        env = _pysisyphus_environment(pysis, nwchem, method)
        result = run_command(
            [pysis, input_path.name],
            cwd=workdir,
            timeout_s=int(method.get("dimer_timeout_s", method.get("timeout_s", 172_800))),
            env=env,
            stdout_name="pysis_dimer.out",
            stderr_name="pysis_dimer.err",
        )
        ts_guess_path = workdir / "ts_opt.xyz"
        dimer_evidence = _dimer_convergence_evidence(workdir)
        optimizer_log = workdir / "ts_optimizer.log"
        dimer_converged = bool(dimer_evidence["accepted"])
        runtime_evidence = _nwchem_qcengine_evidence(
            workdir / "qm_calcs",
            str(method["required_program_version"]),
            expected_basis=str(method.get("basis", "def2-svp")),
            require_spherical=True,
        )
        attempt_id = "saddle_attempt_" + fingerprint_dict(
            {
                "reaction": reaction.artifact_id,
                "strategy": "endpoint_local_saddle_search",
                "endpoint": endpoint.artifact_id,
                "input": _file_hash(input_path),
            }
        )
        attempt = Artifact(
            artifact_id=attempt_id,
            artifact_type="saddle_attempt",
            parents=[reaction.artifact_id, endpoint.artifact_id],
            paths={
                "input": str(input_path),
                "output": result.stdout_path,
                "stderr": result.stderr_path,
                "optimizer_log": str(optimizer_log) if optimizer_log.is_file() else "",
                "ts_guess_xyz": str(ts_guess_path) if ts_guess_path.is_file() else "",
                "dimer_seed_xyz": seed["path"],
            },
            data={
                "reaction_id": _reaction_id(reaction),
                "strategy": "endpoint_local_saddle_search",
                "engine": self.name,
                "endpoint_label": endpoint_label,
                "endpoint_artifact_id": endpoint.artifact_id,
                "dimer_seed_displacement_A": seed["displacement_A"],
                "diagnosis": "saddle_candidate" if dimer_converged else "saddle_optimizer_failed",
                "next_action": "validate_saddle_frequency" if dimer_converged else "adaptive_path_refinement",
            },
            method=_method_record(method, state, task="endpoint_dimer"),
            qc={
                "real_ts_search_executed": True,
                "command_ok": result.ok,
                "optimizer_converged": dimer_evidence["optimizer_converged"],
                "dimer_rotation_converged": dimer_evidence[
                    "rotation_converged"
                ],
                "negative_curvature": dimer_evidence["negative_curvature"],
                "dimer_candidate_accepted": dimer_converged,
                "method_evidence_validated": runtime_evidence["accepted"],
                "fallback_dummy": False,
            },
            provenance={
                "created_by": "PysisyphusEngine",
                "command": asdict(result),
                "input_sha256": _file_hash(input_path),
                "dimer_seed": seed,
                "optimizer_log_sha256": _file_hash(optimizer_log),
                "nwchem_qcengine_evidence": runtime_evidence,
                "dimer_convergence_evidence": dimer_evidence,
            },
        )
        if (
            not result.ok
            or not dimer_converged
            or not runtime_evidence["accepted"]
        ):
            fail = Artifact.failure(
                f"dimer_failed_{_reaction_id(reaction)}_{endpoint_label}",
                "ts_result",
                "Endpoint-local dimer did not produce a validated saddle candidate",
                category=(
                    "pysisyphus_dimer_method_evidence_invalid"
                    if result.ok
                    and dimer_converged
                    and not runtime_evidence["accepted"]
                    else "pysisyphus_dimer_optimizer_failed"
                ),
                parents=[attempt.artifact_id],
                dimer_convergence=dimer_evidence,
                method_evidence=runtime_evidence,
            )
            return TSSearchResult(
                [attempt, fail], record=fail.model_dump(), success=False
            )

        reaction_id = _reaction_id(reaction)
        ts_id = (
            reaction_id.replace("rxn_", "ts_", 1)
            if reaction_id.startswith("rxn_")
            else f"ts_{reaction_id}"
        )
        ts_species = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            ts_guess_path,
            ts_id,
            self.name,
            {
                "ts_guess_strategy": "endpoint_local_saddle_search",
                "dimer_optimizer_converged": True,
                "real_ts_search_executed": True,
            },
        )
        ts_species.provenance.update(
            {
                "created_by": "PysisyphusEngine",
                "electronic_state": state,
                "saddle_attempt_artifact_id": attempt.artifact_id,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
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
            str(workdir / "saddle_freq"),
        )
        calc.method = {
            **(calc.method or {}),
            "backend": self.name,
            "stage": "ts-search",
            "backend_step": "dimer_saddle_freq",
        }
        assessment = validate_ts_frequency_calculation(
            ts_species,
            calc,
            method,
            backend=self.name,
            search_method_evidence_validated=runtime_evidence["accepted"],
        )
        validated = bool(assessment["validated"])
        attempt.data.update(
            {
                "diagnosis": (
                    "saddle_validated" if validated else "saddle_frequency_rejected"
                ),
                "next_action": "validate_connectivity" if validated else "adaptive_path_refinement",
                "saddle_frequency_calculation_id": calc.artifact_id,
            }
        )
        attempt.qc.update(
            {
                "ts_validated_by_frequency": validated,
                "n_imag": calc.data.get("n_imag"),
                "mode_overlap_score": assessment["mode_overlap_score"],
            }
        )
        if not validated:
            fail = Artifact.failure(
                f"dimer_frequency_failed_{reaction_id}_{endpoint_label}",
                "ts_result",
                "Dimer candidate failed the one-imaginary-mode reaction-coordinate gate",
                category="dimer_saddle_frequency_rejected",
                parents=[attempt.artifact_id, ts_species.artifact_id, calc.artifact_id],
                n_imag=calc.data.get("n_imag"),
                imag_freq_cm1=calc.data.get("imag_freq_cm1"),
                mode_overlap_score=assessment["mode_overlap_score"],
            )
            return TSSearchResult(
                [attempt, ts_species, calc, fail],
                record=fail.model_dump(),
                success=False,
            )

        method_record = electronic_structure_method(
            engine="nwchem",
            backend=self.name,
            task="endpoint_dimer_saddle_frequency",
            config=method,
            charge=state["charge"],
            multiplicity=state["multiplicity"],
        )
        method_record["disp_vdw"] = method_record.pop("dispersion", None)
        ts_species.method = {**method_record, "task": "ts_guess"}
        reaction_validated = Artifact(
            artifact_id=reaction.artifact_id + "_with_ts",
            artifact_type="reaction_validated",
            parents=[
                reaction.artifact_id,
                attempt.artifact_id,
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
                "ts_found": True,
                "ts_validated_by_frequency": True,
                "n_imag": calc.data.get("n_imag"),
                "imag_freq_cm1": calc.data.get("imag_freq_cm1"),
                "mode_overlap_score": assessment["mode_overlap_score"],
                "real_neb_executed": False,
                "real_ts_search_executed": True,
                "real_qm_executed": calc.qc.get("real_qm_executed") is True,
                "method_evidence_validated": True,
                "dispersion_applied": calc.qc.get("dispersion_applied"),
                "fallback_dummy": False,
            },
            method=method_record,
            provenance={
                "created_by": "PysisyphusEngine",
                "electronic_state": state,
                "saddle_attempt_artifact_id": attempt.artifact_id,
                "saddle_frequency_calculation_id": calc.artifact_id,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "method_evidence": {
                    "dimer": runtime_evidence,
                    "saddle_frequency": True,
                },
            },
        )
        return TSSearchResult(
            [attempt, ts_species, calc, reaction_validated],
            record=reaction_validated.data,
            success=True,
            ts_species_id=ts_species.artifact_id,
            ts_calc_id=calc.artifact_id,
        )

    def _search_ts_with_endpoint_dimers(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict[str, Any],
        workdir: Path,
    ) -> TSSearchResult:
        artifacts: list[Artifact] = []
        last_record: Any = None
        for endpoint_label, endpoint in (
            ("reactant", reactant),
            ("product", product),
        ):
            result = self._run_endpoint_dimer(
                reaction,
                reactant,
                product,
                endpoint,
                endpoint_label,
                method,
                workdir / endpoint_label,
            )
            artifacts.extend(result.artifacts)
            last_record = result.record
            if result.success:
                return TSSearchResult(
                    artifacts,
                    record=result.record,
                    success=True,
                    ts_species_id=result.ts_species_id,
                    ts_calc_id=result.ts_calc_id,
                )
        return TSSearchResult(
            artifacts, record=last_record, success=False
        )

    def render_irc_input(
        self,
        ts_species: Artifact,
        method_in: dict[str, Any],
        path: str | Path,
    ) -> Path:
        method = self._method(method_in)
        program = str(method.get("program", method.get("calculator_program", "nwchem"))).lower()
        if program != "nwchem":
            raise ValueError("pysisyphus IRC currently supports only NWChem QCEngine")
        if str(method.get("required_program_version") or "") != "7.2.3":
            raise ValueError("pysisyphus/NWChem IRC requires required_program_version=7.2.3")
        if method.get("require_optimized_irc_endpoints", False) and not method.get(
            "optimize_irc_endpoints", False
        ):
            raise ValueError(
                "require_optimized_irc_endpoints requires optimize_irc_endpoints=true"
            )
        model = {
            "method": str(method.get("functional", method.get("xc", "pbe0"))),
            "basis": str(method.get("basis", "def2-svp")),
        }
        ts_xyz = species_xyz_path(ts_species).resolve()
        structure = read_xyz(ts_xyz)
        state = resolve_electronic_state(ts_species.data, method, structure.symbols)
        config: dict[str, Any] = {
            "geom": {"type": "cart", "fn": str(ts_xyz)},
            "calc": {
                "type": "qcengine",
                "program": program,
                "model": model,
                "keywords": _qcengine_keywords(method),
                "charge": state["charge"],
                "mult": state["multiplicity"],
                "pal": max(1, int(method.get("ncores", 1))),
            },
            "irc": {
                "type": str(method.get("irc_integrator", "eulerpc")),
                "max_cycles": int(method.get("irc_maxiter", 100)),
                "step_length": float(
                    method.get("irc_step_length", method.get("irc_step", 0.1))
                ),
                "rms_grad_thresh": float(method.get("irc_rms_grad_threshold", 1.0e-3)),
                "forward": True,
                "backward": True,
                "dump_every": int(method.get("irc_dump_every", 5)),
            },
        }
        endpoint_optimizer = str(
            method.get("endpoint_optimizer", "nwchem")
        ).lower()
        if endpoint_optimizer not in {"nwchem", "pysisyphus"}:
            raise ValueError(
                "endpoint_optimizer must be 'nwchem' or 'pysisyphus'"
            )
        if (
            method.get("optimize_irc_endpoints", False)
            and endpoint_optimizer == "pysisyphus"
        ):
            config["endopt"] = {
                "do_hess": False,
                # Redundant internals are singular for linear and near-linear
                # molecules. Cartesian endpoint optimization is the general
                # fail-closed choice after an IRC displacement.
                "geom": {"type": "cart"},
                "thresh": method.get("endpoint_opt_threshold", "gau"),
                "max_cycles": int(method.get("endpoint_opt_maxiter", 200)),
            }
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return output

    def run_irc(
        self,
        reaction: Artifact,
        ts_species: Artifact,
        reactant: Artifact,
        product: Artifact,
        method_in: dict[str, Any],
        workdir: str | Path,
    ) -> IRCResult:
        method = self._method(method_in)
        wd = Path(workdir)
        try:
            ts_structure = read_xyz(species_xyz_path(ts_species))
            state = resolve_electronic_state(ts_species.data, method, ts_structure.symbols)
            endpoint_states = [
                resolve_electronic_state(
                    endpoint.data,
                    method,
                    read_xyz(species_xyz_path(endpoint)).symbols,
                )
                for endpoint in (reactant, product)
            ]
            state_keys = ("charge", "multiplicity", "electron_count")
            if any(
                endpoint_state.get(key) != state.get(key)
                for endpoint_state in endpoint_states
                for key in state_keys
            ):
                raise ValueError("TS and validated endpoints have different electronic states")
            input_path = self.render_irc_input(
                ts_species, method, wd / "pysis_irc.yaml"
            )
        except (KeyError, OSError, TypeError, ValueError, yaml.YAMLError) as exc:
            fail = Artifact.failure(
                f"irc_failed_{reaction.artifact_id}",
                "irc",
                str(exc),
                category="pysisyphus_irc_input_invalid",
                parents=[
                    reaction.artifact_id,
                    ts_species.artifact_id,
                    reactant.artifact_id,
                    product.artifact_id,
                ],
                recommended_fallback=(
                    "match endpoint/TS state and PBE0/def2-SVPD/D3(0)/NWChem 7.2.3"
                ),
            )
            return IRCResult([fail], record=fail.model_dump(), success=False)
        method_record = _method_record(method, state)
        attempt = Artifact(
            artifact_id="irc_attempt_"
            + fingerprint_dict({"reaction": reaction.artifact_id, "engine": self.name}),
            artifact_type="irc_attempt",
            parents=[reaction.artifact_id, ts_species.artifact_id],
            paths={"input": str(input_path)},
            method=method_record,
            data={
                "reaction_id": _reaction_id(reaction),
                "ts_species_id": ts_species.artifact_id,
                "resolved_charge": state["charge"],
                "resolved_multiplicity": state["multiplicity"],
                "electron_count": state["electron_count"],
            },
            qc={"input_rendered": True, "real_irc_executed": False, "fallback_dummy": False},
            provenance={
                "created_by": "PysisyphusEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "input_sha256": _file_hash(input_path),
            },
        )
        pysis = resolve_executable("pysis", method.get("pysis_executable"))
        nwchem = resolve_executable("nwchem", method.get("executable"))
        runnable = allow_pysisyphus_subprocess(method) and not method.get("dry_run", False)
        if not runnable or pysis is None or nwchem is None:
            reason = (
                "dry_run=True"
                if method.get("dry_run", False)
                else "pysisyphus execution disabled"
                if not allow_pysisyphus_subprocess(method)
                else "pysis executable not found"
                if pysis is None
                else "NWChem executable not found"
            )
            fail = Artifact.failure(
                f"irc_failed_{reaction.artifact_id}",
                "irc",
                reason,
                category="pysisyphus_irc_not_run",
                parents=[attempt.artifact_id],
                input=str(input_path),
                recommended_fallback="install pysisyphus and QCEngine, configure NWChem, and enable subprocesses",
            )
            fail.qc = {"real_irc_executed": False, "irc_validated": False, "fallback_dummy": False}
            return IRCResult([attempt, fail], record=fail.model_dump(), success=False)
        env = _pysisyphus_environment(pysis, nwchem, method)
        result = run_command(
            [pysis, input_path.name],
            cwd=wd,
            timeout_s=int(method.get("timeout_s", 172_800)),
            env=env,
            stdout_name="pysis_irc.out",
            stderr_name="pysis_irc.err",
        )
        raw_first = wd / "finished_first.xyz"
        raw_last = wd / "finished_last.xyz"
        optimize_endpoints = bool(method.get("optimize_irc_endpoints", False))
        require_optimized_endpoints = bool(
            method.get("require_optimized_irc_endpoints", False)
        )
        endpoint_optimizer = str(
            method.get("endpoint_optimizer", "nwchem")
        ).lower()
        endpoint_calculations: list[Artifact] = []
        if (
            optimize_endpoints
            and endpoint_optimizer == "nwchem"
            and result.ok
            and raw_first.is_file()
            and raw_last.is_file()
        ):
            endpoint_optimization, endpoint_calculations = (
                _optimize_irc_endpoints_with_nwchem(
                    {"forward": raw_first, "backward": raw_last},
                    reaction_id=_reaction_id(reaction),
                    state=state,
                    method=method,
                    workdir=wd,
                )
            )
        elif optimize_endpoints and endpoint_optimizer == "pysisyphus":
            endpoint_optimization = _endpoint_optimization_evidence(wd)
        else:
            endpoint_optimization = {
                "accepted": not optimize_endpoints,
                "backend": None,
                "branches": {},
            }
        if optimize_endpoints and endpoint_optimization["accepted"]:
            first = Path(
                endpoint_optimization["branches"]["forward"]["geometry_path"]
            )
            last = Path(
                endpoint_optimization["branches"]["backward"]["geometry_path"]
            )
        else:
            first, last = raw_first, raw_last
        runtime_evidence = _nwchem_qcengine_evidence(
            wd,
            str(method["required_program_version"]),
            expected_basis=str(method.get("basis", "def2-svp")),
            require_spherical=True,
        )
        attempt.paths.update({"output": result.stdout_path, "stderr": result.stderr_path})
        attempt.data["program_version"] = (
            runtime_evidence["raw_outputs"][0]["program_version"]
            if runtime_evidence["raw_outputs"]
            else None
        )
        attempt.qc.update(
            {
                "real_irc_executed": True,
                "command_ok": result.ok,
                "method_evidence_validated": runtime_evidence["accepted"],
                "dispersion_applied": bool(
                    runtime_evidence["raw_outputs"]
                    and all(
                        item["dft_d3_applied"] is True
                        for item in runtime_evidence["raw_outputs"]
                    )
                ),
                "endpoint_optimization_requested": optimize_endpoints,
                "endpoint_optimization_converged": endpoint_optimization[
                    "accepted"
                ],
            }
        )
        attempt.provenance.update(
            {
                "command": asdict(result),
                "output_sha256": _file_hash(result.stdout_path),
                "stderr_sha256": _file_hash(result.stderr_path),
                "nwchem_qcengine_evidence": runtime_evidence,
                "endpoint_optimization_evidence": endpoint_optimization,
            }
        )
        if (
            not result.ok
            or not first.exists()
            or not last.exists()
            or (require_optimized_endpoints and not endpoint_optimization["accepted"])
        ):
            endpoint_reason = (
                "required endpoint optimizations did not both converge"
                if require_optimized_endpoints
                and not endpoint_optimization["accepted"]
                else f"endpoints_found={first.exists() and last.exists()}"
            )
            fail = Artifact.failure(
                f"irc_failed_{reaction.artifact_id}",
                "irc",
                f"pysisyphus IRC failed: returncode={result.returncode}, {endpoint_reason}",
                category=(
                    "pysisyphus_irc_endpoint_optimization_failed"
                    if require_optimized_endpoints
                    and not endpoint_optimization["accepted"]
                    else "pysisyphus_irc_failed"
                ),
                parents=[attempt.artifact_id],
                output=result.stdout_path,
                recommended_fallback=(
                    "inspect both endpoint optimizer logs and repair the IRC basin assignment"
                ),
                endpoint_optimization=endpoint_optimization,
            )
            fail.qc = {
                "real_irc_executed": True,
                "irc_validated": False,
                "endpoint_optimization_converged": endpoint_optimization[
                    "accepted"
                ],
                "fallback_dummy": False,
            }
            return IRCResult(
                [attempt, *endpoint_calculations, fail],
                record=fail.model_dump(),
                success=False,
            )
        if not runtime_evidence["accepted"]:
            fail = Artifact.failure(
                f"irc_failed_{reaction.artifact_id}",
                "irc",
                "pysisyphus completed without complete NWChem 7.2.3/D3(0) raw evidence",
                category="pysisyphus_nwchem_method_evidence_invalid",
                parents=[attempt.artifact_id],
                output=result.stdout_path,
                recommended_fallback=(
                    "preserve every qce_nwchem_stdout file and verify version, normal "
                    "termination, and D3 correction for every IRC calculation"
                ),
                method_evidence=runtime_evidence,
            )
            fail.qc = {
                "real_irc_executed": True,
                "irc_validated": False,
                "method_evidence_validated": False,
                "fallback_dummy": False,
            }
            return IRCResult(
                [attempt, *endpoint_calculations, fail],
                record=fail.model_dump(),
                success=False,
            )
        endpoint_qc = endpoint_pair_match_qc(
            first,
            last,
            species_xyz_path(reactant),
            species_xyz_path(product),
            ts_species.data,
            float(method.get("endpoint_rmsd_threshold_A", 0.75)),
            q_tolerance_A=float(method.get("endpoint_q_tolerance_A", 0.30)),
            require_identity_invariant_geometry=bool(
                method.get("require_identity_invariant_endpoint_match", True)
            ),
            permutation_rmsd_threshold_A=float(
                method.get("endpoint_permutation_rmsd_threshold_A", 0.20)
            ),
            distance_spectrum_threshold_A=float(
                method.get("endpoint_distance_spectrum_threshold_A", 0.08)
            ),
        )
        qc = {
            **endpoint_qc,
            "irc_validated": bool(endpoint_qc.get("irc_validated")),
            "real_irc_executed": True,
            "method_evidence_validated": True,
            "dispersion_applied": True,
            "fallback_dummy": False,
            "command_ok": True,
            "endpoint_optimization_requested": optimize_endpoints,
            "endpoint_optimization_converged": endpoint_optimization[
                "accepted"
            ],
            "optimized_endpoints_required": require_optimized_endpoints,
        }
        record = {
            "reaction_id": _reaction_id(reaction),
            "ts_species_id": ts_species.artifact_id,
            "irc_backend": self.name,
            "calculator_program": "nwchem",
            "forward_endpoint_xyz": str(first),
            "backward_endpoint_xyz": str(last),
            "resolved_charge": state["charge"],
            "resolved_multiplicity": state["multiplicity"],
            "electron_count": state["electron_count"],
            "program_version": runtime_evidence["raw_outputs"][0]["program_version"],
            "endpoint_geometry_source": (
                "converged_endpoint_optimizations"
                if endpoint_optimization["accepted"]
                else "raw_irc_termini"
            ),
            **qc,
        }
        artifact = Artifact(
            artifact_id=str(reaction.artifact_id).replace("_with_ts", "_irc"),
            artifact_type="irc",
            parents=[attempt.artifact_id, reaction.artifact_id, ts_species.artifact_id],
            paths={
                "input": str(input_path),
                "output": result.stdout_path,
                "forward_xyz": str(first),
                "backward_xyz": str(last),
                "forward_optimizer_log": endpoint_optimization["branches"].get(
                    "forward", {}
                ).get("optimizer_log_path")
                or "",
                "backward_optimizer_log": endpoint_optimization["branches"].get(
                    "backward", {}
                ).get("optimizer_log_path")
                or "",
            },
            method=method_record,
            data=record,
            qc=qc,
            provenance={
                "command": asdict(result),
                "created_by": "PysisyphusEngine",
                "electronic_state": state,
                "endpoint_artifact_ids": {
                    "reactant": reactant.artifact_id,
                    "product": product.artifact_id,
                },
                "input_sha256": _file_hash(input_path),
                "output_sha256": _file_hash(result.stdout_path),
                "stderr_sha256": _file_hash(result.stderr_path),
                "forward_xyz_sha256": _file_hash(first),
                "backward_xyz_sha256": _file_hash(last),
                "nwchem_qcengine_evidence": runtime_evidence,
                "endpoint_optimization_evidence": endpoint_optimization,
            },
        )
        return IRCResult(
            [attempt, *endpoint_calculations, artifact],
            record=record,
            success=bool(qc["irc_validated"]),
        )
