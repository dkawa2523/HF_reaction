"""Hessian-based pysisyphus saddle refinement for path-evidenced seeds."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from hfauto.backends.qm.nwchem import NWChemEngine, species_xyz_path
from hfauto.backends.ts.base import (
    TSSearchResult,
    make_ts_species_artifact,
    validate_saddle_seed_hessian,
    validate_ts_frequency_calculation,
)
from hfauto.backends.ts.pysisyphus import (
    _file_hash,
    _nwchem_qcengine_evidence,
    _optimizer_converged,
    _pysisyphus_environment,
    _qcengine_keywords,
    allow_pysisyphus_subprocess,
)
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.methods import electronic_structure_method
from hfauto.chemistry.xyz import read_xyz
from hfauto.chemistry.xyz_trajectory import endpoint_mode_reference
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


_BOHR_PER_ANGSTROM = 1.8897261254578281


def _validated_hessian_initialization(
    seed_xyz: Path, method: dict[str, Any]
) -> tuple[str, dict[str, Any]]:
    """Return an exact-Hessian source, validating any retry cache fail-closed."""

    cache_value = method.get("ts_hessian_init_h5")
    if not cache_value:
        return "calc", {"reused": False, "accepted": True}
    cache = Path(str(cache_value)).resolve()
    evidence_dir = Path(str(method.get("ts_hessian_evidence_dir") or ""))
    if cache.suffix.lower() != ".h5" or not cache.is_file():
        raise ValueError("ts_hessian_init_h5 must be a readable pysisyphus HDF5 file")
    if not evidence_dir.is_dir():
        raise ValueError("ts_hessian_evidence_dir is required for Hessian cache reuse")

    try:
        import h5py

        with h5py.File(cache, "r") as handle:
            cached_coords = np.asarray(handle["coords3d"][:], dtype=float)
            cached_hessian = np.asarray(handle["hessian"][:], dtype=float)
    except (ImportError, KeyError, OSError, ValueError) as exc:
        raise ValueError(f"invalid pysisyphus Hessian cache: {exc}") from exc

    structure = read_xyz(seed_xyz)
    expected_coords = structure.coords * _BOHR_PER_ANGSTROM
    coordinate_delta = (
        float(np.max(np.abs(cached_coords - expected_coords)))
        if cached_coords.shape == expected_coords.shape
        else None
    )
    expected_hessian_shape = (3 * len(structure.symbols),) * 2
    runtime_evidence = _nwchem_qcengine_evidence(
        evidence_dir,
        str(method["required_program_version"]),
        expected_basis=str(method.get("basis", "def2-svp")),
        require_spherical=True,
    )
    accepted = bool(
        coordinate_delta is not None
        and coordinate_delta <= 1.0e-6
        and cached_hessian.shape == expected_hessian_shape
        and np.isfinite(cached_hessian).all()
        and runtime_evidence["accepted"]
    )
    evidence = {
        "reused": True,
        "accepted": accepted,
        "path": str(cache),
        "sha256": _file_hash(cache),
        "coordinate_max_abs_delta_bohr": coordinate_delta,
        "hessian_shape": list(cached_hessian.shape),
        "nwchem_evidence": runtime_evidence,
    }
    if not accepted:
        raise ValueError("Hessian cache failed coordinate, shape, or PES validation")
    return str(cache), evidence


class PysisyphusSaddleEngine:
    """Refine one validated path seed with an exact-Hessian RS-I-RFO search."""

    name = "pysisyphus_saddle"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def render_input(
        self,
        seed: Artifact,
        method_in: dict[str, Any],
        path: str | Path,
    ) -> Path:
        """Render a Cartesian RS-I-RFO job backed by QCEngine/NWChem."""

        method = {**self.config, **(method_in or {})}
        program = str(
            method.get("program", method.get("calculator_program", "nwchem"))
        ).lower()
        if program != "nwchem":
            raise ValueError("pysisyphus saddle refinement supports NWChem only")
        if str(method.get("required_program_version") or "") != "7.2.3":
            raise ValueError(
                "pysisyphus/NWChem saddle refinement requires version 7.2.3"
            )
        # Cartesian coordinates are the backend-neutral default: they do not
        # depend on a particular fragment or internal-coordinate perception.
        # Callers may opt into redundant/DLC/TRIC coordinates per system.
        coordinate_type = str(method.get("ts_coordinate_type", "cart")).lower()
        if coordinate_type not in {"cart", "redund", "dlc", "tric"}:
            raise ValueError(
                "ts_coordinate_type must be cart, redund, dlc, or tric"
            )
        seed_xyz = species_xyz_path(seed).resolve()
        structure = read_xyz(seed_xyz)
        state = resolve_electronic_state(seed.data, method, structure.symbols)
        hessian_init = method.get("_validated_ts_hessian_init")
        if not hessian_init:
            hessian_init, _ = _validated_hessian_initialization(seed_xyz, method)
        config = {
            "geom": {"type": coordinate_type, "fn": str(seed_xyz)},
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
            "tsopt": {
                "type": "rsirfo",
                "hessian_init": str(hessian_init),
                "hessian_update": "bofill",
                "root": 0,
                "assert_neg_eigval": True,
                "trust_radius": float(method.get("ts_trust_radius", 0.10)),
                "trust_min": float(method.get("ts_trust_min", 0.005)),
                "trust_max": float(method.get("ts_trust_max", 0.20)),
                "max_cycles": int(method.get("ts_maxiter", 100)),
                "thresh": str(method.get("ts_opt_threshold", "gau_tight")),
                "dump": True,
            },
        }
        hessian_recalc = method.get("ts_hessian_recalc")
        if hessian_recalc is not None:
            config["tsopt"]["hessian_recalc"] = max(1, int(hessian_recalc))
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return output

    @staticmethod
    def _failure(
        artifacts: list[Artifact],
        reaction: Artifact,
        attempt: Artifact,
        message: str,
        category: str,
        **data: Any,
    ) -> TSSearchResult:
        fail = Artifact.failure(
            f"pysis_saddle_failed_{_reaction_id(reaction)}",
            "ts_result",
            message,
            category=category,
            parents=[attempt.artifact_id],
            **data,
        )
        return TSSearchResult(
            [*artifacts, fail], record=fail.model_dump(), success=False
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
        candidate = dict(method.get("saddle_seed_candidate") or {})
        seed_path = Path(str(candidate.get("xyz_path") or ""))
        token = fingerprint_dict(
            {"reaction": reaction_id, "seed": str(seed_path), "engine": self.name}
        )
        attempt = Artifact(
            artifact_id=f"saddle_attempt_{token}",
            artifact_type="saddle_attempt",
            parents=[reaction.artifact_id],
            paths={"seed_xyz": str(seed_path)},
            data={
                "reaction_id": reaction_id,
                "strategy": "exact_hessian_rsirfo_refinement",
                "engine": self.name,
                "diagnosis": "input_pending",
                "next_action": "validate_seed_hessian",
                "candidate": candidate,
            },
            qc={"fallback_dummy": False},
            provenance={"created_by": "PysisyphusSaddleEngine"},
        )
        artifacts = [attempt]
        if (
            method.get("path_strategy") != "saddle_refinement"
            or not method.get("saddle_seed_evidence_validated")
            or not seed_path.is_file()
        ):
            return self._failure(
                artifacts,
                reaction,
                attempt,
                "RS-I-RFO refinement requires one readable path-evidenced seed",
                "pysisyphus_saddle_seed_invalid",
            )

        structure = read_xyz(seed_path)
        state = resolve_electronic_state(reactant.data, method, structure.symbols)
        ts_id = (
            reaction_id.replace("rxn_", "ts_", 1)
            if reaction_id.startswith("rxn_")
            else f"ts_{reaction_id}"
        )
        seed = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            seed_path,
            f"{ts_id}_seed",
            self.name,
            {
                "ts_guess_strategy": "converged_double_ended_path",
                "scientific_role": "mode_following_seed_only",
            },
        )
        mode_reference = candidate.get("reaction_mode_reference")
        if not isinstance(mode_reference, dict):
            mode_reference = endpoint_mode_reference(
                read_xyz(species_xyz_path(reactant)),
                read_xyz(species_xyz_path(product)),
            )
        seed.data["reaction_mode_reference"] = mode_reference
        supplied_hessian = method.get("saddle_seed_hessian_evidence")
        seed_hessian = (
            supplied_hessian
            if isinstance(supplied_hessian, Artifact)
            else Artifact.model_validate(supplied_hessian)
            if isinstance(supplied_hessian, dict)
            else None
        )
        if seed_hessian is None:
            if method.get("require_supplied_seed_hessian", False):
                artifacts.append(seed)
                return self._failure(
                    artifacts,
                    reaction,
                    attempt,
                    "An independently calculated seed Hessian is required",
                    "pysisyphus_saddle_seed_hessian_missing",
                )
            seed_hessian = NWChemEngine().frequency(
                seed,
                {
                    **method,
                    "method_id": method.get(
                        "seed_hessian_method_id",
                        method.get("method_id", "nwchem_saddle_seed_hessian"),
                    ),
                },
                str(workdir / "seed_hessian"),
            )
        seed_assessment = validate_saddle_seed_hessian(
            seed, seed_hessian, method, backend=self.name
        )
        artifacts.extend([seed, seed_hessian])
        unique_mode = seed_hessian.data.get("n_imag") == 1
        seed_accepted = bool(seed_assessment["accepted"] and unique_mode)
        attempt.qc.update(
            {
                "saddle_seed_hessian_accepted": seed_accepted,
                "seed_n_imag": seed_hessian.data.get("n_imag"),
                "seed_mode_overlap_score": seed_assessment["mode_overlap_score"],
                "seed_driver_mode_mapping_validated": unique_mode,
                "seed_hessian_reused": supplied_hessian is not None,
            }
        )
        if not seed_accepted:
            attempt.data.update(
                {
                    "diagnosis": "saddle_seed_hessian_rejected",
                    "next_action": "adaptive_path_refinement",
                }
            )
            return self._failure(
                artifacts,
                reaction,
                attempt,
                "Seed Hessian lacks one unique path-aligned negative mode",
                "pysisyphus_saddle_seed_hessian_rejected",
            )

        try:
            hessian_init, hessian_cache_evidence = _validated_hessian_initialization(
                seed_path.resolve(), method
            )
        except (OSError, TypeError, ValueError) as exc:
            return self._failure(
                artifacts,
                reaction,
                attempt,
                str(exc),
                "pysisyphus_saddle_hessian_cache_invalid",
            )
        method = {**method, "_validated_ts_hessian_init": hessian_init}
        attempt.qc["saddle_hessian_cache_validated"] = bool(
            hessian_cache_evidence["accepted"]
        )
        attempt.provenance["saddle_hessian_cache"] = hessian_cache_evidence

        try:
            input_path = self.render_input(seed, method, workdir / "pysis_saddle.yaml")
        except (OSError, TypeError, ValueError, yaml.YAMLError) as exc:
            return self._failure(
                artifacts,
                reaction,
                attempt,
                str(exc),
                "pysisyphus_saddle_input_invalid",
            )
        attempt.paths["input"] = str(input_path)
        pysis = resolve_executable("pysis", method.get("pysis_executable"))
        nwchem = resolve_executable("nwchem", method.get("executable"))
        runnable = allow_pysisyphus_subprocess(method) and not method.get(
            "dry_run", False
        )
        if not runnable or pysis is None or nwchem is None:
            reason = (
                "execution disabled"
                if not runnable
                else "pysis executable not found"
                if pysis is None
                else "NWChem executable not found"
            )
            return self._failure(
                artifacts,
                reaction,
                attempt,
                reason,
                "pysisyphus_saddle_not_run",
            )

        result = run_command(
            [pysis, input_path.name],
            cwd=workdir,
            timeout_s=int(method.get("timeout_s", 172_800)),
            env=_pysisyphus_environment(pysis, nwchem, method),
            stdout_name="pysis_saddle.out",
            stderr_name="pysis_saddle.err",
        )
        final_xyz = workdir / "ts_opt.xyz"
        optimizer_log = workdir / "ts_optimizer.log"
        optimizer_converged = _optimizer_converged(optimizer_log, final_xyz)
        runtime_evidence = _nwchem_qcengine_evidence(
            workdir,
            str(method["required_program_version"]),
            expected_basis=str(method.get("basis", "def2-svp")),
            require_spherical=True,
        )
        search_validated = bool(
            result.ok and optimizer_converged and runtime_evidence["accepted"]
        )
        attempt.paths.update(
            {
                "output": result.stdout_path,
                "stderr": result.stderr_path,
                "optimizer_log": str(optimizer_log) if optimizer_log.is_file() else "",
                "ts_guess_xyz": str(final_xyz) if final_xyz.is_file() else "",
            }
        )
        attempt.qc.update(
            {
                "real_ts_search_executed": True,
                "command_ok": result.ok,
                "optimizer_converged": optimizer_converged,
                "method_evidence_validated": runtime_evidence["accepted"],
            }
        )
        attempt.provenance.update(
            {
                "command": asdict(result),
                "input_sha256": _file_hash(input_path),
                "optimizer_log_sha256": _file_hash(optimizer_log),
                "nwchem_qcengine_evidence": runtime_evidence,
            }
        )
        if not search_validated:
            attempt.data.update(
                {
                    "diagnosis": "saddle_optimizer_failed",
                    "next_action": "adaptive_path_refinement",
                }
            )
            return self._failure(
                artifacts,
                reaction,
                attempt,
                "RS-I-RFO did not converge with complete NWChem/D3 evidence",
                "pysisyphus_rsirfo_failed",
                optimizer_converged=optimizer_converged,
                method_evidence=runtime_evidence,
            )

        ts_species = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            final_xyz,
            ts_id,
            self.name,
            {
                "ts_guess_strategy": "exact_hessian_rsirfo_refinement",
                "reaction_mode_reference": mode_reference,
                "optimizer_converged": True,
                "real_ts_search_executed": True,
            },
        )
        ts_species.provenance.update(
            {
                "created_by": "PysisyphusSaddleEngine",
                "electronic_state": state,
                "saddle_attempt_artifact_id": attempt.artifact_id,
            }
        )
        calculation = NWChemEngine().frequency(
            ts_species,
            {
                **method,
                "method_id": method.get(
                    "saddle_method_id", method.get("method_id", "nwchem_saddle")
                ),
            },
            str(workdir / "final_frequency"),
        )
        calculation.method = {
            **(calculation.method or {}),
            "backend": self.name,
            "stage": "ts-search",
            "backend_step": "rsirfo_final_fixed_geometry_frequency",
        }
        assessment = validate_ts_frequency_calculation(
            ts_species,
            calculation,
            method,
            backend=self.name,
            search_method_evidence_validated=search_validated,
        )
        validated = bool(assessment["validated"])
        attempt.data.update(
            {
                "diagnosis": (
                    "saddle_validated" if validated else "saddle_frequency_rejected"
                ),
                "next_action": (
                    "validate_connectivity"
                    if validated
                    else "adaptive_path_refinement"
                ),
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
        artifacts.extend([ts_species, calculation])
        if not validated:
            return self._failure(
                artifacts,
                reaction,
                attempt,
                "RS-I-RFO candidate failed the final first-order-saddle gate",
                "pysisyphus_saddle_frequency_rejected",
                n_imag=calculation.data.get("n_imag"),
                mode_overlap_score=assessment["mode_overlap_score"],
            )

        method_record = electronic_structure_method(
            engine="nwchem",
            backend=self.name,
            task="rsirfo_saddle_frequency",
            config=method,
            charge=state["charge"],
            multiplicity=state["multiplicity"],
        )
        method_record["disp_vdw"] = method_record.pop("dispersion", None)
        ts_species.method = {**method_record, "task": "transition_state"}
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
                "ts_energy_above_endpoints": assessment.get(
                    "ts_energy_above_endpoints"
                ),
                "real_ts_search_executed": True,
                "real_qm_executed": True,
                "method_evidence_validated": True,
                "dispersion_applied": calculation.qc.get("dispersion_applied"),
                "fallback_dummy": False,
            },
            method=method_record,
            provenance={
                "created_by": "PysisyphusSaddleEngine",
                "electronic_state": state,
                "saddle_attempt_artifact_id": attempt.artifact_id,
                "saddle_frequency_calculation_id": calculation.artifact_id,
                "nwchem_qcengine_evidence": runtime_evidence,
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
