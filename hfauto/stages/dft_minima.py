from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.chemistry.basin_identity import permutation_invariant_graph_rmsd
from hfauto.chemistry.method_lineage import path_numerical_settings
from hfauto.chemistry.nwchem_evidence import (
    exact_int,
    parse_nwchem_input_evidence,
    paths_equal,
    source_file_hashes,
)
from hfauto.chemistry.proton_transfer import proton_transfer_metrics
from hfauto.chemistry.reactions import is_proton_transfer_state
from hfauto.chemistry.stoichiometry import molecular_surface_key
from hfauto.core.artifacts import (
    artifact_data_matches,
    canonical_species_id,
    preferred_species,
    preferred_species_by_id,
    species_xyz_path,
)
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import ensure_dir, read_manifest, write_jsonl, write_manifest
from hfauto.core.qc import minimum_promotion_gate, minimum_qc_from_freq
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _minimum_request_fingerprint(
    species: Artifact,
    *,
    engine: str,
    method: dict[str, Any],
) -> str:
    """Identify one minimum calculation without depending on its work directory."""

    return fingerprint_dict(
        {
            "species_id": canonical_species_id(species),
            "source_geometry_artifact_id": species.artifact_id,
            "input_xyz_sha256": sha256_file(species_xyz_path(species)),
            "engine": str(engine).lower(),
            "method": method,
        }
    )


def _checkpoint_files_are_unchanged(
    checkpoint: Manifest,
    calculation: Artifact,
) -> bool:
    """Require hashes made after parsing before reusing an external result."""

    for optimized in checkpoint.latest_artifacts("species_optimized"):
        if optimized.qc.get("dft_calc_id") != calculation.artifact_id:
            continue
        integrity = optimized.provenance.get("validated_source_integrity", {}) or {}
        files = integrity.get("files", {}) or {}
        if integrity.get("accepted") is not True or not files:
            return False
        try:
            return all(
                Path(str(evidence["path"])).is_file()
                and sha256_file(str(evidence["path"])) == evidence["sha256"]
                for evidence in files.values()
            )
        except (KeyError, OSError, TypeError):
            return False
    return False


def _reusable_minimum_calculation(
    checkpoint: Manifest | None,
    species: Artifact,
    *,
    engine: str,
    method: dict[str, Any],
    require_real_qm: bool,
) -> Artifact | None:
    """Return an exact, integrity-checked completed job from a stage checkpoint."""

    if checkpoint is None:
        return None
    expected_fingerprint = _minimum_request_fingerprint(
        species,
        engine=engine,
        method=method,
    )
    expected_input_hash = sha256_file(species_xyz_path(species))
    for calculation in reversed(checkpoint.latest_artifacts("calculation")):
        if (
            calculation.status.status != "success"
            or calculation.data.get("species_id") != canonical_species_id(species)
            or calculation.data.get("source_geometry_artifact_id")
            != species.artifact_id
            or str((calculation.method or {}).get("engine", "")).lower()
            != str(engine).lower()
            or any(
                (calculation.method or {}).get(key) != value
                for key, value in method.items()
            )
            or calculation.provenance.get("input_xyz_sha256")
            != expected_input_hash
            or calculation.qc.get("minimum_accepted") is not True
            or (
                require_real_qm
                and (
                    calculation.qc.get("real_qm_executed") is not True
                    or calculation.qc.get("fallback_dummy") is not False
                )
            )
        ):
            continue
        stored_fingerprint = calculation.provenance.get(
            "dft_minima_request_fingerprint"
        )
        if stored_fingerprint not in (None, expected_fingerprint):
            continue
        command = calculation.provenance.get("command", {}) or {}
        if command.get("returncode") != 0 or command.get("timed_out") is not False:
            continue
        if not _checkpoint_files_are_unchanged(checkpoint, calculation):
            continue
        reused = calculation.model_copy(deep=True)
        reused.provenance["dft_minima_request_fingerprint"] = expected_fingerprint
        reused.provenance["reused_from_checkpoint"] = True
        return reused
    return None


def _next_minimum_workdir(base: Path) -> Path:
    """Preserve an incomplete/failed external attempt when a job is retried."""

    if not base.exists() or not any(base.iterdir()):
        return base
    attempt = 2
    while True:
        candidate = base / f"attempt_{attempt:02d}"
        if not candidate.exists() or not any(candidate.iterdir()):
            return candidate
        attempt += 1


def _proton_state_gate(species: Artifact, calc: Artifact) -> dict[str, Any]:
    """Reject a converged minimum when RC/IP optimisation changed proton state."""
    state = str(species.data.get("state") or "")
    if state not in {"reactant_complex", "ion_pair"} or not is_proton_transfer_state(species):
        return {"accepted": True, "reasons": [], "expected_state": state}
    coordinate = species.data.get("reaction_coordinate", {}) or {}
    atoms = coordinate.get("atoms", {}) or {}
    final_xyz = calc.paths.get("final_xyz")
    try:
        metrics = proton_transfer_metrics(final_xyz, **atoms)
    except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
        return {
            "accepted": False,
            "reasons": [f"proton_state_unreadable:{exc}"],
            "expected_state": state,
        }
    expected_class = "neutral_complex" if state == "reactant_complex" else "ion_pair"
    accepted = metrics["endpoint_class"] == expected_class
    return {
        "accepted": accepted,
        "reasons": []
        if accepted
        else [f"expected_{expected_class}_got_{metrics['endpoint_class']}"],
        "expected_state": state,
        "expected_class": expected_class,
        "metrics": metrics,
    }


def _preferred_species_inputs(
    manifest: Manifest,
    states: set[str],
    *,
    include_optimized: bool = False,
) -> list[Artifact]:
    """Prefer preoptimized geometry artifacts while preserving canonical species IDs."""

    artifact_types = ["species", "species_preopt"]
    if include_optimized:
        artifact_types.append("species_optimized")
    return preferred_species(
        manifest,
        states=states,
        artifact_types=artifact_types,
    )


def _candidate_endpoint_species_ids(
    manifest: Manifest,
    maximum_candidates: int,
    *,
    duplicate_rmsd_A: float = 0.05,
) -> tuple[set[str], list[str]]:
    """Select bounded pairs without comparing energies across compositions.

    Candidates are ordered by low-level energy only within one molecular PES,
    then selected round-robin across composition/charge/multiplicity keys.  A
    larger system therefore cannot consume the DFT budget merely because its
    absolute electronic energy is more negative.
    """

    if maximum_candidates < 1:
        raise ValueError("dft-minima candidate_endpoint_limit must be positive")
    species = preferred_species_by_id(manifest)
    candidates = [
        candidate
        for candidate in manifest.latest_artifacts("reaction_candidate")
        if candidate.status.status == "success"
        and candidate.qc.get("structural_change_detected") is True
    ]

    def discovery_energy(candidate: Artifact) -> float:
        product = species.get(str(candidate.data.get("product_species_id")))
        value = (
            product.data.get("discovery_electronic_energy_hartree")
            if product is not None
            else None
        )
        return float(value) if value is not None else float("inf")

    queues: dict[str, list[Artifact]] = {}
    for candidate in candidates:
        reactant = species.get(str(candidate.data.get("reactant_species_id")))
        product = species.get(str(candidate.data.get("product_species_id")))
        surface_species = reactant or product
        if product is None or surface_species is None:
            continue
        queues.setdefault(molecular_surface_key(surface_species), []).append(
            candidate
        )
    for queue in queues.values():
        queue.sort(
            key=lambda candidate: (
                discovery_energy(candidate) == float("inf"),
                discovery_energy(candidate),
                str(
                    candidate.data.get("candidate_id")
                    or candidate.artifact_id
                ),
            )
        )

    selected: list[Artifact] = []
    selected_product_paths: dict[str, list[Path]] = {
        surface: [] for surface in queues
    }
    while queues and len(selected) < maximum_candidates:
        progress = False
        for surface in sorted(queues):
            queue = queues[surface]
            accepted: Artifact | None = None
            while queue:
                candidate = queue.pop(0)
                product = species.get(
                    str(candidate.data.get("product_species_id"))
                )
                if product is None:
                    continue
                product_path = species_xyz_path(product)
                duplicate = any(
                    (
                        rmsd := permutation_invariant_graph_rmsd(
                            product_path, other_path
                        )
                    )
                    is not None
                    and rmsd <= float(duplicate_rmsd_A)
                    for other_path in selected_product_paths[surface]
                )
                if duplicate:
                    continue
                accepted = candidate
                selected_product_paths[surface].append(product_path)
                break
            if not queue:
                queues.pop(surface, None)
            if accepted is None:
                continue
            selected.append(accepted)
            progress = True
            if len(selected) >= maximum_candidates:
                break
        if not progress:
            break
    endpoint_ids = {
        str(species_id)
        for candidate in selected
        for species_id in (
            candidate.data.get("reactant_species_id"),
            candidate.data.get("product_species_id"),
        )
        if species_id is not None
    }
    return endpoint_ids, [
        str(candidate.data.get("candidate_id") or candidate.artifact_id)
        for candidate in selected
    ]


def _anchor_minimum_species_ids(
    species_inputs: list[Artifact],
    limits: dict[str, int],
    *,
    duplicate_rmsd_A: float = 0.05,
) -> set[str]:
    """Select a bounded, composition-balanced NCI ensemble.

    Low-level and DFT potential-energy surfaces need not have identical basin
    topology.  These anchors let minimum-registry discover that discrepancy
    without weakening the reaction-candidate gate.  Budgets are distributed
    round-robin across molecular surfaces, so absolute energies of different
    formulas are never used to choose between them.
    """

    if any(limit < 1 for limit in limits.values()):
        raise ValueError("dft-minima anchor_minima_per_state limits must be positive")
    selected: set[str] = set()
    for state, limit in limits.items():
        queues: dict[str, list[Artifact]] = {}
        for species in species_inputs:
            if str(species.data.get("state") or "") != state:
                continue
            queues.setdefault(molecular_surface_key(species), []).append(
                species
            )
        selected_paths: dict[str, list[Path]] = {
            surface: [] for surface in queues
        }
        selected_in_state = 0
        while queues and selected_in_state < limit:
            progress = False
            for surface in sorted(queues):
                queue = queues[surface]
                accepted: Artifact | None = None
                while queue:
                    candidate = queue.pop(0)
                    try:
                        path = species_xyz_path(candidate)
                    except (FileNotFoundError, ValueError):
                        accepted = candidate
                        break
                    duplicate = any(
                        other.is_file()
                        and (
                            (
                                rmsd := permutation_invariant_graph_rmsd(
                                    path, other
                                )
                            )
                            is not None
                            and rmsd <= float(duplicate_rmsd_A)
                        )
                        for other in selected_paths[surface]
                    )
                    if duplicate:
                        continue
                    accepted = candidate
                    selected_paths[surface].append(path)
                    break
                if not queue:
                    queues.pop(surface, None)
                if accepted is None:
                    continue
                selected.add(canonical_species_id(accepted))
                selected_in_state += 1
                progress = True
                if selected_in_state >= limit:
                    break
            if not progress:
                break
    return selected


def _validated_minimum_provenance(
    calc: Artifact,
) -> tuple[dict[str, Any], bool, bool]:
    """Normalize the real minimum method/state evidence for downstream paths."""

    method = calc.method or {}
    data = calc.data or {}
    validated_method = {
        "engine": method.get("engine"),
        "functional": method.get("functional"),
        "basis": method.get("basis"),
        "disp_vdw": method.get("disp_vdw"),
        "program_version": data.get("program_version"),
    }
    electronic_state = {
        "charge": data.get("resolved_charge"),
        "multiplicity": data.get("resolved_multiplicity"),
        "electron_count": data.get("electron_count"),
    }
    required_version = method.get("required_program_version")
    evidence_validated = bool(
        calc.qc.get("real_qm_executed") is True
        and calc.qc.get("fallback_dummy") is False
        and calc.qc.get("minimum_accepted") is True
        and all(
            validated_method.get(key) is not None
            for key in ("engine", "functional", "basis", "program_version")
        )
        and all(value is not None for value in electronic_state.values())
        and (method.get("disp_vdw") is None or calc.qc.get("dispersion_applied") is True)
        and (required_version is None or str(data.get("program_version")) == str(required_version))
    )
    numerical_settings = path_numerical_settings(method)
    numerical_evidence_validated = True
    numerical_input_evidence: dict[str, Any] = {}
    validated_source_files: dict[str, dict[str, str]] = {}
    input_method_evidence: dict[str, Any] = {
        "accepted": True,
        "reasons": [],
    }
    if str(method.get("engine", "")).lower() == "nwchem":
        input_path = str(calc.paths.get("input") or "")
        try:
            raw_input = parse_nwchem_input_evidence(input_path)
            numerical_input_evidence = {
                "path": input_path,
                "sha256": sha256_file(input_path),
                "grid": raw_input.get("grid"),
                "scf_energy_tolerance": raw_input.get(
                    "scf_energy_tolerance"
                ),
            }
        except (OSError, TypeError, UnicodeError, ValueError):
            raw_input = {}
        expected_input = {
            "charge": exact_int(data.get("resolved_charge")),
            "multiplicity": exact_int(data.get("resolved_multiplicity")),
            "functional": str(method.get("functional") or "").lower() or None,
            "basis": str(method.get("basis") or "").lower() or None,
            "disp_vdw": exact_int(method.get("disp_vdw")),
            **numerical_settings,
        }
        input_reasons = [
            f"nwchem_input_{key}_does_not_match_calculation"
            for key, expected in expected_input.items()
            if raw_input.get(key) != expected
        ]
        input_method_evidence = {
            "accepted": not input_reasons,
            "reasons": input_reasons,
            "expected": expected_input,
            "observed": {
                key: raw_input.get(key) for key in expected_input
            },
            "path": input_path,
            "sha256": numerical_input_evidence.get("sha256"),
        }
        numerical_evidence_validated = bool(
            raw_input.get("grid") == numerical_settings.get("grid")
            and raw_input.get("scf_energy_tolerance")
            == numerical_settings.get("scf_energy_tolerance")
        )
        evidence_validated = bool(
            evidence_validated and input_method_evidence["accepted"]
        )
        source_hashes, source_reasons = source_file_hashes(calc)
        command = (calc.provenance or {}).get("command")
        command_reasons: list[str] = []
        if not isinstance(command, dict):
            command = {}
            command_reasons.append("source_command_evidence_missing")
        if exact_int(command.get("returncode")) != 0:
            command_reasons.append("source_command_returncode_is_not_zero")
        if command.get("timed_out") is not False:
            command_reasons.append("source_command_timed_out_is_not_false")
        for command_key, path_key in (
            ("stdout_path", "output"),
            ("stderr_path", "stderr"),
        ):
            if not paths_equal(
                command.get(command_key), calc.paths.get(path_key)
            ):
                command_reasons.append(
                    f"source_command_{command_key}_path_mismatch"
                )
        validated_source_files = {
            key: {"path": str(calc.paths.get(key)), "sha256": digest}
            for key, digest in source_hashes.items()
        }
        source_integrity = {
            "accepted": not source_reasons and not command_reasons,
            "reasons": [*source_reasons, *command_reasons],
            "files": validated_source_files,
        }
        evidence_validated = bool(
            evidence_validated and source_integrity["accepted"]
        )
    else:
        source_integrity = {
            "accepted": evidence_validated,
            "reasons": [],
            "files": {},
        }
    provenance = {
        **(calc.provenance or {}),
        "validated_qm_method": validated_method,
        "validated_electronic_state": electronic_state,
        "validated_numerical_settings": (
            numerical_settings if numerical_evidence_validated else {}
        ),
        "validated_numerical_input": numerical_input_evidence,
        "validated_input_method": input_method_evidence,
        "validated_source_integrity": source_integrity,
    }
    if str(method.get("engine", "")).lower() == "nwchem":
        provenance["validated_nwchem_method"] = validated_method
    return (
        provenance,
        bool(evidence_validated and numerical_evidence_validated),
        numerical_evidence_validated,
    )


def _optimized_species_artifacts(
    species: Artifact, calc: Artifact, update_species: bool
) -> list[Artifact]:
    final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
    canonical_id = canonical_species_id(species)
    data = dict(species.data)
    data["species_id"] = canonical_id
    data.setdefault("source_species_id", canonical_id)
    data["minimum_input_species_id"] = canonical_id
    data["source_geometry_artifact_id"] = species.artifact_id
    data["xyz_path"] = str(final_xyz)
    for key in ("resolved_charge", "resolved_multiplicity", "electron_count", "program_version"):
        if calc.data.get(key) is not None:
            data[key] = calc.data[key]
    data["dft_minima"] = {
        "calc_id": calc.artifact_id,
        "engine": (calc.method or {}).get("engine"),
        "method_id": (calc.method or {}).get("method_id"),
        "real_orca_executed": calc.qc.get("real_orca_executed", False),
        "real_qm_executed": calc.qc.get(
            "real_qm_executed", calc.qc.get("real_orca_executed", False)
        ),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "n_imag": calc.data.get("n_imag"),
        "minimum_accepted": calc.qc.get("minimum_accepted"),
        "proton_state_accepted": calc.qc.get("proton_state_accepted"),
    }
    (
        provenance,
        method_evidence_validated,
        numerical_evidence_validated,
    ) = _validated_minimum_provenance(calc)
    qc = {
        "optimized_by_dft_minima": True,
        "dft_calc_id": calc.artifact_id,
        "method_id": (calc.method or {}).get("method_id"),
        "engine": (calc.method or {}).get("engine"),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "real_orca_executed": calc.qc.get("real_orca_executed", False),
        "real_qm_executed": calc.qc.get(
            "real_qm_executed", calc.qc.get("real_orca_executed", False)
        ),
        "scf_converged": calc.qc.get("scf_converged"),
        "geometry_converged": calc.qc.get("geometry_converged"),
        "n_imag": calc.data.get("n_imag"),
        "is_minimum": calc.qc.get("is_minimum"),
        "minimum_accepted": calc.qc.get("minimum_accepted"),
        "proton_state_accepted": calc.qc.get("proton_state_accepted"),
        "proton_state_gate": calc.qc.get("proton_state_gate"),
        "geometry_sane": calc.qc.get("geometry_sane")
        or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
        "geometry_qc": calc.qc.get("geometry_qc", {}),
        "state_method_evidence_validated": method_evidence_validated,
        "numerical_method_evidence_validated": numerical_evidence_validated,
    }
    artifacts = [
        Artifact(
            artifact_id=f"opt_{canonical_id}",
            artifact_type="species_optimized",
            parents=[species.artifact_id, calc.artifact_id],
            paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
            data=data,
            method={"stage": "dft-minima", **(calc.method or {})},
            qc=qc,
            provenance=provenance,
        )
    ]
    if update_species:
        artifacts.append(
            Artifact(
                artifact_id=canonical_id,
                artifact_type="species",
                parents=[species.artifact_id, calc.artifact_id],
                paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
                data=data,
                method={"stage": "dft-minima", **(calc.method or {})},
                qc=qc,
                provenance=provenance,
            )
        )
    return artifacts


class DFTMinimaStage(Stage):
    name = "dft-minima"

    def run(
        self, manifest: Manifest | None, config: dict[str, Any], context: StageContext
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine = get_qm_engine(
            config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {})
        )
        method = {
            "method_id": config.get("method", "r2scan3c"),
            **config.get("settings", {}),
            **config.get("method_settings", {}),
        }
        states = set(
            config.get(
                "states",
                [
                    "candidate",
                    "isolated_candidate",
                    "bare_candidate",
                    "hf_cluster",
                    "encounter_complex",
                    "reactant",
                    "reactant_complex",
                    "product",
                    "product_complex",
                    "adsorbed_reactant",
                    "adsorbed_product",
                    "ion_pair",
                    "reaction_intermediate_candidate",
                    "endpoint_minimum_candidate",
                ],
            )
        )
        update_species = bool(config.get("update_species", True))
        production_mode = str(context.global_config.get("mode", "")).lower() == "production"
        require_real_qm = bool(config.get("require_real_qm", production_mode))
        checkpoint_each_job = bool(config.get("checkpoint_each_job", production_mode))
        resume_completed_jobs = bool(
            config.get("resume_completed_jobs", production_mode)
        )
        checkpoint: Manifest | None = None
        checkpoint_path = out_dir / "manifest.json"
        if resume_completed_jobs and checkpoint_path.is_file():
            checkpoint = read_manifest(checkpoint_path)
        elif resume_completed_jobs and manifest.stage == self.name:
            checkpoint = manifest
        resumed_jobs = 0
        data_filters = dict(config.get("data_filters") or {})
        selection_artifact_type = config.get("selection_artifact_type")
        anchor_limits = {
            str(state): int(limit)
            for state, limit in dict(
                config.get("anchor_minima_per_state") or {}
            ).items()
        }
        if selection_artifact_type and anchor_limits:
            raise ValueError(
                "dft-minima cannot combine selection_artifact_type and "
                "anchor_minima_per_state"
            )
        selected_species_ids: set[str] | None = None
        candidate_endpoint_limit = config.get("candidate_endpoint_limit")
        if candidate_endpoint_limit is not None:
            if selection_artifact_type:
                raise ValueError(
                    "dft-minima cannot combine candidate_endpoint_limit and "
                    "selection_artifact_type"
                )
            selected_species_ids, selected_candidate_ids = (
                _candidate_endpoint_species_ids(
                    manifest,
                    int(candidate_endpoint_limit),
                    duplicate_rmsd_A=float(
                        config.get("candidate_duplicate_rmsd_A", 0.05)
                    ),
                )
            )
            out.metadata["dft_minima_selected_candidate_ids"] = (
                selected_candidate_ids
            )
            out.metadata["dft_minima_candidate_budget_distribution"] = (
                "round_robin_by_molecular_surface"
            )
        if selection_artifact_type:
            selected_species_ids = {
                str(species_id)
                for selection in manifest.latest_artifacts(
                    str(selection_artifact_type)
                )
                if selection.status.status == "success"
                for species_id in selection.data.get("selected_species_ids", [])
            }
            out.metadata["dft_minima_selection_artifact_type"] = str(
                selection_artifact_type
            )
        records: list[dict] = []
        optimized_records: list[dict] = []
        species_inputs = _preferred_species_inputs(
            manifest,
            states,
            include_optimized=bool(config.get("include_optimized_inputs", False)),
        )
        species_inputs.sort(
            key=lambda species: (
                str(species.data.get("state") or ""),
                species.data.get("discovery_electronic_energy_hartree") is None
                and species.data.get("relative_energy_kcal_mol") is None,
                float(
                    species.data.get(
                        "discovery_electronic_energy_hartree",
                        species.data.get(
                            "relative_energy_kcal_mol", float("inf")
                        ),
                    )
                ),
                canonical_species_id(species),
            )
        )
        species_inputs = [
            species
            for species in species_inputs
            if artifact_data_matches(species, data_filters)
        ]
        if anchor_limits:
            anchor_ids = _anchor_minimum_species_ids(
                species_inputs,
                anchor_limits,
                duplicate_rmsd_A=float(
                    config.get("anchor_duplicate_rmsd_A", 0.05)
                ),
            )
            selected_species_ids = (selected_species_ids or set()) | anchor_ids
            out.metadata["dft_minima_anchor_minima_per_state"] = anchor_limits
            out.metadata["dft_minima_anchor_species_ids"] = sorted(anchor_ids)
            out.metadata["dft_minima_anchor_budget_distribution"] = (
                "round_robin_by_molecular_surface"
            )
        if selected_species_ids is not None:
            species_inputs = [
                species
                for species in species_inputs
                if canonical_species_id(species) in selected_species_ids
            ]
            out.metadata["dft_minima_selected_species_count"] = len(
                species_inputs
            )
        per_state_limits = {
            str(state): int(limit)
            for state, limit in dict(
                config.get("max_species_per_state") or {}
            ).items()
        }
        if per_state_limits:
            state_counts: dict[str, int] = {}
            selected_inputs: list[Artifact] = []
            for species in species_inputs:
                state = str(species.data.get("state") or "")
                limit = per_state_limits.get(state)
                if limit is not None and state_counts.get(state, 0) >= limit:
                    continue
                selected_inputs.append(species)
                state_counts[state] = state_counts.get(state, 0) + 1
            species_inputs = selected_inputs
            out.metadata["dft_minima_max_species_per_state"] = per_state_limits
        maximum_species = config.get("max_species")
        if maximum_species is not None:
            maximum = int(maximum_species)
            if maximum < 1:
                raise ValueError("dft-minima max_species must be positive")
            species_inputs = species_inputs[:maximum]
            out.metadata["dft_minima_max_species"] = maximum
        for species in species_inputs:
            canonical_id = canonical_species_id(species)
            workdir = out_dir / canonical_id
            calc = _reusable_minimum_calculation(
                checkpoint,
                species,
                engine=str(config.get("engine", "dummy")),
                method=method,
                require_real_qm=require_real_qm,
            )
            if calc is None:
                calc = engine.optimize_frequency(
                    species,
                    method,
                    str(_next_minimum_workdir(workdir)),
                )
            else:
                resumed_jobs += 1
            calc.method = {**method, **(calc.method or {}), "stage": self.name}
            calc.provenance["dft_minima_request_fingerprint"] = (
                _minimum_request_fingerprint(
                    species,
                    engine=str(config.get("engine", "dummy")),
                    method=method,
                )
            )
            calc.qc.update(minimum_qc_from_freq(calc.data.get("n_imag"), species.data.get("state")))
            accepted, rejection_reasons = minimum_promotion_gate(
                calc, require_real_qm=require_real_qm
            )
            proton_state = (
                _proton_state_gate(species, calc)
                if accepted
                else {
                    "accepted": False,
                    "reasons": ["minimum_gate_failed_before_proton_state_check"],
                    "expected_state": species.data.get("state"),
                }
            )
            if accepted and not proton_state["accepted"]:
                accepted = False
                rejection_reasons = [*rejection_reasons, *proton_state["reasons"]]
            calc.qc.update(
                {
                    "minimum_accepted": accepted,
                    "minimum_rejection_reasons": rejection_reasons,
                    "proton_state_accepted": proton_state["accepted"],
                    "proton_state_gate": proton_state,
                    "require_real_qm": require_real_qm,
                }
            )
            calc.data["species_id"] = canonical_id
            calc.data["source_geometry_artifact_id"] = species.artifact_id
            out.add_artifact(calc)
            records.append(calc.model_dump())
            if accepted:
                final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
                optimized_records.append(
                    {
                        "species_id": canonical_id,
                        "state": species.data.get("state"),
                        "mol_id": species.data.get("mol_id"),
                        "site_id": species.data.get("site_id"),
                        "hf_n": species.data.get("hf_n"),
                        "source_geometry_artifact_id": species.artifact_id,
                        "optimized_xyz_path": final_xyz,
                        "dft_calc_id": calc.artifact_id,
                        "engine": (calc.method or {}).get("engine"),
                        "method_id": (calc.method or {}).get("method_id"),
                        "electronic_energy_hartree": calc.data.get("electronic_energy_hartree"),
                        "gibbs_298K_hartree": calc.data.get("gibbs_298K_hartree"),
                        "n_imag": calc.data.get("n_imag"),
                        "minimum_accepted": calc.qc.get("minimum_accepted"),
                        "proton_state_accepted": calc.qc.get("proton_state_accepted"),
                        "geometry_sane": calc.qc.get("geometry_sane"),
                        "fallback_dummy": calc.qc.get("fallback_dummy", False),
                        "real_orca_executed": calc.qc.get("real_orca_executed", False),
                        "real_qm_executed": calc.qc.get(
                            "real_qm_executed", calc.qc.get("real_orca_executed", False)
                        ),
                    }
                )
                for art in _optimized_species_artifacts(
                    species, calc, update_species=update_species
                ):
                    out.add_artifact(art)
            if checkpoint_each_job:
                out.metadata["dft_minima_resumed_jobs"] = resumed_jobs
                write_jsonl(records, out_dir / "dft_minima_calculation_records.jsonl")
                write_jsonl(optimized_records, out_dir / "optimized_species_records.jsonl")
                write_manifest(out, out_dir)
        out.metadata["dft_minima_resumed_jobs"] = resumed_jobs
        write_jsonl(records, out_dir / "dft_minima_calculation_records.jsonl")
        write_jsonl(optimized_records, out_dir / "optimized_species_records.jsonl")
        return out
