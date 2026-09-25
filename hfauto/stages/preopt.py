from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.chemistry.connectivity import connectivity_changes
from hfauto.chemistry.proton_transfer import (
    artifact_xyz_path,
    hf_endpoint_metrics_from_xyz,
    infer_spectator_hf_pairs,
    reaction_coordinate_atoms,
)
from hfauto.chemistry.reactions import is_proton_transfer_state
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import artifact_data_matches, canonical_species_id
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import (
    ensure_dir,
    read_manifest,
    write_jsonl,
    write_manifest,
)
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _preopt_request_fingerprint(
    species: Artifact, engine: str, method: dict[str, Any]
) -> str:
    return fingerprint_dict(
        {
            "species_id": canonical_species_id(species),
            "input_xyz_sha256": sha256_file(artifact_xyz_path(species)),
            "engine": engine,
            "method": method,
        }
    )


def _reusable_preopt_calculation(
    checkpoint: Manifest | None,
    request_fingerprint: str,
    *,
    require_real_qm: bool,
) -> Artifact | None:
    if checkpoint is None:
        return None
    for calculation in reversed(checkpoint.latest_artifacts("calculation")):
        if (
            calculation.status.status != "success"
            or calculation.provenance.get("preopt_request_fingerprint")
            != request_fingerprint
            or calculation.qc.get("fallback_dummy") is not False
            or (
                require_real_qm
                and calculation.qc.get("real_qm_executed") is not True
            )
        ):
            continue
        final_path = Path(str(calculation.paths.get("final_xyz") or ""))
        expected_hash = calculation.provenance.get("preopt_final_xyz_sha256")
        if (
            not final_path.is_file()
            or not expected_hash
            or sha256_file(final_path) != expected_hash
        ):
            continue
        reused = calculation.model_copy(deep=True)
        reused.provenance["reused_from_checkpoint"] = True
        return reused
    return None


def _next_preopt_workdir(base: Path) -> Path:
    attempt = 0
    while (candidate := base / f"attempt_{attempt:02d}").exists():
        attempt += 1
    return ensure_dir(candidate)


def _state_retention_gate(species: Artifact, calc: Artifact) -> tuple[bool, list[str]]:
    """Keep execution success separate from preservation of the requested state."""

    geometry_qc = calc.qc.get("geometry_qc", {}) or {}
    reasons: list[str] = []
    if calc.qc.get("geometry_sane") is False or geometry_qc.get("geometry_sane") is False:
        reasons.append("geometry_not_sane")
    state = str(species.data.get("state") or "")
    if (
        state in {"reactant_complex", "ion_pair"}
        and is_proton_transfer_state(species)
        and not reasons
    ):
        try:
            atoms = reaction_coordinate_atoms(species)
            spectator_pairs = infer_spectator_hf_pairs(artifact_xyz_path(species), atoms)
            final_xyz = calc.paths.get("final_xyz")
            if not final_xyz:
                raise FileNotFoundError("preoptimization final geometry is missing")
            endpoint = hf_endpoint_metrics_from_xyz(read_xyz(final_xyz), atoms, spectator_pairs)
            expected = "neutral_complex" if state == "reactant_complex" else "ion_pair"
            if endpoint["endpoint_class"] != expected:
                reasons.append(f"{state}_optimized_as_{endpoint['endpoint_class']}")
        except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
            if state == "reactant_complex" and geometry_qc.get(
                "proton_transferred_unintentionally",
                calc.qc.get("proton_transferred_unintentionally"),
            ):
                reasons.append("reactant_proton_state_not_retained")
            elif (
                state == "ion_pair"
                and geometry_qc.get(
                    "proton_transfer_confirmed", calc.qc.get("proton_transfer_confirmed")
                )
                is not True
            ):
                reasons.append("ion_pair_collapsed_during_preoptimization")
            else:
                reasons.append(f"proton_state_unverified:{exc}")
    return not reasons, reasons


def _revised_species_artifacts(
    species: Artifact, calc: Artifact, update_species: bool
) -> list[Artifact]:
    final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
    data = dict(species.data)
    data["species_id"] = canonical_species_id(species)
    data.setdefault("source_species_id", canonical_species_id(species))
    data["preopt_input_species_id"] = canonical_species_id(species)
    data["source_species_artifact_id"] = species.artifact_id
    data["xyz_path"] = str(final_xyz)
    data["preopt"] = {
        "calc_id": calc.artifact_id,
        "engine": (calc.method or {}).get("engine"),
        "method_id": (calc.method or {}).get("method_id"),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
    }
    qc = {
        "preoptimized": True,
        "preopt_calc_id": calc.artifact_id,
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "geometry_sane": calc.qc.get("geometry_sane")
        or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
        "hf_dissociated": calc.qc.get("hf_dissociated")
        or (calc.qc.get("geometry_qc", {}) or {}).get("hf_dissociated"),
        "proton_transferred_unintentionally": calc.qc.get("proton_transferred_unintentionally")
        or (calc.qc.get("geometry_qc", {}) or {}).get("proton_transferred_unintentionally"),
        "geometry_qc": calc.qc.get("geometry_qc", {}),
        "connectivity_change": calc.qc.get("connectivity_change", {}),
    }
    out = [
        Artifact(
            artifact_id=f"preopt_{species.artifact_id}",
            artifact_type="species_preopt",
            parents=[species.artifact_id, calc.artifact_id],
            paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
            data=data,
            method={"stage": "preopt", **(calc.method or {})},
            qc=qc,
        )
    ]
    if update_species:
        # Same artifact_id, later revision. Manifest.latest_artifacts("species")
        # will return this preoptimized geometry while original lineage remains.
        out.append(
            Artifact(
                artifact_id=species.artifact_id,
                artifact_type="species",
                parents=[species.artifact_id, calc.artifact_id],
                paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
                data=data,
                method={"stage": "preopt", **(calc.method or {})},
                qc=qc,
            )
        )
    return out


class PreoptStage(Stage):
    """Low-level geometry preoptimization.

    The stage stores both a dedicated ``species_preopt`` Artifact and, when
    ``update_species`` is true, a latest-revision ``species`` Artifact with the
    same ID. This supports both explicit provenance and simple downstream use.
    """

    name = "preopt"

    def run(
        self, manifest: Manifest | None, config: dict[str, Any], context: StageContext
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        engine = get_qm_engine(
            config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {})
        )
        method = {
            "method_id": config.get("method", "gfn2-xtb"),
            **config.get("settings", {}),
            **config.get("method_settings", {}),
        }
        update_species = bool(config.get("update_species", True))
        production_mode = (
            str(context.global_config.get("mode", "")).lower() == "production"
        )
        require_real_qm = bool(config.get("require_real_qm", production_mode))
        checkpoint_each_job = bool(config.get("checkpoint_each_job", production_mode))
        resume_completed_jobs = bool(
            config.get("resume_completed_jobs", production_mode)
        )
        checkpoint_path = out_dir / "manifest.json"
        checkpoint = (
            read_manifest(checkpoint_path)
            if resume_completed_jobs and checkpoint_path.is_file()
            else None
        )
        working = (
            Manifest.merge(
                [manifest, checkpoint],
                run_id=context.run_id,
                stage=f"{self.name}-resume",
                metadata={**manifest.metadata, **checkpoint.metadata},
            )
            if checkpoint is not None
            else manifest
        )
        out = working.carry_forward(self.name)
        require_connectivity_retention = bool(
            config.get("require_connectivity_retention", production_mode)
        )
        states = {str(value) for value in config.get("states", [])}
        data_filters = dict(config.get("data_filters") or {})
        required_source_keys = {
            str(key) for key in config.get("required_source_data_keys", [])
        }
        records: list[dict] = []
        preopt_species_records: list[dict] = []
        resumed_jobs = 0
        # Always fingerprint the upstream geometries. Checkpoint revisions are
        # outputs and must not recursively become new preopt inputs.
        species_inputs = manifest.latest_artifacts("species")
        for species in species_inputs:
            if states and str(species.data.get("state")) not in states:
                continue
            if not artifact_data_matches(species, data_filters):
                continue
            if not all(
                species.data.get(key) is not None
                for key in required_source_keys
            ):
                continue
            request_fingerprint = _preopt_request_fingerprint(
                species, str(config.get("engine", "dummy")), method
            )
            calc = _reusable_preopt_calculation(
                checkpoint,
                request_fingerprint,
                require_real_qm=require_real_qm,
            )
            if calc is None:
                calc = engine.optimize_frequency(
                    species,
                    method,
                    str(_next_preopt_workdir(out_dir / species.artifact_id)),
                )
            else:
                resumed_jobs += 1
            calc.method = {**(calc.method or {}), "stage": self.name}
            calc.provenance["preopt_request_fingerprint"] = request_fingerprint
            final_for_hash = Path(str(calc.paths.get("final_xyz") or ""))
            if calc.status.status == "success" and final_for_hash.is_file():
                calc.provenance["preopt_final_xyz_sha256"] = sha256_file(
                    final_for_hash
                )
            if calc.status.status == "success":
                state_retained, state_reasons = _state_retention_gate(species, calc)
                if require_real_qm and (
                    calc.qc.get("real_qm_executed") is not True
                    or calc.qc.get("fallback_dummy") is not False
                ):
                    state_retained = False
                    state_reasons.append("real_preoptimization_required")
                try:
                    final_xyz = calc.paths.get("final_xyz")
                    if not final_xyz:
                        raise FileNotFoundError(
                            "preoptimization final geometry is missing"
                        )
                    connectivity = connectivity_changes(
                        read_xyz(artifact_xyz_path(species)), read_xyz(final_xyz)
                    )
                except (OSError, TypeError, ValueError) as exc:
                    connectivity = {
                        "accepted": False,
                        "topology_changed": None,
                        "reason": f"connectivity_assessment_failed:{exc}",
                    }
                if require_connectivity_retention and (
                    connectivity.get("accepted") is not True
                    or connectivity.get("topology_changed") is True
                ):
                    state_retained = False
                    state_reasons.append(str(connectivity["reason"]))
                calc.qc.update(
                    {
                        "preopt_state_retained": state_retained,
                        "preopt_promotion_accepted": state_retained,
                        "preopt_promotion_reasons": state_reasons,
                        "connectivity_change": connectivity,
                        "require_connectivity_retention": (
                            require_connectivity_retention
                        ),
                    }
                )
                final_xyz = calc.paths.get("final_xyz") or species.data.get(
                    "xyz_path"
                )
                rec = {
                    "species_id": species.artifact_id,
                    "state": species.data.get("state"),
                    "mol_id": species.data.get("mol_id"),
                    "site_id": species.data.get("site_id"),
                    "hf_n": species.data.get("hf_n"),
                    "source_species_xyz": species.data.get("xyz_path"),
                    "preopt_xyz_path": final_xyz,
                    "preopt_calc_id": calc.artifact_id,
                    "engine": calc.method.get("engine") if calc.method else None,
                    "method_id": calc.method.get("method_id") if calc.method else None,
                    "electronic_energy_hartree": calc.data.get("electronic_energy_hartree"),
                    "geometry_sane": calc.qc.get("geometry_sane")
                    or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
                    "hf_dissociated": calc.qc.get("hf_dissociated")
                    or (calc.qc.get("geometry_qc", {}) or {}).get("hf_dissociated"),
                    "proton_transferred_unintentionally": calc.qc.get(
                        "proton_transferred_unintentionally"
                    )
                    or (calc.qc.get("geometry_qc", {}) or {}).get(
                        "proton_transferred_unintentionally"
                    ),
                    "fallback_dummy": calc.qc.get("fallback_dummy", False),
                    "state_retained": state_retained,
                    "promotion_accepted": state_retained,
                    "promotion_reasons": state_reasons,
                    "connectivity_change": connectivity,
                }
                preopt_species_records.append(rec)
                if state_retained:
                    for art in _revised_species_artifacts(
                        species, calc, update_species=update_species
                    ):
                        out.add_artifact(art)
                out.add_artifact(
                    Artifact(
                        artifact_id=f"preopt_geom_{species.artifact_id}",
                        artifact_type="preopt_geometry",
                        parents=[species.artifact_id, calc.artifact_id],
                        paths={"xyz": str(final_xyz)},
                        data=rec,
                        qc={
                            "preoptimized": True,
                            "state_retained": state_retained,
                            "promotion_accepted": state_retained,
                            "promotion_reasons": state_reasons,
                            "geometry_sane": rec["geometry_sane"],
                            "hf_dissociated": rec["hf_dissociated"],
                            "proton_transferred_unintentionally": rec[
                                "proton_transferred_unintentionally"
                            ],
                            "fallback_dummy": rec["fallback_dummy"],
                            "connectivity_change": connectivity,
                        },
                        status=ArtifactStatus(
                            status="success" if state_retained else "partial",
                            category=(
                                None
                                if state_retained
                                else "connectivity_changed"
                                if connectivity.get("topology_changed") is True
                                else "state_not_retained"
                            ),
                            reason=None if state_retained else "; ".join(state_reasons),
                        ),
                    )
                )
            out.add_artifact(calc)
            records.append(calc.model_dump())
            if checkpoint_each_job:
                out.metadata["preopt_resumed_jobs"] = resumed_jobs
                write_jsonl(records, out_dir / "preopt_calculation_records.jsonl")
                write_jsonl(
                    preopt_species_records,
                    out_dir / "preopt_species_records.jsonl",
                )
                write_manifest(out, out_dir)
        out.metadata["preopt_resumed_jobs"] = resumed_jobs
        out.metadata["preopt_required_source_data_keys"] = sorted(
            required_source_keys
        )
        write_jsonl(records, out_dir / "preopt_calculation_records.jsonl")
        write_jsonl(preopt_species_records, out_dir / "preopt_species_records.jsonl")
        write_manifest(out, out_dir)
        return out
