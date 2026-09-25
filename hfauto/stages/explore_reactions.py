"""Execute bounded single-ended trials and publish low-level candidates."""

from __future__ import annotations

import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.backends.reaction_discovery.base import DiscoveryResult
from hfauto.backends.registry import get_reaction_discovery_engine
from hfauto.chemistry.basin_identity import (
    permutation_invariant_graph_rmsd,
)
from hfauto.chemistry.connectivity import covalent_adjacency
from hfauto.chemistry.reactions import (
    validate_reaction_coordinate_between_geometries,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import (
    ensure_dir,
    read_manifest,
    write_json,
    write_jsonl,
    write_manifest,
)
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.chemistry import ReactionCandidateRecord, ReactionTrialRecord
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _endpoint_score(
    trial: ReactionTrialRecord, source_path: str, endpoint_path: str
) -> dict[str, Any]:
    source = read_xyz(source_path)
    endpoint = read_xyz(endpoint_path)
    if Counter(source.symbols) != Counter(endpoint.symbols):
        return {
            "composition_preserved": False,
            "changed_graph": False,
            "requested_changes_satisfied": 0,
            "rmsd_A": None,
        }
    if source.symbols != endpoint.symbols:
        return {
            "composition_preserved": True,
            "changed_graph": False,
            "requested_changes_satisfied": 0,
            "rmsd_A": None,
            "reason": "atom_order_changed",
        }
    source_graph = covalent_adjacency(source)
    endpoint_graph = covalent_adjacency(endpoint)
    changed_edges = int(np.sum(source_graph != endpoint_graph) // 2)
    graph_rmsd = permutation_invariant_graph_rmsd(source_path, endpoint_path)
    atom_mapped_graph_changed = changed_edges > 0
    graph_isomorphic = graph_rmsd is not None
    associations_satisfied = sum(
        bool(endpoint_graph[pair.atoms]) and not bool(source_graph[pair.atoms])
        for pair in trial.associations
    )
    dissociations_satisfied = sum(
        bool(source_graph[pair.atoms]) and not bool(endpoint_graph[pair.atoms])
        for pair in trial.dissociations
    )
    coordinate_evidence = None
    selected_coordinate = None
    if trial.reaction_coordinate.terms:
        declared_coordinate = trial.reaction_coordinate.model_dump(mode="json")
        coordinate_evidence = validate_reaction_coordinate_between_geometries(
            {"reaction_coordinate": declared_coordinate},
            source_path,
            endpoint_path,
        )
        selected_coordinate = declared_coordinate
        coordinate_evidence["direction"] = "forward"
        if coordinate_evidence["accepted"] is not True:
            reversed_coordinate = {
                **declared_coordinate,
                "terms": [
                    {
                        **term,
                        "coefficient": -float(term.get("coefficient", 1.0)),
                    }
                    for term in declared_coordinate["terms"]
                ],
            }
            reverse_evidence = validate_reaction_coordinate_between_geometries(
                {"reaction_coordinate": reversed_coordinate},
                source_path,
                endpoint_path,
            )
            if reverse_evidence["accepted"] is True:
                coordinate_evidence = {**reverse_evidence, "direction": "reverse"}
                selected_coordinate = reversed_coordinate
    return {
        "composition_preserved": True,
        # A mapped edge difference caused only by exchanging equivalent atoms
        # is not a chemically new connectivity graph.
        "changed_graph": atom_mapped_graph_changed and not graph_isomorphic,
        "atom_mapped_graph_changed": atom_mapped_graph_changed,
        "permutation_invariant_graph_isomorphic": graph_isomorphic,
        "changed_edge_count": changed_edges,
        "requested_changes_satisfied": associations_satisfied
        + dissociations_satisfied,
        "requested_change_count": (
            len(trial.associations)
            + len(trial.dissociations)
            + len(trial.reaction_coordinate.terms)
        ),
        "reaction_coordinate_evidence": coordinate_evidence,
        "selected_reaction_coordinate": selected_coordinate,
        "rmsd_A": graph_rmsd,
    }


def _select_product_endpoint(
    trial: ReactionTrialRecord,
    source_path: str,
    result: DiscoveryResult,
    *,
    same_graph_product_rmsd_A: float,
) -> tuple[str | None, dict[str, Any]]:
    paths = [
        path
        for path in (
            result.product_xyz_path,
            result.reactant_endpoint_xyz_path,
        )
        if path
    ]
    scored = [(path, _endpoint_score(trial, source_path, path)) for path in paths]
    scored.sort(
        key=lambda item: (
            int(item[1].get("composition_preserved") is True),
            int(item[1].get("changed_graph") is True),
            int(item[1].get("requested_changes_satisfied", 0)),
            float(item[1].get("rmsd_A") or 0.0),
        ),
        reverse=True,
    )
    if not scored:
        return None, {"reason": "discovery_result_has_no_endpoint"}
    path, evidence = scored[0]
    requested_change_count = (
        len(trial.associations)
        + len(trial.dissociations)
        + len(trial.reaction_coordinate.terms)
    )
    graph_changed = evidence.get("changed_graph") is True
    declared_coordinate_changed = bool(
        (evidence.get("reaction_coordinate_evidence") or {}).get("accepted")
        is True
    )
    undeclared_geometry_changed = bool(
        requested_change_count == 0
        and float(evidence.get("rmsd_A") or 0.0)
        >= float(same_graph_product_rmsd_A)
    )
    distinct = bool(
        evidence.get("composition_preserved") is True
        and (
            graph_changed
            or declared_coordinate_changed
            or undeclared_geometry_changed
        )
    )
    rejection_reason = None
    if (
        not distinct
        and requested_change_count > 0
        and evidence.get("atom_mapped_graph_changed") is True
        and evidence.get("permutation_invariant_graph_isomorphic") is True
    ):
        rejection_reason = "equivalent_atom_permutation_only"
    elif not distinct and trial.reaction_coordinate.terms:
        rejection_reason = "requested_reaction_coordinate_not_realized"
    elif not distinct and requested_change_count > 0:
        rejection_reason = "requested_bond_changes_not_realized"
    elif not distinct:
        rejection_reason = "same_graph_change_below_threshold"
    return (path if distinct else None), {
        **evidence,
        **({"reason": rejection_reason} if rejection_reason else {}),
        "low_level_endpoint_distinct": distinct,
        "structural_change_detected": distinct,
        "structural_change_basis": (
            "covalent_graph_change"
            if graph_changed
            else "declared_reaction_coordinate_change"
            if declared_coordinate_changed
            else "undeclared_same_graph_geometry_change"
            if undeclared_geometry_changed
            else None
        ),
        "all_endpoint_scores": [
            {"path": candidate_path, **candidate_evidence}
            for candidate_path, candidate_evidence in scored
        ],
    }


def _hash_paths(paths: dict[str, str]) -> dict[str, str]:
    return {
        key: sha256_file(path)
        for key, path in paths.items()
        if Path(path).is_file()
    }


def _attempt_request_fingerprint(
    trial: ReactionTrialRecord,
    source_path: str,
    *,
    driver: str,
    engine: str,
    settings: dict[str, Any],
) -> str:
    return fingerprint_dict(
        {
            "trial": trial.model_dump(mode="json"),
            "source_xyz_sha256": sha256_file(source_path),
            "driver": driver,
            "engine": engine,
            "settings": settings,
        }
    )


def _attempt_files_unchanged(attempt: Artifact) -> bool:
    expected = dict(attempt.data.get("path_hashes") or {})
    for name, digest in expected.items():
        path = attempt.paths.get(name)
        if not path or not Path(path).is_file() or sha256_file(path) != digest:
            return False
    return True


class ExploreReactionsStage(Stage):
    """Run real discovery backends; no product means no reaction candidate."""

    name = "explore-reactions"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("explore-reactions requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        production_mode = (
            str(context.global_config.get("mode", "")).lower() == "production"
        )
        checkpoint_each_attempt = bool(
            config.get("checkpoint_each_attempt", production_mode)
        )
        resume_completed_attempts = bool(
            config.get("resume_completed_attempts", production_mode)
        )
        checkpoint_path = out_dir / "manifest.json"
        checkpoint = (
            read_manifest(checkpoint_path)
            if resume_completed_attempts and checkpoint_path.is_file()
            else None
        )
        working = (
            Manifest.merge(
                [manifest, checkpoint],
                run_id=manifest.run_id,
                stage=f"{self.name}-resume",
                metadata={**manifest.metadata, **checkpoint.metadata},
            )
            if checkpoint is not None
            else manifest
        )
        out = working.carry_forward(self.name)
        engine = get_reaction_discovery_engine(
            config.get("engine", "readuct"),
            **dict(config.get("engine_settings") or {}),
        )
        source_by_id = preferred_species_by_id(
            working, artifact_types=("species", "species_preopt")
        )
        attempts = [
            dict(item.data)
            for item in working.latest_artifacts("reaction_discovery_attempt")
        ]
        candidates = [
            dict(item.data)
            for item in working.latest_artifacts("reaction_candidate")
        ]
        attempt_outcomes: Counter[str] = Counter()
        for item in working.latest_artifacts("reaction_discovery_attempt"):
            attempt_outcomes[
                str(item.status.reason or item.status.category or item.status.status)
            ] += 1
        max_total_attempts = int(config.get("max_total_attempts", 100))
        max_stage_walltime_s = float(
            config.get("max_stage_walltime_s", float("inf"))
        )
        if max_total_attempts < 1 or max_stage_walltime_s <= 0.0:
            raise ValueError("reaction-discovery execution budgets must be positive")
        stage_started = time.monotonic()
        budget_exhausted_reason: str | None = None
        stop_after_first = bool(config.get("stop_after_first_candidate", True))
        same_graph_rmsd = float(config.get("same_graph_product_rmsd_A", 0.20))
        settings = dict(config.get("settings") or {})
        completed_requests = {
            str(item.data.get("request_fingerprint"))
            for item in working.latest_artifacts("reaction_discovery_attempt")
            if item.data.get("request_fingerprint")
            and _attempt_files_unchanged(item)
        }
        candidate_trials = {
            str(item.data.get("trial_id"))
            for item in working.latest_artifacts("reaction_candidate")
        }

        def checkpoint_now() -> None:
            out.metadata.update(
                {
                    "reaction_discovery_attempt_count": len(attempts),
                    "reaction_discovery_candidate_count": len(candidates),
                    "reaction_discovery_budget": {
                        "max_total_attempts": max_total_attempts,
                        "max_stage_walltime_s": (
                            None
                            if max_stage_walltime_s == float("inf")
                            else max_stage_walltime_s
                        ),
                    },
                }
            )
            write_jsonl(
                attempts, out_dir / "reaction_discovery_attempts.jsonl"
            )
            write_jsonl(candidates, out_dir / "reaction_candidates.jsonl")
            write_manifest(out, out_dir)

        for trial_artifact in working.latest_artifacts("reaction_trial"):
            if len(attempts) >= max_total_attempts:
                budget_exhausted_reason = "max_total_attempts_reached"
                break
            if time.monotonic() - stage_started >= max_stage_walltime_s:
                budget_exhausted_reason = "max_stage_walltime_reached"
                break
            if trial_artifact.status.status != "success":
                continue
            if trial_artifact.data.get("discovery_complete") is True:
                continue
            trial = ReactionTrialRecord.model_validate(trial_artifact.data)
            if stop_after_first and trial.trial_id in candidate_trials:
                continue
            source = source_by_id.get(trial.source_species_id)
            if source is None:
                out.add_artifact(
                    Artifact.failure(
                        f"attempt_missing_source_{trial.trial_id}",
                        "reaction_discovery_attempt",
                        f"source species not found: {trial.source_species_id}",
                        category="discovery_source_missing",
                        parents=[trial_artifact.artifact_id],
                    )
                )
                continue
            source_path = str(source.data.get("xyz_path") or source.paths["xyz"])
            for attempt_index, driver in enumerate(
                trial.driver_order[: trial.max_attempts]
            ):
                if len(attempts) >= max_total_attempts:
                    budget_exhausted_reason = "max_total_attempts_reached"
                    break
                if time.monotonic() - stage_started >= max_stage_walltime_s:
                    budget_exhausted_reason = "max_stage_walltime_reached"
                    break
                attempt_id = f"attempt_{trial.trial_id}_{driver}_{attempt_index}"
                request_fingerprint = _attempt_request_fingerprint(
                    trial,
                    source_path,
                    driver=driver,
                    engine=engine.name,
                    settings=settings,
                )
                if request_fingerprint in completed_requests:
                    continue
                result = engine.explore(
                    trial,
                    source,
                    settings,
                    out_dir / trial.trial_id / f"{attempt_index:02d}_{driver}",
                    driver=driver,
                )
                product_path, endpoint_evidence = _select_product_endpoint(
                    trial,
                    source_path,
                    result,
                    same_graph_product_rmsd_A=same_graph_rmsd,
                )
                successful_candidate = bool(result.success and product_path)
                path_hashes = _hash_paths(result.paths)
                attempt_data = {
                    "attempt_id": attempt_id,
                    "trial_id": trial.trial_id,
                    "source_species_id": trial.source_species_id,
                    "driver": driver,
                    "backend": engine.name,
                    "success": result.success,
                    "candidate_emitted": successful_candidate,
                    "failure_reason": result.failure_reason,
                    "endpoint_evidence": endpoint_evidence,
                    "low_level_ts_validated": result.low_level_ts_validated,
                    "low_level_irc_connected": result.low_level_irc_connected,
                    "imaginary_mode_count": result.imaginary_mode_count,
                    "electronic_energy_hartree": result.electronic_energy_hartree,
                    "biased_energy_used_as_barrier": result.biased_energy_used_as_barrier,
                    "path_hashes": path_hashes,
                    "request_fingerprint": request_fingerprint,
                    **result.data,
                }
                attempt = Artifact(
                    artifact_id=attempt_id,
                    artifact_type="reaction_discovery_attempt",
                    parents=[trial_artifact.artifact_id, source.artifact_id],
                    paths=dict(result.paths),
                    data=attempt_data,
                    method={
                        "engine": engine.name,
                        "driver": driver,
                        "program": "XTB",
                        "method_family": str(
                            config.get("settings", {}).get("method_family", "GFN2")
                        ),
                    },
                    qc={
                        "real_discovery_executed": result.success,
                        "composition_preserved": endpoint_evidence.get(
                            "composition_preserved"
                        ),
                        "low_level_endpoint_distinct": endpoint_evidence.get(
                            "low_level_endpoint_distinct"
                        ),
                        "structural_change_detected": endpoint_evidence.get(
                            "structural_change_detected"
                        ),
                        "biased_energy_used_as_barrier": False,
                    },
                    status=ArtifactStatus(
                        status="success" if successful_candidate else "partial",
                        category=(
                            None
                            if successful_candidate
                            else "no_distinct_low_level_product"
                            if result.success
                            else "reaction_discovery_failed"
                        ),
                        reason=result.failure_reason or endpoint_evidence.get("reason"),
                    ),
                )
                out.add_artifact(attempt)
                attempts.append(attempt_data)
                attempt_outcomes[
                    str(
                        attempt.status.reason
                        or attempt.status.category
                        or attempt.status.status
                    )
                ] += 1
                if not successful_candidate or product_path is None:
                    if checkpoint_each_attempt:
                        checkpoint_now()
                    continue

                token = fingerprint_dict(
                    {
                        "trial_id": trial.trial_id,
                        "driver": driver,
                        "product_sha256": sha256_file(product_path),
                    }
                )[:20]
                product_id = f"spc_discovered_product_{token}"
                product_data = {
                    **source.data,
                    "species_id": product_id,
                    "source_species_id": product_id,
                    "discovery_reactant_species_id": trial.source_species_id,
                    "reaction_trial_id": trial.trial_id,
                    "reaction_discovery_attempt_id": attempt_id,
                    "discovery_electronic_energy_hartree": (
                        result.electronic_energy_hartree
                    ),
                    "state": "product_candidate",
                    "xyz_path": product_path,
                    "preopt": {
                        "engine": engine.name,
                        "program": "XTB",
                        "method_family": "GFN2",
                        "unbiased_endpoint_optimization": True,
                    },
                }
                product = Artifact(
                    artifact_id=product_id,
                    artifact_type="species_preopt",
                    parents=[source.artifact_id, trial_artifact.artifact_id, attempt_id],
                    paths={"xyz": product_path, "source_xyz": source_path},
                    data=product_data,
                    method={"stage": self.name, "engine": engine.name, "program": "XTB"},
                    qc={
                        "preoptimized": True,
                        "unbiased_endpoint_optimization": True,
                        "low_level_endpoint_distinct": True,
                        "structural_change_detected": True,
                        "composition_preserved": True,
                        "biased_energy_used_as_barrier": False,
                    },
                )
                candidate_id = f"candidate_{token}"
                record = ReactionCandidateRecord(
                    candidate_id=candidate_id,
                    trial_id=trial.trial_id,
                    reactant_species_id=trial.source_species_id,
                    product_species_id=product_id,
                    driver=driver,
                    bond_changes=trial.bond_changes(),
                    composition_preserved=True,
                    charge=trial.charge,
                    multiplicity=trial.multiplicity,
                    low_level_endpoint_distinct=True,
                    low_level_ts_validated=result.low_level_ts_validated,
                    low_level_irc_connected=result.low_level_irc_connected,
                    evidence_artifact_id=attempt_id,
                )
                candidate_data = record.model_dump(mode="json")
                candidate_data.update(
                    {
                        "mechanism_family": trial.mechanism_hint,
                        "rationale": trial.rationale,
                        "reaction_coordinate": (
                            endpoint_evidence.get("selected_reaction_coordinate")
                            if trial.reaction_coordinate.terms
                            else {
                                "terms": [
                                    {
                                        "kind": "distance",
                                        "atoms": list(pair.atoms),
                                        "coefficient": coefficient,
                                        "label": pair.label,
                                    }
                                    for pairs, coefficient in (
                                        (trial.dissociations, 1.0),
                                        (trial.associations, -1.0),
                                    )
                                    for pair in pairs
                                ]
                            }
                        ),
                    }
                )
                candidate = Artifact(
                    artifact_id=candidate_id,
                    artifact_type="reaction_candidate",
                    parents=[trial_artifact.artifact_id, attempt_id, source.artifact_id, product_id],
                    data=candidate_data,
                    qc={
                        "composition_preserved": True,
                        "low_level_endpoint_distinct": True,
                        "structural_change_detected": True,
                        "structural_change_basis": endpoint_evidence.get(
                            "structural_change_basis"
                        ),
                        "dft_minima_validated": False,
                        "reaction_promoted": False,
                    },
                )
                out.extend([product, candidate])
                candidates.append(candidate_data)
                candidate_trials.add(trial.trial_id)
                if checkpoint_each_attempt:
                    checkpoint_now()
                if stop_after_first:
                    break

        summary = {
            "reaction_trial_count": len(
                [
                    item
                    for item in working.latest_artifacts("reaction_trial")
                    if item.status.status == "success"
                ]
            ),
            "discovery_attempt_count": len(attempts),
            "reaction_candidate_count": len(candidates),
            "max_total_attempts": max_total_attempts,
            "max_stage_walltime_s": (
                None
                if max_stage_walltime_s == float("inf")
                else max_stage_walltime_s
            ),
            "budget_exhausted_reason": budget_exhausted_reason,
            "biased_energy_accepted_as_barrier": False,
            "outcome": (
                "reaction_candidates_found"
                if candidates
                else "no_distinct_low_level_product"
            ),
            "path_search_eligible": bool(candidates),
            "attempt_outcomes": dict(sorted(attempt_outcomes.items())),
        }
        summary_path = write_json(out_dir / "reaction_discovery_summary.json", summary)
        out.add_artifact(
            Artifact(
                artifact_id="reaction_discovery_summary",
                artifact_type="table",
                paths={"json": str(summary_path)},
                data=summary,
                qc={
                    "bounded_execution": len(attempts) <= max_total_attempts,
                    "execution_budget_exhausted": (
                        budget_exhausted_reason is not None
                    ),
                    "path_search_eligible": bool(candidates),
                    "scientific_negative_result": bool(attempts and not candidates),
                },
            )
        )
        write_jsonl(attempts, out_dir / "reaction_discovery_attempts.jsonl")
        write_jsonl(candidates, out_dir / "reaction_candidates.jsonl")
        write_manifest(out, out_dir)
        return out
