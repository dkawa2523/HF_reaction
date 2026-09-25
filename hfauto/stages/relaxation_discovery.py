"""Promote unbiased connectivity-changing relaxations to discovery candidates."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from hfauto.backends.conformer.crest import parse_crest_topology_endpoint
from hfauto.chemistry.connectivity import connectivity_changes
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.artifacts import canonical_species_id
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.chemistry import (
    BondChangeRecord,
    DrivingAtomPairRecord,
    ReactionCandidateRecord,
    ReactionCoordinateRecord,
    ReactionCoordinateTermRecord,
    ReactionTrialRecord,
)
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _existing_path(*values: object) -> Path | None:
    for value in values:
        if not value:
            continue
        candidate = Path(str(value))
        if candidate.is_file():
            return candidate
    return None


def _preopt_endpoints(manifest: Manifest) -> Iterable[dict[str, Any]]:
    for geometry in manifest.latest_artifacts("preopt_geometry"):
        if not geometry.parents:
            continue
        source = manifest.find(str(geometry.parents[0]))
        final_path = _existing_path(geometry.paths.get("xyz"))
        if source is None or final_path is None:
            continue
        initial_path = _existing_path(
            geometry.paths.get("source_xyz"),
            source.data.get("xyz_path"),
            source.paths.get("xyz"),
        )
        if initial_path is None:
            continue
        yield {
            "source": source,
            "initial_path": initial_path,
            "final_xyz": read_xyz(final_path),
            "evidence": geometry,
            "origin": "preopt",
            "electronic_energy_hartree": geometry.data.get(
                "electronic_energy_hartree"
            ),
        }


def _crest_topology_stop_endpoints(
    manifest: Manifest,
) -> Iterable[dict[str, Any]]:
    for failure in manifest.latest_artifacts("conformer"):
        if (
            failure.status.status != "failed"
            or failure.status.category != "crest_failed"
            or not failure.parents
        ):
            continue
        parsed = parse_crest_topology_endpoint(failure)
        source = manifest.find(str(failure.parents[0]))
        if source is None or parsed is None:
            continue
        initial_path = _existing_path(
            source.data.get("xyz_path"), source.paths.get("xyz")
        )
        if initial_path is None:
            continue
        yield {
            "source": source,
            "initial_path": initial_path,
            "final_xyz": parsed.xyz,
            "evidence": failure,
            "origin": "crest_initial_optimization",
            "electronic_energy_hartree": parsed.electronic_energy_hartree,
            "optimization_trajectory": parsed.trajectory_path,
        }


def _reaction_coordinate(
    formed_bonds: list[list[int]], broken_bonds: list[list[int]]
) -> ReactionCoordinateRecord:
    return ReactionCoordinateRecord(
        terms=[
            *(
                ReactionCoordinateTermRecord(
                    kind="distance",
                    atoms=tuple(pair),
                    coefficient=1.0,
                    label="relaxed_bond_breaking",
                )
                for pair in broken_bonds
            ),
            *(
                ReactionCoordinateTermRecord(
                    kind="distance",
                    atoms=tuple(pair),
                    coefficient=-1.0,
                    label="relaxed_bond_forming",
                )
                for pair in formed_bonds
            ),
        ]
    )


class RelaxationDiscoveryStage(Stage):
    """Normalize unbiased topology changes without claiming a barrier."""

    name = "relaxation-discovery"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("relaxation-discovery requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        endpoints = list(_preopt_endpoints(manifest))
        if bool(config.get("include_crest_topology_stops", True)):
            endpoints.extend(_crest_topology_stop_endpoints(manifest))
        records: list[dict[str, Any]] = []
        rejected: list[dict[str, Any]] = []
        observed_tokens: set[str] = set()

        for endpoint in endpoints:
            source: Artifact = endpoint["source"]
            initial_path: Path = endpoint["initial_path"]
            final_xyz: XYZ = endpoint["final_xyz"]
            connectivity = connectivity_changes(read_xyz(initial_path), final_xyz)
            if (
                connectivity.get("accepted") is not True
                or connectivity.get("topology_changed") is not True
            ):
                rejected.append(
                    {
                        "source_species_id": canonical_species_id(source),
                        "evidence_artifact_id": endpoint["evidence"].artifact_id,
                        "origin": endpoint["origin"],
                        "reason": connectivity.get("reason"),
                    }
                )
                continue
            token = fingerprint_dict(
                {
                    "source_species_id": canonical_species_id(source),
                    "initial_sha256": sha256_file(initial_path),
                    "final_symbols": final_xyz.symbols,
                    "final_coordinates": final_xyz.coords.tolist(),
                    "connectivity": connectivity,
                }
            )
            if token in observed_tokens:
                continue
            observed_tokens.add(token)
            final_path = write_xyz(
                final_xyz, out_dir / "endpoints" / f"relaxation_{token}.xyz"
            )
            source_id = canonical_species_id(source)
            trial_id = f"trial_relaxation_{token}"
            attempt_id = f"attempt_relaxation_{token}"
            product_id = f"spc_relaxation_product_{token}"
            candidate_id = f"candidate_relaxation_{token}"
            formed = list(connectivity["formed_bonds"])
            broken = list(connectivity["broken_bonds"])
            coordinate = _reaction_coordinate(formed, broken)
            trial = ReactionTrialRecord(
                trial_id=trial_id,
                source_species_id=source_id,
                associations=[
                    DrivingAtomPairRecord(atoms=tuple(pair)) for pair in formed
                ],
                dissociations=[
                    DrivingAtomPairRecord(atoms=tuple(pair)) for pair in broken
                ],
                reaction_coordinate=coordinate,
                driver_order=["unbiased_preopt"],
                mechanism_hint="unbiased_connectivity_relaxation",
                rationale=[
                    "An unconstrained low-level optimization changed the covalent graph."
                ],
                component_ids=[
                    str(item.get("component_id"))
                    for item in source.data.get("components", [])
                    if item.get("component_id")
                ],
                charge=int(source.data.get("charge", 0) or 0),
                multiplicity=int(source.data.get("multiplicity", 1) or 1),
                max_attempts=1,
            )
            trial_data = trial.model_dump(mode="json")
            trial_data.update(
                {
                    "source_xyz_path": str(initial_path),
                    "product_geometry_claimed": False,
                    "barrier_claimed": False,
                    "discovery_complete": True,
                }
            )
            trial_artifact = Artifact(
                artifact_id=trial_id,
                artifact_type="reaction_trial",
                parents=[source.artifact_id],
                paths={"reactant": str(initial_path)},
                data=trial_data,
                qc={"bounded_trial": True, "product_geometry_claimed": False},
            )
            attempt_data = {
                "attempt_id": attempt_id,
                "trial_id": trial_id,
                "source_species_id": source_id,
                "driver": "unbiased_preopt",
                "backend": endpoint["origin"],
                "success": True,
                "candidate_emitted": True,
                "endpoint_evidence": connectivity,
                "low_level_ts_validated": False,
                "low_level_irc_connected": False,
                "electronic_energy_hartree": endpoint.get(
                    "electronic_energy_hartree"
                ),
                "biased_energy_used_as_barrier": False,
                "barrier_claimed": False,
            }
            attempt = Artifact(
                artifact_id=attempt_id,
                artifact_type="reaction_discovery_attempt",
                parents=[trial_id, endpoint["evidence"].artifact_id],
                paths={
                    "reactant": str(initial_path),
                    "product": str(final_path),
                    **(
                        {
                            "optimization_trajectory": str(
                                endpoint["optimization_trajectory"]
                            )
                        }
                        if endpoint.get("optimization_trajectory")
                        else {}
                    ),
                },
                data=attempt_data,
                method={
                    "engine": endpoint["origin"],
                    "driver": "unbiased_preopt",
                },
                qc={
                    "composition_preserved": True,
                    "low_level_endpoint_distinct": True,
                    "structural_change_detected": True,
                    "biased_energy_used_as_barrier": False,
                },
            )
            product_data = {
                **source.data,
                "species_id": product_id,
                "source_species_id": product_id,
                "discovery_reactant_species_id": source_id,
                "reaction_trial_id": trial_id,
                "reaction_discovery_attempt_id": attempt_id,
                "discovery_electronic_energy_hartree": endpoint.get(
                    "electronic_energy_hartree"
                ),
                "state": "product_candidate",
                "xyz_path": str(final_path),
            }
            product = Artifact(
                artifact_id=product_id,
                artifact_type="species_preopt",
                parents=[source.artifact_id, attempt_id],
                paths={"xyz": str(final_path), "source_xyz": str(initial_path)},
                data=product_data,
                method={"stage": self.name, "engine": endpoint["origin"]},
                qc={
                    "preoptimized": True,
                    "unbiased_endpoint_optimization": True,
                    "low_level_endpoint_distinct": True,
                    "structural_change_detected": True,
                    "composition_preserved": True,
                    "biased_energy_used_as_barrier": False,
                },
            )
            candidate_record = ReactionCandidateRecord(
                candidate_id=candidate_id,
                trial_id=trial_id,
                reactant_species_id=source_id,
                product_species_id=product_id,
                driver="unbiased_preopt",
                bond_changes=[
                    *(
                        BondChangeRecord(kind="form", atoms=tuple(pair))
                        for pair in formed
                    ),
                    *(
                        BondChangeRecord(kind="break", atoms=tuple(pair))
                        for pair in broken
                    ),
                ],
                composition_preserved=True,
                charge=int(source.data.get("charge", 0) or 0),
                multiplicity=int(source.data.get("multiplicity", 1) or 1),
                low_level_endpoint_distinct=True,
                evidence_artifact_id=attempt_id,
            )
            candidate_data = candidate_record.model_dump(mode="json")
            candidate_data.update(
                {
                    "mechanism_family": "unbiased_connectivity_relaxation",
                    "reaction_coordinate": coordinate.model_dump(mode="json"),
                    "effectively_barrierless_low_level_hypothesis": True,
                    "barrier_claimed": False,
                }
            )
            candidate = Artifact(
                artifact_id=candidate_id,
                artifact_type="reaction_candidate",
                parents=[trial_id, attempt_id, source.artifact_id, product_id],
                data=candidate_data,
                qc={
                    "composition_preserved": True,
                    "low_level_endpoint_distinct": True,
                    "structural_change_detected": True,
                    "dft_minima_validated": False,
                    "reaction_promoted": False,
                    "activation_barrier_validated": False,
                },
            )
            out.extend([trial_artifact, attempt, product, candidate])
            records.append(candidate_data)

        summary = {
            "unbiased_relaxation_count": len(endpoints),
            "connectivity_change_candidate_count": len(records),
            "rejected_or_same_connectivity_count": len(rejected),
            "interpretation": (
                "low_level_effectively_barrierless_hypothesis_only"
                if records
                else "no_connectivity_changing_relaxation_observed"
            ),
            "activation_barrier_claimed": False,
        }
        summary_path = write_json(
            out_dir / "relaxation_discovery_summary.json", summary
        )
        out.add_artifact(
            Artifact(
                artifact_id="relaxation_discovery_summary",
                artifact_type="table",
                paths={"json": str(summary_path)},
                data=summary,
                qc={"activation_barrier_claimed": False},
                status=ArtifactStatus(status="success"),
            )
        )
        write_jsonl(records, out_dir / "relaxation_candidates.jsonl")
        write_jsonl(rejected, out_dir / "relaxation_rejected.jsonl")
        return out
