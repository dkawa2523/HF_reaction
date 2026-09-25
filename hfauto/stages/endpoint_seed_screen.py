"""Screen endpoint orientation seeds before expensive opt/freq calculations."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from hfauto.chemistry.complex_seed_ensemble import select_diverse_seed_records
from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import calculations_by_species


class EndpointSeedScreenStage(Stage):
    """Select DFT inputs; no preoptimization energy is reported as chemistry."""

    name = "endpoint-seed-screen"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("endpoint-seed-screen requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        calculations = calculations_by_species(manifest)
        candidates = preferred_species_by_id(
            manifest,
            states=("endpoint_minimum_candidate",),
            artifact_types=("species_preopt", "species"),
        )
        grouped: defaultdict[tuple[str, str, str], list[dict[str, Any]]] = (
            defaultdict(list)
        )
        for species_id, species in candidates.items():
            if int(species.data.get("endpoint_seed_generation", 0) or 0) != 1:
                continue
            calculation = calculations.get(species_id)
            if calculation is None or (calculation.method or {}).get("stage") != "preopt":
                continue
            reaction_ids = species.data.get("reaction_ids") or []
            for reaction_id in reaction_ids:
                grouped[
                    (
                        str(reaction_id),
                        str(species.data.get("endpoint_role") or "unknown"),
                        str(species.data.get("source_species_id") or ""),
                    )
                ].append(
                    {
                        "species_id": species_id,
                        "species_artifact_id": species.artifact_id,
                        "calculation_artifact_id": calculation.artifact_id,
                        "xyz_path": species.data.get("xyz_path")
                        or species.paths.get("xyz"),
                        "energy_hartree": calculation.data.get(
                            "electronic_energy_hartree"
                        ),
                        "accepted": bool(
                            calculation.status.status == "success"
                            and calculation.qc.get("geometry_converged") is True
                            and calculation.qc.get("geometry_sane") is True
                            and calculation.qc.get("fallback_dummy") is False
                            and calculation.qc.get("preopt_promotion_accepted") is True
                        ),
                    }
                )
        artifacts: list[dict[str, Any]] = []
        for (reaction_id, role, root_species_id), records in sorted(grouped.items()):
            selected = select_diverse_seed_records(
                records,
                maximum_selected=int(config.get("maximum_selected_per_role", 2)),
                energy_window_kcal_mol=float(
                    config.get("energy_window_kcal_mol", 3.0)
                ),
                minimum_distance_rmsd_A=float(
                    config.get("minimum_distance_rmsd_A", 0.03)
                ),
            )
            identifier = "endpoint_seed_selection_" + fingerprint_dict(
                {
                    "reaction_id": reaction_id,
                    "role": role,
                    "root_species_id": root_species_id,
                    "candidate_ids": sorted(
                        record["species_id"] for record in records
                    ),
                }
            )
            data = {
                "reaction_id": reaction_id,
                "endpoint_role": role,
                "root_species_id": root_species_id,
                "candidate_records": records,
                "selected_species_ids": [
                    record["species_id"] for record in selected
                ],
                "selected_species_artifact_ids": [
                    record["species_artifact_id"] for record in selected
                ],
                "selection_method": (
                    "preopt_energy_window_then_maximin_permutation_invariant_"
                    "distance_rmsd"
                ),
                "preopt_energy_is_final_chemistry": False,
                "thresholds": {
                    "maximum_selected_per_role": int(
                        config.get("maximum_selected_per_role", 2)
                    ),
                    "energy_window_kcal_mol": float(
                        config.get("energy_window_kcal_mol", 3.0)
                    ),
                    "minimum_distance_rmsd_A": float(
                        config.get("minimum_distance_rmsd_A", 0.03)
                    ),
                },
            }
            artifact = Artifact(
                artifact_id=identifier,
                artifact_type="endpoint_seed_selection",
                parents=[
                    *[record["species_artifact_id"] for record in records],
                    *[record["calculation_artifact_id"] for record in records],
                ],
                data=data,
                qc={
                    "selection_ready": bool(selected),
                    "selected_count": len(selected),
                    "dft_minimum_claimed": False,
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if selected else "partial",
                    category=None if selected else "endpoint_seed_screen_empty",
                    reason=None if selected else "no accepted preoptimized seeds",
                ),
            )
            out.add_artifact(artifact)
            artifacts.append(artifact.model_dump(mode="json"))
        write_jsonl(artifacts, out_dir / "endpoint_seed_selections.jsonl")
        out.metadata["endpoint_seed_selection_count"] = len(artifacts)
        out.metadata["endpoint_seed_selected_count"] = sum(
            artifact["qc"]["selected_count"] for artifact in artifacts
        )
        return out
