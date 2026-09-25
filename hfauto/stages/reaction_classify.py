"""Publish the final evidence-based class of each molecular reaction."""

from __future__ import annotations

from typing import Any

from hfauto.chemistry.minima import is_accepted_optimized_minimum
from hfauto.chemistry.reaction_classification import classify_reaction
from hfauto.core.artifacts import canonical_species_id, preferred_species_by_id
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.workflow.reaction_evidence import (
    artifacts_by_reaction,
    attempts_by_reaction,
    latest_minimum_attempts_by_species,
)


class ReactionClassifyStage(Stage):
    """Combine basin, path, saddle, IRC, and intermediate evidence."""

    name = "reaction-classify"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        production = str(context.global_config.get("mode", "")).lower() == "production"
        require_real_qm = bool(config.get("require_real_qm", production))
        species = preferred_species_by_id(
            manifest,
            artifact_types=("species_optimized",),
            accept=is_accepted_optimized_minimum,
        )
        minimum_attempts = latest_minimum_attempts_by_species(manifest)
        attempts = attempts_by_reaction(manifest)
        paths = artifacts_by_reaction(manifest, "reaction_path")
        ts = artifacts_by_reaction(manifest, "reaction_validated")
        irc = artifacts_by_reaction(manifest, "reaction_path_validated")
        path_ensembles = artifacts_by_reaction(
            manifest, "path_ensemble_assessment"
        )
        mode_assessments = artifacts_by_reaction(
            manifest, "minimum_mode_following_assessment"
        )
        minima = artifacts_by_reaction(manifest, "species_optimized")
        basin_by_candidate = {
            str(item.data.get("candidate_id")): item
            for item in manifest.latest_artifacts("basin_pair_assessment")
            if item.data.get("candidate_id")
        }
        records: list[dict[str, Any]] = []

        for reaction in manifest.latest_artifacts("reaction"):
            reaction_id = str(reaction.data.get("reaction_id") or reaction.artifact_id)
            reactant = species.get(str(reaction.data.get("reactant_species_id") or ""))
            product = species.get(str(reaction.data.get("product_species_id") or ""))
            endpoint_refinement_reasons: list[str] = []
            for role, endpoint in (
                ("reactant", reactant),
                ("product", product),
            ):
                if endpoint is None:
                    continue
                canonical_id = canonical_species_id(endpoint)
                latest_attempt = minimum_attempts.get(canonical_id)
                promoted_calc_id = endpoint.qc.get("dft_calc_id") or (
                    endpoint.data.get("dft_minima") or {}
                ).get("calc_id")
                if (
                    latest_attempt is not None
                    and latest_attempt.artifact_id != promoted_calc_id
                ):
                    suffix = (
                        "latest_minimum_refinement_rejected"
                        if latest_attempt.qc.get("minimum_accepted") is not True
                        else "latest_minimum_refinement_not_promoted"
                    )
                    endpoint_refinement_reasons.append(
                        f"{role}_{suffix}"
                    )
            candidate_id = str(reaction.data.get("candidate_id") or "")
            bound_basin = basin_by_candidate.get(candidate_id)
            if endpoint_refinement_reasons:
                basin = {
                    "status": "invalid_endpoints",
                    "accepted": False,
                    "same_basin": False,
                    "distinct_basin": False,
                    "reaction_id": reaction_id,
                    "reasons": endpoint_refinement_reasons,
                }
                basin_parent = None
            elif bound_basin is not None:
                basin = dict(bound_basin.data)
                basin_parent = bound_basin.artifact_id
            elif reaction.data.get("basin_assessment"):
                basin = dict(reaction.data["basin_assessment"])
                basin_parent = None
            else:
                basin = {
                    "status": "unresolved",
                    "accepted": False,
                    "reasons": ["minimum_registry_assessment_missing"],
                }
                basin_parent = None

            intermediate_minima = [
                item
                for item in minima.get(reaction_id, [])
                if item.data.get("state") == "reaction_intermediate_candidate"
            ]
            result = classify_reaction(
                reaction,
                basin,
                path_attempts=attempts.get(reaction_id, []),
                reaction_paths=paths.get(reaction_id, []),
                ts_artifacts=ts.get(reaction_id, []),
                irc_artifacts=irc.get(reaction_id, []),
                intermediate_minima=intermediate_minima,
                path_ensemble_assessments=path_ensembles.get(reaction_id, []),
                minimum_mode_following_assessments=mode_assessments.get(
                    reaction_id, []
                ),
                reactant=reactant,
                product=product,
                require_real_qm=require_real_qm,
                barrier_resolution_kcal_mol=float(
                    config.get("barrier_resolution_kcal_mol", 0.05)
                ),
                intermediate_basin_match_rmsd_A=float(
                    config.get("intermediate_basin_match_rmsd_A", 0.20)
                ),
            )
            supported = result["scientific_conclusion_supported"] is True
            parent_ids = [
                reaction.artifact_id,
                *([basin_parent] if basin_parent else []),
                *([reactant.artifact_id] if reactant else []),
                *([product.artifact_id] if product else []),
                *[item.artifact_id for item in paths.get(reaction_id, [])],
                *[item.artifact_id for item in attempts.get(reaction_id, [])],
                *[item.artifact_id for item in ts.get(reaction_id, [])],
                *[item.artifact_id for item in irc.get(reaction_id, [])],
                *[
                    item.artifact_id
                    for item in path_ensembles.get(reaction_id, [])
                ],
                *[item.artifact_id for item in intermediate_minima],
                *[
                    item.artifact_id
                    for item in mode_assessments.get(reaction_id, [])
                ],
            ]
            artifact = Artifact(
                artifact_id=f"reaction_classification_{path_token(reaction_id, 48)}",
                artifact_type="reaction_classification",
                parents=list(dict.fromkeys(parent_ids)),
                data=result,
                qc={
                    "scientific_conclusion_supported": supported,
                    "classification": result["classification"],
                    "require_real_qm": require_real_qm,
                },
                provenance={"created_by": self.name},
                status=ArtifactStatus(
                    status="success" if supported else "partial",
                    category=None if supported else "reaction_classification_unresolved",
                    reason="; ".join(result["reasons"]),
                ),
            )
            out.add_artifact(artifact)
            records.append(artifact.model_dump(mode="json"))

        write_jsonl(records, out_dir / "reaction_classifications.jsonl")
        out.metadata["reaction_classification_count"] = len(records)
        out.metadata["reaction_classification_supported_count"] = sum(
            record["qc"]["scientific_conclusion_supported"] is True
            for record in records
        )
        return out
