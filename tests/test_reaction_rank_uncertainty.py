from pathlib import Path

import pandas as pd
import pytest

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.reaction_rank import ReactionRankStage


def test_reaction_rank_can_use_conservative_method_bound(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="rank", stage="method-panel")
    manifest.extend(
        [
            Artifact(
                artifact_id="thermo_r1",
                artifact_type="thermo",
                data={
                    "reaction_id": "r1",
                    "stoichiometry": {"reactants": [], "products": []},
                    "T_K": 298.15,
                    "delta_G_activation_standard_kcal_mol": 10.0,
                    "production_rank_eligible": True,
                },
            ),
            Artifact(
                artifact_id="u_r1",
                artifact_type="method_uncertainty",
                data={
                    "reaction_id": "r1",
                    "activation_max_abs_shift_from_reference_kcal_mol": 2.5,
                    "method_uncertainty_characterized": True,
                },
            ),
        ]
    )

    result = ReactionRankStage().run(
        manifest,
        {
            "evidence_policy": "production",
            "uncertainty_policy": "conservative_bound",
        },
        StageContext(out_dir=tmp_path, run_id="rank", global_config={}),
    )

    ranking = result.latest_artifacts("ranking")[-1]
    row = pd.read_csv(ranking.paths["csv"]).iloc[0]
    assert row["ranking_point_value"] == 10.0
    assert row["ranking_uncertainty_kcal_mol"] == 2.5
    assert row["ranking_value"] == 12.5


def test_reaction_rank_can_apply_conditional_reactant_population(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="rank", stage="basin-populations")
    manifest.extend(
        [
            Artifact(
                artifact_id="reaction_r1",
                artifact_type="reaction_validated",
                data={
                    "reaction_id": "r1",
                    "stoichiometry": {
                        "reactants": [
                            {"species_id": "spc_reactant", "coefficient": 1.0}
                        ],
                        "products": [
                            {"species_id": "spc_product", "coefficient": 1.0}
                        ],
                    },
                },
            ),
            Artifact(
                artifact_id="thermo_r1",
                artifact_type="thermo",
                data={
                    "reaction_id": "r1",
                    "stoichiometry": {
                        "reactants": [
                            {"species_id": "spc_reactant", "coefficient": 1.0}
                        ],
                        "products": [
                            {"species_id": "spc_product", "coefficient": 1.0}
                        ],
                    },
                    "T_K": 298.15,
                    "delta_G_activation_standard_kcal_mol": 10.0,
                    "production_rank_eligible": True,
                },
            ),
            Artifact(
                artifact_id="population_basin_reactant_298K",
                artifact_type="basin_population",
                data={
                    "species_id": "spc_reactant",
                    "T_K": 298.15,
                    "conditional_population": 0.1,
                    "ensemble_complete": False,
                },
            ),
        ]
    )

    result = ReactionRankStage().run(
        manifest,
        {
            "evidence_policy": "production",
            "population_policy": "conditional_adjustment",
        },
        StageContext(out_dir=tmp_path, run_id="rank", global_config={}),
    )

    row = pd.read_csv(result.latest_artifacts("ranking")[-1].paths["csv"]).iloc[0]
    assert row["reactant_conditional_population"] == 0.1
    assert row["population_adjustment_applied"]
    assert row["population_adjustment_kcal_mol"] > 1.3
    assert row["ranking_value"] > 11.3


def test_population_adjustment_rejects_nonactivation_metric(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="rank", stage="basin-populations")
    manifest.add_artifact(
        Artifact(
            artifact_id="thermo_r1",
            artifact_type="thermo",
            data={
                "reaction_id": "r1",
                "stoichiometry": {"reactants": [], "products": []},
                "T_K": 298.15,
                "delta_G_reaction_standard_kcal_mol": -1.0,
                "production_rank_eligible": True,
            },
        )
    )

    with pytest.raises(ValueError, match="activation free-energy metric"):
        ReactionRankStage().run(
            manifest,
            {
                "metric": "delta_G_reaction_standard_kcal_mol",
                "population_policy": "conditional_adjustment",
            },
            StageContext(out_dir=tmp_path, run_id="rank", global_config={}),
        )
