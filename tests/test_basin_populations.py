from math import exp

import pytest

from hfauto.chemistry.populations import (
    R_KCAL_MOL_K,
    conditional_boltzmann_populations,
)
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.basin_populations import (
    BasinPopulationsStage,
    build_basin_population_rows,
)


def _basin(basin_id: str, species_id: str, method: str = "pbe0") -> Artifact:
    return Artifact(
        artifact_id=basin_id,
        artifact_type="minimum_basin",
        data={
            "basin_id": basin_id,
            "representative_species_id": species_id,
            "member_species_ids": [species_id],
            "signature": {
                "composition": {"C": 1, "H": 1, "N": 1},
                "charge": 0,
                "multiplicity": 1,
                "method": {"method": {"functional": method}},
            },
        },
    )


def _thermo(
    species_id: str,
    gibbs_hartree: float,
    *,
    production_ready: bool = True,
) -> Artifact:
    return Artifact(
        artifact_id=f"thermo_{species_id}",
        artifact_type="species_thermo",
        data={
            "species_id": species_id,
            "T_K": 298.15,
            "G_standard_hartree": gibbs_hartree,
            "production_thermo_ready": production_ready,
        },
    )


def test_conditional_populations_are_normalized_and_use_degeneracy() -> None:
    rows = conditional_boltzmann_populations(
        [
            {"species_id": "low", "G_standard_hartree": -10.0},
            {
                "species_id": "high",
                "G_standard_hartree": -10.0 + 1.0 / HARTREE_TO_KCAL_MOL,
                "degeneracy": 2.0,
            },
        ],
        temperature_K=298.15,
    )

    expected_high_weight = 2.0 * exp(-1.0 / (R_KCAL_MOL_K * 298.15))
    assert sum(row["conditional_population"] for row in rows) == pytest.approx(1.0)
    assert rows[1]["conditional_population"] == pytest.approx(
        expected_high_weight / (1.0 + expected_high_weight)
    )
    assert all(row["population_is_conditional"] for row in rows)


def test_population_groups_require_same_pes_signature() -> None:
    manifest = Manifest.new(run_id="run", stage="thermo")
    manifest.extend(
        [
            _basin("basin_a", "spc_a", "pbe0"),
            _basin("basin_b", "spc_b", "pbe0"),
            _basin("basin_c", "spc_c", "b3lyp"),
            _thermo("spc_a", -10.0),
            _thermo("spc_b", -10.0 + 1.0 / HARTREE_TO_KCAL_MOL),
            _thermo("spc_c", -10.5),
        ]
    )

    rows = build_basin_population_rows(manifest)
    by_species = {row["species_id"]: row for row in rows}

    assert by_species["spc_a"]["observed_basin_count"] == 2
    assert by_species["spc_b"]["observed_basin_count"] == 2
    assert by_species["spc_c"]["observed_basin_count"] == 1
    assert by_species["spc_c"]["conditional_population"] == pytest.approx(1.0)
    assert not any(row["ensemble_complete"] for row in rows)


def test_stage_excludes_nonproduction_thermo_and_writes_audit_table(tmp_path) -> None:
    manifest = Manifest.new(run_id="run", stage="thermo")
    manifest.extend(
        [
            _basin("basin_a", "spc_a"),
            _basin("basin_b", "spc_b"),
            _thermo("spc_a", -10.0),
            _thermo("spc_b", -10.1, production_ready=False),
        ]
    )

    result = BasinPopulationsStage().run(
        manifest,
        {},
        StageContext(out_dir=tmp_path, run_id="run", global_config={}),
    )

    populations = result.latest_artifacts("basin_population")
    assert [item.data["species_id"] for item in populations] == ["spc_a"]
    assert populations[0].data["conditional_population"] == pytest.approx(1.0)
    table = result.find("basin_population_table")
    assert table is not None
    assert table.data["exhaustive_claim_allowed"] is False
    assert (tmp_path / "basin_populations.csv").exists()
