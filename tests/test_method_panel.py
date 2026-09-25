import pytest

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.method_panel import (
    build_electronic_method_panel,
    summarize_method_uncertainty,
)


def _sp(species_id: str, energy: float) -> Artifact:
    return Artifact(
        artifact_id=f"sp_{species_id}",
        artifact_type="calculation",
        method={"method_id": "method_a", "task": "single_point"},
        data={
            "species_id": species_id,
            "electronic_energy_hartree": energy,
        },
        qc={"real_qm_executed": True, "fallback_dummy": False},
    )


def test_method_panel_builds_fixed_geometry_barrier() -> None:
    manifest = Manifest.new(run_id="panel", stage="sp")
    manifest.add_artifact(
        Artifact(
            artifact_id="reaction_with_ts",
            artifact_type="reaction_validated",
            data={
                "reaction_id": "reaction",
                "reactant_species_id": "reactant",
                "product_species_id": "product",
                "ts_species_id": "ts",
            },
        )
    )
    manifest.extend(
        [_sp("reactant", -10.0), _sp("product", -9.99), _sp("ts", -9.95)]
    )

    rows = build_electronic_method_panel(
        manifest,
        reference_method_id="method_a",
    )

    assert len(rows) == 1
    assert rows[0]["complete"] is True
    assert rows[0]["delta_electronic_activation_hartree"] == pytest.approx(0.05)
    assert rows[0]["fixed_geometry_only"] is True
    assert rows[0]["stationary_points_revalidated_on_method"] is False
    assert rows[0]["activation_shift_from_reference_kcal_mol"] == 0.0


def test_method_panel_rejects_dummy_energy() -> None:
    manifest = Manifest.new(run_id="panel", stage="sp")
    manifest.add_artifact(
        Artifact(
            artifact_id="reaction_with_ts",
            artifact_type="reaction_validated",
            data={
                "reaction_id": "reaction",
                "reactant_species_id": "reactant",
                "product_species_id": "product",
                "ts_species_id": "ts",
            },
        )
    )
    calculations = [_sp("reactant", -10.0), _sp("product", -9.99), _sp("ts", -9.95)]
    calculations[-1].qc["fallback_dummy"] = True
    manifest.extend(calculations)

    rows = build_electronic_method_panel(manifest)

    assert len(rows) == 1
    assert rows[0]["complete"] is False
    assert rows[0]["missing_species_ids"] == "ts"


def test_method_uncertainty_reports_reference_shift_and_fixed_geometry_limit() -> None:
    rows = [
        {
            "reaction_id": "reaction",
            "method_id": "reference",
            "complete": True,
            "reference_method_id": "reference",
            "delta_electronic_activation_kcal_mol": 10.0,
            "delta_electronic_reaction_kcal_mol": -2.0,
        },
        {
            "reaction_id": "reaction",
            "method_id": "alternative",
            "complete": True,
            "reference_method_id": "reference",
            "delta_electronic_activation_kcal_mol": 12.5,
            "delta_electronic_reaction_kcal_mol": -1.0,
        },
    ]

    summary = summarize_method_uncertainty(rows)[0]

    assert summary["method_uncertainty_characterized"] is True
    assert summary["activation_range_kcal_mol"] == 2.5
    assert summary["activation_max_abs_shift_from_reference_kcal_mol"] == 2.5
    assert summary["fixed_geometry_only"] is True
    assert summary["stationary_points_revalidated_on_all_methods"] is False
