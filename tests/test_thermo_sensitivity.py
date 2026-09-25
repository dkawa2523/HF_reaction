import csv
from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.thermo_sensitivity import ThermoSensitivityStage


def _calculation(
    tmp_path: Path,
    species_id: str,
    energy: float,
    *,
    transition_state: bool = False,
) -> Artifact:
    output = tmp_path / f"{species_id}.out"
    output.write_text("raw frequency output\n", encoding="utf-8")
    return Artifact(
        artifact_id=f"calc_{species_id}",
        artifact_type="calculation",
        paths={"output": str(output)},
        method={"task": "saddle_freq" if transition_state else "opt_freq"},
        data={
            "species_id": species_id,
            "electronic_energy_hartree": energy,
            "n_imag": 1 if transition_state else 0,
        },
        qc={"real_qm_executed": True, "fallback_dummy": False},
    )


def test_thermo_sensitivity_requires_external_complete_rows(
    tmp_path: Path, monkeypatch
) -> None:
    manifest = Manifest.new(run_id="thermo_panel", stage="thermo")
    manifest.extend(
        [
            _calculation(tmp_path, "reactant", -10.0),
            _calculation(tmp_path, "product", -9.9),
            _calculation(
                tmp_path, "transition_state", -9.9, transition_state=True
            ),
            Artifact(
                artifact_id="validated_reaction",
                artifact_type="reaction_validated",
                data={
                    "reaction_id": "rxn",
                    "ts_species_id": "transition_state",
                    "stoichiometry": {
                        "reactants": [
                            {"species_id": "reactant", "coefficient": 1.0}
                        ],
                        "products": [
                            {"species_id": "product", "coefficient": 1.0}
                        ],
                    },
                },
                qc={"ts_validated_by_frequency": True},
            ),
        ]
    )

    class FakeGoodVibes:
        def __init__(self, cutoff: float):
            self.cutoff = cutoff

        def species_thermo(self, *, species_id, source_calc, **_kwargs):
            coefficient = {
                "reactant": 0.0,
                "product": 1.0e-6,
                "transition_state": 2.0e-6,
            }[species_id]
            return {
                "G_standard_hartree": source_calc.data[
                    "electronic_energy_hartree"
                ]
                + coefficient * self.cutoff,
                "production_thermo_ready": True,
                "goodvibes_csv_path": f"cutoff_{self.cutoff:g}.csv",
            }

    monkeypatch.setattr(
        "hfauto.stages.thermo_sensitivity.get_thermo_engine",
        lambda _name, **settings: FakeGoodVibes(
            float(settings["quasi_rrho_cutoff_cm1"])
        ),
    )
    output = ThermoSensitivityStage().run(
        manifest,
        {"cutoffs_cm1": [50, 100]},
        StageContext(
            out_dir=tmp_path / "sensitivity",
            run_id="thermo_panel",
            global_config={"temperature_K": [298.15]},
        ),
    )

    artifact = output.latest_artifacts("thermo_sensitivity")[0]
    with Path(artifact.paths["csv"]).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert artifact.qc["all_rows_complete"] is True
    assert artifact.data["complete_row_count"] == 2
    assert artifact.data["activation_range_kcal_mol"] > 0.0
    assert {float(row["cutoff_cm1"]) for row in rows} == {50.0, 100.0}
