"""External GoodVibes sensitivity across declared low-frequency cutoffs."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from hfauto.backends.registry import get_thermo_engine
from hfauto.chemistry.stoichiometry import reaction_side_terms
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.units import hartree_to_kcal_mol
from hfauto.stages.base import Stage, StageContext

_MINIMUM_TASKS = {"opt_freq", "minimum_opt_freq"}
_TS_TASKS = {"saddle_freq", "frequency", "ts_frequency"}


def _calculation_task(calculation: Artifact) -> str:
    return str(
        (calculation.method or {}).get("task")
        or calculation.data.get("task")
        or ""
    )


def _frequency_source(
    manifest: Manifest, species_id: str, *, transition_state: bool
) -> Artifact | None:
    """Select one real frequency calculation with the required saddle order."""

    expected_imag = 1 if transition_state else 0
    tasks = _TS_TASKS if transition_state else _MINIMUM_TASKS
    candidates = [
        calculation
        for calculation in manifest.latest_artifacts("calculation")
        if str(calculation.data.get("species_id") or "") == species_id
        and calculation.status.status == "success"
        and calculation.data.get("electronic_energy_hartree") is not None
        and calculation.data.get("n_imag") == expected_imag
        and _calculation_task(calculation) in tasks
        and calculation.qc.get("real_qm_executed") is True
        and calculation.qc.get("fallback_dummy") is False
        and any(
            Path(str(calculation.paths.get(key) or "")).is_file()
            for key in ("output", "out", "log")
        )
    ]
    if not candidates:
        return None
    priority = {"saddle_freq": 3, "opt_freq": 3, "frequency": 2}
    return max(candidates, key=lambda item: priority.get(_calculation_task(item), 1))


def _validated_reactions(manifest: Manifest) -> list[Artifact]:
    return [
        reaction
        for reaction in manifest.latest_artifacts("reaction_validated")
        if reaction.status.status == "success"
        and reaction.qc.get("ts_validated_by_frequency") is True
    ]


def _side_gibbs(
    reaction: Artifact,
    side: str,
    records: dict[str, dict[str, Any]],
) -> float:
    return sum(
        float(term.coefficient)
        * float(records[term.species_id]["G_standard_hartree"])
        for term in reaction_side_terms(reaction, side)
    )


class ThermoSensitivityStage(Stage):
    """Re-run external thermochemistry without re-running quantum chemistry."""

    name = "thermo-sensitivity"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("thermo-sensitivity requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        cutoffs = sorted(
            {float(value) for value in config.get("cutoffs_cm1", [50, 100, 150])}
        )
        if not cutoffs or min(cutoffs) <= 0:
            raise ValueError("thermo-sensitivity cutoffs must be positive")
        temperature = float(
            config.get(
                "temperature_K",
                (context.global_config.get("temperature_K") or [298.15])[0],
            )
        )
        scale = float(config.get("frequency_scale_factor", 0.985))
        pressure = float(context.global_config.get("standard_pressure_bar", 1.0))
        base_settings = dict(config.get("settings") or {})
        rows: list[dict[str, Any]] = []

        for reaction in _validated_reactions(manifest):
            reaction_id = str(
                reaction.data.get("reaction_id") or reaction.artifact_id
            )
            reactants = reaction_side_terms(reaction, "reactants")
            products = reaction_side_terms(reaction, "products")
            ts_species_id = str(reaction.data.get("ts_species_id") or "")
            minimum_ids = {
                term.species_id for term in [*reactants, *products]
            }
            sources = {
                species_id: _frequency_source(
                    manifest, species_id, transition_state=False
                )
                for species_id in minimum_ids
            }
            sources[ts_species_id] = _frequency_source(
                manifest, ts_species_id, transition_state=True
            )

            for cutoff in cutoffs:
                settings = {
                    **base_settings,
                    "allow_subprocess": True,
                    "use_external_goodvibes": True,
                    "prefer_external_goodvibes_g": True,
                    "apply_internal_quasi_rrho": False,
                    "quasi_rrho_cutoff_cm1": cutoff,
                    "frequency_scale_factor": scale,
                    "work_root": str(out_dir / f"cutoff_{cutoff:g}_cm-1"),
                    "extra_args": [
                        "-q",
                        "-f",
                        f"{cutoff:g}",
                        "-v",
                        f"{scale:g}",
                        "--temp",
                        f"{temperature:g}",
                        "--csv",
                        "Goodvibes.csv",
                    ],
                }
                engine = get_thermo_engine("goodvibes", **settings)
                species_records: dict[str, dict[str, Any]] = {}
                missing = sorted(
                    species_id
                    for species_id, source in sources.items()
                    if not species_id or source is None
                )
                if not missing:
                    for species_id, source in sources.items():
                        assert source is not None
                        species_records[species_id] = engine.species_thermo(
                            species_id=species_id,
                            source_calc=source,
                            T_K=temperature,
                            p_bar=None,
                            p_standard_bar=pressure,
                        )
                production_ready = bool(
                    species_records
                    and all(
                        record.get("production_thermo_ready") is True
                        for record in species_records.values()
                    )
                )
                row: dict[str, Any] = {
                    "reaction_id": reaction_id,
                    "cutoff_cm1": cutoff,
                    "frequency_scale_factor": scale,
                    "temperature_K": temperature,
                    "complete": production_ready,
                    "missing_species_ids": ";".join(missing),
                    "external_goodvibes_required": True,
                    "source_calculation_ids": ";".join(
                        sorted(
                            source.artifact_id
                            for source in sources.values()
                            if source is not None
                        )
                    ),
                    "goodvibes_csv_paths": ";".join(
                        sorted(
                            str(record.get("goodvibes_csv_path") or "")
                            for record in species_records.values()
                            if record.get("goodvibes_csv_path")
                        )
                    ),
                }
                if production_ready:
                    reactant_g = _side_gibbs(
                        reaction, "reactants", species_records
                    )
                    product_g = _side_gibbs(
                        reaction, "products", species_records
                    )
                    ts_g = float(
                        species_records[ts_species_id]["G_standard_hartree"]
                    )
                    row.update(
                        {
                            "delta_G_activation_kcal_mol": hartree_to_kcal_mol(
                                ts_g - reactant_g
                            ),
                            "delta_G_reaction_kcal_mol": hartree_to_kcal_mol(
                                product_g - reactant_g
                            ),
                        }
                    )
                rows.append(row)

        complete_rows = [row for row in rows if row["complete"] is True]
        activation_values = [
            float(row["delta_G_activation_kcal_mol"])
            for row in complete_rows
        ]
        summary = {
            "row_count": len(rows),
            "complete_row_count": len(complete_rows),
            "cutoffs_cm1": cutoffs,
            "external_goodvibes_only": True,
            "activation_min_kcal_mol": min(activation_values)
            if activation_values
            else None,
            "activation_max_kcal_mol": max(activation_values)
            if activation_values
            else None,
            "activation_range_kcal_mol": (
                max(activation_values) - min(activation_values)
                if activation_values
                else None
            ),
        }
        jsonl_path = write_jsonl(rows, out_dir / "thermo_sensitivity.jsonl")
        csv_path = out_dir / "thermo_sensitivity.csv"
        fields = sorted({key for row in rows for key in row})
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if fields:
                writer.writeheader()
                writer.writerows(rows)
        out.add_artifact(
            Artifact(
                artifact_id="thermo_sensitivity",
                artifact_type="thermo_sensitivity",
                parents=sorted(
                    {
                        calculation.artifact_id
                        for calculation in manifest.latest_artifacts(
                            "calculation"
                        )
                        if calculation.data.get("species_id")
                    }
                ),
                paths={"csv": str(csv_path), "jsonl": str(jsonl_path)},
                data=summary,
                qc={
                    "all_rows_complete": bool(rows)
                    and len(complete_rows) == len(rows),
                    "external_goodvibes_only": True,
                },
                provenance={"created_by": self.name},
            )
        )
        return out
