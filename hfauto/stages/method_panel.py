"""Compare fixed-geometry electronic energies across configured QM methods."""

from __future__ import annotations

import csv
from pathlib import Path
from statistics import pstdev
from typing import Any

from hfauto.chemistry.stoichiometry import reaction_side_terms
from hfauto.core.ids import path_token
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.units import hartree_to_kcal_mol
from hfauto.stages.base import Stage, StageContext


def _stoichiometric_species(
    reaction: Artifact, side: str
) -> list[tuple[str, float]]:
    return [
        (term.species_id, float(term.coefficient))
        for term in reaction_side_terms(reaction, side)
    ]


def _validated_reactions(manifest: Manifest) -> list[Artifact]:
    by_id: dict[str, Artifact] = {}
    for artifact_type in ("reaction_validated", "reaction_path_validated"):
        for artifact in manifest.latest_artifacts(artifact_type):
            if artifact.status.status != "success":
                continue
            reaction_id = str(
                artifact.data.get("reaction_id") or artifact.artifact_id
            )
            previous = by_id.get(reaction_id)
            if previous is None:
                by_id[reaction_id] = artifact
                continue
            merged = artifact.model_copy(deep=True)
            merged.data = {**previous.data, **artifact.data}
            by_id[reaction_id] = merged
    return [by_id[key] for key in sorted(by_id)]


def build_electronic_method_panel(
    manifest: Manifest,
    *,
    method_ids: set[str] | None = None,
    require_real_qm: bool = True,
    reference_method_id: str | None = None,
) -> list[dict[str, Any]]:
    """Build comparable reaction/activation energies from single points only."""

    energies: dict[tuple[str, str], tuple[float, str]] = {}
    method_metadata: dict[str, dict[str, Any]] = {}
    for calculation in manifest.latest_artifacts("calculation"):
        method_id = str((calculation.method or {}).get("method_id") or "")
        species_id = str(calculation.data.get("species_id") or "")
        task = str(
            (calculation.method or {}).get("task")
            or calculation.data.get("task")
            or ""
        )
        energy = calculation.data.get("electronic_energy_hartree")
        if (
            calculation.status.status != "success"
            or task not in {"single_point", "sp"}
            or not method_id
            or not species_id
            or energy is None
            or (method_ids is not None and method_id not in method_ids)
            or (
                require_real_qm
                and (
                    calculation.qc.get("real_qm_executed") is not True
                    or calculation.qc.get("fallback_dummy") is not False
                )
            )
        ):
            continue
        energies[(method_id, species_id)] = (
            float(energy),
            calculation.artifact_id,
        )
        method_metadata.setdefault(method_id, dict(calculation.method or {}))

    methods = sorted({method_id for method_id, _species_id in energies})
    rows: list[dict[str, Any]] = []
    for reaction in _validated_reactions(manifest):
        reaction_id = str(
            reaction.data.get("reaction_id") or reaction.artifact_id
        )
        reactants = _stoichiometric_species(reaction, "reactants")
        products = _stoichiometric_species(reaction, "products")
        ts_species_id = str(reaction.data.get("ts_species_id") or "")
        for method_id in methods:
            required = [
                *(species_id for species_id, _coefficient in reactants),
                *(species_id for species_id, _coefficient in products),
                ts_species_id,
            ]
            missing = sorted(
                species_id
                for species_id in required
                if species_id and (method_id, species_id) not in energies
            )
            complete = bool(reactants and products and ts_species_id and not missing)
            row: dict[str, Any] = {
                "reaction_id": reaction_id,
                "method_id": method_id,
                "functional": method_metadata[method_id].get("functional"),
                "basis": method_metadata[method_id].get("basis"),
                "disp_vdw": method_metadata[method_id].get("disp_vdw"),
                "grid": method_metadata[method_id].get("grid"),
                "scf_energy_tolerance": method_metadata[method_id].get(
                    "scf_energy_tolerance"
                ),
                "complete": complete,
                "missing_species_ids": ";".join(missing),
                "fixed_geometry_only": True,
                "stationary_points_revalidated_on_method": False,
            }
            if complete:
                reactant_energy = sum(
                    coefficient * energies[(method_id, species_id)][0]
                    for species_id, coefficient in reactants
                )
                product_energy = sum(
                    coefficient * energies[(method_id, species_id)][0]
                    for species_id, coefficient in products
                )
                ts_energy = energies[(method_id, ts_species_id)][0]
                delta_reaction = product_energy - reactant_energy
                delta_activation = ts_energy - reactant_energy
                row.update(
                    {
                        "reactant_electronic_energy_hartree": reactant_energy,
                        "product_electronic_energy_hartree": product_energy,
                        "ts_electronic_energy_hartree": ts_energy,
                        "delta_electronic_reaction_hartree": delta_reaction,
                        "delta_electronic_activation_hartree": delta_activation,
                        "delta_electronic_reaction_kcal_mol": hartree_to_kcal_mol(
                            delta_reaction
                        ),
                        "delta_electronic_activation_kcal_mol": hartree_to_kcal_mol(
                            delta_activation
                        ),
                        "calculation_ids": ";".join(
                            energies[(method_id, species_id)][1]
                            for species_id in required
                        ),
                    }
                )
            rows.append(row)
    if reference_method_id:
        references = {
            row["reaction_id"]: row
            for row in rows
            if row["method_id"] == reference_method_id
            and row["complete"] is True
        }
        for row in rows:
            reference = references.get(row["reaction_id"])
            row["reference_method_id"] = reference_method_id
            if row["complete"] is True and reference is not None:
                row["activation_shift_from_reference_kcal_mol"] = float(
                    row["delta_electronic_activation_kcal_mol"]
                ) - float(reference["delta_electronic_activation_kcal_mol"])
                row["reaction_shift_from_reference_kcal_mol"] = float(
                    row["delta_electronic_reaction_kcal_mol"]
                ) - float(reference["delta_electronic_reaction_kcal_mol"])
    return rows


def summarize_method_uncertainty(
    rows: list[dict[str, Any]],
    *,
    minimum_complete_methods: int = 2,
) -> list[dict[str, Any]]:
    """Summarize method spread while retaining its fixed-geometry limitation."""

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("complete") is True:
            grouped.setdefault(str(row["reaction_id"]), []).append(row)
    summaries: list[dict[str, Any]] = []
    for reaction_id in sorted(grouped):
        complete = grouped[reaction_id]
        activation = [
            float(row["delta_electronic_activation_kcal_mol"])
            for row in complete
        ]
        reaction = [
            float(row["delta_electronic_reaction_kcal_mol"])
            for row in complete
        ]
        reference_rows = [
            row
            for row in complete
            if row.get("reference_method_id") == row.get("method_id")
        ]
        reference = reference_rows[0] if reference_rows else None
        reference_activation = (
            float(reference["delta_electronic_activation_kcal_mol"])
            if reference
            else None
        )
        reference_reaction = (
            float(reference["delta_electronic_reaction_kcal_mol"])
            if reference
            else None
        )
        characterized = len(complete) >= int(minimum_complete_methods)
        summaries.append(
            {
                "reaction_id": reaction_id,
                "complete_method_count": len(complete),
                "method_ids": ";".join(sorted(str(row["method_id"]) for row in complete)),
                "minimum_complete_methods": int(minimum_complete_methods),
                "method_uncertainty_characterized": characterized,
                "activation_min_kcal_mol": min(activation),
                "activation_max_kcal_mol": max(activation),
                "activation_range_kcal_mol": max(activation) - min(activation),
                "activation_half_range_kcal_mol": (max(activation) - min(activation)) / 2.0,
                "activation_population_std_kcal_mol": pstdev(activation),
                "reaction_min_kcal_mol": min(reaction),
                "reaction_max_kcal_mol": max(reaction),
                "reaction_range_kcal_mol": max(reaction) - min(reaction),
                "reaction_population_std_kcal_mol": pstdev(reaction),
                "reference_method_id": (
                    str(reference["method_id"]) if reference else None
                ),
                "reference_activation_kcal_mol": reference_activation,
                "reference_reaction_kcal_mol": reference_reaction,
                "activation_max_abs_shift_from_reference_kcal_mol": (
                    max(abs(value - reference_activation) for value in activation)
                    if reference_activation is not None
                    else None
                ),
                "reaction_max_abs_shift_from_reference_kcal_mol": (
                    max(abs(value - reference_reaction) for value in reaction)
                    if reference_reaction is not None
                    else None
                ),
                "fixed_geometry_only": True,
                "stationary_points_revalidated_on_all_methods": False,
                "interpretation": (
                    "Electronic-energy sensitivity at common geometries; this does not "
                    "prove stationary points or IRC connectivity on alternative PESs."
                ),
            }
        )
    return summaries


class MethodPanelStage(Stage):
    """Publish a fixed-geometry method panel without overstating PES evidence."""

    name = "method-panel"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("method-panel requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        configured = config.get("method_ids")
        rows = build_electronic_method_panel(
            manifest,
            method_ids={str(value) for value in configured} if configured else None,
            require_real_qm=bool(config.get("require_real_qm", True)),
            reference_method_id=(
                str(config["reference_method_id"])
                if config.get("reference_method_id")
                else None
            ),
        )
        uncertainty_rows = summarize_method_uncertainty(
            rows,
            minimum_complete_methods=int(
                config.get("minimum_complete_methods", 2)
            ),
        )
        jsonl_path = write_jsonl(rows, out_dir / "electronic_method_panel.jsonl")
        csv_path = Path(out_dir) / "electronic_method_panel.csv"
        fields = sorted({key for row in rows for key in row})
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        uncertainty_jsonl = write_jsonl(
            uncertainty_rows, out_dir / "method_uncertainty.jsonl"
        )
        uncertainty_csv = Path(out_dir) / "method_uncertainty.csv"
        uncertainty_fields = sorted(
            {key for row in uncertainty_rows for key in row}
        )
        with uncertainty_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=uncertainty_fields)
            if uncertainty_fields:
                writer.writeheader()
                writer.writerows(uncertainty_rows)
        complete_count = sum(row["complete"] is True for row in rows)
        out.add_artifact(
            Artifact(
                artifact_id="electronic_method_panel",
                artifact_type="table",
                parents=[
                    calculation.artifact_id
                    for calculation in manifest.latest_artifacts("calculation")
                    if str((calculation.method or {}).get("task")) == "single_point"
                ],
                paths={"csv": str(csv_path), "jsonl": str(jsonl_path)},
                data={
                    "row_count": len(rows),
                    "complete_row_count": complete_count,
                    "fixed_geometry_only": True,
                    "stationary_points_revalidated_on_method": False,
                },
                qc={"all_rows_complete": bool(rows) and complete_count == len(rows)},
                provenance={"created_by": self.name},
            )
        )
        for row in uncertainty_rows:
            out.add_artifact(
                Artifact(
                    artifact_id=(
                        "method_uncertainty_"
                        + path_token(str(row["reaction_id"]), max_length=48)
                    ),
                    artifact_type="method_uncertainty",
                    parents=["electronic_method_panel"],
                    paths={
                        "csv": str(uncertainty_csv),
                        "jsonl": str(uncertainty_jsonl),
                    },
                    data=row,
                    qc={
                        "fixed_geometry_only": True,
                        "method_uncertainty_characterized": row[
                            "method_uncertainty_characterized"
                        ],
                        "stationary_points_revalidated_on_all_methods": False,
                    },
                    provenance={"created_by": self.name},
                )
            )
        return out
