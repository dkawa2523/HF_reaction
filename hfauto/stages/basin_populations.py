"""Compute conditional populations over observed frequency-validated basins."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from hashlib import sha256
from typing import Any

from hfauto.chemistry.populations import conditional_boltzmann_populations
from hfauto.core.io import ensure_dir, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _signature_key(basin: Artifact, temperature_K: float) -> tuple[str, str]:
    signature = dict(basin.data.get("signature") or {})
    payload = {
        "composition": signature.get("composition"),
        "charge": signature.get("charge"),
        "multiplicity": signature.get("multiplicity"),
        "method": signature.get("method"),
        "T_K": float(temperature_K),
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return sha256(serialized.encode("utf-8")).hexdigest()[:16], serialized


def build_basin_population_rows(
    manifest: Manifest,
    *,
    require_production_thermo: bool = True,
    degeneracy_by_species: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    """Join minimum basins to species thermochemistry and normalize by PES."""

    degeneracies = degeneracy_by_species or {}
    thermo_by_species_temperature: dict[tuple[str, float], Artifact] = {}
    for artifact in manifest.latest_artifacts("species_thermo"):
        data = artifact.data
        species_id = str(data.get("species_id") or "")
        temperature = data.get("T_K")
        gibbs = data.get("G_standard_hartree")
        if (
            artifact.status.status != "success"
            or not species_id
            or temperature is None
            or gibbs is None
            or (
                require_production_thermo
                and data.get("production_thermo_ready") is not True
            )
        ):
            continue
        thermo_by_species_temperature[(species_id, float(temperature))] = artifact

    grouped: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    signatures: dict[str, str] = {}
    for basin in manifest.latest_artifacts("minimum_basin"):
        if basin.status.status != "success":
            continue
        basin_id = str(basin.data.get("basin_id") or basin.artifact_id)
        species_ids = list(
            dict.fromkeys(
                [
                    str(basin.data.get("representative_species_id") or ""),
                    *(str(item) for item in basin.data.get("member_species_ids") or []),
                ]
            )
        )
        candidates = [
            ((species_id, temperature), thermo)
            for (species_id, temperature), thermo in thermo_by_species_temperature.items()
            if species_id in species_ids
        ]
        by_temperature: dict[float, list[tuple[str, Artifact]]] = defaultdict(list)
        for (species_id, temperature), thermo in candidates:
            by_temperature[temperature].append((species_id, thermo))
        for temperature, items in by_temperature.items():
            species_id, thermo = min(
                items,
                key=lambda item: float(item[1].data["G_standard_hartree"]),
            )
            ensemble_hash, signature_json = _signature_key(basin, temperature)
            signatures[ensemble_hash] = signature_json
            grouped[(ensemble_hash, temperature)].append(
                {
                    "ensemble_id": f"ensemble_{ensemble_hash}",
                    "basin_id": basin_id,
                    "species_id": species_id,
                    "thermo_artifact_id": thermo.artifact_id,
                    "G_standard_hartree": float(
                        thermo.data["G_standard_hartree"]
                    ),
                    "degeneracy": float(degeneracies.get(species_id, 1.0)),
                }
            )

    rows: list[dict[str, Any]] = []
    for (ensemble_hash, temperature), entries in sorted(grouped.items()):
        population_rows = conditional_boltzmann_populations(
            entries,
            temperature_K=temperature,
        )
        basin_count = len(population_rows)
        for row in population_rows:
            row.update(
                {
                    "observed_basin_count": basin_count,
                    "ensemble_signature_json": signatures[ensemble_hash],
                    "ensemble_complete": False,
                    "coverage_scope": "observed_frequency_validated_basins_only",
                }
            )
            rows.append(row)
    return rows


class BasinPopulationsStage(Stage):
    """Publish finite-ensemble populations without an exhaustiveness claim."""

    name = "basin-populations"

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("basin-populations requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        rows = build_basin_population_rows(
            manifest,
            require_production_thermo=bool(
                config.get("require_production_thermo", True)
            ),
            degeneracy_by_species={
                str(key): float(value)
                for key, value in (config.get("degeneracy_by_species") or {}).items()
            },
        )
        jsonl_path = write_jsonl(rows, out_dir / "basin_populations.jsonl")
        csv_path = out_dir / "basin_populations.csv"
        fields = sorted({key for row in rows for key in row})
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if fields:
                writer.writeheader()
                writer.writerows(rows)
        summary = {
            "ensemble_count": len({row["ensemble_id"] for row in rows}),
            "population_row_count": len(rows),
            "all_populations_conditional": True,
            "exhaustive_claim_allowed": False,
            "coverage_scope": "observed_frequency_validated_basins_only",
        }
        summary_path = write_json(out_dir / "basin_population_summary.json", summary)

        for row in rows:
            out.add_artifact(
                Artifact(
                    artifact_id=f"population_{row['basin_id']}_{int(row['T_K'])}K",
                    artifact_type="basin_population",
                    parents=[row["basin_id"], row["thermo_artifact_id"]],
                    data=dict(row),
                    qc={
                        "population_is_conditional": True,
                        "ensemble_complete": False,
                    },
                    provenance={"created_by": self.name},
                )
            )
        out.add_artifact(
            Artifact(
                artifact_id="basin_population_table",
                artifact_type="table",
                parents=[row["basin_id"] for row in rows],
                paths={
                    "csv": str(csv_path),
                    "jsonl": str(jsonl_path),
                    "summary_json": str(summary_path),
                },
                data={"table_type": "basin_populations", **summary},
                qc={
                    "population_is_conditional": True,
                    "exhaustive_claim_allowed": False,
                },
                provenance={"created_by": self.name},
            )
        )
        return out
