from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _revised_species_artifacts(species: Artifact, calc: Artifact, update_species: bool) -> list[Artifact]:
    final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
    data = dict(species.data)
    data["species_id"] = _canonical_species_id(species)
    data["source_species_id"] = _canonical_species_id(species)
    data["source_species_artifact_id"] = species.artifact_id
    data["xyz_path"] = str(final_xyz)
    data["preopt"] = {
        "calc_id": calc.artifact_id,
        "engine": (calc.method or {}).get("engine"),
        "method_id": (calc.method or {}).get("method_id"),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
    }
    qc = {
        "preoptimized": True,
        "preopt_calc_id": calc.artifact_id,
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "geometry_sane": calc.qc.get("geometry_sane") or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
        "hf_dissociated": calc.qc.get("hf_dissociated") or (calc.qc.get("geometry_qc", {}) or {}).get("hf_dissociated"),
        "proton_transferred_unintentionally": calc.qc.get("proton_transferred_unintentionally") or (calc.qc.get("geometry_qc", {}) or {}).get("proton_transferred_unintentionally"),
        "geometry_qc": calc.qc.get("geometry_qc", {}),
    }
    out = [
        Artifact(
            artifact_id=f"preopt_{species.artifact_id}",
            artifact_type="species_preopt",
            parents=[species.artifact_id, calc.artifact_id],
            paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
            data=data,
            method={"stage": "preopt", **(calc.method or {})},
            qc=qc,
        )
    ]
    if update_species:
        # Same artifact_id, later revision. Manifest.latest_artifacts("species")
        # will return this preoptimized geometry while original lineage remains.
        out.append(
            Artifact(
                artifact_id=species.artifact_id,
                artifact_type="species",
                parents=[species.artifact_id, calc.artifact_id],
                paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
                data=data,
                method={"stage": "preopt", **(calc.method or {})},
                qc=qc,
            )
        )
    return out


class PreoptStage(Stage):
    """Low-level geometry preoptimization.

    The stage stores both a dedicated ``species_preopt`` Artifact and, when
    ``update_species`` is true, a latest-revision ``species`` Artifact with the
    same ID. This supports both explicit provenance and simple downstream use.
    """

    name = "preopt"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine = get_qm_engine(config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {}))
        method = {
            "method_id": config.get("method", "gfn2-xtb"),
            **config.get("settings", {}),
            **config.get("method_settings", {}),
        }
        update_species = bool(config.get("update_species", True))
        records: list[dict] = []
        preopt_species_records: list[dict] = []
        # Use latest species revisions if preopt is resumed/re-run.
        species_inputs = manifest.latest_artifacts("species") if hasattr(manifest, "latest_artifacts") else list(manifest.iter_artifacts("species"))
        for species in species_inputs:
            workdir = out_dir / species.artifact_id
            calc = engine.optimize_frequency(species, method, str(workdir))
            calc.method = {**(calc.method or {}), "stage": self.name}
            out.add_artifact(calc)
            records.append(calc.model_dump())
            if calc.status.status == "success":
                final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
                rec = {
                    "species_id": species.artifact_id,
                    "state": species.data.get("state"),
                    "mol_id": species.data.get("mol_id"),
                    "site_id": species.data.get("site_id"),
                    "hf_n": species.data.get("hf_n"),
                    "source_species_xyz": species.data.get("xyz_path"),
                    "preopt_xyz_path": final_xyz,
                    "preopt_calc_id": calc.artifact_id,
                    "engine": calc.method.get("engine") if calc.method else None,
                    "method_id": calc.method.get("method_id") if calc.method else None,
                    "electronic_energy_hartree": calc.data.get("electronic_energy_hartree"),
                    "geometry_sane": calc.qc.get("geometry_sane") or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
                    "hf_dissociated": calc.qc.get("hf_dissociated") or (calc.qc.get("geometry_qc", {}) or {}).get("hf_dissociated"),
                    "proton_transferred_unintentionally": calc.qc.get("proton_transferred_unintentionally") or (calc.qc.get("geometry_qc", {}) or {}).get("proton_transferred_unintentionally"),
                    "fallback_dummy": calc.qc.get("fallback_dummy", False),
                }
                preopt_species_records.append(rec)
                for art in _revised_species_artifacts(species, calc, update_species=update_species):
                    out.add_artifact(art)
                out.add_artifact(
                    Artifact(
                        artifact_id=f"preopt_geom_{species.artifact_id}",
                        artifact_type="preopt_geometry",
                        parents=[species.artifact_id, calc.artifact_id],
                        paths={"xyz": str(final_xyz)},
                        data=rec,
                        qc={
                            "preoptimized": True,
                            "geometry_sane": rec["geometry_sane"],
                            "hf_dissociated": rec["hf_dissociated"],
                            "proton_transferred_unintentionally": rec["proton_transferred_unintentionally"],
                            "fallback_dummy": rec["fallback_dummy"],
                        },
                    )
                )
        write_jsonl(records, out_dir / "preopt_calculation_records.jsonl")
        write_jsonl(preopt_species_records, out_dir / "preopt_species_records.jsonl")
        return out
