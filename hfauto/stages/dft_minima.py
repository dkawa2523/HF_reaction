from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_qm_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.qc import minimum_qc_from_freq
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _preferred_species_inputs(manifest: Manifest, states: set[str]) -> list[Artifact]:
    """Prefer preoptimized geometry artifacts while preserving canonical species IDs."""
    preopt_by_id: dict[str, Artifact] = {}
    for art in manifest.iter_artifacts("species_preopt"):
        if art.status.status == "success" and art.data.get("state") in states:
            preopt_by_id[_canonical_species_id(art)] = art
    selected: list[Artifact] = []
    seen: set[str] = set()
    for species in manifest.latest_artifacts("species"):
        if species.data.get("state") not in states:
            continue
        cid = _canonical_species_id(species)
        if cid in seen:
            continue
        seen.add(cid)
        selected.append(preopt_by_id.get(cid, species))
    return selected


def _optimized_species_artifacts(species: Artifact, calc: Artifact, update_species: bool) -> list[Artifact]:
    final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
    canonical_id = _canonical_species_id(species)
    data = dict(species.data)
    data["species_id"] = canonical_id
    data["source_species_id"] = canonical_id
    data["source_geometry_artifact_id"] = species.artifact_id
    data["xyz_path"] = str(final_xyz)
    data["dft_minima"] = {
        "calc_id": calc.artifact_id,
        "engine": (calc.method or {}).get("engine"),
        "method_id": (calc.method or {}).get("method_id"),
        "real_orca_executed": calc.qc.get("real_orca_executed", False),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "n_imag": calc.data.get("n_imag"),
    }
    qc = {
        "optimized_by_dft_minima": True,
        "dft_calc_id": calc.artifact_id,
        "method_id": (calc.method or {}).get("method_id"),
        "engine": (calc.method or {}).get("engine"),
        "fallback_dummy": calc.qc.get("fallback_dummy", False),
        "real_orca_executed": calc.qc.get("real_orca_executed", False),
        "scf_converged": calc.qc.get("scf_converged"),
        "geometry_converged": calc.qc.get("geometry_converged"),
        "n_imag": calc.data.get("n_imag"),
        "is_minimum": calc.qc.get("is_minimum"),
        "geometry_sane": calc.qc.get("geometry_sane") or (calc.qc.get("geometry_qc", {}) or {}).get("geometry_sane"),
        "geometry_qc": calc.qc.get("geometry_qc", {}),
    }
    artifacts = [
        Artifact(
            artifact_id=f"opt_{canonical_id}",
            artifact_type="species_optimized",
            parents=[species.artifact_id, calc.artifact_id],
            paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
            data=data,
            method={"stage": "dft-minima", **(calc.method or {})},
            qc=qc,
        )
    ]
    if update_species:
        artifacts.append(
            Artifact(
                artifact_id=canonical_id,
                artifact_type="species",
                parents=[species.artifact_id, calc.artifact_id],
                paths={"xyz": str(final_xyz), "source_xyz": species.data.get("xyz_path", "")},
                data=data,
                method={"stage": "dft-minima", **(calc.method or {})},
                qc=qc,
            )
        )
    return artifacts


class DFTMinimaStage(Stage):
    name = "dft-minima"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine = get_qm_engine(config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {}))
        method = {
            "method_id": config.get("method", "r2scan3c"),
            **config.get("settings", {}),
            **config.get("method_settings", {}),
        }
        states = set(config.get("states", ["candidate", "isolated_candidate", "bare_candidate", "hf_cluster", "reactant_complex", "ion_pair"]))
        update_species = bool(config.get("update_species", True))
        records: list[dict] = []
        optimized_records: list[dict] = []
        for species in _preferred_species_inputs(manifest, states):
            canonical_id = _canonical_species_id(species)
            workdir = out_dir / canonical_id
            calc = engine.optimize_frequency(species, method, str(workdir))
            calc.method = {**(calc.method or {}), "stage": self.name}
            calc.qc.update(minimum_qc_from_freq(calc.data.get("n_imag"), species.data.get("state")))
            calc.data["species_id"] = canonical_id
            calc.data["source_geometry_artifact_id"] = species.artifact_id
            out.add_artifact(calc)
            records.append(calc.model_dump())
            if calc.status.status == "success":
                final_xyz = calc.paths.get("final_xyz") or species.data.get("xyz_path")
                optimized_records.append(
                    {
                        "species_id": canonical_id,
                        "state": species.data.get("state"),
                        "mol_id": species.data.get("mol_id"),
                        "site_id": species.data.get("site_id"),
                        "hf_n": species.data.get("hf_n"),
                        "source_geometry_artifact_id": species.artifact_id,
                        "optimized_xyz_path": final_xyz,
                        "dft_calc_id": calc.artifact_id,
                        "engine": (calc.method or {}).get("engine"),
                        "method_id": (calc.method or {}).get("method_id"),
                        "electronic_energy_hartree": calc.data.get("electronic_energy_hartree"),
                        "gibbs_298K_hartree": calc.data.get("gibbs_298K_hartree"),
                        "n_imag": calc.data.get("n_imag"),
                        "geometry_sane": calc.qc.get("geometry_sane"),
                        "fallback_dummy": calc.qc.get("fallback_dummy", False),
                        "real_orca_executed": calc.qc.get("real_orca_executed", False),
                    }
                )
                for art in _optimized_species_artifacts(species, calc, update_species=update_species):
                    out.add_artifact(art)
        write_jsonl(records, out_dir / "dft_minima_calculation_records.jsonl")
        write_jsonl(optimized_records, out_dir / "optimized_species_records.jsonl")
        return out
