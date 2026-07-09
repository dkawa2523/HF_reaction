from __future__ import annotations

from typing import Any

import pandas as pd

from hfauto.backends.kinetics.cantera import CanteraEngine
from hfauto.backends.thermo.arkane import ArkaneEngine
from hfauto.backends.registry import get_kinetics_engine
from hfauto.core.io import ensure_dir, read_json, write_json, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


class KineticsStage(Stage):
    """Compute TST kinetics and export production-connector artifacts.

    Phase 10 still calculates core rates with the internal auditable TST engine,
    but it now writes machine-readable connector QC for Cantera and Arkane.  This
    lets production runs distinguish between:
    * internal screening kinetics;
    * validated Cantera mechanism smoke tests;
    * reviewed Arkane execution; and
    * skeleton-only exports.
    """

    name = "kinetics"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)

        engine_name = config.get("engine", config.get("backend", "tst"))
        if engine_name in {"cantera", "arkane"}:
            # Cantera/Arkane remain exporters/validators; TST is the rate source.
            engine_name = "tst"
        settings = {**(config.get("settings", {}) or {}), **{k: v for k, v in config.items() if k not in {"name", "enabled", "engine", "backend", "settings", "cantera", "arkane"}}}
        kinetics_engine = get_kinetics_engine(engine_name, **settings)

        records: list[dict[str, Any]] = []
        for thermo in manifest.iter_artifacts("thermo"):
            rec = kinetics_engine.rate_record(thermo.data)
            if rec is None:
                continue
            rec["quality_tier"] = thermo.data.get("quality_tier")
            rec["confidence_score"] = thermo.data.get("confidence_score")
            rec["thermo_artifact_id"] = thermo.artifact_id
            # Carry Phase 10 thermochemistry readiness flags into kinetics.
            for key in ["production_thermo_ready", "external_goodvibes_executed", "external_goodvibes_status"]:
                if key in thermo.data:
                    rec[key] = thermo.data.get(key)
            records.append(rec)
            out.add_artifact(
                Artifact(
                    artifact_id=f"kin_{thermo.data.get('reaction_id')}_{int(float(thermo.data.get('T_K')))}K",
                    artifact_type="kinetics",
                    parents=[thermo.artifact_id],
                    data=rec,
                    qc={
                        "kinetics_status": "success",
                        "kinetics_quality": rec.get("kinetics_quality"),
                        "tunneling_model": rec.get("tunneling_model"),
                        "transmission_coefficient": rec.get("transmission_coefficient"),
                        "production_thermo_ready": rec.get("production_thermo_ready"),
                    },
                )
            )

        write_jsonl(records, out_dir / "kinetics_records.jsonl")
        pd.DataFrame(records).to_csv(out_dir / "kinetics_records.csv", index=False)

        cantera_cfg = {**(config.get("cantera", {}) or {})}
        cantera_cfg.setdefault("residence_times_s", settings.get("residence_times_s", [0.001, 0.01, 0.1, 1.0]))
        cantera = CanteraEngine(**cantera_cfg)
        mechanism_paths = cantera.export(records, out_dir)
        validation = {}
        try:
            validation = read_json(mechanism_paths.get("cantera_validation_json"))
        except Exception:
            validation = {"cantera_validation_status": "missing_validation_json"}
        out.add_artifact(
            Artifact(
                artifact_id="hfauto_kinetics_export",
                artifact_type="mechanism",
                paths=mechanism_paths,
                data={
                    "format": "hfauto_yaml_and_cantera_connector",
                    "status": "success",
                    "n_reactions": int(len(records)),
                    "production_cantera_ready": bool(validation.get("production_cantera_ready", False)),
                    "cantera_validation_status": validation.get("cantera_validation_status"),
                    "pseudo_species": True,
                    "reason_not_production_ready": "species NASA thermo and true gas-phase mechanism equations require Arkane/curated thermo unless validation says production-ready",
                },
                qc={"cantera_validation_status": validation.get("cantera_validation_status"), "production_cantera_ready": validation.get("production_cantera_ready", False)},
            )
        )
        out.add_artifact(
            Artifact(
                artifact_id="cantera_validation",
                artifact_type="connector_validation",
                paths={"json": mechanism_paths.get("cantera_validation_json"), "reactor_csv": mechanism_paths.get("cantera_reactor_results_csv")},
                data=validation,
                qc={"connector": "cantera", "status": validation.get("cantera_validation_status"), "production_ready": validation.get("production_cantera_ready", False)},
            )
        )

        if bool(config.get("write_arkane_skeleton", config.get("arkane_export", True))):
            arkane_cfg = {**(config.get("arkane", {}) or {})}
            arkane = ArkaneEngine(**arkane_cfg)
            arkane_meta = arkane.export_and_maybe_run(records, out_dir)
            out.add_artifact(
                Artifact(
                    artifact_id="arkane_tst_connector",
                    artifact_type="arkane_input",
                    paths={"python": arkane_meta.get("skeleton_path"), "summary_json": str(out_dir / "arkane" / "arkane_export_summary.json")},
                    data=arkane_meta,
                    qc={"connector": "arkane", "arkane_status": arkane_meta.get("arkane_status"), "production_ready": arkane_meta.get("production_arkane_ready", False)},
                )
            )
        out.add_artifact(Artifact(artifact_id="kinetics_table", artifact_type="table", paths={"csv": str(out_dir / "kinetics_records.csv"), "jsonl": str(out_dir / "kinetics_records.jsonl")}, data={"n_rows": int(len(records)), "table_type": "kinetics"}))
        return out
