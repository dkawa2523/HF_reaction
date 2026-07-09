from __future__ import annotations

from typing import Any

from hfauto.backends.registry import get_ts_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _preferred_species_by_id(manifest: Manifest) -> dict[str, Artifact]:
    selected: dict[str, Artifact] = {}
    for art_type in ["species", "species_preopt", "species_optimized"]:
        for art in manifest.iter_artifacts(art_type):
            if art.status.status != "success":
                continue
            selected[_canonical_species_id(art)] = art
    return selected


def _result_artifacts_record(result):
    if hasattr(result, "artifacts"):
        return list(result.artifacts), getattr(result, "record", None)
    if isinstance(result, dict):
        return list(result.get("artifacts", []) or []), result.get("record")
    return list(result or []), None


class IRCStage(Stage):
    name = "irc"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        species_by_id = _preferred_species_by_id(manifest)
        engine = get_ts_engine(config.get("engine", "dummy"), **(config.get("engine_settings", {}) or {}))
        method = {
            "method_id": config.get("method", "irc"),
            **(config.get("settings", {}) or {}),
            **(config.get("method_settings", {}) or {}),
        }
        records: list[dict] = []
        reactions = list(manifest.latest_artifacts("reaction_validated"))
        for rxn in reactions:
            # Preserve recoverable data when TS frequency is not validated, but do not
            # attempt IRC unless explicitly requested.
            if not rxn.qc.get("ts_validated_by_frequency", False) and not config.get("run_unvalidated_ts", False):
                fail = Artifact.failure(
                    f"irc_skipped_{rxn.artifact_id}",
                    "irc",
                    "TS frequency/mode QC failed; IRC skipped",
                    category="ts_not_validated",
                    parents=[rxn.artifact_id],
                    recoverable=True,
                    recommended_fallback="inspect_ts_or_retry_ts_search_with_scan_optts",
                    reaction_id=rxn.data.get("reaction_id"),
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue
            ts_id = rxn.data.get("ts_species_id")
            rc_id = rxn.data.get("reactant_species_id")
            ip_id = rxn.data.get("product_species_id")
            ts_species = species_by_id.get(str(ts_id)) if ts_id else None
            rc = species_by_id.get(str(rc_id)) if rc_id else None
            ip = species_by_id.get(str(ip_id)) if ip_id else None
            if ts_species is None or rc is None or ip is None:
                fail = Artifact.failure(
                    f"irc_failed_{rxn.artifact_id}",
                    "irc",
                    "missing_ts_or_endpoint_species",
                    category="missing_input",
                    parents=[rxn.artifact_id],
                    recommended_fallback="rerun_ts_search_or_build_hf",
                    reaction_id=rxn.data.get("reaction_id"),
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue
            workdir = out_dir / str(rxn.data.get("reaction_id", rxn.artifact_id))
            result = engine.run_irc(rxn, ts_species, rc, ip, method, workdir)
            irc_artifacts, record = _result_artifacts_record(result)
            for art in irc_artifacts:
                out.add_artifact(art)
            if record is not None:
                records.append(record)
            primary_irc = next((a for a in irc_artifacts if a.artifact_type == "irc"), None)
            if primary_irc is not None and primary_irc.status.status == "success":
                path_validated = Artifact(
                    artifact_id=f"path_validated_{rxn.artifact_id}",
                    artifact_type="reaction_path_validated",
                    parents=[rxn.artifact_id, primary_irc.artifact_id],
                    data={
                        **rxn.data,
                        "irc_artifact_id": primary_irc.artifact_id,
                        "irc_validated": primary_irc.qc.get("irc_validated"),
                        "real_irc_executed": primary_irc.qc.get("real_irc_executed", False),
                    },
                    qc={
                        **rxn.qc,
                        "irc_validated": primary_irc.qc.get("irc_validated"),
                        "real_irc_executed": primary_irc.qc.get("real_irc_executed", False),
                        "fallback_dummy": primary_irc.qc.get("fallback_dummy", False),
                    },
                )
                out.add_artifact(path_validated)
                revised = Artifact(
                    artifact_id=rxn.artifact_id,
                    artifact_type="reaction_validated",
                    parents=[rxn.artifact_id, primary_irc.artifact_id],
                    data={**rxn.data, "irc_artifact_id": primary_irc.artifact_id},
                    qc={**rxn.qc, "irc_validated": primary_irc.qc.get("irc_validated"), "real_irc_executed": primary_irc.qc.get("real_irc_executed", False), "irc_artifact_id": primary_irc.artifact_id},
                )
                out.add_artifact(revised)
        write_jsonl(records, out_dir / "irc_records.jsonl")
        return out
