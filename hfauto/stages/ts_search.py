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
    """Return canonical species id -> best available geometry artifact."""
    selected: dict[str, Artifact] = {}
    # Base species first, then preopt, then DFT optimized to allow upgrades.
    for art_type in ["species", "species_preopt", "species_optimized"]:
        for art in manifest.iter_artifacts(art_type):
            if art.status.status != "success":
                continue
            cid = _canonical_species_id(art)
            selected[cid] = art
    return selected


def _engine_specs(config: dict[str, Any]) -> list[str | dict[str, Any]]:
    if config.get("engine_order"):
        return list(config["engine_order"])
    return [config.get("engine", "dummy")]


def _result_artifacts_record(result):
    if hasattr(result, "artifacts"):
        return list(result.artifacts), getattr(result, "record", None), bool(getattr(result, "success", False) or any(a.artifact_type == "reaction_validated" and a.qc.get("ts_validated_by_frequency") for a in result.artifacts))
    if isinstance(result, dict):
        artifacts = list(result.get("artifacts", []) or [])
        return artifacts, result.get("record"), bool(result.get("success") or any(a.artifact_type == "reaction_validated" and a.qc.get("ts_validated_by_frequency") for a in artifacts))
    artifacts = list(result or [])
    return artifacts, None, bool(any(a.artifact_type == "reaction_validated" and a.qc.get("ts_validated_by_frequency") for a in artifacts))


class TSSearchStage(Stage):
    name = "ts-search"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        species_by_id = _preferred_species_by_id(manifest)
        method_base = {
            "method_id": config.get("method", "ts-search"),
            **(config.get("settings", {}) or {}),
            **(config.get("method_settings", {}) or {}),
        }
        engine_settings = config.get("engine_settings", {}) or {}
        records: list[dict] = []

        for rxn in manifest.latest_artifacts("reaction"):
            rc_id = rxn.data.get("reactant_species_id")
            ip_id = rxn.data.get("product_species_id")
            rc = species_by_id.get(str(rc_id)) if rc_id else None
            ip = species_by_id.get(str(ip_id)) if ip_id else None
            if rc is None or ip is None:
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "missing_endpoint",
                    category="missing_input",
                    parents=[rxn.artifact_id],
                    recommended_fallback="run_build_hf_preopt_dft_minima_before_ts_search",
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
                continue

            success = False
            last_failure: Artifact | None = None
            for eng_idx, eng_spec in enumerate(_engine_specs(config)):
                # String engine specs share global engine_settings; dict specs can override.
                if isinstance(eng_spec, dict):
                    engine = get_ts_engine(eng_spec)
                    method = {**method_base, **(eng_spec.get("settings", {}) or {})}
                    eng_name = eng_spec.get("engine") or eng_spec.get("name")
                else:
                    engine = get_ts_engine(eng_spec, **engine_settings)
                    method = dict(method_base)
                    eng_name = eng_spec
                workdir = out_dir / str(rxn.data.get("reaction_id", rxn.artifact_id)) / f"attempt_{eng_idx:02d}_{str(eng_name).replace('-', '_')}"
                result = engine.search_ts(rxn, rc, ip, method, workdir)
                artifacts, record, result_success = _result_artifacts_record(result)
                for art in artifacts:
                    out.add_artifact(art)
                if record is not None:
                    records.append(record)
                failures = [a for a in artifacts if a.status.status == "failed"]
                if result_success:
                    success = True
                    break
                if failures:
                    last_failure = failures[-1]
            if not success and last_failure is None:
                fail = Artifact.failure(
                    f"ts_failed_{rxn.artifact_id}",
                    "ts_result",
                    "all_ts_backends_failed_without_diagnostics",
                    category="ts_search_failed",
                    parents=[rxn.artifact_id, rc.artifact_id, ip.artifact_id],
                    recommended_fallback="inspect_stage_logs_or_reduce_candidate_set",
                )
                out.add_artifact(fail)
                records.append(fail.model_dump())
        write_jsonl(records, out_dir / "ts_search_records.jsonl")
        return out
