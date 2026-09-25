from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.core.io import write_json
from hfauto.core.schemas.artifact import Artifact
from hfauto_ops.core.array_jobs import write_array_plan
from hfauto_ops.core.autotune import write_autotune_report
from hfauto_ops.core.backend_compare import write_backend_comparison
from hfauto_ops.core.compare import compare_runs
from hfauto_ops.core.qcarchive import write_qcarchive_bundle
from hfauto_ops.core.resource import build_resource_plan, write_resource_plan
from hfauto_ops.core.retry import write_retry_plan
from hfauto_ops.core.reuse import write_reuse_bundle
from hfauto_ops.core.run_index import write_index_bundle
from hfauto_ops.qcarchive.connector import write_qcarchive_plan
from hfauto_ops.qcarchive.export import write_qcarchive_export
from hfauto_ops.reports.operations_report import write_operations_report
from hfauto_ops.schedulers.controller import scheduler_status
from hfauto_ops.schedulers.production import write_multi_scheduler_bundle


def build_operations_bundle(
    run_dir: str | Path,
    out_dir: str | Path,
    pipeline_config: str | Path | None = None,
    ops_config: dict[str, Any] | None = None,
) -> dict[str, Path]:
    ops_config = ops_config or {}
    run = Path(run_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    paths.update(write_index_bundle(run, out))

    resource_csv: Path | None = None
    array_jsonl: Path | None = None
    if pipeline_config:
        plan = build_resource_plan(pipeline_config, run_dir=run, ops_config=ops_config)
        resource_paths = write_resource_plan(plan, out)
        paths.update({f"resource_{k}": v for k, v in resource_paths.items()})
        resource_csv = resource_paths.get("csv")

        array_paths = write_array_plan(
            run,
            out,
            pipeline_config=pipeline_config,
            ops_config=ops_config,
            stages=ops_config.get("artifact_level_stages") or ops_config.get("array_stages"),
        )
        paths.update({f"array_{k}": v for k, v in array_paths.items()})
        array_jsonl = array_paths.get("array_jsonl")

        schedulers = ops_config.get("schedulers") or [ops_config.get("scheduler", "slurm")]
        scheduler_paths = write_multi_scheduler_bundle(
            resource_csv,
            array_jsonl,
            out / "scheduler",
            schedulers=[str(s) for s in schedulers],
            project_root=ops_config.get("project_root", "$PWD"),
        )
        paths.update({f"scheduler_{k}": v for k, v in scheduler_paths.items()})
        paths.update({f"scheduler_status_{k}": v for k, v in scheduler_status(out, scheduler=str(schedulers[0]), dry_run=True).items()})

    paths.update({f"retry_{k}": v for k, v in write_retry_plan(run, out, pipeline_config).items()})

    registry_db = ops_config.get("registry_db") or str(out / "calculation_registry.sqlite")
    paths.update(write_reuse_bundle(
        run,
        out,
        job_plan_csv=paths.get("array_array_csv"),
        registry_db=registry_db,
        cache_runs=ops_config.get("cache_roots") or ops_config.get("cache_runs") or ops_config.get("reuse_from") or [],
    ))

    qca_paths = write_qcarchive_export(run, out / "qcarchive", dataset_name=ops_config.get("qcarchive_dataset_name") or ops_config.get("qcarchive_dataset") or f"hfauto_{run.name}")
    paths.update(qca_paths)
    # Root-level import plan used by operations reports.
    paths.update(write_qcarchive_plan(run, out, dataset_name=ops_config.get("qcarchive_dataset_name") or ops_config.get("qcarchive_dataset") or f"hfauto_{run.name}"))
    if paths.get("artifact_job_csv") or paths.get("array_array_csv"):
        # Legacy/review payload keyed from artifact/array plan for external QCArchive mapping.
        qca_input = paths.get("artifact_job_csv") or paths.get("array_array_csv")
        qca_legacy = write_qcarchive_bundle(qca_input, out, allow_submit=bool(ops_config.get("qcarchive_allow_submit", False)))
        paths.update({f"qcarchive_{k}": v for k, v in qca_legacy.items()})

    paths.update(write_backend_comparison(run, out))

    autotune_dirs = [run, out] + [Path(p) for p in ops_config.get("autotune_history_dirs", [])]
    paths.update(write_autotune_report(autotune_dirs, out, safety_factor=float(ops_config.get("autotune_safety_factor", 1.25))))

    compare_to = ops_config.get("compare_to")
    if compare_to:
        cmp_dir = out / "run_comparison"
        cmp_paths = compare_runs(compare_to, run, cmp_dir)
        paths.update({f"compare_{k}": v for k, v in cmp_paths.items()})

    report = write_operations_report(out, title="hfauto production operations report")
    paths["operations_report"] = report
    summary = {k: str(v) for k, v in paths.items()}
    paths["operations_manifest"] = write_json(out / "operations_manifest.json", summary)
    return paths


def operations_artifact(run_dir: str | Path, out_dir: str | Path, paths: dict[str, Path]) -> Artifact:
    return Artifact(
        artifact_id="operations_bundle",
        artifact_type="operations_bundle",
        paths={k: str(v) for k, v in paths.items()},
        data={
            "operations_package": "hfauto_ops",
            "run_dir": str(run_dir),
            "n_outputs": len(paths),
            "capabilities": [
                "artifact_index",
                "stage_index",
                "resource_plan",
                "artifact_level_array_plan",
                "slurm_pbs_lsf_local_array_scripts",
                "scheduler_submit_status_cancel_dry_run",
                "retry_plan",
                "deduplication_reuse_plan",
                "qcarchive_payload_export",
                "backend_comparison_dashboard",
                "resource_autotune_proposal",
                "operations_report",
            ],
        },
        provenance={"created_at": Artifact.now_iso()},
        qc={"read_only_consumer": True, "scientific_results_modified": False, "production_hpc_integration_plan": True},
    )
