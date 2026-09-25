#!/usr/bin/env python3
"""Refine one path-evidenced saddle seed with pysisyphus RS-I-RFO."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.backends.ts.pysisyphus_saddle import PysisyphusSaddleEngine
from hfauto.core.io import read_manifest, write_manifest
from hfauto.core.schemas.manifest import Manifest
from hfauto.workflow.reaction_inputs import (
    basin_assessment,
    reaction_and_endpoints,
    saddle_seed_candidate,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--basin-manifest", type=Path, required=True)
    parser.add_argument("--seed-hessian-manifest", type=Path, required=True)
    parser.add_argument("--reaction-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--functional", default="pbe0")
    parser.add_argument("--basis", default="def2-svpd")
    parser.add_argument("--disp-vdw", type=int, default=3)
    parser.add_argument("--required-version", default="7.2.3")
    parser.add_argument("--pysis-executable")
    parser.add_argument("--nwchem-executable")
    parser.add_argument("--method-id", default="pysis_rsirfo_nwchem_saddle")
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument(
        "--coordinate-type",
        choices=("cart", "redund", "dlc", "tric"),
        default="tric",
    )
    parser.add_argument("--trust-radius", type=float, default=0.05)
    parser.add_argument("--trust-min", type=float, default=0.005)
    parser.add_argument("--trust-max", type=float, default=0.10)
    parser.add_argument("--hessian-init-h5", type=Path)
    parser.add_argument("--hessian-evidence-dir", type=Path)
    parser.add_argument("--hessian-recalc", type=int)
    parser.add_argument("--mode-overlap", type=float, default=0.50)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=4000)
    parser.add_argument("--timeout-s", type=float, default=172800.0)
    parser.add_argument("--temperature-k", type=float, default=298.15)
    return parser


def _frequency_calculation(manifest: Manifest):
    calculations = [
        artifact
        for artifact in manifest.latest_artifacts("calculation")
        if artifact.status.status == "success"
        and artifact.data.get("frequency_count_complete") is True
    ]
    if len(calculations) != 1:
        raise ValueError(
            "Seed Hessian manifest must contain exactly one complete calculation"
        )
    return calculations[0]


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    source_path = args.manifest.resolve()
    basin_path = args.basin_manifest.resolve()
    hessian_path = args.seed_hessian_manifest.resolve()
    output_dir = args.output_dir.resolve()
    source = read_manifest(source_path)
    basin_source = read_manifest(basin_path)
    hessian_source = read_manifest(hessian_path)
    reaction, reactant, product = reaction_and_endpoints(
        source, str(args.reaction_id)
    )
    candidate, source_attempt = saddle_seed_candidate(
        source, str(args.reaction_id)
    )
    seed_hessian = _frequency_calculation(hessian_source)
    method = {
        "path_strategy": "saddle_refinement",
        "saddle_seed_candidate": candidate,
        "saddle_seed_evidence_validated": True,
        "saddle_seed_hessian_evidence": seed_hessian,
        "saddle_seed_mode_overlap_threshold": float(args.mode_overlap),
        "basin_assessment": basin_assessment(
            basin_source, str(args.reaction_id)
        ),
        "ts_endpoint_energy_tolerance_hartree": 1.0e-5,
        "program": "nwchem",
        "functional": str(args.functional),
        "basis": str(args.basis),
        "disp_vdw": int(args.disp_vdw),
        "required_program_version": str(args.required_version),
        "method_id": str(args.method_id),
        "ts_maxiter": int(args.maxiter),
        "ts_coordinate_type": str(args.coordinate_type),
        "ts_trust_radius": float(args.trust_radius),
        "ts_trust_min": float(args.trust_min),
        "ts_trust_max": float(args.trust_max),
        "ts_opt_threshold": "gau_tight",
        "allow_subprocess": True,
        "fallback_to_dummy": False,
        "ncores": int(args.ncores),
        "memory_mb": int(args.memory_mb),
        "timeout_s": float(args.timeout_s),
        "temperature_K": float(args.temperature_k),
    }
    if args.hessian_init_h5:
        method["ts_hessian_init_h5"] = str(args.hessian_init_h5.resolve())
    if args.hessian_evidence_dir:
        method["ts_hessian_evidence_dir"] = str(args.hessian_evidence_dir.resolve())
    if args.hessian_recalc:
        method["ts_hessian_recalc"] = int(args.hessian_recalc)
    if args.pysis_executable:
        method["pysis_executable"] = str(args.pysis_executable)
    if args.nwchem_executable:
        method["executable"] = str(args.nwchem_executable)
    result = PysisyphusSaddleEngine().search_ts(
        reaction, reactant, product, method, output_dir / "work"
    )
    metadata_method = {
        key: value
        for key, value in method.items()
        if key != "saddle_seed_hessian_evidence"
    }
    manifest = Manifest.merge(
        [source, basin_source, hessian_source],
        run_id=output_dir.name,
        stage="exact-hessian-rsirfo-saddle",
        metadata={
            "source_manifest": str(source_path),
            "basin_manifest": str(basin_path),
            "seed_hessian_manifest": str(hessian_path),
            "reaction_id": str(args.reaction_id),
            "source_saddle_attempt_id": source_attempt.artifact_id,
            "seed_hessian_calculation_id": seed_hessian.artifact_id,
            "method": metadata_method,
        },
    )
    manifest.extend(result.artifacts)
    result_path = write_manifest(manifest, output_dir)
    print(
        f"success={result.success}; artifacts={len(result.artifacts)}; "
        f"manifest={result_path}"
    )
    return 0 if result.success else 1


if __name__ == "__main__":
    raise SystemExit(main())
