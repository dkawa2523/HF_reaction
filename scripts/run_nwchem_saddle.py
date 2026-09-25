#!/usr/bin/env python3
"""Bracket or refine one molecular saddle candidate with NWChem."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.backends.ts.nwchem_saddle import NWChemSaddleEngine
from hfauto.chemistry.xyz import write_xyz
from hfauto.chemistry.xyz_trajectory import (
    endpoint_mode_reference,
    read_xyz_trajectory,
)
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
    parser.add_argument("--reaction-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed-hessian-manifest", type=Path)
    parser.add_argument("--seed-xyz", type=Path)
    parser.add_argument(
        "--strategy",
        choices=("saddle_refinement", "bracketed_saddle_search"),
        default="saddle_refinement",
    )
    parser.add_argument("--bracket-only", action="store_true")
    parser.add_argument(
        "--bracket-fractions",
        type=float,
        nargs="+",
        default=[0.0, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5, 1.0],
    )
    parser.add_argument("--bracket-minimum-prominence", type=float, default=1.0e-5)
    parser.add_argument("--segment-path", type=Path)
    parser.add_argument("--segment-start-index", type=int, default=0)
    parser.add_argument(
        "--hessian-restart-count", type=int, choices=(0, 1), default=0
    )
    parser.add_argument("--functional", default="pbe0")
    parser.add_argument("--basis", default="def2-svpd")
    parser.add_argument("--disp-vdw", type=int)
    parser.add_argument("--required-version")
    parser.add_argument("--method-id", default="nwchem_mode_following_saddle")
    parser.add_argument("--geometry-maxiter", type=int, default=100)
    parser.add_argument("--driver-trust", type=float, default=0.10)
    parser.add_argument("--saddle-step", type=float, default=0.03)
    parser.add_argument("--mode-overlap", type=float, default=0.50)
    parser.add_argument("--mode-margin", type=float, default=0.10)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=4000)
    parser.add_argument("--timeout-s", type=float, default=21600.0)
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
    output_dir = args.output_dir.resolve()
    source = read_manifest(source_path)
    basin_source = read_manifest(basin_path)
    reaction, reactant, product = reaction_and_endpoints(
        source, str(args.reaction_id)
    )
    segment = None
    if args.segment_path is not None:
        segment_path = args.segment_path.resolve()
        images = read_xyz_trajectory(segment_path)
        start_index = int(args.segment_start_index)
        end_index = start_index + 1
        if not 0 <= start_index < len(images) - 1:
            raise ValueError("segment-start-index must select two adjacent path images")
        path_artifacts = [
            artifact
            for artifact in source.latest_artifacts("reaction_path")
            if artifact.status.status == "success"
            and str(artifact.data.get("reaction_id") or "") == str(args.reaction_id)
            and artifact.data.get("reaction_path", {}).get("converged") is True
        ]
        if len(path_artifacts) != 1:
            raise ValueError("source manifest must contain one converged reaction path")
        path_artifact = path_artifacts[0]
        path_images = path_artifact.data["reaction_path"]["images"]
        if len(path_images) != len(images):
            raise ValueError("segment path image count differs from reaction-path evidence")
        segment_dir = output_dir / "segment_inputs"
        start_path = write_xyz(images[start_index], segment_dir / "start.xyz")
        end_path = write_xyz(images[end_index], segment_dir / "end.xyz")
        segment = {
            "accepted": True,
            "source_reaction_path_artifact_id": path_artifact.artifact_id,
            "source_path_xyz": str(segment_path),
            "start_index": start_index,
            "end_index": end_index,
            "start_label": f"path_image_{start_index}",
            "target_label": f"path_image_{end_index}",
            "start_xyz": str(start_path.resolve()),
            "target_xyz": str(end_path.resolve()),
            "start_energy_hartree": float(path_images[start_index]["energy_hartree"]),
            "target_energy_hartree": float(path_images[end_index]["energy_hartree"]),
        }
    candidate: dict = {}
    source_attempt = None
    if args.strategy == "saddle_refinement" and args.seed_xyz is None:
        candidate, source_attempt = saddle_seed_candidate(
            source, str(args.reaction_id)
        )
    if args.seed_xyz is not None:
        seed_xyz = args.seed_xyz.resolve()
        if not seed_xyz.is_file():
            raise FileNotFoundError(seed_xyz)
        candidate = {
            **candidate,
            "xyz_path": str(seed_xyz),
            "source": "explicit_path_evidenced_seed",
        }
        if args.segment_path is not None:
            candidate["reaction_mode_reference"] = endpoint_mode_reference(
                images[start_index],
                images[end_index],
                source=(
                    f"{segment_path}:segment={start_index}-{end_index}"
                ),
            )
    method = {
        "path_strategy": str(args.strategy),
        "endpoint_basin_status": "distinct_basin",
        "basin_assessment": basin_assessment(
            basin_source, str(args.reaction_id)
        ),
        "saddle_seed_candidate": candidate,
        "saddle_seed_evidence_validated": True,
        "saddle_bracket_only": bool(args.bracket_only),
        "saddle_bracket_fractions": [float(value) for value in args.bracket_fractions],
        "saddle_bracket_minimum_prominence_hartree": float(
            args.bracket_minimum_prominence
        ),
        "saddle_seed_hessian_precheck": True,
        "saddle_initial_hessian": True,
        "saddle_seed_mode_overlap_threshold": float(args.mode_overlap),
        "ts_mode_overlap_threshold": float(args.mode_overlap),
        "saddle_seed_mode_overlap_margin": float(args.mode_margin),
        "driver_trust": float(args.driver_trust),
        "driver_saddle_step": float(args.saddle_step),
        "geometry_maxiter": int(args.geometry_maxiter),
        "optimization_convergence": "tight",
        "saddle_max_hessian_restarts": 1,
        "saddle_hessian_restart_count": int(args.hessian_restart_count),
        "ts_endpoint_energy_tolerance_hartree": 1.0e-5,
        "allow_subprocess": True,
        "fallback_to_dummy": False,
        "method_id": str(args.method_id),
        "functional": str(args.functional),
        "basis": str(args.basis),
        "ncores": int(args.ncores),
        "memory_mb": int(args.memory_mb),
        "timeout_s": float(args.timeout_s),
        "temperature_K": float(args.temperature_k),
    }
    if segment is not None:
        method["saddle_bracket_segment"] = segment
    evidence_manifests = [source, basin_source]
    seed_hessian_id = None
    if args.seed_hessian_manifest is not None:
        hessian_source = read_manifest(args.seed_hessian_manifest.resolve())
        seed_hessian = _frequency_calculation(hessian_source)
        method["saddle_seed_hessian_evidence"] = seed_hessian
        seed_hessian_id = seed_hessian.artifact_id
        evidence_manifests.append(hessian_source)
    if args.disp_vdw is not None:
        method["disp_vdw"] = int(args.disp_vdw)
    if args.required_version:
        method["required_program_version"] = str(args.required_version)

    result = NWChemSaddleEngine().search_ts(
        reaction,
        reactant,
        product,
        method,
        output_dir / "work",
    )
    metadata_method = {
        key: value
        for key, value in method.items()
        if key != "saddle_seed_hessian_evidence"
    }
    manifest = Manifest.merge(
        evidence_manifests,
        run_id=output_dir.name,
        stage="mode-following-saddle",
        metadata={
            "source_manifest": str(source_path),
            "basin_manifest": str(basin_path),
            "reaction_id": str(args.reaction_id),
            "source_saddle_attempt_id": (
                source_attempt.artifact_id if source_attempt is not None else None
            ),
            "restart_seed_xyz": (
                str(args.seed_xyz.resolve()) if args.seed_xyz is not None else None
            ),
            "seed_hessian_calculation_id": seed_hessian_id,
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
