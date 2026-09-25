#!/usr/bin/env python3
"""Run a real bidirectional pysisyphus/NWChem IRC from a validated TS."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.backends.ts.pysisyphus import PysisyphusEngine
from hfauto.core.io import read_manifest, write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.workflow.reaction_inputs import validated_ts_and_endpoints


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reaction-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--functional", default="pbe0")
    parser.add_argument("--basis", default="def2-svpd")
    parser.add_argument("--disp-vdw", type=int, default=3)
    parser.add_argument("--required-version", default="7.2.3")
    parser.add_argument("--pysis-executable")
    parser.add_argument("--nwchem-executable")
    parser.add_argument("--integrator", default="eulerpc")
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--step-length", type=float, default=0.08)
    parser.add_argument("--rms-gradient", type=float, default=5.0e-4)
    parser.add_argument("--endpoint-maxiter", type=int, default=250)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=4000)
    parser.add_argument("--timeout-s", type=float, default=172800.0)
    return parser


def assemble_result_manifest(
    source: Manifest,
    *,
    run_id: str,
    source_path: Path,
    reaction_id: str,
    validated_reaction_id: str,
    method: dict[str, object],
    artifacts: list[Artifact],
) -> Manifest:
    """Append IRC evidence without dropping upstream calculation evidence."""

    manifest = source.carry_forward("bidirectional-irc")
    manifest.run_id = run_id
    manifest.metadata.update(
        {
            "source_manifest": str(source_path),
            "reaction_id": reaction_id,
            "validated_reaction_id": validated_reaction_id,
            "method": method,
        }
    )
    manifest.extend(artifacts)
    return manifest


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    source_path = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    source = read_manifest(source_path)
    reaction, ts_species, reactant, product, validated = (
        validated_ts_and_endpoints(source, str(args.reaction_id))
    )
    method = {
        "allow_subprocess": True,
        "fallback_to_dummy": False,
        "program": "nwchem",
        "functional": str(args.functional),
        "basis": str(args.basis),
        "disp_vdw": int(args.disp_vdw),
        "required_program_version": str(args.required_version),
        "irc_integrator": str(args.integrator),
        "irc_maxiter": int(args.maxiter),
        "irc_step_length": float(args.step_length),
        "irc_rms_grad_threshold": float(args.rms_gradient),
        "irc_dump_every": 5,
        "optimize_irc_endpoints": True,
        "require_optimized_irc_endpoints": True,
        "endpoint_optimizer": "nwchem",
        "endpoint_opt_threshold": "gau",
        "endpoint_opt_maxiter": int(args.endpoint_maxiter),
        "endpoint_rmsd_threshold_A": 0.35,
        "endpoint_permutation_rmsd_threshold_A": 0.20,
        "endpoint_q_tolerance_A": 0.08,
        "endpoint_distance_spectrum_threshold_A": 0.08,
        "require_identity_invariant_endpoint_match": True,
        "ncores": int(args.ncores),
        "memory_mb": int(args.memory_mb),
        "timeout_s": float(args.timeout_s),
    }
    if args.pysis_executable:
        method["pysis_executable"] = str(args.pysis_executable)
    if args.nwchem_executable:
        method["executable"] = str(args.nwchem_executable)
    result = PysisyphusEngine().run_irc(
        reaction,
        ts_species,
        reactant,
        product,
        method,
        output_dir / "work",
    )
    # The IRC backend contributes only new path evidence; endpoint and TS
    # calculations remain the single source of truth for downstream thermo.
    manifest = assemble_result_manifest(
        source,
        run_id=output_dir.name,
        source_path=source_path,
        reaction_id=str(args.reaction_id),
        validated_reaction_id=validated.artifact_id,
        method=method,
        artifacts=result.artifacts,
    )
    result_path = write_manifest(manifest, output_dir)
    print(
        f"success={result.success}; artifacts={len(result.artifacts)}; "
        f"manifest={result_path}"
    )
    return 0 if result.success else 1


if __name__ == "__main__":
    raise SystemExit(main())
