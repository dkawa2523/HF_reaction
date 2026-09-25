#!/usr/bin/env python3
"""Run one auditable NWChem zero-temperature string from a reaction manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.backends.ts.nwchem_string import NWChemStringEngine
from hfauto.core.io import read_manifest, write_manifest
from hfauto.core.schemas.manifest import Manifest
from hfauto.workflow.reaction_inputs import reaction_and_endpoints


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reaction-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--initial-path", type=Path)
    parser.add_argument("--resample-initial-path", action="store_true")
    parser.add_argument("--functional", default="pbe0")
    parser.add_argument("--basis", default="def2-svpd")
    parser.add_argument("--disp-vdw", type=int)
    parser.add_argument("--required-version")
    parser.add_argument("--method-id", default="nwchem_zero_temperature_string")
    parser.add_argument("--nbeads", type=int, default=7)
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--stepsize", type=float, default=0.05)
    parser.add_argument("--tolerance", type=float, default=0.00045)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=4000)
    parser.add_argument("--timeout-s", type=float, default=21600.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    source_path = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    source = read_manifest(source_path)
    reaction, reactant, product = reaction_and_endpoints(
        source, str(args.reaction_id)
    )
    method = {
        "path_strategy": "reparameterized_double_ended_path",
        "endpoint_basin_status": "distinct_basin",
        "allow_subprocess": True,
        "fallback_to_dummy": False,
        "method_id": str(args.method_id),
        "functional": str(args.functional),
        "basis": str(args.basis),
        "string_nbeads": int(args.nbeads),
        "string_maxiter": int(args.maxiter),
        "string_stepsize": float(args.stepsize),
        "string_tolerance": float(args.tolerance),
        "ncores": int(args.ncores),
        "memory_mb": int(args.memory_mb),
        "timeout_s": float(args.timeout_s),
    }
    if args.initial_path is not None:
        method["initial_path_xyz"] = str(args.initial_path.resolve())
        method["resample_initial_path"] = bool(args.resample_initial_path)
    if args.disp_vdw is not None:
        method["disp_vdw"] = int(args.disp_vdw)
    if args.required_version:
        method["required_program_version"] = str(args.required_version)
    result = NWChemStringEngine().search_ts(
        reaction,
        reactant,
        product,
        method,
        output_dir / "work",
    )
    manifest = Manifest.merge(
        [source],
        run_id=output_dir.name,
        stage="zero-temperature-string",
        metadata={
            "source_manifest": str(source_path),
            "reaction_id": str(args.reaction_id),
            "method": method,
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
