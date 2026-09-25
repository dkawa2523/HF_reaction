#!/usr/bin/env python3
"""Run one auditable fixed-geometry NWChem Hessian calculation."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.backends.qm.nwchem import NWChemEngine
from hfauto.core.io import write_json_atomic, write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--species-id", default="fixed_geometry_hessian")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--functional", default="pbe0")
    parser.add_argument("--basis", default="def2-svpd")
    parser.add_argument("--disp-vdw", type=int)
    parser.add_argument("--required-version")
    parser.add_argument("--method-id", default="nwchem_fixed_geometry_hessian")
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=4000)
    parser.add_argument("--timeout-s", type=float, default=7200.0)
    parser.add_argument("--temperature-k", type=float, default=298.15)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    xyz = args.xyz.resolve()
    if not xyz.is_file():
        raise FileNotFoundError(xyz)
    output_dir = args.output_dir.resolve()
    result_path = args.result.resolve()
    species = Artifact(
        artifact_id=str(args.species_id),
        artifact_type="species",
        paths={"xyz": str(xyz)},
        data={
            "species_id": str(args.species_id),
            "state": "fixed_geometry_hessian",
            "xyz_path": str(xyz),
            "charge": int(args.charge),
            "multiplicity": int(args.multiplicity),
        },
        qc={"scientific_role": "fixed_geometry_hessian_target"},
        provenance={"created_by": "scripts.run_nwchem_frequency"},
    )
    method = {
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
    if args.disp_vdw is not None:
        method["disp_vdw"] = int(args.disp_vdw)
    if args.required_version:
        method["required_program_version"] = str(args.required_version)
    calculation = NWChemEngine().frequency(species, method, str(output_dir))
    write_json_atomic(result_path, calculation.model_dump())
    manifest = Manifest.new(
        run_id=result_path.parent.name,
        stage="fixed-geometry-hessian",
        metadata={"method": method, "target_xyz": str(xyz)},
    )
    manifest.extend([species, calculation])
    manifest_path = write_manifest(manifest, result_path.parent)
    print(
        f"{calculation.status.status}: {calculation.artifact_id}; "
        f"n_imag={calculation.data.get('n_imag')}; manifest={manifest_path}"
    )
    return 0 if calculation.status.status == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
