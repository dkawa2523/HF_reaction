#!/usr/bin/env python3
"""Merge independent calculation-evidence manifests for downstream stages."""

from __future__ import annotations

import argparse
from pathlib import Path

from hfauto.core.io import read_manifest, write_manifest
from hfauto.core.schemas.manifest import Manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", action="append", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id")
    parser.add_argument("--stage", default="evidence-merge")
    args = parser.parse_args(argv)

    sources = [read_manifest(path.resolve()) for path in args.manifest]
    output_dir = args.output_dir.resolve()
    merged = Manifest.merge(
        sources,
        run_id=args.run_id or output_dir.name,
        stage=str(args.stage),
        metadata={
            "source_manifests": [str(path.resolve()) for path in args.manifest]
        },
    )
    path = write_manifest(merged, output_dir)
    print(f"sources={len(sources)}; artifacts={len(merged.artifacts)}; manifest={path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
