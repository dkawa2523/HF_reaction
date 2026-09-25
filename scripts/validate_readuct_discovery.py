"""Run one real ReaDuct discovery attempt and save its evidence summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hfauto.backends.reaction_discovery.readuct import ReaDuctDiscoveryBackend
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import (
    DrivingAtomPairRecord,
    ReactionTrialRecord,
)


def _pair(value: str) -> DrivingAtomPairRecord:
    first, second = value.split(",", maxsplit=1)
    return DrivingAtomPairRecord(atoms=(int(first), int(second)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("xyz", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--associate", action="append", default=[])
    parser.add_argument("--dissociate", action="append", default=[])
    parser.add_argument("--driver", choices=("nt2", "afir"), default="nt2")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    args = parser.parse_args()

    source = Artifact(
        artifact_id="validation_source",
        artifact_type="species_preopt",
        paths={"xyz": str(args.xyz.resolve())},
        data={
            "species_id": "validation_source",
            "xyz_path": str(args.xyz.resolve()),
            "charge": args.charge,
            "multiplicity": args.multiplicity,
        },
    )
    trial = ReactionTrialRecord(
        trial_id="validation_trial",
        source_species_id="validation_source",
        associations=[_pair(value) for value in args.associate],
        dissociations=[_pair(value) for value in args.dissociate],
        driver_order=[args.driver],
        charge=args.charge,
        multiplicity=args.multiplicity,
        max_attempts=1,
    )
    backend = ReaDuctDiscoveryBackend()
    result = backend.explore(
        trial,
        source,
        {
            "method_family": "GFN2",
            "nt_max_iterations": 250,
            "ts_max_iterations": 250,
            "irc_max_iterations": 200,
            "endpoint_max_iterations": 250,
        },
        args.out,
        driver=args.driver,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {
        **result.__dict__,
        "capability": backend.capability(),
    }
    (args.out / "discovery_result.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if result.success else 2


if __name__ == "__main__":
    raise SystemExit(main())
