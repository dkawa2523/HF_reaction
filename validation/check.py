"""Compare run directories with their validation cases (validation/cases.yaml).

usage: python validation/check.py RUNS_DIR [NAME ...] [--require-complete]

Every case (or the NAMEs) is compared with its requested stages and their records. PASS means
the expectations matched without a known reference contradiction; it is not a general accuracy
certificate. An unfinished pipeline is INCOMPLETE, a mismatch is FAIL, and a BH76 outcome
contradiction is DEVIATION. All three exit 1. A missing run also fails unless superseded.

A twice case's second fresh run is checked when present. --require-complete also fails when
it is missing; use it for the final fresh validation. Without it an absent second run is
explicitly NOT_CHECKED, allowing comparison of single-run archives and cache replays.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping
from pathlib import Path

import yaml

from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.core.evidence import Geometry
from hfauto.pipeline.layout import RunLayout
from hfauto.reporting.validation import (
    Case,
    Reference,
    compare,
    deviations,
    load_cases,
    load_references,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def _execution_problems(layout: RunLayout, case: Case) -> list[str]:
    """Required stages, respecting each chain step's --from/--to, must have finished."""
    if case.superseded:  # historical evidence; its replaced pipeline need not exist any more
        return []
    states = {s.stage_id: s for s in layout.read_state()}
    required: dict[str, None] = {}
    for step in case.chain:
        path = REPO / (step.pipeline if "/" in step.pipeline else
                       f"configs/pipelines/{step.pipeline}.yaml")
        ids = [s["id"] for s in yaml.safe_load(path.read_text(encoding="utf-8"))["stages"]]
        first = (ids.index(step.args[step.args.index("--from") + 1])
                 if "--from" in step.args else 0)
        last = (ids.index(step.args[step.args.index("--to") + 1])
                if "--to" in step.args else len(ids) - 1)
        required.update(dict.fromkeys(ids[first:last + 1]))
    problems = []
    for sid in required:
        state = states.get(sid)
        if state is None or state.status != "done" or state.unfinished():
            label = "not run" if state is None else state.status
            if state is not None and state.unfinished():
                label += " (work remains)"
            problems.append(f"{sid}: {label}")
    return problems


def check(run: Path, name: str, case: Case, references: Mapping[str, Reference]) -> bool:
    layout = RunLayout(run)
    unfinished = _execution_problems(layout, case)
    if unfinished:
        print(f"INCOMPLETE {name}")
        for line in unfinished:
            print(f"    {line}")
        return False
    view = layout.view()

    def load(geometry: Geometry) -> XYZ:
        return read_xyz(layout.run_dir / geometry.file.path)

    problems = compare(case, view, load)
    reference_problems = deviations(case, view, load, references)
    note = f"  (superseded: {' '.join(case.superseded.split())})" if case.superseded else ""
    verdict = "FAIL" if problems else "DEVIATION" if reference_problems else "PASS"
    if reference_problems and not problems:
        note += "  (case expectations matched; reference contradicts the outcome)"
    print(f"{verdict} {name}{note}")
    for line in problems:
        print(f"    {line}")
    for line in reference_problems:
        print(f"    DEVIATION {line}")
    return not problems and not reference_problems


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", type=Path)
    ap.add_argument("names", nargs="*")
    ap.add_argument("--require-complete", action="store_true",
                    help="require a second fresh run for every twice case")
    args = ap.parse_args(argv)
    runs, names = args.runs, args.names
    cases = load_cases(HERE / "cases.yaml")
    references = load_references(HERE / "bh76" / "subset.yaml")
    unknown = sorted(set(names) - set(cases))
    if unknown:
        print(f"no such case: {', '.join(unknown)}")
        return 2
    ok = True
    for name in names or cases:
        case = cases[name]
        if not (runs / name / "run_state.json").is_file():
            print(f"MISSING {name}")
            ok = ok and bool(case.superseded)
            continue
        ok = check(runs / name, name, case, references) and ok
        second = runs / f"{name}.2"  # the second fresh run of a twice case
        if case.twice:
            if (second / "run_state.json").is_file():
                ok = check(second, second.name, case, references) and ok
            else:
                print(f"{'MISSING' if args.require_complete else 'NOT_CHECKED'} {second.name}"
                      " (second fresh run)")
                ok = ok and not args.require_complete
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
