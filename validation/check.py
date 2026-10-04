"""Compare run directories with their validation cases (validation/cases.yaml).

usage: python validation/check.py RUNS_DIR [NAME ...]

Every case (or the NAMEs) with a run dir RUNS_DIR/<name>, and RUNS_DIR/<name>.2 when a twice
case has one, is compared through the view of its done stages: PASS, or FAIL with what it does not
meet. DEVIATION lines name outcomes that contradict a BH76 reference; they never count as a
pass. A case without a run dir is MISSING, which fails unless the case is superseded. Exit 1
when a run fails or is missing.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path

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


def check(run: Path, name: str, case: Case, references: Mapping[str, Reference]) -> bool:
    layout = RunLayout(run)
    view = layout.view()

    def load(geometry: Geometry) -> XYZ:
        return read_xyz(layout.run_dir / geometry.file.path)

    problems = compare(case, view, load)
    note = f"  (superseded: {' '.join(case.superseded.split())})" if case.superseded else ""
    print(f"{'FAIL' if problems else 'PASS'} {name}{note}")
    for line in problems:
        print(f"    {line}")
    for line in deviations(case, view, load, references):
        print(f"    DEVIATION {line}")
    return not problems


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    runs, names = Path(argv[0]), argv[1:]
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
        if case.twice and (second / "run_state.json").is_file():
            ok = check(second, second.name, case, references) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
