"""Single points of the BH76 table (table.py) through hfauto's NWChem engine.

usage: python validation/bh76/geometry_sp.py RUNS_DIR OUT

Every energy-layer method (table.CURRENT and table.CANDIDATES) is computed at the GMTKN55
geometries (xyz/), at the PBE0 structures of the chain points of the cases linked to a BH76 row,
and at the open-shell weak complexes (COMPLEXES). RUNS_DIR/<case> are the cases' run dirs, as
for validation/check.py. Every job starts from the atomic guess and climbs the engine's SCF
ladder (W2: attempt 1 is rung 1, attempt 2 rungs 2-3). Jobs run one at a time under the QM lock
and are cached in OUT/jobs. OUT/measured.json is the table's data, written after every job:
copy it to validation/bh76/ and run table.py.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import yaml
from table import CANDIDATES, CURRENT, HERE, Structure, links, rows, run_points

from hfauto.backends.engines import create
from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.topology import fragments
from hfauto.chemistry.xyz import XYZ, Molecule, hill_formula, read_xyz
from hfauto.core.evidence import Failure
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord
from hfauto.execution.jobs import JobRunner
from hfauto.execution.jobstore import JobStore
from hfauto.pipeline.config import load_method, load_site
from hfauto.pipeline.layout import RunLayout

REPO = HERE.parents[1]
LOCK = Path("/home/user/hfauto_r10/.qm.lock")  # the QM lock of validation/replay.py
SITE = REPO / "configs" / "sites" / "wsl_local.yaml"
# the open-shell weak complexes: case -> the fragment formulas of its lowest such PBE0 minimum
COMPLEXES = {"h_c2h4": "C2H4 + H", "s5_ch3_o2": "CH3 + O2", "s6_oh_ch4": "CH4 + HO"}
_DIIS = re.compile(r"^ d=\s*\d+,ls=", re.MULTILINE)  # one SCF iteration line of NWChem 7.2.3
_CGMIN = re.compile(r"^\s+\d+\s+-\d+\.\d+\s+\d\.\d+D[-+]\d+\s+\d\.\d+D[-+]\d+\s+\d+\.\d+\s*$",
                    re.MULTILINE)  # one cgmin iteration line


def gmtkn55() -> dict[str, Structure]:
    """The GMTKN55 structures of subset.yaml, with their charge and multiplicity."""
    raw = yaml.safe_load((HERE / "subset.yaml").read_text(encoding="utf-8"))["geometries"]
    return {f"bh76/{name}": {"file": f"validation/bh76/xyz/{name}.xyz", "charge": g["q"],
                             "multiplicity": g["m"]} for name, g in raw.items()}


def run_info(run: Path) -> dict[str, str]:
    resolved = yaml.safe_load((run / "resolved_config.yaml").read_text(encoding="utf-8"))
    return {"path": str(run.resolve()), "code_version": next(iter(resolved.values()))[
        "code_version"]}


def complex_structure(run: Path, formula: str) -> tuple[str, Structure] | None:
    """The lowest open-shell PBE0 minimum of the run whose fragments make ``formula``, with
    its freq <S2> (the state the pipeline's layer single point starts from)."""
    layout = RunLayout(run)
    view = layout.view()
    found = []
    for m in view.records(ArtifactType.MINIMUM, MinimumRecord):
        ev = view.evidence(m.freq_calc)
        if m.tier != "dft" or ev.level.multiplicity == 1:
            continue
        x = read_xyz(layout.run_dir / ev.final.file.path)
        parts = sorted(hill_formula([x.symbols[i] for i in g])
                       for g in fragments(x.symbols, x.coords))
        if " + ".join(parts) == formula:
            found.append((ev.energy_hartree, m.minimum_id, ev))
    if not found:
        return None
    _, minimum_id, ev = min(found, key=lambda f: f[0])
    return f"{run.name}/{minimum_id}", {
        "file": str(layout.run_dir / ev.final.file.path), "charge": ev.level.charge,
        "multiplicity": ev.level.multiplicity, "freq_s2": ev.s2}


def single_point(engine: QMEngine, store: JobStore, structure: Structure,
                 method: MethodSpec) -> dict[str, Any]:
    """Energy and <S2>, or the failure, with the job's attempts, SCF iterations and wall
    seconds (every attempt counted)."""
    x = read_xyz(REPO / structure["file"])  # an absolute file stays itself
    mol = Molecule(XYZ([s.capitalize() for s in x.symbols], x.coords), structure["charge"],
                   structure["multiplicity"])
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        result = engine.energy(mol, method)
    key = result.job_key
    attempts = sorted(store.job_dir(key).glob("attempt_*")) if key else []
    texts = [(a / "stdout.txt").read_text(encoding="utf-8", errors="replace")
             for a in attempts if (a / "stdout.txt").is_file()]
    seconds = sum(json.loads((a / "command_result.json").read_text(encoding="utf-8"))
                  ["duration_s"] for a in attempts if (a / "command_result.json").is_file())
    stats = {"attempts": len(attempts), "seconds": round(seconds, 1),
             "iterations": sum(len(_DIIS.findall(t)) + len(_CGMIN.findall(t)) for t in texts)}
    if isinstance(result, Failure):
        return {**stats, "failure": f"{result.kind}: {result.reason}"}
    return {**stats, "energy_hartree": result.energy_hartree, "s2": result.s2}


def collect(runs: Path) -> dict[str, Any]:
    """measured.json without its single points: runs, structures, the points the rows read
    and complexes."""
    structures, points, info, complexes = gmtkn55(), {}, {}, {}
    read: dict[str, set[str]] = {}
    for row in rows():
        read.setdefault(row.reaction, {"ts"}).add(row.end)
    for reaction, linked in links().items():
        for case, equation in linked:
            run = runs / case
            found = (run_points(run, equation) if (run / "run_state.json").is_file()
                     else "no run")
            if isinstance(found, str):
                points.setdefault(reaction, {})[case] = found
                continue
            named = {k: p for k, p in found[0].items() if k in read[reaction]}
            used = {i for p in named.values() for state in p for i in state}
            points.setdefault(reaction, {})[case] = named
            structures |= {k: s for k, s in found[1].items() if k in used}
            info[case] = run_info(run)
    for case, formula in COMPLEXES.items():
        found = complex_structure(runs / case, formula)
        if found is not None:
            complexes[formula] = found[0]
            structures[found[0]] = {**found[1], **structures.get(found[0], {})}
            info[case] = run_info(runs / case)
    for s in structures.values():
        s["sha256"] = hashlib.sha256((REPO / s["file"]).read_bytes()).hexdigest()
    return {"runs": info, "structures": structures, "points": points, "complexes": complexes}


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    runs, out = Path(argv[0]), Path(argv[1])
    measured: dict[str, Any] = {**collect(runs), "sp": {}}
    site = load_site(SITE)
    store = JobStore(out / "jobs")
    engine = create(Capability.QM, "nwchem", jobs=JobRunner(store, cores=site.cores),
                    site=site.engines["nwchem"])
    assert isinstance(engine, QMEngine)
    for method_id in (CURRENT, *CANDIDATES):
        method = load_method(REPO / "configs" / "methods" / f"{method_id}.yaml")
        done = measured["sp"].setdefault(method_id, {})
        for name, structure in measured["structures"].items():
            done[name] = single_point(engine, store, structure, method)
            (out / "measured.json").write_text(json.dumps(measured, indent=1) + "\n",
                                               encoding="utf-8")
            print(method_id, name, done[name], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
