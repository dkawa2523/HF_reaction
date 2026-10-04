"""Run the validation cases (validation/cases.yaml): replays of recorded runs, or fresh runs.

usage: python validation/replay.py OUT [--src DIR] [--runs A,B] [--strict] [--cap SECONDS]

OUT is a directory under /home/user/hfauto_r10 (a relative name is taken there). Each case's
chain runs in OUT/<name> with the working tree, one hfauto run at a time under the QM lock
/home/user/hfauto_r10/.qm.lock; a pipeline still running after --cap seconds gets SIGTERM (its
stage stays incomplete). Superseded cases are not run.

--src DIR  replay: DIR/<name> is copied, the chain's first pipeline is re-staged from its first
           stage and the later pipelines resume. A replay passes when every re-staged stage ran
           no new job (misses 0), hit the cache where the source stage ran jobs, ends in the
           source's status, and leaves the same records (the artifacts except calculations;
           reason strings and report.html are not compared). --strict: exit 1 on a failure.
no --src   fresh runs, a twice case also in OUT/<name>.2; validation/check.py compares them.

Writes OUT/summary.json (merged over calls into one OUT) and OUT/logs/<name>.log. An existing
OUT without summary.json is refused.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess  # noqa: TID251 (a harness outside hfauto: it starts hfauto itself)
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from hfauto.reporting.validation import Case, load_cases

REPO = Path(__file__).resolve().parents[1]
ROOT = Path("/home/user/hfauto_r10")
LOCK = ROOT / ".qm.lock"
SITE = "configs/sites/wsl_local.yaml"
_NOT_COMPARED = frozenset({"reason", "reasons", "report.html"})


def _config(name: str, kind: str) -> str:
    """A path relative to the repository: a name is a file under configs/<kind>."""
    return name if "/" in name else f"configs/{kind}/{name}.yaml"


def run_chain(case: Case, run: Path, log: Path, cap: int, restage: bool) -> list[int]:
    """Run the chain's pipelines in turn in ``run``; their exit codes."""
    rcs: list[int] = []
    for i, step in enumerate(case.chain):
        pipeline = _config(step.pipeline, "pipelines")
        argv = ["flock", str(LOCK), "timeout", "-k", "60", str(cap), sys.executable, "-m",
                "hfauto.cli.main", "run", pipeline, "--system", _config(step.system, "systems"),
                "--site", SITE, "--run-dir", str(run), *step.args]
        if restage and i == 0:
            first = yaml.safe_load((REPO / pipeline).read_text(encoding="utf-8"))["stages"][0]
            argv += ["--from", first["id"]]
        with log.open("a", encoding="utf-8") as out:
            out.write(f"$ {' '.join(argv)}\n")
            out.flush()
            rcs.append(subprocess.run(argv, cwd=REPO, stdout=out, stderr=subprocess.STDOUT,
                                      check=False).returncode)
        if rcs[-1] not in (0, 1):  # 2: a stage failed or the input was refused; else stopped
            break
    return rcs


def _states(run: Path) -> dict[str, dict[str, Any]]:
    return {s["stage_id"]: s for s in json.loads((run / "run_state.json").read_text())}


def stage_problems(src: Path, dest: Path, t0: datetime) -> tuple[list[str], list[str]]:
    """What breaks the replay gate per stage, and the stages re-staged after ``t0``."""
    before, after = _states(src), _states(dest)
    problems, replayed = [], []
    for sid in [*before, *(s for s in after if s not in before)]:
        b, a = before.get(sid), after.get(sid)
        ran = b["jobs"]["hits"] + b["jobs"]["misses"] if b else 0
        if not (a and a["started"] and datetime.fromisoformat(a["started"]) >= t0):
            if ran:
                problems.append(f"{sid}: ran {ran} jobs in the source but was not re-staged")
            continue
        replayed.append(sid)
        if b is None or a["status"] != b["status"]:
            problems.append(f"{sid}: status {b and b['status']} -> {a['status']}")
        if a["jobs"]["misses"]:
            problems.append(f"{sid}: {a['jobs']['misses']} new jobs")
        if ran and not a["jobs"]["hits"]:
            problems.append(f"{sid}: no cache hit although the source stage ran {ran} jobs")
    return problems, replayed


def _strip(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _strip(v) for k, v in value.items() if k not in _NOT_COMPARED}
    if isinstance(value, list):
        return [_strip(v) for v in value]
    return value


def records(run: Path, stage_id: str) -> dict[str, Any]:
    """The stage's artifacts except calculations (covered by the misses), by artifact id."""
    path = run / stage_id / "manifest.json"
    data = json.loads(path.read_text()) if path.is_file() else {"artifacts": []}
    return {a["artifact_id"]: _strip(a) for a in data["artifacts"] if a["type"] != "calculation"}


def diff(a: Any, b: Any, path: str) -> list[str]:
    """Field differences of two JSON values as 'path: old -> new'. A field the source lacks,
    now at an empty default, is a schema addition, not a changed conclusion."""
    if isinstance(a, dict) and isinstance(b, dict):
        added = {k for k in b.keys() - a.keys() if b[k] in ([], None)}
        return [d for k in sorted((a.keys() | b.keys()) - added)
                for d in diff(a.get(k, "<absent>"), b.get(k, "<absent>"), f"{path}.{k}")]
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [d for i, (x, y) in enumerate(zip(a, b, strict=True))
                for d in diff(x, y, f"{path}[{i}]")]
    return [] if a == b else [f"{path}: {json.dumps(a)[:60]} -> {json.dumps(b)[:60]}"]


def record_diffs(src: Path, dest: Path, stage_ids: list[str]) -> list[str]:
    out: list[str] = []
    for sid in stage_ids:
        a, b = records(src, sid), records(dest, sid)
        for aid in sorted(a.keys() | b.keys()):
            out += ([f"{sid}/{aid}: {'added' if aid in b else 'removed'}"]
                    if aid not in a or aid not in b else diff(a[aid], b[aid], f"{sid}/{aid}"))
    return out


def replay(case: Case, src: Path, dest: Path, log: Path, cap: int) -> dict[str, Any]:
    shutil.copytree(src, dest, symlinks=True)
    t0 = datetime.now(UTC).replace(microsecond=0)  # run_state times have whole seconds
    rcs = run_chain(case, dest, log, cap, restage=True)
    problems, replayed = stage_problems(src, dest, t0)
    diffs = record_diffs(src, dest, replayed)
    return {"verdict": "FAIL" if problems or diffs else "PASS", "rcs": rcs,
            "problems": problems, "diffs": diffs}


def _plan(args: argparse.Namespace, cases: dict[str, Case]) -> dict[str, list[Path]]:
    """The run dirs of each case to run; SystemExit naming every case that cannot run."""
    names = args.runs.split(",") if args.runs else [n for n, c in cases.items() if not c.superseded]
    unknown = [n for n in names if n not in cases]
    if unknown:
        raise SystemExit(f"no such case: {', '.join(unknown)}")
    out = ROOT / args.out
    plan = {n: [out / n, *([out / f"{n}.2"] if cases[n].twice and not args.src else [])]
            for n in names}
    errors = _refusals(out, plan, cases, args.src)
    if errors:
        raise SystemExit("\n".join(errors))
    return plan


def _refusals(out: Path, plan: dict[str, list[Path]], cases: dict[str, Case],
              src: Path | None) -> list[str]:
    errors = [f"{n}: superseded" for n in plan if cases[n].superseded]
    if out.exists() and not (out / "summary.json").is_file():
        errors.append(f"{out} exists and was not made by this tool")
    errors += [f"{d} exists" for dirs in plan.values() for d in dirs if d.exists()]
    if src:
        errors += [f"{n}: no source run in {src}" for n in plan
                   if not (src / n / "run_state.json").is_file()]
    return errors


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("--src", type=Path)
    ap.add_argument("--runs")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--cap", type=int, default=7200)
    args = ap.parse_args(argv)
    cases = load_cases(REPO / "validation" / "cases.yaml")
    plan = _plan(args, cases)
    out = ROOT / args.out
    (out / "logs").mkdir(parents=True, exist_ok=True)
    git = ["git", "-C", str(REPO), "describe", "--always", "--dirty"]
    code = subprocess.run(git, capture_output=True, text=True, check=False).stdout.strip()
    path = out / "summary.json"  # merged over the calls into one OUT
    summary: dict[str, Any] = json.loads(path.read_text()) if path.is_file() else {}
    failed = False
    for name, dirs in plan.items():
        log = out / "logs" / f"{name}.log"
        result = (replay(cases[name], args.src / name, dirs[0], log, args.cap) if args.src else
                  {"verdict": "RAN", "rcs": [run_chain(cases[name], d, log, args.cap, False)
                                             for d in dirs]})
        summary[name] = {"code": code, "src": str(args.src or ""), **result}
        path.write_text(json.dumps(summary, indent=1) + "\n")
        print(f"{result['verdict']:5s} {name} rc {result['rcs']}", flush=True)
        for line in [*result.get("problems", []), *result.get("diffs", [])[:15]]:
            print(f"    {line}")
        failed = failed or result["verdict"] == "FAIL"
    return 1 if failed and args.strict else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
