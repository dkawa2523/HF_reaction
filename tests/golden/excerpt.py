"""Regenerate tests/golden/data and SOURCES.json from runs/ (design §10.2; runs/ is only read).
Usage: python tests/golden/excerpt.py [RUNS_DIR]. Each file's `method` names a rule below.
Absolute sources (G25/G26, G27/G28 and G29: WSL job dirs under /home/user/hfauto_r6/m0/golden,
/home/user/hfauto_r6/m1/golden and /home/user/hfauto_r6/SA/sac_golden) are read as they are;
regenerate those on WSL."""

import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_NW_LINE = re.compile(
    r"Northwest Computational Chemistry Package|nwchem branch|Grid used for XC"
    r"|Convergence on energy requested|\((spherical|cartesian)\)|DFT-D3|<S2>"
    r"|Total DFT energy|Charge +:|Spin multiplicity:|^ string:|^@"
)
_NW_STOP = re.compile("Optimization converged|Failed to converge")  # kept with 6 lines around
_NW_BLOCKS = (  # (start, end): keep from a start line through the next end line
    ("echo of input deck", r"^={20,}\s*$"),
    (r'Summary of "ao basis"|XC Information', r"^\s*$"),
    ("Vibrational analysis via the FX method", "Task  times|NWChem Input Module"),
)


def _end(lines: list[str], i: int, end: str) -> int:
    return next((k for k in range(i + 1, len(lines)) if re.search(end, lines[k])), len(lines) - 1)


def nwchem(text: str) -> str:
    """Keep level, energy, geometry, vibration and convergence lines; drop SCF and gradients."""
    lines = text.splitlines(keepends=True)
    keep = {i for i, s in enumerate(lines) if _NW_LINE.search(s)}
    spans = [(len(lines) - 50, len(lines))]
    for start, end in _NW_BLOCKS:
        spans += [(i, _end(lines, i, end)) for i, s in enumerate(lines) if re.search(start, s)]
    geoms = [i for i, s in enumerate(lines) if re.match(r'\s+Geometry "', s)]
    spans += [(i, _end(lines, i, "Atomic Mass")) for i in {geoms[0], geoms[-1]}] if geoms else []
    spans += [(i - 6, i + 6) for i, s in enumerate(lines) if _NW_STOP.search(s)]
    keep.update(k for a, b in spans for k in range(max(a, 0), min(b + 1, len(lines))))
    return "".join(lines[k] for k in sorted(keep))


def readuct_vib_blocks(text: str) -> str:
    lines = text.splitlines(keepends=True)
    starts = [i for i, s in enumerate(lines) if "Vib. Frequencies:" in s]
    return "".join("".join(lines[i : _end(lines, i + 2, r"^\s*$") + 1]) for i in starts)


def jsonl_imaginary_rows(text: str) -> str:
    rows = [r for r in text.splitlines(keepends=True) if r.strip()]
    return "".join(r for r in rows if json.loads(r).get("imaginary_mode_count") is not None)


def head400_tail200(text: str) -> str:
    lines = text.splitlines(keepends=True)
    return "".join(lines[:400] + lines[max(400, len(lines) - 200) :])


_RULES = (nwchem, readuct_vib_blocks, jsonl_imaginary_rows, head400_tail200)
METHODS = {"copy": None} | {f.__name__: f for f in _RULES}


def main(runs: Path) -> None:
    sources = HERE / "SOURCES.json"
    entries = json.loads(sources.read_text(encoding="utf-8"))
    for f in (f for e in entries for f in e["files"] if not f.get("missing")):
        raw = (runs / f["source"]).read_bytes()
        rule = METHODS[f["method"]]  # None: verbatim copy
        out = rule(raw.decode("utf-8", "replace")).encode() if rule else raw
        name = f["golden"].removesuffix(".gz")
        if len(out) > 50 * 1024:  # stored gzipped; mtime 0 keeps the bytes reproducible
            out, name = gzip.compress(out, mtime=0), name + ".gz"
        (HERE / "data" / name).parent.mkdir(parents=True, exist_ok=True)
        (HERE / "data" / name).write_bytes(out)
        f.update(golden=name, source_sha256=hashlib.sha256(raw).hexdigest())
    sources.write_text(json.dumps(entries, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1] / "runs")
