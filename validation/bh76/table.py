"""The BH76 table of the energy-layer choice (P0e; docs/roadmap.md changes 11 and 15).

usage: python validation/bh76/table.py [MEASURED_JSON]   (writes validation/bh76/table.md)

Every P0e row of subset.yaml is an electronic barrier dE = E(TS) - E(zero), with no ZPE, thermal
or spin-orbit term, computed for each energy-layer method at two structure sets: the pipeline's
PBE0-D3BJ/def2-SVPD stationary structures, which are the chain points (thermo.reaction_points)
of the case reaction that cases.yaml links to the row, and the GMTKN55 geometries. g =
|dE//PBE0 - dE//GMTKN55| is the geometry share of the error. A point sums its states, and a
state takes the lowest energy of its minima (fast conformer equilibria, as for the state G).

``run_points`` reads a run. ``render`` turns measured.json (from geometry_sp.py) into the table
and applies the rule (RULE), which was declared before any value was computed. Only the
selection rows choose; the hold-out rows are guidance only.
"""

from __future__ import annotations

import json
import re
import sys
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any, NamedTuple

import yaml

from hfauto.chemistry import thermo as th
from hfauto.chemistry.topology import fragments
from hfauto.chemistry.xyz import XYZ, hill_formula, read_xyz
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Geometry
from hfauto.core.manifest import Manifest
from hfauto.core.method import MethodSpec, level_mismatches
from hfauto.core.records import ArtifactType, MinimumRecord, ReactionRecord, SpeciesRecord
from hfauto.core.system import SystemConfig
from hfauto.pipeline.layout import RunLayout
from hfauto.reporting.validation import load_cases, load_references

HERE = Path(__file__).resolve().parent
CURRENT = "m06-2x-d3_def2-tzvpd"  # energy_method of configs/pipelines
CANDIDATES = ("revm06_def2-tzvpd", "wb97x-d3_def2-tzvpd")
COST_MAX = 2.0
RULE = ("A candidate replaces the current layer only if (1) its selection MUE is <= the current "
        "layer's in every class (HTBH38, NHTBH38), (2) every open-shell weak complex gets its "
        "energy (no scf_unavailable), and (3) its cost per point is <= 2x the current layer's. "
        "Otherwise the trade-off is stated and the current layer stays. Hold-out rows never "
        "choose.")

Point = list[list[str]]  # states, each the structure ids of its minima
Structure = dict[str, Any]  # file, charge, multiplicity (+ sha256, layer_hartree)
SP = Mapping[str, Mapping[str, Any]]  # structure id -> one method's single point


class Row(NamedTuple):
    key: str  # <reaction>.<direction> of subset.yaml
    gmtkn: str  # GMTKN55 row number(s)
    ref: float  # GMTKN55 reference, kcal/mol
    minnesota: float  # Minnesota/08 REF1
    cls: str  # HTBH38 or NHTBH38
    selection: bool
    end: str  # the zero: R, R_sep, P or P_sep
    bh76: dict[str, Point]  # ts and zero at the GMTKN55 geometries

    @property
    def reaction(self) -> str:
        return self.key.partition(".")[0]


def rows(subset: Path = HERE / "subset.yaml") -> list[Row]:
    raw = yaml.safe_load(subset.read_text(encoding="utf-8"))
    references, selection = load_references(subset), set(raw["split"]["selection"])
    out = []
    for name, reaction in raw["reactions"].items():
        directions = {d: r for d, r in reaction.items() if isinstance(r, dict) and "ref_kcal" in r}
        numbers = [n for r in directions.values() for n in r.get("rows", [r.get("row")])]
        cls = re.search(r"\((N?HTBH38)\)", reaction["type"])
        for direction, r in directions.items():
            verbatim = r["verbatim"] if isinstance(r["verbatim"], str) else r["verbatim"][0]
            systems, coefficients = (part.split() for part in verbatim.split("|")[1:3])
            stoich = [(f"bh76/{s}", int(c)) for s, c in zip(systems, coefficients, strict=True)]
            ref = references[f"{name}.{direction}"]
            first = r.get("rows", [r.get("row")])[0]
            out.append(Row(
                f"{name}.{direction}", ref.row, ref.ref_kcal,
                reaction["minnesota_08"]["REF1"][numbers.index(first)],
                cls[1] if cls else "?", name in selection,
                ("P" if direction == "reverse" else "R")
                + ("_sep" if ref.zero.startswith("separated") else ""),
                {"ts": [[s] for s, c in stoich for _ in range(c)],
                 "zero": [[s] for s, c in stoich for _ in range(-c)]}))
    return out


def links(cases: Path = HERE.parent / "cases.yaml") -> dict[str, list[tuple[str, str]]]:
    """BH76 reaction -> the (case, equation) pairs that cases.yaml links to it, in file order."""
    out: dict[str, list[tuple[str, str]]] = {}
    for name, case in load_cases(cases).items():
        for e in case.reactions:
            if e.bh76:
                out.setdefault(e.bh76.partition(".")[0], []).append((name, e.equation))
    return out


def _equation(rx: ReactionRecord, minima: Mapping[str, MinimumRecord], view: Manifest,
              load: Callable[[Geometry], XYZ]) -> str:
    """The case signature of reporting.validation: fragment formulas of the optimized ends."""
    def side(ids: Iterable[str]) -> str:
        found = []
        for i in ids:
            x = load(view.evidence(minima[i].opt_calc).final)
            found += [hill_formula([x.symbols[j] for j in g]) for g in fragments(x.symbols, x.coords)]
        return " + ".join(sorted(found))

    return f"{side(rx.monomers or rx.minima[:1])} -> {side(rx.minima[1:])}"


def _layer(resolved: Mapping[str, Any], view: Manifest) -> dict[str, float]:
    """Subject -> the run's own energy-layer single point (thermo energy_method), in Hartree."""
    method_id = next((s.get("energy_method") for s in resolved["pipeline"]["stages"]
                      if s["stage"] == "thermo"), None)
    if method_id is None:
        return {}
    method = MethodSpec.model_validate(resolved["methods"][method_id])
    return {subject: ev.energy_hartree for (subject, _), (_, ev)
            in th.single_points(view.of(ArtifactType.CALCULATION)).items()
            if not level_mismatches(method, ev.level, version_pin=ev.level.version)}


def run_points(run: Path, equation: str) -> tuple[dict[str, Point], dict[str, Structure]] | str:
    """The chain points (ts, R, P, and R_sep / P_sep where an end's fragments are declared
    monomers) of the run's reaction with ``equation`` and a saddle, with their PBE0 structures
    (freq final geometries) by id ``<case>/<minimum or freq calc id>``; why not, as a string."""
    layout = RunLayout(run)
    view = layout.view()
    minima = {m.minimum_id: m for m in view.records(ArtifactType.MINIMUM, MinimumRecord)
              if m.tier == "dft"}
    freq = {i: view.evidence(m.freq_calc) for i, m in minima.items()}

    def load(geometry: Geometry) -> XYZ:
        return read_xyz(layout.run_dir / geometry.file.path)

    rx = next((r for r in view.records(ArtifactType.REACTION, ReactionRecord)
               if r.saddle is not None and _equation(r, minima, view, load) == equation), None)
    if rx is None:
        return f"no {equation} reaction with a saddle"
    resolved = next(iter(yaml.safe_load(
        layout.resolved_config_path.read_text(encoding="utf-8")).values()))  # first pipeline
    system = SystemConfig.model_validate(resolved["system"])
    monomers = th.declared_monomers(view.records(ArtifactType.SPECIES, SpeciesRecord),
                                    system.compositions)
    state_of = {i: (m.composition_id, m.state_label) for i, m in minima.items()}
    points = th.reaction_points(rx, state_of, {i: th.separated_states(
        load(f.final), f.level.charge, monomers) for i, f in freq.items()})
    if not points.chain:
        return f"{equation}: {rx.outcome}, no chain"
    t = next(i for i, pt in enumerate(points.chain) if not pt.states)
    before, after = points.chain[:t], points.chain[t + 1:]
    named = {"ts": points.chain[t], "R_sep": points.separated, "P": after[0],
             "R": points.complex if points.separated else (before[-1] if before else None),
             "P_sep": after[1] if len(after) > 1 else None}
    ids = {name: [[f"{run.name}/{i}" for i, s in state_of.items() if s == state]
                  for state in pt.states] or [[f"{run.name}/{pt.subject}"]]
           for name, pt in named.items() if pt is not None}
    used = {i for point in ids.values() for state in point for i in state}
    layer = _layer(resolved, view)
    evidence = {**freq, rx.saddle.freq_calc: view.evidence(rx.saddle.freq_calc)}
    return ids, {key: {"file": str(layout.run_dir / ev.final.file.path),
                       "charge": ev.level.charge, "multiplicity": ev.level.multiplicity,
                       "layer_hartree": layer.get(i)}
                 for i, ev in evidence.items() if (key := f"{run.name}/{i}") in used}


def energy(point: Point, sp: SP) -> float | None:
    """Hartree: the sum over states of the lowest energy of their structures; None when one
    structure lacks its energy (it may be the lowest)."""
    total = 0.0
    for state in point:
        E = [sp.get(s, {}).get("energy_hartree") for s in state]
        if not E or None in E:
            return None
        total += min(E)
    return total


def barrier(points: Mapping[str, Point], sp: SP) -> float | None:
    """kcal/mol from the zero to the TS."""
    ts, zero = energy(points["ts"], sp), energy(points["zero"], sp)
    return None if ts is None or zero is None else (ts - zero) * HARTREE_TO_KCAL_MOL


# ---- values and the rule --------------------------------------------------------------------


class Value(NamedTuple):
    row: Row
    case: str  # the case whose PBE0 structures give ``pbe0``
    counted: bool  # the row's value in the MUE: its first case with points
    pbe0: float | None  # kcal/mol at the PBE0 structures
    geometry: float | None  # kcal/mol at the GMTKN55 geometries
    note: str  # why there is no pbe0 value


def values(measured: Mapping[str, Any], table: list[Row], method: str) -> list[Value]:
    sp, out = measured["sp"][method], []
    for row in table:
        geometry, found, first = barrier(row.bh76, sp), [], None
        for case, named in measured["points"].get(row.reaction, {}).items():
            if isinstance(named, str) or row.end not in named:
                note = named if isinstance(named, str) else f"no {row.end} point"
                found.append(Value(row, case, False, None, geometry, note))
                continue
            first = len(found) if first is None else first
            dE = barrier({"ts": named["ts"], "zero": named[row.end]}, sp)
            found.append(Value(row, case, False, dE, geometry,
                               "" if dE is not None else "single point failed"))
        found = found or [Value(row, "", False, None, geometry, "no linked case")]
        found[first or 0] = found[first or 0]._replace(counted=True)
        out += found
    return out


def mue(found: Iterable[Value], *, selection: bool, cls: str | None = None,
        at_geometry: bool = False) -> tuple[float | None, int, int]:
    """(MUE in kcal/mol, rows with a value, rows) over the counted values of one set (and
    class), at the PBE0 structures or at the GMTKN55 geometries."""
    picked = [v for v in found if v.counted and v.row.selection == selection
              and cls in (None, v.row.cls)]
    errors = [abs(x - v.row.ref) for v in picked
              if (x := v.geometry if at_geometry else v.pbe0) is not None]
    return (sum(errors) / len(errors) if errors else None), len(errors), len(picked)


def verdict(v: Value, tolerance: Mapping[str, float]) -> str:
    """pass: |error| <= per_row_kcal and g <= geometry_share_max_kcal (subset.yaml); a larger g
    is geometry-limited, never a pass."""
    if v.pbe0 is None:
        return v.note
    if v.geometry is None:
        return "no g"
    if abs(v.pbe0 - v.geometry) > tolerance["geometry_share_max_kcal"]:
        return "geometry-limited"
    return "pass" if abs(v.pbe0 - v.row.ref) <= tolerance["per_row_kcal"] else "fail"


def cost(measured: Mapping[str, Any], method: str) -> tuple[float | None, float, float, int]:
    """(summed wall-time ratio, lowest and highest per-point ratio, points) against the current
    layer over the structures both converged."""
    a, b = measured["sp"][method], measured["sp"][CURRENT]
    common = [s for s in a if "energy_hartree" in a[s] and "energy_hartree" in b.get(s, {})]
    if not common:
        return None, 0.0, 0.0, 0
    ratios = [a[s]["seconds"] / b[s]["seconds"] for s in common]
    total = sum(a[s]["seconds"] for s in common) / sum(b[s]["seconds"] for s in common)
    return total, min(ratios), max(ratios), len(common)


def _held(measured: Mapping[str, Any], method: str) -> int:
    """The open-shell weak complexes the method gives an energy."""
    return sum("energy_hartree" in measured["sp"][method].get(s, {})
               for s in measured["complexes"].values())


def decide(measured: Mapping[str, Any], table: list[Row]) -> tuple[str, list[str]]:
    """The layer RULE chooses, and each candidate's criteria with its and the current layer's
    values: the trade-off when it is not adopted."""
    current = values(measured, table, CURRENT)
    classes = sorted({r.cls for r in table if r.selection})
    n, adopted, lines = len(measured["complexes"]), [], []
    for method in CANDIDATES:
        found = values(measured, table, method)
        mues = [(c, mue(found, selection=True, cls=c), mue(current, selection=True, cls=c))
                for c in classes]
        ratio = cost(measured, method)[0]
        met = (all(a[0] is not None and b[0] is not None and a[1] == a[2] == b[1]
                   and a[0] <= b[0] for _, a, b in mues),
               _held(measured, method) == n, ratio is not None and ratio <= COST_MAX)
        verdicts = ["met" if m else "not met" for m in met]
        lines.append(
            f"- `{method}`: **{'adopted' if all(met) else 'not adopted'}**. (1) {verdicts[0]}:"
            " selection MUE " + ", ".join(f"{c} {_f(a[0])} (current {_f(b[0])})"
                                          for c, a, b in mues)
            + f"; (2) {verdicts[1]}: {_held(measured, method)}/{n} complexes (current"
            f" {_held(measured, CURRENT)}/{n}); (3) {verdicts[2]}: {_f(ratio)}x.")
        if all(met):
            adopted.append((mue(found, selection=True)[0] or 0.0, method))
    return (min(adopted)[1] if adopted else CURRENT), lines


# ---- the table ------------------------------------------------------------------------------


def _f(value: float | None, digits: int = 2) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


def _row_table(found: list[Value], tolerance: Mapping[str, float]) -> list[str]:
    out = [("| row | GMTKN55 | set | class | case | ref | Minn/08 | s_ref | dE//PBE0 | error |"
            " dE//geom | error//geom | g | verdict |"), "|" + "---|" * 14]
    for v in found:
        r = v.row
        err = None if v.pbe0 is None else v.pbe0 - r.ref
        geom_err = None if v.geometry is None else v.geometry - r.ref
        g = None if v.pbe0 is None or v.geometry is None else abs(v.pbe0 - v.geometry)
        case = v.case if v.counted else f"{v.case} (not counted)"
        out.append(f"| {r.key} | {r.gmtkn} | {'selection' if r.selection else 'hold-out'} |"
                   f" {r.cls} | {case} | {r.ref:.1f} | {r.minnesota:.2f} |"
                   f" {abs(r.ref - r.minnesota):.2f} | {_f(v.pbe0)} | {_f(err)} |"
                   f" {_f(v.geometry)} | {_f(geom_err)} | {_f(g)} | {verdict(v, tolerance)} |")
    return out


def _mue_table(measured: Mapping[str, Any], table: list[Row]) -> list[str]:
    classes = sorted({r.cls for r in table if r.selection})

    def cell(m: tuple[float | None, int, int]) -> str:
        return f"{_f(m[0])} ({m[1]}/{m[2]})"

    out = ["| method | " + " | ".join(f"selection {c}" for c in classes)
           + " | selection | selection //geom | hold-out | hold-out //geom |",
           "|" + "---|" * (len(classes) + 5)]
    for method in (CURRENT, *CANDIDATES):
        found = values(measured, table, method)
        cells = [mue(found, selection=True, cls=c) for c in classes] + [
            mue(found, selection=True), mue(found, selection=True, at_geometry=True),
            mue(found, selection=False), mue(found, selection=False, at_geometry=True)]
        out.append(f"| `{method}` | " + " | ".join(cell(m) for m in cells) + " |")
    return out


def _scf_tables(measured: Mapping[str, Any]) -> list[str]:
    out = ["| method | complex | PBE0 <S2> | attempts | iterations | <S2> | wall s | result |",
           "|---|---|---|---|---|---|---|---|"]
    for method in (CURRENT, *CANDIDATES):
        for label, structure in measured["complexes"].items():
            s = measured["sp"][method].get(structure, {})
            result = (_f(s["energy_hartree"], 6) if "energy_hartree" in s
                      else s.get("failure", "not computed"))
            freq_s2 = measured["structures"][structure].get("freq_s2")
            out.append(f"| `{method}` | {label} | {_f(freq_s2, 4)} | {s.get('attempts', '—')} |"
                       f" {s.get('iterations', '—')} | {_f(s.get('s2'), 4)} |"
                       f" {_f(s.get('seconds'), 0)} | {result} |")
    out += ["", "| method | single points | rung 1 | rungs 2-3 | failed |",
            "|---|---|---|---|---|"]
    for method in (CURRENT, *CANDIDATES):
        sp = list(measured["sp"][method].values())
        ok = [s for s in sp if "energy_hartree" in s]
        out.append(f"| `{method}` | {len(sp)} | {sum(s['attempts'] == 1 for s in ok)} |"
                   f" {sum(s['attempts'] > 1 for s in ok)} | {len(sp) - len(ok)} |")
    return out


def _cost_table(measured: Mapping[str, Any]) -> list[str]:
    out = ["| method | points | summed wall time / current | per point |", "|---|---|---|---|"]
    for method in CANDIDATES:
        total, low, high, n = cost(measured, method)
        out.append(f"| `{method}` | {n} | {_f(total)} | {_f(low)}–{_f(high)} |")
    return out


def _layer_agreement(measured: Mapping[str, Any]) -> str:
    """The current layer from the atomic guess against the runs' own layer single points."""
    sp = measured["sp"][CURRENT]
    diffs = [abs(sp[i]["energy_hartree"] - s["layer_hartree"]) * HARTREE_TO_KCAL_MOL
             for i, s in measured["structures"].items()
             if s.get("layer_hartree") is not None and "energy_hartree" in sp.get(i, {})]
    return (f"`{CURRENT}` from the atomic guess against the runs' own layer single points (from"
            f" the freq vectors): max |ΔE| {_f(max(diffs, default=None), 3)} kcal/mol over"
            f" {len(diffs)} structures.")


_HEADER = """\
# BH76 table: the energy-layer choice

Generated by `validation/bh76/table.py` from `measured.json` (written by `geometry_sp.py`); do
not edit.

- References: GMTKN55 BH76 and Minnesota/08 REF1 as copied into `subset.yaml` (P0e; {sources}).
- dE = E(TS) − E(zero) in kcal/mol, electronic only (no ZPE, thermal or spin-orbit term), at the
  pipeline's PBE0-D3BJ/def2-SVPD stationary structures (dE//PBE0) and at the GMTKN55 geometries
  (dE//geom); g = |dE//PBE0 − dE//geom|.
- Tolerance (subset.yaml): |error| ≤ {per_row_kcal} against GMTKN55 and
  g ≤ {geometry_share_max_kcal}; a larger g is geometry-limited, never a pass.
  s_ref = |GMTKN55 − Minnesota/08|. A row's MUE value is its first linked case with points.
- Single points: `geometry_sp.py`, hfauto's NWChem engine from the atomic guess with the SCF
  ladder (attempts 1: rung 1; 2: rungs 2-3), one job at a time on the site wsl_local (4 MPI
  ranks).
- Runs: {runs}.
- {agreement}
"""


def render(measured: Mapping[str, Any]) -> str:
    """table.md: the rows of subset.yaml with measured.json's single points, and the decision."""
    sources = json.loads((HERE / "SOURCES.json").read_text(encoding="utf-8"))["archives"]
    table, raw = rows(), yaml.safe_load((HERE / "subset.yaml").read_text(encoding="utf-8"))
    tolerance = {"per_row_kcal": raw["tolerance"]["per_row_kcal"], "geometry_share_max_kcal":
                 raw["tolerance"]["composition"]["geometry_share_max_kcal"]}
    chosen, lines = decide(measured, table)
    out = [_HEADER.format(
        sources=", ".join(f"{a['name']} {a['sha256'][:16]}" for a in sources),
        runs="; ".join(f"{c} `{r['path']}` ({r['code_version'][:12]})"
                       for c, r in measured["runs"].items()),
        agreement=_layer_agreement(measured), **tolerance)]
    for method in (CURRENT, *CANDIDATES):
        out += [f"## `{method}`", "", *_row_table(values(measured, table, method), tolerance), ""]
    out += ["## MUE (kcal/mol; rows with a value / rows)", "", *_mue_table(measured, table), "",
            "## SCF on the open-shell weak complexes and on every single point", "",
            *_scf_tables(measured), "", "## Cost against the current layer", "",
            *_cost_table(measured), "", "## Decision", "",
            f"Rule (declared before any value): {RULE}", "", *lines, "",
            f"Energy layer: **`{chosen}`**" + (" (unchanged)." if chosen == CURRENT else ".")]
    return "\n".join(out) + "\n"


def main(argv: list[str]) -> int:
    path = Path(argv[0]) if argv else HERE / "measured.json"
    text = render(json.loads(path.read_text(encoding="utf-8")))
    (HERE / "table.md").write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
