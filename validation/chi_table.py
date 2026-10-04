"""The reaction-mode character χ of every DFT saddle a reaction case checked (QM-free).

usage: python validation/chi_table.py DIR [DIR ...] [--out FILE]

DIR is a run dir or a directory of run dirs. For every done reaction-paths stage and every case
it drove: the case ends as driver.open_case builds them, each DFT saddle job its seeds started
(the case folder's seeds, the low-level NEB's TS, the record's low-level TSs, a stalled search's
last frame) or its ts_calc, and the DFT freq on that saddle, read from the JobStore. χ is
gates.reaction_mode_chi of the freq's first imaginary mode on the case's bond change
(Ctx.change), as validate_ts measures it. A row is connected (C) when a case of the stage claims
that freq with a connected or multi-step outcome.

The summary gives the χ range of the connected TSs; the other first-order saddles at or below
0.065, inside the gap (0.065, 0.69) and from 0.69 up; and each case whose χ rejections differ
from its log's (``ts_rejected:not_reaction_mode``), i.e. a conclusion χ changes. Records of the
older runs (r7, M3, W5) load with two adaptations (``_load``); a run that still does not load is
listed as unreadable.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from pydantic import ValidationError

from hfauto.chemistry import gates
from hfauto.chemistry.geometry import declared_coordinate_gradient
from hfauto.chemistry.identity import member_coords
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz, read_xyz_trajectory
from hfauto.core.evidence import Evidence, Geometry
from hfauto.core.hashing import sha256_file
from hfauto.core.manifest import Manifest
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    CaseOutcome,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.core.records import ArtifactType as T
from hfauto.drivers.minimum import calc_id
from hfauto.drivers.reaction_case.actions import Ctx
from hfauto.drivers.reaction_case.state import CaseRules
from hfauto.pipeline import layout as run_layout
from hfauto.pipeline.layout import RunLayout

CONNECTED = {*CONNECTED_OUTCOMES.values(), CaseOutcome.MULTI_STEP}
LOW, HIGH = 0.065, 0.69  # the empty gap the review declared (U6-P1); a row inside is listed


def _load(path: Path) -> Manifest:
    """A manifest, records of older runs adapted: MinimumRecord.chiral is dropped and a single
    ReactionRecord.low_level_ts becomes a tuple."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    for artifact in data.get("artifacts", []):
        payload = artifact.get("payload") or {}
        if payload.get("kind") == "minimum":
            payload.pop("chiral", None)
        ts = payload.get("low_level_ts")
        if payload.get("kind") == "reaction" and not isinstance(ts, list):
            payload["low_level_ts"] = [] if ts is None else [ts]
    return Manifest.model_validate(data)


run_layout.load_manifest = _load  # RunLayout.view reads every input manifest through it


class Jobs:
    """The run's JobStore: saddle jobs, DFT freqs by start fingerprint, NEB TSs by initial path."""

    def __init__(self, run: Path) -> None:
        self.saddles: list[tuple[str, Geometry | None, Evidence | None]] = []
        self.freqs: dict[str, Evidence] = {}
        self.nebs: dict[str, Geometry] = {}
        for job_json in sorted(run.glob("jobs/*/*/job.json")):
            result_json = job_json.parent / "result.json"
            if not result_json.is_file():
                continue
            job, result = json.loads(job_json.read_text()), json.loads(result_json.read_text())
            kind, data, ok = job["kind"], result["data"], result["kind"] == "calculation"
            if kind == "saddle":
                final = Geometry.model_validate(data["final"]) if data.get("final") else None
                ev = Evidence.model_validate(data) if ok else None
                self.saddles.append((job["key_payload"]["molecule"], final, ev))
            elif kind == "frequencies" and ok and job["key_payload"]["method"]["kind"] == "dft":
                ev = Evidence.model_validate(data)
                self.freqs.setdefault(ev.start.fingerprint, ev)
            elif kind == "neb" and data.get("ts"):  # a path result
                self.nebs[job["key_payload"]["initial_path"]] = Geometry.model_validate(data["ts"])


class Case:
    """One driven case: its context (driver.open_case) and the saddles its seeds reached."""

    def __init__(self, run: Path, out: Manifest, folder: Path, case: ReactionRecord) -> None:
        self.run, self.case, self.folder = run, case, folder
        minima = {m.minimum_id: m for m in out.records(T.MINIMUM, MinimumRecord)}
        species = {s.species_id: s for s in out.records(T.SPECIES, SpeciesRecord)}
        first, ids = species[case.endpoints[0]], self._driven(case, minima)

        def end(m: str, s: str) -> np.ndarray:  # driver._endpoint
            basin = self.load(out.evidence(minima[m].opt_calc).final)
            own = self.load(species[s].geometry).coords
            return member_coords(basin.symbols, basin.coords, minima[m].species_id, s, own)

        collapsed = bool(case.monomers) and ids[0] == ids[1]
        a = self.load(first.geometry).coords if collapsed else end(ids[0], case.endpoints[0])
        b = end(ids[1], case.endpoints[1])
        self.ctx = Ctx(case=case, rt=None, rules=CaseRules(), folder=folder,  # type: ignore[arg-type]
                       symbols=list(first.geometry.symbols), charge=first.charge,
                       multiplicity=first.multiplicity, raw=(a, b), ends=(a, align_mapped(a, b)),
                       energies=(0.0, 0.0), log=print)

    @staticmethod
    def _driven(case: ReactionRecord, minima: dict[str, MinimumRecord]) -> tuple[str, str]:
        """The minima the case was driven between: a reassigned record holds its TS's sides, so
        its endpoint species' own DFT basins (classification.finalize)."""
        if case.source != "reassigned":
            return case.minima
        own = [next(m.minimum_id for m in minima.values() if m.tier == "dft"
                    and (m.species_id == s or s in m.members)) for s in case.endpoints]
        return own[0], own[1]

    def load(self, geometry: Geometry) -> XYZ:
        return read_xyz(self.run / geometry.file.path)

    def key(self, coords: np.ndarray) -> str:
        """The Molecule.fingerprint a job of this case keys a structure by."""
        xyz = XYZ(self.ctx.symbols, np.asarray(coords, dtype=float))
        return Molecule(xyz, self.ctx.charge, self.ctx.multiplicity).fingerprint()

    def saddles(self, jobs: Jobs) -> list[Evidence]:
        """The converged saddles: the ts_calc, then the searches its seeds started, following
        each stalled search to the restart from its last frame."""
        frontier = {self.key(self.load(ts).coords) for ts in self.case.low_level_ts}
        for path in self.folder.glob("*.xyz"):
            frames = read_xyz_trajectory(path)
            if len(frames) == 1:
                frontier.add(self.key(frames[0].coords))
            if (ts := jobs.nebs.get(sha256_file(path))) is not None:
                frontier.add(self.key(self.load(ts).coords))
        found, pending = [], list(jobs.saddles)
        while hits := [s for s in pending if s[0] in frontier]:
            pending = [s for s in pending if s not in hits]
            found += [ev for _, _, ev in hits if ev is not None]
            frontier |= {self.key(self.load(final).coords) for _, final, ev in hits
                         if ev is None and final is not None}
        return found


def rows_of(run: Path) -> list[dict]:
    layout, jobs = RunLayout(run), Jobs(run)
    done = {s.stage_id for s in layout.read_state() if s.status == "done"}
    rows = []
    for config in (yaml.safe_load(layout.resolved_config_path.read_text()) or {}).values():
        policy = gates.Policy(**(config["pipeline"].get("gates") or {}))
        for stage in config["pipeline"]["stages"]:
            if stage["stage"] != "reaction-paths" or stage["id"] not in done:
                continue
            own = _load(layout.manifest_path(stage["id"]))
            out = Manifest.union([layout.view(stage["id"]), own], run_id=run.name,
                                 stage_id=stage["id"])
            cases = own.records(T.REACTION, ReactionRecord)
            claimed = {c.saddle.freq_calc for c in cases if c.saddle and c.outcome in CONNECTED}
            for record in cases:
                if record.log and "" not in record.minima:
                    folder = layout.stage_dir(stage["id"]) / Path(record.log).parent
                    rows += case_rows(Case(run, out, folder, record), out, jobs, policy, claimed)
    return rows


def case_rows(case: Case, out: Manifest, jobs: Jobs, policy: gates.Policy, claimed: set[str]
              ) -> list[dict]:
    record, ctx = case.case, case.ctx
    saddles = ([out.evidence(record.ts_calc)] if record.ts_calc else []) + case.saddles(jobs)
    logged = (case.folder / "log.jsonl").read_text().count('"ts_rejected:not_reaction_mode')
    formed, broken = ctx.change()
    rows = []
    for saddle in saddles:
        freq = jobs.freqs.get(saddle.final.fingerprint)
        if freq is None or not freq.imaginary_modes:
            continue
        x, terms = np.asarray(case.load(freq.start).coords, dtype=float), record.coordinate
        chi = gates.reaction_mode_chi(ctx.symbols, freq.imaginary_modes[0], x, formed | broken,
                                      declared_coordinate_gradient(terms, x) if terms else None)
        nus = sorted(freq.frequencies_cm1 or ())
        gate = gates.is_first_order_saddle(freq, saddle=saddle, policy=policy)
        rows.append({"case": record.reaction_id.removeprefix("rxn_"), "nu1": nus[0],
                     "nu2": nus[1] if len(nus) > 1 else None, "chi": chi, "logged": logged,
                     "gate": "ok" if gate else ",".join(gate.reasons),
                     "connected": calc_id(freq) in claimed,
                     "bonds": " ".join([*(f"+{ctx.symbols[i]}{i}-{ctx.symbols[j]}{j}"
                                          for i, j in sorted(formed)),
                                        *(f"-{ctx.symbols[i]}{i}-{ctx.symbols[j]}{j}"
                                          for i, j in sorted(broken))]) or "none"})
    return rows


def _fmt(r: dict) -> str:
    chi = "-" if r["chi"] is None else f"{r['chi']:.3f}"
    nu2 = "-" if r["nu2"] is None else f"{r['nu2']:.1f}"
    return (f"{r['run']:<32} {r['case']:<34} {r['nu1']:>8.1f} {nu2:>8} {r['gate']:<17} "
            f"{'C' if r['connected'] else '-'} {chi:>6}  {r['bonds']}")


def summary(rows: list[dict]) -> list[str]:
    first = [r for r in rows if r["gate"] == "ok" and r["chi"] is not None]
    conn = [r["chi"] for r in first if r["connected"]]
    other = [r for r in first if not r["connected"]]
    bands = {"<= 0.065": [r for r in other if r["chi"] <= LOW],
             "in the gap": [r for r in other if LOW < r["chi"] < HIGH],
             ">= 0.69": [r for r in other if r["chi"] >= HIGH]}
    out = ["", f"connected first-order saddles: {len(conn)}"
           + (f", chi {min(conn):.3f}..{max(conn):.3f}" if conn else ""),
           "other first-order saddles: " + ", ".join(f"{k} {len(v)}" for k, v in bands.items())]
    out += [f"  in the gap: {_fmt(r)}" for r in bands["in the gap"]]
    out += [f"  connected below {gates.REACTION_MODE_MIN}: {_fmt(r)}" for r in first
            if r["connected"] and r["chi"] < gates.REACTION_MODE_MIN]
    now: dict[tuple[str, str], list[int]] = {}
    for r in rows:
        low = r["gate"] == "ok" and r["chi"] is not None and r["chi"] < gates.REACTION_MODE_MIN
        now.setdefault((r["run"], r["case"]), [r["logged"], 0])[1] += low
    out += [f"  chi rejections {old} in the log, {new} now: {run} {case}"
            for (run, case), (old, new) in now.items() if old != new]
    return out


def runs_in(path: Path) -> list[Path]:
    if (path / "run_state.json").is_file():
        return [path]
    return sorted(p for p in path.iterdir() if (p / "run_state.json").is_file())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("dirs", nargs="+", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    rows, unreadable = [], []
    for run in (r for d in args.dirs for r in runs_in(d)):
        try:
            rows += [{**r, "run": f"{run.parent.name}/{run.name}"} for r in rows_of(run)]
        except ValidationError as exc:  # records of an older schema
            unreadable.append(f"unreadable {run}: {str(exc).splitlines()[0]}")
    rows.sort(key=lambda r: (-1.0 if r["chi"] is None else r["chi"], r["run"], r["case"]))
    head = (f"{'run':<32} {'case':<34} {'nu1':>8} {'nu2':>8} {'gate':<17} C {'chi':>6}  "
            "bond change of the case ends")
    text = "\n".join([head, *map(_fmt, rows), *summary(rows), *unreadable])
    print(text)
    if args.out:
        args.out.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
