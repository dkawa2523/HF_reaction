"""Validation cases (validation/cases.yaml) and their comparison with a run's records.

A case states what one run must conclude. A reaction is signed by the equation of its end
fragments (``CH4 + HO -> CH3 + H2O``: ``topology.fragments`` of the optimized end minima; an
association starts from its separated monomers), its outcome and its dG_eff. Runs that stop
before reactions are signed by their minima, refusals and discovery products. ``compare`` lists
what a run does not meet. ``deviations`` lists outcomes that contradict a literature reference
(validation/bh76): they are reported next to the verdict and never count as a pass. The
numeric reference comparison (electronic barriers) is not made here: a dG_eff is never compared
with a reference dE.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict

from hfauto.chemistry.topology import fragments
from hfauto.chemistry.xyz import XYZ, hill_formula
from hfauto.core.evidence import Geometry
from hfauto.core.manifest import Manifest
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    StandardState,
)

LoadXYZ = Callable[[Geometry], XYZ]
Value = float | tuple[float, float] | None  # a value, a [low, high] band, or no value
_FROZEN = ConfigDict(frozen=True, extra="forbid")
_WITH_BARRIER = frozenset({CaseOutcome.ELEMENTARY_STEP, CaseOutcome.DEGENERATE})


class Step(BaseModel):
    """One pipeline of a run's chain; later steps resume the same run dir."""

    model_config = _FROZEN
    pipeline: str  # a name under configs/pipelines, or a path relative to the repository
    system: str  # a name under configs/systems, or a path relative to the repository
    args: tuple[str, ...] = ()


class Expected(BaseModel):
    model_config = _FROZEN
    equation: str
    outcome: CaseOutcome
    dG_eff: Value  # kcal/mol at the case's (T_K, standard_state)
    bh76: str | None = None  # "<reaction>.<direction>" of validation/bh76/subset.yaml


class Case(BaseModel):
    model_config = _FROZEN
    chain: tuple[Step, ...]
    reactions: tuple[Expected, ...] = ()
    minima: tuple[str, ...] = ()  # fragment formulas of minima the run holds
    refusals: tuple[str, ...] = ()  # reason categories of refused records
    products: int | None = None  # discovery products
    T_K: float = 298.15
    standard_state: StandardState = "1atm"
    tolerance: float = 0.01  # kcal/mol around a dG_eff value
    twice: bool = False  # a boundary system: a fresh validation runs it twice
    superseded: str | None = None  # evidence only, no new run: what replaces it


class Reference(BaseModel):
    """A literature barrier: electronic energy only, measured from ``zero``."""

    model_config = _FROZEN
    row: str
    ref_kcal: float
    zero: str


def load_cases(path: Path) -> dict[str, Case]:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return {name: Case.model_validate(case) for name, case in raw.items()}


def load_references(path: Path) -> dict[str, Reference]:
    """The rows of a reference subset file (validation/bh76/subset.yaml) by
    ``<reaction>.<direction>``."""
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return {
        f"{name}.{direction}": Reference(
            row="/".join(str(r) for r in row.get("rows", [row.get("row")])),
            ref_kcal=row["ref_kcal"], zero=row["zero"])
        for name, reaction in raw["reactions"].items()
        for direction, row in reaction.items()
        if isinstance(row, dict) and "ref_kcal" in row
    }


def _formulas(minima: Iterable[MinimumRecord], view: Manifest, load_xyz: LoadXYZ) -> str:
    """Sorted fragment formulas of the minima's optimized structures, joined by ' + '."""
    out = []
    for minimum in minima:
        xyz = load_xyz(view.evidence(minimum.opt_calc).final)
        out += [hill_formula([xyz.symbols[i] for i in group])
                for group in fragments(xyz.symbols, xyz.coords)]
    return " + ".join(sorted(out))


def _observed(case: Case, view: Manifest,
              load_xyz: LoadXYZ) -> list[tuple[str, CaseOutcome, float | None]]:
    """(equation, outcome, dG_eff) of every reaction with an outcome."""
    minima = {m.minimum_id: m for m in view.records(ArtifactType.MINIMUM, MinimumRecord)}
    dG = {t.reaction_id: t.dG_eff_kcal
          for t in view.records(ArtifactType.REACTION_THERMO, ReactionThermo)
          if (t.T_K, t.standard_state) == (case.T_K, case.standard_state)}

    def side(ids: Iterable[str]) -> str:
        return _formulas((minima[i] for i in ids), view, load_xyz)

    return [(f"{side(r.monomers or r.minima[:1])} -> {side(r.minima[1:])}", r.outcome,
             dG.get(r.reaction_id))
            for r in view.records(ArtifactType.REACTION, ReactionRecord) if r.outcome is not None]


def _key(value: Value) -> float:
    if value is None:
        return -math.inf
    return sum(value) / 2 if isinstance(value, tuple) else value


def _meets(got: float | None, want: Value, tolerance: float) -> bool:
    if want is None or got is None:
        return want is got
    if isinstance(want, tuple):
        return want[0] <= got <= want[1]
    return abs(got - want) <= tolerance


def _reaction_problems(case: Case, observed: list[tuple[str, CaseOutcome, float | None]]
                       ) -> list[str]:
    """Pairs expected and observed reactions by (equation, outcome), then by dG_eff order."""
    want: defaultdict[tuple[str, CaseOutcome], list[Value]] = defaultdict(list)
    got: defaultdict[tuple[str, CaseOutcome], list[float | None]] = defaultdict(list)
    for e in case.reactions:
        want[e.equation, e.outcome].append(e.dG_eff)
    for equation, outcome, value in observed:
        got[equation, outcome].append(value)
    problems = []
    for key in sorted(want.keys() | got.keys()):
        label = f"{key[0]} ({key[1]})"
        expected, seen = sorted(want[key], key=_key), sorted(got[key], key=_key)
        if len(expected) != len(seen):
            problems.append(f"{label}: {len(seen)} reaction(s), expected {len(expected)}")
            continue
        problems += [f"{label}: dG_eff {g}, expected {w}" for w, g in zip(expected, seen)
                     if not _meets(g, w, case.tolerance)]
    return problems


def compare(case: Case, view: Manifest, load_xyz: LoadXYZ) -> list[str]:
    """What the run (the view of its done stages) does not meet; empty when it passes."""
    problems = _reaction_problems(case, _observed(case, view, load_xyz))
    if case.minima:
        held = {_formulas([m], view, load_xyz)
                for m in view.records(ArtifactType.MINIMUM, MinimumRecord)}
        problems += [f"minimum {f}: absent" for f in case.minima if f not in held]
    refused = {a.failure.reason.partition(":")[0] for a in view.artifacts
               if a.failure is not None and a.type != ArtifactType.CALCULATION}
    problems += [f"refusal {r}: absent" for r in case.refusals if r not in refused]
    if case.products is not None:
        n = sum(d.outcome == "product"
                for d in view.records(ArtifactType.DISCOVERY, DiscoveryRecord))
        if n != case.products:
            problems.append(f"discovery products: {n}, expected {case.products}")
    return problems


def deviations(case: Case, view: Manifest, load_xyz: LoadXYZ,
               references: Mapping[str, Reference]) -> list[str]:
    """Referenced reactions whose observed outcome has no barrier: the reference has one."""
    outcomes: defaultdict[str, set[CaseOutcome]] = defaultdict(set)
    for equation, outcome, _ in _observed(case, view, load_xyz):
        outcomes[equation].add(outcome)
    out = []
    for e in case.reactions:
        seen = outcomes.get(e.equation, set())
        if e.bh76 is None or not seen or seen & _WITH_BARRIER:
            continue
        ref = references[e.bh76]
        out.append(f"{e.equation}: {', '.join(sorted(seen))}; BH76 row {ref.row} has a barrier"
                   f" of {ref.ref_kcal} kcal/mol (zero: {ref.zero})")
    return out
