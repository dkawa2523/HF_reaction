"""Report tables (design §8.2 report): ranking by dG_eff with ties, discovery coverage,
method panel.

Imports only hfauto.core, hfauto.chemistry.gates and the standard library. Everything is a
pure function over typed records except ``write_tables``. No confidence score is produced
(AR-27): a reaction is either rankable by ``gates.rankable`` or listed with its blockers; the
method panel's dE_act spread is shown as columns, never as a blocker.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import astuple, dataclass, fields
from pathlib import Path

from hfauto.chemistry.gates import rankable, reaction_tier
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Level
from hfauto.core.manifest import Artifact
from hfauto.core.records import (
    DiscoveryRecord,
    MinimumRecord,
    RankRow,
    ReactionRecord,
    ReactionThermo,
)

_FROM_THERMO = {"dG_eff_kcal", "dG_act_kcal", "dG_rxn_kcal", "dG_act_vs_separated_kcal",
                "band_kcal", "notes"}  # RankRow fields copied from the ReactionThermo


@dataclass(frozen=True)
class CoverageRow:
    metric: str  # attempts | products | negatives | failed | negative_reason | failure_kind
    key: str  # mechanism, reason category (text before ':') or FailureKind
    count: int


@dataclass(frozen=True)
class PanelRow:
    """One reaction at one level of theory; the min/max columns span the reaction's levels."""

    reaction_id: str
    level_key: str  # Level.full_key(): grid and scf_tol are part of the key (CH-25)
    level: str
    dE_rxn_kcal: float
    dE_act_kcal: float | None
    dE_rxn_min_kcal: float
    dE_rxn_max_kcal: float
    dE_act_min_kcal: float | None
    dE_act_max_kcal: float | None


def _tie_ranks(intervals: Mapping[str, tuple[float, float]]) -> dict[str, int]:
    """Competition ranks (1, 1, 3) where overlapping intervals share a rank.

    Overlap is closed transitively: sorted by lower bound, a reaction joins the current
    group while its lower bound does not exceed the group's highest upper bound.
    """
    ranks: dict[str, int] = {}
    rank, reach = 0, -math.inf
    ordered = sorted(intervals.items(), key=lambda item: (item[1], item[0]))
    for position, (reaction_id, (low, high)) in enumerate(ordered, start=1):
        if low > reach:
            rank = position
        reach = max(reach, high)
        ranks[reaction_id] = rank
    return ranks


def _row(reaction: ReactionRecord, thermo: ReactionThermo | None, T: float, state: str
         ) -> RankRow:
    gate = rankable(reaction, thermo)
    values = thermo.model_dump(include=_FROM_THERMO) if thermo is not None else {}
    return RankRow(reaction_id=reaction.reaction_id, outcome=reaction.outcome,
                   tier=reaction_tier(reaction), rankable=gate.ok, rank=None,
                   blockers=gate.reasons, T_K=T, standard_state=state,
                   torsional=reaction.torsional, **values)


def rank_rows(reactions: Iterable[ReactionRecord], reaction_thermo: Iterable[ReactionThermo],
              T: float, state: str) -> list[RankRow]:
    """Rank the rankable reactions by dG_eff at (T, state), sharing a rank where sensitivity
    bands overlap; the rest keep rank None. Reactions without an outcome are unprocessed
    hypotheses and are not listed."""
    thermo = {t.reaction_id: t for t in reaction_thermo
              if math.isclose(t.T_K, T, abs_tol=1e-6) and t.standard_state == state}
    rows = [_row(r, thermo.get(r.reaction_id), T, state) for r in reactions
            if r.outcome is not None]
    ranks = _tie_ranks({r.reaction_id: r.band_kcal or (value, value) for r in rows
                        if r.rankable and (value := r.dG_eff_kcal) is not None})
    rows = [row.model_copy(update={"rank": ranks.get(row.reaction_id)}) for row in rows]
    return sorted(rows, key=lambda r: (
        (r.rank, r.dG_eff_kcal or 0.0) if r.rank else (math.inf, 0.0), r.reaction_id))


def coverage(
    discoveries: Iterable[DiscoveryRecord], failed_artifacts: Iterable[Artifact]
) -> list[CoverageRow]:
    """Attempts, products, negatives and failures per mechanism; negatives per reason
    category; failed artifacts per FailureKind."""
    counts: Counter[tuple[str, str]] = Counter()
    outcome_metric = {"product": "products", "negative": "negatives", "failed": "failed"}
    for discovery in discoveries:
        counts["attempts", discovery.mechanism] += 1
        counts[outcome_metric[discovery.outcome], discovery.mechanism] += 1
        if discovery.outcome == "negative":
            counts["negative_reason", (discovery.reason or "none").partition(":")[0]] += 1
    for artifact in failed_artifacts:
        if artifact.status == "failed" and artifact.failure is not None:
            counts["failure_kind", artifact.failure.kind.value] += 1
    return [CoverageRow(metric, key, n) for (metric, key), n in sorted(counts.items())]


def _label(level: Level) -> str:
    parts = (level.program, level.version, level.method, level.basis, level.dispersion)
    text = " ".join(p for p in (*parts, level.solvation, level.grid) if p)
    return text if level.scf_tol is None else f"{text} scf={level.scf_tol:g}"


def _stationary_points(
    reaction: ReactionRecord,
    calculations: Mapping[str, Evidence],
    minima: Mapping[str, MinimumRecord],
) -> tuple[str, str, str | None] | None:
    """Geometry fingerprints of (reactant minimum, product minimum, TS or None)."""
    try:
        start, end = (calculations[minima[m].opt_calc].final.fingerprint for m in reaction.minima)
    except KeyError:  # missing minimum (placeholder '') or calculation outside the view
        return None
    saddle = calculations.get(reaction.saddle.saddle_calc) if reaction.saddle else None
    return start, end, saddle.final.fingerprint if saddle is not None else None


def _deltas(
    energies: Mapping[str, float], points: tuple[str, str, str | None]
) -> tuple[float, float | None] | None:
    start, end, ts = points
    if start not in energies or end not in energies:
        return None
    e_ts = energies.get(ts) if ts is not None else None
    act = (e_ts - energies[start]) * HARTREE_TO_KCAL_MOL if e_ts is not None else None
    return (energies[end] - energies[start]) * HARTREE_TO_KCAL_MOL, act


def _panel_rows(
    reaction_id: str, found: list[tuple[str, Level, float, float | None]]
) -> list[PanelRow]:
    rxn = [f[2] for f in found]
    act = [f[3] for f in found if f[3] is not None]
    spread = {
        "dE_rxn_min_kcal": min(rxn, default=0.0),
        "dE_rxn_max_kcal": max(rxn, default=0.0),
        "dE_act_min_kcal": min(act, default=None),
        "dE_act_max_kcal": max(act, default=None),
    }
    return [
        PanelRow(reaction_id, key, _label(level), d_rxn, d_act, **spread)
        for key, level, d_rxn, d_act in found
    ]


def method_panel(
    calculations: Mapping[str, Evidence],
    reactions: Iterable[ReactionRecord],
    *,
    minima: Mapping[str, MinimumRecord],
) -> list[PanelRow]:
    """dE_rxn and dE_act of each reaction at every level with energies at its stationary
    points, keyed by ``Level.full_key`` (CH-25).

    A calculation contributes its energy at ``final`` to the stationary point with the same
    geometry fingerprint: fixed-geometry single points, and the opt / saddle / freq jobs that
    produced the claim itself (the reference level).
    """
    energies: dict[str, dict[str, float]] = defaultdict(dict)
    levels: dict[str, Level] = {}
    for ev in calculations.values():
        key = ev.level.full_key()
        levels.setdefault(key, ev.level)
        energies[key].setdefault(ev.final.fingerprint, ev.energy_hartree)
    rows: list[PanelRow] = []
    for reaction in reactions:
        points = _stationary_points(reaction, calculations, minima)
        if points is None:
            continue
        deltas = {key: _deltas(energies[key], points) for key in sorted(energies)}
        found = [(key, levels[key], *d) for key, d in deltas.items() if d is not None]
        rows += _panel_rows(reaction.reaction_id, found)
    return rows


_RANK_HEADER = (
    "rank", "reaction_id", "outcome", "tier", "rankable", "T_K", "standard_state", "dG_eff_kcal",
    "band_low_kcal", "band_high_kcal", "dG_act_kcal", "dG_rxn_kcal", "dG_act_vs_separated_kcal",
    "torsional", "blockers", "notes",
)
_PANEL_COLUMNS = ("dE_act_panel_min_kcal", "dE_act_panel_max_kcal")


def _rank_cells(row: RankRow) -> tuple[object, ...]:
    low, high = row.band_kcal or (None, None)
    return (row.rank, row.reaction_id, row.outcome.value, row.tier, row.rankable, row.T_K,
            row.standard_state, row.dG_eff_kcal, low, high, row.dG_act_kcal, row.dG_rxn_kcal,
            row.dG_act_vs_separated_kcal, row.torsional, ";".join(row.blockers),
            ";".join(row.notes))


def _ranking(rows: Sequence[RankRow], panel_rows: Sequence[PanelRow]
             ) -> tuple[tuple[str, ...], list[tuple[object, ...]]]:
    """ranking.csv; with a method panel (energies at two levels or more) each row also gets
    the dE_act min and max over the levels."""
    if len({p.level_key for p in panel_rows}) < 2:
        return _RANK_HEADER, [_rank_cells(row) for row in rows]
    spread = {p.reaction_id: (p.dE_act_min_kcal, p.dE_act_max_kcal) for p in panel_rows}
    return _RANK_HEADER + _PANEL_COLUMNS, [
        (*_rank_cells(row), *spread.get(row.reaction_id, (None, None))) for row in rows]


def write_tables(
    out_dir: Path,
    rows: Sequence[RankRow],
    coverage_rows: Sequence[CoverageRow],
    panel_rows: Sequence[PanelRow],
) -> dict[str, Path]:
    """Write ranking.csv, coverage.csv and method_panel.csv; None becomes an empty cell."""
    tables = {
        "ranking.csv": _ranking(rows, panel_rows),
        "coverage.csv": ([f.name for f in fields(CoverageRow)], map(astuple, coverage_rows)),
        "method_panel.csv": ([f.name for f in fields(PanelRow)], map(astuple, panel_rows)),
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, (header, body) in tables.items():
        paths[name] = out_dir / name
        with paths[name].open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(body)
    return paths
