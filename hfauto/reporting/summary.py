"""Report tables (design §8.2 report): ranking with ties, discovery coverage, method panel.

Imports only hfauto.core, hfauto.chemistry.gates and the standard library. Everything is a
pure function over typed records except ``write_tables``. No confidence score is produced
(AR-27): a reaction is either rankable by ``gates.rankable`` or listed with its blockers.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import astuple, dataclass, fields
from pathlib import Path
from typing import Literal

from hfauto.chemistry.gates import Policy, rankable, reaction_tier
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

Metric = Literal["dG_act", "dG_rxn"]
SIGN_DISAGREEMENT = "method_sign_disagreement"
_SIGN_DEADBAND_KCAL = 1.0  # a dE_rxn within +-1 kcal/mol has no sign to disagree on
_DEFAULT = Policy()


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
    sign_disagreement: bool


def _value(thermo: ReactionThermo, metric: Metric) -> float | None:
    return thermo.dG_act_kcal if metric == "dG_act" else thermo.dG_rxn_kcal


def _interval(thermo: ReactionThermo, value: float, metric: Metric) -> tuple[float, float]:
    """The sensitivity band (absolute dG_act interval) or the point value itself."""
    if metric == "dG_act" and thermo.band_kcal is not None:
        return thermo.band_kcal
    return value, value


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


def rank_rows(
    reactions: Iterable[ReactionRecord],
    reaction_thermo: Iterable[ReactionThermo],
    participant_notes: Mapping[str, Sequence[str]],
    T: float,
    state: str,
    *,
    metric: Metric = "dG_act",
    policy: Policy = _DEFAULT,
) -> list[RankRow]:
    """Rank the rankable reactions by ``metric`` at (T, state); the rest keep rank None.

    Reactions without an outcome are unprocessed hypotheses and are not listed.
    """
    thermo = {
        t.reaction_id: t
        for t in reaction_thermo
        if math.isclose(t.T_K, T, abs_tol=1e-6) and t.standard_state == state
    }
    rows: list[RankRow] = []
    values: dict[str, float] = {}
    intervals: dict[str, tuple[float, float]] = {}
    for reaction in reactions:
        if reaction.outcome is None:
            continue
        rid, t = reaction.reaction_id, thermo.get(reaction.reaction_id)
        notes = participant_notes.get(rid, ())
        gate = rankable(reaction, t, participant_notes=notes, policy=policy)
        value = _value(t, metric) if t is not None else None
        if gate and t is not None and value is not None:
            values[rid], intervals[rid] = value, _interval(t, value, metric)
        rows.append(
            RankRow(
                reaction_id=rid,
                outcome=reaction.outcome,
                tier=reaction_tier(reaction),
                rankable=rid in values,
                rank=None,
                dG_act_kcal=t.dG_act_kcal if t is not None else None,
                band_kcal=t.band_kcal if t is not None else None,
                blockers=gate.reasons,
                T_K=T,
                standard_state=state,
                dG_rxn_kcal=t.dG_rxn_kcal if t is not None else None,
            )
        )
    ranks = _tie_ranks(intervals)
    order = {rid: (ranks[rid], values[rid]) for rid in ranks}
    rows = [row.model_copy(update={"rank": ranks.get(row.reaction_id)}) for row in rows]
    return sorted(rows, key=lambda r: (order.get(r.reaction_id, (math.inf, 0.0)), r.reaction_id))


def participant_notes(
    reactions: Iterable[ReactionRecord],
    minima: Mapping[str, MinimumRecord],
    disagreements: Collection[str] = (),
) -> dict[str, tuple[str, ...]]:
    """Notes of each reaction's minima and saddle, plus the method-panel sign blocker."""
    out: dict[str, tuple[str, ...]] = {}
    for reaction in reactions:
        notes = [n for m in dict.fromkeys(reaction.minima) if m in minima for n in minima[m].notes]
        notes += reaction.saddle.notes if reaction.saddle is not None else ()
        if reaction.reaction_id in disagreements:
            notes.append(SIGN_DISAGREEMENT)
        out[reaction.reaction_id] = tuple(dict.fromkeys(notes))
    return out


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
        "sign_disagreement": any(v > _SIGN_DEADBAND_KCAL for v in rxn)
        and any(v < -_SIGN_DEADBAND_KCAL for v in rxn),
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
) -> tuple[list[PanelRow], frozenset[str]]:
    """dE_rxn and dE_act of each reaction at every level with energies at its stationary
    points, keyed by ``Level.full_key`` (CH-25), and the reactions whose dE_rxn changes sign
    between levels: above +1 kcal/mol at one level and below -1 at another
    (``_SIGN_DEADBAND_KCAL``), so near-zero reaction energies never block ranking.

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
    disagreements = frozenset(r.reaction_id for r in rows if r.sign_disagreement)
    return rows, disagreements


_RANK_HEADER = (
    "rank",
    "reaction_id",
    "outcome",
    "tier",
    "rankable",
    "T_K",
    "standard_state",
    "dG_act_kcal",
    "band_low_kcal",
    "band_high_kcal",
    "dG_rxn_kcal",
    "blockers",
)


def _rank_cells(row: RankRow) -> tuple[object, ...]:
    low, high = row.band_kcal or (None, None)
    cells = (row.rank, row.reaction_id, row.outcome.value, row.tier, row.rankable, row.T_K,
             row.standard_state)
    return (*cells, row.dG_act_kcal, low, high, row.dG_rxn_kcal, ";".join(row.blockers))


def write_tables(
    out_dir: Path,
    rows: Sequence[RankRow],
    coverage_rows: Sequence[CoverageRow],
    panel_rows: Sequence[PanelRow],
) -> dict[str, Path]:
    """Write ranking.csv, coverage.csv and method_panel.csv; None becomes an empty cell."""
    tables = {
        "ranking.csv": (_RANK_HEADER, map(_rank_cells, rows)),
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
