"""Report tables (design §8.2 report): ranking by dG_eff with ties, discovery coverage,
method panel.

Imports only hfauto.core, hfauto.chemistry and the standard library. Everything is a pure
function over typed records except ``write_tables``. No confidence score is produced (design
§5.4): a reaction is either rankable by ``gates.rankable`` or listed with its blockers (a
barrierless outcome: capture-limited, by its dG_rxn); the method panel's dE_act spread is shown
as columns, never as a blocker. The panel reads the reaction's own points (thermo.participants,
the ones reaction_points gives dE from). A spin-contaminated panel energy
(``gates.energy_spin_ok``, the thermo stage's definition) stays in its row, noted, and out of
the spread.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import astuple, dataclass, fields
from pathlib import Path

from hfauto.chemistry.gates import Policy, energy_spin_ok, rankable, reaction_tier
from hfauto.chemistry.thermo import Sides, participants, single_points
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

_FROM_THERMO = {"dG_eff_kcal", "reference", "dG_act_kcal", "dG_rxn_kcal",
                "dG_act_vs_separated_kcal", "band_kcal", "energy_level",
                "notes"}  # RankRow fields from the ReactionThermo


@dataclass(frozen=True)
class CoverageRow:
    metric: str  # attempts | products | negatives | failed | negative_reason | failure_kind
    key: str  # mechanism, reason category (text before ':') or FailureKind
    count: int


@dataclass(frozen=True)
class PanelRow:
    """One reaction at one level of theory; the min/max columns span the reaction's levels
    except the values named in ``notes``."""

    reaction_id: str
    level_key: str  # Level.full_key(): grid and scf_tol are part of the key
    level: str
    dE_rxn_kcal: float
    dE_act_kcal: float | None
    dE_rxn_min_kcal: float | None
    dE_rxn_max_kcal: float | None
    dE_act_min_kcal: float | None
    dE_act_max_kcal: float | None
    notes: str = ""  # spin_contaminated:dE_rxn / :dE_act, ';'-joined: left out of the spread


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


def _row(reaction: ReactionRecord, thermo: ReactionThermo | None, T: float | None,
         state: str | None) -> RankRow:
    gate = rankable(reaction, thermo)
    values = thermo.model_dump(include=_FROM_THERMO) if thermo is not None else {}
    return RankRow(reaction_id=reaction.reaction_id, outcome=reaction.outcome,
                   tier=reaction_tier(reaction), rankable=gate.ok, rank=None,
                   blockers=gate.reasons, T_K=T, standard_state=state,
                   torsional=reaction.torsional, **values)


def rank_rows(reactions: Iterable[ReactionRecord], reaction_thermo: Iterable[ReactionThermo],
              T: float | None, state: str | None) -> list[RankRow]:
    """Rank the rankable reactions by dG_eff at (T, state), sharing a rank where sensitivity
    bands overlap; the rest keep rank None (a barrierless step: capture-limited, with its
    dG_rxn). Reactions without an outcome are unprocessed hypotheses and are not listed."""
    thermo = {t.reaction_id: t for t in reaction_thermo if (t.T_K, t.standard_state) == (T, state)}
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
    outcome_metric = {"product": "products", "negative": "negatives", "failed": "failed",
                      "unconnected": "unconnected", "not_attempted": "not_attempted"}
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
    text = " ".join(p for p in (*parts, level.grid) if p)
    return text if level.scf_tol is None else f"{text} scf={level.scf_tol:g}"


def _deltas(at: Mapping[str, Evidence], freq: Mapping[str, Evidence], sides: Sides,
            policy: Policy) -> tuple[float, float | None, tuple[str, ...]] | None:
    """dE_rxn, dE_act (kcal/mol) over the own points and the names of those with a
    spin-contaminated energy: a point's freq fails gates.energy_spin_ok."""
    reactant, product, ts = sides
    if not all(i in at for i in (*reactant, *product)):
        return None
    pairs = {"dE_rxn": product, **({"dE_act": ts} if ts and all(t in at for t in ts) else {})}
    start = sum(at[i].energy_hartree for i in reactant)
    values = {name: (sum(at[i].energy_hartree for i in ids) - start) * HARTREE_TO_KCAL_MOL
              for name, ids in pairs.items()}
    contaminated = tuple(name for name, ids in pairs.items()
                         if not all(energy_spin_ok(at[p], freq[p], policy)
                                    for p in (*reactant, *ids)))
    return values["dE_rxn"], values.get("dE_act"), contaminated


def _panel_rows(
    reaction_id: str, found: list[tuple[str, Level, float, float | None, tuple[str, ...]]]
) -> list[PanelRow]:
    rxn = [f[2] for f in found if "dE_rxn" not in f[4]]
    act = [f[3] for f in found if f[3] is not None and "dE_act" not in f[4]]
    spread = {
        "dE_rxn_min_kcal": min(rxn, default=None),
        "dE_rxn_max_kcal": max(rxn, default=None),
        "dE_act_min_kcal": min(act, default=None),
        "dE_act_max_kcal": max(act, default=None),
    }
    return [
        PanelRow(reaction_id, key, _label(level), d_rxn, d_act, **spread,
                 notes=";".join(f"spin_contaminated:{name}" for name in contaminated))
        for key, level, d_rxn, d_act, contaminated in found
    ]


def _lot(level: Level) -> str:
    """A level of theory as a key, charge and multiplicity aside: separated monomers differ in
    them from their complex by design (gates.same_pes(state=False))."""
    return level.model_copy(update={"charge": 0, "multiplicity": 1}).full_key()


def method_panel(
    calculations: Sequence[Artifact],
    reactions: Sequence[ReactionRecord],
    *,
    minima: Mapping[str, MinimumRecord],
    policy: Policy,
) -> list[PanelRow]:
    """dE_rxn and dE_act of each reaction at every level with energies at its own points
    (thermo.participants), keyed by the product's ``Level.full_key``.

    A point is a subject of the sp stage (a minimum_id, or SaddleClaim.freq_calc for the TS):
    its reference energy comes from its freq calculation, the other levels from the sp
    calculations whose parents name it (thermo.single_points: one per subject and level). A
    value with a spin-contaminated energy (_deltas) is noted and left out of the min/max.
    """
    evidence = {a.artifact_id: a.payload for a in calculations if isinstance(a.payload, Evidence)}
    calcs = {m.minimum_id: m.freq_calc for m in minima.values()} | {
        r.saddle.freq_calc: r.saddle.freq_calc for r in reactions if r.saddle is not None}
    freq = {s: evidence[c] for s, c in calcs.items() if c in evidence}
    energies: dict[str, dict[str, Evidence]] = defaultdict(dict)  # _lot -> subject -> ev
    for (subject, _), (_, ev) in single_points(calculations).items():
        if subject in freq:
            energies[_lot(ev.level)][subject] = ev
    for subject, ev in freq.items():  # the reference level, unless an sp is on it
        energies[_lot(ev.level)].setdefault(subject, ev)
    rows: list[PanelRow] = []
    for reaction in reactions:
        sides = participants(reaction)
        product = sides[1][0]
        found = [(at[product].level.full_key(), at[product].level, *d)
                 for at in energies.values() if (d := _deltas(at, freq, sides, policy))]
        rows += _panel_rows(reaction.reaction_id, sorted(found, key=lambda f: f[0]))
    return rows


_RANK_HEADER = (
    "rank", "reaction_id", "outcome", "tier", "rankable", "T_K", "standard_state",
    "energy_level", "dG_eff_kcal", "reference", "band_low_kcal", "band_high_kcal",
    "dG_act_kcal", "dG_rxn_kcal", "dG_act_vs_separated_kcal", "torsional", "blockers", "notes",
)
_PANEL_COLUMNS = ("dE_act_panel_min_kcal", "dE_act_panel_max_kcal")


def _rank_cells(row: RankRow) -> tuple[object, ...]:
    low, high = row.band_kcal or (None, None)
    return (row.rank, row.reaction_id, row.outcome.value, row.tier, row.rankable, row.T_K,
            row.standard_state, row.energy_level, row.dG_eff_kcal, row.reference, low, high,
            row.dG_act_kcal, row.dG_rxn_kcal, row.dG_act_vs_separated_kcal, row.torsional,
            ";".join(row.blockers), ";".join(row.notes))


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
