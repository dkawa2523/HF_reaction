"""Static HTML report (design §4.1 #8, §7.4 ``report``): ``render(view, out_dir, load_xyz=)``.

One self-contained page per view: a summary row per reaction (δG_eff, ΔE‡, ΔG‡, ν_imag,
connection, outcome, blockers, the association quantities and the energy level), an inline-SVG
energy diagram per reaction and the TS imaginary mode animated by 3Dmol.js. The module does not
know the run layout: geometries are reached only through ``load_xyz``. Standard
library, numpy and hfauto.core only; no plotting or templating package.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Sequence
from html import escape
from pathlib import Path

import numpy as np

from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import Geometry
from hfauto.core.manifest import Manifest
from hfauto.core.records import ArtifactType, RankRow, ReactionRecord, ReactionThermo, ReportRecord

THREEDMOL_JS = "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"
MODE_AMPLITUDE_A = 0.5  # largest atomic displacement of the animated imaginary mode (Å)
_DASH = "—"

DiagramLevel = tuple[str, float | None, float | None]  # (label, ΔE, ΔG) kcal/mol

_CSS = """
:root{--ink:#1d2330;--muted:#5b6475;--line:#d8dce3;--bg:#fff;--card:#f5f6f8;
--de:#8a94a6;--dg:#1f6feb}
@media (prefers-color-scheme:dark){:root{--ink:#e6e8ec;--muted:#9aa3b2;--line:#363c47;
--bg:#15181e;--card:#1d2129;--de:#7d8799;--dg:#5ea1ff}}
body{margin:0 auto;max-width:1100px;padding:24px 16px;font:14px/1.5 system-ui,sans-serif;
color:var(--ink);background:var(--bg)}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:0 0 4px}
a{color:var(--dg)}.meta{color:var(--muted);margin:0 0 12px}
.wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;
font-variant-numeric:tabular-nums}th,td{border-bottom:1px solid var(--line);
padding:4px 8px;text-align:left;vertical-align:top}th{font-weight:600}td.n{text-align:right}
section{border-top:1px solid var(--line);margin-top:28px;padding-top:16px}
.figs{display:flex;flex-wrap:wrap;gap:16px;margin-bottom:12px}
figure{margin:0;flex:1 1 360px;max-width:460px}figcaption{color:var(--muted);font-size:12px}
svg.energy{width:100%;height:auto;background:var(--card)}
svg text{fill:var(--ink);font-size:11px}svg .de{stroke:var(--de)}svg .dg{stroke:var(--dg)}
svg .bar{stroke-width:3}svg .link{stroke-width:1;stroke-dasharray:3 3}
.mol{position:relative;width:100%;height:300px;background:#fff}
"""

# 3Dmol reads the 5th-7th xyz columns as per-atom displacements; vibrate() turns them into
# frames between -1 and +1 times the displacement and animate() plays them back and forth.
_VIEWER_JS = """
document.querySelectorAll(".mol").forEach((el) => {
  if (typeof $3Dmol === "undefined") { el.textContent = "3Dmol.js could not be loaded."; return; }
  const viewer = $3Dmol.createViewer(el, {backgroundColor: "white"});
  viewer.addModel(el.dataset.xyz, "xyz");
  viewer.setStyle({}, {stick: {radius: 0.14}, sphere: {scale: 0.25}});
  viewer.vibrate(10, 1, true);
  viewer.animate({loop: "backAndForth", interval: 60});
  viewer.zoomTo();
  viewer.render();
});
"""


def _num(value: float | None, digits: int = 2) -> str:
    return _DASH if value is None else f"{value:.{digits}f}"


def _imag(cm1: float) -> str:
    return f"{abs(cm1):.1f}i"


def _table(header: Sequence[str], rows: Iterable[Sequence[str]], numeric: set[int]) -> str:
    """An HTML table of already escaped cells; columns in ``numeric`` are right-aligned."""
    head = "".join(f"<th>{h}</th>" for h in header)
    body = "".join(
        "<tr>"
        + "".join(f"<td class='n'>{c}</td>" if i in numeric else f"<td>{c}</td>"
                  for i, c in enumerate(row))
        + "</tr>"
        for row in rows
    )
    return f"<div class='wrap'><table><tr>{head}</tr>{body}</table></div>"


def _levels(thermo: ReactionThermo) -> list[DiagramLevel]:
    """Energy-diagram levels relative to the reactant; separated monomers sit at −ΔG_assoc."""
    g_known = any(v is not None for v in (thermo.dG_act_kcal, thermo.dG_rxn_kcal,
                                          thermo.dG_assoc_kcal))
    e_known = thermo.dE_act_kcal is not None or thermo.dE_rxn_kcal is not None
    out: list[DiagramLevel] = []
    if thermo.dG_assoc_kcal is not None:
        out.append(("separated", None, -thermo.dG_assoc_kcal))
    out.append(("reactant", 0.0 if e_known else None, 0.0 if g_known else None))
    for label, e, g in (("TS", thermo.dE_act_kcal, thermo.dG_act_kcal),
                        ("product", thermo.dE_rxn_kcal, thermo.dG_rxn_kcal)):
        if e is not None or g is not None:
            out.append((label, e, g))
    return out


def _energy_svg(points: Sequence[DiagramLevel], caption: str) -> str:
    """Inline SVG: one slot per level, ΔE bars on the left and ΔG bars on the right half."""
    values = [v for _, e, g in points for v in (e, g) if v is not None]
    if not values:
        return f"<p class='meta'>No energies for {escape(caption)}.</p>"
    width, height, top, bottom, left = 460, 240, 38, 200, 12
    high, low = max(values), min(values)
    span = (high - low) or 1.0
    slot = (width - 2 * left) / len(points)

    def y(value: float) -> float:
        return top + (high - value) / span * (bottom - top)

    legend = ("<tspan style='fill:var(--de)'>■ ΔE</tspan> "
              "<tspan style='fill:var(--dg)'>■ ΔG</tspan>")
    parts = [(f"<svg class='energy' viewBox='0 0 {width} {height}' role='img' aria-label='"
              f"energy diagram, kcal/mol'><text x='{left}' y='16'>kcal/mol · "
              f"{escape(caption)} · {legend}</text>")]
    series = {"de": (0.08, [e for _, e, _ in points]), "dg": (0.54, [g for _, _, g in points])}
    for css, (offset, energies) in series.items():
        previous: tuple[float, float] | None = None
        for i, value in enumerate(energies):
            if value is None:
                continue
            x0, yv = left + slot * (i + offset), y(value)
            x1 = x0 + slot * 0.38
            if previous is not None:
                parts.append(f"<line class='link {css}' x1='{previous[0]:.1f}' "
                             f"y1='{previous[1]:.1f}' x2='{x0:.1f}' y2='{yv:.1f}'/>")
            parts.append(f"<line class='bar {css}' x1='{x0:.1f}' y1='{yv:.1f}' "
                         f"x2='{x1:.1f}' y2='{yv:.1f}'/><text x='{(x0 + x1) / 2:.1f}' "
                         f"y='{yv - 5:.1f}' text-anchor='middle'>{value:.1f}</text>")
            previous = (x1, yv)
    parts += [f"<text x='{left + slot * (i + 0.5):.1f}' y='{height - 8}' "
              f"text-anchor='middle'>{escape(label)}</text>"
              for i, (label, _, _) in enumerate(points)]
    return "".join(parts) + "</svg>"


def _mode_xyz(view: Manifest, reaction: ReactionRecord,
             load_xyz: Callable[[Geometry], XYZ]) -> str | None:
    """The TS of ``reaction`` with its first imaginary mode as xyz displacement columns,
    scaled so that the largest atomic displacement is MODE_AMPLITUDE_A; None without one."""
    if reaction.saddle is None:
        return None
    try:
        evidence = view.evidence(reaction.saddle.freq_calc)
    except (KeyError, ValueError):  # the freq calculation is not in this view
        return None
    if not evidence.imaginary_modes:
        return None
    xyz = load_xyz(evidence.final)
    mode = np.asarray(evidence.imaginary_modes[0], dtype=float).reshape(-1, 3)
    mode = mode * MODE_AMPLITUDE_A / max(float(np.linalg.norm(mode, axis=1).max()), 1e-12)
    lines = [str(len(xyz.symbols)), f"{reaction.reaction_id} TS"]
    lines += [f"{s} {x:.6f} {y:.6f} {z:.6f} {dx:.6f} {dy:.6f} {dz:.6f}"
              for s, (x, y, z), (dx, dy, dz) in zip(xyz.symbols, xyz.coords, mode, strict=True)]
    return "\n".join(lines) + "\n"


def _connection(reaction: ReactionRecord) -> str:
    claim = reaction.connection
    if claim is None:
        return "not confirmed" if reaction.saddle is not None else _DASH
    return escape(f"QRC: {claim.minima[0]} ↔ {claim.minima[1]}")


def _band(thermo: ReactionThermo) -> str:
    return _DASH if thermo.band_kcal is None else "{:.2f}–{:.2f}".format(*thermo.band_kcal)


def _thermo_for(records: Sequence[ReactionThermo], T: float | None, state: str | None
                ) -> ReactionThermo | None:
    """The record at the ranking's (T, standard state), the one summary.rank_rows ranks by;
    without a ranking condition (T None) the first: the thermo stage's first T and state."""
    return next((t for t in records if T is None or (t.T_K, t.standard_state) == (T, state)), None)


def _summary_cells(reaction: ReactionRecord, t: ReactionThermo | None,
                   row: RankRow | None) -> list[str]:
    rid, saddle = escape(reaction.reaction_id), reaction.saddle
    blockers = row.blockers if row is not None else t.blockers if t else ()
    rank = row.rank if row is not None else None
    dG_eff, dE_act, dG_act, assoc, vs_separated = map(_num, (None,) * 5 if t is None else (
        t.dG_eff_kcal, t.dE_act_kcal, t.dG_act_kcal, t.dG_assoc_kcal, t.dG_act_vs_separated_kcal))
    return [
        f"<a href='#rxn-{rid}'>{rid}</a>",
        escape(reaction.outcome.value if reaction.outcome else _DASH),
        _DASH if rank is None else str(rank), dG_eff, dE_act, dG_act,
        _imag(saddle.imag_cm1) if saddle is not None else _DASH,
        _connection(reaction), assoc, vs_separated,
        escape(", ".join(blockers)) or _DASH,
        escape(t.energy_level or _DASH) if t else _DASH,
    ]


_SUMMARY_HEADER = ("reaction", "outcome", "rank", "δG<sub>eff</sub>", "ΔE‡", "ΔG‡",
                   "ν<sub>imag</sub> (cm⁻¹)", "connection", "ΔG<sub>assoc</sub>",
                   "ΔG‡ vs separated", "blockers", "energy level")
_THERMO_HEADER = ("T (K)", "state", "ΔE‡", "ΔE<sub>rxn</sub>", "δG<sub>eff</sub>", "ΔG‡",
                  "ΔG<sub>rxn</sub>", "δG<sub>eff</sub> band", "ΔG<sub>assoc</sub>",
                  "ΔG‡ vs separated", "blockers", "notes", "energy level")


def _thermo_cells(t: ReactionThermo) -> list[str]:
    values = (t.dE_act_kcal, t.dE_rxn_kcal, t.dG_eff_kcal, t.dG_act_kcal, t.dG_rxn_kcal)
    return [f"{t.T_K:g}", t.standard_state, *map(_num, values), _band(t),
            _num(t.dG_assoc_kcal), _num(t.dG_act_vs_separated_kcal),
            escape(", ".join(t.blockers)) or _DASH, escape(", ".join(t.notes)) or _DASH,
            escape(t.energy_level or _DASH)]


def _section(reaction: ReactionRecord, thermo: Sequence[ReactionThermo],
             shown: ReactionThermo | None, xyz: str | None) -> str:
    """``thermo``: every record of the reaction, for the table; ``shown``: the diagram's."""
    rid = escape(reaction.reaction_id)
    meta = [f"outcome {escape(reaction.outcome.value if reaction.outcome else _DASH)}",
            f"source {escape(reaction.source)}", f"minima {escape(' → '.join(reaction.minima))}"]
    if reaction.reasons:
        meta.append(f"reasons {escape(', '.join(reaction.reasons))}")
    figures = []
    if shown is not None:
        caption = f"{shown.T_K:g} K, {shown.standard_state}"
        figures.append(f"<figure>{_energy_svg(_levels(shown), caption)}"
                       f"<figcaption>Energy diagram ({escape(caption)})</figcaption></figure>")
    if xyz is not None and reaction.saddle is not None:
        figures.append(f"<figure><div class='mol' data-xyz='{escape(xyz)}'></div><figcaption>"
                       f"TS imaginary mode, ν = {_imag(reaction.saddle.imag_cm1)} cm⁻¹"
                       "</figcaption></figure>")
    table = (_table(_THERMO_HEADER, map(_thermo_cells, thermo), set(range(2, 10))) if thermo
             else "<p class='meta'>No reaction thermochemistry in this view.</p>")
    return (f"<section id='rxn-{rid}'><h2>{rid}</h2><p class='meta'>{' · '.join(meta)}</p>"
            f"<div class='figs'>{''.join(figures)}</div>{table}</section>")


def render(view: Manifest, out_dir: Path, *, load_xyz: Callable[[Geometry], XYZ]) -> Path:
    """Write ``out_dir/report.html`` for the reactions with an outcome and return its path.

    Rows follow the last ReportRecord in the view (ranked first), and so does the thermo shown
    for a reaction: the ReactionThermo at the (T, standard state) of those rows (_thermo_for).
    """
    reactions = [r for r in view.records(ArtifactType.REACTION, ReactionRecord)
                 if r.outcome is not None]
    thermo: dict[str, list[ReactionThermo]] = defaultdict(list)
    for t in view.records(ArtifactType.REACTION_THERMO, ReactionThermo):
        thermo[t.reaction_id].append(t)
    reports = view.records(ArtifactType.REPORT, ReportRecord)
    rows = {row.reaction_id: row for row in reports[-1].rows} if reports else {}
    T, state = next(((r.T_K, r.standard_state) for r in rows.values()), (None, None))
    order = {rid: i for i, rid in enumerate(rows)}
    reactions.sort(key=lambda r: order.get(r.reaction_id, len(order)))
    shown = {r.reaction_id: _thermo_for(thermo[r.reaction_id], T, state) for r in reactions}
    modes = {r.reaction_id: _mode_xyz(view, r, load_xyz) for r in reactions}

    summary = [_summary_cells(r, shown[r.reaction_id], rows.get(r.reaction_id)) for r in reactions]
    run = escape(view.run_id)
    intro = (f"<h1>hfauto report</h1><p class='meta'>run {run} · view of stage "
             f"{escape(view.stage_id)} · {len(reactions)} reaction(s) · energies in kcal/mol</p>")
    body = [intro, _table(_SUMMARY_HEADER, summary, {2, 3, 4, 5, 6, 8, 9}) if reactions
            else "<p>No reactions with an outcome in this view.</p>"]
    body += [_section(r, thermo[r.reaction_id], shown[r.reaction_id], modes[r.reaction_id])
             for r in reactions]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.html"
    path.write_text(_page(run, body, with_mol=any(x is not None for x in modes.values())),
                    encoding="utf-8")
    return path


def _page(run: str, body: Sequence[str], *, with_mol: bool) -> str:
    """The document around ``body``; 3Dmol.js only when a TS mode is shown."""
    script = f"<script src='{THREEDMOL_JS}'></script>" if with_mol else ""
    head = ("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            f"content='width=device-width,initial-scale=1'><title>hfauto report · {run}</title>"
            f"<style>{_CSS}</style>{script}</head>")
    tail = f"<script>{_VIEWER_JS}</script>" if with_mol else ""
    return f"{head}<body>{''.join(body)}{tail}</body></html>\n"
