"""Generic reaction-energy profile rendering."""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_viz.core.html import write_html


def _number(row: dict[str, Any], *keys: str) -> float | None:
    for key in keys:
        value = row.get(key)
        if value in {None, ""}:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def energy_points(
    row: dict[str, Any], *, process: bool = False
) -> tuple[list[str], list[float]]:
    """Return available reactant, TS, and product Gibbs energies.

    Values are relative to the declared reactants. Missing TS or product values
    are omitted instead of being rendered as artificial zero-energy points.
    """

    activation_keys = (
        ("delta_G_activation_process_kcal_mol", "delta_G_act_process_kcal_mol")
        if process
        else (
            "delta_G_activation_standard_kcal_mol",
            "delta_G_act_kcal_mol",
        )
    )
    reaction_keys = (
        ("delta_G_reaction_process_kcal_mol",)
        if process
        else (
            "delta_G_reaction_standard_kcal_mol",
            "delta_G_reaction_kcal_mol",
        )
    )
    labels = ["Reactants"]
    values = [0.0]
    activation = _number(row, *activation_keys)
    if activation is not None:
        labels.append("TS")
        values.append(activation)
    reaction = _number(row, *reaction_keys)
    if reaction is not None:
        labels.append("Products")
        values.append(reaction)
    return labels, values


def _fallback_svg(labels: list[str], values: list[float]) -> str:
    width, height = 760, 330
    lower, upper = min(values + [0.0]), max(values + [0.0])
    if lower == upper:
        lower -= 1.0
        upper += 1.0

    def x(index: int) -> float:
        return 70 + index * (width - 140) / max(len(values) - 1, 1)

    def y(value: float) -> float:
        return height - 50 - (value - lower) * (height - 100) / (upper - lower)

    parts = [
        (
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' "
            f"height='{height}'><rect width='100%' height='100%' fill='white'/>"
        )
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        xpos, ypos = x(index), y(value)
        parts.append(
            f"<line x1='{xpos - 30}' x2='{xpos + 30}' y1='{ypos}' y2='{ypos}' "
            "stroke='#145da0' stroke-width='4'/>"
            f"<text x='{xpos}' y='{ypos - 9}' text-anchor='middle' "
            f"font-size='12'>{value:.2f}</text>"
            f"<text x='{xpos}' y='{height - 20}' text-anchor='middle' "
            f"font-size='11'>{escape(label)}</text>"
        )
        if index < len(labels) - 1:
            parts.append(
                f"<line x1='{xpos + 30}' y1='{ypos}' x2='{x(index + 1) - 30}' "
                f"y2='{y(values[index + 1])}' stroke='#999' stroke-dasharray='4 4'/>"
            )
    return "".join(parts) + "</svg>"


def plotly_energy_div(row: dict[str, Any], title: str = "Energy profile") -> str:
    standard_labels, standard_values = energy_points(row)
    process_labels, process_values = energy_points(row, process=True)
    try:
        import plotly.graph_objects as go
        import plotly.io as pio

        figure = go.Figure()
        figure.add_trace(
            go.Scatter(
                x=standard_labels,
                y=standard_values,
                mode="lines+markers+text",
                name="standard",
                text=[f"{value:.2f}" for value in standard_values],
                textposition="top center",
            )
        )
        if (process_labels, process_values) != (
            standard_labels,
            standard_values,
        ):
            figure.add_trace(
                go.Scatter(
                    x=process_labels,
                    y=process_values,
                    mode="lines+markers",
                    name="process",
                )
            )
        figure.update_layout(
            title=title,
            yaxis_title="Relative Gibbs energy / kcal mol^-1",
            template="plotly_white",
        )
        return pio.to_html(figure, include_plotlyjs="cdn", full_html=False)
    except Exception:
        return _fallback_svg(standard_labels, standard_values)


def write_energy_profile(
    row: dict[str, Any], out_path: str | Path, title: str | None = None
) -> Path:
    title = title or f"Energy profile: {row.get('reaction_id', 'reaction')}"
    fields = [
        "reaction_id",
        "reaction_type",
        "T_K",
        "quality_tier",
        "delta_G_reaction_standard_kcal_mol",
        "delta_G_activation_standard_kcal_mol",
        "delta_E_activation_kcal_mol",
        "activation_connectivity_validated",
    ]
    table = '<table class="data-table">' + "".join(
        f"<tr><th>{escape(key)}</th><td>{escape(str(row.get(key, '')))}</td></tr>"
        for key in fields
    ) + "</table>"
    body = (
        f"<div class='card'>{plotly_energy_div(row, title)}</div>"
        f"<div class='card'><h2>Key values</h2>{table}</div>"
    )
    return write_html(out_path, title, body)


def write_energy_comparison(
    dataframe: pd.DataFrame,
    out_path: str | Path,
    title: str = "Top candidate energy profiles",
    top_n: int = 12,
) -> Path:
    if dataframe.empty:
        return write_html(
            out_path, title, "<div class='warn'>No reaction table available.</div>"
        )
    rows = dataframe.head(top_n).to_dict("records")
    try:
        import plotly.graph_objects as go
        import plotly.io as pio

        figure = go.Figure()
        for row in rows:
            labels, values = energy_points(row)
            figure.add_trace(
                go.Scatter(
                    x=labels,
                    y=values,
                    mode="lines+markers",
                    name=str(row.get("reaction_id") or row.get("candidate_id")),
                )
            )
        figure.update_layout(
            title=title,
            yaxis_title="Relative Gibbs energy / kcal mol^-1",
            template="plotly_white",
        )
        chart = pio.to_html(figure, include_plotlyjs="cdn", full_html=False)
    except Exception:
        labels, values = energy_points(rows[0])
        chart = _fallback_svg(labels, values)
    body = (
        f"<div class='card'>{chart}</div>"
        f"<div class='card'>{dataframe.head(top_n).to_html(index=False, escape=True)}</div>"
    )
    return write_html(out_path, title, body)


def write_energy_profile_static_svg(
    row: dict[str, Any], out_path: str | Path, process: bool = False
) -> dict[str, Path]:
    """Write a static SVG and its exact plotted values as CSV."""

    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    labels, values = energy_points(row, process=process)
    path.write_text(_fallback_svg(labels, values), encoding="utf-8")
    csv_path = path.with_suffix(".csv")
    pd.DataFrame(
        {"state": labels, "relative_G_kcal_mol": values}
    ).to_csv(csv_path, index=False)
    return {"energy_profile_svg": path, "energy_profile_csv": csv_path}
