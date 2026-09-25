"""Generate the reaction-case figures expected by the hfauto dossier design.

This is a stopped-run adapter.  It uses the reactions and atom mappings from
the preoptimization manifest, then reads any available NWChem minima output.
Missing TS/IRC data are shown as missing; no synthetic barrier or path is
created.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from report_partial_nwchem_run import (
    HARTREE_TO_KJ_MOL,
    SHORT_NAME,
    STATUS_COLOR,
    calculation_rows,
    configure_plotting,
    endpoint_rows,
    load_species,
    markdown_table,
    read_xyz,
    save_figure,
)

from hfauto.backends.qm.nwchem import parse_nwchem_output
from hfauto_viz.core.html import write_html
from hfauto_viz.renderers.mol3d_3dmol import write_animation_html, write_structure_html
from hfauto_viz.structure.bond_changes import annotation_shapes_for_frame
from hfauto_viz.structure.xyz import (
    Frame,
    frames_to_xyz,
    interpolate_frames,
    read_xyz_frames,
)


def project_path(project_root: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    if path.is_absolute():
        return path
    return project_root / path


def load_case_inputs(
    manifest_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, Path], dict[str, str]]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    project_root = manifest_path.resolve().parents[3]
    molecule_names = {
        str(artifact["data"]["mol_id"]): str(artifact["data"].get("name"))
        for artifact in payload["artifacts"]
        if artifact.get("artifact_type") == "molecule"
    }
    preopt_xyz: dict[str, Path] = {}
    for artifact in payload["artifacts"]:
        if artifact.get("artifact_type") != "species_preopt":
            continue
        data = artifact.get("data", {})
        species_id = str(
            data.get("species_id")
            or data.get("source_species_id")
            or str(artifact["artifact_id"]).removeprefix("preopt_")
        )
        raw = data.get("xyz_path") or artifact.get("paths", {}).get("xyz")
        path = project_path(project_root, str(raw) if raw else None)
        if path is not None:
            preopt_xyz[species_id] = path
    reactions = [
        dict(artifact.get("data", {}))
        for artifact in payload["artifacts"]
        if artifact.get("artifact_type") == "reaction"
    ]
    reactions.sort(key=lambda row: (str(row.get("mol_id")), int(row.get("hf_n", 0))))
    return reactions, preopt_xyz, molecule_names


def case_label(reaction: dict[str, Any], molecule_names: dict[str, str]) -> str:
    molecule = molecule_names[str(reaction["mol_id"])]
    short = SHORT_NAME.get(molecule, molecule)
    hf_n = int(reaction["hf_n"])
    return f"{short}-{'HF' if hf_n == 1 else f'(HF){hf_n}'}"


def final_or_latest(minima_dir: Path, species_id: str) -> tuple[Path | None, str]:
    workdir = minima_dir / species_id
    final = workdir / "final.xyz"
    if final.exists():
        output = workdir / "nwchem.out"
        completed = False
        if output.exists():
            parsed = parse_nwchem_output(output.read_text(encoding="utf-8", errors="ignore"))
            completed = bool(parsed.get("normal_termination") and parsed.get("geometry_converged"))
        return final, "NWChem final" if completed else "NWChem interrupted"
    candidates = sorted(workdir.glob("final-*.xyz"), key=lambda path: path.stat().st_mtime)
    if candidates:
        return candidates[-1], "NWChem interrupted"
    return None, "unavailable"


def first_frame(path: Path | None) -> Frame | None:
    if path is None:
        return None
    frames = read_xyz_frames(path)
    return frames[0] if frames else None


def reaction_annotations(frame: Frame | None, reaction: dict[str, Any]) -> list[dict[str, Any]]:
    return annotation_shapes_for_frame(
        frame,
        reaction,
        include_labels=True,
        include_bond_lines=True,
    )


def framed_annotations(
    frames: list[tuple[str, Frame]], reaction: dict[str, Any]
) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    for index, (_, frame) in enumerate(frames):
        annotations.extend(
            {**annotation, "frame": index}
            for annotation in reaction_annotations(frame, reaction)
        )
    return annotations


def reaction_coordinates(path: Path, atoms: dict[str, Any]) -> dict[str, float]:
    _, coords = read_xyz(path)
    base = int(atoms["base_atom"])
    transfer_h = int(atoms["transfer_h"])
    leaving_f = int(atoms["leaving_f"])
    r_bh = float(np.linalg.norm(coords[base] - coords[transfer_h]))
    r_hf = float(np.linalg.norm(coords[transfer_h] - coords[leaving_f]))
    return {"r_BH_A": r_bh, "r_HF_A": r_hf, "q_BH_minus_HF_A": r_bh - r_hf}


def optimization_trace(output_path: Path) -> pd.DataFrame:
    if not output_path.exists():
        return pd.DataFrame()
    text = output_path.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(
        r"^@\s+(\d+)\s+(-?\d+\.\d+)\s+\S+\s+([0-9.]+)\s+([0-9.]+)"
        r"\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$",
        flags=re.MULTILINE,
    )
    rows: dict[int, dict[str, float | int]] = {}
    for match in pattern.finditer(text):
        step = int(match.group(1))
        rows[step] = {
            "step": step,
            "energy_hartree": float(match.group(2)),
            "gmax": float(match.group(3)),
            "grms": float(match.group(4)),
            "xrms": float(match.group(5)),
            "xmax": float(match.group(6)),
            "wall_s": float(match.group(7)),
        }
    frame = pd.DataFrame(list(rows.values())).sort_values("step") if rows else pd.DataFrame()
    if not frame.empty:
        frame["relative_energy_kj_mol"] = (
            frame["energy_hartree"] - frame["energy_hartree"].min()
        ) * HARTREE_TO_KJ_MOL
    return frame


def parsed_output(minima_dir: Path, species_id: str) -> dict[str, Any] | None:
    path = minima_dir / species_id / "nwchem.out"
    if not path.exists():
        return None
    return parse_nwchem_output(path.read_text(encoding="utf-8", errors="ignore"))


def energy_for(calculations: pd.DataFrame, species_id: str) -> float | None:
    selected = calculations[calculations["species_id"] == species_id]
    if len(selected) != 1 or not bool(selected.iloc[0]["normal_termination"]):
        return None
    value = selected.iloc[0]["electronic_energy_hartree"]
    return None if pd.isna(value) else float(value)


def association_energy(
    calculations: pd.DataFrame,
    reaction: dict[str, Any],
    endpoint_species_id: str,
) -> float | None:
    endpoint = energy_for(calculations, endpoint_species_id)
    candidate = energy_for(calculations, str(reaction["candidate_species_id"]))
    hf = energy_for(calculations, "spc_hf_cluster_hf1")
    if endpoint is None or candidate is None or hf is None:
        return None
    return (endpoint - candidate - int(reaction["hf_n"]) * hf) * HARTREE_TO_KJ_MOL


def endpoint_status(calculations: pd.DataFrame, species_id: str) -> str:
    selected = calculations[calculations["species_id"] == species_id]
    return str(selected.iloc[0]["status"]) if len(selected) == 1 else "not_started"


def plot_energy_profile(
    label: str,
    reaction: dict[str, Any],
    calculations: pd.DataFrame,
    output_path: Path,
) -> tuple[float | None, float | None]:
    rc_id = str(reaction["reactant_species_id"])
    ip_id = str(reaction["product_species_id"])
    rc_energy = association_energy(calculations, reaction, rc_id)
    ip_energy = association_energy(calculations, reaction, ip_id)
    values = [0.0 if energy_for(calculations, str(reaction["candidate_species_id"])) is not None else None,
              rc_energy, None, ip_energy]
    states = ["A + nHF", "RC", "TS", "IP"]
    colors = ["#4C78A8", STATUS_COLOR[endpoint_status(calculations, rc_id)],
              "#A7ADB4", STATUS_COLOR[endpoint_status(calculations, ip_id)]]
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    for index, (state, value, color) in enumerate(zip(states, values, colors)):
        if value is None:
            ax.scatter(index, 0, marker="x", s=70, color="#A7ADB4")
            ax.text(index, 2, "not computed", ha="center", va="bottom", color="#666666", fontsize=8)
            continue
        ax.hlines(value, index - 0.28, index + 0.28, color=color, linewidth=5)
        ax.text(index, value + 2, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
    finite = [value for value in values if value is not None]
    low = min(finite + [-10.0]) - 12
    high = max(finite + [0.0]) + 12
    ax.set_ylim(low, high)
    ax.set_xticks(range(4), states)
    ax.set_ylabel("Relative electronic energy / kJ mol$^{-1}$")
    ax.set_title(f"{label}: available-state energy profile")
    ax.text(
        0.99,
        0.02,
        "No TS/IRC calculation; levels are not a connected path",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#8C3B2A",
    )
    save_figure(fig, output_path)
    return rc_energy, ip_energy


def coordinate_records(
    reaction: dict[str, Any],
    preopt_xyz: dict[str, Path],
    minima_dir: Path,
) -> pd.DataFrame:
    rc_id = str(reaction["reactant_species_id"])
    ip_id = str(reaction["product_species_id"])
    atoms = dict(reaction["reaction_coordinate"]["atoms"])
    records: list[dict[str, Any]] = []
    candidates: list[tuple[str, Path | None, str]] = [
        ("RC preopt", preopt_xyz.get(rc_id), "preoptimized"),
    ]
    rc_path, rc_source = final_or_latest(minima_dir, rc_id)
    ip_path, ip_source = final_or_latest(minima_dir, ip_id)
    candidates.extend(
        [
            ("RC NWChem", rc_path, rc_source),
            ("IP NWChem", ip_path, ip_source),
            ("IP preopt", preopt_xyz.get(ip_id), "preoptimized"),
        ]
    )
    for stage, path, source in candidates:
        if path is None or not path.exists():
            continue
        records.append({"stage": stage, "source": source, **reaction_coordinates(path, atoms)})
    return pd.DataFrame(records)


def plot_coordinates(label: str, coordinates: pd.DataFrame, output_path: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.3))
    x = np.arange(len(coordinates))
    width = 0.34
    left.bar(x - width / 2, coordinates["r_BH_A"], width, label="r(B-H)", color="#F2A541")
    left.bar(x + width / 2, coordinates["r_HF_A"], width, label="r(H-F)", color="#4C78A8")
    left.set_xticks(x, coordinates["stage"], rotation=40, ha="right")
    left.set_ylabel("Distance / A")
    left.set_title("Mapped proton-contact distances")
    left.legend()
    colors = ["#C7473A" if "interrupted" in value else "#59A14F" for value in coordinates["source"]]
    right.bar(coordinates["stage"], coordinates["q_BH_minus_HF_A"], color=colors)
    right.axhline(0, color="black", linewidth=0.8)
    right.set_ylabel("q = r(B-H) - r(H-F) / A")
    right.set_title("q > 0: H remains closer to F")
    right.tick_params(axis="x", rotation=40)
    fig.suptitle(f"{label}: endpoint reaction-coordinate diagnostics\n(not an IRC path)")
    save_figure(fig, output_path)


def inferred_bonds(frame: Frame) -> list[tuple[int, int]]:
    radii = {"H": 0.31, "C": 0.76, "N": 0.71, "F": 0.57, "O": 0.66}
    bonds: list[tuple[int, int]] = []
    for left in range(len(frame)):
        for right in range(left + 1, len(frame)):
            if frame[left].symbol == frame[right].symbol == "H":
                continue
            distance = np.linalg.norm(
                np.array([frame[left].x, frame[left].y, frame[left].z])
                - np.array([frame[right].x, frame[right].y, frame[right].z])
            )
            cutoff = 1.25 * (
                radii.get(frame[left].symbol, 0.75)
                + radii.get(frame[right].symbol, 0.75)
            )
            if 0.35 < distance <= cutoff:
                bonds.append((left, right))
    return bonds


def draw_frame_3d(
    axis: Any,
    frame: Frame | None,
    title: str,
    atoms: dict[str, Any],
) -> None:
    axis.set_facecolor("white")
    for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
        pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
    if not frame:
        axis.text2D(0.5, 0.5, "Structure unavailable", transform=axis.transAxes,
                    ha="center", va="center", color="#666666")
        axis.set_title(title)
        axis.set_axis_off()
        return
    colors = {"H": "#E8E8E8", "C": "#555555", "N": "#3569B0", "F": "#55A868", "O": "#C7473A"}
    sizes = {"H": 45, "C": 105, "N": 125, "F": 120, "O": 120}
    highlight = {
        int(atoms["base_atom"]): "B",
        int(atoms["transfer_h"]): "H*",
        int(atoms["leaving_f"]): "F*",
    }
    for left, right in inferred_bonds(frame):
        axis.plot(
            [frame[left].x, frame[right].x],
            [frame[left].y, frame[right].y],
            [frame[left].z, frame[right].z],
            color="#888888",
            linewidth=1.5,
            zorder=1,
        )
    for index, atom in enumerate(frame):
        axis.scatter(
            atom.x,
            atom.y,
            atom.z,
            s=sizes.get(atom.symbol, 90),
            color=colors.get(atom.symbol, "#999999"),
            edgecolor="#C7473A" if index in highlight else "#333333",
            linewidth=2.0 if index in highlight else 0.5,
            depthshade=True,
            zorder=2,
        )
        if index in highlight:
            axis.text(atom.x, atom.y, atom.z, f" {highlight[index]}", color="#8C2D23", fontsize=9)
    xyz = np.array([[atom.x, atom.y, atom.z] for atom in frame])
    center = xyz.mean(axis=0)
    radius = max(float(np.ptp(xyz, axis=0).max()) / 2, 1.0)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))
    axis.view_init(elev=22, azim=-58)
    axis.set_title(title)
    axis.set_axis_off()


def write_case_3d_outputs(
    label: str,
    reaction: dict[str, Any],
    preopt_xyz: dict[str, Path],
    minima_dir: Path,
    case_dir: Path,
) -> dict[str, Path]:
    rc_id = str(reaction["reactant_species_id"])
    ip_id = str(reaction["product_species_id"])
    rc_path, rc_source = final_or_latest(minima_dir, rc_id)
    ip_path, ip_source = final_or_latest(minima_dir, ip_id)
    rc_display = rc_path or preopt_xyz.get(rc_id)
    ip_display = ip_path or preopt_xyz.get(ip_id)
    rc_frame = first_frame(rc_display)
    ip_frame = first_frame(ip_display)
    atoms = dict(reaction["reaction_coordinate"]["atoms"])
    fig = plt.figure(figsize=(9.5, 4.6))
    left = fig.add_subplot(1, 2, 1, projection="3d")
    right = fig.add_subplot(1, 2, 2, projection="3d")
    draw_frame_3d(left, rc_frame, f"RC: {rc_source if rc_path else 'preoptimized'}", atoms)
    draw_frame_3d(right, ip_frame, f"IP: {ip_source if ip_path else 'preoptimized'}", atoms)
    fig.suptitle(f"{label}: available endpoint structures")
    save_figure(fig, case_dir / "endpoint_structures_3d.png")
    outputs = {
        "static_3d": case_dir / "endpoint_structures_3d.png",
        "rc_3d": write_structure_html(
            rc_frame,
            case_dir / "rc_3d.html",
            f"{label}: RC 3D ({rc_source if rc_path else 'preoptimized'})",
            library_mode="cdn",
            annotations=reaction_annotations(rc_frame, reaction),
        ),
        "ip_3d": write_structure_html(
            ip_frame,
            case_dir / "ip_3d.html",
            f"{label}: IP 3D ({ip_source if ip_path else 'preoptimized'})",
            library_mode="cdn",
            annotations=reaction_annotations(ip_frame, reaction),
        ),
    }
    comparison_frames: list[tuple[str, Frame]] = []
    for name, path in [
        ("RC preopt", preopt_xyz.get(rc_id)),
        ("RC NWChem", rc_path),
        ("IP NWChem", ip_path),
        ("IP preopt", preopt_xyz.get(ip_id)),
    ]:
        frame = first_frame(path)
        if frame:
            comparison_frames.append((name, frame))
    outputs["endpoint_animation"] = write_animation_html(
        comparison_frames,
        case_dir / "endpoint_comparison_3d.html",
        title=f"{label}: endpoint comparison (not a reaction path)",
        library_mode="cdn",
        annotations=framed_annotations(comparison_frames, reaction),
    )
    preopt_rc = first_frame(preopt_xyz.get(rc_id))
    preopt_ip = first_frame(preopt_xyz.get(ip_id))
    illustrative_frames = interpolate_frames(preopt_rc or [], preopt_ip or [], n=21, prefix="illustrative")
    animation_path = write_animation_html(
        illustrative_frames,
        case_dir / "reaction_animation_illustrative.html",
        title=f"{label}: illustrative endpoint interpolation",
        library_mode="cdn",
        annotations=framed_annotations(illustrative_frames, reaction),
    )
    warning = (
        "<div class='warn'><b>Scientific warning:</b> This is a linear interpolation "
        "between preoptimized RC/IP inputs. It is not an NEB trajectory, transition-state "
        "path, or IRC and must not be interpreted as molecular dynamics.</div>"
    )
    html = animation_path.read_text(encoding="utf-8")
    animation_path.write_text(html.replace("<main>", "<main>" + warning, 1), encoding="utf-8")
    outputs["reaction_animation"] = animation_path
    (case_dir / "reaction_animation_illustrative.xyz").write_text(
        frames_to_xyz(illustrative_frames), encoding="utf-8"
    )
    (case_dir / "path_metadata.json").write_text(
        json.dumps(
            {
                "source": "linear_interpolation_between_preoptimized_endpoints",
                "real_path_frames": False,
                "n_frames": len(illustrative_frames),
                "warning": "Not NEB, TS, IRC, or molecular dynamics.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return outputs


def plot_optimization(
    label: str,
    reaction: dict[str, Any],
    minima_dir: Path,
    calculations: pd.DataFrame,
    output_path: Path,
) -> dict[str, pd.DataFrame]:
    traces: dict[str, pd.DataFrame] = {}
    for endpoint, species_key in [("RC", "reactant_species_id"), ("IP", "product_species_id")]:
        species_id = str(reaction[species_key])
        trace = optimization_trace(minima_dir / species_id / "nwchem.out")
        if not trace.empty:
            traces[endpoint] = trace
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.2))
    endpoint_color = {"RC": "#4C78A8", "IP": "#F2A541"}
    for endpoint, trace in traces.items():
        species_id = str(reaction["reactant_species_id" if endpoint == "RC" else "product_species_id"])
        status = endpoint_status(calculations, species_id)
        line_style = "--" if status == "interrupted" else "-"
        left.plot(trace["step"], trace["relative_energy_kj_mol"], marker="o", markersize=3,
                  color=endpoint_color[endpoint], linestyle=line_style, label=f"{endpoint} ({status})")
        right.plot(trace["step"], trace["grms"], marker="o", markersize=3,
                   color=endpoint_color[endpoint], linestyle=line_style, label=endpoint)
    if not traces:
        for axis in (left, right):
            axis.text(0.5, 0.5, "No NWChem optimization trace", transform=axis.transAxes,
                      ha="center", va="center", color="#666666")
    left.set_xlabel("Optimization step")
    left.set_ylabel("E - lowest sampled E / kJ mol$^{-1}$")
    left.set_title("Electronic-energy relaxation")
    right.set_xlabel("Optimization step")
    right.set_ylabel("Gradient RMS / a.u.")
    right.set_yscale("log")
    right.set_title("Gradient convergence")
    if traces:
        left.legend()
        right.legend()
    fig.suptitle(f"{label}: endpoint optimization histories")
    save_figure(fig, output_path)
    return traces


def plot_frequencies(
    label: str,
    reaction: dict[str, Any],
    minima_dir: Path,
    output_path: Path,
) -> dict[str, list[float]]:
    spectra: dict[str, list[float]] = {}
    for endpoint, species_key in [("RC", "reactant_species_id"), ("IP", "product_species_id")]:
        parsed = parsed_output(minima_dir, str(reaction[species_key]))
        if parsed and parsed.get("normal_termination"):
            spectra[endpoint] = [float(value) for value in parsed.get("frequencies_cm1", [])]
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    markers = {"RC": "o", "IP": "x"}
    colors = {"RC": "#4C78A8", "IP": "#F2A541"}
    for endpoint, frequencies in spectra.items():
        modes = np.arange(1, len(frequencies) + 1)
        ax.scatter(modes, frequencies, marker=markers[endpoint], color=colors[endpoint],
                   s=24, alpha=0.8, label=endpoint)
        negative = [(index + 1, value) for index, value in enumerate(frequencies) if value < 0]
        if negative:
            ax.scatter([row[0] for row in negative], [row[1] for row in negative],
                       marker="o", facecolors="none", edgecolors="#C7473A", s=75, linewidths=1.5)
    if not spectra:
        ax.text(0.5, 0.5, "No completed frequency calculation", transform=ax.transAxes,
                ha="center", va="center", color="#666666")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Vibrational mode index")
    ax.set_ylabel("Frequency / cm$^{-1}$")
    ax.set_title(f"{label}: endpoint vibrational spectra")
    if spectra:
        ax.legend()
    save_figure(fig, output_path)
    return spectra


def status_summary(calculations: pd.DataFrame, species_id: str) -> tuple[str, str]:
    selected = calculations[calculations["species_id"] == species_id]
    if len(selected) != 1:
        return "not started", "—"
    row = selected.iloc[0]
    n_imag = "—" if pd.isna(row["n_imag"]) else str(int(row["n_imag"]))
    return str(row["status"]), n_imag


def network_node_color(status: str) -> str:
    return STATUS_COLOR.get(status, "#A7ADB4")


def write_reaction_network(
    output_dir: Path,
    reactions: list[dict[str, Any]],
    molecule_names: dict[str, str],
    calculations: pd.DataFrame,
) -> dict[str, Path]:
    fig, axis = plt.subplots(figsize=(12.0, 8.2))
    x_positions = {"separated": 0.0, "RC": 1.8, "TS": 3.6, "IP": 5.4}
    dot_lines = [
        "digraph stopped_hfauto_network {",
        "rankdir=LR;",
        'graph [bgcolor="white", nodesep=0.35, ranksep=0.65];',
        'node [shape=box, style="rounded,filled", fontname="Helvetica"];',
    ]
    for row_index, reaction in enumerate(reactions):
        label = case_label(reaction, molecule_names)
        y = len(reactions) - row_index - 1
        rc_id = str(reaction["reactant_species_id"])
        ip_id = str(reaction["product_species_id"])
        rc_status = endpoint_status(calculations, rc_id)
        ip_status = endpoint_status(calculations, ip_id)
        candidate_ok = energy_for(calculations, str(reaction["candidate_species_id"])) is not None
        nodes = {
            "separated": ("A + nHF", "minimum" if candidate_ok else "not_started"),
            "RC": ("RC", rc_status),
            "TS": ("TS\nnot computed", "not_started"),
            "IP": ("IP", ip_status),
        }
        for stage, (node_label, status) in nodes.items():
            x = x_positions[stage]
            box = FancyBboxPatch(
                (x - 0.46, y - 0.27),
                0.92,
                0.54,
                boxstyle="round,pad=0.04",
                facecolor=network_node_color(status),
                edgecolor="#555555",
                linewidth=1.0,
            )
            axis.add_patch(box)
            axis.text(x, y, node_label, ha="center", va="center", color="white", fontsize=8)
            node_id = f"case_{row_index}_{stage}"
            dot_lines.append(
                f'{node_id} [label="{label}\\n{node_label}", '
                f'fillcolor="{network_node_color(status)}"];'
            )
        for start, end, real in [
            ("separated", "RC", rc_status in {"minimum", "imaginary_modes"}),
            ("RC", "TS", False),
            ("TS", "IP", False),
        ]:
            arrow = FancyArrowPatch(
                (x_positions[start] + 0.47, y),
                (x_positions[end] - 0.47, y),
                arrowstyle="-|>",
                mutation_scale=11,
                linewidth=1.4,
                linestyle="-" if real else "--",
                color="#555555" if real else "#A7ADB4",
            )
            axis.add_patch(arrow)
            style = "solid" if real else "dashed"
            dot_lines.append(
                f"case_{row_index}_{start} -> case_{row_index}_{end} "
                f'[style="{style}", color="{("#555555" if real else "#A7ADB4")}"];'
            )
        axis.text(-0.65, y, label, ha="right", va="center", fontsize=9, weight="bold")
    dot_lines.append("}")
    axis.set_xlim(-2.0, 6.1)
    axis.set_ylim(-0.7, len(reactions) - 0.3)
    axis.set_xticks(list(x_positions.values()), ["Separated", "RC", "TS", "IP"])
    axis.set_yticks([])
    axis.set_title("Amine-HF reaction-case network at manual stop")
    axis.legend(
        handles=[
            Patch(color=STATUS_COLOR["minimum"], label="true minimum / available reference"),
            Patch(color=STATUS_COLOR["imaginary_modes"], label="imaginary mode(s)"),
            Patch(color=STATUS_COLOR["interrupted"], label="interrupted"),
            Patch(color=STATUS_COLOR["not_started"], label="not computed"),
        ],
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.13),
    )
    axis.spines[:].set_visible(False)
    network_png = output_dir / "reaction_network.png"
    svg_path = output_dir / "reaction_network.svg"
    fig.tight_layout()
    fig.savefig(network_png, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    dot_path = output_dir / "reaction_network.dot"
    dot_path.write_text("\n".join(dot_lines), encoding="utf-8")
    dot_executable = shutil.which("dot")
    if dot_executable:
        subprocess.run(
            [dot_executable, "-Tsvg", str(dot_path), "-o", str(svg_path)],
            check=True,
            timeout=30,
        )
    network_html = write_html(
        output_dir / "reaction_network.html",
        "Amine-HF reaction network",
        "<div class='warn'>Dashed grey edges and nodes are not computed. TS/IRC data are absent.</div>"
        "<div class='card'><img src='reaction_network.png' style='max-width:100%'></div>"
        "<div class='card'><a href='reaction_network.dot'>Graphviz DOT</a> | "
        "<a href='reaction_network.svg'>Graphviz SVG</a></div>",
    )
    return {"png": network_png, "dot": dot_path, "svg": svg_path, "html": network_html}


def write_case_report(
    output_dir: Path,
    reactions: list[dict[str, Any]],
    molecule_names: dict[str, str],
    calculations: pd.DataFrame,
    endpoints: pd.DataFrame,
    preopt_xyz: dict[str, Path],
    minima_dir: Path,
) -> None:
    sections: list[str] = []
    index_rows: list[dict[str, Any]] = []
    for position, reaction in enumerate(reactions, start=1):
        label = case_label(reaction, molecule_names)
        slug = f"{position:02d}_{label.replace('(', '').replace(')', '').replace('-', '_')}"
        case_dir = output_dir / "case_reports" / slug
        case_dir.mkdir(parents=True, exist_ok=True)
        rc_energy, ip_energy = plot_energy_profile(
            label, reaction, calculations, case_dir / "energy_profile.png"
        )
        coordinates = coordinate_records(reaction, preopt_xyz, minima_dir)
        coordinates.to_csv(case_dir / "reaction_coordinates.csv", index=False)
        plot_coordinates(label, coordinates, case_dir / "reaction_coordinates.png")
        write_case_3d_outputs(label, reaction, preopt_xyz, minima_dir, case_dir)
        traces = plot_optimization(
            label, reaction, minima_dir, calculations, case_dir / "optimization_trace.png"
        )
        for endpoint, trace in traces.items():
            trace.to_csv(case_dir / f"optimization_trace_{endpoint.lower()}.csv", index=False)
        spectra = plot_frequencies(label, reaction, minima_dir, case_dir / "frequency_spectrum.png")
        for endpoint, frequencies in spectra.items():
            pd.DataFrame({"mode": range(1, len(frequencies) + 1), "frequency_cm1": frequencies}).to_csv(
                case_dir / f"frequencies_{endpoint.lower()}.csv", index=False
            )
        rc_id = str(reaction["reactant_species_id"])
        ip_id = str(reaction["product_species_id"])
        rc_status, rc_imag = status_summary(calculations, rc_id)
        ip_status, ip_imag = status_summary(calculations, ip_id)
        molecule = molecule_names[str(reaction["mol_id"])]
        match = endpoints[
            (endpoints["molecule"] == molecule) & (endpoints["hf_n"] == int(reaction["hf_n"]))
        ]
        collapse = "not evaluable"
        endpoint_note = "RC/IPの両方が正常終了していないため、独立性を判定できない。"
        if len(match) == 1:
            row = match.iloc[0]
            collapse = "same endpoint" if bool(row["same_endpoint"]) else "distinct"
            endpoint_note = (
                f"RC/IP RMSD={row['rmsd_A']:.6g} Å、|ΔE|="
                f"{abs(row['delta_e_ip_minus_rc_kj_mol']):.6g} kJ mol⁻¹で、"
                "同一端点へ収束した。"
            )
        optimized = coordinates[coordinates["stage"].isin(["RC NWChem", "IP NWChem"])]
        q_note = "NWChem端点座標なし。"
        if not optimized.empty:
            q_values = ", ".join(
                f"{row.stage}: q={row.q_BH_minus_HF_A:.3f} Å" for row in optimized.itertuples()
            )
            q_note = f"{q_values}。q>0なのでHはF側に近い。"
        energy_note = "正常終了した端点エネルギーなし。"
        if rc_energy is not None:
            energy_note = f"RC電子会合エネルギーは{rc_energy:.1f} kJ mol⁻¹。"
        if ip_energy is not None:
            energy_note += f" IPは{ip_energy:.1f} kJ mol⁻¹。"
        quality_note = (
            f"RC={rc_status}（虚振動{rc_imag}）、IP={ip_status}（虚振動{ip_imag}）。"
        )
        conclusion = (
            "TS/IRCが未計算であるため、障壁と速度定数は評価不能。"
            if collapse == "not evaluable"
            else "RC/IPが同一端点なので、現端点を用いたTS/IRC探索は不適切。"
        )
        pd.DataFrame(
            [
                {
                    "case": label,
                    "reaction_id": reaction["reaction_id"],
                    "rc_status": rc_status,
                    "rc_n_imag": rc_imag,
                    "ip_status": ip_status,
                    "ip_n_imag": ip_imag,
                    "rc_association_energy_kj_mol": rc_energy,
                    "ip_association_energy_kj_mol": ip_energy,
                    "endpoint_classification": collapse,
                }
            ]
        ).to_csv(case_dir / "case_summary.csv", index=False)
        relative = f"case_reports/{slug}"
        sections.append(
            f"""## {position}. {label}

### 利用可能状態のエネルギープロファイル

![{label} energy]({relative}/energy_profile.png)

{energy_note} TSは未計算なので、RC–TS–IPを線で接続していない。

### 端点の反応座標

![{label} coordinates]({relative}/reaction_coordinates.png)

{q_note} この図はRC→TS→IPの実経路ではなく、preoptとNWChem端点の診断である。

### 3D分子構造とアニメーション

![{label} 3D endpoints]({relative}/endpoint_structures_3d.png)

- [RCインタラクティブ3D]({relative}/rc_3d.html)
- [IPインタラクティブ3D]({relative}/ip_3d.html)
- [preopt/NWChem端点比較アニメーション]({relative}/endpoint_comparison_3d.html)
- [説明用RC→IP補間アニメーション]({relative}/reaction_animation_illustrative.html)
- [補間XYZ]({relative}/reaction_animation_illustrative.xyz)
- [経路メタデータ]({relative}/path_metadata.json)

説明用アニメーションはpreopt RC/IPの線形補間であり、NEB、TS、IRC、分子動力学ではない。実反応運動の証拠には使用できない。

### 最適化履歴

![{label} optimization]({relative}/optimization_trace.png)

RCは{len(traces.get('RC', []))}ステップ、IPは{len(traces.get('IP', []))}ステップの履歴を取得した。中断軌跡は破線で示す。

### 振動スペクトル

![{label} frequencies]({relative}/frequency_spectrum.png)

{quality_note} 赤い中空円は虚振動を示す。

### ケース判定

{endpoint_note} {conclusion}
"""
        )
        index_rows.append(
            {
                "case": label,
                "reaction_id": reaction["reaction_id"],
                "RC": rc_status,
                "IP": ip_status,
                "endpoint": collapse,
            }
        )
    index = pd.DataFrame(index_rows)
    index.to_csv(output_dir / "case_index.csv", index=False)
    write_reaction_network(output_dir, reactions, molecule_names, calculations)
    table = markdown_table(
        index,
        ["case", "reaction_id", "RC", "IP", "endpoint"],
        ["case", "reaction_id", "RC", "IP", "endpoint"],
    )
    report = f"""# アミン–HF反応ケース別グラフ・評価

## このレポートの位置づけ

本レポートは、元のhfauto反応ドシエが想定するreaction ID単位で図を整理したものである。対象は4分子×HF₁/HF₂の8ケースである。

本runではTS、NEB、IRC、熱化学ステージまで到達していない。そのため、存在しないTS障壁やIRC経路を補間・捏造せず、各図では未計算として表示した。エネルギー図は利用可能な電子エネルギー準位、反応座標図はpreoptと最適化端点の比較であり、実反応経路ではない。

## 正しく計算できているか

| 判定層 | 結論 | 根拠 |
|---|---|---|
| NWChemプログラム実行 | 正常 | 15構造で正常終了、SCF・幾何最適化・振動解析出力を取得 |
| 極小点探索 | 部分成功 | 15正常終了中、虚振動0は10構造、5構造は虚振動あり |
| RC/IP端点構築 | 不成立 | 正常終了した5対がRMSD・電子エネルギーとも同一端点へ収束 |
| TS・反応経路 | 未計算 | NEB/TS/IRCステージ未実行 |
| 自由エネルギー・速度定数 | 未成立 | TS/IRC、単一点補正、quasi-RRHO熱化学が未完了 |

したがって、電子状態計算プログラムは正しく動いているが、目的としたプロトン移動反応を正しく完了できたとは言えない。

## 反応ネットワーク

![Reaction network](reaction_network.png)

- [ネットワークHTML](reaction_network.html)
- [Graphviz SVG](reaction_network.svg)
- [Graphviz DOT](reaction_network.dot)

緑は真の最小点、橙は虚振動あり、赤は中断、灰色は未計算を示す。破線は実計算経路が存在しない接続である。

## ケース一覧

{table}

## 全ケースに共通する結論

- 正常終了した5つのRC/IP対はすべて同一端点へ収束した。
- 得られたNWChem端点ではq=r(B–H)−r(H–F)>0で、移動Hは窒素よりフッ素に近い。
- TMA–(HF)₂とアニリン–HFは虚振動を持ち、安定最小点として採用できない。
- アニリン–(HF)₂はRCが中断、IPは開始直後に停止した。
- ピリジンはpreopt端点だけ存在し、NWChemエネルギー・最適化・振動解析はない。
- 独立RC/IP、TS、IRCがないため、反応障壁・速度定数・本来の4点自由エネルギープロファイルは未確定である。

{"".join(sections)}

## 次に必要な計算

1. 拘束付きN–H/H–F緩和スキャンから独立RC/IPを構築する。
2. 虚振動構造を負モード±方向へ変位して再最適化する。
3. 端点独立性ゲート通過後にのみNEB/TSを実行する。
4. TSの虚振動モードとIRC両端接続を検証してから自由エネルギープロファイルを完成させる。
"""
    (output_dir / "CASE_REPORT.md").write_text(report, encoding="utf-8")
    summary_path = output_dir / "REPORT.md"
    if summary_path.exists():
        summary = summary_path.read_text(encoding="utf-8")
        case_link = (
            "\n> **ケース別グラフ:** "
            "[8反応ケースのエネルギー・反応座標・最適化・振動解析](CASE_REPORT.md)\n"
        )
        if case_link.strip() not in summary:
            first_break = summary.find("\n")
            summary = summary[: first_break + 1] + case_link + summary[first_break + 1 :]
            summary_path.write_text(summary, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    output_dir = (args.output_dir or run_dir / "analysis_report").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / "05_preopt" / "manifest.json"
    species, _ = load_species(manifest_path)
    calculations = calculation_rows(species, run_dir / "06_dft-minima")
    endpoints = endpoint_rows(calculations, run_dir / "06_dft-minima")
    reactions, preopt_xyz, molecule_names = load_case_inputs(manifest_path)
    configure_plotting()
    write_case_report(
        output_dir,
        reactions,
        molecule_names,
        calculations,
        endpoints,
        preopt_xyz,
        run_dir / "06_dft-minima",
    )
    print(f"wrote {len(reactions)} case dossiers to {output_dir / 'case_reports'}")


if __name__ == "__main__":
    main()
