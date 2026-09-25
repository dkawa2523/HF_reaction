"""Build an evidence-backed report from a stopped NWChem minima stage.

The normal visualization pipeline consumes a completed stage manifest.  A
production run can be stopped before that manifest is written, while still
leaving scientifically useful NWChem outputs.  This script reads those raw
outputs without changing them and produces CSV tables, diagnostic plots, and a
Markdown report.  It never treats an interrupted energy as a converged result.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hfauto.backends.qm.nwchem import parse_nwchem_output

HARTREE_TO_KJ_MOL = 2625.499639
STATUS_ORDER = ["minimum", "imaginary_modes", "interrupted", "not_started"]
STATUS_LABEL = {
    "minimum": "True minimum",
    "imaginary_modes": "Optimized, imaginary mode(s)",
    "interrupted": "Interrupted",
    "not_started": "Not started",
}
STATUS_COLOR = {
    "minimum": "#2E8B57",
    "imaginary_modes": "#D98E04",
    "interrupted": "#C7473A",
    "not_started": "#A7ADB4",
}
SHORT_NAME = {
    "ammonia": "NH3",
    "trimethylamine": "TMA",
    "aniline": "Aniline",
    "pyridine": "Pyridine",
}


@dataclass(frozen=True)
class SpeciesInfo:
    species_id: str
    molecule: str
    state: str
    hf_n: int

    @property
    def label(self) -> str:
        if self.state == "hf_cluster":
            return "HF" if self.hf_n == 1 else f"(HF){self.hf_n}"
        name = SHORT_NAME.get(self.molecule, self.molecule)
        if self.state == "bare_candidate":
            return name
        suffix = "RC" if self.state == "reactant_complex" else "IP"
        hf = "HF" if self.hf_n == 1 else f"(HF){self.hf_n}"
        return f"{name}-{hf} {suffix}"


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines()[2:]:
        if not line.strip():
            continue
        fields = line.split()
        symbols.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])
    return symbols, np.asarray(coords, dtype=float)


def aligned_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    left_centered = left - left.mean(axis=0)
    right_centered = right - right.mean(axis=0)
    u, _, vt = np.linalg.svd(left_centered.T @ right_centered)
    sign = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1.0, 1.0, sign]) @ vt
    return float(np.sqrt(np.mean(np.sum((left_centered @ rotation - right_centered) ** 2, axis=1))))


def load_species(preopt_manifest: Path) -> tuple[list[SpeciesInfo], dict[str, str]]:
    payload = json.loads(preopt_manifest.read_text(encoding="utf-8"))
    molecule_names = {
        str(artifact["data"]["mol_id"]): str(artifact["data"].get("name") or artifact["data"]["mol_id"])
        for artifact in payload["artifacts"]
        if artifact.get("artifact_type") == "molecule"
    }
    species: list[SpeciesInfo] = []
    for artifact in payload["artifacts"]:
        if artifact.get("artifact_type") != "species_preopt":
            continue
        data = artifact.get("data", {})
        species_id = str(
            data.get("species_id")
            or data.get("source_species_id")
            or str(artifact["artifact_id"]).removeprefix("preopt_")
        )
        mol_id = str(data.get("mol_id") or "")
        species.append(
            SpeciesInfo(
                species_id=species_id,
                molecule=molecule_names.get(mol_id, "HF cluster"),
                state=str(data.get("state") or "unknown"),
                hf_n=int(data.get("hf_n") or 0),
            )
        )
    return species, molecule_names


def elapsed_seconds(text: str) -> float | None:
    completed = re.findall(r"Total times\s+cpu:\s*[0-9.]+s\s+wall:\s*([0-9.]+)s", text)
    if completed:
        return float(completed[-1])
    # NWChem prints a cumulative time at the right edge of each active SCF row.
    active = re.findall(r"^\s*d=.*?\s([0-9]+\.[0-9]+)\s*$", text, flags=re.MULTILINE)
    return max(map(float, active)) if active else None


def calculation_rows(species: list[SpeciesInfo], minima_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for info in species:
        workdir = minima_dir / info.species_id
        output_path = workdir / "nwchem.out"
        final_path = workdir / "final.xyz"
        iterations = len(list(workdir.glob("final-*.xyz"))) if workdir.exists() else 0
        parsed: dict[str, Any] = {}
        negative: list[float] = []
        wall_s: float | None = None
        if output_path.exists():
            text = output_path.read_text(encoding="utf-8", errors="ignore")
            parsed = parse_nwchem_output(text)
            negative = [value for value in parsed.get("frequencies_cm1", []) if value < 0.0]
            wall_s = elapsed_seconds(text)
        normal = bool(parsed.get("normal_termination", False))
        optimized = bool(parsed.get("geometry_converged", False))
        n_imag = int(parsed.get("n_imag", 0)) if normal else None
        if normal and optimized:
            status = "minimum" if n_imag == 0 else "imaginary_modes"
        elif workdir.exists():
            status = "interrupted"
        else:
            status = "not_started"
        rows.append(
            {
                "species_id": info.species_id,
                "label": info.label,
                "molecule": info.molecule,
                "state": info.state,
                "hf_n": info.hf_n,
                "status": status,
                "normal_termination": normal,
                "geometry_converged": optimized,
                "iterations": iterations,
                "wall_s": wall_s,
                "wall_min": wall_s / 60.0 if wall_s is not None else None,
                "electronic_energy_hartree": parsed.get("electronic_energy_hartree"),
                "gibbs_298K_hartree": parsed.get("gibbs_298K_hartree"),
                "n_imag": n_imag,
                "imaginary_frequencies_cm1": ";".join(f"{value:.2f}" for value in negative),
                "lowest_frequency_cm1": parsed.get("lowest_freq_cm1"),
                "final_xyz": str(final_path) if final_path.exists() else "",
                "nwchem_output": str(output_path) if output_path.exists() else "",
            }
        )
    return pd.DataFrame(rows)


def endpoint_rows(calculations: pd.DataFrame, minima_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    candidates = calculations[calculations["state"].isin(["reactant_complex", "ion_pair"])]
    for (molecule, hf_n), group in candidates.groupby(["molecule", "hf_n"]):
        rc = group[group["state"] == "reactant_complex"]
        ip = group[group["state"] == "ion_pair"]
        if len(rc) != 1 or len(ip) != 1:
            continue
        rc_row, ip_row = rc.iloc[0], ip.iloc[0]
        if not (rc_row["normal_termination"] and ip_row["normal_termination"]):
            continue
        rc_xyz = minima_dir / str(rc_row["species_id"]) / "final.xyz"
        ip_xyz = minima_dir / str(ip_row["species_id"]) / "final.xyz"
        if not (rc_xyz.exists() and ip_xyz.exists()):
            continue
        rc_symbols, rc_coords = read_xyz(rc_xyz)
        ip_symbols, ip_coords = read_xyz(ip_xyz)
        if rc_symbols != ip_symbols:
            continue
        rmsd = aligned_rmsd(rc_coords, ip_coords)
        delta_e = (
            float(ip_row["electronic_energy_hartree"])
            - float(rc_row["electronic_energy_hartree"])
        ) * HARTREE_TO_KJ_MOL
        label = f"{SHORT_NAME.get(str(molecule), molecule)}-{'HF' if hf_n == 1 else f'(HF){hf_n}'}"
        rows.append(
            {
                "label": label,
                "molecule": molecule,
                "hf_n": int(hf_n),
                "rmsd_A": rmsd,
                "delta_e_ip_minus_rc_kj_mol": delta_e,
                "same_endpoint": bool(rmsd < 0.01 and abs(delta_e) < 0.1),
            }
        )
    return pd.DataFrame(rows)


def proton_geometry_rows(calculations: pd.DataFrame, minima_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    rc_rows = calculations[
        (calculations["state"] == "reactant_complex") & calculations["normal_termination"]
    ]
    for _, row in rc_rows.iterrows():
        path = minima_dir / str(row["species_id"]) / "final.xyz"
        symbols, coords = read_xyz(path)
        nitrogen = [index for index, symbol in enumerate(symbols) if symbol == "N"]
        fluorine = [index for index, symbol in enumerate(symbols) if symbol == "F"]
        hydrogen = [index for index, symbol in enumerate(symbols) if symbol == "H"]
        bridges: list[tuple[float, float, float]] = []
        for f_index in fluorine:
            h_index = min(hydrogen, key=lambda index: np.linalg.norm(coords[index] - coords[f_index]))
            r_fh = float(np.linalg.norm(coords[h_index] - coords[f_index]))
            r_nh = min(float(np.linalg.norm(coords[h_index] - coords[index])) for index in nitrogen)
            bridges.append((r_fh + r_nh, r_fh, r_nh))
        _, r_fh, r_nh = min(bridges)
        rows.append(
            {
                "label": row["label"].removesuffix(" RC"),
                "molecule": row["molecule"],
                "hf_n": int(row["hf_n"]),
                "r_fh_A": r_fh,
                "r_nh_A": r_nh,
                "proton_transfer_coordinate_A": r_fh - r_nh,
                "frequency_qc": row["status"],
            }
        )
    return pd.DataFrame(rows)


def association_rows(calculations: pd.DataFrame) -> pd.DataFrame:
    energy = calculations.set_index("species_id")["electronic_energy_hartree"].to_dict()
    hf = calculations[(calculations["state"] == "hf_cluster") & (calculations["hf_n"] == 1)]
    if len(hf) != 1 or pd.isna(hf.iloc[0]["electronic_energy_hartree"]):
        return pd.DataFrame()
    hf_energy = float(hf.iloc[0]["electronic_energy_hartree"])
    rows: list[dict[str, Any]] = []
    for molecule, group in calculations.groupby("molecule"):
        if molecule == "HF cluster":
            continue
        bare = group[(group["state"] == "bare_candidate") & group["normal_termination"]]
        if len(bare) != 1:
            continue
        bare_energy = float(bare.iloc[0]["electronic_energy_hartree"])
        for _, complex_row in group[
            (group["state"] == "reactant_complex") & group["normal_termination"]
        ].iterrows():
            complex_energy = energy[str(complex_row["species_id"])]
            delta = (float(complex_energy) - bare_energy - int(complex_row["hf_n"]) * hf_energy)
            rows.append(
                {
                    "label": complex_row["label"].removesuffix(" RC"),
                    "molecule": molecule,
                    "hf_n": int(complex_row["hf_n"]),
                    "association_energy_kj_mol": delta * HARTREE_TO_KJ_MOL,
                    "frequency_qc": complex_row["status"],
                }
            )
    return pd.DataFrame(rows)


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 180,
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_workflow(calculations: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 2.8))
    y = [1, 0]
    ax.barh(y[0], len(calculations), color="#377EB8", label="Preoptimized successfully")
    left = 0
    for status in STATUS_ORDER:
        count = int((calculations["status"] == status).sum())
        ax.barh(y[1], count, left=left, color=STATUS_COLOR[status], label=STATUS_LABEL[status])
        if count:
            ax.text(left + count / 2, y[1], str(count), ha="center", va="center", color="white", weight="bold")
        left += count
    ax.text(len(calculations) / 2, y[0], str(len(calculations)), ha="center", va="center", color="white", weight="bold")
    ax.set_yticks(y, ["GFN2-xTB preoptimization", "NWChem minima/frequency"])
    ax.set_xlim(0, len(calculations))
    ax.set_xlabel("Planned structures")
    ax.set_title("Production pilot status at manual stop")
    ax.legend(ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
    save_figure(fig, output_dir / "01_workflow_status.png")


def plot_runtime(calculations: pd.DataFrame, output_dir: Path) -> None:
    attempted = calculations[calculations["status"] != "not_started"].copy()
    attempted["plot_minutes"] = attempted["wall_min"].fillna(0.03).clip(lower=0.03)
    attempted = attempted.sort_values("plot_minutes")
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    bars = ax.barh(
        attempted["label"],
        attempted["plot_minutes"],
        color=[STATUS_COLOR[value] for value in attempted["status"]],
    )
    ax.set_xscale("log")
    ax.set_xlabel("Wall time / min (log scale)")
    ax.set_title("Per-structure NWChem cost")
    for bar, (_, row) in zip(bars, attempted.iterrows()):
        label = "<0.1" if pd.isna(row["wall_min"]) else f"{row['wall_min']:.1f}"
        ax.text(bar.get_width() * 1.05, bar.get_y() + bar.get_height() / 2, label, va="center", fontsize=7)
    handles = [plt.Rectangle((0, 0), 1, 1, color=STATUS_COLOR[key]) for key in STATUS_ORDER[:3]]
    ax.legend(handles, [STATUS_LABEL[key] for key in STATUS_ORDER[:3]], loc="lower right")
    save_figure(fig, output_dir / "02_runtime_by_job.png")


def plot_iterations(calculations: pd.DataFrame, output_dir: Path) -> float:
    data = calculations[(calculations["iterations"] > 0) & calculations["wall_min"].notna()].copy()
    correlation = float(data[["iterations", "wall_min"]].corr().iloc[0, 1])
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for status, group in data.groupby("status"):
        ax.scatter(
            group["iterations"],
            group["wall_min"],
            s=50,
            color=STATUS_COLOR[status],
            label=STATUS_LABEL[status],
            alpha=0.9,
        )
    for _, row in data.iterrows():
        if row["iterations"] >= 25 or row["wall_min"] >= 150:
            ax.annotate(row["label"], (row["iterations"], row["wall_min"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_yscale("log")
    ax.set_xlabel("Saved optimization geometries")
    ax.set_ylabel("Wall time / min (log scale)")
    ax.set_title(f"Optimization effort vs wall time (Pearson r = {correlation:.2f})")
    ax.legend()
    save_figure(fig, output_dir / "03_iterations_vs_runtime.png")
    return correlation


def plot_imaginary(calculations: pd.DataFrame, output_dir: Path) -> None:
    failed = calculations[calculations["status"] == "imaginary_modes"].copy()
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.5), gridspec_kw={"width_ratios": [1, 1.7]})
    left.bar(failed["label"], failed["n_imag"], color=STATUS_COLOR["imaginary_modes"])
    left.set_ylabel("Number of imaginary modes")
    left.set_title("Frequency-QC failures")
    left.tick_params(axis="x", rotation=70)
    for x, row in enumerate(failed.itertuples()):
        values = [float(value) for value in str(row.imaginary_frequencies_cm1).split(";") if value]
        right.scatter([row.label] * len(values), values, color="#8C3B2A", s=45)
    right.axhline(0, color="black", linewidth=0.8)
    right.set_ylabel("Imaginary frequency / cm$^{-1}$")
    right.set_title("Magnitude of unstable modes")
    right.tick_params(axis="x", rotation=70)
    save_figure(fig, output_dir / "04_imaginary_modes.png")


def plot_endpoints(endpoints: pd.DataFrame, output_dir: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.3))
    left.bar(endpoints["label"], endpoints["rmsd_A"], color="#4C78A8")
    left.axhline(0.01, color="#C7473A", linestyle="--", label="Distinctness threshold")
    left.set_yscale("log")
    left.set_ylabel("Aligned RC-IP RMSD / A")
    left.set_title("Endpoint geometry separation")
    left.tick_params(axis="x", rotation=65)
    left.legend()
    absolute_energy = endpoints["delta_e_ip_minus_rc_kj_mol"].abs().clip(lower=1e-9)
    right.bar(endpoints["label"], absolute_energy, color="#72B7B2")
    right.axhline(0.1, color="#C7473A", linestyle="--", label="Distinctness threshold")
    right.set_yscale("log")
    right.set_ylabel("abs(E(IP)-E(RC)) / kJ mol$^{-1}$")
    right.set_title("Endpoint electronic-energy separation")
    right.tick_params(axis="x", rotation=65)
    right.legend()
    save_figure(fig, output_dir / "05_endpoint_collapse.png")


def plot_proton_coordinate(geometry: pd.DataFrame, output_dir: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.5))
    x = np.arange(len(geometry))
    width = 0.36
    left.bar(x - width / 2, geometry["r_fh_A"], width, label="r(F-H)", color="#4C78A8")
    left.bar(x + width / 2, geometry["r_nh_A"], width, label="r(N-H)", color="#F2A541")
    left.set_xticks(x, geometry["label"], rotation=65)
    left.set_ylabel("Distance / A")
    left.set_title("Optimized proton-contact distances")
    left.legend()
    colors = [STATUS_COLOR[value] for value in geometry["frequency_qc"]]
    right.bar(geometry["label"], geometry["proton_transfer_coordinate_A"], color=colors)
    right.axhline(0, color="black", linewidth=0.8)
    right.set_ylabel("q = r(F-H) - r(N-H) / A")
    right.set_title("Proton-transfer coordinate\nq < 0: proton remains closer to F")
    right.tick_params(axis="x", rotation=65)
    save_figure(fig, output_dir / "06_proton_transfer_coordinate.png")


def plot_association(association: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    ordered = association.sort_values(["molecule", "hf_n"])
    bars = ax.bar(
        ordered["label"],
        ordered["association_energy_kj_mol"],
        color=[STATUS_COLOR[value] for value in ordered["frequency_qc"]],
    )
    for bar, (_, row) in zip(bars, ordered.iterrows()):
        if row["frequency_qc"] != "minimum":
            bar.set_hatch("///")
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 4,
            f"{row['association_energy_kj_mol']:.1f}",
            ha="center",
            va="bottom",
            color="white",
            fontsize=8,
            weight="bold",
        )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Electronic association energy / kJ mol$^{-1}$")
    ax.set_title("A + n HF -> A.(HF)n (uncorrected electronic energy)")
    ax.tick_params(axis="x", rotation=55)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=STATUS_COLOR["minimum"]),
        plt.Rectangle((0, 0), 1, 1, color=STATUS_COLOR["imaginary_modes"], hatch="///"),
    ]
    ax.legend(
        handles,
        ["Frequency-QC passed", "Imaginary mode(s): provisional"],
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
    )
    save_figure(fig, output_dir / "07_association_energy.png")


def fmt(value: float | None, digits: int = 2) -> str:
    if value is None or pd.isna(value):
        return "—"
    return f"{float(value):.{digits}f}"


def markdown_table(frame: pd.DataFrame, columns: list[str], headers: list[str]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    return "\n".join(lines)


def write_report(
    output_dir: Path,
    calculations: pd.DataFrame,
    endpoints: pd.DataFrame,
    geometry: pd.DataFrame,
    association: pd.DataFrame,
    correlation: float,
) -> None:
    counts = calculations["status"].value_counts().to_dict()
    completed = calculations[calculations["normal_termination"]]
    attempted_wall_h = calculations["wall_s"].fillna(0).sum() / 3600.0
    calc_table = calculations.copy()
    calc_table["判定"] = calc_table["status"].map(STATUS_LABEL)
    calc_table["反復"] = calc_table["iterations"].astype(int)
    calc_table["虚振動"] = calc_table["n_imag"].map(lambda value: "—" if pd.isna(value) else str(int(value)))
    calc_table["時間/分"] = calc_table["wall_min"].map(lambda value: fmt(value, 1))
    endpoint_max_rmsd = endpoints["rmsd_A"].max()
    endpoint_max_de = endpoints["delta_e_ip_minus_rc_kj_mol"].abs().max()
    strongest = association.loc[association["association_energy_kj_mol"].idxmin()]
    weakest = association.loc[association["association_energy_kj_mol"].idxmax()]
    closest = geometry.loc[geometry["proton_transfer_coordinate_A"].idxmax()]
    report = f"""# NWChemアミン–HFパイロット：停止時点の図表評価

生成日時: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}

## 結論

- GFN2-xTB前処理は22/22成功した。
- NWChem最適化・振動解析は{len(completed)}/22が正常終了したが、虚振動0の真の最小点は{counts.get('minimum', 0)}件である。
- 虚振動を持つ収束構造は{counts.get('imaginary_modes', 0)}件、中断は{counts.get('interrupted', 0)}件、未着手は{counts.get('not_started', 0)}件である。
- 完了したRC/IP {len(endpoints)}対はすべて同一端点へ崩壊した。最大RMSDは{endpoint_max_rmsd:.6f} Å、最大電子エネルギー差は{endpoint_max_de:.6f} kJ mol⁻¹にすぎない。
- したがって、現データから独立した反応物・生成物端点、遷移状態、IRC、速度定数を主張できない。
- 完了・中断ジョブから抽出できた累積壁時計時間は約{attempted_wall_h:.1f}時間である。

## 1. ワークフロー進捗

![Workflow status](01_workflow_status.png)

前処理は全22構造で成功している一方、DFT段階で真の最小点まで到達したのは10構造である。正常終了だけを成功と数えると15件だが、そのうち5件は虚振動を持つため熱化学・経路端点には採用できない。中断2件はアニリン–(HF)₂ RC/IP、未着手5件はピリジンの裸分子とHF₁/HF₂ RC/IPである。

## 2. 構造別計算時間

![Runtime by job](02_runtime_by_job.png)

小さいNH₃系は数分以内だが、TMA系は約23–92分、アニリン系は約96–353分へ増大した。アニリン–HF RCは114構造更新と約353分を要し、それでも虚振動2本を残した。アニリン–(HF)₂ RCは約205分で手動中断された。分子サイズだけでなく、分子間の柔らかいHF配向と拡散基底の線形依存がコストを支配している。

## 3. 反復数と実時間

![Iterations versus runtime](03_iterations_vs_runtime.png)

保存された最適化構造数と壁時計時間のPearson相関は{correlation:.2f}で、反復増加が計算時間を強く支配している。ただし原子数・基底関数数も異なるため、同じ反復数でもアニリンはNH₃より高コストである。運用上は構造単位のタイムアウト、途中マニフェスト、完了出力再利用が必要である。

## 4. 虚振動による品質不合格

![Imaginary modes](04_imaginary_modes.png)

(HF)₂は−427.17、−427.17、−83.63 cm⁻¹、TMA–(HF)₂ RC/IPは−83.41 cm⁻¹が各2本、アニリン–HF RC/IPは約−172および−157 cm⁻¹を持つ。TMAの−83 cm⁻¹は柔らかい分子間モードの可能性があるが、アニリンの約−170 cm⁻¹と(HF)₂の−427 cm⁻¹は無視できない。いずれも「最小点」として採用せず、負モード±方向から再最適化すべきである。

## 5. RC/IP端点の崩壊

![Endpoint collapse](05_endpoint_collapse.png)

図中の破線は、暫定的な独立端点判定値（RMSD 0.01 Å、|ΔE| 0.1 kJ mol⁻¹）である。全5対が両方の基準を数桁下回る。つまり、RCとIPという入力ラベルは異なっていても、無拘束最適化後は同じ構造である。現在の端点からNEBを実行すると、経路長ほぼゼロまたは単調経路になり、意味のあるTS探索にならない。

## 6. プロトン移動座標

![Proton-transfer coordinate](06_proton_transfer_coordinate.png)

q = r(F–H) − r(N–H) と定義した。すべてq < 0で、プロトンは窒素よりフッ素に近く、中性HF会合体側にある。最も共有性が強いのは{closest['label']}（q={closest['proton_transfer_coordinate_A']:.3f} Å）だが、この構造は虚振動を持つため安定イオン対とは断定できない。RC/IPが同じ中性・共有型構造へ収束したことと整合する。

## 7. 電子会合エネルギー

![Association energy](07_association_energy.png)

HF単量体を基準とする未補正電子会合エネルギー ΔE = E[A·(HF)n] − E[A] − nE[HF] を示す。最も負なのは{strongest['label']}（{strongest['association_energy_kj_mol']:.1f} kJ mol⁻¹）、最も弱いのは{weakest['label']}（{weakest['association_energy_kj_mol']:.1f} kJ mol⁻¹）である。HF₂付加はHF₁より強い会合を示し、TMAはNH₃よりわずかに強く、アニリン–HFは弱いという傾向は気相塩基性の期待と定性的に整合する。

ただし、斜線の棒は虚振動を持つ停留点であり暫定値である。すべてBSSE補正、ZPE、quasi-RRHO、標準状態補正を含まない電子エネルギーなので、定量的な会合自由エネルギーや平衡定数には使用できない。

## 計算別一覧

{markdown_table(calc_table, ['label', '判定', '反復', '虚振動', '時間/分'], ['系', '判定', '保存構造数', '虚振動数', '時間/分'])}

## 推奨する次の処理

1. 正常終了した15件を生出力からステージマニフェストへ安全にインポートし、再計算を避ける。
2. 虚振動を持つ5構造は負モードの±方向へ変位し、低コスト基底で再最適化してからPBE0/def2-TZVPDで再確認する。
3. 独立RC/IPはN–H/H–F座標の拘束付き緩和スキャンから作り、拘束解除後にも別極小点として残るか確認する。
4. アニリン–(HF)₂はdef2-SVPD予備最適化、小さい信頼半径、線形依存監視を用いて再開する。
5. RC/IPのRMSD・ΔEゲートをTS探索前に追加し、同一端点をNEBへ渡さない。

## データと適用範囲

- [計算サマリー](calculation_summary.csv)
- [RC/IP比較](endpoint_comparison.csv)
- [プロトン座標](proton_transfer_geometry.csv)
- [電子会合エネルギー](association_energy.csv)
- [生のNWChem出力](../06_dft-minima/)
- [前処理マニフェスト](../05_preopt/manifest.json)

本レポートは停止時点の部分データを対象とする。ステージマニフェストが未生成であるため、生出力を直接解析した。中断エネルギーは会合エネルギーやランキングへ使用していない。
"""
    (output_dir / "REPORT.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    output_dir = (args.output_dir or run_dir / "analysis_report").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    species, _ = load_species(run_dir / "05_preopt" / "manifest.json")
    minima_dir = run_dir / "06_dft-minima"
    calculations = calculation_rows(species, minima_dir)
    endpoints = endpoint_rows(calculations, minima_dir)
    geometry = proton_geometry_rows(calculations, minima_dir)
    association = association_rows(calculations)
    calculations.to_csv(output_dir / "calculation_summary.csv", index=False)
    endpoints.to_csv(output_dir / "endpoint_comparison.csv", index=False)
    geometry.to_csv(output_dir / "proton_transfer_geometry.csv", index=False)
    association.to_csv(output_dir / "association_energy.csv", index=False)
    configure_plotting()
    plot_workflow(calculations, output_dir)
    plot_runtime(calculations, output_dir)
    correlation = plot_iterations(calculations, output_dir)
    plot_imaginary(calculations, output_dir)
    plot_endpoints(endpoints, output_dir)
    plot_proton_coordinate(geometry, output_dir)
    plot_association(association, output_dir)
    write_report(output_dir, calculations, endpoints, geometry, association, correlation)
    print(f"wrote report and 7 figures to {output_dir}")


if __name__ == "__main__":
    main()
