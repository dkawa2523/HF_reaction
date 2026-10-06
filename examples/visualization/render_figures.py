"""Render inspectable scientific figures from saved hfauto observations.

This renderer never starts quantum chemistry. Energy plots show discrete stored
thermochemical points or separately labelled sampled electronic profiles. Bonds
in the ball-and-stick panels are distance-based depiction, not assigned chemistry.
"""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from mpl_toolkits.mplot3d import proj3d

HERE = Path(__file__).resolve().parent
COLORS = {
    "ink": "#17334C", "muted": "#617588", "grid": "#DDE5EB", "paper": "#FFFFFF",
    "reactant": "#2A7590", "ts": "#C37126", "product": "#75549A",
    "pass": "#287664", "hold": "#A36727", "fail": "#A54D4F",
}
ATOM_COLORS = {
    "H": "#F5F7FA", "C": "#6B7683", "N": "#3C6BB6", "O": "#D14B4B",
    "F": "#67B478", "B": "#D8A58B",
}
RADII = {"H": .31, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57}
TITLES = {
    "nh3_inversion": ("01", "NH₃ の反転", "鞍点から両側の極小へ接続を確認"),
    "nh3_hf_exchange": ("02", "NH₃・HF の水素交換", "未知候補の探索から DFT 認証・熱化学比較へ"),
    "tma_hf2_same_basin": ("03", "TMA・(HF)₂ の端点判定", "別の入力候補が同じ極小へ戻る"),
    "bh3_nh3_association": ("04", "BH₃ ＋ NH₃ の会合", "保存された電子エネルギー走査で経路を読む"),
    "oh_ch4_abstraction": ("05", "OH ＋ CH₄ の水素引き抜き", "開殻反応の停留点・接続・比較基準を確認"),
}
ROLE_NAMES = {
    "reactant": "反応物側", "product": "生成物側", "ts": "遷移状態 (TS)",
    "declared_reactant": "宣言した反応物入力", "declared_product": "宣言した生成物入力",
    "minimum": "収束した同じ極小", "complex": "錯体", "separated": "分離した単量体",
}
ENERGY_LABELS = {"R": "R\n反応物側", "P": "P\n生成物側", "TS": "TS",
                 "R separated": "分離した R", "R (separated)": "分離した R"}


def configure() -> None:
    """Use a font shipped on this Windows host; retain vector text in SVG."""
    candidates = [Path("/mnt/c/Windows/Fonts/meiryo.ttc"), Path("C:/Windows/Fonts/meiryo.ttc")]
    for path in candidates:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
            bold_path = path.with_name("meiryob.ttc")
            if bold_path.exists():
                font_manager.fontManager.addfont(str(bold_path))
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(path)).get_name()
            break
    plt.rcParams.update({
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 10.5,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "text.color": COLORS["ink"], "axes.labelcolor": COLORS["ink"],
        "xtick.color": COLORS["muted"], "ytick.color": COLORS["muted"],
        "axes.edgecolor": COLORS["grid"], "svg.fonttype": "none",
        "savefig.facecolor": "white", "axes.unicode_minus": False,
    })


def wrapped(value: object, width: int = 40) -> str:
    return "\n".join(part for line in str(value).splitlines()
                     for part in textwrap.wrap(line, width=width, break_long_words=True))


def numeric(value: object) -> float | None:
    return float(value) if isinstance(value, (int, float)) else None


def panel_title(ax, title: str) -> None:
    ax.set_title(title, loc="left", fontweight="bold", pad=12)


def clean_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=COLORS["grid"], lw=.7, zorder=0)


def energy_diagram(ax, example: dict) -> None:
    """Show stored free-energy states; dashed joins only establish point order."""
    free = example.get("free_energy")
    panel_title(ax, "A  自由エネルギーの離散点")
    if not free or not free.get("points"):
        ax.set_axis_off()
        notes = {
            "tma_hf2_same_basin": "同じ極小への合流\n認証した TS はないため、障壁図は作らない。",
            "bh3_nh3_association": "会合は通常の TS 障壁順位と分ける。\n右の走査は電子エネルギー E であり、G ではない。",
        }
        ax.text(.04, .74, notes.get(example["id"], "この例の順位に使える G は未登録。"),
                transform=ax.transAxes, fontsize=12.5, va="top", linespacing=1.7)
        ax.text(.04, .28, "未観測の点を補間したり、仮の障壁を追加したりしない。",
                transform=ax.transAxes, fontsize=9.5, color=COLORS["muted"], va="top")
        return
    points = [p for p in free["points"] if numeric(p.get("relative_kcal")) is not None]
    xs = np.arange(len(points))
    ys = [p["relative_kcal"] for p in points]
    for i, p in enumerate(points):
        color = COLORS.get(p.get("role", ""), COLORS["reactant"])
        ax.hlines(ys[i], i - .19, i + .19, color=color, lw=3.5, zorder=3)
        offset = 10 if ys[i] >= 0 else -17
        ax.annotate(f"{ys[i]:+.2f}" if abs(ys[i]) > .005 else "0.00", (i, ys[i]),
                    xytext=(0, offset), textcoords="offset points", ha="center", fontsize=10.5,
                    color=color, fontweight="bold")
    if len(points) > 1:
        ax.plot(xs, ys, "--", color=COLORS["muted"], lw=1, zorder=1, alpha=.7)
    span = max(max(ys) - min(ys), 1.0)
    ax.set_ylim(min(ys) - .24 * span, max(ys) + .38 * span)
    ax.set_xlim(-.55, len(points) - .45)
    labels = [ENERGY_LABELS.get(p.get("label"),
              p.get("label", ROLE_NAMES.get(p.get("role"), p.get("role", "")))) for p in points]
    ax.set_xticks(xs, labels)
    ax.set_ylabel("相対 G / kcal mol⁻¹")
    clean_axis(ax)
    zero_label = "登録 R 極小"
    if example["id"] in {"nh3_hf_exchange", "oh_ch4_abstraction"}:
        zero_label = "R 錯体"
    if free.get("zero") == "declared separated monomers":
        zero_label = "分離した単量体"
    ax.text(.02, .95, "G = 0: " + zero_label, transform=ax.transAxes, fontsize=9,
            color=COLORS["muted"], va="top")
    ax.text(0, -.24, "破線は離散点の順序。連続したポテンシャル曲線ではない。",
            transform=ax.transAxes, fontsize=8.5, color=COLORS["muted"])
    value = free.get("dG_eff_kcal")
    if numeric(value) is not None:
        ax.text(.98, .95, f"δG_eff = {value:.2f} kcal mol⁻¹", ha="right", va="top",
                transform=ax.transAxes, color=COLORS["ink"], fontsize=10.5,
                bbox={"boxstyle": "round,pad=.4", "fc": "#F0F5F8", "ec": "none"})


def electronic_profile(ax, example: dict) -> None:
    profiles = example.get("electronic_profiles") or []
    title = "B  低レベル候補の保存 IRC（未認証）" if example["id"] == "tma_hf2_same_basin" \
        else "B  計算済み電子エネルギー走査"
    panel_title(ax, title)
    if not profiles:
        ax.set_axis_off()
        ax.text(.02, .72, "この図で採用した連続走査はない。\n保存済みの停留点を A に示す。",
                transform=ax.transAxes, va="top", fontsize=11.5, linespacing=1.7)
        ax.text(.02, .26, "A の点と構造を合わせて、\n停留性・接続・同一状態の判定を読む。",
                transform=ax.transAxes, va="top", fontsize=9.5, color=COLORS["muted"],
                linespacing=1.6)
        return
    for index, profile in enumerate(profiles):
        points = [p for p in profile.get("points", [])
                  if numeric(p.get("x")) is not None and numeric(p.get("relative_kcal")) is not None]
        if not points:
            continue
        ax.plot([p["x"] for p in points], [p["relative_kcal"] for p in points], "o-",
                ms=4.2, lw=1.1, color=[COLORS["reactant"], COLORS["product"]][index % 2],
                label=f"走査 {index + 1}" if len(profiles) > 1 else "保存された計算点")
    profile = profiles[0]
    coordinate = profile.get("coordinate_kind", "走査座標")
    coordinate = ("保存 IRC 反復数（符号付き、距離の単位ではない）" if "IRC iteration" in coordinate
                  else "B–N 距離 / Å" if "B-N distance" in coordinate else coordinate)
    ax.set_xlabel(wrapped(coordinate, 32))
    ax.set_ylabel("相対 E / kcal mol⁻¹")
    clean_axis(ax)
    zero = {"last saved backward IRC frame": "最終 backward IRC 点",
            "optimized adduct": "最適化した付加体"}.get(profile.get("zero"), str(profile.get("zero", "")))
    ax.text(.02, .95, "E = 0: " + zero, transform=ax.transAxes, fontsize=9,
            color=COLORS["muted"], va="top")
    if len(profiles) > 1:
        ax.legend(loc="best", frameon=False, fontsize=9)
    ax.text(0, -.24, "点間の線は計算点の結線。スプライン補間はしていない。",
            transform=ax.transAxes, fontsize=8.5, color=COLORS["muted"])


def decision_lines(example: dict) -> list[tuple[str, str]]:
    reaction, decisions = example.get("reaction", {}), example.get("decisions", {})
    outcome = reaction.get("outcome", "未登録")
    labels = {
        "ts_valid": "TS を認証", "same_basin": "同じ極小へ合流",
        "barrierless": "走査の分解能内で障壁なし", "accepted": "経路を認証",
        "capture_limited": "捕獲律速として扱う", "certified": "経路を認証",
        "degenerate_rearrangement": "同じ化学状態を結ぶ TS を認証",
        "elementary_step": "異なる化学状態を結ぶ TS を認証",
        "barrierless_at_resolution": "走査の分解能内で障壁なし",
    }
    lines = [("判定", labels.get(outcome, str(outcome)))]
    imag = decisions.get("imag_cm1")
    if isinstance(imag, list):
        values = [v for v in imag if numeric(v) is not None and v < 0]
        if values:
            lines.append(("虚振動", ", ".join(f"{v:.1f}" for v in values) + " cm⁻¹"))
    elif numeric(imag) is not None:
        lines.append(("虚振動", f"{imag:.1f} cm⁻¹"))
    connection = decisions.get("connection_label")
    connection = {"degenerate_rearrangement": "両側から同じ登録極小へ接続",
                  "elementary_step": "両側から別の登録極小へ接続",
                  "same_basin": "宣言端点が同じ極小へ収束",
                  "barrierless_at_resolution": "拘束走査を確認。TS は採用しない"}.get(connection, connection)
    if not connection and decisions.get("connection"):
        minima = decisions["connection"].get("minima", [])
        connection = "両側から同じ登録極小へ接続" if len(set(minima)) == 1 else "両側から別の登録極小へ接続"
    if connection:
        lines.append(("接続", str(connection)))
    if reaction.get("degenerate"):
        lines.append(("化学状態", "両側は同じ状態。原子対応・配置は保持。"))
    free = example.get("free_energy")
    if free:
        ref = free.get("reference")
        if ref:
            lines.append(("比較基準", {"complex": "反応物錯体", "separated": "分離した反応物"}.get(ref, str(ref))))
        if numeric(free.get("dG_rxn_kcal")) is not None:
            heading = "反応自由エネルギー（錯体間）" if example["id"] == "oh_ch4_abstraction" \
                else "反応自由エネルギー"
            lines.append((heading, f"{free['dG_rxn_kcal']:+.2f} kcal mol⁻¹"))
    blockers = decisions.get("blockers") or (free or {}).get("blockers") or []
    if blockers:
        labels = {"outcome:same_basin": "同じ極小への合流", "thermo_unavailable": "順位に使う熱化学なし",
                  "outcome:barrierless_at_resolution": "障壁順位から除外（捕獲律速）"}
        lines.append(("順位の扱い", "; ".join(labels.get(x, str(x)) for x in blockers)))
    reasons = reaction.get("reasons") or []
    if reasons and not blockers and not reaction.get("degenerate"):
        reasons = ["両端の化学状態を照合し、接続を認証" if x == "connection:elementary" else str(x)
                   for x in reasons]
        lines.append(("根拠", "; ".join(reasons)))
    return lines[:6]


def decision_panel(ax, example: dict) -> None:
    panel_title(ax, "C  反応判断の根拠")
    ax.set_axis_off()
    y = .95
    for heading, value in decision_lines(example):
        ax.text(.02, y, heading, transform=ax.transAxes, fontsize=9, fontweight="bold",
                color=COLORS["muted"], va="top")
        lines = wrapped(value, 25)
        ax.text(.02, y - .07, lines, transform=ax.transAxes, fontsize=10.5, va="top",
                linespacing=1.5)
        y -= .12 + .065 * len(lines.splitlines())
    if y > .14:
        ax.text(.02, .04, "計算の終了 ≠ 科学的認証\n順位 ≠ 初期系の生成物比",
                transform=ax.transAxes, fontsize=9, color=COLORS["muted"], va="bottom",
                linespacing=1.6)


def atom_coords(structure: dict) -> np.ndarray:
    return np.array([[atom["x"], atom["y"], atom["z"]] for atom in structure["atoms"]], dtype=float)


def depiction_bonds(structure: dict) -> set[tuple[int, int]]:
    """Distance criterion for drawing only; not a bond order or state assignment."""
    xyz, atoms = atom_coords(structure), structure["atoms"]
    bonds = set()
    for i in range(len(atoms)):
        for j in range(i):
            limit = 1.20 * (RADII.get(atoms[i]["element"], .8)
                            + RADII.get(atoms[j]["element"], .8))
            distance = np.linalg.norm(xyz[i] - xyz[j])
            if .35 < distance < limit:
                bonds.add((j, i))
    return bonds


def aligned_coords(structure: dict, reference: dict) -> np.ndarray:
    """Proper rigid rotation preserves chirality and original atom indices."""
    mobile, fixed = atom_coords(structure), atom_coords(reference)
    mobile = mobile - mobile.mean(axis=0)
    if mobile.shape != fixed.shape:
        return mobile
    if [a["element"] for a in structure["atoms"]] != [a["element"] for a in reference["atoms"]]:
        return mobile
    fixed = fixed - fixed.mean(axis=0)
    u, _, vt = np.linalg.svd(mobile.T @ fixed)
    rotation = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    return mobile @ rotation


def selected_structures(example: dict) -> list[dict]:
    structures = example.get("structures", [])
    roles = (["declared_reactant", "minimum", "declared_product"]
             if example["id"] == "tma_hf2_same_basin" else ["reactant", "ts", "product"])
    selected = [next((s for s in structures if s["role"] == role), None) for role in roles]
    return [s for s in selected if s is not None]


def atom_labels(ax, xyz: np.ndarray, atoms: list[dict], changed: set) -> None:
    """Place index labels outside the projected spheres, keeping them legible."""
    projected = np.asarray(proj3d.proj_transform(*xyz.T, ax.get_proj())[:2]).T
    display = ax.transData.transform(projected)
    pixel_per_pt = ax.figure.dpi / 72
    used = []
    offsets = [(0, 17), (17, 8), (-17, 8), (17, -12), (-17, -12), (0, -19)]
    for i, atom in enumerate(atoms):
        if len(atoms) > 10 and atom["element"] == "H" and not any(i in pair for pair in changed):
            continue
        label = f"{atom['element']}{atom.get('index', i) + 1}"
        width, height = len(label) * 5.5 * pixel_per_pt, 11 * pixel_per_pt
        candidates = []
        for dx, dy in offsets:
            center = display[i] + np.array([dx, dy]) * pixel_per_pt
            bounds = (center[0] - width / 2, center[1] - height / 2,
                      center[0] + width / 2, center[1] + height / 2)
            collisions = sum(not (bounds[2] < b[0] or bounds[0] > b[2]
                                  or bounds[3] < b[1] or bounds[1] > b[3]) for b in used)
            nearest = np.min(np.linalg.norm(display - center, axis=1))
            candidates.append((collisions * 1000 - nearest, dx, dy, bounds))
        _, dx, dy, bounds = min(candidates)
        used.append(bounds)
        ax.annotate(label, projected[i], xytext=(dx, dy), textcoords="offset points",
                    ha="center", va="center", fontsize=8.5, color=COLORS["ink"],
                    arrowprops={"arrowstyle": "-", "color": "#98A6B1", "lw": .55}, zorder=100)


def draw_molecule(ax, structure: dict, reference: dict, changed: set, extent: float) -> None:
    xyz = aligned_coords(structure, reference)
    atoms = structure["atoms"]
    displayed = depiction_bonds(structure) | (changed if structure["role"] == "ts" else set())
    for i, j in sorted(displayed):
        line = xyz[[i, j]]
        is_changed = (i, j) in changed and structure["role"] == "ts"
        ax.plot(*line.T, color=COLORS["ts"] if is_changed else "#718393",
                lw=1.5 if is_changed else 3.8, ls="--" if is_changed else "-", alpha=.9)
    sizes = [300 if a["element"] == "H" else 570 for a in atoms]
    ax.scatter(*xyz.T, s=sizes, c=[ATOM_COLORS.get(a["element"], "#A8A3BB") for a in atoms],
               edgecolors="#657583", linewidths=.65, depthshade=True, zorder=4)
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_zlim(-extent, extent)
    ax.set_box_aspect((1, 1, 1), zoom=1.42)
    ax.view_init(elev=19, azim=-43)
    ax.set_axis_off()
    atom_labels(ax, xyz, atoms, changed)
    raw_label = structure.get("label", "")
    labels = {"reactant": "R 側の極小", "product": "P 側の極小", "certified TS": "認証した TS",
              "QRC side +": "R: +側の接続確認構造", "QRC side -": "P: −側の接続確認構造",
              "QRC side - (symmetry image)": "P: QRC+ の対称像\n（独立 QM ではない）",
              "neutral": "R 入力: neutral", "shared_proton": "P 入力: shared_proton",
              "common DFT minimum": "収束した同じ DFT 極小",
              "declared entrance complex": "R 入力: 会合の入口配置",
              "optimized adduct": "P: 最適化した付加体"}
    label = labels.get(raw_label, raw_label or ROLE_NAMES.get(structure["role"], structure["role"]))
    ax.set_title(wrapped(label, 18), fontsize=10.5, pad=-2, fontweight="bold")
    observed = structure.get("observed_role", "")
    if observed:
        observed = {"saved DFT mode-follow minimum": "保存された DFT mode-follow 極小",
                    "saved DFT QRC endpoint": "保存された DFT 接続確認端点",
                    "saved DFT saddle frequency geometry": "DFT 鞍点の振動計算構造",
                    "saved optimized geometry": "保存された DFT 最適化構造",
                    "declared input; not a certified minimum": "宣言入力（認証済み極小ではない）",
                    "declared input; not separated gas monomers": "G の分離基準とは別の入口配置",
                    "symmetry_image_of_qrc_plus (not an independent QM endpoint)":
                    "QRC+ と等長な対称像。独立 QM 端点ではない。"}.get(observed, observed)
        ax.text2D(.5, -.015, wrapped(observed, 37), transform=ax.transAxes, fontsize=8.5,
                  color=COLORS["muted"], ha="center", va="top")


def structures_panel(fig, grid, example: dict) -> None:
    selected = selected_structures(example)
    if not selected:
        return
    reference = next((s for s in selected if s["role"] == "ts"), selected[0])
    endpoints = [s for s in selected if s["role"] in {"reactant", "product"}]
    changed = (depiction_bonds(endpoints[0]) ^ depiction_bonds(endpoints[-1])) if len(endpoints) == 2 else set()
    extent = max(max(np.max(np.abs(aligned_coords(s, reference))) for s in selected) + .5, 1.4)
    for i, structure in enumerate(selected):
        ax = fig.add_subplot(grid[1, i], projection="3d")
        draw_molecule(ax, structure, reference, changed, extent)
    if len(selected) < 3:
        ax = fig.add_subplot(grid[1, len(selected)])
        ax.set_axis_off()
        ax.text(.08, .55, "認証した TS 構造はない。" if example["id"] == "bh3_nh3_association"
                else "この役割の構造は保存されていない。", transform=ax.transAxes,
                fontsize=11, color=COLORS["muted"])


def conditions_text(free: dict | None) -> str:
    if not free:
        return "自由エネルギーを採用していない例題"
    conditions = free.get("conditions", {})
    if isinstance(conditions, str):
        return conditions
    return f"気相 / {conditions.get('T_K', '未記録')} K / " + str(conditions.get("standard_state", "未記録"))


def save(fig, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.png", dpi=300)
    svg = output / f"{stem}.svg"
    fig.savefig(svg)
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())
                   + "\n", encoding="utf-8")
    plt.close(fig)


def example_figure(example: dict, output: Path) -> None:
    index, title, purpose = TITLES.get(example["id"], ("", example["title"], ""))
    fig = plt.figure(figsize=(15.6, 9.6))
    grid = fig.add_gridspec(2, 3, left=.055, right=.975, bottom=.125, top=.815,
                           height_ratios=[1, 1.08], wspace=.30, hspace=.54)
    fig.text(.055, .95, f"{index}  {title}", fontsize=22, fontweight="bold")
    fig.text(.055, .905, purpose, fontsize=13.5, color=COLORS["muted"])
    free = example.get("free_energy")
    method = (free or {}).get("method", "")
    profiles = example.get("electronic_profiles") or []
    methods = ["G: " + str(method) + " + 熱補正 [GoodVibes]"] if method else []
    if free and free.get("geometry_method"):
        methods += ["構造・振動: " + free["geometry_method"]]
    if profiles:
        methods += ["走査 E: " + str(profiles[0].get("method", ""))]
    fig.text(.055, .865, wrapped(" | ".join(methods), 130), fontsize=9.5,
             color=COLORS["muted"], va="top")
    energy_diagram(fig.add_subplot(grid[0, 0]), example)
    electronic_profile(fig.add_subplot(grid[0, 1]), example)
    decision_panel(fig.add_subplot(grid[0, 2]), example)
    structures_panel(fig, grid, example)
    elements = sorted({a["element"] for s in selected_structures(example) for a in s["atoms"]})
    handles = [Line2D([], [], marker="o", ls="", ms=7,
                      markerfacecolor=ATOM_COLORS.get(e, "#A8A3BB"),
                      markeredgecolor="#657583", label=e) for e in elements]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(.975, .093),
               ncol=len(elements), frameon=False, fontsize=9, handletextpad=.4, columnspacing=1)
    structure_title = "D  保存された計算座標の 3D 表示"
    if example["id"] == "nh3_hf_exchange":
        structure_title = "D  保存座標と明示した対称像の 3D 表示"
    fig.text(.055, .077, structure_title, fontsize=10.5, fontweight="bold")
    fig.text(.055, .052, "線は距離基準による描画。TS の橙破線は端点間で変化する接触。"
             "番号は XYZ 順序（1 始まり、重い原子・反応中心を表示）。", fontsize=8.8,
             color=COLORS["muted"])
    runs = example.get("provenance", {}).get("runs", [])
    run_labels = [r["run_path"].replace("/home/user/hfauto_r10/", "") for r in runs]
    fig.text(.055, .030, "計算記録: " + " / ".join(run_labels), fontsize=8.5,
             color=COLORS["muted"])
    fig.text(.055, .009, conditions_text(free) + "  |  原データ: data/data.json・data/xyz/"
             "  |  図: references.json [matplotlib] / 計算法の文献は資料に記載",
             fontsize=8.8, color=COLORS["muted"])
    save(fig, output, example["id"])


def rounded_box(ax, xy, size, title: str, body: str, color: str) -> None:
    x, y = xy
    width, height = size
    ax.add_patch(FancyBboxPatch((x, y), width, height,
                 boxstyle="round,pad=.012,rounding_size=.025", fc="#F2F6F9", ec=color, lw=1.5))
    ax.text(x + .025, y + height - .035, title, va="top", fontsize=12,
            fontweight="bold", color=color)
    ax.text(x + .025, y + height - .1, body, va="top", fontsize=10.2, linespacing=1.7)


def arrow(ax, start, end, color=COLORS["muted"]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                               color=color, linewidth=1.5))


def decision_chart(examples: list[dict], output: Path) -> None:
    """Explain the actual distinctions without pretending this is a kinetic network."""
    fig, ax = plt.subplots(figsize=(15.6, 8.8))
    fig.subplots_adjust(left=.035, right=.985, bottom=.08, top=.86)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.text(.055, .935, "反応候補から、何を根拠に判断するか", fontsize=22, fontweight="bold")
    fig.text(.055, .885, "終了した計算を観測に分け、停留点・接続・比較基準を確認する", fontsize=13,
             color=COLORS["muted"])
    rounded_box(ax, (.025, .68), (.255, .265), "1  候補を作る / 与える",
                "discover: 結合変更から仮説を生成\nknown_endpoints: 両側の候補を宣言\n候補数だけでは反応の証明にならない", COLORS["reactant"])
    rounded_box(ax, (.365, .68), (.255, .265), "2  停留性を確認",
                "極小: 不適切な虚振動がないか\nTS: 虚振動と残る勾配を確認\n電子状態・計算手法の整合を確認", COLORS["ts"])
    rounded_box(ax, (.705, .68), (.255, .265), "3  経路の観測を集める",
                "TS: 虚モードと QRC の両側接続\n会合: 拘束走査の形状\n端点: 原子対応と極小の照合", COLORS["product"])
    arrow(ax, (.287, .81), (.355, .81))
    arrow(ax, (.627, .81), (.695, .81))
    branches = [
        (.025, "同じ極小", "TMA・(HF)₂\n別名の端点が同じ basin に収束\n→ 新しい反応の障壁は描かない", COLORS["muted"]),
        (.365, "TS を持つ経路", "NH₃ 反転 / NH₃・HF / OH＋CH₄\n停留性と両側接続を満たす\n→ 同じ基準の G で比較", COLORS["pass"]),
        (.705, "TS を採用しない会合", "BH₃＋NH₃\n保存走査では解像度内で障壁なし\n→ 捕獲律速を別扱い", COLORS["hold"]),
    ]
    for x, title, body, color in branches:
        rounded_box(ax, (x, .245), (.255, .28), title, body, color)
    arrow(ax, (.83, .668), (.83, .555))
    arrow(ax, (.81, .595), (.49, .595))
    arrow(ax, (.49, .595), (.49, .555))
    arrow(ax, (.49, .595), (.15, .595))
    arrow(ax, (.15, .595), (.15, .555))
    ax.text(.025, .145, "比較の条件", fontsize=11, fontweight="bold")
    ax.text(.18, .145, "理論・温度・標準状態・分離 / 錯体の基準をそろえる。未達の点は順位に入れない。",
            fontsize=11)
    ax.text(.025, .075, "次の開発段階", fontsize=11, fontweight="bold")
    ax.text(.18, .075, "異なる出発状態の δG_eff の順位は、初期系の主要生成物・収率・反応速度ではない。",
            fontsize=11)
    fig.text(.055, .033, "対応する計算例と観測値: data/data.json  |  現行仕様: docs/design.md  |  "
             "将来の反応網・速度評価: docs/improvement-plan.md", fontsize=9, color=COLORS["muted"])
    save(fig, output, "reaction_decision_chart")


def overview(examples: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(15.6, 8.4))
    fig.subplots_adjust(left=.045, right=.98, bottom=.12, top=.8)
    ax.set_axis_off()
    fig.text(.055, .935, "5 つの例題で読む、現在の hfauto", fontsize=22, fontweight="bold")
    fig.text(.055, .87, "探索・認証・比較の異なる観測を並べる。例題間の障壁の大小で主要生成物を決めない。",
             fontsize=12.5, color=COLORS["muted"])
    rows = []
    for example in examples:
        index, title, purpose = TITLES.get(example["id"], ("", example["title"], ""))
        lines = decision_lines(example)
        free = example.get("free_energy")
        value = (free or {}).get("dG_eff_kcal")
        energy = f"δG_eff {value:.2f} kcal mol⁻¹" if numeric(value) is not None else "TS 障壁の順位なし"
        rows.append([index + "  " + title, purpose, lines[0][1] + "\n" + energy])
    table = ax.table(cellText=rows, colLabels=["計算例", "説明できる処理", "観測した判断"],
                     colWidths=[.29, .40, .31], cellLoc="left", colLoc="left", bbox=(0, 0, 1, 1))
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor(COLORS["grid"])
        cell.PAD = .10
        cell.set_facecolor("#17334C" if row == 0 else ("#F1F6F9" if row % 2 else "white"))
        cell.get_text().set_color("white" if row == 0 else COLORS["ink"])
        if row == 0:
            cell.get_text().set_fontweight("bold")
    fig.text(.055, .055, "実測構造・理論・基準・判定ログ・出典 run は各例題の図と data/data.json に記録。"
             "過去の計算を現在の処理の説明に用い、コード全体の科学的検証とは区別する。",
             fontsize=9.5, color=COLORS["muted"])
    save(fig, output, "examples_overview")


def coverage_outcomes(example: dict) -> dict:
    coverage = example.get("coverage", {})
    counts = coverage.get("counts", {})
    if isinstance(counts, list):
        counts = next((row for row in counts if row.get("mechanism") == "nt2"), {})
    outcomes = coverage.get("outcomes", counts.get("outcomes", counts))
    outcomes = {key.removesuffix("s") if key in {"negatives", "products"} else key: value
                for key, value in outcomes.items()}
    required = {"negative", "product", "unconnected", "not_attempted"}
    return outcomes if required.issubset(outcomes) else {}


def exploration_coverage(examples: list[dict], output: Path) -> bool:
    """Count report-time discovery outcomes, including DFT admission decisions."""
    selected = [e for e in examples if e["id"] in {"nh3_hf_exchange", "tma_hf2_same_basin"}
                and coverage_outcomes(e)]
    if not selected:
        return False
    fig = plt.figure(figsize=(15.6, 7.5))
    grid = fig.add_gridspec(1, 2, left=.14, right=.97, bottom=.31, top=.69,
                          width_ratios=[2.4, 1], wspace=.25)
    ax = fig.add_subplot(grid[0, 0])
    fig.text(.055, .935, "探索クエリの内訳から、被覆を読む", fontsize=22, fontweight="bold")
    fig.text(.055, .872, "保存された分類、未試行、report 時点に残った候補を区別する", fontsize=13,
             color=COLORS["muted"])
    order = ["negative", "not_attempted", "unconnected", "product"]
    labels = {"negative": "negative（探索・入口で採用外）", "not_attempted": "未試行",
              "unconnected": "unconnected（接続未確認）", "product": "product（保存候補辺）"}
    colors = {"negative": "#91A7B7", "not_attempted": "#E1E8ED",
              "unconnected": COLORS["ts"], "product": COLORS["pass"]}
    totals, rows = [], []
    for y, example in zip(range(len(selected) - 1, -1, -1), selected, strict=True):
        counts = coverage_outcomes(example)
        total = sum(counts[k] for k in order)
        totals.append(total)
        left = 0
        for key in order:
            value = counts[key]
            ax.barh(y, value, left=left, height=.38, color=colors[key],
                    label=labels[key] if y == len(selected) - 1 else None)
            if value / total > .045:
                ax.text(left + value / 2, y, f"{value:,}", ha="center", va="center",
                        fontsize=11, color="white" if key == "negative" else COLORS["ink"])
            left += value
        ax.text(total + 45, y, f"保存分類 {total:,}", va="center", fontsize=10.5)
        ax.text(0, y - .29, f"試行済み分類 {total - counts['not_attempted']:,} / 未試行 "
                f"{counts['not_attempted']:,}", fontsize=9.5, color=COLORS["muted"])
        rows.append(["S10 pilot 全体\nNH₃・HF / TMA・HF" if example["id"] == "nh3_hf_exchange"
                     else "TMA・(HF)₂",
                     str(counts["product"]), str(counts["unconnected"])])
    ax.set_xlim(0, max(totals) * 1.19)
    ax.set_ylim(-.6, len(selected) - .5)
    ax.set_yticks(list(range(len(selected) - 1, -1, -1)), [row[0] for row in rows])
    ax.set_xlabel("保存された探索クエリの分類件数（計算失敗は別集計）", labelpad=15)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=15)
    ax.grid(axis="x", color=COLORS["grid"], lw=.7)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(-.13, 1.31), ncol=2, frameon=False, fontsize=10)
    right = fig.add_subplot(grid[0, 1])
    right.set_axis_off()
    right.text(0, 1.02, "小さな区分の実数", fontsize=12, fontweight="bold")
    table = right.table(cellText=rows, colLabels=["系", "product", "unconnected"],
                        colWidths=[.44, .25, .31], cellLoc="center", bbox=(0, .3, 1, .58))
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor(COLORS["grid"])
        cell.set_facecolor("#F1F6F9" if row == 0 else "white")
        cell.get_text().set_color(COLORS["ink"])
    right.text(0, .13, "product は保存された候補辺。\nDFT で認証した TS 数とは別。", fontsize=9.5,
               color=COLORS["muted"], linespacing=1.7)
    fig.text(.055, .185, "S10 の集計は NH₃・HF と TMA・HF を含む pilot 全体。"
             "本資料で紹介する認証経路は NH₃・HF。", fontsize=10.5)
    fig.text(.055, .135, "report 時点に保存された DiscoveryRecord の分類。DFT 入口の選抜結果を含む。"
             "計算失敗は別集計。", fontsize=10.5)
    fig.text(.055, .087, "同じ negative の理由や failure_kind は別の集計。"
             "この積上げ棒へ加えない。数の大小だけで探索の汎用性・化学的網羅性を保証しない。",
             fontsize=9.5, color=COLORS["muted"])
    fig.text(.055, .035, "原データ: data/data.json coverage / 各 run の report/coverage.csv  "
             "|  候補探索の手法: references.json [nt2]  |  作図: [matplotlib]", fontsize=9,
             color=COLORS["muted"])
    save(fig, output, "exploration_coverage")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "data" / "data.json")
    parser.add_argument("--output", type=Path, default=HERE / "figures")
    args = parser.parse_args()
    document = json.loads(args.data.read_text(encoding="utf-8"))
    configure()
    for example in document["examples"]:
        example_figure(example, args.output)
    decision_chart(document["examples"], args.output)
    overview(document["examples"], args.output)
    has_coverage = exploration_coverage(document["examples"], args.output)
    print(f"Rendered {len(document['examples']) + 2 + has_coverage} figures "
          "as PNG (300 dpi) and SVG.")


if __name__ == "__main__":
    main()
