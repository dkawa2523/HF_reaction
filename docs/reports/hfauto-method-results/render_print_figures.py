"""Render compact, legible report figures from the same saved example data.

The landscape example figures remain untouched. This print composition reuses
their energy points, structure depiction and decision definitions at A4 text width.
"""

from __future__ import annotations

import argparse
import json
import runpy
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.text import Text
from mpl_toolkits.mplot3d import proj3d

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE = runpy.run_path(str(REPO / "examples" / "visualization" / "render_figures.py"))
INK = BASE["COLORS"]["ink"]
MUTED = BASE["COLORS"]["muted"]
WIDTH, HEIGHT = 6.7, 6.2


def wrap_print(value: str, max_units: int) -> str:
    """Wrap full-width Japanese and half-width Latin text at print column width."""
    lines, current, units = [], "", 0
    for char in value:
        size = 2 if unicodedata.east_asian_width(char) in {"F", "W"} else 1
        if char == "\n" or units + size > max_units:
            lines.append(current.strip())
            current, units = "", 0
        if char != "\n":
            current += char
            units += size
    if current:
        lines.append(current.strip())
    return "\n".join(lines)


def compact_axis(ax) -> None:
    """Remove duplicated plot captions and keep all point/tick text at 9 pt."""
    for artist in list(ax.texts):
        value = artist.get_text()
        if value.startswith(("破線は", "点間の線は")):
            artist.remove()
        elif value.startswith(("G = 0:", "δG_eff =")):
            artist.remove()
        elif value == "E = 0: 最終 backward IRC 点":
            artist.set_text("E = 0: IRC 左端")
    for artist in ax.findobj(match=Text):
        artist.set_fontsize(9)
    ax.title.set_fontsize(9.5)
    ax.tick_params(axis="both", pad=2)
    ax.xaxis.labelpad = 4
    ax.yaxis.labelpad = 4


def concise_decisions(example: dict) -> list[str]:
    lines = []
    free = example.get("free_energy") or {}
    if free.get("dG_eff_kcal") is not None:
        lines.append(f"δG_eff: {free['dG_eff_kcal']:.2f} kcal mol⁻¹")
    for heading, value in BASE["decision_lines"](example):
        if heading in {"根拠", "化学状態"}:
            continue
        heading = {"反応自由エネルギー": "ΔG_R→P",
                   "反応自由エネルギー（錯体間）": "ΔG_R→P（錯体間）",
                   "比較基準": "δG_eff 基準"}.get(heading, heading)
        value = value.replace("同じ化学状態を結ぶ TS を認証", "同じ化学状態の TS を認証")
        value = value.replace("異なる化学状態を結ぶ TS を認証", "別の化学状態の TS を認証")
        lines.append(f"{heading}: {value}")
    return lines


def profile_evidence(example: dict) -> str:
    free = example.get("free_energy") or {}
    if example["id"] == "nh3_hf_exchange":
        return ("同じ化学状態の TS を認証。両側から同じ登録極小へ接続。\n"
                f"虚振動 {example['decisions']['imag_cm1']:.1f} cm⁻¹；"
                f"δG_eff {free['dG_eff_kcal']:.2f} kcal mol⁻¹（R 錯体基準）。")
    if example["id"] == "tma_hf2_same_basin":
        return "宣言端点が同じ DFT 極小へ収束。認証した TS・G 順位はない。"
    return ("拘束走査の分解能内で障壁なし。TS 障壁の順位から除外（捕獲律速）。\n"
            f"ΔG_R→P {free['dG_rxn_kcal']:+.2f} kcal mol⁻¹（分離した単量体基準）。")


def evidence_axis(ax, example: dict) -> None:
    ax.set_axis_off()
    ax.set_title("C  判断の根拠", loc="left", fontweight="bold", fontsize=9.5, pad=10)
    text = "\n".join(wrap_print(line, 36) for line in concise_decisions(example))
    ax.text(0, .98, text, transform=ax.transAxes, va="top", fontsize=9, linespacing=1.48)


def compact_atom_labels(ax, position: tuple, structure: dict, reference: dict,
                        changed: set) -> None:
    """Place 9 pt index labels after the 3D projection has its final axes size."""
    ax.figure.canvas.draw()
    xyz = BASE["aligned_coords"](structure, reference)
    projected = np.asarray(proj3d.proj_transform(*xyz.T, ax.get_proj())[:2]).T
    display = ax.transData.transform(projected)
    scale = ax.figure.dpi / 72
    panel = ax.figure.transFigure.transform([
        [position[0], position[1]], [position[0] + position[2], position[1] + position[3]]])
    used = []
    offsets = [(radius * np.cos(angle), radius * np.sin(angle))
               for radius in [16, 23, 30, 37]
               for angle in np.linspace(0, 2 * np.pi, 16, endpoint=False)]
    for i, atom in enumerate(structure["atoms"]):
        if len(structure["atoms"]) > 10 and atom["element"] == "H" \
                and not any(i in pair for pair in changed):
            continue
        label = f"{atom['element']}{atom['index'] + 1}"
        width, height = len(label) * 5.7 * scale, 14 * scale
        candidates = []
        for dx, dy in offsets:
            center = display[i] + np.array([dx, dy]) * scale
            bounds = (center[0] - width / 2, center[1] - height / 2,
                      center[0] + width / 2, center[1] + height / 2)
            collisions = sum(not (bounds[2] + 4 * scale < b[0] or bounds[0] - 4 * scale > b[2]
                                  or bounds[3] + 4 * scale < b[1] or bounds[1] - 4 * scale > b[3])
                             for b in used)
            outside = sum([bounds[0] < panel[0, 0], bounds[1] < panel[0, 1],
                           bounds[2] > panel[1, 0], bounds[3] > panel[1, 1]])
            nearest = np.min(np.linalg.norm(display - center, axis=1)) / scale
            atom_overlap = max(0, 12 - nearest)
            distance = np.hypot(dx, dy)
            candidates.append((outside * 10000 + collisions * 1000 + atom_overlap * 20 + distance,
                               dx, dy, bounds))
        _, dx, dy, bounds = min(candidates)
        used.append(bounds)
        ax.annotate(label, projected[i], xytext=(dx, dy), textcoords="offset points",
                    ha="center", va="center", fontsize=9, color=INK,
                    arrowprops={"arrowstyle": "-", "color": "#98A6B1", "lw": .55}, zorder=100)


def small_structure(fig, position, structure: dict, reference: dict, changed: set,
                    extent: float, example_id: str) -> None:
    ax = fig.add_axes(position, projection="3d")
    BASE["draw_molecule"](ax, structure, reference, changed, extent)
    ax.collections[-1].set_sizes(np.array([85 if a["element"] == "H" else 190
                                          for a in structure["atoms"]]))
    for line in ax.lines:
        line.set_linewidth(1.0 if line.get_linestyle() == "--" else 1.6)
    for text in list(ax.texts):
        text.remove()
    titles = {"reactant": "R 側の保存極小", "ts": "認証した TS", "product": "P 側の保存極小",
              "declared_reactant": "R 入力: neutral", "declared_product": "P 入力: shared_proton",
              "minimum": "同じ DFT 極小"}
    if example_id == "nh3_hf_exchange":
        titles.update({"reactant": "QRC+ 保存端点",
                       "product": "QRC+ の対称像\n（独立 QM ではない）"})
    if example_id == "bh3_nh3_association":
        titles.update({"reactant": "入口配置", "product": "最適化した付加体"})
    ax.set_title(titles[structure["role"]], fontsize=9, fontweight="bold", pad=-4)
    compact_atom_labels(ax, position, structure, reference, changed)


def print_example(example: dict, output: Path) -> None:
    index, title, _ = BASE["TITLES"][example["id"]]
    fig = plt.figure(figsize=(WIDTH, HEIGHT))
    fig.text(.055, .975, f"{index}  {title}", va="top", fontsize=13, fontweight="bold")
    free = example.get("free_energy")
    profiles = example.get("electronic_profiles") or []
    method_parts = []
    if free:
        method_parts.append("G: " + free["method"] + " + 熱補正 (GoodVibes)")
        method_parts.append("構造・振動: " + free["geometry_method"] + " / 298.15 K, 1 atm")
    if profiles:
        method_parts.append("走査 E: " + profiles[0]["method"])
    fig.text(.055, .918, "\n".join(method_parts), va="top", fontsize=9,
             color=MUTED, linespacing=1.25)
    left = fig.add_axes((.095, .525, .36, .245))
    BASE["energy_diagram"](left, example)
    compact_axis(left)
    if free:
        labels = ["R_sep" if p["role"] == "reactant_separated" else "TS" if p["role"] == "ts"
                  else "R" if p["role"] == "reactant" else "P" for p in free["points"]]
        left.set_xticks(range(len(labels)), labels, fontsize=9)
        zero = "分離単量体=0" if example["id"] == "bh3_nh3_association" \
            else "R 錯体=0" if example["id"] in {"nh3_hf_exchange", "oh_ch4_abstraction"} \
            else "R 極小=0"
        left.set_title("A  離散 G（" + zero + "）", loc="left", fontsize=9.5,
                       fontweight="bold", pad=10)
    if not free:
        left.clear()
        left.set_axis_off()
        left.set_title("A  端点の判断", loc="left", fontsize=9.5, fontweight="bold", pad=10)
        left.text(.04, .8, "端点は同じ極小へ合流。\n認証した TS はない。\nG・TS 障壁の順位は付けない。",
                  transform=left.transAxes, va="top", fontsize=9, linespacing=1.7)
    right = fig.add_axes((.585, .525, .375, .245))
    if profiles:
        BASE["electronic_profile"](right, example)
        compact_axis(right)
        right.set_title("B  保存計算点の E" + ("（低レベル・未認証）"
                        if example["id"] == "tma_hf2_same_basin" else ""),
                        loc="left", fontsize=9.5, fontweight="bold", pad=10)
        coordinate = profiles[0]["coordinate_kind"]
        right.set_xlabel("保存 IRC 反復数（符号付き）" if "IRC iteration" in coordinate
                         else "B–N 距離 / Å", fontsize=9)
        fig.text(.055, .413, "C  " + profile_evidence(example), fontsize=9,
                 va="top", linespacing=1.4)
    else:
        evidence_axis(right, example)
    selected = BASE["selected_structures"](example)
    reference = next((s for s in selected if s["role"] == "ts"), selected[0])
    endpoints = [s for s in selected if s["role"] in {"reactant", "product"}]
    changed = (BASE["depiction_bonds"](endpoints[0]) ^ BASE["depiction_bonds"](endpoints[-1])) \
        if len(endpoints) == 2 else set()
    extent = max(max(np.max(np.abs(BASE["aligned_coords"](s, reference)))
                     for s in selected) + .5, 1.4)
    fig.text(.055, .330, "D  保存構造と原子対応", fontsize=9.5, fontweight="bold")
    for i, structure in enumerate(selected):
        small_structure(fig, (.025 + i * .325, .085, .31, .190), structure, reference,
                        changed, extent, example["id"])
    if len(selected) < 3:
        fig.text(.705, .203, "認証した TS\n構造はない。", fontsize=9, color=MUTED,
                 ha="center", linespacing=1.5)
    foot = "G は停留点の離散値。E は保存計算点の結線（G と別の量）。" if profiles \
        else "G は停留点の離散値。点間の線は連続した経路の曲線ではない。"
    if not free:
        foot = "低レベルの候補経路 E を示す。認証した TS 障壁・G 順位はない。"
    fig.text(.055, .041, foot, fontsize=9, color=MUTED)
    fig.text(.055, .011, "構造は距離に基づく描画。原子番号は XYZ 順序（1 始まり）。", fontsize=9,
             color=MUTED)
    BASE["save"](fig, output, example["id"])


def print_overview(examples: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(WIDTH, HEIGHT))
    fig.subplots_adjust(left=.045, right=.965, bottom=.15, top=.85)
    ax.set_axis_off()
    fig.text(.055, .975, "5 つの例題で読む、現在の hfauto", fontsize=13,
             fontweight="bold", va="top")
    fig.text(.055, .915, "探索・認証・比較の異なる観測を並べる。", fontsize=9, color=MUTED)
    rows = []
    for example in examples:
        index, title, purpose = BASE["TITLES"][example["id"]]
        outcome = BASE["decision_lines"](example)[0][1]
        free = example.get("free_energy") or {}
        value = free.get("dG_eff_kcal")
        energy = f"δG_eff {value:.2f} kcal mol⁻¹" if value is not None else "TS 障壁の順位なし"
        rows.append([wrap_print(index + " " + title, 20), wrap_print(purpose, 30),
                     wrap_print(outcome, 28) + "\n" + energy])
    table = ax.table(cellText=rows, colLabels=["計算例", "説明できる処理", "観測した判断"],
                     colWidths=[.26, .38, .36], cellLoc="left", colLoc="left", bbox=(0, 0, 1, 1))
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, _), cell in table.get_celld().items():
        cell.PAD = .06
        cell.set_edgecolor(BASE["COLORS"]["grid"])
        cell.set_facecolor(INK if row == 0 else "#F1F6F9" if row % 2 else "white")
        cell.get_text().set_color("white" if row == 0 else INK)
        if row == 0:
            cell.get_text().set_fontweight("bold")
    note = "例題間の δG_eff の大小から、初期系の主要生成物・収率を決めない。\n数値・座標・条件・判定ログは同じ原データに保存。"
    fig.text(.055, .095, note, fontsize=9, color=MUTED, va="top", linespacing=1.5)
    BASE["save"](fig, output, "examples_overview")


def print_decision_chart(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(WIDTH, HEIGHT))
    fig.subplots_adjust(left=.03, right=.985, bottom=.07, top=.86)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.text(.055, .975, "反応候補から、何を根拠に判断するか", fontsize=13,
             fontweight="bold", va="top")
    fig.text(.055, .918, "候補・停留点・接続・比較基準を分けて確認する。", fontsize=9, color=MUTED)
    items = [
        (.025, .62, "1 候補", "discover: 結合変更\nknown_endpoints:\n両端を宣言\n候補だけで証明しない", "reactant"),
        (.36, .62, "2 停留性", "極小と TS の振動\n残る勾配を確認\n計算手法と\n電子状態を確認", "ts"),
        (.695, .62, "3 経路の観測", "TS: 虚モード / QRC\n会合: 拘束走査\n端点: 極小の照合\n原子対応を保持", "product"),
        (.025, .17, "同じ極小", "TMA・(HF)₂\n宣言端点が合流\n→ 新しい TS 障壁\nを描かない", "muted"),
        (.36, .17, "TS を認証", "NH₃ 反転 / NH₃・HF\nOH＋CH₄\n停留性と接続を確認\n→ 同じ条件の G", "pass"),
        (.695, .17, "TS を採用しない", "BH₃＋NH₃\n走査の分解能内で\n障壁なし\n→ 捕獲律速を別扱い", "hold"),
    ]
    for x, y, title, body, color in items:
        BASE["rounded_box"](ax, (x, y), (.275, .31), title, body, BASE["COLORS"][color])
    for start, end in [((.303, .78), (.35, .78)), ((.638, .78), (.685, .78)),
                       ((.835, .61), (.835, .51)), ((.835, .51), (.16, .51)),
                       ((.16, .51), (.16, .49)), ((.5, .51), (.5, .49)),
                       ((.835, .51), (.835, .49))]:
        BASE["arrow"](ax, start, end)
    for text in ax.findobj(match=Text):
        text.set_fontsize(9)
    ax.text(.025, .07, "理論・温度・標準状態・分離/錯体の基準をそろえる。", fontsize=9)
    fig.text(.055, .04, "計算終了・科学的認証・生成物比の予測は、それぞれ異なる判断。", fontsize=9,
             color=MUTED)
    BASE["save"](fig, output, "reaction_decision_chart")


def print_coverage(examples: list[dict], output: Path) -> None:
    selected = []
    for example in examples:
        if example["id"] not in {"nh3_hf_exchange", "tma_hf2_same_basin", "oh_ch4_abstraction"}:
            continue
        for record in example["coverage"]["counts"]:
            if record["mechanism"] not in {"nt2", "relaxation"}:
                continue
            label = "S10 pilot 全体\nNH₃・HF / TMA・HF" if example["id"] == "nh3_hf_exchange" \
                else "TMA・(HF)₂" if example["id"] == "tma_hf2_same_basin" \
                else "S6 OH＋CH₄\n" + record["mechanism"]
            outcomes = {key.removesuffix("s") if key in {"negatives", "products"} else key: value
                        for key, value in record["outcomes"].items()}
            selected.append((label, outcomes))
    fig = plt.figure(figsize=(WIDTH, 7.2))
    fig.text(.055, .975, "探索クエリの内訳から、被覆を読む", fontsize=13,
             fontweight="bold", va="top")
    fig.text(.055, .919, "保存された分類・未試行・report 時点の候補を区別する。", fontsize=9, color=MUTED)
    ax = fig.add_axes((.27, .42, .69, .38))
    order = ["negative", "not_attempted", "unconnected", "product"]
    labels = {"negative": "negative（探索・入口で採用外）", "not_attempted": "未試行",
              "unconnected": "unconnected（接続未確認）", "product": "product（保存候補辺）"}
    colors = {"negative": "#91A7B7", "not_attempted": "#E1E8ED",
              "unconnected": BASE["COLORS"]["ts"], "product": BASE["COLORS"]["pass"]}
    rows, totals = [], []
    largest_total = max(sum(counts.values()) for _, counts in selected)
    for y, (label, outcomes) in zip(range(len(selected) - 1, -1, -1), selected, strict=True):
        counts = {key: outcomes.get(key, 0) for key in order}
        total = sum(counts[key] for key in order)
        totals.append(total)
        left = 0
        for key in order:
            value = counts[key]
            ax.barh(y, value, left=left, height=.40, color=colors[key],
                    label=labels[key] if y == len(selected) - 1 else None)
            if value > largest_total * .08:
                ax.text(left + value / 2, y, f"{value:,}", fontsize=9, ha="center", va="center",
                        color="white" if key == "negative" else INK)
            left += value
        ax.text(total + 50, y, f"保存分類 {total:,}", fontsize=9, va="center")
        ax.text(0, y - .40, f"試行済み分類 {total - counts['not_attempted']:,} / "
                f"未試行 {counts['not_attempted']:,}", fontsize=9, color=MUTED, va="center")
        rows.append([label.replace("\nNH₃・HF / TMA・HF", "").replace("\n", " / "),
                     str(counts["negative"]), str(counts["product"]),
                     str(counts["unconnected"]), str(counts["not_attempted"])])
    ax.set_xlim(0, max(totals) * 1.40)
    ax.set_ylim(-.6, len(selected) - .4)
    ax.set_yticks(range(len(selected) - 1, -1, -1), [label for label, _ in selected], fontsize=9)
    ax.set_xlabel("保存分類件数（計算失敗は別集計）", fontsize=9, labelpad=5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.grid(axis="x", color=BASE["COLORS"]["grid"], lw=.6)
    ax.set_axisbelow(True)
    handles, legend_labels = ax.get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper left", bbox_to_anchor=(.04, .886),
               ncol=2, frameon=False, fontsize=9, handlelength=1.3, columnspacing=1.1)
    table_axis = fig.add_axes((.055, .17, .905, .18))
    table_axis.set_axis_off()
    table = table_axis.table(cellText=rows, colLabels=["保存分類の内訳", "採用外", "product", "未接続", "未試行"],
                            colWidths=[.32, .18, .15, .20, .15], cellLoc="center", bbox=(0, 0, 1, 1))
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor(BASE["COLORS"]["grid"])
        cell.set_facecolor("#F1F6F9" if row == 0 else "white")
        cell.get_text().set_color(INK)
    notes = ["report 時点に保存された DiscoveryRecord の分類。DFT 入口の選抜結果を含む。",
             "計算失敗は別集計。product は認証した TS の数とは別。",
             "S10 は NH₃・HF と TMA・HF を含む pilot 全体。紹介する認証経路は NH₃・HF。",
             "同じ negative の理由や failure_kind は別集計。この積上げに加えない。"]
    fig.text(.055, .135, "\n".join(wrap_print(note, 88) for note in notes), fontsize=9,
             color=MUTED, va="top", linespacing=1.5)
    BASE["save"](fig, output, "exploration_coverage")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=REPO / "examples" / "visualization" / "data" / "data.json")
    parser.add_argument("--output", type=Path, default=HERE / "figures")
    args = parser.parse_args()
    BASE["configure"]()
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5,
                         "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9})
    document = json.loads(args.data.read_text(encoding="utf-8"))
    for example in document["examples"]:
        print_example(example, args.output)
    print_overview(document["examples"], args.output)
    print_decision_chart(args.output)
    print_coverage(document["examples"], args.output)
    print(f"Rendered {len(document['examples']) + 3} report figures at {WIDTH} inches wide "
          f"({HEIGHT} inches high; coverage 7.2 inches high).")


if __name__ == "__main__":
    main()
