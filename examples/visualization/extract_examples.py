"""Extract a small, auditable illustration bundle from completed, read-only saved runs.

WSL: /home/user/.venvs/hfauto-prod/bin/python examples/visualization/extract_examples.py
No jobs, reports, thermal corrections or geometries are recomputed. The one derived
geometry is the exact symmetry image that the recorded NH3-HF QRC used instead of a
second optimization; its transformation is reproduced with the existing identity code.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from hfauto.chemistry.identity import IMAGE_A, carry
from hfauto.core.constants import HARTREE_TO_KCAL_MOL


STAGES = ("structures", "conformers", "screen", "explore", "dft", "paths", "sp", "thermo")
SPECS = (
    ("nh3_inversion", "NH3の反転", "NH3 inversion", "W5/fresh/nh3_planar_seed", None,
     "nh3_inversion"),
    ("nh3_hf_exchange", "NH3・HFの水素交換", "NH3 + HF", "W4/fresh/s10_amine_pilot2",
     "W5/thermo_replay/s10_amine_pilot2", "rxn_discovery_db5f945102"),
    ("tma_hf2_same_basin", "TMA・(HF)2の極小への合流", "(CH3)3N + 2 HF",
     "W4/fresh/s19_tma_hf2", None, "neutral_to_shared_proton"),
    ("bh3_nh3_association", "BH3 + NH3の会合", "BH3 + NH3", "W3/fresh/bh3_nh3", None,
     "bh3_nh3_association"),
    ("oh_ch4_abstraction", "OH + CH4の水素引き抜き", "OH + CH4", "W4/fresh/s6_oh_ch4",
     "W5/thermo_replay/s6_oh_ch4", "rxn_discovery_6494d6a583"),
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def frames(text: str) -> list[tuple[str, list[str], np.ndarray]]:
    """Parse saved XYZ frames without introducing coordinates or changing atom order."""
    lines, out, i = text.splitlines(), [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        comment = lines[i + 1].strip()
        rows = [row.split() for row in lines[i + 2:i + 2 + n]]
        if len(rows) != n:
            raise ValueError("truncated XYZ")
        out.append((comment, [row[0] for row in rows],
                    np.array([[float(v) for v in row[1:4]] for row in rows])))
        i += n + 2
    return out


def xyz_text(symbols: list[str], coords: np.ndarray, comment: str) -> str:
    return f"{len(symbols)}\n{comment}\n" + "".join(
        f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}\n"
        for symbol, (x, y, z) in zip(symbols, coords, strict=True))


class Run:
    def __init__(self, path: Path):
        self.path = path
        self.files: dict[str, dict] = {}
        self.used_calcs: dict[str, dict] = {}
        self.artifacts: dict[str, dict] = {}
        self.by_stage: dict[str, list[dict]] = {}
        config = self.read("resolved_config.yaml").decode()
        self.code_version = re.search(r"code_version:\s*(\S+)", config).group(1)
        self.state = json.loads(self.read("run_state.json"))
        for stage in STAGES:
            file = f"{stage}/manifest.json"
            if not (path / file).exists():
                continue
            artifacts = json.loads(self.read(file))["artifacts"]
            self.by_stage[stage] = artifacts
            self.artifacts.update((a["artifact_id"], a) for a in artifacts)

    def read(self, file: str, expected_sha: str | None = None) -> bytes:
        data = (self.path / file).read_bytes()
        digest = sha(data)
        if expected_sha and expected_sha != digest:
            raise ValueError(f"source hash mismatch: {self.path / file}")
        self.files[file] = {"path": file, "sha256": digest, "bytes": len(data)}
        return data

    def payload(self, artifact: str) -> dict:
        return self.artifacts[artifact]["payload"]

    def calc(self, calc_id: str) -> dict:
        value = self.payload(calc_id)
        self.used_calcs[calc_id] = {
            "id": calc_id, "job_key": value.get("job_key"), "task": value["task"],
            "level": value["level"], "energy_hartree": value["energy_hartree"],
            "final_file": value["final"]["file"], "s2": value.get("s2"),
        }
        return value

    def geometry(self, ref: dict) -> tuple[list[str], np.ndarray, str]:
        file = ref.get("file", ref)
        raw = self.read(file["path"], file.get("sha256"))
        _, symbols, coords = frames(raw.decode())[-1]
        return symbols, coords, file["path"]

    def provenance(self) -> dict:
        return {"run_path": str(self.path), "code_version": self.code_version,
                "source_files": sorted(self.files.values(), key=lambda row: row["path"]),
                "calculations": list(self.used_calcs.values()), "stages": self.state}


def save_structure(out: Path, run: Run, example: str, role: str, label: str, ref: dict,
                   observed_role: str = "saved optimized geometry", **extra) -> dict:
    symbols, coords, source = run.geometry(ref)
    file = f"xyz/{example}/{role}.xyz"
    dest = out / file
    dest.parent.mkdir(parents=True, exist_ok=True)
    raw = run.read(source)
    dest.write_bytes(raw)
    return {"role": role, "label": label, "xyz_file": file,
            "atoms": [{"index": i, "element": s, "x": float(x), "y": float(y), "z": float(z)}
                      for i, (s, (x, y, z)) in enumerate(zip(symbols, coords, strict=True))],
            "source_run": str(run.path), "source_path": source, "sha256": sha(raw),
            "observed_role": observed_role, **extra}


def geometry_of_minimum(run: Run, minimum: str) -> dict:
    return run.calc(run.payload(minimum)["opt_calc"])["final"]


def thermochemistry(run: Run, reaction: dict) -> tuple[dict | None, str | None]:
    rows = [a["payload"] for a in run.by_stage.get("thermo", [])]
    row = next(p for p in rows if p["kind"] == "reaction_thermo"
               and p["reaction_id"] == reaction["reaction_id"] and p["T_K"] == 298.15
               and p["standard_state"] == "1atm")
    if row["energy_level"] is None or row["dG_rxn_kcal"] is None:
        return None, "認証した別の端点・TSがなく、共通基準の反応Gは保存されていない。"
    species = {p["subject"]: p for p in rows if p["kind"] == "species_thermo"
               and p["T_K"] == row["T_K"]}
    reactants = reaction["monomers"] or reaction["minima"][:1]
    products = reaction["minima"][1:]
    G_R = sum(species[s]["G_hartree"] for s in reactants)
    G_P = sum(species[s]["G_hartree"] for s in products)
    points = [{"role": "reactant", "label": "R (separated)" if reaction["monomers"] else "R",
               "G_hartree": G_R, "relative_kcal": 0.0, "subjects": reactants,
               "geometry_role": None if reaction["monomers"] else "reactant"}]
    if (assoc := row.get("dG_assoc_kcal")) is not None and not reaction["monomers"]:
        points.insert(0, {"role": "reactant_separated", "label": "R separated",
                          "G_hartree": None, "relative_kcal": -assoc,
                          "subjects": [], "geometry_role": None,
                          "derived_from": "保存されたdG_assoc = G_complex - G_separated"})
    if reaction["saddle"] is not None:
        sub = reaction["saddle"]["freq_calc"]
        ts = species[sub]
        g = ts["G_hartree"]
        delta = (g - G_R) * HARTREE_TO_KCAL_MOL
        if not math.isclose(delta, row["dG_act_kcal"], abs_tol=1e-7):
            raise ValueError(f"TS common zero inconsistent: {reaction['reaction_id']}")
        points.append({"role": "ts", "label": "TS", "G_hartree": g,
                       "relative_kcal": delta, "subjects": [sub], "geometry_role": "ts"})
    dG = (G_P - G_R) * HARTREE_TO_KCAL_MOL
    if not math.isclose(dG, row["dG_rxn_kcal"], abs_tol=1e-7):
        raise ValueError(f"endpoint common zero inconsistent: {reaction['reaction_id']}")
    points.append({"role": "product", "label": "P", "G_hartree": G_P,
                   "relative_kcal": dG, "subjects": products, "geometry_role": "product"})
    subjects = [*reactants, *products]
    if reaction["saddle"] is not None:
        subjects.append(reaction["saddle"]["freq_calc"])
    for subject in dict.fromkeys(subjects):
        for field in ("freq_calc", "energy_calc"):
            if (calc_id := species[subject].get(field)) is not None:
                run.calc(calc_id)
    return {"method": row["energy_level"], "geometry_method": "pbe0-d3bj/def2-svpd",
            "units": "kcal/mol", "conditions": {"phase": "gas", "T_K": row["T_K"],
            "standard_state": row["standard_state"], "settings_sha": species[products[0]]["settings_sha"]},
            "zero": "declared separated monomers" if reaction["monomers"] else "registered reactant minimum",
            "zero_G_hartree": G_R, "points": points,
            "dG_act_kcal": row["dG_act_kcal"], "dG_rxn_kcal": row["dG_rxn_kcal"],
            "dG_eff_kcal": row["dG_eff_kcal"], "reference": row["reference"],
            "dG_act_vs_separated_kcal": row["dG_act_vs_separated_kcal"],
            "band_kcal": row["band_kcal"], "blockers": row["blockers"],
            "species_records": [{**species[s], "m_semantics": "mirror image count; not spin multiplicity"}
                                for s in dict.fromkeys(subjects)],
            "notes": ["停留点間の離散的なG図。線を描いても連続IRCのGではない。",
                      "3D表示のQRC端点は同じ登録basinに割り当てられた構造で、Gの振動計算構造とは別の場合がある。"]}, None


def coverage(run: Run) -> dict:
    file = "report/coverage.csv"
    rows = list(csv.DictReader(run.read(file).decode().splitlines()))
    for row in rows:
        row["count"] = int(row["count"])
    outcomes = {"products", "negatives", "failed", "unconnected", "not_attempted"}
    counts: dict[str, Counter] = {}
    for row in rows:
        if row["metric"] in outcomes:
            counts.setdefault(row["key"], Counter())[row["metric"]] += row["count"]
    summaries = [{"mechanism": key, "registered_queries": sum(value.values()),
                  "executed_queries": sum(value.values()) - value["not_attempted"],
                  "outcomes": dict(value)} for key, value in counts.items()]
    return {"grain": "report時点の成功artifactに保存されたDiscoveryRecordの最終分類（DFT入口選抜の更新と未試行を含む）", "rows": rows,
            "counts": summaries,
            "notes": ["CSVのattemptsとregistered_queriesはcoverage分類総数。executed_queriesは未試行を除く試行済み分類数。",
                      "productsは後段の選抜後も保持された探索辺。negativesは探索で候補化しなかった分類とDFT入口の非採用を含む。negative_reasonはその内訳。",
                      "計算失敗artifactはfailure_kindに別集計されるため、NT2の実試行総数は試行済み分類数に該当する計算失敗を加えて求める。",
                      "有限予算・代表クエリの探索結果であり、全反応空間の網羅率ではない。"]}


def irc_profile(out: Path, run: Run, example: str, reaction: dict) -> list[dict]:
    if not reaction["low_level_ts"]:
        return []
    attempt = Path(reaction["low_level_ts"][0]["file"]["path"]).parent
    paths = [attempt / f"irc_{side}/irc_{side}.irc.{direction}.trj.xyz"
             for side, direction in (("b", "backward"), ("f", "forward"))]
    if not all((run.path / p).exists() for p in paths):
        return []
    trajectories = [frames(run.read(str(p)).decode()) for p in paths]
    job = json.loads(run.read(str(attempt.parent / "job.json")))
    result = json.loads(run.read(str(attempt / "result.json")))
    zero = float(trajectories[0][-1][0])
    points = []
    for side, trajectory in enumerate(trajectories):
        stride = max(1, (len(trajectory) - 1) // 20)
        indices = sorted(set(range(0, len(trajectory), stride)) | {len(trajectory) - 1})
        for i in indices:
            if side == 1 and i == 0:
                continue
            comment, symbols, coords = trajectory[i]
            energy = float(comment)
            x = (-1 if side == 0 else 1) * i
            file = f"xyz/{example}/irc_{x:+05d}.xyz"
            dest = out / file
            dest.parent.mkdir(parents=True, exist_ok=True)
            raw = xyz_text(symbols, coords, f"saved IRC frame {i}; E={energy} Eh").encode()
            dest.write_bytes(raw)
            points.append({"x": x, "energy_hartree": energy,
                           "relative_kcal": (energy - zero) * HARTREE_TO_KCAL_MOL,
                           "source_path": str(paths[side]), "source_frame": i,
                           "geometry_file": file, "geometry_sha256": sha(raw)})
    return [{"id": "low_level_irc", "method": "GFN2-xTB (SCINE ReaDuct)",
             "electronic_temperature_K": job["key_payload"]["settings"]["electronic_temperature_K"],
             "settings": job["key_payload"]["settings"], "trial": job["key_payload"]["trial"],
             "version": job["version_pin"], "low_level_outcome": result["outcome"],
             "dft_outcome": reaction["outcome"], "units": "kcal/mol",
             "zero": "last saved backward IRC frame", "zero_energy_hartree": zero,
             "coordinate_kind": "signed saved IRC iteration (not mass-weighted distance)",
             "coordinate_label": "Saved IRC iteration (signed; no distance unit)",
             "points": sorted(points, key=lambda p: p["x"]),
             "saved_frame_counts": [len(t) for t in trajectories],
             "notes": ["描画は保存IRCから選んだ約40点。各点は実測値、補間した点はない。",
                       "GFN2-xTBの電子エネルギー。DFTの停留点G図とは別の量・別の理論。",
                       "この低レベル候補だけではDFT認証を意味しない。TMA候補はDFT端点がsame_basinで反応として採用されない。"]}]


def scan_profile(out: Path, run: Run, example: str, reaction: dict) -> list[dict]:
    minimum = run.payload(reaction["minima"][0])
    calc = run.calc(minimum["opt_calc"])
    items = [(minimum["opt_calc"], calc)] + [(a["artifact_id"], a["payload"])
             for a in run.by_stage["paths"] if a["type"] == "calculation"]
    zero = minimum["energy_hartree"]
    points = []
    for i, (calc_id, calc) in enumerate(items):
        if calc["level"]["method"] != "pbe0":
            continue
        run.calc(calc_id)
        symbols, coords, _ = run.geometry(calc["final"])
        distance = float(np.linalg.norm(coords[0] - coords[4]))
        geometry = save_structure(out, run, example, f"scan_{i:02d}",
                                  f"B-N {distance:.3f} A", calc["final"],
                                  observed_role="saved constrained scan / midpoint single point")
        points.append({"x": distance, "energy_hartree": calc["energy_hartree"],
                       "relative_kcal": (calc["energy_hartree"] - zero) * HARTREE_TO_KCAL_MOL,
                       "source_calc": calc_id, "task": calc["task"],
                       "geometry_file": geometry["xyz_file"], "geometry_sha256": geometry["sha256"]})
    return [{"id": "association_scan", "method": "pbe0-d3bj/def2-svpd",
             "units": "kcal/mol", "zero": "optimized adduct", "zero_energy_hartree": zero,
             "coordinate_kind": "B-N distance / angstrom", "points": sorted(points, key=lambda p: p["x"]),
             "coordinate_label": "B-N distance / angstrom",
             "notes": ["8点の拘束走査と保存された中点SP。現在の分解能内でbarrierless。",
                       "電子Eの走査であり、Gの曲線・捕獲速度・厳密な無障壁の証明ではない。"]}]


def extract(base: Path, out: Path) -> dict:
    examples = []
    for eid, title, system, relative, thermo_relative, reaction_id in SPECS:
        run = Run(base / relative)
        thermal = run if thermo_relative is None else Run(base / thermo_relative)
        reaction = run.payload(reaction_id)
        structures = []
        if eid == "tma_hf2_same_basin":
            for role, endpoint in zip(("declared_reactant", "declared_product"), reaction["endpoints"], strict=True):
                structures.append(save_structure(out, run, eid, role, endpoint,
                                  run.payload(f"species_{endpoint}")["geometry"], "declared input; not a certified minimum"))
            structures.append(save_structure(out, run, eid, "minimum", "common DFT minimum",
                                              geometry_of_minimum(run, reaction["minima"][0])))
        elif eid == "nh3_inversion":
            for role, species in (("reactant", "species_planar_nh3_mf1"), ("product", "species_planar_nh3_mf2")):
                structures.append(save_structure(out, run, eid, role, role,
                                                  run.payload(species)["geometry"], "saved DFT mode-follow minimum"))
        elif eid == "bh3_nh3_association":
            structures.append(save_structure(out, run, eid, "reactant", "declared entrance complex",
                                              run.payload("species_complex")["geometry"], "declared input; not separated gas monomers"))
            structures.append(save_structure(out, run, eid, "product", "optimized adduct",
                                              geometry_of_minimum(run, reaction["minima"][1])))
            for i, monomer in enumerate(reaction["monomers"]):
                structures.append(save_structure(out, run, eid, f"monomer_{i}", run.payload(monomer)["species_id"],
                                                  geometry_of_minimum(run, monomer)))
        else:
            side_ids = reaction["connection"]["side_calcs"]
            plus = run.calc(side_ids[0])
            structures.append(save_structure(out, run, eid, "reactant", "QRC side +", plus["final"], "saved DFT QRC endpoint"))
            if side_ids[0] != side_ids[1]:
                structures.append(save_structure(out, run, eid, "product", "QRC side -",
                                                  run.calc(side_ids[1])["final"], "saved DFT QRC endpoint"))
            else:
                freq = run.calc(reaction["saddle"]["freq_calc"])
                symbols, x, _ = run.geometry(freq["final"])
                _, start, _ = run.geometry(plus["start"])
                _, final, _ = run.geometry(plus["final"])
                minus_start = 2 * x - start
                rmsd, image = carry(symbols, minus_start, start, final)
                if rmsd > IMAGE_A:
                    raise ValueError(f"recorded QRC symmetry image not reproduced: {rmsd}")
                raw = xyz_text(symbols, image, "Recorded QRC minus: exact symmetry image of saved plus optimum; no second QM job").encode()
                file = f"xyz/{eid}/product.xyz"
                (out / file).write_bytes(raw)
                structures.append({"role": "product", "label": "QRC side - (symmetry image)", "xyz_file": file,
                    "atoms": [{"index": i, "element": s, "x": float(a), "y": float(b), "z": float(c)}
                              for i, (s, (a, b, c)) in enumerate(zip(symbols, image, strict=True))],
                    "source_run": str(run.path), "source_path": plus["final"]["file"]["path"], "sha256": sha(raw),
                    "observed_role": "symmetry_image_of_qrc_plus (not an independent QM endpoint)",
                    "transform_provenance": {"function": "hfauto.chemistry.identity.carry", "fit_rmsd_A": rmsd,
                        "threshold_A": IMAGE_A, "plus_calc": side_ids[0], "TS_freq_calc": reaction["saddle"]["freq_calc"],
                        "minus_start": "2 * saved TS coordinates - saved plus start coordinates"}})
        if reaction["saddle"] is not None:
            freq = run.calc(reaction["saddle"]["freq_calc"])
            structures.insert(1, save_structure(out, run, eid, "ts", "certified TS", freq["final"],
                                               "saved DFT saddle frequency geometry", imaginary_mode=freq["imaginary_modes"][0]))
        g, unavailable = thermochemistry(thermal, reaction)
        log_path = f"paths/{reaction['log']}"
        log = [json.loads(line) for line in run.read(log_path).decode().splitlines() if line.strip()]
        row_csv = list(csv.DictReader(thermal.read("report/ranking.csv").decode().splitlines()))
        ranking = next(row for row in row_csv if row["reaction_id"] == reaction_id)
        profiles = scan_profile(out, run, eid, reaction) if eid == "bh3_nh3_association" else irc_profile(out, run, eid, reaction)
        state_level = run.calc(run.payload(reaction["minima"][0])["freq_calc"])["level"]
        examples.append({"id": eid, "title": title, "system": system,
            "composition": reaction["reactants"], "reaction": reaction,
            "charge": state_level["charge"], "multiplicity": state_level["multiplicity"],
            "structures": structures, "free_energy": g, "free_energy_unavailable_reason": unavailable,
            "electronic_profiles": profiles, "decisions": {"log": log, "imag_cm1": None if reaction["saddle"] is None else reaction["saddle"]["imag_cm1"],
                "chi": None, "chi_note": "共通のchi記録欄がないため数値を作らない。",
                "reasons": reaction["reasons"], "connection_label": reaction["outcome"],
                "connection": reaction["connection"], "blockers": ranking["blockers"].split(";") if ranking["blockers"] else [],
                "ranking_row": ranking}, "coverage": coverage(run),
            "provenance": {"runs": [run.provenance()] if thermal is run else [run.provenance(), thermal.provenance()]}})
    return {"schema_version": "hfauto.visualization.v1", "generated_at": datetime.now(timezone.utc).isoformat(),
            "atom_index_base": 0, "coordinate_units": "angstrom", "energy_conversion": {"hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL},
            "atom_correspondence": "XYZ row index preserved within a reaction case; isolated monomers have their own indexing",
            "examples": examples}


def check_bundle(out: Path) -> None:
    """Check the shipped scientific quantities, atom mapping and portable geometry files."""
    data = json.loads((out / "data.json").read_text(encoding="utf-8"))
    assert [e["id"] for e in data["examples"]] == [s[0] for s in SPECS]
    for example in data["examples"]:
        mapping = []
        for structure in example["structures"]:
            raw = (out / structure["xyz_file"]).read_bytes()
            assert sha(raw) == structure["sha256"]
            _, symbols, coords = frames(raw.decode())[-1]
            atoms = structure["atoms"]
            assert [a["element"] for a in atoms] == symbols
            assert [a["index"] for a in atoms] == list(range(len(atoms)))
            assert np.allclose(coords, [[a["x"], a["y"], a["z"]] for a in atoms], atol=1e-8)
            if structure["role"] in ("reactant", "ts", "product"):
                mapping.append(symbols)
        assert not mapping or all(symbols == mapping[0] for symbols in mapping)
        if (g := example["free_energy"]) is not None:
            for point in g["points"]:
                if point["G_hartree"] is not None:
                    assert math.isclose((point["G_hartree"] - g["zero_G_hartree"]) * HARTREE_TO_KCAL_MOL,
                                        point["relative_kcal"], abs_tol=1e-7)
        for profile in example["electronic_profiles"]:
            for point in profile["points"]:
                assert math.isclose((point["energy_hartree"] - profile["zero_energy_hartree"]) * HARTREE_TO_KCAL_MOL,
                                    point["relative_kcal"], abs_tol=1e-7)
                raw = (out / point["geometry_file"]).read_bytes()
                assert sha(raw) == point["geometry_sha256"]
        for count in example["coverage"]["counts"]:
            assert count["registered_queries"] == sum(count["outcomes"].values())
            assert count["executed_queries"] == count["registered_queries"] - count["outcomes"].get("not_attempted", 0)
        if example["id"] == "tma_hf2_same_basin":
            assert example["reaction"]["outcome"] == "same_basin"
            assert example["free_energy"] is None
            assert all(s["role"] != "ts" for s in example["structures"])
        if example["id"] == "nh3_inversion":
            assert example["reaction"]["minima"][0] == example["reaction"]["minima"][1]
            assert example["reaction"]["degenerate"]
        if example["id"] == "nh3_hf_exchange":
            product = next(s for s in example["structures"] if s["role"] == "product")
            assert "not an independent QM" in product["observed_role"]
    print("Bundle check passed: 5 examples; geometry hashes, atom mapping, G/E zeros and coverage grains consistent.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", type=Path, default=Path("/home/user/hfauto_r10"))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument("--check", action="store_true", help="Check the existing portable bundle without reading source runs")
    args = parser.parse_args()
    if args.check:
        check_bundle(args.output)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    data = extract(args.runs_root, args.output)
    (args.output / "data.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    check_bundle(args.output)
    print(f"Extracted {len(data['examples'])} examples to {args.output / 'data.json'}; source runs read only.")


if __name__ == "__main__":
    main()
