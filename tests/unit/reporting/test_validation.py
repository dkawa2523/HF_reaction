"""reporting.validation on the golden manifest excerpts of VAL9/R2 runs (tests/golden G36-G39)
and the integrity of validation/: the cases, their inputs and the BH76 references."""

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from hfauto.chemistry.xyz import read_xyz
from hfauto.core.manifest import Manifest
from hfauto.reporting.validation import compare, deviations, load_cases, load_references

REPO = Path(__file__).resolve().parents[3]
VALIDATION = REPO / "validation"
GOLDEN = REPO / "tests" / "golden" / "data" / "manifests"
CASES = load_cases(VALIDATION / "cases.yaml")
REFERENCES = load_references(VALIDATION / "bh76" / "subset.yaml")
RUNS = ("s6_oh_ch4", "s5_ch3_o2", "oxalic_two_step", "h_c2h4")


def _run(name: str):
    """The view of the excerpt's dft, paths and thermo manifests, and its xyz loader."""
    manifests = []
    for stage in ("dft", "paths", "thermo"):
        path = GOLDEN / name / stage / "manifest.json"
        text = (path.read_text(encoding="utf-8") if path.is_file()
                else gzip.decompress(path.with_suffix(".json.gz").read_bytes()).decode())
        manifests.append(Manifest.model_validate_json(text))
    view = Manifest.union(manifests, run_id=name, stage_id="view")
    return view, lambda geometry: read_xyz(GOLDEN / name / geometry.file.path)


@pytest.mark.parametrize("name", RUNS)
def test_the_golden_runs_meet_their_cases(name):
    assert compare(CASES[name], *_run(name)) == []


def test_a_barrierless_outcome_against_a_bh76_barrier_is_a_deviation_not_a_pass():
    view, load = _run("h_c2h4")
    assert compare(CASES["h_c2h4"], view, load) == []
    (line,) = deviations(CASES["h_c2h4"], view, load, REFERENCES)
    assert line.startswith("C2H4 + H -> C2H5: barrierless_at_resolution") and "row 33" in line
    assert "2.0 kcal/mol (zero: separated H + C2H4)" in line
    assert deviations(CASES["s6_oh_ch4"], *_run("s6_oh_ch4"), REFERENCES) == []


def _with(name: str, index: int, **changes):
    case = CASES[name]
    reactions = list(case.reactions)
    reactions[index] = type(reactions[index]).model_validate(
        {**reactions[index].model_dump(), **changes})
    return case.model_copy(update={"reactions": tuple(reactions)})


def test_the_signature_value_and_outcome_each_decide():
    view, load = _run("s6_oh_ch4")
    abstraction = "CH4 + HO -> CH3 + H2O (elementary_step)"
    wrong = compare(_with("s6_oh_ch4", 0, equation="CH4 + HO -> CH4O + H"), view, load)
    assert f"{abstraction}: 1 reaction(s), expected 0" in wrong
    assert compare(_with("s6_oh_ch4", 0, dG_eff=(5.0, 5.5)), view, load) == [
        f"{abstraction}: dG_eff 5.9292872080974375, expected (5.0, 5.5)"]
    assert compare(_with("s6_oh_ch4", 0, outcome="degenerate_rearrangement"), view, load)
    assert compare(_with("s6_oh_ch4", 1, dG_eff=0.0), view, load)  # the parent has no value


def test_two_steps_with_one_equation_pair_by_value():
    view, load = _run("oxalic_two_step")
    parent, split1, split2 = CASES["oxalic_two_step"].reactions
    listed = CASES["oxalic_two_step"].model_copy(update={"reactions": (split2, parent, split1)})
    assert compare(listed, view, load) == []  # the listing order is free
    assert compare(_with("oxalic_two_step", 1, dG_eff=12.71), view, load) == [
        "C2H2O4 -> C2H2O4 (elementary_step): dG_eff 11.812262816616082, expected 12.71"]


def test_minima_refusals_and_products_of_runs_without_reactions():
    view, load = _run("h_c2h4")  # minima: H, C2H4, the H···C2H4 complex and C2H5
    base = CASES["h_c2h4"]
    assert compare(base.model_copy(update={"minima": ("C2H4 + H", "C2H5", "H")}), view,
                   load) == []
    assert compare(base.model_copy(update={"minima": ("CH4",), "refusals": ("declare_multiplicity",),
                                           "products": 1}), view, load) == [
        "minimum CH4: absent", "refusal declare_multiplicity: absent",
        "discovery products: 0, expected 1"]


def test_every_case_input_and_reference_exists():
    for name, case in CASES.items():
        for e in case.reactions:
            assert e.bh76 is None or e.bh76 in REFERENCES, (name, e.bh76)
        if case.superseded:
            continue
        for step in case.chain:
            pipeline = (REPO / step.pipeline if "/" in step.pipeline
                        else REPO / "configs" / "pipelines" / f"{step.pipeline}.yaml")
            system = (REPO / step.system if "/" in step.system
                      else REPO / "configs" / "systems" / f"{step.system}.yaml")
            assert pipeline.is_file() and system.is_file(), (name, step)


def test_the_bh76_copies_are_verbatim():
    sources = json.loads((VALIDATION / "bh76" / "SOURCES.json").read_text(encoding="utf-8"))
    stored = {p.relative_to(VALIDATION / "bh76").as_posix()
              for p in (VALIDATION / "bh76").rglob("*") if p.is_file()}
    assert stored == {f["path"] for f in sources["files"]} | {"SOURCES.json"}
    for f in sources["files"]:
        data = (VALIDATION / "bh76" / f["path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == f["sha256"], f["path"]


def test_every_bh76_reference_is_cited_and_the_split_is_disjoint():
    raw = yaml.safe_load((VALIDATION / "bh76" / "subset.yaml").read_text(encoding="utf-8"))
    assert all(source["cite"] for source in raw["sources"].values())
    for name, reaction in raw["reactions"].items():
        rows = [r for r in reaction.values() if isinstance(r, dict) and "ref_kcal" in r]
        assert rows and all(r["zero"] and r["verbatim"] for r in rows), name
        assert reaction["minnesota_08"]["REF1"], name
    selection, holdout = set(raw["split"]["selection"]), set(raw["split"]["holdout"])
    assert not selection & holdout and selection | holdout == set(raw["reactions"])
    assert {key.partition(".")[0] for key in REFERENCES} == set(raw["reactions"])


def test_check_fails_a_missing_run_unless_its_case_is_superseded(tmp_path, capsys):
    spec = importlib.util.spec_from_file_location("check", VALIDATION / "check.py")
    check = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(check)
    assert check.main([str(tmp_path), "hcn"]) == 1
    assert check.main([str(tmp_path), "hono_walltime"]) == 0
    assert capsys.readouterr().out.splitlines() == ["MISSING hcn", "MISSING hono_walltime"]
