import dataclasses
import re

import pytest
from pydantic import ValidationError, create_model

from hfauto.chemistry.gates import Policy
from hfauto.drivers.reaction_case.state import ReactionPathsPolicy
from hfauto.pipeline.config import PipelineConfig, load
from hfauto.stages.reaction_paths import ReactionPathsConfig
from hfauto.stages.spec import StageConfig

FILES = {  # pipeline, system and site first: the arguments of load()
    "pipelines/demo.yaml": "{pipeline_id: demo, gates: {noise_cm1: 12.0}, stages: [{id: screen,"
    " stage: minima, engine: xtb, method: gfn2, init_hessian: {method: gfn2}}]}",
    "systems/h2.yaml": "{system_id: h2, species: [{id: h2, xyz: xyz/h2.xyz, multiplicity: 1}]}",
    "sites/local.yaml": "{site: local, scratch_root: /home/u/s, cores: 4}",
    "pipelines/methods/gfn2.yaml": "{id: gfn2, kind: xtb, gfn: 2}",  # next to the pipeline
    "systems/xyz/h2.xyz": "2\n\nH 0 0 0\nH 0 0 0.74\n",
}


def _write(root):  # -> the pipeline, system and site paths
    for name, text in FILES.items():
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text(text, encoding="utf-8")
    return [root / name for name in list(FILES)[:3]]


def test_load_reads_the_four_layers(tmp_path):
    resolved = load(*_write(tmp_path))
    assert list(resolved.methods) == ["gfn2"] and resolved.policy().noise_cm1 == 12.0
    assert resolved.system.species[0].xyz == (tmp_path / "systems/xyz/h2.xyz").resolve()
    assert resolved.pipeline.stages[0].settings()["init_hessian"] == {"method": "gfn2"}
    assert re.fullmatch(r"[0-9a-f]{40}(-dirty)?|unknown", resolved.code_version)


def test_methods_next_to_the_pipeline_come_first_then_configs_methods(tmp_path, monkeypatch):
    pipeline, *rest = _write(tmp_path)
    shared = tmp_path / "configs/methods/gfn2.yaml"  # the CLI's configs/ of the working dir
    shared.parent.mkdir(parents=True)
    shared.write_text("{id: gfn2, kind: xtb, gfn: 1}", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    other = tmp_path / "elsewhere/demo.yaml"  # a pipeline outside configs/ (P1c)
    other.parent.mkdir()
    other.write_text(FILES["pipelines/demo.yaml"], encoding="utf-8")
    assert load(pipeline, *rest).methods["gfn2"].gfn == 2
    assert load(other, *rest).methods["gfn2"].gfn == 1


def test_unknown_keys_are_errors():
    stage = {"id": "a", "stage": "structures"}
    for bad in ({"merge_with": "base"}, {"gates": {"no_such_gate": 1.0}}, {"stages": [stage] * 2},
                {"conditions": {"temperatures_K": [298.15]}}):  # temperatures belong to thermo
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"pipeline_id": "p", "stages": [stage], **bad})
    config = create_model("MinimaConfig", __base__=StageConfig, level=(str, "screen"))
    assert config.model_validate({"level": "dft"}).level == "dft"
    with pytest.raises(ValidationError):
        config.model_validate({"level": "dft", "no_such_key": True})


def test_seven_knobs_five_gates_and_a_typed_reaction_paths_budget():
    """U9-P2: the gates are the chemical thresholds only; the reaction-paths policy is typed and
    counts only (no wall clock)."""
    assert {f.name for f in dataclasses.fields(Policy)} == {
        "noise_cm1", "saddle_cm1", "resolution_kcal", "reaction_window_kcal", "spin_tol"}
    assert set(ReactionPathsPolicy.model_fields) == {"max_saddle_attempts", "max_split_depth"}
    stage = {"id": "a", "stage": "structures"}
    for gone in ("spin_contamination_tol", "qrc_min_drop_hartree", "thermo_zpe_tol_hartree"):
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"pipeline_id": "p", "stages": [stage],
                                           "gates": {gone: 1.0}})
    base = {"method": "m", "engines": {"qm": "q", "saddle": "s", "path": "p"}}
    config = ReactionPathsConfig.model_validate({**base, "policy": {"max_saddle_attempts": 3}})
    assert config.policy == ReactionPathsPolicy(max_saddle_attempts=3)
    assert ReactionPathsConfig.model_validate(base).policy.max_saddle_attempts == 2
    for gone in ("screen", "string_beads", "qrc_bounds_A", "walltime_s", "walltime_h"):
        with pytest.raises(ValidationError):
            ReactionPathsConfig.model_validate({**base, "policy": {gone: 1}})


def test_stationary_point_stages_must_share_one_method():
    stages = [{"id": "screen", "stage": "minima", "level": "screen", "method": "gfn2"},
              {"id": "dft", "stage": "minima", "level": "dft", "method": "a"},
              {"id": "paths", "stage": "reaction-paths", "method": "b"}]
    with pytest.raises(ValidationError, match=r"different methods: \{'dft': 'a', 'paths': 'b'\}"):
        PipelineConfig.model_validate({"pipeline_id": "p", "stages": stages})
    stages[2]["method"] = "a"  # the xTB screen level is not a stationary-point level
    assert len(PipelineConfig.model_validate({"pipeline_id": "p", "stages": stages}).stages) == 3
