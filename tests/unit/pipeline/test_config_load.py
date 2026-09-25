import re

import pytest
from pydantic import ValidationError, create_model

from hfauto.pipeline.config import PipelineConfig, load
from hfauto.stages.spec import StageConfig

FILES = {  # pipeline, system and site first: the arguments of load()
    "pipelines/demo.yaml": "{pipeline_id: demo, gates: {noise_cm1: 12.0}, stages: [{id: screen,"
    " stage: minima, engine: xtb, method: gfn2, init_hessian: {method: gfn2}}]}",
    "systems/h2.yaml": "{system_id: h2, species: [{id: h2, xyz: xyz/h2.xyz}]}",
    "sites/local.yaml": "{site: local, scratch_root: /home/u/s, cores: 4, memory_mb: 11000}",
    "methods/gfn2.yaml": "{id: gfn2, kind: xtb, gfn: 2}",
    "systems/xyz/h2.xyz": "2\n\nH 0 0 0\nH 0 0 0.74\n",
}


def test_load_reads_the_four_layers(tmp_path):
    for name, text in FILES.items():
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text(text, encoding="utf-8")
    resolved = load(*(tmp_path / name for name in list(FILES)[:3]))
    assert list(resolved.methods) == ["gfn2"] and resolved.policy().noise_cm1 == 12.0
    assert resolved.system.species[0].xyz == (tmp_path / "systems/xyz/h2.xyz").resolve()
    assert resolved.pipeline.stages[0].settings()["init_hessian"] == {"method": "gfn2"}
    assert re.fullmatch(r"[0-9a-f]{40}(-dirty)?|unknown", resolved.code_version)


def test_unknown_keys_are_errors():
    stage = {"id": "a", "stage": "structures"}
    for bad in ({"merge_with": "base"}, {"gates": {"no_such_gate": 1.0}}, {"stages": [stage] * 2}):
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"pipeline_id": "p", "stages": [stage], **bad})
    config = create_model("MinimaConfig", __base__=StageConfig, level=(str, "screen"))
    assert config.model_validate({"level": "dft"}).level == "dft"
    with pytest.raises(ValidationError):
        config.model_validate({"level": "dft", "fallback_to_dummy": True})
