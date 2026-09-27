"""configs/ (design §7.5): every file validates, and the four layers keep their invariants."""

from pathlib import Path, PurePosixPath, PureWindowsPath

import pytest
import yaml

from hfauto.chemistry.electronic_state import check_electronic_state
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.system import load_system
from hfauto.pipeline import config
from hfauto.stages import catalog

CONFIGS = Path(__file__).resolve().parents[2] / "configs"
FORBIDDEN = {"fallback", "dummy", "require_real_qm", "extra_args", "engine_order", "mode", "xfine"}


def files(layer):
    return sorted((CONFIGS / layer).glob("*.yaml"))


def scalars(node):  # every key and leaf value of a YAML tree
    if isinstance(node, dict):
        node = [*node, *node.values()]
    return [s for item in node for s in scalars(item)] if isinstance(node, list) else [node]


@pytest.mark.parametrize("path", files("pipelines"), ids=lambda p: p.stem)
def test_pipeline_stages_and_their_methods_validate(path):
    pipeline = config.PipelineConfig.model_validate(yaml.safe_load(path.read_text("utf-8")))
    for entry in pipeline.stages:
        catalog.get(entry.stage).spec.config.model_validate(entry.settings())
    assert config.method_ids([e.settings() for e in pipeline.stages]) <= {
        p.stem for p in files("methods")}


def test_sites_and_methods_validate():
    assert all(config.load_site(path).engines for path in files("sites"))
    assert all(config.load_method(CONFIGS / "methods", p.stem) for p in files("methods"))


@pytest.mark.parametrize("path", files("systems"), ids=lambda p: p.stem)
def test_system_loads_and_its_xyz_fit_the_state_and_the_reaction_atom_order(path):
    system = load_system(path)
    xyz = {s.id: read_xyz(s.xyz) for s in system.species if s.xyz is not None}
    for species in system.species:  # elements, charge and multiplicity (xyz default 1)
        if species.id in xyz:
            m = species.multiplicity
            check_electronic_state(xyz[species.id].symbols, species.charge, 1 if m is None else m,
                                   declared=m is not None)
    for reaction in system.reactions:
        assert xyz[reaction.reactant].symbols == xyz[reaction.product].symbols


def test_no_forbidden_keys_and_absolute_paths_only_in_sites():
    for path in sorted(CONFIGS.glob("*/*.yaml")):
        text = [v for v in scalars(yaml.safe_load(path.read_text("utf-8"))) if isinstance(v, str)]
        assert not FORBIDDEN & set(text), path
        absolute = [v for v in text if PurePosixPath(v).is_absolute() or PureWindowsPath(v).drive]
        assert path.parent.name == "sites" or not absolute, path
