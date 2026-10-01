"""configs/ (design §7.5): every file validates, and the four layers keep their invariants."""

from pathlib import Path, PurePosixPath, PureWindowsPath

import numpy as np
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
    assert all(config.load_method(p) for p in files("methods"))


def pipeline(name):
    path = CONFIGS / "pipelines" / f"{name}.yaml"
    return config.PipelineConfig.model_validate(yaml.safe_load(path.read_text("utf-8")))


def method(method_id):
    return config.load_method(CONFIGS / "methods" / f"{method_id}.yaml")


@pytest.mark.parametrize("name", ["discover", "known_endpoints"])
def test_default_ranking_is_m06_2x_single_points_on_the_pbe0_stationary_points(name):
    """One energy layer after the paths: the sp method is thermo's energy_method, on the
    stationary points' numerics but another functional, basis and D3 (not PBE0 reused)."""
    cfg = pipeline(name)
    assert [e.id for e in cfg.stages][-4:] == ["paths", "sp", "thermo", "report"]
    sp, thermo = (cfg.entry(i).settings() for i in ("sp", "thermo"))
    assert sp["methods"] == [thermo["energy_method"]] == ["m06-2x-d3_def2-tzvpd"]
    layer, stationary = method(thermo["energy_method"]), method(
        cfg.entry("paths").settings()["method"])
    assert (layer.kind, layer.functional, layer.basis, layer.dispersion) == (
        "dft", "m06-2x", "def2-tzvpd", "d3zero")
    assert (layer.grid, layer.scf_energy_tol) == (stationary.grid, stationary.scf_energy_tol)
    assert (stationary.functional, stationary.basis, stationary.dispersion) == (
        "pbe0", "def2-svpd", "d3bj")


def test_method_panel_shows_levels_and_defines_no_second_energy_layer():
    assert [(e.id, e.stage) for e in pipeline("method_panel").stages] == [
        ("panel_sp", "sp"), ("panel_report", "report")]


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


_C3 = np.array([[-0.5, -np.sqrt(0.75), 0.0], [np.sqrt(0.75), -0.5, 0.0], [0.0, 0.0, 1.0]])
SADDLE_SEEDS = {  # a seed off its symmetry by more than noise would relax to a minimum instead
    "nh3_planar_seed/d3h.xyz": (np.diag([1.0, 1.0, -1.0]), _C3),
    "dme_c2v_seed/c2v_eclipsed.xyz": (np.diag([-1.0, 1.0, 1.0]), np.diag([1.0, 1.0, -1.0])),
}


@pytest.mark.parametrize("name", SADDLE_SEEDS)
def test_saddle_seeds_keep_their_point_group_exactly(name):
    xyz = read_xyz(CONFIGS / "systems" / "xyz" / name)
    x, symbols = np.asarray(xyz.coords), np.asarray(xyz.symbols)
    for op in SADDLE_SEEDS[name]:  # every image atom lands on a like atom
        d = np.linalg.norm(x[:, None] - (x @ op.T)[None], axis=2)
        d[symbols[:, None] != symbols[None]] = np.inf
        assert d.min(axis=0).max() < 1e-6, name


def test_no_forbidden_keys_and_absolute_paths_only_in_sites():
    for path in sorted(CONFIGS.glob("*/*.yaml")):
        text = [v for v in scalars(yaml.safe_load(path.read_text("utf-8"))) if isinstance(v, str)]
        assert not FORBIDDEN & set(text), path
        absolute = [v for v in text if PurePosixPath(v).is_absolute() or PureWindowsPath(v).drive]
        assert path.parent.name == "sites" or not absolute, path
