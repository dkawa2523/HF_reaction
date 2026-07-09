from pathlib import Path

import numpy as np

from hfauto.backends.ts.base import endpoints_compatible, reaction_coordinate_progress_score
from hfauto.backends.ts.orca_nebts import ORCANEBTSEngine, ORCANEBTSInputRenderer
from hfauto.chemistry.hf_builder import XYZ, write_xyz
from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def _endpoint_artifacts(tmp_path: Path):
    rc_xyz = tmp_path / "rc.xyz"
    ip_xyz = tmp_path / "ip.xyz"
    # N···H-F -> N-H···F toy endpoints, atom order identical.
    write_xyz(XYZ(["N", "H", "F"], np.array([[0.0, 0.0, 0.0], [1.60, 0.0, 0.0], [2.53, 0.0, 0.0]]), "rc"), rc_xyz)
    write_xyz(XYZ(["N", "H", "F"], np.array([[0.0, 0.0, 0.0], [1.05, 0.0, 0.0], [2.45, 0.0, 0.0]]), "ip"), ip_xyz)
    rc_data = {
        "species_id": "spc_rc",
        "state": "reactant_complex",
        "charge": 0,
        "multiplicity": 1,
        "hf_n": 1,
        "xyz_path": str(rc_xyz),
        "atom_order_key": "map1",
        "components": [{"name": "B", "atom_indices": [0]}, {"name": "HF_1", "atom_indices": [1, 2]}],
        "reaction_coordinate": {"type": "distance_difference", "atoms": {"base_atom": 0, "transfer_h": 1, "leaving_f": 2}},
    }
    ip_data = {**rc_data, "species_id": "spc_ip", "state": "ion_pair", "xyz_path": str(ip_xyz)}
    rxn_data = {
        "reaction_id": "rxn_toy_PT",
        "reaction_type": "proton_transfer",
        "mol_id": "toy",
        "site_id": "toy_N0",
        "site_type": "aliphatic_amine",
        "hf_n": 1,
        "reactant_species_id": "spc_rc",
        "product_species_id": "spc_ip",
        "atom_order_key": "map1",
        "reaction_coordinate": rc_data["reaction_coordinate"],
    }
    rc = Artifact(artifact_id="spc_rc", artifact_type="species", paths={"xyz": str(rc_xyz)}, data=rc_data)
    ip = Artifact(artifact_id="spc_ip", artifact_type="species", paths={"xyz": str(ip_xyz)}, data=ip_data)
    rxn = Artifact(artifact_id="rxn_toy_PT", artifact_type="reaction", parents=["spc_rc", "spc_ip"], data=rxn_data)
    return rxn, rc, ip


def test_orca_nebts_input_renderer_contains_required_blocks(tmp_path):
    _rxn, rc, ip = _endpoint_artifacts(tmp_path)
    text = ORCANEBTSInputRenderer.render_nebts(rc, ip, {"method_id": "r2scan3c_orca_nebts", "n_images": 6, "ncores": 4, "base_keywords": ["r2SCAN-3c"]}, product_xyz_name="product.xyz")
    assert "NEB-TS" in text
    assert "%neb" in text
    assert 'NEB_END_XYZFILE "product.xyz"' in text
    assert "NImages 6" in text
    assert "* xyzfile 0 1" in text


def test_orca_nebts_offline_fallback_returns_ts_and_path(tmp_path):
    rxn, rc, ip = _endpoint_artifacts(tmp_path)
    ok, compat = endpoints_compatible(rc, ip)
    assert ok, compat
    result = ORCANEBTSEngine().search_ts(
        rxn,
        rc,
        ip,
        {"method_id": "r2scan3c_orca_nebts", "allow_subprocess": False, "fallback_to_dummy": True},
        tmp_path / "ts_work",
    )
    types = {a.artifact_type for a in result.artifacts}
    assert {"species", "calculation", "ts_path", "reaction_validated"} <= types
    assert result.ts_species_id == "ts_toy_PT"
    calc = next(a for a in result.artifacts if a.artifact_type == "calculation")
    assert calc.qc["fallback_dummy"] is True
    assert calc.qc["ts_validated_by_frequency"] is True
    assert Path(next(a for a in result.artifacts if a.artifact_type == "ts_path").paths["nebts_input"]).exists()


def test_reaction_coordinate_progress_score(tmp_path):
    rxn, rc, ip = _endpoint_artifacts(tmp_path)
    result = ORCANEBTSEngine().search_ts(rxn, rc, ip, {"allow_subprocess": False, "fallback_to_dummy": True}, tmp_path / "ts_work")
    ts = next(a for a in result.artifacts if a.artifact_type == "species")
    score = reaction_coordinate_progress_score(rc, ip, ts.data["xyz_path"], rxn.data["reaction_coordinate"])
    assert score["reaction_coordinate_between_endpoints"] is True
    assert score["reaction_coordinate_progress_score"] >= 0.5


def test_phase5_orca_ts_offline_pipeline(tmp_path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase5_ts_irc_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(root / "examples" / "candidates.sdf")
    run_dir = run_pipeline(cfg, "phase5_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "rank"
    assert any(a.artifact_type == "ts_path" for a in manifest.artifacts)
    assert any(a.artifact_type == "irc" for a in manifest.artifacts)
    assert any(a.artifact_type == "reaction_validated" and a.qc.get("ts_engine") == "orca_nebts" for a in manifest.artifacts)
    ts_calcs = [a for a in manifest.iter_artifacts("calculation") if (a.method or {}).get("engine") == "orca_nebts"]
    assert ts_calcs
    assert any(a.qc.get("fallback_dummy") for a in ts_calcs)
