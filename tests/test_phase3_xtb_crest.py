from pathlib import Path
from types import SimpleNamespace

from hfauto.backends.conformer.crest import assign_relative_energies, parse_crest_energies, read_multixyz
from hfauto.backends.qm.xtb import XTBEngine, parse_xtb_total_energy
from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_xtb_energy_parser():
    assert parse_xtb_total_energy("| TOTAL ENERGY       -5.123456 Eh") == -5.123456
    assert parse_xtb_total_energy("total E = -10.5") == -10.5
    assert parse_xtb_total_energy("no energy") is None


def test_crest_multixyz_and_energy_parser(tmp_path):
    xyz = tmp_path / "crest_conformers.xyz"
    xyz.write_text(
        "2\nconf1\nH 0 0 0\nF 0 0 0.93\n"
        "2\nconf2\nH 0 0 0\nF 0 0 0.95\n",
        encoding="utf-8",
    )
    blocks = read_multixyz(xyz)
    assert len(blocks) == 2
    assert blocks[0].symbols == ["H", "F"]
    energies = tmp_path / "crest.energies"
    energies.write_text("1 -100.000000\n2 -99.999000\n", encoding="utf-8")
    raw = parse_crest_energies(tmp_path, blocks)
    rel = assign_relative_energies(raw, len(blocks))
    assert rel[0] == 0.0
    assert rel[1] > 0.0


def test_xtb_backend_subprocess_mock(monkeypatch, tmp_path):
    xyz = tmp_path / "input_species.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = Artifact(
        artifact_id="spc_hf1",
        artifact_type="species",
        paths={"xyz": str(xyz)},
        data={
            "species_id": "spc_hf1",
            "state": "hf_cluster",
            "hf_n": 1,
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(xyz),
            "components": [{"name": "HF_1", "atom_indices": [0, 1]}],
        },
    )

    monkeypatch.setattr("hfauto.backends.qm.xtb.shutil.which", lambda exe: "/usr/bin/xtb")

    def fake_run(cmd, cwd, text, capture_output, timeout):
        Path(cwd, "xtbopt.xyz").write_text("2\nopt\nH 0 0 0\nF 0 0 0.94\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="TOTAL ENERGY      -1.234567\n", stderr="")

    monkeypatch.setattr("hfauto.backends.qm.xtb.subprocess.run", fake_run)
    calc = XTBEngine().optimize_frequency(
        species,
        {"method_id": "gfn2-xtb", "allow_subprocess": True, "gfn": 2},
        str(tmp_path / "work"),
    )
    assert calc.status.status == "success"
    assert calc.data["electronic_energy_hartree"] == -1.234567
    assert calc.qc["fallback_dummy"] is False
    assert calc.qc["geometry_sane"] is True
    assert Path(calc.paths["final_xyz"]).exists()


def test_phase3_xtb_crest_pipeline_safe_fallback(tmp_path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase3_xtb_crest.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(root / "examples" / "candidates.sdf")
    run_dir = run_pipeline(cfg, "phase3_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "rank"
    assert any(a.artifact_type == "conformer" and a.qc.get("fallback_reason") for a in manifest.artifacts)
    assert any(a.artifact_type == "calculation" and a.method and a.method.get("engine") == "xtb" for a in manifest.artifacts)
    assert any(a.artifact_type == "preopt_geometry" for a in manifest.artifacts)
    assert any(a.artifact_type == "ranking" for a in manifest.artifacts)
