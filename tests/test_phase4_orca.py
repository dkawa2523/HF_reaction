from pathlib import Path
from types import SimpleNamespace

from hfauto.backends.qm.orca import OrcaEngine, OrcaInputRenderer, parse_orca_cartesian_blocks, parse_orca_output
from hfauto.core.config import load_yaml
from hfauto.core.executables import CommandResult
from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


SAMPLE_ORCA_OUTPUT = """
Program Version 6.1.0
FINAL SINGLE POINT ENERGY     -100.123456789
SCF CONVERGED AFTER 12 CYCLES
THE OPTIMIZATION HAS CONVERGED

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
H       0.000000    0.000000    0.000000
F       0.000000    0.000000    0.940000

VIBRATIONAL FREQUENCIES
-----------------------
   0:          50.00 cm**-1
   1:        1120.00 cm**-1
   2:        3920.00 cm**-1

Zero point energy                ...      0.012345 Eh
Total Enthalpy                   ...    -100.101000 Eh
Final Gibbs free energy          ...    -100.110000 Eh
Magnitude (Debye)                :      1.234
ORCA TERMINATED NORMALLY
"""


def test_orca_parser_core_fields():
    parsed = parse_orca_output(SAMPLE_ORCA_OUTPUT)
    assert parsed["electronic_energy_hartree"] == -100.123456789
    assert parsed["zpe_hartree"] == 0.012345
    assert parsed["gibbs_298K_hartree"] == -100.110000
    assert parsed["thermal_correction_gibbs_hartree"] is not None
    assert parsed["n_imag"] == 0
    assert parsed["hf_stretch_cm1"] == 3920.0
    assert parsed["dipole_D"] == 1.234
    assert parsed["scf_converged"] is True
    assert parsed["geometry_converged"] is True
    assert parsed["normal_termination"] is True


def test_orca_cartesian_block_parser():
    blocks = parse_orca_cartesian_blocks(SAMPLE_ORCA_OUTPUT)
    assert len(blocks) == 1
    assert blocks[0].symbols == ["H", "F"]
    assert abs(blocks[0].coords[1][2] - 0.94) < 1e-12


def test_orca_input_renderer_keywords(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = Artifact(
        artifact_id="spc_hf",
        artifact_type="species",
        data={"species_id": "spc_hf", "state": "hf_cluster", "charge": 0, "multiplicity": 1, "xyz_path": str(xyz)},
    )
    inp = OrcaInputRenderer.render(species, {"method_id": "r2scan3c_orca", "ncores": 8, "memory_mb": 3000}, "opt_freq")
    assert "r2SCAN-3c" in inp
    assert "Opt" in inp and "Freq" in inp
    assert "%pal nprocs 8 end" in inp
    assert "* xyzfile 0 1" in inp


def test_orca_backend_subprocess_mock(monkeypatch, tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = Artifact(
        artifact_id="spc_hf",
        artifact_type="species",
        data={
            "species_id": "spc_hf",
            "state": "hf_cluster",
            "hf_n": 1,
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(xyz),
            "components": [{"name": "HF_1", "atom_indices": [0, 1]}],
        },
    )

    monkeypatch.setattr("hfauto.backends.qm.orca.resolve_executable", lambda name, explicit=None: "/usr/bin/orca")

    def fake_run(command, cwd, timeout_s=3600, env=None, stdout_name="stdout.txt", stderr_name="stderr.txt"):
        cwd = Path(cwd)
        (cwd / stdout_name).write_text(SAMPLE_ORCA_OUTPUT, encoding="utf-8")
        (cwd / stderr_name).write_text("", encoding="utf-8")
        (cwd / "orca.xyz").write_text("2\nopt\nH 0 0 0\nF 0 0 0.94\n", encoding="utf-8")
        return CommandResult(command=list(command), cwd=str(cwd), returncode=0, stdout_path=str(cwd / stdout_name), stderr_path=str(cwd / stderr_name), duration_s=0.1, executable="/usr/bin/orca")

    monkeypatch.setattr("hfauto.backends.qm.orca.run_command", fake_run)
    calc = OrcaEngine().optimize_frequency(
        species,
        {"method_id": "r2scan3c_orca", "allow_subprocess": True, "fallback_to_dummy": False},
        str(tmp_path / "work"),
    )
    assert calc.status.status == "success"
    assert calc.data["electronic_energy_hartree"] == -100.123456789
    assert calc.data["gibbs_298K_hartree"] == -100.110000
    assert calc.qc["real_orca_executed"] is True
    assert calc.qc["fallback_dummy"] is False
    assert Path(calc.paths["final_xyz"]).exists()


def test_phase4_orca_offline_pipeline(tmp_path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase4_orca_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(root / "examples" / "candidates.sdf")
    run_dir = run_pipeline(cfg, "phase4_offline_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "rank"
    assert any(a.artifact_type == "calculation" and a.method and a.method.get("engine") == "orca" for a in manifest.artifacts)
    assert any(a.artifact_type == "calculation" and a.qc.get("fallback_dummy") for a in manifest.artifacts)
    assert any(a.artifact_type == "species_optimized" for a in manifest.artifacts)
    assert any(a.artifact_type == "ranking" for a in manifest.artifacts)
