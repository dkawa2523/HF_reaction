import os
from pathlib import Path

import pytest

from hfauto.backends.qm.nwchem import (
    NWCHEM_BACKEND_SCHEMA_VERSION,
    NWChemEngine,
    NWChemInputRenderer,
    extract_nwchem_projected_imaginary_mode,
    extract_nwchem_projected_imaginary_modes,
    parse_nwchem_output,
)
from hfauto.backends.registry import get_qm_engine, get_ts_engine
from hfauto.backends.ts.nwchem_neb import (
    NWChemNEBEngine,
    classify_neb_energy_profile,
    nwchem_neb_converged,
    select_neb_ts_guess,
)
from hfauto.backends.ts.pysisyphus import (
    PysisyphusEngine,
    _nwchem_qcengine_evidence,
    _optimize_irc_endpoints_with_nwchem,
    _pysisyphus_environment,
)
from hfauto.core.executables import CommandResult
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact

SAMPLE_OUTPUT = """
 Northwest Computational Chemistry Package (NWChem) 7.2.3
 Total DFT energy = -100.2583655705
 Optimization converged
 P.Frequency       -523.10       -0.00        0.00
 P.Frequency          0.00       82.40     4159.25
 Temperature = 298.15K
 Zero-Point correction to Energy = 5.943 kcal/mol ( 0.009471 au)
 Thermal correction to Enthalpy = 7.100 kcal/mol ( 0.011314 au)
 Total Entropy = 45.000 cal/mol-K
 Total times  cpu: 1.0s wall: 1.1s
"""


MULTIBLOCK_PROJECTED_MODES = """
 Northwest Computational Chemistry Package (NWChem) 7.2.3
 Total DFT energy = -100.2583655705
 Optimization converged
          -------------------------------------------------
          NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES
          -------------------------------------------------
             (Projected Frequencies expressed in cm-1)

                    1           2           3           4           5           6

 P.Frequency        0.00        0.00        0.00       10.00       20.00       30.00

           1        0.01        0.11        0.21        0.31        0.41        0.51
           2        0.02        0.12        0.22        0.32        0.42        0.52
           3        0.03        0.13        0.23        0.33        0.43        0.53
           4        0.04        0.14        0.24        0.34        0.44        0.54
           5        0.05        0.15        0.25        0.35        0.45        0.55
           6        0.06        0.16        0.26        0.36        0.46        0.56
           7        0.07        0.17        0.27        0.37        0.47        0.57
           8        0.08        0.18        0.28        0.38        0.48        0.58
           9        0.09        0.19        0.29        0.39        0.49        0.59

                    7           8           9

 P.Frequency       40.00     -250.00       50.00

           1        0.71        0.81        0.91
           2        0.72        0.82        0.92
           3        0.73        0.83        0.93
           4        0.74        0.84        0.94
           5        0.75        0.85        0.95
           6        0.76        0.86        0.96
           7        0.77        0.87        0.97
           8        0.78        0.88        0.98
           9        0.79        0.89        0.99

 ----------------------------------------------------------------------------
 Total times  cpu: 1.0s wall: 1.1s
"""


def _species(xyz: Path, artifact_id: str = "spc_hf") -> Artifact:
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="species",
        data={
            "species_id": artifact_id,
            "state": "hf_cluster",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(xyz),
        },
    )


def _assert_backend_schema_recorded(artifact: Artifact) -> None:
    field = "nwchem_backend_schema_version"
    assert artifact.data[field] == NWCHEM_BACKEND_SCHEMA_VERSION
    assert artifact.method is not None
    assert artifact.method[field] == NWCHEM_BACKEND_SCHEMA_VERSION
    assert artifact.provenance[field] == NWCHEM_BACKEND_SCHEMA_VERSION


def test_nwchem_backend_schema_version_is_explicit():
    assert NWCHEM_BACKEND_SCHEMA_VERSION == 1


def test_nwchem_parser_core_fields():
    parsed = parse_nwchem_output(SAMPLE_OUTPUT)
    assert parsed["program_version"] == "7.2.3"
    assert parsed["electronic_energy_hartree"] == -100.2583655705
    assert parsed["frequencies_cm1"] == [-523.1, 82.4, 4159.25]
    assert parsed["thermochemistry_frequencies_cm1"] == [82.4, 4159.25]
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -523.1
    assert parsed["zpe_hartree"] == 0.009471
    assert parsed["gibbs_298K_hartree"] is not None
    assert parsed["normal_termination"] is True
    assert parsed["dft_d3_applied"] is False
    assert parsed["frequency_analysis_present"] is True
    assert parsed["raw_frequency_count"] == 6


def test_nwchem_parser_reads_sole_imaginary_mode_from_multiple_projected_blocks():
    parsed = parse_nwchem_output(MULTIBLOCK_PROJECTED_MODES, atom_count=3)

    assert parsed["raw_frequency_count"] == 9
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -250.0
    assert parsed["imaginary_mode_displacements"] == [
        [0.81, 0.82, 0.83],
        [0.84, 0.85, 0.86],
        [0.87, 0.88, 0.89],
    ]
    mode = extract_nwchem_projected_imaginary_mode(
        MULTIBLOCK_PROJECTED_MODES,
        atom_count=3,
    )
    assert mode is not None
    assert mode["mode_number"] == 8
    assert mode["frequency_cm1"] == -250.0
    assert mode["cartesian_displacements"] == parsed[
        "imaginary_mode_displacements"
    ]
    assert mode["component_units"] == "amu^-1/2"
    assert "do not divide by sqrt(atomic_mass) again" in mode[
        "mass_weighting_handling"
    ]
    assert parsed["projected_imaginary_mode"] == mode


@pytest.mark.parametrize(
    "malformed_output",
    [
        MULTIBLOCK_PROJECTED_MODES.replace(
            "           9        0.79        0.89        0.99\n", ""
        ),
        MULTIBLOCK_PROJECTED_MODES.replace(
            "                    7           8           9",
            "                    1           2           3",
        ),
        MULTIBLOCK_PROJECTED_MODES.replace(
            "        0.75        0.85        0.95",
            "        0.75         NaN        0.95",
        ),
        MULTIBLOCK_PROJECTED_MODES.replace(
            "        0.00        0.00        0.00       10.00",
            "        0.00        0.00        0.00      -20.00",
        ),
    ],
    ids=["missing_cartesian_row", "duplicate_modes", "nonfinite", "two_imaginary"],
)
def test_nwchem_projected_mode_parser_fails_closed(malformed_output):
    parsed = parse_nwchem_output(malformed_output, atom_count=3)

    assert parsed["imaginary_mode_displacements"] is None


def test_nwchem_projected_mode_parser_rejects_3n_mismatch():
    parsed = parse_nwchem_output(MULTIBLOCK_PROJECTED_MODES, atom_count=2)

    assert parsed["imaginary_mode_displacements"] is None


def test_nwchem_projected_mode_extractor_requires_exactly_one_significant_mode():
    two_imaginary = MULTIBLOCK_PROJECTED_MODES.replace(
        "       10.00       20.00       30.00",
        "      -20.00       20.00       30.00",
    )

    assert (
        extract_nwchem_projected_imaginary_mode(two_imaginary, atom_count=3)
        is None
    )
    assert extract_nwchem_projected_imaginary_mode(
        MULTIBLOCK_PROJECTED_MODES,
        atom_count=2,
    ) is None

    modes = extract_nwchem_projected_imaginary_modes(
        two_imaginary,
        atom_count=3,
    )
    assert modes is not None
    assert [mode["mode_number"] for mode in modes] == [4, 8]


def test_nwchem_parser_does_not_treat_missing_frequency_output_as_zero_imaginary_modes():
    parsed = parse_nwchem_output(
        """
 Northwest Computational Chemistry Package (NWChem) 7.2.3
 Total DFT energy = -100.2583655705
 Optimization converged
 Total times  cpu: 1.0s wall: 1.1s
"""
    )

    assert parsed["frequency_analysis_present"] is False
    assert parsed["raw_frequency_count"] == 0
    assert parsed["frequencies_cm1"] == []
    assert parsed["n_imag"] is None
    assert parsed["n_imag_raw"] is None


def test_nwchem_parser_records_d3_evidence():
    parsed = parse_nwchem_output(
        SAMPLE_OUTPUT
        + "\n DFT-D3 Model\n Dispersion correction = -0.001234500000\n"
    )
    assert parsed["dft_d3_applied"] is True
    assert parsed["dispersion_correction_hartree"] == -0.0012345


def test_nwchem_parser_reads_property_mulliken_population_table():
    parsed = parse_nwchem_output(
        SAMPLE_OUTPUT
        + """
 ---------------- Total gross population on atoms ----------------
    1 N       7.0000       7.2145
    2 H       1.0000       0.8420
    3 F       9.0000       8.9435
 ---------------- Bond indices ----------------
"""
    )

    assert parsed["population_analysis_present"] is True
    assert parsed["mulliken_population_analysis"] == [
        {
            "atom_index": 0,
            "element": "N",
            "nuclear_charge": 7,
            "electron_population": 7.2145,
            "partial_charge_e": pytest.approx(-0.2145),
        },
        {
            "atom_index": 1,
            "element": "H",
            "nuclear_charge": 1,
            "electron_population": 0.842,
            "partial_charge_e": pytest.approx(0.158),
        },
        {
            "atom_index": 2,
            "element": "F",
            "nuclear_charge": 9,
            "electron_population": 8.9435,
            "partial_charge_e": pytest.approx(0.0565),
        },
    ]


def test_nwchem_parser_counts_a_soft_negative_mode_as_imaginary():
    parsed = parse_nwchem_output(
        """
 P.Frequency         -5.00        0.00       18.40
 Total times  cpu: 1.0s wall: 1.1s
"""
    )

    assert parsed["frequencies_cm1"] == [18.4]
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -5.0


def test_nwchem_parser_keeps_internal_imaginary_before_six_projected_zeros():
    parsed = parse_nwchem_output(
        """
 P.Frequency         -3.08       -0.00       -0.00       -0.00       -0.00        0.00
 P.Frequency          0.00       31.90       55.94       87.21      211.21      259.55
 P.Frequency        329.36      381.79      415.04      462.94      513.56      540.77
 P.Frequency        629.11      712.32      732.92      753.88      814.14      827.61
 P.Frequency        836.65      895.08      955.52      986.73     1004.13     1013.37
 P.Frequency       1058.60     1086.69     1131.58     1169.97     1190.06     1193.17
 P.Frequency       1275.07     1307.24     1344.53     1402.54     1504.20     1536.44
 P.Frequency       1621.84     1676.09     1688.10     2481.22     3187.27     3197.23
 P.Frequency       3205.17     3213.28     3224.40     3524.58     3588.69     3629.55
 Total times  cpu: 1.0s wall: 1.1s
"""
    )

    # A non-linear 18-atom system has six projected rigid-body zero modes and
    # 48 internal modes.  Here the six zeros are modes 2--7, so mode 1 is an
    # internal negative curvature and must not be discarded by position.
    assert parsed["raw_frequency_count"] == 54
    assert parsed["raw_frequencies_cm1"][:8] == [
        -3.08,
        -0.0,
        -0.0,
        -0.0,
        -0.0,
        0.0,
        0.0,
        31.9,
    ]
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -3.08


def test_nwchem_parser_retains_positive_soft_mode_for_thermochemistry():
    parsed = parse_nwchem_output(
        """
 P.Frequency         -0.80        0.00        0.75        6.83       18.40
 Total times  cpu: 1.0s wall: 1.1s
"""
    )

    # The established spectral field remains filtered at 10 cm-1.
    assert parsed["frequencies_cm1"] == [18.4]
    # Projected near-zero modes are removed, while the +6.83 cm-1 physical mode
    # is available to quasi-RRHO thermochemistry.
    assert parsed["thermochemistry_frequencies_cm1"] == [6.83, 18.4]
    assert parsed["n_imag"] == 0
    assert parsed["imaginary_frequency_cutoff_cm1"] == -1.0


def test_nwchem_parser_uses_configured_imaginary_noise_cutoff():
    parsed = parse_nwchem_output(
        """
 P.Frequency        -25.00       -5.00       -0.80       18.40
 Total times  cpu: 1.0s wall: 1.1s
""",
        imaginary_frequency_cutoff_cm1=-10.0,
    )

    assert parsed["imaginary_frequency_cutoff_cm1"] == -10.0
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -25.0


def test_nwchem_parser_uses_only_final_projected_frequency_analysis():
    parsed = parse_nwchem_output(
        """
 (Projected Frequencies expressed in cm-1)
 P.Frequency       -500.00       20.00       30.00
 (Projected Frequencies expressed in cm-1)
 P.Frequency        -80.00       40.00       50.00
 Total times  cpu: 1.0s wall: 1.1s
"""
    )

    assert parsed["raw_frequencies_cm1"] == [-80.0, 40.0, 50.0]
    assert parsed["raw_frequency_count"] == 3
    assert parsed["n_imag"] == 1
    assert parsed["imag_freq_cm1"] == -80.0


def test_nwchem_parser_detects_thermochemistry_temperature_mismatch():
    parsed = parse_nwchem_output(SAMPLE_OUTPUT, temperature_K=373.15)

    assert parsed["requested_temperature_K"] == 373.15
    assert parsed["thermochemistry_temperature_K"] == 298.15
    assert parsed["temperature_consistent"] is False


def test_nwchem_input_renderer_is_readable(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    rendered = NWChemInputRenderer.render(
        _species(xyz), {"functional": "pbe0", "basis": "def2-svp"}, "opt_freq"
    )
    assert "geometry units angstrom" in rendered
    assert "* library def2-svp" in rendered
    assert "xc pbe0" in rendered
    assert "task dft optimize" in rendered
    assert "task dft frequencies" in rendered


def test_nwchem_optimize_task_does_not_request_or_claim_frequencies(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    rendered = NWChemInputRenderer.render(
        _species(xyz), {"functional": "pbe0", "basis": "def2-svp"}, "optimize"
    )

    assert "task dft optimize" in rendered
    assert "task dft frequencies" not in rendered


def test_nwchem_fixed_geometry_frequency_requests_only_one_hessian(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    rendered = NWChemInputRenderer.render(
        _species(xyz), {"functional": "pbe0", "basis": "def2-svp"}, "frequency"
    )

    assert rendered.count("task dft frequencies") == 1
    assert "task dft optimize" not in rendered
    assert "task dft saddle" not in rendered


def test_nwchem_saddle_can_build_exact_initial_hessian(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    rendered = NWChemInputRenderer.render(
        _species(xyz),
        {
            "functional": "pbe0",
            "basis": "def2-svp",
            "saddle_initial_hessian": True,
        },
        "saddle_freq",
    )

    assert rendered.count("task dft frequencies") == 2
    assert rendered.index("task dft frequencies") < rendered.index("driver")
    assert "inhess 2" in rendered
    assert "firstneg" in rendered
    assert rendered.index("task dft saddle") < rendered.rindex(
        "task dft frequencies"
    )


def test_nwchem_saddle_can_follow_selected_hessian_mode(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    rendered = NWChemInputRenderer.render(
        _species(xyz),
        {
            "functional": "pbe0",
            "basis": "def2-svp",
            "saddle_initial_hessian": True,
            "saddle_mode_number": 2,
            "saddle_follow_first_negative": False,
            "driver_trust": 0.15,
            "driver_saddle_step": 0.05,
        },
        "saddle_freq",
    )

    assert "inhess 2" in rendered
    assert "moddir 2" in rendered
    assert "nofirstneg" in rendered
    assert "trust 0.15" in rendered
    assert "sadstp 0.05" in rendered


def test_nwchem_saddle_skips_duplicate_frequency_after_seed_precheck(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    rendered = NWChemInputRenderer.render(
        _species(xyz),
        {
            "functional": "pbe0",
            "basis": "def2-svp",
            "saddle_initial_hessian": True,
            "saddle_initial_hessian_only": True,
            "saddle_mode_number": 2,
            "saddle_follow_first_negative": False,
        },
        "saddle_freq",
    )

    assert rendered.count("task dft hessian") == 1
    assert rendered.count("task dft frequencies") == 1
    assert rendered.index("task dft hessian") < rendered.index("driver")


def test_nwchem_renderer_rejects_conflicting_species_charge_fields(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = _species(xyz)
    species.data["formal_charge"] = 1

    with pytest.raises(ValueError, match="Conflicting charge specifications"):
        NWChemInputRenderer.render(
            species, {"functional": "pbe0", "basis": "def2-svp"}, "single_point"
        )


def test_nwchem_renderer_rejects_zero_multiplicity(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = _species(xyz)
    species.data["multiplicity"] = 0

    with pytest.raises(ValueError, match="Multiplicity must be positive"):
        NWChemInputRenderer.render(
            species, {"functional": "pbe0", "basis": "def2-svp"}, "single_point"
        )


def test_nwchem_renderer_rejects_species_method_multiplicity_conflict(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Conflicting multiplicity specifications"):
        NWChemInputRenderer.render(
            _species(xyz),
            {"functional": "pbe0", "basis": "def2-svp", "multiplicity": 3},
            "single_point",
        )


def test_nwchem_calculation_id_includes_input_geometry_hash(tmp_path):
    first_xyz = tmp_path / "first.xyz"
    second_xyz = tmp_path / "second.xyz"
    first_xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    second_xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 1.03\n", encoding="utf-8")
    engine = NWChemEngine()
    method = {"allow_subprocess": False}

    first = engine.single_point(
        _species(first_xyz), method, str(tmp_path / "first_work")
    )
    second = engine.single_point(
        _species(second_xyz), method, str(tmp_path / "second_work")
    )

    assert first.status.category == "nwchem_not_run"
    assert second.status.category == "nwchem_not_run"
    assert first.artifact_id != second.artifact_id


def test_nwchem_calculation_id_includes_resolved_electronic_state(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    neutral = _species(xyz)
    cation = _species(xyz)
    cation.data.update({"charge": 1, "multiplicity": 2})
    engine = NWChemEngine()

    neutral_result = engine.single_point(
        neutral, {"allow_subprocess": False}, str(tmp_path / "neutral_work")
    )
    cation_result = engine.single_point(
        cation, {"allow_subprocess": False}, str(tmp_path / "cation_work")
    )

    assert neutral_result.status.category == "nwchem_not_run"
    assert cation_result.status.category == "nwchem_not_run"
    assert neutral_result.artifact_id != cation_result.artifact_id


def test_nwchem_backend_schema_is_recorded_on_failed_calculation(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    artifact = NWChemEngine().single_point(
        _species(xyz), {"allow_subprocess": False}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert artifact.data["resolved_charge"] == 0
    assert artifact.data["resolved_multiplicity"] == 1
    assert artifact.data["electron_count"] == 10
    _assert_backend_schema_recorded(artifact)


def test_nwchem_backend_schema_is_recorded_on_successful_calculation(
    tmp_path, monkeypatch
):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")

    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable", lambda *_args, **_kwargs: "nwchem"
    )

    def fake_run(input_path, _method):
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT,
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().single_point(
        _species(xyz),
        {
            "allow_subprocess": True,
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
        },
        str(tmp_path / "work"),
    )

    assert artifact.status.status == "success"
    assert artifact.data["resolved_charge"] == 0
    assert artifact.data["resolved_multiplicity"] == 1
    assert artifact.data["electron_count"] == 10
    assert artifact.method["charge"] == 0
    assert artifact.method["multiplicity"] == 1
    assert artifact.method["grid"] == "xfine"
    assert artifact.method["scf_energy_tolerance"] == 1.0e-8
    _assert_backend_schema_recorded(artifact)


def test_nwchem_timeout_failure_preserves_complete_execution_provenance(
    tmp_path, monkeypatch
):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    workdir = tmp_path / "work"
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )

    captured: dict[str, CommandResult] = {}

    def fake_run(input_path, _method):
        output_path = input_path.parent / "nwchem.out"
        stderr_path = input_path.parent / "nwchem.err"
        output_path.write_text(
            " Northwest Computational Chemistry Package (NWChem) 7.2.3\n"
            " DFT-D3 Model\n"
            " Dispersion correction = -0.001234500000\n"
            " Total DFT energy = -100.2583655705\n",
            encoding="utf-8",
        )
        stderr_path.write_text("", encoding="utf-8")
        (input_path.parent / "final-000.xyz").write_text(
            xyz.read_text(encoding="utf-8"), encoding="utf-8"
        )
        result = CommandResult(
            command=["nwchem", input_path.name],
            cwd=str(input_path.parent),
            returncode=124,
            stdout_path=str(output_path),
            stderr_path=str(stderr_path),
            duration_s=7200.0,
            timed_out=True,
            executable="nwchem",
        )
        captured["result"] = result
        return result, output_path.read_text(encoding="utf-8")

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().optimize(
        _species(xyz),
        {
            "allow_subprocess": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
        },
        str(workdir),
    )

    assert artifact.status.status == "failed"
    assert artifact.status.category == "nwchem_failed"
    assert artifact.provenance["command"] == captured["result"].to_dict()
    assert artifact.provenance["created_by"] == "NWChemEngine"
    assert artifact.provenance["input_xyz_sha256"] == sha256_file(xyz)
    assert artifact.provenance["electronic_state"] == {
        "charge": 0,
        "multiplicity": 1,
        "uhf": 0,
        "electron_count": 10,
    }
    _assert_backend_schema_recorded(artifact)


def test_optfreq_rejects_duplicate_or_extra_frequency_blocks(tmp_path, monkeypatch):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )

    def fake_run(input_path, _method):
        (input_path.parent / "final-000.xyz").write_text(
            xyz.read_text(encoding="utf-8"), encoding="utf-8"
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT + "\n P.Frequency 100.0 200.0 300.0\n",
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().optimize_frequency(
        _species(xyz), {"allow_subprocess": True}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert artifact.data["raw_frequency_count"] == 9
    assert artifact.data["expected_raw_frequency_count"] == 6
    assert artifact.data["frequency_count_complete"] is False
    assert artifact.data["n_imag"] is None


def test_failed_saddle_does_not_expose_initial_hessian_as_final_ts_frequency(
    tmp_path, monkeypatch
):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )
    initial_hessian_output = SAMPLE_OUTPUT.replace(
        "Optimization converged", "Initial Hessian evaluated"
    ) + "\n Failed to converge in maximum number of steps\n"

    def fake_run(input_path, _method):
        (input_path.parent / "final-001.xyz").write_text(
            xyz.read_text(encoding="utf-8"), encoding="utf-8"
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=1,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            initial_hessian_output,
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().saddle_frequency(
        _species(xyz), {"allow_subprocess": True}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert artifact.data["geometry_converged"] is False
    assert artifact.data["frequency_analysis_present"] is False
    assert artifact.data["frequency_count_complete"] is False
    assert artifact.data["n_imag"] is None
    assert artifact.data["zpe_hartree"] is None
    initial = artifact.data["preoptimization_hessian_analysis"]
    assert initial["frequency_analysis_present"] is True
    assert initial["n_imag"] == 1
    assert initial["imag_freq_cm1"] == -523.1


def test_nwchem_optimization_rejects_final_atom_order_change(tmp_path, monkeypatch):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )

    def fake_run(input_path, _method):
        (input_path.parent / "final-000.xyz").write_text(
            "2\nreordered\nF 0 0 0.93\nH 0 0 0\n", encoding="utf-8"
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT,
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().optimize(
        _species(xyz), {"allow_subprocess": True}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert "atom_order_ok=False" in str(artifact.status.reason)
    assert "atom_order_or_count_mismatch" in str(artifact.status.reason)


def test_nwchem_final_xyz_uses_highest_numeric_suffix_not_mtime(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    species = _species(xyz)
    engine = NWChemEngine()
    workdir = tmp_path / "work"
    input_path = engine.render_input(species, {}, workdir, "optimize")
    lower = workdir / "final-002.xyz"
    higher = workdir / "final-010.xyz"
    lower.write_text("2\nlower\nH 0.1 0 0\nF 0 0 1.03\n", encoding="utf-8")
    higher.write_text("2\nhigher\nH 0.2 0 0\nF 0 0 1.13\n", encoding="utf-8")
    lower_text = lower.read_text(encoding="utf-8")
    higher_text = higher.read_text(encoding="utf-8")
    os.utime(lower, (2_000_000_000.0, 2_000_000_000.0))
    os.utime(higher, (1_000_000_000.0, 1_000_000_000.0))

    final, error = engine._final_xyz(workdir, species, input_path, "optimize")

    assert error is None
    assert final == workdir / "final.xyz"
    assert final.read_text(encoding="utf-8") == higher_text
    assert lower.read_text(encoding="utf-8") == lower_text
    assert higher.read_text(encoding="utf-8") == higher_text

    os.utime(lower, (1_000_000_000.0, 1_000_000_000.0))
    os.utime(higher, (2_000_000_000.0, 2_000_000_000.0))
    repeated, repeated_error = engine._final_xyz(
        workdir, species, input_path, "optimize"
    )
    assert repeated_error is None
    assert repeated.read_text(encoding="utf-8") == higher_text


@pytest.mark.parametrize(
    ("highest_frame", "expected_error"),
    [
        ("not an xyz file\n", "highest_numeric_frame_unreadable"),
        (
            "2\nreordered\nF 0 0 0.93\nH 0 0 0\n",
            "atom_order_or_count_mismatch",
        ),
        ("2\nnonfinite\nH nan 0 0\nF 0 0 0.93\n", "nonfinite_coordinates"),
    ],
)
def test_nwchem_invalid_highest_numeric_frame_returns_failed_artifact(
    tmp_path, monkeypatch, highest_frame, expected_error
):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )

    def fake_run(input_path, _method):
        (input_path.parent / "final-001.xyz").write_text(
            xyz.read_text(encoding="utf-8"), encoding="utf-8"
        )
        highest = input_path.parent / "final-002.xyz"
        highest.write_text(highest_frame, encoding="utf-8")
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT,
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().optimize(
        _species(xyz), {"allow_subprocess": True}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert expected_error in str(artifact.status.reason)
    assert expected_error in str(artifact.data["final_xyz_error"])
    assert artifact.paths["final_xyz"] == ""
    assert not (tmp_path / "work" / "final.xyz").exists()
    assert (tmp_path / "work" / "final-002.xyz").read_text(
        encoding="utf-8"
    ) == highest_frame


def test_nwchem_missing_numeric_frame_returns_failed_artifact(tmp_path, monkeypatch):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    monkeypatch.setattr(
        "hfauto.backends.qm.nwchem.resolve_executable",
        lambda *_args, **_kwargs: "nwchem",
    )

    def fake_run(input_path, _method):
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(input_path.parent / "nwchem.out"),
                stderr_path=str(input_path.parent / "nwchem.err"),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT,
        )

    monkeypatch.setattr("hfauto.backends.qm.nwchem.run_nwchem_input", fake_run)
    artifact = NWChemEngine().optimize(
        _species(xyz), {"allow_subprocess": True}, str(tmp_path / "work")
    )

    assert artifact.status.status == "failed"
    assert artifact.data["final_xyz_error"] == "final_xyz_numeric_frame_missing"


def test_nwchem_neb_renderer_and_registry(tmp_path):
    start = tmp_path / "start.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    engine = NWChemNEBEngine()
    input_path = engine.render_neb_input(
        _species(start, "reactant"),
        _species(end, "product"),
        {
            "path_image_count": 7,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "neb.nw",
    )
    rendered = input_path.read_text(encoding="utf-8")
    assert "geometry endgeom" in rendered
    assert "nbeads 7" in rendered
    assert "disp vdw 3" in rendered
    assert "task dft neb ignore" in rendered
    assert isinstance(get_qm_engine("nwchem"), NWChemEngine)
    assert isinstance(get_ts_engine("nwchem_neb"), NWChemNEBEngine)


def test_nwchem_neb_path_ids_include_the_rendered_discretization(tmp_path):
    start = tmp_path / "start.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    reaction = Artifact(
        artifact_id="rxn_hf",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hf"},
    )
    method = {
        "allow_subprocess": False,
        "functional": "pbe0",
        "basis": "def2-svpd",
        "disp_vdw": 3,
        "required_program_version": "7.2.3",
    }
    engine = NWChemNEBEngine()
    paths = []
    for image_count in (5, 7):
        result = engine.search_ts(
            reaction,
            _species(start, "reactant"),
            _species(end, "product"),
            {**method, "path_image_count": image_count},
            tmp_path / f"neb_{image_count}",
        )
        paths.append(
            next(
                artifact
                for artifact in result.artifacts
                if artifact.artifact_type == "reaction_path"
            )
        )

    assert paths[0].artifact_id != paths[1].artifact_id


def test_nwchem_neb_convergence_parser_accepts_runtime_spacing_only():
    assert nwchem_neb_converged("@neb   NEB calculation converged\n") is True
    assert nwchem_neb_converged("@neb NEB calculation not converged\n") is False


def test_neb_guess_selects_highest_internal_energy(tmp_path):
    path = tmp_path / "path.xyz"
    path.write_text(
        "2\nenergy -10.0\nH 0 0 0\nF 0 0 1\n"
        "2\nenergy -9.0\nH 0 0 0.1\nF 0 0 1.1\n"
        "2\nenergy -8.0\nH 0 0 0.2\nF 0 0 1.2\n"
        "2\nenergy -9.5\nH 0 0 0.3\nF 0 0 1.3\n"
        "2\nenergy -10.2\nH 0 0 0.4\nF 0 0 1.4\n",
        encoding="utf-8",
    )
    xyz, index, reason = select_neb_ts_guess(path)
    assert index == 2
    assert xyz.coords[0][2] == 0.2
    assert reason == "highest_energy_internal_bead"


def test_neb_profile_does_not_force_a_saddle_on_monotonic_path(tmp_path):
    path = tmp_path / "monotonic.xyz"
    path.write_text(
        "1\nenergy -10.0\nH 0 0 0\n"
        "1\nenergy -9.9\nH 0 0 0.1\n"
        "1\nenergy -9.8\nH 0 0 0.2\n",
        encoding="utf-8",
    )
    profile = classify_neb_energy_profile(path)
    assert profile["neb_profile_classification"] == "monotonic_no_internal_maximum"

    barrier = tmp_path / "barrier.xyz"
    barrier.write_text(
        "1\nenergy -10.0\nH 0 0 0\n"
        "1\nenergy -9.0\nH 0 0 0.1\n"
        "1\nenergy -10.2\nH 0 0 0.2\n",
        encoding="utf-8",
    )
    profile = classify_neb_energy_profile(barrier)
    assert profile["neb_profile_classification"] == "resolved_internal_maximum"

    unresolved = tmp_path / "unresolved.xyz"
    unresolved.write_text(
        "1\nenergy -10.0\nH 0 0 0\n"
        "1\nenergy -9.2\nH 0 0 0.1\n"
        "1\nenergy -9.6\nH 0 0 0.2\n"
        "1\nenergy -9.0\nH 0 0 0.3\n",
        encoding="utf-8",
    )
    profile = classify_neb_energy_profile(unresolved)
    assert profile["neb_profile_classification"] == "unresolved_internal_profile"


def test_nwchem_neb_rejects_missing_runtime_d3_evidence_before_saddle(
    tmp_path, monkeypatch
):
    start = tmp_path / "start.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    reaction = Artifact(
        artifact_id="rxn_hf",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hf"},
    )

    def fake_run(input_path, method, output_name):
        output_path = input_path.parent / output_name
        error_path = input_path.parent / "nwchem_neb.err"
        output_path.write_text(SAMPLE_OUTPUT, encoding="utf-8")
        error_path.write_text("", encoding="utf-8")
        (input_path.parent / "hfauto.nebpath.xyz").write_text(
            "2\nenergy -10.0\nH 0 0 0\nF 0 0 0.93\n"
            "2\nenergy -9.0\nH 0 0 0.1\nF 0 0 1.03\n"
            "2\nenergy -10.0\nH 0 0 0.2\nF 0 0 1.13\n",
            encoding="utf-8",
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(output_path),
                stderr_path=str(error_path),
                duration_s=0.1,
                executable="nwchem",
            ),
            SAMPLE_OUTPUT + "\n@neb NEB calculation converged\n",
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.resolve_executable", lambda *args: "nwchem"
    )
    monkeypatch.setattr("hfauto.backends.ts.nwchem_neb.run_nwchem_input", fake_run)
    result = NWChemNEBEngine().search_ts(
        reaction,
        _species(start, "reactant"),
        _species(end, "product"),
        {
            "allow_subprocess": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "neb_work",
    )

    failure = next(
        artifact for artifact in result.artifacts if artifact.status.status == "failed"
    )
    path = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "reaction_path"
    )
    assert failure.status.category == "nwchem_neb_method_evidence_invalid"
    assert path.qc["dispersion_in_input"] is True
    assert path.qc["method_evidence_validated"] is False


def test_nwchem_neb_success_preserves_complete_method_state_and_provenance(
    tmp_path, monkeypatch
):
    start = tmp_path / "start.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    reaction = Artifact(
        artifact_id="rxn_hf",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hf"},
    )
    raw_output = (
        SAMPLE_OUTPUT
        + "\n DFT-D3 Model\n Dispersion correction = -0.001234500000\n"
        + "@neb NEB calculation converged\n"
    )

    def fake_run(input_path, method, output_name):
        output_path = input_path.parent / output_name
        error_path = input_path.parent / "nwchem_neb.err"
        output_path.write_text(raw_output, encoding="utf-8")
        error_path.write_text("", encoding="utf-8")
        (input_path.parent / "hfauto.nebpath.xyz").write_text(
            "2\nenergy -10.0\nH 0 0 0\nF 0 0 0.93\n"
            "2\nenergy -9.0\nH 0 0 0.1\nF 0 0 1.03\n"
            "2\nenergy -10.0\nH 0 0 0.2\nF 0 0 1.13\n",
            encoding="utf-8",
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=0,
                stdout_path=str(output_path),
                stderr_path=str(error_path),
                duration_s=0.1,
                executable="nwchem",
            ),
            raw_output,
        )

    def fake_saddle(self, species, method, workdir):
        output = Path(workdir) / "nwchem_saddle.out"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(raw_output, encoding="utf-8")
        return Artifact(
            artifact_id="calc_ts",
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"output": str(output), "final_xyz": species.data["xyz_path"]},
            method={"engine": "nwchem", "task": "saddle_freq"},
            data={
                "program_version": "7.2.3",
                "n_imag": 1,
                "imag_freq_cm1": -500.0,
            },
            qc={
                "real_qm_executed": True,
                "dispersion_applied": True,
                "fallback_dummy": False,
            },
            provenance={"created_by": "NWChemEngine"},
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.resolve_executable", lambda *args: "nwchem"
    )
    monkeypatch.setattr("hfauto.backends.ts.nwchem_neb.run_nwchem_input", fake_run)
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.NWChemEngine.saddle_frequency", fake_saddle
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.estimate_reaction_mode_overlap",
        lambda *args, **kwargs: (0.9, "cartesian_imaginary_mode"),
    )
    result = NWChemNEBEngine().search_ts(
        reaction,
        _species(start, "reactant"),
        _species(end, "product"),
        {
            "allow_subprocess": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "neb_success",
    )

    validated = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "reaction_validated"
    )
    assert result.success is True
    assert validated.method == {
        "engine": "nwchem",
        "backend": "nwchem_neb",
        "task": "ts_search",
        "functional": "pbe0",
        "basis": "def2-svpd",
        "disp_vdw": 3,
        "required_program_version": "7.2.3",
        "charge": 0,
        "multiplicity": 1,
    }
    assert validated.data["electron_count"] == 10
    assert validated.qc["method_evidence_validated"] is True
    assert validated.provenance["method_evidence"] == {
        "neb": True,
        "neb_method_identity": True,
        "neb_seed_sufficient": True,
        "saddle_frequency": True,
    }


def test_interrupted_neb_is_only_a_seed_for_independent_saddle_validation(
    tmp_path, monkeypatch
):
    start = tmp_path / "start.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    reaction = Artifact(
        artifact_id="rxn_hf",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hf"},
    )
    interrupted_output = (
        " Northwest Computational Chemistry Package (NWChem) 7.2.3\n"
        " Total DFT energy = -10.0\n"
        " DFT-D3 Model\n Dispersion correction = -0.001\n"
    )
    completed_saddle_output = SAMPLE_OUTPUT + (
        "\n DFT-D3 Model\n Dispersion correction = -0.001\n"
    )

    def fake_run(input_path, method, output_name):
        output_path = input_path.parent / output_name
        error_path = input_path.parent / "nwchem_neb.err"
        output_path.write_text(interrupted_output, encoding="utf-8")
        error_path.write_text("", encoding="utf-8")
        (input_path.parent / "hfauto.nebpath_000004.xyz").write_text(
            "2\nenergy -10.0\nH 0 0 0\nF 0 0 0.93\n"
            "2\nenergy -9.0\nH 0 0 0.1\nF 0 0 1.03\n"
            "2\nenergy -10.0\nH 0 0 0.2\nF 0 0 1.13\n",
            encoding="utf-8",
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=124,
                stdout_path=str(output_path),
                stderr_path=str(error_path),
                duration_s=10.0,
                timed_out=True,
                executable="nwchem",
            ),
            interrupted_output,
        )

    def fake_saddle(self, species, method, workdir):
        output = Path(workdir) / "nwchem_saddle.out"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(completed_saddle_output, encoding="utf-8")
        return Artifact(
            artifact_id="calc_ts",
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"output": str(output), "final_xyz": species.data["xyz_path"]},
            method={"engine": "nwchem", "task": "saddle_freq"},
            data={
                "program_version": "7.2.3",
                "n_imag": 1,
                "imag_freq_cm1": -500.0,
            },
            qc={
                "real_qm_executed": True,
                "dispersion_applied": True,
                "fallback_dummy": False,
            },
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.resolve_executable", lambda *args: "nwchem"
    )
    monkeypatch.setattr("hfauto.backends.ts.nwchem_neb.run_nwchem_input", fake_run)
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.NWChemEngine.saddle_frequency", fake_saddle
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_neb.estimate_reaction_mode_overlap",
        lambda *args, **kwargs: (0.9, "cartesian_imaginary_mode"),
    )

    result = NWChemNEBEngine().search_ts(
        reaction,
        _species(start, "reactant"),
        _species(end, "product"),
        {
            "allow_subprocess": True,
            "allow_interrupted_path_seed": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "interrupted_neb",
    )

    path = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "reaction_path"
    )
    validated = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "reaction_validated"
    )
    assert path.qc["path_energy_publishable"] is False
    assert path.qc["method_evidence_validated"] is False
    assert path.qc["method_identity_validated"] is True
    assert path.qc["unconverged_path_used_as_seed_only"] is True
    assert validated.qc["ts_validated_by_frequency"] is True
    assert validated.provenance["method_evidence"] == {
        "neb": False,
        "neb_method_identity": True,
        "neb_seed_sufficient": True,
        "saddle_frequency": True,
    }


def test_nwchem_no_implicit_dummy_fallback(tmp_path):
    xyz = tmp_path / "hf.xyz"
    xyz.write_text("2\nHF\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    calc = NWChemEngine().single_point(_species(xyz), {"allow_subprocess": False}, str(tmp_path / "work"))
    assert calc.status.status == "failed"
    assert calc.status.category == "nwchem_not_run"


def test_pysisyphus_nwchem_irc_input_and_no_dummy(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nTS\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")
    species = _species(xyz, "ts_hf")
    engine = PysisyphusEngine()
    input_path = engine.render_irc_input(
        species,
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "irc_step": 0.08,
        },
        tmp_path / "irc.yaml",
    )
    rendered = input_path.read_text(encoding="utf-8")
    assert "type: qcengine" in rendered
    assert "program: nwchem" in rendered
    assert "type: eulerpc" in rendered
    assert "dft__disp: vdw 3" in rendered
    assert "step_length: 0.08" in rendered
    assert f"fn: {xyz.resolve().as_posix()}" in rendered.replace("\\", "/")
    assert "endopt:" not in rendered

    reaction = Artifact(artifact_id="rxn_hf", artifact_type="reaction", data={})
    result = engine.run_irc(
        reaction,
        species,
        _species(xyz, "reactant"),
        _species(xyz, "product"),
        {
            "allow_subprocess": False,
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "not_run",
    )
    failure = next(artifact for artifact in result.artifacts if artifact.status.status == "failed")
    assert failure.status.category == "pysisyphus_irc_not_run"
    assert failure.qc["fallback_dummy"] is False


def test_pysisyphus_environment_bounds_qcengine_resources() -> None:
    env = _pysisyphus_environment(
        "/opt/pysis/bin/pysis",
        "/opt/nwchem/bin/nwchem",
        {"ncores": 4, "memory_mb": 4096},
    )

    assert env["QCENGINE_NCORES"] == "4"
    assert env["QCENGINE_MEMORY"] == "4.0"
    assert env["PATH"].split(os.pathsep)[0].replace("\\", "/") == (
        "/opt/pysis/bin"
    )


def test_pysisyphus_environment_enables_matching_nwchem_mpi(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "nwchem" / "bin"
    bin_dir.mkdir(parents=True)
    nwchem = bin_dir / "nwchem"
    prterun = bin_dir / "prterun"
    nwchem.write_text("", encoding="utf-8")
    prterun.write_text("", encoding="utf-8")

    env = _pysisyphus_environment(
        str(tmp_path / "pysis"),
        str(nwchem),
        {"ncores": 4, "memory_mb": 4096},
    )

    assert env["QCENGINE_USE_MPIEXEC"] == "true"
    assert env["QCENGINE_MPIEXEC_COMMAND"].endswith(
        "prterun -np {total_ranks}"
    )


def test_pysisyphus_prefers_canonical_irc_step_length(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nTS\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")
    path = PysisyphusEngine().render_irc_input(
        _species(xyz, "ts_hf"),
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "irc_step": 0.08,
            "irc_step_length": 0.05,
        },
        tmp_path / "canonical_step.yaml",
    )

    assert "step_length: 0.05" in path.read_text(encoding="utf-8")


def test_pysisyphus_internal_endpoint_optimizer_is_explicit(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nTS\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")
    path = PysisyphusEngine().render_irc_input(
        _species(xyz, "ts_hf"),
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "optimize_irc_endpoints": True,
            "endpoint_optimizer": "pysisyphus",
        },
        tmp_path / "with_endopt.yaml",
    )

    assert "endopt:" in path.read_text(encoding="utf-8")


def test_nwchem_endpoint_optimizer_is_decoupled_from_irc_input(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nTS\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")
    path = PysisyphusEngine().render_irc_input(
        _species(xyz, "ts_hf"),
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "optimize_irc_endpoints": True,
            "endpoint_optimizer": "nwchem",
        },
        tmp_path / "separated_endopt.yaml",
    )

    assert "endopt:" not in path.read_text(encoding="utf-8")


def test_nwchem_endpoint_optimizer_requires_two_real_convergences(
    tmp_path, monkeypatch
):
    raw_endpoints = {}
    for direction in ("forward", "backward"):
        path = tmp_path / f"{direction}.xyz"
        path.write_text("1\nend\nH 0 0 0\n", encoding="utf-8")
        raw_endpoints[direction] = path

    def fake_optimize(self, species, method, workdir):
        final = Path(workdir) / "final.xyz"
        final.parent.mkdir(parents=True, exist_ok=True)
        final.write_text("1\nend\nH 0 0 0\n", encoding="utf-8")
        return Artifact(
            artifact_id=f"calc_{species.artifact_id}",
            artifact_type="calculation",
            paths={"final_xyz": str(final)},
            qc={
                "real_qm_executed": True,
                "geometry_converged": True,
                "normal_termination": True,
            },
        )

    monkeypatch.setattr(NWChemEngine, "optimize", fake_optimize)
    evidence, calculations = _optimize_irc_endpoints_with_nwchem(
        raw_endpoints,
        reaction_id="rxn",
        state={"charge": 0, "multiplicity": 1},
        method={"functional": "pbe0", "basis": "def2-svpd", "disp_vdw": 3},
        workdir=tmp_path / "relax",
    )

    assert evidence["accepted"] is True
    assert evidence["backend"] == "nwchem"
    assert len(calculations) == 2


def test_pysisyphus_rejects_conflicting_qcengine_dispersion_keyword(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nTS\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="conflicts"):
        PysisyphusEngine().render_irc_input(
            _species(xyz, "ts_hf"),
            {
                "program": "nwchem",
                "functional": "pbe0",
                "basis": "def2-svpd",
                "disp_vdw": 3,
                "required_program_version": "7.2.3",
                "qcengine_keywords": {"dft__disp": "vdw 4"},
            },
            tmp_path / "bad_dispersion.yaml",
        )


def test_pysisyphus_nwchem_evidence_is_fail_closed(tmp_path):
    raw = tmp_path / "forward_001.000.qce_nwchem_stdout"
    raw.write_text(
        SAMPLE_OUTPUT
        + "\n DFT-D3 Model\n Dispersion correction = -0.001234500000\n",
        encoding="utf-8",
    )
    evidence = _nwchem_qcengine_evidence(tmp_path, "7.2.3")
    assert evidence["accepted"] is True
    assert evidence["raw_output_count"] == 1
    assert evidence["raw_outputs"][0]["sha256"]

    (tmp_path / "backward_002.000.qce_nwchem_stdout").write_text(
        SAMPLE_OUTPUT, encoding="utf-8"
    )
    assert _nwchem_qcengine_evidence(tmp_path, "7.2.3")["accepted"] is False


def test_pysisyphus_nwchem_evidence_enforces_basis_representation(tmp_path):
    raw = tmp_path / "calculator_000.001.qce_nwchem_stdout"
    common = (
        SAMPLE_OUTPUT
        + "\n C library def2-svpd"
        + "\n DFT-D3 Model\n Dispersion correction = -0.001234500000\n"
    )
    raw.write_text(
        common + '\n Basis "ao basis" -> "" (cartesian)\n', encoding="utf-8"
    )
    cartesian = _nwchem_qcengine_evidence(
        tmp_path,
        "7.2.3",
        expected_basis="def2-svpd",
        require_spherical=True,
    )
    assert cartesian["accepted"] is False
    assert cartesian["raw_outputs"][0]["basis_representations"] == ["cartesian"]

    raw.write_text(
        common + '\n Basis "ao basis" -> "" (spherical)\n', encoding="utf-8"
    )
    spherical = _nwchem_qcengine_evidence(
        tmp_path,
        "7.2.3",
        expected_basis="def2-svpd",
        require_spherical=True,
    )
    assert spherical["accepted"] is True
    assert spherical["raw_outputs"][0]["basis_name_matched"] is True


def test_pysisyphus_success_preserves_method_state_and_raw_evidence(
    tmp_path, monkeypatch
):
    xyz = tmp_path / "state.xyz"
    xyz.write_text("2\nstate\nH 0 0 0\nF 0 0 1.1\n", encoding="utf-8")
    ts_species = _species(xyz, "ts_hf")
    reactant = _species(xyz, "reactant")
    product = _species(xyz, "product")
    reaction = Artifact(
        artifact_id="rxn_hf_with_ts",
        artifact_type="reaction_validated",
        data={"reaction_id": "rxn_hf"},
    )

    def fake_command(command, *, cwd, **kwargs):
        workdir = Path(cwd)
        first = workdir / "finished_first.xyz"
        last = workdir / "finished_last.xyz"
        first.write_text(xyz.read_text(encoding="utf-8"), encoding="utf-8")
        last.write_text(xyz.read_text(encoding="utf-8"), encoding="utf-8")
        (workdir / "qce_nwchem_stdout_000").write_text(
            SAMPLE_OUTPUT
            + "\n C library def2-svpd"
            + '\n Basis "ao basis" -> "" (spherical)'
            + "\n DFT-D3 Model\n Dispersion correction = -0.001234500000\n",
            encoding="utf-8",
        )
        output = workdir / kwargs["stdout_name"]
        error = workdir / kwargs["stderr_name"]
        output.write_text("pysisyphus completed\n", encoding="utf-8")
        error.write_text("", encoding="utf-8")
        return CommandResult(
            command=list(command),
            cwd=str(workdir),
            returncode=0,
            stdout_path=str(output),
            stderr_path=str(error),
            duration_s=0.1,
            executable=str(command[0]),
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.resolve_executable",
        lambda *args, **kwargs: "/usr/bin/true",
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.nwchem_environment",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr("hfauto.backends.ts.pysisyphus.run_command", fake_command)
    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.endpoint_pair_match_qc",
        lambda *args, **kwargs: {"irc_validated": True},
    )
    result = PysisyphusEngine().run_irc(
        reaction,
        ts_species,
        reactant,
        product,
        {
            "allow_subprocess": True,
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "irc_success",
    )

    irc = next(
        artifact for artifact in result.artifacts if artifact.artifact_type == "irc"
    )
    assert result.success is True
    assert irc.method["engine"] == "nwchem"
    assert irc.method["backend"] == "pysisyphus"
    assert irc.method["functional"] == "pbe0"
    assert irc.method["basis"] == "def2-svpd"
    assert irc.method["disp_vdw"] == 3
    assert irc.data["electron_count"] == 10
    assert irc.qc["method_evidence_validated"] is True
    evidence = irc.provenance["nwchem_qcengine_evidence"]
    assert evidence["accepted"] is True
    assert evidence["raw_outputs"][0]["sha256"]
