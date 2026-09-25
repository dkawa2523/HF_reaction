from pathlib import Path

import pytest

from hfauto.backends.registry import get_ts_engine
from hfauto.backends.ts.nwchem_string import (
    NWChemStringEngine,
    nwchem_string_converged,
    parse_string_path_energies,
)
from hfauto.chemistry.path_diagnostics import parse_string_optimization_history
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.executables import CommandResult
from hfauto.core.schemas.artifact import Artifact


def _species(path: Path, artifact_id: str) -> Artifact:
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="species_optimized",
        paths={"xyz": str(path)},
        data={
            "species_id": artifact_id,
            "xyz_path": str(path),
            "charge": 0,
            "multiplicity": 1,
        },
    )


def _endpoints(tmp_path: Path) -> tuple[Artifact, Artifact, Artifact]:
    start = tmp_path / "start.xyz"
    middle = tmp_path / "middle.xyz"
    end = tmp_path / "end.xyz"
    start.write_text("2\nstart\nH 0 0 0\nF 0 0 0.93\n", encoding="utf-8")
    middle.write_text("2\nmiddle\nH 0 0 0.2\nF 0 0 1.13\n", encoding="utf-8")
    end.write_text("2\nend\nH 0 0 0.4\nF 0 0 1.33\n", encoding="utf-8")
    reaction = Artifact(
        artifact_id="rxn_hf",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hf"},
    )
    return reaction, _species(start, "reactant"), _species(end, "product")


def _string_run_result(
    raw_output: str, *, returncode: int, final: bool
):
    """Return a small deterministic NWChem runner for path-contract tests."""

    def run(input_path, method, output_name):
        output = input_path.parent / output_name
        error = input_path.parent / "nwchem_string.err"
        output.write_text(raw_output, encoding="utf-8")
        error.write_text("", encoding="utf-8")
        path_name = (
            "hfauto.stringpath_final.xyz"
            if final
            else "hfauto.stringpath_000004.xyz"
        )
        energy_name = (
            "hfauto.string_final_epath"
            if final
            else "hfauto.string_epath"
        )
        (input_path.parent / path_name).write_text(
            "2\ngeometry\nH 0 0 0\nF 0 0 0.93\n"
            "2\ngeometry\nH 0 0 0.1\nF 0 0 1.03\n"
            "2\ngeometry\nH 0 0 0.3\nF 0 0 1.23\n"
            "2\ngeometry\nH 0 0 0.4\nF 0 0 1.33\n",
            encoding="utf-8",
        )
        (input_path.parent / energy_name).write_text(
            "0.0 -10.0\n0.33 -8.5\n0.67 -9.0\n1.0 -10.1\n",
            encoding="utf-8",
        )
        return (
            CommandResult(
                command=["nwchem", input_path.name],
                cwd=str(input_path.parent),
                returncode=returncode,
                stdout_path=str(output),
                stderr_path=str(error),
                duration_s=0.1,
                executable="nwchem",
            ),
            raw_output,
        )

    return run


def test_string_renderer_uses_fixed_endpoints_and_evidenced_middle(
    tmp_path: Path,
) -> None:
    reaction, reactant, product = _endpoints(tmp_path)
    middle = tmp_path / "middle.xyz"
    engine = NWChemStringEngine()

    path = engine.render_string_input(
        reactant,
        product,
        {
            "path_strategy": "adaptive_double_ended_path",
            "adaptive_path_seed_candidate": {"xyz_path": str(middle)},
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "string_nbeads": 9,
            "string_interpol": 3,
        },
        tmp_path / "string.nw",
    )

    rendered = path.read_text(encoding="utf-8")
    assert reaction.artifact_id == "rxn_hf"
    assert "geometry midgeom" in rendered
    assert "nbeads 9" in rendered
    assert "interpol 3" in rendered
    assert "freeze1 .true." in rendered
    assert "freezeN .true." in rendered
    assert "hasmiddle" in rendered
    assert "task dft string ignore" in rendered
    assert isinstance(get_ts_engine("nwchem_string"), NWChemStringEngine)


def test_string_convergence_and_step_parser_are_fail_closed() -> None:
    converged_output = (
        "@zts  1 0.120000 0.500000 -10 -9 -10 -9 0\n"
        "@zts  2 0.010000 0.040000 -10 -9 -10 -9 0\n"
        "@zts The string calculation converged\n"
    )
    failed_output = "@zts The string calculation failed to converge\n"

    assert nwchem_string_converged(converged_output) is True
    assert nwchem_string_converged(failed_output) is False
    history = parse_string_optimization_history(
        converged_output, converged=True
    )
    assert history.iterations == 2
    assert history.rms_displacement == 0.01
    assert history.max_displacement == 0.04
    assert history.max_gradient is None


def test_string_stagnation_uses_reported_tolerance_and_two_windows() -> None:
    stalled = (
        "@zts Covergence Tolerance = 0.00045\n"
        "@zts  1 0.0110 0.32 -10 -9 -10 -9 0\n"
        "@zts  2 0.0039 0.09 -10 -9 -10 -9 0\n"
        "@zts  3 0.0150 0.39 -10 -9 -10 -9 0\n"
        "@zts  4 0.0037 0.09 -10 -9 -10 -9 0\n"
        "@zts  5 0.0140 0.35 -10 -9 -10 -9 0\n"
        "@zts  6 0.0040 0.10 -10 -9 -10 -9 0\n"
    )
    improving = (
        "@zts Convergence Tolerance = 0.00045\n"
        "@zts  1 0.1000 0.32 -10 -9 -10 -9 0\n"
        "@zts  2 0.0800 0.20 -10 -9 -10 -9 0\n"
        "@zts  3 0.0600 0.15 -10 -9 -10 -9 0\n"
        "@zts  4 0.0300 0.10 -10 -9 -10 -9 0\n"
        "@zts  5 0.0200 0.08 -10 -9 -10 -9 0\n"
        "@zts  6 0.0100 0.05 -10 -9 -10 -9 0\n"
    )

    assert parse_string_optimization_history(
        stalled, converged=False
    ).stagnant is True
    assert parse_string_optimization_history(
        improving, converged=False
    ).stagnant is False


def test_string_energy_parser_reads_final_image_block(tmp_path: Path) -> None:
    energy_path = tmp_path / "job.string_final_epath"
    energy_path.write_text(
        "# String path\n"
        " 0.0 -10.0\n"
        " 0.5 -9.0\n"
        " 1.0 -10.2\n",
        encoding="utf-8",
    )

    assert parse_string_path_energies(energy_path, 3) == [-10.0, -9.0, -10.2]
    assert parse_string_path_energies(energy_path, 4) == []


def test_string_renderer_rejects_initial_path_with_wrong_endpoint(
    tmp_path: Path,
) -> None:
    _reaction, reactant, product = _endpoints(tmp_path)
    initial = tmp_path / "initial.xyz"
    initial.write_text(
        "2\nbad start\nH 0 0 0\nF 0 0 2.0\n"
        "2\nmiddle\nH 0 0 0.2\nF 0 0 1.13\n"
        "2\nend\nH 0 0 0.4\nF 0 0 1.33\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="does not start at the reactant"):
        NWChemStringEngine().render_string_input(
            reactant,
            product,
            {
                "initial_path_xyz": str(initial),
                "functional": "pbe0",
                "basis": "def2-svpd",
            },
            tmp_path / "string.nw",
        )


def test_string_renderer_rejects_missing_requested_initial_path(
    tmp_path: Path,
) -> None:
    _reaction, reactant, product = _endpoints(tmp_path)

    with pytest.raises(ValueError, match="initial_path_xyz is not readable"):
        NWChemStringEngine().render_string_input(
            reactant,
            product,
            {
                "initial_path_xyz": str(tmp_path / "missing.xyz"),
                "functional": "pbe0",
                "basis": "def2-svpd",
            },
            tmp_path / "string.nw",
        )


def test_string_renderer_can_resample_a_validated_initial_path(
    tmp_path: Path,
) -> None:
    _reaction, reactant, product = _endpoints(tmp_path)
    initial = tmp_path / "initial.xyz"
    initial.write_text(
        "2\nstart\nH 0 0 0\nF 0 0 0.95\n"
        "2\nmiddle\nH 0 0 0.2\nF 0 0 1.13\n"
        "2\nend\nH 0 0 0.4\nF 0 0 1.30\n",
        encoding="utf-8",
    )

    rendered_path = NWChemStringEngine().render_string_input(
        reactant,
        product,
        {
            "initial_path_xyz": str(initial),
            "resample_initial_path": True,
            "string_nbeads": 5,
            "functional": "pbe0",
            "basis": "def2-svpd",
        },
        tmp_path / "string.nw",
    )

    rendered = rendered_path.read_text(encoding="utf-8")
    images = read_xyz_trajectory(tmp_path / "initial_path.xyz")
    assert "nbeads 5" in rendered
    assert len(images) == 5
    assert images[0].coords.tolist() == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.93]]
    assert images[-1].coords.tolist() == [[0.0, 0.0, 0.4], [0.0, 0.0, 1.33]]


def test_string_path_ids_include_the_rendered_discretization(
    tmp_path: Path, monkeypatch
) -> None:
    reaction, reactant, product = _endpoints(tmp_path)
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_string.resolve_executable",
        lambda *args: None,
    )
    engine = NWChemStringEngine()

    results = [
        engine.search_ts(
            reaction,
            reactant,
            product,
            {
                "path_strategy": "reparameterized_double_ended_path",
                "string_nbeads": image_count,
                "functional": "pbe0",
                "basis": "def2-svpd",
            },
            tmp_path / f"images_{image_count}",
        )
        for image_count in (5, 9)
    ]
    path_ids = [
        next(
            artifact.artifact_id
            for artifact in result.artifacts
            if artifact.artifact_type == "reaction_path"
        )
        for result in results
    ]

    assert len(set(path_ids)) == 2


def test_converged_string_publishes_seed_without_claiming_ts(
    tmp_path: Path, monkeypatch
) -> None:
    reaction, reactant, product = _endpoints(tmp_path)
    raw_output = """
 Northwest Computational Chemistry Package (NWChem) 7.2.3
 Total DFT energy = -10.0
 DFT-D3 Model
 Dispersion correction = -0.001
 @zts  1 0.010000 0.040000 -10 -9 -10 -9 0
 @zts The string calculation converged
 Total times  cpu: 1.0s wall: 1.1s
"""

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_string.resolve_executable",
        lambda *args: "nwchem",
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_string.run_nwchem_input",
        _string_run_result(raw_output, returncode=0, final=True),
    )
    result = NWChemStringEngine().search_ts(
        reaction,
        reactant,
        product,
        {
            "path_strategy": "adaptive_double_ended_path",
            "endpoint_basin_status": "distinct_basin",
            "allow_subprocess": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "work",
    )

    saddle_attempt = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    seed = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "species"
    )
    assert result.success is False
    assert saddle_attempt.data["diagnosis"] == "resolved_saddle_candidate"
    assert saddle_attempt.data["search_evidence_validated"] is True
    assert saddle_attempt.data["candidate"]["path_image_index"] == 1
    assert seed.data["state"] == "saddle_seed"
    assert seed.qc["selection"] == "highest_energy_internal_image"
    assert not any(
        artifact.artifact_type == "reaction_validated"
        for artifact in result.artifacts
    )
    assert result.artifacts[-1].status.category == "saddle_refinement_required"


def test_interrupted_string_can_supply_geometry_seed_but_not_barrier(
    tmp_path: Path, monkeypatch
) -> None:
    reaction, reactant, product = _endpoints(tmp_path)
    raw_output = """
 Northwest Computational Chemistry Package (NWChem) 7.2.3
 Total DFT energy = -10.0
 DFT-D3 Model
 Dispersion correction = -0.001
 @zts  4 0.003600 0.090000 -10 -9 -10 -9 100
"""
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_string.resolve_executable",
        lambda *args: "nwchem",
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_string.run_nwchem_input",
        _string_run_result(raw_output, returncode=124, final=False),
    )

    result = NWChemStringEngine().search_ts(
        reaction,
        reactant,
        product,
        {
            "path_strategy": "adaptive_double_ended_path",
            "endpoint_basin_status": "distinct_basin",
            "allow_subprocess": True,
            "allow_interrupted_path_seed": True,
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "interrupted",
    )

    path = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "reaction_path"
    )
    seed = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "species"
    )
    assert path.qc["path_energy_publishable"] is False
    assert path.qc["unconverged_path_used_as_seed_only"] is True
    assert path.qc["method_evidence_validated"] is False
    assert path.qc["method_identity_validated"] is True
    assert seed.qc["unconverged_path_used_as_seed_only"] is True
    assert seed.qc["method_evidence_validated"] is False
    assert seed.qc["saddle_seed_method_evidence_sufficient"] is True
    assert seed.data["state"] == "saddle_seed"
    assert not any(
        artifact.artifact_type == "reaction_validated"
        for artifact in result.artifacts
    )
