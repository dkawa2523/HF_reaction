from pathlib import Path

import yaml

from hfauto.chemistry.reactions import validate_reaction_coordinate_between_geometries

ROOT = Path(__file__).parents[1]
CONFIG = ROOT / "configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml"
REACTANT = ROOT / "examples/m3_trimethylamine_hf2/neutral.xyz"
PRODUCT = ROOT / "examples/m3_trimethylamine_hf2/shared_proton.xyz"
RECOVERY_CONFIG = ROOT / "configs/pipelines/tma_hf2_recovered_path.yaml"


def _reaction_coordinate() -> dict:
    return {
        "min_change": 0.05,
        "terms": [
            {
                "kind": "distance",
                "atoms": [13, 14],
                "coefficient": 1.0,
                "label": "transfer_H_to_accepting_F",
            },
            {
                "kind": "distance",
                "atoms": [1, 13],
                "coefficient": -1.0,
                "label": "base_N_to_transfer_H",
            },
        ],
    }


def test_tma_hf2_reviewed_minima_have_declared_coordinate_progress() -> None:
    qc = validate_reaction_coordinate_between_geometries(
        {"reaction_coordinate": _reaction_coordinate()}, REACTANT, PRODUCT
    )

    assert qc["accepted"] is True
    assert qc["progress"] > 0.05


def test_tma_hf2_pipeline_keeps_one_pes_and_validates_endpoint_basins() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    stages = {stage["name"]: stage for stage in config["stages"]}
    minimum = stages["dft-minima"]["settings"]
    neb = stages["ts-search"]["settings"]
    irc = stages["irc"]["settings"]

    for settings in (minimum, neb, irc):
        assert settings["functional"] == "pbe0"
        assert settings["basis"] == "def2-svpd"
        assert settings["disp_vdw"] == 3
        assert str(settings["required_program_version"]) == "7.2.3"
    assert neb["nbeads"] == 9
    assert stages["ts-search"]["engine_order"] == [
        "nwchem_neb",
        "nwchem_string",
        "nwchem_saddle",
    ]
    assert neb["string_nbeads"] == 9
    assert neb["ts_mode_overlap_threshold"] == 0.50
    assert stages["ts-search"]["require_reaction_case"] is True
    assert stages["irc"]["require_reaction_case"] is True
    assert irc["optimize_irc_endpoints"] is True
    assert irc["require_optimized_irc_endpoints"] is True
    assert stages["generate-reactions"]["max_total_trials"] == 12
    assert stages["explore-reactions"]["engine"] == "readuct"
    assert stages["minimum-registry"]["require_real_qm"] is True


def test_tma_hf2_recovery_is_bounded_and_interrupted_paths_are_seed_only() -> None:
    config = yaml.safe_load(RECOVERY_CONFIG.read_text(encoding="utf-8"))
    stages = {stage["name"]: stage for stage in config["stages"]}
    ts = stages["ts-search"]

    assert ts["max_total_attempts"] == 4
    assert ts["max_stage_walltime_s"] == 14400
    assert ts["settings"]["string_maxiter"] == 12
    assert ts["settings"]["allow_interrupted_path_seed"] is True
    assert stages["irc"]["require_reaction_case"] is True
