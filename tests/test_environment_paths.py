from hfauto.core.environment import (
    _configured_external_paths,
    _configured_external_requirements,
    _configured_python_requirements,
)


def test_production_preflight_uses_stage_owned_explicit_executables() -> None:
    config = {
        "stages": [
            {
                "name": "ts-search",
                "engine_order": ["nwchem_saddle"],
                "settings": {
                    "executable": "/opt/nwchem",
                    "mpi_executable": "/opt/mpirun",
                },
            },
            {
                "name": "irc",
                "engine": "pysisyphus",
                "settings": {
                    "program": "nwchem",
                    "pysis_executable": "/opt/pysis",
                },
            },
        ]
    }

    assert _configured_external_paths(config) == {
        "nwchem": "/opt/nwchem",
        "mpirun": "/opt/mpirun",
        "pysis": "/opt/pysis",
    }


def test_production_preflight_detects_crest_nci_requirements() -> None:
    config = {
        "stages": [
            {
                "name": "build-complexes",
                "crest_nci_settings": {
                    "executable": "/opt/crest",
                    "xtb_executable": "/opt/xtb",
                },
                "compositions": [
                    {
                        "components": [
                            {"selector": "base", "count": 1},
                            {"selector": "acid", "count": 2},
                        ]
                    }
                ],
            }
        ]
    }

    assert _configured_external_requirements(config) == {"crest", "xtb"}
    assert _configured_external_paths(config) == {
        "crest": "/opt/crest",
        "xtb": "/opt/xtb",
    }


def test_conformer_preflight_reads_stage_owned_crest_paths() -> None:
    config = {
        "stages": [
            {
                "name": "conformers",
                "engine": "crest",
                "executable": "/opt/crest",
                "xtb_executable": "/opt/xtb",
            }
        ]
    }

    assert _configured_external_requirements(config) == {"crest", "xtb"}
    assert _configured_external_paths(config) == {
        "crest": "/opt/crest",
        "xtb": "/opt/xtb",
    }


def test_reaction_discovery_preflight_declares_python_dependencies() -> None:
    config = {
        "stages": [
            {"name": "build-complexes"},
            {"name": "explore-reactions", "engine": "readuct"},
        ]
    }

    assert _configured_python_requirements(config) == {
        "hfauto",
        "numpy",
        "pandas",
        "pydantic",
        "rdkit",
        "scine_readuct",
        "scine_xtb_wrapper",
    }

    config["stages"][1]["engine"] = "future_backend"
    assert "scine_readuct" not in _configured_python_requirements(config)
    assert "scine_xtb_wrapper" not in _configured_python_requirements(config)
