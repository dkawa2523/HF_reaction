import json
from pathlib import Path

from rdkit import Chem

from hfauto.chemistry.site_detection import detect_basic_sites_from_smiles
from hfauto.core.config import load_yaml

ROOT = Path(__file__).parents[1]


def test_amine_panel_definition_and_generated_sdfs():
    panel = json.loads((ROOT / "examples/production/amine_hf_panel.json").read_text(encoding="utf-8"))
    entries = panel["entries"]
    assert len(entries) == 10
    assert len({entry["candidate_id"] for entry in entries}) == 10
    assert sum("pilot" in entry["subsets"] for entry in entries) == 4

    for subset, expected in [("pilot", 4), ("full", 10)]:
        supplier = Chem.SDMolSupplier(
            str(ROOT / f"examples/production/amine_hf_{subset}.sdf"), removeHs=False
        )
        molecules = [mol for mol in supplier if mol is not None]
        assert len(molecules) == expected
        assert all(mol.GetNumConformers() == 1 and mol.GetConformer().Is3D() for mol in molecules)
        assert all(mol.HasProp("candidate_id") and mol.GetProp("panel_subset") == subset for mol in molecules)


def test_panel_exercises_site_selectivity_and_sterics():
    imidazole = detect_basic_sites_from_smiles("c1ncc[nH]1")
    assert sum(not site["excluded"] for site in imidazole) == 1
    assert any(site["excluded"] and site["site_type"] == "pyrrole_N_excluded" for site in imidazole)

    dabco = detect_basic_sites_from_smiles("C1CN2CCN1CC2")
    assert sum(not site["excluded"] for site in dabco) == 2
    triethylamine = detect_basic_sites_from_smiles("CCN(CC)CC")
    assert triethylamine[0]["local_environment"]["steric_proxy"] > 0.7


def test_tma_hf2_study_is_bounded_and_uses_real_diffuse_dft():
    config = load_yaml(
        ROOT / "configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml"
    )
    stages = {stage["name"]: stage for stage in config["stages"]}
    assert stages["dft-minima"]["settings"]["basis"].endswith("d")
    assert stages["dft-minima"]["settings"]["optimization_convergence"] == "default"
    assert stages["ts-search"]["settings"]["basis"] == "def2-svpd"
    assert stages["irc"]["engine"] == "pysisyphus"
    assert stages["generate-reactions"]["max_total_trials"] == 12
    assert stages["explore-reactions"]["max_total_attempts"] == 24
    assert stages["generate-reactions"]["required_source_data_keys"] == [
        "crest_nci_conformer_id"
    ]
    assert stages["preopt"]["require_real_qm"] is True
    assert stages["preopt"]["require_connectivity_retention"] is True
    assert stages["preopt"]["required_source_data_keys"] == [
        "crest_nci_conformer_id"
    ]
    assert stages["preopt"]["checkpoint_each_job"] is True
    assert stages["preopt"]["resume_completed_jobs"] is True
    assert stages["relaxation-discovery"]["include_crest_topology_stops"] is True
    assert stages["dft-minima"]["candidate_endpoint_limit"] == 3
    assert stages["dft-minima"]["candidate_duplicate_rmsd_A"] == 0.05
    assert all(
        not (stage.get("settings", {}) or {}).get("fallback_to_dummy", False)
        for stage in config["stages"]
    )


def test_tma_hf2_keeps_minima_ts_and_irc_on_one_pes():
    config = load_yaml(
        ROOT / "configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml"
    )
    stages = {stage["name"]: stage for stage in config["stages"]}
    expected = {
        "functional": "pbe0",
        "basis": "def2-svpd",
        "disp_vdw": 3,
        "required_program_version": "7.2.3",
    }
    for stage_name in ("dft-minima", "ts-search", "irc"):
        settings = stages[stage_name]["settings"]
        observed = {key: settings[key] for key in expected}
        observed["required_program_version"] = str(
            observed["required_program_version"]
        )
        assert observed == expected
    assert stages["minimum-registry"]["require_real_qm"] is True
    assert stages["reaction-plan"]["max_path_attempts"] == 3


def test_standard_discovery_pipelines_normalize_unbiased_relaxations() -> None:
    for config_name in (
        "molecular_reaction_discovery_nwchem.yaml",
        "m3_trimethylamine_hf2_nwchem.yaml",
    ):
        config = load_yaml(ROOT / "configs/pipelines" / config_name)
        names = [stage["name"] for stage in config["stages"]]
        preopt = names.index("preopt")
        relaxation = names.index("relaxation-discovery")
        generation = names.index("generate-reactions")
        assert preopt < relaxation < generation

        stages = {stage["name"]: stage for stage in config["stages"]}
        assert stages["preopt"]["require_real_qm"] is True
        assert stages["preopt"]["require_connectivity_retention"] is True
        assert stages["preopt"]["checkpoint_each_job"] is True
        assert stages["preopt"]["resume_completed_jobs"] is True
        assert stages["relaxation-discovery"]["include_crest_topology_stops"] is True


def test_production_nwchem_paths_preserve_interrupted_geometry_as_seed_only() -> None:
    config_names = (
        "molecular_reaction_discovery_nwchem.yaml",
        "m3_trimethylamine_hf2_nwchem.yaml",
        "tma_hf2_validated_minima_path.yaml",
        "tma_hf2_recovered_path.yaml",
        "hcn_alternative_pes_validation.yaml",
        "hono_isomerization_nwchem.yaml",
        "formaldehyde_hydroxymethylene_nwchem.yaml",
    )
    for config_name in config_names:
        config = load_yaml(ROOT / "configs/pipelines" / config_name)
        ts = next(stage for stage in config["stages"] if stage["name"] == "ts-search")
        assert ts["settings"]["allow_interrupted_path_seed"] is True


def test_pilot_campaign_covers_four_chemistries_and_hf1_to_hf3() -> None:
    config = load_yaml(ROOT / "configs/pipelines/amine_hf_pilot_hf1_hf3.yaml")
    complex_stages = [
        stage for stage in config["stages"] if stage["name"] == "build-complexes"
    ]
    campaign = complex_stages[-1]
    compositions = campaign["compositions"]

    assert campaign["max_components"] == 5
    assert len(compositions) == 15
    assert {item["composition_id"].split("_")[0] for item in compositions} == {
        "NH3",
        "TMA",
        "ANL",
        "PYR",
    }
    assert {
        sum(int(component.get("count", 1)) for component in item["components"])
        for item in compositions
    } == {2, 3, 4, 5}
    water_series = [
        item for item in compositions if item["composition_id"].endswith("_H2O")
    ]
    assert {item["composition_id"] for item in water_series} == {
        "TMA_HF1_H2O",
        "TMA_HF2_H2O",
        "TMA_HF3_H2O",
    }
    stages = {stage["name"]: stage for stage in config["stages"]}
    assert stages["preopt"]["require_connectivity_retention"] is True
    assert stages["preopt"]["required_source_data_keys"] == [
        "crest_nci_conformer_id"
    ]
    assert stages["relaxation-discovery"]["include_crest_topology_stops"] is True
    assert stages["generate-reactions"]["max_total_trials"] == 30
    assert stages["explore-reactions"]["max_total_attempts"] == 60
    assert campaign["checkpoint_each_definition"] is True
    assert campaign["resume_completed_definitions"] is True


def test_pilot_dft_followup_is_composition_balanced_and_end_to_end() -> None:
    config = load_yaml(
        ROOT / "configs/pipelines/amine_hf_pilot_dft_followup.yaml"
    )
    stages = {stage["name"]: stage for stage in config["stages"]}
    names = [stage["name"] for stage in config["stages"]]

    assert names == [
        "dft-minima",
        "minimum-registry",
        "connect-minima",
        "reaction-plan",
        "ts-search",
        "irc",
        "reaction-classify",
        "thermo",
        "basin-populations",
        "reaction-rank",
        "discovery-audit",
    ]
    dft = stages["dft-minima"]
    assert dft["candidate_endpoint_limit"] == 4
    assert dft["max_species"] == 8
    assert dft["checkpoint_each_job"] is True
    assert dft["resume_completed_jobs"] is True
    assert stages["ts-search"]["max_reactions"] == 4
    assert stages["ts-search"]["max_total_attempts"] == 12
    expected_surface = {
        "functional": "pbe0",
        "basis": "def2-svpd",
        "disp_vdw": 3,
        "required_program_version": "7.2.3",
        "grid": "xfine",
        "scf_energy_tolerance": 1.0e-8,
    }
    for stage_name in ("dft-minima", "ts-search", "irc"):
        observed = {
            key: str(stages[stage_name]["settings"][key])
            if key == "required_program_version"
            else stages[stage_name]["settings"][key]
            for key in expected_surface
        }
        assert observed == expected_surface


def test_water_extension_contains_only_incremental_microsolvated_compositions() -> None:
    config = load_yaml(ROOT / "configs/pipelines/amine_hf_water_extension.yaml")
    assert [stage["name"] for stage in config["stages"]] == ["build-complexes"]
    campaign = config["stages"][0]
    compositions = campaign["compositions"]

    assert {item["composition_id"] for item in compositions} == {
        "TMA_HF1_H2O",
        "TMA_HF2_H2O",
        "TMA_HF3_H2O",
    }
    assert all(
        sum(component.get("role") == "relay" for component in item["components"])
        == 1
        for item in compositions
    )
    assert campaign["checkpoint_each_definition"] is True
