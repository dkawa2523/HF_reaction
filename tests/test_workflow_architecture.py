from pathlib import Path

from hfauto.chemistry.electronic_state import (
    coupled_multiplicities,
    select_coupled_multiplicities,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.reaction_plan import make_reaction_case_artifact
from hfauto.stages.ts_search import TSSearchStage
from hfauto.workflow.reaction_state import decide_reaction_case
from hfauto.workflow.ts_execution import (
    configured_ts_attempts,
    next_authorized_attempt,
    ts_attempt_key,
)

ROOT = Path(__file__).parents[1]
NEUTRAL = ROOT / "examples/m3_trimethylamine_hf2/neutral.xyz"
SHARED_PROTON = ROOT / "examples/m3_trimethylamine_hf2/shared_proton.xyz"


def _minimum(identifier: str, path: Path) -> Artifact:
    return Artifact(
        artifact_id=identifier,
        artifact_type="species_optimized",
        paths={"xyz": str(path)},
        data={
            "species_id": identifier.removeprefix("opt_"),
            "xyz_path": str(path),
            "charge": 0,
            "multiplicity": 1,
            "n_imag": 0,
        },
        qc={
            "minimum_accepted": True,
            "is_minimum": True,
            "scf_converged": True,
            "geometry_converged": True,
            "geometry_sane": True,
            "fallback_dummy": False,
        },
    )


def _generic_reaction() -> Artifact:
    return Artifact(
        artifact_id="rxn_generic",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn_generic",
            "reactant_species_id": "tma_hf2_neutral",
            "product_species_id": "tma_hf2_shared_proton",
            "mechanism_family": "generic_rearrangement",
            "reaction_coordinate": {
                "min_change": 0.05,
                "terms": [
                    {"kind": "distance", "atoms": [13, 14], "coefficient": 1.0},
                    {"kind": "distance", "atoms": [1, 13], "coefficient": -1.0},
                ],
            },
            "basin_assessment": {
                "status": "distinct_basin",
                "accepted": True,
            },
        },
        qc={"distinct_registry_basins": True},
    )


def test_spin_coupling_enumerates_valid_surfaces() -> None:
    assert coupled_multiplicities([1, 1]) == (1,)
    assert coupled_multiplicities([2, 2]) == (1, 3)
    assert coupled_multiplicities([3, 2]) == (2, 4)
    assert coupled_multiplicities([2, 2, 2]) == (2, 4)
    assert select_coupled_multiplicities([2, 2], policy="lowest") == (1,)


def test_ts_attempt_identity_includes_engine_settings() -> None:
    attempts = configured_ts_attempts(
        {
            "engine_order": [
                {"engine": "nwchem_saddle", "settings": {"driver_trust": 0.1}},
                {"engine": "nwchem_saddle", "settings": {"driver_trust": 0.05}},
            ]
        },
        {"functional": "pbe0"},
    )
    assert ts_attempt_key(attempts[0], None) != ts_attempt_key(attempts[1], None)
    first = next_authorized_attempt(attempts, None, set())
    assert first is not None
    second = next_authorized_attempt(attempts, None, {first[2]})
    assert second is not None
    assert second[1].method["driver_trust"] == 0.05


def test_ts_stage_advances_neb_string_and_saddle_in_one_run(
    tmp_path: Path, monkeypatch
) -> None:
    reaction = _generic_reaction()
    case = decide_reaction_case(
        reaction_id="rxn_generic",
        basin_assessment={"status": "distinct_basin", "accepted": True},
    )
    manifest = Manifest.new(run_id="test", stage="reaction-plan")
    manifest.extend(
        [
            reaction,
            _minimum("opt_tma_hf2_neutral", NEUTRAL),
            _minimum("opt_tma_hf2_shared_proton", SHARED_PROTON),
            make_reaction_case_artifact(
                case, parents=[reaction.artifact_id], created_by="test"
            ),
        ]
    )
    calls: list[tuple[str, str]] = []

    class Engine:
        def __init__(self, name: str) -> None:
            self.name = name

        def search_ts(self, reaction, reactant, product, method, workdir):
            strategy = str(method["path_strategy"])
            calls.append((self.name, strategy))
            if self.name == "nwchem_neb":
                diagnosis = "path_unresolved"
                artifact_type = "path_attempt"
                candidate = None
            elif self.name == "nwchem_string":
                diagnosis = "resolved_saddle_candidate"
                artifact_type = "path_attempt"
                candidate = {"xyz_path": str(SHARED_PROTON)}
            else:
                return [
                    Artifact(
                        artifact_id="validated_ts",
                        artifact_type="reaction_validated",
                        data={"reaction_id": "rxn_generic"},
                        qc={"ts_validated_by_frequency": True},
                    )
                ]
            return [
                Artifact(
                    artifact_id=f"{self.name}_attempt",
                    artifact_type=artifact_type,
                    data={
                        "reaction_id": "rxn_generic",
                        "strategy": strategy,
                        "diagnosis": diagnosis,
                        "candidate": candidate,
                    },
                )
            ]

    monkeypatch.setattr(
        "hfauto.stages.ts_search.get_ts_engine",
        lambda spec, **settings: Engine(
            str(spec.get("engine") if isinstance(spec, dict) else spec)
        ),
    )
    result = TSSearchStage().run(
        manifest,
        {
            "engine_order": [
                "nwchem_neb",
                "nwchem_string",
                "pysisyphus_saddle",
            ],
            "require_reaction_case": True,
            "require_validated_dft_minima": False,
        },
        StageContext(
            out_dir=tmp_path / "ts",
            run_id="test",
            global_config={"mode": "development"},
        ),
    )
    assert calls == [
        ("nwchem_neb", "double_ended_path"),
        ("nwchem_string", "adaptive_double_ended_path"),
        ("pysisyphus_saddle", "saddle_refinement"),
    ]
    final_case = result.latest_artifacts("reaction_case")[-1]
    assert final_case.data["state"] == "ts_validated"
    assert final_case.data["irc_allowed"] is True
