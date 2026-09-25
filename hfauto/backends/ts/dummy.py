from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.backends.ts.base import (
    IRCResult,
    TSSearchResult,
    make_ts_species_artifact,
    midpoint_ts_xyz,
)
from hfauto.chemistry.reaction_path_qc import (
    endpoint_pair_match_qc,
    estimate_reaction_mode_overlap,
    reaction_coordinate_vector,
)
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.frequency_qc import imaginary_frequency_cutoff_from_method
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.ids import path_token
from hfauto.core.qc import ts_qc
from hfauto.core.schemas.artifact import Artifact


class DummyTSEngine:
    """Deterministic TS/IRC backend for software-contract tests.

    All artifacts are explicitly labelled dummy/fallback so scientific quality
    tiers are not upgraded by offline pipeline runs.
    """

    name = "dummy"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def search_ts(self, reaction: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> TSSearchResult:
        imaginary_cutoff = imaginary_frequency_cutoff_from_method(method)
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        ts_id = str(reaction.data.get("reaction_id") or reaction.artifact_id).replace("rxn_", "ts_", 1)
        ts_xyz = midpoint_ts_xyz(
            reactant,
            product,
            wd / f"{path_token(ts_id, max_length=28)}.xyz",
        )
        ts_species = make_ts_species_artifact(
            reaction,
            reactant,
            product,
            ts_xyz,
            ts_id=ts_id,
            source=self.name,
            extra_qc={"endpoint_midpoint": True, "engine_is_dummy": True, "fallback_dummy": True, "scientific_use": "software_test_only"},
        )
        calc = DummyQMEngine().optimize_frequency(ts_species, {**method, "method_id": method.get("method_id", "dummy-ts")}, str(wd / "dummy_qm"))
        calc.method = {**(calc.method or {}), "engine": self.name, "stage": "ts-search", "task": "opt_freq"}
        # The deterministic dummy uses the exact reaction-coordinate vector as
        # synthetic mode input.  It remains explicitly non-scientific, while the
        # shared production gate exercises a literal vector projection.
        synthetic_mode = reaction_coordinate_vector(ts_species.data, ts_xyz)
        overlap, overlap_method = estimate_reaction_mode_overlap(
            ts_species.data,
            ts_xyz,
            calc.data.get("n_imag"),
            calc.data.get("imag_freq_cm1"),
            mode_displacements=synthetic_mode,
            imaginary_frequency_cutoff_cm1=imaginary_cutoff,
        )
        calc.qc.update({"mode_overlap_score": overlap, "mode_overlap_method": overlap_method, "imag_mode_matches_reaction_coordinate": overlap >= 0.7})
        calc.qc.update(
            ts_qc(
                calc.data.get("n_imag"),
                calc.data.get("imag_freq_cm1"),
                overlap,
                imaginary_frequency_cutoff_cm1=imaginary_cutoff,
            )
        )
        calc.qc.update({"engine_is_dummy": True, "fallback_dummy": True, "scientific_use": "software_test_only"})
        calc.data.update({"mode_overlap_score": overlap, "mode_overlap_method": overlap_method, "ts_backend": self.name})
        rv = Artifact(
            artifact_id=reaction.artifact_id + "_with_ts",
            artifact_type="reaction_validated",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id, ts_species.artifact_id, calc.artifact_id],
            data={**reaction.data, "ts_species_id": ts_species.artifact_id, "ts_calc_id": calc.artifact_id, "ts_backend": self.name},
            qc={
                "ts_engine": self.name,
                "ts_found": True,
                "ts_validated_by_frequency": bool(calc.qc.get("ts_validated_by_frequency")),
                "n_imag": calc.data.get("n_imag"),
                "imag_freq_cm1": calc.data.get("imag_freq_cm1"),
                "mode_overlap_score": overlap,
                "mode_overlap_method": overlap_method,
                "ts_backend_is_dummy": True,
                "fallback_dummy": True,
                "scientific_use": "software_test_only",
            },
        )
        path_art = Artifact(
            artifact_id="reaction_path_"
            + fingerprint_dict({"reaction": reaction.artifact_id, "engine": self.name}),
            artifact_type="reaction_path",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            paths={"ts_xyz": str(ts_xyz)},
            method={"engine": self.name, "task": "ts_guess"},
            data={
                "reaction_id": reaction.data.get("reaction_id", reaction.artifact_id),
                "ts_species_id": ts_species.artifact_id,
                "reaction_path": analyze_reaction_path(
                    reaction_id=str(
                        reaction.data.get("reaction_id", reaction.artifact_id)
                    ),
                    engine=self.name,
                    comments=[],
                    converged=False,
                    energy_source="dummy_midpoint_has_no_path_energies",
                ).model_dump(),
            },
            qc={
                "neb_profile_classification": "missing_profile",
                "fallback_dummy": True,
                "scientific_use": "software_test_only",
            },
        )
        path_validated = Artifact(
            artifact_id="path_validated_" + fingerprint_dict({"reaction": reaction.artifact_id, "engine": self.name}),
            artifact_type="reaction_path_validated",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id, ts_species.artifact_id, calc.artifact_id],
            paths={"ts_xyz": str(ts_xyz), "reactant_xyz": reactant.data.get("xyz_path") or reactant.paths.get("xyz") or "", "product_xyz": product.data.get("xyz_path") or product.paths.get("xyz") or ""},
            method={"engine": self.name, "task": "path_validation"},
            data={"reaction_id": reaction.data.get("reaction_id", reaction.artifact_id), "ts_species_id": ts_species.artifact_id, "ts_calc_id": calc.artifact_id},
            qc={"path_validated": bool(calc.qc.get("ts_validated_by_frequency")), "fallback_dummy": True, "scientific_use": "software_test_only"},
        )
        artifacts = [path_art, path_validated, ts_species, calc, rv]
        return TSSearchResult(artifacts=artifacts, ts_species_id=ts_species.artifact_id, ts_calc_id=calc.artifact_id, record=rv.data)

    def run_irc(self, reaction: Artifact, ts_species: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> IRCResult:
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        forward = wd / "irc_forward_endpoint.xyz"
        backward = wd / "irc_backward_endpoint.xyz"
        prod = read_xyz(product.data.get("xyz_path") or product.paths.get("xyz"))
        reac = read_xyz(reactant.data.get("xyz_path") or reactant.paths.get("xyz"))
        write_xyz(XYZ(list(prod.symbols), prod.coords.copy(), "state=irc_forward dummy_expected_product"), forward)
        write_xyz(XYZ(list(reac.symbols), reac.coords.copy(), "state=irc_backward dummy_expected_reactant"), backward)
        qc = endpoint_pair_match_qc(
            forward,
            backward,
            reactant.data.get("xyz_path") or reactant.paths.get("xyz"),
            product.data.get("xyz_path") or product.paths.get("xyz"),
            ts_species.data,
            float(method.get("endpoint_rmsd_threshold_A", 0.75)),
            q_tolerance_A=float(method.get("endpoint_q_tolerance_A", 0.30)),
            require_identity_invariant_geometry=bool(
                method.get("require_identity_invariant_endpoint_match", True)
            ),
            permutation_rmsd_threshold_A=float(
                method.get("endpoint_permutation_rmsd_threshold_A", 0.20)
            ),
            distance_spectrum_threshold_A=float(
                method.get("endpoint_distance_spectrum_threshold_A", 0.08)
            ),
        )
        rec = {
            "reaction_id": reaction.data.get("reaction_id", reaction.artifact_id.replace("_with_ts", "")),
            "ts_species_id": ts_species.artifact_id,
            "forward_endpoint_xyz": str(forward),
            "backward_endpoint_xyz": str(backward),
            "irc_backend": self.name,
            "real_irc_executed": False,
            "fallback_dummy": True,
            "scientific_use": "software_test_only_not_irc",
            **qc,
        }
        artifact = Artifact(
            artifact_id=str(reaction.artifact_id).replace("_with_ts", "_irc"),
            artifact_type="irc",
            parents=[reaction.artifact_id, ts_species.artifact_id, reactant.artifact_id, product.artifact_id],
            paths={"forward_xyz": str(forward), "backward_xyz": str(backward), "endpoint_a_xyz": str(forward), "endpoint_b_xyz": str(backward)},
            data=rec,
            method={"engine": self.name, "task": "irc"},
            qc={"irc_validated": bool(qc.get("irc_validated")), "irc_backend_is_dummy": True, "fallback_dummy": True, "real_irc_executed": False},
        )
        return IRCResult(artifacts=[artifact], record=rec)
