
from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from hfauto.core.artifacts import canonical_species_id
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.units import kcal_mol_to_hartree


class DummyQMEngine:
    """Deterministic fake QM backend for testing DAGs and data contracts.

    The energies are deliberately constructed so association, proton transfer and TS
    outputs are non-trivial but stable. They are not chemically predictive.
    """

    name = "dummy"

    SITE_STRENGTH: ClassVar[dict[str, float]] = {
        "amidine_like": 1.15,
        "aliphatic_amine": 1.00,
        "pyridine_like": 0.78,
        "imine": 0.65,
        "sulfide": 0.45,
        "ether_oxygen": 0.32,
        "carbonyl_oxygen": 0.25,
    }

    def __init__(self, **kwargs):
        self.config = kwargs

    def _strength(self, species: Artifact) -> float:
        return float(self.SITE_STRENGTH.get(species.data.get("site_type"), 0.5))

    def _assoc_kcal(self, species: Artifact) -> float:
        hf_n = int(species.data.get("hf_n", 0) or 0)
        strength = self._strength(species)
        return -(4.0 + 4.5 * strength) * max(1, hf_n) ** 0.75

    def _barrier_kcal(self, species: Artifact) -> float:
        hf_n = int(species.data.get("hf_n", 0) or 0)
        strength = self._strength(species)
        return max(0.5, 12.0 - 7.0 * strength - 2.0 * max(0, hf_n - 1))

    def _ionpair_delta_kcal(self, species: Artifact) -> float:
        hf_n = int(species.data.get("hf_n", 0) or 0)
        strength = self._strength(species)
        return 5.0 - 8.0 * strength - 1.5 * max(0, hf_n - 1)

    def _fake_energy(self, species: Artifact, method: dict, task: str) -> float:
        state = species.data.get("state", "")
        hf_n = int(species.data.get("hf_n", 0) or 0)
        # Reference species have additivity; complexes are stabilized so ΔG_assoc is finite.
        if state in {"bare_candidate", "candidate", "isolated_candidate"}:
            base = -100.0
        elif state == "hf_cluster":
            base = -0.010 * hf_n
        elif state == "reactant_complex":
            base = -100.0 - 0.010 * hf_n + kcal_mol_to_hartree(self._assoc_kcal(species))
        elif state == "ion_pair":
            rc = -100.0 - 0.010 * hf_n + kcal_mol_to_hartree(self._assoc_kcal(species))
            base = rc + kcal_mol_to_hartree(self._ionpair_delta_kcal(species))
        elif state == "transition_state":
            rc = -100.0 - 0.010 * hf_n + kcal_mol_to_hartree(self._assoc_kcal(species))
            base = rc + kcal_mol_to_hartree(self._barrier_kcal(species))
        else:
            base = -100.0 - 0.01 * hf_n
        if task == "single_point":
            base -= 0.001
        return base

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        Path(workdir).mkdir(parents=True, exist_ok=True)
        task = str(method.get("task", "opt_freq"))
        energy = self._fake_energy(species, method, task)
        calc_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": task})
        out_path = str(Path(workdir) / "dummy_qm.out")
        Path(out_path).write_text(f"Dummy QM output\nE={energy}\n", encoding="utf-8")
        n_imag = 1 if species.data.get("state") == "transition_state" else 0
        imag = -1000.0 if n_imag == 1 else None
        hf_n = int(species.data.get("hf_n", 0) or 0)
        strength = self._strength(species)
        hf_stretch = None
        if species.data.get("state") in {"reactant_complex", "ion_pair", "transition_state"}:
            hf_stretch = 3961.0 - 45.0 * hf_n - 120.0 * strength
        return Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"output": out_path, "final_xyz": species.data.get("xyz_path", "")},
            method={"engine": self.name, "method_id": method.get("method_id", "dummy"), "task": task},
            data={
                "calc_id": calc_id,
                "species_id": canonical_species_id(species),
                "task": task,
                "engine": self.name,
                "method_id": method.get("method_id", "dummy"),
                "electronic_energy_hartree": energy,
                "zpe_hartree": 0.01 + 0.001 * hf_n,
                "gibbs_298K_hartree": energy + 0.02 + 0.001 * hf_n,
                "enthalpy_298K_hartree": energy + 0.015 + 0.001 * hf_n,
                "n_imag": n_imag,
                "imag_freq_cm1": imag,
                "frequencies_cm1": ([imag] if imag is not None else []) + [35.0, 75.0, 140.0, 520.0, 1100.0] + ([hf_stretch] if hf_stretch is not None else [3600.0]),
                "lowest_freq_cm1": imag if imag is not None else 35.0,
                "hf_stretch_cm1": hf_stretch,
                "dummy_assoc_kcal_mol": self._assoc_kcal(species) if species.data.get("state") == "reactant_complex" else None,
                "dummy_barrier_kcal_mol": self._barrier_kcal(species) if species.data.get("state") == "transition_state" else None,
            },
            qc={"scf_converged": True, "geometry_converged": True, "n_imag": n_imag, "engine_is_dummy": True, "scientific_use": "software_test_only"},
        )

    def single_point(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        Path(workdir).mkdir(parents=True, exist_ok=True)
        task = "single_point"
        energy = self._fake_energy(species, method, task)
        calc_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": task})
        out_path = str(Path(workdir) / "dummy_sp.out")
        Path(out_path).write_text(f"Dummy high-level SP output\nE={energy}\n", encoding="utf-8")
        return Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"output": out_path, "geometry": species.data.get("xyz_path", "")},
            method={"engine": self.name, "method_id": method.get("method_id", "dummy"), "task": task},
            data={
                "calc_id": calc_id,
                "species_id": canonical_species_id(species),
                "task": task,
                "engine": self.name,
                "method_id": method.get("method_id", "dummy"),
                "electronic_energy_hartree": energy,
            },
            qc={"scf_converged": True, "engine_is_dummy": True},
        )
