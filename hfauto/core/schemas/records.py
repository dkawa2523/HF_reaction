from __future__ import annotations

"""Typed record contracts used inside Artifact.data.

The workflow still passes generic Artifacts between stages, but these models document
and validate the canonical payloads that reviewers should expect. Stages may write
additional backend-specific fields under ``extras`` without changing the core schema.
"""

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


class MoleculeRecord(BaseModel):
    mol_id: str
    source_sdf_index: int
    name: str
    canonical_smiles: Optional[str] = None
    isomeric_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    formula: Optional[str] = None
    exact_mw: Optional[float] = None
    formal_charge: Optional[int] = 0
    multiplicity: int = 1
    num_atoms: Optional[int] = None
    num_fragments: int = 1
    has_3d: bool = False
    sdf_props: dict[str, Any] = Field(default_factory=dict)
    ingest_qc: dict[str, Any] = Field(default_factory=dict)
    extras: dict[str, Any] = Field(default_factory=dict)


class SiteRecord(BaseModel):
    site_id: str
    mol_id: str
    atom_index: int
    site_type: str
    site_smarts: str
    priority: int
    site_confidence: float
    excluded: bool = False
    exclude_reason: Optional[str] = None
    local_environment: dict[str, Any] = Field(default_factory=dict)
    basicity_proxy: dict[str, Any] = Field(default_factory=dict)


class ConformerRecord(BaseModel):
    conformer_id: str
    mol_id: str
    source: str
    relative_energy_kcal_mol: float
    boltzmann_weight_298K: float
    xyz_path: str
    selected_for_hf_build: bool = True
    rmsd_cluster_id: Optional[str] = None
    extras: dict[str, Any] = Field(default_factory=dict)


class SpeciesRecord(BaseModel):
    species_id: str
    mol_id: str
    site_id: Optional[str] = None
    conformer_id: Optional[str] = None
    state: Literal[
        "molecule",
        "candidate",
        "hf_cluster",
        "reactant_complex",
        "ion_pair",
        "transition_state",
        "probe_reactant",
        "probe_product",
    ]
    hf_n: int = 0
    charge: int = 0
    multiplicity: int = 1
    xyz_path: str
    atom_order_key: Optional[str] = None
    components: list[dict[str, Any]] = Field(default_factory=list)


class ReactionRecord(BaseModel):
    reaction_id: str
    reaction_type: str
    mol_id: str
    site_id: Optional[str] = None
    hf_n: int = 1
    reactant_species_id: str
    product_species_id: str
    ts_species_id: Optional[str] = None
    reaction_coordinate: dict[str, Any] = Field(default_factory=dict)


class CalculationRecord(BaseModel):
    calc_id: str
    species_id: str
    task: str
    engine: str
    method_id: str
    status: str = "success"
    energies: dict[str, Optional[float]] = Field(default_factory=dict)
    frequencies: dict[str, Optional[float]] = Field(default_factory=dict)
    qc: dict[str, Any] = Field(default_factory=dict)
    paths: dict[str, str] = Field(default_factory=dict)


class ThermoRecord(BaseModel):
    reaction_id: str
    mol_id: str
    site_id: Optional[str] = None
    hf_n: int
    reaction_type: str
    T_K: float = 298.15
    p_standard_bar: float = 1.0
    delta_G_assoc_kcal_mol: Optional[float] = None
    delta_G_assoc_standard_kcal_mol: Optional[float] = None
    delta_G_assoc_pressure_corrected_kcal_mol: Optional[float] = None
    delta_G_ionpair_kcal_mol: Optional[float] = None
    delta_G_act_kcal_mol: Optional[float] = None
    K_assoc_standard: Optional[float] = None
    K_assoc_process_adjusted: Optional[float] = None
    K_ionpair: Optional[float] = None
    thermo_backend: Optional[str] = None
    quasi_rrho_applied: bool = False
    frequency_scale_factor: Optional[float] = None
    quality_tier: str = "Q0"
    confidence_score: Optional[float] = None
    extras: dict[str, Any] = Field(default_factory=dict)


class KineticsRecord(BaseModel):
    reaction_id: str
    mol_id: str
    site_id: Optional[str] = None
    hf_n: int
    T_K: float
    delta_G_act_kcal_mol: Optional[float] = None
    imag_freq_cm1: Optional[float] = None
    k_TST_s_inv: Optional[float] = Field(default=None, alias="k_TST_s-1")
    k_corrected_s_inv: Optional[float] = Field(default=None, alias="k_corrected_s-1")
    tunneling_model: str = "none"
    transmission_coefficient: float = 1.0
    kinetics_quality: Optional[str] = None
    extras: dict[str, Any] = Field(default_factory=dict)


class DescriptorRecord(BaseModel):
    species_id: str
    mol_id: str
    site_id: Optional[str] = None
    hf_n: int
    r_HF_A: Optional[float] = None
    delta_r_HF_A: Optional[float] = None
    nu_HF_cm1: Optional[float] = None
    delta_nu_HF_cm1: Optional[float] = None
    B_H_distance_A: Optional[float] = None
    B_H_F_angle_deg: Optional[float] = None
    descriptor_quality: str = "geometry_only"


class SpeciesThermoRecord(BaseModel):
    species_id: str
    state: Optional[str] = None
    mol_id: Optional[str] = None
    site_id: Optional[str] = None
    hf_n: Optional[int] = None
    T_K: float
    G_hartree: float
    G_standard_hartree: Optional[float] = None
    pressure_bar: Optional[float] = None
    pressure_correction_hartree: Optional[float] = None
    thermal_backend: str = "simple"
    thermal_correction_source: Optional[str] = None
    quasi_rrho_applied: bool = False
    low_frequency_count: Optional[int] = None
    extras: dict[str, Any] = Field(default_factory=dict)


class KineticsRecord(BaseModel):
    reaction_id: str
    mol_id: Optional[str] = None
    site_id: Optional[str] = None
    hf_n: Optional[int] = None
    T_K: float
    delta_G_act_kcal_mol: float
    k_TST_s_1: Optional[float] = Field(default=None, alias="k_TST_s-1")
    k_corrected_s_1: Optional[float] = Field(default=None, alias="k_corrected_s-1")
    tunneling_model: str = "none"
    tunneling_factor: float = 1.0
    kinetics_quality: str = "unknown"
    extras: dict[str, Any] = Field(default_factory=dict)
