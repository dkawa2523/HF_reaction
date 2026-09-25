"""Artifact payload records (design §5.3). Payloads are a union discriminated by ``kind``."""

from __future__ import annotations

from enum import StrEnum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field

from hfauto.core.evidence import Evidence, FileRef, Geometry

_FROZEN = ConfigDict(frozen=True, extra="forbid")


class ArtifactType(StrEnum):
    SPECIES = "species"
    CALCULATION = "calculation"
    MINIMUM = "minimum"
    DISCOVERY = "discovery"
    REACTION = "reaction"
    SPECIES_THERMO = "species_thermo"
    REACTION_THERMO = "reaction_thermo"
    REPORT = "report"


class SpeciesRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["species"] = "species"
    species_id: str
    composition_id: str  # chemistry.xyz.composition_key
    formula: str
    charge: int
    multiplicity: int
    geometry: Geometry
    source: Literal[
        "input", "conformer", "placement", "crest_topology", "discovery",
        "mode_follow", "connection", "intermediate",
    ]
    state_label: str  # topology.state_label
    energy_hartree: float | None = None  # low-level reference energy
    level_key: str | None = None


class MinimumRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["minimum"] = "minimum"
    minimum_id: str
    basin_id: str
    composition_id: str
    species_id: str  # species of the representative structure
    tier: Literal["screen", "dft"]
    level_key: str  # Level.full_key()
    opt_calc: str  # calculation artifact ids
    freq_calc: str
    energy_hartree: float
    state_label: str
    n_fragments: int
    members: tuple[str, ...] = ()  # species that fell into this basin (collapsed seeds included)
    notes: tuple[str, ...] = ()


class ReactionTrial(BaseModel):
    model_config = _FROZEN
    trial_id: str
    source_minimum: str
    kind: Literal["polar_h", "h_shift", "heavy_bond", "association"]
    mechanism: Literal["nt2", "afir"]
    associations: tuple[tuple[int, int], ...] = ()
    dissociations: tuple[tuple[int, int], ...] = ()
    perturbed: bool = False  # linear molecule bent or randomly displaced


class DiscoveryRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["discovery"] = "discovery"
    discovery_id: str
    source_minimum: str
    mechanism: Literal["nt2", "afir", "relaxation", "mode_follow"]
    trial: ReactionTrial | None = None
    outcome: Literal["product", "negative", "failed"]
    reason: str | None = None
    product_species: str | None = None
    ts: Geometry | None = None  # low-level TS with one projected imaginary mode
    ts_imag_cm1: float | None = None
    barrier_kj_mol: float | None = None  # evaluated at 300 K
    reaction_kj_mol: float | None = None
    electronic_temperature_K: float = 300.0


class StoichTerm(BaseModel):
    model_config = _FROZEN
    composition_id: str
    coefficient: int


class CoordinateTerm(BaseModel):
    model_config = _FROZEN
    kind: Literal["distance", "angle", "dihedral"]
    atoms: tuple[int, ...]
    coefficient: float = 1.0


class BarrierVerdict(BaseModel):
    model_config = _FROZEN
    verdict: Literal["proceed", "barrierless", "negative_evidence", "unavailable"]
    max_rel_low_kcal: float | None = None
    max_rel_dft_kcal: float | None = None
    n_dft_points: int = 0
    max_node_spacing_A: float | None = None  # resolution of the barrierless claim
    below_zpe: bool = False
    seed: Geometry | None = None  # proceed only: low-level TS or DFT HEI node
    reasons: tuple[str, ...] = ()


class SaddleClaim(BaseModel):
    model_config = _FROZEN
    saddle_calc: str
    freq_calc: str
    imag_cm1: float
    energy_hartree: float
    notes: tuple[str, ...] = ()


class ConnectionClaim(BaseModel):
    model_config = _FROZEN
    method: Literal["qrc"] = "qrc"
    side_calcs: tuple[str, str]
    minima: tuple[str, str]
    amplitude_A: float
    notes: tuple[str, ...] = ()


class CaseOutcome(StrEnum):
    ELEMENTARY_STEP = "elementary_step"
    DEGENERATE = "degenerate_rearrangement"
    REASSIGNED = "reassigned_step"
    MULTI_STEP = "multi_step"
    BARRIERLESS = "barrierless_at_resolution"
    SAME_BASIN = "same_basin"
    NO_PRODUCT = "no_product_basin"
    OUT_OF_WINDOW = "out_of_window"
    UNRESOLVED = "unresolved_within_budget"
    BLOCKED = "blocked_upstream"


class ReactionRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["reaction"] = "reaction"
    reaction_id: str
    parent_id: str | None = None
    reactants: tuple[StoichTerm, ...]
    products: tuple[StoichTerm, ...]
    minima: tuple[str, str]  # (reactant side, product side); equal for degenerate reactions
    endpoints: tuple[str, str]  # species ids used for path calculations
    degenerate: bool = False
    source: Literal["declared", "discovery", "mode_follow", "conformer", "split", "reassigned"]
    coordinate: tuple[CoordinateTerm, ...] = ()
    torsional: bool = False
    n_h_transferred: int = 0
    low_level_ts: Geometry | None = None
    negative_evidence: tuple[str, ...] = ()
    barrier: BarrierVerdict | None = None
    saddle: SaddleClaim | None = None
    connection: ConnectionClaim | None = None
    outcome: CaseOutcome | None = None
    reasons: tuple[str, ...] = ()
    log: str | None = None  # cases/<reaction_id>/log.jsonl (write-only)


class SpeciesThermo(BaseModel):
    model_config = _FROZEN
    kind: Literal["species_thermo"] = "species_thermo"
    subject: str  # minimum_id or SaddleClaim.freq_calc
    freq_calc: str
    energy_calc: str | None = None  # sp used for the composite energy
    T_K: float
    G_hartree: float | None  # None means thermo_unavailable
    H_hartree: float | None
    zpe_hartree: float | None
    settings_sha: str
    population: float | None = None  # Boltzmann weight within the composition
    notes: tuple[str, ...] = ()


class ReactionThermo(BaseModel):
    model_config = _FROZEN
    kind: Literal["reaction_thermo"] = "reaction_thermo"
    reaction_id: str
    T_K: float
    standard_state: Literal["1atm", "1bar", "1M"]
    dE_act_kcal: float | None
    dE_rxn_kcal: float | None
    dzpe_act_kcal: float | None  # dE0(act) - dE(act), used by rankable
    dG_act_kcal: float | None
    dG_rxn_kcal: float | None
    dG_assoc_kcal: float | None = None
    dG_act_vs_separated_kcal: float | None = None
    band_kcal: tuple[float, float] | None = None
    blockers: tuple[str, ...] = ()


class RankRow(BaseModel):
    model_config = _FROZEN
    reaction_id: str
    outcome: CaseOutcome
    tier: Literal["screening", "minima", "saddle", "connected"]
    rankable: bool
    rank: int | None
    dG_act_kcal: float | None
    band_kcal: tuple[float, float] | None
    blockers: tuple[str, ...]


class ReportRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["report"] = "report"
    rows: tuple[RankRow, ...]
    tables: dict[str, FileRef]  # ranking.csv, coverage.csv, method_panel.csv, report.html


Payload = Annotated[
    Evidence | SpeciesRecord | MinimumRecord | DiscoveryRecord | ReactionRecord
    | SpeciesThermo | ReactionThermo | ReportRecord,
    Field(discriminator="kind"),
]
