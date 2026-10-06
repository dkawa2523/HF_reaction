"""Artifact payload records (design §5.3). Payloads are a union discriminated by ``kind``."""

from __future__ import annotations

from enum import StrEnum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field

from hfauto.core.evidence import Evidence, FileRef, Geometry

_FROZEN = ConfigDict(frozen=True, extra="forbid")

StandardState = Literal["1atm", "1bar", "1M"]
ConnectionLabel = Literal["elementary", "degenerate", "reassigned", "failed"]  # gates.connection
ProfileSource = Literal["screen", "string", "scan"]  # scan: an association's relaxed scan


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
    members: tuple[str, ...] = ()  # species that fell into this basin (collapsed seeds included)
    notes: tuple[str, ...] = ()  # *_imaginary_mode, spin_contaminated, endpoint_was_saddle


class ReactionTrial(BaseModel):
    """One edit class of a state (chemistry.trials), as realised on one of its structures."""

    model_config = _FROZEN
    trial_id: str  # of the state's composition and the class only
    kind: str  # the edit: f<formed>b<broken>
    associations: tuple[tuple[int, int], ...] = ()
    dissociations: tuple[tuple[int, int], ...] = ()
    source_minimum: str | None = None  # read by nothing; kept so that older manifests parse


class DiscoveryRecord(BaseModel):
    """An edge between two states, low-level (explore) or a DFT-tier saddle's sides
    (mode_follow), or an explore attempt that added none. An edge (``product``) runs from
    ``source_species`` to ``product_species``, two species in one atom order: its ends."""

    model_config = _FROZEN
    kind: Literal["discovery"] = "discovery"
    discovery_id: str
    mechanism: Literal["nt2", "relaxation", "mode_follow"]
    trial: ReactionTrial | None = None
    # product: an edge reached from the start states; unconnected: an edge not reached (no
    # hypothesis); not_attempted: a class the budget cut (without a trial: a state left)
    outcome: Literal["product", "unconnected", "negative", "failed", "not_attempted"]
    reason: str | None = None
    source_species: str | None = None  # an attempt's: the structure it started from
    product_species: str | None = None
    ts: Geometry | None = None  # None: an edge without a TS (a relaxation)
    dE_act_kcal: float | None = None  # low level, from the source end
    dE_rxn_kcal: float | None = None
    generation: int | None = None  # explore: edges from the start states to its source end + 1
    ts_calc: str | None = None  # mode_follow: the DFT saddle's opt calculation
    # read by nothing; kept so that older manifests parse
    source_minimum: str | None = None
    ts_imag_cm1: float | None = None
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
    """Class of the latest DFT profile (profile.judge), the DFT minima energies at its ends
    (an association's scan: the sum of the separated monomers' and the adduct's)."""

    model_config = _FROZEN
    verdict: Literal["barrierless", "single", "intermediate", "unavailable"]
    source: ProfileSource
    reasons: tuple[str, ...] = ()
    points: tuple[str, ...] = ()  # the calc ids of the profile's nodes, in path order


class SaddleClaim(BaseModel):
    model_config = _FROZEN
    saddle_calc: str
    freq_calc: str
    imag_cm1: float
    energy_hartree: float
    notes: tuple[str, ...] = ()


class ConnectionClaim(BaseModel):
    """The QRC sides of a validated TS and the minima they reached."""

    model_config = _FROZEN
    side_calcs: tuple[str, str]
    minima: tuple[str, str]


class CaseOutcome(StrEnum):
    ELEMENTARY_STEP = "elementary_step"
    DEGENERATE = "degenerate_rearrangement"
    REASSIGNED = "reassigned_step"
    MULTI_STEP = "multi_step"
    BARRIERLESS = "barrierless_at_resolution"
    SAME_BASIN = "same_basin"
    OUT_OF_WINDOW = "out_of_window"
    UNRESOLVED = "unresolved_within_budget"
    BLOCKED = "blocked_upstream"


# the outcome of each connected label: a validated TS whose QRC sides reached assigned minima
CONNECTED_OUTCOMES: dict[ConnectionLabel, CaseOutcome] = {
    "elementary": CaseOutcome.ELEMENTARY_STEP,
    "degenerate": CaseOutcome.DEGENERATE,
    "reassigned": CaseOutcome.REASSIGNED,
}


class ReactionRecord(BaseModel):
    model_config = _FROZEN
    kind: Literal["reaction"] = "reaction"
    reaction_id: str  # a split child is "<parent reaction_id>_split<n>"
    reactants: tuple[StoichTerm, ...]
    products: tuple[StoichTerm, ...]
    minima: tuple[str, str]  # (reactant side, product side); equal for degenerate reactions
    endpoints: tuple[str, str]  # species ids used for path calculations
    # an association (chemistry.hypotheses): the DFT minima of the separated monomers, repeated
    # by count, its reactant side; minima[0], the complex, identifies it but ends no path
    monomers: tuple[str, ...] = ()
    degenerate: bool = False
    source: Literal["declared", "discovery", "mode_follow", "split", "reassigned"]
    coordinate: tuple[CoordinateTerm, ...] = ()
    torsional: bool = False
    # every distinct low-level TS of the hypothesis's case key (chemistry.hypotheses), in
    # priority order: the TS seeds the case tries in turn
    low_level_ts: tuple[Geometry, ...] = ()
    # a DFT stationary point at the case Level (a mode-follow saddle or a split parent's TS):
    # validated first instead of searched again
    ts_calc: str | None = None
    barrier: BarrierVerdict | None = None
    saddle: SaddleClaim | None = None
    connection: ConnectionClaim | None = None
    outcome: CaseOutcome | None = None
    reasons: tuple[str, ...] = ()
    log: str | None = None  # cases/<reaction_id>/log.jsonl (write-only)
    steps: tuple[str, ...] = ()  # a multi-step parent: its child reaction ids, in path order


class SpeciesThermo(BaseModel):
    model_config = _FROZEN
    kind: Literal["species_thermo"] = "species_thermo"
    subject: str  # minimum_id or SaddleClaim.freq_calc
    freq_calc: str
    energy_calc: str | None = None  # sp used for the composite energy
    T_K: float
    G_hartree: float | None  # None: unavailable (notes say why, e.g. energy_layer_missing)
    H_hartree: float | None
    zpe_hartree: float | None
    settings_sha: str
    # the point group the G counts (chemistry.symmetry); None in manifests written before it
    point_group: str | None = None
    sigma: int | None = None
    m: int | None = None
    notes: tuple[str, ...] = ()


class ReactionThermo(BaseModel):
    model_config = _FROZEN
    kind: Literal["reaction_thermo"] = "reaction_thermo"
    reaction_id: str
    T_K: float
    standard_state: StandardState
    dE_act_kcal: float | None
    dE_rxn_kcal: float | None
    dG_act_kcal: float | None  # TST barrier seen from the reaction's own reactant minimum
    dG_rxn_kcal: float | None
    dG_assoc_kcal: float | None = None  # the precursor complex against its separated monomers
    dG_act_vs_separated_kcal: float | None = None
    band_kcal: tuple[float, float] | None = None  # dG_eff over the qs x cutoff variants
    blockers: tuple[str, ...] = ()
    # the ranking quantity of a connected outcome (chemistry.thermo.chain_barrier); None for a
    # barrierless one (capture-limited: no TST barrier)
    dG_eff_kcal: float | None = None
    reference: Literal["complex", "separated"] | None = None  # dG_eff's zero, with R_sep only
    energy_level: str | None = None  # method/basis of the energies, e.g. m06-2x-d3/def2-tzvpd
    # the profile judged again on the energy layer, and its highest point above the reactant
    layer_verdict: Literal["barrierless", "single", "unavailable"] | None = None
    dE_act_path_kcal: float | None = None
    notes: tuple[str, ...] = ()  # submerged_barrier


class RankRow(BaseModel):
    model_config = _FROZEN
    reaction_id: str
    outcome: CaseOutcome
    tier: Literal["screening", "minima", "saddle", "connected"]
    rankable: bool
    rank: int | None
    blockers: tuple[str, ...]
    # rank_rows always sets these; the defaults keep reports of older manifests valid.
    # html.render shows the ReactionThermo of the rows' (T_K, standard_state); T_K None: the first.
    dG_act_kcal: float | None = None
    band_kcal: tuple[float, float] | None = None
    T_K: float | None = None
    standard_state: str | None = None
    dG_rxn_kcal: float | None = None
    dG_eff_kcal: float | None = None
    reference: Literal["complex", "separated"] | None = None
    dG_act_vs_separated_kcal: float | None = None
    torsional: bool = False
    energy_level: str | None = None
    notes: tuple[str, ...] = ()


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
