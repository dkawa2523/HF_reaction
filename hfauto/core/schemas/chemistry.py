"""Small, chemistry-neutral records for molecular reaction workflows.

The records describe *what* is being calculated.  Molecule-specific reaction
rules belong in chemistry modules, not in these transport models.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

Phase = Literal["gas", "condensed", "solution", "unknown"]
BondChangeKind = Literal["form", "break", "increase_order", "decrease_order", "transfer"]
CoordinateKind = Literal["distance", "angle", "dihedral"]
DiscoveryDriver = Literal[
    "nt2", "afir", "unbiased_preopt", "dft_minimum_ensemble"
]
EndpointEvidenceLevel = Literal["low_level", "dft_minimum_ensemble"]


class ComponentRecord(BaseModel):
    """One molecular or material fragment inside a chemical state."""

    component_id: str
    atom_indices: list[int]
    label: str | None = None
    role: str | None = None
    source_artifact_id: str | None = None
    element_counts: dict[str, int] = Field(default_factory=dict)
    charge: int = 0
    multiplicity: int = 1

    @field_validator("atom_indices")
    @classmethod
    def validate_atom_indices(cls, value: list[int]) -> list[int]:
        indices = [int(index) for index in value]
        if not indices:
            raise ValueError("a component must contain at least one atom")
        if min(indices) < 0 or len(indices) != len(set(indices)):
            raise ValueError("component atom indices must be unique and non-negative")
        return indices


class ChemicalStateRecord(BaseModel):
    """Composition, electronic state, and environment for one geometry."""

    state_id: str
    state_type: str
    components: list[ComponentRecord]
    atom_count: int | None = None
    charge: int = 0
    multiplicity: int = 1
    phase: Phase = "unknown"
    environment: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_component_partition(self) -> ChemicalStateRecord:
        indices = [index for component in self.components for index in component.atom_indices]
        if len(indices) != len(set(indices)):
            raise ValueError("components must not share atom indices")
        if self.atom_count is not None:
            expected = set(range(int(self.atom_count)))
            if set(indices) != expected:
                raise ValueError("components must partition every atom exactly once")
        return self


class BondChangeRecord(BaseModel):
    """An atom-mapped change used to define and validate a reaction path."""

    kind: BondChangeKind
    atoms: tuple[int, int]
    label: str | None = None
    reactant_order: float | None = None
    product_order: float | None = None
    min_distance_change_A: float = 0.15
    weight: float = 1.0

    @field_validator("atoms")
    @classmethod
    def validate_atoms(cls, value: tuple[int, int]) -> tuple[int, int]:
        first, second = (int(value[0]), int(value[1]))
        if first < 0 or second < 0 or first == second:
            raise ValueError("bond-change atoms must be distinct non-negative indices")
        return first, second

    @field_validator("min_distance_change_A")
    @classmethod
    def validate_distance_change(cls, value: float) -> float:
        if float(value) < 0.0:
            raise ValueError("min_distance_change_A must be non-negative")
        return float(value)

    @field_validator("weight")
    @classmethod
    def validate_weight(cls, value: float) -> float:
        if float(value) <= 0.0:
            raise ValueError("bond-change weight must be positive")
        return float(value)


class ReactionCoordinateTermRecord(BaseModel):
    """One internal-coordinate term; angular values are evaluated in radians."""

    kind: CoordinateKind
    atoms: tuple[int, ...]
    coefficient: float = 1.0
    label: str | None = None

    @model_validator(mode="after")
    def validate_atom_count(self) -> ReactionCoordinateTermRecord:
        expected = {"distance": 2, "angle": 3, "dihedral": 4}[self.kind]
        if len(self.atoms) != expected:
            raise ValueError(f"{self.kind} coordinate requires {expected} atom indices")
        if min(self.atoms, default=-1) < 0 or len(set(self.atoms)) != len(self.atoms):
            raise ValueError("coordinate atoms must be unique and non-negative")
        return self


class ReactionCoordinateRecord(BaseModel):
    """Weighted internal-coordinate definition for path and mode validation."""

    terms: list[ReactionCoordinateTermRecord] = Field(default_factory=list)
    min_change: float = 0.05

    @field_validator("min_change")
    @classmethod
    def validate_min_change(cls, value: float) -> float:
        if float(value) < 0.0:
            raise ValueError("reaction-coordinate min_change must be non-negative")
        return float(value)


class StoichiometricTermRecord(BaseModel):
    """One positive-coefficient species on one side of a reaction equation."""

    species_id: str
    coefficient: float = 1.0
    role: str | None = None

    @field_validator("coefficient")
    @classmethod
    def validate_coefficient(cls, value: float) -> float:
        if float(value) <= 0.0:
            raise ValueError("stoichiometric coefficient must be positive")
        return float(value)


class ReactionHypothesisRecord(BaseModel):
    """A mechanism hypothesis independent of a particular path-search engine."""

    reaction_id: str
    reactant_state_id: str
    product_state_id: str
    mechanism_family: str
    bond_changes: list[BondChangeRecord]
    reactants: list[StoichiometricTermRecord] = Field(default_factory=list)
    products: list[StoichiometricTermRecord] = Field(default_factory=list)
    reaction_coordinate: ReactionCoordinateRecord = Field(
        default_factory=ReactionCoordinateRecord
    )
    rationale: list[str] = Field(default_factory=list)
    environment: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def populate_endpoint_stoichiometry(self) -> ReactionHypothesisRecord:
        if not self.reactants:
            self.reactants = [StoichiometricTermRecord(species_id=self.reactant_state_id)]
        if not self.products:
            self.products = [StoichiometricTermRecord(species_id=self.product_state_id)]
        return self


class DrivingAtomPairRecord(BaseModel):
    """One unordered atom pair driven together or apart during discovery."""

    atoms: tuple[int, int]
    label: str | None = None

    @field_validator("atoms")
    @classmethod
    def validate_atoms(cls, value: tuple[int, int]) -> tuple[int, int]:
        first, second = sorted((int(value[0]), int(value[1])))
        if first < 0 or first == second:
            raise ValueError("driving atoms must be distinct non-negative indices")
        return first, second


class ReactionTrialRecord(BaseModel):
    """A bounded single-ended search request with no assumed product geometry."""

    trial_id: str
    source_species_id: str
    associations: list[DrivingAtomPairRecord] = Field(default_factory=list)
    dissociations: list[DrivingAtomPairRecord] = Field(default_factory=list)
    reaction_coordinate: ReactionCoordinateRecord = Field(
        default_factory=ReactionCoordinateRecord
    )
    driver_order: list[DiscoveryDriver] = Field(default_factory=lambda: ["nt2", "afir"])
    mechanism_hint: str = "unspecified"
    priority: float = 0.0
    rationale: list[str] = Field(default_factory=list)
    component_ids: list[str] = Field(default_factory=list)
    charge: int = 0
    multiplicity: int = 1
    max_attempts: int = 2
    settings: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_trial(self) -> ReactionTrialRecord:
        if (
            not self.associations
            and not self.dissociations
            and not self.reaction_coordinate.terms
        ):
            raise ValueError(
                "a reaction trial requires a driving pair or reaction-coordinate term"
            )
        associated = {pair.atoms for pair in self.associations}
        dissociated = {pair.atoms for pair in self.dissociations}
        if associated & dissociated:
            raise ValueError("one atom pair cannot be both associated and dissociated")
        if not self.driver_order or len(self.driver_order) != len(set(self.driver_order)):
            raise ValueError("driver_order must contain unique discovery drivers")
        if self.multiplicity < 1 or self.max_attempts < 1:
            raise ValueError("multiplicity and max_attempts must be positive")
        return self

    def bond_changes(self) -> list[BondChangeRecord]:
        """Translate the driving definition into the shared path contract."""

        return [
            *(
                BondChangeRecord(kind="form", atoms=pair.atoms, label=pair.label)
                for pair in self.associations
            ),
            *(
                BondChangeRecord(kind="break", atoms=pair.atoms, label=pair.label)
                for pair in self.dissociations
            ),
        ]


class ReactionCandidateRecord(BaseModel):
    """A proposed endpoint pair with an explicit level of structural evidence."""

    candidate_id: str
    trial_id: str
    reactant_species_id: str
    product_species_id: str
    driver: DiscoveryDriver
    bond_changes: list[BondChangeRecord]
    composition_preserved: bool
    charge: int
    multiplicity: int
    endpoint_evidence_level: EndpointEvidenceLevel = "low_level"
    low_level_endpoint_distinct: bool | None = None
    low_level_ts_validated: bool = False
    low_level_irc_connected: bool = False
    evidence_artifact_id: str
