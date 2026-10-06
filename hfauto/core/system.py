"""System file models (design §5.2); xyz paths are relative to the system file."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, PositiveInt, model_validator

from hfauto.core.records import CoordinateTerm

_ATOMS_PER_TERM = {"distance": 2, "angle": 3, "dihedral": 4}


class SpeciesInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    xyz: Path | None = None
    smiles: str | None = None
    charge: int = 0
    multiplicity: PositiveInt  # 2S + 1, always declared: no structure implies its spin state
    role: Literal["monomer", "endpoint"] = "monomer"
    # (degeneracy, cm⁻¹) of the ground term's levels for a linear radical with Λ > 0 (OH, SH,
    # NO, ClO): its electronic term and spin-orbit lowering; others get the spin multiplet (an
    # isolated atom its NIST term)
    electronic_levels: tuple[tuple[PositiveInt, float], ...] | None = None

    @model_validator(mode="after")
    def _one_source(self) -> SpeciesInput:
        if (self.xyz is None) == (self.smiles is None):
            raise ValueError(f"species {self.id}: give exactly one of xyz or smiles")
        if self.electronic_levels and sum(g for g, _ in self.electronic_levels) % self.multiplicity:
            raise ValueError(f"species {self.id}: the levels' degeneracies must add up to a "
                             "multiple of the multiplicity")
        return self


class CompositionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    components: dict[str, PositiveInt]  # species id -> count; the charge is their sum
    multiplicity: int | None = None  # None: the only one spin coupling allows


class ReactionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    reactant: str  # species id with role=endpoint; atom order must match the product
    product: str
    coordinate: list[CoordinateTerm] = []
    torsional: bool | None = None  # None means True when no bond changes

    @model_validator(mode="after")
    def _check_coordinate(self) -> ReactionInput:
        for term in self.coordinate:
            atoms, n = term.atoms, _ATOMS_PER_TERM[term.kind]
            if len(atoms) != n or min(atoms) < 0 or len(set(atoms)) != n:
                raise ValueError(f"reaction {self.id}: a {term.kind} needs {n} distinct atom "
                                 f"indices >= 0, got {list(atoms)}")
        return self


def _duplicate_ids(items: Sequence[SpeciesInput | CompositionInput | ReactionInput]) -> list[str]:
    return sorted(k for k, n in Counter(item.id for item in items).items() if n > 1)


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system_id: str
    species: list[SpeciesInput]
    compositions: list[CompositionInput] = []
    reactions: list[ReactionInput] = []

    @model_validator(mode="after")
    def _check_references(self) -> SystemConfig:
        named = [*self.species, *self.compositions]  # conformers keys artifacts by these ids
        for what, items in (("species/composition", named), ("reaction", self.reactions)):
            duplicated = _duplicate_ids(items)
            if duplicated:
                raise ValueError(f"duplicate {what} ids: {duplicated}")
        species = {s.id: s for s in self.species}
        for composition in self.compositions:
            unknown = sorted(set(composition.components) - set(species))
            if unknown:
                raise ValueError(f"composition {composition.id}: unknown species {unknown}")
        for reaction in self.reactions:
            for species_id in (reaction.reactant, reaction.product):
                endpoint = species.get(species_id)
                if endpoint is None or endpoint.role != "endpoint":
                    raise ValueError(
                        f"reaction {reaction.id}: {species_id!r} is not an endpoint species"
                    )
                if endpoint.xyz is None:
                    raise ValueError(f"reaction {reaction.id}: endpoint {species_id!r} must be "
                                     "given as xyz; SMILES cannot fix the atom mapping or the "
                                     "conformer")
        return self


def load_system(path: Path) -> SystemConfig:
    """Read a system YAML file and resolve xyz paths against its directory."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    config = SystemConfig.model_validate(raw)
    base = path.resolve().parent
    species = [
        s.model_copy(update={"xyz": (base / s.xyz).resolve()}) if s.xyz is not None else s
        for s in config.species
    ]
    return config.model_copy(update={"species": species})
