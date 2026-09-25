"""System file models (design §5.2); xyz paths are relative to the system file."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, model_validator

from hfauto.core.records import CoordinateTerm


class SpeciesInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    xyz: Path | None = None
    smiles: str | None = None
    charge: int = 0
    multiplicity: int = 1
    role: Literal["monomer", "endpoint"] = "monomer"


class CompositionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    components: dict[str, int]  # species id -> count; charge and multiplicity follow from them


class ReactionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    reactant: str  # species id with role=endpoint; atom order must match the product
    product: str
    coordinate: list[CoordinateTerm] = []
    torsional: bool | None = None  # None means True when no bond changes


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system_id: str
    species: list[SpeciesInput]
    compositions: list[CompositionInput] = []
    reactions: list[ReactionInput] = []

    @model_validator(mode="after")
    def _check_references(self) -> SystemConfig:
        counts = Counter(s.id for s in self.species)
        duplicated = sorted(k for k, n in counts.items() if n > 1)
        if duplicated:
            raise ValueError(f"duplicate species ids: {duplicated}")
        roles = {s.id: s.role for s in self.species}
        for composition in self.compositions:
            unknown = sorted(set(composition.components) - set(roles))
            if unknown:
                raise ValueError(f"composition {composition.id}: unknown species {unknown}")
        for reaction in self.reactions:
            for species_id in (reaction.reactant, reaction.product):
                if roles.get(species_id) != "endpoint":
                    raise ValueError(
                        f"reaction {reaction.id}: {species_id!r} is not an endpoint species"
                    )
        return self


class Conditions(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    temperatures_K: tuple[float, ...] = (298.15,)
    standard_states: tuple[Literal["1atm", "1bar", "1M"], ...] = ("1atm",)


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
