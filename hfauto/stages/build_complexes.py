"""Config-driven assembly of arbitrary molecular encounter complexes."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.conformer.crest import CRESTConformerBackend
from hfauto.chemistry.complexes import append_fragment, element_counts, write_complex
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.ids import species_id as make_species_id
from hfauto.core.io import (
    ensure_dir,
    read_jsonl,
    read_manifest,
    write_jsonl,
    write_manifest,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import (
    ChemicalStateRecord,
    ComponentRecord,
)
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext

SCHEMA_VERSION = "hfauto.build_complexes.v1"
SOURCE_ARTIFACT_TYPES = {
    "molecule",
    "conformer",
    "species",
    "species_preopt",
    "species_optimized",
}

_ACCEPTOR_PRIORITY = ("N", "O", "S", "P", "F", "Cl", "Br", "I")


def _next_attempt_root(base: Path) -> Path:
    """Return a fresh calculation root without overwriting raw CREST evidence."""

    ensure_dir(base)
    attempt = 0
    while (candidate := base / f"attempt_{attempt:02d}").exists():
        attempt += 1
    return ensure_dir(candidate)


def _successful_nci_counts(
    manifest: Manifest, valid_request_fingerprints: set[str]
) -> dict[str, int]:
    """Count NCI species proven to match completed build requests."""

    counts: dict[str, int] = {}
    for species in manifest.latest_artifacts("species"):
        composition_id = species.data.get("composition_id")
        if (
            species.status.status != "success"
            or composition_id is None
            or not species.qc.get("crest_nci_searched")
            or species.provenance.get("build_complex_request_fingerprint")
            not in valid_request_fingerprints
        ):
            continue
        key = str(composition_id)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _definition_request_fingerprint(
    manifest: Manifest,
    definition: dict[str, Any],
    config: dict[str, Any],
) -> str:
    """Identify one build request including its exact geometry inputs."""

    source_hashes: list[dict[str, str]] = []
    if definition.get("xyz_path"):
        path = Path(str(definition["xyz_path"]))
        source_hashes.append({"path": str(path), "sha256": sha256_file(path)})
    else:
        for fragment in list(definition.get("fragments") or []):
            source = _resolve_source(manifest, fragment["source"])
            geometry = species_xyz_path(source)
            source_hashes.append(
                {
                    "artifact_id": source.artifact_id,
                    "sha256": sha256_file(geometry),
                }
            )
    nci_settings = (
        {
            **dict(config.get("crest_nci_settings") or {}),
            **dict(definition.get("nci_settings") or {}),
        }
        if definition.get("crest_nci")
        else {}
    )
    return fingerprint_dict(
        {
            "schema_version": SCHEMA_VERSION,
            "definition": definition,
            "source_geometries": source_hashes,
            "nci_settings": nci_settings,
        }
    )


def _source_artifacts(manifest: Manifest) -> list[Artifact]:
    return [
        artifact
        for artifact in manifest.latest_artifacts()
        if artifact.status.status == "success" and artifact.artifact_type in SOURCE_ARTIFACT_TYPES
    ]


def _require_usable_source(artifact: Artifact) -> Artifact:
    if artifact.status.status != "success":
        raise ValueError(f"source artifact is not successful: {artifact.artifact_id}")
    if artifact.artifact_type not in SOURCE_ARTIFACT_TYPES:
        raise ValueError(
            f"unsupported source artifact type {artifact.artifact_type!r}: {artifact.artifact_id}"
        )
    return artifact


def _resolve_source(manifest: Manifest, source: str | dict[str, Any]) -> Artifact:
    if isinstance(source, str):
        artifact = manifest.find(source)
        if artifact is None:
            raise ValueError(f"source artifact not found: {source}")
        return _require_usable_source(artifact)
    if not isinstance(source, dict):
        raise TypeError("fragment source must be an artifact id or selector mapping")
    artifact_id = source.get("artifact_id")
    if artifact_id:
        artifact = manifest.find(str(artifact_id))
        if artifact is None:
            raise ValueError(f"source artifact not found: {artifact_id}")
        return _require_usable_source(artifact)
    criteria = {key: value for key, value in source.items() if key != "artifact_type"}
    matches = [
        artifact
        for artifact in _source_artifacts(manifest)
        if (
            source.get("artifact_type") is None or artifact.artifact_type == source["artifact_type"]
        )
        and all(artifact.data.get(key) == value for key, value in criteria.items())
    ]
    if len(matches) != 1:
        raise ValueError(
            f"source selector must resolve exactly one artifact; selector={source!r}, matches={len(matches)}"
        )
    return matches[0]


def _fragment_xyz(source: Artifact) -> XYZ:
    return read_xyz(species_xyz_path(source))


def _first_index(symbols: list[str], priority: tuple[str, ...]) -> int:
    for symbol in priority:
        if symbol in symbols:
            return symbols.index(symbol)
    return next(
        (index for index, symbol in enumerate(symbols) if symbol != "H"),
        0,
    )


def _host_anchor(source: Artifact) -> int:
    return _first_index(_fragment_xyz(source).symbols, _ACCEPTOR_PRIORITY)


def _guest_anchor(source: Artifact) -> int:
    symbols = _fragment_xyz(source).symbols
    if "H" in symbols and any(symbol in _ACCEPTOR_PRIORITY for symbol in symbols):
        return symbols.index("H")
    return _first_index(symbols, _ACCEPTOR_PRIORITY)


def _guest_orientation_atom(source: Artifact, anchor: int) -> int | None:
    """Return the heavy atom behind a donor H, when one is present."""

    symbols = _fragment_xyz(source).symbols
    if symbols[anchor] != "H":
        return None
    return next(
        (
            index
            for symbol in _ACCEPTOR_PRIORITY
            for index, observed in enumerate(symbols)
            if observed == symbol and index != anchor
        ),
        None,
    )


def _composition_definitions(
    manifest: Manifest,
    compositions: list[dict[str, Any]],
    *,
    max_components: int,
) -> list[dict[str, Any]]:
    """Expand selector/count compositions into a bounded set of geometry seeds."""

    definitions: list[dict[str, Any]] = []
    for composition in compositions:
        composition_id = str(
            composition.get("composition_id")
            or composition.get("system_id")
            or ""
        )
        if not composition_id:
            raise ValueError("composition_id must be non-empty")
        expanded: list[tuple[Artifact, dict[str, Any], int]] = []
        for component in list(composition.get("components") or []):
            if "selector" not in component:
                raise ValueError(
                    f"composition {composition_id!r} component requires selector"
                )
            count = int(component.get("count", 1))
            if count < 1:
                raise ValueError("component count must be positive")
            source = _resolve_source(manifest, component["selector"])
            expanded.extend(
                (source, component, copy_index) for copy_index in range(count)
            )
        if not 1 <= len(expanded) <= max_components:
            raise ValueError(
                f"composition {composition_id!r} expands to {len(expanded)} "
                f"components; allowed range is 1..{max_components}"
            )
        max_seeds = max(1, int(composition.get("max_seeds", 4)))
        distances = [
            float(value)
            for value in (composition.get("seed_distances_A") or [1.7, 2.1])
        ]
        if not distances or any(value <= 0.0 for value in distances):
            raise ValueError("seed_distances_A must contain positive distances")
        default_topologies = ["chain", "star"] if len(expanded) > 2 else ["star"]
        topologies = list(
            composition.get("placement_topologies") or default_topologies
        )
        if not topologies or any(item not in {"chain", "star"} for item in topologies):
            raise ValueError("placement_topologies must contain 'chain' or 'star'")
        for seed_index in range(max_seeds):
            topology = str(topologies[seed_index % len(topologies)])
            fragments: list[dict[str, Any]] = []
            for component_index, (source, component, copy_index) in enumerate(expanded):
                component_id = f"{composition_id}_component_{component_index:02d}"
                fragment = {
                    "source": source.artifact_id,
                    "component_id": component_id,
                    "label": component.get("label"),
                    "role": component.get("role"),
                }
                if component_index:
                    host_index = 0 if topology == "star" else component_index - 1
                    host_source = expanded[host_index][0]
                    guest_anchor = _guest_anchor(source)
                    placement = {
                        "host_component_id": (
                            f"{composition_id}_component_{host_index:02d}"
                        ),
                        "host_atom_index": _host_anchor(host_source),
                        "guest_atom_index": guest_anchor,
                        "distance_A": distances[
                            (seed_index + component_index - 1) % len(distances)
                        ],
                        "twist_degrees": float((120 * seed_index) % 360),
                        "min_interfragment_distance_A": float(
                            composition.get("min_interfragment_distance_A", 0.55)
                        ),
                        "placement_topology": topology,
                    }
                    orientation_atom = _guest_orientation_atom(
                        source, guest_anchor
                    )
                    if orientation_atom is not None:
                        placement.update(
                            {
                                "guest_orientation_atom_index": orientation_atom,
                                "orientation": "away_from_host",
                            }
                        )
                    fragment["placement"] = placement
                fragment["copy_index"] = copy_index
                fragments.append(fragment)
            definition = {
                "system_id": f"{composition_id}_seed{seed_index:02d}",
                "composition_id": composition_id,
                "state": composition.get("state", "encounter_complex"),
                "phase": composition.get("phase", "gas"),
                "environment": dict(composition.get("environment") or {}),
                "fragments": fragments,
                "crest_nci": bool(
                    composition.get("crest_nci", len(expanded) > 1)
                ),
                "nci_settings": dict(composition.get("nci_settings") or {}),
                "max_nci_complexes": int(
                    composition.get("max_nci_complexes", 4)
                ),
                "placement_topology": topology,
            }
            for key in ("charge", "multiplicity"):
                if composition.get(key) is not None:
                    definition[key] = composition[key]
            definitions.append(definition)
    return definitions


def _electronic_state(
    definition: dict[str, Any], sources: list[Artifact], symbols: list[str]
) -> tuple[int, int]:
    charge = int(
        definition.get(
            "charge",
            sum(
                int(item.data.get("charge", item.data.get("formal_charge", 0)) or 0)
                for item in sources
            ),
        )
    )
    multiplicities = [int(item.data.get("multiplicity", 1) or 1) for item in sources]
    if "multiplicity" not in definition and any(value != 1 for value in multiplicities):
        raise ValueError("system multiplicity is required when an input fragment is open-shell")
    state = resolve_electronic_state(
        {
            "charge": charge,
            "multiplicity": int(definition.get("multiplicity", 1)),
        },
        symbols=symbols,
    )
    return int(state["charge"]), int(state["multiplicity"])


def _build_system(
    manifest: Manifest,
    definition: dict[str, Any],
    out_dir: Path,
) -> Artifact:
    system_id = str(definition["system_id"])
    source_xyz = definition.get("xyz_path")
    fragment_definitions = list(definition.get("fragments", []) or [])
    if source_xyz and fragment_definitions:
        raise ValueError(f"system {system_id!r} cannot define both xyz_path and fragments")
    if source_xyz:
        return _build_xyz_system(definition, out_dir)
    if not fragment_definitions:
        raise ValueError(f"system {system_id!r} requires xyz_path or fragments")

    sources = [_resolve_source(manifest, fragment["source"]) for fragment in fragment_definitions]
    first_xyz = _fragment_xyz(sources[0])
    assembled = XYZ(list(first_xyz.symbols), first_xyz.coords.copy(), first_xyz.comment)
    components: list[ComponentRecord] = []
    placements: list[dict[str, Any]] = []

    first_definition = fragment_definitions[0]
    first_component_id = str(first_definition.get("component_id") or sources[0].artifact_id)
    components.append(
        ComponentRecord(
            component_id=first_component_id,
            label=first_definition.get("label"),
            role=first_definition.get("role"),
            source_artifact_id=sources[0].artifact_id,
            atom_indices=list(range(len(first_xyz.symbols))),
            element_counts=element_counts(first_xyz.symbols),
            charge=int(sources[0].data.get("charge", sources[0].data.get("formal_charge", 0)) or 0),
            multiplicity=int(sources[0].data.get("multiplicity", 1) or 1),
        )
    )

    for fragment_definition, source in zip(fragment_definitions[1:], sources[1:]):
        guest = _fragment_xyz(source)
        placement = dict(fragment_definition.get("placement", {}) or {})
        host_component_id = str(placement.get("host_component_id") or components[0].component_id)
        host_component = next(
            (component for component in components if component.component_id == host_component_id),
            None,
        )
        if host_component is None:
            raise ValueError(
                f"unknown host component {host_component_id!r} in system {system_id!r}"
            )
        host_local_index = int(placement["host_atom_index"])
        if not 0 <= host_local_index < len(host_component.atom_indices):
            raise IndexError("host_atom_index is outside the selected host component")
        offset = len(assembled.symbols)
        minimum_allowed = float(
            placement.get("min_interfragment_distance_A", 0.60)
        )
        assembled, placement_qc = append_fragment(
            assembled,
            guest,
            host_anchor_index=host_component.atom_indices[host_local_index],
            guest_anchor_index=int(placement["guest_atom_index"]),
            distance_A=float(placement["distance_A"]),
            direction=placement.get("direction"),
            guest_orientation_atom_index=placement.get("guest_orientation_atom_index"),
            orientation=str(placement.get("orientation", "toward_host")),
            twist_degrees=float(placement.get("twist_degrees", 0.0)),
            minimum_interfragment_distance_A=minimum_allowed,
        )
        placement_qc.update(
            {
                "component_id": str(fragment_definition.get("component_id") or source.artifact_id),
                "host_component_id": host_component_id,
                "minimum_allowed_interfragment_distance_A": minimum_allowed,
            }
        )
        if placement_qc["minimum_interfragment_distance_A"] < minimum_allowed:
            raise ValueError(
                f"fragment collision in {system_id}: minimum distance "
                f"{placement_qc['minimum_interfragment_distance_A']:.3f} A < {minimum_allowed:.3f} A"
            )
        placements.append(placement_qc)
        component_id = str(fragment_definition.get("component_id") or source.artifact_id)
        if component_id in {component.component_id for component in components}:
            raise ValueError(f"duplicate component_id {component_id!r} in system {system_id!r}")
        components.append(
            ComponentRecord(
                component_id=component_id,
                label=fragment_definition.get("label"),
                role=fragment_definition.get("role"),
                source_artifact_id=source.artifact_id,
                atom_indices=list(range(offset, offset + len(guest.symbols))),
                element_counts=element_counts(guest.symbols),
                charge=int(source.data.get("charge", source.data.get("formal_charge", 0)) or 0),
                multiplicity=int(source.data.get("multiplicity", 1) or 1),
            )
        )

    charge, multiplicity = _electronic_state(definition, sources, assembled.symbols)
    state_type = str(definition.get("state", "encounter_complex"))
    species_id = make_species_id(system_id)
    xyz_path = write_complex(
        assembled, out_dir / "systems" / f"{system_id}.xyz", system_id=system_id
    )
    chemical_state = ChemicalStateRecord(
        state_id=species_id,
        state_type=state_type,
        components=components,
        atom_count=len(assembled.symbols),
        charge=charge,
        multiplicity=multiplicity,
        phase=str(definition.get("phase", "unknown")),
        environment=dict(definition.get("environment", {}) or {}),
    )
    data = {
        "species_id": species_id,
        "system_id": system_id,
        "state": state_type,
        "charge": charge,
        "multiplicity": multiplicity,
        "xyz_path": str(xyz_path),
        "atom_order_key": str(
            definition.get("atom_order_key")
            or "components:"
            + "|".join(
                f"{component.source_artifact_id}:{len(component.atom_indices)}"
                for component in components
            )
        ),
        "components": [component.model_dump() for component in components],
        "chemical_state": chemical_state.model_dump(),
        "element_counts": element_counts(assembled.symbols),
        "environment": chemical_state.environment,
        "schema_version": SCHEMA_VERSION,
    }
    if definition.get("composition_id") is not None:
        data["composition_id"] = str(definition["composition_id"])
    # A multicomponent state is not a revision of any one molecular input.
    # Retaining one fragment's mol_id makes downstream conformer IDs collide
    # with that fragment's own conformers.
    if len(sources) == 1 and sources[0].data.get("mol_id"):
        data["mol_id"] = str(sources[0].data["mol_id"])
    return Artifact(
        artifact_id=species_id,
        artifact_type="species",
        parents=list(dict.fromkeys(source.artifact_id for source in sources)),
        paths={"xyz": str(xyz_path)},
        data=data,
        qc={
            "assembled": True,
            "atom_count": len(assembled.symbols),
            "component_count": len(components),
            "placements": placements,
        },
    )


def _nci_species(seed: Artifact, conformer: Artifact, rank: int) -> Artifact:
    """Promote one CREST NCI geometry to a normal species workflow input."""

    identifier = make_species_id(seed.data["system_id"], "nci", rank)
    xyz_path = str(conformer.paths["xyz"])
    data = {
        **seed.data,
        "species_id": identifier,
        "source_species_id": identifier,
        "source_complex_seed_id": seed.data["species_id"],
        "system_id": identifier.removeprefix("spc_"),
        "state": "encounter_complex",
        "xyz_path": xyz_path,
        "crest_nci_conformer_id": conformer.artifact_id,
        "relative_energy_kcal_mol": conformer.data.get(
            "relative_energy_kcal_mol"
        ),
    }
    return Artifact(
        artifact_id=identifier,
        artifact_type="species",
        parents=[seed.artifact_id, conformer.artifact_id],
        paths={"xyz": xyz_path},
        data=data,
        method={"stage": "build-complexes", **(conformer.method or {})},
        qc={
            "assembled": True,
            "crest_nci_searched": True,
            "nci_search_fallback": False,
            "component_count": len(data.get("components") or []),
            "relative_energy_kcal_mol": conformer.data.get(
                "relative_energy_kcal_mol"
            ),
        },
    )


def _build_xyz_system(definition: dict[str, Any], out_dir: Path) -> Artifact:
    """Import one reviewed, atom-mapped XYZ as a configured chemical state."""

    system_id = str(definition["system_id"])
    source = Path(str(definition["xyz_path"]))
    structure = read_xyz(source)
    state = resolve_electronic_state(
        {
            "charge": int(definition.get("charge", 0)),
            "multiplicity": int(definition.get("multiplicity", 1)),
        },
        symbols=structure.symbols,
    )
    component = ComponentRecord(
        component_id=str(definition.get("component_id", system_id)),
        label=definition.get("label"),
        role=definition.get("role", "molecule"),
        atom_indices=list(range(len(structure.symbols))),
        element_counts=element_counts(structure.symbols),
        charge=int(state["charge"]),
        multiplicity=int(state["multiplicity"]),
    )
    state_type = str(definition.get("state", "candidate"))
    identifier = make_species_id(system_id)
    xyz_path = write_complex(
        structure, out_dir / "systems" / f"{system_id}.xyz", system_id=system_id
    )
    chemical_state = ChemicalStateRecord(
        state_id=identifier,
        state_type=state_type,
        components=[component],
        atom_count=len(structure.symbols),
        charge=int(state["charge"]),
        multiplicity=int(state["multiplicity"]),
        phase=str(definition.get("phase", "gas")),
        environment=dict(definition.get("environment", {}) or {}),
    )
    data = {
        "species_id": identifier,
        "system_id": system_id,
        "state": state_type,
        "charge": int(state["charge"]),
        "multiplicity": int(state["multiplicity"]),
        "xyz_path": str(xyz_path),
        "atom_order_key": str(definition.get("atom_order_key") or f"xyz:{system_id}"),
        "components": [component.model_dump()],
        "chemical_state": chemical_state.model_dump(),
        "element_counts": element_counts(structure.symbols),
        "environment": chemical_state.environment,
        "source_xyz_path": str(source),
        "source_xyz_sha256": sha256_file(source),
        "schema_version": SCHEMA_VERSION,
    }
    return Artifact(
        artifact_id=identifier,
        artifact_type="species",
        paths={"xyz": str(xyz_path), "source_xyz": str(source)},
        data=data,
        qc={"assembled": True, "atom_count": len(structure.symbols), "component_count": 1},
        provenance={"created_by": "BuildComplexesStage", "source_xyz_sha256": sha256_file(source)},
    )


class BuildComplexesStage(Stage):
    """Build configured chemical states; reaction definitions are a later stage."""

    name = "build-complexes"
    accepts_empty_manifest = True

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        out_dir = ensure_dir(context.out_dir)
        production_mode = (
            str(context.global_config.get("mode", "")).lower() == "production"
        )
        checkpoint_each_definition = bool(
            config.get("checkpoint_each_definition", production_mode)
        )
        resume_completed_definitions = bool(
            config.get("resume_completed_definitions", production_mode)
        )
        checkpoint_path = out_dir / "manifest.json"
        checkpoint = (
            read_manifest(checkpoint_path)
            if resume_completed_definitions and checkpoint_path.is_file()
            else None
        )
        source = manifest or Manifest.new(
            run_id=context.run_id, stage=f"{self.name}-input"
        )
        working = (
            Manifest.merge(
                [source, checkpoint],
                run_id=context.run_id,
                stage=f"{self.name}-resume",
                metadata={**source.metadata, **checkpoint.metadata},
            )
            if checkpoint is not None
            else source
        )
        out = working.carry_forward(self.name)
        records_path = out_dir / "chemical_state_records.jsonl"
        records: list[dict[str, Any]] = (
            read_jsonl(records_path) if checkpoint is not None else []
        )
        system_definitions = list(config.get("systems", []) or [])
        compositions = list(config.get("compositions", []) or [])
        if not system_definitions and not compositions:
            raise ValueError(
                "build-complexes requires systems or compositions"
            )

        completed_ids = {
            str(value)
            for value in working.metadata.get(
                "build_complex_completed_definition_ids", []
            )
        }
        completed_requests = {
            str(key): str(value)
            for key, value in dict(
                working.metadata.get("build_complex_completed_requests") or {}
            ).items()
        }
        nci_counts = _successful_nci_counts(
            working, set(completed_requests.values())
        )
        failed_ids = {
            str(value)
            for value in working.metadata.get(
                "build_complex_failed_definition_ids", []
            )
        }
        resumed_count = 0
        crest = CRESTConformerBackend()

        def checkpoint_now() -> None:
            out.metadata.update(
                {
                    "build_complex_completed_definition_ids": sorted(
                        completed_ids
                    ),
                    "build_complex_completed_requests": dict(
                        sorted(completed_requests.items())
                    ),
                    "build_complex_failed_definition_ids": sorted(failed_ids),
                    "build_complex_resumed_definition_count": resumed_count,
                    "build_complex_nci_counts": dict(sorted(nci_counts.items())),
                }
            )
            write_jsonl(records, records_path)
            write_manifest(out, out_dir)

        def process_definition(
            definition: dict[str, Any], request_fingerprint: str
        ) -> bool:
            system_id = str(definition.get("system_id") or "unnamed")
            try:
                artifact = _build_system(out, definition, out_dir)
            except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
                failure = Artifact.failure(
                    f"build_complex_failed_{system_id}",
                    "chemical_state",
                    str(exc),
                    category="complex_definition_invalid",
                    recoverable=True,
                    definition=definition,
                )
                failure.provenance["build_complex_request_fingerprint"] = (
                    request_fingerprint
                )
                out.add_artifact(failure)
                records.append(failure.model_dump())
                return False
            out.add_artifact(artifact)
            artifact.provenance["build_complex_request_fingerprint"] = (
                request_fingerprint
            )
            records.append(artifact.model_dump())
            if not definition.get("crest_nci"):
                return True
            composition_id = str(
                definition.get("composition_id") or definition["system_id"]
            )
            maximum = max(1, int(definition.get("max_nci_complexes", 4)))
            remaining = maximum - nci_counts.get(composition_id, 0)
            if remaining <= 0:
                return True
            nci_config = {
                **dict(config.get("crest_nci_settings") or {}),
                **dict(definition.get("nci_settings") or {}),
                "nci": True,
                "fallback_to_rdkit": False,
                "selected_max": remaining,
            }
            conformers = crest.generate(
                artifact,
                nci_config,
                str(
                    _next_attempt_root(
                        out_dir / "crest_nci" / composition_id / system_id
                    )
                ),
            )
            promoted = 0
            for conformer in conformers:
                existing = out.find(conformer.artifact_id)
                if (
                    existing is not None
                    and conformer.status.status == "success"
                    and set(existing.parents) != set(conformer.parents)
                ):
                    collision = Artifact.failure(
                        f"build_complex_identity_collision_{system_id}",
                        "conformer",
                        "CREST output artifact id already belongs to a different lineage",
                        category="artifact_identity_collision",
                        parents=[artifact.artifact_id, existing.artifact_id],
                        recoverable=False,
                        colliding_artifact_id=conformer.artifact_id,
                        existing_parents=existing.parents,
                        proposed_parents=conformer.parents,
                    )
                    out.add_artifact(collision)
                    records.append(collision.model_dump())
                    continue
                out.add_artifact(conformer)
                records.append(conformer.model_dump())
                if conformer.status.status != "success":
                    continue
                rank = nci_counts.get(composition_id, 0)
                if rank >= maximum:
                    break
                species = _nci_species(artifact, conformer, rank)
                species.provenance["build_complex_request_fingerprint"] = (
                    request_fingerprint
                )
                out.add_artifact(species)
                records.append(species.model_dump())
                nci_counts[composition_id] = rank + 1
                promoted += 1
            return promoted > 0

        def run_definition(definition: dict[str, Any]) -> None:
            nonlocal resumed_count
            definition_id = str(definition["system_id"])
            request_fingerprint = _definition_request_fingerprint(
                out, definition, config
            )
            previous_fingerprint = completed_requests.get(definition_id)
            if (
                resume_completed_definitions
                and previous_fingerprint == request_fingerprint
            ):
                resumed_count += 1
                return
            if (
                previous_fingerprint is not None
                and previous_fingerprint != request_fingerprint
                and definition.get("crest_nci")
            ):
                nci_counts[str(definition.get("composition_id") or definition_id)] = 0
            completed = process_definition(definition, request_fingerprint)
            if completed:
                completed_ids.add(definition_id)
                completed_requests[definition_id] = request_fingerprint
                failed_ids.discard(definition_id)
            else:
                completed_ids.discard(definition_id)
                completed_requests.pop(definition_id, None)
                failed_ids.add(definition_id)
            if checkpoint_each_definition:
                checkpoint_now()

        configured_ids = [
            str(item.get("system_id") or "") for item in system_definitions
        ]
        if "" in configured_ids or len(configured_ids) != len(set(configured_ids)):
            raise ValueError("build-complexes system_id values must be non-empty and unique")
        for definition in system_definitions:
            run_definition(definition)

        # Compositions may select systems created immediately above. Keeping
        # this order in one stage removes duplicate build-complexes stages and
        # gives partial pipeline restarts one unambiguous boundary.
        composition_definitions = _composition_definitions(
            out,
            compositions,
            max_components=int(config.get("max_components", 3)),
        )
        generated_ids = [
            str(item.get("system_id") or "") for item in composition_definitions
        ]
        all_ids = [*configured_ids, *generated_ids]
        if "" in generated_ids or len(all_ids) != len(set(all_ids)):
            raise ValueError("expanded build-complexes system_id values must be unique")
        for definition in composition_definitions:
            run_definition(definition)

        checkpoint_now()
        return out
