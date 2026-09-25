"""Molecule-level read-only dossier generation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hfauto_viz.core.html import dataframe_to_html, esc, write_html
from hfauto_viz.data.loaders import RunData
from hfauto_viz.renderers.mol3d_3dmol import html_for_structure
from hfauto_viz.renderers.rdkit_2d import write_molecule_svg
from hfauto_viz.structure.xyz import read_xyz_frames


def _first_geometry_for_molecule(run: RunData, mol_id: str):
    for species in run.species_index().values():
        if str(species.data.get("mol_id") or "") != str(mol_id):
            continue
        path = run.artifact_xyz_path(species)
        if path is None:
            continue
        frames = read_xyz_frames(path)
        if frames:
            return frames[0]
    return None


def _reaction_links(
    run: RunData, mol_id: str, *, max_reactions: int
) -> str:
    table = run.read_table("reaction_results")
    if table.empty or "mol_id" not in table.columns:
        return ""
    rows = table[table["mol_id"].astype(str) == str(mol_id)].head(
        max_reactions
    )
    items: list[str] = []
    for _, row in rows.iterrows():
        reaction_id = str(row.get("reaction_id"))
        barrier = row.get(
            "delta_G_activation_standard_kcal_mol",
            row.get("delta_G_act_kcal_mol", ""),
        )
        items.append(
            f"<li><a href='../../reactions/{esc(reaction_id)}/dossier.html'>"
            f"{esc(reaction_id)}</a> activation Gibbs energy={esc(barrier)}</li>"
        )
    if not items:
        return ""
    return (
        "<div class='card'><h2>Linked reactions</h2><ul>"
        + "".join(items)
        + "</ul></div>"
    )


def render_molecule_dossier(
    run: RunData,
    mol_id: str,
    out_dir: str | Path,
    render_reactions: bool = True,
    max_reactions: int = 3,
    library_mode: str = "cdn",
) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    molecule = run.molecule_index().get(mol_id)
    data = molecule.data if molecule else {"mol_id": mol_id}
    svg = write_molecule_svg(
        out / "structure_2d.svg",
        data.get("canonical_smiles") or data.get("isomeric_smiles"),
        legend=data.get("name") or mol_id,
    )
    candidates = run.read_table("candidate_summary")
    candidate_rows = (
        candidates[candidates["mol_id"].astype(str) == str(mol_id)]
        if not candidates.empty and "mol_id" in candidates.columns
        else candidates.head(0)
    )
    geometry = _first_geometry_for_molecule(run, mol_id)
    structure_view = write_html(
        out / "best_structure.html",
        "Best available structure",
        html_for_structure(
            geometry,
            "molecule_structure",
            "Best available structure",
            library_mode,
            asset_prefix="../../assets",
        ),
    )
    links = (
        _reaction_links(run, mol_id, max_reactions=max_reactions)
        if render_reactions
        else ""
    )
    identity = dataframe_to_html(pd.DataFrame([data]))
    body = (
        "<div class='grid'><div class='card'>"
        f"<h2>{esc(data.get('name', mol_id))}</h2>"
        f"<div class='svg-box'>{svg.read_text(encoding='utf-8')}</div></div>"
        f"<div class='card'><h2>Identity</h2>{identity}</div></div>"
        "<div class='card'><h2>Candidate summary</h2>"
        f"{dataframe_to_html(candidate_rows)}</div>"
        "<div class='card'><a href='best_structure.html'>"
        f"Best available 3D structure</a></div>{links}"
    )
    dossier = write_html(
        out / "dossier.html", f"Molecule dossier: {mol_id}", body
    )
    return {
        "dossier": dossier,
        "structure_2d": svg,
        "best_structure": structure_view,
    }
