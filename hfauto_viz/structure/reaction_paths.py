from __future__ import annotations

from pathlib import Path

from hfauto_viz.structure.xyz import frames_to_xyz, read_xyz_frames


def reaction_frame_xyz(run, reaction_id: str, out_path: str | Path) -> Path | None:
    """Create a simple reactant/TS/product path from run artifacts if available."""
    rxn = run.reaction_index().get(reaction_id)
    if not rxn:
        return None
    data = rxn.data
    sp_idx = run.species_index()
    frames = []
    for label, sid in [
        ("reactant", data.get("reactant_species_id")),
        ("TS", data.get("ts_species_id")),
        ("product", data.get("product_species_id")),
    ]:
        sp = sp_idx.get(str(sid)) if sid else None
        if sp:
            fs = read_xyz_frames(run.artifact_xyz_path(sp))
            if fs:
                frames.append((label, fs[0]))
    if not frames:
        return None
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(frames_to_xyz(frames), encoding="utf-8")
    return out
