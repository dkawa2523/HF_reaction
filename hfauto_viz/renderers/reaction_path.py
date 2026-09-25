from __future__ import annotations

from pathlib import Path

from hfauto_viz.structure.xyz import frames_to_xyz, parse_xyz_frames, reaction_path_frames


def build_reactant_ts_product_path(reactant_xyz: str, ts_xyz: str, product_xyz: str, out_path: str | Path, n_each: int = 8) -> Path:
    """Build a conservative reactant/TS/product multi-frame XYZ for visualization."""
    r = parse_xyz_frames(reactant_xyz)
    t = parse_xyz_frames(ts_xyz)
    p = parse_xyz_frames(product_xyz)
    frames = reaction_path_frames(r[0] if r else None, t[0] if t else None, p[0] if p else None)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(frames_to_xyz(frames), encoding="utf-8")
    return out
