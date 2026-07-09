from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.core.schemas.artifact import Artifact
from hfauto_viz.data.loaders import RunData
from hfauto_viz.structure.xyz import Frame, read_xyz_frames, interpolate_frames


def _path_from_artifact(run: RunData, art: Artifact | None, keys: list[str]) -> Path | None:
    if art is None:
        return None
    for key in keys:
        raw = art.paths.get(key) or art.data.get(key)
        p = run.resolve_path(raw) if raw else None
        if p and p.exists():
            return p
    return None


def reaction_path_frames_from_run(
    run: RunData,
    reaction_id: str,
    reactant: Frame | None,
    ts: Frame | None,
    product: Frame | None,
    interpolation_points: int = 12,
) -> tuple[list[tuple[str, Frame]], dict[str, Any]]:
    """Prefer real IRC/NEB/scan frames, fallback to smooth interpolation."""
    source: dict[str, Any] = {"source": "interpolated_reactant_ts_product", "real_path_frames": False}
    # Search IRC/path artifacts first.
    for art_type in ["irc", "irc_attempt", "reaction_path_validated", "ts_path"]:
        for art in run.latest_artifacts(art_type):
            if str(art.data.get("reaction_id") or "") != str(reaction_id):
                continue
            p = _path_from_artifact(run, art, [
                "path_xyz", "irc_path_xyz", "irc_xyz", "trajectory_xyz", "neb_path_xyz", "scan_xyz", "xyz",
            ])
            if p:
                frames = read_xyz_frames(p)
                if frames:
                    source = {
                        "source": f"artifact:{art.artifact_type}",
                        "artifact_id": art.artifact_id,
                        "path": str(p),
                        "real_path_frames": bool(art.qc.get("real_irc_executed") or art.qc.get("real_orca_executed") or art.data.get("real_irc_executed")),
                        "n_frames": len(frames),
                    }
                    return [(f"{art_type}_{i:03d}", fr) for i, fr in enumerate(frames)], source
    frames: list[tuple[str, Frame]] = []
    if reactant and ts and len(reactant) == len(ts):
        frames.extend(interpolate_frames(reactant, ts, n=interpolation_points, prefix="reactant_to_ts"))
    elif reactant:
        frames.append(("reactant", reactant))
    if ts and product and len(ts) == len(product):
        second = interpolate_frames(ts, product, n=interpolation_points, prefix="ts_to_product")
        frames.extend(second[1:] if frames else second)
    else:
        if ts and not any(label.startswith("ts") or "ts" in label for label, _ in frames):
            frames.append(("transition_state", ts))
        if product:
            frames.append(("product", product))
    if not frames:
        for label, frame in [("reactant", reactant), ("transition_state", ts), ("product", product)]:
            if frame:
                frames.append((label, frame))
    source["n_frames"] = len(frames)
    return frames, source
