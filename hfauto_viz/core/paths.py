from __future__ import annotations
from pathlib import Path


def safe_name(value: str, max_len: int = 220) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))[:max_len]


def ensure_dir(path: str | Path) -> Path:
    p = Path(path); p.mkdir(parents=True, exist_ok=True); return p


def latest_manifest_path(run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    pointer = run_dir / "manifest.path"
    if pointer.exists():
        raw = pointer.read_text(encoding="utf-8").strip()
        p = Path(raw)
        candidates = [p, run_dir / raw]
        if not p.is_absolute():
            candidates.append(run_dir.parent.parent / p)
        for c in candidates:
            if c.exists():
                return c.resolve()
    manifests = sorted(run_dir.glob("*/manifest.json"), key=lambda p: p.stat().st_mtime)
    if not manifests:
        raise FileNotFoundError(f"No manifest found under {run_dir}")
    return manifests[-1].resolve()


def resolve_artifact_path(path_value: str | Path | None, run_dir: str | Path) -> Path | None:
    if not path_value:
        return None
    p = Path(path_value)
    if p.is_absolute():
        return p
    run_dir = Path(run_dir)
    candidates = [Path.cwd() / p, run_dir / p, run_dir.parent.parent / p]
    if run_dir.name in p.parts:
        idx = p.parts.index(run_dir.name)
        candidates.append(run_dir.joinpath(*p.parts[idx + 1:]))
    for c in candidates:
        if c.exists():
            return c.resolve()
    return candidates[0]
