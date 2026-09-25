"""pysisyphus side of ``pysis_gs``; runs in ``python -m hfauto.execution.worker``.

pysisyphus is imported only inside ``run_growing_string``, so only the worker process loads
it. The returned dict becomes the attempt's ``result.json``; coordinates are in Å.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from hfauto.core.constants import BOHR_TO_ANGSTROM

_XTB_VERSION = re.compile(r"xtb version (\S+)")


def run_growing_string(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """Run the pysisyphus input ``job["run_dict"]`` (cos + optional tsopt) in ``workdir``."""
    from pysisyphus.run import run_from_dict

    return summarize(run_from_dict(job["run_dict"], cwd=workdir), workdir)


def summarize(result: Any, workdir: Path) -> dict[str, Any]:
    """Images, energies, optimizer history and the TS of a pysisyphus ``RunResult``.

    The TS is reported only when its optimization converged; pysisyphus leaves it unset
    when the splined HEI is the first or last image.
    """
    cos, opt = result.cos, result.cos_opt
    ts = None
    if result.ts_geom is not None and result.ts_opt is not None and result.ts_opt.is_converged:
        ts = {"coords": _angstrom(result.ts_geom.cart_coords),
              "energy": float(result.ts_geom.energy)}
    climbing = int(cos.get_hei_index()) if getattr(cos, "started_climbing", False) else None
    return {
        "images": [_angstrom(image.cart_coords) for image in cos.images],
        "energies": [float(image.energy) for image in cos.images],
        "max_forces": [float(f) for f in opt.max_forces],
        "converged": bool(opt.is_converged),
        "climbing_image": climbing,
        "ts": ts,
        "xtb_version": _xtb_version(workdir),
    }


def _angstrom(coords: Any) -> list[list[float]]:
    flat = [float(x) * BOHR_TO_ANGSTROM for x in coords]
    return [flat[i : i + 3] for i in range(0, len(flat), 3)]


def _xtb_version(workdir: Path) -> str | None:
    """The xTB version printed in the calculator outputs pysisyphus keeps in qm_calcs/."""
    for path in sorted((workdir / "qm_calcs").glob("*.out")):
        found = _XTB_VERSION.search(path.read_text(encoding="utf-8", errors="replace"))
        if found:
            return found.group(1)
    return None
