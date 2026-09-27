"""pysisyphus side of ``pysis_neb``; runs in ``python -m hfauto.execution.worker``.

pysisyphus is imported only inside ``run_neb``, so only the worker process loads it. The
returned dict becomes the attempt's ``result.json``; the NEB's images are pysisyphus's own
trajectories (Å): ``final_geometries.trj``, written whether or not the NEB converged, else the
last cycle's ``current_geometries.trj`` when an error stopped it.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from hfauto.core.constants import BOHR_TO_ANGSTROM

IMAGES = ("final_geometries.trj", "current_geometries.trj")
_XTB_VERSION = re.compile(r"xtb version (\S+)")


def neb_images(workdir: Path) -> Path | None:
    """The NEB's final images, else those of its last cycle; None before its first cycle."""
    return next((workdir / name for name in IMAGES if (workdir / name).is_file()), None)


def run_neb(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """Run ``job["run_dict"]`` (CI-NEB, then the TS optimization from its climbing image).

    When the run raises after a NEB cycle (mostly the TS optimization), the images still stand:
    ``ts`` is None and ``error`` says why. A run that raises before any images propagates.
    """
    from pysisyphus.run import run_from_dict

    try:
        result = run_from_dict(job["run_dict"], cwd=workdir)
    except Exception as exc:
        if neb_images(workdir) is None:
            raise
        return {"ts": None, "error": f"{type(exc).__name__}: {exc}",
                "xtb_version": _xtb_version(workdir)}
    return {"ts": ts_coords(result), "xtb_version": _xtb_version(workdir)}


def ts_coords(result: Any) -> list[list[float]] | None:
    """The TS (Å) of a pysisyphus ``RunResult`` when its optimization converged; pysisyphus
    leaves it unset when the splined HEI is the first or last image."""
    if result.ts_geom is None or result.ts_opt is None or not result.ts_opt.is_converged:
        return None
    flat = [float(x) * BOHR_TO_ANGSTROM for x in result.ts_geom.cart_coords]
    return [flat[i : i + 3] for i in range(0, len(flat), 3)]


def _xtb_version(workdir: Path) -> str | None:
    """The xTB version printed in the calculator outputs pysisyphus keeps in qm_calcs/."""
    for path in sorted((workdir / "qm_calcs").glob("*.out")):
        found = _XTB_VERSION.search(path.read_text(encoding="utf-8", errors="replace"))
        if found:
            return found.group(1)
    return None
