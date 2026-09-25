from __future__ import annotations

import platform
import shutil
import sys
from collections.abc import Mapping
from importlib import import_module, metadata
from pathlib import Path
from typing import Any

from hfauto.core.config import load_yaml
from hfauto.core.executables import resolve_executable

PYTHON_DISTRIBUTIONS = {
    "hfauto": "hfauto",
    "pydantic": "pydantic",
    "pandas": "pandas",
    "numpy": "numpy",
    "rdkit": "rdkit",
    "requests": "requests",
    "bs4": "beautifulsoup4",
    "plotly": "plotly",
    "graphviz": "graphviz",
    "cantera": "cantera",
    "goodvibes": "goodvibes",
    "pysisyphus": "pysisyphus",
    "qcengine": "qcengine",
    "scine_readuct": "scine-readuct",
    "scine_xtb_wrapper": "scine-xtb-wrapper",
}


def _distribution_status(import_name: str, distribution_name: str) -> dict[str, Any]:
    try:
        version = metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        return {
            "name": import_name,
            "distribution": distribution_name,
            "available": False,
            "version": None,
            "error": "distribution_not_installed",
        }
    try:
        import_module(import_name)
    except Exception as exc:
        return {
            "name": import_name,
            "distribution": distribution_name,
            "available": False,
            "version": version,
            "error": f"{type(exc).__name__}: {exc}",
        }
    return {
        "name": import_name,
        "distribution": distribution_name,
        "available": True,
        "version": version,
        "error": None,
    }


def _configured_external_requirements(config: dict[str, Any]) -> set[str]:
    required: set[str] = set()
    for stage in config.get("stages", []):
        if not stage.get("enabled", True):
            continue
        name = stage.get("name")
        engine = stage.get("engine")
        if name == "conformers" and engine == "crest":
            required.update({"crest", "xtb"})
        if name == "build-complexes":
            compositions = list(stage.get("compositions") or [])
            uses_crest_nci = any(
                composition.get("crest_nci", _component_count(composition) > 1)
                for composition in compositions
            )
            if uses_crest_nci:
                required.update({"crest", "xtb"})
        if name == "preopt" and engine == "xtb":
            required.add("xtb")
        if name in {"dft-minima", "sp"} and engine in {"nwchem", "orca"}:
            required.add(str(engine))
            if engine == "nwchem" and int((stage.get("settings", {}) or {}).get("ncores", 1)) > 1:
                required.add("mpirun")
        if name in {"ts-search", "irc"}:
            specs = stage.get("engine_order") or [engine]
            if any("orca" in str(spec) for spec in specs):
                required.add("orca")
            if any("nwchem" in str(spec) for spec in specs):
                required.add("nwchem")
                settings = stage.get("settings", {}) or {}
                if int(settings.get("ncores", 1)) > 1:
                    required.add("mpirun")
            if any("pysisyphus" in str(spec) for spec in specs):
                required.add("pysis")
        if name == "thermo" and engine == "goodvibes":
            settings = stage.get("settings", {}) or {}
            if settings.get("use_external_goodvibes") or settings.get("run_external"):
                required.add("goodvibes")
    return required


def _component_count(composition: Mapping[str, Any]) -> int:
    return sum(
        max(0, int(component.get("count", 1)))
        for component in (composition.get("components") or [])
    )


def _configured_python_requirements(config: Mapping[str, Any]) -> set[str]:
    """Return import names required by the enabled workflow stages."""

    required = {"hfauto", "pydantic", "pandas", "numpy"}
    stage_names = {
        str(stage.get("name") or "").replace("_", "-")
        for stage in config.get("stages", [])
        if stage.get("enabled", True)
    }
    if stage_names.intersection(
        {"ingest", "enrich", "detect-sites", "conformers", "build-complexes"}
    ):
        required.add("rdkit")
    uses_readuct = any(
        str(stage.get("name") or "").replace("_", "-")
        == "explore-reactions"
        and str(stage.get("engine") or "readuct").lower() == "readuct"
        for stage in config.get("stages", [])
        if stage.get("enabled", True)
    )
    if uses_readuct:
        required.update({"scine_readuct", "scine_xtb_wrapper"})

    required_external = _configured_external_requirements(dict(config))
    if "pysis" in required_external:
        required.update({"pysisyphus", "qcengine"})
    return required


def _configured_external_paths(config: dict[str, Any]) -> dict[str, str]:
    """Collect explicit executable paths from the stages that own them."""

    paths: dict[str, str] = {}
    for stage in config.get("stages", []):
        if not stage.get("enabled", True):
            continue
        name = str(stage.get("name") or "")
        engine = str(stage.get("engine") or "")
        specs = stage.get("engine_order") or [engine]
        settings = {
            **(stage.get("settings") or {}),
            **(stage.get("method_settings") or {}),
        }
        uses_nwchem = bool(
            engine == "nwchem"
            or any("nwchem" in str(spec) for spec in specs)
            or (
                name == "irc"
                and str(settings.get("program") or "").lower() == "nwchem"
            )
        )
        if uses_nwchem and settings.get("executable"):
            paths["nwchem"] = str(settings["executable"])
        if uses_nwchem and settings.get("mpi_executable"):
            paths["mpirun"] = str(settings["mpi_executable"])
        if any("pysisyphus" in str(spec) for spec in specs) and settings.get(
            "pysis_executable"
        ):
            paths["pysis"] = str(settings["pysis_executable"])
        if name == "thermo" and engine == "goodvibes" and settings.get(
            "executable"
        ):
            paths["goodvibes"] = str(settings["executable"])
        if engine in {"crest", "xtb"}:
            executable = settings.get("executable") or stage.get("executable")
            if executable:
                paths[engine] = str(executable)
        if name == "conformers" and engine == "crest":
            xtb_executable = settings.get("xtb_executable") or stage.get(
                "xtb_executable"
            )
            if xtb_executable:
                paths["xtb"] = str(xtb_executable)
        if name == "build-complexes":
            nci_settings = dict(stage.get("crest_nci_settings") or {})
            if nci_settings.get("executable"):
                paths["crest"] = str(nci_settings["executable"])
            if nci_settings.get("xtb_executable"):
                paths["xtb"] = str(nci_settings["xtb_executable"])
    return paths


def environment_report(config_path: str | Path | Mapping[str, Any] | None = None) -> dict[str, Any]:
    config = dict(config_path) if isinstance(config_path, Mapping) else (load_yaml(config_path) if config_path else {})
    packages = [_distribution_status(k, v) for k, v in PYTHON_DISTRIBUTIONS.items()]
    configured_paths = _configured_external_paths(config)
    executables = []
    for name in ["orca", "nwchem", "pysis", "mpirun", "xtb", "crest", "goodvibes", "dot", "sbatch", "qsub"]:
        path = resolve_executable(name, configured_paths.get(name))
        if name == "mpirun" and path is None:
            nwchem = resolve_executable("nwchem")
            sibling = Path(nwchem).resolve().parent / "mpirun" if nwchem else None
            path = str(sibling) if sibling and sibling.exists() else None
        executables.append({"name": name, "available": path is not None, "path": path})

    required_external = _configured_external_requirements(config)
    required_python = _configured_python_requirements(config)
    package_map = {item["name"]: item for item in packages}
    executable_map = {item["name"]: item for item in executables}
    missing_python = sorted(name for name in required_python if not package_map[name]["available"])
    missing_external = sorted(name for name in required_external if not executable_map[name]["available"])
    if "goodvibes" in required_external and package_map["goodvibes"]["available"]:
        missing_external = [name for name in missing_external if name != "goodvibes"]

    return {
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "config": str(config_path) if config_path and not isinstance(config_path, Mapping) else None,
        "packages": packages,
        "executables": executables,
        "configured_external_paths": configured_paths,
        "graphviz_dot": shutil.which("dot"),
        "required_python": sorted(required_python),
        "required_external": sorted(required_external),
        "missing_python": missing_python,
        "missing_external": missing_external,
        "ready": not missing_python and not missing_external,
    }
