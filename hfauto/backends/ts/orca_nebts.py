from __future__ import annotations

"""ORCA NEB-TS / OptTS / IRC backend for Phase 5.

The backend is intentionally auditable:
- NEB-TS, OptTS, and IRC inputs are rendered even in offline mode;
- ORCA subprocess execution is opt-in;
- fallback artifacts are explicitly marked as dummy/software-test-only;
- real ORCA outputs are parsed into common Calculation/IRC artifacts.
"""

from dataclasses import asdict
from pathlib import Path
from typing import Any

from hfauto.backends.qm.orca import OrcaEngine, allow_orca_subprocess, parse_orca_cartesian_blocks, parse_orca_output
from hfauto.backends.ts.base import IRCResult, TSSearchResult, canonical_species_id, make_ts_species_artifact, midpoint_ts_xyz
from hfauto.backends.ts.dummy import DummyTSEngine
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.hf_builder import XYZ, read_xyz, write_xyz
from hfauto.chemistry.reaction_path_qc import endpoint_pair_match_qc, estimate_pt_mode_overlap
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.qc import ts_qc
from hfauto.core.schemas.artifact import Artifact


def _reaction_id(reaction: Artifact) -> str:
    return str(reaction.data.get("reaction_id") or reaction.artifact_id)


def _ts_id_for_reaction(reaction: Artifact) -> str:
    rid = _reaction_id(reaction)
    return rid.replace("rxn_", "ts_", 1) if rid.startswith("rxn_") else f"ts_{rid}"


def _base_keywords(method: dict[str, Any]) -> list[str]:
    kw = list(method.get("base_keywords") or method.get("orca_keywords") or [])
    if not kw:
        mid = str(method.get("method_id", "r2SCAN-3c"))
        if mid.lower() in {"r2scan3c", "r2scan-3c", "r2scan3c_orca", "orca_nebts_r2scan3c"}:
            kw = ["r2SCAN-3c"]
        else:
            kw = [mid]
    forbidden = {"neb-ts", "fast-neb-ts", "tight-neb-ts", "zoom-neb-ts", "optts", "irc", "opt", "freq"}
    return [str(x) for x in kw if str(x).lower() not in forbidden]


def _job_keywords(method: dict[str, Any], job: str) -> list[str]:
    if job == "neb_ts" and method.get("neb_keywords"):
        return list(method["neb_keywords"])
    if job == "optts" and method.get("optts_keywords"):
        return list(method["optts_keywords"])
    if job == "irc" and method.get("irc_keywords"):
        return list(method["irc_keywords"])
    base = _base_keywords(method)
    if job == "neb_ts":
        base += [str(method.get("neb_job", method.get("neb_keyword", "NEB-TS")))]
        if method.get("freq_in_neb", False):
            base += ["Freq"]
    elif job == "optts":
        base += ["OptTS", "Freq"]
    elif job == "irc":
        base += ["IRC"]
    if "TightSCF" not in base and "tightscf" not in [x.lower() for x in base]:
        base += ["TightSCF"]
    return list(dict.fromkeys(base))


class ORCANEBTSInputRenderer:
    """Stateless input renderers used by backend and tests."""

    @staticmethod
    def render_nebts(reactant: Artifact, product: Artifact | None, method: dict[str, Any], product_xyz_name: str = "product.xyz") -> str:
        rc_path = method.get("reactant_xyz_name") or Path(str(reactant.data.get("xyz_path") or reactant.paths.get("xyz"))).name
        charge = int(reactant.data.get("charge", 0) or 0)
        mult = int(reactant.data.get("multiplicity", 1) or 1)
        nprocs = int(method.get("ncores", method.get("nprocs", 1)) or 1)
        mem_mb = int(method.get("memory_mb", method.get("maxcore_mb", 2000)) or 2000)
        nimages = method.get("n_images", method.get("nimages"))
        lines = [
            "# generated_by=hfauto",
            "# task=neb_ts",
            "! " + " ".join(_job_keywords(method, "neb_ts")),
            "",
            f"%pal nprocs {nprocs} end",
            f"%maxcore {mem_mb}",
            "%neb",
            f"  NEB_END_XYZFILE \"{product_xyz_name}\"",
            f"  PREOPT_ENDS {str(bool(method.get('preopt_ends', True))).upper()}",
        ]
        if nimages is not None:
            lines.append(f"  NImages {int(nimages)}")
        for k, v in (method.get("neb_settings", {}) or {}).items():
            lines.append(f"  {k} {v}")
        lines += ["END", ""]
        for block in method.get("neb_extra_blocks", []) or []:
            lines += [str(block).rstrip(), ""]
        lines += [f"* xyzfile {charge} {mult} {rc_path}", ""]
        return "\n".join(lines)

    @staticmethod
    def render_optts(ts_species: Artifact, method: dict[str, Any]) -> str:
        xyz_path = Path(str(ts_species.data.get("xyz_path") or ts_species.paths.get("xyz"))).resolve()
        charge = int(ts_species.data.get("charge", 0) or 0)
        mult = int(ts_species.data.get("multiplicity", 1) or 1)
        nprocs = int(method.get("ncores", method.get("nprocs", 1)) or 1)
        mem_mb = int(method.get("memory_mb", method.get("maxcore_mb", 2000)) or 2000)
        lines = [
            "# generated_by=hfauto",
            "# task=optts",
            "! " + " ".join(_job_keywords(method, "optts")),
            "",
            f"%pal nprocs {nprocs} end",
            f"%maxcore {mem_mb}",
        ]
        if method.get("calc_hess", True):
            lines += ["%geom", "  Calc_Hess true", "END"]
        lines += ["", f"* xyzfile {charge} {mult} {xyz_path}", ""]
        return "\n".join(lines)

    @staticmethod
    def render_irc(ts_species: Artifact, method: dict[str, Any]) -> str:
        xyz_path = Path(str(ts_species.data.get("xyz_path") or ts_species.paths.get("xyz"))).resolve()
        charge = int(ts_species.data.get("charge", 0) or 0)
        mult = int(ts_species.data.get("multiplicity", 1) or 1)
        nprocs = int(method.get("ncores", method.get("nprocs", 1)) or 1)
        mem_mb = int(method.get("memory_mb", method.get("maxcore_mb", 2000)) or 2000)
        max_iter = int(method.get("irc_maxiter", method.get("irc_max_iter", 80)))
        step = float(method.get("irc_step", method.get("irc_step_size", 0.10)))
        lines = [
            "# generated_by=hfauto",
            "# task=irc",
            "! " + " ".join(_job_keywords(method, "irc")),
            "",
            f"%pal nprocs {nprocs} end",
            f"%maxcore {mem_mb}",
            "%irc",
            f"  MaxIter {max_iter}",
            f"  StepSize {step}",
        ]
        for k, v in (method.get("irc_settings", {}) or {}).items():
            lines.append(f"  {k} {v}")
        lines += ["END", "", f"* xyzfile {charge} {mult} {xyz_path}", ""]
        return "\n".join(lines)


# Backward-compatible alias used by earlier tests/docs.
OrcaTSInputRenderer = ORCANEBTSInputRenderer


class ORCANEBTSEngine:
    name = "orca_nebts"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        return {**self.config, **(method or {})}

    def _executable(self, method: dict[str, Any]) -> str | None:
        return resolve_executable("orca", method.get("executable"))

    def render_neb_input(self, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> Path:
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        # Keep local copies for review/reproducibility; render uses absolute RC path and local product filename.
        write_xyz(read_xyz(reactant.data.get("xyz_path") or reactant.paths.get("xyz")), wd / "reactant.xyz")
        write_xyz(read_xyz(product.data.get("xyz_path") or product.paths.get("xyz")), wd / "product.xyz")
        inp = wd / "orca_nebts.inp"
        inp.write_text(ORCANEBTSInputRenderer.render_nebts(reactant, product, {**self._method(method), "reactant_xyz_name": "reactant.xyz"}, product_xyz_name="product.xyz"), encoding="utf-8")
        return inp

    def render_irc_input(self, ts_species: Artifact, method: dict[str, Any], workdir: str | Path) -> Path:
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        inp = wd / "orca_irc.inp"
        inp.write_text(ORCANEBTSInputRenderer.render_irc(ts_species, self._method(method)), encoding="utf-8")
        return inp

    def _path_artifact(self, reaction: Artifact, inp: Path, reactant: Artifact, product: Artifact, stage: str) -> Artifact:
        return Artifact(
            artifact_id=f"ts_path_{fingerprint_dict({'reaction': reaction.artifact_id, 'stage': stage, 'engine': self.name})}",
            artifact_type="ts_path",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            paths={"nebts_input": str(inp), "reactant_xyz": str(inp.parent / "reactant.xyz"), "product_xyz": str(inp.parent / "product.xyz")},
            method={"engine": self.name, "task": "neb_ts", "keywords": _job_keywords(self._method({}), "neb_ts")},
            data={"reaction_id": _reaction_id(reaction)},
            qc={"input_rendered": True},
        )

    def _dummy_fallback(self, reaction: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: Path, reason: str, path_art: Artifact) -> TSSearchResult:
        res = DummyTSEngine().search_ts(reaction, reactant, product, {**method, "method_id": method.get("method_id", "orca_nebts_fallback")}, workdir / "dummy_fallback")
        out = [path_art]
        for art in res.artifacts:
            art.parents = [path_art.artifact_id] + art.parents
            art.qc.update({"fallback_dummy": True, "fallback_reason": reason, "intended_engine": self.name, "scientific_use": "software_test_only_not_ts"})
            art.method = {**(art.method or {}), "engine": self.name, "intended_engine": self.name, "backend_fallback": "dummy"}
            if art.artifact_type == "reaction_validated":
                art.qc["ts_engine"] = self.name
                art.data["ts_backend"] = self.name
            if art.artifact_type == "calculation":
                art.data.update({"requested_engine": self.name, "calculation_level": "dummy_fallback_for_orca_nebts", "real_orca_executed": False})
        out.extend(res.artifacts)
        return TSSearchResult(artifacts=out, ts_species_id=res.ts_species_id, ts_calc_id=res.ts_calc_id, record={**res.record, "fallback_reason": reason})

    def _run_orca(self, inp: Path, method: dict[str, Any], stdout_name: str) -> tuple[Any, str]:
        exe = self._executable(method)
        if exe is None:
            raise FileNotFoundError("ORCA executable not found")
        result = run_command([str(exe), str(inp.name)], cwd=inp.parent, timeout_s=int(method.get("timeout_s", 86400)), env=method.get("env"), stdout_name=stdout_name, stderr_name=stdout_name.replace(".out", ".err"))
        out_text = Path(result.stdout_path).read_text(encoding="utf-8", errors="ignore")
        err_text = Path(result.stderr_path).read_text(encoding="utf-8", errors="ignore") if Path(result.stderr_path).exists() else ""
        if err_text.strip():
            out_text += "\nSTDERR\n" + err_text
        return result, out_text

    def _find_ts_xyz_from_output(self, output_text: str, reactant: Artifact, workdir: Path) -> Path | None:
        rc_symbols = read_xyz(reactant.data.get("xyz_path") or reactant.paths.get("xyz")).symbols
        candidates = []
        for pat in ["*TS*.xyz", "*ts*.xyz", "*NEB*.xyz", "*neb*.xyz", "*.xyz"]:
            candidates.extend(workdir.glob(pat))
        seen = set()
        for p in sorted(candidates, key=lambda x: x.stat().st_mtime if x.exists() else 0.0, reverse=True):
            if p.name in {"reactant.xyz", "product.xyz"} or p in seen:
                continue
            seen.add(p)
            try:
                if read_xyz(p).symbols == rc_symbols:
                    return p
            except Exception:
                pass
        blocks = parse_orca_cartesian_blocks(output_text)
        if blocks:
            out = workdir / "ts_from_orca_output.xyz"
            write_xyz(blocks[-1], out)
            return out
        return None

    def search_ts(self, reaction: Artifact, reactant: Artifact, product: Artifact, method_in: dict[str, Any], workdir: str | Path) -> TSSearchResult:
        method = self._method(method_in)
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        inp = self.render_neb_input(reactant, product, method, wd / "neb")
        path_art = self._path_artifact(reaction, inp, reactant, product, "ts-search")
        # Endpoint gate required for NEB/GSM.
        if read_xyz(reactant.data.get("xyz_path") or reactant.paths.get("xyz")).symbols != read_xyz(product.data.get("xyz_path") or product.paths.get("xyz")).symbols:
            fail = Artifact.failure("ts_failed_" + _reaction_id(reaction), "ts_result", "reactant/product atom order mismatch", category="endpoint_mismatch", parents=[path_art.artifact_id], recommended_fallback="rebuild_hf_endpoints")
            return TSSearchResult(artifacts=[path_art, fail], record=fail.model_dump())
        should_run = allow_orca_subprocess(method) and not method.get("dry_run", False) and self._executable(method) is not None
        if not should_run:
            reason = "dry_run=True" if method.get("dry_run", False) else ("ORCA execution disabled" if not allow_orca_subprocess(method) else "ORCA executable not found")
            if method.get("fallback_to_dummy", False):
                return self._dummy_fallback(reaction, reactant, product, method, wd, reason, path_art)
            fail = Artifact.failure("ts_failed_" + _reaction_id(reaction), "ts_result", reason + "; NEB-TS input was written", category="orca_nebts_not_run", parents=[path_art.artifact_id], recommended_fallback="set allow_subprocess=true, configure executable, or enable fallback_to_dummy", input=str(inp))
            return TSSearchResult(artifacts=[path_art, fail], record=fail.model_dump())
        result, output_text = self._run_orca(inp, method, "orca_nebts.out")
        path_art.paths.update({"nebts_output": result.stdout_path, "nebts_stderr": result.stderr_path, "command_result": str(Path(result.cwd) / "command_result.json")})
        path_art.provenance["neb_command"] = asdict(result)
        ts_xyz_guess = self._find_ts_xyz_from_output(output_text, reactant, Path(result.cwd))
        if not result.ok or ts_xyz_guess is None:
            fail = Artifact.failure("ts_failed_" + _reaction_id(reaction), "ts_result", f"ORCA NEB-TS failed or no TS geometry was found: returncode={result.returncode}", category="orca_nebts_failed", parents=[path_art.artifact_id], recommended_fallback="scan_optts or inspect ORCA NEB output", output=result.stdout_path)
            return TSSearchResult(artifacts=[path_art, fail], record=fail.model_dump())
        ts_id = _ts_id_for_reaction(reaction)
        ts_xyz = wd / f"{ts_id}_neb_guess.xyz"
        write_xyz(read_xyz(ts_xyz_guess), ts_xyz)
        ts_species = make_ts_species_artifact(reaction, reactant, product, ts_xyz, ts_id=ts_id, source=self.name, extra_qc={"real_neb_executed": True})
        opt_method = {**method, "orca_keywords": _job_keywords(method, "optts"), "method_id": method.get("optts_method_id", method.get("method_id", "orca_optts"))}
        calc = OrcaEngine().optimize_frequency(ts_species, opt_method, str(wd / "optts"))
        calc.method = {**(calc.method or {}), "engine": self.name, "stage": "ts-search", "task": "opt_freq", "backend_step": "optts_freq"}
        final_xyz = calc.paths.get("final_xyz") or ts_species.data["xyz_path"]
        ts_species.data["xyz_path"] = str(final_xyz)
        ts_species.paths["xyz"] = str(final_xyz)
        geom_qc = geometry_qc_from_xyz(ts_species.data, final_xyz)
        overlap, overlap_method = estimate_pt_mode_overlap(ts_species.data, final_xyz, calc.data.get("n_imag"), calc.data.get("imag_freq_cm1"), Path(calc.paths.get("output", "")).read_text(encoding="utf-8", errors="ignore") if calc.paths.get("output") and Path(calc.paths.get("output")).exists() else None)
        calc.qc.update({"mode_overlap_score": overlap, "mode_overlap_method": overlap_method, "imag_mode_matches_reaction_coordinate": overlap >= float(method.get("mode_overlap_threshold", 0.7)), "geometry_qc": geom_qc})
        calc.qc.update(ts_qc(calc.data.get("n_imag"), calc.data.get("imag_freq_cm1"), overlap))
        calc.data.update({"mode_overlap_score": overlap, "mode_overlap_method": overlap_method, "ts_backend": self.name, **{k: v for k, v in geom_qc.items() if v is not None}})
        rv = Artifact(
            artifact_id=reaction.artifact_id + "_with_ts",
            artifact_type="reaction_validated",
            parents=[reaction.artifact_id, path_art.artifact_id, ts_species.artifact_id, calc.artifact_id],
            data={**reaction.data, "ts_species_id": ts_species.artifact_id, "ts_calc_id": calc.artifact_id, "ts_backend": self.name},
            qc={"ts_engine": self.name, "ts_found": calc.status.status == "success", "ts_validated_by_frequency": bool(calc.qc.get("ts_validated_by_frequency")), "n_imag": calc.data.get("n_imag"), "imag_freq_cm1": calc.data.get("imag_freq_cm1"), "mode_overlap_score": overlap, "real_ts_search_executed": True, "fallback_dummy": bool(calc.qc.get("fallback_dummy", False))},
        )
        return TSSearchResult(artifacts=[path_art, ts_species, calc, rv], ts_species_id=ts_species.artifact_id, ts_calc_id=calc.artifact_id, record=rv.data)

    def _dummy_irc_fallback(self, reaction: Artifact, ts_species: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: Path, reason: str, path_art: Artifact) -> IRCResult:
        res = DummyTSEngine().run_irc(reaction, ts_species, reactant, product, method, workdir / "dummy_irc")
        for art in res.artifacts:
            art.parents = [path_art.artifact_id] + art.parents
            art.qc.update({"fallback_dummy": True, "fallback_reason": reason, "intended_engine": self.name, "scientific_use": "software_test_only_not_irc"})
            art.method = {**(art.method or {}), "engine": self.name, "intended_engine": self.name, "backend_fallback": "dummy"}
            art.data.update({"requested_engine": self.name, "fallback_reason": reason})
        return IRCResult(artifacts=[path_art] + res.artifacts, record={**res.record, "fallback_reason": reason})

    def _find_irc_endpoints(self, workdir: Path, output_text: str) -> tuple[Path | None, Path | None]:
        # Prefer explicit XYZ files if ORCA emitted them, otherwise use first/last coordinate block.
        xyzs = sorted([p for p in workdir.glob("*.xyz") if "endpoint" in p.name.lower() or "irc" in p.name.lower()], key=lambda p: p.stat().st_mtime if p.exists() else 0.0)
        if len(xyzs) >= 2:
            return xyzs[0], xyzs[-1]
        blocks = parse_orca_cartesian_blocks(output_text)
        if len(blocks) >= 2:
            a = workdir / "irc_endpoint_a.xyz"
            b = workdir / "irc_endpoint_b.xyz"
            write_xyz(blocks[0], a)
            write_xyz(blocks[-1], b)
            return a, b
        return None, None

    def run_irc(self, reaction: Artifact, ts_species: Artifact, reactant: Artifact, product: Artifact, method_in: dict[str, Any], workdir: str | Path) -> IRCResult:
        method = self._method(method_in)
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        inp = self.render_irc_input(ts_species, method, wd)
        path_art = Artifact(
            artifact_id="irc_attempt_" + fingerprint_dict({"reaction": reaction.artifact_id, "ts": ts_species.artifact_id, "engine": self.name}),
            artifact_type="irc_attempt",
            parents=[reaction.artifact_id, ts_species.artifact_id],
            paths={"irc_input": str(inp)},
            method={"engine": self.name, "task": "irc", "keywords": _job_keywords(method, "irc")},
            data={"reaction_id": reaction.data.get("reaction_id", reaction.artifact_id), "ts_species_id": ts_species.artifact_id},
            qc={"input_rendered": True},
        )
        should_run = allow_orca_subprocess(method) and not method.get("dry_run", False) and self._executable(method) is not None
        if not should_run:
            reason = "dry_run=True" if method.get("dry_run", False) else ("ORCA execution disabled" if not allow_orca_subprocess(method) else "ORCA executable not found")
            if method.get("fallback_to_dummy", False):
                return self._dummy_irc_fallback(reaction, ts_species, reactant, product, method, wd, reason, path_art)
            fail = Artifact.failure("irc_failed_" + reaction.artifact_id, "irc", reason + "; IRC input was written", category="orca_irc_not_run", parents=[path_art.artifact_id], recommended_fallback="set allow_subprocess=true, configure executable, or enable fallback_to_dummy", input=str(inp))
            return IRCResult(artifacts=[path_art, fail], record=fail.model_dump())
        result, output_text = self._run_orca(inp, method, "orca_irc.out")
        path_art.paths.update({"irc_output": result.stdout_path, "irc_stderr": result.stderr_path, "command_result": str(Path(result.cwd) / "command_result.json")})
        parsed = parse_orca_output(output_text)
        fwd, bwd = self._find_irc_endpoints(Path(result.cwd), output_text)
        qc = endpoint_pair_match_qc(fwd, bwd, reactant.data.get("xyz_path") or reactant.paths.get("xyz"), product.data.get("xyz_path") or product.paths.get("xyz"), ts_species.data, float(method.get("endpoint_rmsd_threshold_A", 0.75)))
        qc.update({"normal_termination": bool(parsed.get("normal_termination")), "real_irc_executed": bool(result.ok), "fallback_dummy": False, "command_ok": bool(result.ok)})
        qc["irc_validated"] = bool(qc.get("irc_validated") and result.ok)
        rec = {"reaction_id": reaction.data.get("reaction_id", reaction.artifact_id.replace("_with_ts", "")), "ts_species_id": ts_species.artifact_id, "irc_backend": self.name, "forward_endpoint_xyz": str(fwd) if fwd else None, "backward_endpoint_xyz": str(bwd) if bwd else None, "output": result.stdout_path, **qc}
        irc_id = str(reaction.artifact_id).replace("_with_ts", "_irc") if "_with_ts" in reaction.artifact_id else f"irc_{reaction.artifact_id}"
        artifact = Artifact(
            artifact_id=irc_id,
            artifact_type="irc",
            parents=[reaction.artifact_id, ts_species.artifact_id, reactant.artifact_id, product.artifact_id, path_art.artifact_id],
            paths={"input": str(inp), "output": result.stdout_path, "forward_xyz": str(fwd) if fwd else "", "backward_xyz": str(bwd) if bwd else "", "command_result": str(Path(result.cwd) / "command_result.json")},
            data=rec,
            method={"engine": self.name, "stage": "irc", "method_id": method.get("method_id", "orca_irc"), "task": "irc", "keywords": _job_keywords(method, "irc")},
            qc=qc,
            provenance={"command": asdict(result), "created_by": "ORCANEBTSEngine"},
        )
        return IRCResult(artifacts=[path_art, artifact], record=rec)
