"""SCINE ReaDuct 6.1 adapter for bounded AFIR and NT2 exploration."""

from __future__ import annotations

from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.backends.reaction_discovery.base import DiscoveryResult
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.chemistry.xyz_trajectory import write_xyz_trajectory
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import ReactionTrialRecord


class _TrajectoryCollector:
    def __init__(self) -> None:
        self.images: list[XYZ] = []
        self.energies: list[float | None] = []

    def collect(self, cycle: int, atoms, results, tag: str) -> None:
        import scine_utilities as utils

        symbols = [utils.ElementInfo.symbol(element) for element in atoms.elements]
        energy_value = getattr(results, "energy", None)
        energy = float(energy_value) if energy_value is not None else None
        self.images.append(
            XYZ(
                symbols=symbols,
                coords=np.asarray(atoms.positions, dtype=float)
                / float(utils.BOHR_PER_ANGSTROM),
                comment=f"readuct cycle={cycle} tag={tag} energy_hartree={energy}",
            )
        )
        self.energies.append(energy)


def _calculator_energy(calculator) -> float | None:
    if not calculator.has_results():
        return None
    results = calculator.get_results()
    value = getattr(results, "energy", None)
    return float(value) if value is not None else None


def _write_system(systems: dict, name: str, path: Path) -> str:
    import scine_utilities as utils

    path.parent.mkdir(parents=True, exist_ok=True)
    utils.io.write(str(path), systems[name].structure, f"readuct_system={name}")
    return str(path)


def _write_collector(collector: _TrajectoryCollector, path: Path) -> str | None:
    if not collector.images:
        return None
    if len(collector.images) == 1:
        collector.images.append(
            XYZ(
                symbols=list(collector.images[0].symbols),
                coords=collector.images[0].coords.copy(),
                comment=collector.images[0].comment + " duplicate_final_frame",
            )
        )
    return str(write_xyz_trajectory(collector.images, path))


def _flat_pairs(pairs) -> list[int]:
    return [atom for pair in pairs for atom in pair.atoms]


def _afir_bias_target(
    trial: ReactionTrialRecord,
) -> tuple[Any | None, bool | None, dict[str, Any]]:
    """Choose one unambiguous AFIR force; other trial terms remain validation targets."""

    if len(trial.associations) == 1:
        return (
            trial.associations[0],
            True,
            {
                "afir_bias_target": "association",
                "directly_biased_pair": list(trial.associations[0].atoms),
                "unbiased_dissociation_pair_count": len(trial.dissociations),
            },
        )
    if not trial.associations and len(trial.dissociations) == 1:
        return (
            trial.dissociations[0],
            False,
            {
                "afir_bias_target": "dissociation",
                "directly_biased_pair": list(trial.dissociations[0].atoms),
                "unbiased_dissociation_pair_count": 0,
            },
        )
    return (
        None,
        None,
        {
            "reason": "afir_requires_one_unambiguous_bias_pair",
            "association_pair_count": len(trial.associations),
            "dissociation_pair_count": len(trial.dissociations),
        },
    )


class ReaDuctDiscoveryBackend:
    """Run real GFN2-xTB discovery; never supplies a dummy fallback."""

    name = "readuct"

    def __init__(self, **defaults: Any) -> None:
        self.defaults = defaults

    @staticmethod
    def capability() -> dict[str, Any]:
        try:
            import scine_readuct  # noqa: F401
            import scine_utilities  # noqa: F401
            import scine_xtb_wrapper  # noqa: F401
        except (ImportError, ModuleNotFoundError) as exc:
            return {"available": False, "reason": str(exc)}
        return {
            "available": True,
            "readuct_version": metadata.version("scine-readuct"),
            "utilities_version": metadata.version("scine-utilities"),
            "xtb_wrapper_version": metadata.version("scine-xtb-wrapper"),
        }

    def _load_system(self, source: Artifact, settings: dict[str, Any]):
        path = str(source.data.get("xyz_path") or source.paths["xyz"])
        return self._load_xyz(
            path,
            charge=int(source.data.get("charge", 0) or 0),
            multiplicity=int(source.data.get("multiplicity", 1) or 1),
            settings=settings,
        )

    @staticmethod
    def _load_xyz(
        path: str | Path,
        *,
        charge: int,
        multiplicity: int,
        settings: dict[str, Any],
    ):
        import scine_utilities as utils
        import scine_xtb_wrapper  # noqa: F401 - registers the calculator module

        return utils.core.load_system_into_calculator(
            str(path),
            str(settings.get("method_family", "GFN2")),
            program="XTB",
            molecular_charge=int(charge),
            spin_multiplicity=int(multiplicity),
        )

    def _relax_nt2_driven_endpoint(
        self,
        collector: _TrajectoryCollector,
        trial: ReactionTrialRecord,
        cfg: dict[str, Any],
        directory: Path,
    ) -> tuple[str | None, bool, dict[str, Any]]:
        """Remove the drive and relax the final NT2 frame as a product candidate."""

        import scine_readuct as readuct

        if not collector.images:
            return None, False, {"reason": "nt2_trajectory_has_no_frames"}
        seed_path = write_xyz(
            collector.images[-1], directory / "nt2_driven_endpoint_seed.xyz"
        )
        systems = {
            "driven_endpoint": self._load_xyz(
                seed_path,
                charge=trial.charge,
                multiplicity=trial.multiplicity,
                settings=cfg,
            )
        }
        systems, converged = readuct.run_opt_task(
            systems,
            ["driven_endpoint"],
            output=["product_unbiased"],
            convergence_max_iterations=int(
                cfg.get("endpoint_max_iterations", 200)
            ),
            stop_on_error=False,
        )
        if "product_unbiased" not in systems:
            return str(seed_path), False, {
                "reason": "nt2_driven_endpoint_unbiased_optimization_failed",
                "unbiased_opt_converged": bool(converged),
            }
        product_path = _write_system(
            systems,
            "product_unbiased",
            directory / "nt2_product_unbiased.xyz",
        )
        return product_path, bool(converged), {
            "driven_endpoint_seed": str(seed_path),
            "unbiased_opt_converged": bool(converged),
            "electronic_energy_hartree": _calculator_energy(
                systems["product_unbiased"]
            ),
        }

    def explore(
        self,
        trial: ReactionTrialRecord,
        source: Artifact,
        settings: dict[str, Any],
        workdir: str | Path,
        *,
        driver: str,
    ) -> DiscoveryResult:
        cfg = {**self.defaults, **settings, **trial.settings}
        capability = self.capability()
        if not capability["available"]:
            return DiscoveryResult(
                driver=driver,
                success=False,
                failure_reason=f"readuct_unavailable:{capability['reason']}",
                data={"capability": capability},
            )
        directory = Path(workdir)
        directory.mkdir(parents=True, exist_ok=True)
        try:
            if driver == "nt2":
                return self._run_nt2(trial, source, cfg, directory, capability)
            if driver == "afir":
                return self._run_afir(trial, source, cfg, directory, capability)
            raise ValueError(f"unsupported ReaDuct discovery driver: {driver}")
        except Exception as exc:
            return DiscoveryResult(
                driver=driver,
                success=False,
                failure_reason=f"{type(exc).__name__}:{exc}",
                data={"capability": capability},
            )

    def _run_nt2(
        self,
        trial: ReactionTrialRecord,
        source: Artifact,
        cfg: dict[str, Any],
        directory: Path,
        capability: dict[str, Any],
    ) -> DiscoveryResult:
        import scine_readuct as readuct

        systems = {"start": self._load_system(source, cfg)}
        collector = _TrajectoryCollector()
        systems, nt_ok = readuct.run_nt2_task(
            systems,
            ["start"],
            False,
            [collector.collect],
            output=["nt2_guess"],
            nt_associations=_flat_pairs(trial.associations),
            nt_dissociations=_flat_pairs(trial.dissociations),
            nt_total_force_norm=float(cfg.get("nt_total_force_norm", 0.1)),
            convergence_max_iterations=int(cfg.get("nt_max_iterations", 200)),
            stop_on_error=False,
        )
        trajectory = _write_collector(collector, directory / "nt2_trajectory.xyz")
        if "nt2_guess" not in systems:
            return DiscoveryResult(
                driver="nt2",
                success=False,
                trajectory_path=trajectory,
                failure_reason="nt2_did_not_produce_a_ts_guess",
                data={"nt2_converged": nt_ok, "capability": capability},
            )
        guess_path = _write_system(
            systems, "nt2_guess", directory / "nt2_guess.xyz"
        )
        configured_order = cfg.get("ts_optimizer_order")
        if configured_order is not None:
            optimizer_order = list(configured_order)
        elif cfg.get("ts_optimizer") is not None:
            optimizer_order = [str(cfg["ts_optimizer"])]
        else:
            optimizer_order = ["bofill", "evf", "dimer"]
        ts_attempts: list[dict[str, Any]] = []
        ts_name: str | None = None
        ts_ok = False
        for optimizer in dict.fromkeys(str(value) for value in optimizer_order):
            output_name = f"ts_{optimizer}"
            systems, converged = readuct.run_tsopt_task(
                systems,
                ["nt2_guess"],
                output=[output_name],
                optimizer=optimizer,
                convergence_max_iterations=int(
                    cfg.get("ts_max_iterations", 200)
                ),
                stop_on_error=False,
            )
            produced = output_name in systems
            ts_attempts.append(
                {
                    "optimizer": optimizer,
                    "converged": bool(converged),
                    "system_produced": produced,
                }
            )
            if converged and produced:
                ts_name = output_name
                ts_ok = True
                break
        if ts_name is None:
            product_path, product_ok, recovery = (
                self._relax_nt2_driven_endpoint(
                    collector, trial, cfg, directory
                )
            )
            if product_ok and product_path is not None:
                return DiscoveryResult(
                    driver="nt2",
                    success=True,
                    product_xyz_path=product_path,
                    trajectory_path=trajectory,
                    electronic_energy_hartree=recovery.get(
                        "electronic_energy_hartree"
                    ),
                    low_level_ts_validated=False,
                    low_level_irc_connected=False,
                    biased_energy_used_as_barrier=False,
                    paths={
                        "product": product_path,
                        "nt2_guess": guess_path,
                        "driven_endpoint_seed": str(
                            recovery["driven_endpoint_seed"]
                        ),
                        **({"trajectory": trajectory} if trajectory else {}),
                    },
                    data={
                        "nt2_converged": nt_ok,
                        "ts_converged": False,
                        "ts_optimizer_attempts": ts_attempts,
                        "endpoint_recovered_without_ts_claim": True,
                        **recovery,
                    },
                )
            return DiscoveryResult(
                driver="nt2",
                success=False,
                trajectory_path=trajectory,
                failure_reason="nt2_ts_optimization_failed",
                paths={
                    "nt2_guess": guess_path,
                    **({"trajectory": trajectory} if trajectory else {}),
                },
                data={
                    "nt2_converged": nt_ok,
                    "ts_converged": False,
                    "ts_optimizer_attempts": ts_attempts,
                    "endpoint_recovery": recovery,
                },
            )
        ts_path = _write_system(systems, ts_name, directory / "ts.xyz")
        systems, hessian_ok = readuct.run_hessian_task(
            systems,
            [ts_name],
            output=["ts_hessian"],
            stop_on_error=False,
        )
        hessian_name = "ts_hessian" if "ts_hessian" in systems else "ts"
        imaginary_count: int | None = None
        if hessian_ok and systems[hessian_name].has_results():
            results = systems[hessian_name].get_results()
            hessian = getattr(results, "hessian", None)
            if hessian is not None:
                eigenvalues = np.linalg.eigvalsh(
                    np.asarray(hessian, dtype=float)
                )
                imaginary_count = int(
                    np.sum(eigenvalues < float(cfg.get("hessian_negative_threshold", -1.0e-6)))
                )
        first_order = bool(ts_ok and hessian_ok and imaginary_count == 1)
        if not first_order:
            product_path, product_ok, recovery = (
                self._relax_nt2_driven_endpoint(
                    collector, trial, cfg, directory
                )
            )
            if product_ok and product_path is not None:
                return DiscoveryResult(
                    driver="nt2",
                    success=True,
                    product_xyz_path=product_path,
                    ts_xyz_path=ts_path,
                    trajectory_path=trajectory,
                    electronic_energy_hartree=recovery.get(
                        "electronic_energy_hartree"
                    ),
                    low_level_ts_validated=False,
                    low_level_irc_connected=False,
                    imaginary_mode_count=imaginary_count,
                    biased_energy_used_as_barrier=False,
                    paths={
                        "product": product_path,
                        "nt2_guess": guess_path,
                        "ts_rejected": ts_path,
                        "driven_endpoint_seed": str(
                            recovery["driven_endpoint_seed"]
                        ),
                        **({"trajectory": trajectory} if trajectory else {}),
                    },
                    data={
                        "nt2_converged": nt_ok,
                        "ts_converged": ts_ok,
                        "ts_optimizer_attempts": ts_attempts,
                        "hessian_converged": hessian_ok,
                        "endpoint_recovered_without_ts_claim": True,
                        **recovery,
                    },
                )
            return DiscoveryResult(
                driver="nt2",
                success=False,
                ts_xyz_path=ts_path,
                trajectory_path=trajectory,
                low_level_ts_validated=False,
                imaginary_mode_count=imaginary_count,
                failure_reason="nt2_guess_is_not_a_frequency_validated_first_order_saddle",
                paths={
                    "nt2_guess": guess_path,
                    "ts": ts_path,
                    **({"trajectory": trajectory} if trajectory else {}),
                },
                data={
                    "nt2_converged": nt_ok,
                    "ts_converged": ts_ok,
                    "ts_optimizer_attempts": ts_attempts,
                    "hessian_converged": hessian_ok,
                    "endpoint_recovery": recovery,
                },
            )
        systems, irc_ok = readuct.run_irc_task(
            systems,
            [hessian_name],
            output=["irc_forward", "irc_backward"],
            convergence_max_iterations=int(cfg.get("irc_max_iterations", 150)),
            irc_initial_step_size=float(cfg.get("irc_initial_step_size", 0.3)),
            stop_on_error=False,
        )
        if not {"irc_forward", "irc_backward"}.issubset(systems):
            return DiscoveryResult(
                driver="nt2",
                success=False,
                ts_xyz_path=ts_path,
                trajectory_path=trajectory,
                low_level_ts_validated=True,
                imaginary_mode_count=imaginary_count,
                failure_reason="low_level_irc_endpoints_missing",
                data={"irc_converged": irc_ok},
            )
        systems, forward_ok = readuct.run_opt_task(
            systems,
            ["irc_forward"],
            output=["endpoint_forward"],
            convergence_max_iterations=int(cfg.get("endpoint_max_iterations", 200)),
            stop_on_error=False,
        )
        systems, backward_ok = readuct.run_opt_task(
            systems,
            ["irc_backward"],
            output=["endpoint_backward"],
            convergence_max_iterations=int(cfg.get("endpoint_max_iterations", 200)),
            stop_on_error=False,
        )
        if not {"endpoint_forward", "endpoint_backward"}.issubset(systems):
            return DiscoveryResult(
                driver="nt2",
                success=False,
                ts_xyz_path=ts_path,
                trajectory_path=trajectory,
                low_level_ts_validated=True,
                imaginary_mode_count=imaginary_count,
                failure_reason="unbiased_low_level_endpoint_optimization_failed",
                data={"forward_converged": forward_ok, "backward_converged": backward_ok},
            )
        forward_path = _write_system(systems, "endpoint_forward", directory / "endpoint_forward.xyz")
        backward_path = _write_system(systems, "endpoint_backward", directory / "endpoint_backward.xyz")
        return DiscoveryResult(
            driver="nt2",
            success=bool(forward_ok and backward_ok),
            product_xyz_path=forward_path,
            reactant_endpoint_xyz_path=backward_path,
            ts_xyz_path=ts_path,
            trajectory_path=trajectory,
            electronic_energy_hartree=_calculator_energy(systems["endpoint_forward"]),
            low_level_ts_validated=True,
            low_level_irc_connected=bool(irc_ok),
            imaginary_mode_count=imaginary_count,
            biased_energy_used_as_barrier=False,
            paths={
                "product": forward_path,
                "reactant_endpoint": backward_path,
                "nt2_guess": guess_path,
                "ts": ts_path,
                **({"trajectory": trajectory} if trajectory else {}),
            },
            data={
                "capability": capability,
                "nt2_converged": nt_ok,
                "ts_converged": ts_ok,
                "ts_optimizer_attempts": ts_attempts,
                "hessian_converged": hessian_ok,
                "irc_converged": irc_ok,
                "forward_converged": forward_ok,
                "backward_converged": backward_ok,
            },
        )

    def _run_afir(
        self,
        trial: ReactionTrialRecord,
        source: Artifact,
        cfg: dict[str, Any],
        directory: Path,
        capability: dict[str, Any],
    ) -> DiscoveryResult:
        import scine_readuct as readuct

        pair, attractive, bias_evidence = _afir_bias_target(trial)
        if pair is None or attractive is None:
            return DiscoveryResult(
                driver="afir",
                success=False,
                failure_reason=str(bias_evidence["reason"]),
                data={"afir_bias_evidence": bias_evidence},
            )
        systems = {"start": self._load_system(source, cfg)}
        collector = _TrajectoryCollector()
        systems, afir_ok = readuct.run_afir_task(
            systems,
            ["start"],
            False,
            [collector.collect],
            output=["afir_biased"],
            afir_lhs_list=[pair.atoms[0]],
            afir_rhs_list=[pair.atoms[1]],
            afir_attractive=attractive,
            afir_energy_allowance=float(cfg.get("afir_energy_allowance_kj_mol", 300.0)),
            convergence_max_iterations=int(cfg.get("afir_max_iterations", 200)),
            stop_on_error=False,
        )
        trajectory = _write_collector(collector, directory / "afir_biased_trajectory.xyz")
        if "afir_biased" not in systems:
            return DiscoveryResult(
                driver="afir",
                success=False,
                trajectory_path=trajectory,
                failure_reason="afir_did_not_produce_a_biased_endpoint",
                data={
                    "afir_converged": afir_ok,
                    "capability": capability,
                    "afir_bias_evidence": bias_evidence,
                },
            )
        biased_path = _write_system(systems, "afir_biased", directory / "afir_biased.xyz")
        systems, opt_ok = readuct.run_opt_task(
            systems,
            ["afir_biased"],
            output=["product_unbiased"],
            convergence_max_iterations=int(cfg.get("endpoint_max_iterations", 200)),
            stop_on_error=False,
        )
        if "product_unbiased" not in systems:
            return DiscoveryResult(
                driver="afir",
                success=False,
                trajectory_path=trajectory,
                failure_reason="afir_unbiased_endpoint_optimization_failed",
                paths={"biased_endpoint": biased_path},
                data={
                    "afir_converged": afir_ok,
                    "unbiased_opt_converged": opt_ok,
                    "afir_bias_evidence": bias_evidence,
                },
            )
        product_path = _write_system(systems, "product_unbiased", directory / "product_unbiased.xyz")
        return DiscoveryResult(
            driver="afir",
            success=bool(opt_ok),
            product_xyz_path=product_path,
            trajectory_path=trajectory,
            electronic_energy_hartree=_calculator_energy(systems["product_unbiased"]),
            biased_energy_used_as_barrier=False,
            paths={
                "product": product_path,
                "biased_endpoint": biased_path,
                **({"trajectory": trajectory} if trajectory else {}),
            },
            data={
                "capability": capability,
                "afir_converged": afir_ok,
                "unbiased_opt_converged": opt_ok,
                "afir_energy_is_biased_and_excluded_from_barriers": True,
                "afir_bias_evidence": bias_evidence,
            },
        )
