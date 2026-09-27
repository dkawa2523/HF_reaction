"""GoodVibes 4.3.0 worker (design §6.3): ``goodvibes.api.compute_thermo(qcdata=...)``.

Run as ``python -m hfauto.execution.worker hfauto.backends.goodvibes.worker:compute job.json``;
goodvibes is imported only inside ``compute``, so only the worker process loads it.

QCData fields that ``calc_bbe`` (goodvibes/thermo.py of 4.3.0, checked in the WSL install)
actually reads for the numbers:

- ``scf_energy``: E in H and G; ``multiplicity``: S_elec; ``molecular_mass`` (amu): S_trans.
- ``frequency_wn`` / ``im_frequency_wn``: only ``frequency_wn`` enters ZPE, U_vib and S_vib.
  The engine sends ``chemistry.thermo.thermo_frequencies`` (no negative mode) and ``invert``
  stays None.
- ``rotemp`` (K): S_rot and the calculation gate (empty -> no thermochemistry);
  ``linear_mol`` selects the linear formula, which uses ``rotemp[0]``; ``symmno``: sigma in
  S_rot (kept 1: ``symm=True`` adds -R ln sigma from pymsym itself). An atom's rotemp is
  (0, 0, 0): not empty, and with no frequency every rotational and vibrational term is 0, so
  G holds the translational, electronic and pV terms only.
- ``zero_point_corr``: only a gate (None -> no thermochemistry) and the monatomic test
  (== 0.0); its value is never used.
- ``atom_nums`` and ``cartesians``: pymsym's point group when ``symm=True``.
- ``roconst`` (GHz): only with ``inertia != "global"`` (not used here).
- ``file``: must not exist on disk, otherwise calc_bbe re-parses it (``parse_data``); it is
  also the label of the symmetry correction. ``job_type`` matters only for
  ``invert="auto"`` (not used). ``charge``, ``atom_types``, ``program`` are carried along.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.elements import atomic_number, mass
from hfauto.chemistry.gates import zpe_hartree
from hfauto.chemistry.vibrations import rotational_constants_ghz

_QCDATA_FILE = "hfauto-qcdata-without-file"  # never created: see the module docstring
_K_PER_GHZ = 0.0479924307  # h / k_B in K per GHz (CODATA 2018)


def qcdata_fields(job: dict[str, Any]) -> dict[str, Any]:
    """Keyword arguments of ``goodvibes.io.QCData`` for the freq job described by ``job``.

    ``job`` holds symbols, coords (Å, input frame), energy_hartree, charge, multiplicity,
    frequencies_cm1 (projected; imaginary < 0 if any) and n_external of one freq Evidence.
    """
    symbols: Sequence[str] = job["symbols"]
    coords = np.asarray(job["coords"], dtype=float).reshape(-1, 3)
    freqs = [float(nu) for nu in job["frequencies_cm1"]]
    linear = job["n_external"] == 5
    constants = rotational_constants_ghz(symbols, coords)
    if linear:  # (0, B, B): GoodVibes' linear formula reads the first, non-zero value
        constants = tuple(b for b in constants if b > 0.0)
    return {
        "file": _QCDATA_FILE,
        "program": "hfauto",
        "scf_energy": float(job["energy_hartree"]),
        "charge": int(job["charge"]),
        "multiplicity": int(job["multiplicity"]),
        "atom_types": list(symbols),
        "atom_nums": [atomic_number(s) for s in symbols],
        "cartesians": coords.tolist(),
        "frequency_wn": [nu for nu in freqs if nu > 0.0],
        "im_frequency_wn": [nu for nu in freqs if nu < 0.0],
        "linear_mol": linear,
        "molecular_mass": sum(mass(s) for s in symbols),
        "symmno": 1,
        "roconst": list(constants),
        "rotemp": [b * _K_PER_GHZ for b in constants],
        "zero_point_corr": zpe_hartree(freqs),
    }


def _error(kind: str, reason: str) -> dict[str, Any]:
    return {"error": {"kind": kind, "reason": reason}}


def compute(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """One result per (settings, temperature); vib_scale scales the frequencies and the ZPE."""
    import goodvibes
    from goodvibes.api import compute_thermo
    from goodvibes.constants import J_TO_AU
    from goodvibes.io import QCData
    from goodvibes.thermo import calc_rotational_entropy

    if goodvibes.__version__ != job["version_pin"]:
        found, pin = goodvibes.__version__, job["version_pin"]
        return _error("method_mismatch", f"goodvibes {found} != version pin {pin}")
    fields = qcdata_fields(job)
    atom = len(fields["atom_nums"]) == 1
    results = []
    for s in job["settings"]:
        for T in job["temperatures_K"]:
            s_rot = calc_rotational_entropy(T, fields["rotemp"], symmno=1, monatomic=atom,
                                            linear=fields["linear_mol"]) / J_TO_AU
            if not (atom or s_rot > 0.0):
                return _error("input_invalid", f"nonpositive_rotational_entropy:{s_rot}")
            r = compute_thermo(
                qcdata=QCData(**fields), QS=s["qs"], s_freq_cutoff=s["cutoff_cm1"],
                temperature=T, freq_scale_factor=s["vib_scale"],
                zpe_scale_factor=s["vib_scale"], invert=None, symm=s["symmetry"],
            )
            results.append({
                "settings_sha": s["sha"], "T_K": T, "G_hartree": r.qh_gibbs_free_energy,
                "H_hartree": r.enthalpy, "E_hartree": r.scf_energy, "zpe_hartree": r.zpe,
                "n_real": len(r.frequency_wn or ()), "S_rot": s_rot,
            })
    return {"results": results}
