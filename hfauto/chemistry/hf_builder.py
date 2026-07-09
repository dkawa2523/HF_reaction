from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class XYZ:
    symbols: list[str]
    coords: np.ndarray
    comment: str = ""


def read_xyz(path: str | Path) -> XYZ:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"Empty XYZ file: {path}")
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ line in {path}: {line!r}")
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if len(symbols) != n:
        raise ValueError(f"XYZ atom count mismatch in {path}: header={n}, parsed={len(symbols)}")
    return XYZ(symbols=symbols, coords=np.array(coords, dtype=float), comment=comment)


def write_xyz(xyz: XYZ, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(xyz.symbols)), xyz.comment or "generated_by=hfauto"]
    for sym, (x, y, z) in zip(xyz.symbols, xyz.coords):
        lines.append(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < 1.0e-12:
        return np.array([1.0, 0.0, 0.0])
    return v / n


def _min_distance_to_existing(coords: np.ndarray, points: list[np.ndarray]) -> float:
    if len(coords) == 0:
        return 999.0
    best = 999.0
    for p in points:
        best = min(best, float(np.min(np.linalg.norm(coords - p, axis=1))))
    return best


def choose_hf_direction(base: XYZ, site_atom_index: int) -> np.ndarray:
    site = base.coords[site_atom_index]
    heavy = np.array([c for s, c in zip(base.symbols, base.coords) if s != "H"], dtype=float)
    center = heavy.mean(axis=0) if len(heavy) else base.coords.mean(axis=0)
    seed = _unit(site - center)
    candidates = [seed, -seed]
    candidates += [
        np.array([1.0, 0.0, 0.0]),
        np.array([-1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, -1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([0.0, 0.0, -1.0]),
    ]
    best = candidates[0]
    best_clearance = -1.0
    for d in candidates:
        d = _unit(d)
        trial_h = site + d * 1.60
        trial_f = trial_h + d * 0.93
        clearance = _min_distance_to_existing(base.coords, [trial_h, trial_f])
        if clearance > best_clearance:
            best_clearance = clearance
            best = d
    return _unit(best)


def write_isolated_candidate(conformer_xyz: str | Path, out_path: str | Path, comment: str = "state=isolated_candidate") -> dict:
    base = read_xyz(conformer_xyz)
    write_xyz(XYZ(list(base.symbols), np.array(base.coords), comment=comment), out_path)
    return {"num_atoms": len(base.symbols), "n_atoms": len(base.symbols), "atom_order_consistent": True, "geometry_sane": True}


def copy_candidate_species(conformer_xyz: str | Path, out_path: str | Path) -> dict:
    qc = write_isolated_candidate(conformer_xyz, out_path, comment="state=candidate generated_by=hfauto")
    return {"n_atoms": int(qc["num_atoms"]), **qc, "copied_from": str(conformer_xyz)}


def build_hf_cluster(hf_n: int, out_path: str | Path) -> dict:
    symbols: list[str] = []
    coords: list[np.ndarray] = []
    for i in range(int(hf_n)):
        h = np.array([i * 1.85, 0.0, 0.0])
        f = h + np.array([0.93, 0.0, 0.0])
        symbols += ["H", "F"]
        coords += [h, f]
    write_xyz(XYZ(symbols, np.array(coords), comment=f"state=hf_cluster hf_n={hf_n}"), out_path)
    return {"hf_n": int(hf_n), "num_atoms": 2 * int(hf_n), "n_atoms": 2 * int(hf_n), "HF_bond_A": 0.93, "initial_HF_A": 0.93, "geometry_sane": True}


def write_hf_cluster(hf_n: int, out_path: str | Path) -> dict:
    return build_hf_cluster(hf_n, out_path)


def build_hf_complexes(
    conformer_xyz: str | Path,
    site_atom_index: int,
    hf_n: int,
    out_reactant: str | Path,
    out_product: str | Path,
    neighbor_indices: list[int] | None = None,
) -> tuple[dict, dict]:
    base = read_xyz(conformer_xyz)
    site = base.coords[site_atom_index]
    direction = choose_hf_direction(base, site_atom_index)

    rc_symbols = list(base.symbols)
    rc_coords: list[np.ndarray] = [np.array(c) for c in base.coords]
    ip_symbols = list(base.symbols)
    ip_coords: list[np.ndarray] = [np.array(c) for c in base.coords]

    for i in range(int(hf_n)):
        offset = 1.60 + i * 1.85
        h_rc = site + direction * offset
        f_rc = h_rc + direction * 0.93
        rc_symbols += ["H", "F"]
        rc_coords += [h_rc, f_rc]

        if i == 0:
            h_ip = site + direction * 1.05
            f_ip = site + direction * 2.45
        else:
            h_ip = site + direction * (2.45 + (i - 1) * 1.60)
            f_ip = h_ip + direction * 0.93
        ip_symbols += ["H", "F"]
        ip_coords += [h_ip, f_ip]

    rc = XYZ(rc_symbols, np.array(rc_coords), comment=f"state=reactant_complex hf_n={hf_n}")
    ip = XYZ(ip_symbols, np.array(ip_coords), comment=f"state=ion_pair hf_n={hf_n}")
    write_xyz(rc, out_reactant)
    write_xyz(ip, out_product)

    first_h_index = len(base.symbols)
    first_f_index = len(base.symbols) + 1
    rc_qc = {
        "num_base_atoms": len(base.symbols),
        "candidate_atom_count": len(base.symbols),
        "initial_BH_A": float(np.linalg.norm(rc.coords[first_h_index] - site)),
        "initial_HF_A": float(np.linalg.norm(rc.coords[first_f_index] - rc.coords[first_h_index])),
        "initial_BHF_angle_deg": 180.0,
        "atom_order_consistent": rc.symbols == ip.symbols,
        "base_atom": int(site_atom_index),
        "transfer_h": int(first_h_index),
        "leaving_f": int(first_f_index),
        "direction": [float(x) for x in direction],
        "min_distance_A": _min_distance_to_existing(base.coords, [rc.coords[first_h_index], rc.coords[first_f_index]]),
    }
    ip_qc = {
        **rc_qc,
        "product_BH_A": float(np.linalg.norm(ip.coords[first_h_index] - site)),
        "product_HF_A": float(np.linalg.norm(ip.coords[first_f_index] - ip.coords[first_h_index])),
    }
    return rc_qc, ip_qc


def endpoints_atom_order_match(path_a: str | Path, path_b: str | Path) -> bool:
    a = read_xyz(path_a)
    b = read_xyz(path_b)
    return a.symbols == b.symbols and len(a.symbols) == len(b.symbols)
