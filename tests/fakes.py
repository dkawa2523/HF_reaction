"""Test doubles for the capability Protocols (design §6.4); the API is frozen after Wave 3.
Files go below ``root`` (the run directory); FileRefs are relative to it. Behaviour one test
needs (a collapsing saddle, a soft mode...) belongs in a subclass inside that test."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import minimize
from scipy.spatial.distance import pdist

from hfauto.backends import protocols as bp
from hfauto.chemistry.profile import hei
from hfauto.chemistry.thermo import Thermal
from hfauto.chemistry.vibrations import projected_frequencies, shape_hessian, to_canonical_npy
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    read_xyz,
    read_xyz_trajectory,
    write_xyz,
    write_xyz_trajectory,
    written_geometry,
)
from hfauto.core.constants import BOHR_TO_ANGSTROM as BOHR
from hfauto.core.constants import CM1_TO_HARTREE
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry, Level, PathProfile
from hfauto.core.evidence import FailureKind as Kind
from hfauto.core.hashing import sha256_file
from hfauto.core.method import MethodSpec

GMAX_TOL = 1.0e-5  # Eh/bohr: largest gradient component at convergence
_FD = 1.0e-4  # central-difference step (Å)
_R_AB, _K_AB, _K_BEND, _HALF = 2.8, 0.5, 0.05, 0.4  # line triatomic A–H–B (Å, Eh/Å²)


def _fd(f: Callable[[np.ndarray], Any], x: np.ndarray) -> np.ndarray:
    return np.array([f(x + d) - f(x - d) for d in np.eye(x.size) * _FD]) / (2 * _FD)


@dataclass(frozen=True)
class PES:
    """energy(coords Å) -> Eh; gradient (flat, Eh/Å) and Hessian (Eh/bohr²) by differences."""

    symbols: tuple[str, ...]
    energy: Callable[[np.ndarray], float]
    points: Mapping[str, np.ndarray]  # named stationary points and starting structures

    def gradient(self, coords: np.ndarray) -> np.ndarray:
        return _fd(self.energy, np.asarray(coords, dtype=float).ravel())

    def hessian(self, coords: np.ndarray) -> np.ndarray:
        h = _fd(self.gradient, np.asarray(coords, dtype=float).ravel())
        return 0.5 * (h + h.T) * BOHR**2

    def molecule(self, point: str, charge: int = 0, multiplicity: int = 1) -> Molecule:
        return Molecule(XYZ(list(self.symbols), np.array(self.points[point], dtype=float)),
                        charge, multiplicity)


def harmonic() -> PES:
    """Water-like triangle held by pairwise springs: a single minimum."""
    ref = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])
    d0 = pdist(ref)
    return PES(("O", "H", "H"), lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - d0) ** 2)),
               {"minimum": ref, "start": ref + [[0, 0, 0], [0.05, 0.02, 0], [-0.03, 0, 0.01]]})


def _wells(symbols: tuple[str, ...], crit: Sequence[float], height: float = 0.01,
           ref: tuple[int, int] = (0, 1), names: Sequence[str] = ("reactant", "ts", "product"),
           r_ab: float = _R_AB) -> PES:
    """Collinear A–H–B, A···B r_ab apart; V'(x) ∝ Π(x − crit) with x the H position along A–B
    over _HALF, so minima and maxima alternate; V(crit[ref[1]]) − V(crit[ref[0]]) = height."""
    v = Polynomial.fromroots(crit).integ()
    v = v * (height / (v(crit[ref[1]]) - v(crit[ref[0]])))

    def energy(x: np.ndarray) -> float:
        a, h, b = np.asarray(x, dtype=float).reshape(3, 3)
        r, rel = float(np.linalg.norm(b - a)), h - 0.5 * (a + b)
        s = float(rel @ (b - a)) / r
        return float(v(s / _HALF) + 0.5 * _K_AB * (r - r_ab) ** 2
                     + 0.5 * _K_BEND * (rel @ rel - s * s))

    points = {n: np.array([[-r_ab / 2, 0, 0], [c * _HALF, 0, 0], [r_ab / 2, 0, 0]])
              for n, c in zip(names, crit, strict=True)}
    return PES(symbols, energy, points)


# ``height``: the barrier above the reactant (flat_uphill: the rise). double_well's reactant
# (H at N) lies below its product; symmetric_double_well's minima are permutation images;
# triple_well has an intermediate between "ts1" and "ts2" (H bonded to both N and O: a chemical
# state of its own, as a well that splits a bond-changing case must be); flat_uphill's product sits
# ~3e-6 Eh below a shoulder (barrierless at any practical resolution).
double_well = partial(_wells, ("N", "H", "O"), (-1, 0.1, 1))
symmetric_double_well = partial(_wells, ("F", "H", "F"), (-1, 0, 1))
triple_well = partial(_wells, ("N", "H", "O"), (-1, -0.5, 0, 0.5, 1),
                      names=("reactant", "ts1", "intermediate", "ts2", "product"), r_ab=2.6)
flat_uphill = partial(_wells, ("N", "H", "O"), (-1, 0.9, 1), ref=(0, 2))

# NH3·HF with H3 pointing at F: the double H exchange (amine_pilot2, −1265i) is the source
# relabelled, H3 ↔ H4 (formed (0,4),(3,5), broken (0,3),(4,5)), one state label and one basin.
NH3_HF_SYMBOLS = ("N", "H", "H", "H", "H", "F")
NH3_HF = np.array([[0.0, 0.0, 0.0], [-0.37, 0.47, 0.83], [-0.37, 0.47, -0.83], [0.55, 0.85, 0.0],
                   [1.694, 0.391, 0.0], [2.6, 0.6, 0.0]])
NH3_HF_EXCHANGED = NH3_HF[[0, 1, 2, 4, 3, 5]]


def _ref(root: Path, path: Path) -> FileRef:
    return FileRef(path=path.relative_to(root).as_posix(), sha256=sha256_file(path))


def write_hessian(root: Path, path: str | Path, hessian: np.ndarray) -> FileRef:
    """Write a canonical Hessian .npy (``path`` relative to root) and return its FileRef."""
    return _ref(Path(root), to_canonical_npy(hessian, Path(root) / path))


def write_geometry(root: Path, path: str | Path, symbols: Sequence[str], coords: np.ndarray
                   ) -> Geometry:
    """Write an xyz (``path`` relative to root); the fingerprint is of the file as written."""
    target = write_xyz(XYZ(list(symbols), np.asarray(coords, dtype=float).reshape(-1, 3)),
                       Path(root) / path)
    return written_geometry(target, partial(_ref, Path(root)))


def xyz_loader(root: Path) -> Callable[[Geometry], XYZ]:
    return lambda geometry: read_xyz(Path(root) / geometry.file.path)


def fake_level(method: MethodSpec, mol: Molecule) -> Level:
    same = method.model_dump(include={"basis", "dispersion", "grid",
                                      "electronic_temperature_K"})
    return Level(program="fake", version="0", charge=mol.charge, multiplicity=mol.multiplicity,
                 method=method.functional or method.wft_method or f"gfn{method.gfn}",
                 scf_tol=method.scf_energy_tol, **same)


def _cap(step: np.ndarray, limit: float) -> np.ndarray:
    largest = float(np.linalg.norm(step.reshape(-1, 3), axis=1).max())
    return step * min(1.0, limit / largest) if largest > 0 else step


def _minimize(pes: PES, x0: np.ndarray, h0: np.ndarray | None, maxiter: int = 500
              ) -> tuple[np.ndarray, bool]:
    """BFGS with Armijo backtracking; h0 (Eh/bohr²) seeds the inverse Hessian."""
    x, inv = np.asarray(x0, dtype=float).ravel(), np.eye(np.size(x0)) * 2.0
    if h0 is not None:
        w, vecs = np.linalg.eigh(h0 / BOHR**2)
        inv = vecs @ np.diag(1.0 / np.maximum(np.abs(w), 0.05)) @ vecs.T
    g, e = pes.gradient(x), pes.energy(x)
    for _ in range(maxiter):
        if np.abs(g).max() * BOHR < GMAX_TOL:
            return x.reshape(-1, 3), True
        step, t = _cap(-inv @ g, 0.2), 1.0
        while pes.energy(x + t * step) > e + 1e-4 * t * (g @ step) and t > 1e-6:
            t *= 0.5
        x, s = x + t * step, t * step
        g_old, g, e = g, pes.gradient(x), pes.energy(x)
        y = g - g_old
        if y @ s > 1e-12:
            m = np.eye(x.size) - np.outer(s, y) / (y @ s)
            inv = m @ inv @ m.T + np.outer(s, s) / (y @ s)
    return x.reshape(-1, 3), False


def _minimize_fixed(pes: PES, x0: np.ndarray, bond: tuple[int, int, float]
                    ) -> tuple[np.ndarray, bool]:
    """SLSQP with atoms i and j held r Å apart (NWChem's frozen zcoord bond)."""
    i, j, r = bond

    def stretch(x: np.ndarray) -> float:
        return float(np.linalg.norm(x[3 * j:3 * j + 3] - x[3 * i:3 * i + 3]) - r)

    res = minimize(pes.energy, np.ravel(x0), jac=pes.gradient, method="SLSQP",
                   constraints=[{"type": "eq", "fun": stretch}],
                   options={"ftol": 1e-12, "maxiter": 500})
    return res.x.reshape(-1, 3), bool(res.success)


def _saddle(pes: PES, x0: np.ndarray, mode: np.ndarray | None, maxiter: int = 200
            ) -> tuple[np.ndarray, bool]:
    """Eigenvector following: climb along the tracked mode (initially ``mode``, else the
    lowest one) and descend along all others."""
    x, v = np.asarray(x0, dtype=float).ravel(), mode
    for _ in range(maxiter):
        g = pes.gradient(x)
        if np.abs(g).max() * BOHR < GMAX_TOL:
            return x.reshape(-1, 3), True
        w, vecs = np.linalg.eigh(pes.hessian(x) / BOHR**2)
        w, vecs = w[np.abs(w) > 1e-6], vecs[:, np.abs(w) > 1e-6]  # drop external modes
        i = 0 if v is None else int(np.argmax(np.abs(vecs.T @ np.ravel(v))))
        v, dq = vecs[:, i], -(vecs.T @ g) / np.abs(w)
        dq[i] = -dq[i]
        x = x + _cap(vecs @ dq, 0.1)
    return x.reshape(-1, 3), False


class _Fake:
    name: ClassVar[str] = "fake"

    def __init__(self) -> None:
        self.calls: list[Any] = []

    @classmethod
    def requirements(cls) -> bp.Requirements:
        return bp.Requirements()

    def supports(self, method: MethodSpec) -> bool:
        return True


class _Surface(_Fake):
    def __init__(self, root: Path, pes: PES) -> None:
        super().__init__()
        self.root, self.pes = Path(root), pes

    def _key(self, *parts: object) -> str:
        text = json.dumps([type(self).__name__, *parts], default=str)
        return hashlib.sha256(text.encode()).hexdigest()

    def _start(self, mol: Molecule, key: str) -> Geometry:
        return write_geometry(self.root, f"fake/{key[:16]}/start.xyz", mol.xyz.symbols,
                              mol.xyz.coords)

    def _hessian(self, ev: Evidence, start: Geometry, key: str, near_A: float = 0.0
                 ) -> np.ndarray | Failure:
        """The freq Hessian at ``start``, or at most ``near_A`` per atom from it (same frame)."""
        at, here = (read_xyz(self.root / g.file.path) for g in (ev.final, start))
        same_atoms = list(at.symbols) == list(here.symbols)
        shift = np.linalg.norm(at.coords - here.coords, axis=1).max() if same_atoms else np.inf
        if ev.hessian is None or (ev.final.fingerprint != start.fingerprint and shift > near_A):
            return Failure(kind=Kind.INPUT_INVALID, reason="hessian_geometry_mismatch", job_key=key)
        return np.load(self.root / ev.hessian.path)

    def _evidence(self, task: Any, mol: Molecule, method: MethodSpec, key: str, start: Geometry,
                  final: np.ndarray | None = None, **extra: Any) -> Evidence:
        end = start if final is None else write_geometry(
            self.root, f"fake/{key[:16]}/final.xyz", mol.xyz.symbols, final)
        energy = self.pes.energy(mol.xyz.coords if final is None else final)
        return Evidence(engine=self.name, task=task, level=fake_level(method, mol), start=start,
                        final=end, energy_hartree=energy, output=end.file, job_key=key, **extra)


class FakeQM(_Surface):
    """calls: "energy" (+"+scf_guess"), "optimize" (+"+init_hessian", +"+fixed_bond"),
    "frequencies"."""

    def energy(self, mol, method, *, scf_guess: Evidence | None = None, deadline=None
               ) -> Evidence | Failure:
        """scf_guess only enters the key and the call, when given (a PES has no SCF branches)."""
        self.calls.append("energy+scf_guess" if scf_guess else "energy")
        guess = (scf_guess.job_key,) if scf_guess else ()
        key = self._key("sp", mol.fingerprint(), method.signature(), *guess)
        return self._evidence("sp", mol, method, key, self._start(mol, key))

    def optimize(self, mol, method, *, init_hessian: Evidence | None = None,
                 fixed_bond: tuple[int, int, float] | None = None,
                 scf_guess: Evidence | None = None, deadline=None) -> Evidence | Failure:
        """Like NWChem: init_hessian within 0.5 Å; fixed_bond (i, j, r) at mol within 1e-4 Å
        and then held exactly; scf_guess only enters the key (a PES has no SCF branches)."""
        self.calls.append("+".join(["optimize", *(["init_hessian"] if init_hessian else []),
                                    *(["fixed_bond"] if fixed_bond else [])]))
        key = self._key("opt", mol.fingerprint(), method.signature(),
                        init_hessian and init_hessian.job_key, fixed_bond,
                        scf_guess and scf_guess.job_key)
        start = self._start(mol, key)
        h0 = None if init_hessian is None else self._hessian(init_hessian, start, key, 0.5)
        if isinstance(h0, Failure):
            return h0
        if fixed_bond is None:
            x, ok = _minimize(self.pes, mol.xyz.coords, h0)
        else:
            i, j, r = fixed_bond
            if abs(np.linalg.norm(mol.xyz.coords[j] - mol.xyz.coords[i]) - r) > 1e-4:
                return Failure(kind=Kind.INPUT_INVALID, reason="fixed_bond_mismatch", job_key=key)
            x, ok = _minimize_fixed(self.pes, mol.xyz.coords, fixed_bond)
        if not ok:
            return Failure(kind=Kind.GEOMETRY_MAXITER, reason="maxiter", job_key=key)
        return self._evidence("opt", mol, method, key, start, x)

    def frequencies(self, mol, method, *, scf_guess=None, deadline=None) -> Evidence | Failure:
        self.calls.append("frequencies")
        key = self._key("freq", mol.fingerprint(), method.signature())
        start, h = self._start(mol, key), self.pes.hessian(mol.xyz.coords)
        freqs, modes, n_external = projected_frequencies(h, mol.xyz.symbols, mol.xyz.coords)
        npy = to_canonical_npy(h, self.root / f"fake/{key[:16]}/hessian.npy")
        imaginary = tuple(tuple(map(float, m)) for f, m in zip(freqs, modes, strict=True) if f < 0)
        return self._evidence("freq", mol, method, key, start, hessian=_ref(self.root, npy),
                              frequencies_cm1=tuple(map(float, freqs)), n_external=n_external,
                              imaginary_modes=imaginary)


class FakeSaddle(_Surface):  # calls: "refine"; like NWChemSaddle: a Hessian within 0.5 Å,
    def refine(self, seed, method, *, hessian: Evidence, mode, deadline=None  # the only
               ) -> Evidence | Failure:  # negative mode followed, the last frame at maxiter
        self.calls.append("refine")
        key = self._key(seed.fingerprint(), method.signature(), hessian.job_key,
                        np.round(np.ravel(mode), 6).tolist())
        start = self._start(seed, key)
        h = self._hessian(hessian, start, key, 0.5)
        if isinstance(h, Failure):
            return h
        _, vectors = np.linalg.eigh(shape_hessian(h, seed.xyz.coords, mode))
        x, ok = _saddle(self.pes, seed.xyz.coords, vectors[:, 0])
        if not ok:
            last = write_geometry(self.root, f"fake/{key[:16]}/last.xyz", seed.xyz.symbols, x)
            return Failure(kind=Kind.GEOMETRY_MAXITER, reason="maxiter", final=last, job_key=key)
        return self._evidence("saddle", seed, method, key, start, x)


class FakePath(_Surface):
    """The initial path as it stands (no relaxation) with its PES energies. Each call takes the
    next ``script`` entry: "pes" (also once empty), "barrierless" / "single" / "intermediate"
    (scripted energies) or "failed". With ``tsopt`` (the low-level NEB) a TS is optimized from
    the highest interior image, like pysis_neb's TSOpt from the climbing image."""

    def __init__(self, root: Path, pes: PES, script: Sequence[str] = (), *,
                 tsopt: bool = False) -> None:
        super().__init__(root, pes)
        self.script, self.tsopt = list(script), tsopt

    def find_path(self, start, end, method, *, images: int, initial_path: FileRef,
                  deadline=None) -> PathProfile | Failure:
        shape = self.script.pop(0) if self.script else "pes"
        self.calls.append(f"find_path:{shape}")
        key = self._key(start.fingerprint(), end.fingerprint(), method.signature(), images,
                        initial_path.sha256, shape)
        if shape == "failed":
            return Failure(kind=Kind.NONZERO_EXIT, reason="scripted", job_key=key)
        frames = [i.coords for i in read_xyz_trajectory(self.root / initial_path.path)]
        assert len(frames) == images, "the initial path holds the images"
        t, e = np.linspace(0.0, 1.0, images), np.array([self.pes.energy(f) for f in frames])
        bumps = {"barrierless": 0 * t, "single": np.sin(np.pi * t),
                 "intermediate": np.sin(2 * np.pi * t) ** 2}
        e = e[0] + (e[-1] - e[0]) * t + 0.01 * bumps[shape] if shape in bumps else e
        folder = self.root / "fake" / key[:16]
        xyz = write_xyz_trajectory([XYZ(list(start.xyz.symbols), f) for f in frames],
                                   folder / "images.xyz")
        x, ok = (_saddle(self.pes, hei(frames, e, 1 + int(np.argmax(e[1:-1])))[2], None)
                 if self.tsopt else (None, False))
        ts = write_geometry(self.root, folder / "ts.xyz", start.xyz.symbols, x) if ok else None
        return PathProfile(engine=self.name, level=fake_level(method, start), ts=ts,
                           images=_ref(self.root, xyz), energies_hartree=tuple(e.tolist()))


class _Scripted(_Fake):  # replays the results in order, or calls a callable script
    def __init__(self, script: Callable[..., Any] | Sequence[Any]) -> None:
        super().__init__()
        self._script = script if callable(script) else list(script)

    def _next(self, *args: Any) -> Any:
        self.calls.append(args)
        if callable(self._script):
            return self._script(*args)
        if not self._script:
            raise AssertionError(f"{type(self).__name__}: script exhausted")
        return self._script.pop(0)


class FakeConformers(_Scripted):
    def search(self, mol, method, settings) -> bp.ConformerEnsemble | Failure:
        return self._next(mol, method, settings)


class FakeDiscovery(_Scripted):
    def explore(self, source, trial, method, settings) -> bp.DiscoveryResult | Failure:
        return self._next(source, trial, method, settings)


def fake_species_thermo(symbols, sym, modes_cm1, *, multiplicity, settings, T) -> Thermal:
    """chemistry.thermo.species_thermo without GoodVibes (monkeypatch it in): H = ZPE and
    G = ZPE - 1e-6 Eh x cutoff x n_modes, so the cutoff variants spread the band."""
    zpe = 0.5 * settings.vib_scale * sum(modes_cm1) * CM1_TO_HARTREE
    return Thermal(G=zpe - 1e-6 * settings.cutoff_cm1 * len(modes_cm1), H=zpe, zpe=zpe)
