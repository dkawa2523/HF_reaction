"""placement.seeds and is_small_rigid (§8.2): rigid fragments at Σ r_vdW + 0.5 Å, pairwise
distinct, reproducible and independent of the input orientation."""

import numpy as np
import pytest
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from hfauto.chemistry import placement
from hfauto.chemistry.elements import vdw_radius
from hfauto.chemistry.identity import permutation_invariant_rmsd as pirmsd
from hfauto.chemistry.xyz import XYZ


def _xyz(symbols: str, *rows) -> XYZ:
    return XYZ(symbols.split(), np.array(rows, dtype=float))


NH3 = _xyz("N H H H", [0, 0, .1], [.94, 0, -.25], [-.47, .814, -.25], [-.47, -.814, -.25])
HF, AR = _xyz("H F", [0, 0, 0], [0, 0, .92]), _xyz("Ar", [0, 0, 0])
CH4 = _xyz("C H H H H", [0, 0, 0], [.63, .63, .63], [-.63, -.63, .63], [-.63, .63, -.63],
           [.63, -.63, -.63])
MEOH = _xyz("C O H H H H", [0, 0, 0], [1.43, 0, 0], [1.75, .9, 0], [-.36, 1.03, 0],
            [-.36, -.51, .89], [-.36, -.51, -.89])
HONO = _xyz("H O N O", [.95, .3, 0], [0, 0, 0], [-.5, 1.3, 0], [-1.7, 1.3, 0])
FECL3 = _xyz("Fe Cl Cl Cl", [0, 0, 0], [2.13, 0, 0], [-1.065, 1.845, 0], [-1.065, -1.845, 0])
CASES = [(NH3, [HF]), (CH4, [AR]), (NH3, [HF, HF]), (FECL3, [CH4]), (MEOH, [NH3])]


def _gap(a: XYZ, b: np.ndarray, b_symbols) -> float:
    radii = np.add.outer([vdw_radius(s) for s in a.symbols], [vdw_radius(s) for s in b_symbols])
    return float((cdist(a.coords, b) - radii).min())


@pytest.mark.parametrize("host,guests", CASES)
def test_seeds_are_rigid_at_contact_and_distinct(host, guests):
    seeds = placement.seeds(host, guests)
    assert len(seeds) == 6
    for seed in seeds:
        start = len(host.symbols)
        np.testing.assert_allclose(cdist(seed.coords[:start], seed.coords[:start]),
                                   cdist(host.coords, host.coords), atol=1e-8)
        for guest in guests:  # rigid, at contact with every atom placed before it
            end = start + len(guest.symbols)
            moved = seed.coords[start:end]
            np.testing.assert_allclose(cdist(moved, moved), cdist(guest.coords, guest.coords),
                                       atol=1e-8)
            placed = XYZ(seed.symbols[:start], seed.coords[:start])
            assert _gap(placed, moved, guest.symbols) == pytest.approx(0.5, abs=1e-6)
            start = end
    assert min(pirmsd(a.symbols, a.coords, b.coords)[0]
               for i, a in enumerate(seeds) for b in seeds[i + 1:]) >= 0.1


@pytest.mark.parametrize("host,guests", CASES)
def test_seeds_are_reproducible_and_do_not_depend_on_the_input_orientation(host, guests):
    seeds = placement.seeds(host, guests)
    assert all(np.array_equal(a.coords, b.coords)
               for a, b in zip(seeds, placement.seeds(host, guests), strict=True))
    q = [Rotation.random(random_state=k).as_matrix() for k in range(len(guests) + 1)]
    turned = [XYZ(m.symbols, m.coords @ r.T) for m, r in zip((host, *guests), q, strict=True)]
    again = placement.seeds(turned[0], turned[1:])
    for a, b in zip(seeds, again, strict=True):  # the whole complex turns with the host
        np.testing.assert_allclose(b.coords, a.coords @ q[0].T, atol=1e-6)


@pytest.mark.parametrize("mol,rigid", [(NH3, True), (HF, True), (AR, True), (HONO, False),
                                       (MEOH, False)])  # H–O–N=O and C–O rotate
def test_small_rigid_molecules_skip_the_search(mol, rigid):
    assert placement.is_small_rigid(mol.symbols, mol.coords) == rigid
