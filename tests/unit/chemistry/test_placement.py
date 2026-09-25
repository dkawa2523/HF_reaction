"""placement.seeds: H-bond cone, rigid fallback, collisions, determinism (§8.2, CH-19)."""

import numpy as np
import pytest
from scipy.spatial.distance import cdist

from hfauto.chemistry import placement
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


@pytest.mark.parametrize("host,guests", [(NH3, [HF]), (CH4, [AR]), (NH3, [HF, HF])])
def test_seeds_are_rigid_collision_free_and_distinct(host, guests):
    seeds, sizes = placement.seeds(host, guests), [len(m.symbols) for m in (host, *guests)]
    owner, bounds = np.repeat(np.arange(len(sizes)), sizes), np.cumsum([0, *sizes])
    heavy = np.array([s != "H" for s in seeds[0].symbols])
    inter, both = owner[:, None] != owner, np.outer(heavy, heavy)
    assert len(seeds) == 6  # CH4·Ar has no acceptor: rigid random fallback
    for seed in seeds:
        d = cdist(seed.coords, seed.coords)
        for m, s, e in zip((host, *guests), bounds[:-1], bounds[1:], strict=True):
            np.testing.assert_allclose(d[s:e, s:e], cdist(m.coords, m.coords), atol=1e-8)  # rigid
        assert d[inter & both].min() >= 2.2 and d[inter & ~both].min() >= 1.2
    assert min(pirmsd(a.symbols, a.coords, b.coords)[0]
               for i, a in enumerate(seeds) for b in seeds[i + 1:]) >= 0.1


def test_only_polar_hydrogens_donate():
    for seed in placement.seeds(MEOH, [NH3]):
        d = cdist(seed.coords[:6], seed.coords[6:])
        closest = tuple(int(i) for i in np.unravel_index(np.argmin(d), d.shape))
        assert closest in {(1, 1), (1, 2), (1, 3), (2, 0)}  # N–H···O or O–H···N, never C–H


def test_seeds_are_reproducible_and_depend_on_the_rng_seed():
    first = placement.seeds(NH3, [HF])
    assert all(np.array_equal(a.coords, b.coords)
               for a, b in zip(first, placement.seeds(NH3, [HF]), strict=True))
    assert not np.allclose(placement.seeds(NH3, [HF], rng_seed=1)[0].coords, first[0].coords)
