"""The element table: supported range, Z and masses against RDKit, today's values kept, the
maximum coordination of reaction trials."""

import pytest

from hfauto.chemistry import elements as el


def test_supported_range_is_z_1_to_57_and_72_to_86():
    assert {e.z for e in el.ELEMENTS.values()} == {*range(1, 58), *range(72, 87)}
    assert all(symbol == e.symbol for symbol, e in el.ELEMENTS.items())
    for unsupported in ("Ce", "Lu", "Fr", "Xx", "cl"):
        with pytest.raises(KeyError):
            el.atomic_number(unsupported)


def test_z_and_masses_agree_with_rdkit():
    chem = pytest.importorskip("rdkit.Chem")
    table = chem.GetPeriodicTable()
    for symbol, e in el.ELEMENTS.items():
        assert table.GetAtomicNumber(symbol) == e.z
        assert e.mass == pytest.approx(table.GetMostCommonIsotopeMass(e.z), abs=1e-5)
        assert 0.2 < e.covalent_radius < e.vdw_radius < 3.5


def test_values_of_h_to_kr_and_i_are_unchanged():
    assert (el.mass("H"), el.mass("C"), el.mass("Kr"), el.mass("I")) == (
        1.007825, 12.0, 83.911498, 126.904472)
    assert (el.covalent_radius("V"), el.covalent_radius("Fe"), el.covalent_radius("I")) == (
        1.53, 1.32, 1.39)
    assert (el.vdw_radius("C"), el.vdw_radius("N"), el.vdw_radius("Ni")) == (1.70, 1.55, 1.63)
    assert (el.atomic_number("Te"), el.vdw_radius("Fe"), el.vdw_radius("Rn")) == (52, 2.44, 2.20)


def test_max_coordination_rule():
    for n, symbols in ((1, "H F Cl I At"), (3, "O"), (4, "B C N"), (0, "He Ne"),
                       (6, "Al Si P S Se Te Xe Pb Bi"), (9, "Li Na Mg Fe Pt La Hg")):
        assert {el.max_coordination(s) for s in symbols.split()} == {n}
