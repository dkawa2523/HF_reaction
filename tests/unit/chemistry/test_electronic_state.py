import pytest

from hfauto.chemistry.electronic_state import (
    check_electronic_state,
    composition_multiplicity,
    coupled_multiplicities,
    low_spin_coupled,
)


def test_elements_and_parity_of_the_declared_multiplicity():
    for symbols, charge, multiplicity in ((["N", "H", "H", "H"], 0, 1), (["Fe"], 2, 5),
                                          (["F", "C", "F"], 0, 1), (["F", "C", "F"], 0, 3)):
        check_electronic_state(symbols, charge, multiplicity)  # the declaration is trusted
    for symbols in (["Ce", "Cl", "Cl", "Cl"], ["Fr"], ["D", "H"]):
        with pytest.raises(ValueError, match=f"unsupported_element:{symbols[0]}"):
            check_electronic_state(symbols, 0, 2)
    for symbols, charge, multiplicity in ((["N", "H", "H", "H"], 0, 2),
                                          (["N", "H", "H", "H"], 1, 1),
                                          (["C", "H", "H", "H"], 0, 1), (["H"], 0, 4)):
        with pytest.raises(ValueError, match="parity"):
            check_electronic_state(symbols, charge, multiplicity)


def test_coupled_multiplicities():
    assert coupled_multiplicities([2, 2]) == (1, 3)
    assert coupled_multiplicities([3, 2]) == (2, 4)
    with pytest.raises(ValueError):
        coupled_multiplicities([])


@pytest.mark.parametrize(("components", "multiplicity", "expected"), [
    ([2, 3], 2, True),  # CH3. + O2 doublet
    ([2, 3], 4, False),  # its quartet is the high-spin coupling
    ([2, 2], 1, True),  # OH. + OH. singlet
    ([2, 2], 3, False),
    ([3, 3, 1], 3, True),  # O2 + O2 + N2 triplet
    ([2, 1], 2, False),  # one open-shell component
    ([1, 1], 1, False),  # H2O + H2O
])
def test_low_spin_coupled(components, multiplicity, expected):
    assert low_spin_coupled(components, multiplicity) is expected


@pytest.mark.parametrize(("components", "declared", "expected"), [
    ([1, 1], None, 1),  # HF + HF
    ([2, 1], None, 2),  # radical + closed shell
    ([2, 2], 3, 3),  # OH. + OH. declared triplet
    ([2, 3], 2, 2),  # CH3. + O2 declared doublet: low-spin coupled, accepted
    ([2, 3], 4, 4),  # CH3. + O2 declared quartet
    ([3, 3], 3, 3),  # O2 + O2 declared triplet: low-spin coupled, accepted
    ([6, 1], 6, 6),  # FeCl3 + CH4
])
def test_composition_multiplicity(components, declared, expected):
    assert composition_multiplicity(components, declared) == expected


@pytest.mark.parametrize(("components", "declared", "match"), [
    ([2, 2], None, r"declare_multiplicity: candidates \(1, 3\)"),
    ([2, 3], None, r"declare_multiplicity: candidates \(2, 4\)"),
    ([2, 2], 1, "low_spin_singlet_unsupported"),
    ([3, 3], 1, "low_spin_singlet_unsupported"),
    ([2, 2], 5, "not among"),
    ([2, 1], 1, "not among"),
])
def test_composition_multiplicity_rejects(components, declared, match):
    with pytest.raises(ValueError, match=match):
        composition_multiplicity(components, declared)
