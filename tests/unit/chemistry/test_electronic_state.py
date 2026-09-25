import pytest

from hfauto.chemistry.electronic_state import check_electronic_state, coupled_multiplicities


def test_electron_parity_and_spin_coupling():
    check_electronic_state(["N", "H", "H", "H"], 0, 1)
    check_electronic_state(["Xx", "H"], 0, 2)  # unknown element: no parity check
    for charge, multiplicity in ((0, 2), (1, 1), (0, 0)):
        with pytest.raises(ValueError):
            check_electronic_state(["N", "H", "H", "H"], charge, multiplicity)
    assert coupled_multiplicities([2, 2]) == (1, 3)
    assert coupled_multiplicities([3, 2]) == (2, 4)
    with pytest.raises(ValueError):
        coupled_multiplicities([])
