"""Electronic-state checks: electron-count parity and spin coupling of components."""

from __future__ import annotations

from collections.abc import Iterable

_ELEMENTS = """
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La
Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po
At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg
Cn Nh Fl Mc Lv Ts Og
""".split()  # noqa: SIM905 - compact periodic-table data is clearer as grouped text
_ATOMIC_NUMBERS = {symbol.upper(): number for number, symbol in enumerate(_ELEMENTS, 1)}


def check_electronic_state(symbols: Iterable[str], charge: int, multiplicity: int) -> None:
    """Raise ValueError for a non-positive multiplicity or one the electron count cannot have.

    Elements outside the table skip the parity check rather than assume an atomic number.
    """

    if multiplicity < 1:
        raise ValueError(f"Multiplicity must be positive, got {multiplicity}")
    normalized = [symbol.strip().upper() for symbol in symbols]
    if not all(symbol in _ATOMIC_NUMBERS for symbol in normalized):
        return
    electron_count = sum(_ATOMIC_NUMBERS[symbol] for symbol in normalized) - charge
    unpaired = multiplicity - 1
    if electron_count < unpaired or (electron_count - unpaired) % 2:
        raise ValueError(
            "Electron-count parity is incompatible with multiplicity: "
            f"electrons={electron_count}, multiplicity={multiplicity}"
        )


def coupled_multiplicities(component_multiplicities: Iterable[int]) -> tuple[int, ...]:
    """Return every spin multiplicity allowed by angular-momentum coupling.

    Multiplicity is represented as ``2S + 1``.  Working with ``2S`` keeps the
    calculation exact for both integer and half-integer fragment spins.
    """

    values = list(component_multiplicities)
    if not values or any(value < 1 for value in values):
        raise ValueError("component multiplicities must be positive")
    coupled_twice_spins = {values[0] - 1}
    for multiplicity in values[1:]:
        fragment_twice_spin = multiplicity - 1
        coupled_twice_spins = {
            total_twice_spin
            for current in coupled_twice_spins
            for total_twice_spin in range(
                abs(current - fragment_twice_spin),
                current + fragment_twice_spin + 1,
                2,
            )
        }
    return tuple(sorted(value + 1 for value in coupled_twice_spins))
