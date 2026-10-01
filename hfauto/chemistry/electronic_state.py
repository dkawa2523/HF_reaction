"""Electronic state of the input: supported elements, the spin declaration, electron-count
parity, and the multiplicity of a composition from its components."""

from __future__ import annotations

from collections.abc import Iterable

from hfauto.chemistry.elements import ELEMENTS

# d-block (Sc–Zn, Y–Cd, La, Hf–Hg): neither RDKit radicals nor the xyz default give the spin.
_D_BLOCK = frozenset((*range(21, 31), *range(39, 49), 57, *range(72, 81)))


def check_electronic_state(symbols: Iterable[str], charge: int, multiplicity: int, *,
                           declared: bool = True) -> None:
    """Raise ValueError for an element outside the table (``unsupported_element``), a species
    with a d-block element whose multiplicity was not declared (``declare_multiplicity``), or a
    multiplicity the electron count cannot have."""

    names = list(symbols)
    unsupported = sorted({name for name in names if name not in ELEMENTS})
    if unsupported:
        raise ValueError(f"unsupported_element:{','.join(unsupported)}")
    metals = sorted({name for name in names if ELEMENTS[name].z in _D_BLOCK})
    if metals and not declared:
        raise ValueError(f"declare_multiplicity: d-block element(s) {','.join(metals)}")
    if multiplicity < 1:
        raise ValueError(f"Multiplicity must be positive, got {multiplicity}")
    electron_count = sum(ELEMENTS[name].z for name in names) - charge
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


def low_spin_coupled(component_multiplicities: Iterable[int], multiplicity: int) -> bool:
    """Whether at least two open-shell components couple below their high-spin multiplicity
    Σ(m_i − 1) + 1 (CH3· + O2 as a doublet; its quartet is high-spin).

    UKS describes the bound isomers of such a composition, but where the components separate
    its solution is broken-symmetry and spin-contaminated.
    """

    unpaired = [m - 1 for m in component_multiplicities if m > 1]
    return len(unpaired) >= 2 and multiplicity < sum(unpaired) + 1


def composition_multiplicity(component_multiplicities: Iterable[int],
                             declared: int | None) -> int:
    """The declared multiplicity of a composition, or the only one coupling allows.

    A low-spin coupling is accepted. Raise ValueError when several multiplicities are allowed
    and none is declared, when the declared one is not allowed, or for a low-spin-coupled
    singlet (``low_spin_singlet_unsupported``: a singlet runs as RKS, which has no
    broken-symmetry solution; declare a recombination product as a monomer instead).
    """

    components = list(component_multiplicities)
    allowed = coupled_multiplicities(components)
    listed = ", ".join(map(str, allowed))
    if declared is None:
        if len(allowed) > 1:
            raise ValueError(f"declare_multiplicity: candidates ({listed})")
        return allowed[0]
    if declared not in allowed:
        raise ValueError(f"multiplicity {declared} is not among the coupled ones ({listed})")
    if declared == 1 and low_spin_coupled(components, declared):
        raise ValueError("low_spin_singlet_unsupported: RKS has no broken-symmetry solution; "
                         "declare a recombination product as a monomer")
    return declared
