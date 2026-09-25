"""Small, shared electronic-state validation for external QM backends."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

_ELEMENTS = """
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La
Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po
At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg
Cn Nh Fl Mc Lv Ts Og
""".split()  # noqa: SIM905 - compact periodic-table data is clearer as grouped text
_ATOMIC_NUMBERS = {symbol.upper(): number for number, symbol in enumerate(_ELEMENTS, 1)}


def _exact_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be an integer, got boolean {value}")
    try:
        parsed = int(value)
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be an integer, got {value!r}") from exc
    if numeric != float(parsed):
        raise ValueError(f"{label} must be an integer, got {value!r}")
    return parsed


def resolve_electronic_state(
    species_data: dict[str, Any] | None,
    method: dict[str, Any] | None = None,
    symbols: Iterable[str] | None = None,
) -> dict[str, int | None]:
    """Resolve charge/multiplicity once and reject contradictory inputs.

    Electron-parity validation is performed when every element is in the small
    built-in table.  Unknown elements leave ``electron_count`` unset rather than
    silently applying an incorrect atomic number.
    """

    data = species_data or {}
    settings = method or {}
    charge_values = {
        label: _exact_int(value, label)
        for label, value in (
            ("species.charge", data.get("charge")),
            ("species.formal_charge", data.get("formal_charge")),
            ("method.charge", settings.get("charge")),
            ("method.formal_charge", settings.get("formal_charge")),
        )
        if value is not None
    }
    if len(set(charge_values.values())) > 1:
        rendered = ", ".join(f"{key}={value}" for key, value in charge_values.items())
        raise ValueError(f"Conflicting charge specifications: {rendered}")
    charge = next(iter(charge_values.values()), 0)

    multiplicity_values = {
        label: _exact_int(value, label)
        for label, value in (
            ("species.multiplicity", data.get("multiplicity")),
            ("method.multiplicity", settings.get("multiplicity")),
        )
        if value is not None
    }
    if len(set(multiplicity_values.values())) > 1:
        rendered = ", ".join(f"{key}={value}" for key, value in multiplicity_values.items())
        raise ValueError(f"Conflicting multiplicity specifications: {rendered}")
    multiplicity = next(iter(multiplicity_values.values()), 1)
    if multiplicity < 1:
        raise ValueError(f"Multiplicity must be positive, got {multiplicity}")

    electron_count: int | None = None
    if symbols is not None:
        normalized = [str(symbol).strip().upper() for symbol in symbols]
        if all(symbol in _ATOMIC_NUMBERS for symbol in normalized):
            electron_count = sum(_ATOMIC_NUMBERS[symbol] for symbol in normalized) - charge
            unpaired = multiplicity - 1
            if electron_count < unpaired or (electron_count - unpaired) % 2:
                raise ValueError(
                    "Electron-count parity is incompatible with multiplicity: "
                    f"electrons={electron_count}, multiplicity={multiplicity}"
                )
    return {
        "charge": charge,
        "multiplicity": multiplicity,
        "uhf": multiplicity - 1,
        "electron_count": electron_count,
    }


def coupled_multiplicities(
    component_multiplicities: Iterable[int],
) -> tuple[int, ...]:
    """Return every spin multiplicity allowed by angular-momentum coupling.

    Multiplicity is represented as ``2S + 1``.  Working with ``2S`` keeps the
    calculation exact for both integer and half-integer fragment spins.
    """

    values = [
        _exact_int(value, "component multiplicity")
        for value in component_multiplicities
    ]
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


def select_coupled_multiplicities(
    component_multiplicities: Iterable[int],
    *,
    policy: str = "all",
    requested: Iterable[int] | None = None,
) -> tuple[int, ...]:
    """Select valid spin surfaces without silently guessing an invalid one."""

    allowed = coupled_multiplicities(component_multiplicities)
    if requested is not None:
        selected = tuple(
            dict.fromkeys(
                _exact_int(value, "requested multiplicity") for value in requested
            )
        )
        invalid = sorted(set(selected) - set(allowed))
        if invalid:
            raise ValueError(
                f"requested multiplicities {invalid} are not spin-coupling allowed; "
                f"allowed: {list(allowed)}"
            )
        if not selected:
            raise ValueError("requested multiplicities must not be empty")
        return selected
    normalized = str(policy).strip().lower()
    if normalized == "all":
        return allowed
    if normalized == "lowest":
        return allowed[:1]
    if normalized == "highest":
        return allowed[-1:]
    raise ValueError("spin_coupling_policy must be 'all', 'lowest', or 'highest'")
