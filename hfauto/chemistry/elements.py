"""The one element table: Z = 1–57 and 72–86 (no Ce–Lu, nothing beyond Rn).

This is the supported range: ``check_electronic_state`` rejects any other symbol at the
structures stage, so everything downstream may index the table directly.

Columns: symbol, Z, most-abundant-isotope mass (amu), covalent radius (Å), van der Waals
radius (Å), maximum coordination. H–Kr and I keep the values hfauto has always used: masses
from NIST AME2016, covalent radii from Cordero et al., Dalton Trans. 2008, van der Waals radii
from Bondi, J. Phys. Chem. 1964 (Be, B, Al, Ca, Ge from Mantina et al., J. Phys. Chem. A 2009).
The other rows were generated once from the RDKit PeriodicTable (masses; its covalent radii are
Cordero 2008), with van der Waals radii from Alvarez, Dalton Trans. 2013 (also Sc–Co; Po, At and
Rn, which Alvarez does not cover, from Mantina 2009). RDKit is not needed at run time.

The maximum coordination is the number of bonded neighbours a reaction trial may leave on an
atom: H 1, B/C/N 4, O 3, halogens 1, the other p-block elements of period ≥ 3 6, metals 9,
He/Ne 0. Hypervalent iodine is out of scope.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Element:
    symbol: str
    z: int
    mass: float  # amu, most abundant isotope
    covalent_radius: float  # Å
    vdw_radius: float  # Å
    max_coordination: int


_TABLE = """
H    1    1.007825  0.31  1.20  1
He   2    4.002603  0.28  1.40  0
Li   3    7.016003  1.28  1.82  9
Be   4    9.012183  0.96  1.53  9
B    5   11.009305  0.84  1.92  4
C    6   12.000000  0.76  1.70  4
N    7   14.003074  0.71  1.55  4
O    8   15.994915  0.66  1.52  3
F    9   18.998403  0.57  1.47  1
Ne  10   19.992440  0.58  1.54  0
Na  11   22.989769  1.66  2.27  9
Mg  12   23.985042  1.41  1.73  9
Al  13   26.981539  1.21  1.84  6
Si  14   27.976927  1.11  2.10  6
P   15   30.973762  1.07  1.80  6
S   16   31.972071  1.05  1.80  6
Cl  17   34.968853  1.02  1.75  1
Ar  18   39.962383  1.06  1.88  6
K   19   38.963706  2.03  2.75  9
Ca  20   39.962591  1.76  2.31  9
Sc  21   44.955908  1.70  2.58  9
Ti  22   47.947942  1.60  2.46  9
V   23   50.943957  1.53  2.42  9
Cr  24   51.940506  1.39  2.45  9
Mn  25   54.938044  1.39  2.45  9
Fe  26   55.934936  1.32  2.44  9
Co  27   58.933194  1.26  2.40  9
Ni  28   57.935342  1.24  1.63  9
Cu  29   62.929598  1.32  1.40  9
Zn  30   63.929142  1.22  1.39  9
Ga  31   68.925574  1.22  1.87  6
Ge  32   73.921178  1.20  2.11  6
As  33   74.921595  1.19  1.85  6
Se  34   79.916522  1.20  1.90  6
Br  35   78.918338  1.20  1.85  1
Kr  36   83.911498  1.16  2.02  6
Rb  37   84.911790  2.20  3.21  9
Sr  38   87.905612  1.95  2.84  9
Y   39   88.905848  1.90  2.75  9
Zr  40   89.904704  1.75  2.52  9
Nb  41   92.906378  1.64  2.56  9
Mo  42   97.905408  1.54  2.45  9
Tc  43   96.906365  1.47  2.44  9
Ru  44  101.904349  1.46  2.46  9
Rh  45  102.905504  1.42  2.44  9
Pd  46  105.903486  1.39  2.15  9
Ag  47  106.905097  1.45  2.53  9
Cd  48  113.903358  1.44  2.49  9
In  49  114.903878  1.42  2.43  6
Sn  50  119.902195  1.39  2.42  6
Sb  51  120.903816  1.39  2.47  6
Te  52  129.906224  1.38  1.99  6
I   53  126.904472  1.39  1.98  1
Xe  54  131.904154  1.40  2.06  6
Cs  55  132.905452  2.44  3.48  9
Ba  56  137.905247  2.15  3.03  9
La  57  138.906353  2.07  2.98  9
Hf  72  179.946550  1.75  2.63  9
Ta  73  180.947996  1.70  2.53  9
W   74  183.950931  1.62  2.57  9
Re  75  186.955753  1.51  2.49  9
Os  76  191.961481  1.44  2.48  9
Ir  77  192.962926  1.41  2.41  9
Pt  78  194.964791  1.36  2.29  9
Au  79  196.966569  1.36  2.32  9
Hg  80  201.970643  1.32  2.45  9
Tl  81  204.974427  1.45  2.47  6
Pb  82  207.976652  1.46  2.60  6
Bi  83  208.980399  1.48  2.54  6
Po  84  208.982430  1.40  1.97  6
At  85  209.987148  1.50  2.02  1
Rn  86  222.017571  1.50  2.20  6
"""

ELEMENTS: dict[str, Element] = {
    row[0]: Element(row[0], int(row[1]), float(row[2]), float(row[3]), float(row[4]),
                    int(row[5]))
    for row in (line.split() for line in _TABLE.strip().splitlines())
}


def atomic_number(symbol: str) -> int:
    return ELEMENTS[symbol].z


def mass(symbol: str) -> float:
    return ELEMENTS[symbol].mass


def covalent_radius(symbol: str) -> float:
    return ELEMENTS[symbol].covalent_radius


def vdw_radius(symbol: str) -> float:
    return ELEMENTS[symbol].vdw_radius


def max_coordination(symbol: str) -> int:
    return ELEMENTS[symbol].max_coordination
