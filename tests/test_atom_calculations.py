from crosssection.atom import Atom
import numpy as np

from crosssection.atom_calculations import calc_Z, calc_B, calc_N, calc_Ni, calc_Mi
from crosssection.constants import HELIUM, NEON


# Helium
def test_calc_Z():
    assert calc_Z(HELIUM) == HELIUM
    assert calc_Z(NEON) == NEON


def test_calc_B():
    value = calc_B(NEON)
    assert list(value) == [21.7, 48.47, 866.9]


def test_calc_N():
    assert list(calc_N(NEON)) == [6, 2, 2]


def test_calc_Ni():
    assert list(calc_Ni(NEON)) == [6.963, 0.7056, 1.686]


def test_calc_Mi():
    assert list(calc_Mi(NEON)) == [1.552, 4.8e-2, 1.642e-2]
