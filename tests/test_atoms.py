from crosssection.atom import Atom
import numpy as np


def test_atom_neon():
    neon = Atom(10)

    assert neon.Z == 10.
    assert neon.B == np.array([21.7, 48.47, 866.9])
    assert neon.N == np.array([6, 2, 2])

    assert neon.Ni == np.array([6.963, 0.7056, 1.686])
    assert neon.Mi == np.array([1.552, 4.8e-2, 1.642e-2])
