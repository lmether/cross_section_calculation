from crosssection.atom_calculations import calc_B, calc_N
from crosssection.constants import c, m_e, e_charge, a_0, R
import numpy as np


class Atom:
    def __init__(self, Ne):
        self.Z = Ne
        self.B = calc_B(Ne)
        self.N = calc_N(Ne)
        self.Ni = 3
        self.Mi = 4
        self.ai_below = 5
        self.ai_above = 6
        self.U_bed = 7


rest_en_eV = (c ** 2) * m_e / e_charge


class AtomEnergyProperty:
    def __init__(self, atom, T):
        """
        Calculate additional properties of an Atom with given Energy

        :param atom : Atom instance
        :param T    : Energy
        """

        self.T_rel = T  # (1. - 1./((T) / self.rest_en_eV) ** 2 ) ** (1. / 2.)
        self.T = rest_en_eV / 2. * (1. - (rest_en_eV / (self.T_rel + rest_en_eV)) ** 2.)
        self.t = self.T / atom.B  # T in units of binding energy(vector with entries for each subshell)
        self.w_max = (self.t - 1) / 2
        self.u = atom.U_bed / atom.B  # Constant factor from derivation
        self.S = 4 * np.pi * a_0 ** 2 * atom.N * (R / atom.B) ** 2
