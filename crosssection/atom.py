from crosssection.constants import c, m_e, e_charge, a_0, R
import numpy as np


class Atom:
    def __init__(self, Ne):
        self.Z = None
        self.B = None
        self.N = None
        self.Ni = None
        self.Mi = None
        self.ai_below = None
        self.ai_above = None
        self.U_bed = None
        self.M_squared = None
        self.C = None
        self.corr_fact = None


rest_en_eV = (c ** 2) * m_e / e_charge


class AtomEnergyProperty:
    def __init__(self, atom, T):
        """
        Calculate additional properties of an Atom with given Energy

        :param atom: Atom instance
        :param T: Energy
        """

        self.T_rel = T  # (1. - 1./((T) / self.rest_en_eV) ** 2 ) ** (1. / 2.)
        self.T = rest_en_eV / 2. * (1. - (rest_en_eV / (self.T_rel + rest_en_eV)) ** 2.)
        self.t = self.T / atom.B  # T in units of binding energy(vector with entries for each subshell)
        self.w_max = (self.t - 1) / 2
        self.u = atom.U_bed / atom.B  # Constant factor from derivation
        self.S = 4 * np.pi * a_0 ** 2 * atom.N * (R / atom.B) ** 2
