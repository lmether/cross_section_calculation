import sys
import numpy as np
from scipy.integrate import quad


# ACCORDING TO CERN VACUUM TECHNICAL NOTE 96-01 (1996) THE CROSS SECTIONS OF 7TeV PROTONS IS EQUAL TO THE ONE OF 3.81GeV ELECTRONS.
# THUS AN ADAPTED VERSION OF THE BINARY ENCOUNTER BETHE VRIENS MODEL IS USED TO DETERMINE THE DIFFERENTIAL CROSS SECTION,
# self.w_maxHOSE ANALYTICAL FORM CAN BE INTEGRATED AS BELOW, TO OBTAIN THE TOTAL CROSS SECTION.


########################################## CONSTANTS


########################################## INPUT VALUES

########################################## ANSATZ 1:	BEBVM


########################################## INTEGRATED VERSION OF BEBVM INTEGRATED FROM A FINAL EJECTED KINETIC ENERGY OF 0 TO w_max = ( (T + U) / B - 1 ) / 2 TO OBTAIN TOTAL CROSS SECTION. ORIGINAL FORM OF THE DIFFERENTIAL CROSS SECTION:
########################################## d_sigma/d_w = f_1 * ( f_2 * (1 / (w + 1) + 1 / (t - w)) + 1 / (w + 1)**2 + 1 / (t - w)**2 + bbar**2 / (1 + tbar / 2)**2 + factor_3 * (1 / (w + 1)**3 + 1 / (t - w)**3) )


############################### COMPOSED QUANTITIES

class Atom:
    def __init__(self, B, N, Ni, Mi, ai_below, ai_above, U_bed):
        self.B = B
        self.N = N
        self.Ni = Ni
        self.Mi = Mi
        self.ai_below = ai_below
        self.ai_above = ai_above
        self.U_bed = U_bed


class AtomFactory(object):
    @staticmethod
    def get_neon():
        B = np.array([21.7, 48.47, 866.9])  # Bound electron binding energy in eV w.r.t the shell
        N = np.array([6, 2, 2])  # Total occupation number of requested atom w.r.t the shell

        Ni = np.array([6.963, 0.7056, 1.686])  # Values of Ni for requested atom, for each subshell
        Mi = np.array([1.552, 4.8e-2, 1.642e-2])  # Values of Mi for requested atom, for each subshell

        ai_below_2s_ion_threshold = np.array(
            [[4.8791, -2.882, -7.4711e-1, 0., 0., 0., 0.],
             [0., 1.7769, 2.8135, -3.151e1, 6.3469e1, -5.2528e1, 1.5982e1],
             [0., 0., 5.2475, -2.8121, 0., 0.,
              0.]])  # Coefficients in a 7-th order fit for dipole oscillator strength for energies below 2s ionization 																											# energy of  (originally scaled via units of bounding energy)

        ai_above_2s_ion_threshold = np.array([[0., -5.8514, 3.2930e2, -1.6788e3, 3.2985e3, -2.3250e3, 0.],
                                              [0., 1.7769, 2.8135, -3.151e1, 6.3469e1, -5.2528e1, 1.5982e1],
                                              [0., 0., 5.2475, -2.8121, 0., 0.,
                                               0.]])  # Coefficients in a 7-th order fit for dipole oscillator strength for energies above 2s ionization 																											# energy (originally scaled via units of bounding energy)

        U_bed = np.array([1.1602e2, 1.4188e2, 1.2591e3])  # Bound electron kinetic energy in eV in the BED model

        return Atom(B, N, Ni, Mi, ai_below_2s_ion_threshold, ai_above_2s_ion_threshold, U_bed)

    @staticmethod
    def get_helium(self):
        # TODO
        pass



class CrossSectionCalc:
    # here die constants

    m_e = 9.11e-31  # Electron mass in kg
    c = 3e8  # Speed of light in m/s
    a_0 = 5.29e-11  # Bohr radius in m
    alpha = 1. / 137  # Fine structure constant
    R = 13.6  # Rydberg constant in eV


    def __init__(self, T, atom=AtomFactory.get_neon()):
        self.atom = atom
        self.f_1_bed = 4 * np.pi * self.a_0 ** 2 * atom.N * (self.R / atom.B) ** 2  # Constant factor from derivation
        self.bbar = atom.B / (self.m_e * self.c ** 2)
        self.ubar = atom.U_bed / (self.m_e * self.c ** 2)
        self.beta_b = 1 - 1 / (1 + self.bbar) ** 2
        self.beta_u = 1 - 1 / (1 + self.ubar) ** 2
        self.T = T
        self.t = T / atom.B  # T in units of binding energy
        self.tbar = T / (self.m_e * self.c ** 2)  # nicht
        self.beta_t = 1 - 1 / (1 + self.tbar) ** 2  # nicht
        self.phi = np.cos(np.sqrt(self.alpha ** 2 / (self.beta_t + self.beta_b) * np.log(self.beta_t / self.beta_b)))  # nicht
        self.f_2 = - self.phi / (self.t + 1) * (1 + 2 * self.tbar) / (1 + self.tbar / 2) ** 2  # nicht
        self.f_3 = np.log(self.beta_t * (1 + self.tbar) ** 2) - self.beta_t - np.log(2 * self.bbar)  # nicht
        self.w_max = ((self.T + atom.U_bed) / self.atom.B - 1) / 2


    def total_cross_section_bebvm(self):
        data_set_len = len(self.atom.U_bed)
        total_cross_sec = np.zeros(data_set_len)
        for i in range(data_set_len):
            f_1 = 4 * np.pi * self.alpha ** 4 * self.a_0 ** 2 * self.atom.N / (
                    2 * self.bbar * (self.beta_t + (self.beta_b + self.beta_u[i])) / 2)
            total_cross_sec[i] = np.sum(f_1 * (
                    self.f_2 * (np.log((self.w_max + 1) / (self.t - self.w_max) * self.t)) - 1 / (self.w_max + 1) + 1 + 1 / (
                    self.t - self.w_max) - 1 / self.t + self.bbar ** 2 * self.w_max / (
                            1 + self.tbar / 2) ** 2 + self.f_3 / 2 * (
                            1 / (self.t - self.w_max) ** 2 - 1 / self.t ** 2 - 1 / (self.w_max + 1) ** 2 + 1)))
        return 1.e6 * total_cross_sec

    def differential_cross_section_subshells_bed(self):
        N_shells = len(self.w_max)  # Number of shells (equivalent to number of integrations in the end)
        u = self.atom.U_bed / self.atom.B
        osc_str = np.zeros(N_shells)
        diff_cross_sec_subshells = np.zeros(N_shells)
        for i in range(N_shells):
            ai = self.atom.ai_below if (self.w_max[i] + u[i] - 1.) * self.atom.B[
                i] < 48.47 else self.atom.ai_above  # different coefficients for different engergy regimes, W+U-B heeby corresponds to the
            osc_str[i] = np.dot(ai, np.array([-1 / (self.w_max + u - 1.) ** (i + 2) for i in range(
                6)]))  # Dipole oscillator strength (ionization acc. to QED happens via 																															# exchange of photon)
            
            diff_cross_sec_subshells[i] = self.f_1_bed / (self.t + u + 1.) * (
                    (self.atom.Ni / self.atom.N - 2.) / (self.t + 1.) * (1. / (self.w_max + 1.) + 1. / (self.t - self.w_max)) + (2. - self.atom.Ni / self.atom.N) * (
                    1. / (self.t - self.w_max) ** 2 + 1. / (self.w_max + 1.) ** 2.) + np.log(self.t) / self.N * np.dot(1. / (self.w_max + 1.),
                                                                                                                 osc_str))  # Vector with corresponding cross sections for ionization from every 																															# subshell

        return 1.e4 * diff_cross_sec_subshells  # Converting from m**2 to cm**2

    def total_cross_section_bed(self):
        integrated_cross_sec = np.zeros(len(self.w_max))
        for i in range(len(self.w_max)):
            # total_cross_sec
            pass



