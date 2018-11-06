import sys
import numpy as np
from scipy.integrate import quad

# ACCORDING TO CERN VACUUM TECHNICAL NOTE 96-01 (1996) THE CROSS SECTIONS OF 7TeV PROTONS IS EQUAL TO THE ONE OF 3.81GeV ELECTRONS.
# THUS AN ADAPTED VERSION OF THE BINARY ENCOUNTER BETHE VRIENS MODEL IS USED TO DETERMINE THE DIFFERENTIAL CROSS SECTION,
# WHOSE ANALYTICAL FORM CAN BE INTEGRATED AS BELOW, TO OBTAIN THE TOTAL CROSS SECTION.


########################################## CONSTANTS


########################################## INPUT VALUES

########################################## ANSATZ 1:	BEBVM


########################################## INTEGRATED VERSION OF BEBVM INTEGRATED FROM A FINAL EJECTED KINETIC ENERGY OF 0 TO w_max = ( (T + U) / B - 1 ) / 2 TO OBTAIN TOTAL CROSS SECTION. ORIGINAL FORM OF THE DIFFERENTIAL CROSS SECTION:
########################################## d_sigma/d_w = f_1 * ( f_2 * (1 / (w + 1) + 1 / (t - w)) + 1 / (w + 1)**2 + 1 / (t - w)**2 + bbar**2 / (1 + tbar / 2)**2 + factor_3 * (1 / (w + 1)**3 + 1 / (t - w)**3) )


############################### COMPOSED QUANTITIES


class NatureConstants:
    m_e = 9.11e-31  # Electron mass in kg
    c = 3e8  # Speed of light in m/s
    a_0 = 5.29e-11  # Bohr radius in m
    alpha = 1. / 137  # Fine structure constant
    R = 13.6  # Rydberg constant in eV

    B = np.array([21.7, 48.47, 866.9])  # Bound electron binding energy in eV w.r.t the shell
    N = np.array([6, 2, 2])  # Total occupation number of requested atom w.r.t the shell

    U_bebvm = np.linspace(0, 800, num=8, endpoint=True)  # Bound electron kinetic energy in eV of the BEBVM approximation (currently for testing)
    bbar = B / (m_e * c ** 2)
    ubar = U_bebvm / (m_e * c ** 2)
    beta_b = 1 - 1 / (1 + bbar) ** 2
    beta_u = 1 - 1 / (1 + ubar) ** 2

    Ni = np.array([6.963, 0.7056, 1.686])  # Values of Ni for requested atom, for each subshell
    Mi = np.array([1.552, 4.8e-2, 1.642e-2])  # Values of Mi for requested atom, for each subshell

    ai_below_2s_ion_threshold = np.array(
        [[4.8791, -2.882, -7.4711e-1, 0., 0., 0., 0.], [0., 1.7769, 2.8135, -3.151e1, 6.3469e1, -5.2528e1, 1.5982e1],
         [0., 0., 5.2475, -2.8121, 0., 0.,
          0.]])  # Coefficients in a 7-th order fit for dipole oscillator strength for energies below 2s ionization 																											# energy of  (originally scaled via units of bounding energy)

    ai_above_2s_ion_threshold = np.array([[0., -5.8514, 3.2930e2, -1.6788e3, 3.2985e3, -2.3250e3, 0.],
                                          [0., 1.7769, 2.8135, -3.151e1, 6.3469e1, -5.2528e1, 1.5982e1],
                                          [0., 0., 5.2475, -2.8121, 0., 0.,
                                           0.]])  # Coefficients in a 7-th order fit for dipole oscillator strength for energies above 2s ionization 																											# energy (originally scaled via units of bounding energy)

    U_bed = np.array([1.1602e2, 1.4188e2, 1.2591e3])  # Bound electron kinetic energy in eV in the BED model
    f_1_bed = 4 * np.pi * a_0 ** 2 * N * (R / B) ** 2  # Constant factor from derivation

nc = NatureConstants

class IrgendeinGebilde:
    def __init__(self, T):
        self.T = T
        self.t = T / nc.B  # T in units of binding energy
        self.tbar = T / (nc.m_e * nc.c ** 2)         # nicht
        self.beta_t = 1 - 1 / (1 + self.tbar) ** 2          # nicht
        self.phi = np.cos(np.sqrt(nc.alpha ** 2 / (self.beta_t + nc.beta_b) * np.log(self.beta_t / nc.beta_b)))      # nicht
        self.f_2 = - self.phi / (self.t + 1) * (1 + 2 * self.tbar) / (1 + self.tbar / 2) ** 2  # nicht
        self.f_3 = np.log(self.beta_t * (1 + self.tbar) ** 2) - self.beta_t - np.log(2 * nc.bbar)  # nicht


# w_max = ( (T + U) / B - 1 ) / 2			#Devided by 2 since the outgoing electrons are indistinguishable and the exchange term is already included
# f_1 = 4 * np.pi * alpha**4 * a_0**2 * N / (2 * bbar * ( beta_t + (beta_b + beta_u) / 2 ))

#T = float(sys.argv[1]) * 1.e6  # Incident electron kinetic energy in MeV

#irgendein = IrgendeinGebilde(T)

def total_cross_section_bebvm(w, ig = IrgendeinGebilde(3810)) :
    data_set_len = len(nc.U_bebvm)
    total_cross_sec = np.zeros(data_set_len)
    for i in range(data_set_len):
        f_1 = 4 * np.pi * nc.alpha ** 4 * nc.a_0 ** 2 * nc.N / (2 * nc.bbar * (ig.beta_t + (nc.beta_b + nc.beta_u[i])) / 2)
        total_cross_sec[i] = np.sum(f_1 * (
                    ig.f_2 * (np.log((w + 1) / (ig.t - w) * ig.t)) - 1 / (w + 1) + 1 + 1 / (ig.t - w) - 1 / ig.t + nc.bbar ** 2 * w / (
                        1 + ig.tbar / 2) ** 2 + ig.f_3 / 2 * (1 / (ig.t - w) ** 2 - 1 / ig.t ** 2 - 1 / (w + 1) ** 2 + 1)))
    return 1.e6 * total_cross_sec

'''
for i in range(len(nc.U_bebvm)):
    w_max = ((T + nc.U_bebvm[i]) / B - 1) / 2  # Devided by 2 since the outgoing electrons are indistinguishable
    print(np.log(total_cross_section_bebvm(w_max)))
'''
########################################## ANSATZ 1:	BED

########################################## INTEGRATE THE VERSION OF BED FROM A FINAL EJECTED KINETIC ENERGY OF 0 TO w_max = ( (T + U) / B - 1 ) / 2 TO OBTAIN TOTAL CROSS SECTION. ORIGINAL FORM OF THE DIFFERENTIAL CROSS SECTION:
########################################## d_sigma/d_W = f_1 / (t + u + 1) * ( (Ni/N - 2) / (t + 1) * (1 / (w + 1) + 1 / (t - w) ) + (2 - Ni / N) * (1 / (t - w)**2 + 1 / (w + 1)**2) + np.log(t) / (N * (w + 1)) * osc_str)





def differential_cross_section_subshells_bed(w, ig = IrgendeinGebilde(3810)):
    N_shells = len(w)  # Number of shells (equivalent to number of integrations in the end)
    u = nc.U_bed / nc.B
    osc_str = np.zeros(N_shells)
    for i in range(N_shells):
        ai = nc.ai_below if (w[i] + u[i] - 1.) * B[
            i] < 48.47 else nc.ai_above  # different coefficients for different engergy regimes, W+U-B heeby corresponds to the
        osc_str = np.dot(ai, np.array([-1 / (w + u - 1.) ** (i + 2) for i in range(
            6)]))  # Dipole oscillator strength (ionization acc. to QED happens via 																															# exchange of photon)

        diff_cross_sec_subshells = nc.f_1_bed / (ig.t + u + 1.) * (
                    (nc.Ni / nc.N - 2.) / (ig.t + 1.) * (1. / (w + 1.) + 1. / (ig.t - w)) + (2. - nc.Ni / nc.N) * (
                        1. / (ig.t - w) ** 2 + 1. / (w + 1.) ** 2.) + np.log(ig.t) / nc.N * np.dot(1. / (w + 1.),
                                                                                          osc_str))  # Vector with corresponding cross sections for ionization from every 																															# subshell

    return 1.e4 * diff_cross_sec_subshells  # Converting from m**2 to cm**2


def total_cross_section_bed(ig = IrgendeinGebilde(3810)):
    w_max = ((ig.T + nc.U) / nc.B - 1) / 2
    integrated_cross_sec = np.zeros(len(w_max))
    for i in range(len(w_max)):
        #total_cross_sec
        print(np.log(total_cross_section_bed(w_max)))


