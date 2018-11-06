
import sys
import numpy as np


# ACCORDING TO CERN VACUUM TECHNICAL NOTE 96-01 (1996) THE CROSS SECTIONS OF 7TeV PROTONS IS EQUAL TO THE ONE OF 3.81GeV ELECTRONS. 
# THUS AN ADAPTED VERSION OF THE BINARY ENCOUNTER BETHE VRIENS MODEL IS USED TO DETERMINE THE DIFFERENTIAL CROSS SECTION, 
# WHOSE ANALYTICAL FORM CAN BE INTGRATED AS BELOW, TO OBTAIN THE TOTAL CROSS SECTION.


########################################## CONSTANTS

m_e = 9.11e-31		# Electron mass in kg
c = 3e8			# Speed of light in m/s
a_0 = 5.29e-11		# Bohr radius in m
alpha = 1./137 		# Fine structure constant

########################################## INPUT VALUES

B = np.array([870.2, 48.5, 21.7]) 									# Bound electron binding energy in eV
T = float(sys.argv[2])*1.e6				# Incident electron kinetic energy in MeV
#W = float(sys.argv[3])*1.e6			# Ejected electron kinetic energy in MeV
U = float(sys.argv[3])					# Bound electron kinetic energy in eV
N = np.array([2, 2, 6])							# Total occupation number of requested atom

# For Neon:
# B ~ 8.032 MeV; T ~ 3810 MeV ; U ~ 5 MeV ; N ~ 8

########################################## COMPOSED QUANTITIES

t = T/B
#w = W/B
u = U/B
tbar = T/(m_e * c**2)
bbar = B/(m_e * c**2)
ubar = U/(m_e * c**2)
beta_t = 1 - 1/(1 + tbar)**2
beta_b = 1 - 1/(1 + bbar)**2
beta_u = 1 - 1/(1 + ubar)**2
phi = np.cos( np.sqrt( alpha**2 / (beta_t + beta_b) * np.log(beta_t / beta_b) ) )
w_max = ( (T + U) / B - 1 ) / 2			# Devided by 2 since the outgoing electrons are indistinguishable and the exchange term is already included

factor_1 = 4 * np.pi * alpha**4 * a_0**2 * N / (2 * bbar * ( beta_t + (beta_b + beta_u) / 2 ))

factor_2 = - phi / (t + 1) * (1 + 2 * tbar) / (1 + tbar / 2)**2

factor_3 = np.log( beta_t * (1 + tbar)**2 ) - beta_t - np.log(2 * bbar)

print(t, tbar, beta_t, factor_1, factor_2, factor_3)

########################################## INTEGRATED VERSION OF BEBVM INTEGRATED FROM A FINAL EJECTED KINETIC ENERGY OF 0 TO w_max = (T + U) / B - 1 TO OBTAIN TOTAL CROSS SECTION. ORIGINAL FORM OF THE DIFFERENTIAL CROSS SECTION:
########################################## d_sigma/d_w = factor_1 * ( factor_2 * (1 / (w + 1) + 1 / (t - w)) + 1 / (w + 1)**2 + 1 / (t - w)**2 + bbar**2 / (1 + tbar / 2)**2 + factor_3 * (1 / (w + 1)**3 + 1 / (t - w)**3) )

def total_cross_section(w, t = t, tbar = tbar, bbar = bbar, factor_1 = factor_1, factor_2 = factor_2, factor_3 = factor_3):
	total_cross_sec = np.sum(factor_1 * ( factor_2 * (np.log((w + 1) / (t - w) * t )) - 1 / (w + 1) + 1 + 1 / (t - w) - 1 / t  + bbar**2 * w / (1 + tbar / 2)**2 + factor_3 / 2 * (1 / (t - w)**2 - 1 / t**2 - 1 / (w + 1)**2 + 1 ) ) )
	return 1e6 * total_cross_sec
		

print(total_cross_section(w_max))
