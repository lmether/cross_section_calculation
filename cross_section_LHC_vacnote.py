
import math



def crosssection(energy, M_squared=0.695, C=8.115, Z=1., correction_factor=1.2):
	
	en_joule = energy * 1.6e-19
	m0 = 1.67e-27
	c = 3e8
	rest_en = m0 * c**2.	
	beta = (1. - (rest_en / en_joule )**2.)**(1./2.)
	gamma = beta / (1. - beta**2.)**(1./2.)
	x = 2. * math.log(gamma) - beta**2.
	omega = 1.874e-24 * (Z/beta)**2. * (M_squared * x + C) * 1.e4		# from m^2 to cm^2
	omega = omega * correction_factor						# Correction factor extrapolated from measurements at 26GeV

	print('Beta: ' + str(beta) + '  Cross section: ' + str(omega))


print('H2:')
crosssection(en)
print('')
print('Helium:')
crosssection(en, M_squared=0.752, C=7.571, Z=1.)
print('')
print('Neon:')
crosssection(en, M_squared=2.02, C=18.17, Z=1., correction_factor = 1.5)
print('')
print('Nitrogen:')
crosssection(en, M_squared=3.74, C=34.84, Z=1., correction_factor = 1.5)
