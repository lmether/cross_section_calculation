#from decimal import *
import sys
import math
import numpy as np
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt



def crosssection(energy, M_squared=0.695, C=8.115, Z=1., correction_factor=1.2):
	
	en_joule = energy * 1.6e-19
	m0 = 1.67e-27
	c = 3e8
	rest_en = m0 * c**2.	
	beta = (1. - (rest_en / en_joule )**2.)**(1./2.)
	gamma = beta / (1. - beta**2.)**(1./2.)
	x = 2. * math.log(gamma) - beta**2.
	omega = 1.874e-20 * (Z/beta)**2. * (M_squared * x + C) 				# from m^2 to cm^2
	omega = omega * correction_factor						# Correction due to flaws of the measurement, extrapolated from measurmenets at 26GeV and stupidly assumed to be constant because FUCK ANY THEORETICAL THOUGHTS ON THIS, but I'm gonna take this anyways cause I used it 
											# already before finding out about this fauxpas and I really don't wanna start over again. 
	
	print('Beta: ' + str(beta) + '  Cross section: ' + str(omega))


en = float(sys.argv[1])


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
