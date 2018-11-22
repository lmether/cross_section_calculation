import unittest

from cross_section import CrossSectionCalcBed, AtomFactory
import matplotlib.pyplot as plt
import numpy as np

"""
class TestStringMethods(unittest.TestCase):


    def test_total_cross_section_bebvm_LHC(self):
        parameters = [0, 3810, 9999]

        for parameter in parameters:
            cross = CrossSectionCalc(parameter)
            print(cross.total_cross_section_bebvm())
            #cross.differential_cross_section_subshells_bed()
            #cross.total_cross_section_bed()

        assert True


    def test_integ(self):
        parameters = [0, 3810, 9999]
        cross = CrossSectionCalcBed(3.81)

        assert True
"""

class TestCrossSectionCalc(unittest.TestCase):

    '''
    def test_init(self):
        cross = CrossSectionCalcBed(0, atom=AtomFactory.get_hydrogen())
        for w_max in cross.w_max:
            self.assertEquals(w_max, 0)
        for elem in cross.t:
            self.assertEquals(elem, 0)
        crossH = CrossSectionCalcBed(13.6057, atom=AtomFactory.get_hydrogen())
        for w_max in crossH.w_max:
            self.assertEquals(w_max, 0.5e9)
        for elem in crossH.t:
            self.assertEquals(elem, 1)
        crossH.calculate()
        '''

    def test_behaviour(self):

        x = np.logspace(2.,9., num = 1000)
        y = np.zeros(len(x))
        for i, T in enumerate(x):
            cross = CrossSectionCalcBed(T, atom=AtomFactory.get_neon())
            y[i] = cross.calculate()
        plt.semilogx(x,y)
        plt.show()


    def test_osc_strength(self):
        cross = CrossSectionCalcBed(1e1,atom=AtomFactory.get_helium())
        coeff1 = cross.atom.ai_below[0]
        coeff2 = cross.atom.ai_above[0]
        x_values = np.logspace(0.0, 1.0, num = 100000)
        x_values_powermatrix = np.array([x_values**i for i in range(1,8)])
        func1_vals = np.dot(coeff1, x_values_powermatrix)
        func2_vals = np.dot(coeff2, x_values_powermatrix)
        plt.semilogx(x_values, func1_vals)
        plt.show()
        #assert all(func1_val > 0 for func1_val in func1_vals) and all(func2_val > 0 for func2_val in func2_vals)




if __name__ == '__main__':
    unittest.main()