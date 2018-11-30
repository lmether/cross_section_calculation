import unittest

from cross_section import CrossSectionCalcBed, AtomFactory
import matplotlib.pyplot as plt
import numpy as np


class TestCrossSectionCalc(unittest.TestCase):


    def test_Mi_calc(self):
        log_boundary_a = 4.
        log_boundary_b = 11.
        x = np.logspace(log_boundary_a, log_boundary_b, num=10)
        y = np.zeros(len(x))
        z = np.zeros(len(x))
        for i, T in enumerate(x):
            cross = CrossSectionCalcBed(3.81e9, atom=AtomFactory.get_helium())
            y[i] = cross.Mi_calculation(0, T, 0)
            sample_array = np.linspace(0, T, num = 1000000)
            for j in range(len(sample_array) - 1):
                z[i] += cross.Mi_calculation(sample_array[j], sample_array[j + 1], 0)
            print(z - y)
        plt.semilogx(x, y)
        plt.semilogx(x, z)
        plt.show()


    """
    
       def test_behaviour_full_bed(self):

        x = np.logspace(1.3, 3.3, num=100)
        y = np.zeros(len(x))
        for i, T in enumerate(x):
            cross = CrossSectionCalcBed(T, atom=AtomFactory.get_hydrogen())
            y[i] = cross.calculate()
        plt.semilogx(x, y)
        plt.show()

    
        def test_mod_osc_str_behaviour(self):
        log_boundary_a = 5.
        log_boundary_b = 6.
        x = np.logspace(log_boundary_a, log_boundary_b, num=1000)
        y = np.zeros(len(x))
        for i, T in enumerate(x):
            cross = CrossSectionCalcBed(3.81e9, atom=AtomFactory.get_neon())
            y[i] = cross.calculate_modified_oscillator_strength(T, 2)
        plt.semilogx(x, y)
        plt.show()


    
    
    def test_asymptotic_behaviour(self):
        x = np.logspace(8., 10., num=100)
        full_bed = np.zeros(len(x))
        bethe = np.zeros(len(x))
        for i, T in enumerate(x):
            cross = CrossSectionCalcBed(T, atom=AtomFactory.get_hydrogen())
            full_bed[i] = cross.calculate()
            bethe[i] = cross.bethe_asymptotic()
        plt.semilogx(x,bethe - full_bed)
        plt.show()
    """


if __name__ == '__main__':
    unittest.main()
