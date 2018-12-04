'''
import sys
sys.path.append("")
import matplotlib
matplotlib.use('agg')
'''
import unittest
from cross_section import CrossSectionCalcBed, CrossSectionCalcBeb, AtomFactory
import matplotlib.pyplot as plt
import numpy as np


class TestCrossSectionCalc(unittest.TestCase):


    def test_behaviour(self):
        log_boundary_a = 9.
        log_boundary_b = 9.7
        x = np.logspace(log_boundary_a, log_boundary_b, num=100)
        y_bed = np.zeros(len(x))
        y_beb = np.zeros(len(x))
        y_bethe = np.zeros(len(x))
        Atom = AtomFactory.get_helium()
        for i, T in enumerate(x):
            crossBed = CrossSectionCalcBed(T, atom = Atom)
            crossBeb = CrossSectionCalcBeb(T, atom = Atom)
            y_bed[i] = crossBed.calculate()
            y_beb[i] = crossBeb.calculate()
            y_bethe[i] = crossBed.bethe_asymptotic()
        plt.semilogx(x, y_bed, label = 'BED model')
        plt.semilogx(x, y_beb, label = 'BEB model')
        plt.semilogx(x, y_bethe, label = 'Asymptotic behaviour (Bethe theory)')
        plt.legend()
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
