
import unittest
import matplotlib.pyplot as plt
import numpy as np

from crosssection.cross_section import CrossSectionCalcBeb, CrossSectionCalcBed, CrossSectionCalcVacNote, AtomFactory
"""

class TestCrossSectionCalc(unittest.TestCase):


    def test_behaviour(self):
        Atom = AtomFactory.get_helium()
        log_boundary_a = 4.
        log_boundary_b = 12.
        x = np.logspace(log_boundary_a, log_boundary_b, num = 50)
        y_bed = np.zeros(len(x))
        y_beb = np.zeros(len(x))
        y_vac_note = np.zeros(len(x))
        y_bethe = np.zeros(len(x))
        for i, T in enumerate(x):
            crossBed = CrossSectionCalcBed(T, atom = Atom)
            crossBeb = CrossSectionCalcBeb(T, atom = Atom)
            crossVacNote = CrossSectionCalcVacNote(T, atom = Atom)
            y_bed[i] = crossBed.calculate()
            y_beb[i] = crossBeb.calculate()
            y_bethe[i] = crossBed.bethe_asymptotic()
            y_vac_note[i] = crossVacNote.calculate()
        #plt.semilogy(x, np.repeat(2.64e-19, len(x)))
        #plt.semilogy(np.repeat(1.4e7, 100), np.linspace(1e-19,2.64e-19, 100))
        plt.semilogx(x, y_bed, label = 'BED model')
        plt.semilogx(x, y_beb, label = 'BEB model')
        plt.semilogx(x, y_vac_note, label = 'Vacuum Note Calculation', linestyle = ':')
        plt.semilogx(x, y_bethe, label = 'Asymptotic behaviour (Bethe theory)')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    unittest.main()
"""