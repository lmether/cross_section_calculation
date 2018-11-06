import unittest

from cross_section import CrossSectionCalc


class TestStringMethods(unittest.TestCase):


    def test_total_cross_section_bebvm_LHC(self):
        parameters = [0, 3810, 9999]

        for parameter in parameters:
            cross = CrossSectionCalc(parameter)
            print(cross.total_cross_section_bebvm())
            #cross.differential_cross_section_subshells_bed()
            #cross.total_cross_section_bed()

        assert True




if __name__ == '__main__':
    unittest.main()