import unittest
import numpy as np
from cross_section import NatureConstants, total_cross_section_bebvm, IrgendeinGebilde



class TestIrgendeinGebilde(unittest.TestCase):
    def test_irgendwas(self):
        irgendwas = IrgendeinGebilde(0)
        print(irgendwas.t)
        assert 0. == irgendwas.t.all()
        assert 0. == irgendwas.tbar
        assert 0. == irgendwas.beta_t
        print(irgendwas.f_2.all())
        print(irgendwas.f_3.all())
        print(irgendwas.phi)
        assert 0. == irgendwas.T

class TestStringMethods(unittest.TestCase):


    def test_total_cross_section_bebvm_LHC(self):
        nc = NatureConstants
        ig = IrgendeinGebilde(3810)
        w_max = ((ig.T + nc.U_bebvm[0:3]) / nc.B - 1) / 2  # Devided by 2 since the outgoing electrons are indistinguishable
        print(total_cross_section_bebvm(w_max))

        assert True



    def test_total_cross_section_bebvm_0(self):
        nc = NatureConstants
        ig = IrgendeinGebilde(0)

        for i in range(len(nc.U_bebvm)):
            w_max = ((ig.T + nc.U_bebvm[i]) / nc.B - 1) / 2  # Devided by 2 since the outgoing electrons are indistinguishable
            print(total_cross_section_bebvm(w_max))

        assert True


    def test_total_cross_section_bebvm_negativ(self):
        nc = NatureConstants
        ig = IrgendeinGebilde(-1)
        for i in range(len(nc.U_bebvm)):
            w_max = ((ig.T + nc.U_bebvm[i]) / nc.B - 1) / 2  # Devided by 2 since the outgoing electrons are indistinguishable
            print(total_cross_section_bebvm(w_max))

        assert True


if __name__ == '__main__':
    unittest.main()