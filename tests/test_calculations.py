from crosssection.atom import Atom
from crosssection.calculations import calculate_crosssection
from crosssection.cross_section import AtomFactory


def test_cacluate_cross_section_neon():
    neon = AtomFactory.get_neon()

    value = calculate_crosssection(neon, 15)
    expected_value = 83421
    assert expected_value == value


def test_cacluate_cross_section_custom():
    atom = Atom(8)

    value = calculate_crosssection(atom, 15)
    expected_value = 83421
    assert expected_value == value
