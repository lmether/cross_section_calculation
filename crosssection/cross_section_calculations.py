import numpy as np

from crosssection.atom import AtomEnergyProperty


def calculate_crosssection(atom, energy):
    cross_section = None  # TODO

    bla = atom.Z

    return cross_section



def calculate_modified_oscillator_strength(atom, T, w, n_shell):
    atom_energy_properties = AtomEnergyProperty(T, atom)

    threshold = ((w) + 1.0) * atom_energy_properties.atom.B[n_shell]
    coefficient_matrix = (
        atom_energy_properties.atom.ai_below if threshold < (48.47) else atom_energy_properties.atom.ai_above
    )
    quotient = [1 / ((w) + (1.0)) ** i for i in range(1, 8)]
    osc_str = np.reshape(np.dot(coefficient_matrix, np.array(quotient)), (len(atom_energy_properties.atom.B)))
    return 1 / (w + 1) * osc_str[n_shell]
