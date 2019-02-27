import matplotlib.pyplot as plt
import numpy as np

import crosssection.save_cross_section as scs

import argparse

from crosssection.cross_section import AtomFactory, CrossSectionCalcBed, CrossSectionCalcBeb, CrossSectionCalcVacNote

atom_library = {'h': AtomFactory.get_hydrogen(), 'h2': AtomFactory.get_h2(), 'he': AtomFactory.get_helium(),
                   'ne': AtomFactory.get_neon(), 'n2': AtomFactory.get_nitrogen()}

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Description of optional Parameters."
    )

    parser.add_argument(
        "-v",
        "--version",
        help="Version number",
        action="store_true",
    )

    parser.add_argument(
        "--example",
        help="Example usage",
        action="store_true",
    )

    parser.add_argument(
        "-e",
        "--energy",
        help="Single energy to calculate the ionization cross section for.",
        type=float,
    )

    parser.add_argument(
        "-el",
        "--energylower",
        help="Lower boundary of the range of energies to use.",
        type=float,
    )

    parser.add_argument(
        "-eu",
        "--energyupper",
        help="Upper boundary of the range of energies to use.",
        type=float,
    )

    parser.add_argument(
        "-ebin",
        "--energybinning",
        help="Energy binning to use for the range of energies.",
        type=float,
    )

    parser.add_argument(
        "-s",
        "--species",
        help="Particle species to use. Until now, hydrogen(both atomic and molecular), helium, neon and molecular nitrogen are implemented. \n"
             "Nomenclature: \n"
             "Hydrogen: h / h2 \n"
             "Helium: he \n"
             "Neon: ne \n"
             "Nitrogen n2 \n",
        type=str,
    )

    return parser.parse_args()

def init_calculation(method, Atom):

    if method == 'bed':
        return lambda T: CrossSectionCalcBed(T, atom = Atom)

    if method == 'beb':
        return lambda T: CrossSectionCalcBeb(T, atom = Atom)

    raise Exception('Calculation method has not been assigned yet. Work in progress,')


def example():
    Atom = AtomFactory.get_h2()
    log_boundary_a = 7.
    log_boundary_b = 13.
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
    plt.semilogx(x, y_bed, label = 'BED model')
    plt.semilogx(x, y_beb, label = 'BEB model')
    plt.semilogx(x, y_vac_note, label = 'Vacuum Note Calculation', linestyle = ':')
    plt.semilogx(x, y_bethe, label = 'Asymptotic behaviour (Bethe theory)')
    plt.legend()
    plt.show()


def calc_CS(energy,species):

    try:
        Atom = atom_library[species]
        cross = init_calculation(energy, Atom.calc_method)


    except Exception as err:
        print(err)

    return cross.calculate()



def main():

    args = parse_arguments()

    try:
        if args.version:
            print('0.1 alpha')

        if args.example:
            print('Example curve for molecular hydrogen in a region 1e7 to 1e13')
            example()

        if args.energy and not args.species:
            print('No species defined. Add one via -s or --species .')

        if args.species and not args.energy:
            print('No energy defined. Add one via -e or --energy .')

        if args.energy and args.species:
            return calc_CS(args.energy, args.species)

        if (args.energylower or args.energyupper) and not (args.energyupper and args.energylower):
            raise Exception('You are missing a boundary.')

        if (args.species and args.energylower and args.energyupper and args.energybinning):
            number_datapoints = (args.energyupper - args.energylower)//args.energybinning
            dataset = []
            for energy in np.linspace(args.energylower, args.energyupper, number_datapoints):
                dataset.append([energy, calc_CS(energy, args.species)])
            scs.save_as_mat(dataset, args.energylower, args.energyupper)

    except Exception as err:
        print(err)

if __name__ == "__main__":
    main()
