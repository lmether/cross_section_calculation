import matplotlib.pyplot as plt
import numpy as np


import argparse

from crosssection.cross_section import AtomFactory, CrossSectionCalcBed, CrossSectionCalcBeb, CrossSectionCalcVacNote

atom_collection = {'h': AtomFactory.get_hydrogen(), 'h2': AtomFactory.get_h2(), 'he': AtomFactory.get_helium(),
                   'ne': AtomFactory.get_nitrogen(), 'n2': AtomFactory.get_nitrogen()}

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
        help="Energy to use.",
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

    pass


def main():
    args = parse_arguments()
    if args.version:
        print('0.1 alpha')
    if args.energy and args.species:
        calc_CS(args.energy, args.species)
    if args.energy:
        print(args.energy)
    if args.example:
        print('Example curve for molecular hydrogen in a region 1e7 to 1e13')
        example()

if __name__ == "__main__":
    main()
