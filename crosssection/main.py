import matplotlib.pyplot as plt
import numpy as np


import argparse

from crosssection.cross_section import AtomFactory, CrossSectionCalcBed, CrossSectionCalcBeb, CrossSectionCalcVacNote


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Ich bin eins dumma."
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
        type=int,
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

def main():
    args = parse_arguments()
    if args.version:
        print('0.1 alpha')
    if args.energy:
        print(args.energy)
    if args.example:
        example()

if __name__ == "__main__":
    main()
