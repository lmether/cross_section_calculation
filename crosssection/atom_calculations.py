
import numpy as np
#import sympy



def singular_initial_states(Ne):
    in_state0 = np.zeros(Ne)
    for electron in range(Ne):
        in_state0[electron] = 3
    return in_state0



def calc_B(Ne):
    # TODO
    return np.array([21.7, 48.47, 866.9])

def calc_N(Ne):
    # TODO
    return Ne//3



def calc_Z(Ne):
    # TODO
    return Ne


def calc_Mi(Ne):
    # TODO
    return Ne


def calc_Ni(Ne):
    # TODO
    return Ne


def calc_(Ne):
    # TODO
    return Ne