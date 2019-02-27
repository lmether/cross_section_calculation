import scipy.io as sio
import os

def save_as_mat(array, energy1 = '', energy2 = ''):

    cross_sec_dict = {'Ionization cross section': array}
    path = os.getcwd()
    sio.savemat(path + '/cross_sections' + energy1 + '_' + energy2, mdict = cross_sec_dict)

    return None