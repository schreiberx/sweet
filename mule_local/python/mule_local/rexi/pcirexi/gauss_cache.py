from os import path
from sympy.integrals.quadrature import gauss_legendre as g_le, gauss_lobatto as g_lo
import pickle

from mule_local.rexi.pcirexi.section.arbitrary_spline_basis_generator import generate_basis_polynomials


class GaussCache:
    lobatto_dict: dict
    legendre_dict: dict
    spline_dict: dict
    file_path: str

    def __init__(self, file_path='gauss_cache.pkl'):
        self.file_path = file_path
        if path.isfile(file_path):
            file = open(self.file_path, "rb")
            pickle_load = pickle.load(file)
            if len(pickle_load) == 2:
                self.spline_dict = {}
                self.legendre_dict, self.lobatto_dict = pickle_load
            else:
                self.legendre_dict, self.lobatto_dict, self.spline_dict = pickle_load
            file.close()
        else:
            self.legendre_dict = {}
            self.lobatto_dict = {}
            self.spline_dict = {}

    def gauss_legendre(self, points, accuracy):
        if (points, accuracy) in self.legendre_dict:
            return self.legendre_dict[(points, accuracy)]
        else:
            print("run " + str(points) + " ac: " + str(accuracy))
            value = g_le(points, accuracy)
            self.legendre_dict[(points, accuracy)] = value
            self.pickle_dicts()
            return value

    def gauss_lobatto(self, points, accuracy):
        if (points, accuracy) in self.lobatto_dict:
            return self.lobatto_dict[(points, accuracy)]
        else:
            print("run " + str(points) + " ac: " + str(accuracy))
            value = g_lo(points, accuracy)
            self.lobatto_dict[(points, accuracy)] = value
            self.pickle_dicts()
            return value

    def spline_basis_polynomials(self, samples, n, precision):
        tuple = (samples, n, precision)
        if tuple in self.spline_dict:
            return self.spline_dict[tuple]
        else:
            print("run ")
            value = generate_basis_polynomials(samples, n, precision)
            self.spline_dict[tuple] = value
            self.pickle_dicts()
            return value

    def pickle_dicts(self):
        file = open(self.file_path, "wb")
        pickle.dump((self.legendre_dict, self.lobatto_dict, self.spline_dict), file)
        file.close()
