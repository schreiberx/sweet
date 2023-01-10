#! /usr/bin/env python3

import os
import sys

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d)
from EFloat import *
from Functions import *
sys.path.pop()


import math
import struct


class REXICoefficients:

    def __init__(self):

        self.efloat = None

        self.function_name = None

        self.alphas = []
        self.betas = []
        self.gamma = None

        self.unique_id_string = None
        
        self.normalized_steady_state = False
        return

    def len(self):
        return len(self.alphas)


    def toFloat(self):
        c = REXICoefficients()

        c.efloat = EFloat("float")
        c.function_name = self.function_name

        if self.efloat != None:
            if self.efloat.abs(self.efloat.im(self.gamma)) > 1e-10:
                print("WARNING")
                print("WARNING")
                print("WARNING")
                print("WARNING: Imaginary value "+str(self.efloat.im(self.gamma))+" should be close to zero")
                print("WARNING")
                print("WARNING")
                print("WARNING")

            c.gamma = float(self.efloat.re(self.gamma))
        else:
            c.gamma = float(self.gamma.real)

        c.alphas = [complex(a) for a in self.alphas]
        c.betas = [complex(b) for b in self.betas]

        c.unique_id_string = self.unique_id_string

        return c



    def writestr(self, f, string):
        f.write(string.encode('ascii'))


    def writeefloat(self, f, value):
        float_value = float(value)
        s = struct.pack('d', float_value)
        f.write(s)

    def writeefloatcplx(self, f, value, binary):
        if binary:
            self.writeefloat(f, self.efloat.re(value))
            self.writeefloat(f, self.efloat.im(value))
        else:
            self.writestr(f, self.efloat.floatToStr(self.efloat.re(value))+"\t")
            self.writestr(f, self.efloat.floatToStr(self.efloat.im(value))+"\n")
        


    def write_file(self, filename, binary=True):
        f = open(filename, "wb")

        self.writestr(f, "# N "+str(len(self.alphas))+"\n")

        if binary:
            self.writestr(f, "# binary 1\n")
        else:
            self.writestr(f, "# binary 0\n")

        self.writestr(f, "# function_name "+str(self.function_name)+"\n")
        if self.gamma != None:
            self.writestr(f, "# gamma\n")
            self.writeefloatcplx(f, self.gamma, binary)

        self.writestr(f, "# alphas\n")
        for a in self.alphas:
            self.writeefloatcplx(f, a, binary)

        self.writestr(f, "# betas\n")
        for b in self.betas:
            self.writeefloatcplx(f, b, binary)



    def symmetric_reduction(self):
        """
        Exploit a symmetry of the poles in case that they are complex conjugate symmetric
        """
        
        import numpy as np

        def print_coeffs():
            for i in range(len(self.alphas)):
                print(f"alpha: {self.alphas[i]}\t\tbeta: {self.betas[i]}")
            print("")

        print_coeffs()
        print("Number of coefficients before symmetric reduction: "+str(len(self.alphas)))

        for i in range(len(self.alphas)):
            if i >= len(self.alphas):
                break

            # Search for corresponding conjugate one
            a1 = self.alphas[i]
            for j in range(i+1, len(self.alphas)):
                if j >= len(self.alphas):
                    break

                a2 = self.alphas[j]
                if np.isclose(np.real(a1), np.real(a2)) and np.isclose(np.imag(a1), -np.imag(a2)):
                    del self.alphas[j]
                    self.betas[i] *= 2.0

        print_coeffs()
        print("Number of coefficients after symmetric reduction: "+str(len(self.alphas)))
        
        self.symmetric_reduction_applied = True
        if not '_symred' in self.unique_id_string:
                self.unique_id_string += '_symred'
        
        print("*"*80)
        print("TODO: Not yet tested!")
        print("*"*80)


    def beta_filter(
        self,
        beta_threshold
    ):
        """
        Filter out beta coefficients which amplitude is below a given threshold
        """

        print("WARNING: THIS IS JUST A DUMMY IMPLEMENTATION (beta filter)")

        new_alphas = []
        new_betas = []

        self.symmetric_reduction_applied = True
        if not '_betafilter' in self.unique_id_string:
                self.unique_id_string += '_betafilter'

        return


    def normalize_steady_state(
        self
    ):
        """
        Rescale beta coefficients so that they the solution converges to 1 for dt -> 0
        """
        
        if self.normalized_steady_state:
            raise Exception("Normalization already triggered")
        
        if self.gamma != 0:
            raise Exception("Normalization is not supported for beta != 0")
        
        if self.function_name == None:
            raise Exception("Function name required for normalization")
            
        
        x = 0
        val = 0
        for i in range(len(self.alphas)):
            val += self.betas[i] / (x - self.alphas[i])

        fun = Functions(self.function_name, efloat_mode='float')     
        
        target_value = fun.eval(0)
        
        rescale = target_value/val
        
        self.betas = [i*rescale for i in self.betas]
        
        self.unique_id_string += "_nrm"
        

    def eval(self, x):
        retval = self.gamma

        for i in range(len(self.alphas)):
            retval += self.betas[i] / (x - self.alphas[i])

        return retval


if __name__ == "__main__":

    r = REXICoefficients()
    r.efloat = EFloat()
    r.gamma = 1e-12
    r.alphas.append(1j+2.0)
    r.betas.append(3j+4.0)

    r.write_file("/tmp/test.txt")

    r2 = r.toFloat()
