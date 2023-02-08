#! /usr/bin/env python3

#
# Author: Martin Schreiber
# Email: schreiberx@gmail.com
# Date: 2017-06-17
#


import sys
import math
import cmath
import os

d = os.path.dirname(os.path.realpath(__file__))

sys.path.append(d+'/..')
import EFloat as ef
from trexi.TREXI_GaussianCoefficients import *
from trexi.TREXI_GaussianPhi0 import *

sys.path.pop()


class TREXI_GaussianPhi0:

    def __init__(
            self,
                        gaussphi0_N,    # required argument
                        gaussphi0_basis_function_spacing,    # required argument

            floatmode = None
    ):
        self.efloat = ef.EFloat(floatmode)

        self.h = self.efloat.to(gaussphi0_basis_function_spacing)
        self.M = int(gaussphi0_N)

        self.b = [self.efloat.exp(self.h*self.h)*self.efloat.exp(-1j*(float(m)*self.h)) for m in range(-self.M, self.M+1)]

        # Generate dummy Gaussian function
        fafcoeffs = TREXI_GaussianCoefficients()
        fafcoeffs.function_name = 'gaussianorig'
        fafcoeffs.function_scaling = self.efloat.to(1.0)


    def output(self):
        for i in self.b:
            print(i)


    def fun(
        self,
        i_x
    ):
        return self.efloat.exp(self.efloat.i*i_x)


