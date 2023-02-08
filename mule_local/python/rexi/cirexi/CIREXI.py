#! /usr/bin/env python3
#
# Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
# Date: 2017-08-16
#


import math
import cmath
import numpy as np
import sys

from mule.rexi.EFloat import *
from mule.rexi.Functions import *
from mule.rexi.REXICoefficients import *



class CIREXI:
    """
    Cauchy Contour Integration REXI method
    """

    #
    # Constructor
    # See setup(...) for documentation on parameters
    # 
    def __init__(
        self,
        efloat_mode = "float"
    ):
        self.efloat_mode = efloat_mode
        self.efloat = EFloat(efloat_mode)

        self.alphas = []
        self.betas = []
        self.unique_id_string = ""



    def setup(
        self,
        function_name,
        N: int,
        R: float = None,
        lambda_shift: complex = None,

        lambda_max_real: float = None,
        lambda_include_imag: float = None,
        half_shifted: bool = True
    ):
        self.alphas = []
        self.betas = []
        self.function_name = ""

        if lambda_shift != None:
            self.setup_shifted_circle(function_name, N=N, R=R, lambda_shift=lambda_shift, half_shifted=half_shifted)

        elif lambda_max_real != None:
            self.setup_limited_circle(function_name, N=N, lambda_max_real=lambda_max_real, lambda_include_imag=lambda_include_imag, half_shifted=half_shifted)

        elif lambda_include_imag != None:
            self.setup_circle(function_name, N=N, R=lambda_include_imag, half_shifted=half_shifted)

        elif R != None:
            self.setup_shifted_circle(function_name, N=N, R=R, lambda_shift=0)

        else:
            raise Exception("Can't calculate circle radius, please provide imag/real limits.")

        coeffs = REXICoefficients()
        coeffs.alphas = self.alphas[:]
        coeffs.betas = self.betas[:]
        coeffs.gamma = 0
        coeffs.function_name = self.function_name

        coeffs.unique_id_string = self.getUniqueId()

        return coeffs



    def setup_shifted_circle(
        self,
        function_name,
        N: int,
        R: float,
        lambda_shift: complex,
        half_shifted: bool = True
    ):
        """
        Setup coefficients for Cauchy contour integral coefficients
        using a shifted circle with the shift given by lambda_shift
        which fulfills the constraints given by lambda_*

        Args:
            N (int):
                Number of points on Cauchy contour.

            R (float):
                Radius of circle for Cauchy contour

            lambda_shift (complex):
    `            Shift of center of circle on complex plane
        """

        self.function_name = function_name
        self.lambda_shift = lambda_shift

        self.fun = Functions(function_name, efloat_mode = self.efloat_mode)

        for j in range(N):
            if half_shifted:
                theta_j = self.efloat.pi2*(j+self.efloat.to(0.5))/N
            else:
                theta_j = self.efloat.pi2*j/N

            # sampling position of support point
            pos = R*self.efloat.exp(self.efloat.i*theta_j)

            # shifted position
            alpha = pos + lambda_shift

            self.betas.append(-self.fun.eval(alpha)*pos / N)
            self.alphas.append(alpha)


        self.unique_id_string = "shic"
        self.unique_id_string += "_"+function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(self.efloat.re(lambda_shift))
        self.unique_id_string += "_"+str(self.efloat.im(lambda_shift))



    def setup_circle(
        self,
        function_name,
        N: int,
        R: float,
        half_shifted: bool = True
    ):
        """
        Setup circle contour integral centered at origin

        Args:
            N (int):
                Number of quadrature poles

            R (float):
                Radius of circle
        """

        self.setup_shifted_circle(function_name, N, R, 0.0, half_shifted)

        self.unique_id_string = ""
        #self.unique_id_string = "circ"
        self.unique_id_string += function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(R)



    def setup_limited_circle(
        self,
        function_name: str,
        N: int,
        lambda_max_real: float,
        lambda_include_imag: float,
        half_shifted: bool = True
    ):
        """
        Setup coefficients for Cauchy contour integral coefficients circle
        which fulfills the constraints given by lambda_*

        Args:
            N (int):
                Number of points on Cauchy contour.

            lambda_max_real (float):
                Maximum allowed real number on contour.

            lambda_include_imag (float):
                Include at least this imaginary value as part of the contour.
                Even, if the contour has to be enlarged
        """

        # If the maximal real number is larger than the max to be included imaginary number, set a lower value for the max real number
        if lambda_max_real > lambda_include_imag:
            lambda_max_real = lambda_include_imag

        x0 = lambda_max_real
        xm = lambda_include_imag
        r = (x0*x0 + xm*xm)/(2.0*x0)
        center = lambda_max_real-r

        self.setup_shifted_circle(function_name, N, r, center, half_shifted=half_shifted)

        self.unique_id_string = ""
        #self.unique_id_string = "limc"
        self.unique_id_string += function_name
        self.unique_id_string += "_"+str(N)
        self.unique_id_string += "_"+str(lambda_max_real)
        self.unique_id_string += "_"+str(lambda_include_imag)
        self.unique_id_string += "_h"+str(int(half_shifted))



    def getUniqueId(self):
        return "CIREXI_"+self.unique_id_string
