#! /usr/bin/env python3

#
# Author: Martin Schreiber
# Email: schreiberx@gmail.com
# Date: 2017-06-16
#

import sys
import math
import numpy
import cmath
import traceback
import os

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+'/..')

from REXICoefficients import *
from trexi.TREXI_GaussianPhi0 import *
from trexi.TREXI_GaussianCoefficients import *

sys.path.pop()


#
# This class offers the coefficients for the
#
# [T ] erry Haut et al. based
#
# [R ] ational approximation of an
# [Ex] ponential
# [I ] ntegrator
#
class TREXI:
    
    def __init__(
        self,
        efloat_mode = None,
    ):
        self.efloat = ef.EFloat(efloat_mode)



    def setup(
        self,
        # gaussian approximation of phi function

        M = None,        # Number of basis functions
        h = None,        # spacing of basisfunctions

        test_min = None,    # Test range
        test_max = None,

        function_name = "phi0",    # function to approximate

        # parameters which are specific to basis_function_name = "gaussianorig"
        ratgauss_N = None,
        ratgauss_basis_function_spacing = 1.0,
        ratgauss_basis_function_rat_shift = None,

        reduce_to_half = False,

        floatmode = None
    ):
        self.M = M
        self.h = h

        self.test_min = test_min
        self.test_max = test_max

        self.function_name = function_name

        self.ratgauss_N = ratgauss_N
        self.ratgauss_basis_function_spacing = ratgauss_basis_function_spacing
        self.ratgauss_basis_function_rat_shift = ratgauss_basis_function_rat_shift

        self.reduce_to_half = reduce_to_half

        if self.function_name == "phi0":
            self.setup_rexi_gaussianorig_phi0(floatmode = floatmode)

            # setup range of approximation interval
            self.test_max = self.h*self.M-self.efloat.pi
            self.test_min = -self.test_max

            #    
            # \return \f$ cos(x) + i*sin(x) \f$
            #
            def eval_fun(i_x):
                return self.efloat.re(self.efloat.exp(self.efloat.i*i_x))
            self.eval_fun = eval_fun

            def eval_exp(i_x, i_u0):
                return self.efloat.re(self.efloat.exp(self.efloat.i*i_x)*i_u0)
            self.eval_exp = eval_exp


            #    
            # \return \f$ cos(x) + i*sin(x) \f$
            #
            def eval_fun_cplx(i_x):
                return self.efloat.exp(self.efloat.i*i_x)
            self.eval_fun_cplx = eval_fun_cplx

            def eval_exp_cplx(i_x, i_u0):
                return self.efloat.exp(self.efloat.i*i_x)*i_u0
            self.eval_exp_cplx = eval_exp_cplx

        else:
            raise Exception("TODO")


        self.unique_id_string = ""
        self.unique_id_string += "_"+function_name
        self.unique_id_string += "_"+str(M)
        self.unique_id_string += "_"+str(h)


        coeffs = REXICoefficients()
        coeffs.function_name = function_name
        coeffs.efloat = self.efloat
        coeffs.alphas = self.alpha
        coeffs.betas = self.beta
        coeffs.gamma = 0

        coeffs.unique_id_string = self.getUniqueId()

        return coeffs



    #
    # Setup phi0(x) = Re(exp(ix))
    #
    # The coefficients are then made available in
    # self.alpha_reim
    # self.beta_re
    # self.beta_im
    #
    def setup_rexi_gaussianorig_phi0(self, floatmode = None):

        self.ratgauss_coeffs = TREXI_GaussianCoefficients(floatmode = floatmode)

        #
        # Step 1) Load rational approximation of gaussian basis functions
        #
        self.ratgauss_coeffs.load_orig_ratgaussian_poles()

        # load parameters from faf coefficients
        self.ratgauss_N = self.ratgauss_coeffs.N
        self.ratgauss_basis_function_spacing = self.ratgauss_coeffs.basis_function_spacing
        self.ratgauss_basis_function_rat_shift = self.ratgauss_coeffs.basis_function_rat_shift


        #
        # Step 2) Load gaussian approximation of phi0
        #
        self.phi0 = TREXI_GaussianPhi0(
            gaussphi0_N = self.M,
            gaussphi0_basis_function_spacing = self.h
        )

        #
        # Step 3) Merge
        #

        # Use variable naming from original REXI paper
        L = self.ratgauss_coeffs.N//2
        h = self.h
        M = self.M
        N = M+L

        self.alpha_reim = [0 for i in range(2*N+1)]
        self.beta_re = [0 for i in range(2*N+1)]
        self.beta_im = [0 for i in range(2*N+1)]

        for l in range(-L, L+1):
            for m in range(-M, M+1):
                n = l+m
                w = self.ratgauss_coeffs.weights_cplx[l+L]
                # NOTE: We use the conjugate here!
                w = self.efloat.conj(w)

                self.alpha_reim[n+N] = h*(self.ratgauss_coeffs.basis_function_rat_shift + self.efloat.i*n)

                self.beta_re[n+N] += self.efloat.re(self.phi0.b[m+M])*h*w
                self.beta_im[n+N] += self.efloat.im(self.phi0.b[m+M])*h*w

        #
        # Merge real and imaginary approximation together
        #
        alpha_new = []
        beta_new = []
        M = len(self.alpha_reim)
        for i in range(M):
            beta_new.append(0.5*(self.beta_re[i] + 1.j*self.beta_im[i]))
            alpha_new.append(self.alpha_reim[i])

        for i in range(M):
            beta_new.append(-0.5*(self.beta_re[i] - 1.j*self.beta_im[i]))
            alpha_new.append(-self.alpha_reim[i])

        self.alpha = [-a for a in alpha_new]
        self.beta = beta_new

        if self.reduce_to_half:
            # reduce the computational amount to its half,
            # see understanding REXI in the documentation folder

            alpha_new = []
            beta_new = []

            N = len(self.alpha)//2
            alpha_new = self.alpha[0:N//2+1] + self.alpha[N:N+N//2+1]
            beta_new = self.beta[0:N//2+1] + self.beta[N:N+N//2+1]

            for i in range(N//2):
                beta_new[i] *= 2.0
                beta_new[N//2+1+i] *= 2.0

            self.alpha = alpha_new
            self.beta = beta_new


    def getUniqueId(self):
        return "TREXI"+self.unique_id_string



    #
    # approx with linear operator, input: complex value, output: complex value
    #
    def approx_fun_cplx_linop(self, i_x, i_u0):
        # Use linear system to avoid complex value (halving doesn't work with this)
        retval = numpy.array([0, 0], dtype=complex)
        u0 = numpy.array([i_u0.real, i_u0.imag], dtype=complex)

        if i_x.imag != 0.0:
            raise Exception("Imaginary value for i_x")

        # Convert to linear solver
        L = numpy.array([[0, -i_x], [i_x, 0]], dtype=complex)

        N = len(self.alpha)
        for n in range(N):
            retval += self.beta[n]*numpy.linalg.solve(L - numpy.eye(2, dtype=complex)*self.alpha[n], u0)

        return retval[0].real + retval[1].real*1.j



    #
    # approx
    #
    def approx_fun_cplx(self, i_x):

        retval = 0
        # Split computation into real part of \f$ cos(x) \f$ and imaginary part \f$ sin(x) \f$
        for n in range(len(self.alpha)):
            denom = (self.efloat.i*i_x - self.alpha[n])
            retval += (self.beta[n] / denom)

        return retval



    #
    # approx
    #
    def approx_fun_cplx(self, i_x, i_u0):

        retval = 0
        # Split computation into real part of \f$ cos(x) \f$ and imaginary part \f$ sin(x) \f$
        for n in range(len(self.alpha)):
            denom = (self.efloat.i*i_x - self.alpha[n])
            retval += (self.beta[n] / denom * i_u0)

        return retval



    def runTests(self):

        h = self.h
        N = self.M

        maxerror = 0
        d = int((self.test_max-self.test_min)*22)

        for u0 in [1.0, (1.0 + 2.0j)/numpy.sqrt(5.0)]:

            print("+"*80)
            print("+"*80)
            print("N: "+str(len(self.alpha)))
            print("u0: "+str(u0))
            print("+"*80)
            print("+"*80)

            if True:
                print("*"*40)
                print(">>> ODE TEST (complex) <<<")

                maxerror = 0
                for x in self.efloat.linspace(self.test_min, self.test_max, d):
                    #a = self.eval_fun_cplx(x)
                    #b = self.approx_fun_cplx(x)

                    a = self.eval_exp_cplx(x, u0)
                    b = self.approx_fun_cplx(x, u0)

                    e = abs(a-b)
                    maxerror = max(maxerror, e)

                print("max error: "+str(maxerror))

                if maxerror > 1e-9:
                    if self.reduce_to_half:
                        print("\t\tERROR ignored (reduce_to_half), real-only violated")
                    else:
                        raise Exception("Error threshold exceeded")
                else:
                    print("\t\tOK")


            if True:
                print("*"*40)
                print(">>> PDE TEST (real-valued PDE) <<<")

                maxerror = 0
                for x in self.efloat.linspace(self.test_min, self.test_max, d):
                    a = self.eval_exp_cplx(x, u0)
                    b = self.approx_fun_cplx_linop(x, u0)

                    e = abs(a-b)
                    maxerror = max(maxerror, e)

                print("max cplx error: "+str(maxerror))

                if maxerror > 1e-9:
                    if self.reduce_to_half and not self.merge:
                        print("\t\tERROR ignored (reduce to half)")
                    else:
                        raise Exception("Error threshold exceeded")
                else:
                    print("\t\tOK")


def main():
    # double precision
    floatmode = 'float'

    # multiprecision floating point numbers
    #floatmode = 'mpfloat'


    # load float handling library
    efloat = ef.EFloat(floatmode)

    h = efloat.to(0.2)
    M = 64
    print("M: "+str(M))
    print("h: "+str(h))
    print("")

    for half in [False, True]:

        print("*"*80)
        print("* Original REXI coefficients")

        rexi = TREXI()
        rexi.setup(
            M = M,
            h = h,

            function_name = "phi0",

            reduce_to_half = half,
            floatmode = floatmode,
        )

        print("")
        print("*"*80)
        print("Running tests...")
        rexi.runTests()

        print("*"*80)
        print("* REXI coefficients based on original REXI mu")
        print("*"*80)

        rexi = TREXI()
        rexi.setup(
            M = M,
            h = h,

            ratgauss_basis_function_spacing = 1.0,
            ratgauss_basis_function_rat_shift = efloat.to('-4.31532151087502402475593044073320925235748291015625'),

            reduce_to_half = half,

            floatmode = floatmode,
        )
        print("*"*80)
        print("Rational approx. of Gaussian data:")
        #rexi.ratgauss_coeffs.print_coefficients()

        print("")
        print("*"*80)
        print("Running tests...")
        rexi.runTests()




if __name__ == "__main__":

    try:
        main()

    except SystemExit as e:
        sys.exit(e)

    except Exception as e:
        print("* ERROR:")
        print(str(e))
        traceback.print_exc()
        sys.exit(1)

