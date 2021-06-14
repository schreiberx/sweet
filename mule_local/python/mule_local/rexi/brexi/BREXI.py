#! /usr/bin/env python3


import sys
import numpy
import mule_local.rexi.brexi.rk_co as rk_co

import mule_local.rexi.EFloat as ef
from mule_local.rexi.REXICoefficients import *
from mule_local.rexi.Functions import *




class BREXI:

    def __init__(
        self,
        efloat_mode = None,
    ):
        self.efloat = ef.EFloat(efloat_mode)


    def setup(
        self,
        N,
        quadrature_method = "gauss_legendre",
        tanhsinh_h = None
    ):
        self.N = N
        self.quadrature_method = quadrature_method
        self.quadrature_method_short = quadrature_method[0].upper()


        co = rk_co.rk_co(N, quadrature_method)
        
        """
        Step 1) Diagonalization
        """

        """
        Compute Eigendecomposition with Eigenvalues in D
        """
        D, E = numpy.linalg.eig(co.A)

        """
        Compute diagonal-only W so that
            W * E^-1 * ONES = ONES
        with ONES a vector with only ones
        """
        #
        # E^-1 * ONES = eta
        # <=> E*eta = ONE
        # <=> eta = solve(E, ONE)
        #
        eta = numpy.linalg.solve(E, numpy.ones(N))

        #
        # W * eta = ONES
        # diag-only W
        # diag(W) = eta^-1
        # diag(W_inv) = eta
        #
        W_inv = numpy.diag(eta)
        b_tilde = co.b.T.dot(E.dot(W_inv))

        # Create modified/unified REXI form
        """
        Step 2) Get REXI formulation
        """

        gamma = 1.0 - numpy.sum(b_tilde/D)
        alphas = 1/D
        betas = -b_tilde/(D*D)

        coeffs = REXICoefficients()
        coeffs.function_name = "phi0"
        coeffs.efloat = self.efloat
        coeffs.alphas = alphas
        coeffs.betas = betas
        coeffs.gamma = gamma

        coeffs.unique_id_string = self.getUniqueId()

        return coeffs



    def getUniqueId(self):
        return "BREXI_"+self.quadrature_method_short+"_"+str(self.N)



if __name__ == "__main__":
    numpy.set_printoptions(precision=20)

    for method in ['gauss_legendre', 'gauss_chebyshev_u']:

        brexi = BREXI()

        N=16

        coeffs = brexi.setup(N=N, quadrature_method=method)
        filepre = 'brexi_'+method


        # Convert to floating point
        coeffs = coeffs.toFloat()


        print("Alphas:")
        for i in coeffs.alphas:
            print(" + "+str(i))
        print("")

        print("Betas:")
        for i in coeffs.betas:
            print(" + "+str(i))
        print("")

        print("Gamma:")
        print(coeffs.gamma)
        print("")


        if True:
        #if False:
            import matplotlib
            matplotlib.use('Agg')
            import numpy as np
            import matplotlib.pyplot as plt

            x = np.real(coeffs.alphas)
            y = np.imag(coeffs.alphas)

            plt.clf()
            plt.scatter(x, y)
            plt.savefig('output_'+filepre+'_alphas.pdf')

            plt.clf()
            plt.scatter(x, y, s=np.log(np.absolute(coeffs.betas))+1.0)
            plt.savefig('output_'+filepre+'_alphas_scaled.pdf')

            x = np.real(coeffs.betas)
            y = np.imag(coeffs.betas)

            plt.clf()
            plt.scatter(x, y)
            plt.savefig('output_'+filepre+'_betas.pdf')


            #
            # Error plot
            #
            plt.clf()
            function = Functions(
                function_name = "phi0",
                efloat_mode = "float"
            )

            test_range = [-N, N]
            num_test_samples = 4096

            max_error = 0
            xpts = np.linspace(test_range[0], test_range[1], num_test_samples)
            yvals = np.zeros(num_test_samples)

            for i in range(num_test_samples):
                x = xpts[i]
                lam = 1j*x

                y = function.eval(lam)
                yn = coeffs.eval(lam)

                err = np.abs(y-yn)

                yvals[i] = err

            plt.plot(xpts, yvals, 'r-')

            plt.xlim(test_range[0], test_range[1])
            plt.ylim(1e-12, 10)
            plt.yscale("log")

            plt.savefig('output_'+filepre+'_error_oscillatory.pdf')
