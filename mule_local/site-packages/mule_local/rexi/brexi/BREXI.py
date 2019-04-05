#! /usr/bin/env python3

"""
Copyright 2017, Jed Brown

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
Author: Jed Brown

Other contributions by: Matthew Normile, Martin Schreiber
"""


import numpy
import mule_local.rexi.brexi.rkanalysis as rkanalysis

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
        quadrature_method = "gauss",
        tanhsinh_h = None
    ):
        self.N = N
        self.quadrature_method = quadrature_method

        if quadrature_method == "gauss":
            A, b, c = rkanalysis.gauss(N)
            self.quadrature_method_short = "G"

        elif quadrature_method == "chebyshev":
            A, b, c = rkanalysis.chebyshev(N)
            self.quadrature_method_short = "C"

        elif quadrature_method == "radau":
            A, b, c = rkanalysis.radau(N)
            self.quadrature_method_short = "R"

        elif quadrature_method == "tanhsinh":
            A, b, c = rkanalysis.tanhsinh(N, h=tanhsinh_h)
            self.quadrature_method_short = "T"

        else:
            raise Exception("Quadrature method '"+quadrature_method+"' not supported")

        # Create modified/unified REXI form
        A, b, c = rkanalysis.butcher_diag(A, b)
        alphas, betas = rkanalysis.butcher2rexi(A, b)
        gamma = self.efloat.to(1.0) - numpy.sum(betas/alphas)

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


if __name__ == "__main__" and False:

    b = BREXI()

    # TANHSINH
    # 16: h = 0.1
    # 32: h = ???
    # 64: h = ???
    # 128: h = ???

    filepre = 'brexi_tanhsinh'

    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import matplotlib.pyplot as plt
    plt.clf()

    N=32

    for h in [0.01*1.3**i for i in range(10)]:
        print("*"*80)
        print("h="+str(h))

        #try:
        if True:
            brexi = BREXI()
            coeffs = brexi.setup(N=N, quadrature_method="tanhsinh", tanhsinh_h=h)

        #except Exception as e:
        #    print(e)
        #    continue

        # Convert to floating point
        coeffs = coeffs.toFloat()

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

        plt.plot(xpts, yvals, label="h="+str(h))

    plt.legend()
    plt.xlim(test_range[0], test_range[1])
    plt.ylim(1e-12, 10)
    plt.yscale("log")

    plt.savefig('output_'+filepre+'_error_oscillatory.pdf')

    sys.exit(1)


if __name__ == "__main__":

    b = BREXI()

    #coeffs = b.setup(16, "gauss")

    # TANHSINH
    # 16: h = 0.1
    # 32: h = ???
    # 64: h = ???
    # 128: h = ???

    brexi = BREXI()
    N=32
    h=0.04
    coeffs = brexi.setup(N=N, quadrature_method="tanhsinh", tanhsinh_h=h)

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

        filepre = 'brexi_tanhsinh'

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
