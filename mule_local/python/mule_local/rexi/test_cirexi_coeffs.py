#! /usr/bin/env python3

import sys
import numpy as np

from Functions import *
from trexi.TREXI import *
from cirexi.CIREXI import *
from brexi.BREXI import *
from elrexi.ELREXI import *


for rexi_method in [
        "cirexi",
    ]:

    for function_name in ["phi0"]: #, "ups1", "ups2", "ups3"]:

        # CI-REXI: Number of quadrature poles
        #N = 256
        N = 32

        # CI-REXI: Value on imaginary axis to be included
        #same used for ELREXI as max imag of ellipse
        lambda_include_imag = 20

        # CI-REXI: Maximum value of quadrature pole
        #same used for ELREXI as max real for ellipse
        lambda_max_real = 10


        # Butcher-REXI
        butcherOrder = 10


        # Testing: number of samples
        num_test_samples = 12345

        # Testing: Range (start, end)
        test_range = [None, None]

        # Error to test for
        error_eps = 1e-8

        # efloat_mode
        efloat_mode = "float"
        #efloat_mode = "mpfloat"


        cirexi = CIREXI(efloat_mode=efloat_mode)

        coeffs = cirexi.setup(
            function_name = function_name,
            N = N,
            lambda_max_real = lambda_max_real,
            lambda_include_imag = lambda_include_imag
        )

        # Convert to floating point
        coeffs = coeffs.toFloat()

        unique_id_string = cirexi.getUniqueId()

        #test_range = [-lambda_include_imag*0.5, lambda_include_imag*0.5]
        test_range = [-1.0, 1.0]


        function = Functions(
            function_name = function_name,
            efloat_mode = "float"
        )


        print("")
        print("COEFF_NR\tALPHAS\t\t\tBETAS")
        for i in range(len(coeffs.alphas)):
            print(str(i)+"\t"+str(coeffs.alphas[i])+"\t"+str(coeffs.betas[i]))

        max_error = 0
        verbosity = 0
        for x in np.linspace(test_range[0], test_range[1], num_test_samples):
            lam = 1j*x

            y = function.eval(lam)
            yn = coeffs.eval(lam)

            err = np.abs(y-yn)

            if verbosity > 0:
                #if True:
                if False:
                    print("x="+str(lam)+"\t\terror="+str(err))
                else:
                    print("Lambda: "+str(lam))
                    print(" +  exact: "+str(y))
                    print(" + approx: "+str(yn))
                    print(" + Error: "+str(err))
                    print("")

            max_error = max(max_error, err)


        if verbosity == 0:
            print(" + test_range: ["+str(test_range[0])+", "+str(test_range[1])+"]")
            print(" + Error: "+str(max_error))

        if max_error > error_eps:
            raise Exception("Error threshold "+str(error_eps)+" exceeded")

