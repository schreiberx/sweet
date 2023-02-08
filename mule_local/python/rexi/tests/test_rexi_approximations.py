#! /usr/bin/env python3

import sys
import os
import numpy as np

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+"/..")

from Functions import *
from trexi.TREXI import *
from cirexi.CIREXI import *
from brexi.BREXI import *
from elrexi.ELREXI import *

sys.path.pop()

for rexi_method in [
        #"trexi",
        #"cirexi",
        #"elrexi",
        #"brexi_tanhsinh",
        "brexi_gauss",
        "brexi_chebyshev",
        "brexi_jacobi",
    ]:

    function = Functions(
        function_name = "phi3",
        efloat_mode = "mpfloat"
    )

    """
    x = 1e-16j
    x = 0.6j+0.5
    for i in range(6):
        val = function.phiNRec(i, x)
        print(str(float(val.real))+"\t"+str(float(val.imag)))
        val = function.phiNSeries(i, x)
        print(str(float(val.real))+"\t"+str(float(val.imag)))
        print("")

    sys.exit(1)
    """

    for function_name in ["phi0", "phi1", "phi2", "phi3"]: #, "ups1", "ups2", "ups3"]:
    #for function_name in ["phi0"]:

        # T-REXI: Number of Gaussian basis functions
        M = 256

        # T-REXI: Spacing between Gaussian basis functions
        h = 0.2


        # CI-REXI: Number of quadrature poles
        #N = M
        N = 256

        # CI-REXI: Radius
        #R = M*h

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
        test_range_real = [None, None]
        test_range_imag = [None, None]

        # Error to test for
        error_eps = 1e-8

        # verbose
        verbosity = 10
        verbosity = 0

        # efloat_mode
        efloat_mode = "float"
        #efloat_mode = "mpfloat"


        coeffs = None

        print("REXI method: "+rexi_method)

        if rexi_method == "trexi":

            if function_name != "phi0":
                continue

            trexi = TREXI(efloat_mode=efloat_mode)
            coeffs = trexi.setup(
                function_name = function_name,
                M = N,
                h = h
            )

            # Convert to floating point
            coeffs = coeffs.toFloat()
            unique_id_string = trexi.getUniqueId()

            test_range_real = None
            test_range_imag = [-N*h*0.95, N*h*0.95]


        elif rexi_method == "cirexi":

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
            test_range_real = [-1.0, 1.0]
            test_range_imag = [-1.0, 1.0]

        elif rexi_method == "elrexi":

            elrexi = ELREXI(efloat_mode=efloat_mode)

            coeffs = elrexi.setup(
                function_name = function_name,
                N = N,
                lambda_max_real = lambda_max_real,
                lambda_max_imag = lambda_include_imag
            )

            # Convert to floating point
            coeffs = coeffs.toFloat()

            unique_id_string = elrexi.getUniqueId()

            #test_range = [-lambda_include_imag*0.5, lambda_include_imag*0.5]
            test_range_real = [-1.0, 1.0]
            test_range_imag = [-1.0, 1.0]

        elif rexi_method == "brexi_gauss":

            if function_name != "phi0":
                continue

            brexi = BREXI(efloat_mode=efloat_mode)
            coeffs = brexi.setup(N=butcherOrder, quadrature_method="gauss_legendre")

            # Convert to floating point
            coeffs = coeffs.toFloat()
            unique_id_string = brexi.getUniqueId()
            test_range_real = [butcherOrder*0.25, butcherOrder*0.25]
            test_range_imag = [butcherOrder*0.5, butcherOrder*0.5]



        elif rexi_method == "brexi_jacobi":

            if function_name != "phi0":
                continue

            brexi = BREXI(efloat_mode=efloat_mode)
            coeffs = brexi.setup(N=butcherOrder, quadrature_method="gauss_jacobi")

            # Convert to floating point
            coeffs = coeffs.toFloat()
            unique_id_string = brexi.getUniqueId()
            test_range_real = [butcherOrder*0.25, butcherOrder*0.25]
            test_range_imag = [butcherOrder*0.5, butcherOrder*0.5]


        elif rexi_method == "brexi_chebyshev":

            if function_name != "phi0":
                continue

            brexi = BREXI(efloat_mode=efloat_mode)
            coeffs = brexi.setup(N=butcherOrder, quadrature_method="gauss_chebyshev_u")

            # Convert to floating point
            coeffs = coeffs.toFloat()
            unique_id_string = brexi.getUniqueId()
            test_range_real = [butcherOrder*0.15, butcherOrder*0.15]
            test_range_imag = [butcherOrder*0.25, butcherOrder*0.25]


        else:
            raise Exception("Unsupported REXI method")

        function = Functions(
            function_name = function_name,
            efloat_mode = "float"
        )


        print("")
        print(unique_id_string)
        print(" + function_name: "+function_name)

        if test_range_real != None:
            max_error = 0
            for x in np.linspace(test_range_real[0], test_range_real[1], num_test_samples):
                lam = x

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
                print(" + test_range_real: ["+str(test_range_real[0])+", "+str(test_range_real[1])+"]")
                print(" + Error: "+str(max_error))

            if max_error > error_eps:
                raise Exception("Error threshold "+str(error_eps)+" exceeded")



        if test_range_imag != None:
            max_error = 0
            for x in np.linspace(test_range_imag[0], test_range_imag[1], num_test_samples):
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
                print(" + test_range_imag: ["+str(test_range_imag[0])+", "+str(test_range_imag[1])+"]")
                print(" + Error: "+str(max_error))

            if max_error > error_eps:
                raise Exception("Error threshold "+str(error_eps)+" exceeded")

        #coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_txt.txt", False)
        #coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_bin.txt", True)

