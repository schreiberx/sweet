#! /usr/bin/env python3

import sys
import os
import numpy as np

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+"/..")

from Functions import *
from cirexi.CIREXI import *


sys.path.pop()


for rexi_method in ["cirexi"]:

    for function_name in ["phi0", "phi1", "phi2"]: #, "ups1", "ups2", "ups3"]:
        
        # CI-REXI: Number of quadrature poles
        for N in [32, 128, 256]:
            print()
            print("function: "+str(function_name))
            print("N: "+str(N))
            
            # CI-REXI: Value on imaginary axis to be included
            #same used for ELREXI as max imag of ellipse
            lambda_include_imag = 10
            
            # CI-REXI: Maximum value of quadrature pole
            #same used for ELREXI as max real for ellipse
            lambda_max_real = 20
            
            # Error to test for
            error_eps = 1e-10
            
            # efloat_mode
            efloat_mode = "float"
            #efloat_mode = "mpfloat"
            
            cirexi = CIREXI(efloat_mode=efloat_mode)
            
            coeffs = cirexi.setup(
                function_name = function_name,
                N = N,
                lambda_max_real = lambda_max_real,
                lambda_include_imag = lambda_include_imag,
                half_shifted = False
            )
            
            coeffs.normalize_steady_state()
            
            # Convert to floating point
            coeffs = coeffs.toFloat()
            
            unique_id_string = cirexi.getUniqueId()
            
            function = Functions(
                function_name = function_name,
                efloat_mode = "float"
            )
            
            verbosity = 1
            
            if verbosity > 1:
                print("")
                print("COEFF_NR\tALPHAS\t\t\tBETAS")
                for i in range(len(coeffs.alphas)):
                    print(str(i)+"\t"+str(coeffs.alphas[i])+"\t"+str(coeffs.betas[i]))
            
            max_error = 0
            
            lam = 0
            
            y = function.eval(lam)
            
            if function_name in ['phi0', 'phi1']:
                """
                Should be close to 1.0
                """
                assert np.isclose(y, 1.0)
            
            yn = coeffs.eval(lam)
            
            err = np.abs(y-yn)
            
            if verbosity > 0:
                print("error="+str(err))
            
            if err > error_eps:
                raise Exception("Error threshold "+str(error_eps)+" exceeded")
            
