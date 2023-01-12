#! /usr/bin/env python3

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

from Functions import *
from trexi.TREXI import *
from cirexi.CIREXI import *
from brexi.BREXI import *
from elrexi.ELREXI import *
#from fcfrexi.FCFREXI import *

# verbose
#verbosity = 100
verbosity = 10

rexi_methods = [
        "trexi",
        "cirexi",
        #"fcfrexi",
        "elrexi",
        #"brexi_tanhsinh",
        #"brexi_gauss",
        #"brexi_radau",
        "brexi_chebyshev"
    ]
n_met = len(rexi_methods)
m_met = int(n_met/2)

#Most rexi functions only work for phi0
#functions = ["phi0", "phi1", "phi2", "phi3"]: #, "ups1", "ups2", "ups3"]:
functions = ["phi0"] #, "phi1", "phi2", "phi3"]: #, "ups1", "ups2", "ups3"]:

errors = [None]*n_met
met_names = [None]*n_met

for imet, rexi_method in enumerate(rexi_methods):

    function = Functions(
        function_name = "phi3",
        efloat_mode = "mpfloat"
    )

    for function_name in functions:
        print("Function_name: "+function_name)

        # CI-REXI: Value on imaginary axis to be included
        #same used for ELREXI as max imag of ellipse
        lambda_include_imag = 20

        # CI-REXI: Maximum value of quadrature pole
        #same used for ELREXI as max real for ellipse
        lambda_max_real = 8
            
        #Parameter to be looped over
        nend = 128
        nini = 8
        K_tests = np.arange(nini, nend, 4)

        # Testing: Range (start, end)
        test_range = [-25,25]
        print("Test_range: ["+str(test_range[0])+", "+str(test_range[1])+"]")
        
        # Testing: number of samples
        num_test_samples = 1000 #12345
        x_points = np.linspace(test_range[0], test_range[1], num_test_samples)

        # Error to test for
        error_eps = 1e-8       

        # efloat_mode
        efloat_mode = "float"
        #efloat_mode = "mpfloat"

        coeffs = None

        #Error matrices
        errors[imet]=np.ones((len(K_tests), num_test_samples))

        for j, k in enumerate(K_tests):

            if rexi_method == "trexi":
                # T-REXI: Number of Gaussian basis functions
                M = k #number of poles 256

                # T-REXI: Spacing between Gaussian basis functions
                h = 0.2

                if function_name != "phi0":
                    continue

                trexi = TREXI(efloat_mode=efloat_mode)
                coeffs = trexi.setup(
                    function_name = function_name,
                    M = M,
                    h = h
                )

                unique_id_string = trexi.getUniqueId()

            elif rexi_method == "cirexi":
                # CI-REXI: Number of quadrature poles
                #N = M
                N = k #256

                cirexi = CIREXI(efloat_mode=efloat_mode)

                coeffs = cirexi.setup(
                    function_name = function_name,
                    N = k,
                    lambda_max_real = lambda_max_real,
                    lambda_include_imag = lambda_include_imag
                )

                unique_id_string = cirexi.getUniqueId()

            elif rexi_method == "elrexi":
                # CI-REXI: Number of quadrature poles
                #N = M
                M = k #256

                elrexi = ELREXI(efloat_mode=efloat_mode)

                coeffs = elrexi.setup(
                    function_name = function_name,
                    N = M,
                    lambda_max_real = lambda_max_real,
                    lambda_max_imag = lambda_include_imag
                )

                unique_id_string = elrexi.getUniqueId()

            elif rexi_method == "fcfrexi":
                # CI-REXI: Number of quadrature poles
                #N = M
                M = k #256

                if k > 70:
                    continue

                if function_name != "phi0":
                    continue
                
                
                fcfrexi = FCFREXI(efloat_mode=efloat_mode)
                
                try:
                    coeffs = fcfrexi.setup(
                        function_name = function_name,
                        K = k,
                        rat_num = int(k-1),
                        rat_denom = int(k),
                        lambda_max_imag = lambda_include_imag
                    )
                except:
                    continue

                unique_id_string = fcfrexi.getUniqueId()

            elif rexi_method == "brexi_gauss":
                
                # Butcher-REXI
                butcherOrder = k

                if function_name != "phi0":
                    continue

                brexi = BREXI(efloat_mode=efloat_mode)
                coeffs = brexi.setup(N=butcherOrder, quadrature_method="gauss")

                unique_id_string = brexi.getUniqueId()
                
            elif rexi_method == "brexi_radau":
                
                # Butcher-REXI
                butcherOrder = k

                if function_name != "phi0":
                    continue

                brexi = BREXI(efloat_mode=efloat_mode)
                coeffs = brexi.setup(N=butcherOrder, quadrature_method="radau")

                unique_id_string = brexi.getUniqueId()
            

            elif rexi_method == "brexi_chebyshev":
                # Butcher-REXI
                butcherOrder = k
                if butcherOrder > 60:
                    continue

                if function_name != "phi0":
                    continue

                brexi = BREXI(efloat_mode=efloat_mode)
                coeffs = brexi.setup(N=butcherOrder, quadrature_method="chebyshev")

                # Convert to floating point
                coeffs = coeffs.toFloat()
                unique_id_string = brexi.getUniqueId()

            elif rexi_method == "brexi_tanhsinh":

                # Butcher-REXI
                butcherOrder = k
            
                if function_name != "phi0":
                    continue

                brexi = BREXI(efloat_mode=efloat_mode)
                coeffs = brexi.setup(N=butcherOrder, quadrature_method="tanhsinh")

                unique_id_string = brexi.getUniqueId()

            else:
                raise Exception("Unsupported REXI method")

            # Convert to floating point
            coeffs = coeffs.toFloat()

            function = Functions(
                function_name = function_name,
                efloat_mode = "float"
            )
           
            max_error = 0
            
            for i, x in enumerate(x_points):
                lam = 1j*x
                
                y = function.eval(lam)
                yn = coeffs.eval(lam)

                err = np.abs(y-yn)
                max_error = max(max_error, err)

                #For plotting purposes
                if err > 10.0:
                    err = np.nan
                if err < 1e-15:
                    err = 1e-16
                errors[imet][j,i]=err

                if verbosity > 10:
                        print(unique_id_string+" + Lambda: "+str(lam)+" +  exact: "+str(y)+" + approx: "+str(yn)+
                        " + Error: "+str(errors[imet][j,i]), imet, i, j)
        
            if verbosity > 0:
                print(unique_id_string+" + Error: "+str(max_error))

            #if max_error > error_eps:
            #    raise Exception("Error threshold "+str(error_eps)+" exceeded")
            
            coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_txt.txt", False)
            coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_bin.txt", True)

    #Save last names used
    met_names[imet]=unique_id_string


fig, axes = plt.subplots(m_met,2, figsize=(12,12))

plt.suptitle(r'$e^{ix}$'+" reconstruction error")
levels = np.geomspace(1e-16, 1, num=17)
cmap=cm.viridis

X, Y = np.meshgrid(x_points, K_tests )
axes1d=axes.reshape(-1)
for i, rexi_method in enumerate(rexi_methods):
    cs = axes1d[i].contourf(X, Y, errors[i], levels, locator=ticker.LogLocator(), cmap=cmap)
    axes1d[i].set(xlabel="x", ylabel="Poles", title=str(met_names[i]))    

fig.tight_layout(pad=4.0)
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(top=0.9)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar=fig.colorbar(cs, cax=cbar_ax)
cbar.set_label('Error in exp(ix)')
plt.figtext(0.3,0.01, "lambda_include_imag="+str(lambda_include_imag)+" lambda_max_real="+str(lambda_max_real) )
plt.savefig("expi_recon.png")

plt.show()
