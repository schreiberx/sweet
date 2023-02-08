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

met_coefs = [None]*n_met
met_names = [None]*n_met

for imet, rexi_method in enumerate(rexi_methods):

    function = Functions(
        function_name = "phi3",
        efloat_mode = "mpfloat"
    )
    
    for function_name in functions:
        #print("Function_name: "+function_name)
        
        # CI-REXI: Value on imaginary axis to be included
        #same used for ELREXI as max imag of ellipse
        lambda_include_imag = 10

        # CI-REXI: Maximum value of quadrature pole
        #same used for ELREXI as max real for ellipse
        lambda_max_real = 9
            
        # Error to test for
        error_eps = 1e-8       

        # efloat_mode
        efloat_mode = "float"
        #efloat_mode = "mpfloat"

        coeffs = None

        if rexi_method == "trexi":
            # T-REXI: Number of Gaussian basis functions
            M = 64

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
            N = 64

            cirexi = CIREXI(efloat_mode=efloat_mode)

            coeffs = cirexi.setup(
                function_name = function_name,
                N = N,
                lambda_max_real = lambda_max_real,
                lambda_include_imag = lambda_include_imag
            )

            unique_id_string = cirexi.getUniqueId()

        elif rexi_method == "elrexi":
            # CI-REXI: Number of quadrature poles
            #N = M
            M = 40

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
            k = 40

            if function_name != "phi0":
                continue

            fcfrexi = FCFREXI(efloat_mode=efloat_mode)

            coeffs = fcfrexi.setup(
                function_name = function_name,
                K = k,
                rat_num = int(k-1),
                rat_denom = int(k),
                lambda_max_imag = lambda_include_imag
            )

            unique_id_string = fcfrexi.getUniqueId()

        elif rexi_method == "brexi_gauss":
            
            # Butcher-REXI
            butcherOrder = 16

            if function_name != "phi0":
                continue

            brexi = BREXI(efloat_mode=efloat_mode)
            coeffs = brexi.setup(N=butcherOrder, quadrature_method="gauss")

            unique_id_string = brexi.getUniqueId()
            
        elif rexi_method == "brexi_radau":
            
            # Butcher-REXI
            butcherOrder = 16

            if function_name != "phi0":
                continue

            brexi = BREXI(efloat_mode=efloat_mode)
            coeffs = brexi.setup(N=butcherOrder, quadrature_method="radau")

            unique_id_string = brexi.getUniqueId()
        

        elif rexi_method == "brexi_chebyshev":
            # Butcher-REXI
            butcherOrder = 16
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
            butcherOrder = 16
        
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
        
        met_coefs[imet]=coeffs
        #print(coeffs.alphas)
        coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_txt.txt", False)
        coeffs.write_file("/tmp/REXI_"+rexi_method+"_"+unique_id_string+"_bin.txt", True)

    #Save last names used
    met_names[imet]=unique_id_string

fig, axes = plt.subplots(m_met,2, figsize=(12,10))

plt.suptitle(" Poles for different schemes")
levels = np.geomspace(1e-16, 1, num=17)
cmap=cm.viridis


axes1d=axes.reshape(-1)
for i, rexi_method in enumerate(rexi_methods):
    poles=np.asarray(met_coefs[i].alphas)
    cs = axes1d[i].scatter(poles.real,poles.imag, label=str(met_names[i]), marker=".")
    axes1d[i].set(xlabel="real", ylabel="imag", title=str(met_names[i]))
    axes1d[i].axis('equal')  
    
fig.tight_layout(pad=10.0)
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(top=0.9)

plt.savefig("poles.png")
plt.show()


fig, axes = plt.subplots(m_met,2, figsize=(12,10))

plt.suptitle(" Beta coefficients for different schemes")
levels = np.geomspace(1e-16, 1, num=17)
cmap=cm.viridis


axes1d=axes.reshape(-1)
for i, rexi_method in enumerate(rexi_methods):
    poles=np.asarray(met_coefs[i].betas)
    cs = axes1d[i].scatter(poles.real,poles.imag, label=str(met_names[i]), marker=".")
    axes1d[i].set(xlabel="real", ylabel="imag", title=str(met_names[i]))
    axes1d[i].axis('equal')  
    
fig.tight_layout(pad=10.0)
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(top=0.9)

plt.savefig("betas.png")
plt.show()