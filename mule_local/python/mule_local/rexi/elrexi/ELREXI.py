#! /usr/bin/env python3
#
# Author: Pedro Peixoto <ppeixoto@usp.br>
# Date 2020-03-13
# 
# based on 
# CIREXI.py by
# Martin Schreiber <M.Schreiber@exeter.ac.uk>
# Date: 2017-08-16
#

import math
import cmath
import numpy as np
import sys

from mule_local.rexi.EFloat import *
from mule_local.rexi.Functions import *
from mule_local.rexi.REXICoefficients import *

class ELREXI:

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
    	lambda_max_real: float = None,
    	lambda_max_imag: float = None
    ):
    	"""
    	Setup coefficients for Cauchy contour integral coefficients
    	using an ellipse with height given by lambda_max_imag
    	and width given by lambda_max_real
    	See doc/rexi/rexi_wit_laplace for details

    	Args:
    		N (int):
    			Number of points on Cauchy contour.

    	"""

    	self.alphas = []
    	self.betas = []
    	self.function_name = function_name
    	self.N = N
    	gamma = lambda_max_imag
    	delta = lambda_max_real
    	r = np.sqrt((1-delta/gamma) / (1+delta/gamma) )
    	d = r + (1.0/r)
    	self.fun = Functions(function_name, efloat_mode = self.efloat_mode)
    	
    	#print("Setting up ellipse:")
    	#print("gamma: ", gamma)
    	#print("delta: ", delta)
    	#print("r:     ", r)

    	if lambda_max_real == None or lambda_max_imag == None:
    		raise Exception("ELREXI: please provide max imag and real values of contour")

    	c1 = (gamma*self.efloat.i/d)
    	c2 = -(gamma/d)
    	for j in range(N):
    		theta_j = self.efloat.pi2*(j+self.efloat.to(0.5))/N

    		# sampling position of support point
    		sn = c1*(r*self.efloat.exp(self.efloat.i*theta_j)+(1.0/r)*self.efloat.exp(-self.efloat.i*theta_j))

    		#metric term
    		sn_prime = c2*(r*self.efloat.exp(self.efloat.i*theta_j)-(1.0/r)*self.efloat.exp(-self.efloat.i*theta_j))


    		self.betas.append(-self.efloat.i*self.fun.eval(sn)*sn_prime / N)
    		self.alphas.append(-sn)

    	coeffs = REXICoefficients()
    	coeffs.alphas = self.alphas[:]
    	coeffs.betas = self.betas[:]
    	coeffs.gamma = 0
    	coeffs.function_name = self.function_name

    	self.unique_id_string = ""
    	self.unique_id_string += self.function_name
    	self.unique_id_string += "_N"+str(self.N)
    	self.unique_id_string += "_im"+str(self.efloat.re(lambda_max_imag))
    	self.unique_id_string += "_re"+str(self.efloat.re(lambda_max_real))

    	coeffs.unique_id_string = self.getUniqueId()

    	return coeffs

    def getUniqueId(self):
    	return "ELREXI_"+self.unique_id_string
