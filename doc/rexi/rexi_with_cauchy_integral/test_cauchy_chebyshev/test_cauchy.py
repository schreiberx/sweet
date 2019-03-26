#! /usr/bin/env python3
#
# Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
# Date: 2017-08-16
#


import math
import cmath
import numpy as np
import sys
from CauchyPhiQuadrature import CauchyPhiQuadrature
import matplotlib.pyplot as plt



#
# Cauchy quadrature starts here
#

#
# Contour boundary
#
contour_boundary = {
	'shape': 'circle',
	'R': 10,		# radius
	'mu' : 1.0,		# shift
}


#
# Contour integration method
#
cont_int_method = 'trapezoidal'
#cont_int_method = 'chebyshev'

if cont_int_method == 'trapezoidal':
	contour_int_method = {
		'method': cont_int_method,	# contour integral method
		'N': 64,		# number of poles
	}
	testR = 10

elif cont_int_method == 'chebyshev':
	contour_int_method = {
		'method': cont_int_method,	# contour integral method
		'N': 20,		# number of poles
	}
	testR = 10


#for phiN in range(5+1):
for phiN in range(1):
	#
	# Phi function name
	#
	phistr = "phi"+str(phiN)
	print("*"*80)
	print(phistr)

	#
	# Use only upper half
	#
	#half = True
	half = False

	cq = CauchyPhiQuadrature(phiN, contour_boundary, contour_int_method)
	cq.plot("cauchy_"+phistr+"_quadrature_range.pdf")



	for pde_id in range(5):
		#
		# Setup linear operator
		#
		# S: size
		# L: linear operator
		#
		if pde_id == 0:
			S = 1
			L = np.array([[-1.0j]])
			desc = "ODE Oscillation, lambda="+str(L[0][0])

		if pde_id == 1:
			S = 1
			L = np.array([[1.0j]])
			desc = "ODE Oscillation, lambda="+str(L[0][0])

		elif pde_id == 2:
			S = 1
			L = np.array([[-1.0]])
			desc = "ODE Diffusion, lambda="+str(L[0][0])

		elif pde_id == 3:
			S = 1
			L = np.array([[1.0]])
			desc = "ODE Antidiffusion, lambda="+str(L[0][0])

		elif pde_id == 4:
			S = 2
			# PDE Oscillator
			L = np.array([[0.0, 1.0], [-1.0, 0]])
			desc = "PDE Oscillator, lambda=+-1j"



		print("PDE ID "+str(pde_id))

		#
		# Setup CauchyPhiQuadrature
		#
		cq = CauchyPhiQuadrature(phiN, contour_boundary, contour_int_method)


		#if True:
		if False:
			print("ALPHA:")
			for i in range(len(cq.alpha)):
				print(str(i)+"\t"+str(cq.alpha[i]))

			print("BETA:")
			for i in range(len(cq.beta)):
				print(str(i)+"\t"+str(cq.beta[i]))

		start = -testR*1.3
		end = -start
		dt = (end-start)/30

		U0 = np.zeros(S)
		U0[0] = 1.666



		if True:
		#if False:
			#
			# Treat everything as PDE (also diagonal L's)
			#
			t = start
			maxerror = 0
			x = []
			y = []
			while t < end:
				anal = cq.analytical_phi_pde(t*L, U0)
				approx = cq.approx_phi_pde(t*L, U0)
				err = sum(abs(approx-anal))
				maxerror = max(maxerror, err)

				output = ""
				output += str(round(t, 2))
				#output += "\tanal = "+str(anal)
				#output += "\tnumerical= "+str(approx)

				output += "\terror = "+str(err)

				print(output)

				x.append(t)
				y.append(err)
				t += dt

			print("max error = "+str(maxerror))

			plt.clf()
			plt.yscale('log')
			plt.title(desc)
			plt.plot(x, y, "-bo")
			plt.savefig("cauchy_"+phistr+"_errors_pde_"+str(pde_id)+".pdf")

		else:
			if S != 1:
				raise Exception("Only 1x1 Linear operators are supported")

			L = L[0][0]
			U0 = U0[0]

			t = start+dt
			maxerror = 0
			while t < end:
				anal = cq.analytical_phi_ode(t*L, U0)
				approx = cq.approx_phi_ode(t*L, U0)
				err = abs(approx-anal)
				maxerror = max(maxerror, err)

				output = ""
				output += str(round(t, 2))
				#output += "\tanal = "+str(anal)
				#output += "\tnumerical= "+str(approx)
				output += "\terror = "+str(err)

				print(output)

				t += dt

			print("max error = "+str(maxerror))

