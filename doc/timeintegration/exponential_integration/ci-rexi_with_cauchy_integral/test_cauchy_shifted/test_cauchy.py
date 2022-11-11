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
def testPhi(
		P,
		R,
		Rl,
		z0,
		phiN,
		pde_id,
		output_console = True,
		output_plot = True
	):

	#
	# Phi function name
	#
	phistr = "phi"+str(phiN)

	cq = CauchyPhiQuadrature(phiN, P, R, z0, False, Rl)
	if output_plot:
		cq.plot("cauchy_"+phistr+"_quadrature_range.pdf")



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


	if output_console:
		print("PDE ID "+str(pde_id))
		print("Description: "+desc)


	#
	# Setup CauchyPhiQuadrature
	#
	cq = CauchyPhiQuadrature(phiN, P, R, z0, False, Rl)


	#if True:
	if False:
		if output_console:
			print("ALPHA:")
			for i in range(len(cq.alpha)):
				print(str(i)+"\t"+str(cq.alpha[i]))

			print("BETA:")
			for i in range(len(cq.beta)):
				print(str(i)+"\t"+str(cq.beta[i]))

	if True:
	#if False:
		start = -100
		end = -start
		dt = (end-start)/100
	else:
		start = -R*1.3
		end = -start
		dt = (end-start)/30


	eps = 1e-10
	U0 = np.zeros(S)
	U0[0] = 1.666


	t = start
	globalmaxerror = 0
	x = []
	y = []
	hits = 0
	while t < end:
		anal = cq.analytical_phi_pde(t*L, U0)
		approx = cq.approx_phi_pde(t*L, U0)
		localmaxerr = max(abs(approx-anal))
		globalmaxerror = max(globalmaxerror, localmaxerr)

		if output_console:
			output = ""
			output += str(round(t, 2))
			#output += "\tanal = "+str(anal)
			#output += "\tnumerical= "+str(approx)

			output += "\tlocal_maxerror = "+str(localmaxerr)

			print(output)

		x.append(t)
		y.append(localmaxerr)
		t += dt

		if localmaxerr < eps:
			hits = hits + 1

	if output_console:
		print("max error = "+str(globalmaxerror))
		print("hits = "+str(hits))

	if output_plot:
		plt.clf()
		plt.yscale('log')
		plt.title(desc)
		plt.plot(x, y, "-bo")
		plt.savefig("cauchy_"+phistr+"_errors_pde_"+str(pde_id)+".pdf")

	return hits


def testPhiShifted(
		P,
		max_real_evalue,
		max_imag_evalue,
		phiN,
		pde_id,
		output_console = True,
		output_plot = True
	):
		x0 = max_real_evalue
		xm = max_imag_evalue
		r = (x0*x0 + xm*xm)/(2.0*x0)
		center = max_real_evalue-r

		testPhi(P, r, [], center, phiN, pde_id, output_console, output_plot)

		if output_console:
			print("max_real_evalue: "+str(max_real_evalue))
			print("max_imag_evalue: "+str(max_imag_evalue))
			print("radius: "+str(r))
			print("center: "+str(center))



P = 128
R = 25

hits = testPhiShifted(P, 10, R, 0, 1, True, True)
hits = testPhi(P, R, [], 0, 0, 1, True, True)

print(hits)
sys.exit(1)


# Number of quadrature poles
print("P\tR\tz0\thits")
for P in [64*p for p in range(1, 10)]:
	for R in [10*r for r in range(1, 10+1)]:
		Rl = []
		# Shift
		for z0 in [5*z for z in range(-20, 20+1)]:
			#for phiN in range(5+1):
			for phiN in range(1):
				#for pde_id in range(3):
				for pde_id in range(1):
					hits = testPhi(P, R, Rl, z0, phiN, pde_id, False, False)
					print(str(P)+"\t"+str(R)+"\t"+str(z0)+"\t"+str(hits))
					sys.stdout.flush()



