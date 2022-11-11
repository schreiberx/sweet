#! /usr/bin/env python3
#
# Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
# Date: 2017-08-16
#


import math
import cmath
import numpy as np
import sys

import matplotlib.pyplot as plt


class CauchyPhiQuadrature:
	alpha = []
	beta = []

	#
	# Phi 0-N functions
	#
	def phi(self, n, z):
		if n == 0:
			return cmath.exp(z)

		if n != 0:
			if abs(z) < 1e-8:
				return 1.0/math.factorial(n)

				raise Exception("Z close to zero, not yet supported for phi "+str(n)+" !!!")


		return (self.phi(n-1, z) - 1.0/math.factorial(n-1))/z

		raise Exception("Phi function not supported yet")


	#
	# Constructor
	# See setup(...) for documentation on parameters
	# 
	def __init__(self, phiN, contour_boundary, contour_int_method):
		if phiN == -1:
			return

		self.setup(phiN, contour_boundary, contour_int_method)



	def setup(
		self,
		phiN,			# phi function id, use 0 for exp(lambda) standard ODE integration
		contour_boundary,	# boundary description
		contour_int_method	# integration method
	):
		self.contour_boundary = contour_boundary
		self.contour_int_method = contour_int_method

		if self.contour_boundary['shape'] == 'circle':
			self.R = contour_boundary['R']
			self.mu = contour_boundary['mu']

			if self.contour_int_method['method'] == 'trapezoidal':
				self.phiN = phiN
				self.N = contour_int_method['N']


				#
				# Compute support points of quadrature
				#
				self.coords = []
				for j in range(self.N):
					theta_j = 2.0*math.pi*(j+0.5)/self.N
					gamma_j = self.R*cmath.exp(1j*theta_j) + self.mu
					self.coords.append(gamma_j)


				self.alpha = []
				self.beta = []
				for j in range(self.N):
					theta_j = 2.0*math.pi*(j+0.5)/self.N
					gamma_j = self.R*cmath.exp(1j*theta_j) + self.mu

					k = self.R*cmath.exp(1j*theta_j)

					beta = -self.phi(phiN, gamma_j)*k
					beta /= self.N
					alpha = -(k + self.mu)

					self.alpha.append(alpha)
					self.beta.append(beta)

				return


			elif self.contour_int_method['method'] == 'chebyshev':
				self.phiN = phiN
				self.N = contour_int_method['N']


		raise Exception("Unsupported combination of contour_boundary "+contour_boundary['shape']+" and integration method "+contour_int_method['method'])


	def plot(self, filename = None):
		points_re = []
		points_im = []
		for j in range(self.N):
			points_re.append(self.coords[j].real)
			points_im.append(self.coords[j].imag)

		points_re.append(points_re[0])
		points_im.append(points_im[0])

		plt.plot(points_re, points_im,  '-bo')

		if filename != None:
			plt.savefig(filename)
		else:
			plt.show()



	def approx_phi_pde(self, dt_L, U):
		S = len(dt_L)

		accum = np.array([0.j, 0.j])
		for j in range(len(self.alpha)):
			M_inv = np.linalg.inv(dt_L + np.identity(S)*self.alpha[j])
			accum += self.beta[j] * np.dot(M_inv, U)

		return accum

	def approx_phi_ode(self, dt_L, U):
		accum = 0.0
		for j in range(len(self.alpha)):
			M_inv = 1.0/(dt_L + self.alpha[j])
			accum += self.beta[j] * M_inv * U

		return accum


	def analytical_phi_pde(self, dt_L, U):
		S = len(dt_L)

		# Setup eigenvalues and Eigenvectors for analytical solution
		LEvals, LEvecs = np.linalg.eig(dt_L)
		LEvecs_inv = np.linalg.inv(LEvecs)

		if True:
			error = np.sum(np.absolute(dt_L - np.dot(np.dot(LEvecs, np.diag(LEvals)), LEvecs_inv)))
			if error > 1e-10:
				raise Exception("Error "+str(error)+" too large")

		Uwave = np.dot(LEvecs_inv, U)
		tmp = np.array([self.phi(self.phiN, LEvals[i])*Uwave[i] for i in range(S)])
		return np.dot(LEvecs, tmp)

	def analytical_phi_ode(self, dt_L, U):
		# Setup eigenvalues and Eigenvectors for analytical solution
		return self.phi(self.phiN, dt_L)*U


