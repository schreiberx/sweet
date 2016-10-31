#! /usr/bin/env python3

#
# GaussianApproximation.py
#
#  Created on: 2 Aug 2015
#      Author: Martin Schreiber <schreiberx@gmail.com>
#
# Changelog:
# 	2016-05-14: Converted to Python
#                   Source: https://github.com/schreiberx/sweet
#

import math
import cmath
import sys

#
# This class provides the weights and coefficients for the
# approximation of a Gaussian
#
# \f$
# 	  exp(-(x*x)/(4*h*h))/sqrt(4 \pi)
# \f$
#
# with a sum over complex rational functions.
#
# See e.g. Near optimal rational approximations of large data sets, Damle et. al.
#

class GaussianApproximation:

	def __init__(self):
		"""	
		mu and a coefficients from
		"A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
		"""

		self.mu = -4.315321510875024 + 1j*0

		self.L = 11

		self.a = [
					-1.0845749544592896e-7 + 1j*				2.77075431662228e-8		,
					1.858753344202957e-8 + 1j*				-9.105375434750162e-7	,
					3.6743713227243024e-6 + 1j*				7.073284346322969e-7	,
					-2.7990058083347696e-6 + 1j*				0.0000112564827639346	,
					0.000014918577548849352 + 1j*			-0.0000316278486761932	,
					-0.0010751767283285608 + 1j*				-0.00047282220513073084	,
					0.003816465653840016 + 1j*				0.017839810396560574	,
					0.12124105653274578 + 1j*				-0.12327042473830248	,
					-0.9774980792734348 + 1j*				-0.1877130220537587		,
					1.3432866123333178 + 1j*					3.2034715228495942	,
					4.072408546157305 + 1j*					-6.123755543580666		,
					-9.442699917778205 + 1j*					0.			,
					4.072408620272648 + 1j*					6.123755841848161		,
					1.3432860877712938 + 1j*					-3.2034712658530275	,
					-0.9774985292598916 + 1j*				0.18771238018072134		,
					0.1212417070363373 + 1j*					0.12326987628935386	,
					0.0038169724770333343 + 1j*				-0.017839242222443888	,
					-0.0010756025812659208 + 1j*				0.0004731874917343858	,
					0.000014713754789095218 + 1j*			0.000031358475831136815	,
					-2.659323898804944e-6 + 1j*				-0.000011341571201752273	,
					3.6970377676364553e-6 + 1j*				-6.517457477594937e-7	,
					3.883933649142257e-9 + 1j*				9.128496023863376e-7	,
					-1.0816457995911385e-7 + 1j*				-2.954309729192276e-8	,
		]


	def evalGaussian(	
			self,
			x,		# x-coefficient for Gaussian basis function
			h		# h-coefficient for Gaussian basis function
	):
		return math.exp(-(x*x)/(4.0*h*h))/math.sqrt(4.0*math.pi)


	"""
	evaluate approximation of Gaussian basis function
	with sum of complex rational functions
	"""
	def approxGaussian(	
			self,
			x,		# x-coefficient for Gaussian basis function
			h		# h-coefficient for Gaussian basis function
	):
		# scale x, since it depends linearly on h:
		# x^2 ~ h^2
		x /= h

		sum = 0

		for l in range(0, len(self.a)):
			j = l-self.L

			# WORKS with max error 7.15344e-13
			sum += (self.a[l]/(1j*x + self.mu + 1j*j)).real

		return sum

	def print():
		print("mu: "+str(self.mu))

		for l in range(0, len(self.a)):
			print(self.a(l))


#
# This class computes an approximation of an exponential \f$ e^{ix} \f$ by
# a combination of Gaussians.
#
# The approximation is given by
#
# \f$
# 		e^{ix} - \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)}
# \f$
#
# see eq. (3.1) in
# "A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
#
# Here, the coefficents b_m are given by the coefficents \f$ c_m \f$ (Yes, it's the \f$ c_m \f$ here!)
# in the equation below equation (3.4):
#
# \f$
#    c_m := \int_{-1/(2h)}^{1/2h}  exp(-2 \pi i m h \xi) F(\xi)/\psi_h(\xi) d \xi
# \f$
#
# with F and \f$ \psi \f$ the functions f and \f$ \psi \f$ in Fourier space
#
class ExponentialApproximation:

	def __init__(
			self,
			i_h,
			i_M
	):
		self.ga = GaussianApproximation()

		self.h = i_h
		self.M = i_M

		 #
		 # See Section 3.1 for approx. of general function
		 #
		 # Note, that in the following we use the Fourier transformation
		 #  \f$ F(f(x)) := (xi) -> \int_{-\inf}^\inf {  exp(-i*2*\pi*x*\xi)  } dx \f$
		 # hence, with 2*Pi in the exponent
		 #
		 # STEP 1: specialize on \f$ f(x) = e^{i x} \f$
		 #
		 # In Fourier space, the function f(x) is given by
		 #
		 # \f$
		 # 		F(\xi) = \delta( \xi-1.0/(2 \pi) )
		 # \f$
		 #
		 # with \delta the Kronecker delta
		 #
		 # This simplifies the equation
		 #
		 # \f$
		 #    c_m := h * int_{-1/(2 h)}^{1/(2 h)}  exp(-2 \pi i m h \xi) F(\xi)/Psi_h(\xi)   d \xi
		 # \f$
		 #
		 # to
		 #
		 # \f$
		 #    c_m := h * exp(-i m h) / Psi_h( 1/(2*\pi) )
		 # \f$
		 #
		 #
		 # STEP 2: Compute \f$ Psi_h( 1/(2*\pi) ) \f$
		 #
		 # Furthermore, \psi_h(\xi) is given by
		 #
		 # \f$
		 #    \psi_h(xi) := 1/sqrt(2) * e^{-(2*\pi*h*\xi)^2}
		 # \f$
		 #
		 # and restricting it to xi=1/(2*pi) (see Kronecker delta above), yields
		 #
		 # \f$
		 #    \psi_h(1/2 \pi) := 1/\sqrt{2} * e^{-h^2}
		 # \f$
		 #
		 # URL:
		 # http://www.wolframalpha.com/input/?i=FourierTransform[1%2Fsqrt%284*pi%29*exp%28-x^2%2F%284*h^2%29%29%2C+x%2C+\[Omega]%2C+FourierParameters+-%3E+{0%2C+2+Pi}]
		 # RESULT: \f$ exp(-4 h^2 \pi^2 \omega^2)/\sqrt(1/h^2)
		 #         = h * exp(-(2 h \pi \omega)^2) \f$
		 # and with \f$ \xi = 1/(2 \pi) \f$:
		 # 			\f$ h * exp(-h^2) \f$
		 #
		 # STEP 3: combine (1) and (2):
		 #
		 # This simplifies c_m to
		 #
		 # \f$
		 #    c_m := e^{h^2} * e^{-i*m*h}
		 # \f$
		 #
		 # Let's hope, that these equations are right.
		 #

		self.b = [math.exp(self.h*self.h)*cmath.exp(-1j*(float(m)*self.h)) for m in range(-self.M, self.M+1)]

	def print(self):
		for i in self.b:
			print(i)


	def eval_e_ix(
		self,
		i_x
	):
		return cmath.exp(1j*i_x)


	def approx_e_ix(
		self,
		i_x
	):
		sum = 0

		# \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$
		for m in range(-self.M, self.M+1):
			sum += self.b[m+self.M] * self.ga.approxGaussian(i_x+float(m)*h, h)

		return sum


class REXI:
	
	def __init__(
		self,
		i_h,
		i_M,
		i_reduce_to_half = True
	):
		self.ga = GaussianApproximation();
		self.ea = ExponentialApproximation(i_h, i_M);

		self.L = self.ga.L;
		self.N = i_M+self.ga.L;
		self.M = i_M;

		M = self.M
		N = self.N
		L = self.L

		self.alpha = [0 for i in range(0, 2*self.N+1)]
		self.beta_re = [0 for i in range(0, 2*self.N+1)]
		self.beta_im = [0 for i in range(0, 2*self.N+1)]

		for l in range(-self.L, L+1):
			for m in range(-M, M+1):
				n = l+m;
				self.alpha[n+N] = i_h*(self.ga.mu + 1j*n);

				self.beta_re[n+N] += self.ea.b[m+M].real*i_h*self.ga.a[l+L];
				self.beta_im[n+N] += self.ea.b[m+M].imag*i_h*self.ga.a[l+L];

		#for (int n = -N; n < N+1; n++)
		#{
		#	alpha[n+N] = i_h*(ga.mu + complex(0, n));
		#
		#	int L1 = std::max(-L, n-M);
		#	int L2 = std::min(L, n+M);
		#
		#	beta_re[n+N] = 0;
		#	for (int k = L1; k < L2; k++)
		#	{
		#		assert(k+L >= 0);
		#		assert(k+L < 2*L+1);
		#		assert(n-k+M >= 0);
		#		assert(n-k+M < 2*M+1);
		#
		#		beta_re[n+N] += ga.a[k+L]*ea.b[n-k+M].real();
		#	}
		#
		#	beta_re[n+N] *= i_h;
		#}


		if i_reduce_to_half:
			# reduce the computational amount to its half,
			# see understanding REXI in the documentation folder

			self.alpha = self.alpha[:N+1]
			self.beta_re = self.beta_re[:N+1]
			self.beta_im = self.beta_im[:N+1]

			# N+1 contains the pole and we don't rescale this one by 2 but all the other ones
			for i in range(0, N):
				self.beta_re[i] *= 2.0
				self.beta_im[i] *= 2.0

	#	
	# \return \f$ cos(x) + i*sin(x) \f$
	#
	def eval_e_ix(self, i_x):
		return cmath.exp(1j*i_x)

	#
	# approx
	#
	def approx_e_ix(self, i_x):
		sum_re = 0;
		sum_im = 0;

		S = len(self.alpha)

		# Split computation into real part of \f$ cos(x) \f$ and imaginary part \f$ sin(x) \f$
		for n in range(0, S):
			denom = (1j*i_x + self.alpha[n]);
			sum_re += (self.beta_re[n] / denom).real
			sum_im += (self.beta_im[n] / denom).real

		return sum_re + 1j*sum_im;


##################################################
##################################################

h = 0.2
M = 32

##################################################
##################################################


print("GaussianApproximation")

g = GaussianApproximation()
h = 1
for x in range(0, 10):
	a = g.evalGaussian(x, h)
	b = g.approxGaussian(x, h)

	print(str(a)+"\t"+str(b)+"\t"+str(abs(a-b)))

##################################################
##################################################

print()

print("ExponentialApproximation")
ea = ExponentialApproximation(h, M)

for x in range(-int(h*M)+1, int(h*M)):
	a = ea.eval_e_ix(x)
	b = ea.approx_e_ix(x)

	print(str(a)+"\t"+str(b)+"\t"+str(abs(a-b)))

##################################################
##################################################

print()

print("REXI")
h = 0.2
M = 64
rexi = REXI(h, M)

for x in range(-int(h*M)+1, int(h*M)):
	a = rexi.eval_e_ix(x)
	b = rexi.approx_e_ix(x)

	print(str(a)+"\t"+str(b)+"\t"+str(abs(a-b)))


