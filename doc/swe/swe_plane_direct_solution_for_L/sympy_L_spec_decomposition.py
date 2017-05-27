#! /usr/bin/python2

#
# Compute analytical solution for linear operator in spectral space
#
# See Embid/Madja, Averaging over fast Gravity waves for geophysical flows with arbitrary potential vorticity, 1996
#
# This file is part of the SWEET development
# Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
#
#

from sympy.matrices import Matrix
from sympy import *
from sympy.printing import print_ccode
from sympy.abc import u, v, h
import re as re_exp
import sys


####################################
# SETUP START
####################################
#
# Choose for which L operator to compute the EV decomposition in spectral space:
# 	1) Embid and Madja formulation with U=(u,v,h)
# 	2) Similar to 1, but with U=(h,u,v)
#	3) Full-dimensional with U=(h,u,v)
#	4) Full-dimensional with U=(h,u,v), but with f=0
#
LNr=3

####################################
# SETUP END
####################################

# use real assumptions for complex valued optimizations (conjugate)

# discretization in 2D Cartesian space
x = symbols('x0:2', real=True, positive=True)

# discretization in 2D spectral space
k = symbols('k0:2', real=True, positive=True)

# Domain size
# TODO: Not yet working
s = symbols('s0:2', real=True, positive=True)
#s = [1,1]
#s = [1/s[0],1/s[1]]

# do the normalisation?
# If we use the inverse of the Eigenvector matrix (instead of the Hermitian of the Eigenvalue matrix),
# we don't need to do any normalisation
normalise = False
#normalise = True

# Setup some variables used in the respective formulation
if LNr in [1,2]:
	F = symbols('F', real=True, positive=True)
	Fsinv = 1/sqrt(F)
	#Fsinv = symbols('F^-1/2', real=True, positive=True)

elif LNr in [3,4]:
	H0 = symbols('H0', real=True, positive=True)
	g = symbols('g', real=True, positive=True)
	f = symbols('f', real=True, positive=True)


if LNr == 1:
	#
	# U=(u, v, h)
	# f = 1
	#
	U = Matrix([u, v, h])
	U_spec = Matrix([symbols('u_spec', real=True), symbols('v_spec', real=True), symbols('h_spec', real=True)])

	def getL(U):
		return Matrix(
			[
				[0*U[0],			-1*U[1],		Fsinv*diff(U[2], x[0])	],
				[1*U[0],			0*U[1],			Fsinv*diff(U[2], x[1])	],
				[Fsinv*diff(U[0], x[0]),	Fsinv*diff(U[1], x[1]),	0*U[2]			]
			]
		)

elif LNr == 2:
	#
	# U=(h, u, v)
	# f = 1
	#
	U = Matrix([h, u, v])
	U_spec = Matrix([symbols('h_spec', real=True), symbols('u_spec', real=True), symbols('v_spec', real=True)])

	def getL(U):
		return Matrix(
			[
				[0*U[0],			Fsinv*diff(U[1], x[0]),		Fsinv*diff(U[2], x[1])	],
				[Fsinv*diff(U[0], x[0]),	0*U[1],				-1*U[2]			],
				[Fsinv*diff(U[0], x[1]),	1*U[1],				0*U[2]			]
			]
		)

elif LNr == 3:
	#
	# f = variable
	#
	U = Matrix([h, u, v])
	U_spec = Matrix([symbols('h_spec', real=True), symbols('u_spec', real=True), symbols('v_spec', real=True)])

	def getL(U):
		return Matrix(
			[
				[0*U[0],		H0*diff(U[1], x[0]),		H0*diff(U[2], x[1])	],
				[g*diff(U[0], x[0]),	0*U[1],				-f*U[2]			],
				[g*diff(U[0], x[1]),	f*U[1],				0*U[2]			]
			]
		)

elif LNr == 4:
	#
	# f = 0
	#
	U = Matrix([h, u, v])
	U_spec = Matrix([symbols('h_spec', real=True), symbols('u_spec', real=True), symbols('v_spec', real=True)])

	def getL(U):
		return Matrix(
			[
				[0*U[0],		H0*diff(U[1], x[0]),		H0*diff(U[2], x[1])	],
				[g*diff(U[0], x[0]),	0*U[1],				0*U[2]			],
				[g*diff(U[0], x[1]),	0*U[1],				0*U[2]			]
			]
		)


# Substitutions to make equations nicer
eq_subs = {}
inv_eq_subs = {}


# Substitutions
if LNr in [1,2]:
	w = symbols('w', real=True, positive=True)
	omega_2 = 1 + Fsinv**2*pi*pi*k[0]*k[0] + Fsinv**2*pi*pi*k[1]*k[1]
	eq_subs[omega_2] = w*w
	eq_subs[-omega_2] = -w*w
	eq_subs[expand(omega_2*F)] = w*w*F
	eq_subs[expand(-omega_2*F)] = -w*w*F
	inv_eq_subs[w] = sqrt(omega_2)

elif LNr in [3,4]:
	w = symbols('w', real=True, positive=True)
	omega_2 = H0*g*2*2*pi*pi*k[0]*k[0]*s[1]*s[1]+H0*g*2*2*pi*pi*k[1]*k[1]*s[0]*s[0]+f*f*s[0]*s[0]*s[1]*s[1]
	eq_subs[omega_2] = w*w
	eq_subs[-omega_2] = -w*w
	inv_eq_subs[w] = sqrt(omega_2)

	wg = symbols('wg', real=True, positive=True)
	omegag_2 = g*g*2*2*pi*pi*k[0]*k[0]*s[1]*s[1]+g*g*2*2*pi*pi*k[1]*k[1]*s[0]*s[0]+f*f*s[0]*s[0]*s[1]*s[1]
	eq_subs[omegag_2] = wg*wg
	eq_subs[-omegag_2] = -wg*wg
	inv_eq_subs[wg] = sqrt(omegag_2)


# replace k0*k0 + k1*k1 by 'K2' symbol
K2 = symbols('K2', real=True, positive=True)
asdf = pi*pi*k[0]*k[0]+pi*pi*k[1]*k[1]
eq_subs[asdf] = K2
inv_eq_subs[K2] = asdf


# Solution expressed with variables in spectral space
U_spec = U_spec.transpose()

#
# Solution given in Cartesian space, but expressed in spectral space
# Remove 2*pi for "2pi-sized" domain
#
U_fromSpec = exp(I*2*pi*(k[0]*x[0]/s[0]+k[1]*x[1]/s[1])) * U_spec


print
print
print("L(U):")
M = getL(U_fromSpec)
pprint(M)

# Remove spectral basis components (This should be possible in general)
Lik = Matrix.zeros(3,3)
Lik[:,0] = simplify(M[:,0]/U_fromSpec[0])
Lik[:,1] = simplify(M[:,1]/U_fromSpec[1])
Lik[:,2] = simplify(M[:,2]/U_fromSpec[2])


if True:
	print
	print
	print("Lik(U):")
	pprint(Lik)


print
print
print("Using omega(k)^2:")
pprint(omega_2)
if LNr in [3,4]:
	print("Using omegag(k)^2:")
	pprint(omegag_2)
print


#
# Compute the entire eigenvector decomposition
#
def getEVStuff(Lik):

	# compute Eigenvalue/Eigenvector decomposition
	e_vals_vects = Lik.eigenvects()
	if True:
		print
		print
		print("Eigenvectors:")
		pprint(e_vals_vects)

	# check for multiplicity of 2
#	if 2 in e_vals_vects[:][1]:
#		print("Multiplicity of 2 not yet supported!")
#		sys.exit(1)

	# special handler for 3 times the same EV (special case for f=0)
	if e_vals_vects[0][1] == 3:
		# Multiplicity of 3 given, probably this is the solution without the Coriolis frequency
		evals = [simplify(e_vals_vects[0][0].subs(eq_subs)) for i in range(0,3)]
		multi = [simplify(e_vals_vects[0][1]) for i in range(0,3)]
		evects = [Matrix([simplify(e_vals_vects[0][2][i][j].subs(eq_subs)).subs(eq_subs) for j in range(0,3)]) for i in range(0,3)]

		if normalise:
			norm_facs = [simplify(sqrt(evects[i].dot(evects[0].conjugate()).subs(eq_subs)).subs(eq_subs)).subs(eq_subs) for i in range(0,3)]
		else:
			norm_facs = [1,1,1]

		nevects = [simplify(evects[i]/norm_facs[0]) for i in range(0,3)]

	else:
		evals = [simplify(e_vals_vects[i][0].subs(eq_subs)) for i in range(0,3)]
		multi = [simplify(e_vals_vects[i][1]) for i in range(0,3)]
		evects = [Matrix([simplify(e_vals_vects[i][2][0][j].subs(eq_subs)).subs(eq_subs) for j in range(0,3)]) for i in range(0,3)]

		if normalise:
			norm_facs = [simplify(sqrt(evects[i].dot(evects[i].conjugate()).subs(eq_subs)).subs(eq_subs)).subs(eq_subs) for i in range(0,3)]
		else:
			norm_facs = [1,1,1]

		nevects = [simplify(evects[i]/norm_facs[i]) for i in range(0,3)]

	for i in range(0,3):
		print
		print
		print("*"*80)
		print("Eigenvector/-values "+str(i)+":")
		print("*"*80)
		print
		print("Multiplicity:")
		print(pretty(multi[i]))
		print
		print("Eigenvalues:")
		print(pretty(evals[i]))
		if False:
			print
			print("Eigenvector:")
			print(pretty(evects[i]))
			if normalise:
				print
				print("Normalization denominator value (1/value):")
				print(pretty(norm_facs[i]))
		if normalise:
			print
			print("Orthonormal Eigenvector:")
			print(pretty(nevects[i]))


	print("*"*80)
	print("VALIDATION")
	print("*"*80)


	print("Eigenvector property Lik*V-lambda*V:")
	for i in range(0,3):
		vec = simplify(expand((Lik*evects[i]-evects[i]*evals[i]).subs(inv_eq_subs)))
		pprint(vec.transpose())

		if vec.dot(vec) != 0:
			print("EV FAILED or too complex!")
			sys.exit(1)

		print ("OK")


	if normalise:
		print("Eigenvector orthonormality, only for LNr=1 and LNr=2:")
		if LNr in [1,2]:
			# assure orthonormality of EVs
			for i in range(0,3):
				for j in range(0,3):
					va = nevects[i]
					vb = nevects[j]

					dot = va.dot(vb.conjugate())
					# note, we have to expand the complex stuff here, since the conjugate of a real-valued square root is not determined
					dot = dot.subs(inv_eq_subs).expand(complex=True)
					dot = simplify(dot)
					print(str(i)+", "+str(j)+": "+str(dot))

					dot_result = 1 if i == j else 0

					# don't check for orthonormality if the Eigenvectors cannot be orthonormal
					if LNr in [3,4] and dot_result != 0:
						continue

					if simplify(dot-dot_result) != 0:
						print "VALIDATION FAILED (orthonormality of EVs)"
						sys.exit(1)


	print("Eigenvalues - imaginary only:")
	# assure imaginary-only valued eigenvalues
	for i in range(0,3):
		pprint(re(evals[i]))
		if re(evals[i]) != 0:
			print "VALIDATION FAILED (imaginary-only eigenvalues)"
			sys.exit(1)


	print("*"*80)
	print("VALIDATION SUCCESSFUL!")
	print("*"*80)

	return (evals, nevects)


k0zero_subs = {k[0]:0, K2:pi*pi*k[1]*k[1]}
k1zero_subs = {k[1]:0, K2:pi*pi*k[0]*k[0]}
kzero_subs = {k[0]:0, k[1]:0, K2:0}


(evals, nevects) = getEVStuff(Lik)
(evals0, nevects0) = getEVStuff(Lik.subs(kzero_subs))
(evalsk0zero, nevectsk0zero) = getEVStuff(Lik.subs(k0zero_subs))
(evalsk1zero, nevectsk1zero) = getEVStuff(Lik.subs(k1zero_subs))


#
# function to make code nice looking
#
def code_nice_replace(code):
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 2.0\)' , r'\1*\1', code)
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 3.0\)' , r'\1*\1*\1', code)
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 4.0\)' , r'\1*\1*\1*\1', code)
	code = re_exp.sub(r'sqrt\(' , r'std::sqrt((complex)', code)
	code = re_exp.sub(r'1\.0L' , r'1.0', code)
	code = re_exp.sub(r'2\.0L' , r'2.0', code)
	code = re_exp.sub(r'^([0-9][0-9]*)([*])' , r'\1.\2', code)
	return code


def _ccode(expr):
	for i in range(0, 30):
		expr = expr.subs(i, float(i))
	return code_nice_replace(ccode(expr))


print "*"*80
print "*"*80
print "*"*80

#
# Output C-Code
#
if True:
	print
	print("*"*80)
	print("* INFO: s0/s1: Domain size in X/Y direction")
	print("*"*80)
	print
	print("*"*80)
	print(" C CODE OUTPUT for LNr "+str(LNr))
	print("*"*80)
	print
	print("std::complex<double> eigenvalues[3];")
	print("std::complex<double> eigenvectors[3][3];")
	print
	print("if (k0 == 0 && k1 == 0)")
	print("{")

	if LNr == 3:
		print("\tcomplex wg = std::sqrt((complex)"+_ccode(simplify(omegag_2.subs(kzero_subs)))+");")

	print
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evals0[i]))+";")
	print

	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevects0[i][j])+";")

	print("}")
	print("else if (k0 == 0)")
	print("{")

	if LNr == 3:
		print("\tcomplex wg = std::sqrt((complex)"+_ccode(simplify(omegag_2.subs(k0zero_subs)))+");")

	print
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evalsk0zero[i]))+";")
	print

	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevectsk0zero[i][j])+";")

	print("}")
	print("else if (k1 == 0)")
	print("{")

	if LNr == 3:
		print("\tcomplex wg = std::sqrt((complex)"+_ccode(simplify(omegag_2.subs(k1zero_subs)))+");")

	print
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evalsk1zero[i]))+";")
	print

	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevectsk1zero[i][j])+";")

	print("}")
	print("else")
	print("{")
	print("\tcomplex K2 = "+_ccode(inv_eq_subs[K2])+";")
	print("\tcomplex w = std::sqrt((complex)"+_ccode(omega_2)+");")
	print

	if LNr == 3:
		print("\tcomplex wg = std::sqrt((complex)"+_ccode(omegag_2)+");")
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evals[i]))+";")

	print
	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevects[i][j])+";")

	print("}")
	print

