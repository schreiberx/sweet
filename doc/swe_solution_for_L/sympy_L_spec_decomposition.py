#! /usr/bin/python2

from sympy.matrices import Matrix
from sympy import *
import sys
from sympy.printing import print_ccode
from sympy.abc import u, v, h
import re as re_exp

# use real assumptions for complex valued optimizations (conjugate)
x = symbols('x1:3', real=True, positive=True)
k = symbols('k1:3', real=True, positive=True)


LNr=3

if LNr == 1 or LNr == 2:
	F = symbols('F', real=True, positive=True)
	Fsinv = 1/sqrt(F)
	#Fsinv = symbols('F^-1/2', real=True, positive=True)

elif LNr == 3:
	H0 = symbols('H0', real=True, positive=True)
	g = symbols('g', real=True, positive=True)
	f = symbols('f', real=True, positive=True)

K2 = symbols('K2', real=True, positive=True)


if LNr == 1:
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


# Substitutions
if LNr == 1:
	eq_subs = {}
	inv_eq_subs = {}

	w = symbols('w', real=True, positive=True)
	omega_2 = 1 + Fsinv**2*k[0]*k[0] + Fsinv**2*k[1]*k[1]
	eq_subs[omega_2] = w*w
	eq_subs[-omega_2] = -w*w
	eq_subs[expand(omega_2*F)] = w*w*F
	eq_subs[expand(-omega_2*F)] = -w*w*F
	inv_eq_subs[w] = sqrt(omega_2)

elif LNr == 2:
	eq_subs = {}
	inv_eq_subs = {}

	w = symbols('w', real=True, positive=True)
	omega_2 = 1 + Fsinv**2*k[0]*k[0] + Fsinv**2*k[1]*k[1]
	eq_subs[omega_2] = w*w
	eq_subs[-omega_2] = -w*w
	eq_subs[expand(omega_2*F)] = w*w*F
	eq_subs[expand(-omega_2*F)] = -w*w*F
	inv_eq_subs[w] = sqrt(omega_2)
else:
	eq_subs = {}
	inv_eq_subs = {}

	w = symbols('w', real=True, positive=True)
	omega_2 = H0*g*k[0]*k[0]+H0*g*k[1]*k[1]+f*f
	eq_subs[omega_2] = w*w
	eq_subs[-omega_2] = -w*w
	inv_eq_subs[w] = sqrt(omega_2)

	wg = symbols('wg', real=True, positive=True)
	omegag_2 = g*g*k[0]*k[0]+g*g*k[1]*k[1]+f*f
	eq_subs[omegag_2] = wg*wg
	eq_subs[-omegag_2] = -wg*wg
	inv_eq_subs[wg] = sqrt(omegag_2)



# replace k1*k1 + k2*k2 by 'K2' symbol
eq_subs[k[0]*k[0]+k[1]*k[1]] = K2
inv_eq_subs[K2] = k[0]*k[0]+k[1]*k[1]


U_spec = U_spec.transpose()

U_fromSpec = exp(I*(k[0]*x[0]+k[1]*x[1])) * U_spec


print
print
print("L(U):")
M = getL(U_fromSpec)
pprint(M)


#Lik = Matrix[
Lik = Matrix.zeros(3,3)
Lik[:,0] = M[:,0]/U_fromSpec[0]
Lik[:,1] = M[:,1]/U_fromSpec[1]
Lik[:,2] = M[:,2]/U_fromSpec[2]

if True:
	print
	print
	print("Lik(U):")
	pprint(Lik)

print
print("Using omega(k)^2:")
pprint(omega_2)
if LNr == 3:
	print("Using omegag(k)^2:")
	pprint(omegag_2)
print


def getEVStuff(Lik):
	e_vals_vects = Lik.eigenvects()
	if False:
		print
		print
		print("Eigenvectors:")
		pprint(e_vals_vects)


	evals = [simplify(e_vals_vects[i][0].subs(eq_subs)) for i in range(0,3)]
	multi = [simplify(e_vals_vects[i][1]) for i in range(0,3)]
	evects = [Matrix([simplify(e_vals_vects[i][2][0][j].subs(eq_subs)).subs(eq_subs) for j in range(0,3)]) for i in range(0,3)]

	norm_facs = [simplify(sqrt(evects[i].dot(evects[i].conjugate()).subs(eq_subs)).subs(eq_subs)).subs(eq_subs) for i in range(0,3)]
	nevects = [simplify(evects[i]/norm_facs[i]) for i in range(0,3)]

	#pprint(expand(norm_facs[1]))
	#sys.exit(1)

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
			print
			print("Normalization denominator value (1/value):")
			print(pretty(norm_facs[i]))
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
				if simplify(dot-dot_result) != 0:
					print "VALIDATION FAILED (orthonormality of EVs)"
					if LNr != 3:
						sys.exit(1)

					print "*"*80
					print "* FALLBACK TEST with (g=2, H0=1)"
					dot = simplify(dot.subs({g:1,H0:1,k[0]:3, k[1]:4, f:1}))
					if simplify(dot-dot_result) != 0:
						print "VALIDATION FAILED (orthonormality of EVs)"
						pprint(simplify(expand(dot*dot)))
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


kzero_subs = {k[0]:0, k[1]:0, K2:0}

(evals, nevects) = getEVStuff(Lik)
(evals0, nevects0) = getEVStuff(Lik.subs(kzero_subs))


#
# function to make code nice looking
#
def code_nice_replace(code):
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 2\)' , r'\1*\1', code)
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 3\)' , r'\1*\1*\1', code)
	code = re_exp.sub(r'pow\(([_a-zA-Z0-9]*), 4\)' , r'\1*\1*\1*\1', code)
	code = re_exp.sub(r'sqrt\(([^\)])\)' , r'sqrt((double)\1)', code)
	code = re_exp.sub(r'1\.0L' , r'1.0', code)
	code = re_exp.sub(r'2\.0L' , r'2.0', code)
	return code

def _ccode(expr):
	return code_nice_replace(ccode(expr))


if True:
	print("*"*80)
	print(" C CODE OUTPUT")
	print("*"*80)
	print
	print("std::complex<double> eigenvalues[3];")
	print("std::complex<double> eigenvectors[3][3];")
	print("if (k1 != 0 || k2 != 0)")
	print("{")

	print("\tdouble K2 = k1*k1+k2*k2;")
	print("\tdouble w = std::sqrt("+_ccode(omega_2)+");")
	print
	if LNr == 3:
		print("double omegag = std::sqrt("+_ccode(omegag_2)+");")
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evals[i]))+";")
	print
	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevects[i][j])+";")
	print("}")
	print("else")
	print("{")

	if LNr == 3:
		print("double omegag0 = std::sqrt("+_ccode(omegag_2.subs(kzero_subs))+");")
	print
	for i in range(0,3):
		print("\teigenvalues["+str(i)+"] = "+_ccode(im(evals0[i]))+";")
	print
	for i in range(0,3):
		for j in range(0,3):
			print("\teigenvectors["+str(i)+"]["+sstr(j)+"] = "+_ccode(nevects0[i][j])+";")
	print("}")
	print

