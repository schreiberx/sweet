#! /usr/bin/env python3

from sympy.matrices import Matrix
from sympy import *
from sympy.printing import print_ccode
from sympy.abc import u, v, h
import re as re_exp
import sys

a11 = symbols('a11')
a12 = symbols('a12')
a13 = symbols('a13')
a21 = symbols('a21')
a22 = symbols('a22')
a23 = symbols('a23')
a31 = symbols('a31')
a32 = symbols('a32')
a33 = symbols('a33')


M = Matrix([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])

A = M
print(A)
pprint(A)
print()
evals = A.eigenvals()
print(evals)
pprint(evals)

ev=list(evals.keys())
print(ev)
#evects = A.eigenvects()
#print(evects)

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



print("*"*80)
print(" C CODE OUTPUT ")
print("*"*80)
print()
print("std::complex<double> eigenvalues[3];")
print()
	
for i in range(0,3):
	print("\t eigenvalues["+str(i)+"] = "+str(ev[i])+";")

print()



