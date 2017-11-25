#! /usr/bin/env python3

from sympy.matrices import Matrix
from sympy import *
from sympy.printing import print_ccode
from sympy.abc import u, v, h
import re as re_exp
import sys

a00 = symbols('a00')
a01 = symbols('a01')
a02 = symbols('a02')
a10 = symbols('a10')
a11 = symbols('a11')
a12 = symbols('a12')
a20 = symbols('a20')
a21 = symbols('a21')
a22 = symbols('a22')


M = Matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])

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



print("*"*80)
print(" C CODE OUTPUT ")
print("*"*80)
print()
print("std::complex<double> eigenvalues[3];")
print()
	
for i in range(0,3):
	print("\t o_eval["+str(i)+"] = "+str(ccode(simplify(ev[i])))+";")

print()



