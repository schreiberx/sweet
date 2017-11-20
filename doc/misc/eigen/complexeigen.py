#! /usr/bin/env python3

from sympy import *

a11, b11 = symbols('a11 b11', real=True)
a12, b12 = symbols('a12 b12', real=True)
a21, b21 = symbols('a21 b21', real=True)
a22, b22 = symbols('a22 b22', real=True)

M = Matrix([[a11+I*b11, a12+I*b12], [a21+I*b21, a22+I*b22]])
A = M
print(A)
print()
pprint(A)
print()
evals = A.eigenvals()
print(evals)
pprint(evals)
print()
evects = A.eigenvects()
print(evects)
