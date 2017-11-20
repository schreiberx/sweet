#! /usr/bin/env python3

from sympy.matrices import Matrix
from sympy import *
from sympy.printing import print_ccode
from sympy.abc import u, v, h
import re as re_exp
import sys

a11, b11 = symbols('a11 b11', real=True)
a12, b12 = symbols('a12 b12', real=True)
a13, b13 = symbols('a13 b13', real=True)
a21, b21 = symbols('a21 b21', real=True)
a22, b22 = symbols('a22 b22', real=True)
a23, b23 = symbols('a23 b23', real=True)
a31, b31 = symbols('a31 b31', real=True)
a32, b32 = symbols('a32 b32', real=True)
a33, b33 = symbols('a33 b33', real=True)


M = Matrix([[a11+I*b11, a12+I*b12, a13+I*b13], [a21+I*b21, a22+I*b22, a23+I*b23], [a31+I*b31, a32+I*b32, a33+I*b33]])

A = M
print(A)
pprint(A)
print()
evals = A.eigenvals()
print(evals)

#evects = A.eigenvects()
#print(evects)
