#! /usr/bin/env python3

from sympy import *

a, b, c, d = symbols('a b c d')
M = Matrix([[a, b], [c, d]])
A = M

evals = A.eigenvals()
print(evals)

evects = A.eigenvects()
print(evects)
