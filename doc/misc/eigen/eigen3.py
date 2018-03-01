#! /usr/bin/env python3

from sympy import *

a, b, c, d, e, f, g, h, i = symbols('a b c d e f g h i')
M = Matrix([[a, b, c], [d, e, f], [g, h, i]])
A = M

evals = A.eigenvals()
print(evals)

evects = A.eigenvects()
print(evects)
