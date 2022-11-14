#! /usr/bin/env python3

from sympy import *

init_printing()

D = Symbol('D')
G = Symbol('G')

M = Matrix([[0,G],[D,0]])
print("")
print("M:")
pprint(M)

P,D = M.diagonalize()
Pinv = P**-1

P = simplify(P)
D = simplify(D)
Pinv = simplify(Pinv)

print("")
print("P:")
pprint(P)

print("")
print("D:")
pprint(D)

print("")
print("Pinv:")
pprint(P**-1)

print("")
print("R (should be 0):")
R = P*D*Pinv - M
pprint(simplify(R))
