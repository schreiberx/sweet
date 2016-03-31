#! /usr/bin/env python3

import math

a_list = [	-0.86304 + 0j,
		-0.86304 + 1j,
		-0.86304 + 10j,
		-0.86304 + 100j,
		-0.86304 + 1000j
	]


f = 1.0
g = 1.0

eta0_hat = 1.0
u0_hat = 1.0
v0_hat = 0.5

eta_bar = 1.0

for a in a_list:
	for k in [0.01, 0.1, 1.0]:
#	for k in [0.1]:
		lhs = (a*a+f*f)+g*eta_bar*(2.0*math.pi*k)**2
		rhs =	(f*f+a*a)/a * eta0_hat - eta_bar *1.0j*2.0*math.pi*k*(u0_hat + v0_hat) - (f * eta_bar / a) *1.0j*2.0*math.pi*k*(v0_hat - u0_hat)

		eta_hat = rhs / lhs
		print("a="+str(a)+", k="+str(k)+"    > "+str(eta_hat))
