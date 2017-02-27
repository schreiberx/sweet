#!/usr/bin/env python2

import sys
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt


if len(sys.argv) < 2:
	print("Execute with [executable] [input]")
	sys.exit(1)

filename = sys.argv[1]
print("Loading data from "+filename)
data = np.loadtxt(filename, skiprows=0).transpose()

if 1:
	print("Computing EV decomposition with numpy")
	w, v = np.linalg.eig(data)
else:
	print("Computing EV decomposition with multiprecision")
	mp.prec = 53	# default: 53
	w, v = mp.eig(mp.matrix(data))
	w = np.matrix(w)
	v = np.matrix(v)	# TODO: fix this conversion

f = filename+"_evalues_real.csv"
print("Writing data to "+f)
np.savetxt(f, w.real, delimiter='\t')

f = filename+"_evalues_imag.csv"
print("Writing data to "+f)
np.savetxt(f, w.imag, delimiter='\t')

f = filename+"_evectors_real.csv"
print("Writing data to "+f)
np.savetxt(f, v.real, delimiter='\t')

f = filename+"_evectors_imag.csv"
print("Writing data to "+f)
np.savetxt(f, v.imag, delimiter='\t')
