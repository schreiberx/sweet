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

rows,cols = data.shape
print("Loaded "+str(rows)+" rows")
print("Loaded "+str(cols)+" cols")

if data.shape[0] != data.shape[1]:
	print("Assuming complex values => splitting them!")
	data = data[:rows/2] + data[rows/2:]*1j

if data.shape[0] != data.shape[1]:
	print("Fatal error, non-square sized matrix!")

if 1:
	print("Computing EV decomposition with numpy")
	w, v = np.linalg.eig(data)
else:
	print("Computing EV decomposition with multiprecision")
	mp.prec = 53	# default: 53
	w, v = mp.eig(mp.matrix(data))
	w = np.matrix(w)
	v = np.matrix(v)	# TODO: fix this conversion


f = filename+"_evalues_complex.csv"
print("Writing data to "+f)
# simply append the imaginary parts as columns at the end of the real parts!
np.savetxt(f, np.column_stack([w.real, w.imag]), delimiter='\t')

if False:
	f = filename+"_evectors_complex.csv"
	print("Writing data to "+f)
	np.savetxt(f, np.column_stack([v.real, v.imag]), delimiter='\t')

