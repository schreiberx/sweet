#!/usr/bin/env python2

import sys
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt


if len(sys.argv) < 2:
	print("Execute with [executable] [input]")
	sys.exit(1)

filename = sys.argv[1]


time = -1
g = -1
h = -1
r = -1
f = -1


#
# Load meta data
#

new_headers = []
if True:
	fileh=open(filename)

	while True:
		headerstr = fileh.readline()

		if headerstr[0] != '#':
			break

		if headerstr[2] == 't':
			time = float(headerstr[4:])

		elif headerstr[2] == 'g':
			g = float(headerstr[4:])

		elif headerstr[2] == 'h':
			h = float(headerstr[4:])

		elif headerstr[2] == 'r':
			r = float(headerstr[4:])

		elif headerstr[2] == 'f':
			f = float(headerstr[4:])

		else:
			print("ERROR: Unknown tag "+headerstr[2])
			sys.exit(1)

		if headerstr[-1] == '\n':
			headerstr = headerstr[1:-1]

		new_headers.append(headerstr[1:])


	fileh.close()

	if time == -1:
		print("time meta information not found")
		sys.exit(1)

	if g == -1:
		print("g meta information not found")
		sys.exit(1)

	if h == -1:
		print("h meta information not found")
		sys.exit(1)

	if r == -1:
		print("r meta information not found")
		sys.exit(1)

	if f == -1:
		print("f meta information not found")
		sys.exit(1)

print("Time: "+str(time))
print("g: "+str(g))
print("h: "+str(h))
print("r: "+str(r))
print("f: "+str(f))
	
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


f = filename+"_evalues_orig_complex.csv"
print("Writing data to "+f)
# simply append the imaginary parts as columns at the end of the real parts!
np.savetxt(f, np.column_stack([w.real, w.imag]), delimiter='\t', header='\n'.join(new_headers))


we = np.log(w)/time


print("Headers: "+str(new_headers))

f = filename+"_evalues_complex.csv"
print("Writing data to "+f)
# simply append the imaginary parts as columns at the end of the real parts!
np.savetxt(f, np.column_stack([we.real, we.imag]), delimiter='\t', header='\n'.join(new_headers))


if False:
	f = filename+"_evectors_complex.csv"
	print("Writing data to "+f)
	np.savetxt(f, np.column_stack([v.real, v.imag]), delimiter='\t', header='\n'.join(new_headers))

