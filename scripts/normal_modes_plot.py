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
data = np.loadtxt(filename, skiprows=0)

rows,cols = data.shape
print("Loaded "+str(rows)+" rows")
print("Loaded "+str(cols)+" cols")

if cols > 2:
	print("Fatal error, cols > 2!")

if cols == 2:
	print("Assuming complex values => splitting them!")
	data = data[:,0] + data[:,1]*1j


#print("EValues:")
#print(data)

print("Plotting EValues")

fig, ax = plt.subplots()
ax.grid(True)
ax.scatter(
	data.real,
	data.imag,
	c='blue',
	s=3,
	alpha=0.3,
	edgecolors='blue'
)


if len(sys.argv) >= 3:
	plt.title("NORMAL MODES - "+sys.argv[2], fontsize=12)

plt.savefig("output.png", dpi=300)
plt.close()


