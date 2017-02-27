#!/usr/bin/env python2

import sys
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt


if len(sys.argv) < 2:
	print("Execute with [executable] [input real] [input imag]")
#	print("Execute with [executable] [input] [output_prefix]")
	sys.exit(1)

filename = sys.argv[1]
print("Loading data from "+filename)
data = np.loadtxt(filename, skiprows=0)
print(data)

filename = sys.argv[2]
print("Loading data from "+filename)
data2 = np.loadtxt(filename, skiprows=0)
print(data2)

data = data + data2*1j

print("EValues:")
print(data)

print("Plotting EValues")

fig, ax = plt.subplots()
ax.grid(True)
ax.scatter(
	data.real,
	data.imag,
	c='black',
	s=3,
	alpha=1.0, edgecolors='black'
)

plt.savefig("output.png", dpi=300)
plt.close()
