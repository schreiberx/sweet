#!/usr/bin/env python3
#
#   normal_mode_analysis
#
# Read normal mode (mode wise) output from sweet
#
#  Pedro Peixoto (pedrosp@ime.usp.br)
#
#--------------------------------------------------

import sys
import numpy as np
from mpmath import mp
import math
import re
import os
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#Check arguments
if len(sys.argv) < 2:
	print("Execute with [executable] [input]")
	sys.exit(1)


#File with the complex modes 
filename = sys.argv[1]
name, ext = os.path.splitext(filename)
print("Loading data from "+filename)

#load data
data = np.loadtxt(filename, skiprows=0)

#split real and imag parts
#print(data)
data_real = data[:, 0::2]
#print(data_real)
data_imag = data[:, 1::2]
#print(data_imag)

#get dimensions
rows,cols = data_real.shape
print("Loaded "+str(rows)+" rows")
print("Loaded "+str(cols)+" cols")

#mode indexes
#x = np.arange(0, 1, 1/cols)
#y = np.arange(0, 1, 1/rows)
x = np.arange(0, cols, 1)
y = np.arange(0, rows, 1)
X, Y = np.meshgrid(x, y)
print(X)
#plot

plt.figure()
CS = plt.contour(X, Y, data_imag)
plt.clabel(CS, inline=1, fontsize=10)
plt.title(name)
plt.xlabel('x wavenumber')
plt.ylabel('y wavenumber')

plt.savefig(name+str("imag.png"), dpi=300)


quit()

alim=[-0.000001,0.000001]

if True:
	print("Plotting Real parts of omega")

	datap = data[:]

	fig, ax = plt.subplots()
	ax.grid(True)
	ax.scatter(
		datap.real,
		datap.imag,
		c='blue',
		s=3,
		alpha=0.3,
		edgecolors='none'
	)
	ax.set_xlabel('Re(lambda)')
	ax.set_ylabel('Im(lambda)')

	minx = min(datap.real)*1.1
	maxx = max(datap.real)*1.1

	miny = min(datap.imag)*1.1
	maxy = max(datap.imag)*1.1

	if abs(minx-maxx) > 1e-9:
		plt.xlim([minx, maxx])

	if abs(miny-maxy) > 1e-9:
		plt.ylim([miny, maxy])

	if len(sys.argv) >= 3:
		plt.title("Eigenvalues lambda - "+sys.argv[2], fontsize=10)

	plt.savefig("output_lambda.png", dpi=300)




if True:
	print("Plotting exp(lambda)")
	#
	# INFORMATION
	#
	# If the eigenvalues are too small, nothing will be plotted!
	#

	datap = np.exp(data)

#	for i in range(0, len(datap)):
#		print(str(datap[i].real)+"\t"+str(datap[i].imag))

	fig, ax = plt.subplots()

	ax.set_xlabel('Re(exp(lambda))')
	ax.set_ylabel('Im(exp(lambda))')

	minx = min(datap.real)*1.1
	maxx = max(datap.real)*1.1

	miny = min(datap.imag)*1.1
	maxy = max(datap.imag)*1.1

	if abs(minx-maxx) > 1e-9:
		plt.xlim([minx, maxx])
	else:
		plt.xlim(alim)

	if abs(miny-maxy) > 1e-9:
		plt.ylim([miny, maxy])
	else:
		plt.ylim(alim)

	if len(sys.argv) >= 3:
		plt.title("Eigenvalues exp(lambda) - "+sys.argv[2], fontsize=10)

	ax.scatter(
		datap.real,
		datap.imag,
		c='blue',
		s=3,
		alpha=0.3,
		edgecolors='none'
	)
	ax.grid(True)

	plt.savefig("output_exp_lambda.png", dpi=300)



num_threshold=1e-10
#lowhigh_threshold=1e-4
lowhigh_threshold=0

#
# Compute artificial solution for f-plane
#
f_plane_modes = []

# determine truncation
k=len(data)
print(k)
k=k/3
print(k)

N = 0
t = 0
while t < k:
	t += N+1
#	print("T: "+str(t))
	N = N+1

if t != k:
	print("ERROR: Unable to compute truncation")
	sys.exit(1)

non_geostrophic_modes=[]
for n in range(0, N):
	# compute w in side brackets (equation 47) in John's paper
	f = 0.00014584
	a = 6371220.0
	h = 100000.0

	w = 0
	if True:
		w = math.sqrt(f*f + n*(n+1)*h/(a*a))
	else:
		from sympy import Symbol, solve
		w = Symbol("w")
		s = solve(w*(w*w-f*f-n*(n+1)*h/(a*a)))
#	print(str(n)+str(": ")+str(w))

	for i in range(0, n+1):
		non_geostrophic_modes.append(w)
		non_geostrophic_modes.append(-w)

non_geostrophic_modes.sort()



stat_modes=[]
high_freq_modes=[]
low_freq_modes=[]

datai = data.imag
datai.sort()
#print(datai)

for d in datai:
	ad = abs(d)
	if ad < num_threshold:
		stat_modes.append(d)
	elif ad < lowhigh_threshold:
		low_freq_modes.append(d)
	else:
		high_freq_modes.append(d)



if True:

	print("Stationary modes: "+str(len(stat_modes)))
	#print(stat_modes)


if lowhigh_threshold != 0:
	print("Low frequency modes: "+str(len(low_freq_modes)))
	#print(low_freq_modes)

	datap = low_freq_modes

	fig, ax = plt.subplots()
	ax.grid(True)
	ax.scatter(
		range(-len(datap)/2, len(datap)/2),
		datap,
		c='blue',
		s=3,
		alpha=0.3,
		edgecolors='none'
	)
	ax.set_xlabel('Modes')
	ax.set_ylabel('Eigenvalue (imaginary)')
	if len(datap) > 0:
		plt.ylim([min(datap)*1.1, max(datap)*1.1])


	if len(sys.argv) >= 3:
		plt.title("Low Frequencies - "+sys.argv[2], fontsize=10)

	plt.savefig("output_freq_low.png", dpi=300)


if True:
	print("High frequency modes: "+str(len(high_freq_modes)))
	print("High frequency modes (computed): "+str(len(non_geostrophic_modes)))
	#print(high_freq_modes)

	datap = high_freq_modes

	fig, ax = plt.subplots()
	ax.grid(True)


	ax.scatter(
		range(-len(non_geostrophic_modes)/2, len(non_geostrophic_modes)/2),
		non_geostrophic_modes,
		c='blue',
		s=5,
		alpha=0.3,
		edgecolors='none'
	)

	ax.scatter(
		range(-len(datap)/2, len(datap)/2),
		datap,
		c='red',
		s=1,
		alpha=1,
		edgecolors='none'
	)
	ax.set_xlabel('Modes')
	ax.set_ylabel('Eigenvalue (imaginary)')
	if len(datap) > 0:
		plt.ylim([min(datap)*1.1, max(datap)*1.1])


	if len(sys.argv) >= 3:
		plt.title("Low Frequencies - "+sys.argv[2], fontsize=10)

	plt.savefig("output_freq_high.png", dpi=300)


