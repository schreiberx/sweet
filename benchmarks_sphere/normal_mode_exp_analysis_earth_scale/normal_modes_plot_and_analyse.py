#!/usr/bin/env python2

import sys
import os
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
import math


curid = os.path.basename(os.getcwd())
print("*"*80)
print(curid)
print("*"*80)


if len(sys.argv) < 2:
	print("Execute with [executable] [input]")
	sys.exit(1)

ext = 'png'
if len(sys.argv) > 2:
	ext = sys.argv[3]


filename = sys.argv[1]


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


print("Time: "+str(time))
print("g: "+str(g))
print("h: "+str(h))
print("r: "+str(r))
print("f: "+str(f))
	

print("Loading data from "+filename)
data = np.loadtxt(filename, skiprows=0)

rows,cols = data.shape
print("Loaded "+str(rows)+" rows")
print("Loaded "+str(cols)+" cols")

alim=[-0.000001,0.000001]

if cols > 2:
	print("Fatal error, cols > 2!")

if cols == 2:
	print("Assuming complex values => splitting them!")
	data = data[:,0] + data[:,1]*1j

#if True:
if False:
	print("Plotting lambdas")

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

	plt.savefig("output_lambda."+ext, dpi=300)




#if True:
if False:
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

	plt.savefig("output_exp_lambda."+ext, dpi=300)



#num_threshold=1e-10
num_threshold=0
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


swe_modes=[]
for n in range(0, N):
	# compute w in side brackets (equation 47) in John's paper

	w = 0
	if True:
		w = math.sqrt(f*f + n*(n+1)*h/(r*r))
	else:
		from sympy import Symbol, solve
		w = Symbol("w")
		s = solve(w*(w*w-f*f-n*(n+1)*h/(r*r)))
#	print(str(n)+str(": ")+str(w))

	for i in range(0, n+1):
		swe_modes.append(w)
		swe_modes.append(-w)
#		print(w)

	for i in range(0, n+1):
		swe_modes.append(0)

swe_modes.sort()
print("Min/max SWE modes on the f-sphere")
print(min(swe_modes))
print(max(swe_modes))



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

	plt.savefig("output_freq_low."+ext, dpi=300)



if True:
	print("High frequency modes: "+str(len(high_freq_modes)))
	print("High frequency modes (computed): "+str(len(swe_modes)))

	datap = high_freq_modes[:]

	s = 0.7
	fig, ax = plt.subplots(figsize=(10.0*s, 5.0*s))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15)

	ax.grid(linestyle='-', linewidth='0.5', color='grey')

	ax.scatter(
		range(-len(swe_modes)/2, len(swe_modes)/2),
		swe_modes,
		c='black',
		s=7.5,
		alpha=0.3,
		edgecolors='none',
		marker="P"
	)

	ax.scatter(
		range(-len(datap)/2, len(datap)/2),
		datap,
		c='red',
		s=3,
		alpha=1,
		edgecolors='none'
	)

	ax.set_xlabel('Eigenmode ID')
	ax.set_ylabel('Eigenvalue (imaginary)')

	if len(datap) > 0:
		plt.ylim([min(min(datap), min(swe_modes))*1.1, max(max(datap), max(swe_modes))*1.1])


#	if len(sys.argv) >= 2:
#		plt.title("High Frequencies - "+sys.argv[2], fontsize=10)

	legend_labels=['Analytical Eigenvalues for f-sphere', 'Numerical Eigenvalues']
	leg = plt.legend(legend_labels, ncol=1, loc='lower right')
	leg.get_frame().set_alpha(1) 

	plt.savefig("output_freq_high."+ext, dpi=300)



#if True:
if False:
	datap = high_freq_modes[:]

	print("modes a: "+str(len(datap)))
	print("modes b: "+str(len(swe_modes)))

	s = 0.7
	fig, ax = plt.subplots(figsize=(10.0*s, 5.0*s))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15)

	ax.grid(linestyle='-', linewidth='0.5', color='grey')

	datap = [swe_modes[i]-datap[i] for i in range(0, len(swe_modes))]
	ax.scatter(
		range(-len(datap)/2, len(datap)/2),
		datap,
		#datap,
		c='red',
		s=3,
		alpha=1,
		edgecolors='none'
	)

	ax.set_xlabel('Eigenmode ID')
	ax.set_ylabel('Eigenvalue (imaginary)')

	if len(datap) > 0:
		plt.ylim([min(datap)*1.1, max(datap)*1.1])

	if len(sys.argv) >= 3:
		plt.title("DIFF "+sys.argv[2], fontsize=10)

	legend_labels=['Difference of Eigenvalues for f-sphere']
	leg = plt.legend(legend_labels, ncol=1, loc='lower right')
	leg.get_frame().set_alpha(1) 

	plt.savefig("output_freq_high_diff."+ext, dpi=300)

