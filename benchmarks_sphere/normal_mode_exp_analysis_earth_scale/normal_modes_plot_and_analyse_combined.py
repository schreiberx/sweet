#!/usr/bin/env python2

import sys
import os
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math


curid = os.path.basename(os.getcwd())
print("*"*80)
print(curid)
print("*"*80)


if len(sys.argv) < 3:
	print("Execute with [executable] [output] [reference] [comparison]")
	sys.exit(1)

output_file = sys.argv[1]
ref_file = sys.argv[2]
comp_files = sys.argv[3:]

ext = 'png'
if 'pdf' in output_file:
	ext = 'pdf'

class Params:
	time = -1
	g = -1
	f = -1
	h = -1
	r = -1


def load_data(filename):
	p = Params()

	fileh=open(filename)

	while True:
		headerstr = fileh.readline()

		if headerstr[0] != '#':
			break

		if headerstr[2] == 't':
			p.time = float(headerstr[4:])

		elif headerstr[2] == 'g':
			p.g = float(headerstr[4:])

		elif headerstr[2] == 'h':
			p.h = float(headerstr[4:])

		elif headerstr[2] == 'r':
			p.r = float(headerstr[4:])

		elif headerstr[2] == 'f':
			p.f = float(headerstr[4:])

		else:
			print("ERROR: Unknown tag "+headerstr[2])
			sys.exit(1)

	fileh.close()

	if p.time == -1:
		print("Warning: time meta information not found")
		print(filename)
		sys.exit(1)

	if p.g == -1:
		print("Warning: g meta information not found")
		print(filename)
		sys.exit(1)

	if p.h == -1:
		print("Warning: h meta information not found")
		print(filename)
		sys.exit(1)

	if p.r == -1:
		print("Warning: r meta information not found")
		print(filename)
		sys.exit(1)

	if p.f == -1:
		print("Warning: f meta information not found")
		print(filename)
		sys.exit(1)

	data_ref = np.loadtxt(filename, skiprows=0)

	rows,cols = data_ref.shape
	print("Loaded "+str(rows)+" rows")
	print("Loaded "+str(cols)+" cols")

	if cols > 2:
		print("Fatal error, cols > 2!")

	if cols == 2:
		print("Assuming complex values => splitting them!")
		data_ref = data_ref[:,0] + data_ref[:,1]*1j
	
	return (data_ref, p)


(data_ref, p) = load_data(ref_file)
data_ref = data_ref.imag

data_ref.sort()

if 'fsphere1' in ref_file:
	# compute analytical Eigenvalues for fsphere

	# determine truncation
	k=len(data_ref)
	print(k)
	k=k/3
	print(k)

	N = 0
	t = 0
	while t < k:
		t += N+1
		N = N+1

	if t != k:
		print("ERROR: Unable to compute truncation")
		sys.exit(1)

	f = p.f
	r = p.r
	h = p.h

	data_ref=[]
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
			data_ref.append(w)
			data_ref.append(-w)
	#		print(w)

		for i in range(0, n+1):
			data_ref.append(0)

	data_ref.sort()


print("Time: "+str(p.time))
print("g: "+str(p.g))
print("h: "+str(p.h))
print("r: "+str(p.r))
print("f: "+str(p.f))


# markers
markers = []
for m in Line2D.markers:
    try:
        if m != ' ':
            markers.append(m)
    except TypeError:
        pass

print markers


if True:
	s = 0.7
	fig, ax = plt.subplots(figsize=(10.0*s, 5.0*s))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15)

	ax.grid(linestyle='-', linewidth='0.5', color='grey')

	num_plots = len(comp_files)+1

	colormap = plt.cm.gist_ncar
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

	legend_labels=[]
	min_val = 0
	max_val = 0

	if False:
		plt.plot(
			range(-len(data_ref)/2, len(data_ref)/2),
			data_ref,
			linewidth=0.5,
			linestyle='-'#,
	#			marker=marker
		)
		legend_labels.append('reference')

	plt.yscale("log", nonposy='clip')

	ik = 4
	for comp_file in comp_files:
		ik = ik+1

		print("Loading "+comp_file)
		(data_cmp, p_cmp) = load_data(comp_file)
		data_cmp = data_cmp.imag
		data_cmp.sort()

		if len(data_ref) != len(data_cmp):
			print("Number of Eigenvalues don't match")
			sys.exit(1) 

		datap = [data_ref[i]-data_cmp[i] for i in range(0, len(data_ref))]
#		datap = [data_cmp[i] for i in range(0, len(data_ref))]

		min_val = min(min_val, min(datap))
		max_val = max(max_val, max(datap))

		tmp = comp_file[:]

		tmp = tmp.replace('_rexim00000000', '')
		tmp = tmp.replace('_rexim0000000', '_M')
		tmp = tmp.replace('_rexim000000', '_M')
		tmp = tmp.replace('_rexim00000', '_M')
		tmp = tmp.replace('_rexim0000', '_M')

		tmp = tmp.replace('script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere1_', '')
		tmp = tmp.replace('script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere0_', '')
		tmp = tmp.replace('_t-0000001_o000.0001', '')
		tmp = tmp.replace('C0000', 'C')
		tmp = tmp.replace('_Tsm00', '_TM')
		tmp = tmp.replace('_tso0', '')
		tmp = tmp.replace('_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv', '')

		legend_labels.append(tmp)


		if True:
			marker = markers[(ik-1) % len(markers)]
			plt.plot(
				range(-len(datap)/2, len(datap)/2),
				[abs(i) for i in datap],
				linewidth=1,
				linestyle='-',
				marker=marker,
				markevery=30,
				markersize=4
			)



		if False:
			ax.scatter(
				range(-len(datap)/2, len(datap)/2),
				datap,
				s=2,
				alpha=1,
				edgecolors='none'
			)

	ax.set_xlabel('Eigenmode ID')
	ax.set_ylabel('Max. error on Eigenvalues')

	if len(datap) > 0:
		plt.ylim([min_val*2.0, max_val*2.0])

	leg = plt.legend(legend_labels, ncol=1, loc='lower right', fontsize=6)
	leg.get_frame().set_alpha(1) 

	plt.savefig(output_file, dpi=300)

