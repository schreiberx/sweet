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


if len(sys.argv) < 2:
	print("Execute with [executable] [reference_file] [output] [dispersion files]")
	sys.exit(1)


plot_reference = True

if 'dont_plot_reference' == sys.argv[1]:
	plot_reference = False


output_file = sys.argv[2]
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
		sys.exit(1)

	if cols == 2:
		print("Assuming complex values => splitting them!")
		data_ref = data_ref[:,0] + data_ref[:,1]*1j

	return (data_ref, p)



def getTruncation(
	k	# number of reference Eigenvalues (3x number of DoFs)
)	:
	k=k/3

	N = 0
	t = 0
	while t < k:
		t += N+1
		N = N+1

	if t != k:
		print("ERROR: Unable to compute truncation")
		sys.exit(1)

	return N



# get truncation
(data_dummy, p) = load_data(comp_files[0])
data_dummy = data_dummy.imag

num_evalues = len(data_dummy)
N = getTruncation(num_evalues)


# compute reference data
data_ref=[]
if 'fsphere1' in comp_files[0]:
	# compute analytical Eigenvalues for fsphere

	f = p.f
	r = p.r
	h = p.h

	for n in range(0, N):
		# compute w in side brackets (equation 47) in John's paper

		w = 0
		if True:
			w = math.sqrt(f*f + n*(n+1)*h/(r*r))
		else:
			from sympy import Symbol, solve
			w = Symbol("w")
			s = solve(w*(w*w-f*f-n*(n+1)*h/(r*r)))

		for i in range(0, n+1):
			data_ref.append(w)
			data_ref.append(-w)

		for i in range(0, n+1):
			data_ref.append(0)

	data_ref.sort()

else:
	print("TODO")
	sys.exit(1)

##################################################################
##################################################################


# markers
markers = []
for m in Line2D.markers:
    try:
        if m != ' ':
            markers.append(m)
    except TypeError:
        pass

print markers



##################################################################
##################################################################

print("Time: "+str(p.time))
print("g: "+str(p.g))
print("h: "+str(p.h))
print("r: "+str(p.r))
print("f: "+str(p.f))

##################################################################
##################################################################
##################################################################


epsilon=1e-14

#
# Relative dispersion error for single stage two-level schemes
# (See Durran, Numerical Methods or Fluid Dynamics, page 46)
# parameter a selects the type of two-level scheme:
#   0: forward euler
#   1: backward euler
#   0.5: CN (trapezoidal)
#
def RSingle(	omega,	# wave speed
		dt,	# time step size
		a	# alpha coefficient
):
	if abs(omega) < epsilon:
		return 1.0

	retval = 1.0/(omega*dt)*math.atan( (omega*dt) / (1.0-a*(1.0-a)*((omega*dt)**2.0)))

	return retval


def RSingleRK2(	omega,	# wave speed
		dt	# time step size
):
	if abs(omega) < epsilon:
		return 1.0

	retval = 1.0/(omega*dt)*math.atan( (omega*dt) / (1.0-0.5*((omega*dt)**2.0)))

	return retval


if True:
	s = 0.7
	fig, ax = plt.subplots(figsize=(10.0*s, 5.0*s))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15)

	ax.grid(linestyle='-', linewidth='0.5', color='grey')

	alpha_vals = [0.0, 1.0, 0.5]	# alpha values for relative dispersion shifts (see RSingle for a description)

	# number of total plots
	num_plots = len(alpha_vals)

	colormap = plt.cm.gist_ncar
	#plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
	colors = [colormap(i) for i in np.linspace(0, 0.9, num_plots)]

	legend_labels=[]
#	plt.yscale("log", nonposy='clip')

	# timestep numbers
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

		tmp = comp_file[:]

		tmp = tmp.replace('_rexim00000000', '')
		tmp = tmp.replace('_rexim0000000', '_M')
		tmp = tmp.replace('_rexim000000', '_M')
		tmp = tmp.replace('_rexim00000', '_M')
		tmp = tmp.replace('_rexim0000', '_M')

		tmp = tmp.replace('script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid1_fsphere1_', '')
		tmp = tmp.replace('script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid1_fsphere0_', '')
		tmp = tmp.replace('_t-0000001_o000.0001', '')
		tmp = tmp.replace('C0000', 'C')
		tmp = tmp.replace('_tsm00', '_TM')
		tmp = tmp.replace('_tsm', '_TM')
		tmp = tmp.replace('_tso0', '')
		tmp = tmp.replace('_rexih0.15', '')
		tmp = tmp.replace('_rexihalf0', '')
		tmp = tmp.replace('_rexihalf1', '')
		tmp = tmp.replace('_rexiextmodes02', '')
		tmp = tmp.replace('_rexiextmodes04', '')
		tmp = tmp.replace('_rexiprealloc0', '')
		tmp = tmp.replace('/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv', '')
		tmp = tmp.replace('/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv', '')

		legend_labels.append(tmp)


		datax = []
		datap = []

		for i in range(len(data_cmp)):
			if 'TM_l_rexi' in tmp:
				if abs(data_cmp[i]) < epsilon:
					continue

			datax.append(data_ref[i])

			rel_diff = (data_cmp[i]-data_ref[i])/data_ref[i]+1.0

			# Filter out this point
			if not 'TM_l_rexi' in tmp:
				if abs(data_cmp[i]) < epsilon:
					rel_diff = 1.0

			datap.append(rel_diff)

		marker = markers[(ik-1) % len(markers)]
		color = colors[(ik-1) % len(colors)]
		plt.plot(
			datax,
			datap,
			linewidth=1,
			linestyle='-',
			marker=marker,
			markevery=10,
			markersize=4,
			color = color
		)



	ref_linewidth = 0.5
	ref_linestyle = '--'
	ref_color = 'k'
	ref_marker = ''
	ref_markevery = 10
	ref_markersize = 2

	if plot_reference:
		legend_labels.append("Analytical dispersion errors for RK1/2 and CN")
		for alpha in alpha_vals:
			ik = ik+1

			print("Plotting relative dispersion for alpha "+str(alpha))
			#legend_labels.append("Dispersion with alpha="+str(alpha))

			#
			# generate graph (datax, datap) for
			# omega vs. omega+relative error
			#
			datax = data_ref

			dt = p.time
			datap = [RSingle(i, dt, alpha) for i in data_ref]

#			ref_marker = markers[(ik-1) % len(markers)]
			plt.plot(
				datax,
				datap,
				linewidth=ref_linewidth,
				linestyle=ref_linestyle,
				color=ref_color,
				marker=ref_marker,
				markevery=ref_markevery,
				markersize=ref_markersize
			)
		
		if True:
			ik = ik+1

			print("Plotting relative dispersion for alpha "+str(alpha))
			#legend_labels.append("Dispersion for RK2")

			#
			# generate graph (datax, datap) for
			# omega vs. omega+relative error
			#
			datax = data_ref

			dt = p.time
			datap = [RSingleRK2(i, dt) for i in data_ref]

			plt.plot(
				datax,
				datap,
				linewidth=ref_linewidth,
				linestyle=ref_linestyle,
				color=ref_color,
				marker=ref_marker,
				markevery=ref_markevery,
				markersize=ref_markersize
			)





	ax.set_xlabel('Dispersion speed')
	ax.set_ylabel('Error in dispersion')

	leg = plt.legend(legend_labels, ncol=1, loc='lower right', fontsize=6)
	leg.get_frame().set_alpha(1) 

	plt.savefig(output_file, dpi=300)

