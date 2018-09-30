#! /usr/bin/env python3

import sys
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from postprocessing_read_output_benchmark import *

#
# First, use
#   ./postprocessing.py > postprocessing_output.txt
# to generate the .txt file
#


# Mode: wallclocktime, dt
metric = None
if len(sys.argv) > 1:
	metric = sys.argv[1]

if metric == None:
	raise Exception("Please specify metric ('dt' or 'wallclocktime') as 1st argument")

# input filename
input_filename = None
if len(sys.argv) > 2:
	input_filename = sys.argv[2]

if input_filename == None:
	raise Exception("Please specify input file name as 2nd argument")

# output filename
output_filename = None
if len(sys.argv) > 3:
	output_filename = sys.argv[3]

if output_filename == None:
	raise Exception("Please specify input file name as 3rd argument")

print("INFO: "+"*"*70)
print("INFO: Metric: "+metric)
print("INFO: Input filename: "+input_filename)
print("INFO: Output filename: "+output_filename)
print("INFO: "+"*"*70)


#################################################################################
#################################################################################
#################################################################################


linestyles = ['-', '--', ':', '-.']

# axis
global ax
ax = None

def plot(x, y, marker, linestyle, label):
	global ax

	# plot values and prev_name
	print("INFO: plot.label "+label)
	print("INFO: plot.x: "+str(x))
	print("INFO: plot.y: "+str(y))

	if len(x) == 0:
		return

	ax.plot(x, y, marker=marker, linestyle=linestyle, label=label)

	px = x[:]
	py = y[:]

	px = px[0::3]
	py = py[0::3]

	if len(px) % 2 == 0:
		px.append(x[-1])
		py.append(y[-1])

	for i, txt in enumerate(px):
		text = "%.1f/%.1f" % (py[i], px[i])

		if metric == 'dt':
			#ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=8)
			ax.annotate(px[i], (px[i]*1.03, py[i]*0.92), fontsize=8)
		elif metric == 'wallclocktime':
			ax.annotate(text, (px[i]*1.03, py[i]*1.03), fontsize=8)
		else:
			raise Exception("Unknown metric "+metric)


print("INFO: Loading data")
solver_groups = postprocessing_read_output_benchmark(input_filename)

values_y = []
values_x = []

fig, ax = plt.subplots(figsize=(10,7))

ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

markers = []
for m in Line2D.markers:
	try:
		if len(m) == 1 and m != ' ' and m != '':
			markers.append(m)

	except TypeError:
		pass

for key, solver_group in solver_groups.items():
	print("INFO: solver_group.name: "+solver_group['name'])


	# setup accumulation buffers
	values_x = []
	values_y = []
	prev_ts_name = ''
	c = 2
	for simdata in solver_group['simdata']:
		if prev_ts_name == '':
			# Initialize
			prev_ts_name = simdata['ts_name']

		elif prev_ts_name != simdata['ts_name']:

			# plot data
			if len(values_x) > 0:
				plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], simdata['ts_name'])

				# empty accumulation buffers
				values_x = []
				values_y = []

			prev_ts_name = simdata['ts_name']

		if metric == 'wallclocktime':
			values_x.append(simdata['wallclocktime'])

		elif metric == 'dt':
			values_x.append(simdata['dt'])
		else:
			raise Exception("Unknown metric")

		values_y.append(simdata['l1'])


		c += 1

	if metric == 'wallclocktime':
		plt.xlabel("Wallclock time")
	elif metric == 'dt':
		plt.xlabel("Timestep size")
	else:
		raise Exception("Unknown metric")

	# plot data
	if len(values_x) > 0:
		plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], simdata['ts_name'])

	plt.legend()
	plt.savefig(output_filename)

sys.exit(1)

c = 2
for l in lines:
	if l[-1] == '\n':
		l = l[:-1]

	# Plot and restart
	if l == 'Running tests for new group:':
		print(" + New group detected")
		if len(values_x) > 0:
			plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)

		prev_name = d[0]
		values_y = []
		values_x = []
		c = c+1
		continue

	d = l.split("\t")

	# Only consider lines starting with "script_" to contain valid data
	line_start = "script_"
	if l[:len(line_start)] != line_start:
		continue

	if len(d) != 5:
		raise Exception("ERROR: should have 5 columns")

	# skip invalid nan's
	if d[1] == 'nan':
		print("NaN value detected, skipping")
		continue

	# Start with empty dictionary
	simdata = {}

	prev_name = d[0]
	# make names more pretty
	# script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_rexi_lc_n_erk_ver1_tso2_tsob2_C000060_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space007_time128
	prev_name = prev_name.replace('script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_', '')
	prev_name = prev_name.replace('tso2_tsob2_', '')
	prev_name = re.sub(r"C[0-9]*_", "", prev_name)
	prev_name = prev_name.replace('REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_', '')

	if prev_name[-1] == '_':
		prev_name = prev_name[:-1]

	print("+ Processing "+prev_name)

	# Setup some data
	simdata['name'] = d[0],
	simdata['l1'] = float(d[1])
	simdata['l2'] = float(d[2])
	simdata['linf'] = float(d[3])
	simdata['wallclocktime'] = float(d[4])

	m = re.search('_C([0-9]*)', d[0])
	simdata['dt'] = float(m.group(1))

	# Always use max error
	values_y.append(simdata['l1'])

	if metric == 'wallclocktime':
		values_x.append(simdata['wallclocktime'])
		plt.xlabel("Wallclock time")

	elif metric == 'dt':
		values_x.append(simdata['dt'])
		plt.xlabel("Timestep size")
	else:
		raise Exception("Unknown metric")


	plt.ylabel("Error")


if len(values_x) > 0:
	plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)

#if metric == 'dt':
	#ax.xaxis.set_ticks([2**i for i in range(0, 10)])
	#ax.yaxis.set_ticks([2**i for i in range(-5, 5)])

plt.legend()

if output_filename != '':
	plt.savefig(output_filename)
else:
	plt.show()

