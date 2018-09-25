#! /usr/bin/env python3

import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D
import os


datafiles = sys.argv[1:]



#################################################################################
#################################################################################
#################################################################################

fig, ax = plt.subplots(figsize=(8,5))


#ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
ax.set_ylim(1e-1, 3000)


output_filename = "output_rexi_performance_breakdown.pdf"


#
# All identifiers from the output and their placeholders
#
timings_identifiers = {
	'main': ' + SimulationBenchmarkTimings.main: ',
	'main_setup': ' + SimulationBenchmarkTimings.main_setup: ',
	'main_simulationloop': ' + SimulationBenchmarkTimings.main_simulationloop: ',
	'rexi': ' + SimulationBenchmarkTimings.rexi: ',
	'rexi_setup': ' + SimulationBenchmarkTimings.rexi_setup: ',
	'rexi_shutdown': ' + SimulationBenchmarkTimings.rexi_shutdown: ',
	'rexi_timestepping': ' + SimulationBenchmarkTimings.rexi_timestepping: ',
	'rexi_timestepping_solver': ' + SimulationBenchmarkTimings.rexi_timestepping_solver: ',
	'rexi_timestepping_broadcast': ' + SimulationBenchmarkTimings.rexi_timestepping_broadcast: ',
	'rexi_timestepping_reduce': ' + SimulationBenchmarkTimings.rexi_timestepping_reduce: ',
	'rexi_timestepping_miscprocessing': ' + SimulationBenchmarkTimings.rexi_timestepping_miscprocessing: ',
}

#
# Bars to plot
#
plot_bars_and_labels = {
	'main_simulationloop': 'Total wallclock',
#	'main_nonlinearities': 'Non-linearities',	# This is a placeholder and handled in a special way
	'rexi_timestepping': 'Total REXI',
	'rexi_timestepping_miscprocessing': 'REXI misc',
	'rexi_timestepping_broadcast': 'REXI broadcast',
	'rexi_timestepping_solver': 'REXI solver',
	'rexi_timestepping_reduce': 'REXI reduce',
}


# 'name' contains the pretty name of the run
bar_values = {
	'name':  []	# name for this run
}
for key in plot_bars_and_labels:
	bar_values[key] = []


if len(datafiles) == 0:
	print("Please specify datafiles")
	sys.exit(1)


c = 2
for i in datafiles:
	
	filename = i+'/output.out'
	if not os.path.isfile(filename):
		print("File '"+filename+"' doesn't exist")
		continue

	with open(filename) as f:
		lines = f.readlines()

	#
	# prepare temporary timing list
	# this helps to setup dummy data in case that some timing doesn't exist
	#
	timings = {key: -1.0 for key in timings_identifiers}

	valid_found = False
	for l in lines:
		for key, prefix in timings_identifiers.items():
			if l[0:len(prefix)] == prefix:
				timings[key] = float(l[len(prefix):])
				valid_found = True
				break

	if not valid_found:
		continue

	#
	# Append to bar values
	#
	for key in timings:
		# only if bar should be plotted
		if key in bar_values:
			bar_values[key].append(timings[key])


	#
	# Some info output
	#
	name = i
	name = name.replace('script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_rexi_lc_n_erk_ver1_tso2_tsob2_C000360_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time', '')
	name = re.sub(r"C[0-9]*_", "", name)
#	name = name.replace('REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128_', '')
#	name = name.replace('tso2_tsob2_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_', '')
#	name = name.replace('M0128_MPI_space01_', '')
	if name[-1] == '_':
		name = name[:-1]

	name = "MPI ranks = "+str(int(name))
	bar_values['name'].append(name)

	print([bar_values[key][-1] for key in bar_values])

labels_timings = [plot_bars_and_labels[key] for key in plot_bars_and_labels]
values_timings = [bar_values[key] for key in plot_bars_and_labels]

#print(labels_timings)
#print(values_timings)

# number of different runs
num_runs = len(bar_values['name'])

if num_runs == 0:
	print("No valid output found from one run")
	sys.exit(1)


# number of timings
num_timings = len(plot_bars_and_labels)

# adjust depending on number of labels_timings
width = 1.0/num_timings*0.9
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

ax.set_xticks([i-width/2 for i in range(num_timings)])
ax.set_xticklabels(bar_values['name'])
ax.set_ylabel('Wallclock time (seconds)')



for tick in ax.get_xticklabels():
    tick.set_rotation(45)



rects = []

# iterate over all runs
for j in range(num_timings):
	# scatter plot bars for each run in different colors 'colors'
	rects_ = ax.bar([i - width*(num_runs/2-j) for i in range(num_runs)], values_timings[j], width, color=colors[j])
	rects.append(rects_)


for j in range(num_timings):
	for rect, val in zip(rects[j], values_timings[j]):
		height = rect.get_height()
		text =  "%.2f" % height
		ax.text(rect.get_x() + rect.get_width()/2, height+math.log(1.0+0.1*height), text, ha='center', va='bottom', size=8, rotation=90)


ax.legend([rect[0] for rect in rects], labels_timings)
plt.tight_layout()

if output_filename != '':
	plt.savefig(output_filename)
else:
	plt.show()

