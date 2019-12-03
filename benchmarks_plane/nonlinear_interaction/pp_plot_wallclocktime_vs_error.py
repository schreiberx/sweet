#! /usr/bin/env python3

import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D


from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

if len(sys.argv) > 1:
	muletag = sys.argv[1]
	output_filename = sys.argv[2]
else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [jobdata mule tag for y axis] [output_filename.pdf] [jobdir1] [jobdir2] ... [jobdirN]")
	print("")
	sys.exit(1)


if len(sys.argv) > 3:
	# Load Jobs specified via program parameters
	jd = JobsData(job_dirs=sys.argv[3:])
else:
	# Load all Jobs
	jd = JobsData()


# Consolidate data...
jdc = JobsDataConsolidate(jd)

# ... which belongs to the same time integration method
jdc_groups = jdc.create_groups(['runtime.timestepping_method'])

#
# Filter to exclude data which indicates instabilities
#
def data_filter(x, y, jd):
	if y == None:
		return True

	if 'runtime.max_simulation_time' in jd:
		if jd['runtime.max_simulation_time'] <= 24*60*60:
			if y > 100:
				return True
		elif jd['runtime.max_simulation_time'] <= 10*24*60*60:
			if y > 1000:
				return True

	return False



# Exctract data suitable for plotting
jdc_groups_data = JobsData_GroupsPlottingScattered(
			jdc_groups,
#			'runtime.timestep_size',
			'output.simulation_benchmark_timings.main_timestepping',
			muletag,
			data_filter=data_filter
		)

data = jdc_groups_data.get_data()


def label(d):
	val = d['runtime.timestepping_method'].replace('_', '\\_')+', $\Delta t = '+str(d['runtime.timestep_size'])+'$'
	#val = d['runtime.timestepping_method'].replace('_', '\\_')+', $\Delta t = '+str(d['output.simulation_benchmark_timings.main_timestepping'])+'$'
	return val



##########################################################
# Plotting starts here
##########################################################

print("*"*80)
print("*"*80)
print("*"*80)


fontsize=18
figsize=(10, 10)

fig, ax = plt.subplots(figsize=(10,10))
#plt.rc('text', usetex=True)

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

linestyles = ['-', '--', ':', '-.']


c = 0

title = ''
for key, d in data.items():
	x = d['x_values']
	y = d['y_values']
	l = key.replace('_', '\\_')

	x = [float(a) for a in x]
	y = [float(a) for a in y]

	xo = []
	yo = []
	for i in range(len(x)):
		if np.isnan(y[i]) or np.isnan(x[i]):
			continue
		xo.append(x[i])
		yo.append(y[i])

	x = xo[:]
	y = yo[:]


	print(" + "+l)
	print("\tx: "+str(x))
	print("\ty: "+str(y))
	ax.plot(x, y, marker=markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=l)

	c = c + 1


if title != '':
	plt.title(title, fontsize=fontsize)

plt.xlabel("Wallclock time (sec)", fontsize=fontsize)

#
# Name of data
#
dataname = "TODO"
if 'prog_h' in muletag:
	dataname = "surface height $h$"

#
# Norm
#
if 'linf' in muletag:
	norm = "$L_\infty$"
else:
	norm = "$L_{TODO}$"


plt.ylabel(norm+" error on "+dataname, fontsize=fontsize)



plt.legend(fontsize=15)

plt.savefig(output_filename, transparent=True, bbox_inches='tight', pad_inches=0)

plt.close()

