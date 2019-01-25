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
	output_filename = sys.argv[1]
else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [output_filename.pdf] [jobdir1] [jobdir2] ... [jobdirN]")
	print("")
	sys.exit(1)


if len(sys.argv) > 2:
	# Load Jobs specified via program parameters
	jd = JobsData(job_dirs=sys.argv[2:])
else:
	# Load all Jobs
	jd = JobsData()


# Consolidate data...
jdc = JobsDataConsolidate(jd)

# ... which belongs to the same time integration method
jdc_groups = jdc.create_groups(['runtime.timestepping_method'])

# Exctract data suitable for plotting
jdc_groups_data = JobsData_GroupsPlottingScattered(
			jdc_groups,
			'runtime.timestep_size',
			'plane_data_diff_prog_h_pert.norm_linf',
		)

data = jdc_groups_data.get_data()


def label(d):
	val = d['runtime.timestepping_method'].replace('_', '\\_')+', $\Delta t = '+str(d['runtime.timestep_size'])+'$'
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
plt.rc('text', usetex=True)

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

	print(" + "+l)
	print(x)
	print(y)
	ax.plot(x, y, marker=markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=l)

	c = c + 1


if title != '':
	plt.title(title, fontsize=fontsize)

plt.xlabel("Timestep size $\Delta t$ (sec)", fontsize=fontsize)
plt.ylabel("$L_\infty$ error on surface height $h$", fontsize=fontsize)

plt.legend(fontsize=15)

plt.savefig(output_filename, transparent=True, bbox_inches='tight', pad_inches=0)

plt.close()

