#! /usr/bin/env python3

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

import numpy as np


datafiles = sys.argv[1:]



#################################################################################
#################################################################################
#################################################################################

fig, ax = plt.subplots(figsize=(8,5))


#ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
ax.set_ylim(4e-1, 300)


output_filename = "output_rexi_performance_breakdown.pdf"



bar_name = []
bar_preprocess = []
bar_broadcast = []
bar_reduce = []
bar_solve = []
bar_wallclock_time = []



if len(datafiles) == 0:
	print("Please specify datafiles")
	sys.exit(1)


c = 2
for i in datafiles:
	
	with open(i+'/output.out') as f:
		lines = f.readlines()

	wallclock_time = -1.0
	rexi_preprocess = -1.0
	rexi_broadcast = -1.0
	rexi_reduce = -1.0
	rexi_solve = -1.0

	for l in lines:
		tag = "REXI STOPWATCH preprocessing: "
		if l[0:len(tag)] == tag:
			rexi_preprocess = float(l[len(tag):])
			continue

		tag = "REXI STOPWATCH reduce: "
		if l[0:len(tag)] == tag:
			rexi_reduce = float(l[len(tag):])
			continue

		tag = "REXI STOPWATCH solve_rexi_terms: "
		if l[0:len(tag)] == tag:
			rexi_solve = float(l[len(tag):])
			continue

		tag = "REXI STOPWATCH broadcast: "
		if l[0:len(tag)] == tag:
			rexi_broadcast = float(l[len(tag):])
			continue

		tag = "Wallclock time (seconds): "
		if l[0:len(tag)] == tag:
			wallclock_time = float(l[len(tag):])
			continue


	if rexi_broadcast < 0:
		continue


	prev_name = i
	prev_name = prev_name.replace('script_ln2_g9.81_h10000_f7.2921e-05_a6371220_u0.0_U0_fsph0_tsm_', '')
	#prev_name = prev_name.replace('tso2_tsob2_', '')
	prev_name = re.sub(r"C[0-9]*_", "", prev_name)
	prev_name = prev_name.replace('REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128_', '')
	prev_name = prev_name.replace('tso2_tsob2_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_', '')
	prev_name = prev_name.replace('M0128_MPI_space01_', '')


	if prev_name[-1] == '_':
		prev_name = prev_name[:-1]

	print(prev_name)


	bar_name.append(prev_name)
	bar_preprocess.append(rexi_preprocess)
	bar_broadcast.append(rexi_broadcast)
	bar_solve.append(rexi_solve)
	bar_reduce.append(rexi_reduce)
	bar_wallclock_time.append(wallclock_time)

	print("")
	print(prev_name)
	print(wallclock_time)
	#print(rexi_preprocess)
	print(rexi_broadcast)
	print(rexi_solve)
	print(rexi_reduce)
	


if len(bar_name) == 0:
	print("No valid data found")
	sys.exit(1)



N = len(bar_name)
ind = range(N)
width = 0.2
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']


#rects0 = ax.bar([ind[i] - width*2 for i in range(N)], bar_preprocess, width, color=colors[0])
rects0 = ax.bar([ind[i] - width*2 for i in range(N)], bar_wallclock_time, width, color=colors[0])
rects1 = ax.bar([ind[i] - width*1 for i in range(N)], bar_broadcast, width, color=colors[1])
rects2 = ax.bar([ind[i] + width*0 for i in range(N)], bar_solve, width, color=colors[2])
rects3 = ax.bar([ind[i] + width*1 for i in range(N)], bar_reduce, width, color=colors[3])
ax.set_xticklabels(bar_name)


ax.set_ylabel('Wallclock time (seconds)')

#ax.legend((rects0[0], rects1[0], rects2[0], rects3[0]), ('preprocess', 'broadcast', 'solve', 'reduce'))
ax.legend((rects0[0], rects1[0], rects2[0], rects3[0]), ('wallclock time', 'broadcast', 'solve', 'reduce'))

ax.set_xticks([ind[i]-width/2 for i in range(N)])
ax.set_xticklabels(bar_name)


for tick in ax.get_xticklabels():
    tick.set_rotation(45)

plt.legend()
plt.tight_layout()

if output_filename != '':
	plt.savefig(output_filename)
else:
	plt.show()

