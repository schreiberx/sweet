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
	
	filename = i+'/output.out'
	if not os.path.isfile(filename):
		print("File '"+filename+"' doesn't exist")
		continue

	with open(filename) as f:
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
	print("prev_name: "+str(prev_name))
	print("wallclock_time: "+str(wallclock_time))
	print("rexi_preprocess: "+str(rexi_preprocess))
	print("rexi_broadcast: "+str(rexi_broadcast))
	print("rexi_solve: "+str(rexi_solve))
	print("rexi_reduce: "+str(rexi_reduce))
	

labels = ['Total wallclock', 'REXI misc', 'REXI broadcast', 'REXI solve', 'REXI reduce']
values = [bar_wallclock_time, bar_preprocess, bar_broadcast, bar_solve, bar_reduce]
M = len(values)

if len(bar_name) == 0:
	print("No valid data found")
	sys.exit(1)



N = len(bar_name)
ind = range(N)
width = 0.18
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

ax.set_xticklabels(bar_name)
ax.set_ylabel('Wallclock time (seconds)')

ax.set_xticks([ind[i]-width/2 for i in range(N)])
ax.set_xticklabels(bar_name)


for tick in ax.get_xticklabels():
    tick.set_rotation(45)



rects = []

for j in range(M):
	rects_ = ax.bar([i - width*(M/2-j) for i in range(N)], values[j], width, color=colors[j])
	rects.append(rects_)


for j in range(M):
	for rect, val in zip(rects[j], values[j]):
		height = rect.get_height()
		text =  "%.2f" % height
		ax.text(rect.get_x() + rect.get_width()/2, height+math.log(1.0+0.1*height), text, ha='center', va='bottom', size=8, rotation=90)


ax.legend([rect[0] for rect in rects], labels)


plt.legend()
plt.tight_layout()

if output_filename != '':
	plt.savefig(output_filename)
else:
	plt.show()

