#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import sys
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

#
# First, use
#   ./postprocessing.py > postprocessing_output.txt
# to generate the .txt file
#

fig, ax = plt.subplots(figsize=(8,5))
#fig, ax = plt.subplots(figsize=(13,7))


ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

#mode = 'simtime'
mode = 'dt'


with open('./postprocessing_h_compare_ref_output_text.txt') as f:
	lines = f.readlines()


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ' and m != '':
            markers.append(m)

    except TypeError:
        pass



linestyles = ['-', '--', ':', '-.']



if len(sys.argv) > 1:
	output_filename = sys.argv[1]
else:
	output_filename = "./postprocessing_h_plot_err_vs_dt.pdf"

if len(sys.argv) > 2:
	plot_set = sys.argv[2:]
else:
	plot_set = []



def plot(x, y, marker, linestyle, label):

	# plot values and prev_name
	print(label)
	#print(values_err)
	#print(values_time)
	#print("")

	if len(x) == 0:
		return


	if len(plot_set) != 0:
		if prev_name not in plot_set:
			return

	ax.plot(x, y, marker=marker, linestyle=linestyle, label=label)



prev_name = ''
values_err = []
values_time = []
c = 2
for l in lines:
	if l[-1] == '\n':
		l = l[0:-1]

	d = l.split("\t")

	if d[0] == 'Running tests for new group:' or len(d) != 5:
		if len(values_time) == 0:
			continue

		plot(values_time, values_err, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)
#		for i, txt in enumerate(values_time):
#			ax.annotate("%.1f" % txt, (values_time[i]*1.03, values_err[i]*1.03))

		prev_name = d[0]
		values_err = []
		values_time = []
		c = c+1
		continue


	if d[0] == 'SIMNAME':
		continue

	prev_name = d[0]
	prev_name = prev_name.replace('script_g9.80616_h10000_f7.292e-05_a6371220_fsph0_u0_U0_tsm_', '')
	prev_name = prev_name.replace('h0.15_nrm1_hlf0_bf0_ext00_M0128_MPI_space01_', '')
	print(prev_name)

	prev_name = prev_name.replace('tsob1_', '')

#	prev_name = prev_name.replace('C0100_', '')
#	prev_name = prev_name.replace('C0200_', '')
#	prev_name = prev_name.replace('C0400_', '')
#	prev_name = prev_name.replace('C0800_', '')
#	prev_name = prev_name.replace('C1600_', '')
#	prev_name = prev_name.replace('C3200_', '')
#	prev_name = prev_name.replace('C129600_', '')

	prev_name = prev_name.replace('_time512', '')
	prev_name = prev_name.replace('l_rexi_tso0_REXITER', 'l_rexi')

	prev_name = prev_name.replace('REXITER_m00000256_h0.15_nrm1_hlf0_bf0_ext00_M0128', '')

#	prev_name = prev_name.replace('_mr10.0_mi30.0', '')
#	prev_name = prev_name.replace('_n0064_sx50.0_sy50.0', '')
#	prev_name = prev_name.replace('_n0064', '')
#	prev_name = prev_name.replace('_sx50.0_sy50.0', '')

	prev_name = re.sub(r"_mu.*", "", prev_name)
	prev_name = re.sub(r"0000", "", prev_name)

	values_err.append(float(d[1]))



	if mode == 'simtime':
		#
		# SIMTIME
		#
		values_time.append(float(d[4]))
		plt.xlabel("simulation time")

	elif mode == 'dt':
		#
		# DT
		#
		m = re.search('_C([0-9]*)', d[0])
		dt = float(m.group(1))
		values_time.append(dt)
		plt.xlabel("Timestep size")

	plt.ylabel("Error")


plot(values_time, values_err, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)

plt.legend()

plt.savefig(output_filename)
#plt.show()

