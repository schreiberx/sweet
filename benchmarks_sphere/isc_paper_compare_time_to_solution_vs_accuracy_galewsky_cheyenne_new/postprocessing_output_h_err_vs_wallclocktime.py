#! /usr/bin/env python3

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

#
# First, use
#   ./postprocessing.py > postprocessing_output.txt
# to generate the .txt file
#

fig, ax = plt.subplots(figsize=(10,7))


ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

mode = 'wallclocktime'
#mode = 'dt'


with open('postprocessing_output_h.txt') as f:
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
	output_filename = "./postprocessing_output_h_err_vs_"+mode+".pdf"

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

	if d[0] == 'Running tests for new group:':
		plot(values_time, values_err, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)
		for i, txt in enumerate(values_time):
			ax.annotate("%.1f" % txt, (values_time[i]*1.03, values_err[i]*1.03))

		prev_name = d[0]
		values_err = []
		values_time = []
		c = c+1
		continue

	if len(d) != 5:
		continue

	if d[0] == 'SIMNAME':
		continue

	prev_name = d[0]
	prev_name = prev_name.replace('script_ln2_g9.81_h10000_f7.2921e-05_a6371220_u0.0_U0_fsph0_tsm_', '')
	prev_name = prev_name.replace('tso2_tsob2_', '')
	prev_name = prev_name.replace('REXICI_n00000128_sx50.0_sy50.0_mu0.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf1_bf1e-16_', '')
	prev_name = prev_name.replace('REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd1.3000E+02_gfe2.0000E+00_nrm0_hlf0_bf1e-16_', '')
	prev_name = prev_name.replace('ext00_M0128_MPI_space01_time001', '')
	prev_name = prev_name.replace('ext00_M0128_MPI_space01_time128', '')

	prev_name = prev_name.replace('C000060_', '')
	prev_name = prev_name.replace('C000120_', '')
	prev_name = prev_name.replace('C000180_', '')
	prev_name = prev_name.replace('C000360_', '')
	prev_name = prev_name.replace('C000720_', '')

	if prev_name[-1]:
		prev_name = prev_name[0:-1]

	#prev_name = re.sub(r"_mu.*", "", prev_name)
	#prev_name = re.sub(r"0000", "", prev_name)

	values_err.append(float(d[1]))



	if mode == 'wallclocktime':
		#
		# SIMTIME
		#
		values_time.append(float(d[4]))
		plt.xlabel("Wallclock time")

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

