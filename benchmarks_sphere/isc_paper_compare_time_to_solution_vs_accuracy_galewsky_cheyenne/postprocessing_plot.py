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


# Mode: wallclocktime, dt
mode = 'dt'
if len(sys.argv) > 1:
	mode = sys.argv[1]

# input filename
input_filename = 'postprocessing_output_h.txt'
if len(sys.argv) > 2:
	input_filename = sys.argv[2]

# output filename
#output_filename = "./postprocessing_output_h_err_vs_"+mode+".pdf"
output_filename = ""
if len(sys.argv) > 3:
	output_filename = sys.argv[3]


#################################################################################
#################################################################################
#################################################################################

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



linestyles = ['-', '--', ':', '-.']


with open(input_filename) as f:
	lines = f.readlines()



def plot(x, y, marker, linestyle, label):

	# plot values and prev_name
	print(label)

	if len(x) == 0:
		return

	ax.plot(x, y, marker=marker, linestyle=linestyle, label=label)

	px = x[:]
	py = y[:]

	px = px[0::2]
	py = py[0::2]

	if len(px) % 2 == 0:
		px.append(x[-1])
		py.append(y[-1])

	for i, txt in enumerate(px):
		text = "%.1f/%.1f" % (py[i], px[i])

		if mode == 'dt':
			ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=8)
		elif mode == 'wallclocktime':
			ax.annotate(text, (px[i]*1.03, py[i]*1.03), fontsize=8)



prev_name = ''
values_y = []
values_x = []
c = 2
for l in lines:
	if l[-1] == '\n':
		l = l[0:-1]

	d = l.split("\t")

	if d[0] == 'Running tests for new group:':
		plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)

		prev_name = d[0]
		values_y = []
		values_x = []
		c = c+1
		continue

	if len(d) != 5:
		continue

	if d[0] == 'SIMNAME':
		continue

	prev_name = d[0]
	prev_name = prev_name.replace('script_ln2_g9.81_h10000_f7.2921e-05_a6371220_u0.0_U0_fsph0_tsm_', '')
	#prev_name = prev_name.replace('tso2_tsob2_', '')
	prev_name = re.sub(r"C[0-9]*_", "", prev_name)
	prev_name = prev_name.replace('REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128', '')
	prev_name = prev_name.replace('M0128_MPI_space01_time128', '')
	prev_name = prev_name.replace('M0128_MPI_space01_time001', '')

	if prev_name[-1] == '_':
		prev_name = prev_name[:-1]

	# skip invalid nan's
	if d[1] == 'nan':
		continue

	values_y.append(float(d[1]))



	if mode == 'wallclocktime':
		#
		# SIMTIME
		#
		values_x.append(float(d[4]))
		plt.xlabel("Wallclock time")

	elif mode == 'dt':
		#
		# DT
		#
		m = re.search('_C([0-9]*)', d[0])
		dt = float(m.group(1))
		values_x.append(dt)
		plt.xlabel("Timestep size")

	plt.ylabel("Error")


plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)



plt.legend()

if output_filename != '':
	plt.savefig(output_filename)
else:
	plt.show()

