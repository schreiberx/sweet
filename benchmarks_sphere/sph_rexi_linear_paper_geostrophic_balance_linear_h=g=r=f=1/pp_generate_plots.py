#! /usr/bin/env python2

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import sys
import os

plt.figure(figsize=(9, 5))
fontsize = 12

# 2: height
# 3: velocity u
# 4: velocity v
column=2

files = sys.argv[1:]


output_filename = 'output.pdf'
if files[0] == 'output':
	output_filename = files[1]
	files = files[2:]

num_plots = len(files)

# set color map
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

# markers
markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass

# Title
plt.title("Geostrphic balance benchmark")

#plt.xticks(labelsx, fontsize=fontsize)
plt.xlabel("Simulation time", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("Lmax error in height", fontsize=fontsize)

legend_labels=[]
for i in range(0, len(files)):
	f = files[i]

	data = np.loadtxt(f, skiprows=0, usecols=(1,2), delimiter='\t')
	data = data.transpose()

	marker = markers[i % len(markers)]
	plt.plot(data[0], data[1], linewidth=1.0, linestyle='-', marker=marker)

	f = f.replace('script_g1_h1_f1_a1_u0_U0_fsph0_tsm_l_rexi_tso1_tsob1_C0000.1_REXITER_', '')
	f = f.replace('_h0.15', '')
	f = f.replace('_bf0_ext02_M0064', '')

	#f = f.replace('script_modes128_bench10_nonlin0_g1_h1_f1_a1_u0_pdeid1_', '')
	#f = f.replace('_rexiextmodes02_rexipar1_C00000.01_t000000.1_o00000.01_robert1', '')

	legend_labels.append(f)



plt.legend(legend_labels, ncol=1, loc='upper center'
#           bbox_to_anchor=[0.5, 1.1], 
#           columnspacing=1.0, labelspacing=0.0,
#           handletextpad=0.0, handlelength=1.5,
#           fancybox=True, shadow=True
	)

print(output_filename)

if '.png' in output_filename:
	plt.savefig(output_filename, dpi=200)
else:
	plt.savefig(output_filename)
plt.close()

