#! /usr/bin/env python3

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

#
# First, use
# postprocess_generrors.sh > errors.txt
# to generate the .txt file
#


# Mode: wallclocktime, dt
mode = 'dt'
if len(sys.argv) > 1:
	mode = sys.argv[1]

# input filename
input_filename = 'errors.txt'
if len(sys.argv) > 2:
	input_filename = sys.argv[2]

# output filename
output_filename = "./postprocessing_output_h_err_vs_"+mode+".pdf"
if len(sys.argv) > 3:
	output_filename = sys.argv[3]

print(output_filename)
#################################################################################
#################################################################################
#################################################################################


#Plot style and definitions
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
	x, y = (list(t) for t in zip(*sorted(zip(x,y))))
	
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

		if mode == 'dt':
			#ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=8)
			ax.annotate(px[i], (px[i]*1.03, py[i]*0.92), fontsize=8)
		elif mode == 'wallclocktime':
			ax.annotate(text, (px[i]*1.03, py[i]*1.03), fontsize=8)



prev_name = ''
values_y = []
values_x = []
c = 2
prev_variable = ''
prev_method = ''

head = lines[0]
#print(head)
if head[-1] == '\n':
	head = head[0:-1]
head = head.split(" ")
#print(head)

imethod1 = head.index("Method1Paper")
imethod2 = head.index("Method2Paper")
il1er = head.index("L1Error")
il2er = head.index("RMSError")
imax = head.index("MaxError")
idt1 = head[:imethod2].index("dt")
idt2 = imethod2+head[imethod2:].index("dt")
itime1 = head[:imethod2].index("Time")
itime2 = imethod2+head[imethod2:].index("Time")

#Choose what error to plot
#ierr = il2er
ierr = imax
idt = idt2
imethod = imethod2
prev_name=''

for l in lines[1:]:
	if l[-1] == '\n': #get rid of \n
		l = l[0:-1]

	d = l.split(" ") #line data
	
	name = d[imethod]
	if prev_name == '':
		prev_name = name
		
	if  name != prev_name:
		
		plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], prev_name)

		prev_name = name
		values_y = []
		values_x = []
		c = c+1
		continue

	# skip invalid nan's
	if d[ierr] == 'nan':
		continue

	values_y.append(float(d[ierr]))
	values_x.append(float(d[idt]))
	plt.xlabel("Timestep size (sec)")
	plt.ylabel(head[imax])
	
plot(values_x, values_y, markers[c % len(markers)], linestyles[c % len(linestyles)], name)


#if mode == 'dt':
	#ax.xaxis.set_ticks([2**i for i in range(0, 10)])
	#ax.yaxis.set_ticks([2**i for i in range(-5, 5)])

plt.legend()

#if output_filename != '':
print(output_filename)
plt.savefig(output_filename)
#else:
#	plt.show()

