#! /usr/bin/python2

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import sys
import numpy
import csv


if len(sys.argv) < 2:
	print("Usage: "+sys.argv[0]+" [infile.csv] [outfile.whatever]")
	sys.exit(1)

infile = sys.argv[1]


if len(sys.argv) > 2:
	outfile = sys.argv[2]
else:
	outfile = ""


legend_min = -1
legend_max = 1

if len(sys.argv) > 3:
	legend_min = float(sys.argv[3])

if len(sys.argv) > 4:
	legend_max = float(sys.argv[4])

legend_norm = Normalize(legend_min, legend_max)
#print "Legend: "+str(legend_min)+", "+str(legend_max)

col_title = ''
row_title = ''
title = ''
subtitle = ''

print "Processing "+infile

col_labels = []
row_labels = []
Z = []
with open(infile, 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	c = 0
	for row in spamreader:
		some_data = [float(d) for d in row]

		# read row labels
#		row_labels.append(some_data[0])

		# append data
		Z.append(some_data[:])

Z = numpy.array(Z)



#
# PLOTTING
#
fig, ax = plt.subplots()

row_title = 'y'
col_title = 'x'
title = sys.argv[1]

print
print "Plotting..."
print " + title: "+title
print " + subtitle: "+subtitle
print " + row_title: "+row_title
print " + col_title: "+col_title
print " + array size: "+str(Z.shape)
print

heatmap = ax.pcolor(Z, cmap=cm.rainbow, norm=legend_norm)
#heatmap = ax.pcolor(Z, cmap=cm.rainbow)

#plt.title(title, y=1.08, fontweight='bold')
plt.title(title, fontweight='bold')
#plt.suptitle(subtitle, y=0.95)

#legend
cbar = plt.colorbar(heatmap)
#cbar.set_label('', rotation=270)

#fig.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.3)

ax.set_xlim(0, Z.shape[0])
ax.set_xlabel(col_title)

ax.set_ylim(0, Z.shape[1])
ax.set_ylabel(row_title)


if outfile != '':
	plt.savefig(outfile, format=outfile[-3:])
else:
	plt.show()

