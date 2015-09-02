#! /usr/bin/python2

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import sys
import numpy
import csv


if len(sys.argv) < 2:
	print("Usage: "+sys.argv[0]+" [infile.csv] [outfile.pdf]")
	sys.exit(1)

infile = sys.argv[1]

if len(sys.argv) > 2:
	outfile = sys.argv[2]
else:
	outfile = ""


column_labels = []
row_labels = []
Z = []
with open(infile, 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	c = 0

	# read column headers
	row = next(spamreader)

	# axis descriptions in the format [Y]/[X]
	XY_desc = row[0]
	column_labels = [float(d) for d in row[1:]]

	# convert to int if integer
	column_labels = [int(c) if round(c) == c else c for c in column_labels]

	for row in spamreader:
		float_data = [float(d) for d in row]
		# read row headers
		row_labels.append(float_data[0])

		# read data
		Z.append(float_data[1:])

Z = numpy.array(Z)

print "Array size: "+str(Z.shape)

#
# PLOTTING
#
fig, ax = plt.subplots()


desc = XY_desc.split('\\')

heatmap = ax.pcolor(Z, cmap=cm.rainbow, norm=LogNorm(10e-12, vmax=1))

plt.title(infile)

#legend
cbar = plt.colorbar(heatmap)
cbar.set_label('RMS error in height', rotation=270)

fig.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.3)

# put the major ticks at the middle of each cell
ax.set_xticks(numpy.arange(len(column_labels))+0.5, minor=False)
ax.set_xticklabels(column_labels, minor=False, rotation=45)
ax.set_xlim(0, len(column_labels))
if len(desc) == 2:
	ax.set_xlabel(desc[1])

ax.set_yticks(numpy.arange(len(row_labels))+0.5, minor=False)
ax.set_yticklabels(row_labels, minor=False)
ax.set_ylim(0, len(row_labels))
if len(desc) == 2:
	ax.set_ylabel(desc[0])

# want a more natural, table-like display
ax.invert_yaxis()


if outfile != '':
	plt.savefig(outfile, format='pdf')
else:
	plt.show()

