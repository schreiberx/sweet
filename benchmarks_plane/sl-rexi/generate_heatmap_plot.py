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


legend_min = 1e-11
legend_max = 1e-0

if len(sys.argv) > 3:
	legend_min = float(sys.argv[3])

if len(sys.argv) > 4:
	legend_max = float(sys.argv[4])

legend_norm = LogNorm(legend_min, legend_max)

print "Legend: "+str(legend_min)+", "+str(legend_max)

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

	# read col headers
	while True:
		row = next(spamreader)
		if row[0][0:4] == '#TI ':
			title = row[0][4:]
			continue

		if row[0][0:4] == '#TS ':
			subtitle = row[0][4:]
			continue

		if row[0][0:4] == '#TX ':
			col_title = row[0][4:]
			continue

		if row[0][0:4] == '#TY ':
			row_title = row[0][4:]
			continue

		# skip first lines starting with #
		if row[0][0] != '#':
			break

	# axis descriptions in the format [Y]/[X]
	XY_desc = row[0]
	col_labels = [float(d) for d in row[1:]]

	# convert to int if integer
	col_labels = [int(c) if round(c) == c else c for c in col_labels]

	for row in spamreader:
		some_data = [float(d) for d in row]

		# read row labels
		if type(some_data[0]) == float:
			some_data[0] = '{0:.4G}'.format(some_data[0])
			print some_data[0]

		row_labels.append(some_data[0])

		# read data
		Z.append(some_data[1:])

Z = numpy.array(Z)



#
# PLOTTING
#
fig, ax = plt.subplots()

desc = XY_desc.split('\\')

if row_title == '':
	row_title = desc[0]
if col_title == '':
	col_title = desc[1]
if title == '':
	title = infile

print "Plotting..."
print " + title: "+title
print " + subtitle: "+subtitle
print " + row_title: "+row_title
print " + col_title: "+col_title
print " + array size: "+str(Z.shape)

heatmap = ax.pcolor(Z, cmap=cm.rainbow, norm=legend_norm)

plt.title(title, y=1.08, fontweight='bold')
plt.suptitle(subtitle, y=0.95)

#legend
cbar = plt.colorbar(heatmap)
cbar.set_label('RMS error in height', rotation=270)

fig.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.3)

# put the major ticks at the middle of each cell
ax.set_xticks(numpy.arange(len(col_labels))+0.5, minor=False)
ax.set_xticklabels(col_labels, minor=False, rotation=45)
ax.set_xlim(0, len(col_labels))
if len(desc) == 2:
	ax.set_xlabel(col_title)

ax.set_yticks(numpy.arange(len(row_labels))+0.5, minor=False)
ax.set_yticklabels(row_labels, minor=False)
ax.set_ylim(0, len(row_labels))
if len(desc) == 2:
	ax.set_ylabel(row_title)

# want a more natural, table-like display
ax.invert_yaxis()


if outfile != '':
	plt.savefig(outfile, format='pdf')
else:
	plt.show()

