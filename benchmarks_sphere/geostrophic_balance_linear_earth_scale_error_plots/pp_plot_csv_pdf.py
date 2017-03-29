#! /usr/bin/python3

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys

first = True

zoom_lat = True
zoom_lat = False

zoom_lat = 'eta' in sys.argv[1]

fontsize=8

figsize=(9, 3)

files = sys.argv[1:]

refdataavailable = False

if files[0] == 'reference':
	reffilename = files[1]
	files = files[2:]

	print("Loading reference solution from '"+reffilename+"'")
	refdata = np.loadtxt(reffilename, skiprows=3)
	refdata = refdata[1:,1:]

	refdataavailable = True


#for filename in sys.argv[1:]:
for filename in files:

	print(filename)
	data = np.loadtxt(filename, skiprows=3)

	labelsx = data[0,1:]
	labelsy = data[1:,0]

	data = data[1:,1:]

	if zoom_lat:
		while labelsy[1] < 10:
			labelsy = labelsy[1:]
			data = data[1:]

		while labelsy[-2] > 80:
			labelsy = labelsy[0:-2]
			data = data[0:-2]

	if first:
		lon_min = labelsx[0]
		lon_max = labelsx[-1]

		lat_min = labelsy[0]
		lat_max = labelsy[-1]

		new_labelsx = np.linspace(lon_min, lon_max, 7)
		new_labelsy = np.linspace(lat_min, lat_max, 7)

	labelsx = np.interp(new_labelsx, labelsx, labelsx)
	labelsy = np.interp(new_labelsy, labelsy, labelsy)


	if first:
		cmin = np.amin(data)
		cmax = np.amax(data)

	extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])

	plt.figure(figsize=figsize)
	plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('rainbow'))
	plt.clim(cmin, cmax)

	cbar = plt.colorbar()
	cbar.ax.tick_params(labelsize=fontsize) 

	plt.title(filename, fontsize=fontsize)

	ax = plt.gca()
	ax.xaxis.set_label_coords(0.5, -0.075)

	plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("Longitude", fontsize=fontsize)

	plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("Latitude", fontsize=fontsize)

	#plt.show()
	outfilename = filename.replace('.csv', '.pdf')
	print(outfilename)

	#plt.savefig(outfilename, dpi=200)
	plt.savefig(outfilename)
	plt.close()

	first = False
