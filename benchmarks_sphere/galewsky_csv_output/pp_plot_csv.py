#! /usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys

first = True

s = 2e-5
vort_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))

zoom_lat = True
zoom_lat = False

zoom_lat = 'vort' in sys.argv[1]

fontsize=8

figsize=(9, 3)

for filename in sys.argv[1:]:

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


#		while labelsx[1] < 90:
#			tmplabelsx = labelsx[0]
#			labelsx[0:-1] = labelsx[1:]
#			labelsx[-1] = tmplabelsx
#
#			tmpdata = data[:,0]
#			data[:,0:-1] = data[:,1:]
#			data[:,-1] = tmpdata

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

		if 'vort' in filename:
			cmin *= 1.2
			cmax *= 1.2

	extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])

	plt.figure(figsize=figsize)

	plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto')

	plt.clim(cmin, cmax)
	cbar = plt.colorbar()
	cbar.ax.tick_params(labelsize=fontsize) 

	plt.title(filename, fontsize=fontsize)

	if 'vort' in filename:
		plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=vort_contour_levels, linewidths=0.5)
	else:
		if cmin != cmax:
			plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)
		


	ax = plt.gca()
	ax.xaxis.set_label_coords(0.5, -0.075)

	plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("Longitude", fontsize=fontsize)

	plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("Latitude", fontsize=fontsize)

	#plt.show()
	outfilename = filename.replace('.csv', '.png')
	print(outfilename)

	plt.savefig(outfilename, dpi=200)
	plt.close()

	first = False
