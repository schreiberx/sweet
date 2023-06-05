#! /usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys

first = True

s = 2e-5
eta_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))
hs = 5
h_contour_levels = np.append(np.arange(900, 1000-hs, hs), np.arange(1000+hs, 1100, hs))

zoom_lat = True
zoom_lat = False

zoom_lat = 'eta' in sys.argv[1]

fontsize=8

figsize=(9, 3)

for filename in sys.argv[1:]:

    print(filename)
    #data = np.loadtxt(filename, skiprows=3)
    data = np.loadtxt(filename)

#    labelsx = data[0,1:]
#    labelsy = data[1:,0]

    data = data[1:,1:]

    if first:
    	cmin = np.amin(data)
    	cmax = np.amax(data)

    	if 'eta' in filename:
    		cmin = 1e-4
    		cmax = -1e-4
    		#cmin *= 1.2
    		#cmax *= 1.2

    plt.figure(figsize=figsize)


    plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto')

    plt.clim(cmin, cmax)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=fontsize) 

    plt.title(filename, fontsize=fontsize)


    if 'prog_eta' in filename:
    	plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
    elif 'prog_h' in filename:
    	plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
#    elif '_u' in filename:
#    	hs = 0.001
#    	h_contour_levels = np.append(np.arange(-0.1, 0-hs, hs), np.arange(hs, 0.1, hs))
#    	plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
    else:
    	if cmin != cmax:
    		pass
    		#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)

    ax = plt.gca()
    ax.xaxis.set_label_coords(0.5, -0.075)

    plt.xlabel("Longitude", fontsize=fontsize)

    plt.ylabel("Latitude", fontsize=fontsize)

    #plt.show()
    outfilename = filename.replace('.csv', '.png')
    print(outfilename)

    plt.savefig(outfilename, dpi=200)
    plt.close()

    first = False
