#! /usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys

first = True

s = 200
h_contour_levels_ = [np.arange(8000, 10000, s), np.arange(10000, 12000, s)]

s = 2000
phi_contour_levels_ = [np.arange(80000, 100000, s), np.arange(100000, 120000, s)]

s = 4.0 #*9.80616
ht0diff_contour_levels_ = [np.arange(-s*30, -s+s/2, s), np.arange(s, s*30, s)]
print(ht0diff_contour_levels_)

s = 2e-5
vort_contour_levels_ = [np.arange(-s*30, -s+s/2, s), np.arange(s, s*30, s)]
s = 4e-7
div_contour_levels_ = [np.arange(-s*30, -s+s/2, s), np.arange(s, s*30, s)]
print(div_contour_levels_)

print("h_contour_levels:")
print(h_contour_levels_)
print("")

print("phi_contour_levels:")
print(phi_contour_levels_)
print("")

print("ht0diff_contour_levels:")
print(ht0diff_contour_levels_)
print("")

print("vort_contour_levels:")
print(vort_contour_levels_)
print("")

print("div_contour_levels:")
print(vort_contour_levels_)
print("")


fontsize = 8
figsize = (9, 3)


for filename in sys.argv[1:]:

    print("")
    print("*"*80)
    print("* "+filename)
    print("*"*80)
    #data = np.loadtxt(filename, skiprows=3)
    data = np.loadtxt(filename)

    with open(filename, 'r') as f:
        line = f.readline()
        if line[0:3] == "#TI":
            # eliminate lat/lon coordinates
            data = data[1:,1:]

    # Detect SWEET Output file headers
    #TI
    #TX Longitude
    #TY Latitude


    data_ = np.zeros_like(data)
    print(data.shape)
    s = data.shape[1]//2
    data_[:,s:] = data[:,:s]
    data_[:,:s] = data[:,s:]
    #data = data_


    print("+ min(data): "+str(np.min(data)))
    print("+ max(data): "+str(np.max(data)))

    if 'prog_vort_' in filename:
        contour_levels = vort_contour_levels_

    elif 'prog_phi_' in filename:
        contour_levels = phi_contour_levels_

    elif 'prog_h_' in filename:
        contour_levels = h_contour_levels_

    elif 'prog_ht0diff_' in filename:
        contour_levels = ht0diff_contour_levels_

    elif 'prog_div_' in filename:
        contour_levels = div_contour_levels_

    else:
        contour_levels = None
        cmin = np.min(data)
        cmax = np.max(data)

    # min/max from contour
    if contour_levels != None:
        cmin = np.min(contour_levels[0])
        cmax = np.max(contour_levels[1])


    plt.figure(figsize=figsize)

    plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto')

    plt.clim(cmin, cmax)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=fontsize) 

    plt.title(filename, fontsize=fontsize)

    if contour_levels != None:
        plt.contour(data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=contour_levels[0], linewidths=0.5, linestyles='solid')
        plt.contour(data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=contour_levels[1], linewidths=0.5, linestyles='dashed')

    ax = plt.gca()
    ax.xaxis.set_label_coords(0.5, -0.075)

    plt.xlabel("Longitude", fontsize=fontsize)

    plt.ylabel("Latitude", fontsize=fontsize)

    outfilename = filename.replace('.csv', '.png')
    print(outfilename)

    plt.tight_layout()
    plt.savefig(outfilename, dpi=200)
    plt.close()

    first = False
