#! /usr/bin/env python3
#
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mule.postprocessing.JobData import *

import numpy as np
import sys
import os

import scipy
from scipy.interpolate import RectBivariateSpline

import matplotlib.pylab as pylab
params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18}
pylab.rcParams.update(params)

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": [r'\usepackage{amsmath}']
})


Omega = 7.292e-5
g = 9.80616
H = 1e4
######def computePV(lon, lat, vrt, phi_pert):
######    return vrt * phi_pert / g

def computePV(lon, lat, vrt, phi):
    return (vrt + 2. * Omega * np.transpose(np.vstack([np.sin(lat * np.pi / 180.)] * vrt.shape[1]))) / (phi / g)

## interpolate from fine to coarse grid
def interpolateSolution(data_ref, data_cmp):

    (ny_ref, nx_ref) = data_ref.shape
    (ny_cmp, nx_cmp) = data_cmp.shape

    multiplier_j = (ny_ref)/(ny_cmp)
    multiplier_i = (nx_ref)/(nx_cmp)

    #Comparison via interpolation
    #print("Interpolation")
    # A-grid REFERENCE (file1) - sweet outputs only A grids physical space
    dx_ref=1.0/(nx_ref)
    dy_ref=1.0/(ny_ref)

    x_ref = np.arange(0, 1, dx_ref)
    x_ref = np.linspace(0, 1, nx_ref, endpoint=False)

    y_ref = np.arange(0, 1, dy_ref)
    y_ref = np.linspace(0, 1, ny_ref, endpoint=False)

    x_ref += dx_ref/2
    y_ref += dy_ref/2
    X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

    #Creat cubic interpolation of reference file
    print(y_ref.shape, x_ref.shape, data_ref.shape)
    interp_spline = RectBivariateSpline(y_ref, x_ref, data_ref)

    #A-grid cmp file (file2)
    dx_cmp=1.0/nx_cmp
    dy_cmp=1.0/ny_cmp

    x_cmp = np.arange(0, 1, dx_cmp)
    x_cmp = np.linspace(0, 1, nx_cmp, endpoint=False)

    y_cmp = np.arange(0, 1, dy_cmp)
    y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)

    x_cmp += dx_cmp/2
    y_cmp += dy_cmp/2
    X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

    #Get reduced reference resolution
    data_ref_low = interp_spline(y_cmp, x_cmp)

    return data_ref_low


def computeDiff(data_ref, data_cmp, relative = True):

    (ny_ref, nx_ref) = data_ref.shape
    (ny_cmp, nx_cmp) = data_cmp.shape

    multiplier_j = (ny_ref)/(ny_cmp)
    multiplier_i = (nx_ref)/(nx_cmp)

    print("MULTIPLIERS", multiplier_i, multiplier_j, nx_ref, ny_ref, nx_cmp, ny_cmp)

    if multiplier_i == 1 and multiplier_j == 1: #Grids are the same

        data = data_cmp - data_ref

        if relative:
            data = data / data_ref

    else:
        data_ref_low = interpolateSolution(data_ref, data_cmp)

        data = data_cmp - data_ref_low

        if relative:
            data = data / data_ref_low

    return data;


## command line parameters
geometry = sys.argv[2];                     ## plane / sphere
plot_type = sys.argv[3];                    ## map, contour lines or both
extent_type = sys.argv[4];                  ## full domain or zoom
compare_solution = int(sys.argv[5]);        ## 0: plot solution; 1: plot diff (sol - ref); 2: plot diff unstable jet paper
filename_ref = sys.argv[6];                 ## path containing ref solution
roll = int(sys.argv[7])                     ## 1: sum 180^o in longitude
ref_type = sys.argv[8]                      ## fine or ref
lim_ref = int(sys.argv[9])                  ## 0: limits of map legend defined from solution; 1: limits of map legend defined from ref solution
input_lim = None
if len(sys.argv) > 10:
    input_lim = float(sys.argv[10]);        ## limits of map legend: (-input_lim, input_lim)


assert geometry in ["plane", "sphere"]
assert plot_type in ["map", "contour", "both"]
assert extent_type in ["full", "zoom_gaussian_bump", "zoom_unstable_jet", "zoom_unstable_jet_2"]
assert compare_solution <= 2
assert ref_type in ["fine", "ref"]

#Figure definitions
fontsize=26
figsize=(9, 7)
if extent_type == "zoom_unstable_jet":
    figsize=(10, 3)

## get ref solution
jd = JobData(os.path.split(filename_ref)[0])
jobdata = jd.get_flattened_data()
try:
    data_ref = np.loadtxt(filename_ref);
    ref_ok = True
except:
    ref_ok = False

## get all solutions
filenames = [sys.argv[1]]
for filename in filenames:

    jd = JobData(os.path.split(filename)[0])
    jobdata = jd.get_flattened_data()

    filename_orig = filename
    if 'potvrt' in filename:
        filename = filename.replace('potvrt', 'vrt')

    #Load data
    print(filename)
    data = np.loadtxt(filename)

    if 'potvrt' in filename_orig:
        filename = filename_orig

    if not ref_ok:
        data_ref = data.copy()

    full_data = data;
    full_data_ref = data_ref;

    if geometry == "sphere":
        ## remove first line and colum (lat/lon values)
        data = data[1:, 1:]
        data_ref = data_ref[1:, 1:]

    data_orig = data.copy();
    data_ref_low = interpolateSolution(data_ref, data);
    ## compute difference between solution and ref
    ## only if plotting map!
    if compare_solution and plot_type == "map":
        relative = not ("vrt" in filename)  ## do not compute relative error for vorticity (div by zero!)
        data = computeDiff(data_ref, data, relative = relative);

    #Set physical grid for axis
    if geometry == "plane":
        x_min = 0
        x_max = 4.00316e+07/1000/1000
        y_min = 0
        y_max = 4.00316e+07/1000/1000
        n = data.shape[0]
        x = np.linspace(x_min, x_max, n)
        y = np.linspace(y_min, y_max, n)
    else:
        x = full_data[0, 1:];
        y = full_data[1:, 0];
        x_ref = full_data_ref[0, 1:];
        y_ref = full_data_ref[1:, 0];

    if ('prog_h' in filename or 'prog_phi' in filename) and not compare_solution: ## define in km
        data = data / 1000;
        full_data = full_data / 1000;
        data_ref = data_ref / 1000;
        full_data_ref = full_data_ref / 1000;
        data_ref_low = data_ref_low / 1000;


    if 'potvrt' in filename:
        #####data_phi_pert = np.loadtxt(filename.replace("potvrt", "phi_pert"));
        #####data_phi_pert = data_phi_pert[1:, 1:];
        #####data = computePV(x, y, data, data_phi_pert);
        data_phi = np.loadtxt(filename.replace("potvrt", "phi"));
        data_phi = data_phi[1:, 1:];
        data = computePV(x, y, data, data_phi);
        filename = filename_orig

    ## load solution from paper of scott et al.
    ##print(float(filename.split("_t")[-1].split(".csv")[0]) / 24)
    if compare_solution == 1:
        if plot_type == "map" or plot_type == "both" or extent_type == "full" or ("prog_potvrt" not in filename):
            sys.exit()
        if "iter" in filename:
            time = float(filename.split("_t")[-1].split("_iter")[0]) / 24
        else:
            time = float(filename.split("_t")[-1].split(".csv")[0]) / 24
        print("TIME", time)
        if time == 5:
            x_ref, y_ref, data_ref = np.loadtxt('/home/jcaldass/Development/IME/PostDocIMEUSP/Biblio/SWE/UnstableJet/qj2667-sup-0001-files1.dat',usecols=(0,1,2),unpack=True);
        elif time == 6:
            x_ref, y_ref, data_ref = np.loadtxt('/home/jcaldass/Development/IME/PostDocIMEUSP/Biblio/SWE/UnstableJet/qj2667-sup-0002-files2.dat',usecols=(0,1,2),unpack=True);
        else:
            sys.exit()
        data_ref = data_ref * 1e-8
        print(data_ref, np.min(data_ref), np.max(data_ref))


    #Get max/min
    cmin = np.amin(data)
    cmax = np.amax(data)
    print(filename + " cmin cmax diff " + str(cmin) + " " + str(cmax));


    if geometry == "plane":
        #Labels
        labelsx = np.linspace(x_min, x_max, 7)
        labelsy = np.linspace(y_min, y_max, 7)

        #Contour levels for fields
        extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])

    else:
                ## TODO: check this!
        len_x = x.size;
        len_x_ref = x_ref.size;
        if roll:
            data = np.roll(data, len_x // 2, axis = 1);
            if compare_solution != 1:
                data_ref = np.roll(data_ref, len_x_ref // 2, axis = 1);
                data_ref_low = np.roll(data_ref_low, len_x // 2, axis = 1);


        aspect = 'auto';
        if extent_type == "full":
            if not roll:
                extent = (0, 360., -90., 90.)
            else:
                extent = (-180., 180., -90., 90.)

        else:
            if extent_type == "zoom_gaussian_bump":
                lon_min = -50. + 180.
                lon_max = 50. + 180.
                lat_min = -75.
                lat_max = -15.
            elif extent_type == "zoom_unstable_jet":
                lon_min = 0.
                lon_max = 360.
                lat_min = 0.
                lat_max = 90.
                aspect = 1.75;
            elif extent_type == "zoom_unstable_jet_2":
                lon_min = 240.
                lon_max = 300.
                lat_min = 15.
                lat_max = 75.


            extent = (lon_min, lon_max, lat_min, lat_max)

            ## select data corresponding to extent region
            idx_x = np.argwhere( (x >= lon_min) & (x <= lon_max) ).reshape(-1)
            idx_y = np.argwhere( (y >= lat_min) & (y <= lat_max) ).reshape(-1)
            idx_x_ref = np.argwhere( (x_ref >= lon_min) & (x_ref <= lon_max) ).reshape(-1)
            idx_y_ref = np.argwhere( (y_ref >= lat_min) & (y_ref <= lat_max) ).reshape(-1)
            x = x[idx_x]
            y = y[idx_y]
            x_ref = x_ref[idx_x_ref]
            y_ref = y_ref[idx_y_ref]
            data = data[np.ix_(idx_y, idx_x)]
            data_orig = data_orig[np.ix_(idx_y, idx_x)]
            if compare_solution != 1:
                data_ref = data_ref[np.ix_(idx_y_ref, idx_x_ref)]
                data_ref_low = data_ref_low[np.ix_(idx_y, idx_x)]

    #Start plotting
    plt.figure(figsize=figsize)

    #Colorbar
    origin = "lower";

    ## without path
    filename_file = os.path.basename(filename);

    print("FILENAME", filename)
    if plot_type == "map" or plot_type == "both":

        ## set colormap
        if compare_solution or 'vrt' in filename_file:
            cmap = 'seismic' ## blue: -; red: +
        else:
            cmap = 'jet'

        ## plot map
        if plot_type == "map":
            data_to_plot = data
        elif plot_type == "both":
            data_to_plot = data_ref_low
        plt.imshow(data_to_plot, interpolation='nearest', extent=extent, origin=origin, aspect=aspect, cmap=plt.get_cmap(cmap))

        ## set colorbar limits
        plt.clim(cmin, cmax)
        cref=max(abs(cmin),abs(cmax))
        if not input_lim is None:  ## user defined colorbar limits
            cref = input_lim
            cmax = cref
        if compare_solution or ('vrt' in filename_file and not lim_ref):   ## symmetric colorbar
            plt.clim(-cref, +cref)
        elif lim_ref: ## limits obtained from reference solution
            data_min_ref = np.min(data_ref.reshape(-1));
            data_max_ref = np.max(data_ref.reshape(-1));
            if 'vrt' in filename_file:  ## symmetric colorbar
                cref = max(abs(data_min_ref), abs(data_max_ref))
                plt.clim(-cref, cref)
            else:
                plt.clim(data_min_ref, data_max_ref)
        else:
            plt.clim(cmin, cmax)

        if aspect == 'auto':
            cbar = plt.colorbar()
        else:
            cbar = plt.colorbar(fraction = 0.017, pad = 0.04)


    if plot_type == "contour" or plot_type == "both":

        ## set contour levels
        if not compare_solution:
            data_min = np.min(data.reshape(-1));
            data_max = np.max(data.reshape(-1));
        else:
            data_min = np.min(data.reshape(-1));
            data_max = np.max(data.reshape(-1));
            data_min_ref = np.min(data_ref.reshape(-1));
            data_max_ref = np.max(data_ref.reshape(-1));
            data_min = min(data_min, data_min_ref);
            data_max = max(data_max, data_max_ref);
        dlevels = (data_max - data_min) / 10.;

        ###if compare_solution == 1:
        ###    dlevels = 0.2 * Omega / H

        levels = np.arange(data_min, data_max + dlevels / 2, dlevels);

        ## contour lines: solution (black) and reference (red)
        plt.contour(x, y, data, colors="black", origin='lower', extent=extent, vmin=data_min, vmax=data_max, linewidths=0.5, levels = levels)
        if compare_solution == 1:
            print(levels)
            plt.tricontour(x_ref, y_ref, data_ref, colors="red", origin='lower', extent=extent, vmin=data_min, vmax=data_max, linewidths=0.25, levels = levels)
        elif compare_solution == 2:
            plt.contour(x, y, data_ref_low, colors="red", origin='lower', extent=extent, vmin=data_min, vmax=data_max, linewidths=0.25, levels = levels)


    #Set tittle
    title=""
    if 'diag_vort' in filename_file:
        ##title += "Vorticity (1/t)"
        if (plot_type == "map" or plot_type == "both") and not compare_solution:
            cbar.set_label('1/s', rotation=270, size=fontsize, labelpad = 30)
    if 'prog_h' in filename_file:
        ##title += "Depth (km) "
        if (plot_type == "map" or plot_type == "both") and not compare_solution:
            cbar.set_label('km', rotation=270, size=fontsize, labelpad = 30)
    if 'prog_u' in filename_file:
        ##title += "U-Velocity (m/s) "
        if (plot_type == "map" or plot_type == "both") and not compare_solution:
            cbar.set_label('m/s', rotation=270, size=fontsize, labelpad = 30)
    if 'prog_v' in filename_file:
        ##title += "V-Velocity (m/s) "
        if (plot_type == "map" or plot_type == "both") and not compare_solution:
            cbar.set_label('m/s', rotation=270, size=fontsize, labelpad = 30)
    if 'prog_phi' in filename_file:
        ##title += "Depth (km) "
        if (plot_type == "map" or plot_type == "both") and not compare_solution:
            cbar.set_label(r'$\mathrm{m}^2/\mathrm{s}^2$', rotation=270, size=fontsize, labelpad = 30)

    if plot_type == "map" or plot_type == "both":
        cbar.ax.tick_params(labelsize=fontsize) 

    #Axis
    ax = plt.gca()
    if not extent_type == "zoom_unstable_jet":
        ax.xaxis.set_label_coords(0.5, -0.075)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    #plt.xticks(labelsx, fontsize=fontsize)
    if geometry == "plane":
        plt.xlabel("x (1000 km)", fontsize=fontsize)
    elif geometry == "sphere":
        plt.xlabel("Longitude (degrees)", fontsize=fontsize)
        if extent_type == "zoom_unstable_jet":
            plt.xlabel("Longitude (degrees)", fontsize=fontsize * .9)

    #plt.yticks(labelsy, fontsize=fontsize)
    if geometry == "plane":
        plt.ylabel("y (1000 km)", fontsize=fontsize)
    elif geometry == "sphere":
        plt.ylabel("Latitude (degrees)", fontsize=fontsize)
        if extent_type == "zoom_unstable_jet":
            plt.ylabel("Latitude (degrees)", fontsize=fontsize * .9)

    # Save file as eps
    if not compare_solution:
        outfilename = filename
    if compare_solution:
        outfilename = filename.replace('.csv', '_diff.pdf')
    if compare_solution:
        outfilename = outfilename[:-9] + "_" + plot_type + outfilename[-9:];
        outfilename = outfilename[:-9] + "_" + extent_type + outfilename[-9:];
        outfilename = outfilename[:-9] + "_comparison_" + ref_type + outfilename[-9:];
    else:
        outfilename = outfilename[:-4] + "_" + plot_type + outfilename[-4:];
        outfilename = outfilename[:-4] + "_" + extent_type + outfilename[-4:];

    outfilename = outfilename.replace('.csv', '.pdf')
    outfilename = outfilename.replace('output', 'plot')
    print(outfilename)

    plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)

    plt.close()
