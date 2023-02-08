#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def get_data_phys(filename):
    if filename.endswith(".sweet"):
        from mule.postprocessing.SphereDataSpectral import SphereDataSpectral

        if 0:
            from mule.postprocessing.SphereDataOperators import SphereDataOperators

            sphere_data = SphereDataSpectral(filename, setup_physical=False)
            ops = SphereDataOperators(file_info=sphere_data.file_info)
            data_phys = ops.spec2phys(sphere_data.data_spectral)
            file_info = sphere_data.file_info

            file_info['lons'] = ops.lons
            file_info['lats'] = ops.lats


        else:
            sphere_data = SphereDataSpectral(filename, setup_physical=True)
            data_phys = sphere_data.data_physical
            file_info = sphere_data.file_info

        # Convert to degree
        file_info['lons'] = file_info['lons']/(2*np.pi)*360
        file_info['lats'] = file_info['lats']/(2*np.pi)*360


    elif filename.endswith(".csv"):
        data = np.loadtxt(filename, skiprows=3)

        data_phys = data[1:,1:]

        file_info = {}

        # Lat/lon given in degree
        file_info['lons'] = data[0,1:]
        file_info['lats'] = data[1:,0]

        data_phys = np.flip(data_phys, axis=0)
        file_info['lats'] = np.flip(file_info['lats'])


    else:
        raise Exception("Unknown file ending")

    return (data_phys, file_info)




print(sys.argv)
if len(sys.argv) <= 1:
    print("Usage:")
    print(" "+sys.argv[0]+" [file.csv|file.sweet]")
    sys.exit(1)


input_filepath = sys.argv[1]
print("Loading data from '"+input_filepath+"'")

input_filepath_noext = os.path.splitext(input_filepath)[0]
input_filename_noext = os.path.basename(input_filepath_noext)

(data_phys, file_info) = get_data_phys(input_filepath)


# Clear plot
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(13, 8))

# Locations of ticks
xtickslocs = np.arange(data_phys.shape[1]) + 0.5
ytickslocs = np.arange(data_phys.shape[0]) + 0.5

# Labels of ticks
xticklabels = file_info['lons']
yticklabels = file_info['lats']

xticklabels = np.array([round(_, 1) for _ in xticklabels])
yticklabels = np.array([round(_, 1) for _ in yticklabels])

assert len(xtickslocs) == len(xticklabels)
assert len(ytickslocs) == len(yticklabels)

if True:
    """
    Cleanup ticks so that there are only Nx ticks
    """
    Nx = 16
    N = len(xticklabels)
    step = max(1, N//Nx)
    r = np.arange(Nx, dtype=int)*step
    xtickslocs = xtickslocs[r]
    xticklabels = xticklabels[r]

    Ny = 8
    N = len(yticklabels)
    step = max(1, N//Ny)
    r = np.arange(Ny, dtype=int)*step
    ytickslocs = ytickslocs[r]
    yticklabels = yticklabels[r]


# Make pixel centered around integer coordinates
extent = [-0.5, data_phys.shape[1]-0.5, data_phys.shape[0]-0.5, -0.5]
imhandler = ax.imshow(data_phys, cmap="viridis", extent=extent)

if 'vrt' in input_filename_noext:
    e=2e-5
    ax.contour(data_phys, levels=np.arange(e, e*50, e), linestyles='solid', linewidths=0.2, colors='black')
    ax.contour(data_phys, levels=np.arange(-e*50, 0, e), linestyles='dashed', linewidths=0.2, colors='black')
else:
    e=2e-5
    ax.contour(data_phys, colors="black", origin='lower', linewidths=0.5)


# Fontsize
fontsize = 12

# Colorbar
cbar = fig.colorbar(imhandler, ax=ax)
cbar.ax.tick_params(labelsize=fontsize) 


# Axis labels
#ax = fig.gca()

ax.set_xticks(xtickslocs)
ax.set_xticklabels([round(_, 1) for _ in xticklabels], fontsize=fontsize)
ax.set_xlabel("Longitude", fontsize=fontsize)

ax.set_yticks(ytickslocs)
ax.set_yticklabels([round(_, 1) for _ in yticklabels], fontsize=fontsize)
ax.set_ylabel("Latitude", fontsize=fontsize)


ax.set_title(input_filename_noext)

output_filepath = input_filepath_noext+".png"

fig.tight_layout()

plt.show()

#print("Writing to "+str(output_filepath))
#fig.savefig(output_filepath, dpi=300)

