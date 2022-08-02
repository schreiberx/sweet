#! /usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import sys
import mule.plot_config as pc


import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams.update({'font.size': 12})


matplotlib.rcParams["figure.facecolor"] = "ffffff"
matplotlib.rcParams["axes.facecolor"] = "ffffff"



"""
Load first file to get grid layout
"""
data_all = np.genfromtxt(sys.argv[1])
print(data_all.shape)

xs = data_all[0,1:]
ys = data_all[1:,0]
data = data_all[1:,1:]


X, Y = np.meshgrid(xs, ys)

for file in sys.argv[1:]:

    # Plot
    from matplotlib import pyplot, ticker, cm, axis

    # Plot
    fig, ax = pc.setup(square_figsize=True)
    ax.set_aspect(1.0)

    ps = pc.PlotStyles()

    plt.axhline(0, color='black', linewidth=0.5, linestyle="dashed")
    plt.axvline(0, color='black', linewidth=0.5, linestyle="dashed")


    c = 0
    colors = ps.colors

    legend_lines = []
    legend_labels = []


    """
    Load first file to get grid layout
    """
    data_all = np.genfromtxt(file)

    xs = data_all[0,1:]
    ys = data_all[1:,0]
    data = data_all[1:,1:]

    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]

    extent = (xs[0]-0.5*dx, xs[-1]+0.5*dx, ys[0]-0.5*dy, ys[-1]+0.5*dy)

    X, Y = np.meshgrid(xs, ys)

    Z = data

    aspect = xs[-1]/ys[-1]
    #aspect = 1
    fmt = '%1.1e'
    if "stability" in file:
        contour_range = [1e-8]
        contour_range = [0.5]

        vmin = -1
        vmax = +1

        Z[Z<vmin] = np.inf
        Z[Z>vmax] = np.inf

        ax.set_title("Stability plots")

        cs = plt.contour(X, Y, Z, contour_range, colors=["white" for _ in contour_range], linewidths=[0.5 for _ in contour_range], extent=extent)
        pyplot.clabel(cs, cs.levels, fontsize=4, fmt=fmt)
        ax.clabel(cs, cs.levels, fontsize=4)

        im = plt.imshow(np.flip(Z, 0), extent=extent, aspect=aspect)

        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


    elif "errors" in file:
        contour_range = [-1e-5, -1e-10, 1e-10, 1e-5]
        contour_range = [1e-5]

        vmin = 0
        vmax = +5

        Z[Z<vmin] = np.inf
        Z[Z>vmax] = np.inf

        ax.set_title("Error plots")

        cs = plt.contour(X, Y, Z, contour_range, colors=["white" for _ in contour_range], linewidths=[0.5 for _ in contour_range], extent=extent)
        pyplot.clabel(cs, cs.levels, fontsize=4, fmt=fmt)
        ax.clabel(cs, cs.levels, fontsize=4)

        if 1:
            from matplotlib.colors import LogNorm
            im = plt.imshow(np.flip(Z, 0), extent=extent, aspect=aspect, norm=LogNorm(vmin=1e-12, vmax=vmax))
        else:
            im = plt.imshow(np.flip(Z, 0), extent=extent, aspect=aspect)

        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


    #l = cs.legend_elements()
    #legend_lines += l[0]
    #legend_labels += [file]

    c = c+1

    ax.set_xlabel("Timestep size $\Delta t$")
    ax.set_ylabel("Imaginary value of $\lambda_2$")

    #plt.legend(legend_lines, legend_labels)
    #plt.show()

    fig.tight_layout()

    plot_file = file.replace(".csv", ".pdf")
    if file == plot_file:
        raise Exception("Filename is strange, no .csv file ending?")

    pc.savefig(plot_file)

