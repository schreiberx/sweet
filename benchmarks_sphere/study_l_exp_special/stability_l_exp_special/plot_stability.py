#! /usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import sys
import mule.plot_config as pc


import matplotlib
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

    extent = (xs[0], xs[-1], ys[0], ys[-1])

    X, Y = np.meshgrid(xs, ys)

    Z = data

    if 1:
        contour_range = [0.0]

        im = plt.imshow(Z, extent=extent)
        fig.colorbar(im, ax=ax)
        
        cs = plt.contour(X, Y, Z, contour_range, colors="black", linewidth=0.1)
        #pyplot.clabel(cs, cs.levels, fontsize=8, fmt={0: method})
        #ax.clabel(cs, cs.levels, fontsize=8, fmt={1: file})

        l = cs.legend_elements()
        legend_lines += l[0]
        legend_labels += [file]

    c = c+1

    #plt.xlabel("Re(z)")
    #plt.ylabel("Im(z)", labelpad=-5)

    plt.legend(legend_lines, legend_labels)
    plt.tight_layout()
    #plt.show()


    plot_file = file.replace(".csv", ".pdf")
    if file == plot_file:
        raise Exception("Filename is strange, no .csv file ending?")

    pc.savefig(plot_file)

