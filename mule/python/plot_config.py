#! /usr/bin/env python3

#
# Author: Martin Schreiber <martin.schreiber@tum.de>
#
# Creation Date: Somewhen 2019
#
# Changelog:
# 2020-12-09: Changes for priorizing colors suitable for color vision defected people
# 2019-01-14: Minor updates & Cleanups (Martin)
# 2019-02-01: Added markevery option
#

#
# There is an example given below which plots a picture if this file is executed directly
#

import matplotlib.pyplot as plt

def setup(
        scale = 1.0,    # Scaling factor for default image size
                        # Use only this to resize your image.
        nrows = 1,      # Number of rows in plot
        ncols = 1,      # Number of colums in plot
        figsize = None, # Size of figure
        square_figsize = False,
        title = True,
        uselatex = False
    ):
    """
    Setup the plotting with default parameters which are optimized for PDF DinA4 plots

    Use the 'scale' parameter to enlarge or shrink the picture if required
    """

    from matplotlib import rc

    default_figsize = (3, 2)
    if figsize == None:
        figsize = (default_figsize[0]*ncols, default_figsize[1]*nrows)

    if square_figsize:
        _ = max(figsize)
        figsize = (_, _)

    if uselatex:
        rc('text', usetex=True)
        rc('pgf', texsystem='pdflatex')
        rc('font', family='serif')
        rc('font', size=10)
        rc('axes', labelsize=10)
        rc('axes', titlesize=10)
        rc('figure', titlesize=12)

    #
    # Important: Optimize this only for PDF output!
    #
    if 1:
        rc('figure', dpi = 300)
        rc('font', size = 6)
        rc('legend', fontsize = 5)
        rc('axes', titlesize=6)
        rc('figure', titlesize=7)

    if 1:
        # PNG output dpi
        rc('savefig', dpi = 300)


    figsize = (figsize[0]*scale, figsize[1]*scale)

    # If there's no title, decrease by 10%
    if title == False:
        figsize = (figsize[0], figsize[1]*0.9)

    # Start new plot
    plt.close()
    return plt.subplots(nrows, ncols, figsize=figsize)


def savefig(
        filename,
        fig = None,
        overwrite_file = True,
        verbose = 0,
        **kwargs
    ):
    """
    Homebrew method to save figure

    Has some nice features such as removing particular tags in PDFs to ensure bytewise reproducibility
    """


    if not overwrite_file:
        import os
        if os.path.exists(filename):
            raise Exception("Output file '"+filename+"' already exists")

    if ".pdf" in filename:
        # Remove all metadata from PDF output file
        if 'metadata' in kwargs:
            raise Exception("TODO")
        else:
            kwargs['metadata'] = {'Creator': None, 'Producer': None, 'CreationDate': None}

    if verbose > 0:
        print("Saving image to figure '"+filename+"'")

    if fig == None:
        plt.savefig(filename, **kwargs)
    else:
        plt.savefig(filename, **kwargs)


def show():
    plt.show()


class PlotStyles:
    """
    Class which provides a variety of plot styles

    There are three different styles, for lines, markers and colors

    Using different lines and markers is in particular helpful for colorblind people.
    """

    def __init__(self):

        self.markers = []

        # We start with all variants of crosses 
        self.markers += ['+', 'x', '1', '2', '3', '4']

        # dot
        self.markers += ['.']

        # Triangles
        self.markers += [4, 5, 6, 7, 8, 9, 10, 11]

        self.linestyles = [
            #(?, (line, spacing)
            (0, (1, 1)),

            (0, (8, 1)),
            (0, (5, 1)),

            (0, (3, 2, 5, 2)),
            (0, (3, 1, 1, 1)),

            (0, (6, 1, 6, 1, 6, 1)),
        ]

        #self.colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

        # New order of colors which has two colorblind suitable colors first
        self.colors = ['#b4b400', '#00b4b4', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

        self.reset()

    def reset(self):
        """
        Reset all style counters
        """
        self.c_colors = 0
        self.c_markers = 0
        self.c_linestyles = 0

    def set_color_counter(self, c):
        """
        Set style for color counter

        This is e.g. helpful to group some lines to have the same color
        """
        self.c_colors = c

    def set_marker_counter(self, c):
        """
        Set the marker style counter
        """
        self.c_markers = c

    def set_linestyle_counter(self, c):
        """
        Set the line style counter
        """
        self.c_linestyles = c


    def getNextStyle(self, num_points=None, num_markers=15):
        """
        Return a set of parameters which are ready to be used
        e.g. for plt.plot(...) with the ** Python feature, e.g.
        plt.plot(x, y, **ps.getNextStyle(), label="f(x) = cos(x)")
        """
        retval = {}

        retval['color'] = self.colors[self.c_colors % len(self.colors)]
        retval['marker'] = self.markers[self.c_markers % len(self.markers)]
        retval['linestyle'] = self.linestyles[self.c_linestyles % len(self.linestyles)]

        if num_points != None:
            retval['markevery'] = num_points//num_markers
            if retval['markevery'] == 0:
                retval['markevery'] = 1

        self.c_colors += 1
        self.c_markers += 1
        self.c_linestyles += 1

        return retval



if __name__ == "__main__":
    import numpy as np

    # Import plot_config
    import plot_config as pc

    #
    # Call setup() for each new plot
    # Using scaling factor in case you need larger plots
    #
    # The default scaling of 1.0 is optimized for smaller
    # plots on presentations if you store them as a .pdf file
    #
    fig, ax = pc.setup(scale=1.5)

    #
    # Get a handler to different plotting styles
    #
    ps = pc.PlotStyles()


    x = np.linspace(0, 1, 80)

    y = np.sin(x*10)


    #
    # The next plotting style can be loaded with ps.getNextStyle()
    # This returns a dictionary of parameters and values
    #
    # In order to use it as a parameter, we use the ** prefix
    #
    ax.plot(x, y, **ps.getNextStyle(), label="f(x) = sin(x*10)")

    y = np.cos(x*10)
    ax.plot(x, y, **ps.getNextStyle(), label="f(x) = cos(x*10)")

    for i in range(10):
        s = 0.3
        y = np.cos(x*10+i*s)

        if 1:
            # Directly use plot style
            plt.plot(x, y, **ps.getNextStyle(), label="f(x) = cos(x*10+"+str(i)+"*"+str(s)+")")
        else:
            # Change plot styles
            pstyle = ps.getNextStyle()
            # E.g. use no markers
            pstyle['marker'] = ''
            plt.plot(x, y, **pstyle, label="f(x) = cos(x*10+"+str(i)+"*"+str(s)+")")

    y = x*1
    ax.plot(x[-len(x)//7:], y[-len(x)//7:], color='black', linestyle="dashed", linewidth=1, label="ref. order 1")

    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.set_title("plot_config example")
    ax.legend()
    fig.tight_layout()
    #outfile = "/tmp/asdf.pdf"
    #print("Saving to '"+outfile+"'")
    #plt.savefig(outfile)
    pc.show()

