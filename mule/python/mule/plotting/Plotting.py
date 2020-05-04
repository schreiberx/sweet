#! /usr/bin/env python3

import math
import sys


import matplotlib
matplotlib.use('Agg')

from mule.InfoError import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def mule_plotting_usetex(usetex = True):
    import matplotlib
    matplotlib.rcParams['text.usetex'] = usetex


print("Warning: mule.plotting activates 'tex' per default")
mule_plotting_usetex(True)


class Plotting(InfoError):

    def __init__(self):
        InfoError.__init__(self, 'Plotting')

        self.reset()


    def reset(self):
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        self.markers = []
        for m in Line2D.markers:
            try:
                # Filter out some markers
                if m in [' ', ',', '.', 'None', None]:
                    continue

                self.markers.append(m)

            except TypeError:
                pass

        self.linestyles = ['-', '--', ':', '-.']
        self.fillstyle = 'full'
        self.fillstyle = 'none'

        plt.clf()
        plt.cla()
        plt.close()


    def _get_color(self, i):
        return self.colors[i % len(self.colors)]

    def _get_marker(self, i):
        return self.markers[i % len(self.markers)]

    def _get_linestyle(self, i):
        return self.linestyles[i % len(self.linestyles)]


    def plot_setup(
            self,

            figsize = (8,8/1.61803398875),

            title = None,
            subtitle = None,

            xlabel = None,
            xscale = "linear",
            xlim = None,

            yscale = "linear",
            ylabel = None,
            ylim = None,

            annotate = False,
            annotate_each_nth_value = 3,
            annotate_fontsize = False,
            annotate_text_template = None,

            bars_annotate_with_values = False,
            bars_annotate_with_labels = False,
            bars_filled = False,

            legend = True,
            legend_fontsize = None,

            grid = False,

            tight_layout = True,
            outfile = None,

            lambda_fun = None,
        ):

        if isinstance(figsize, float) or isinstance(figsize, int):
            self.figsize = (figsize,figsize/1.61803398875)
        else:
            self.figsize = figsize

        self.title = title
        self.subtitle = subtitle

        self.xlabel = xlabel
        self.xscale = xscale
        self.xlim = xlim

        self.ylabel = ylabel
        self.yscale = yscale
        self.ylim = ylim

        self.annotate = annotate
        self.annotate_each_nth_value = annotate_each_nth_value
        self.annotate_fontsize = annotate_fontsize
        self.annotate_text_template = annotate_text_template

        self.bars_annotate_with_values = bars_annotate_with_values
        self.bars_annotate_with_labels = bars_annotate_with_labels
        self.bars_filled = bars_filled


        self.legend = legend
        self.legend_fontsize = legend_fontsize

        self.grid = grid

        self.tight_layout = tight_layout

        self.outfile = outfile

        self.lambda_fun = lambda_fun


    def plot_start(self):

        self.reset()

        self.fig, self.ax = plt.subplots(figsize=self.figsize)


        if self.title != None:
            plt.title(self.title)

        if self.subtitle != None:
            plt.suptitle(self.subtitle)

        if self.xlabel != None:
            plt.xlabel(self.xlabel)
        if self.xscale != None:
            self.ax.set_xscale(self.xscale, nonposx='clip')
        if self.xlim != None:
            self.ax.set_xlim(self.xlim[0], self.xlim[1])

        if self.ylabel != None:
            plt.ylabel(self.ylabel)
        if self.yscale != None:
            self.ax.set_yscale(self.yscale, nonposy='clip')
        if self.ylim != None:
            self.ax.set_ylim(self.ylim[0], self.ylim[1])

        if self.grid != False:
            #plt.grid(True, which="both", ls="-", color='0.65')
            plt.grid(True, ls="-", color='0.5')

    
        if self.lambda_fun != None:
            self.lambda_fun(self)


    def plot_finish(self):
        if self.legend:
            plt.legend(fontsize=self.legend_fontsize)

        if self.tight_layout:
            plt.tight_layout()

        if self.outfile != None:
            self.info("Plotting to '"+self.outfile+"'")
            plt.savefig(self.outfile)
        else:
            plt.show()




class Plotting_ScatteredData(Plotting):
    def __init__(self):
        Plotting.__init__(self)

    def plot_data(
            self,
            data_plotting,
        ):
        """
        Parameters:
        data_plotting: dict
            key:    string which is also used to sort the data
                'label': Label to use for legend
                'x_values': x coordinates
                'y_values': y coordinates
        """

        c = 0
        for key in sorted(data_plotting.keys()):
            data = data_plotting[key]

            marker = self._get_marker(c)
            linestyle = self._get_linestyle(c)
            color = self._get_color(c)

            if 'label' in data:
                label = data['label']
            else:
                label = key

            x_values = data['x_values']
            y_values = data['y_values']

            self.ax.plot(
                    x_values,
                    y_values,
                    marker=marker,
                    linestyle=linestyle,
                    label=label,
                    fillstyle=self.fillstyle,
                )

            if self.annotate:
                x = x_values
                y = y_values

                px = x[:]
                py = y[:]

                px = px[0::self.annotate_each_nth_value]
                py = py[0::self.annotate_each_nth_value]

                if len(px) % 2 == 0:
                    px.append(x[-1])
                    py.append(y[-1])

                for i, txt in enumerate(px):
                    if self.annotate_text_template != None:
                        text = self.annotate_text_template.format(px[i], py[i])
                        self.ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=self.annotate_fontsize)
                    else:
                        text = "{:.1f}" % (px[i])
                        self.ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=self.annotate_fontsize)

            c += 1
    

    def plot(
            self,
            data_plotting,
            **kwargs
        ):
        self.plot_setup(**kwargs)
        self.plot_start()
        self.plot_data(data_plotting)
        self.plot_finish()



class Plotting_Bars(Plotting):
    def __init__(self):
        Plotting.__init__(self)

    def plot_data_from_tabledata(
            self,
            data_table,
        ):

        #self.fig, self.ax = plt.subplots(figsize=figsize)

        # Size of table data
        sx = len(data_table[0])-1
        sy = len(data_table)-1

        data = [[data_table[iy+1][ix+1] for ix in range(sx)] for iy in range(sy)]

        bar_names = data_table[0][1:]
        print("Bar labels: "+str(bar_names))

        group_names = [data_table[j][0] for j in range(1, sy+1)]
        print("Group names: "+str(group_names))

        # Number of rows
        print("Number of groups: "+str(len(group_names)))

        num_bars_in_group = len(bar_names)
        if num_bars_in_group == 0:
            print("No groups found")
            sys.exit(1)

        print("Number of bars in each group: "+str(num_bars_in_group))
        
        # adjust depending on number of bar_names
        group_width = 1.0/len(group_names)

        if self.bars_filled == True:
            group_width *= 0.9
        else:
            # More space if there's no filling in bars
            group_width *= 0.8

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        bar_width = group_width/len(bar_names)

        # relative position of first bar
        if len(bar_names) % 2:
            # odd
            group_first_bar_rel_position = -(len(bar_names)-1)//2*bar_width
        else:
            # even
            group_first_bar_rel_position = -(len(bar_names)-1)//2*bar_width + bar_width*0.5

        # Center of bars in each group
        group_bars_center = [(i+0.5)/len(group_names) for i in range(len(group_names))]

        # Squeeze together
        group_bars_center = [x for x in group_bars_center]

        self.ax.set_xticks(group_bars_center)
        self.ax.set_xticklabels(group_names)

        self.ax.set_ylabel('Wallclock time (seconds)')

        for tick in self.ax.get_xticklabels():
            tick.set_rotation(45)

        rects = []
        # iterate over all runs
        for j in range(len(bar_names)):

            # scatter plot bars for each run in different colors 'colors'
            bars_x = group_bars_center[:]
            bars_x = [x + group_first_bar_rel_position for x in bars_x]
            bars_x = [x + bar_width*j for x in bars_x]
            bars_height = [data[i][j] for i in range(len(group_names))]

            if self.bars_filled:
                rects_ = self.ax.bar(
                        bars_x,
                        bars_height,
                        bar_width,
                        color=colors[j % len(colors)],
                    )
            else:
                rects_ = self.ax.bar(
                        bars_x,
                        bars_height,
                        bar_width,
                        color=colors[j % len(colors)],
                        fill=False
                    )

            rects.append(rects_)


        if self.bars_annotate_with_values or self.bars_annotate_with_labels:
            for j in range(len(bar_names)):
                for rect in rects[j]:
                    pos_x = rect.get_x() + rect.get_width()/2

                    if self.bars_annotate_with_values:
                        height = rect.get_height()
                        pos_y = height+math.log(1.0+0.1*height)
                        text = "%.2f" % height
                        self.ax.text(pos_x, pos_y, text, ha='center', va='bottom', size=8, rotation=90)

                    if self.bars_annotate_with_labels:
                        pos_y = 1e-10
                        if self.ylim != None:
                            pos_y = self.ylim[0]*1.5
                        text = bar_names[j]
                        self.ax.text(pos_x, pos_y, text, ha='center', va='bottom', size=8, rotation=90)




    

    def gen_plot_from_tabledata(
            self,
            data_plotting,
            outfile=None,
            **kwargs
        ):
        self.plot_setup(outfile=outfile, **kwargs)
        self.plot_start()
        self.plot_data_from_tabledata(data_plotting)
        self.plot_finish()





if __name__ == "__main__":

    p = Plotting()

