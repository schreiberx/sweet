#! /usr/bin/env python3

import math
import sys

from mule.InfoError import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import matplotlib
print("Warning: mule.plotting activates 'tex' per default")
matplotlib.rcParams['text.usetex'] = True


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
			xlabel = None,
			xscale = "linear",
			yscale = "linear",
			ylabel = None,

			annotate = False,

			legend = True,
			legend_fontsize = None,

			tight_layout = True,
			outfile = None,
		):
		self.figsize = figsize

		self.title = title

		self.xlabel = xlabel
		self.xscale = xscale

		self.ylabel = ylabel
		self.yscale = yscale

		self.annotate = annotate,

		self.legend = legend
		self.legend_fontsize = legend_fontsize

		self.tight_layout = tight_layout

		self.outfile = outfile

	def plot_start(self):

		self.reset()

		self.fig, self.ax = plt.subplots(figsize=self.figsize)

		self.ax.set_xscale(self.xscale, nonposx='clip')
		self.ax.set_yscale(self.yscale, nonposy='clip')

		if self.title != None:
			plt.title(self.title)

		if self.xlabel != None:
			plt.xlabel(self.xlabel)

		if self.ylabel != None:
			plt.ylabel(self.ylabel)


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
			annotation = False,
		):
		"""
		Parameters:
		data_plotting: dict
			key:	string which is also used to sort the data
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
					label=label
				)

			if annotation:
				x = x_values
				y = y_values

				px = x[:]
				py = y[:]

				px = px[0::annotation_render_each_nth_value]
				py = py[0::annotation_render_each_nth_value]

				if len(px) % 2 == 0:
					px.append(x[-1])
					py.append(y[-1])

				for i, txt in enumerate(px):

					text = "%.1f" % (px[i])
					self.ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=8)
					if False:
						if mode == 'dt':
							#self.ax.annotate(text, (px[i]*1.03, py[i]*0.92), fontsize=8)
							self.ax.annotate(px[i], (px[i]*1.03, py[i]*0.92), fontsize=8)
						elif mode == 'wallclocktime':
							text = "%.1f/%.1f" % (py[i], px[i])
							self.ax.annotate(text, (px[i]*1.03, py[i]*1.03), fontsize=8)

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
			xlabel = None,
			ylabel = None,
			title = None,
			subtitle = None,
			xscale = None,
			yscale = None,
			annotate_bars_with_values = False,
			annotate_bars_with_labels = False,
			legend = True,
			filled_bars = True,
			ylim = None,
			figsize=(8,8/1.61803398875)
		):

		self.fig, self.ax = plt.subplots(figsize=figsize)

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

		if filled_bars == True:
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

			if filled_bars:
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


		if annotate_bars_with_values or annotate_bars_with_labels:
			for j in range(len(bar_names)):
				for rect in rects[j]:
					pos_x = rect.get_x() + rect.get_width()/2

					if annotate_bars_with_values:
						height = rect.get_height()
						pos_y = height+math.log(1.0+0.1*height)
						text = "%.2f" % height
						self.ax.text(pos_x, pos_y, text, ha='center', va='bottom', size=8, rotation=90)

					if annotate_bars_with_labels:
						pos_y = 1e-10
						if ylim != None:
							pos_y = ylim[0]*1.5
						text = bar_names[j]
						self.ax.text(pos_x, pos_y, text, ha='center', va='bottom', size=8, rotation=90)

		if ylim != None:
			self.ax.set_ylim(ylim[0], ylim[1])

		if legend:
			self.ax.legend([rect[0] for rect in rects], bar_names)

		if xscale != None:
			self.ax.set_xscale(xscale, nonposx='clip')

		if yscale != None:
			self.ax.set_yscale(yscale, nonposy='clip')

		if title != None:
			plt.title(title)

		if xlabel != None:
			plt.xlabel(xlabel)

		# Doesn't work :-(
		#if subtitle != None:
		#	plt.suptitle(subtitle)
		#	plt.subplots_adjust(top=0.85)

		if ylabel != None:
			plt.ylabel(ylabel)

		self.fig.tight_layout()



	def plot_data_annotated(
			self,
			data_plotting,
			xlabel = None,
			ylabel = None,
			title = None,
			xscale = "linear",
			yscale = "linear",
			annotation_render_each_nth_value = 3
		):
		raise Exception("TODO")

	

	def gen_plot_from_tabledata(
			self,
			data_plotting,
			outfile=None,
			**kwargs
		):
		self.plot_start(**kwargs)
		self.plot_data_from_tabledata(data_plotting, **kwargs)
		self.plot_finish(outfile, **kwargs)





if __name__ == "__main__":

	p = Plotting()

