#! /usr/bin/env python3

from SWEET import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

class SWEETPostprocessingPlotting(InfoError):

	def __init__(self):
		InfoError.__init__(self, 'SWEETPostprocessingPlotting')

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


	def __get_color(self, i):
		return self.colors[i % len(self.colors)]

	def __get_marker(self, i):
		return self.markers[i % len(self.markers)]

	def __get_linestyle(self, i):
		return self.linestyles[i % len(self.linestyles)]


	def plot_start(self):
		self.reset()


	def plot_finish(
			self,
			outfile = None,
		):

		if outfile != None:
			self.info("Plotting to '"+outfile+"'")
			plt.savefig(outfile)
		else:
			plt.show()

	def plot_scattered_data(
			self,
			data_plotting,
			xlabel = None,
			ylabel = None,
			title = None,
			xscale = "linear",
			yscale = "linear",
		):

		self.fig, self.ax = plt.subplots(figsize=(10,7))

		self.ax.set_xscale(xscale, nonposx='clip')
		self.ax.set_yscale(yscale, nonposy='clip')

		if title != None:
			plt.title(title)

		if xlabel != None:
			plt.xlabel(xlabel)

		if ylabel != None:
			plt.ylabel(ylabel)


		c = 0
		for key, values in data_plotting.items():
			marker = self.__get_marker(c)
			linestyle = self.__get_linestyle(c)
			color = self.__get_color(c)

			label = key
			x_values = values['x_values']
			y_values = values['y_values']

			self.ax.plot(x_values, y_values, marker=marker, linestyle=linestyle, label=label)
			c += 1

		plt.legend()

		

	def plot_scattered(
			self,
			data_plotting,
			outfile=None,
			**kwargs
		):
		self.reset()
		self.plot_scattered_data(data_plotting, **kwargs)
		self.plot_finish(outfile)



if __name__ == "__main__":

	p = SWEETPostprocessingPlotting()

