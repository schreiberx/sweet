#! /usr/bin/env python3

from SWEET import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

class SWEETPlotting:

	def __init__(self):
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


	def get_color(self, i):
		return self.colors[i % len(self.colors)]

	def get_marker(self, i):
		return self.markers[i % len(self.markers)]

	def get_linestyle(self, i):
		return self.linestyles[i % len(self.linestyles)]


	def plot_scattered(	
			self,
			data_plotting,
			x_label = None,
			y_label = None,
			title = None,
			outfile = None
		):

		plt.clf()
		plt.cla()
		plt.close()

		fig, ax = plt.subplots(figsize=(8,6))

		ax.set_xscale("linear")
		ax.set_yscale("linear")

		if title != None:
			plt.title(title)

		if x_label != None:
			plt.xlabel(x_label)

		if y_label != None:
			plt.ylabel(y_label)


		c = 0
		for key, values in data_plotting.items():
			marker = self.get_marker(c)
			linestyle = self.get_linestyle(c)
			color = self.get_color(c)

			label = key
			x_values = values['x_values']
			y_values = values['y_values']

			ax.plot(x_values, y_values, marker=marker, linestyle=linestyle, label=label)
			c += 1

		plt.legend()

		if outfile != None:
			plt.savefig(outfile)
		else:
			plt.show()




if __name__ == "__main__":

	p = SWEETPlotting()

