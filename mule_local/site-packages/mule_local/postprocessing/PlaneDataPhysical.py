#! /usr/bin/env python3

import numpy as np


class PlaneDataPhysical:

	def __init__(self, filename = None):

		if filename != None:
			self.read_file(filename)

		pass

	def read_file(self, filename):
		"""
		Load plane data stored in physical space from .csv file
		"""
		print("Loading file: "+filename)

		try:
			data = np.loadtxt(
					filename,
					skiprows=0,	# Don't skip any data rows
					ndmin=2		# Generate at least 2 dimensional array
				)
		except Exception as e:
			raise e

		# First row and col are longitude and latitude coordinates
		self.labelsx = []#data[0,0:]
		self.labelsy = []#data[0:,0]
		self.data = data #[1:,1:]

