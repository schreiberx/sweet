#! /usr/bin/env python3

import numpy as np
import sys
import math

from sweet.postprocessing.SphereDataPhysical import *



class SphereDataPhysicalDiff:

	def __init__(self, filename_a = None, filename_b = None):

		if filename_b != None:
			self.compute_diff(filename_a, filename_b)

		pass



	def compute_diff(
			self,
			filename_a,
			filename_b
		):

		file_a = SphereDataPhysical(filename_a)
		file_b = SphereDataPhysical(filename_b)

		self.norm_l1_value = 0.0
		self.norm_l2_value = 0.0
		self.norm_linf_value = 0.0
		self.norm_rms_value = 0.0

		size_ref_j = len(file_a.data)
		size_ref_i = len(file_a.data[0])
		size_cmp_j = len(file_b.data)
		size_cmp_i = len(file_b.data[0])

		multiplier_j = (size_ref_j+1)/(size_cmp_j+1)
		multiplier_i = (size_ref_i+1)/(size_cmp_i+1)


		print ("Dimensions of reference solution: ", size_ref_i, size_ref_j)
		print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
		if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
			print ("Grids are not aligned")
			print ("Try to use (TODO) interpolation script")
			print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
			print ("Multipliers: ", multiplier_i, multiplier_j)
			raise Exception("Grids not properly aligned")

		multiplier_j = int(multiplier_j)
		multiplier_i = int(multiplier_i)

		print("Using multipliers (int): ", multiplier_i, multiplier_j)

		for j in range(0, size_cmp_j):
			for i in range(0, size_cmp_i):
				value = file_b.data[j,i]-file_a.data[j*multiplier_j,i*multiplier_i]

				# http://mathworld.wolfram.com/L1-Norm.html
				self.norm_l1_value += abs(value)
				# http://mathworld.wolfram.com/L2-Norm.html
				self.norm_l2_value += value*value
				# http://mathworld.wolfram.com/L-Infinity-Norm.html
				self.norm_linf_value = max(abs(value), self.norm_linf_value)

				# http://mathworld.wolfram.com/Root-Mean-Square.html
				self.norm_rms_value += value*value

		# Compute sqrt() for l2 norm
		self.norm_l2_value  = math.sqrt(self.norm_l2_value)

		# Divide by 1/sqrt(N)
		self.norm_rms_value /= math.sqrt(size_cmp_i*size_cmp_j)

		#
		# Warning! We normalize here with the number of samples!
		# This doesn't really follow the definition of the Ln norms
		#
		self.norm_l1_value /= (size_cmp_i*size_cmp_j)
		self.norm_l2_value /= (size_cmp_i*size_cmp_j)



	def print(self):

		print("")
		print(" + Warning: L1 and L2 norm are normalized here to be able to compare different resolutions!")
		print(" + norm l1: "+str(self.norm_l1_value))
		print(" + norm l2: "+str(self.norm_l2_value))
		print(" + norm linf: "+str(self.norm_linf_value))
		print(" + norm rms: "+str(self.norm_rms_value))



	def write_file(
			self,
			picklefile,
			tagname = None
		):

		#
		# If picklefile is specified, write norm data to pickle file.
		# This can be later on further postprocessed!
		#
		if picklefile != None:
			import pickle

			if tagname != None:
				tagname += '.'
			else:
				tagname = ''

			pickle_data = {
				tagname+'norm_l1' : self.norm_l1_value,
				tagname+'norm_l2' : self.norm_l2_value,
				tagname+'norm_linf' : self.norm_linf_value,
				tagname+'norm_rms' : self.norm_rms_value,
			}

			print(" + picklefile: "+str(picklefile))

			with open(picklefile, 'wb') as f:
				# Pickle the 'data' dictionary using the highest protocol available.
				pickle.dump(pickle_data, f)

		print("")



if __name__ == "__main__":
	if len(sys.argv) <= 3:
		print("")
		print("Usage:")
		print("	"+sys.argv[0]+" [infile A] [infile B] [picklefile output file (optional)] [reference tagname]")
		print("")
		print("	infile A:")
		print("		First input .csv file with physical space data on the sphere")
		print("")
		print("	infile B:")
		print("		Second input .csv file with physical space data on the sphere")
		print("")
		print("	picklefile:")
		print("		If given, output is pickled into this file")
		print("		diff.error_l1")
		print("		diff.error_l2")
		print("		diff.error_linf")
		print("		diff.error_rms")
		print("")
		print(" reference tagname:")
		print("		How to name value in .pickle file")
		print("")
		sys.exit(1)



	filename_a = sys.argv[1]
	filename_b = sys.argv[2]

	picklefile = None
	if len(sys.argv) >= 3:
		picklefile = sys.argv[3]

	tagname = None
	if len(sys.argv) > 4:
		tagname = sys.argv[4]

	s = SphereDataPhysicalDiff()
	s.compute_diff(filename_a, filename_b)
	s.print()
	s.write_file(picklefile, tagname)
