#! /usr/bin/env python3

import numpy as np
import sys
import math


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




def loadSphereDataFromFile(filename):
	"""
	Load sphere data stored in physical space from .csv file

	Return:
	-------
		Tuple with (data, longitude angles, latitude angles)
	"""
	print("Loading file: "+filename_a)

	try:
		data = np.loadtxt(filename, skiprows=0)
	except:
		raise Exception("Unable to open '"+filename+"'")
		sys.exit(1)

	# First row and col are longitude and latitude coordinates
	labelsx = data[0,0:]
	labelsy = data[0:,0]
	data = data[1:,1:]
	return data, labelsx, labelsy



filename_a = sys.argv[1]
filename_b = sys.argv[2]

picklefile = None
if len(sys.argv) >= 3:
	picklefile = sys.argv[3]

tagname = ''
if len(sys.argv) > 4:
	tagname = sys.argv[4]

file_a_data, file_a_labelsx, file_a_labelsy = loadSphereDataFromFile(filename_a)
file_b_data, file_b_labelsx, file_b_labelsy = loadSphereDataFromFile(filename_b)

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0
norm_rms_value = 0.0

size_ref_j = len(file_a_data)
size_ref_i = len(file_a_data[0])
size_cmp_j = len(file_b_data)
size_cmp_i = len(file_b_data[0])

multiplier_j = (size_ref_j+1)/(size_cmp_j+1)
multiplier_i = (size_ref_i+1)/(size_cmp_i+1)


if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
	print ("Dimensions of reference solution: ", size_ref_i, size_ref_j)
	print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
	print ("Multipliers: ", multiplier_i, multiplier_j)
	print ("Grids are not aligned")
	print ("Try to use (TODO) interpolation script")
	sys.exit(1)

multiplier_j = int(multiplier_j)
multiplier_i = int(multiplier_i)

print("Using multipliers (int): ", multiplier_i, multiplier_j)

for j in range(0, size_cmp_j):
	for i in range(0, size_cmp_i):
		value = file_b_data[j,i]-file_a_data[j*multiplier_j,i*multiplier_i]

		# http://mathworld.wolfram.com/L1-Norm.html
		norm_l1_value += abs(value)
		# http://mathworld.wolfram.com/L2-Norm.html
		norm_l2_value += value*value
		# http://mathworld.wolfram.com/L-Infinity-Norm.html
		norm_linf_value = max(abs(value), norm_linf_value)

		# http://mathworld.wolfram.com/Root-Mean-Square.html
		norm_rms_value += value*value

# Compute sqrt() for l2 norm
norm_l2_value  = math.sqrt(norm_l2_value)

# Divide by 1/sqrt(N)
norm_rms_value /= math.sqrt(size_cmp_i*size_cmp_j)

#
# Warning! We normalize here with the number of samples!
# This doesn't really follow the definition of the Ln norms
#
norm_l1_value /= (size_cmp_i*size_cmp_j)
norm_l2_value /= (size_cmp_i*size_cmp_j)


print("")
print(" + Warning: L1 and L2 norm are normalized here to be able to compare different resolutions!")
print(" + norm l1: "+str(norm_l1_value))
print(" + norm l2: "+str(norm_l2_value))
print(" + norm linf: "+str(norm_linf_value))
print(" + norm rms: "+str(norm_rms_value))

#
# If picklefile is specified, write norm data to pickle file.
# This can be later on further postprocessed!
#
if picklefile != None:
	import pickle

	if tagname != '':
		tagname += '.'

	pickle_data = {
		tagname+'norm_l1' : norm_l1_value,
		tagname+'norm_l2' : norm_l2_value,
		tagname+'norm_linf' : norm_linf_value,
		tagname+'norm_rms' : norm_rms_value,
	}

	print(" + picklefile: "+str(picklefile))

	with open(picklefile, 'wb') as f:
		# Pickle the 'data' dictionary using the highest protocol available.
		pickle.dump(pickle_data, f)

print("")
