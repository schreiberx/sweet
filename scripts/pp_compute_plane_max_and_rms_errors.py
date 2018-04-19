#! /usr/bin/env python3

import numpy as np
import sys
import math

sys.stdout.write('.')

def loadDataFromFile(filename):
	global prefix

	try:
		data = np.loadtxt(filename, skiprows=0, ndmin=2)
	except:
		prefix = filename if len(sys.argv) <= 3 else sys.argv[3]
		print(prefix+": UNABLE TO OPEN '"+filename+"'")
		sys.exit(1)

	#labelsx = data[0,0:]
	#labelsy = data[0:,0]
	#data = data[1:,1:]
	return data

ref_data = loadDataFromFile(sys.argv[1])
cmp_data = loadDataFromFile(sys.argv[2])
prefix = sys.argv[3]

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0

(size_ref_j, size_ref_i) = ref_data.shape
(size_cmp_j, size_cmp_i) = cmp_data.shape

multiplier_j = (size_ref_j)/(size_cmp_j)
multiplier_i = (size_ref_i)/(size_cmp_i)


if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
	print ("Dimensions of reference solution: ", size_ref_i, size_ref_j)
	print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
	print ("Multipliers: ", multiplier_i, multiplier_j)
	print ("Grids are not aligned")
	sys.exit(1)

multiplier_j = int(multiplier_j)
multiplier_i = int(multiplier_i)
#print ("Multipliers (int): ", multiplier_i, multiplier_j)

for j in range(0, size_cmp_j):
	for i in range(0, size_cmp_i):
		#print("(",i,",",j,",", i*multiplier_i,",", j*multiplier_j,")", end="")

		value = cmp_data[j,i]-ref_data[j*multiplier_j,i*multiplier_i]

		norm_l1_value += abs(value)
		norm_l2_value += value*value
		norm_linf_value = max(norm_linf_value, abs(value))

norm_l1_value = norm_l1_value/(size_cmp_i*size_cmp_j)
norm_l2_value = math.sqrt(norm_l2_value/(size_cmp_i*size_cmp_j))

#
# L1, L2, Linf
#
print(prefix+"\t"+str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value))
