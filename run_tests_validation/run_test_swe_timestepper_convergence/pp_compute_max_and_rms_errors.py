#! /usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys
import math


def loadDataFromFile(filename):
	global prefix

	try:
		data = np.loadtxt(filename, skiprows=0)
	except:
		prefix = filename if len(sys.argv) <= 3 else sys.argv[3]
		print(prefix+": UNABLE TO OPEN "+filename)
		sys.exit(1)

	labelsx = data[0,1:]
	labelsy = data[1:,0]
	data = data[1:,1:]
	return data

ref_data = loadDataFromFile(sys.argv[1])
cmp_data = loadDataFromFile(sys.argv[2])
prefix = sys.argv[3]

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0

size_j = len(cmp_data)
size_i = len(cmp_data[0])
for j in range(0, size_j):
	for i in range(0, size_i):
		value = cmp_data[j,i]-ref_data[j,i]

		norm_l1_value += abs(value)
		norm_l2_value += value*value
		norm_linf_value = max(norm_linf_value, abs(value))

norm_l1_value = norm_l1_value/(size_i*size_j)
norm_l2_value = math.sqrt(norm_l2_value/(size_i*size_j))

#
# L1, L2, Linf
#
print(prefix+"\t"+str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value))
