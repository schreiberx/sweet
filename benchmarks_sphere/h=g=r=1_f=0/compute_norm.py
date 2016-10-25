#! /usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys


def loadDataFromFile(filename):
	data = np.loadtxt(filename, skiprows=3)
	labelsx = data[0,1:]
	labelsy = data[1:,0]
	data = data[1:,1:]
	return data

ref_data = loadDataFromFile(sys.argv[1])

for filename in sys.argv[2:]:

	cmp_data = loadDataFromFile(filename)

	max_norm_value = -1.0
	for j in range(0, len(cmp_data)):
		for i in range(0, len(cmp_data[0])):
			max_norm_value = max(max_norm_value, abs(cmp_data[j,i]-ref_data[j,i]))

	print(filename+"\t"+str(max_norm_value))
