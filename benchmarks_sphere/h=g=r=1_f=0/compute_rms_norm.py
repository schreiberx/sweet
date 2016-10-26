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

	norm_value = 0.0

	size_j = len(cmp_data)
	size_i = len(cmp_data[0])
	for j in range(0, size_j):
		for i in range(0, size_i):
			value = cmp_data[j,i]-ref_data[j,i]
			norm_value += value*value

	norm_value /= sqrt(size_i*size_j)

	print(filename+"\t"+str(norm_value))
