#! /usr/bin/env python3

import glob
import subprocess
import os
import sys
import math

from math import log
from subprocess import Popen, PIPE


def extract_errors(filename):

	file = open(filename, "r")
	output=file.read()
	match_list = [
		'DIAGNOSTICS BENCHMARK DIFF H:',
		'DIAGNOSTICS BENCHMARK DIFF U:',
		'DIAGNOSTICS BENCHMARK DIFF V:'
	]

	vals = ["x" for i in range(3)]

	ol = output.splitlines(True)
	for o in ol:
		o = o.replace('\n', '')
		o = o.replace('\r', '')
		for i in range(0, len(match_list)):
			m = match_list[i]
			if o[0:len(m)] == m:
				vals[i] = o[len(m)+1:]

	return vals


datafile="output_prog_h_pert_t00000000000.10000000.csv"


groups = [
	#['l1', 1],	# Group name / convergence order
	#['l2', 2],

	#['ln1', 1],
	#['ln2', 2],

	#['ln1test', 1],
	#['ln2test', 2],

        ['ln2space',2]
]

for group_info in groups:
	group = group_info[0]
	conv_order = group_info[1]
	print("")
	print("Group to be tested: "+str(group)+" "+str(conv_order))
	print("")
	#List of output files (all methods in group)
	outputs = glob.glob("script_"+group+"*.out")
	outputs.sort()
			
	#
	# Determine test groups, each group has the same TS method
	#
	test_group_methods = []
	prev_test_name = ""
	for output in outputs:
		# Reset convergence test?
		pos = output.find("_N")
		test_name = output[0:pos]
		if test_name != prev_test_name:
			test_group_methods.append(test_name)
                        
		prev_test_name = test_name

		#test_group_methods[-1].append(test_name)

	print("")
	print("Methods to be analysed:"+str(test_group_methods))
	print("")

	for method in test_group_methods:
		# List of outputs for this method
		outputs = glob.glob(method+"*.out")
		outputs.sort()

		conv_test = []
		prev_conv_value = 0.0
		prev_h_error = 0.0
		pos = outputs[1].find("_N")
		test_name = outputs[1][0:pos]
		print("----------------------")
		print("Method: "+test_name)
		print("Resolution    MaxErrorH    MaxErrorU  MaxErrorV  RatioH")

		i = -1
		n = len(outputs)
		ratios = [0 for j in range(n)]

		for output in outputs:
			i = i+1 #output index
			test_res = output[pos+2:len(output)-4]
			result = extract_errors(output)

			if result[0] != "x" and float(result[0]) != 0 :
				ratios[i] = float(prev_h_error)/float(result[0])
			else:
				print("This method could not be evaluated: "+str(result[0]))
				sys.exit(1)

			prev_h_error=result[0]

			print(str(test_res)+"\t"+str(result[0])+"\t"+str(result[1])+"\t"+str(result[2])+"\t"+str(ratios[i]))

		conv_test = math.log(max(ratios), 2)
		print("Expected convergence: "+str(conv_order))
		print("Measured convergence: "+str(conv_test))

		# test these first convergence tests

		a = abs(1.0-abs(conv_test)/conv_order )
		#print a
		if a > 0.05 and conv_test < conv_order:
			print("ERROR: Convergence order not matching expected")
			print(a)
			sys.exit(1)
		else:
			print("*Convergence seems correct or higher than expected!!*")

	print("-----------------------")
