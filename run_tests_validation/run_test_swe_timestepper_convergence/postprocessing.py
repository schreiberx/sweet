#! /usr/bin/env python2

import glob
import subprocess
import os
import sys
from subprocess import Popen, PIPE



datafile="output_prog_h_pert_t00000000000.10000000.csv"


groups = [
	['l1', 1],	# Group name / convergence order
	['l2', 2],

	['ln1', 1],
	['ln2', 2],

	['ln1test', 1],
	['ln2test', 2]
]

for group_info in groups:
	group = group_info[0]
	conv_order = group_info[1]

	ref_dir=""
	for i in glob.glob("script_"+group+"_ref*"):
		if os.path.isdir(i):
			ref_dir=i
			break

	print("")
	print("Using reference: "+ref_dir)

	prev_conv_value = 0.0

	print("SIMNAME\tL1\tL2\tLinf\tCONV")
	directories = glob.glob("script_"+group+"_*")
	directories.sort()

	prev_test_name = ""

	#
	# Determine test groups, each group has the same TS method
	#
	test_group_dirs = []
	for rundir in directories:
		if not os.path.isdir(rundir):
			continue

		# Skip reference solution
		if rundir.find('_ref') >= 0:
			continue

		# Reset convergence test?
		pos = rundir.find("_C")
		test_name = rundir[0:pos]

		if test_name != prev_test_name:
			test_group_dirs.append([])
		prev_test_name = test_name

		test_group_dirs[-1].append(rundir)



	for g in test_group_dirs:

		conv_test = []
		prev_conv_value = 0.0

		print("")
		print("Running tests for new group")
		for rundir in g:
			progparams = ['./pp_compute_max_and_rms_errors.py', ref_dir+"/"+datafile, rundir+"/"+datafile, rundir]

			p = Popen(progparams, stdout=PIPE, stderr=PIPE)
			output, error = p.communicate()

			if p.returncode != 0:
				print("*"*80)
				print("Exit code "+str(p.returncode))
				print("EXEC: "+(" ".join(progparams)))
				print("STDOUT: "+str(output))
				print("STDERR: "+str(error))
				continue

			result = output

			# last line contains convergence info
			if result[-1] == '\n':
				result = result[0:-1]

			last_conv_value = float(result.split("\t")[-1])

			if prev_conv_value == 0.0:
				conv_test.append(0.0)
			else:
				conv_test.append(float(last_conv_value)/float(prev_conv_value))

			print(result+"\t"+str(conv_test[-1]))
			prev_conv_value = last_conv_value


		print("Expected convergence: "+str(conv_order))
		print("Measured convergence: "+str(conv_test))

		# test these first convergence tests
		for i in range(1,4):
			if conv_test[i] == 0.0:
				print("Invalid convergence of 0")
				sys.exit(1)

			a = abs(1.0-abs(conv_test[i]-conv_order)/conv_order )
			print("Testing convergence "+str(conv_test[i])+": "+str(a))

			if a > 0.05:
				print("ERROR: First tests should converge, but was not the case")
				sys.exit(1)

