#! /usr/bin/env python3

import glob
import subprocess
import os
import sys
from subprocess import Popen, PIPE


# Load simulation time to get reference file
import jobs_create as jc
simtime = jc.p.runtime.max_simulation_time
print(simtime)
t = ("%8.8f" % simtime).zfill(20)

datafile="output_prog_h_pert_t"+t+".csv"




groups = [
	#['l1', 1],	# Group name / convergence order
	#['l2', 2],

	#['ln1', 1],
	#['ln2', 2],
	#['ln4', 4],

	#['ln1test', 1],
	['ln2space', 2]
]

for group_info in groups:
	group = group_info[0]
	conv_order = group_info[1]

	ref_dir=""
	for i in glob.glob("script_"+group+"_ref*"):
		if os.path.isdir(i):
			ref_dir=i
			break

	print("Group:"+str(group_info))
	print("Using reference: "+str(ref_dir))

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
		print("Running tests for new group:")
		for rundir in g:
			progparams = ['./pp_compute_max_and_rms_errors.py', ref_dir+"/"+datafile, rundir+"/"+datafile, rundir]
			#print (progparams)
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

			result = result.decode()

			# last line contains convergence info
			if result[-1] == '\n':
				result = result[0:-1]

			last_conv_value = float(result.split('\t')[-1])

			if prev_conv_value == 0.0:
				conv_test.append(0.0)
			else:
				conv_test.append(float(last_conv_value)/float(prev_conv_value))

			print(result+"\t"+str(conv_test[-1]))
			prev_conv_value = last_conv_value


		print("Expected convergence: "+str(conv_order))
		print("Measured convergence: "+str(conv_test))
                
		# test these first convergence tests
		if 'ln4' in group_info[0]:
			# Comparing RK4 with 4-th order methods requires a more relaxed test
			# This works and was empirically determined
			test_range = range(2,5)
			max_error_rate = 0.5
		elif 'space' in group_info[0]:
			test_range = range(2,5)
			max_error_rate = 10.00
		elif 'ln2' in group_info[0]:
			test_range = range(1,3)
			max_error_rate = 0.5
		else:
			test_range = range(1,3)
			max_error_rate = 0.05

		print(test_range)

		for i in test_range:
			if conv_test[i] == 0.0:
				if 'ln4' not in group_info[0]:
					print("Invalid convergence of 0")
					sys.exit(1)

			conv_rate = 2**conv_order
			a = abs(conv_test[i]-conv_rate)/conv_rate
			print("Testing convergence "+str(conv_test[i])+": "+str(a))

			if a > max_error_rate:
				if group_info[0] == "ln2space":
					print("ERROR: First tests should converge, but was not the case")
					continue

				print("ERROR: Convergence rate not given")
				print("ERROR: Tested with convergence rate "+str(conv_test[i]))
				print("ERROR: Expected convergence rate "+str(conv_rate))
				print("ERROR: Relative error: "+str(a))
				print("ERROR: Max error: "+str(max_error_rate))
				sys.exit(1)
