#! /usr/bin/env python3

import glob
import subprocess
import os
import sys
from subprocess import Popen, PIPE


# Load simulation time to get reference file
import jobs_create as jc
simtime = jc.p.runtime.simtime
print(simtime)
t = ("%8.8f" % simtime).zfill(20)

datafile="output_prog_h_t"+t+".csv"




groups = [
	['l1', 1],	# Group name / convergence order
	['l2', 2],

	['ln1', 1],
	['ln2', 2],
	['ln4', 4],

	['ln1test', 1],
	['ln2space', 2]
]



for group_info in groups:
	group = group_info[0]

	ref_dir=""
	for i in glob.glob("script_"+group+"_ref*"):
		if os.path.isdir(i):
			ref_dir=i
			break

	print("Group:"+str(group_info))
	print("Using reference: "+str(ref_dir))

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

			# remove line break
			if result[-1] == '\n':
				result = result[0:-1]

			f=open(rundir+"/output.out")
			tag = "Simulation time (seconds): "
			secs = 0
			for l in f:
				if l[0:len(tag)] == tag:
					secs = l[len(tag):]

			if secs[-1] == '\n':
				secs = secs[0:-1]

			result += "\t"+str(secs)

			print(result)

