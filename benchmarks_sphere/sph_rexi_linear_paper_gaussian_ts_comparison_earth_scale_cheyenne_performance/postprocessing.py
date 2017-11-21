#! /usr/bin/env python3

import glob
import subprocess
import os
import sys
from subprocess import Popen, PIPE


# Load simulation time to get reference file
import jobs_create as jc
simtime = jc.p.runtime.simtime


if jc.p.runtime.bench_id == 100:
	simtimef = simtime/(60.0*60.0)
else:
	simtimef = simtime

t = ("%8.8f" % simtimef).zfill(20)

datafile="output_prog_h_t"+t+".csv"




groups = [
	['', 2]
]



for group_info in groups:
	group = group_info[0]

	ref_dir=""
	if group == '':
		for i in glob.glob("script_ref*"):
			if os.path.isdir(i):
				ref_dir=i
				break
	else:
		for i in glob.glob("script_"+group+"_ref*"):
			if os.path.isdir(i):
				ref_dir=i
				break

	print("Group:"+str(group_info))
	print("Using reference: "+str(ref_dir))

	print("SIMNAME\tL1\tL2\tLinf\tCONV")
	if group == '':
		directories = glob.glob("script_*")
	else:
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
			#print(" ".join(progparams))

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
			tag = "Wallclock time (seconds): "
			secs = '-'
			for l in f:
				if l[0:len(tag)] == tag:
					secs = l[len(tag):]

			if len(secs) > 0:
				if secs[-1] == '\n':
					secs = secs[0:-1]

			result += "\t"+str(secs)

			print(result)

