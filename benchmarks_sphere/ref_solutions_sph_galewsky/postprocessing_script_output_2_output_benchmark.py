#! /usr/bin/env python3

import glob
import subprocess
import os
import io
import sys
from subprocess import Popen, PIPE


# Load simulation time to get reference file
from contextlib import redirect_stdout

# redirect output
with io.StringIO() as buf, redirect_stdout(buf):
	import jobs_create_scripts as jc
	# ignore output
	#output = buf.getvalue()


simtime = jc.p.runtime.simtime


if jc.p.runtime.bench_id == 100:
	simtimef = simtime/(60.0*60.0)
else:
	simtimef = simtime

t = ("%8.8f" % simtimef).zfill(20)


varname = None
if len(sys.argv) < 2:
	print("Usage: "+sys.argv[0]+" ['h' or 'vort' or 'potvort']")

varname = sys.argv[1]

datafile="output_prog_"+varname+"_t"+t+".csv"




solver_groups = [
	['l1', 1],	# Group name / convergence order
	['l2', 2],

	['ln1', 1],
	['ln2', 2],
	['ln4', 4],

	['ln1test', 1],
	['ln2space', 2]
]



for solver_group_info in solver_groups:

	solver_group = {}
	solver_group['name'] = solver_group_info[0]
	solver_group['conv_order'] = solver_group_info[1]

	dirs = glob.glob("script_"+solver_group['name']+"_*")
	dirs.sort()

	# Skip benchmarks without groups
	if len(dirs) == 0:
		print("INFO: No data found for solver_group.name '"+str(solver_group['name'])+"'")
		continue


	ref_dir=""
	for i in glob.glob("script_"+solver_group['name']+"_ref*"):
		if os.path.isdir(i):
			ref_dir=i
			break

	if ref_dir == "":
		raise Exception("Reference solution not found")

	solver_group['ref_solution'] = ref_dir

	print("DATA: solver_group.name: "+str(solver_group['name']))
	print("DATA: solver_group.conv_order: "+str(solver_group['conv_order']))
	print("DATA: solver_group.ref_solution: "+str(solver_group['ref_solution']))
	print("INFO:\tSIMNAME\tL1\tL2\tLinf\tCONV")

	prev_test_name = ""

	#
	# Determine test groups.
	# Each group has the same TS method
	#
	test_group_dirs = []
	for rundir in dirs:
		if not os.path.isdir(rundir):
			continue

		# Skip reference solution
		if rundir.find('_ref') >= 0:
			continue

		# Reset convergence test?
		# detect new group (time integration method)
		pos = rundir.find("_C")
		test_name = rundir[0:pos]

		if test_name != prev_test_name:
			test_group_dirs.append([])
		prev_test_name = test_name

		test_group_dirs[-1].append(rundir)


	for g in test_group_dirs:

		conv_test = []
		prev_conv_value = 0.0

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

			fs = rundir+"/output.out"
			if os.path.exists(fs):
				f=open(fs, 'r')
			else:
				fs = rundir+".out"
				f=open(fs, 'r')

			tag = " + SimulationBenchmarkTimings.main_simulationloop: "
			secs = "-1"
			for l in f:
				if l[0:len(tag)] == tag:
					secs = l[len(tag):]
					if len(secs) == 0:
						raise Exception("len of secs == 0")

			if secs[-1] == '\n':
				secs = secs[0:-1]

			result += "\t"+str(secs)

			print("DATA: "+result)

