#! /usr/bin/env python2

import glob
import subprocess
import os


datafile="output_prog_h_t00000000000.10000000.csv"

for group in ['l1', 'l2', 'ln1', 'ln2']:
	ref_dir=""
	for i in glob.glob("script_"+group+"_ref*"):
		if os.path.isdir(i):
			ref_dir=i
			break

	print("")
	print("Using reference: "+ref_dir)

	prev_value = 0

	print("SIMNAME\tL1\tL2\tLinf\tCONV")
	directories = glob.glob("script_"+group+"*")
	directories.sort()
	for dir in directories:
		if not os.path.isdir(dir):
			continue

		try:
			result = subprocess.check_output(['./pp_compute_max_and_rms_errors.py', ref_dir+"/"+datafile, dir+"/"+datafile, dir])
		except:
			print("Error in processing file "+dir)
			continue

		if result[-1] == '\n':
			result = result[0:-1]

		last_value = result.split("\t")[-1]
		conv = (0.0 if prev_value == 0.0 else float(last_value)/float(prev_value))
		print(result+"\t"+str(conv))

		prev_value = float(last_value)
