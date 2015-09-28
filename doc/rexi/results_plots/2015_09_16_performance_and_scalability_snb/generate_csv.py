#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket

import params

output_file_prefix = 'summary'
if len(sys.argv) > 1:
	output_file_prefix = sys.argv[1]

def extract_errors(output):
	match_list = [
		'DIAGNOSTICS ANALYTICAL RMS H:',
		'DIAGNOSTICS ANALYTICAL RMS U:',
		'DIAGNOSTICS ANALYTICAL RMS V:',
		'DIAGNOSTICS ANALYTICAL MAXABS H:',
		'DIAGNOSTICS ANALYTICAL MAXABS U:',
		'DIAGNOSTICS ANALYTICAL MAXABS V:',
		'Simulation time (seconds):'
	]

	vals = ["x" for i in range(len(match_list))]

	ol = output.splitlines(True)
	for o in ol:
		o = o.replace('\n', '')
		o = o.replace('\r', '')
		for i in range(0, len(match_list)):
			m = match_list[i]
			if o[0:len(m)] == m:
				vals[i] = o[len(m)+1:]

	return vals



for a in params.A_list:
	for n in params.N_list:

		output_csv = output_file_prefix+"_time_a"+str(a)+"_n"+str(n).zfill(4)+".csv"

		print "Writing output to "+output_csv
		fd = open(output_csv, "w")

		fd.write("#TI Threads vs. simulation time N^2="+str(n)+"^2, DT="+str(params.max_time)+", A="+str(a)+"\n")
		fd.write("#TX Simulation type\n")
		fd.write("#TY number of threads (compact)\n")
		fd.write("threads\\SimType\t"+"\t".join([i[0] for i in params.tests])+"\n")
		for threads in params.T_list:
			fd.write(str(threads))
			for i in params.tests:
				id = i[0]
				parameters = i[1]
				title = i[2]
				compile = i[3]

				file_id = id+"_t"+str(threads).zfill(3)+"_n"+str(n).zfill(4)+"_a"+str(a)

				output_filename = "run_"+file_id+".txt"

				try:
					content = open(output_filename, 'r').read()

					vals = extract_errors(content)

					h_rms = vals[0]
					seconds = vals[6]
					fd.write("\t"+str(seconds))
				except:
					fd.write("\tE")
					

			fd.write("\n")



#for threads in [params.T_list[-1]]:
for threads in [params.T_list[0]]:
	for a in params.A_list:
		output_csv = output_file_prefix+"_error_a"+str(a)+"_n"+str(n).zfill(4)+".csv"

		print "Writing output to "+output_csv
		fd = open(output_csv, "w")

		fd.write("#TI Error vs. resolution N^2="+str(n)+"^2\n")
		fd.write("#TX Simulation type\n")
		fd.write("#TY Resolution\n")
		fd.write("N\\SimType\t"+"\t".join([i[0] for i in params.tests])+"\n")

		for n in params.N_list:
			fd.write(str(n))
			for i in params.tests:
				id = i[0]
				parameters = i[1]
				title = i[2]
				compile = i[3]

				file_id = id+"_t"+str(threads).zfill(3)+"_n"+str(n).zfill(4)+"_a"+str(a)

				output_filename = "run_"+file_id+".txt"
#				print output_filename

				try:
					content = open(output_filename, 'r').read()

					vals = extract_errors(content)

					h_rms = float(vals[0])
					fd.write("\t"+str(h_rms))
				except:
					fd.write("\tE")
					

			fd.write("\n")
