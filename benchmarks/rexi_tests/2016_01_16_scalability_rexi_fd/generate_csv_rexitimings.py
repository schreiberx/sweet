#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket

import inspect
dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(dir)

import params

output_file_prefix = 'summary'
if len(sys.argv) > 1:
	output_file_prefix = sys.argv[1]


error_output = {}

def extract_errors(output):
	match_list = [
		'DIAGNOSTICS ANALYTICAL RMS H:',
		'DIAGNOSTICS ANALYTICAL RMS U:',
		'DIAGNOSTICS ANALYTICAL RMS V:',
		'DIAGNOSTICS ANALYTICAL MAXABS H:',
		'DIAGNOSTICS ANALYTICAL MAXABS U:',
		'DIAGNOSTICS ANALYTICAL MAXABS V:',
		'Simulation time (seconds):',
		'STOPWATCH preprocessing:',
		'STOPWATCH broadcast:',
		'STOPWATCH reduce:',
		'STOPWATCH solve_rexi_terms:'
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
	for n in params.res_list:
			output_csv = output_file_prefix+"_time_a"+str(a)+"_n"+str(n).zfill(4)+".csv"

			print "Writing output to "+output_csv
			fd = open(output_csv, "w")

			fd.write("#TI Threads vs. simulation time N^2="+str(n)+"^2, DT="+str(params.max_time)+", A="+str(a)+"\n")
			fd.write("#TX Simulation type\n")
			fd.write("#TY number of threads (compact)\n")
			fd.write("threads\\SimType\t"+"\tseconds\trexi_preproc\trexi_broadcast\trexi_reduce\trexi_solver\n")

			for threads in params.thread_list:
				for r in params.mpi_ranks:

					total_threads = r*threads
					if total_threads > params.max_total_threads:
						continue

					for i in params.tests:

						id = i[0]
						parameters = i[1]
						title = i[2]
						compile = i[3]

						#fd.write(str(total_threads)+" cores ("+str(r)+" ranks / "+str(threads)+" threads) "+timing_name)
						fd.write(i[0]+" ("+str(threads).zfill(5)+" threads / "+str(r).zfill(5)+" ranks) "+str(total_threads)+" cores")

						# don't create mpi-parallelized jobs for non-mpi parallelized tests
						if id[0:3] == 'nr_':
							if r > 1:
								fd.write("\t-")
								continue


						file_id = id+"_t"+str(threads).zfill(3)+"_n"+str(n).zfill(4)+"_r"+str(r).zfill(4)+"_a"+str(a)

						output_filename = "run_"+file_id+".txt"

						try:
							content = open(output_filename, 'r').read()

							vals = extract_errors(content)

							h_rms = vals[0]
							seconds = vals[6]
							time_bcast = vals[7]
							time_preproc = vals[8]
							time_reduce = vals[9]
							time_rexi_term = vals[10]

	#						if len(params.tests) > 1:
	#							error_output[output_filename] = ''

							for timing_id in [6, 7, 8, 9, 10]:
								fd.write("\t"+str(vals[timing_id]))
						except:
							fd.write("\tE")
							if len(params.tests) > 1:
								error_output[output_filename] = ''
							

						fd.write("\n")




if len(error_output) != 0:
	output = "\n".join(error_output)
	print "ERROR OUTPUT"
	print output
	print
	print "ERROR OUTPUT regarding the following script files"
	print output.replace('.txt', '.sh')
	print


