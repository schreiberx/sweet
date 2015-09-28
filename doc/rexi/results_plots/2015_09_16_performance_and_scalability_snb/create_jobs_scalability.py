#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket
import params



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



for i in params.tests:
	id = i[0]
	parameters = i[1]
	title = i[2]
	compile_c = i[3]

	for a in params.A_list:
		for threads in params.T_list:
			for n in params.N_list:
				compile_command = compile_c+" --numa-block-allocator="+str(a)+" "

				file_id = id+"_t"+str(threads).zfill(3)+"_n"+str(n).zfill(4)+"_a"+str(a)
				script_filename = params.curdir_name+"/run_"+file_id+".sh"
				output_filename = params.curdir_name+"/run_"+file_id+".txt"

				command = "./build/"+file_id+" "
				command += params.default_params
				command += ' -C '+str(params.cfl)
				command += ' -N '+str(n)
				command += ' -U '+str(params.hyperviscosity[n])
				command += parameters

				print("Creating script "+script_filename)
				fd = open(script_filename, "w")

				fd.write("""#!/bin/bash
#SBATCH -o """+output_filename+"""
#SBATCH -D """+os.getcwd()+"""
#SBATCH -J """+file_id+"""
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=end
#SBATCH --mail-user=martin.schreiber@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00

source /etc/profile.d/modules.sh
module load fftw/mpi/3.3
make clean


. ~/bin/local_vars.sh
. ~/bin/intel_vars.sh

cd ../../../

"""+compile_command+""" --program-binary-name="""+file_id+"""

declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
# this ignored hyperthreads
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS="""+str(threads)+"""

"""+command+""" #> """+output_filename+"""
""")


