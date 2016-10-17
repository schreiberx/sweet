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

fftw_written = False

fftw_wisdom_fileprefix = "FFTW_WISDOM_nofreq_T"


for i in params.tests:
	id = i[0]
	parameters = i[1]
	title = i[2]
	compile_c = i[3]
	id_c = i[4]

	for mpi_ranks in params.mpi_ranks:
		for a in params.A_list:
			for threads in params.thread_list:
				for n in params.res_list:

					# don't create mpi-parallelized jobs for non-mpi parallelized tests
					if id[0:3] == 'nr_':
						if mpi_ranks > 1:
							continue

					total_threads = mpi_ranks * threads

					if total_threads > params.max_total_threads:
						continue

					compile_command = compile_c+" --numa-block-allocator="+str(a)+" "

					# deactivate SWEET threading if geq 1 thread or REXI Parallel sum is enabled
					no_threading = threads <= 1 or '--rexi-parallel-sum=enable' in compile_command
					if no_threading:
						compile_command += " --threading=off"
					else:
						compile_command += " --threading=omp"

					file_id = id+"_t"+str(threads).zfill(3)+"_n"+str(n).zfill(4)+"_r"+str(mpi_ranks).zfill(4)+"_a"+str(a)
					file_id_compile = id_c+"_t"+("no" if no_threading else "yes")+"_a"+str(a)
					script_filename = params.curdir_name+"/run_"+file_id+".sh"
					compile_script_filename = params.curdir_name+"/compile_"+file_id_compile+".sh"
					output_filename = params.curdir_name+"/run_"+file_id+".txt"
					output_filename_err = params.curdir_name+"/run_"+file_id+".err"

					command = "./build/"+file_id_compile+" "
					command += params.default_params
					command += ' -C '+str(params.cfl)
					command += ' -N '+str(n)
					command += ' -U '+str(params.hyperviscosity[n])
					command += parameters



					script_headers = """#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --ntasks="""+str(mpi_ranks)+"""
#SBATCH --cpus-per-task="""+str(threads)+"""
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --time="""+params.timeout+"""

#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS="""+str(threads)+"""
"""

					env_header = """
. /etc/profile.d/modules.sh

module unload gcc
module unload fftw

module unload python
module load python/2.7_anaconda_nompi


module unload intel
module load intel/16.0

module unload mpi.intel
module load mpi.intel/5.1

module load gcc/5

cd """+os.getcwd()+"""
cd ../../../

. local_software/env_vars.sh

"""


					print("Creating compile script "+compile_script_filename)
					fd = open(compile_script_filename, "w")
					fd.write("#! /bin/bash\n\n"+env_header+"\n\n"+compile_command+" --program-binary-name="+file_id_compile+" || exit 1\n")


					if not fftw_written:
						for t in range(0, max(params.thread_list)+1):
							fd = open("fftw_plan_"+str(t).zfill(3)+".sh", "w")
							fd.write("""#! /bin/bash

#SBATCH -o fftw_plan_"""+str(t).zfill(3)+""".txt
#SBATCH -J fftw_plan_"""+str(t).zfill(3)+"""
"""+script_headers+"""
"""+env_header+"""

mpiexec.hydra -genv OMP_NUM_THREADS """+str(t)+""" -envall -ppn 1 ./fftw_gen_wisdoms_all.sh """+str(t)+" "+fftw_wisdom_fileprefix+str(t)+"\n")
						fftw_written = True


					print("Creating script "+script_filename)
					fd = open(script_filename, "w")


					fd.write("""#! /bin/bash

#SBATCH -o """+output_filename+"""
###SBATCH -e """+output_filename_err+"""
#SBATCH -J """+file_id+"""
"""+script_headers+"""

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo

"""+env_header+"""

# force to use FFTW WISDOM data
declare -x SWEET_FFTW_LOAD_WISDOM_FROM_FILE=\""""+fftw_wisdom_fileprefix+("0" if no_threading else str(threads))+"""\"

time -p mpiexec.hydra -genv OMP_NUM_THREADS """+str(threads)+""" -envall -ppn """+str(28/threads)+""" -n """+str(mpi_ranks)+""" """+command+"""

""")
