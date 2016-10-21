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

fftw_written = False


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


					"""
					http://community.hartree.stfc.ac.uk/wiki/site/admin/idataplex%20phase-2%20-%20further%20info.html
					84 nodes, each node has 2 x 12 core Intel Xeon processors, A total of 2,016 cores.

					Interconnect is Infiniband from Mellanox (FDR Connect-IB 56 GB/s).
					42 of 84 have accelerators - an Intel Phi 5110P.

					http://community.hartree.stfc.ac.uk/wiki/site/admin/jobs2.html
					Jobs will by default go into the q1h32 queue.

					Hyperthreading is deactivated!
					24 cores per node available

					# run 24 MPI tasks per node
					#BSUB -R "span[ptile=24]"

					# phase2 system has 2x12 cores per node

					# number of overall MPI tasks
					#BSUB -n 128
					"""

					script_headers = """#BSUB -R "span[ptile="""+str(24/threads)+"""]"
#BSUB -R "affinity[core("""+str(threads)+"""):cpubind=core:distribute=pack]"
#BSUB -R same[type:model]
###BSUB -R order[hosts]
# All on the same rack - TODO: does this really work?
##BSUB -R "cu[type=rack]"
##BSUB -R "cu[maxcus=1]"
# dedicated network
###BSUB -network "type=sn_all: usage=dedicated"
###BSUB -network "type=sn_single: usage=dedicated"
# exclusive resource
#BSUB -x
#BSUB -q idbq
#BSUB -W """+params.timeout+"""
#BSUB -n """+str(mpi_ranks)+"""

#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS="""+str(threads)+"""
"""

					env_header = """
. /etc/profile.d/modules.sh

module unload gcc
# This FFTW does not support OMP!!!
module unload fftw
module load python/2.7.8
module load intel/15.2.164
module load intel_mpi

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
							fd.write("""#BSUB -o fftw_plan_"""+str(t).zfill(3)+""".out
#BSUB -J fftw_plan_"""+str(t).zfill(3)+"""
"""+script_headers+"""
"""+env_header+"""

./fftw_gen_wisdoms_all.sh """+str(t)+"\n")
						fftw_written = True


					print("Creating script "+script_filename)
					fd = open(script_filename, "w")


					fd.write(
"""#BSUB -o """+output_filename+"""
#BSUB -e """+output_filename_err+"""
#BSUB -J """+file_id+"""
"""+script_headers+"""

echo "LSB_BIND_CPU_LIST"
echo "$LSB_BIND_CPU_LIST"

echo "LSB_BIND_MEM_LIST"
echo "$LSB_BIND_MEM_LIST"

echo "LSB_BIND_MEM_POLICY"
echo "$LSB_BIND_MEM_POLICY"

#    RM_CPUTASKn"
echo "RM_MEM_AFFINITY"
echo "$RM_MEM_AFFINITY"

echo "OMP_NUM_THREADS"
echo "$OMP_NUM_THREADS"

echo

"""+env_header+"""

# force to use FFTW WISDOM data
declare -x SWEET_FFTW_LOAD_WISDOM_FROM_FILE="FFTW_WISDOM_T"""+("0" if no_threading else str(threads))+""""

mpiexec.hydra -genv OMP_NUM_THREADS """+str(threads)+""" -envall  -np """+str(mpi_ranks)+""" """+command+"""

""")
