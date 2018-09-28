

import platform
import socket
import sys
import os



def platform_autodetect():
	fqdn = socket.getfqdn()
	if ".cheyenne" in fqdn:
		return True
	else:
		return False
	

#
# Return an identifier string if system detected, otherwise 'None'
#
def platform_id():
	return 'cheyenne'



#
# Header for scripts
# These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.
#
def job_script_header(jobgeneration, dirname):

	# total number of used MPI ranks
	par_total_cores = self.par_space_cores*self.par_time_cores
	par_mpi_ranks_total = self.par_time_cores
	par_mpi_ranks_per_node = math.floor(self.cores_per_node/self.par_space_cores)

	cwd = os.getcwd()

	for i in range(6):
		if i == 5:
			raise Exception("Unable to find SWEET main directory")

		sweetdir=os.path.normpath(os.getcwd()+("/.."*i))
		if os.path.exists(sweetdir+'/local_software'):
			break

	###########################################################
	# SWEET specific allocation starts here
	###########################################################
	#
	# Number of cores to use for space parallelization
	# -> self.par_space_cores
	#
	# Number of cores to use for time parallelization
	# -> self.par_time_cores
	#

	###########################################################
	# Parallelization model
	###########################################################
	#
	# Number of cores per MPI rank
	# -> self.pm_time_cores_per_mpi_rank
	# -> self.pm_space_cores_per_mpi_rank
	#

	#
	# SETUP variables which are shared by all job scripts
	#

	#if False:
	if True:
		print(" > par_space_cores: "+str(self.par_space_cores))
		print(" > par_time_cores: "+str(self.par_time_cores))

		###########################################################
		# Compact in time
		###########################################################

		pm_cores_per_mpi_rank = self.pm_time_cores_per_mpi_rank * self.pm_space_cores_per_mpi_rank
		print(" + pm_cores_per_mpi_rank: "+str(pm_cores_per_mpi_rank))

		#if self.pm_time_cores_per_mpi_rank != 1:
		#	raise Exception("Not yet supported")

		if self.pm_space_cores_per_mpi_rank > self.cores_per_node:
			raise Exception("Not yet supported")



		###########################################################
		# SPACE
		###########################################################

		# Cores per space group on each node
		space_ranks_per_node = int(math.ceil(self.cores_per_node/pm_cores_per_mpi_rank))
		print(" + space_ranks_per_node: "+str(space_ranks_per_node))

		if space_ranks_per_node == 0:
			raise Exception("Too many cores per MPI rank")

		# Total number of space mpi ranks
		space_num_ranks = int(math.ceil(self.par_space_cores/space_ranks_per_node))
		print(" + space_num_ranks: "+str(space_num_ranks))

		# Total number of nodes INCLUDING IDLING ONES !!!
		space_total_num_nodes = int(math.ceil(space_num_ranks/space_ranks_per_node))
		print(" + space_total_num_nodes: "+str(space_total_num_nodes))

		# Total number of cores in space INCLUDING IDLING ONES !!!
		space_total_num_cores = space_total_num_nodes*self.cores_per_node
		print(" + space_total_num_cores: "+str(space_total_num_cores))



		###########################################################
		# TIME
		###########################################################
		time_ranks_per_node = int(math.ceil(18/self.pm_time_cores_per_mpi_rank))

		time_ranks_per_node = min(time_ranks_per_node*self.pm_time_cores_per_mpi_rank, self.par_time_cores)
		time_num_ranks = self.par_time_cores

		print(" + space_ranks_per_node: "+str(space_ranks_per_node))
		print(" + time_ranks_per_node: "+str(time_ranks_per_node))
		print(" + space_num_ranks: "+str(space_num_ranks))
		print(" + time_num_ranks: "+str(time_num_ranks))



		###########################################################
		###########################################################
		###########################################################

		# Total number of nodes
		par_total_cores = space_num_ranks*self.pm_space_cores_per_mpi_rank * time_num_ranks * self.pm_time_cores_per_mpi_rank
		print(" + par_total_cores: "+str(par_total_cores))

		# Total number of nodes
		num_nodes = int(math.ceil(par_total_cores/self.cores_per_node))
		print(" + num_nodes: "+str(num_nodes))

		# Number of cores (CPUs) per node
		num_cores_per_node = self.cores_per_node

		# Number of OpenMP threads to use per MPI threads
		# TODO: Include compile options and runtime options to determine this number
		num_omp_threads_per_mpi_thread = pm_cores_per_mpi_rank

		# Ranks per node
		num_ranks_per_node = int(math.ceil(self.cores_per_node / pm_cores_per_mpi_rank))

		# Total number of MPI ranks
		mpi_ranks_total = space_num_ranks*time_num_ranks

		print(" + mpi_ranks_total: "+str(mpi_ranks_total))


	else:


		###########################################################
		# Compact in space
		###########################################################
		raise Exception("Not implemented yet")

		###########################################################
		# Compact in time
		###########################################################

		pm_cores_per_mpi_rank = self.pm_time_cores_per_mpi_rank * self.pm_space_cores_per_mpi_rank

		if self.pm_time_cores_per_mpi_rank != 1:
			raise Exception("Not yet supported")

		if self.pm_space_cores_per_mpi_rank > self.cores_per_node:
			raise Exception("Not yet supported")



		###########################################################
		# TIME
		###########################################################

		# Cores per time group on each node
		time_ranks_per_node = int(math.ceil(self.cores_per_node/self.pm_cores_per_mpi_rank))

		if time_ranks_per_node == 0:
			raise Exception("Too many cores per MPI rank")

		# Total number of time mpi ranks
		time_num_ranks = int(math.ceil(self.par_time_cores/time_ranks_per_node))

		# Total number of nodes INCLUDING IDLING ONES !!!
		time_total_num_nodes = int(math.ceil(time_num_ranks/time_ranks_per_node))

		# Total number of cores in time INCLUDING IDLING ONES !!!
		time_total_num_cores = time_total_num_nodes*self.cores_per_node



		###########################################################
		# SPACE
		###########################################################
		space_ranks_per_node = 1
		space_num_cores = 1

		###########################################################

		# Total number of nodes
		par_total_cores = self.par_space_cores*time_total_num_cores

		# Total number of nodes
		num_nodes = int(math.ceil(par_total_cores/self.cores_per_node))

		# Number of cores (CPUs) per node
		num_cores_per_node = self.cores_per_node

		# Number of OpenMP threads to use per MPI threads
		# TODO: Include compile options and runtime options to determine this number
		num_omp_threads_per_mpi_thread = self.pm_time_cores_per_mpi_rank*self.pm_space_cores_per_mpi_rank

		# Total number of MPI ranks
		mpi_ranks_total = par_total_cores



	max_wallclock_seconds_str = str(datetime.timedelta(seconds=self.max_wallclock_seconds)).zfill(8)

	mkl = False
	if compileOptions.mkl == 'enable':
		mkl = True


	#
	# CHEYENNE:
	#  - Dual socket (18 cores / socket)
	#  - 36 cores in total per node
	#
	# 9-D enhanced hypercube topology
	# 100-Gbps link bandwidth — 0.5 μs latency
	# 36 TB/s bisection bandwidth
	#

	content = "#!/bin/bash\n"
	content += "# TARGET MACHINE: "+jobgeneration.cluster.target_machine+"\n"
	content += """#
## project code
#PBS -A NCIS0002
## regular limit: 12 hours
## economy queue
#PBS -q economy
## shared queue
######PBS -q share
## wall-clock time (hrs:mins:secs)
#PBS -l walltime="""+max_wallclock_seconds_str+"""
## select: number of nodes
## ncpus: number of CPUs per node
## mpiprocs: number of ranks per node
#PBS -l select="""+str(num_nodes)+""":ncpus="""+str(num_cores_per_node)+""":mpiprocs="""+str(num_ranks_per_node)+""":ompthreads="""+str(num_omp_threads_per_mpi_thread)+"\n"

	#"default": 2301000 
	#"turbo": 2301000
	#"rated": 2300000
	#"slow": 1200000
	if jobgeneration.cluster.force_turbo_off:
		content += "#PBS -l select=cpufreq=2300000\n"

	content += """#
#PBS -N """+jobid[0:100]+"""
#PBS -o """+cwd+"/"+dirname+"""/output.out
#PBS -e """+cwd+"/"+dirname+"""/output.err

export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""

module load openmpi
"""+("module load mkl" if mkl else "")+"""

"""

	#
	# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs/omplace-and-dplace
	#
	self.tmp_mpi_exec_prefix = "mpirun -n "+str(mpi_ranks_total)+" "

	# TODO: This seems to make trouble
	#mpi_exec_prefix += " omplace -vv "
	#mpi_exec_prefix += " dplace -s 1 "
#	mpi_exec_prefix += " omplace "
#	mpi_exec_prefix += " -nt "+str(num_omp_threads_per_mpi_thread)+" "
#	mpi_exec_prefix += " -vv "
	# -tm pthreads didn't make any difference in performance for single-threaded programs which 1 thread per socket
	#mpi_exec_prefix += " -tm pthreads "

	fn = 'compile_cheyenne.sh'
	f = open(fn, 'w')
	f.write("#! /bin/bash\n")
	f.write("\n")
#	f.write("SWEETROOT=\""+dirpath+"/../../../\"\n")
#	f.write("cd \"$SWEETROOT\"\n")
#	f.write("\n")

#	if jobgeneration.cluster.target_machine == 'cheyenne_impi':
#		# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/intel-mpi-and-open-mpi
#		f.write("module load impi\n")
#		f.write("\n")
#
#	if jobgeneration.cluster.target_machine == 'cheyenne_openmpi':
#		# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/intel-mpi-and-open-mpi
#		f.write("module load openmpi\n")
#		f.write("\n")
#
#	if jobgeneration.compile.mkl == 'enable':
#		f.write("module load mkl\n")
#		f.write("\n")

	f.write("\n")
	f.write("scons "+jobgeneration.compile.getSConsParams()+'\n')
	f.write("\n")
	os.chmod(fn, 0o755)

	return content


def job_script_mpiexec(jobgeneration):
	return self.tmp_mpi_exec_prefix


#
# Prefix before the executable is executed
#
def job_script_exec_prefix(jobgeneration):
	content = ""

	content += """
#
# Cheyenne (job_script_exec_prefix)
#
SCONS="scons """+jobgeneration.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS
"""

	return content


