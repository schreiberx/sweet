#! /usr/bin/env python3

import socket
import math
import sys
import os
import multiprocessing
import datetime


class SWEETClusterOptions:


	def __init__(self, target_machine = ''):
		self.target_machine = target_machine

		# Total number of cores without hyperthreading
		self.total_max_cores = -1

		# Number of physical cores per node
		self.cores_per_node = -1

		# Cores per MPI thread
		self.pm_time_cores_per_mpi_rank = 1
		self.pm_space_cores_per_mpi_rank = 1

		self.setupTargetMachine(target_machine)

		self.exec_prefix = ''

		#
		# Setup default values
		#

		# OpenMP threads in space
		self.par_space_cores = 1

		# Cores in time
		self.par_time_cores = 1

		# max wallclock time, default: 1h
		self.max_wallclock_seconds = 60*60

		# Force disabling of turbo mode (if supported)
		self.force_turbo_off = False



	def setupTargetMachine(
		self,
		target_machine
	):
		self.target_machine = target_machine

		auto = False
		if self.target_machine == 'auto':
			auto = True

			# Autodetect host with FQDN
			#hostname = socket.gethostname()
			fqdn = socket.getfqdn()

			self.target_machine = None

			if ".gw4.metoffice.gov.uk" in fqdn:
				self.target_machine = "isambard"

			elif ".yellowstone" in fqdn:
				self.target_machine = "yellowstone"

			elif ".cheyenne" in fqdn:
				self.target_machine = "cheyenne"

			else:
				self.target_machine = socket.gethostname()


		if self.target_machine == 'isambard':
			self.total_max_cores = 4096
			self.total_max_nodes = 512
			self.total_max_cores = self.cores_per_node*self.total_max_nodes
			raise Exception("TODO")

		elif self.target_machine == "yellowstone":
			self.cores_per_node = 16
			self.total_max_nodes = 512
			self.total_max_cores = self.cores_per_node*self.total_max_nodes
			raise Exception("TODO")

		elif self.target_machine == "cheyenne" or self.target_machine == "cheyenne_impi" or self.target_machine == "cheyenne_openmpi":
			self.cores_per_node = 36

			# REAL number:
			# self.total_max_nodes = 4032

			# Low number to avoid accidentally wasting computing time
			self.total_max_nodes = 128
			self.total_max_cores = self.cores_per_node*self.total_max_nodes

		elif self.target_machine == "martinium":
			self.cores_per_node = 2
			self.total_max_nodes = 1
			self.total_max_cores = self.cores_per_node*self.total_max_nodes

		elif self.target_machine == "coolmuc_mpp2":
			# https://www.lrz.de/services/compute/linux-cluster/overview/
			self.cores_per_node = 28

			# REAL number:
			# self.total_max_nodes = 4032

			# Low number to avoid accidentally wasting computing time
			self.total_max_nodes = 384
			self.total_max_cores = self.cores_per_node*self.total_max_nodes


		else:
			if not auto and self.target_machine != '':
				raise Exception("Invalid target machine '"+self.target_machine+"'")

			#print("Invalid target machine '"+self.target_machine+"'")
			#print("Unknown Target: "+str(self.target_machine))
			#print("Using default values")

			self.total_max_cores = multiprocessing.cpu_count()
			self.cores_per_node = self.total_max_cores


		self.total_max_nodes = self.total_max_cores//self.cores_per_node

		if self.total_max_nodes*self.cores_per_node != self.total_max_cores:
			raise Exception("Inconsistency detected")



	def setup(
		par_space_cores,
		par_time_cores,
	):
		self.par_space_cores = par_space_cores
		self.par_time_cores = par_time_cores



	def getUniqueID(self):
		retval = 'MPI'
		retval += '_space'+str(self.par_space_cores).zfill(2).zfill(2)
		retval += '_time'+str(self.par_time_cores).zfill(3)

		return retval



	##################################################################
	# return header for job execution script and the parallel job execution command
	##################################################################
	#
	# \return header, exec_prefix
	#
	def getScriptHeader(self, jobid, runtimeOptions, compileOptions, dirname):
		#if self.par_mpi_time_threads != 1:
		#	if self.max_cores_per_node % self.par_space_threads != 0:
		#		raise ValueError('Number of cores on node not evenly dividable by space threads')


		print("Target: "+self.target_machine)

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
		# SETUP the following variables:
		#
		# content
		#    e.g. "/bin/bash\n#PBS....."
		#
		# mpi_exec_prefix
		#    e.g. mpirun -n XXX
		#
		if self.target_machine == 'yellowstone':
			#
			# YELLOWSTONE:
			# Each node has 16 cores
			# 8 cores per socket
			# hyperthreading enabled
			#
			# More example job scripts:
			# https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/using-computing-resources/running-jobs/platform-lsf-job-script-examples
			#
			content = "#!/bin/bash\n"
			content += "# TARGET MACHINE: "+self.target_machine+"\n"

			content += """
#
# LSF batch script to run an MPI application
#
# YELLOW STONE SPECIFIC!!!
# https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/
#
#BSUB -P NCIS0002	# project code
#BSUB -W 02:00		# wall-clock time (hrs:mins)
#
#BSUB -n """+str(par_mpi_ranks_total)+"""	 number of tasks in job
#BSUB -R "span[ptile=16]"    # run 16 MPI tasks per node
#
#BSUB -outdir """+dirname+"""
#BSUB -J """+job_id+"""	# job name
#BSUB -o """+dirname+""".out  # output file name in which %J is replaced by the job ID
#BSUB -e """+dirname+""".out  # error file name in which %J is replaced by the job ID
#
## https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/using-computing-resources/queues-and-charges
#BSUB -q small
#

"""
			raise Exception("TODO")
			#
			# AND PAY ATTENTION TO NODE SHARING IF USING ONLY SINGLE CORE!!!!!!!!!!!
			#
			mpi_exec_prefix = ""


		elif self.target_machine == "cheyenne_openmpi":
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
			content += "# TARGET MACHINE: "+self.target_machine+"\n"
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
			if self.force_turbo_off:
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
			mpi_exec_prefix = "mpirun -n "+str(mpi_ranks_total)+" "
			# TODO: This seems to make trouble
			#mpi_exec_prefix += " omplace -vv "
			#mpi_exec_prefix += " dplace -s 1 "
		#	mpi_exec_prefix += " omplace "
		#	mpi_exec_prefix += " -nt "+str(num_omp_threads_per_mpi_thread)+" "
		#	mpi_exec_prefix += " -vv "
			# -tm pthreads didn't make any difference in performance for single-threaded programs which 1 thread per socket
			#mpi_exec_prefix += " -tm pthreads "


		elif self.target_machine == "cheyenne_impi":
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
			content += "# TARGET MACHINE: "+self.target_machine+"\n"
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
			if self.force_turbo_off:
				content += "#PBS -l select=cpufreq=2300000\n"

			content += """#
#PBS -N """+jobid[0:100]+"""
#PBS -o """+cwd+"/"+dirname+"""/output.out
#PBS -e """+cwd+"/"+dirname+"""/output.err

export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""

module load impi
"""+("module load mkl" if mkl else "")+"""

"""

			#
			# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs/omplace-and-dplace
			#
			mpi_exec_prefix = "mpirun -n "+str(mpi_ranks_total)+" "
			# TODO: This seems to make trouble
			#mpi_exec_prefix += " omplace -vv "
			#mpi_exec_prefix += " dplace -s 1 "
		#	mpi_exec_prefix += " omplace "
		#	mpi_exec_prefix += " -nt "+str(num_omp_threads_per_mpi_thread)+" "
		#	mpi_exec_prefix += " -vv "
			# -tm pthreads didn't make any difference in performance for single-threaded programs which 1 thread per socket
			#mpi_exec_prefix += " -tm pthreads "


		elif self.target_machine == "cheyenne":

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
			content += "# TARGET MACHINE: "+self.target_machine+"\n"
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
			if self.force_turbo_off:
				content += "#PBS -l select=cpufreq=2300000\n"

			content += """#
#PBS -N """+jobid[0:100]+"""
#PBS -o """+cwd+"/"+dirname+"""/output.out
#PBS -e """+cwd+"/"+dirname+"""/output.err

export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""

"""+("module load mkl" if mkl else "")+"""

"""

			#
			# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs/omplace-and-dplace
			#
			mpi_exec_prefix = "mpiexec_mpt -n "+str(mpi_ranks_total)+" "
			# TODO: This seems to make trouble
			#mpi_exec_prefix += " omplace -vv "
			#mpi_exec_prefix += " dplace -s 1 "
			mpi_exec_prefix += " omplace "
			mpi_exec_prefix += " -nt "+str(num_omp_threads_per_mpi_thread)+" "
			mpi_exec_prefix += " -vv "
			# -tm pthreads didn't make any difference in performance for single-threaded programs which 1 thread per socket
			#mpi_exec_prefix += " -tm pthreads "


		elif self.target_machine == "coolmuc_mpp2":
			#
			# Coolmuc:
			#  - Dual socket (14 cores / socket)
			#  - 28 physical cores in total per node
			#

			content = """#!/bin/bash
#SBATCH -o """+cwd+"/"+dirname+"""/output.out
#SBATCH -e """+cwd+"/"+dirname+"""/output.err
#SBATCH -D """+cwd+"/"+dirname+"""/
#SBATCH -J """+jobid[0:100]+"""
###SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --ntasks="""+str(num_nodes*num_cores_per_node)+"""
#SBATCH --cpus-per-task="""+str(int(num_cores_per_node/num_ranks_per_node))+"""
# the above is a good match for the
# CooLMUC2 architecture.
#SBATCH --mail-type=end 
#SBATCH --mail-user=schreiberx@gmail.com
#SBATCH --export=NONE 
#SBATCH --time="""+max_wallclock_seconds_str+"""
source /etc/profile.d/modules.sh
#PBS -l select="""+str(num_nodes)+""":ncpus="""+str(num_cores_per_node)+""":mpiprocs="""+str(num_ranks_per_node)+""":ompthreads="""+str(num_omp_threads_per_mpi_thread)+"\n"

			#if self.force_turbo_off:
			#	content += "#PBS -l select=freq=rated\n"

			content += """
export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""
export OMP_PROC_BIND=SCATTER

"""+self.environment_vars

			#
			# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs/omplace-and-dplace
			#
			mpi_exec_prefix = "mpirun -n "+str(mpi_ranks_total)+" "
			# TODO: This seems to make trouble


		elif self.target_machine == 'mac-login-intel':

			content = """#!/bin/bash
#SBATCH -o """+cwd+"/"+dirname+"""/output.out
#SBATCH -e """+cwd+"/"+dirname+"""/output.err
#####SBATCH -D /home/hpc/<project>/<user>
#SBATCH -J """+jobid[0:100]+"""
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=end
#SBATCH --mail-user=schreiberx@gmail.com
#SBATCH --export=NONE
###SBATCH --time=01:30:00

source /etc/profile.d/modules.sh

# optimal values for sph sphere T128 REXI stuff
export OMP_PROC_BIND=COMPACT
export OMP_NUM_THREADS=32

#export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""
"""

			# one process per node, only one node
			mpi_exec_prefix = "mpiexec.hydra -ppn 1 -n 1"


		elif self.target_machine == 'mac-login-amd':

			content = """#!/bin/bash
#SBATCH -o """+cwd+"/"+dirname+"""/output.out
#SBATCH -e """+cwd+"/"+dirname+"""/output.err
#####SBATCH -D /home/hpc/<project>/<user>
#SBATCH -J """+jobid[0:100]+"""
#SBATCH --get-user-env
#SBATCH --partition=bdz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=end
#SBATCH --mail-user=schreiberx@gmail.com
#SBATCH --export=NONE
###SBATCH --time=01:30:00

source /etc/profile.d/modules.sh

# optimal values for sph sphere T128 REXI stuff
export OMP_PROC_BIND=SPREAD
export OMP_NUM_THREADS=16

#export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""
"""

			# one process per node, only one node
			mpi_exec_prefix = "mpiexec.hydra -ppn 1 -n 1"


		else:
			content = ""
			content += "#!/bin/bash\n"
			content += "\n"
			#content += "export OMP_PROC_BIND=CLOSE\n"
			content += "\n"
			mpi_exec_prefix = ""

	#	else:
	#		raise Exception("Invalid target machine "+self.target_machine)
 


		if len(mpi_exec_prefix) > 1:
			if mpi_exec_prefix[-1] != " ":
				mpi_exec_prefix += " "

		if self.exec_prefix != '':
			mpi_exec_prefix += self.exec_prefix

		return content, mpi_exec_prefix


