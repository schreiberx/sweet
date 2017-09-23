#! /usr/bin/env python3

import socket
import math
import sys
import os
import multiprocessing


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

		#
		# Setup default values
		#

		# OpenMP threads in space
		self.par_space_cores = 1

		# Cores in time
		self.par_time_cores = 1



	def setupTargetMachine(
		self,
		target_machine
	):
		self.target_machine = target_machine

		if self.target_machine == '':
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

		elif self.target_machine == "cheyenne":
			self.cores_per_node = 36

			# REAL number:
			# self.total_max_nodes = 4032

			# Low number to avoid accidentally wasting computing time
			self.total_max_nodes = 128
			self.total_max_cores = self.cores_per_node*self.total_max_nodes

		else:
			print("Unknown Target: "+str(self.target_machine))
			print("Using default values")

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
	def getScriptHeader(self, jobid, runtimeOptions, dirname):
		#if self.par_mpi_time_threads != 1:
		#	if self.max_cores_per_node % self.par_space_threads != 0:
		#		raise ValueError('Number of cores on node not evenly dividable by space threads')

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

			###########################################################
			# Compact in time
			###########################################################

			pm_cores_per_mpi_rank = self.pm_time_cores_per_mpi_rank * self.pm_space_cores_per_mpi_rank

			if self.pm_time_cores_per_mpi_rank != 1:
				raise Exception("Not yet supported")

			if self.pm_space_cores_per_mpi_rank > self.cores_per_node:
				raise Exception("Not yet supported")



			###########################################################
			# SPACE
			###########################################################

			# Cores per space group on each node
			space_ranks_per_node = int(math.ceil(self.cores_per_node/pm_cores_per_mpi_rank))

			if space_ranks_per_node == 0:
				raise Exception("Too many cores per MPI rank")

			# Total number of space mpi ranks
			space_num_ranks = int(math.ceil(self.par_space_cores/space_ranks_per_node))

			# Total number of nodes INCLUDING IDLING ONES !!!
			space_total_num_nodes = int(math.ceil(space_num_ranks/space_ranks_per_node))

			# Total number of cores in space INCLUDING IDLING ONES !!!
			space_total_num_cores = space_total_num_nodes*self.cores_per_node



			###########################################################
			# TIME
			###########################################################
			time_ranks_per_node = int(math.ceil(18/self.pm_time_cores_per_mpi_rank))
			time_num_cores = self.par_time_cores
			time_num_ranks = self.par_time_cores

			#time_total_num_nodes = 1
			#time_total_num_cores = 1



			###########################################################
			###########################################################
			###########################################################

			# Total number of nodes
			par_total_cores = space_num_ranks*self.pm_space_cores_per_mpi_rank * time_num_ranks*self.pm_time_cores_per_mpi_rank

			# Total number of nodes
			num_nodes = space_ranks_per_node * time_ranks_per_node

			# Number of cores (CPUs) per node
			num_cores_per_node = self.cores_per_node

			# Number of OpenMP threads to use per MPI threads
			# TODO: Include compile options and runtime options to determine this number
			num_omp_threads_per_mpi_thread = pm_cores_per_mpi_rank


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


		elif self.target_machine == 'cheyenne':

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
			content += "# TARGET MACHINE: "+self.target_machine
			content += """#
## project code
#PBS -A NCIS0002
## regular limit: 12 hours
## economy queue
#PBS -q economy
## shared queue
######PBS -q share
## wall-clock time (hrs:mins:secs)
#PBS -l walltime=00:30:00
## select one chunk with one CPU in it
#PBS -l select="""+str(num_nodes)+""":ncpus="""+str(num_cores_per_node)+""":mpiprocs="""+str(num_cores_per_node)+""":ompthreads="""+str(num_omp_threads_per_mpi_thread)+"""
#
#PBS -N """+jobid[0:100]+"""
#PBS -o """+cwd+"/"+dirname+"""/output.out
#PBS -e """+cwd+"/"+dirname+"""/output.err

export OMP_NUM_THREADS="""+str(num_omp_threads_per_mpi_thread)+"""

"""

			#
			# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs/omplace-and-dplace
			#
			mpi_exec_prefix = "mpiexec_mpt -n "+str(par_total_cores)+" "
			# TODO: This seems to make trouble
			#mpi_exec_prefix += " omplace -vv "
			mpi_exec_prefix += " dplace -s 1 "

		else:
			content = ""
			mpi_exec_prefix = ""

		return content, mpi_exec_prefix


