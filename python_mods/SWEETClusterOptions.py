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
			raise Exception("TODO")

		elif self.target_machine == "yellowstone":
			self.total_max_cores = 4096
			self.cores_per_node = 16
			raise Exception("TODO")

		elif self.target_machine == "cheyenne":
			self.total_max_cores = 4096
			self.cores_per_node = 36

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
		retval += '_space'+str(self.par_space_cores)
		retval += '_time'+str(self.par_time_cores)

		return retval



	def getScriptHeader(self):
		#if self.par_mpi_time_threads != 1:
		#	if self.max_cores_per_node % self.par_space_threads != 0:
		#		raise ValueError('Number of cores on node not evenly dividable by space threads')

		real_time_threads = self.par_time_cores

		# total number of used MPI ranks
		total_cores = self.par_space_cores*self.par_time_cores
		mpi_ranks_total = self.par_time_cores
		mpi_ranks_per_node = math.floor(self.cores_per_node/self.par_space_cores)

		for i in range(6):
			if i == 5:
				raise Exception("Unable to find SWEET main directory")

			sweetdir=os.path.normpath(os.getcwd()+("/.."*i))
			if os.path.exists(sweetdir+'/local_software'):
				break

		content = "#!/bin/bash\n"
		content += "# TARGET MACHINE: "+self.target_machine

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
#BSUB -n """+str(mpi_ranks_total)+"""	 number of tasks in job
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

		elif self.target_machine == 'cheyenne':

			content += """#
## project code
#PBS -A NCIS0002
## regular limit: 12 hours
## economy queue
#PBS -q economy
## shared queue
##PBS -q share
## wall-clock time (hrs:mins:secs)
#PBS -l walltime=01:00:00
## select one chunk with one CPU in it
#PBS -l select=1:ncpus=1
#
#PBS -N """+jobid+"""
#PBS -o """+cwd+"/"+jobid+""".out
#PBS -e """+cwd+"/"+jobid+""".err

export MPI_USE_ARRAY=false

"""


		return content
