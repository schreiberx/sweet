#! /usr/bin/env python3

import math
import sys
import os


class SWEETClusterOptions:
	def __init__(self):
		self.target_machine = ''
		#self.target_machine = 'yellowstone'

		self.total_max_cores = 4096
		self.max_cores_per_node = 16
		self.total_max_nodes = self.total_max_cores/self.max_cores_per_node


		# OpenMP threads in space
		self.par_space_threads = 1

		# MPI threads i time
		self.par_mpi_time_threads = 1



	def getScriptHeader(self):
		if self.par_mpi_time_threads != 1:
			if self.max_cores_per_node % self.par_space_threads != 0:
				raise ValueError('Number of cores on node not evenly dividable by space threads')

		real_time_threads = self.par_mpi_time_threads

		# total number of used MPI ranks
		total_cores = self.par_space_threads*self.par_mpi_time_threads
		mpi_ranks_total = self.par_mpi_time_threads
		mpi_ranks_per_node = math.floor(self.max_cores_per_node/self.par_space_threads)

		sweetdir=os.path.normpath(os.getcwd()+'/../../')

		content = "#!/bin/bash\n"

		if self.target_machine == '':
			content += "\n"
		elif self.target_machine == 'yellowstone':
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
			raise Exception("TODO")

		else:
			print("Target machine "+str(self.target_machine)+" not supported")
			sys.exit(1)

		return content
