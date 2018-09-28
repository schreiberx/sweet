#! /usr/bin/env python3

import socket
import math
import sys
import os
import multiprocessing
import datetime

from SWEETPlatforms import *
from InfoError import *



class SWEETParallelizationDimOptions(InfoError):
	"""
	Parallelization information for one dimension
	"""


	def __init__(self, dim_name = ''):
		"""
		Constructor

		Parameters
		----------
		dim_name: string
			Information on this dimension for error output
		"""

		InfoError.__init__(self, "SWEETParallelizationDimOptions."+dim_name)

		# Name of dimension for debugging purpose
		self.dim_name = dim_name

		# overall number of cores
		self.num_cores = None

		# number of cores per MPI rank
		self.num_cores_per_rank = None

		# number of cores per MPI rank
		self.num_threads_per_rank = None

		# Number of ranks
		self.num_ranks = None



	def print(self):
		self.hline()
		self.info("dim_name: "+str(self.dim_name))
		self.info("num_cores: "+str(self.num_cores))
		self.info("num_cores_per_rank: "+str(self.num_cores_per_rank))
		self.info("num_threads_per_rank: "+str(self.num_threads_per_rank))
		self.info("num_ranks: "+str(self.num_ranks))



	def setup(self, dim_id = None):
		"""
		Setup and infer some values automatically.

		In case of some values impossible to infer, generate error
		"""

		if self.dim_name == '':
			if dim_id != None:
				self.dim_name = "DIM"+str(dim_id)
				InfoError.setup(self, "SWEETParallelizationDimOptions."+self.dim_name)

		if self.num_threads_per_rank == None:
			# Use same number of thraeds as there are cores per rank per default
			self.num_threads_per_rank = self.num_cores_per_rank

		if self.num_cores == None:
			if self.num_cores_per_rank != None and self.num_ranks != None:
				self.num_cores = self.num_cores_per_rank*self.num_ranks

		if self.num_cores == None:
			self.error("num_cores not specified")

		if self.num_cores_per_rank == None:
			self.error("num_cores_per_rank not specified")

		if self.num_threads_per_rank == None:
			self.p_error("num_threads_per_rank not specified")

		if self.num_ranks == None:
			self.p_error("num_ranks not specified")
		
