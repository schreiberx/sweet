import math
from functools import reduce
import operator

from InfoError import *
from SWEETPlatformResources import *
from SWEETParallelizationDimOptions import *

__all__ = ['SWEETParallelization']

def _prod(iterable):
	return reduce(operator.mul, iterable, 1)


class SWEETParallelization(InfoError):

	"""
	This class stores information on how each program should be executed on a platform

	The application has to initialize the variables in a sufficient way so that the
	final configuration to execute the program on a cluster can be inferred from this.

	Terminology:
	------------

	'num_threads':
		Number of running threads within one MPI rank

	'num_cores':
		Number of physical processing cores

	'rank':
		MPI rank

	"""
	def __init__(self):
		InfoError.__init__(self, "SWEETParallelization")

		self.reset()

		#
		# WARNIING:
		# Leave these variables here to ensure being not influenced by reset()
		#

		#
		# Disable utilization of `mpiexec` to run job
		# This is required to run e.g. the validation scripts should be
		# (currently) executed on a single node and without MPI support
		self.mpiexec_disabled = False

		# Force disabling of turbo mode (if supported)
		self.force_turbo_off = False


	def reset(self):
		"""
		Reset functionality for a fresh configuration in case that a new setup is triggered
		"""

		# Number of cores per rank
		self.num_cores_per_rank : int = None 

		# Number of threads per rank
		self.num_threads_per_rank = None

		# Number of ranks per node
		self.num_ranks_per_node = None

		# Number of cores per node
		self.num_cores_per_node : int = None 

		# Number of total ranks
		self.num_ranks = None

		# Number of total nodes
		self.num_nodes = None

		# Number of total cores
		self.num_cores = None

		# max wallclock time, default: 1h
		self.max_wallclock_seconds = 60*60

		# number of OpenMP threads for each MPI rank
		self.omp_threads_per_mpi_rank = None


		# List with parallelization information in each dimension
		# Note, that space dimension can and should be treated as a single dimension
		self.pardims = None




	def print(self):
		if self.pardims == None:
			print("No dimension-wise parallelization information specified")
		else:
			for i in self.pardims:
				i.hline()
				i.print()

		self.hline()
		self.hline()
		self.info("num_cores_per_rank: "+str(self.num_cores_per_rank))
		self.info("num_threads_per_rank: "+str(self.num_threads_per_rank))
		self.info("num_ranks_per_node: "+str(self.num_ranks_per_node))
		self.info("num_cores_per_node: "+str(self.num_cores_per_node))
		self.info("num_ranks: "+str(self.num_ranks))
		self.info("num_nodes: "+str(self.num_nodes))
		self.info("num_cores: "+str(self.num_cores))
		self.info("max_wallclock_seconds: "+str(self.max_wallclock_seconds))
		self.info("force_turbo_off: "+str(self.force_turbo_off))



	def dummy_setup_if_no_setup(self, platform_resources : SWEETPlatformResources):
		"""
		Setup a dummy parallelization dimension to use one rank on one node and all cores on the node
		"""

		if self.pardims == None:
			dummy = SWEETParallelizationDimOptions("dummy")
			dummy.num_cores = platform_resources.num_cores_per_node
			dummy.num_cores_per_rank = dummy.num_cores
			dummy.num_threads_per_rank = dummy.num_cores
			dummy.num_ranks = 1
			dummy.omp_threads_per_mpi_rank = dummy.num_threads_per_rank
			self.setup([dummy], platform_resources)
			self.print()



	def setup(self, list_pardims, platform_resources : SWEETPlatformResources):
		"""
		Setup data which is required by the platform specific scripts to
		generate the job scripts

		Parameters
		----------
		list_pardims:	SWEETParallelizationDimOptions
			List with options for parallelization in each dimension

		platform_resources
			reference to jobgeneration class

		#mode : string
		#	'serial': No parallelization
		"""

		self.reset()

		self.pardims = list_pardims

#		# Support space-only parallelization without list
		if not isinstance(self.pardims, list):
			self.pardims = [self.pardims]


		# First, we setup each dimension
		# This also runs a validation checker over it
		dim_id = 0
		for i in self.pardims:
			i.setup(dim_id)
			dim_id += 1

		# Compute total number of resources over all dimensions
		self.num_cores_per_rank = _prod(i.num_cores_per_rank for i in self.pardims)

		# Check if number of cores per rank exceeds the available number of cores per node
		if self.num_cores_per_rank > platform_resources.num_cores_per_node:
			self.print()
			self.error("Invalid config for parallelization: self.num_cores_per_rank >= platform_resources.num_cores_per_node")

		# Number of ranks
		self.num_ranks = _prod(i.num_ranks for i in self.pardims)
		if self.num_ranks <= 0:
			self.error("self.num_ranks <= 0")

		# Check how many ranks we can run on each node
		self.num_ranks_per_node = int(math.ceil(platform_resources.num_cores_per_node // self.num_cores_per_rank))
		if self.num_ranks_per_node <= 0:
			self.error("self.num_ranks_per_node <= 0")

		# Reduce ranks per node if only a single node is used with all ranks on this particular node
		if self.num_ranks_per_node > self.num_ranks:
			self.num_ranks_per_node = self.num_ranks

		# Compute number of cores per node
		if self.num_cores_per_node == None:
			self.num_cores_per_node = self.num_cores_per_rank*self.num_ranks_per_node
			
		#
		# Compute raw numbers and compare to new number
		# The new number must be always \leq than the raw number
		# due to additional restrictions
		#
		# We do this mainly for debugging restrictions
		#
		# VALIDATION for inconsistencies
		raw_num_ranks = _prod(i.num_ranks for i in self.pardims)
		if self.num_ranks < raw_num_ranks:
			self.print()
			self.error("Internal error: self.num_ranks < raw_num_ranks")

		# Number of nodes
		self.num_nodes = int(math.ceil(self.num_ranks / self.num_ranks_per_node))
		if self.num_nodes <= 0:
			self.error("self.num_nodes <= 0")

		# VALIDATION for inconsistencies
		if self.num_nodes * self.num_ranks_per_node != self.num_ranks:
			self.print()
			self.error("Error: self.num_nodes * self.num_ranks_per_node != self.num_ranks\n******* Please change your job settings to avoid this *******")

		self.num_cores = self.num_nodes * platform_resources.num_cores_per_node

		#
		# VALIDATION for hardware restrictions
		#

		# Enough computing cores?
		if self.num_ranks*self.num_cores_per_rank > platform_resources.num_cores:
			self.print()
			self.error("Invalid config for parallelization: self.num_ranks*self.num_cores_per_rank >= platform_resources.num_cores")

		if self.num_cores > platform_resources.num_cores:
			self.print()
			self.error("Invalid config for parallelization: self.num_cores >= platform_resources.num_cores")


		#
		# Finally, setup variables without any restrictions
		#

		# Number of total (e.g. OpenMP) threads per rank (There are no restrictions for logical threading)
		self.num_threads_per_rank = _prod(i.num_threads_per_rank for i in self.pardims)



	def getUniqueID(self, i_filter):
		"""
		Return a unique ID including *all* string and number attributes of this class

		i_filter:
			list of filter names to filter out from unique ID generation
		"""
		retval = ''
		if not 'parallelization' in i_filter:
			retval += 'PAR'

			# cores per rank
			retval += "_cpr"+str(self.num_cores_per_rank)

			# mpi ranks
			retval += "_rks"+str(self.num_ranks)

			for i in self.pardims:
				retval += '_'+i.dim_name+str(i.num_cores).zfill(3)

		return retval


if __name__ == "__main__":
	p = SWEETParallelization()
	s = p.getUniqueID()
	p.info(s)
	p.print()

	p.info("FIN")
