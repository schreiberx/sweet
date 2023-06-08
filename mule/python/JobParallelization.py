import math
import sys
from functools import reduce
import operator

from mule.InfoError import *
from mule.JobPlatformResources import *
from mule.JobParallelizationDimOptions import *

__all__ = ['JobParallelization']

def _prod(iterable):
    return reduce(operator.mul, iterable, 1)


class JobParallelization(InfoError):

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
    def __init__(self, dummy_init = False):

        self.init_phase = True

        InfoError.__init__(self, "JobParallelization")

        self.reset(dummy_init)

        #
        # WARNING:
        # Leave these variables here to ensure being not influenced by reset()
        #

        #
        # Disable utilization of `mpiexec` to run job
        # This is required to run e.g. the validation scripts should be
        # (currently) executed on a single node and without MPI support
        self.mpiexec_disabled = False

        # Force disabling of turbo mode (if supported)
        self.force_turbo_off = False

        # Qualitative settings

        # Allow oversubscription (aka Hyperthreading)
        self.core_oversubscription = False

        # affinities:
        # compact, scatter
        self.core_affinity = None

        # max wallclock time, default: 1h
        self.max_wallclock_seconds = 60*60
        
        self.init_phase = False


    def __setattr__(self, name, value):

        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value


    def reset(self, dummy_init = False):
        """
        Reset functionality for a fresh configuration in case that a new setup is triggered
        """

        # Number of cores per rank
        self.num_cores_per_rank : int = None 

        # Number of threads per rank
        self.num_threads_per_rank : int = None

        # Number of ranks per node
        self.num_ranks_per_node : int = None

        # Number of cores per node
        self.num_cores_per_node : int = None 

        # Number of total ranks
        self.num_ranks : int = None

        # Number of total nodes
        self.num_nodes : int = None

        # Number of total cores
        self.num_cores : int = None


        # String for OMP_NUM_THREADS
        self.omp_num_threads : str = None


        # List with parallelization information in each dimension
        # Note, that space dimension can and should be treated as a single dimension
        self.pardims = None

        self.pardims_dict = {}

    def get_max_wallclock_seconds_hh_mm_ss(self):
        """
        Return properly formatted self.max_wallclock_seconds usable for job scripts
        """
        secs = self.max_wallclock_seconds
        # seconds
        s = int(secs)
        m = s // 60
        s = s % 60
        h = m // 60
        m = m % 60

        stest = h*60*60 + m*60 + s
        if int(secs) != stest:
            print(secs)
            print(stest)
            raise Exception("Internal error!")

        return str(h).zfill(2)+":"+str(m).zfill(2)+":"+str(s).zfill(2)




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
        self.info("mpiexec_disabled: "+str(self.mpiexec_disabled))
        self.info("force_turbo_off: "+str(self.force_turbo_off))



    def dummy_setup_if_no_setup(self, platform_resources : JobPlatformResources):
        """
        Setup a dummy parallelization dimension to use one rank on one node and all cores on the node
        """

        if self.pardims == None:
            dummy = JobParallelizationDimOptions("dummy")
            dummy.num_cores = platform_resources.num_cores_per_node
            dummy.num_cores_per_rank = dummy.num_cores
            dummy.num_threads_per_rank = dummy.num_cores
            dummy.num_ranks = 1
            self.setup([dummy], platform_resources)
            self.print()


    def setup(self, list_pardims, platform_resources : JobPlatformResources, override_insufficient_resources : bool = False):
        """
        Setup data which is required by the platform specific scripts to
        generate the job scripts

        Parameters
        ----------
        list_pardims:    JobParallelizationDimOptions
            List with options for parallelization in each dimension

        platform_resources:
            reference to jobgeneration class

        override_insufficient_resources:
            If True, insufficient resources will not lead to abortion of this test case.
            This is useful for CI systems which have less cores than required.

        #mode : string
        #    'serial': No parallelization
        """

        self.reset()

        self.pardims = list_pardims

#        # Support space-only parallelization without list
        if not isinstance(self.pardims, list):
            self.pardims = [self.pardims]

        # First, we setup each dimension
        # This also runs a validation checker over it
        dim_id = 0
        self.pardims_dict = {}
        for i in self.pardims:
            i.setup(dim_id)
            dim_id += 1

            self.pardims_dict[i.dim_name] = i

        # Compute total number of resources over all dimensions
        self.num_cores_per_rank = _prod(i.num_cores_per_rank for i in self.pardims)

        # Check if number of cores per rank exceeds the available number of cores per node
        if self.num_cores_per_rank > platform_resources.num_cores_per_node:
            self.print()
            self.error("Invalid config for parallelization: self.num_cores_per_rank > platform_resources.num_cores_per_node")

        # Number of total MPI ranks
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

        # Enough computing nodes?
        if self.num_nodes > platform_resources.num_nodes:
            if not override_insufficient_resources:
                self.print()
                self.error("Invalid config for parallelization: self.num_nodes > platform_resources.num_nodes")
            else:
                self.num_nodes = platform_resources.num_nodes


        # VALIDATION for inconsistencies
        if self.num_nodes * self.num_ranks_per_node != self.num_ranks:
            if not override_insufficient_resources:
                self.print()
                self.error("Error: self.num_nodes * self.num_ranks_per_node != self.num_ranks\n******* Please change your job settings to avoid this *******")
            else:
                if self.num_nodes != 1:
                    self.error("Error: override_insufficient_resources only works for one node")
                self.num_ranks_per_node = self.num_ranks

        self.num_cores = self.num_nodes * platform_resources.num_cores_per_node

        #
        # VALIDATION for hardware restrictions
        #

        # Enough computing cores?
        if self.num_ranks*self.num_cores_per_rank > platform_resources.num_cores:
            if not override_insufficient_resources:
                self.print()
                self.error("Invalid config for parallelization: self.num_ranks*self.num_cores_per_rank > platform_resources.num_cores")
            else:
                print("Ignoring insufficient cores due to override_insufficient_resources=True")

        if self.num_cores > platform_resources.num_cores:
            if not override_insufficient_resources:
                self.print()
                self.error("Invalid config for parallelization: self.num_cores > platform_resources.num_cores")
            else:
                print("Ignoring insufficient cores due to override_insufficient_resources=True")


        #
        # Finally, setup variables without any restrictions
        #
        self.num_threads_per_rank = _prod(i.num_threads_per_rank for i in self.pardims)

        # Prepare OpenMP nested OMP_NUM_THREADS environment variable
        lThreads = []
        for p in self.pardims:
            if p.threading:
                lThreads.append(p.num_threads_per_rank)

        if len(lThreads) > 0:
            self.omp_num_threads = ','.join([str(n) for n in lThreads])



    def getUniqueID(self, i_filters):
        """
        Return a unique ID including *all* string and number attributes of this class

        i_filter:
            list of filter names to filter out from unique ID generation
        """
        retval = ''
        if not 'parallelization' in i_filters:

            if not 'parallelization.mpi_ranks' in i_filters:
                # mpi ranks
                retval += "_r"+str(self.num_ranks).zfill(5)

            if not 'parallelization.cores_per_rank' in i_filters:
                # cores per rank
                retval += "_cpr"+str(self.num_cores_per_rank).zfill(3)

            if not 'parallelization.threads_per_rank' in i_filters:
                # threads per rank
                retval += "_tpr"+str(self.num_threads_per_rank).zfill(3)

            if not 'parallelization.dims' in i_filters:
                retval += "_DIMS"
                for i in self.pardims:
                    retval += '_'+i.dim_name+str(i.num_cores).zfill(3)

            if retval != '':
                retval = 'PAR'+retval

        return retval


if __name__ == "__main__":
    p = JobParallelization()
    s = p.getUniqueID()
    p.info(s)
    p.print()

    p.info("FIN")
