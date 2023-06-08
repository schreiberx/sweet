import math
from mule.JobMule import *



"""
Parallelization parameters
"""

def setupParallelization(jg):

    pint = False

    # Number of threads to use in time and space dimension
    if 'REXI(' in jg.runtime.timestepping_method:

        pint = True

        # Number of ranks in time
        if jg.compile.sweet_mpi == "enable":
            # Use simply 2 ranks for testing things out
            num_ranks_time = 2
        else:
            num_ranks_time = 1

        divval = 1
        if jg.platform_resources.num_nodes == 1:
            divval = num_ranks_time

        # Number of cores for time parallelization
        if jg.compile.rexi_thread_parallel_sum == "enable":
            # Split between space and time
            num_threads_time = max(1, int(math.sqrt(jg.platform_resources.num_cores_per_socket//divval)))
            num_cores_time = num_threads_time
            num_threads_space = max(1, int(jg.platform_resources.num_cores_per_socket//divval//num_threads_time))
            threading_time = True

        else:
            # Parallelization is only in spatial dimension
            num_threads_time = 1
            num_threads_space = jg.platform_resources.num_cores_per_socket//divval
            threading_time = False

        
        # Number of ranks in space
        num_ranks_space = 1

    else:
        # Number of cores
        num_threads_time = 1
        num_threads_space = jg.platform_resources.num_cores_per_socket
        threading_time = False

        # Number of ranks
        num_ranks_time = 1
        num_ranks_space = 1

    assert num_threads_time*num_threads_space <= jg.platform_resources.num_cores_per_socket

    ptime = JobParallelizationDimOptions('time')
    ptime.num_cores_per_rank = num_threads_time
    ptime.num_threads_per_rank = num_threads_time
    ptime.num_ranks = num_ranks_time
    ptime.threading = jg.compile.rexi_thread_parallel_sum == "enable"

    pspace = JobParallelizationDimOptions('space')
    pspace.num_cores_per_rank = num_threads_space
    pspace.num_threads_per_rank = num_threads_space
    pspace.num_ranks = num_ranks_space
    pspace.threading = jg.compile.threading != "off"

    # Setup parallelization
    if pint:
        jg.setup_parallelization([ptime,pspace])
    else:
        jg.setup_parallelization([pspace])

    jg.runtime.sh_setup_num_threads = num_threads_space

