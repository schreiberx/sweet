#! /usr/bin/env python3

import sys
from itertools import product

from mule.JobGeneration import JobGeneration
from mule.JobParallelizationDimOptions import JobParallelizationDimOptions
from mule.sdc import getSDCSetup

p = JobGeneration()
verbose = True

p.runtime.paramsSDC = getSDCSetup(
    nNodes=3,
    nIter=3,
    nodeType='RADAU-RIGHT', 
    qDeltaImplicit='OPT-QMQD-0', 
    qDeltaExplicit='PIC', 
    qDeltaInitial='BEPAR',
    diagonal=True,
    initialSweepType="QDELTA",
    useEndUpdate=False
)

p.compile.mode = 'release'
p.compile.gui = 'enable'
p.compile.mode = 'debug'

#
# Mode and Physical resolution
#
p.runtime.space_res_spectral = 128
p.runtime.space_res_physical = -1

p.parallelization.core_oversubscription = False
p.parallelization.core_affinity = 'compact'

p.compile.threading = 'omp'
p.compile.rexi_thread_parallel_sum = 'disable'

gen_reference_solution = False
p.runtime.benchmark_name = "galewsky"

p.runtime.max_simulation_time = 60*60*24*8    # 8 days

p.runtime.output_timestep_size = 60*60  # Generate output every 1 hour
p.runtime.output_file_mode = 'bin'

params_timestep_size_reference = 30.0
base_timestep_size = 128/p.runtime.space_res_spectral*1200.0

# Parallelization
nSpacePar = int(sys.argv[1]) if len(sys.argv) > 1 else p.platform_resources.num_cores_per_socket
params_pspace_num_cores_per_rank = [nSpacePar]
params_pspace_num_threads_per_rank = [nSpacePar]

unique_id_filter = []
unique_id_filter.append('compile')
unique_id_filter.append('runtime.max_simulation_time')

p.unique_id_filter = unique_id_filter

def estimateWallclockTime(p):
    return 12*60*60


p.compile.lapack = 'enable'
p.compile.mkl = 'disable'
p.compilecommand_in_jobscript = False


#
# Run simulation on plane or sphere
#
p.compile.program = 'programs/pde_sweSphere'

p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.compile.benchmark_timings = 'enable'
p.compile.quadmath = 'disable'

p.runtime.verbosity = 0

# Leave instability checks activated
p.runtime.instability_checks = 0
# Don't activate them for wallclock time studies since they are pretty costly!!!
#p.runtime.instability_checks = 0

p.runtime.viscosity = 0.0

#
# allow including this file
#
if __name__ == "__main__":

    ts_methods = [
        # REFERENCE METHOD, not computed by default
        ['ln_erk',        4,    4],

        ###########
        # IMEX Euler
        ###########
        # ['l_irk',        1,    1],

        ###########
        # IMEX SDC
        ###########
        ['ln_imex_sdc',        1,    1],
        # ['ln_erk',        4,    4],
    ]


    #
    # Reference solution
    #
    p.reference_job_unique_id = None

    if gen_reference_solution:
        tsm = ts_methods[0]

        p.runtime.timestep_size  = params_timestep_size_reference

        p.runtime.timestepping_method = tsm[0]
        p.runtime.timestepping_order = tsm[1]
        p.runtime.timestepping_order2 = tsm[2]

        pspace = JobParallelizationDimOptions('space')
        pspace.num_cores_per_rank = 1
        pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
        pspace.num_ranks = 1

        # Setup parallelization
        p.setup_parallelization([pspace])

        if verbose:
            pspace.print()
            p.parallelization.print()

        p.parallelization.max_wallclock_seconds = estimateWallclockTime(p)

        p.gen_jobscript_directory('job_benchref_'+p.getUniqueID())

        # Use this as a reference job
        p.reference_job_unique_id = p.job_unique_id


    for tsm in ts_methods[1:]:
        p.runtime.timestepping_method = tsm[0]
        p.runtime.timestepping_order = tsm[1]
        p.runtime.timestepping_order2 = tsm[2]

        tsm_name = tsm[0]
        params_timestep_sizes = [base_timestep_size]

        for (
                pspace_num_cores_per_rank,
                pspace_num_threads_per_rank,
                p.runtime.timestep_size
            ) in product(
                params_pspace_num_cores_per_rank,
                params_pspace_num_threads_per_rank,
                params_timestep_sizes
            ):
            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = pspace_num_cores_per_rank
            pspace.num_threads_per_rank = pspace_num_threads_per_rank
            pspace.num_ranks = 1
            pspace.setup()

            p.setup_parallelization([pspace])

            if verbose:
                pspace.print()
                p.parallelization.print()

            p.parallelization.max_wallclock_seconds = estimateWallclockTime(p)

            p.gen_jobscript_directory(f'job_bench_{p.getUniqueID()}')

    p.write_compilecommands()
