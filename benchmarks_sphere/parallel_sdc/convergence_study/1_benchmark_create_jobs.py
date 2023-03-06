#! /usr/bin/env python3
import sys
from itertools import product

from mule import JobGeneration, JobParallelizationDimOptions
from mule.sdc import getSDCSetup
jg = JobGeneration()

verbose = False
#verbose = True


##################################################
##################################################
##################################################

jg.compile.mode = 'release'

jg.parallelization.core_oversubscription = False
jg.parallelization.core_affinity = 'compact'

jg.compile.threading = 'omp'

gen_reference_solution = True

jg.runtime.max_simulation_time = 60*60*24*8    # 8 days

#params_timestep_sizes_explicit = [30]
params_timestep_sizes_explicit = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360]

#params_timestep_sizes_implicit = [30]
params_timestep_sizes_implicit = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360, 480, 600, 720, 960]

#params_timestep_sizes_rexi = [30]
params_timestep_sizes_exp = [15, 30, 60, 120, 180, 240, 300, 360, 480, 600, 720, 960]



# Parallelization
params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
params_ptime_num_cores_per_rank = [1]

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# No output
#jg.runtime.output_filename = "-"

##################################################


#
# Force deactivating Turbo mode
#
jg.parallelization.force_turbo_off = True


def estimateWallclockTime(jg):
    # 1 Hour max
    return 60*60

jg.compile.lapack = 'enable'
jg.compile.mkl = 'disable'

# Request dedicated compile script
jg.compilecommand_in_jobscript = False


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_sphere'

jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

jg.compile.benchmark_timings = 'enable'
jg.compile.quadmath = 'disable'


#
# Activate Fortran source
#
jg.compile.fortran_source = 'enable'


# Verbosity mode
jg.runtime.verbosity = 0

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 128
jg.runtime.space_res_physical = -1

#
# Benchmark
#
jg.runtime.benchmark_name = "galewsky"

#
# Compute error
#
jg.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 1

# Leave instability checks activated
# Don't activate them since they are pretty costly!!!
jg.runtime.instability_checks = 0
jg.runtime.viscosity = 0.0

#
# This allows including this file
#
sdcBaseParameters = {
    'nNodes': 3,
    'nodeType': 'RADAU-RIGHT',
    'nIter': 3,
    'qDeltaImplicit': 'BE',
    'qDeltaExplicit': 'FE',
    'qDeltaInitial': 'BEPAR',
    'initialSweepType': 'COPY',
    'diagonal': False,
    'useEndUpdate': False,
}

if __name__ == "__main__":

    sdcVariablePars = {'nIter': [1, 2, 3]}

    #
    # Reference solution
    #
    if gen_reference_solution:
        jg.runtime.timestep_size  = params_timestep_sizes_explicit[0]

        jg.runtime.timestepping_method = 'ln_erk'
        jg.runtime.timestepping_order = 4
        jg.runtime.timestepping_order2 = 4

        # Update TIME parallelization
        ptime = JobParallelizationDimOptions('time')
        ptime.num_cores_per_rank = 1
        ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
        ptime.num_ranks = 1

        pspace = JobParallelizationDimOptions('space')
        pspace.num_cores_per_rank = params_pspace_num_cores_per_rank[0]
        pspace.num_threads_per_rank = params_pspace_num_threads_per_rank[0]
        pspace.num_ranks = params_ptime_num_cores_per_rank[0]

        # Setup parallelization
        jg.setup_parallelization([pspace, ptime])

        if verbose:
            pspace.print()
            ptime.print()
            jg.parallelization.print()

        jg.reference_job = True
        jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

        jg.gen_jobscript_directory(f'job_benchref_RK4_dt{jg.runtime.timestep_size:06.2f}')
        jg.reference_job = False

        jg.reference_job_unique_id = jg.job_unique_id



    #
    # Create job scripts for SDC
    #
    for dt in params_timestep_sizes_implicit:
        for val in sdcVariablePars['nIter']:

            jg.runtime.timestepping_method = 'ln_imex_sdc'
            jg.runtime.timestepping_order = 1
            jg.runtime.timestepping_order2 = 1
            jg.runtime.timestep_size = dt
            
            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = params_pspace_num_cores_per_rank[0]
            pspace.num_threads_per_rank = params_pspace_num_threads_per_rank[0]
            pspace.num_ranks = params_ptime_num_cores_per_rank[0]
            pspace.setup()

            # Update TIME parallelization
            ptime = JobParallelizationDimOptions('time')
            ptime.num_cores_per_rank = 1
            ptime.num_threads_per_rank = 1
            ptime.num_ranks = 1
            ptime.setup()

            jg.setup_parallelization([pspace, ptime])

            if verbose:
                pspace.print()
                ptime.print()
                jg.parallelization.print()

            jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)
            jg.gen_jobscript_directory(f'job_bench_SDC_dt{dt:06.2f}_K{val}')

    # Write compile script
    jg.write_compilecommands("./compile_platform_"+jg.platforms.platform_id+".sh")


