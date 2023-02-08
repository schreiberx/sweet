#! /usr/bin/env python3

import os
import sys
import math
import numpy as np

from itertools import product

from mule.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

p = JobGeneration()

verbose = False
#verbose = True

##################################################

##################################################

p.compile.mode = 'release'
#p.compile.sweet_mpi = 'disable'


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

#params_timestep_sizes_explicit_ = [15*(2**i) for i in range(0, 4)]
#params_timestep_sizes_explicit_ = [60]

base_timestep_size = 128/p.runtime.space_res_spectral*300.0
params_timestep_sizes_explicit_ = [base_timestep_size]

params_timestep_sizes_implicit_ = [15*(2**i) for i in range(2, 6)]
params_timestep_sizes_sl_ = [15*(2**i) for i in range(2, 6)]


# Parallelization
params_pspace_num_cores_per_rank = [p.platform_resources.num_cores_per_socket]
params_pspace_num_threads_per_rank = [p.platform_resources.num_cores_per_socket]

unique_id_filter = []
unique_id_filter.append('compile')
#unique_id_filter.append('runtime.galewsky_params')
unique_id_filter.append('runtime.rexi')
unique_id_filter.append('runtime.benchmark')
unique_id_filter.append('runtime.max_simulation_time')

p.unique_id_filter = unique_id_filter

#p.runtime.output_timestep_size = p.runtime.max_simulation_time


##########################################################################
##########################################################################
##########################################################################

def estimateWallclockTime(p):
    return 12*60*60


p.compile.lapack = 'enable'
p.compile.mkl = 'disable'

p.compilecommand_in_jobscript = False


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_sphere'

p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.compile.benchmark_timings = 'enable'

p.compile.quadmath = 'disable'


#
# Activate Fortran source
#
#p.compile.fortran_source = 'enable'


# Verbosity mode
p.runtime.verbosity = 0


#
# Compute error
#
p.runtime.compute_error = 0


#
# Preallocate the REXI matrices
#
#p.runtime.rexi_sphere_preallocation = 1

# Leave instability checks activated
p.runtime.instability_checks = 1
# Don't activate them for wallclock time studies since they are pretty costly!!!
#p.runtime.instability_checks = 0

p.runtime.viscosity = 0.0



#
# allow including this file
#
if __name__ == "__main__":

    ts_methods = [
        # REFERENCE METHOD, not computed by default
        ['ln_erk',        4,    4,    0],


        ###########
        # Runge-Kutta
        ###########
        ['ln_erk',        4,    4,    0],
        #['ln_erk',        2,    2,    0],

        ###########
        # CN
        ###########
        #['lg_irk_lc_n_erk_ver0',    2,    2,    0],
        #['lg_irk_lc_n_erk_ver1',    2,    2,    0],

        #['l_irk_na_sl_nd_settls_ver1',    2,    2,    0],
        #['l_irk_na_sl_nd_settls_ver2',    2,    2,    0],

        #['lg_irk_na_sl_lc_nd_settls_ver1',    2,    2,    0],
        #['lg_irk_na_sl_lc_nd_settls_ver2',    2,    2,    0],
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

        if len(tsm) > 4:
            s = tsm[4]
            p.runtime.load_from_dict(tsm[4])

        tsm_name = tsm[0]
        if 'ln_erk' in tsm_name:
            params_timestep_sizes = params_timestep_sizes_explicit_
        elif 'l_erk' in tsm_name or 'lg_erk' in tsm_name:
            params_timestep_sizes = params_timestep_sizes_explicit_
        elif 'l_irk' in tsm_name or 'lg_irk' in tsm_name:
            params_timestep_sizes = params_timestep_sizes_implicit_
        elif '_sl' in tsm_name:
            params_timestep_sizes = params_timestep_sizes_sl_
        else:
            print("Unable to identify time stepping method "+tsm_name)
            sys.exit(1)

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

            p.gen_jobscript_directory('job_bench_'+p.getUniqueID())

    p.write_compilecommands()
