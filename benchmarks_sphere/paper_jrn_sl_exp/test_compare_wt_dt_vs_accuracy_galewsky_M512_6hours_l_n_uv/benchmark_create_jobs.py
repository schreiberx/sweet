#! /usr/bin/env python3

import os
import sys
import math

from itertools import product

# REXI
from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *

efloat_mode = "float"
#efloat_mode = "mpfloat"


from mule.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *
jg = JobGeneration()

verbose = False
#verbose = True

##################################################

##################################################

jg.compile.mode = 'release'
if '_gcc' in os.getenv('MULE_PLATFORM_ID'):
    jg.compile.compiler = 'gcc'
else:
    jg.compile.compiler = 'intel'
jg.compile.sweet_mpi = 'enable'


jg.parallelization.core_oversubscription = False
jg.parallelization.core_affinity = 'compact'

jg.compile.threading = 'omp'
jg.compile.rexi_thread_parallel_sum = 'disable'

gen_reference_solution = True

jg.runtime.max_simulation_time = 60*60*6
#jg.runtime.max_simulation_time = 30*16

jg.runtime.max_wallclock_time = 30*60       # 30 minutes max

#space_res_spectral_ = [64, 128, 256]
space_res_spectral_ = [512]


# Reference time step size
timestep_size_reference = 5


params_timestep_sizes_explicit = [15/2*(2**i) for i in range(10)]
params_timestep_sizes_implicit = [15/2*(2**i) for i in range(10)]
params_timestep_sizes_exp = [15/2*(2**i) for i in range(9)]


# Parallelization
params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
params_ptime_num_cores_per_rank = [1]

unique_id_filter = []
#unique_id_filter.append('simparams')
unique_id_filter.append('compile')
unique_id_filter.append('runtime.disc_space')
unique_id_filter.append('runtime.timestep_order')
#unique_id_filter.append('timestep_size')
unique_id_filter.append('runtime.rexi')
unique_id_filter.append('runtime.benchmark')

unique_id_filter.append('parallelization')

jg.unique_id_filter = unique_id_filter


jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# No output
#jg.runtime.output_filename = "-"


#
# Force deactivating Turbo mode
#
jg.parallelization.force_turbo_off = True




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
# Benchmark
#
jg.runtime.benchmark_name = "galewsky"

#
# Binary output
#
jg.runtime.output_file_mode = "bin"

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


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
#
jg.runtime.rexi_method = ''

jg.runtime.viscosity = 0.0

jg.runtime.rexi_method = 'direct'



def estimateWallclockTime(jg):

    if jg.reference_job:
        return 2*24*60*60

    return 1*60*60




#
# allow including this file
#
if __name__ == "__main__":

    ts_methods = [
        ['ln_erk_split_uv',        4,    4,    0],

        ###########
        # Runge-Kutta
        ###########
        ['ln_erk_split_aa_uv',     2,    2,    0],
        ['ln_erk_split_uv',        2,    2,    0],


        ###########
        # SETTLS variants
        ###########
        ['l_irk_na_sl_nr_settls_uv_only',    2,    2,    0],

        ['l_irk_na_sl_nr_settls_ver0_uv',    2,    2,    0],
        ['l_irk_na_sl_nr_settls_ver1_uv',    2,    2,    0],

        ['lg_irk_na_sl_lc_nr_settls_ver0_uv',    2,    2,    0],
        ['lg_irk_na_sl_lc_nr_settls_ver1_uv',    2,    2,    0],

        ###########
        # EXP variants
        ###########
        ['lg_exp_na_sl_lc_nr_settls_ver0_uv',    2,    2,    0],
        ['lg_exp_na_sl_lc_nr_settls_ver1_uv',    2,    2,    0],
    ]

    for space_res_spectral in space_res_spectral_:
        jg.runtime.space_res_spectral = space_res_spectral
        #jg.runtime.reuse_plans = 2    # enforce using plans (todo, enforcing not yet implemented)!

        #
        # Reference solution
        #
        if gen_reference_solution:
            tsm = ts_methods[0]
            jg.runtime.timestep_size  = timestep_size_reference

            jg.runtime.timestepping_method = tsm[0]
            jg.runtime.timestepping_order = tsm[1]
            jg.runtime.timestepping_order2 = tsm[2]

            # Update TIME parallelization
            ptime = JobParallelizationDimOptions('time')
            ptime.num_cores_per_rank = 1
            ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
            ptime.num_ranks = 1

            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = 1
            pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
            pspace.num_ranks = 1

            # Setup parallelization
            jg.setup_parallelization([pspace, ptime])

            if verbose:
                pspace.print()
                ptime.print()
                jg.parallelization.print()

            if len(tsm) > 4:
                s = tsm[4]
                jg.load_from_dict(tsm[4])

            jg.reference_job = True
            jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

            _a = jg.runtime.max_wallclock_time
            jg.runtime.max_wallclock_time = 2*24*60*60       # 30 minutes max

            jg.reference_job_unique_id = None
            jg.gen_jobscript_directory('job_benchref_'+jg.getUniqueID())
            jg.reference_job = False

            jg.runtime.max_wallclock_time = _a

            jg.reference_job_unique_id = jg.job_unique_id



        #
        # Create job scripts
        #
        for tsm in ts_methods[1:]:

            jg.runtime.timestepping_method = tsm[0]
            jg.runtime.timestepping_order = tsm[1]
            jg.runtime.timestepping_order2 = tsm[2]

            if len(tsm) > 4:
                s = tsm[4]
                jg.runtime.load_from_dict(tsm[4])

            exp_integrator = False

            tsm_name = tsm[0]
            if 'l_erk' in tsm_name or 'lg_erk' in tsm_name or 'ln_erk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_explicit
            elif 'l_na_erk' in tsm_name or 'ln_erk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_explicit
            elif '_irk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_implicit
            elif '_rexi' in tsm_name or '_exp' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_exp
                exp_integrator = True
            else:
                raise Exception("Unable to identify time stepping method "+tsm_name)


            for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, params_timestep_sizes):
                pspace = JobParallelizationDimOptions('space')
                pspace.num_cores_per_rank = pspace_num_cores_per_rank
                pspace.num_threads_per_rank = pspace_num_threads_per_rank
                pspace.num_ranks = 1
                pspace.setup()

                #if not exp_integrator or 'lg_' in tsm_name:
                if True:

                    # Always use direct REXI method if no parallel-in-time
                    jg.runtime.rexi_method = 'direct'

                    # Update TIME parallelization
                    ptime = JobParallelizationDimOptions('time')
                    ptime.num_cores_per_rank = 1
                    ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
                    ptime.num_ranks = 1
                    ptime.setup()

                    jg.setup_parallelization([pspace, ptime])

                    if verbose:
                        pspace.print()
                        ptime.print()
                        jg.parallelization.print()

                    jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

                    jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

                else:

                    raise Exception("This branch shouldn't be taken, yet")
    #
    # SHTNS plan generation scripts
    #
    #jg.runtime.reuse_plans = 1    # search for awesome plans and store them

    #
    # Create dummy scripts to be used for SHTNS script generation
    #

    # No parallelization in time
    ptime = JobParallelizationDimOptions('time')
    ptime.num_cores_per_rank = 1
    ptime.num_threads_per_rank = 1
    ptime.num_ranks = 1
    ptime.setup()

    for tsm in ts_methods[1:2]:

        jg.runtime.timestepping_method = tsm[0]
        jg.runtime.timestepping_order = tsm[1]
        jg.runtime.timestepping_order2 = tsm[2]

        if not '_rexi' in jg.runtime.timestepping_method:
            jg.runtime.rexi_method = ''
        else:
            jg.runtime.rexi_method = 'ci'

        if len(tsm) > 4:
            s = tsm[4]
            jg.runtime.load_from_dict(tsm[4])


        for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, [params_timestep_sizes_explicit[0]]):
            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = pspace_num_cores_per_rank
            pspace.num_threads_per_rank = pspace_num_threads_per_rank
            pspace.num_ranks = 1
            pspace.setup()

            jg.setup_parallelization([pspace, ptime])

            # Use 10 minutes per default to generate plans
            jg.parallelization.max_wallclock_seconds = 60*10

            # Set simtime to 0
            #jg.runtime.max_simulation_time = 0

            # No output
            jg.runtime.output_timestep_size = -1
            jg.runtime.output_filename = "-"

            jobdir = 'job_plan_'+jg.getUniqueID()
            jg.gen_jobscript_directory(jobdir)



    # Write compile script
    jg.write_compilecommands("./compile_platform_"+jg.platforms.platform_id+".sh")


    print("")
    print("Timestepping methods:")
    m = [i[0] for i in ts_methods[1:]]
    m.sort()

    for i in m:
        print(i)
    print("")


